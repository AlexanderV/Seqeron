// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S107, S1121, S3776, S1117, S125, S2325, S1066, S3358, S4136, S2692, IDE0059
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The raw <c>deflate</c> of zlib 1.3 (<c>deflate.c</c> <c>deflate_slow</c>, <c>trees.c</c>) for level 6 (<c>Z_DEFAULT_COMPRESSION</c>), <c>windowBits</c> -15, <c>memLevel</c> 8, <c>Z_DEFAULT_STRATEGY</c>, one call with <c>Z_FINISH</c>:
/// the BAM writer of EMBOSS 6.6.0 (<c>bamDeflateBlock</c> of <c>ajseqbam.c</c>) compresses its BGZF blocks so, and the bytes of the file are the bytes of the installed library (A5-8c3h, B09 WP59). The library of the .NET
/// runtime (zlib-ng) makes other, equally valid, streams; the ones of this class are byte for byte those of zlib.
/// </summary>
public sealed class EmbossZlibDeflate
{
    private const int WBits = 15, WSize = 1 << WBits, WMask = WSize - 1, MemLevel = 8, HashBits = MemLevel + 7, HashSize = 1 << HashBits, HashMask = HashSize - 1, HashShift = (HashBits + 3 - 1) / 3;
    private const int MinMatch = 3, MaxMatch = 258, MinLookahead = MaxMatch + MinMatch + 1, MaxDist = WSize - MinLookahead, TooFar = 4096;
    private const int LitBufSize = 1 << (MemLevel + 6), SymEnd = (LitBufSize - 1) * 3;
    private const int GoodLength = 8, MaxLazy = 16, NiceLength = 128, MaxChain = 128;
    private const int LengthCodes = 29, Literals = 256, LCodes = Literals + 1 + LengthCodes, DCodes = 30, BlCodes = 19, HeapSize = 2 * LCodes + 1, MaxBits = 15, MaxBlBits = 7, EndBlock = 256;
    private const int Rep3_6 = 16, Repz3_10 = 17, Repz11_138 = 18;

    private static readonly int[] ExtraLBits = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
    private static readonly int[] ExtraDBits = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
    private static readonly int[] ExtraBlBits = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 7];
    private static readonly byte[] BlOrder = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

    private static readonly ushort[] StaticLCode = new ushort[LCodes + 2], StaticDCode = new ushort[DCodes];
    private static readonly byte[] StaticLLen = new byte[LCodes + 2], StaticDLen = new byte[DCodes];
    private static readonly byte[] DistCode = new byte[512], LengthCode = new byte[MaxMatch - MinMatch + 1];
    private static readonly int[] BaseLength = new int[LengthCodes], BaseDist = new int[DCodes];

    static EmbossZlibDeflate()
    {
        // tr_static_init
        int length = 0, code;
        for (code = 0; code < LengthCodes - 1; code++)
        {
            BaseLength[code] = length;
            for (int n = 0; n < 1 << ExtraLBits[code]; n++)
                LengthCode[length++] = (byte)code;
        }

        LengthCode[length - 1] = (byte)code;
        int dist = 0;
        for (code = 0; code < 16; code++)
        {
            BaseDist[code] = dist;
            for (int n = 0; n < 1 << ExtraDBits[code]; n++)
                DistCode[dist++] = (byte)code;
        }

        dist >>= 7;
        for (; code < DCodes; code++)
        {
            BaseDist[code] = dist << 7;
            for (int n = 0; n < 1 << (ExtraDBits[code] - 7); n++)
                DistCode[256 + dist++] = (byte)code;
        }

        var blCount = new ushort[MaxBits + 1];
        int nn = 0;
        while (nn <= 143) { StaticLLen[nn++] = 8; blCount[8]++; }
        while (nn <= 255) { StaticLLen[nn++] = 9; blCount[9]++; }
        while (nn <= 279) { StaticLLen[nn++] = 7; blCount[7]++; }
        while (nn <= 287) { StaticLLen[nn++] = 8; blCount[8]++; }
        GenCodes(StaticLCode, StaticLLen, LCodes + 1, blCount);
        for (nn = 0; nn < DCodes; nn++)
        {
            StaticDLen[nn] = 5;
            StaticDCode[nn] = (ushort)BiReverse(nn, 5);
        }
    }

    // a tree: Freq / Code share a field, Dad / Len share a field in zlib (ct_data); separate arrays here, the sharing never matters for the result
    private sealed class Tree(int size)
    {
        public readonly ushort[] Freq = new ushort[size], Code = new ushort[size], Dad = new ushort[size];
        public readonly ushort[] Len = new ushort[size];
    }

    private sealed class Desc(Tree dyn, int[] extraBits, int extraBase, int elems, int maxLength, byte[]? staticLen)
    {
        public readonly Tree Dyn = dyn;
        public int MaxCode;
        public readonly int[] ExtraBits = extraBits;
        public readonly int ExtraBase = extraBase, Elems = elems, MaxLength = maxLength;
        public readonly byte[]? StaticLen = staticLen;
    }

    private readonly byte[] _input;
    private int _inPos;
    private readonly byte[] _window = new byte[2 * WSize];
    private readonly ushort[] _prev = new ushort[WSize], _head = new ushort[HashSize];
    private int _insH, _blockStart, _matchLength, _matchAvailable, _strstart, _matchStart, _lookahead, _prevLength, _insert;
    private readonly byte[] _symBuf = new byte[LitBufSize * 3 + 8];
    private int _symNext;
    private readonly Tree _dynL = new(HeapSize), _dynD = new(2 * DCodes + 1), _blTree = new(2 * BlCodes + 1);
    private readonly Desc _lDesc, _dDesc, _blDesc;
    private readonly ushort[] _blCount = new ushort[MaxBits + 1];
    private readonly int[] _heap = new int[2 * LCodes + 1];
    private int _heapLen, _heapMax;
    private readonly byte[] _depth = new byte[2 * LCodes + 1];
    private long _optLen, _staticLen;
    private readonly List<byte> _out = [];
    private uint _biBuf;
    private int _biValid;

    private EmbossZlibDeflate(byte[] input)
    {
        _input = input;
        _lDesc = new Desc(_dynL, ExtraLBits, Literals + 1, LCodes, MaxBits, StaticLLen);
        _dDesc = new Desc(_dynD, ExtraDBits, 0, DCodes, MaxBits, StaticDLen);
        _blDesc = new Desc(_blTree, ExtraBlBits, 0, BlCodes, MaxBlBits, null);
        _matchLength = _prevLength = MinMatch - 1;
    }

    /// <summary>Compresses <paramref name="input"/> as one raw deflate stream (<c>deflateInit2(level 6, -15, 8)</c>, <c>deflate(Z_FINISH)</c>).</summary>
    internal static byte[] Compress(ReadOnlySpan<byte> input)
    {
        var d = new EmbossZlibDeflate(input.ToArray());
        d.InitBlock();
        d.DeflateSlow();
        return [.. d._out];
    }

    private void InitBlock()
    {
        Array.Clear(_dynL.Freq, 0, LCodes);
        Array.Clear(_dynD.Freq, 0, DCodes);
        Array.Clear(_blTree.Freq, 0, BlCodes);
        _dynL.Freq[EndBlock] = 1;
        _optLen = _staticLen = 0;
        _symNext = 0;
    }

    // ---- deflate.c

    private int ReadBuf(int dest, int size)
    {
        int len = Math.Min(_input.Length - _inPos, size);
        if (len == 0)
            return 0;
        Array.Copy(_input, _inPos, _window, dest, len);
        _inPos += len;
        return len;
    }

    private void FillWindow()
    {
        do
        {
            int more = 2 * WSize - _lookahead - _strstart;
            if (_strstart >= WSize + MaxDist)
            {
                Array.Copy(_window, WSize, _window, 0, WSize - more);
                _matchStart -= WSize;
                _strstart -= WSize;
                _blockStart -= WSize;
                if (_insert > _strstart)
                    _insert = _strstart;
                SlideHash();
                more += WSize;
            }

            if (_inPos >= _input.Length)
                break;
            int n = ReadBuf(_strstart + _lookahead, more);
            _lookahead += n;
            if (_lookahead + _insert >= MinMatch)
            {
                int str = _strstart - _insert;
                _insH = _window[str];
                _insH = ((_insH << HashShift) ^ _window[str + 1]) & HashMask;
                while (_insert > 0)
                {
                    _insH = ((_insH << HashShift) ^ _window[str + MinMatch - 1]) & HashMask;
                    _prev[str & WMask] = _head[_insH];
                    _head[_insH] = (ushort)str;
                    str++;
                    _insert--;
                    if (_lookahead + _insert < MinMatch)
                        break;
                }
            }
        }
        while (_lookahead < MinLookahead && _inPos < _input.Length);
    }

    private void SlideHash()
    {
        for (int n = 0; n < HashSize; n++)
        {
            int m = _head[n];
            _head[n] = (ushort)(m >= WSize ? m - WSize : 0);
        }

        for (int n = 0; n < WSize; n++)
        {
            int m = _prev[n];
            _prev[n] = (ushort)(m >= WSize ? m - WSize : 0);
        }
    }

    private int InsertString(int str)
    {
        _insH = ((_insH << HashShift) ^ _window[str + MinMatch - 1]) & HashMask;
        int head = _head[_insH];
        _prev[str & WMask] = (ushort)head;
        _head[_insH] = (ushort)str;
        return head;
    }

    private int LongestMatch(int curMatch)
    {
        int chainLength = MaxChain;
        int scan = _strstart;
        int bestLen = _prevLength;
        int niceMatch = NiceLength;
        int limit = _strstart > MaxDist ? _strstart - MaxDist : 0;
        int strend = _strstart + MaxMatch;
        byte scanEnd1 = _window[scan + bestLen - 1];
        byte scanEnd = _window[scan + bestLen];
        if (_prevLength >= GoodLength)
            chainLength >>= 2;
        if (niceMatch > _lookahead)
            niceMatch = _lookahead;
        do
        {
            int match = curMatch;
            if (_window[match + bestLen] != scanEnd || _window[match + bestLen - 1] != scanEnd1 || _window[match] != _window[scan] || _window[match + 1] != _window[scan + 1])
                continue;
            int s = scan + 2;
            int m = match + 2;
            // scan[2] == match[2] (the hash keys are equal): the comparison goes on from index 3
            do
            {
                s++;
                m++;
            }
            while (s < strend && _window[s] == _window[m]);
            int len = MaxMatch - (strend - s);
            if (len > bestLen)
            {
                _matchStart = curMatch;
                bestLen = len;
                if (len >= niceMatch)
                    break;
                scanEnd1 = _window[scan + bestLen - 1];
                scanEnd = _window[scan + bestLen];
            }
        }
        while ((curMatch = _prev[curMatch & WMask]) > limit && --chainLength != 0);
        return bestLen <= _lookahead ? bestLen : _lookahead;
    }

    private void DeflateSlow()
    {
        while (true)
        {
            if (_lookahead < MinLookahead)
            {
                FillWindow();
                if (_lookahead == 0)
                    break;
            }

            int hashHead = 0;
            if (_lookahead >= MinMatch)
                hashHead = InsertString(_strstart);
            _prevLength = _matchLength;
            int prevMatch = _matchStart;
            _matchLength = MinMatch - 1;
            if (hashHead != 0 && _prevLength < MaxLazy && _strstart - hashHead <= MaxDist)
            {
                _matchLength = LongestMatch(hashHead);
                if (_matchLength <= 5 && _matchLength == MinMatch && _strstart - _matchStart > TooFar)
                    _matchLength = MinMatch - 1;
            }

            if (_prevLength >= MinMatch && _matchLength <= _prevLength)
            {
                int maxInsert = _strstart + _lookahead - MinMatch;
                bool bflush = TallyDist(_strstart - 1 - prevMatch, _prevLength - MinMatch);
                _lookahead -= _prevLength - 1;
                _prevLength -= 2;
                do
                {
                    if (++_strstart <= maxInsert)
                        InsertString(_strstart);
                }
                while (--_prevLength != 0);
                _matchAvailable = 0;
                _matchLength = MinMatch - 1;
                _strstart++;
                if (bflush)
                    FlushBlock(false);
            }
            else if (_matchAvailable != 0)
            {
                bool bflush = TallyLit(_window[_strstart - 1]);
                if (bflush)
                    FlushBlock(false);
                _strstart++;
                _lookahead--;
            }
            else
            {
                _matchAvailable = 1;
                _strstart++;
                _lookahead--;
            }
        }

        if (_matchAvailable != 0)
        {
            TallyLit(_window[_strstart - 1]);
            _matchAvailable = 0;
        }

        _insert = _strstart < MinMatch - 1 ? _strstart : MinMatch - 1;
        FlushBlock(true);
    }

    private void FlushBlock(bool last)
    {
        TrFlushBlock(_blockStart >= 0 ? _blockStart : -1, _strstart - _blockStart, last);
        _blockStart = _strstart;
    }

    // ---- trees.c

    private bool TallyLit(int c)
    {
        _symBuf[_symNext++] = 0;
        _symBuf[_symNext++] = 0;
        _symBuf[_symNext++] = (byte)c;
        _dynL.Freq[c]++;
        return _symNext == SymEnd;
    }

    private static int DCode(int dist) => dist < 256 ? DistCode[dist] : DistCode[256 + (dist >> 7)];

    private bool TallyDist(int dist, int lc)
    {
        _symBuf[_symNext++] = (byte)dist;
        _symBuf[_symNext++] = (byte)(dist >> 8);
        _symBuf[_symNext++] = (byte)lc;
        dist--;
        _dynL.Freq[LengthCode[lc] + Literals + 1]++;
        _dynD.Freq[DCode(dist)]++;
        return _symNext == SymEnd;
    }

    private void SendBits(int value, int length)
    {
        _biBuf |= (uint)value << _biValid;
        _biValid += length;
        while (_biValid >= 16)
        {
            _out.Add((byte)_biBuf);
            _out.Add((byte)(_biBuf >> 8));
            _biBuf >>= 16;
            _biValid -= 16;
        }
    }

    private void BiWindup()
    {
        if (_biValid > 8)
        {
            _out.Add((byte)_biBuf);
            _out.Add((byte)(_biBuf >> 8));
        }
        else if (_biValid > 0)
        {
            _out.Add((byte)_biBuf);
        }

        _biBuf = 0;
        _biValid = 0;
    }

    private static int BiReverse(int code, int len)
    {
        int res = 0;
        do
        {
            res |= code & 1;
            code >>= 1;
            res <<= 1;
        }
        while (--len > 0);
        return res >> 1;
    }

    private static void GenCodes(ushort[] code, byte[] len, int maxCode, ushort[] blCount)
    {
        var nextCode = new ushort[MaxBits + 1];
        int c = 0;
        for (int bits = 1; bits <= MaxBits; bits++)
        {
            c = (c + blCount[bits - 1]) << 1;
            nextCode[bits] = (ushort)c;
        }

        for (int n = 0; n <= maxCode; n++)
        {
            int l = len[n];
            if (l == 0)
                continue;
            code[n] = (ushort)BiReverse(nextCode[l]++, l);
        }
    }

    private bool Smaller(Tree tree, int n, int m) => tree.Freq[n] < tree.Freq[m] || (tree.Freq[n] == tree.Freq[m] && _depth[n] <= _depth[m]);

    private void PqDownHeap(Tree tree, int k)
    {
        int v = _heap[k];
        int j = k << 1;
        while (j <= _heapLen)
        {
            if (j < _heapLen && Smaller(tree, _heap[j + 1], _heap[j]))
                j++;
            if (Smaller(tree, v, _heap[j]))
                break;
            _heap[k] = _heap[j];
            k = j;
            j <<= 1;
        }

        _heap[k] = v;
    }

    private void GenBitlen(Desc desc)
    {
        var tree = desc.Dyn;
        int maxCode = desc.MaxCode;
        var stree = desc.StaticLen;
        int[] extra = desc.ExtraBits;
        int baseIdx = desc.ExtraBase, maxLength = desc.MaxLength;
        int overflow = 0;
        Array.Clear(_blCount);
        tree.Len[_heap[_heapMax]] = 0;
        int h;
        for (h = _heapMax + 1; h < HeapSize; h++)
        {
            int n = _heap[h];
            int bits = tree.Len[tree.Dad[n]] + 1;
            if (bits > maxLength)
            {
                bits = maxLength;
                overflow++;
            }

            tree.Len[n] = (ushort)bits;
            if (n > maxCode)
                continue;
            _blCount[bits]++;
            int xbits = 0;
            if (n >= baseIdx)
                xbits = extra[n - baseIdx];
            int f = tree.Freq[n];
            _optLen += (long)f * (bits + xbits);
            if (stree is not null)
                _staticLen += (long)f * (stree[n] + xbits);
        }

        if (overflow == 0)
            return;
        do
        {
            int bits = maxLength - 1;
            while (_blCount[bits] == 0)
                bits--;
            _blCount[bits]--;
            _blCount[bits + 1] += 2;
            _blCount[maxLength]--;
            overflow -= 2;
        }
        while (overflow > 0);
        for (int bits = maxLength; bits != 0; bits--)
        {
            int n = _blCount[bits];
            while (n != 0)
            {
                int m = _heap[--h];
                if (m > maxCode)
                    continue;
                if (tree.Len[m] != bits)
                {
                    _optLen += ((long)bits - tree.Len[m]) * tree.Freq[m];
                    tree.Len[m] = (ushort)bits;
                }

                n--;
            }
        }
    }

    private void BuildTree(Desc desc)
    {
        var tree = desc.Dyn;
        var stree = desc.StaticLen;
        int elems = desc.Elems;
        int maxCode = -1;
        _heapLen = 0;
        _heapMax = HeapSize;
        for (int n = 0; n < elems; n++)
        {
            if (tree.Freq[n] != 0)
            {
                _heap[++_heapLen] = maxCode = n;
                _depth[n] = 0;
            }
            else
            {
                tree.Len[n] = 0;
            }
        }

        while (_heapLen < 2)
        {
            int node = _heap[++_heapLen] = maxCode < 2 ? ++maxCode : 0;
            tree.Freq[node] = 1;
            _depth[node] = 0;
            _optLen--;
            if (stree is not null)
                _staticLen -= stree[node];
        }

        desc.MaxCode = maxCode;
        for (int n = _heapLen / 2; n >= 1; n--)
            PqDownHeap(tree, n);
        int nodeIdx = elems;
        do
        {
            int n = _heap[1];
            _heap[1] = _heap[_heapLen--];
            PqDownHeap(tree, 1);
            int m = _heap[1];
            _heap[--_heapMax] = n;
            _heap[--_heapMax] = m;
            tree.Freq[nodeIdx] = (ushort)(tree.Freq[n] + tree.Freq[m]);
            _depth[nodeIdx] = (byte)(Math.Max(_depth[n], _depth[m]) + 1);
            tree.Dad[n] = tree.Dad[m] = (ushort)nodeIdx;
            _heap[1] = nodeIdx++;
            PqDownHeap(tree, 1);
        }
        while (_heapLen >= 2);
        _heap[--_heapMax] = _heap[1];
        GenBitlen(desc);
        // gen_codes: the codes of the dynamic tree
        var nextCode = new ushort[MaxBits + 1];
        int code = 0;
        for (int bits = 1; bits <= MaxBits; bits++)
        {
            code = (code + _blCount[bits - 1]) << 1;
            nextCode[bits] = (ushort)code;
        }

        for (int n = 0; n <= maxCode; n++)
        {
            int len = tree.Len[n];
            if (len == 0)
                continue;
            tree.Code[n] = (ushort)BiReverse(nextCode[len]++, len);
        }
    }

    private static void ScanTree(Tree tree, int maxCode, Tree bl)
    {
        int prevlen = -1, nextlen = tree.Len[0], count = 0, maxCount = 7, minCount = 4;
        if (nextlen == 0)
        {
            maxCount = 138;
            minCount = 3;
        }

        tree.Len[maxCode + 1] = 0xffff;
        for (int n = 0; n <= maxCode; n++)
        {
            int curlen = nextlen;
            nextlen = tree.Len[n + 1];
            if (++count < maxCount && curlen == nextlen)
                continue;
            if (count < minCount)
                bl.Freq[curlen] += (ushort)count;
            else if (curlen != 0)
            {
                if (curlen != prevlen)
                    bl.Freq[curlen]++;
                bl.Freq[Rep3_6]++;
            }
            else if (count <= 10)
                bl.Freq[Repz3_10]++;
            else
                bl.Freq[Repz11_138]++;
            count = 0;
            prevlen = curlen;
            if (nextlen == 0)
            {
                maxCount = 138;
                minCount = 3;
            }
            else if (curlen == nextlen)
            {
                maxCount = 6;
                minCount = 3;
            }
            else
            {
                maxCount = 7;
                minCount = 4;
            }
        }
    }

    private void SendCode(int c, Tree tree) => SendBits(tree.Code[c], tree.Len[c]);

    private void SendTree(Tree tree, int maxCode)
    {
        int prevlen = -1, nextlen = tree.Len[0], count = 0, maxCount = 7, minCount = 4;
        if (nextlen == 0)
        {
            maxCount = 138;
            minCount = 3;
        }

        for (int n = 0; n <= maxCode; n++)
        {
            int curlen = nextlen;
            nextlen = tree.Len[n + 1];
            if (++count < maxCount && curlen == nextlen)
                continue;
            if (count < minCount)
            {
                do
                {
                    SendCode(curlen, _blTree);
                }
                while (--count != 0);
            }
            else if (curlen != 0)
            {
                if (curlen != prevlen)
                {
                    SendCode(curlen, _blTree);
                    count--;
                }

                SendCode(Rep3_6, _blTree);
                SendBits(count - 3, 2);
            }
            else if (count <= 10)
            {
                SendCode(Repz3_10, _blTree);
                SendBits(count - 3, 3);
            }
            else
            {
                SendCode(Repz11_138, _blTree);
                SendBits(count - 11, 7);
            }

            count = 0;
            prevlen = curlen;
            if (nextlen == 0)
            {
                maxCount = 138;
                minCount = 3;
            }
            else if (curlen == nextlen)
            {
                maxCount = 6;
                minCount = 3;
            }
            else
            {
                maxCount = 7;
                minCount = 4;
            }
        }
    }

    private int BuildBlTree()
    {
        ScanTree(_dynL, _lDesc.MaxCode, _blTree);
        ScanTree(_dynD, _dDesc.MaxCode, _blTree);
        BuildTree(_blDesc);
        int maxBlindex;
        for (maxBlindex = BlCodes - 1; maxBlindex >= 3; maxBlindex--)
        {
            if (_blTree.Len[BlOrder[maxBlindex]] != 0)
                break;
        }

        _optLen += 3L * (maxBlindex + 1) + 5 + 5 + 4;
        return maxBlindex;
    }

    private void SendAllTrees(int lcodes, int dcodes, int blcodes)
    {
        SendBits(lcodes - 257, 5);
        SendBits(dcodes - 1, 5);
        SendBits(blcodes - 4, 4);
        for (int rank = 0; rank < blcodes; rank++)
            SendBits(_blTree.Len[BlOrder[rank]], 3);
        SendTree(_dynL, lcodes - 1);
        SendTree(_dynD, dcodes - 1);
    }

    private void CompressBlock(ushort[] lcode, byte[] llen, ushort[] dcode, byte[] dlen)
    {
        int sx = 0;
        while (sx < _symNext)
        {
            int dist = _symBuf[sx++] & 0xff;
            dist += (_symBuf[sx++] & 0xff) << 8;
            int lc = _symBuf[sx++];
            if (dist == 0)
            {
                SendBits(lcode[lc], llen[lc]);
            }
            else
            {
                int code = LengthCode[lc];
                SendBits(lcode[code + Literals + 1], llen[code + Literals + 1]);
                int extra = ExtraLBits[code];
                if (extra != 0)
                {
                    lc -= BaseLength[code];
                    SendBits(lc, extra);
                }

                dist--;
                code = DCode(dist);
                SendBits(dcode[code], dlen[code]);
                extra = ExtraDBits[code];
                if (extra != 0)
                {
                    dist -= BaseDist[code];
                    SendBits(dist, extra);
                }
            }
        }

        SendBits(lcode[EndBlock], llen[EndBlock]);
    }

    private void TrFlushBlock(int bufStart, int storedLen, bool last)
    {
        BuildTree(_lDesc);
        BuildTree(_dDesc);
        int maxBlindex = BuildBlTree();
        long optLenb = (_optLen + 3 + 7) >> 3;
        long staticLenb = (_staticLen + 3 + 7) >> 3;
        if (staticLenb <= optLenb)
            optLenb = staticLenb;
        if (storedLen + 4 <= optLenb && bufStart >= 0)
        {
            SendBits(last ? 1 : 0, 3);
            BiWindup();
            _out.Add((byte)storedLen);
            _out.Add((byte)(storedLen >> 8));
            _out.Add((byte)~storedLen);
            _out.Add((byte)(~storedLen >> 8));
            for (int i = 0; i < storedLen; i++)
                _out.Add(_window[bufStart + i]);
        }
        else if (staticLenb == optLenb)
        {
            SendBits((1 << 1) + (last ? 1 : 0), 3);
            CompressBlock(StaticLCode, StaticLLen, StaticDCode, StaticDLen);
        }
        else
        {
            SendBits((2 << 1) + (last ? 1 : 0), 3);
            SendAllTrees(_lDesc.MaxCode + 1, _dDesc.MaxCode + 1, maxBlindex + 1);
            var lLen = new byte[HeapSize];
            var dLen = new byte[2 * DCodes + 1];
            for (int i = 0; i < lLen.Length; i++)
                lLen[i] = (byte)_dynL.Len[i];
            for (int i = 0; i < dLen.Length; i++)
                dLen[i] = (byte)_dynD.Len[i];
            CompressBlock(_dynL.Code, lLen, _dynD.Code, dLen);
        }

        InitBlock();
        if (last)
            BiWindup();
    }
}
