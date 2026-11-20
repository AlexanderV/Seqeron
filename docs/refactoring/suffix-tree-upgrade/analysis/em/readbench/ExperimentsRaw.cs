#if EXPERIMENTS
using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;
using SuffixTree.Esa;
using SuffixTree.Esa.External;

namespace ReadBench;

// Prototype of the specialised mapped descent of 1b.9a: pointers for every array, the checks, the child layout and the text access as struct type arguments, so that the JIT emits one tight loop per combination.

internal interface IFlag
{
    static abstract bool On { get; }
}

internal struct Yes : IFlag
{
    public static bool On => true;
}

internal struct No : IFlag
{
    public static bool On => false;
}

/// <summary>The text as the descent needs it: a symbol at a position inside the text, and the equality of a stretch of the text with a stretch of the pattern.</summary>
internal interface IRawText
{
    static abstract bool NeedsPackedPattern { get; }

    int Symbol(int position);

    /// <summary>text[start .. start + length) equals pattern[from .. from + length) (the caller checked that the stretch lies inside the text).</summary>
    bool Equals(int start, int from, int length, ReadOnlySpan<byte> pattern, ReadOnlySpan<ulong> packed);
}

internal readonly unsafe struct BytesText(byte* p) : IRawText
{
    public static bool NeedsPackedPattern => false;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Symbol(int position) => p[position];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(int start, int from, int length, ReadOnlySpan<byte> pattern, ReadOnlySpan<ulong> packed)
        => new ReadOnlySpan<byte>(p + start, length).SequenceEqual(pattern.Slice(from, length));
}

/// <summary>The packed text of the container (TEXT codec 4, big-endian bit stream, 64 zero bytes of padding), read without a bounds test: the caller passes positions inside the text.</summary>
internal readonly unsafe struct PackedText(byte* p, int bits) : IRawText
{
    public static bool NeedsPackedPattern => true;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ulong Window(long i)
    {
        long bit = i * bits;
        byte* q = p + (bit >> 3);
        int s = (int)(bit & 7);
        ulong w = BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q));
        // branch free: for s = 0 the shift of the carry byte by 8 leaves nothing
        return (w << s) | ((ulong)q[8] >> (8 - s));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Symbol(int position) => (int)(Window(position) >> (64 - bits));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Equals(int start, int from, int length, ReadOnlySpan<byte> pattern, ReadOnlySpan<ulong> packed)
    {
        int ks = 64 / bits;
        int done = 0;
        while (done < length)
        {
            int rest = length - done;
            ulong x = Window(start + done) ^ PatternWindow(packed, from + done);
            if (rest >= ks)
            {
                if ((x >> (64 - (ks * bits))) != 0) return false;
                done += ks;
            }
            else
            {
                if ((x >> (64 - (rest * bits))) != 0) return false;
                return true;
            }
        }
        return true;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ulong PatternWindow(ReadOnlySpan<ulong> packed, int symbol)
    {
        long bit = (long)symbol * bits;
        int w = (int)(bit >> 6);
        int s = (int)(bit & 63);
        ulong hi = packed[w] << s;
        // the second word exists (the pattern is stored with one word of padding); for s = 0 the shift by 64 must give 0
        ulong lo = s == 0 ? 0 : packed[w + 1] >> (64 - s);
        return hi | lo;
    }

    /// <summary>The pattern as a big-endian bit stream of 64-bit words, the first symbol in the top bits (one zero word of padding follows).</summary>
    internal static void Pack(ReadOnlySpan<byte> pattern, Span<ulong> words, int bits)
    {
        words.Clear();
        long bit = 0;
        for (int i = 0; i < pattern.Length; i++, bit += bits)
        {
            int w = (int)(bit >> 6);
            int s = (int)(bit & 63);
            ulong v = pattern[i];
            int room = 64 - s - bits;
            if (room >= 0) words[w] |= v << room;
            else
            {
                words[w] |= v >> -room;
                words[w + 1] |= v << (64 + room);
            }
        }
    }
}

internal interface IRawChild
{
    int First(int lb, int rb, bool lastChild);

    /// <summary>The l-index after k in the interval that ends at rb, or -1.</summary>
    int Next(int k, int rb);

    /// <summary>Asks for the cache line of the cell that Next(k) reads.</summary>
    void Prefetch(int k);
}

internal interface IRawLcp
{
    int Get(int rank);
}

internal readonly unsafe struct Lcp8(byte* p, CompactLcp fallback) : IRawLcp
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Get(int rank)
    {
        int v = p[rank];
        return v != 255 ? v : fallback[rank];
    }
}

internal readonly unsafe struct RawPackedChild<TLcp, TCheck>(int* cld, TLcp lcp, int n) : IRawChild
    where TLcp : struct, IRawLcp
    where TCheck : struct, IFlag
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int First(int lb, int rb, bool lastChild) => cld[lastChild ? lb : rb];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Prefetch(int k)
    {
        if (Sse.IsSupported && (uint)k <= (uint)n) Sse.Prefetch0(cld + k);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Next(int k, int rb)
    {
        if (k >= rb) return -1;
        int lk = lcp.Get(k);
        if (lcp.Get(k + 1) == lk) return k + 1;
        int v = cld[k];
        if (TCheck.On && (uint)v > (uint)n) throw new InvalidDataException("damaged");
        return lcp.Get(v) == lk ? v : -1;
    }
}

internal readonly unsafe struct RawTaggedChild(uint* cld) : IRawChild
{
    private const uint ValueMask = (1u << 30) - 1;
    private const uint NextTag = 1u << 30;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int First(int lb, int rb, bool lastChild) => (int)(cld[lastChild ? lb : rb] & ValueMask);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Prefetch(int k)
    {
        if (Sse.IsSupported) Sse.Prefetch0(cld + k);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Next(int k, int rb)
    {
        uint c = cld[k];
        return (c & ~ValueMask) == NextTag ? (int)(c & ValueMask) : -1;
    }
}

internal static unsafe class RawDescent
{
    /// <summary>The descent of DescendMapped with raw pointers; TCheck decides at compile time whether the walk is bounded and every value read from the file is validated.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static bool Run<TCheck, TPf, TChild, TLcp, TText>(TChild child, TLcp lcp, TText text, int* sa, int n, ReadOnlySpan<byte> pat, ReadOnlySpan<ulong> packed, out int lo, out int hi)
        where TCheck : struct, IFlag
        where TPf : struct, IFlag
        where TChild : struct, IRawChild
        where TLcp : struct, IRawLcp
        where TText : struct, IRawText
    {
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = child.First(lb, rb, isLast);
            if (TCheck.On && boundary >= 0 && (boundary <= a || boundary > rb)) Damaged();
            if (TPf.On && boundary >= 0 && boundary <= n)
            {
                if (Sse.IsSupported) Sse.Prefetch0(sa + boundary);
                child.Prefetch(boundary);
            }
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = (uint)p < (uint)n ? text.Symbol(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0)
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                a = boundary;
                boundary = child.Next(boundary, rb);
                if (TCheck.On && boundary >= 0 && (boundary <= a || boundary > rb)) Damaged();
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                int start = sa[a] + c;
                if ((uint)start > (uint)n || m - c > n - start || !text.Equals(start, c, m - c, pat, packed))
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                lo = a;
                hi = a;
                return true;
            }

            int first = child.First(a, b, childIsLast);
            if (TCheck.On && (uint)first > (uint)n) Damaged();
            if (TPf.On && first <= n)
            {
                if (Sse.IsSupported) Sse.Prefetch0(sa + first);
                child.Prefetch(first);
            }
            int childEll = lcp.Get(first);
            if (TCheck.On && (childEll <= ell || childEll > n)) Damaged();
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1 && ((uint)from >= (uint)n || !text.Equals(from + 1, c + 1, limit - c - 1, pat, packed)))
            {
                lo = 0;
                hi = -1;
                return false;
            }
            c = limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        lo = lb;
        hi = rb;
        return true;
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void Damaged() => throw new InvalidDataException("damaged");
}

/// <summary>Query objects over the prototype descent: the data are the mapped arrays of the file (or copies of them) and a text of bytes or the packed text of the file.</summary>
internal static unsafe class ExperimentsRaw
{
    private static readonly List<object> s_keep = [];

    internal static Program.Q Create(string variant, string file)
    {
        // ix-rx-<arrays>-<text>-<child>-<check>: arrays m (mapped file) | r (resident copies), text b (bytes) | p (packed), child p (packed, the file's) | t (tagged), check c | t (trusted)
        string[] parts = variant.Split('-');
        char arrays = parts[2][0], text = parts[3][0], child = parts[4][0], check = parts[5][0];
        SesaFile f = SesaFile.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Checked });
        EsaIndex index = EsaIndex.FromFile(f);
        int n = index.Length;
        int* sa = index.Sa.Pointer;
        int* cld = (int*)(f.Data + f.Parsed.Child.Offset);
        byte* lcpBase = index.Lcp.MappedBase;
        byte* textPacked = f.TextReader_Pointer();
        CompactLcp lcp = index.Lcp;
        if (arrays == 'r')
        {
            int[] saCopy = GC.AllocateUninitializedArray<int>(n + 1, pinned: true);
            new ReadOnlySpan<int>(sa, n + 1).CopyTo(saCopy);
            int[] cldCopy = GC.AllocateUninitializedArray<int>(n + 1, pinned: true);
            new ReadOnlySpan<int>(cld, n + 1).CopyTo(cldCopy);
            byte[] lcpCopy = GC.AllocateUninitializedArray<byte>(n + 1, pinned: true);
            new ReadOnlySpan<byte>(lcpBase, n + 1).CopyTo(lcpCopy);
            long tl = PackedTextFormat.TotalBytes(n, f.TextReader.Bits);
            byte[] textCopy = GC.AllocateUninitializedArray<byte>((int)tl, pinned: true);
            new ReadOnlySpan<byte>(textPacked, (int)tl).CopyTo(textCopy);
            s_keep.AddRange([saCopy, cldCopy, lcpCopy, textCopy]);
            sa = (int*)Marshal.UnsafeAddrOfPinnedArrayElement(saCopy, 0);
            cld = (int*)Marshal.UnsafeAddrOfPinnedArrayElement(cldCopy, 0);
            lcpBase = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(lcpCopy, 0);
            textPacked = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(textCopy, 0);
        }
        byte* textBytes = null;
        if (text == 'b')
        {
            byte[] bytes = GC.AllocateUninitializedArray<byte>(n + 64, pinned: true);
            f.TextReader.ReadSymbols(0, bytes.AsSpan(0, n));
            bytes.AsSpan(n).Clear();
            s_keep.Add(bytes);
            textBytes = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(bytes, 0);
        }
        uint* tagged = null;
        if (child == 't')
        {
            ChildTable t = ChildTable.Build(ChildLayoutKind.Tagged, n, lcp, out _);
            object boxed = t.Tagged;
            var cells = (uint[])typeof(TaggedChild).GetField("_cld", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!.GetValue(boxed)!;
            s_keep.Add(cells);
            tagged = (uint*)GCHandle.Alloc(cells, GCHandleType.Pinned).AddrOfPinnedObject();
        }
        return new RawQ(index, n, sa, cld, tagged, lcpBase, lcp, textPacked, textBytes, f.TextReader.Bits, text, child, check, parts.Length > 6 && parts[6][0] == 'f', f);
    }

    private static byte* TextReader_Pointer(this SesaFile f) => f.Data + f.Parsed.Text.Offset;

    private sealed class RawQ : Program.Q
    {
        private readonly int[] _map;
        private readonly int _n;
        private readonly int* _sa;
        private readonly int* _cld;
        private readonly uint* _tagged;
        private readonly byte* _lcp8;
        private readonly CompactLcp _lcp;
        private readonly byte* _textPacked;
        private readonly byte* _textBytes;
        private readonly int _bits;
        private readonly char _text, _child, _check;
        private readonly int _code;
        private readonly EsaIndex _index;
        private readonly SesaFile _file;

        internal RawQ(EsaIndex index, int n, int* sa, int* cld, uint* tagged, byte* lcp8, CompactLcp lcp, byte* textPacked, byte* textBytes, int bits, char text, char child, char check, bool prefetch, SesaFile file)
        {
            _index = index;
            _file = file;
            _n = n;
            _sa = sa;
            _cld = cld;
            _tagged = tagged;
            _lcp8 = lcp8;
            _lcp = lcp;
            _textPacked = textPacked;
            _textBytes = textBytes;
            _bits = bits;
            _text = text;
            _child = child;
            _check = check;
            _code = (check == 'c' ? 0 : 8) + (prefetch ? 4 : 0) + (child == 'p' ? 0 : 2) + (text == 'b' ? 0 : 1);
            char[] alphabet = index.TextView.Alphabet;
            _map = new int[alphabet[^1] + 1];
            Array.Fill(_map, -1);
            for (int s = 0; s < alphabet.Length; s++) _map[alphabet[s]] = s;
        }

        private bool Match(ReadOnlySpan<byte> pat, out int lo, out int hi)
        {
            Span<ulong> packed = stackalloc ulong[(((pat.Length * _bits) + 63) / 64) + 2];
            if (_text == 'p') PackedText.Pack(pat, packed, _bits);
            var l8 = new Lcp8(_lcp8, _lcp);
            switch (_code)
            {
                case 0: return RawDescent.Run<Yes, No, RawPackedChild<Lcp8, Yes>, Lcp8, BytesText>(new RawPackedChild<Lcp8, Yes>(_cld, l8, _n), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 1: return RawDescent.Run<Yes, No, RawPackedChild<Lcp8, Yes>, Lcp8, PackedText>(new RawPackedChild<Lcp8, Yes>(_cld, l8, _n), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 2: return RawDescent.Run<Yes, No, RawTaggedChild, Lcp8, BytesText>(new RawTaggedChild(_tagged), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 3: return RawDescent.Run<Yes, No, RawTaggedChild, Lcp8, PackedText>(new RawTaggedChild(_tagged), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 4: return RawDescent.Run<Yes, Yes, RawPackedChild<Lcp8, Yes>, Lcp8, BytesText>(new RawPackedChild<Lcp8, Yes>(_cld, l8, _n), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 5: return RawDescent.Run<Yes, Yes, RawPackedChild<Lcp8, Yes>, Lcp8, PackedText>(new RawPackedChild<Lcp8, Yes>(_cld, l8, _n), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 6: return RawDescent.Run<Yes, Yes, RawTaggedChild, Lcp8, BytesText>(new RawTaggedChild(_tagged), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 7: return RawDescent.Run<Yes, Yes, RawTaggedChild, Lcp8, PackedText>(new RawTaggedChild(_tagged), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 8: return RawDescent.Run<No, No, RawPackedChild<Lcp8, No>, Lcp8, BytesText>(new RawPackedChild<Lcp8, No>(_cld, l8, _n), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 9: return RawDescent.Run<No, No, RawPackedChild<Lcp8, No>, Lcp8, PackedText>(new RawPackedChild<Lcp8, No>(_cld, l8, _n), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 10: return RawDescent.Run<No, No, RawTaggedChild, Lcp8, BytesText>(new RawTaggedChild(_tagged), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 11: return RawDescent.Run<No, No, RawTaggedChild, Lcp8, PackedText>(new RawTaggedChild(_tagged), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 12: return RawDescent.Run<No, Yes, RawPackedChild<Lcp8, No>, Lcp8, BytesText>(new RawPackedChild<Lcp8, No>(_cld, l8, _n), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 13: return RawDescent.Run<No, Yes, RawPackedChild<Lcp8, No>, Lcp8, PackedText>(new RawPackedChild<Lcp8, No>(_cld, l8, _n), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                case 14: return RawDescent.Run<No, Yes, RawTaggedChild, Lcp8, BytesText>(new RawTaggedChild(_tagged), l8, new BytesText(_textBytes), _sa, _n, pat, packed, out lo, out hi);
                case 15: return RawDescent.Run<No, Yes, RawTaggedChild, Lcp8, PackedText>(new RawTaggedChild(_tagged), l8, new PackedText(_textPacked, _bits), _sa, _n, pat, packed, out lo, out hi);
                default: throw new InvalidOperationException();
            }
        }

        private bool Dense(string pattern, Span<byte> pat)
        {
            for (int i = 0; i < pattern.Length; i++)
            {
                char c = pattern[i];
                int s = c < _map.Length ? _map[c] : -1;
                if (s < 0) return false;
                pat[i] = (byte)s;
            }
            return true;
        }

        internal override int Count(string pattern)
        {
            Span<byte> pat = stackalloc byte[pattern.Length];
            if (!Dense(pattern, pat)) return 0;
            return Match(pat, out int lo, out int hi) ? hi - lo + 1 : 0;
        }

        internal override int Locate(string pattern)
        {
            Span<byte> pat = stackalloc byte[pattern.Length];
            if (!Dense(pattern, pat)) return 0;
            if (!Match(pat, out int lo, out int hi)) return 0;
            var result = new int[hi - lo + 1];
            new ReadOnlySpan<int>(_sa + lo, result.Length).CopyTo(result);
            return result.Length;
        }
    }
}
#endif
