#if EXPERIMENTS
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using SuffixTree;
using SuffixTree.Esa;
using SuffixTree.Esa.External;

namespace ReadBench;

/// <summary>
/// Engine-level experiments of the diagnosis of 1b.9a: the descent is run on the index directly (no facade), so that one factor at a time can be exchanged. A variant is "ix-" followed by factors:
///   ix-memory / ix-memory-packed   the in-memory index (child layout tagged / packed)
///   ix-file                        the shipped mapped arm (EsaIndex.TryMatch over a mapped container)
///   ix-bytes                       the mapped arrays, but the text is a byte per symbol (a pinned copy): the packed decode is gone
///   ix-raw                         pointers for the mapped arrays and the byte text, nothing checked; the packed child table of the file
///   ix-rawtagged                   the same with a tagged child table derived from the file (a managed uint[])
///   ix-res-packed / ix-res-tagged  the sections copied into managed arrays and run through the in-memory generic descent (the engine of the in-memory tree)
/// </summary>
internal static unsafe class ExperimentsIx
{
    private static byte[]? s_textBytes;
    private static readonly List<object> s_keep = [];

    internal static Program.Q Create(string variant, string file, string text)
    {
        if (variant.StartsWith("ix-rx-", StringComparison.Ordinal)) return ExperimentsRaw.Create(variant, file);
        switch (variant)
        {
            case "ix-memory": return new IndexQ(EsaIndex.Build(text, ChildLayoutKind.Tagged));
            case "ix-memory-packed": return new IndexQ(EsaIndex.Build(text, ChildLayoutKind.Packed));
            case "ix-file": return new IndexQ(EsaIndex.FromFile(SesaFile.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Checked })));
            case "ix-file-trusted": return new IndexQ(EsaIndex.FromFile(SesaFile.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Trusted })));
        }

        SesaFile f = SesaFile.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Checked });
        EsaIndex index = EsaIndex.FromFile(f);
        int n = index.Length;
        s_textBytes = GC.AllocateUninitializedArray<byte>(n + 64, pinned: true);
        s_keep.Add(s_textBytes);
        f.TextReader.ReadSymbols(0, s_textBytes.AsSpan(0, n));
        s_textBytes.AsSpan(n).Clear();
        switch (variant)
        {
            case "ix-bytes": return new BytesQ(index, f, s_textBytes);
            case "ix-raw": return new RawQ(index, f, s_textBytes, tagged: false);
            case "ix-rawtagged": return new RawQ(index, f, s_textBytes, tagged: true);
            case "ix-res-packed": return new ResidentQ(index, s_textBytes, ChildLayoutKind.Packed);
            case "ix-res-tagged": return new ResidentQ(index, s_textBytes, ChildLayoutKind.Tagged);
            default: throw new ArgumentException($"unknown variant '{variant}'");
        }
    }

    /// <summary>The dense symbols of a pattern, as the facade maps them (a code unit that is not in the alphabet is no match).</summary>
    private static bool Dense(int[] map, string pattern, Span<byte> pat)
    {
        for (int i = 0; i < pattern.Length; i++)
        {
            char c = pattern[i];
            int s = c < map.Length ? map[c] : -1;
            if (s < 0) return false;
            pat[i] = (byte)s;
        }
        return true;
    }

    private static int[] SymbolMap(char[] alphabet)
    {
        var map = new int[alphabet[^1] + 1];
        Array.Fill(map, -1);
        for (int s = 0; s < alphabet.Length; s++) map[alphabet[s]] = s;
        return map;
    }

    private sealed class IndexQ(EsaIndex index) : Program.Q
    {
        internal override int Count(string pattern) => index.CountOccurrences(pattern.AsSpan());

        internal override int Locate(string pattern) => index.FindAllOccurrences(pattern.AsSpan()).Count;
    }

    private abstract class DenseQ : Program.Q
    {
        private readonly int[] _map;
        protected readonly EsaIndex Index;

        protected DenseQ(EsaIndex index)
        {
            Index = index;
            _map = SymbolMap(index.TextView.Alphabet);
        }

        protected abstract bool Match(ReadOnlySpan<byte> pat, out int lo, out int hi);

        internal override int Count(string pattern)
        {
            Span<byte> pat = stackalloc byte[pattern.Length];
            if (!Dense(_map, pattern, pat)) return 0;
            return Match(pat, out int lo, out int hi) ? hi - lo + 1 : 0;
        }

        internal override int Locate(string pattern)
        {
            Span<byte> pat = stackalloc byte[pattern.Length];
            if (!Dense(_map, pattern, pat)) return 0;
            if (!Match(pat, out int lo, out int hi)) return 0;
            return Index.Sa.ToArray(lo, hi - lo + 1).Length;
        }
    }

    /// <summary>The shipped mapped descent with the text replaced by bytes: a copy of DescendMapped whose text reads are plain loads.</summary>
    private sealed class BytesQ : DenseQ
    {
        private readonly SesaFile _f;
        private readonly byte[] _text;

        internal BytesQ(EsaIndex index, SesaFile f, byte[] text) : base(index)
        {
            _f = f;
            _text = text;
        }

        protected override bool Match(ReadOnlySpan<byte> pat, out int lo, out int hi)
            => DescendBytes(Index.Child.Mapped, Index.Sa, Index.Lcp, _text, Index.Length, pat, _f.IsChecked, out lo, out hi);
    }

    internal static bool DescendBytes(MappedChild child, SaView sa, CompactLcp lcp, byte[] textArray, int n, ReadOnlySpan<byte> pat, bool isChecked, out int lo, out int hi)
    {
        ReadOnlySpan<byte> text = textArray;
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = child.FirstLIndex(lb, rb, isLast);
            if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) throw new InvalidDataException("damaged");
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = (uint)p < (uint)n ? text[p] : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0) { lo = 0; hi = -1; return false; }
                a = boundary;
                boundary = child.NextLIndex(boundary, rb);
                if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) throw new InvalidDataException("damaged");
            }
            bool childIsLast = boundary < 0;
            if (a == b)
            {
                int start = sa[a] + c;
                if ((uint)start > (uint)n || m - c > n - start || !text.Slice(start, m - c).SequenceEqual(pat.Slice(c))) { lo = 0; hi = -1; return false; }
                lo = a; hi = a; return true;
            }
            int childEll = lcp[child.FirstLIndex(a, b, childIsLast)];
            if (isChecked && (childEll <= ell || childEll > n)) throw new InvalidDataException("damaged");
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1 && ((uint)from >= (uint)n || !text.Slice(from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))) { lo = 0; hi = -1; return false; }
            c = limit; lb = a; rb = b; ell = childEll; isLast = childIsLast;
        }
        lo = lb; hi = rb;
        return true;
    }

    /// <summary>Pointers for every array of the mapped file and the byte text, no checks: the descent as tight as the data allow, with the packed child table of the file or a tagged one derived from it.</summary>
    private sealed class RawQ : DenseQ
    {
        private readonly int* _sa;
        private readonly int* _cld;
        private readonly uint* _tagged;
        private readonly byte* _lcp8;
        private readonly CompactLcp _lcp;
        private readonly byte* _text;
        private readonly int _n;
        private readonly bool _isTagged;

        internal RawQ(EsaIndex index, SesaFile f, byte[] text, bool tagged) : base(index)
        {
            _n = index.Length;
            _sa = index.Sa.Pointer;
            _lcp = index.Lcp;
            _lcp8 = _lcp.Kind == LcpKind.U8 ? _lcp.MappedBase : null;
            _text = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(text, 0);
            _isTagged = tagged;
            if (tagged)
            {
                // the tagged table the in-memory tree would build from the same LCP array; its cells are a private uint[] of the struct: read it by reflection and keep it pinned
                ChildTable t = ChildTable.Build(ChildLayoutKind.Tagged, _n, _lcp, out _);
                object boxed = t.Tagged;
                var cells = (uint[])typeof(TaggedChild).GetField("_cld", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!.GetValue(boxed)!;
                s_keep.Add(cells);
                _tagged = (uint*)GCHandle.Alloc(cells, GCHandleType.Pinned).AddrOfPinnedObject();
                _cld = null;
            }
            else
            {
                _cld = (int*)((byte*)f.Data + f.Parsed.Child.Offset);
                _tagged = null;
            }
        }

        protected override bool Match(ReadOnlySpan<byte> pat, out int lo, out int hi)
        {
            return _isTagged ? DescendRawTagged(pat, out lo, out hi) : DescendRawPacked(pat, out lo, out hi);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int Lcp(int rank)
        {
            if (_lcp8 is not null)
            {
                int v = _lcp8[rank];
                if (v != 255) return v;
            }
            return _lcp[rank];
        }

        private bool DescendRawPacked(ReadOnlySpan<byte> pat, out int lo, out int hi)
        {
            int n = _n, m = pat.Length;
            int lb = 0, rb = n, ell = 0, c = 0;
            bool isLast = true;
            while (c < m)
            {
                int symbol = pat[c];
                int a = lb;
                int boundary = _cld[isLast ? lb : rb];
                int b;
                while (true)
                {
                    b = boundary >= 0 ? boundary - 1 : rb;
                    int p = _sa[a] + ell;
                    int key = (uint)p < (uint)n ? _text[p] : -1;
                    if (key == symbol) break;
                    if (key > symbol || boundary < 0) { lo = 0; hi = -1; return false; }
                    a = boundary;
                    // NextLIndex of the packed layout
                    int k = boundary;
                    if (k >= rb) boundary = -1;
                    else
                    {
                        int lk = Lcp(k);
                        if (Lcp(k + 1) == lk) boundary = k + 1;
                        else
                        {
                            int v = _cld[k];
                            boundary = Lcp(v) == lk ? v : -1;
                        }
                    }
                }
                bool childIsLast = boundary < 0;
                if (a == b)
                {
                    int start = _sa[a] + c;
                    if ((uint)start > (uint)n || m - c > n - start || !new ReadOnlySpan<byte>(_text + start, m - c).SequenceEqual(pat.Slice(c))) { lo = 0; hi = -1; return false; }
                    lo = a; hi = a; return true;
                }
                int childEll = Lcp(_cld[childIsLast ? a : b]);
                int limit = Math.Min(childEll, m);
                int from = _sa[a] + c;
                if (limit > c + 1 && ((uint)from >= (uint)n || !new ReadOnlySpan<byte>(_text + from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))) { lo = 0; hi = -1; return false; }
                c = limit; lb = a; rb = b; ell = childEll; isLast = childIsLast;
            }
            lo = lb; hi = rb;
            return true;
        }

        private bool DescendRawTagged(ReadOnlySpan<byte> pat, out int lo, out int hi)
        {
            const uint ValueMask = (1u << 30) - 1;
            const uint NextTag = 1u << 30;
            int n = _n, m = pat.Length;
            int lb = 0, rb = n, ell = 0, c = 0;
            bool isLast = true;
            while (c < m)
            {
                int symbol = pat[c];
                int a = lb;
                int boundary = (int)(_tagged[isLast ? lb : rb] & ValueMask);
                int b;
                while (true)
                {
                    b = boundary >= 0 ? boundary - 1 : rb;
                    int p = _sa[a] + ell;
                    int key = (uint)p < (uint)n ? _text[p] : -1;
                    if (key == symbol) break;
                    if (key > symbol || boundary < 0) { lo = 0; hi = -1; return false; }
                    a = boundary;
                    uint cell = _tagged[boundary];
                    boundary = (cell & ~ValueMask) == NextTag ? (int)(cell & ValueMask) : -1;
                }
                bool childIsLast = boundary < 0;
                if (a == b)
                {
                    int start = _sa[a] + c;
                    if ((uint)start > (uint)n || m - c > n - start || !new ReadOnlySpan<byte>(_text + start, m - c).SequenceEqual(pat.Slice(c))) { lo = 0; hi = -1; return false; }
                    lo = a; hi = a; return true;
                }
                int childEll = Lcp((int)(_tagged[childIsLast ? a : b] & ValueMask));
                int limit = Math.Min(childEll, m);
                int from = _sa[a] + c;
                if (limit > c + 1 && ((uint)from >= (uint)n || !new ReadOnlySpan<byte>(_text + from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))) { lo = 0; hi = -1; return false; }
                c = limit; lb = a; rb = b; ell = childEll; isLast = childIsLast;
            }
            lo = lb; hi = rb;
            return true;
        }
    }

    /// <summary>The sections copied into managed arrays and run through the in-memory generic descent: what a resident open mode would run.</summary>
    private sealed class ResidentQ : DenseQ
    {
        private readonly int[] _sa;
        private readonly CompactLcp _lcp;
        private readonly byte[] _text;
        private readonly ChildLayoutKind _kind;
        private readonly TaggedChild _tagged;
        private readonly PackedChild _packed;

        internal ResidentQ(EsaIndex index, byte[] text, ChildLayoutKind kind) : base(index)
        {
            _text = text;
            _kind = kind;
            _sa = index.Sa.ToArray(0, index.Sa.Length);
            _lcp = CompactLcp.FromDense(index.Lcp.ToArray());
            ChildTable t = ChildTable.Build(kind, index.Length, _lcp, out _);
            _tagged = t.Tagged;
            _packed = t.Packed;
        }

        protected override bool Match(ReadOnlySpan<byte> pat, out int lo, out int hi)
        {
            ReadOnlySpan<byte> text = _text.AsSpan(0, Index.Length);
            return _kind == ChildLayoutKind.Tagged
                ? EsaIndex.Descend(_tagged, _sa, _lcp, text, pat, out lo, out hi)
                : EsaIndex.Descend(_packed, _sa, _lcp, text, pat, out lo, out hi);
        }
    }
}
#endif
