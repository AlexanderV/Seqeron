using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa;

/// <summary>Whether a pattern descent over a container checks what it reads (a compile-time constant of the instantiation, so that the checks of <see cref="EsaAccessMode.Trusted"/> are not in its code).</summary>
internal interface IAccessChecks
{
    /// <summary>True in <see cref="EsaAccessMode.Checked"/> mode.</summary>
    static abstract bool IsChecked { get; }
}

/// <summary><see cref="EsaAccessMode.Checked"/>: every value read from the file that is used as an index is validated, and the walk is bounded.</summary>
internal readonly struct CheckedAccess : IAccessChecks
{
    public static bool IsChecked => true;
}

/// <summary><see cref="EsaAccessMode.Trusted"/>: no checks in the hot loop.</summary>
internal readonly struct TrustedAccess : IAccessChecks
{
    public static bool IsChecked => false;
}

/// <summary>One LCP entry of a container by rank, in the representation of the file (<see cref="LcpKind"/>); escaped entries are resolved through the overflow list of <see cref="CompactLcp"/>.</summary>
internal interface IMappedLcp
{
    int this[int rank] { get; }
}

/// <summary>The LCP base array of 8-bit entries (escape value 255).</summary>
internal readonly unsafe struct Lcp8Reader(byte* values, CompactLcp lcp) : IMappedLcp
{
    public int this[int rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            int v = values[rank];
            return v != byte.MaxValue ? v : lcp[rank];
        }
    }
}

/// <summary>The LCP base array of 16-bit entries (escape value 65,535).</summary>
internal readonly unsafe struct Lcp16Reader(ushort* values, CompactLcp lcp) : IMappedLcp
{
    public int this[int rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            int v = values[rank];
            return v != ushort.MaxValue ? v : lcp[rank];
        }
    }
}

/// <summary>The LCP base array of 32-bit entries; in checked mode a value above the length of the array is damage.</summary>
internal readonly unsafe struct Lcp32Reader<TAccess>(int* values, int length) : IMappedLcp
    where TAccess : struct, IAccessChecks
{
    public int this[int rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            int v = values[rank];
            if (TAccess.IsChecked && (uint)v > (uint)length) MappedSearch.ThrowDamaged();
            return v;
        }
    }
}

/// <summary>The cells of the child table of a container as the descent reads them: the first l-index of an interval and the l-index after one.</summary>
internal interface IMappedCells
{
    /// <summary>The first l-index of the interval <c>[lb..rb]</c> (<paramref name="lastChild"/>: the cell is at its left bound, else at its right bound).</summary>
    int First(int lb, int rb, bool lastChild);

    /// <summary>The l-index after <paramref name="k"/> in the interval that ends at <paramref name="rb"/>, or -1.</summary>
    int Next(int k, int rb);
}

/// <summary>CHILD codec 1 (layout of <see cref="PackedChild"/>): the kind of a field is told by comparing LCP values, so a step to the next l-index reads the LCP array.</summary>
internal readonly unsafe struct PackedCells<TAccess, TLcp>(int* cld, TLcp lcp, int n) : IMappedCells
    where TAccess : struct, IAccessChecks
    where TLcp : struct, IMappedLcp
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int First(int lb, int rb, bool lastChild) => cld[lastChild ? lb : rb];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Next(int k, int rb)
    {
        if (k >= rb) return -1;
        int lk = lcp[k];
        if (lcp[k + 1] == lk) return k + 1;
        int v = cld[k];
        if (TAccess.IsChecked && (uint)v > (uint)n) MappedSearch.ThrowDamaged();
        return lcp[v] == lk ? v : -1;
    }
}

/// <summary>CHILD codec 2 (layout of <see cref="TaggedChild"/>): the kind of a field is in bits 30-31, no LCP read.</summary>
internal readonly unsafe struct TaggedCells(uint* cld) : IMappedCells
{
    private const uint ValueMask = (1u << 30) - 1;
    private const uint KindMask = 3u << 30;
    private const uint NextKind = 1u << 30;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int First(int lb, int rb, bool lastChild) => (int)(cld[lastChild ? lb : rb] & ValueMask);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Next(int k, int rb)
    {
        uint cell = cld[k];
        return (cell & KindMask) == NextKind ? (int)(cell & ValueMask) : -1;
    }
}

/// <summary>
/// The pattern descent over the sections of a container (work package 1b.9a): <see cref="EsaIndex.Descend{TChild, TSym}"/> for the child table (CHILD codec 1 or 2), the suffix array, the LCP array and the packed text (TEXT codec 4)
/// read in place, written as one loop over raw pointers. The arrays are fetched once, at the construction; the access mode, the LCP width and the codec of the child table are fixed at compile time (a struct type argument each), so
/// that the JIT emits a loop of its own for each of the twelve combinations that a file can need and none of them tests a field that is the same for the whole file. The text is compared in 64-bit windows against the pattern packed the same way
/// (a stretch of up to <c>64 / b</c> symbols per step), and the packed pattern is only built when the first stretch has to be compared.
/// <para>
/// The reads are those of the readable reference <see cref="EsaIndex.DescendMapped"/> with the same checks, and the checks are the same: in <see cref="EsaAccessMode.Checked"/> mode every rank that the walk takes from the child table
/// must increase and stay inside its interval, a child must be deeper than its parent, and a rank or a value that was read from a section and is used as an index is validated before it is used, so a damaged file is an
/// <see cref="InvalidDataException"/> and never a read outside a section; the text positions that come from the suffix array are clamped to the text (a position outside it reads nothing), as they are in the reference.
/// </para>
/// </summary>
internal sealed unsafe class MappedSearch
{
    private readonly int* _sa;
    private readonly int* _child;
    private readonly byte* _text;
    private readonly int _bits;
    private readonly int _n;
    private readonly CompactLcp _lcp;
    private readonly int _variant;

    /// <summary>The search over the sections of an opened container.</summary>
    /// <param name="file">The container; it must outlive this object (the facade owns both).</param>
    internal MappedSearch(SesaFile file)
    {
        ArgumentNullException.ThrowIfNull(file);
        _sa = file.Sa.Pointer;
        _child = file.Child.Mapped.Pointer;
        _text = file.TextReader.Pointer;
        _bits = file.TextReader.Bits;
        _n = (int)file.TextReader.Length;
        _lcp = file.Lcp;
        int kind = _lcp.Kind switch { LcpKind.U8 => 0, LcpKind.U16 => 1, _ => 2 };
        _variant = (file.IsChecked ? 0 : 6) + (2 * kind) + (file.Child.Mapped.IsTagged ? 1 : 0);
    }

    /// <summary>
    /// The suffix-array interval <c>[lo..hi]</c> of the NON-EMPTY <paramref name="pattern"/> of dense symbols (each below the alphabet size of the text); false when it does not occur.
    /// </summary>
    /// <exception cref="InvalidDataException">Checked mode: a table of the file is inconsistent.</exception>
    internal bool TryMatch(ReadOnlySpan<byte> pattern, out int lo, out int hi)
    {
        int words = (int)((((long)pattern.Length * _bits) + 63) / 64) + 2;
        Span<ulong> packed = words <= 128 ? stackalloc ulong[words] : new ulong[words];
        switch (_variant)
        {
            case 0: { var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<CheckedAccess, Lcp8Reader, PackedCells<CheckedAccess, Lcp8Reader>>(new PackedCells<CheckedAccess, Lcp8Reader>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            case 1: { var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<CheckedAccess, Lcp8Reader, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
            case 2: { var lcp = new Lcp16Reader((ushort*)_lcp.MappedBase, _lcp); return Run<CheckedAccess, Lcp16Reader, PackedCells<CheckedAccess, Lcp16Reader>>(new PackedCells<CheckedAccess, Lcp16Reader>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            case 3: { var lcp = new Lcp16Reader((ushort*)_lcp.MappedBase, _lcp); return Run<CheckedAccess, Lcp16Reader, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
            case 4: { var lcp = new Lcp32Reader<CheckedAccess>((int*)_lcp.MappedBase, _lcp.Length); return Run<CheckedAccess, Lcp32Reader<CheckedAccess>, PackedCells<CheckedAccess, Lcp32Reader<CheckedAccess>>>(new PackedCells<CheckedAccess, Lcp32Reader<CheckedAccess>>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            case 5: { var lcp = new Lcp32Reader<CheckedAccess>((int*)_lcp.MappedBase, _lcp.Length); return Run<CheckedAccess, Lcp32Reader<CheckedAccess>, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
            case 6: { var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<TrustedAccess, Lcp8Reader, PackedCells<TrustedAccess, Lcp8Reader>>(new PackedCells<TrustedAccess, Lcp8Reader>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            case 7: { var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<TrustedAccess, Lcp8Reader, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
            case 8: { var lcp = new Lcp16Reader((ushort*)_lcp.MappedBase, _lcp); return Run<TrustedAccess, Lcp16Reader, PackedCells<TrustedAccess, Lcp16Reader>>(new PackedCells<TrustedAccess, Lcp16Reader>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            case 9: { var lcp = new Lcp16Reader((ushort*)_lcp.MappedBase, _lcp); return Run<TrustedAccess, Lcp16Reader, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
            case 10: { var lcp = new Lcp32Reader<TrustedAccess>((int*)_lcp.MappedBase, _lcp.Length); return Run<TrustedAccess, Lcp32Reader<TrustedAccess>, PackedCells<TrustedAccess, Lcp32Reader<TrustedAccess>>>(new PackedCells<TrustedAccess, Lcp32Reader<TrustedAccess>>(_child, lcp, _n), lcp, pattern, packed, out lo, out hi); }
            default: { var lcp = new Lcp32Reader<TrustedAccess>((int*)_lcp.MappedBase, _lcp.Length); return Run<TrustedAccess, Lcp32Reader<TrustedAccess>, TaggedCells>(new TaggedCells((uint*)_child), lcp, pattern, packed, out lo, out hi); }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private bool Run<TAccess, TLcp, TCells>(TCells cells, TLcp lcp, ReadOnlySpan<byte> pat, Span<ulong> packed, out int lo, out int hi)
        where TAccess : struct, IAccessChecks
        where TLcp : struct, IMappedLcp
        where TCells : struct, IMappedCells
    {
        int* sa = _sa;
        byte* text = _text;
        int bits = _bits;
        int n = _n;
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        bool patternPacked = false;
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = cells.First(lb, rb, isLast);
            if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamaged();
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = (uint)p < (uint)n ? Symbol(text, bits, p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0)
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                a = boundary;
                boundary = cells.Next(a, rb);
                if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamaged();
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                // a leaf: the rest of the pattern must be a prefix of the rest of the suffix
                int start = sa[a] + c;
                if ((uint)start > (uint)n || m - c > n - start || !Equal(text, bits, start, pat, packed, c, m - c, ref patternPacked))
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                lo = a;
                hi = a;
                return true;
            }

            int first = cells.First(a, b, childIsLast);
            if (TAccess.IsChecked && (uint)first > (uint)n) ThrowDamaged();
            int childEll = lcp[first];
            if (TAccess.IsChecked && (childEll <= ell || childEll > n)) ThrowDamaged();
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1 && ((uint)from >= (uint)n || limit - c - 1 > n - (from + 1) || !Equal(text, bits, from + 1, pat, packed, c + 1, limit - c - 1, ref patternPacked)))
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

    /// <summary>The 64 stream bits of the packed text that start at symbol <paramref name="i"/> (inside the text), the first symbol in the top bits.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong Window(byte* text, int bits, long i)
    {
        long bit = i * bits;
        byte* q = text + (bit >> 3);
        int s = (int)(bit & 7);
        ulong w = BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q));
        // without a test of s: for s = 0 the carry byte shifted by 8 is 0
        return (w << s) | ((ulong)q[8] >> (8 - s));
    }

    /// <summary>The dense symbol at <paramref name="i"/> (inside the text).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int Symbol(byte* text, int bits, long i) => (int)(Window(text, bits, i) >> (64 - bits));

    /// <summary>
    /// Whether the text symbols <c>[start, start + length)</c> (inside the text) equal the pattern symbols <c>[from, from + length)</c>, a stretch of <c>64 / bits</c> symbols per step, the pattern packed on the first use.
    /// </summary>
    private static bool Equal(byte* text, int bits, int start, ReadOnlySpan<byte> pattern, Span<ulong> packed, int from, int length, ref bool patternPacked)
    {
        if (!patternPacked)
        {
            PackPattern(pattern, packed, bits);
            patternPacked = true;
        }
        int per = 64 / bits;
        for (int done = 0; done < length; done += per)
        {
            ulong x = Window(text, bits, start + done) ^ PatternWindow(packed, bits, from + done);
            int take = Math.Min(per, length - done);
            if ((x >> (64 - (take * bits))) != 0) return false;
        }
        return true;
    }

    /// <summary>The 64 bits of the packed pattern that start at symbol <paramref name="symbol"/> (the stream has one zero word more than the pattern needs).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong PatternWindow(Span<ulong> packed, int bits, int symbol)
    {
        long bit = (long)symbol * bits;
        int w = (int)(bit >> 6);
        int s = (int)(bit & 63);
        ulong high = packed[w] << s;
        return s == 0 ? high : high | (packed[w + 1] >> (64 - s));
    }

    /// <summary>The pattern as a big-endian bit stream of 64-bit words (the first symbol in the top bits, <paramref name="bits"/> per symbol), the unused bits and the padding words zero.</summary>
    private static void PackPattern(ReadOnlySpan<byte> pattern, Span<ulong> words, int bits)
    {
        words.Clear();
        long bit = 0;
        for (int i = 0; i < pattern.Length; i++, bit += bits)
        {
            int w = (int)(bit >> 6);
            int room = 64 - (int)(bit & 63) - bits;
            ulong v = pattern[i];
            if (room >= 0)
            {
                words[w] |= v << room;
            }
            else
            {
                words[w] |= v >> -room;
                words[w + 1] |= v << (64 + room);
            }
        }
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    internal static void ThrowDamaged() => throw new InvalidDataException("the child table or the LCP array is inconsistent (the l-indices of an interval must increase and stay inside it, a child must be deeper than its parent, and every rank must lie inside its section): the file is damaged");
}
