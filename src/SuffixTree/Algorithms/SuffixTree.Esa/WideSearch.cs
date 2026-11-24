using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa;

/// <summary>The text of a container for the long descent: the dense symbol at a position inside the text (the reads of an index outside it are the caller's: the descent clamps every position to the text).</summary>
internal interface IWideText
{
    /// <summary>The dense symbol at <paramref name="position"/> (<c>0 &lt;= position &lt; n</c>).</summary>
    int Symbol(long position);
}

/// <summary>The pattern of the long descent: dense symbols at a 64-bit index (a pattern of a test can be longer than a span).</summary>
internal interface IWidePattern
{
    /// <summary>Number of symbols.</summary>
    long Length { get; }

    /// <summary>The dense symbol at <paramref name="index"/> (<c>0 &lt;= index &lt; Length</c>).</summary>
    int this[long index] { get; }
}

/// <summary>The child table and the LCP array as the long descent reads them: the first l-index of an interval, the l-index after one, the LCP value at an l-index.</summary>
internal interface IWideCells
{
    /// <summary>The first l-index of the interval <c>[lb..rb]</c> (<paramref name="lastChild"/>: the cell is at its left bound, else at its right bound).</summary>
    long First(long lb, long rb, bool lastChild);

    /// <summary>The l-index after <paramref name="k"/> in the interval that ends at <paramref name="rb"/>, or -1.</summary>
    long Next(long k, long rb);

    /// <summary>The LCP value at the l-index <paramref name="first"/>: the depth of the interval whose first l-index it is.</summary>
    long Depth(long first);
}

/// <summary>The packed TEXT section as the long descent reads it: the 16-bit load and two shifts of the 32-bit descent (<see cref="MappedSearch.Symbol"/>), at a 64-bit position.</summary>
internal readonly unsafe struct PackedSymbols(byte* text, int bits) : IWideText
{
    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Symbol(long position) => MappedSearch.Symbol(text, bits, position);
}

/// <summary>A pattern of dense symbols in a span (the pattern of every call of the index).</summary>
internal readonly ref struct SpanPattern : IWidePattern
{
    private readonly ReadOnlySpan<byte> _symbols;

    /// <summary>The pattern of <paramref name="symbols"/>.</summary>
    internal SpanPattern(ReadOnlySpan<byte> symbols) => _symbols = symbols;

    /// <inheritdoc />
    public long Length => _symbols.Length;

    /// <inheritdoc />
    public int this[long index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if ((ulong)index >= (ulong)_symbols.Length) ThrowIndex(index);
            return Unsafe.Add(ref MemoryMarshal.GetReference(_symbols), (nint)index);
        }
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void ThrowIndex(long index) => throw new ArgumentOutOfRangeException(nameof(index), index, "the pattern is read outside its symbols");
}

/// <summary>
/// The child table (packed CHILD codec 1 over any stored width, tagged codec 2 over 4 bytes) and the LCP array (a base array of 1, 2, 4 or 8 bytes per entry whose escape, if it has one, is resolved through the overflow list of a
/// <see cref="CompactLcp"/>) of a container, at 64-bit ranks. The reads are those of <see cref="PackedCells{TAccess, TLcp}"/> and <see cref="TaggedCells"/> of the 32-bit descent, with the same checks in Checked mode: a rank
/// that a cell sent the walk to must lie in <c>0..n</c>, and a 32- or 64-bit LCP value must not exceed the text. The type arguments are the strategies that the JIT specialises: the access mode, the stored width of the cells and
/// the kind of the LCP base array (its escape is a constant of the instantiation).
/// </summary>
internal readonly struct WideCells<TAccess, TPos, TLcp>(TPos cells, TLcp lcpBase, CompactLcp escapes, bool tagged, long n) : IWideCells
    where TAccess : struct, IAccessChecks
    where TPos : struct, IPosReader
    where TLcp : struct, IPosReader
{
    private const ulong TaggedValueMask = (1UL << 30) - 1;
    private const ulong TaggedKindMask = 3UL << 30;
    private const ulong TaggedNext = 1UL << 30;

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long First(long lb, long rb, bool lastChild)
    {
        ulong cell = cells.Get(lastChild ? lb : rb);
        return tagged ? (long)(cell & TaggedValueMask) : (long)cell;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Next(long k, long rb)
    {
        if (tagged)
        {
            ulong cell = cells.Get(k);
            return (cell & TaggedKindMask) == TaggedNext ? (long)(cell & TaggedValueMask) : -1;
        }
        if (k >= rb) return -1;
        long lk = LcpAt(k);
        if (LcpAt(k + 1) == lk) return k + 1;
        ulong v = cells.Get(k);
        if (TAccess.IsChecked && v > (ulong)n) MappedSearch.ThrowDamaged();
        return LcpAt((long)v) == lk ? (long)v : -1;
    }

    /// <inheritdoc />
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Depth(long first) => LcpAt(first);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private long LcpAt(long rank)
    {
        ulong v = lcpBase.Get(rank);
        if (TLcp.ElementBytes == 1 ? v == byte.MaxValue : TLcp.ElementBytes == 2 && v == ushort.MaxValue) return escapes[rank];
        if (TAccess.IsChecked && TLcp.ElementBytes >= 4 && v > (ulong)(n + 1)) MappedSearch.ThrowDamaged(); // as the 32-bit reader: a value above the length of the array
        return (long)v;
    }
}

/// <summary>
/// The pattern descent at 64-bit ranks, positions and lengths over any stored width (ADR-005 section 3.3, work package S2.4): the loop of <see cref="MappedSearch"/> with <c>long</c> in every register, the suffix array a reader of
/// <see cref="IPosReader"/> (4 or 5 bytes per entry, or a closed form in a test), the child table and the LCP array a <see cref="IWideCells"/>, the text an <see cref="IWideText"/> and the pattern an <see cref="IWidePattern"/>. The
/// checks of <see cref="EsaAccessMode.Checked"/> mode are those of the 32-bit loop: every rank that the walk takes from the child table must increase and stay inside its interval, a child must be deeper than its parent and not deeper
/// than the text, a rank that came from a section is validated before it is used as an index, and a text position that came from the suffix array is clamped to the text. The text is compared symbol by symbol (the 32-bit loop compares
/// 64-bit windows: a text of more than 2^31 - 2 symbols is not the hot path of S1 and the packed compare is not worth a second packed-pattern buffer here).
/// </summary>
#pragma warning disable IDE0060, S2326 // the stub of the commit of the tests: the loop is written in the next commit
internal static class WideSearch
{
    /// <summary>
    /// The suffix-array interval <c>[lo..hi]</c> of the NON-EMPTY pattern (every symbol a dense symbol of the text); false, with the empty interval <c>[0..-1]</c>, when it does not occur.
    /// </summary>
    /// <param name="sa">The suffix array: <c>n + 1</c> entries.</param>
    /// <param name="cells">The child table and the LCP array.</param>
    /// <param name="text">The text.</param>
    /// <param name="n">The length of the text.</param>
    /// <param name="pat">The pattern.</param>
    /// <param name="lo">The first rank.</param>
    /// <param name="hi">The last rank.</param>
    /// <exception cref="InvalidDataException">Checked mode: a table of the file is inconsistent.</exception>
    internal static bool Run<TAccess, TSa, TCells, TText, TPat>(in TSa sa, in TCells cells, in TText text, long n, in TPat pat, out long lo, out long hi)
        where TAccess : struct, IAccessChecks
        where TSa : struct, IPosReader
        where TCells : struct, IWideCells
        where TText : struct, IWideText
        where TPat : struct, IWidePattern, allows ref struct
    {
        // the commit of the tests (S2.4, step b): the loop is written in the next commit; every test of the long descent is red until then
        throw new NotImplementedException("the long descent (work package S2.4)");
    }

    /// <summary>Whether the text symbols <c>[start, start + length)</c> (inside the text) equal the pattern symbols <c>[from, from + length)</c>.</summary>
    internal static bool Equal<TText, TPat>(in TText text, long start, in TPat pat, long from, long length)
        where TText : struct, IWideText
        where TPat : struct, IWidePattern, allows ref struct
    {
        for (long i = 0; i < length; i++)
        {
            if (text.Symbol(start + i) != pat[from + i]) return false;
        }
        return true;
    }
}

#pragma warning restore IDE0060, S2326

/// <summary>The sections of a container (or of the memory of a test) that the long descent needs, as pointers and counts.</summary>
internal readonly unsafe struct WideSections
{
    /// <summary>The sections.</summary>
    /// <param name="isChecked">Checked access (else Trusted).</param>
    /// <param name="posWidth">The stored width of a suffix-array entry and of a child cell: 4 or 5 bytes.</param>
    /// <param name="lcpKind">The kind of the LCP base array.</param>
    /// <param name="tagged">CHILD codec 2 (4-byte cells only).</param>
    /// <param name="sa">The suffix array: <c>n + 1</c> entries of <paramref name="posWidth"/> bytes.</param>
    /// <param name="child">The child table: <c>n + 1</c> cells of <paramref name="posWidth"/> bytes.</param>
    /// <param name="lcpBase">The LCP base array: <c>n + 1</c> entries of the width of <paramref name="lcpKind"/>.</param>
    /// <param name="lcp">The LCP array (its overflow list resolves the escapes).</param>
    /// <param name="text">The packed text.</param>
    /// <param name="bits">Bits per symbol of the text.</param>
    /// <param name="n">The length of the text.</param>
    internal WideSections(bool isChecked, int posWidth, LcpKind lcpKind, bool tagged, byte* sa, byte* child, byte* lcpBase, CompactLcp lcp, byte* text, int bits, long n)
    {
        IsChecked = isChecked;
        PosWidth = posWidth;
        LcpKind = lcpKind;
        Tagged = tagged;
        Sa = sa;
        Child = child;
        LcpBase = lcpBase;
        Lcp = lcp;
        Text = text;
        Bits = bits;
        N = n;
    }

    internal bool IsChecked { get; }

    internal int PosWidth { get; }

    internal LcpKind LcpKind { get; }

    internal bool Tagged { get; }

    internal byte* Sa { get; }

    internal byte* Child { get; }

    internal byte* LcpBase { get; }

    internal CompactLcp Lcp { get; }

    internal byte* Text { get; }

    internal int Bits { get; }

    internal long N { get; }
}

/// <summary>
/// The long descent of one index, with the strategies chosen once (the access mode, the stored width of the suffix array and the cells, the kind of the LCP base array): <see cref="Create"/> is the one switch of the matrix, one
/// instantiation of <see cref="WideSearch.Run{TAccess, TSa, TCells, TText, TPat}"/> per row (2 access modes x 2 widths x 4 LCP kinds = 16; the tagged codec is a flag of <see cref="WideCells{TAccess, TPos, TLcp}"/>, not a type), of which
/// the JIT compiles the ones that a process uses. <see cref="MappedSearch"/> holds one for every container that its 32-bit loop cannot read.
/// </summary>
internal abstract class WideSearcher
{
    /// <summary>The interval of the NON-EMPTY <paramref name="pattern"/> of dense symbols; false when it does not occur.</summary>
    /// <exception cref="InvalidDataException">Checked mode: a table of the file is inconsistent.</exception>
    internal abstract bool TryMatch(ReadOnlySpan<byte> pattern, out long lo, out long hi);

    /// <summary>The searcher of the sections <paramref name="s"/>.</summary>
    /// <exception cref="ArgumentOutOfRangeException">The stored width is not 4 or 5, or a tagged table over 5-byte cells.</exception>
    internal static WideSearcher Create(in WideSections s)
    {
        if (s.PosWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(s), s.PosWidth, "the stored width of a position is 4 or 5 bytes");
        if (s.Tagged && s.PosWidth == 5) throw new ArgumentException("the tagged codec (CHILD codec 2) exists over 4-byte cells only (ADR-005 section 6.1)", nameof(s));
        return s.IsChecked ? Select<CheckedAccess>(s) : Select<TrustedAccess>(s);
    }

    private static unsafe WideSearcher Select<TAccess>(in WideSections s)
        where TAccess : struct, IAccessChecks
    {
        long count = s.N + 1;
        return s.PosWidth == 5
            ? SelectLcp<TAccess, Pos40<TAccess>>(s, new Pos40<TAccess>(s.Sa, count), new Pos40<TAccess>(s.Child, count))
            : SelectLcp<TAccess, Pos32<TAccess>>(s, new Pos32<TAccess>(s.Sa, count), new Pos32<TAccess>(s.Child, count));
    }

    private static unsafe WideSearcher SelectLcp<TAccess, TPos>(in WideSections s, TPos sa, TPos child)
        where TAccess : struct, IAccessChecks
        where TPos : struct, IPosReader
    {
        long count = s.N + 1;
        return s.LcpKind switch
        {
            LcpKind.U8 => Make<TAccess, TPos, Pos8<TAccess>>(s, sa, child, new Pos8<TAccess>(s.LcpBase, count)),
            LcpKind.U16 => Make<TAccess, TPos, Pos16<TAccess>>(s, sa, child, new Pos16<TAccess>(s.LcpBase, count)),
            LcpKind.U32 => Make<TAccess, TPos, Pos32<TAccess>>(s, sa, child, new Pos32<TAccess>(s.LcpBase, count)),
            _ => Make<TAccess, TPos, Pos64<TAccess>>(s, sa, child, new Pos64<TAccess>(s.LcpBase, count)),
        };
    }

    private static unsafe WideSearcher Make<TAccess, TPos, TLcp>(in WideSections s, TPos sa, TPos child, TLcp lcpBase)
        where TAccess : struct, IAccessChecks
        where TPos : struct, IPosReader
        where TLcp : struct, IPosReader =>
        new WideSearcher<TAccess, TPos, TLcp>(sa, new WideCells<TAccess, TPos, TLcp>(child, lcpBase, s.Lcp, s.Tagged, s.N), new PackedSymbols(s.Text, s.Bits), s.N);
}

/// <inheritdoc cref="WideSearcher"/>
internal sealed class WideSearcher<TAccess, TPos, TLcp>(TPos sa, WideCells<TAccess, TPos, TLcp> cells, PackedSymbols text, long n) : WideSearcher
    where TAccess : struct, IAccessChecks
    where TPos : struct, IPosReader
    where TLcp : struct, IPosReader
{
    /// <inheritdoc />
    internal override bool TryMatch(ReadOnlySpan<byte> pattern, out long lo, out long hi) =>
        WideSearch.Run<TAccess, TPos, WideCells<TAccess, TPos, TLcp>, PackedSymbols, SpanPattern>(sa, cells, text, n, new SpanPattern(pattern), out lo, out hi);
}
