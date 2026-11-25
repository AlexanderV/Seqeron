using System.Text;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// An inclusive interval of suffix-array ranks <c>[First..Last]</c>; <see cref="Empty"/> (<c>[0..-1]</c>) when nothing matches. Rank 0 is the empty suffix, so the ranks of a non-empty pattern are in <c>[1..n]</c>.
/// </summary>
internal readonly record struct RankRange(long First, long Last)
{
    /// <summary>The range of no rank.</summary>
    internal static RankRange Empty { get; } = new(0, -1);

    /// <summary>Number of ranks.</summary>
    internal long Count => Last < First ? 0 : Last - First + 1;

    /// <summary>True when no rank matches.</summary>
    internal bool IsEmpty => Last < First;

    /// <summary>Two empty ranges are equal whatever their bounds.</summary>
    internal bool SameAs(RankRange other) => IsEmpty ? other.IsEmpty : First == other.First && Last == other.Last;

    public override string ToString() => IsEmpty ? "[empty]" : $"[{First}..{Last}]";
}

/// <summary>
/// A text whose suffix array, LCP array and child table have a closed form (work package S2.1, ADR-005 section 3.7 layer 1): a pure function of <c>long</c> arguments that answers every question of the index over
/// a text of any length (2^31, 2^32, 2^33, 2^40 symbols and more) without memory, so that the 64-bit engine can be tested against it where no real text exists. Layout of ADR-002 section 2: ranks <c>0..n</c>,
/// <c>SuffixAt(0) = n</c> (the empty suffix), <c>Lcp(0) = Lcp(1) = 0</c>, the end of the text smaller than every symbol. Every argument and every result is <c>long</c> (a count of internal nodes, a depth, a rank, a position);
/// the only <c>Int128</c> is <see cref="LcpSum"/>, which is of the order of n^2.
/// <para>
/// The child table is given at the semantic level of the reference layout (<c>Reference3Child</c>), the values and not a codec: <see cref="Down"/>, <see cref="Up"/> and <see cref="Next"/> are the three arrays (0 = undefined),
/// <see cref="FirstLIndex"/> and <see cref="NextLIndex"/> are the two reads of the interface <c>IChildTable</c>; <see cref="PackedCell"/> derives the cell of CHILD codec 1 from them.
/// </para>
/// </summary>
public abstract class ClosedFormText
{
    protected ClosedFormText(long length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        Length = length;
    }

    /// <summary>The text length n.</summary>
    internal long Length { get; }

    /// <summary>A short name for messages.</summary>
    internal abstract string Name { get; }

    /// <summary>The period of the infinite word of which the text is a prefix (1 for <c>a^n</c>).</summary>
    internal abstract long Period { get; }

    /// <summary>Whether <see cref="Down"/>, <see cref="Up"/>, <see cref="Next"/> (and <see cref="PackedCell"/>) have a closed form; false for <c>u^m</c> of a de Bruijn word, whose child table is the engine's own scan.</summary>
    internal virtual bool HasChildTable => true;

    /// <summary>The number of distinct letters in the text (0 for the empty text).</summary>
    internal abstract int Sigma { get; }

    /// <summary>How many times the symbol of <paramref name="code"/> occurs in the text.</summary>
    internal abstract long SymbolCount(int code);

    /// <summary>The dense code (0 = 'a', 1 = 'b') of the infinite word at <paramref name="position"/> (any non-negative position, the text length does not matter).</summary>
    internal abstract int InfiniteCodeAt(long position);

    /// <summary>The dense code of the symbol at <paramref name="position"/> in <c>[0, n)</c>.</summary>
    internal int CodeAt(long position)
    {
        if ((ulong)position >= (ulong)Length) throw new ArgumentOutOfRangeException(nameof(position), position, $"the text has {Length} symbols");
        return InfiniteCodeAt(position);
    }

    /// <summary>The letter of a dense code.</summary>
    internal static char LetterOf(int code) => (char)('a' + code);

    /// <summary>The symbol at <paramref name="position"/>.</summary>
    internal char CharAt(long position) => LetterOf(CodeAt(position));

    /// <summary><c>SA[rank]</c> for <c>rank</c> in <c>[0, n]</c> (<c>SA[0] = n</c>).</summary>
    internal abstract long SuffixAt(long rank);

    /// <summary><c>ISA[position]</c> for <c>position</c> in <c>[0, n]</c> (<c>ISA[n] = 0</c>).</summary>
    internal abstract long RankOf(long position);

    /// <summary><c>LCP[rank]</c> for <c>rank</c> in <c>[0, n]</c>: the common prefix of the suffixes at the ranks <c>rank - 1</c> and <c>rank</c>; <c>LCP[0] = LCP[1] = 0</c>.</summary>
    internal abstract long Lcp(long rank);

    /// <summary>The sum of <c>LCP[r]</c> over every rank (about n^2 / 2, so wider than <c>long</c>).</summary>
    internal abstract Int128 LcpSum { get; }

    /// <summary>
    /// <c>min LCP[lo..hi]</c> (inclusive; the range minimum that the LCP pyramid answers). The default reads the range (at most 2^20 ranks); a family with a closed form for it overrides it, so that the levels of the pyramid of a text of
    /// 2^32 symbols can be written without reading its 2^32 LCP values (work package S2.5).
    /// </summary>
    internal virtual long LcpMin(long lo, long hi)
    {
        if (lo < 0 || hi < lo || hi > Length) throw new ArgumentOutOfRangeException(nameof(lo), lo, $"[{lo}, {hi}] is not a range of {Length + 1} ranks");
        if (hi - lo > 1 << 20) throw new NotSupportedException($"{Name} has no closed form for the minimum of {hi - lo + 1} LCP values");
        long min = long.MaxValue;
        for (long r = lo; r <= hi; r++) min = Math.Min(min, Lcp(r));
        return min;
    }

    /// <summary>The number of internal nodes other than the root (the lcp-intervals of value above 0).</summary>
    internal abstract long InternalCount { get; }

    /// <summary><c>NODE_COUNT</c> of the terminated tree: <c>1 + InternalCount + (n + 1)</c>, 1 for the empty text.</summary>
    internal long NodeCount => Length == 0 ? 1 : 1 + InternalCount + Length + 1;

    /// <summary>The deepest leaf (<c>EsaStatistics.MaxDepth</c>): n.</summary>
    internal long MaxDepth => Length;

    /// <summary>Length of the longest repeated substring (0 when no substring occurs twice).</summary>
    internal abstract long LrsDepth { get; }

    /// <summary>The header form of <c>LRS_RANK</c>: the left bound of the lcp-interval of value <see cref="LrsDepth"/> that closes first, -1 when there is no repeat.</summary>
    internal abstract long LrsRank { get; }

    /// <summary>The start position of the longest repeated substring as the header names it (<c>SA[LRS_RANK]</c>), -1 when there is none.</summary>
    internal long LrsPosition => LrsDepth == 0 ? -1 : SuffixAt(LrsRank);

    /// <summary>The ranks of the suffixes that start with <paramref name="pattern"/> (letters other than 'a' and 'b' occur nowhere); the empty pattern is <c>[1..n]</c> (count n, as <c>CountOccurrences("")</c>).</summary>
    internal abstract RankRange RangeOf(ReadOnlySpan<char> pattern);

    /// <summary>The ranks of the suffixes that start with the pattern <c>a^length</c> (a pattern too long for a string).</summary>
    internal abstract RankRange RangeOfRun(long length);

    /// <summary>The ranks of the suffixes that start with <c>text[start .. start + length)</c>, the pattern being a factor of the text itself (a pattern too long for a string).</summary>
    internal abstract RankRange RangeOfFactor(long start, long length);

    /// <summary>The number of occurrences of <paramref name="pattern"/>.</summary>
    internal long CountOf(ReadOnlySpan<char> pattern) => RangeOf(pattern).Count;

    /// <summary><c>Down[i]</c> of the reference layout: the first l-index of the last child interval that starts at rank <c>i</c> (and of the root at 0); 0 = none.</summary>
    internal abstract long Down(long i);

    /// <summary><c>Up[j]</c>: the first l-index of the non-last child interval that ends at rank <c>j</c>; 0 = none.</summary>
    internal abstract long Up(long j);

    /// <summary><c>Next[k]</c>: the l-index after <c>k</c> in its interval; 0 = none.</summary>
    internal abstract long Next(long k);

    /// <summary><c>IChildTable.FirstLIndex</c>: the first l-index of the internal interval <c>[lb..rb]</c>, the field sitting at <c>lb</c> for a last child and at <c>rb</c> otherwise.</summary>
    internal long FirstLIndex(long lb, long rb, bool lastChild) => lastChild ? Down(lb) : Up(rb);

    /// <summary><c>IChildTable.NextLIndex</c>: the next l-index of <c>k</c>, -1 when <c>k</c> is the last.</summary>
    internal long NextLIndex(long k, long rb)
    {
        long v = Next(k);
        return v == 0 ? -1 : v;
    }

    /// <summary>
    /// The cell at <paramref name="rank"/> of CHILD codec 1 (the packed layout: one value per rank, the kind told by comparing LCP values): the next l-index when it is not <c>rank + 1</c>, else the up field of the
    /// non-last child that ends here, else the down field of the last child that starts here (the fields never meet). 0 = empty.
    /// </summary>
    internal long PackedCell(long rank)
    {
        long next = Next(rank);
        if (next != 0 && next != rank + 1) return next;
        long up = Up(rank);
        return up != 0 ? up : Down(rank);
    }

    /// <summary>The text as a string, for the oracles that need one (refused above 2^26 symbols).</summary>
    internal string Materialize()
    {
        if (Length > 1L << 26) throw new InvalidOperationException($"a text of {Length} symbols is not materialized");
        var sb = new StringBuilder((int)Length);
        for (long i = 0; i < Length; i++) sb.Append(CharAt(i));
        return sb.ToString();
    }

    /// <summary>Throws for a rank outside <c>[0, n]</c>.</summary>
    protected void CheckRank(long rank, string name = "rank")
    {
        if ((ulong)rank > (ulong)Length) throw new ArgumentOutOfRangeException(name, rank, $"a rank is in [0, {Length}]");
    }

    public override string ToString() => $"{Name} (n = {Length})";
}
