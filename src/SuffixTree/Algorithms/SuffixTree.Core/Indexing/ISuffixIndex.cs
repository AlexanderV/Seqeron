using System;

namespace SuffixTree.Indexing;

/// <summary>
/// The engine-neutral core of a text index (ADR-001 section 6): count and extend a pattern symbol by symbol, locate lazily, read the text.
/// Symbols are UTF-16 code units <c>0..65535</c>; the terminator is -1 and is not a symbol of the text. Implementations are safe for concurrent
/// readers; an operation an engine cannot do (<see cref="Support"/> reports <see cref="IndexSupport.Unsupported"/>) throws
/// <see cref="NotSupportedException"/> naming the engine and the operation and never returns a wrong or empty answer.
/// </summary>
public interface ISuffixIndex
{
    /// <summary>Engine name: "esa", later "fm", "r-index", "ukkonen".</summary>
    string Engine { get; }

    /// <summary>Number of symbols n of the indexed text.</summary>
    long TextLength { get; }

    /// <summary>The indexed text.</summary>
    IIndexText Text { get; }

    /// <summary>The locus of the empty string; <c>Count(Root) == TextLength</c>.</summary>
    IndexLocus Root { get; }

    /// <summary>How well this engine answers <paramref name="operation"/>.</summary>
    IndexSupport Support(IndexOperation operation);

    /// <summary>
    /// A short statement of the cost of <paramref name="operation"/> (for example "O(width), builds ISA on first use"); never null for
    /// <see cref="IndexSupport.Degraded"/>, may be null otherwise.
    /// </summary>
    string? SupportNote(IndexOperation operation);

    /// <summary>Number of start positions of the string of <paramref name="locus"/> in the text (overlapping occurrences included).</summary>
    long Count(in IndexLocus locus);

    /// <summary>
    /// Appends <paramref name="symbol"/> to the locus string. Returns false, and leaves <paramref name="result"/> default, when the extended string
    /// does not occur in the text or <paramref name="symbol"/> is not a code unit (outside 0..65535, the terminator included).
    /// </summary>
    bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result);

    /// <summary>
    /// Appends <paramref name="pattern"/> to the locus string; an empty pattern returns the locus itself. Returns false, and leaves
    /// <paramref name="result"/> default, when the extended string does not occur in the text.
    /// </summary>
    bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result);

    /// <summary>Prepends <paramref name="symbol"/> to the locus string; same failure convention as <see cref="TryExtendRight(in IndexLocus, int, out IndexLocus)"/>.</summary>
    bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result);

    /// <summary>
    /// Lazy, paged locate: writes up to <paramref name="positions"/>.Length start positions of the locus string, after skipping the first
    /// <paramref name="skip"/> of them in engine order (suffix-array engines: suffix-array order, not ascending), and returns the number written.
    /// Paging with growing <paramref name="skip"/> enumerates every position exactly once.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="skip"/> is negative.</exception>
    int Locate(in IndexLocus locus, long skip, Span<long> positions);
}
