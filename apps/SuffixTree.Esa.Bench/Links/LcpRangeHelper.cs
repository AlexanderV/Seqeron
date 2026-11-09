#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// The RMQ + PSV/NSV helper of link kernel (c) (ADR-001 section 5, prototype P1): a pyramid of block minima over the LCP array (fanout F, a power of two; level k holds the
/// minimum of F^k consecutive entries) answers range-minimum queries with the leftmost minimum position and the previous / next strictly smaller value in O(F * levels)
/// sequential reads, with 4 / (F - 1) bytes per entry. The LCP entries themselves are read through the <see cref="CompactLcp"/> of the index.
/// </summary>
internal sealed class LcpRangeHelper
{
    internal LcpRangeHelper(CompactLcp lcp, int fanout = 16) => throw new NotImplementedException();

    /// <summary>Number of LCP entries (n + 1).</summary>
    internal int Length => throw new NotImplementedException();

    /// <summary>The fanout F.</summary>
    internal int Fanout => throw new NotImplementedException();

    /// <summary>Number of pyramid levels above the LCP array.</summary>
    internal int Levels => throw new NotImplementedException();

    /// <summary>Bytes of the pyramid (the LCP array is not counted).</summary>
    internal long PayloadBytes => throw new NotImplementedException();

    /// <summary>Minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    internal int Min(int lo, int hi) => throw new NotImplementedException();

    /// <summary>Leftmost position of the minimum of <c>LCP[lo..hi]</c> (inclusive).</summary>
    internal int ArgMin(int lo, int hi) => throw new NotImplementedException();

    /// <summary>The largest <c>j &lt;= pos</c> with <c>LCP[j] &lt; bound</c>, or -1.</summary>
    internal int PrevLess(int pos, int bound) => throw new NotImplementedException();

    /// <summary>The smallest <c>j &gt;= pos</c> with <c>LCP[j] &lt; bound</c>, or <see cref="Length"/>.</summary>
    internal int NextLess(int pos, int bound) => throw new NotImplementedException();

    /// <summary>Previous smaller value (strict): the largest <c>j &lt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or -1.</summary>
    internal int PreviousSmaller(int p) => throw new NotImplementedException();

    /// <summary>Next smaller value (strict): the smallest <c>j &gt; p</c> with <c>LCP[j] &lt; LCP[p]</c>, or <see cref="Length"/>.</summary>
    internal int NextSmaller(int p) => throw new NotImplementedException();

    /// <summary><see cref="PrevLess"/> with step accounting (one step per entry read).</summary>
    internal int PrevLess<TCnt>(int pos, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter => throw new NotImplementedException();

    /// <summary><see cref="NextLess"/> with step accounting.</summary>
    internal int NextLess<TCnt>(int pos, int bound, ref TCnt cnt)
        where TCnt : struct, IStepCounter => throw new NotImplementedException();

    /// <summary><see cref="ArgMin"/> with step accounting.</summary>
    internal int ArgMin<TCnt>(int lo, int hi, ref TCnt cnt)
        where TCnt : struct, IStepCounter => throw new NotImplementedException();
}
