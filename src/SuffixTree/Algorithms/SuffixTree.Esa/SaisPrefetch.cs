namespace SuffixTree.Esa;

/// <summary>
/// Software prefetch distances of the SA-IS kernel (iteration 1o.2), in loop iterations. The induction passes and the naming read the text and the types at
/// positions that come from the suffix array, i.e. at random; the position that iteration <c>i + distance</c> will need is already known when iteration <c>i</c> runs
/// (it is an entry of the suffix array that is final by then), so its cache line is requested early. The distances come from a measured sweep (PROGRESS-S1.md, 1o.2):
/// the gain is flat from 16 to 48 iterations (an iteration takes 6-8 ns, an L3 hit about 30-60 ns), the values are in the middle of the plateau. The sites where the sweep
/// showed nothing (compaction of the sorted LMS positions, placement of the sorted LMS suffixes, the gather of the reduced ranks) are not prefetched.
/// </summary>
internal static class SaisPrefetch
{
    /// <summary>Induction passes (L and S, both rounds).</summary>
    internal const int Induce = 32;

    /// <summary>Naming of the LMS substrings.</summary>
    internal const int Name = 16;
}
