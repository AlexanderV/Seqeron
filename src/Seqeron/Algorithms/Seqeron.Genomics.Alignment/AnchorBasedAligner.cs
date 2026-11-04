using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Anchor-based (seed-chain-extend) pairwise aligner: exact-match anchors between a reference and a
/// query are taken from a suffix tree, a colinear chain of anchors is selected, and the regions
/// between the chain's anchors are aligned with the canonical Needleman-Wunsch dynamic program.
/// <para>
/// The anchor/chain/fill scheme is the one of MUMmer (Delcher et al. 1999, Nucleic Acids Res
/// 27:2369: exact matches from a suffix tree, longest colinear subset, DP on the gaps between them)
/// and LAGAN (Brudno et al. 2003). The anchors here are the maximal exact matches (MEMs) of
/// <c>SuffixTree.FindExactMatchAnchors</c>: one MEM per matching diagonal run found by suffix-link
/// streaming (<b>not</b> the full set of MEMs, and not unique matches), so anchors from different
/// diagonals may overlap; the chain step handles overlaps by trimming (see <see cref="ChainAnchors"/>).
/// </para>
/// <para>
/// <b>Optimality:</b> the alignment is forced through the chained anchors, so it is a valid global
/// alignment whose score is <b>at most</b> that of the unconstrained optimum
/// (<c>SequenceAligner.GlobalAlign</c>); it equals the optimum whenever an optimal path contains all
/// chained anchors. Gap scoring is the linear model of <c>GlobalAlign</c>
/// (<see cref="ScoringMatrix.GapExtend"/> per gap position, <see cref="ScoringMatrix.GapOpen"/> unused),
/// so the reported <c>Score</c> is the linear score of the returned alignment.
/// </para>
/// <para>
/// <b>Cost:</b> O(n+m) suffix-tree streaming + O(a²) chaining over a anchors + O(Σδᵢ²) DP over the
/// inter-anchor gaps (very large gaps use the exact linear-space Hirschberg variant).
/// </para>
/// </summary>
internal static class AnchorBasedAligner
{
    /// <summary>
    /// Minimum anchor length to consider. Shorter matches are likely noise
    /// and don't provide reliable alignment anchors.
    /// For DNA, k=8 gives specificity of 4^8 = 65,536.
    /// </summary>
    private const int DefaultMinAnchorLength = 8;

    /// <summary>
    /// Largest inter-anchor segment (cells = (m+1)·(n+1)) aligned with the full O(m·n)-memory
    /// matrix; bigger segments use the exact linear-space Hirschberg recursion.
    /// </summary>
    private const long MaxGapCellsForFullMatrix = 25_000_000;

    /// <summary>
    /// Represents an exact-match anchor between reference and query.
    /// </summary>
    internal readonly record struct Anchor(
        int RefStart,
        int QueryStart,
        int Length) : IComparable<Anchor>
    {
        public int RefEnd => RefStart + Length;
        public int QueryEnd => QueryStart + Length;

        public int CompareTo(Anchor other)
        {
            int cmp = RefStart.CompareTo(other.RefStart);
            return cmp != 0 ? cmp : QueryStart.CompareTo(other.QueryStart);
        }

        public static bool operator <(Anchor left, Anchor right) => left.CompareTo(right) < 0;
        public static bool operator >(Anchor left, Anchor right) => left.CompareTo(right) > 0;
        public static bool operator <=(Anchor left, Anchor right) => left.CompareTo(right) <= 0;
        public static bool operator >=(Anchor left, Anchor right) => left.CompareTo(right) >= 0;
    }

    /// <summary>
    /// Performs anchor-based pairwise alignment of a query against a reference
    /// that has a pre-built suffix tree.
    /// </summary>
    /// <param name="reference">The reference sequence string (uppercase).</param>
    /// <param name="refTree">Pre-built suffix tree for the reference.</param>
    /// <param name="query">The query sequence string (uppercase).</param>
    /// <param name="scoring">Scoring matrix for NW on gaps.</param>
    /// <param name="minAnchorLength">Minimum anchor length (default: 8).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    internal static AlignmentResult AlignWithAnchors(
        string reference,
        SuffixTree.SuffixTree refTree,
        string query,
        ScoringMatrix scoring,
        int minAnchorLength = DefaultMinAnchorLength)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);
        ArgumentNullException.ThrowIfNull(scoring);

        if (reference.Length == 0 || query.Length == 0)
            return StitchAlignment(reference, query, new List<Anchor>(), scoring);

        // Step 1: maximal exact matches (one MEM per matching diagonal run) by O(n+m) suffix-link streaming
        var memTuples = refTree.FindExactMatchAnchors(query, Math.Max(1, minAnchorLength));
        var mems = new List<Anchor>(memTuples.Count);
        foreach (var (refPos, queryPos, len) in memTuples)
            mems.Add(new Anchor(refPos, queryPos, len));

        // Step 2: weighted colinear chain of the (possibly overlapping) anchors
        var chain = ChainAnchors(mems);

        // Step 3: no anchors -> the gap-filler aligns the whole pair (canonical NW)

        // Step 4: Align gaps between anchors using NW, stitch everything together
        return StitchAlignment(reference, query, chain, scoring);
    }

    /// <summary>
    /// Selects the colinear chain of anchors with the largest total number of matched bases
    /// (maximum-weight chain of colinear non-overlapping fragments, weight = length;
    /// Abouelhoda &amp; Ohlebusch 2005, J Discrete Algorithms 3:321-341; MUMmer's longest increasing
    /// subsequence of matches, Delcher et al. 1999, is the unit-weight, overlap-free special case).
    /// </summary>
    /// <remarks>
    /// Anchors from different diagonals may overlap in either sequence (repeats). Instead of
    /// discarding one of two overlapping anchors, the front of the later anchor is trimmed by the
    /// smallest amount t that makes it start at or after the earlier anchor's end in both sequences;
    /// a suffix of an exact match is still an exact match, so the trimmed anchor <c>(RefStart+t,
    /// QueryStart+t, Length−t)</c> remains valid, and the chain gains only the non-overlapping
    /// remainder. The returned anchors are therefore pairwise non-overlapping and strictly increasing
    /// in both coordinates (<c>RefEnd ≤ next.RefStart</c>, <c>QueryEnd ≤ next.QueryStart</c>).
    /// O(a²) dynamic program over a anchors; deterministic (ties: earlier predecessor / earlier anchor).
    /// </remarks>
    internal static List<Anchor> ChainAnchors(List<Anchor> anchors)
    {
        ArgumentNullException.ThrowIfNull(anchors);

        var sorted = anchors
            .Where(a => a.Length > 0 && a.RefStart >= 0 && a.QueryStart >= 0)
            .OrderBy(a => a.RefStart).ThenBy(a => a.QueryStart).ThenByDescending(a => a.Length)
            .ToList();
        int n = sorted.Count;
        if (n == 0)
            return new List<Anchor>();

        // best[i] = most matched bases of a chain ending with (a front-trimmed) sorted[i];
        // trimOf[i] = how much of sorted[i]'s front the chosen predecessor forces us to drop.
        var best = new long[n];
        var prev = new int[n];
        var trimOf = new int[n];

        for (int i = 0; i < n; i++)
        {
            best[i] = sorted[i].Length;
            prev[i] = -1;

            for (int j = 0; j < i; j++)
            {
                int trim = Math.Max(0, Math.Max(
                    sorted[j].RefEnd - sorted[i].RefStart,
                    sorted[j].QueryEnd - sorted[i].QueryStart));
                int gain = sorted[i].Length - trim;
                if (gain <= 0)
                    continue;

                if (best[j] + gain > best[i])
                {
                    best[i] = best[j] + gain;
                    prev[i] = j;
                    trimOf[i] = trim;
                }
            }
        }

        int bestIdx = 0;
        for (int i = 1; i < n; i++)
            if (best[i] > best[bestIdx])
                bestIdx = i;

        var chain = new List<Anchor>();
        for (int idx = bestIdx; idx >= 0; idx = prev[idx])
        {
            var a = sorted[idx];
            int t = prev[idx] >= 0 ? trimOf[idx] : 0;
            chain.Add(new Anchor(a.RefStart + t, a.QueryStart + t, a.Length - t));
        }

        chain.Reverse();
        return chain;
    }

    /// <summary>
    /// Stitches together the final alignment by:
    /// 1. Aligning the gap before the first anchor
    /// 2. For each anchor: emit the exact match
    /// 3. For each gap between anchors: align with NW
    /// 4. Aligning the gap after the last anchor
    /// </summary>
    private static AlignmentResult StitchAlignment(
        string reference,
        string query,
        List<Anchor> chain,
        ScoringMatrix scoring)
    {
        var aligned1 = new StringBuilder();
        var aligned2 = new StringBuilder();
        int totalScore = 0;

        int refPos = 0;
        int queryPos = 0;

        foreach (var anchor in chain)
        {
            // Align the gap before this anchor
            string refGap = reference[refPos..anchor.RefStart];
            string queryGap = query[queryPos..anchor.QueryStart];

            var (gapAligned1, gapAligned2, gapScore) = AlignGap(refGap, queryGap, scoring);
            aligned1.Append(gapAligned1);
            aligned2.Append(gapAligned2);
            totalScore += gapScore;

            // Emit the anchor (exact match — no DP needed)
            string anchorText = reference[anchor.RefStart..anchor.RefEnd];
            aligned1.Append(anchorText);
            aligned2.Append(anchorText);
            totalScore += anchor.Length * scoring.Match;

            refPos = anchor.RefEnd;
            queryPos = anchor.QueryEnd;
        }

        // Align the trailing gap after the last anchor
        string trailingRef = reference[refPos..];
        string trailingQuery = query[queryPos..];

        var (trailAligned1, trailAligned2, trailScore) = AlignGap(trailingRef, trailingQuery, scoring);
        aligned1.Append(trailAligned1);
        aligned2.Append(trailAligned2);
        totalScore += trailScore;

        return new AlignmentResult(
            AlignedSequence1: aligned1.ToString(),
            AlignedSequence2: aligned2.ToString(),
            Score: totalScore,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0,
            StartPosition2: 0,
            EndPosition1: reference.Length - 1,
            EndPosition2: query.Length - 1);
    }

    /// <summary>
    /// Aligns one inter-anchor segment with the linear-gap model of <c>GlobalAlign</c>
    /// (<see cref="ScoringMatrix.GapExtend"/> per gap position). An empty side is an all-gap run;
    /// identical segments are the trivial diagonal; everything else goes through the canonical
    /// Needleman-Wunsch (linear-space Hirschberg variant above <see cref="MaxGapCellsForFullMatrix"/>
    /// cells), so the gap alignment is always optimal.
    /// </summary>
    private static (string Aligned1, string Aligned2, int Score) AlignGap(
        string refGap, string queryGap, ScoringMatrix scoring)
    {
        if (refGap.Length == 0 && queryGap.Length == 0)
            return ("", "", 0);

        if (refGap.Length == 0)
            return (new string('-', queryGap.Length), queryGap, queryGap.Length * scoring.GapExtend);

        if (queryGap.Length == 0)
            return (refGap, new string('-', refGap.Length), refGap.Length * scoring.GapExtend);

        if (refGap == queryGap)
            return (refGap, queryGap, refGap.Length * scoring.Match);

        var result = SequenceAligner.GlobalAlignLinearSpace(refGap, queryGap, scoring, MaxGapCellsForFullMatrix);
        return (result.AlignedSequence1, result.AlignedSequence2, result.Score);
    }
}
