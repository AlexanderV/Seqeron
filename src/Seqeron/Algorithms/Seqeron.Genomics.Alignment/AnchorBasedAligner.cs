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
/// so the reported <c>Score</c> is the linear score of the returned alignment; with <c>affineGaps</c> the inter-anchor segments are
/// aligned with the exact affine (Gotoh) aligner and a gap of length k costs <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>
/// (anchors separate every segment, so no gap run crosses a segment boundary and <c>Score</c> is the affine score of the returned rows).
/// </para>
/// <para>
/// <b>Cost:</b> O(n+m) suffix-tree streaming + O(a log² a) chaining over a anchors (range-maximum sparse dynamic
/// programming, see <see cref="ChainAnchors"/>) + O(Σδᵢ²) DP over the inter-anchor gaps (very large gaps use the exact linear-space Hirschberg variant).
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
    /// <param name="affineGaps">Align the inter-anchor segments with affine gaps (open + k·extend) instead of the linear model.</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    internal static AlignmentResult AlignWithAnchors(
        string reference,
        SuffixTree.SuffixTree refTree,
        string query,
        PairScoring scoring,
        int minAnchorLength = DefaultMinAnchorLength,
        bool affineGaps = false)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);

        if (reference.Length == 0 || query.Length == 0)
            return StitchAlignment(reference, query, new List<Anchor>(), scoring, affineGaps);

        // Step 1: maximal exact matches (one MEM per matching diagonal run) by O(n+m) suffix-link streaming
        var memTuples = refTree.FindExactMatchAnchors(query, Math.Max(1, minAnchorLength));
        var mems = new List<Anchor>(memTuples.Count);
        foreach (var (refPos, queryPos, len) in memTuples)
            mems.Add(new Anchor(refPos, queryPos, len));

        // Step 2: weighted colinear chain of the (possibly overlapping) anchors
        var chain = ChainAnchors(mems);

        // Step 3: no anchors -> the gap-filler aligns the whole pair (canonical NW)

        // Step 4: Align gaps between anchors using NW, stitch everything together
        return StitchAlignment(reference, query, chain, scoring, affineGaps);
    }

    /// <summary>
    /// Selects the colinear chain of anchors with the largest total number of matched bases
    /// (maximum-weight chain of colinear non-overlapping fragments, weight = length;
    /// Abouelhoda &amp; Ohlebusch 2005, J Discrete Algorithms 3:321-341; MUMmer's longest increasing
    /// subsequence of matches, Delcher et al. 1999, is the unit-weight, overlap-free special case).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Anchors from different diagonals may overlap in either sequence (repeats). Instead of
    /// discarding one of two overlapping anchors, the front of the later anchor is trimmed by the
    /// smallest amount t that makes it start at or after the earlier anchor's end in both sequences;
    /// a suffix of an exact match is still an exact match, so the trimmed anchor <c>(RefStart+t,
    /// QueryStart+t, Length−t)</c> remains valid, and the chain gains only the non-overlapping
    /// remainder. The returned anchors are therefore pairwise non-overlapping and strictly increasing
    /// in both coordinates (<c>RefEnd ≤ next.RefStart</c>, <c>QueryEnd ≤ next.QueryStart</c>).
    /// Deterministic (ties: earlier predecessor / earlier anchor).
    /// </para>
    /// <para>
    /// <b>Algorithm (A1-16):</b> the recurrence is best[i] = Lᵢ + max(0, maxⱼ&lt;ᵢ best[j] − trim(j,i)) with
    /// trim(j,i) = max(0, endRef(j) − refStart(i), endQuery(j) − queryStart(i)) over the predecessors j that leave
    /// a positive remainder (trim &lt; Lᵢ ⇔ endRef(j) &lt; endRef(i) and endQuery(j) &lt; endQuery(i)). It is evaluated by
    /// range-maximum queries instead of the all-pairs scan (the sparse dynamic programming of Abouelhoda &amp;
    /// Ohlebusch): with δ = query − ref the diagonal, the predecessors split by δⱼ ≤ δᵢ (the reference overlap is the
    /// larger one) and δⱼ &gt; δᵢ (the query overlap is), and inside each half trim is a function of ONE end coordinate
    /// (0 when it is ≤ the new anchor's start, else end − start), so each half is two orthogonal range-maximum queries
    /// (diagonal prefix × end-coordinate prefix / interval) answered by a Fenwick tree over diagonal ranks whose nodes
    /// hold max segment trees over the end coordinate. Time O(a log² a) (O(a log a) per query set for anchors without
    /// overlaps), memory O(a log a); the chain, its score and every tie decision are identical to the quadratic
    /// scan (kept as the oracle in the test project).
    /// </para>
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

        var refEnd = new long[n];
        var queryEnd = new long[n];
        var diag = new long[n];
        for (int i = 0; i < n; i++)
        {
            refEnd[i] = sorted[i].RefEnd;
            queryEnd[i] = sorted[i].QueryEnd;
            diag[i] = (long)sorted[i].QueryStart - sorted[i].RefStart;
        }

        // Diagonal ranks (1-based): ascending for the "δj <= δi" half, descending for the "δj > δi" half.
        var distinctDiag = diag.Distinct().OrderBy(d => d).ToArray();
        int nd = distinctDiag.Length;
        int[] rankAsc = diag.Select(d => Array.BinarySearch(distinctDiag, d) + 1).ToArray();
        int[] rankDesc = rankAsc.Select(r => nd - r + 1).ToArray();

        // best[i] = most matched bases of a chain ending with (a front-trimmed) sorted[i];
        // trimOf[i] = how much of sorted[i]'s front the chosen predecessor forces us to drop.
        var best = new long[n];
        var prev = new int[n];
        var trimOf = new int[n];

        // value channels of both range structures: 0 = best[j] (end coordinate not overlapped),
        // 1 = best[j] − end (overlapped: trim = end − start of the new anchor)
        var refValue1 = new long[n];
        var queryValue1 = new long[n];
        var refSide = new RangeMaxFenwick(rankAsc, refEnd, nd, best, refValue1);
        var querySide = new RangeMaxFenwick(rankDesc, queryEnd, nd, best, queryValue1);

        for (int i = 0; i < n; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per anchor of the chaining recurrence
            long len = sorted[i].Length;
            long rs = sorted[i].RefStart;
            long qs = sorted[i].QueryStart;
            best[i] = len;
            prev[i] = -1;

            // A predecessor must leave best[j] − trim > 0 to beat "no predecessor"; among equal gains the smallest j wins.
            long bestGain = 0;
            int bestJ = -1;

            void Consider(int j, long gain)
            {
                if (j >= 0 && gain > 0 && (gain > bestGain || (gain == bestGain && j < bestJ)))
                {
                    bestGain = gain;
                    bestJ = j;
                }
            }

            // diag(j) <= diag(i): the reference overlap is the larger one, trim = max(0, refEnd(j) − refStart(i)).
            int j1 = refSide.QueryAtMost(rankAsc[i], rs, 0);
            Consider(j1, j1 < 0 ? 0 : best[j1]);
            int j2 = refSide.QueryBetween(rankAsc[i], rs, refEnd[i], 1);
            Consider(j2, j2 < 0 ? 0 : refValue1[j2] + rs);

            // diag(j) > diag(i): the query overlap is the larger one, trim = max(0, queryEnd(j) − queryStart(i)).
            int above = nd - rankAsc[i];
            int j3 = querySide.QueryAtMost(above, qs, 0);
            Consider(j3, j3 < 0 ? 0 : best[j3]);
            int j4 = querySide.QueryBetween(above, qs, queryEnd[i], 1);
            Consider(j4, j4 < 0 ? 0 : queryValue1[j4] + qs);

            if (bestJ >= 0)
            {
                best[i] += bestGain;
                prev[i] = bestJ;
                trimOf[i] = (int)Math.Max(0, Math.Max(refEnd[bestJ] - rs, queryEnd[bestJ] - qs));
            }

            refValue1[i] = best[i] - refEnd[i];
            queryValue1[i] = best[i] - queryEnd[i];
            refSide.Insert(i);
            querySide.Insert(i);
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
    /// Offline 2-D range-maximum structure of the chaining DP. Points p = 0..n−1 have a fixed x rank (1-based, prefix
    /// dimension) and a fixed y (interval dimension); their two value channels live in caller-owned arrays that the caller
    /// fills before <see cref="Insert"/>. A query returns the inserted point with the largest channel value (smallest
    /// index on ties) with x rank ≤ X and y in the requested range. A Fenwick tree over x whose node k holds the points
    /// routed through k sorted by (y, p) and one iterative max segment tree per channel over them:
    /// O(log² n) per operation, O(n log n) memory.
    /// </summary>
    private sealed class RangeMaxFenwick
    {
        private readonly int[] _xRank;
        private readonly long[] _y;
        private readonly int _xSize;
        private readonly long[][] _channel;   // [0], [1]: caller-owned value arrays indexed by point id
        private readonly int[][] _ids;        // per Fenwick node: point ids sorted by (y, id)
        private readonly int[][][] _trees;    // [channel][node]: max segment tree of point ids (-1 = empty)

        public RangeMaxFenwick(int[] xRank, long[] y, int xSize, long[] value0, long[] value1)
        {
            _xRank = xRank;
            _y = y;
            _xSize = xSize;
            _channel = new[] { value0, value1 };

            int n = y.Length;
            var count = new int[xSize + 1];
            for (int p = 0; p < n; p++)
                for (int k = xRank[p]; k <= xSize; k += k & -k)
                    count[k]++;

            _ids = new int[xSize + 1][];
            _trees = new[] { new int[xSize + 1][], new int[xSize + 1][] };
            for (int k = 1; k <= xSize; k++)
            {
                _ids[k] = new int[count[k]];
                for (int c = 0; c < 2; c++)
                {
                    _trees[c][k] = new int[2 * count[k]];
                    Array.Fill(_trees[c][k], -1);
                }

                count[k] = 0;
            }

            // Points in (y, id) order, so every node's list comes out sorted.
            foreach (int p in Enumerable.Range(0, n).OrderBy(q => y[q]).ThenBy(q => q))
                for (int k = xRank[p]; k <= xSize; k += k & -k)
                    _ids[k][count[k]++] = p;
        }

        /// <summary>Makes point <paramref name="p"/> visible to queries (its channel values must already be set).</summary>
        public void Insert(int p)
        {
            for (int k = _xRank[p]; k <= _xSize; k += k & -k)
            {
                int m = _ids[k].Length;
                int pos = LowerBound(_ids[k], _y[p], p) + m;
                for (int c = 0; c < 2; c++)
                {
                    var tree = _trees[c][k];
                    tree[pos] = p;
                    for (int i = pos >> 1; i >= 1; i >>= 1)
                        tree[i] = Better(tree[2 * i], tree[2 * i + 1], _channel[c]);
                }
            }
        }

        /// <summary>Best inserted point with x rank ≤ <paramref name="xPrefix"/> and y ≤ <paramref name="yMax"/>.</summary>
        public int QueryAtMost(int xPrefix, long yMax, int channel) =>
            Query(xPrefix, channel, ids => 0, ids => UpperBoundY(ids, yMax));

        /// <summary>Best inserted point with x rank ≤ <paramref name="xPrefix"/> and <paramref name="yLowExclusive"/> &lt; y &lt; <paramref name="yHighExclusive"/>.</summary>
        public int QueryBetween(int xPrefix, long yLowExclusive, long yHighExclusive, int channel) =>
            Query(xPrefix, channel, ids => UpperBoundY(ids, yLowExclusive), ids => LowerBoundY(ids, yHighExclusive));

        private int Query(int xPrefix, int channel, Func<int[], int> from, Func<int[], int> to)
        {
            int result = -1;
            for (int k = Math.Min(xPrefix, _xSize); k > 0; k -= k & -k)
            {
                var ids = _ids[k];
                int l = from(ids), r = to(ids);
                if (l >= r)
                    continue;

                var tree = _trees[channel][k];
                int m = ids.Length;
                int lo = l + m, hi = r + m;
                while (lo < hi)
                {
                    if ((lo & 1) == 1)
                        result = Better(result, tree[lo++], _channel[channel]);
                    if ((hi & 1) == 1)
                        result = Better(result, tree[--hi], _channel[channel]);
                    lo >>= 1;
                    hi >>= 1;
                }
            }

            return result;
        }

        private static int Better(int a, int b, long[] value)
        {
            if (a < 0)
                return b;
            if (b < 0)
                return a;
            return value[a] > value[b] || (value[a] == value[b] && a < b) ? a : b;
        }

        /// <summary>First position whose (y, id) is not below (<paramref name="y"/>, <paramref name="id"/>).</summary>
        private int LowerBound(int[] ids, long y, int id)
        {
            int lo = 0, hi = ids.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                int q = ids[mid];
                if (_y[q] < y || (_y[q] == y && q < id))
                    lo = mid + 1;
                else
                    hi = mid;
            }

            return lo;
        }

        /// <summary>Number of points with y &lt; <paramref name="y"/>.</summary>
        private int LowerBoundY(int[] ids, long y)
        {
            int lo = 0, hi = ids.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (_y[ids[mid]] < y)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            return lo;
        }

        /// <summary>Number of points with y ≤ <paramref name="y"/>.</summary>
        private int UpperBoundY(int[] ids, long y)
        {
            int lo = 0, hi = ids.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (_y[ids[mid]] <= y)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            return lo;
        }
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
        PairScoring scoring,
        bool affineGaps)
    {
        var aligned1 = new StringBuilder();
        var aligned2 = new StringBuilder();
        int totalScore = 0;

        int refPos = 0;
        int queryPos = 0;

        foreach (var anchor in chain)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per anchor (each gap is aligned by the polling linear-space DP)
            // Align the gap before this anchor
            string refGap = reference[refPos..anchor.RefStart];
            string queryGap = query[queryPos..anchor.QueryStart];

            var (gapAligned1, gapAligned2, gapScore) = AlignGap(refGap, queryGap, scoring, affineGaps);
            aligned1.Append(gapAligned1);
            aligned2.Append(gapAligned2);
            totalScore += gapScore;

            // Emit the anchor (exact match — no DP needed)
            string anchorText = reference[anchor.RefStart..anchor.RefEnd];
            aligned1.Append(anchorText);
            aligned2.Append(anchorText);
            totalScore += IdenticalScore(anchorText, scoring);

            refPos = anchor.RefEnd;
            queryPos = anchor.QueryEnd;
        }

        // Align the trailing gap after the last anchor
        string trailingRef = reference[refPos..];
        string trailingQuery = query[queryPos..];

        var (trailAligned1, trailAligned2, trailScore) = AlignGap(trailingRef, trailingQuery, scoring, affineGaps);
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
    /// (<see cref="ScoringMatrix.GapExtend"/> per gap position), or, with <paramref name="affineGaps"/>, the affine model of
    /// <c>GlobalAlignAffine</c> (<see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>). An empty side is an
    /// all-gap run; identical segments are the trivial diagonal; everything else goes through the canonical
    /// Needleman-Wunsch / Gotoh (linear-space variant above <see cref="MaxGapCellsForFullMatrix"/>
    /// cells), so the gap alignment is always optimal.
    /// </summary>
    private static (string Aligned1, string Aligned2, int Score) AlignGap(
        string refGap, string queryGap, PairScoring scoring, bool affineGaps)
    {
        if (refGap.Length == 0 && queryGap.Length == 0)
            return ("", "", 0);

        int open = affineGaps ? scoring.GapOpen : 0;

        if (refGap.Length == 0)
            return (new string('-', queryGap.Length), queryGap, open + queryGap.Length * scoring.GapExtend);

        if (queryGap.Length == 0)
            return (refGap, new string('-', refGap.Length), open + refGap.Length * scoring.GapExtend);

        if (refGap == queryGap)
            return (refGap, queryGap, IdenticalScore(refGap, scoring));

        var result = SequenceAligner.GlobalLinearSpace(refGap, queryGap, scoring, affineGaps, MaxGapCellsForFullMatrix);
        return (result.AlignedSequence1, result.AlignedSequence2, result.Score);
    }

    /// <summary>Score of an exact match run: Σ matrix(c, c) (Length · Match for the scalar model).</summary>
    private static int IdenticalScore(string text, PairScoring scoring)
    {
        int total = 0;
        foreach (char c in text)
            total += scoring.Sub(c, c);
        return total;
    }
}
