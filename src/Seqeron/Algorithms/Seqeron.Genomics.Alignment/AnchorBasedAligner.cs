using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Anchor-based sequence aligner that uses a suffix tree to find exact-match anchors
/// between sequences, then applies Needleman-Wunsch only to the gaps between anchors.
/// <para>
/// This approach is inspired by MUMmer (Delcher et al., 1999) and LAGAN (Brudno et al., 2003),
/// which use suffix trees to find Maximal Exact Matches (MEMs) as anchors and align only
/// the intervening regions with dynamic programming.
/// </para>
/// <para>
/// <b>Complexity improvement:</b>
/// <list type="bullet">
/// <item>Standard NW: O(L²) per pair</item>
/// <item>Anchor-based: O(L) for tree build + O(L) for anchor finding + O(Σδᵢ²) for gap NW</item>
/// </list>
/// When anchors cover fraction α of the sequence, the effective cost drops by ~(1-α)².
/// For closely related sequences (α ≈ 0.8), this yields a ~25× speedup.
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
    /// Minimum gap length to trigger NW alignment. Gaps shorter than this
    /// are trivially aligned character-by-character.
    /// </summary>
    private const int TrivialGapThreshold = 3;

    /// <summary>
    /// Maximum gap size for which we run full NW. Larger gaps fall back
    /// to standard alignment to avoid excessive memory.
    /// </summary>
    private const int MaxGapForNW = 5000;

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
        if (string.IsNullOrEmpty(reference) || string.IsNullOrEmpty(query))
            return AlignmentResult.Empty;

        // Step 1: Find all Maximal Exact Matches (MEMs) using O(n+m) suffix-link streaming
        var memTuples = refTree.FindExactMatchAnchors(query, minAnchorLength);
        var mems = new List<Anchor>(memTuples.Count);
        foreach (var (refPos, queryPos, len) in memTuples)
            mems.Add(new Anchor(refPos, queryPos, len));

        // Step 2: Chain anchors — find the longest consistent chain (monotonic in both coords)
        var chain = ChainAnchors(mems);

        // Step 3: If we have too few anchors, fall back to standard NW
        if (chain.Count == 0)
        {
            return SequenceAligner.GlobalAlign(reference, query, scoring);
        }

        // Step 4: Align gaps between anchors using NW, stitch everything together
        return StitchAlignment(reference, query, chain, scoring);
    }

    /// <summary>
    /// Chains anchors using a Longest Increasing Subsequence (LIS) approach
    /// to find the longest chain of non-overlapping anchors consistent in both
    /// reference and query coordinates.
    /// </summary>
    internal static List<Anchor> ChainAnchors(List<Anchor> anchors)
    {
        if (anchors.Count == 0)
            return new List<Anchor>();

        // Sort by reference position
        var sorted = anchors.OrderBy(a => a.RefStart).ThenBy(a => a.QueryStart).ToList();

        // Remove overlapping anchors — keep the longest at each position
        var cleaned = RemoveOverlaps(sorted);

        if (cleaned.Count == 0)
            return new List<Anchor>();

        // LIS on query positions to find longest consistent chain
        return LongestIncreasingSubsequenceChain(cleaned);
    }

    /// <summary>
    /// Removes overlapping anchors, preferring longer ones.
    /// </summary>
    private static List<Anchor> RemoveOverlaps(List<Anchor> sorted)
    {
        var result = new List<Anchor>();

        foreach (var anchor in sorted)
        {
            // Check if this anchor overlaps with the last accepted anchor
            if (result.Count > 0)
            {
                var last = result[^1];
                if (anchor.RefStart < last.RefEnd || anchor.QueryStart < last.QueryEnd)
                {
                    // Overlap — keep the longer one
                    if (anchor.Length > last.Length)
                    {
                        result[^1] = anchor;
                    }
                    continue;
                }
            }

            result.Add(anchor);
        }

        return result;
    }

    /// <summary>
    /// Finds the longest chain of anchors where both RefStart and QueryStart
    /// are strictly increasing (and non-overlapping). Uses patience sorting (LIS).
    /// </summary>
    private static List<Anchor> LongestIncreasingSubsequenceChain(List<Anchor> anchors)
    {
        int n = anchors.Count;
        if (n == 0) return new List<Anchor>();
        if (n == 1) return new List<Anchor> { anchors[0] };

        // dp[i] = length of LIS ending at anchors[i]
        var dp = new int[n];
        var prev = new int[n];

        for (int i = 0; i < n; i++)
        {
            dp[i] = 1;
            prev[i] = -1;
        }

        int bestLen = 1, bestIdx = 0;

        for (int i = 1; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                // Anchor j must end before anchor i starts in both coords
                if (anchors[j].RefEnd <= anchors[i].RefStart &&
                    anchors[j].QueryEnd <= anchors[i].QueryStart &&
                    dp[j] + 1 > dp[i])
                {
                    dp[i] = dp[j] + 1;
                    prev[i] = j;
                }
            }

            if (dp[i] > bestLen)
            {
                bestLen = dp[i];
                bestIdx = i;
            }
        }

        // Reconstruct the chain
        var chain = new List<Anchor>();
        int idx = bestIdx;
        while (idx >= 0)
        {
            chain.Add(anchors[idx]);
            idx = prev[idx];
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
    /// Aligns a gap segment. For very short or equal gaps, uses trivial alignment.
    /// For moderate gaps, uses standard NW (anchors keep gaps small, typically &lt;100bp).
    /// For gaps exceeding MaxGapForNW, uses a simpler heuristic.
    /// </summary>
    private static (string Aligned1, string Aligned2, int Score) AlignGap(
        string refGap, string queryGap, ScoringMatrix scoring)
    {
        // Both empty — nothing to do
        if (refGap.Length == 0 && queryGap.Length == 0)
            return ("", "", 0);

        // One side empty — all gaps
        if (refGap.Length == 0)
        {
            string gaps = new string('-', queryGap.Length);
            int score = scoring.GapOpen + queryGap.Length * scoring.GapExtend;
            return (gaps, queryGap, score);
        }

        if (queryGap.Length == 0)
        {
            string gaps = new string('-', refGap.Length);
            int score = scoring.GapOpen + refGap.Length * scoring.GapExtend;
            return (refGap, gaps, score);
        }

        // Identical gaps — trivial match
        if (refGap == queryGap)
        {
            int score = refGap.Length * scoring.Match;
            return (refGap, queryGap, score);
        }

        // Very large gap — fall back to simple diagonal alignment
        if (refGap.Length > MaxGapForNW || queryGap.Length > MaxGapForNW)
        {
            return SimpleAlignLargeGap(refGap, queryGap, scoring);
        }

        // Standard NW for remaining gaps — anchors keep gaps small (typically <100bp),
        // so full NW is cheap and always finds the optimal alignment.
        var result = SequenceAligner.GlobalAlign(refGap, queryGap, scoring);
        return (result.AlignedSequence1, result.AlignedSequence2, result.Score);
    }

    /// <summary>
    /// Simple alignment for very large gaps where full NW would be too expensive.
    /// Aligns matching prefix, pads the difference with gaps, aligns matching suffix.
    /// </summary>
    private static (string Aligned1, string Aligned2, int Score) SimpleAlignLargeGap(
        string refGap, string queryGap, ScoringMatrix scoring)
    {
        var a1 = new StringBuilder();
        var a2 = new StringBuilder();
        int score = 0;

        int minLen = Math.Min(refGap.Length, queryGap.Length);

        // Align character-by-character up to min length
        for (int i = 0; i < minLen; i++)
        {
            a1.Append(refGap[i]);
            a2.Append(queryGap[i]);
            score += refGap[i] == queryGap[i] ? scoring.Match : scoring.Mismatch;
        }

        // Pad the shorter one with gaps
        if (refGap.Length > minLen)
        {
            string tail = refGap[minLen..];
            a1.Append(tail);
            a2.Append(new string('-', tail.Length));
            score += scoring.GapOpen + tail.Length * scoring.GapExtend;
        }
        else if (queryGap.Length > minLen)
        {
            string tail = queryGap[minLen..];
            a1.Append(new string('-', tail.Length));
            a2.Append(tail);
            score += scoring.GapOpen + tail.Length * scoring.GapExtend;
        }

        return (a1.ToString(), a2.ToString(), score);
    }
}
