using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Performs pairwise sequence alignment using dynamic programming algorithms.
/// Supports global (Needleman-Wunsch), local (Smith-Waterman), and semi-global alignment.
/// </summary>
public static class SequenceAligner
{
    #region Scoring Matrices

    /// <summary>
    /// Simple DNA scoring: +1 match, -1 mismatch.
    /// </summary>
    public static readonly ScoringMatrix SimpleDna = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -2,
        GapExtend: -1);

    /// <summary>
    /// BLAST default DNA scoring: +2 match, -3 mismatch.
    /// </summary>
    public static readonly ScoringMatrix BlastDna = new(
        Match: 2,
        Mismatch: -3,
        GapOpen: -5,
        GapExtend: -2);

    /// <summary>
    /// High identity DNA scoring for closely related sequences.
    /// </summary>
    public static readonly ScoringMatrix HighIdentityDna = new(
        Match: 5,
        Mismatch: -4,
        GapOpen: -10,
        GapExtend: -1);

    #endregion

    #region Global Alignment (Needleman-Wunsch)

    /// <summary>
    /// Performs global alignment using the Needleman-Wunsch algorithm.
    /// Aligns entire sequences end-to-end.
    /// </summary>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    public static AlignmentResult GlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs global alignment on raw sequence strings.
    /// </summary>
    public static AlignmentResult GlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

        return GlobalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs global alignment with cancellation support.
    /// Recommended for aligning long sequences.
    /// </summary>
    /// <param name="sequence1">First sequence string.</param>
    /// <param name="sequence2">Second sequence string.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Alignment result with aligned sequences and score.</returns>
    public static AlignmentResult GlobalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

        var seq1 = sequence1.ToUpperInvariant();
        var seq2 = sequence2.ToUpperInvariant();
        var score = scoring ?? SimpleDna;

        int m = seq1.Length;
        int n = seq2.Length;

        // Initialize scoring matrix
        var matrix = new int[m + 1, n + 1];

        // Initialize first row and column
        for (int i = 0; i <= m; i++)
            matrix[i, 0] = i * score.GapExtend + (i > 0 ? score.GapOpen : 0);
        for (int j = 0; j <= n; j++)
            matrix[0, j] = j * score.GapExtend + (j > 0 ? score.GapOpen : 0);

        // Fill the matrix with cancellation checks
        for (int i = 1; i <= m; i++)
        {
            if (i % 100 == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                progress?.Report((double)i / m * 0.5); // First half is matrix fill
            }

            for (int j = 1; j <= n; j++)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? score.Match : score.Mismatch;

                int diag = matrix[i - 1, j - 1] + matchScore;
                int up = matrix[i - 1, j] + score.GapExtend;
                int left = matrix[i, j - 1] + score.GapExtend;

                matrix[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Traceback
        cancellationToken.ThrowIfCancellationRequested();
        progress?.Report(0.75);

        var aligned1 = new System.Text.StringBuilder();
        var aligned2 = new System.Text.StringBuilder();
        int ii = m, jj = n;

        while (ii > 0 || jj > 0)
        {
            if ((ii + jj) % 200 == 0)
                cancellationToken.ThrowIfCancellationRequested();

            if (ii > 0 && jj > 0)
            {
                int matchScore = seq1[ii - 1] == seq2[jj - 1] ? score.Match : score.Mismatch;
                if (matrix[ii, jj] == matrix[ii - 1, jj - 1] + matchScore)
                {
                    aligned1.Insert(0, seq1[ii - 1]);
                    aligned2.Insert(0, seq2[jj - 1]);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0 && matrix[ii, jj] == matrix[ii - 1, jj] + score.GapExtend)
            {
                aligned1.Insert(0, seq1[ii - 1]);
                aligned2.Insert(0, '-');
                ii--;
            }
            else if (jj > 0)
            {
                aligned1.Insert(0, '-');
                aligned2.Insert(0, seq2[jj - 1]);
                jj--;
            }
            else
            {
                break;
            }
        }

        progress?.Report(1.0);

        return new AlignmentResult(
            AlignedSequence1: aligned1.ToString(),
            AlignedSequence2: aligned2.ToString(),
            Score: matrix[m, n],
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0,
            StartPosition2: 0,
            EndPosition1: seq1.Length - 1,
            EndPosition2: seq2.Length - 1);
    }

    /// <summary>
    /// Performs global alignment on DNA sequences with cancellation support.
    /// </summary>
    public static AlignmentResult GlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return GlobalAlign(sequence1.Sequence, sequence2.Sequence, scoring, cancellationToken, progress);
    }

    private static AlignmentResult GlobalAlignCore(string seq1, string seq2, ScoringMatrix scoring)
    {
        int m = seq1.Length;
        int n = seq2.Length;

        // Initialize scoring matrix
        var score = new int[m + 1, n + 1];

        // Initialize first row and column
        for (int i = 0; i <= m; i++)
            score[i, 0] = i * scoring.GapExtend + (i > 0 ? scoring.GapOpen : 0);
        for (int j = 0; j <= n; j++)
            score[0, j] = j * scoring.GapExtend + (j > 0 ? scoring.GapOpen : 0);

        // Fill the matrix
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                int diag = score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Traceback
        return Traceback(seq1, seq2, score, m, n, scoring, AlignmentType.Global);
    }

    #endregion

    #region Local Alignment (Smith-Waterman)

    /// <summary>
    /// Performs local alignment using the Smith-Waterman algorithm.
    /// Finds the best local alignment between subsequences.
    /// </summary>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result with aligned subsequences and score.</returns>
    public static AlignmentResult LocalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return LocalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    /// <summary>
    /// Performs local alignment on raw sequence strings.
    /// </summary>
    public static AlignmentResult LocalAlign(
        string sequence1,
        string sequence2,
        ScoringMatrix? scoring = null)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            return AlignmentResult.Empty;

        return LocalAlignCore(
            sequence1.ToUpperInvariant(),
            sequence2.ToUpperInvariant(),
            scoring ?? SimpleDna);
    }

    private static AlignmentResult LocalAlignCore(string seq1, string seq2, ScoringMatrix scoring)
    {
        int m = seq1.Length;
        int n = seq2.Length;

        var score = new int[m + 1, n + 1];
        int maxScore = 0;
        int maxI = 0, maxJ = 0;

        // Fill the matrix (with zero floor for local alignment)
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                int diag = score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(0, Math.Max(diag, Math.Max(up, left)));

                if (score[i, j] > maxScore)
                {
                    maxScore = score[i, j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        // Traceback from max score
        return TracebackLocal(seq1, seq2, score, maxI, maxJ, scoring);
    }

    private static AlignmentResult TracebackLocal(
        string seq1, string seq2, int[,] score, int endI, int endJ, ScoringMatrix scoring)
    {
        var aligned1 = new StringBuilder();
        var aligned2 = new StringBuilder();

        int i = endI, j = endJ;
        int startI = endI, startJ = endJ;

        while (i > 0 && j > 0 && score[i, j] > 0)
        {
            startI = i;
            startJ = j;

            int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

            if (score[i, j] == score[i - 1, j - 1] + matchScore)
            {
                aligned1.Insert(0, seq1[i - 1]);
                aligned2.Insert(0, seq2[j - 1]);
                i--; j--;
            }
            else if (score[i, j] == score[i - 1, j] + scoring.GapExtend)
            {
                aligned1.Insert(0, seq1[i - 1]);
                aligned2.Insert(0, '-');
                i--;
            }
            else
            {
                aligned1.Insert(0, '-');
                aligned2.Insert(0, seq2[j - 1]);
                j--;
            }
        }

        return new AlignmentResult(
            AlignedSequence1: aligned1.ToString(),
            AlignedSequence2: aligned2.ToString(),
            Score: score[endI, endJ],
            AlignmentType: AlignmentType.Local,
            StartPosition1: startI - 1,
            StartPosition2: startJ - 1,
            EndPosition1: endI - 1,
            EndPosition2: endJ - 1);
    }

    #endregion

    #region Semi-Global Alignment

    /// <summary>
    /// Performs semi-global alignment (free end gaps).
    /// Useful for aligning a shorter sequence to a longer one.
    /// </summary>
    /// <param name="sequence1">First DNA sequence (typically shorter/query).</param>
    /// <param name="sequence2">Second DNA sequence (typically longer/reference).</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Alignment result.</returns>
    public static AlignmentResult SemiGlobalAlign(
        DnaSequence sequence1,
        DnaSequence sequence2,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        return SemiGlobalAlignCore(sequence1.Sequence, sequence2.Sequence, scoring ?? SimpleDna);
    }

    private static AlignmentResult SemiGlobalAlignCore(string seq1, string seq2, ScoringMatrix scoring)
    {
        int m = seq1.Length;
        int n = seq2.Length;

        var score = new int[m + 1, n + 1];

        // Free gaps at start of seq2 (first row is 0)
        for (int i = 1; i <= m; i++)
            score[i, 0] = i * scoring.GapExtend;

        // Fill the matrix
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                int diag = score[i - 1, j - 1] + matchScore;
                int up = score[i - 1, j] + scoring.GapExtend;
                int left = score[i, j - 1] + scoring.GapExtend;

                score[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Find max in last row (free gaps at end of seq2)
        int maxScore = score[m, 0];
        int maxJ = 0;
        for (int j = 1; j <= n; j++)
        {
            if (score[m, j] > maxScore)
            {
                maxScore = score[m, j];
                maxJ = j;
            }
        }

        return Traceback(seq1, seq2, score, m, maxJ, scoring, AlignmentType.SemiGlobal);
    }

    #endregion

    #region Traceback

    private static AlignmentResult Traceback(
        string seq1, string seq2, int[,] score, int i, int j,
        ScoringMatrix scoring, AlignmentType alignType)
    {
        var aligned1 = new StringBuilder();
        var aligned2 = new StringBuilder();

        // Add trailing gaps for semi-global
        if (alignType == AlignmentType.SemiGlobal)
        {
            for (int k = seq2.Length; k > j; k--)
            {
                aligned1.Insert(0, '-');
                aligned2.Insert(0, seq2[k - 1]);
            }
        }

        while (i > 0 || j > 0)
        {
            if (i > 0 && j > 0)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                if (score[i, j] == score[i - 1, j - 1] + matchScore)
                {
                    aligned1.Insert(0, seq1[i - 1]);
                    aligned2.Insert(0, seq2[j - 1]);
                    i--; j--;
                    continue;
                }
            }

            if (i > 0 && (j == 0 || score[i, j] == score[i - 1, j] + scoring.GapExtend))
            {
                aligned1.Insert(0, seq1[i - 1]);
                aligned2.Insert(0, '-');
                i--;
            }
            else
            {
                aligned1.Insert(0, '-');
                aligned2.Insert(0, seq2[j - 1]);
                j--;
            }
        }

        return new AlignmentResult(
            AlignedSequence1: aligned1.ToString(),
            AlignedSequence2: aligned2.ToString(),
            Score: score[seq1.Length, alignType == AlignmentType.SemiGlobal ? aligned2.ToString().Replace("-", "").Length : seq2.Length],
            AlignmentType: alignType,
            StartPosition1: 0,
            StartPosition2: 0,
            EndPosition1: seq1.Length - 1,
            EndPosition2: seq2.Length - 1);
    }

    #endregion

    #region Alignment Statistics

    /// <summary>
    /// Calculates alignment statistics from an alignment result.
    /// </summary>
    public static AlignmentStatistics CalculateStatistics(AlignmentResult alignment)
    {
        ArgumentNullException.ThrowIfNull(alignment);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return AlignmentStatistics.Empty;

        int matches = 0, mismatches = 0, gaps = 0;
        int alignmentLength = alignment.AlignedSequence1.Length;

        for (int i = 0; i < alignmentLength; i++)
        {
            char c1 = alignment.AlignedSequence1[i];
            char c2 = alignment.AlignedSequence2[i];

            if (c1 == '-' || c2 == '-')
                gaps++;
            else if (c1 == c2)
                matches++;
            else
                mismatches++;
        }

        double identity = alignmentLength > 0 ? (double)matches / alignmentLength * 100 : 0;
        double similarity = alignmentLength > 0 ? (double)(matches + mismatches) / alignmentLength * 100 : 0;
        double gapPercent = alignmentLength > 0 ? (double)gaps / alignmentLength * 100 : 0;

        return new AlignmentStatistics(
            Matches: matches,
            Mismatches: mismatches,
            Gaps: gaps,
            AlignmentLength: alignmentLength,
            Identity: identity,
            Similarity: similarity,
            GapPercent: gapPercent);
    }

    /// <summary>
    /// Generates a visual alignment string showing matches.
    /// </summary>
    public static string FormatAlignment(AlignmentResult alignment, int lineWidth = 60)
    {
        ArgumentNullException.ThrowIfNull(alignment);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return "";

        var sb = new StringBuilder();
        int length = alignment.AlignedSequence1.Length;

        for (int start = 0; start < length; start += lineWidth)
        {
            int end = Math.Min(start + lineWidth, length);

            sb.AppendLine(alignment.AlignedSequence1[start..end]);

            // Match line
            for (int i = start; i < end; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                char c2 = alignment.AlignedSequence2[i];

                if (c1 == c2 && c1 != '-')
                    sb.Append('|');
                else if (c1 == '-' || c2 == '-')
                    sb.Append(' ');
                else
                    sb.Append('.');
            }
            sb.AppendLine();

            sb.AppendLine(alignment.AlignedSequence2[start..end]);
            sb.AppendLine();
        }

        return sb.ToString();
    }

    #endregion

    #region Multiple Sequence Alignment (Anchor-Based)

    /// <summary>
    /// Performs progressive multiple sequence alignment using an anchor-based approach.
    /// <para>
    /// <b>Algorithm:</b>
    /// <list type="number">
    /// <item>Selects the best center sequence using suffix tree LCS distances (O(k²·L))</item>
    /// <item>Builds a suffix tree on the center sequence (O(L))</item>
    /// <item>For each other sequence, finds exact-match anchors via the suffix tree (O(L))</item>
    /// <item>Applies Needleman-Wunsch only to gaps between anchors (O(Σδᵢ²))</item>
    /// </list>
    /// </para>
    /// <para>
    /// For closely related sequences where anchors cover ~80% of length, this yields
    /// approximately 25× speedup over the standard O(k·L²) approach.
    /// </para>
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: SimpleDna).</param>
    /// <returns>Multiple alignment result.</returns>
    public static MultipleAlignmentResult MultipleAlign(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        if (seqList.Count == 1)
        {
            return new MultipleAlignmentResult(
                AlignedSequences: new[] { seqList[0].Sequence },
                Consensus: seqList[0].Sequence,
                TotalScore: 0);
        }

        var effectiveScoring = scoring ?? SimpleDna;

        // For 2 sequences or very short sequences, fall back to standard NW
        // (anchor overhead not worth it for short sequences)
        const int anchorThreshold = 30;
        bool useAnchors = seqList.Count >= 2
            && seqList.All(s => s.Length >= anchorThreshold);

        if (!useAnchors)
        {
            return MultipleAlignClassic(seqList, effectiveScoring);
        }

        // Step 1: Select the best center sequence using LCS-based distances.
        // The center is the sequence that maximizes total similarity to all others.
        int centerIdx = SelectCenterSequence(seqList);

        var centerSeq = seqList[centerIdx];
        string centerStr = centerSeq.Sequence;

        // Step 2: Build suffix tree once on the center sequence — O(L)
        var centerTree = SuffixTree.SuffixTree.Build(centerStr);

        // Step 3: Align all other sequences to center using anchor-based approach
        var aligned = new List<string>(seqList.Count);
        int totalScore = 0;

        // Pre-fill slots
        for (int i = 0; i < seqList.Count; i++)
            aligned.Add("");

        aligned[centerIdx] = centerStr;

        for (int i = 0; i < seqList.Count; i++)
        {
            if (i == centerIdx) continue;

            var result = AnchorBasedAligner.AlignWithAnchors(
                centerStr, centerTree, seqList[i].Sequence, effectiveScoring);

            aligned[centerIdx] = result.AlignedSequence1; // May have gaps inserted
            aligned[i] = result.AlignedSequence2;
            totalScore += result.Score;
        }

        // Step 4: Reconcile multiple pairwise alignments into a single MSA.
        // Since we aligned each seq to center independently, center may have different
        // gap patterns. We merge by re-aligning to the longest center representation.
        var (mergedAligned, mergedScore) = ReconcileAlignments(
            seqList, centerIdx, centerStr, centerTree, effectiveScoring);

        // Pad sequences to same length
        int maxLen = mergedAligned.Max(s => s.Length);
        for (int i = 0; i < mergedAligned.Count; i++)
        {
            if (mergedAligned[i].Length < maxLen)
                mergedAligned[i] = mergedAligned[i].PadRight(maxLen, '-');
        }

        // Generate consensus
        string consensus = BuildConsensus(mergedAligned, maxLen);

        return new MultipleAlignmentResult(
            AlignedSequences: mergedAligned.ToArray(),
            Consensus: consensus,
            TotalScore: mergedScore);
    }

    /// <summary>
    /// Classic star alignment (fallback for short sequences or small sets).
    /// </summary>
    internal static MultipleAlignmentResult MultipleAlignClassic(
        List<DnaSequence> seqList, ScoringMatrix scoring)
    {
        var aligned = new List<string> { seqList[0].Sequence };
        int totalScore = 0;

        for (int i = 1; i < seqList.Count; i++)
        {
            var result = GlobalAlign(seqList[0], seqList[i], scoring);
            aligned.Add(result.AlignedSequence2);
            totalScore += result.Score;
        }

        int maxLen = aligned.Max(s => s.Length);
        for (int i = 0; i < aligned.Count; i++)
        {
            if (aligned[i].Length < maxLen)
                aligned[i] = aligned[i].PadRight(maxLen, '-');
        }

        string consensus = BuildConsensus(aligned, maxLen);

        return new MultipleAlignmentResult(
            AlignedSequences: aligned.ToArray(),
            Consensus: consensus,
            TotalScore: totalScore);
    }

    /// <summary>
    /// Selects the center sequence for star alignment by finding the sequence
    /// with the highest total LCS similarity to all others.
    /// Uses suffix trees for O(L) per-pair comparison instead of O(L²) NW.
    /// </summary>
    private static int SelectCenterSequence(List<DnaSequence> sequences)
    {
        int k = sequences.Count;
        if (k <= 2) return 0;

        var totalSimilarity = new int[k];

        // For each sequence, compute sum of LCS lengths to all others
        for (int i = 0; i < k; i++)
        {
            var tree = sequences[i].SuffixTree;
            for (int j = 0; j < k; j++)
            {
                if (i == j) continue;
                string lcs = tree.LongestCommonSubstring(sequences[j].Sequence);
                totalSimilarity[i] += lcs.Length;
            }
        }

        // Return index of sequence with maximum total similarity
        int bestIdx = 0;
        for (int i = 1; i < k; i++)
        {
            if (totalSimilarity[i] > totalSimilarity[bestIdx])
                bestIdx = i;
        }

        return bestIdx;
    }

    /// <summary>
    /// Reconciles multiple pairwise alignments against the center into a single MSA.
    /// Uses the simple approach of aligning each sequence to center independently
    /// and then merging gap columns.
    /// </summary>
    private static (List<string> Aligned, int TotalScore) ReconcileAlignments(
        List<DnaSequence> sequences,
        int centerIdx,
        string centerStr,
        SuffixTree.SuffixTree centerTree,
        ScoringMatrix scoring)
    {
        int k = sequences.Count;

        // Perform all pairwise alignments to center
        var pairwiseResults = new AlignmentResult[k];
        int totalScore = 0;

        for (int i = 0; i < k; i++)
        {
            if (i == centerIdx)
            {
                pairwiseResults[i] = new AlignmentResult(
                    centerStr, centerStr, 0, AlignmentType.Global, 0, 0,
                    centerStr.Length - 1, centerStr.Length - 1);
                continue;
            }

            pairwiseResults[i] = AnchorBasedAligner.AlignWithAnchors(
                centerStr, centerTree, sequences[i].Sequence, scoring);
            totalScore += pairwiseResults[i].Score;
        }

        // Merge: build a gap map from all center alignments.
        // For each position in the original center sequence, track how many gaps
        // need to be inserted before it (maximum across all pairwise alignments).
        var gapsBefore = new int[centerStr.Length + 1];

        for (int i = 0; i < k; i++)
        {
            if (i == centerIdx) continue;

            string alignedCenter = pairwiseResults[i].AlignedSequence1;
            int origPos = 0;
            int gapCount = 0;

            for (int j = 0; j < alignedCenter.Length; j++)
            {
                if (alignedCenter[j] == '-')
                {
                    gapCount++;
                }
                else
                {
                    gapsBefore[origPos] = Math.Max(gapsBefore[origPos], gapCount);
                    origPos++;
                    gapCount = 0;
                }
            }
            // Trailing gaps
            gapsBefore[centerStr.Length] = Math.Max(gapsBefore[centerStr.Length], gapCount);
        }

        // Now re-project each pairwise alignment into the merged coordinate space
        var merged = new List<string>(k);

        for (int i = 0; i < k; i++)
        {
            string alignedCenter = pairwiseResults[i].AlignedSequence1;
            string alignedOther = pairwiseResults[i].AlignedSequence2;

            var sb = new StringBuilder();
            int alignPos = 0;

            for (int p = 0; p <= centerStr.Length; p++)
            {
                int neededGaps = gapsBefore[p];

                // Count how many gaps this alignment has at this center position
                int thisGaps = 0;
                int tempAlignPos = alignPos;
                while (tempAlignPos < alignedCenter.Length && alignedCenter[tempAlignPos] == '-')
                {
                    thisGaps++;
                    tempAlignPos++;
                }

                // Emit this alignment's gaps (from the other sequence)
                for (int g = 0; g < thisGaps; g++)
                {
                    if (alignPos < alignedOther.Length)
                        sb.Append(alignedOther[alignPos]);
                    else
                        sb.Append('-');
                    alignPos++;
                }

                // Emit additional gaps to match the maximum
                for (int g = thisGaps; g < neededGaps; g++)
                {
                    sb.Append('-');
                }

                // Emit the character at center position p
                if (p < centerStr.Length && alignPos < alignedOther.Length)
                {
                    sb.Append(alignedOther[alignPos]);
                    alignPos++;
                }
            }

            // Emit any remaining characters
            while (alignPos < alignedOther.Length)
            {
                sb.Append(alignedOther[alignPos]);
                alignPos++;
            }

            merged.Add(sb.ToString());
        }

        return (merged, totalScore);
    }

    /// <summary>
    /// Builds a majority-vote consensus from aligned sequences.
    /// </summary>
    private static string BuildConsensus(List<string> aligned, int length)
    {
        var consensus = new StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0, ['-'] = 0 };

            foreach (var seq in aligned)
            {
                if (pos < seq.Length && counts.ContainsKey(seq[pos]))
                    counts[seq[pos]]++;
            }

            char mostCommon = counts.Where(kv => kv.Key != '-')
                                   .OrderByDescending(kv => kv.Value)
                                   .FirstOrDefault().Key;

            consensus.Append(mostCommon == default ? '-' : mostCommon);
        }

        return consensus.ToString();
    }

    #endregion
}

