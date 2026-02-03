using System;
using System.Collections.Generic;
using System.Text;
using System.Threading;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics;

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
        return CancellableOperations.GlobalAlign(sequence1, sequence2, scoring, cancellationToken, progress);
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

    #region Multiple Sequence Alignment (Simple)

    /// <summary>
    /// Performs progressive multiple sequence alignment.
    /// Uses a simple star alignment approach with the first sequence as reference.
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

        // Use first sequence as reference
        var aligned = new List<string> { seqList[0].Sequence };
        int totalScore = 0;

        for (int i = 1; i < seqList.Count; i++)
        {
            var result = GlobalAlign(seqList[0], seqList[i], effectiveScoring);
            aligned.Add(result.AlignedSequence2);
            totalScore += result.Score;
        }

        // Pad sequences to same length
        int maxLen = aligned.Max(s => s.Length);
        for (int i = 0; i < aligned.Count; i++)
        {
            if (aligned[i].Length < maxLen)
                aligned[i] = aligned[i].PadRight(maxLen, '-');
        }

        // Generate consensus
        var consensus = new StringBuilder(maxLen);
        for (int pos = 0; pos < maxLen; pos++)
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

        return new MultipleAlignmentResult(
            AlignedSequences: aligned.ToArray(),
            Consensus: consensus.ToString(),
            TotalScore: totalScore);
    }

    #endregion
}
