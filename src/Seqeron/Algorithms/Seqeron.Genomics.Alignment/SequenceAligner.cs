using System.Buffers;
using System.Text;

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

        // Initialize first row and column with linear gap penalty d = GapExtend.
        // Standard Needleman-Wunsch: F(i,0) = d*i, F(0,j) = d*j
        // Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
        for (int i = 0; i <= m; i++)
            matrix[i, 0] = i * score.GapExtend;
        for (int j = 0; j <= n; j++)
            matrix[0, j] = j * score.GapExtend;

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

        var chars1 = new List<char>();
        var chars2 = new List<char>();
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
                    chars1.Add(seq1[ii - 1]);
                    chars2.Add(seq2[jj - 1]);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0 && (jj == 0 || matrix[ii, jj] == matrix[ii - 1, jj] + score.GapExtend))
            {
                chars1.Add(seq1[ii - 1]);
                chars2.Add('-');
                ii--;
            }
            else if (jj > 0)
            {
                chars1.Add('-');
                chars2.Add(seq2[jj - 1]);
                jj--;
            }
            else
            {
                break;
            }
        }

        chars1.Reverse();
        chars2.Reverse();

        progress?.Report(1.0);

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
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

        // Use pooled flat array instead of 2D array to reduce GC pressure.
        // For anchor-based MSA, this method is called 100+ times per alignment
        // with small gap segments (typically <100bp), so pooling is effective.
        int rows = m + 1;
        int cols = n + 1;
        int totalCells = rows * cols;
        var pool = ArrayPool<int>.Shared;
        int[] buf = pool.Rent(totalCells);

        try
        {
            // Initialize first row and column with linear gap penalty d = GapExtend.
            // Standard Needleman-Wunsch: F(i,0) = d*i, F(0,j) = d*j
            // Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
            for (int i = 0; i <= m; i++)
                buf[i * cols] = i * scoring.GapExtend;
            for (int j = 0; j <= n; j++)
                buf[j] = j * scoring.GapExtend;

            // Fill the matrix
            for (int i = 1; i <= m; i++)
            {
                int rowOff = i * cols;
                int prevRowOff = (i - 1) * cols;
                char c1 = seq1[i - 1];

                for (int j = 1; j <= n; j++)
                {
                    int matchScore = c1 == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                    int diag = buf[prevRowOff + (j - 1)] + matchScore;
                    int up = buf[prevRowOff + j] + scoring.GapExtend;
                    int left = buf[rowOff + (j - 1)] + scoring.GapExtend;

                    buf[rowOff + j] = Math.Max(diag, Math.Max(up, left));
                }
            }

            // Copy to 2D for Traceback (Traceback uses int[,])
            var score = new int[rows, cols];
            Buffer.BlockCopy(buf, 0, score, 0, totalCells * sizeof(int));

            return Traceback(seq1, seq2, score, m, n, scoring, AlignmentType.Global);
        }
        finally
        {
            pool.Return(buf);
        }
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
        var chars1 = new List<char>();
        var chars2 = new List<char>();

        int i = endI, j = endJ;
        int startI = endI, startJ = endJ;

        while (i > 0 && j > 0 && score[i, j] > 0)
        {
            startI = i;
            startJ = j;

            int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

            if (score[i, j] == score[i - 1, j - 1] + matchScore)
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add(seq2[j - 1]);
                i--; j--;
            }
            else if (score[i, j] == score[i - 1, j] + scoring.GapExtend)
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add('-');
                i--;
            }
            else
            {
                chars1.Add('-');
                chars2.Add(seq2[j - 1]);
                j--;
            }
        }

        chars1.Reverse();
        chars2.Reverse();

        return new AlignmentResult(
            AlignedSequence1: new string(chars1.ToArray()),
            AlignedSequence2: new string(chars2.ToArray()),
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
        // Preserve the traceback origin for correct score retrieval.
        // For semi-global (fitting): score[m, maxJ] is the optimal fitting score.
        // For global (NW): score[m, n] is the bottom-right cell.
        int endI = i, endJ = j;

        var chars1 = new List<char>();
        var chars2 = new List<char>();

        // Add trailing gaps for semi-global (appended first, will be at end after reverse).
        // chars2 is reversed in its entirety below, so the unmatched trailing reference
        // suffix seq2[j .. n-1] must be appended in REVERSE order here; otherwise the
        // suffix would itself be reversed in the final output and AlignedSequence2 would
        // no longer reproduce the reference.
        if (alignType == AlignmentType.SemiGlobal)
        {
            for (int k = seq2.Length; k > j; k--)
            {
                chars1.Add('-');
                chars2.Add(seq2[k - 1]);
            }
        }

        while (i > 0 || j > 0)
        {
            if (i > 0 && j > 0)
            {
                int matchScore = seq1[i - 1] == seq2[j - 1] ? scoring.Match : scoring.Mismatch;

                if (score[i, j] == score[i - 1, j - 1] + matchScore)
                {
                    chars1.Add(seq1[i - 1]);
                    chars2.Add(seq2[j - 1]);
                    i--; j--;
                    continue;
                }
            }

            if (i > 0 && (j == 0 || score[i, j] == score[i - 1, j] + scoring.GapExtend))
            {
                chars1.Add(seq1[i - 1]);
                chars2.Add('-');
                i--;
            }
            else
            {
                chars1.Add('-');
                chars2.Add(seq2[j - 1]);
                j--;
            }
        }

        chars1.Reverse();
        chars2.Reverse();

        string aligned1Str = new string(chars1.ToArray());
        string aligned2Str = new string(chars2.ToArray());

        return new AlignmentResult(
            AlignedSequence1: aligned1Str,
            AlignedSequence2: aligned2Str,
            Score: score[endI, endJ],
            AlignmentType: alignType,
            StartPosition1: 0,
            StartPosition2: 0,
            EndPosition1: seq1.Length - 1,
            EndPosition2: seq2.Length - 1);
    }

    #endregion

    #region Alignment Statistics

    // EMBOSS srspair markup legend (Rice, Longden & Bleasby 2000; EMBOSS AlignFormats):
    //   '|' = identical column, ':' = similar column (positive substitution score, not
    //   identical), ' ' = gap or non-positive (mismatch) column.
    //   Source: https://emboss.sourceforge.net/docs/themes/AlignFormats.html
    private const char IdentityMark = '|';
    private const char SimilarityMark = ':';
    private const char GapOrMismatchMark = ' ';

    /// <summary>
    /// Calculates alignment statistics (Identity, Similarity, Gaps) from an alignment
    /// result, following the EMBOSS needle/water convention.
    /// </summary>
    /// <param name="alignment">The pairwise alignment to summarise.</param>
    /// <param name="scoring">
    /// Scoring model used to decide whether a non-identical aligned pair counts as
    /// "similar". A column is similar when its substitution score is positive
    /// (Rice et al. 2000). For the DNA models used by this class (Match &gt; 0,
    /// Mismatch &lt; 0) no non-identical pair scores positively, so Similarity equals
    /// Identity. Defaults to <see cref="SimpleDna"/>.
    /// </param>
    /// <remarks>
    /// Identity = identical columns / alignment length × 100 (length includes gap
    /// columns); Similarity = (identical + similar) columns / length × 100;
    /// Gaps = gap columns / length × 100. Denominator is the full alignment length
    /// including gaps, per the EMBOSS needle definition.
    /// </remarks>
    public static AlignmentStatistics CalculateStatistics(
        AlignmentResult alignment,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return AlignmentStatistics.Empty;

        var score = scoring ?? SimpleDna;

        int matches = 0, mismatches = 0, gaps = 0, similarSubstitutions = 0;
        int alignmentLength = alignment.AlignedSequence1.Length;

        for (int i = 0; i < alignmentLength; i++)
        {
            char c1 = alignment.AlignedSequence1[i];
            char c2 = alignment.AlignedSequence2[i];

            if (c1 == '-' || c2 == '-')
            {
                gaps++;
            }
            else if (c1 == c2)
            {
                matches++;
            }
            else
            {
                mismatches++;
                // A non-identical column is "similar" iff its substitution score is
                // positive (Rice, Longden & Bleasby 2000).
                if (score.Mismatch > 0)
                    similarSubstitutions++;
            }
        }

        // Length (including gaps) is the EMBOSS needle denominator for all three
        // percentages. Source: EMBOSS needle documentation (Rice et al. 2000).
        int similar = matches + similarSubstitutions;
        double identity = (double)matches / alignmentLength * 100;
        double similarity = (double)similar / alignmentLength * 100;
        double gapPercent = (double)gaps / alignmentLength * 100;

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
    /// Generates a visual three-line pairwise alignment using the EMBOSS srspair
    /// markup legend ('|' identity, ':' similarity, ' ' gap/mismatch).
    /// </summary>
    /// <param name="alignment">The pairwise alignment to render.</param>
    /// <param name="lineWidth">Residues per block (must be positive; default 60).</param>
    /// <param name="scoring">
    /// Scoring model used to mark similar (positive-score) columns; defaults to
    /// <see cref="SimpleDna"/>. See <see cref="CalculateStatistics"/>.
    /// </param>
    public static string FormatAlignment(
        AlignmentResult alignment,
        int lineWidth = 60,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(lineWidth);

        if (string.IsNullOrEmpty(alignment.AlignedSequence1))
            return "";

        var score = scoring ?? SimpleDna;
        var sb = new StringBuilder();
        int length = alignment.AlignedSequence1.Length;

        for (int start = 0; start < length; start += lineWidth)
        {
            int end = Math.Min(start + lineWidth, length);

            // Emit '\n' explicitly (not AppendLine/Environment.NewLine) so the formatted alignment
            // is byte-identical across platforms — the EMBOSS-style block must not depend on OS CRLF.
            sb.Append(alignment.AlignedSequence1[start..end]).Append('\n');

            // Markup line (EMBOSS srspair legend).
            for (int i = start; i < end; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                char c2 = alignment.AlignedSequence2[i];

                if (c1 == '-' || c2 == '-')
                    sb.Append(GapOrMismatchMark);
                else if (c1 == c2)
                    sb.Append(IdentityMark);
                else if (score.Mismatch > 0)
                    sb.Append(SimilarityMark);
                else
                    sb.Append(GapOrMismatchMark);
            }
            sb.Append('\n');

            sb.Append(alignment.AlignedSequence2[start..end]).Append('\n');
            sb.Append('\n');
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

        // Step 1: Select the best center sequence using k-mer cosine similarity.
        // The center is the sequence that maximizes total similarity to all others.
        int centerIdx = SelectCenterSequence(seqList);

        var centerSeq = seqList[centerIdx];
        string centerStr = centerSeq.Sequence;

        // Step 2: Build suffix tree once on the center sequence — O(L)
        var centerTree = SuffixTree.SuffixTree.Build(centerStr);

        // Step 3: Reconcile pairwise alignments into a single MSA.
        // ReconcileAlignments performs all pairwise alignments (parallelized for k≥3)
        // and merges gap patterns from independent center-vs-other alignments.
        var (mergedAligned, _) = ReconcileAlignments(
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

        // Compute true sum-of-pairs score on the final aligned sequences.
        // Per Wikipedia MSA: SP = sum over all C(k,2) pairs of column-based scores.
        int spScore = ComputeSumOfPairsScore(mergedAligned, effectiveScoring);

        return new MultipleAlignmentResult(
            AlignedSequences: mergedAligned.ToArray(),
            Consensus: consensus,
            TotalScore: spScore);
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
            TotalScore: ComputeSumOfPairsScore(aligned, scoring));
    }

    /// <summary>
    /// Selects the center sequence for star alignment by finding the sequence
    /// with the highest total similarity to all others.
    /// Uses 4-mer frequency cosine similarity: O(k²·L) with small constant,
    /// instead of O(k²·L) suffix-tree LCS which builds k trees.
    /// </summary>
    private static int SelectCenterSequence(List<DnaSequence> sequences)
    {
        int k = sequences.Count;

        // The center is the star-alignment hub: every other sequence is aligned to it and
        // the final MSA is projected onto the center's coordinate space. An EMPTY center
        // has no positions to carry that projection, so non-empty sequences would lose all
        // their content during reconciliation. Always prefer a non-empty center when one
        // exists; only fall through to the heuristic among the non-empty candidates.
        int firstNonEmpty = -1;
        for (int i = 0; i < k; i++)
        {
            if (sequences[i].Sequence.Length > 0)
            {
                firstNonEmpty = i;
                break;
            }
        }

        // All sequences empty (or none at all): the trivial empty-coordinate MSA is correct.
        if (firstNonEmpty < 0) return 0;

        if (k <= 2) return firstNonEmpty;

        // Build 4-mer frequency vectors for each sequence.
        // 4-mer over ACGT → 256 possible k-mers, stored as int[256].
        const int KmerLen = 4;
        const int VocabSize = 256; // 4^4
        var profiles = new int[k][];
        var norms = new double[k];

        for (int i = 0; i < k; i++)
        {
            profiles[i] = new int[VocabSize];
            string seq = sequences[i].Sequence;
            for (int p = 0; p <= seq.Length - KmerLen; p++)
            {
                int hash = KmerHash(seq, p, KmerLen);
                if (hash >= 0)
                    profiles[i][hash]++;
            }

            double norm = 0;
            for (int v = 0; v < VocabSize; v++)
                norm += (double)profiles[i][v] * profiles[i][v];
            norms[i] = Math.Sqrt(norm);
        }

        // Compute pairwise cosine similarity; pick sequence with max total
        var totalSimilarity = new double[k];

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                double dot = 0;
                for (int v = 0; v < VocabSize; v++)
                    dot += (double)profiles[i][v] * profiles[j][v];

                double sim = (norms[i] > 0 && norms[j] > 0)
                    ? dot / (norms[i] * norms[j])
                    : 0;

                totalSimilarity[i] += sim;
                totalSimilarity[j] += sim;
            }
        }

        // Restrict the center choice to NON-EMPTY candidates (an empty center cannot carry
        // the reconciliation projection — see above). firstNonEmpty is guaranteed valid.
        int bestIdx = firstNonEmpty;
        for (int i = firstNonEmpty + 1; i < k; i++)
        {
            if (sequences[i].Sequence.Length == 0)
                continue;
            if (totalSimilarity[i] > totalSimilarity[bestIdx])
                bestIdx = i;
        }

        return bestIdx;
    }

    /// <summary>
    /// Hashes a 4-mer at the given position into [0..255].
    /// Returns -1 if any character is not A/C/G/T.
    /// </summary>
    private static int KmerHash(string seq, int pos, int len)
    {
        int hash = 0;
        for (int i = 0; i < len; i++)
        {
            int code = seq[pos + i] switch
            {
                'A' or 'a' => 0,
                'C' or 'c' => 1,
                'G' or 'g' => 2,
                'T' or 't' => 3,
                _ => -1
            };
            if (code < 0) return -1;
            hash = (hash << 2) | code;
        }
        return hash;
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

        // Perform all pairwise alignments to center — parallelized for k≥3
        var pairwiseResults = new AlignmentResult[k];
        int totalScore = 0;

        // Center aligns to itself (identity)
        pairwiseResults[centerIdx] = new AlignmentResult(
            centerStr, centerStr, 0, AlignmentType.Global, 0, 0,
            centerStr.Length - 1, centerStr.Length - 1);

        if (k >= 4)
        {
            // Parallel: suffix tree is immutable, each slot is independent
            Parallel.For(0, k, i =>
            {
                if (i == centerIdx) return;

                var result = AnchorBasedAligner.AlignWithAnchors(
                    centerStr, centerTree, sequences[i].Sequence, scoring);
                pairwiseResults[i] = result;
                Interlocked.Add(ref totalScore, result.Score);
            });
        }
        else
        {
            // Sequential for small k (avoid thread pool overhead)
            for (int i = 0; i < k; i++)
            {
                if (i == centerIdx) continue;

                pairwiseResults[i] = AnchorBasedAligner.AlignWithAnchors(
                    centerStr, centerTree, sequences[i].Sequence, scoring);
                totalScore += pairwiseResults[i].Score;
            }
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
    /// Per Wikipedia MSA: consensus derived from aligned columns using majority voting.
    /// Gap characters participate in the vote; on tie, nucleotides are preferred.
    /// </summary>
    private static string BuildConsensus(List<string> aligned, int length)
    {
        var consensus = new StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0, ['-'] = 0 };

            foreach (var seq in aligned)
            {
                if (pos < seq.Length && counts.TryGetValue(seq[pos], out int value))
                    counts[seq[pos]] = ++value;
            }

            // Include all characters (including gaps) in majority vote.
            // On tie between gap and nucleotide, prefer nucleotide.
            char mostCommon = counts.OrderByDescending(kv => kv.Value)
                                   .ThenBy(kv => kv.Key == '-' ? 1 : 0)
                                   .First().Key;

            consensus.Append(mostCommon);
        }

        return consensus.ToString();
    }

    /// <summary>
    /// Computes the column-based sum-of-pairs (SP) score for a multiple sequence alignment.
    /// Per Wikipedia MSA: "sum of all of the pairs of characters at each position in the alignment."
    /// Scores all C(k,2) sequence pairs using column-based scoring:
    /// match/mismatch from scoring matrix, gap-nucleotide penalty per position, gap-gap = 0.
    /// </summary>
    private static int ComputeSumOfPairsScore(List<string> aligned, ScoringMatrix scoring)
    {
        if (aligned.Count < 2) return 0;

        int length = aligned[0].Length;
        int totalScore = 0;

        for (int pos = 0; pos < length; pos++)
        {
            for (int i = 0; i < aligned.Count; i++)
            {
                char ci = pos < aligned[i].Length ? aligned[i][pos] : '-';

                for (int j = i + 1; j < aligned.Count; j++)
                {
                    char cj = pos < aligned[j].Length ? aligned[j][pos] : '-';

                    // Gap-gap is neutral (standard SP convention) and simply contributes nothing.
                    if (ci != '-' || cj != '-')
                    {
                        if (ci == '-' || cj == '-')
                            totalScore += scoring.GapExtend; // Gap-nucleotide penalty
                        else if (ci == cj)
                            totalScore += scoring.Match;
                        else
                            totalScore += scoring.Mismatch;
                    }
                }
            }
        }

        return totalScore;
    }

    #endregion

    #region Multiple Sequence Alignment (Guide-Tree Progressive / Feng-Doolittle)

    // Progressive (guide-tree) multiple sequence alignment, the Feng-Doolittle /
    // Clustal-style method, added as a SECOND aligner alongside the star MSA above.
    // The star MSA (MultipleAlign) is unchanged.
    //
    // Algorithm (Feng & Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    //   1. Compute all C(k,2) pairwise global (Needleman-Wunsch) alignments and derive a
    //      distance matrix. Distance = 1 - fractional identity, where fractional identity =
    //      identical columns / alignment length of the pairwise NW alignment. (Feng-Doolittle
    //      converts pairwise alignment scores to distances; identity-based distance is the
    //      standard, simplest such conversion used by Clustal-style tools.)
    //   2. Build a guide tree by UPGMA over that distance matrix — "an efficient clustering
    //      method such as neighbor-joining or unweighted pair group method with arithmetic
    //      mean (UPGMA)" (Wikipedia "Multiple sequence alignment"). UPGMA is the classic
    //      Feng-Doolittle choice. The UPGMA cluster-merge formula mirrors
    //      PhylogeneticAnalyzer.BuildUPGMA (proportional averaging, Sokal & Michener 1958);
    //      it is reimplemented here because the Phylogenetics project already depends on this
    //      Alignment project, so a reverse reference would be circular.
    //   3. Progressively align along the guide tree from the tips: sequence-sequence, then
    //      sequence-profile, then profile-profile, using the Needleman-Wunsch recurrence over
    //      *columns* (a profile is the set of already-aligned rows). Column-vs-column score is
    //      the average match/mismatch over all cross-profile residue pairs (sum-of-pairs profile
    //      scoring). The "once a gap, always a gap" rule (Feng & Doolittle 1987: "Once a gap is
    //      introduced ... it is preserved within all subsequent fusions") is enforced: gaps
    //      already present in a profile are never removed; merging two profiles only inserts
    //      whole new gap columns, never edits existing columns.
    //
    // Sources retrieved this session:
    //   - Feng & Doolittle (1987) "Progressive sequence alignment as a prerequisite to correct
    //     phylogenetic trees", J Mol Evol 25:351-360 (PubMed 3118049) — progressive method,
    //     guide tree, "once a gap, always a gap".
    //   - Wikipedia "Multiple sequence alignment", § Progressive alignment construction
    //     (https://en.wikipedia.org/wiki/Multiple_sequence_alignment) — guide tree built by
    //     NJ/UPGMA; sequences added sequentially along the guide tree.
    //   - Wikipedia UPGMA / Sokal & Michener (1958) — cluster averaging formula.

    /// <summary>
    /// Performs guide-tree <b>progressive</b> multiple sequence alignment (the Feng-Doolittle /
    /// Clustal-style method), as an additive alternative to the star-alignment
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Feng &amp; Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    /// </para>
    /// <list type="number">
    /// <item>All pairwise Needleman-Wunsch alignments give a distance matrix using
    /// distance = 1 − fractional identity (identical columns / pairwise alignment length).</item>
    /// <item>A UPGMA guide tree is built over that matrix (the classic Feng-Doolittle choice).</item>
    /// <item>Profiles are aligned progressively from the tips with the NW recurrence over
    /// columns and sum-of-pairs profile scoring, enforcing "once a gap, always a gap".</item>
    /// </list>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignProgressive(
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
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Step 1: pairwise NW alignments -> distance matrix (1 - fractional identity).
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        }

        // Step 2: UPGMA guide tree.
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);

        // Step 3: progressive alignment along the tree (post-order: tips first).
        Profile aligned = AlignProfileSubtree(root, rawSeqs, effectiveScoring);

        // Reproject rows back into input order (the tree visited leaves in tree order).
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Pairwise distance for the guide tree: 1 − (identical columns / pairwise NW alignment
    /// length). Two empty sequences (no comparable column) are treated as distance 0.
    /// Source: Feng &amp; Doolittle (1987) convert pairwise alignment to a distance; the
    /// identity-based distance d = 1 − fractional identity is the standard Clustal-style choice.
    /// </summary>
    private static double PairwiseIdentityDistance(string a, string b, ScoringMatrix scoring)
    {
        if (a.Length == 0 && b.Length == 0)
            return 0.0;
        if (a.Length == 0 || b.Length == 0)
            return 1.0; // no shared residues possible

        var pa = GlobalAlignCore(a, b, scoring);
        int length = pa.AlignedSequence1.Length;
        if (length == 0)
            return 0.0;

        int identical = 0;
        for (int i = 0; i < length; i++)
        {
            char c1 = pa.AlignedSequence1[i];
            char c2 = pa.AlignedSequence2[i];
            if (c1 != '-' && c2 != '-' && c1 == c2)
                identical++;
        }

        return 1.0 - (double)identical / length;
    }

    /// <summary>
    /// A node of the progressive-alignment guide tree. Leaves carry a single input-sequence
    /// index; internal nodes carry two children. (Distinct from PhylogeneticAnalyzer.PhyloNode,
    /// which lives in a project that already depends on this one.)
    /// </summary>
    private sealed class ProgressiveGuideNode
    {
        public int LeafIndex = -1; // >= 0 for leaves, -1 for internal nodes
        public ProgressiveGuideNode? Left;
        public ProgressiveGuideNode? Right;
        public bool IsLeaf => Left == null && Right == null;
    }

    /// <summary>
    /// Builds a UPGMA guide tree over the pairwise distance matrix. Uses proportional
    /// (cluster-size weighted) averaging, identical to PhylogeneticAnalyzer.BuildUPGMA
    /// (Sokal &amp; Michener 1958). Tie-breaking is deterministic: the lowest (i, j) index
    /// pair wins, so hand-derived test cases have a unique tree.
    /// </summary>
    private static ProgressiveGuideNode BuildProgressiveGuideTree(int k, double[,] distance)
    {
        var nodes = new Dictionary<int, ProgressiveGuideNode>();
        var sizes = new Dictionary<int, int>();
        var dist = new Dictionary<(int, int), double>();

        for (int i = 0; i < k; i++)
        {
            nodes[i] = new ProgressiveGuideNode { LeafIndex = i };
            sizes[i] = 1;
        }
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                dist[(i, j)] = distance[i, j];
                dist[(j, i)] = distance[i, j];
            }

        // Keep active clusters in a sorted list so tie-breaking is by ascending index.
        var active = Enumerable.Range(0, k).ToList();

        while (active.Count > 1)
        {
            double minDist = double.MaxValue;
            int minI = -1, minJ = -1;
            for (int ii = 0; ii < active.Count; ii++)
            {
                for (int jj = ii + 1; jj < active.Count; jj++)
                {
                    int a = active[ii];
                    int b = active[jj];
                    var key = a < b ? (a, b) : (b, a);
                    double d = dist.GetValueOrDefault(key, 0);
                    if (d < minDist)
                    {
                        minDist = d;
                        minI = a;
                        minJ = b;
                    }
                }
            }

            int newSize = sizes[minI] + sizes[minJ];
            var newNode = new ProgressiveGuideNode
            {
                Left = nodes[minI],
                Right = nodes[minJ]
            };

            foreach (int c in active)
            {
                if (c == minI || c == minJ) continue;
                var keyIK = minI < c ? (minI, c) : (c, minI);
                var keyJK = minJ < c ? (minJ, c) : (c, minJ);
                double dIK = dist.GetValueOrDefault(keyIK, 0);
                double dJK = dist.GetValueOrDefault(keyJK, 0);
                double newDist = (dIK * sizes[minI] + dJK * sizes[minJ]) / newSize;
                var newKey = minI < c ? (minI, c) : (c, minI);
                dist[newKey] = newDist;
            }

            nodes[minI] = newNode;
            sizes[minI] = newSize;
            nodes.Remove(minJ);
            sizes.Remove(minJ);
            active.Remove(minJ);
            // active stays ascending: minJ removed, minI keeps its (smaller) slot.
        }

        return nodes[active[0]];
    }

    /// <summary>
    /// A profile: a block of aligned rows of equal length, each tagged with the input-sequence
    /// index it represents (so the final output can be reordered back to input order).
    /// </summary>
    private sealed class Profile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    /// <summary>
    /// Post-order traversal of the guide tree: align each subtree's profile, combining child
    /// profiles by profile-profile NW alignment. Leaves start as a single-row profile.
    /// </summary>
    private static Profile AlignProfileSubtree(
        ProgressiveGuideNode node, string[] rawSeqs, ScoringMatrix scoring)
    {
        if (node.IsLeaf)
        {
            var leaf = new Profile();
            leaf.Rows.Add(rawSeqs[node.LeafIndex]);
            leaf.SequenceIndices.Add(node.LeafIndex);
            return leaf;
        }

        var left = AlignProfileSubtree(node.Left!, rawSeqs, scoring);
        var right = AlignProfileSubtree(node.Right!, rawSeqs, scoring);
        return AlignProfiles(left, right, scoring);
    }

    /// <summary>
    /// Aligns two profiles with the Needleman-Wunsch recurrence over columns. The column-vs-column
    /// substitution score is the average match/mismatch over all cross-profile residue pairs
    /// (sum-of-pairs profile scoring); a gap column scores GapExtend per residue averaged.
    /// Existing columns are never edited — gaps are inserted as whole new all-gap columns into one
    /// profile, enforcing "once a gap, always a gap" (Feng &amp; Doolittle 1987).
    /// </summary>
    private static Profile AlignProfiles(Profile p1, Profile p2, ScoringMatrix scoring)
    {
        int m = p1.Width;
        int n = p2.Width;

        // F(i,j): best score aligning first i columns of p1 with first j columns of p2.
        var f = new double[m + 1, n + 1];
        for (int i = 1; i <= m; i++)
            f[i, 0] = f[i - 1, 0] + ProfileColumnVsGap(p1, i - 1, p2.Rows.Count, scoring);
        for (int j = 1; j <= n; j++)
            f[0, j] = f[0, j - 1] + ProfileColumnVsGap(p2, j - 1, p1.Rows.Count, scoring);

        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + ProfileColumnScore(p1, i - 1, p2, j - 1, scoring);
                double up = f[i - 1, j] + ProfileColumnVsGap(p1, i - 1, p2.Rows.Count, scoring);
                double left = f[i, j - 1] + ProfileColumnVsGap(p2, j - 1, p1.Rows.Count, scoring);
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Traceback, building column lists for both profiles.
        var cols1 = new List<int>(); // column index in p1, or -1 for a gap column
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + ProfileColumnScore(p1, ii - 1, p2, jj - 1, scoring);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0)
            {
                double up = f[ii - 1, jj] + ProfileColumnVsGap(p1, ii - 1, p2.Rows.Count, scoring);
                if (jj == 0 || AreClose(f[ii, jj], up))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(-1);
                    ii--;
                    continue;
                }
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new Profile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendReprojectedRows(merged, p1, cols1);
        AppendReprojectedRows(merged, p2, cols2);
        return merged;
    }

    /// <summary>
    /// Re-emits each row of <paramref name="source"/> according to <paramref name="columnMap"/>:
    /// a non-negative entry copies that existing column (gaps in it preserved), a -1 inserts a
    /// fresh gap. Existing columns are copied verbatim ("once a gap, always a gap").
    /// </summary>
    private static void AppendReprojectedRows(Profile dest, Profile source, List<int> columnMap)
    {
        foreach (var row in source.Rows)
        {
            var sb = new StringBuilder(columnMap.Count);
            foreach (int col in columnMap)
                sb.Append(col < 0 ? '-' : row[col]);
            dest.Rows.Add(sb.ToString());
        }
    }

    /// <summary>Profile-profile column score: average match/mismatch/gap over all cross pairs.</summary>
    private static double ProfileColumnScore(Profile p1, int col1, Profile p2, int col2, ScoringMatrix scoring)
    {
        double total = 0;
        int pairs = p1.Rows.Count * p2.Rows.Count;
        for (int r1 = 0; r1 < p1.Rows.Count; r1++)
        {
            char a = p1.Rows[r1][col1];
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                char b = p2.Rows[r2][col2];
                if (a == '-' && b == '-')
                    total += 0; // gap-gap neutral
                else if (a == '-' || b == '-')
                    total += scoring.GapExtend;
                else
                    total += a == b ? scoring.Match : scoring.Mismatch;
            }
        }
        return total / pairs;
    }

    /// <summary>
    /// Score for aligning one profile column against an all-gap column of the other profile of
    /// <paramref name="otherRowCount"/> rows: average gap penalty over all cross pairs.
    /// Residue-vs-gap costs GapExtend; gap-vs-gap is neutral.
    /// </summary>
    private static double ProfileColumnVsGap(Profile p, int col, int otherRowCount, ScoringMatrix scoring)
    {
        double total = 0;
        int pairs = p.Rows.Count * otherRowCount;
        for (int r = 0; r < p.Rows.Count; r++)
        {
            char a = p.Rows[r][col];
            // The other side is a gap for all otherRowCount rows.
            if (a != '-')
                total += scoring.GapExtend * otherRowCount;
            // a == '-' : gap-gap, neutral.
        }
        return total / pairs;
    }

    private const double ProfileScoreEpsilon = 1e-9;

    private static bool AreClose(double x, double y) => Math.Abs(x - y) <= ProfileScoreEpsilon;

    #endregion

    #region Multiple Sequence Alignment (Iterative Refinement — MUSCLE tree-dependent restricted partitioning)

    // Iterative refinement of a progressive MSA, added as a THIRD aligner alongside the star
    // (MultipleAlign) and progressive (MultipleAlignProgressive) methods. Both of those remain
    // byte-for-byte unchanged. This removes the single-pass "once a gap, always a gap" limitation
    // of the progressive seed: an early gap-placement error CAN now be corrected, because the
    // alignment is repeatedly re-split and re-aligned and a strictly better arrangement replaces it.
    //
    // Algorithm — MUSCLE Stage 3, "tree-dependent restricted partitioning" (Edgar 2004):
    //   Quoting Edgar (2004), Nucleic Acids Res 32(5):1792-1797, §"Stage 3, Refinement":
    //   3.1 "An edge is chosen from TREE2 (edges are visited in order of decreasing distance from
    //        the root). TREE2 is divided into two subtrees by deleting the edge."
    //   3.2/3.3 "The profile of the multiple alignment in each subtree is computed. A new multiple
    //        alignment is produced by re-aligning the two profiles."
    //   3.4 "If the SP score is improved, the new alignment is kept, otherwise it is discarded."
    //   "Steps 3.1-3.4 are repeated until convergence or until a user-defined limit is reached."
    //
    //   The same scheme also realises the Barton-Sternberg (1987) idea of returning to an existing
    //   alignment and re-aligning a part of it against the rest; here the partition is the
    //   guide-tree edge rather than a single removed sequence, which is exactly Edgar's restricted
    //   partitioning. Re-aligning the two sub-profiles uses the EXISTING profile-profile NW
    //   (AlignProfiles) and is accepted only on a non-decreasing sum-of-pairs (SP) score, so the
    //   refined alignment is provably never worse than the progressive seed.
    //
    //   SP score: "the MSA program optimizes the sum of all of the pairs of characters at each
    //   position in the alignment (the so-called sum of pair score)" (Wikipedia "Multiple sequence
    //   alignment"). Computed by the existing ComputeSumOfPairsScore (column-based: match/mismatch
    //   from the matrix, residue-gap = GapExtend, gap-gap neutral).
    //
    // Determinism: edges are enumerated in a fixed order (decreasing distance from the root, i.e.
    // internal nodes nearest the leaves first, with a deterministic tie-break); profile splitting,
    // re-projection and the accept-on-improvement rule are deterministic. No RNG is used.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Edgar RC (2004) "MUSCLE: multiple sequence alignment with high accuracy and high
    //     throughput", Nucleic Acids Res 32(5):1792-1797.
    //     https://academic.oup.com/nar/article/32/5/1792/2380623 — Stage 3 steps 3.1-3.4 quoted above.
    //   - Barton GJ, Sternberg MJ (1987) "A strategy for the rapid multiple alignment of protein
    //     sequences...", J Mol Biol 198(2):327-337. https://pubmed.ncbi.nlm.nih.gov/3430611/ —
    //     iterative refinement: re-align a part of the alignment to the rest, iterate to a final
    //     alignment.
    //   - Wikipedia "Multiple sequence alignment", §Iterative methods / sum-of-pairs score
    //     https://en.wikipedia.org/wiki/Multiple_sequence_alignment — SP-score definition; iterative
    //     methods "can return to previously calculated ... sub-MSAs ... optimizing a general
    //     objective function".

    // Default cap on full refinement passes over all edges (Edgar's "user-defined limit"). A pass
    // that makes no accepted change means convergence, so this is only an upper bound.
    private const int DefaultMaxRefinementIterations = 16;

    /// <summary>
    /// Performs <b>iterative refinement</b> of a progressive multiple sequence alignment using
    /// MUSCLE-style tree-dependent restricted partitioning (Edgar 2004, Stage 3), removing the
    /// single-pass "once a gap, always a gap" limitation of
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Edgar 2004, Nucleic Acids Res 32(5):1792-1797; Barton &amp; Sternberg 1987):
    /// </para>
    /// <list type="number">
    /// <item>Build the progressive (Feng-Doolittle / UPGMA) seed alignment and keep its guide tree.</item>
    /// <item>Visit each internal guide-tree edge (leaves-first, deterministic order). Deleting the
    /// edge partitions the sequences into two groups.</item>
    /// <item>Project the current alignment onto each group, drop columns that became all-gap, and
    /// realign the two sub-profiles with the existing profile-profile Needleman-Wunsch.</item>
    /// <item>Accept the re-alignment only if the sum-of-pairs (SP) score does not decrease.</item>
    /// <item>Repeat full passes until no edge yields an improvement (convergence) or the iteration
    /// cap is reached.</item>
    /// </list>
    /// <para>
    /// The result is therefore guaranteed to have an SP score <b>≥</b> that of the progressive seed,
    /// and is deterministic (no RNG; fixed edge order). The star
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and progressive
    /// aligners are unchanged.
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="maxIterations">
    /// Upper bound on full refinement passes (Edgar's "user-defined limit"); must be positive.
    /// Defaults to <see cref="DefaultMaxRefinementIterations"/>. Convergence usually stops earlier.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="maxIterations"/> &lt; 1.</exception>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);

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
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Step 1: progressive seed (same pipeline as MultipleAlignProgressive) + retain guide tree.
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }

        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);
        Profile current = AlignProfileSubtree(root, rawSeqs, effectiveScoring);

        // Step 2-5: iterative refinement over the guide-tree edges.
        current = RefineByTreePartitioning(current, root, effectiveScoring, maxIterations);

        // Reproject rows back into input order.
        var orderedRows = new string[k];
        for (int r = 0; r < current.Rows.Count; r++)
            orderedRows[current.SequenceIndices[r]] = current.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Repeatedly splits the current alignment at each internal guide-tree edge into two
    /// sub-profiles, realigns them with profile-profile NW, and keeps the result only when the SP
    /// score does not decrease (Edgar 2004 Stage 3, steps 3.1-3.4). Iterates full passes until a
    /// pass accepts no change (convergence) or <paramref name="maxIterations"/> is reached.
    /// </summary>
    private static Profile RefineByTreePartitioning(
        Profile current, ProgressiveGuideNode root, ScoringMatrix scoring, int maxIterations)
    {
        // 3.1 enumerate the partitions induced by deleting each internal edge of the guide tree.
        // An internal edge is the edge above each non-root internal node; deleting it isolates that
        // node's leaf set as one group and the remaining leaves as the other. Edges are ordered by
        // decreasing distance from the root (deepest internal nodes — nearest the leaves — first),
        // with a deterministic tie-break on the sorted leaf-set, matching Edgar's visiting order.
        var partitions = EnumerateEdgePartitions(root);
        if (partitions.Count == 0)
            return current; // 2 leaves: the single edge is the root split, no internal edge to refine.

        // The SP score the existing pipeline reports for this alignment (input-order independent).
        int CurrentSpScore(Profile p) => ComputeSumOfPairsScore(p.Rows, scoring);

        int bestScore = CurrentSpScore(current);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            bool improvedThisPass = false;

            foreach (var group in partitions)
            {
                // 3.1 split the current alignment into the two leaf groups.
                var groupSet = group;
                var sideA = SplitProfile(current, idx => groupSet.Contains(idx));
                var sideB = SplitProfile(current, idx => !groupSet.Contains(idx));
                if (sideA.Rows.Count == 0 || sideB.Rows.Count == 0)
                    continue;

                // 3.2/3.3 realign the two profiles with the existing profile-profile NW.
                Profile candidate = AlignProfiles(sideA, sideB, scoring);
                int candidateScore = CurrentSpScore(candidate);

                // 3.4 keep only if the SP score does not decrease. Strict ">" would also be valid;
                // ">=" lets an equal-score rearrangement settle, but acceptance still never lowers
                // the score, so the monotonic-non-decreasing guarantee holds. To keep convergence
                // well-defined we only treat a STRICT improvement as "made progress".
                if (candidateScore > bestScore)
                {
                    current = candidate;
                    bestScore = candidateScore;
                    improvedThisPass = true;
                }
            }

            if (!improvedThisPass)
                break; // convergence: a full pass changed nothing.
        }

        return current;
    }

    /// <summary>
    /// Returns, for each internal (non-root) node of the guide tree, the set of input-sequence
    /// indices in its subtree — i.e. one side of the partition obtained by deleting that node's
    /// parent edge. Ordered by decreasing distance from the root (deepest first) with a
    /// deterministic tie-break, per Edgar (2004) "edges are visited in order of decreasing distance
    /// from the root".
    /// </summary>
    private static List<HashSet<int>> EnumerateEdgePartitions(ProgressiveGuideNode root)
    {
        var result = new List<(int Depth, List<int> Sorted, HashSet<int> Set)>();

        void Visit(ProgressiveGuideNode node, int depth, bool isRoot)
        {
            if (node.IsLeaf)
                return;

            // The root edge is not an internal edge (deleting it gives the trivial whole-vs-empty
            // split handled by the seed); only record proper internal nodes.
            if (!isRoot)
            {
                var leaves = new List<int>();
                CollectLeaves(node, leaves);
                leaves.Sort();
                result.Add((depth, leaves, new HashSet<int>(leaves)));
            }

            Visit(node.Left!, depth + 1, false);
            Visit(node.Right!, depth + 1, false);
        }

        Visit(root, 0, true);

        // Decreasing distance from the root = larger depth first; tie-break deterministically on the
        // sorted leaf list so the iteration order is fully reproducible.
        result.Sort((a, b) =>
        {
            int byDepth = b.Depth.CompareTo(a.Depth);
            if (byDepth != 0) return byDepth;
            return CompareIntLists(a.Sorted, b.Sorted);
        });

        return result.Select(r => r.Set).ToList();
    }

    private static int CompareIntLists(List<int> a, List<int> b)
    {
        int n = Math.Min(a.Count, b.Count);
        for (int i = 0; i < n; i++)
        {
            int c = a[i].CompareTo(b[i]);
            if (c != 0) return c;
        }
        return a.Count.CompareTo(b.Count);
    }

    private static void CollectLeaves(ProgressiveGuideNode node, List<int> leaves)
    {
        if (node.IsLeaf)
        {
            leaves.Add(node.LeafIndex);
            return;
        }
        CollectLeaves(node.Left!, leaves);
        CollectLeaves(node.Right!, leaves);
    }

    /// <summary>
    /// Projects the current alignment onto the rows whose input-sequence index satisfies
    /// <paramref name="keep"/>, then drops every column that is all-gap within the selected rows.
    /// The result is a valid profile of the selected subset; gaps inside retained columns are
    /// preserved (no existing column is edited — "once a gap, always a gap" within the sub-profile).
    /// </summary>
    private static Profile SplitProfile(Profile current, Func<int, bool> keep)
    {
        var rowIdx = new List<int>();
        for (int r = 0; r < current.Rows.Count; r++)
            if (keep(current.SequenceIndices[r]))
                rowIdx.Add(r);

        var sub = new Profile();
        if (rowIdx.Count == 0)
            return sub;

        int width = current.Width;
        // Identify columns that are NOT all-gap within the selected rows.
        var keptCols = new List<int>(width);
        for (int c = 0; c < width; c++)
        {
            bool allGap = true;
            foreach (int r in rowIdx)
            {
                if (current.Rows[r][c] != '-')
                {
                    allGap = false;
                    break;
                }
            }
            if (!allGap)
                keptCols.Add(c);
        }

        foreach (int r in rowIdx)
        {
            var sb = new StringBuilder(keptCols.Count);
            foreach (int c in keptCols)
                sb.Append(current.Rows[r][c]);
            sub.Rows.Add(sb.ToString());
            sub.SequenceIndices.Add(current.SequenceIndices[r]);
        }

        return sub;
    }

    #endregion

    #region Multiple Sequence Alignment (Consistency-based — T-Coffee)

    // Consistency-based multiple sequence alignment, the T-Coffee method (Notredame, Higgins &
    // Heringa 2000), added as a FOURTH aligner alongside the star (MultipleAlign), progressive
    // (MultipleAlignProgressive) and iterative (MultipleAlignIterative) methods. All three of those
    // remain byte-for-byte unchanged. This optimises the T-Coffee CONSISTENCY objective — a distinct
    // objective class from the fixed-matrix sum-of-pairs (SP) score the other methods optimise.
    //
    // Algorithm — Notredame, Higgins & Heringa (2000) "T-Coffee: A novel method for fast and accurate
    // multiple sequence alignment", J Mol Biol 302:205-217 (DOI 10.1006/jmbi.2000.4042). Figure 1
    // pipeline: primary library -> extension -> extended library -> progressive alignment.
    //
    //   1. PRIMARY LIBRARY (p.207). For every pair of sequences (i,j) compute pairwise alignments —
    //      a GLOBAL (Needleman-Wunsch, the ClustalW library in the paper) and a LOCAL (Smith-Waterman,
    //      the Lalign library). Each aligned residue pair (residue at position p in Si aligned to the
    //      residue at position q in Sj) "receives a weight equal to percent identity within the
    //      pair-wise alignment it comes from". The global and local libraries are combined by
    //      "signal addition": "If any pair is duplicated between the two libraries, it is merged into
    //      a single entry that has a weight equal to the sum of the two weights." Pairs that never
    //      occur have weight 0.
    //
    //   2. LIBRARY EXTENSION — the triplet consistency transformation (pp.208-209). For each library
    //      pair (Si.p, Sj.q), and for every intermediate sequence Sk, if Si.p aligns to Sk.r (weight
    //      W1) and Sk.r aligns to Sj.q (weight W2), the triplet through Sk contributes min(W1, W2):
    //      "we associate that alignment with a weight equal to the minimum of W1 and W2". The extended
    //      weight is "the sum of all the weights gathered through the examination of all the triplets
    //      involving that pair" PLUS the direct primary weight — the paper's worked example: direct
    //      88, through-C min(77,100)=77, extended = 88 + 77 = 165. Uninformative triplets contribute 0,
    //      so extension never lowers a weight.
    //
    //   3. PROGRESSIVE ALIGNMENT (pp.209-210). A UPGMA guide tree (the existing
    //      BuildProgressiveGuideTree; the paper uses NJ — merge order only, same library and objective)
    //      directs a progressive alignment. Profiles are aligned by dynamic programming whose
    //      position-specific column score is the SUM of the EXTENDED-LIBRARY weights over all
    //      cross-profile residue pairs in the two columns (group-vs-group: "the average library scores
    //      in each column of existing alignment are taken" — here summed over the cross pairs, which
    //      is the same objective up to a per-merge constant). No fixed substitution matrix is used and
    //      "gap-opening penalties and gap-extension penalties [are] set to zero" (p.210): a gap column
    //      scores 0. "Once a gap is introduced ... it cannot be shifted later" is enforced as in the
    //      progressive aligner (whole all-gap columns inserted, existing columns never edited).
    //
    // Determinism: no RNG. Primary alignments, the guide tree, the extension sum and the DP tie-breaks
    // are all deterministic.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Notredame C, Higgins DG, Heringa J (2000) J Mol Biol 302:205-217. Full text fetched/read:
    //     https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf (Figures 1-2, pp.207-211).
    //   - T-Coffee Technical Documentation (min() rule for triplet legs):
    //     https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
    //   - Gotoh O (1982) J Mol Biol 162:705-708 (the DP used in the progressive phase).

    /// <summary>
    /// Performs <b>consistency-based</b> multiple sequence alignment — the T-Coffee method of
    /// Notredame, Higgins &amp; Heringa (2000) — optimising the T-Coffee consistency objective rather
    /// than the fixed-matrix sum-of-pairs score used by
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>,
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, int)"/>
    /// (all of which are unchanged).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Notredame et al. 2000, J Mol Biol 302:205-217):
    /// </para>
    /// <list type="number">
    /// <item>Build a <b>primary library</b> from all pairwise global (Needleman-Wunsch) and local
    /// (Smith-Waterman) alignments; each aligned residue pair is weighted by the pairwise percent
    /// identity, and global+local weights are summed (signal addition).</item>
    /// <item><b>Extend</b> the library by the triplet consistency transformation: each residue pair
    /// accumulates, over every intermediate sequence, the minimum of the two connecting weights,
    /// added to the direct weight.</item>
    /// <item><b>Progressively align</b> along a UPGMA guide tree with dynamic programming whose
    /// column score is the summed extended-library weight (no substitution matrix, zero gap
    /// penalty).</item>
    /// </list>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">
    /// Scoring matrix used only to compute the pairwise alignments that seed the library
    /// (default: <see cref="SimpleDna"/>). The progressive phase uses the library weights, not this
    /// matrix, as the substitution score.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignConsistency(
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
        int k = seqList.Count;
        var rawSeqs = seqList.Select(s => s.Sequence).ToArray();

        // Steps 1-2: primary library (global + local, signal-added) then triplet extension.
        var library = BuildExtendedLibrary(rawSeqs, effectiveScoring);

        // Step 3a: UPGMA guide tree over (1 - fractional identity), reusing the progressive pipeline.
        var distance = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                double d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], effectiveScoring);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance);

        // Step 3b: progressive alignment driven by the EXTENDED-LIBRARY weights.
        ConsistencyProfile aligned = AlignConsistencySubtree(root, rawSeqs, library);

        // Reproject rows back to input order.
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen);
        // TotalScore is reported in the same SP currency as the sibling aligners for comparability.
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    // ----- Library representation -------------------------------------------------------------

    // A residue-pair key (Si.posI, Sj.posJ) with i &lt; j, used to store library weights. Stored
    // canonically with the smaller sequence index first so lookups are symmetric.
    private readonly record struct ResiduePairKey(int SeqA, int PosA, int SeqB, int PosB);

    private static ResiduePairKey MakeKey(int s1, int p1, int s2, int p2) =>
        s1 < s2 || (s1 == s2 && p1 <= p2)
            ? new ResiduePairKey(s1, p1, s2, p2)
            : new ResiduePairKey(s2, p2, s1, p1);

    /// <summary>
    /// Builds the T-Coffee <b>extended</b> library: the primary library (global + local pairwise
    /// alignments, weighted by percent identity and combined by signal addition) transformed by the
    /// triplet consistency extension (Notredame et al. 2000, pp.207-209).
    /// </summary>
    private static Dictionary<ResiduePairKey, double> BuildExtendedLibrary(
        string[] seqs, ScoringMatrix scoring)
    {
        int k = seqs.Length;

        // ---- Step 1: primary library (global + local), keyed per (sequence,position) pair. ----
        var primary = new Dictionary<ResiduePairKey, double>();

        // Per-(i,j) adjacency used by the extension: position p in Si -> list of (Sj, q, weight).
        // Built from the primary library after combination so weights are the signal-added totals.
        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                AddAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                    GlobalAlignCore(seqs[i], seqs[j], scoring));

                if (seqs[i].Length > 0 && seqs[j].Length > 0)
                {
                    var local = LocalAlignCore(seqs[i], seqs[j], scoring);
                    AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
                }
            }
        }

        // ---- Step 2: triplet consistency extension. ----
        // Build, per sequence position, the list of partners (other sequence, position, weight) from
        // the primary library, so triplets can be enumerated in O(neighbours) per residue.
        var neighbours = BuildNeighbourIndex(primary);

        // Extended weight = direct primary weight + sum over intermediates Sk of min(W1, W2).
        var extended = new Dictionary<ResiduePairKey, double>(primary);

        foreach (var entry in primary)
        {
            var key = entry.Key; // (SeqA.PosA, SeqB.PosB), SeqA < SeqB
            int sA = key.SeqA, pA = key.PosA, sB = key.SeqB, pB = key.PosB;

            // Every intermediate residue r in some Sk that aligns to BOTH endpoints contributes
            // min(W(A,k), W(k,B)). We iterate the neighbours of A and intersect on (Sk, r).
            if (!neighbours.TryGetValue((sA, pA), out var partnersOfA))
                continue;
            if (!neighbours.TryGetValue((sB, pB), out var partnersOfB))
                continue;

            // Index B's partners by (seq,pos) for O(1) intersection.
            foreach (var partner in partnersOfA)
            {
                int kSeq = partner.Key.Item1;
                int kPos = partner.Key.Item2;
                double w1 = partner.Value;
                if (kSeq == sA || kSeq == sB) continue; // intermediate must be a third sequence
                if (!partnersOfB.TryGetValue((kSeq, kPos), out double w2)) continue;
                double contribution = Math.Min(w1, w2);
                if (contribution <= 0) continue;
                extended[key] = extended.GetValueOrDefault(key) + contribution;
            }
        }

        return extended;
    }

    /// <summary>
    /// Adds a global pairwise alignment to the primary library. Every aligned residue pair (both
    /// non-gap) gets the alignment's percent identity weight; duplicate pairs are signal-added.
    /// </summary>
    private static void AddAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult pa)
    {
        double weight = PercentIdentity(pa.AlignedSequence1, pa.AlignedSequence2);
        if (weight <= 0) return;

        int posI = 0, posJ = 0;
        string a = pa.AlignedSequence1, b = pa.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    /// <summary>
    /// Adds a local pairwise alignment to the primary library. The local alignment covers a segment
    /// of each sequence starting at <c>StartPosition1</c>/<c>StartPosition2</c>; matched residue pairs
    /// in that segment are weighted by the local segment's percent identity and signal-added.
    /// </summary>
    private static void AddLocalAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult la)
    {
        double weight = PercentIdentity(la.AlignedSequence1, la.AlignedSequence2);
        if (weight <= 0) return;

        int posI = la.StartPosition1, posJ = la.StartPosition2;
        string a = la.AlignedSequence1, b = la.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    private static void AddPairWeight(
        Dictionary<ResiduePairKey, double> lib, int s1, int p1, int s2, int p2, double weight)
    {
        var key = MakeKey(s1, p1, s2, p2);
        lib[key] = lib.GetValueOrDefault(key) + weight;
    }

    /// <summary>
    /// Percent identity of a pairwise alignment = 100 × (columns where both residues are equal and
    /// non-gap) / (alignment length). The percent-identity weight per Notredame et al. (2000) p.207.
    /// </summary>
    private static double PercentIdentity(string a, string b)
    {
        int len = a.Length;
        if (len == 0) return 0;
        int identical = 0;
        for (int c = 0; c < len; c++)
            if (a[c] != '-' && a[c] == b[c])
                identical++;
        // PercentIdentityScale: the weight is a percentage (Notredame et al. 2000, p.207).
        return Math.Round(PercentIdentityScale * identical / len);
    }

    // Weight is expressed as a percentage (0-100), matching the paper's integer "Prim Weight" values.
    private const double PercentIdentityScale = 100.0;

    /// <summary>
    /// Inverts the primary library into a per-residue adjacency: (seq,pos) -&gt; {(otherSeq,otherPos):
    /// weight}. Both directions of each stored pair are emitted so the triplet extension can walk from
    /// either endpoint.
    /// </summary>
    private static Dictionary<(int Seq, int Pos), Dictionary<(int Seq, int Pos), double>>
        BuildNeighbourIndex(Dictionary<ResiduePairKey, double> primary)
    {
        var index = new Dictionary<(int, int), Dictionary<(int, int), double>>();

        void Add(int s1, int p1, int s2, int p2, double w)
        {
            if (!index.TryGetValue((s1, p1), out var inner))
            {
                inner = new Dictionary<(int, int), double>();
                index[(s1, p1)] = inner;
            }
            inner[(s2, p2)] = w;
        }

        foreach (var (key, w) in primary)
        {
            Add(key.SeqA, key.PosA, key.SeqB, key.PosB, w);
            Add(key.SeqB, key.PosB, key.SeqA, key.PosA, w);
        }

        return index;
    }

    // ----- Consistency profile + progressive DP ----------------------------------------------

    /// <summary>
    /// A profile for consistency alignment: like <see cref="Profile"/> but every cell also records the
    /// original residue position in its source sequence (or -1 for a gap), so the library weight for a
    /// residue pair can be looked up during the position-specific DP.
    /// </summary>
    private sealed class ConsistencyProfile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        // PositionMaps[r][col] = original residue index in sequence SequenceIndices[r], or -1 for gap.
        public List<int[]> PositionMaps { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    private static ConsistencyProfile AlignConsistencySubtree(
        ProgressiveGuideNode node, string[] rawSeqs, Dictionary<ResiduePairKey, double> library)
    {
        if (node.IsLeaf)
        {
            int idx = node.LeafIndex;
            string seq = rawSeqs[idx];
            var leaf = new ConsistencyProfile();
            leaf.Rows.Add(seq);
            leaf.SequenceIndices.Add(idx);
            var map = new int[seq.Length];
            for (int p = 0; p < seq.Length; p++) map[p] = p;
            leaf.PositionMaps.Add(map);
            return leaf;
        }

        var left = AlignConsistencySubtree(node.Left!, rawSeqs, library);
        var right = AlignConsistencySubtree(node.Right!, rawSeqs, library);
        return AlignConsistencyProfiles(left, right, library);
    }

    /// <summary>
    /// Aligns two consistency profiles with the Needleman-Wunsch recurrence over columns, scoring a
    /// column pair by the SUM of extended-library weights over all cross-profile residue pairs (gap
    /// columns score 0 — zero gap penalty, Notredame et al. 2000 p.210). Existing columns are never
    /// edited; gaps are inserted as whole all-gap columns ("once a gap, always a gap").
    /// </summary>
    private static ConsistencyProfile AlignConsistencyProfiles(
        ConsistencyProfile p1, ConsistencyProfile p2, Dictionary<ResiduePairKey, double> library)
    {
        int m = p1.Width;
        int n = p2.Width;

        var f = new double[m + 1, n + 1];
        // First row/column: aligning a column against a gap scores 0 (zero gap penalty).
        // (f[i,0] = f[0,j] = 0 for all i,j.)

        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + ColumnLibraryScore(p1, i - 1, p2, j - 1, library);
                double up = f[i - 1, j];     // gap in p2: 0
                double left = f[i, j - 1];   // gap in p1: 0
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        var cols1 = new List<int>();
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + ColumnLibraryScore(p1, ii - 1, p2, jj - 1, library);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0 && (jj == 0 || AreClose(f[ii, jj], f[ii - 1, jj])))
            {
                cols1.Add(ii - 1);
                cols2.Add(-1);
                ii--;
                continue;
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new ConsistencyProfile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendConsistencyRows(merged, p1, cols1);
        AppendConsistencyRows(merged, p2, cols2);
        return merged;
    }

    private static void AppendConsistencyRows(
        ConsistencyProfile dest, ConsistencyProfile source, List<int> columnMap)
    {
        for (int r = 0; r < source.Rows.Count; r++)
        {
            string row = source.Rows[r];
            int[] srcMap = source.PositionMaps[r];
            var sb = new StringBuilder(columnMap.Count);
            var newMap = new int[columnMap.Count];
            for (int c = 0; c < columnMap.Count; c++)
            {
                int col = columnMap[c];
                if (col < 0)
                {
                    sb.Append('-');
                    newMap[c] = -1;
                }
                else
                {
                    sb.Append(row[col]);
                    newMap[c] = srcMap[col];
                }
            }
            dest.Rows.Add(sb.ToString());
            dest.PositionMaps.Add(newMap);
        }
    }

    /// <summary>
    /// Column-vs-column score = sum over all cross-profile residue pairs of their extended-library
    /// weight (0 when a residue is a gap or the pair is absent from the library). This realises the
    /// position-specific library scoring of Notredame et al. (2000): the DP maximises the summed
    /// consistency weight rather than a fixed substitution matrix.
    /// </summary>
    private static double ColumnLibraryScore(
        ConsistencyProfile p1, int col1, ConsistencyProfile p2, int col2,
        Dictionary<ResiduePairKey, double> library)
    {
        double total = 0;
        for (int r1 = 0; r1 < p1.Rows.Count; r1++)
        {
            int posA = p1.PositionMaps[r1][col1];
            if (posA < 0) continue; // gap contributes 0
            int seqA = p1.SequenceIndices[r1];
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                int posB = p2.PositionMaps[r2][col2];
                if (posB < 0) continue;
                int seqB = p2.SequenceIndices[r2];
                var key = MakeKey(seqA, posA, seqB, posB);
                if (library.TryGetValue(key, out double w))
                    total += w;
            }
        }
        return total;
    }

    /// <summary>
    /// Test/diagnostic accessor: builds the T-Coffee extended library for the given sequences and
    /// returns the extended weight of the residue pair (Si.posI, Sj.posJ), and that pair's primary
    /// (pre-extension) weight, so the consistency transformation can be verified directly
    /// (Notredame et al. 2000, p.209).
    /// </summary>
    internal static (double Primary, double Extended) GetLibraryWeights(
        string[] seqs, int seqI, int posI, int seqJ, int posJ, ScoringMatrix? scoring = null)
    {
        var effectiveScoring = scoring ?? SimpleDna;
        var primaryOnly = BuildPrimaryLibraryForDiagnostics(seqs, effectiveScoring);
        var extended = BuildExtendedLibrary(seqs, effectiveScoring);
        var key = MakeKey(seqI, posI, seqJ, posJ);
        return (primaryOnly.GetValueOrDefault(key), extended.GetValueOrDefault(key));
    }

    private static Dictionary<ResiduePairKey, double> BuildPrimaryLibraryForDiagnostics(
        string[] seqs, ScoringMatrix scoring)
    {
        var primary = new Dictionary<ResiduePairKey, double>();
        int k = seqs.Length;
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                AddAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                    GlobalAlignCore(seqs[i], seqs[j], scoring));
                if (seqs[i].Length > 0 && seqs[j].Length > 0)
                    AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                        LocalAlignCore(seqs[i], seqs[j], scoring));
            }
        return primary;
    }

    #endregion
}

