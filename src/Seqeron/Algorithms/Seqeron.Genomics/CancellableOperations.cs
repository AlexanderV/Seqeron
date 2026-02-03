using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics;

/// <summary>
/// Provides cancellable versions of long-running genomic operations.
/// </summary>
public static class CancellableOperations
{
    private const int DefaultCheckInterval = 1000;

    #region K-mer Counting with Cancellation

    /// <summary>
    /// Counts all k-mers in a sequence with cancellation support.
    /// </summary>
    /// <param name="sequence">The sequence to analyze.</param>
    /// <param name="k">The k-mer length.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Dictionary mapping k-mers to their counts.</returns>
    public static Dictionary<string, int> CountKmers(
        string sequence,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, int>();

        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        if (k > sequence.Length)
            return new Dictionary<string, int>();

        var seq = sequence.AsSpan();
        var counts = new Dictionary<string, int>();
        int total = sequence.Length - k + 1;

        for (int i = 0; i <= sequence.Length - k; i++)
        {
            if (i % DefaultCheckInterval == 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                progress?.Report((double)i / total);
            }

            var kmer = new string(seq.Slice(i, k));
            if (!counts.TryAdd(kmer, 1))
                counts[kmer]++;
        }

        progress?.Report(1.0);
        return counts;
    }

    /// <summary>
    /// Counts k-mers asynchronously with cancellation support.
    /// </summary>
    public static Task<Dictionary<string, int>> CountKmersAsync(
        string sequence,
        int k,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        return Task.Run(() => CountKmers(sequence, k, cancellationToken, progress), cancellationToken);
    }

    #endregion

    #region Approximate Matching with Cancellation

    /// <summary>
    /// Finds approximate matches with cancellation support.
    /// </summary>
    public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
        string sequence,
        string pattern,
        int maxMismatches,
        CancellationToken cancellationToken)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
            yield break;

        if (maxMismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches), "Cannot be negative.");

        var seq = sequence.ToUpperInvariant();
        var pat = pattern.ToUpperInvariant();

        if (pat.Length > seq.Length)
            yield break;

        for (int i = 0; i <= seq.Length - pat.Length; i++)
        {
            if (i % DefaultCheckInterval == 0)
                cancellationToken.ThrowIfCancellationRequested();

            int mismatches = 0;
            var positions = new List<int>();

            for (int j = 0; j < pat.Length && mismatches <= maxMismatches; j++)
            {
                if (seq[i + j] != pat[j])
                {
                    mismatches++;
                    positions.Add(j);
                }
            }

            if (mismatches <= maxMismatches)
            {
                yield return new ApproximateMatchResult(
                    i,
                    seq.Substring(i, pat.Length),
                    mismatches,
                    positions.AsReadOnly(),
                    MismatchType.Substitution
                );
            }
        }
    }

    #endregion

    #region Repeat Finding with Cancellation

    /// <summary>
    /// Finds microsatellites with cancellation support.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        var seq = sequence.ToUpperInvariant();
        var reported = new HashSet<(int Start, int End)>();
        int totalPositions = seq.Length * (maxUnitLength - minUnitLength + 1);
        int processed = 0;

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            for (int i = 0; i <= seq.Length - unitLen * minRepeats; i++)
            {
                if (processed % DefaultCheckInterval == 0)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report((double)processed / totalPositions);
                }
                processed++;

                var unit = seq.Substring(i, unitLen);

                if (IsRedundantUnit(unit))
                    continue;

                int repeats = 1;
                int j = i + unitLen;

                while (j + unitLen <= seq.Length && seq.Substring(j, unitLen) == unit)
                {
                    repeats++;
                    j += unitLen;
                }

                if (repeats >= minRepeats)
                {
                    int end = i + (repeats * unitLen) - 1;

                    bool isContained = false;
                    foreach (var r in reported)
                    {
                        if (r.Start <= i && r.End >= end)
                        {
                            isContained = true;
                            break;
                        }
                    }

                    if (!isContained)
                    {
                        reported.Add((i, end));
                        yield return new MicrosatelliteResult(
                            Position: i,
                            RepeatUnit: unit,
                            RepeatCount: repeats,
                            TotalLength: repeats * unitLen,
                            RepeatType: ClassifyRepeatType(unitLen));
                    }
                }
            }
        }

        progress?.Report(1.0);
    }

    private static bool IsRedundantUnit(string unit)
    {
        if (unit.Length <= 1) return false;

        for (int subLen = 1; subLen < unit.Length; subLen++)
        {
            if (unit.Length % subLen != 0) continue;

            var subUnit = unit.AsSpan(0, subLen);
            bool isRedundant = true;

            for (int i = subLen; i < unit.Length; i += subLen)
            {
                if (!unit.AsSpan(i, subLen).SequenceEqual(subUnit))
                {
                    isRedundant = false;
                    break;
                }
            }

            if (isRedundant) return true;
        }

        return false;
    }

    private static RepeatType ClassifyRepeatType(int unitLength) => unitLength switch
    {
        1 => RepeatType.Mononucleotide,
        2 => RepeatType.Dinucleotide,
        3 => RepeatType.Trinucleotide,
        4 => RepeatType.Tetranucleotide,
        5 => RepeatType.Pentanucleotide,
        6 => RepeatType.Hexanucleotide,
        _ => RepeatType.Complex
    };

    #endregion

    #region Sequence Assembly with Cancellation

    /// <summary>
    /// Finds all overlaps between reads with cancellation support.
    /// </summary>
    public static IReadOnlyList<SequenceAssembler.Overlap> FindAllOverlaps(
        IReadOnlyList<string> reads,
        int minOverlap,
        double minIdentity,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        var overlaps = new List<SequenceAssembler.Overlap>();
        int total = reads.Count * reads.Count;
        int processed = 0;

        for (int i = 0; i < reads.Count; i++)
        {
            for (int j = 0; j < reads.Count; j++)
            {
                if (processed % 100 == 0)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report((double)processed / total);
                }
                processed++;

                if (i == j) continue;

                var overlap = SequenceAssembler.FindOverlap(reads[i], reads[j], minOverlap, minIdentity);
                if (overlap.HasValue)
                {
                    overlaps.Add(new SequenceAssembler.Overlap(
                        i, j, overlap.Value.length, overlap.Value.pos1, overlap.Value.pos2));
                }
            }
        }

        progress?.Report(1.0);
        return overlaps;
    }

    #endregion

    #region Motif Finding with Cancellation

    /// <summary>
    /// Finds degenerate motifs with cancellation support.
    /// </summary>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(
        string sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(motif))
            yield break;

        var seq = sequence.ToUpperInvariant();
        var motifUpper = motif.ToUpperInvariant();

        for (int i = 0; i <= seq.Length - motifUpper.Length; i++)
        {
            if (i % DefaultCheckInterval == 0)
                cancellationToken.ThrowIfCancellationRequested();

            bool matches = true;
            for (int j = 0; j < motifUpper.Length && matches; j++)
            {
                char motifChar = motifUpper[j];
                char seqChar = seq[i + j];

                matches = motifChar switch
                {
                    'A' or 'T' or 'G' or 'C' => motifChar == seqChar,
                    'R' => seqChar == 'A' || seqChar == 'G',
                    'Y' => seqChar == 'C' || seqChar == 'T',
                    'S' => seqChar == 'G' || seqChar == 'C',
                    'W' => seqChar == 'A' || seqChar == 'T',
                    'K' => seqChar == 'G' || seqChar == 'T',
                    'M' => seqChar == 'A' || seqChar == 'C',
                    'B' => seqChar != 'A',
                    'D' => seqChar != 'C',
                    'H' => seqChar != 'G',
                    'V' => seqChar != 'T',
                    'N' => true,
                    _ => motifChar == seqChar
                };
            }

            if (matches)
            {
                yield return new MotifMatch(
                    Position: i,
                    MatchedSequence: seq.Substring(i, motifUpper.Length),
                    Pattern: motifUpper,
                    Score: 1.0);
            }
        }
    }

    #endregion

    #region Alignment with Cancellation

    /// <summary>
    /// Performs global alignment with cancellation support.
    /// Useful for very long sequences.
    /// </summary>
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
        var score = scoring ?? SequenceAligner.SimpleDna;

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

    #endregion
}
