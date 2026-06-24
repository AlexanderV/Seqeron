using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Finds various types of repeats in DNA sequences including microsatellites (STRs),
/// minisatellites (VNTRs), inverted repeats, and direct repeats.
/// </summary>
public static class RepeatFinder
{
    #region Microsatellite (STR) Detection

    /// <summary>
    /// Finds microsatellites (Short Tandem Repeats) in a DNA sequence.
    /// STRs are 1-6 bp motifs repeated consecutively.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length (default: 1).</param>
    /// <param name="maxUnitLength">Maximum repeat unit length (default: 6).</param>
    /// <param name="minRepeats">Minimum number of repeats to report (default: 3).</param>
    /// <returns>Collection of microsatellite repeats found.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minUnitLength < 1) throw new ArgumentOutOfRangeException(nameof(minUnitLength));
        if (maxUnitLength < minUnitLength) throw new ArgumentOutOfRangeException(nameof(maxUnitLength));
        if (minRepeats < 2) throw new ArgumentOutOfRangeException(nameof(minRepeats));

        return FindMicrosatellitesCore(sequence.Sequence, minUnitLength, maxUnitLength, minRepeats);
    }

    /// <summary>
    /// Finds microsatellites with cancellation support.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length.</param>
    /// <param name="maxUnitLength">Maximum repeat unit length.</param>
    /// <param name="minRepeats">Minimum number of repeats.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>Collection of microsatellite repeats found.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindMicrosatellitesCancellable(
            sequence.Sequence, minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        if (minUnitLength < 1) throw new ArgumentOutOfRangeException(nameof(minUnitLength));
        if (maxUnitLength < minUnitLength) throw new ArgumentOutOfRangeException(nameof(maxUnitLength));
        if (minRepeats < 2) throw new ArgumentOutOfRangeException(nameof(minRepeats));

        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindMicrosatellitesCore(sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats))
            yield return result;
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string with cancellation support.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        if (minUnitLength < 1) throw new ArgumentOutOfRangeException(nameof(minUnitLength));
        if (maxUnitLength < minUnitLength) throw new ArgumentOutOfRangeException(nameof(maxUnitLength));
        if (minRepeats < 2) throw new ArgumentOutOfRangeException(nameof(minRepeats));

        return FindMicrosatellitesCancellable(
            sequence, minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCancellable(
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
        const int checkInterval = 1000;

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            for (int i = 0; i <= seq.Length - unitLen * minRepeats; i++)
            {
                if (processed % checkInterval == 0)
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
                            RepeatType: ClassifyRepeatType(unit));
                    }
                }
            }
        }

        progress?.Report(1.0);
    }

    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCore(
        string seq,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats)
    {
        var reported = new HashSet<(int Start, int End)>();

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            for (int i = 0; i <= seq.Length - unitLen * minRepeats; i++)
            {
                string unit = seq.Substring(i, unitLen);

                // Skip if unit is just repetition of smaller unit
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

                    // Avoid reporting overlapping/contained repeats
                    if (!reported.Any(r => r.Start <= i && r.End >= end))
                    {
                        reported.Add((i, end));
                        yield return new MicrosatelliteResult(
                            Position: i,
                            RepeatUnit: unit,
                            RepeatCount: repeats,
                            TotalLength: repeats * unitLen,
                            RepeatType: ClassifyRepeatType(unit));
                    }
                }
            }
        }
    }

    private static bool IsRedundantUnit(string unit)
    {
        if (unit.Length <= 1) return false;

        // Check if unit is made of smaller repeating pattern
        for (int subLen = 1; subLen < unit.Length; subLen++)
        {
            if (unit.Length % subLen != 0) continue;

            string subUnit = unit.Substring(0, subLen);
            bool isRedundant = true;

            for (int i = subLen; i < unit.Length; i += subLen)
            {
                if (unit.Substring(i, subLen) != subUnit)
                {
                    isRedundant = false;
                    break;
                }
            }

            if (isRedundant) return true;
        }

        return false;
    }

    private static RepeatType ClassifyRepeatType(string unit)
    {
        return unit.Length switch
        {
            1 => RepeatType.Mononucleotide,
            2 => RepeatType.Dinucleotide,
            3 => RepeatType.Trinucleotide,
            4 => RepeatType.Tetranucleotide,
            5 => RepeatType.Pentanucleotide,
            6 => RepeatType.Hexanucleotide,
            _ => RepeatType.Complex
        };
    }

    #endregion

    #region Approximate (Imperfect/Interrupted) Tandem Repeat Detection — TRF model

    // --- Tandem Repeats Finder (Benson 1999) reported alignment-scoring parameters -------------------
    // Benson G (1999) "Tandem repeats finder: a program to analyze DNA sequences", Nucleic Acids Res
    // 27(2):573-580, https://doi.org/10.1093/nar/27.2.573. The TRF README/usage (Benson-Genomics-Lab/TRF)
    // states the recommended parameter set "2 7 7 80 10 50 500" = Match Mismatch Delta PM PI Minscore
    // MaxPeriod, and: "The recomended values for Match Mismatch and Delta are 2, 7, and 7 respectively."
    // The TRF definitions page gives Match weight "+2 in all options here. Mismatch and indel weights
    // (interpreted as negative numbers) are either 3, 5, or 7." Score is a Smith-Waterman style alignment
    // score (sum of column weights) computed by wraparound dynamic programming; a tandem repeat is
    // reported when its score is at least Minscore (Benson 1999: "Only those repeats scoring at least 50
    // with these parameters are reported").

    /// <summary>Match weight per aligned identical column. Benson (1999) recommended Match = +2.</summary>
    private const int TrfMatchWeight = 2;

    /// <summary>Mismatch penalty per substituted column. Benson (1999) recommended Mismatch = 7 (applied negatively).</summary>
    private const int TrfMismatchPenalty = -7;

    /// <summary>Indel (gap) penalty per gap column. Benson (1999) recommended Delta = 7 (applied negatively); TRF uses a flat per-column indel weight.</summary>
    private const int TrfIndelPenalty = -7;

    /// <summary>
    /// Default minimum alignment score to report a tandem repeat. Benson (1999): "Only those repeats
    /// scoring at least 50 ... are reported"; the recommended TRF parameter set uses Minscore = 50.
    /// </summary>
    public const int DefaultApproximateMinScore = 50;

    /// <summary>
    /// TRF flat-indel scoring matrix: Match +2, Mismatch -7, indel -7 per gap column (Benson 1999,
    /// recommended set "2 7 7"). The library aligner charges <see cref="ScoringMatrix.GapExtend"/> per
    /// gap column with no separate open cost, which matches TRF's flat indel weight.
    /// </summary>
    private static readonly ScoringMatrix TrfScoring = new(
        Match: TrfMatchWeight,
        Mismatch: TrfMismatchPenalty,
        GapOpen: TrfIndelPenalty,
        GapExtend: TrfIndelPenalty);

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats using the Tandem Repeats Finder
    /// alignment model (Benson 1999). Unlike <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>
    /// — which detects only PERFECT (exact) tandem tracts — this opt-in detector tolerates substitutions
    /// and indels within the repeat: a candidate pattern of each period is aligned against tandem copies
    /// of itself across the sequence, the consensus pattern is determined by majority rule, and the
    /// resulting alignment yields the reported statistics (period size, copy number, percent matches,
    /// percent indels, consensus, alignment score). A repeat is reported when its alignment score is at
    /// least <paramref name="minScore"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minPeriod">Minimum period (motif) size to consider (default: 1).</param>
    /// <param name="maxPeriod">Maximum period (motif) size to consider (default: 6).</param>
    /// <param name="minScore">Minimum TRF alignment score to report (default: <see cref="DefaultApproximateMinScore"/> = 50, per Benson 1999).</param>
    /// <returns>Non-overlapping approximate tandem repeats, best alignment score first.</returns>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        DnaSequence sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindApproximateTandemRepeatsCore(sequence.Sequence, minPeriod, maxPeriod, minScore);
    }

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats in a raw sequence string using the
    /// Tandem Repeats Finder alignment model (Benson 1999). See
    /// <see cref="FindApproximateTandemRepeats(DnaSequence,int,int,int)"/>.
    /// </summary>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        string sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<ApproximateTandemRepeatResult>();

        return FindApproximateTandemRepeatsCore(sequence.ToUpperInvariant(), minPeriod, maxPeriod, minScore);
    }

    private static IReadOnlyList<ApproximateTandemRepeatResult> FindApproximateTandemRepeatsCore(
        string seq,
        int minPeriod,
        int maxPeriod,
        int minScore)
    {
        if (minPeriod < 1) throw new ArgumentOutOfRangeException(nameof(minPeriod));
        if (maxPeriod < minPeriod) throw new ArgumentOutOfRangeException(nameof(maxPeriod));

        var candidates = new List<ApproximateTandemRepeatResult>();
        if (string.IsNullOrEmpty(seq))
            return candidates;

        // For every starting position and every period, grow a window of tandem copies and score it
        // against the majority-rule consensus by alignment. This is a deterministic, exhaustive
        // substitute for TRF's probabilistic k-tuple seeding (the honest residual).
        for (int period = minPeriod; period <= maxPeriod; period++)
        {
            // A repeat needs at least two contiguous copies (Benson 1999: "two or more contiguous,
            // approximate copies of a pattern").
            const int MinCopiesForRepeat = 2;
            for (int start = 0; start + period * MinCopiesForRepeat <= seq.Length; start++)
            {
                var best = EvaluateApproximateRepeat(seq, start, period, minScore);
                if (best is not null)
                    candidates.Add(best.Value);
            }
        }

        // Report best (highest-scoring) repeats first, suppressing any whose span is contained in an
        // already-accepted higher-scoring repeat.
        var accepted = new List<ApproximateTandemRepeatResult>();
        foreach (var c in candidates.OrderByDescending(r => r.AlignmentScore).ThenBy(r => r.Start).ThenBy(r => r.Period))
        {
            int cEnd = c.Start + c.SpanLength;
            bool contained = accepted.Any(a => a.Start <= c.Start && a.Start + a.SpanLength >= cEnd);
            if (!contained)
                accepted.Add(c);
        }

        return accepted
            .OrderByDescending(r => r.AlignmentScore)
            .ThenBy(r => r.Start)
            .ToList();
    }

    /// <summary>
    /// Grows the tandem window from <paramref name="start"/> with the given <paramref name="period"/>,
    /// determines the consensus by majority rule, aligns the window against tandem copies of the
    /// consensus with TRF scoring, and returns the repeat statistics if the alignment score reaches
    /// <paramref name="minScore"/>. Returns the longest scoring window for this (start, period).
    /// </summary>
    private static ApproximateTandemRepeatResult? EvaluateApproximateRepeat(
        string seq,
        int start,
        int period,
        int minScore)
    {
        // Candidate pattern is the first copy at the window start (Benson 1999: "An initial candidate
        // pattern P is drawn from the sequence").
        ApproximateTandemRepeatResult? best = null;

        // Extend the window one copy at a time; the window length need not be an exact multiple of the
        // period (the trailing copy may be partial / contain indels), so we extend in single-base steps
        // but only evaluate when at least two copies are spanned.
        for (int spanLen = period * 2; start + spanLen <= seq.Length; spanLen++)
        {
            string window = seq.Substring(start, spanLen);

            // Consensus by majority rule over the period-aligned columns of the window.
            string consensus = MajorityConsensus(window, period);

            // Reference = a WHOLE number of tandem copies of the consensus pattern covering the window
            // (TRF aligns the sequence against tandem copies of the pattern). The copy count is rounded
            // up so the reference is at least as long as the window; tiling to a partial trailing copy
            // would inject a spurious end-gap and understate the match percentage.
            int copies = (spanLen + period - 1) / period;
            string reference = TileTo(consensus, copies * period);

            AlignmentResult alignment = SequenceAligner.GlobalAlign(window, reference, TrfScoring);
            var stats = ComputeTrfStatistics(alignment, period, consensus);

            if (stats.AlignmentScore >= minScore &&
                (best is null || stats.AlignmentScore > best.Value.AlignmentScore))
            {
                best = stats with { Start = start, SpanLength = spanLen };
            }
        }

        return best;
    }

    /// <summary>
    /// Determines the consensus pattern of length <paramref name="period"/> by majority rule over the
    /// period-aligned columns of <paramref name="window"/> (Benson 1999: "we determine a consensus
    /// pattern by majority rule from the alignment"). Ties are broken by first-seen base for determinism.
    /// </summary>
    private static string MajorityConsensus(string window, int period)
    {
        var consensus = new char[period];
        for (int col = 0; col < period; col++)
        {
            var counts = new Dictionary<char, int>();
            var order = new List<char>();
            for (int i = col; i < window.Length; i += period)
            {
                char b = window[i];
                if (!counts.TryGetValue(b, out int n))
                {
                    counts[b] = 1;
                    order.Add(b);
                }
                else
                {
                    counts[b] = n + 1;
                }
            }

            char bestBase = order[0];
            int bestCount = counts[bestBase];
            foreach (char b in order)
            {
                if (counts[b] > bestCount)
                {
                    bestBase = b;
                    bestCount = counts[b];
                }
            }
            consensus[col] = bestBase;
        }
        return new string(consensus);
    }

    /// <summary>Tiles <paramref name="pattern"/> head-to-tail until it reaches <paramref name="length"/> characters.</summary>
    private static string TileTo(string pattern, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = pattern[i % pattern.Length];
        return new string(chars);
    }

    /// <summary>
    /// Reads the TRF reported statistics from a column-by-column alignment of the observed window
    /// (sequence 1) against tandem copies of the consensus (sequence 2). Percent matches and percent
    /// indels are each expressed over the total alignment columns ("between adjacent copies overall",
    /// Benson 1999). The alignment score is the library aligner's column-weight sum.
    /// </summary>
    private static ApproximateTandemRepeatResult ComputeTrfStatistics(
        AlignmentResult alignment,
        int period,
        string consensus)
    {
        string a = alignment.AlignedSequence1;
        string b = alignment.AlignedSequence2;
        int columns = a.Length;

        int matches = 0;
        int mismatches = 0;
        int indels = 0;
        for (int i = 0; i < columns; i++)
        {
            if (a[i] == '-' || b[i] == '-')
                indels++;
            else if (a[i] == b[i])
                matches++;
            else
                mismatches++;
        }

        double percentMatches = columns > 0 ? (double)matches / columns * 100.0 : 0.0;
        double percentIndels = columns > 0 ? (double)indels / columns * 100.0 : 0.0;

        // Copy number = aligned repeat length / period (Benson 1999: "Number of copies aligned with the
        // consensus pattern"). The aligned repeat length is the number of observed (non-gap) bases.
        int observedBases = a.Count(c => c != '-');
        double copyNumber = period > 0 ? (double)observedBases / period : 0.0;

        return new ApproximateTandemRepeatResult(
            Start: 0,
            SpanLength: observedBases,
            Period: period,
            ConsensusSize: consensus.Length,
            Consensus: consensus,
            CopyNumber: copyNumber,
            PercentMatches: percentMatches,
            PercentIndels: percentIndels,
            AlignmentScore: alignment.Score);
    }

    #endregion

    #region Inverted Repeat Detection

    /// <summary>
    /// Finds inverted repeats (sequences that are reverse complements of each other).
    /// These can form hairpin/stem-loop structures.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minArmLength">Minimum length of each arm (default: 4).</param>
    /// <param name="maxLoopLength">Maximum loop length between arms (default: 50).</param>
    /// <param name="minLoopLength">Minimum loop length (default: 3).</param>
    /// <returns>Collection of inverted repeats found.</returns>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        DnaSequence sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minArmLength < 2) throw new ArgumentOutOfRangeException(nameof(minArmLength));
        if (minLoopLength < 0) throw new ArgumentOutOfRangeException(nameof(minLoopLength));

        return FindInvertedRepeatsCore(sequence.Sequence, minArmLength, maxLoopLength, minLoopLength);
    }

    /// <summary>
    /// Finds inverted repeats in a raw sequence string.
    /// </summary>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        string sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        if (minArmLength < 2) throw new ArgumentOutOfRangeException(nameof(minArmLength));
        if (minLoopLength < 0) throw new ArgumentOutOfRangeException(nameof(minLoopLength));

        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindInvertedRepeatsCore(sequence.ToUpperInvariant(), minArmLength, maxLoopLength, minLoopLength))
            yield return result;
    }

    private static IEnumerable<InvertedRepeatResult> FindInvertedRepeatsCore(
        string seq,
        int minArmLength,
        int maxLoopLength,
        int minLoopLength)
    {
        var reported = new HashSet<(int, int, int)>();

        for (int i = 0; i <= seq.Length - 2 * minArmLength - minLoopLength; i++)
        {
            for (int armLen = minArmLength; i + armLen <= seq.Length; armLen++)
            {
                string leftArm = seq.Substring(i, armLen);
                string leftArmRevComp = DnaSequence.GetReverseComplementString(leftArm);

                // Search for right arm
                int minJ = i + armLen + minLoopLength;
                int maxJ = Math.Min(i + armLen + maxLoopLength, seq.Length - armLen);

                for (int j = minJ; j <= maxJ; j++)
                {
                    if (j + armLen > seq.Length) break;

                    string rightArm = seq.Substring(j, armLen);

                    if (rightArm == leftArmRevComp)
                    {
                        int loopLength = j - (i + armLen);
                        string loop = loopLength > 0 ? seq.Substring(i + armLen, loopLength) : "";

                        var key = (i, j, armLen);
                        if (!reported.Contains(key))
                        {
                            reported.Add(key);
                            yield return new InvertedRepeatResult(
                                LeftArmStart: i,
                                RightArmStart: j,
                                ArmLength: armLen,
                                LoopLength: loopLength,
                                LeftArm: leftArm,
                                RightArm: rightArm,
                                Loop: loop,
                                CanFormHairpin: loopLength >= 3);
                        }
                    }
                }
            }
        }
    }

    #endregion

    #region Direct Repeat Detection

    /// <summary>
    /// Finds direct repeats (identical sequences appearing multiple times).
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5).</param>
    /// <param name="maxLength">Maximum repeat length (default: 50).</param>
    /// <param name="minSpacing">Minimum spacing between repeats (default: 1).</param>
    /// <returns>Collection of direct repeats found.</returns>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        DnaSequence sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minLength < 2) throw new ArgumentOutOfRangeException(nameof(minLength));
        if (maxLength < minLength) throw new ArgumentOutOfRangeException(nameof(maxLength));

        return FindDirectRepeatsCore(sequence.Sequence, minLength, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds direct repeats in a raw sequence string.
    /// </summary>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        string sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        // Mirror the DnaSequence overload's numeric validation onto the raw-string surface.
        // A degenerate minLength < 2 (e.g. 0) yields a zero-length candidate whose suffix-tree
        // lookup matches EVERY position, blowing the result set up with O(n^2) spurious
        // empty-/single-base "repeats"; that is undisciplined fuzzing failure, so it is rejected
        // here exactly as the typed overload rejects it. Validation is hoisted into an eager
        // wrapper so the exception surfaces at the call, not only on enumeration.
        if (minLength < 2) throw new ArgumentOutOfRangeException(nameof(minLength));
        if (maxLength < minLength) throw new ArgumentOutOfRangeException(nameof(maxLength));

        return FindDirectRepeatsRaw(sequence, minLength, maxLength, minSpacing);
    }

    private static IEnumerable<DirectRepeatResult> FindDirectRepeatsRaw(
        string sequence,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindDirectRepeatsCore(sequence.ToUpperInvariant(), minLength, maxLength, minSpacing))
            yield return result;
    }

    private static IEnumerable<DirectRepeatResult> FindDirectRepeatsCore(
        string seq,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        // Use SuffixTree for efficient O(m+k) pattern matching instead of O(n) per pattern
        var suffixTree = global::SuffixTree.SuffixTree.Build(seq);
        var reported = new HashSet<(int, int, int)>();

        for (int len = minLength; len <= maxLength; len++)
        {
            for (int i = 0; i <= seq.Length - len * 2 - minSpacing; i++)
            {
                string repeat = seq.Substring(i, len);

                // Use SuffixTree.FindAllOccurrences for O(m+k) lookup
                var occurrences = suffixTree.FindAllOccurrences(repeat);

                foreach (int j in occurrences.Where(p => p > i + len - 1 + minSpacing).OrderBy(p => p))
                {
                    var key = (i, j, len);
                    if (!reported.Contains(key))
                    {
                        reported.Add(key);
                        yield return new DirectRepeatResult(
                            FirstPosition: i,
                            SecondPosition: j,
                            RepeatSequence: repeat,
                            Length: len,
                            Spacing: j - i - len);
                    }
                }
            }
        }
    }

    #endregion

    #region Tandem Repeat Summary

    /// <summary>
    /// Gets a summary of all tandem repeats in a sequence.
    /// </summary>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        DnaSequence sequence,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var microsatellites = FindMicrosatellites(sequence, 1, 6, minRepeats).ToList();

        var byType = microsatellites
            .GroupBy(m => m.RepeatType)
            .ToDictionary(g => g.Key, g => g.ToList());

        int totalBases = microsatellites.Sum(m => m.TotalLength);

        // PercentageOfSequence is the fraction of the sequence that is tandem-repeat, so it must use the
        // DISTINCT bases covered by any microsatellite (the union of their spans), not the raw sum of repeat
        // lengths: overlapping repeats (e.g. a homopolymer run also matched as a dinucleotide repeat) would
        // otherwise double-count bases and push the percentage above 100. Covered bases ≤ sequence length, so
        // the percentage is always in [0, 100]. TotalRepeatBases keeps the (possibly overlapping) repeat content.
        long coveredBases = CountCoveredBases(microsatellites, sequence.Length);
        double percentageOfSequence = sequence.Length > 0
            ? (double)coveredBases / sequence.Length * 100
            : 0;

        return new TandemRepeatSummary(
            TotalRepeats: microsatellites.Count,
            TotalRepeatBases: totalBases,
            PercentageOfSequence: percentageOfSequence,
            MononucleotideRepeats: byType.GetValueOrDefault(RepeatType.Mononucleotide)?.Count ?? 0,
            DinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Dinucleotide)?.Count ?? 0,
            TrinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Trinucleotide)?.Count ?? 0,
            TetranucleotideRepeats: byType.GetValueOrDefault(RepeatType.Tetranucleotide)?.Count ?? 0,
            LongestRepeat: microsatellites.OrderByDescending(m => m.TotalLength).FirstOrDefault(),
            MostFrequentUnit: microsatellites
                .GroupBy(m => m.RepeatUnit)
                .OrderByDescending(g => g.Count())
                .FirstOrDefault()?.Key);
    }

    /// <summary>
    /// Counts the number of distinct sequence bases covered by at least one microsatellite, i.e. the length
    /// of the union of the half-open spans <c>[Position, Position + TotalLength)</c>. Spans are clamped to the
    /// sequence length and merged in start order, so the result never exceeds <paramref name="sequenceLength"/>.
    /// </summary>
    private static long CountCoveredBases(
        IReadOnlyList<MicrosatelliteResult> microsatellites, int sequenceLength)
    {
        if (microsatellites.Count == 0)
        {
            return 0;
        }

        var intervals = microsatellites
            .Select(m => (Start: m.Position, End: Math.Min(sequenceLength, m.Position + m.TotalLength)))
            .Where(iv => iv.End > iv.Start)
            .OrderBy(iv => iv.Start)
            .ToList();

        long covered = 0;
        int currentStart = -1;
        int currentEnd = -1;
        foreach (var iv in intervals)
        {
            if (iv.Start > currentEnd)
            {
                // Disjoint from the current run: close it out and start a new one.
                covered += currentEnd - currentStart;
                currentStart = iv.Start;
                currentEnd = iv.End;
            }
            else
            {
                // Overlapping or adjacent: extend the current run.
                currentEnd = Math.Max(currentEnd, iv.End);
            }
        }

        covered += currentEnd - currentStart;
        return covered;
    }

    #endregion

    #region Palindrome Detection

    /// <summary>
    /// Finds palindromic sequences (sequences that read the same 5' to 3' on both strands).
    /// These are recognition sites for many restriction enzymes.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum palindrome length (default: 4, must be even).</param>
    /// <param name="maxLength">Maximum palindrome length (default: 12).</param>
    /// <returns>Collection of palindromes found.</returns>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        DnaSequence sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minLength < 4 || minLength % 2 != 0)
            throw new ArgumentOutOfRangeException(nameof(minLength), "Must be even and >= 4");
        if (maxLength < minLength)
            throw new ArgumentOutOfRangeException(nameof(maxLength));

        return FindPalindromesCore(sequence.Sequence, minLength, maxLength);
    }

    /// <summary>
    /// Finds palindromes in a raw sequence string.
    /// </summary>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        string sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        if (minLength < 4 || minLength % 2 != 0)
            throw new ArgumentOutOfRangeException(nameof(minLength), "Must be even and >= 4");
        if (maxLength < minLength)
            throw new ArgumentOutOfRangeException(nameof(maxLength));

        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindPalindromesCore(sequence.ToUpperInvariant(), minLength, maxLength))
            yield return result;
    }

    private static IEnumerable<PalindromeResult> FindPalindromesCore(
        string seq,
        int minLength,
        int maxLength)
    {
        for (int len = minLength; len <= maxLength; len += 2) // Palindromes must be even length
        {
            for (int i = 0; i <= seq.Length - len; i++)
            {
                string candidate = seq.Substring(i, len);
                string revComp = DnaSequence.GetReverseComplementString(candidate);

                if (candidate == revComp)
                {
                    yield return new PalindromeResult(
                        Position: i,
                        Sequence: candidate,
                        Length: len);
                }
            }
        }
    }

    #endregion
}

/// <summary>
/// Type of repeat unit.
/// </summary>
public enum RepeatType
{
    Mononucleotide,  // A, T, G, C
    Dinucleotide,    // AT, GC, CA, etc.
    Trinucleotide,   // CAG, CGG, etc.
    Tetranucleotide, // GATA, AAAT, etc.
    Pentanucleotide, // AAAAT, etc.
    Hexanucleotide,  // AAAAAG, etc.
    Complex          // Longer than 6bp
}

/// <summary>
/// Result of microsatellite (STR) detection.
/// </summary>
public readonly record struct MicrosatelliteResult(
    int Position,
    string RepeatUnit,
    int RepeatCount,
    int TotalLength,
    RepeatType RepeatType)
{
    /// <summary>
    /// Gets the full repeat sequence.
    /// </summary>
    public string FullSequence => string.Concat(Enumerable.Repeat(RepeatUnit, RepeatCount));
}

/// <summary>
/// Result of approximate (imperfect / interrupted) tandem-repeat detection, following the statistics
/// reported by Tandem Repeats Finder (Benson 1999): period size, copy number, percent matches, percent
/// indels, consensus pattern/size, and alignment score.
/// </summary>
public readonly record struct ApproximateTandemRepeatResult(
    int Start,
    int SpanLength,
    int Period,
    int ConsensusSize,
    string Consensus,
    double CopyNumber,
    double PercentMatches,
    double PercentIndels,
    int AlignmentScore);

/// <summary>
/// Result of inverted repeat detection.
/// </summary>
public readonly record struct InvertedRepeatResult(
    int LeftArmStart,
    int RightArmStart,
    int ArmLength,
    int LoopLength,
    string LeftArm,
    string RightArm,
    string Loop,
    bool CanFormHairpin)
{
    /// <summary>
    /// Total length of the inverted repeat structure.
    /// </summary>
    public int TotalLength => 2 * ArmLength + LoopLength;
}

/// <summary>
/// Result of direct repeat detection.
/// </summary>
public readonly record struct DirectRepeatResult(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    int Length,
    int Spacing);

/// <summary>
/// Result of palindrome detection.
/// </summary>
public readonly record struct PalindromeResult(
    int Position,
    string Sequence,
    int Length);

/// <summary>
/// Summary of tandem repeats in a sequence.
/// </summary>
public readonly record struct TandemRepeatSummary(
    int TotalRepeats,
    int TotalRepeatBases,
    double PercentageOfSequence,
    int MononucleotideRepeats,
    int DinucleotideRepeats,
    int TrinucleotideRepeats,
    int TetranucleotideRepeats,
    MicrosatelliteResult? LongestRepeat,
    string? MostFrequentUnit);
