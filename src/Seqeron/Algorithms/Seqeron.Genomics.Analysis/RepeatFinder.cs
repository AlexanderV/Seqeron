using Seqeron.Genomics.Alignment;

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
    /// <remarks>
    /// <para>
    /// For every unit length <c>p</c> in <c>[minUnitLength, maxUnitLength]</c> the detector reports each
    /// <b>maximal perfect run</b> of period <c>p</c> exactly once: a run is the longest interval
    /// <c>S[a..e)</c> with <c>S[x] = S[x+p]</c> for all <c>a ≤ x &lt; e−p</c> that cannot be extended to the
    /// left or right (a maximal repetition, Kolpakov &amp; Kucherov 1999). The run is reported at its
    /// left end <c>a</c> with <c>RepeatUnit = S[a..a+p)</c> (the motif phase at the run start) and
    /// <c>RepeatCount = ⌊(e−a)/p⌋</c> complete copies; a trailing partial copy (&lt; p bases) is not
    /// counted, and rotations of the same run (e.g. <c>TA</c> inside <c>ATATATA</c>) are not reported again.
    /// This is the per-motif-size convention of MISA (Thiel et al. 2003; leftmost regex match
    /// <c>([acgt]{p})\1{k-1,}</c>) and pytrf/Krait (Du et al. 2018: seed at the run start, <c>repeat = length / p</c>).
    /// </para>
    /// <para>
    /// A unit that is itself a repetition of a shorter unit (e.g. <c>ATAT</c>, <c>AA</c>) is not reported —
    /// its run is reported at the primitive unit length (MISA "reject false type motifs"). Only units made of
    /// the unambiguous bases A/C/G/T are reported (MISA <c>[acgt]</c>; pytrf skips <c>N</c>), so runs of
    /// <c>N</c> (assembly gaps) or other symbols are never microsatellites. Runs of different unit lengths may
    /// overlap and are reported independently. Input is case-insensitive; positions are 0-based.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length (default: 1).</param>
    /// <param name="maxUnitLength">Maximum repeat unit length (default: 6).</param>
    /// <param name="minRepeats">Minimum number of complete consecutive copies to report (default: 3).</param>
    /// <returns>Collection of microsatellite repeats found, ordered by unit length then position.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        return FindMicrosatellitesCore(
            sequence.Sequence, minUnitLength, maxUnitLength, minRepeats, CancellationToken.None, null);
    }

    /// <summary>
    /// Finds microsatellites with cancellation support. Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>.
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
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        return FindMicrosatellitesCore(
            sequence.Sequence, minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string (case-insensitive). Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>; <c>null</c>/empty input yields no results.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(
            sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats, CancellationToken.None, null);
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
        ValidateMicrosatelliteParameters(minUnitLength, maxUnitLength, minRepeats);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(
            sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats, cancellationToken, progress);
    }

    private static void ValidateMicrosatelliteParameters(int minUnitLength, int maxUnitLength, int minRepeats)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minUnitLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxUnitLength, minUnitLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(minRepeats, 2);
    }

    /// <summary>
    /// Single scan core (upper-case input): for each unit length, visits only run starts
    /// (left-maximal positions), extends the run to its right-maximal end and reports it once.
    /// O(n) character comparisons per unit length.
    /// </summary>
    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCore(
        string seq,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats,
        CancellationToken cancellationToken,
        IProgress<double>? progress)
    {
        int n = seq.Length;
        long unitLengths = (long)maxUnitLength - minUnitLength + 1;
        double totalPositions = Math.Max(1.0, (double)n * unitLengths);
        int sinceCheck = 0;
        const int checkInterval = 1000;

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            if ((long)unitLen * minRepeats > n)
                break; // no longer unit length can fit minRepeats copies either

            int i = 0;
            while (i + unitLen * minRepeats <= n)
            {
                if (++sinceCheck >= checkInterval)
                {
                    sinceCheck = 0;
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report(((double)(unitLen - minUnitLength) * n + i) / totalPositions);
                }

                // Run start (left-maximal): the period-p run cannot be extended one base to the left.
                // Positions inside a run (rotations / suffixes of it) are skipped.
                if (i > 0 && seq[i - 1] == seq[i - 1 + unitLen])
                {
                    i++;
                    continue;
                }

                // Right-maximal end e (exclusive): extend while S[x] == S[x - p].
                int e = i + unitLen;
                while (e < n && seq[e] == seq[e - unitLen])
                    e++;

                int repeats = (e - i) / unitLen;
                if (repeats >= minRepeats && IsReportableUnit(seq, i, unitLen))
                {
                    string unit = seq.Substring(i, unitLen);
                    yield return new MicrosatelliteResult(
                        Position: i,
                        RepeatUnit: unit,
                        RepeatCount: repeats,
                        TotalLength: repeats * unitLen,
                        RepeatType: ClassifyRepeatType(unit));
                }

                // Every position in (i, e - p] lies inside this run (not left-maximal); the next possible
                // run start of the same period is e - p + 1.
                i = e - unitLen + 1;
            }
        }

        cancellationToken.ThrowIfCancellationRequested();
        progress?.Report(1.0);
    }

    /// <summary>
    /// A unit is reportable when it consists only of A/C/G/T and is primitive (not a power of a shorter word).
    /// </summary>
    private static bool IsReportableUnit(string seq, int start, int unitLen)
    {
        for (int k = start; k < start + unitLen; k++)
        {
            char c = seq[k];
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                return false;
        }

        return !IsRedundantUnit(seq.Substring(start, unitLen));
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
        ArgumentOutOfRangeException.ThrowIfLessThan(minPeriod, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxPeriod, minPeriod);

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

    #region TRF Bernoulli statistical-significance scoring (Benson 1999)

    // --- Tandem Repeats Finder probabilistic / Bernoulli model (Benson 1999) -------------------------
    // Benson G (1999) "Tandem repeats finder: a program to analyze DNA sequences", Nucleic Acids Res
    // 27(2):573-580, https://doi.org/10.1093/nar/27.2.573. TRF detailed description / definitions pages
    // (tandem.bu.edu/trf/trf.desc.html, trf.definitions.html; Benson-Genomics-Lab/TRF README), captured
    // VERBATIM 2026-06-24:
    //   * "We model alignment of two tandem copies of a pattern of length n by a sequence of n
    //      independent Bernoulli trials (coin-tosses)."
    //   * "The probability of success, P(Heads), which we also call PM or matching probability,
    //      represents the average percent identity between the copies."
    //   * "A second probability, PI or indel probability, specifies the average percentage of
    //      insertions and deletions between the copies."
    //   * The reported statistics (definitions page items 5-6) are "Percent of matches between adjacent
    //      copies overall" and "Percent of indels between adjacent copies overall", and the alignment
    //      explanation states the statistics refer to "the matches, mismatches and indels overall between
    //      adjacent copies in the sequence, NOT between the sequence and the consensus pattern."
    //   * Default probabilistic data: "PM=80 and PI=10" ("PM = .80 and PI = .10 by default").
    // Faithfully reproducible here: the Bernoulli match/indel PROBABILITY ESTIMATES (PM, PI) computed
    // between ADJACENT COPIES, and the Bernoulli-mean expected matches PM*d over d aligned positions.
    // NOT reproducible without TRF's non-redistributable simulation tables: the percentile cut-offs of
    // R(d,k,PM) (sum-of-heads, "the largest x such that 95% of the time R(d,k,PM) >= x") and the random
    // walk W(d,PI) distance band used for k-tuple SEEDING — that residual is genome-scale performance.

    /// <summary>Benson (1999) default Bernoulli matching probability PM = 0.80 ("PM = .80 by default").</summary>
    public const double TrfDefaultMatchProbability = 0.80;

    /// <summary>Benson (1999) default Bernoulli indel probability PI = 0.10 ("PI = .10 by default").</summary>
    public const double TrfDefaultIndelProbability = 0.10;

    /// <summary>
    /// Computes the Tandem Repeats Finder Bernoulli-model statistical measures (Benson 1999) for a
    /// detected tandem-repeat tract. Benson models the alignment of two adjacent copies of the pattern
    /// as a sequence of independent Bernoulli trials whose success probability P(Heads) = <c>PM</c>
    /// (matching probability) is "the average percent identity between the copies", with a second
    /// probability <c>PI</c> (indel probability) = "the average percentage of insertions and deletions
    /// between the copies". This method estimates <c>PM</c> and <c>PI</c> from the observed tract by
    /// aligning each pair of ADJACENT copies (TRF: statistics are "between adjacent copies in the
    /// sequence, not between the sequence and the consensus pattern") and counting match / mismatch /
    /// indel columns. The Bernoulli-mean expected number of matches over the aligned positions
    /// (<c>PM × columns</c>) is reported as a significance reference, and the estimate is compared to
    /// Benson's default PM (0.80) so callers can judge whether the tract is at least as conserved as a
    /// "significant" random tandem repeat under the model.
    /// </summary>
    /// <remarks>
    /// This is the opt-in probabilistic measure; <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>
    /// and <see cref="FindApproximateTandemRepeats(DnaSequence,int,int,int)"/> are unchanged. Candidate
    /// k-tuple SEEDING (the R(d,k,PM) sum-of-heads percentile cut-off and the W(d,PI) random-walk band)
    /// is NOT reproduced — it depends on TRF's non-redistributable simulation tables and is a
    /// genome-scale performance heuristic, not a per-repeat statistic.
    /// </remarks>
    /// <param name="repeatTract">The observed tandem-repeat tract (≥ 2 copies of the period).</param>
    /// <param name="period">The repeat period (copy length), ≥ 1.</param>
    /// <param name="expectedMatchProbability">
    /// The Bernoulli PM the tract is assessed against (default <see cref="TrfDefaultMatchProbability"/> = 0.80,
    /// Benson 1999). The tract is flagged <see cref="TandemRepeatBernoulliStatistics.MeetsExpectedMatchProbability"/>
    /// when its estimated PM ≥ this value.
    /// </param>
    /// <returns>The Bernoulli-model statistics for the tract.</returns>
    public static TandemRepeatBernoulliStatistics ComputeBernoulliStatistics(
        string repeatTract,
        int period,
        double expectedMatchProbability = TrfDefaultMatchProbability)
    {
        ArgumentNullException.ThrowIfNull(repeatTract);
        ArgumentOutOfRangeException.ThrowIfLessThan(period, 1);
        if (expectedMatchProbability < 0.0 || expectedMatchProbability > 1.0)
            throw new ArgumentOutOfRangeException(nameof(expectedMatchProbability));

        string tract = repeatTract.ToUpperInvariant();
        if (tract.Length < period * 2)
            throw new ArgumentException(
                "A tandem repeat needs at least two contiguous copies of the period.", nameof(repeatTract));

        // Segment the tract into copies of the period (the last copy may be partial) and align each pair
        // of ADJACENT copies. Each Bernoulli trial is one alignment column between the two adjacent
        // copies: heads = match, tails = mismatch or indel (Benson 1999).
        int matches = 0;
        int mismatches = 0;
        int indels = 0;

        int copyCount = (tract.Length + period - 1) / period;
        for (int c = 0; c + 1 < copyCount; c++)
        {
            int leftStart = c * period;
            int rightStart = (c + 1) * period;
            string left = tract.Substring(leftStart, Math.Min(period, tract.Length - leftStart));
            string right = tract.Substring(rightStart, Math.Min(period, tract.Length - rightStart));

            AlignmentResult pair = SequenceAligner.GlobalAlign(left, right, TrfScoring);
            string a = pair.AlignedSequence1;
            string b = pair.AlignedSequence2;

            // GlobalAlign returns AlignmentResult.Empty (no aligned strings) only when an input copy is
            // empty; with ≥ 2 whole copies both adjacent copies are non-empty, so columns are present.
            int columns = a.Length;
            for (int i = 0; i < columns; i++)
            {
                if (a[i] == '-' || b[i] == '-') indels++;
                else if (a[i] == b[i]) matches++;
                else mismatches++;
            }
        }

        int totalColumns = matches + mismatches + indels;

        // PM = matching probability = average percent identity between adjacent copies (Benson 1999).
        // Heads in the Bernoulli model are matches; PM is the fraction of trials that are heads.
        double matchProbability = totalColumns > 0 ? (double)matches / totalColumns : 0.0;

        // PI = indel probability = average percentage of insertions and deletions between the copies.
        double indelProbability = totalColumns > 0 ? (double)indels / totalColumns : 0.0;

        // Bernoulli-mean expected matches over the aligned positions: E[heads] = PM * d for d trials with
        // success probability PM (the mean of a length-d Bernoulli(PM) sequence). Reported as the
        // significance reference (the expected number of matching positions a random tandem repeat with
        // this match probability would show over the same number of trials).
        double expectedMatches = matchProbability * totalColumns;

        return new TandemRepeatBernoulliStatistics(
            Period: period,
            AdjacentCopyPairs: Math.Max(0, copyCount - 1),
            BernoulliTrials: totalColumns,
            Matches: matches,
            Mismatches: mismatches,
            Indels: indels,
            MatchProbability: matchProbability,
            IndelProbability: indelProbability,
            PercentMatches: matchProbability * 100.0,
            PercentIndels: indelProbability * 100.0,
            ExpectedMatches: expectedMatches,
            MeetsExpectedMatchProbability: matchProbability >= expectedMatchProbability);
    }

    #endregion

    #region Inverted Repeat Detection

    /// <summary>
    /// Finds exact (perfect-stem) inverted repeats: a left arm followed, after a loop, by a right arm
    /// equal to the reverse complement of the left arm. Such structures can form hairpin/stem-loop
    /// (single strand) or cruciform (duplex) structures.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Model (exact, maximal stems).</b> A stem is a triple (<c>LeftArmStart = i</c>,
    /// <c>RightArmStart = j</c>, <c>ArmLength = A</c>) such that base <c>s[i+k]</c> is the Watson–Crick
    /// complement of <c>s[j+A−1−k]</c> for every <c>0 ≤ k &lt; A</c>, with
    /// <c>A ≥ minArmLength</c> and <c>minLoopLength ≤ j − (i + A) ≤ maxLoopLength</c>.
    /// Only stems that are <b>not contained, in both arms, in another such stem</b> are reported. This is
    /// the reporting rule of EMBOSS <c>palindrome</c> (M. Faller; <c>palindrome.c</c>: every stem is extended
    /// inward from its outer pair and a stem that is a subset of an already found stem in both halves —
    /// <c>palindrome_AInB</c> — is dropped, default <c>-overlap Y</c>), run with <c>-nummismatches 0</c>.
    /// Consequently every reported stem is maximal: it cannot be extended outward, and it can be extended
    /// inward only by making the loop shorter than <paramref name="minLoopLength"/>; sub-stems and
    /// "slipped" re-pairings lying inside a longer stem are not reported. With
    /// <paramref name="minLoopLength"/> = 0 the result set equals EMBOSS <c>palindrome</c>
    /// (<c>-nummismatches 0 -overlap Y</c>, <c>-gaplimit = maxLoopLength</c>, <c>-minpallen = minArmLength</c>,
    /// unbounded <c>-maxpallen</c>), cross-checked on random sequences.
    /// </para>
    /// <para>
    /// <b>Pairing.</b> Only the unambiguous bases A, C, G, T pair, via the canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/> (A↔T, C↔G). As in EMBOSS <c>einverted</c>
    /// (which scores a match only for a/c/g/t), N and other IUPAC ambiguity codes never form a pair, so a run
    /// of N is not reported as a stem.
    /// </para>
    /// <para>
    /// <b>Not einverted.</b> This is an exact-stem finder: no mismatches, no gaps (bulges), no score threshold.
    /// For imperfect, score-based inverted repeats use EMBOSS <c>einverted</c> (Durbin &amp; Thierry-Mieg
    /// dynamic programming).
    /// </para>
    /// <para>Coordinates are 0-based; results are ordered by <c>LeftArmStart</c>, then <c>RightArmStart</c>.
    /// Parameters are validated eagerly.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minArmLength">Minimum length of each arm (default: 4, must be ≥ 2).</param>
    /// <param name="maxLoopLength">Maximum loop length between arms (default: 50, must be ≥ <paramref name="minLoopLength"/>).</param>
    /// <param name="minLoopLength">Minimum loop length (default: 3, must be ≥ 0).</param>
    /// <returns>The maximal exact inverted repeats.</returns>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        DnaSequence sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength);

        return FindInvertedRepeatsCore(sequence.Sequence, minArmLength, maxLoopLength, minLoopLength);
    }

    /// <summary>
    /// Finds exact maximal inverted repeats in a raw sequence string (case-insensitive).
    /// Same model and validation as <see cref="FindInvertedRepeats(DnaSequence, int, int, int)"/>;
    /// <c>null</c> or empty input yields no results. Characters other than A/C/G/T never pair.
    /// </summary>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        string sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<InvertedRepeatResult>();

        return FindInvertedRepeatsCore(sequence.ToUpperInvariant(), minArmLength, maxLoopLength, minLoopLength);
    }

    private static void ValidateInvertedRepeatParameters(int minArmLength, int maxLoopLength, int minLoopLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minArmLength, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLoopLength, minLoopLength);
    }

    /// <summary>True when <paramref name="a"/> and <paramref name="b"/> form a Watson–Crick pair (ACGT only).</summary>
    private static bool IsWatsonCrickPair(char a, char b) =>
        a is 'A' or 'C' or 'G' or 'T' && SequenceExtensions.GetComplementBase(a) == b;

    /// <summary>
    /// Maximal exact stems. A stem is identified by its innermost pair (<c>iIn</c>, <c>rIn</c>) with loop
    /// <c>L = rIn − iIn − 1</c>; it is inward-maximal when the next inner pair does not pair or would leave a loop
    /// shorter than <c>minLoop</c>. Each inward-maximal stem is extended outward as far as the bases pair
    /// (EMBOSS palindrome's candidate set with 0 mismatches); then stems contained in both arms of another stem
    /// are removed (<c>palindrome_AInB</c>, see <see cref="IsContainedInShiftedStem"/>).
    /// </summary>
    private static List<InvertedRepeatResult> FindInvertedRepeatsCore(
        string seq,
        int minArmLength,
        int maxLoopLength,
        int minLoopLength)
    {
        int n = seq.Length;
        var results = new List<InvertedRepeatResult>();
        if (n < 2 * minArmLength + minLoopLength)
            return results;

        // Largest loop that can still fit two minimum arms; caps a huge maxLoopLength (no overflow / no hang).
        int maxLoop = (int)Math.Min(maxLoopLength, (long)n - 2 * minArmLength);

        for (int iIn = minArmLength - 1; iIn < n; iIn++)
        {
            for (int loop = minLoopLength; loop <= maxLoop; loop++)
            {
                int rIn = iIn + loop + 1;
                if (rIn + minArmLength - 1 >= n)
                    break;
                if (!IsWatsonCrickPair(seq[iIn], seq[rIn]))
                    continue;
                // Inward-maximal: the next inner pair must fail, unless it would make the loop < minLoop.
                if (loop - 2 >= minLoopLength && IsWatsonCrickPair(seq[iIn + 1], seq[rIn - 1]))
                    continue;

                int arm = 1;
                while (iIn - arm >= 0 && rIn + arm < n && IsWatsonCrickPair(seq[iIn - arm], seq[rIn + arm]))
                    arm++;
                if (arm < minArmLength)
                    continue;

                int iOut = iIn - arm + 1;
                if (IsContainedInShiftedStem(seq, iOut, rIn, arm, loop - minLoopLength))
                    continue;

                results.Add(new InvertedRepeatResult(
                    LeftArmStart: iOut,
                    RightArmStart: rIn,
                    ArmLength: arm,
                    LoopLength: loop,
                    LeftArm: seq.Substring(iOut, arm),
                    RightArm: seq.Substring(rIn, arm),
                    Loop: seq.Substring(iIn + 1, loop),
                    CanFormHairpin: loop >= 3));
            }
        }

        results.Sort(static (x, y) =>
        {
            int c = x.LeftArmStart.CompareTo(y.LeftArmStart);
            return c != 0 ? c : x.RightArmStart.CompareTo(y.RightArmStart);
        });
        return results;
    }

    /// <summary>
    /// True when the stem (<paramref name="left"/>, <paramref name="right"/>, <paramref name="arm"/>) lies, in both
    /// arms, inside a stem on another diagonal. Covering both arms of the stem from a diagonal shifted by
    /// <c>m ≥ 1</c> needs a stem with arm <c>arm + m</c> and loop <c>loop − m</c> that either keeps the left start and
    /// the right start (<c>(left, right, arm + m)</c>) or keeps both arm ends (<c>(left − m, right − m, arm + m)</c>);
    /// every containing stem contains one of these, so testing <c>1 ≤ m ≤ loop − minLoop</c> is exact.
    /// </summary>
    private static bool IsContainedInShiftedStem(string seq, int left, int right, int arm, int maxShift)
    {
        int n = seq.Length;
        for (int m = 1; m <= maxShift; m++)
        {
            int len = arm + m;
            if (right + len <= n && IsExactStem(seq, left, right, len))
                return true;
            if (left - m >= 0 && IsExactStem(seq, left - m, right - m, len))
                return true;
        }
        return false;
    }

    /// <summary>True when s[i+k] pairs with s[j+len−1−k] for all k (checked from the innermost pair out).</summary>
    private static bool IsExactStem(string seq, int i, int j, int len)
    {
        for (int k = len - 1; k >= 0; k--)
        {
            if (!IsWatsonCrickPair(seq[i + k], seq[j + len - 1 - k]))
                return false;
        }
        return true;
    }

    #endregion


    #region Direct Repeat Detection

    /// <summary>
    /// Finds exact direct repeats reported as <b>maximal repeated pairs</b> (Gusfield 1997, §7.12;
    /// Kurtz &amp; Schleiermacher 1999 REPuter forward repeats; MUMmer <c>repeat-match -f</c>,
    /// Kurtz et al. 2004).
    /// </summary>
    /// <remarks>
    /// <para>
    /// A pair of 0-based positions <c>i &lt; j</c> with length <c>L</c> is reported when
    /// <c>S[i..i+L) = S[j..j+L)</c>, the pair is <b>right-maximal</b> (<c>L</c> is the full length of the
    /// common prefix of the suffixes at <c>i</c> and <c>j</c>) and <b>left-maximal</b> (<c>i = 0</c> or
    /// <c>S[i−1] ≠ S[j−1]</c>). Every exact direct repeat is therefore reported exactly once, at its full
    /// extent — its shorter sub-copies (nested windows of the same pair of copies) are not repeated.
    /// This is exactly the forward-strand output of MUMmer <c>repeat-match -f -n minLength</c>
    /// (converted to 0-based positions), then filtered by <paramref name="maxLength"/> and
    /// <paramref name="minSpacing"/>.
    /// </para>
    /// <para>
    /// Only A/C/G/T match (case-insensitive); any other symbol (N, IUPAC codes, gaps) never matches and
    /// terminates a repeat — the MUMmer <c>mummer -n</c> convention ("match only the characters a, c, g,
    /// or t"), consistent with the ACGT-only rule of the other <see cref="RepeatFinder"/> methods.
    /// </para>
    /// <para>
    /// Filters: <c>minLength ≤ L ≤ maxLength</c> — a maximal repeat longer than <paramref name="maxLength"/>
    /// is <b>not</b> reported (it is not truncated into sub-windows; pass a larger <paramref name="maxLength"/>
    /// to see it); <c>Spacing = j − i − L ≥ minSpacing</c>. <paramref name="minSpacing"/> may be negative to
    /// admit overlapping copies (e.g. <c>int.MinValue</c> returns every maximal pair, as repeat-match does);
    /// <c>0</c> admits abutting (tandem) copies.
    /// </para>
    /// <para>
    /// Algorithm: suffix array + Kasai LCP array (shared with <see cref="SequenceComplexity"/>), bottom-up
    /// traversal of the lcp-interval tree with per-left-character position lists (Gusfield 1997 §7.12.3;
    /// Abouelhoda, Kurtz &amp; Ohlebusch 2004, maximal repeated pairs on enhanced suffix arrays):
    /// O(n log² n + z) time for z maximal pairs of length within [minLength, maxLength], O(n) extra space.
    /// Results are ordered by (FirstPosition, SecondPosition); each position pair occurs at most once.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5, must be ≥ 2).</param>
    /// <param name="maxLength">Maximum repeat length (default: 50, must be ≥ <paramref name="minLength"/>).</param>
    /// <param name="minSpacing">Minimum number of bases between the copies (default: 1; negative admits overlap).</param>
    /// <returns>Maximal direct-repeat pairs, sorted by (FirstPosition, SecondPosition).</returns>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        DnaSequence sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        return FindDirectRepeatsCore(sequence.Sequence, minLength, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds maximal exact direct-repeat pairs in a raw sequence string (case-insensitive; non-ACGT symbols
    /// never match). <c>null</c> or empty input yields no results. See
    /// <see cref="FindDirectRepeats(DnaSequence,int,int,int)"/> for the reporting convention.
    /// </summary>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        string sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        // Same numeric validation as the DnaSequence overload, eager (at the call, not on enumeration).
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<DirectRepeatResult>();

        return FindDirectRepeatsCore(sequence.ToUpperInvariant(), minLength, maxLength, minSpacing);
    }

    /// <summary>Left-character class for a suffix start with no matchable left neighbour (p = 0 or non-ACGT).</summary>
    private const int UniqueLeftClass = 4;
    private const int LeftClassCount = 5;

    private static int AcgtCode(char c) => c switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1,
    };

    private static List<DirectRepeatResult> FindDirectRepeatsCore(
        string seq,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        var results = new List<DirectRepeatResult>();
        int n = seq.Length;
        if (n <= minLength)
            return results;

        // Symbols: A/C/G/T → 0..3; every other symbol gets a unique code 4 + p, so it never matches.
        var symbols = new int[n];
        var leftClass = new int[n];
        for (int p = 0; p < n; p++)
        {
            int code = AcgtCode(seq[p]);
            symbols[p] = code >= 0 ? code : 4 + p;
            leftClass[p] = UniqueLeftClass;
        }
        for (int p = 1; p < n; p++)
        {
            if (symbols[p - 1] < 4)
                leftClass[p] = symbols[p - 1];
        }

        int[] sa = SequenceComplexity.BuildSuffixArray(symbols);
        int[] lcp = SequenceComplexity.BuildLcpArray(symbols, sa);

        // Bottom-up lcp-interval traversal. Each interval keeps its suffix positions in one linked list
        // per left-character class; merging a child into its parent at string depth ℓ emits every pair
        // (p from the child, q already in the parent) whose left characters differ (or are undefined):
        // those pairs have LCP exactly ℓ (right-maximal) and are left-maximal.
        var next = new int[n];
        var stackLcp = new int[n + 1];
        var stackHead = new int[(n + 1) * LeftClassCount];
        var stackTail = new int[(n + 1) * LeftClassCount];
        var childHead = new int[LeftClassCount];
        var childTail = new int[LeftClassCount];
        int top = 0;
        stackLcp[0] = 0;
        Array.Fill(stackHead, -1, 0, LeftClassCount);

        for (int r = 1; r <= n; r++)
        {
            // Pending child: the leaf for suffix sa[r − 1].
            Array.Fill(childHead, -1);
            int leaf = sa[r - 1];
            next[leaf] = -1;
            childHead[leftClass[leaf]] = leaf;
            childTail[leftClass[leaf]] = leaf;

            int h = r < n ? lcp[r] : 0;
            while (stackLcp[top] > h)
            {
                MergeDirectRepeatLists(seq, top, childHead, childTail, stackLcp, stackHead, stackTail, next,
                    minLength, maxLength, minSpacing, results);
                for (int c = 0; c < LeftClassCount; c++)
                {
                    childHead[c] = stackHead[top * LeftClassCount + c];
                    childTail[c] = stackTail[top * LeftClassCount + c];
                }
                top--;
            }

            if (stackLcp[top] < h)
            {
                top++;
                stackLcp[top] = h;
                Array.Fill(stackHead, -1, top * LeftClassCount, LeftClassCount);
            }

            MergeDirectRepeatLists(seq, top, childHead, childTail, stackLcp, stackHead, stackTail, next,
                minLength, maxLength, minSpacing, results);
        }

        results.Sort(static (a, b) =>
        {
            int c = a.FirstPosition.CompareTo(b.FirstPosition);
            return c != 0 ? c : a.SecondPosition.CompareTo(b.SecondPosition);
        });
        return results;
    }

    /// <summary>
    /// Emits the maximal pairs between a child interval's lists and the lists already accumulated in the
    /// interval at stack slot <paramref name="node"/>, then concatenates the child lists into the node.
    /// </summary>
    private static void MergeDirectRepeatLists(
        string seq,
        int node,
        int[] childHead,
        int[] childTail,
        int[] stackLcp,
        int[] stackHead,
        int[] stackTail,
        int[] next,
        int minLength,
        int maxLength,
        int minSpacing,
        List<DirectRepeatResult> results)
    {
        int length = stackLcp[node];
        int baseIdx = node * LeftClassCount;

        if (length >= minLength && length <= maxLength)
        {
            for (int a = 0; a < LeftClassCount; a++)
            {
                if (childHead[a] < 0) continue;
                for (int b = 0; b < LeftClassCount; b++)
                {
                    if (a == b && a != UniqueLeftClass) continue; // same left character: not left-maximal
                    int nodeListHead = stackHead[baseIdx + b];
                    if (nodeListHead < 0) continue;

                    for (int p = childHead[a]; p >= 0; p = next[p])
                    {
                        for (int q = nodeListHead; q >= 0; q = next[q])
                        {
                            int i = Math.Min(p, q);
                            int j = Math.Max(p, q);
                            long spacing = (long)j - i - length;
                            if (spacing < minSpacing) continue;
                            results.Add(new DirectRepeatResult(
                                FirstPosition: i,
                                SecondPosition: j,
                                RepeatSequence: seq.Substring(i, length),
                                Length: length,
                                Spacing: (int)spacing));
                        }
                    }
                }
            }
        }

        for (int c = 0; c < LeftClassCount; c++)
        {
            if (childHead[c] < 0) continue;
            int idx = baseIdx + c;
            if (stackHead[idx] < 0)
            {
                stackHead[idx] = childHead[c];
            }
            else
            {
                next[stackTail[idx]] = childHead[c];
            }
            stackTail[idx] = childTail[c];
        }
    }

    #endregion

    #region Tandem Repeat Summary

    /// <summary>
    /// Gets a summary of all perfect microsatellites (1–6 bp units) in a sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The summary aggregates exactly the list returned by
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/> with unit lengths 1–6 and the given
    /// <paramref name="minRepeats"/> (one maximal primitive run per unit length; no second scan), following the
    /// statistics conventions of MISA (Thiel et al. 2003, <c>misa.pl</c> <c>.statistics</c> output) and Krait
    /// (Du et al. 2018, <c>statistics.py</c>):
    /// </para>
    /// <list type="bullet">
    /// <item><description><c>TotalRepeats</c> — number of reported microsatellites (MISA "Total number of identified SSRs").</description></item>
    /// <item><description>Per-class counts for all six unit sizes, mono … hexa (MISA "Distribution to different repeat
    /// type classes": one row per unit size; Krait <c>motifTypeStatis</c>: Mono, Di, Tri, Tetra, Penta, Hexa). The six
    /// counts always sum to <c>TotalRepeats</c>.</description></item>
    /// <item><description><c>TotalRepeatBases</c> — sum of the repeat lengths (Krait "Length (bp)" = <c>SUM(length)</c>);
    /// runs of different unit lengths that overlap are each counted in full.</description></item>
    /// <item><description><c>PercentageOfSequence</c> — percent of the sequence's bases covered by at least one reported
    /// microsatellite: |∪ [Position, Position+TotalLength)| / sequence length × 100, so always in [0, 100]. When no two
    /// reported runs overlap this equals <c>TotalRepeatBases</c> / length × 100 (Krait's relative density expressed as a
    /// percentage; DnaSequence contains only A/C/G/T, so Krait's valid-base denominator equals the length).</description></item>
    /// <item><description><c>LongestRepeat</c> — the run with the largest <c>TotalLength</c>; ties go to the shorter unit,
    /// then the leftmost position (the <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/> order).</description></item>
    /// <item><description><c>MostFrequentUnit</c> — the repeat unit string (as reported, i.e. the motif phase at the run
    /// start, not rotation- or strand-canonicalized; MISA table "Frequency of identified SSR motifs") that occurs in the
    /// most runs; ties go to the unit whose first run appears first in that order.</description></item>
    /// </list>
    /// </remarks>
    /// <param name="sequence">DNA sequence to summarize.</param>
    /// <param name="minRepeats">Minimum number of complete copies (≥ 2; default 3), applied to every unit length.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minRepeats"/> is less than 2.</exception>
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
            PentanucleotideRepeats: byType.GetValueOrDefault(RepeatType.Pentanucleotide)?.Count ?? 0,
            HexanucleotideRepeats: byType.GetValueOrDefault(RepeatType.Hexanucleotide)?.Count ?? 0,
            // MicrosatelliteResult is a struct: FirstOrDefault() on the struct sequence would return default(...)
            // (Position 0, RepeatUnit null) instead of null when nothing is found, so lift to the nullable type first.
            LongestRepeat: microsatellites
                .OrderByDescending(m => m.TotalLength)
                .Select(m => (MicrosatelliteResult?)m)
                .FirstOrDefault(),
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
    /// Finds reverse-complement palindromes (sequences that read the same 5'→3' on both strands, i.e. equal
    /// their own reverse complement — the recognition-site form of most Type II restriction enzymes).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Reports <b>every</b> window <c>S[i..i+len)</c> with <c>minLength ≤ len ≤ maxLength</c> that equals its
    /// reverse complement, including windows nested inside longer palindromes and overlapping ones — the
    /// Rosalind REVP ("Locating Restriction Sites") convention: "the position and length of every reverse
    /// palindrome in the string having length between 4 and 12". This differs from
    /// <see cref="FindInvertedRepeats(DnaSequence,int,int,int)"/>, which reports only maximal stems (EMBOSS
    /// <c>palindrome</c>).
    /// </para>
    /// <para>
    /// A reverse-complement palindrome always has even length (the middle base of an odd window would have to be
    /// its own complement), so only even lengths are scanned; an odd <paramref name="maxLength"/> is effectively
    /// rounded down. Only the unambiguous bases A, C, G, T pair (canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/>): a window containing N, another IUPAC ambiguity
    /// code, U, a gap or any other symbol is never reported (as in the REVP alphabet and EMBOSS
    /// <c>einverted</c>; EMBOSS <c>palindrome</c> likewise rejects all-N stems) — e.g. <c>NNNN</c> or <c>SSSS</c>
    /// equal their symbolic reverse complement but need not be palindromic once resolved (<c>SS</c> may be CC).
    /// </para>
    /// <para>
    /// Positions are 0-based (REVP positions are 1-based: add 1). Results are ordered by <c>Position</c>,
    /// then <c>Length</c> (the order of the REVP sample output). Parameters are validated eagerly on both
    /// overloads. Cost O(n·maxLength) comparisons (one capped centre expansion per inter-base gap) plus output.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum palindrome length (default: 4; must be even and ≥ 4).</param>
    /// <param name="maxLength">Maximum palindrome length (default: 12; must be ≥ <paramref name="minLength"/>).</param>
    /// <returns>Every palindromic window, ordered by position then length.</returns>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        DnaSequence sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidatePalindromeParameters(minLength, maxLength);

        return FindPalindromesCore(sequence.Sequence, minLength, maxLength);
    }

    /// <summary>
    /// Finds reverse-complement palindromes in a raw sequence string (case-insensitive).
    /// Same model and validation as <see cref="FindPalindromes(DnaSequence, int, int)"/>; <c>null</c> or empty
    /// input yields no results. Characters other than A/C/G/T never pair.
    /// </summary>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        string sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ValidatePalindromeParameters(minLength, maxLength);

        return string.IsNullOrEmpty(sequence)
            ? Array.Empty<PalindromeResult>()
            : FindPalindromesCore(sequence.ToUpperInvariant(), minLength, maxLength);
    }

    private static void ValidatePalindromeParameters(int minLength, int maxLength)
    {
        if (minLength < 4 || minLength % 2 != 0)
            throw new ArgumentOutOfRangeException(nameof(minLength), minLength, "Must be even and >= 4");
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);
    }

    private static IEnumerable<PalindromeResult> FindPalindromesCore(
        string seq,
        int minLength,
        int maxLength)
    {
        int n = seq.Length;
        int minHalf = minLength / 2;
        int maxHalf = Math.Min(maxLength, n) / 2; // capping by n also avoids int overflow for huge maxLength
        if (maxHalf < minHalf)
            yield break;

        // radius[c] = number of consecutive Watson–Crick pairs (seq[c−1−k], seq[c+k]) around the gap before c,
        // capped at maxHalf. The window (i, 2h) is a palindrome iff radius[i + h] ≥ h.
        var radius = new int[n];
        for (int c = 1; c < n; c++)
        {
            int r = 0;
            while (r < maxHalf && c - r - 1 >= 0 && c + r < n && IsWatsonCrickPair(seq[c - r - 1], seq[c + r]))
                r++;
            radius[c] = r;
        }

        for (int i = 0; i + 2 * minHalf <= n; i++)
        {
            for (int h = minHalf; h <= maxHalf && i + 2 * h <= n; h++)
            {
                if (radius[i + h] >= h)
                    yield return new PalindromeResult(Position: i, Sequence: seq.Substring(i, 2 * h), Length: 2 * h);
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
/// Tandem Repeats Finder Bernoulli-model statistical measures (Benson 1999) for a detected repeat tract.
/// Benson models the alignment of two adjacent copies as a sequence of independent Bernoulli trials;
/// <see cref="MatchProbability"/> = P(Heads) = PM (matching probability) is the average percent identity
/// between adjacent copies, and <see cref="IndelProbability"/> = PI (indel probability) is the average
/// percentage of insertions and deletions between them. The statistics are computed between ADJACENT
/// copies, not between the tract and the consensus pattern.
/// </summary>
public readonly record struct TandemRepeatBernoulliStatistics(
    int Period,
    int AdjacentCopyPairs,
    int BernoulliTrials,
    int Matches,
    int Mismatches,
    int Indels,
    double MatchProbability,
    double IndelProbability,
    double PercentMatches,
    double PercentIndels,
    double ExpectedMatches,
    bool MeetsExpectedMatchProbability);

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
/// Summary of perfect microsatellites (1–6 bp units) in a sequence; see
/// <see cref="RepeatFinder.GetTandemRepeatSummary(DnaSequence,int)"/> for each field's definition.
/// The six per-class counts (mono … hexa) sum to <see cref="TotalRepeats"/>.
/// </summary>
public readonly record struct TandemRepeatSummary(
    int TotalRepeats,
    int TotalRepeatBases,
    double PercentageOfSequence,
    int MononucleotideRepeats,
    int DinucleotideRepeats,
    int TrinucleotideRepeats,
    int TetranucleotideRepeats,
    int PentanucleotideRepeats,
    int HexanucleotideRepeats,
    MicrosatelliteResult? LongestRepeat,
    string? MostFrequentUnit);
