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
            sequence.Sequence, UniformThresholds(sequence.Length, minUnitLength, maxUnitLength, minRepeats), CancellationToken.None, null);
    }

    /// <summary>
    /// Finds microsatellites with cancellation support. Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length.</param>
    /// <param name="maxUnitLength">Maximum repeat unit length.</param>
    /// <param name="minRepeats">Minimum number of repeats.</param>
    /// <param name="cancellationToken">Cancellation token, checked every 1000 visited run starts and once at the end;
    /// cancellation surfaces as <see cref="OperationCanceledException"/> while the result is enumerated.</param>
    /// <param name="progress">Optional progress reporter: non-decreasing values in [0, 1) (fraction of the unit-length ×
    /// position scan space visited) every 1000 visited run starts, then exactly 1.0 when the scan completes. Reported
    /// synchronously during enumeration (the result is lazy).</param>
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
            sequence.Sequence, UniformThresholds(sequence.Length, minUnitLength, maxUnitLength, minRepeats), cancellationToken, progress);
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
            sequence.ToUpperInvariant(), UniformThresholds(sequence.Length, minUnitLength, maxUnitLength, minRepeats), CancellationToken.None, null);
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
            sequence.ToUpperInvariant(), UniformThresholds(sequence.Length, minUnitLength, maxUnitLength, minRepeats), cancellationToken, progress);
    }

    /// <summary>
    /// MISA default microsatellite definition (Thiel et al. 2003; <c>misa.ini</c>
    /// <c>definition(unit_size,min_repeats): 1-10 2-6 3-5 4-5 5-5 6-5</c>): minimum number of complete copies per unit
    /// length — mononucleotide ≥ 10, dinucleotide ≥ 6, tri- to hexanucleotide ≥ 5. Read-only.
    /// </summary>
    public static IReadOnlyDictionary<int, int> MisaDefaultMinRepeats { get; } =
        new System.Collections.ObjectModel.ReadOnlyDictionary<int, int>(new Dictionary<int, int>
        {
            [1] = 10, [2] = 6, [3] = 5, [4] = 5, [5] = 5, [6] = 5,
        });

    /// <summary>
    /// MISA default maximal number of bases interrupting two SSRs in a compound microsatellite
    /// (<c>misa.ini</c> <c>interruptions(max_difference_for_2_SSRs): 100</c>).
    /// </summary>
    public const int MisaDefaultMaxInterruption = 100;

    /// <summary>
    /// Finds microsatellites with a separate minimum number of copies per unit length (MISA-style definition,
    /// e.g. <see cref="MisaDefaultMinRepeats"/> = <c>1-10 2-6 3-5 4-5 5-5 6-5</c>). Only the unit lengths present as
    /// keys are searched. Detection semantics are those of
    /// <see cref="FindMicrosatellites(DnaSequence,int,int,int)"/> (one maximal primitive ACGT run per unit length,
    /// reported at its left end with complete copies); a run is reported when its copy number is at least the
    /// threshold of its unit length. Results are ordered by unit length, then position.
    /// </summary>
    /// <remarks>
    /// Compared with a <c>misa.pl</c> run using the same definition, results are identical except where MISA's
    /// left-to-right regex scan (<c>([acgt]{p})\2{k-1,}</c>, resumed after each match) differs from the maximal-run
    /// convention: (a) a run of the same unit length that overlaps the previous match by fewer than p bases is
    /// truncated (MISA) instead of reported from its true left end; (b) a match whose unit is not primitive
    /// (e.g. <c>ATAT</c>) is rejected by MISA after it has consumed the bases, hiding a primitive run of that unit
    /// length that starts inside it.
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minRepeatsByUnitLength">Unit length (≥ 1) → minimum number of complete copies (≥ 2); at least one entry.</param>
    /// <param name="cancellationToken">Cancellation token, checked every 1000 visited run starts.</param>
    /// <param name="progress">Optional progress reporter (non-decreasing values in [0, 1], final 1.0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="minRepeatsByUnitLength"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="minRepeatsByUnitLength"/> is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A unit length is &lt; 1 or a minimum copy number is &lt; 2.</exception>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength);

        return FindMicrosatellitesCore(sequence.Sequence, thresholds, cancellationToken, progress);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string (case-insensitive; non-ACGT symbols never form a unit) with
    /// a minimum copy number per unit length. Same semantics as
    /// <see cref="FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},CancellationToken,IProgress{double})"/>;
    /// <c>null</c>/empty input yields no results.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(sequence.ToUpperInvariant(), thresholds, cancellationToken, progress);
    }

    /// <summary>
    /// Finds microsatellites with a minimum number of copies per unit length and an explicit scan convention.
    /// <see cref="MicrosatelliteScanMode.MaximalRuns"/> is
    /// <see cref="FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},CancellationToken,IProgress{double})"/>;
    /// <see cref="MicrosatelliteScanMode.MisaRegex"/> reproduces the SSR list of <c>misa.pl</c> v1.0 exactly
    /// (Thiel et al. 2003): for each unit length p in ascending order the leftmost match of
    /// <c>([acgt]{p})\2{t−1,}</c> (case-insensitive, all complete copies) is taken, the scan resumes at its end,
    /// and a non-primitive motif (e.g. <c>ATAT</c>, <c>AA</c>) is rejected only after its bases were consumed. The two
    /// conventions differ exactly in the two cases listed in the remarks of the map-based overload. Results are
    /// ordered by unit length, then position — misa.pl's SSR numbering.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minRepeatsByUnitLength">Unit length (≥ 1) → minimum number of complete copies (≥ 2); at least one entry.</param>
    /// <param name="scanMode">Scan convention.</param>
    /// <param name="cancellationToken">Cancellation token, checked every 1000 visited positions.</param>
    /// <param name="progress">Optional progress reporter (non-decreasing values in [0, 1), final 1.0).</param>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        MicrosatelliteScanMode scanMode,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength);
        ValidateScanMode(scanMode);

        return FindMicrosatellitesCore(sequence.Sequence, thresholds, cancellationToken, progress, scanMode);
    }

    /// <summary>
    /// Raw-string counterpart (case-insensitive; <c>null</c>/empty input yields no results) of
    /// <see cref="FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},MicrosatelliteScanMode,CancellationToken,IProgress{double})"/>.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        MicrosatelliteScanMode scanMode,
        CancellationToken cancellationToken = default,
        IProgress<double>? progress = null)
    {
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength);
        ValidateScanMode(scanMode);

        if (string.IsNullOrEmpty(sequence))
            return [];

        return FindMicrosatellitesCore(sequence.ToUpperInvariant(), thresholds, cancellationToken, progress, scanMode);
    }

    private static void ValidateScanMode(MicrosatelliteScanMode scanMode)
    {
        if (scanMode is not (MicrosatelliteScanMode.MaximalRuns or MicrosatelliteScanMode.MisaRegex))
            throw new ArgumentOutOfRangeException(nameof(scanMode), scanMode, "Unknown scan mode.");
    }

    /// <summary>
    /// Validates a per-unit-length threshold map and returns it as (unit length, minRepeats) pairs in ascending
    /// unit length.
    /// </summary>
    private static (int UnitLength, int MinRepeats)[] ThresholdsFromMap(
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        int maxUnitLength = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(minRepeatsByUnitLength);
        if (minRepeatsByUnitLength.Count == 0)
            throw new ArgumentException("At least one unit length must be given.", nameof(minRepeatsByUnitLength));

        var thresholds = new (int UnitLength, int MinRepeats)[minRepeatsByUnitLength.Count];
        int k = 0;
        foreach (var (unitLength, minRepeats) in minRepeatsByUnitLength)
        {
            if (unitLength < 1 || unitLength > maxUnitLength)
                throw new ArgumentOutOfRangeException(nameof(minRepeatsByUnitLength), unitLength,
                    $"Unit lengths must be in [1, {maxUnitLength}].");
            if (minRepeats < 2)
                throw new ArgumentOutOfRangeException(nameof(minRepeatsByUnitLength), minRepeats,
                    $"The minimum number of copies for unit length {unitLength} must be at least 2.");
            thresholds[k++] = (unitLength, minRepeats);
        }

        Array.Sort(thresholds, (a, b) => a.UnitLength.CompareTo(b.UnitLength));
        return thresholds;
    }

    private static void ValidateMicrosatelliteParameters(int minUnitLength, int maxUnitLength, int minRepeats)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minUnitLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxUnitLength, minUnitLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(minRepeats, 2);
    }

    /// <summary>
    /// Single scan core (upper-case input): for each (unit length, minimum copies) threshold in ascending unit
    /// length, visits only run starts (left-maximal positions), extends the run to its right-maximal end and
    /// reports it once. O(n) character comparisons per unit length.
    /// </summary>
    /// <remarks>
    /// Progress (when a reporter is given) is the fraction of the (unit length × position) scan space already
    /// visited: reported every 1000 visited run starts as <c>(k·n + i) / (K·n)</c> for the k-th of K unit lengths,
    /// so the values are non-decreasing and lie in [0, 1); a final <c>1.0</c> is reported once the scan is
    /// complete. The token is checked at the same points and once more before the final report.
    /// </remarks>
    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCore(
        string seq,
        (int UnitLength, int MinRepeats)[] thresholds,
        CancellationToken cancellationToken,
        IProgress<double>? progress,
        MicrosatelliteScanMode scanMode = MicrosatelliteScanMode.MaximalRuns)
    {
        int n = seq.Length;
        double totalPositions = Math.Max(1.0, (double)n * thresholds.Length);
        int sinceCheck = 0;
        const int checkInterval = 1000;
        int[]? nonAcgtBefore = null; // MISA scan: prefix count of non-ACGT symbols

        for (int k = 0; k < thresholds.Length; k++)
        {
            var (unitLen, minRepeats) = thresholds[k];
            if ((long)unitLen * minRepeats > n)
                continue; // minRepeats copies of this unit cannot fit

            if (scanMode == MicrosatelliteScanMode.MisaRegex)
            {
                if (nonAcgtBefore is null)
                {
                    nonAcgtBefore = new int[n + 1];
                    for (int x = 0; x < n; x++)
                        nonAcgtBefore[x + 1] = nonAcgtBefore[x] + (AcgtCode(seq[x]) < 0 ? 1 : 0);
                }

                // misa.pl: while ($seq =~ /(([acgt]{p})\2{t-1,})/ig) — leftmost match at or after the end of the
                // previous one, greedy (all complete copies); a non-primitive motif is rejected only after the match
                // has consumed its bases ("next if $redundant" keeps pos()).
                int start = 0;
                int runEnd = -1; // S[x] == S[x − p] for every x in [start + p, runEnd)
                while (start + unitLen * minRepeats <= n)
                {
                    if (++sinceCheck >= checkInterval)
                    {
                        sinceCheck = 0;
                        cancellationToken.ThrowIfCancellationRequested();
                        progress?.Report(((double)k * n + start) / totalPositions);
                    }

                    if (runEnd < start + unitLen)
                    {
                        runEnd = start + unitLen;
                        while (runEnd < n && seq[runEnd] == seq[runEnd - unitLen])
                            runEnd++;
                    }

                    int copies = (runEnd - start) / unitLen;
                    if (copies < minRepeats || nonAcgtBefore[start + unitLen] != nonAcgtBefore[start])
                    {
                        start++; // no match here ([acgt]{p} fails, or fewer than t copies)
                        continue;
                    }

                    string motif = seq.Substring(start, unitLen);
                    if (!IsRedundantUnit(motif))
                    {
                        yield return new MicrosatelliteResult(
                            Position: start,
                            RepeatUnit: motif,
                            RepeatCount: copies,
                            TotalLength: copies * unitLen,
                            RepeatType: ClassifyRepeatType(motif));
                    }

                    start += copies * unitLen; // resume the scan at pos()
                }

                continue;
            }

            int i = 0;
            while (i + unitLen * minRepeats <= n)
            {
                if (++sinceCheck >= checkInterval)
                {
                    sinceCheck = 0;
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report(((double)k * n + i) / totalPositions);
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
    /// One (unit length, minRepeats) pair per unit length in [min, max], capped at the longest unit whose
    /// <paramref name="minRepeats"/> copies still fit a sequence of length <paramref name="n"/> (longer units can
    /// never be reported, so a huge <paramref name="maxUnitLength"/> costs nothing).
    /// </summary>
    private static (int UnitLength, int MinRepeats)[] UniformThresholds(
        int n, int minUnitLength, int maxUnitLength, int minRepeats)
    {
        long lastFitting = Math.Min(maxUnitLength, n / minRepeats);
        long count = Math.Max(0, lastFitting - minUnitLength + 1);
        var thresholds = new (int, int)[count];
        for (int k = 0; k < count; k++)
            thresholds[k] = (minUnitLength + k, minRepeats);
        return thresholds;
    }

    /// <summary>
    /// A unit is reportable when it consists only of A/C/G/T and is primitive (not a power of a shorter word).
    /// </summary>
    private static bool IsReportableUnit(string seq, int start, int unitLen)
    {
        for (int k = start; k < start + unitLen; k++)
        {
            if (AcgtCode(seq[k]) < 0)
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

    #region Compound Microsatellites (MISA)

    /// <summary>
    /// Finds compound microsatellites: two or more microsatellites separated by at most
    /// <paramref name="maxInterruption"/> bases (MISA, Thiel et al. 2003, <c>misa.pl</c> types <c>c</c> and <c>c*</c>).
    /// The component SSRs are those of
    /// <see cref="FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},CancellationToken,IProgress{double})"/>
    /// with <paramref name="minRepeatsByUnitLength"/> (default <see cref="MisaDefaultMinRepeats"/>); they are chained
    /// exactly as by <see cref="AssembleCompoundMicrosatellites"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minRepeatsByUnitLength">Unit length → minimum copies; <c>null</c> = MISA default <c>1-10 2-6 3-5 4-5 5-5 6-5</c>.</param>
    /// <param name="maxInterruption">Maximal number of bases between two adjacent SSRs of a compound (≥ 0; MISA default 100).</param>
    /// <returns>Compound microsatellites ordered by start position (only records with ≥ 2 components).</returns>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int>? minRepeatsByUnitLength = null,
        int maxInterruption = MisaDefaultMaxInterruption)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfNegative(maxInterruption);
        var ssrs = FindMicrosatellites(sequence, minRepeatsByUnitLength ?? MisaDefaultMinRepeats);
        return AssembleCompoundsCore(sequence.Sequence, ssrs, maxInterruption);
    }

    /// <summary>
    /// Finds compound microsatellites in a raw sequence string (case-insensitive; interruption strings keep the
    /// input symbols, lower-cased as in MISA). See
    /// <see cref="FindCompoundMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},int)"/>.
    /// </summary>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        string sequence,
        IReadOnlyDictionary<int, int>? minRepeatsByUnitLength = null,
        int maxInterruption = MisaDefaultMaxInterruption)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxInterruption);
        var ssrs = FindMicrosatellites(sequence, minRepeatsByUnitLength ?? MisaDefaultMinRepeats);
        return string.IsNullOrEmpty(sequence) ? [] : AssembleCompoundsCore(sequence, ssrs, maxInterruption);
    }

    /// <summary>
    /// Finds compound microsatellites with an explicit scan convention for the component SSRs. With
    /// <see cref="MicrosatelliteScanMode.MisaRegex"/> the result equals the <c>c</c>/<c>c*</c> rows of <c>misa.pl</c>
    /// v1.0 run with the same <c>misa.ini</c>, except that SSRs sharing a start position are chained in unit-length
    /// order here, while misa.pl orders such ties by Perl's per-process randomised hash order
    /// (<c>sort { $start{$a} &lt;=&gt; $start{$b} } keys %start</c>; different <c>PERL_HASH_SEED</c> values give different
    /// misa.pl outputs). Ties need two primitive runs of different unit lengths starting at one position, which is
    /// impossible with the MISA default thresholds (Fine–Wilf).
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minRepeatsByUnitLength">Unit length → minimum copies; <c>null</c> = MISA default <c>1-10 2-6 3-5 4-5 5-5 6-5</c>.</param>
    /// <param name="maxInterruption">Maximal number of bases between two adjacent SSRs of a compound (≥ 0; MISA default 100).</param>
    /// <param name="scanMode">Scan convention for the component SSRs.</param>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int>? minRepeatsByUnitLength,
        int maxInterruption,
        MicrosatelliteScanMode scanMode)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfNegative(maxInterruption);
        var ssrs = FindMicrosatellites(sequence, minRepeatsByUnitLength ?? MisaDefaultMinRepeats, scanMode);
        return AssembleCompoundsCore(sequence.Sequence, ssrs, maxInterruption);
    }

    /// <summary>
    /// Raw-string counterpart of
    /// <see cref="FindCompoundMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},int,MicrosatelliteScanMode)"/>.
    /// </summary>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        string sequence,
        IReadOnlyDictionary<int, int>? minRepeatsByUnitLength,
        int maxInterruption,
        MicrosatelliteScanMode scanMode)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxInterruption);
        var ssrs = FindMicrosatellites(sequence, minRepeatsByUnitLength ?? MisaDefaultMinRepeats, scanMode);
        return string.IsNullOrEmpty(sequence) ? [] : AssembleCompoundsCore(sequence, ssrs, maxInterruption);
    }

    /// <summary>
    /// Finds compound microsatellites whose components are microsatellites of unit length 1–6 with at least
    /// <paramref name="minRepeats"/> complete copies (one threshold for every unit length).
    /// </summary>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        DnaSequence sequence,
        int minRepeats,
        int maxInterruption = MisaDefaultMaxInterruption)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindCompoundMicrosatellites(sequence, UniformMap(minRepeats), maxInterruption);
    }

    /// <summary>
    /// Raw-string counterpart of <see cref="FindCompoundMicrosatellites(DnaSequence,int,int)"/>.
    /// </summary>
    public static IReadOnlyList<CompoundMicrosatelliteResult> FindCompoundMicrosatellites(
        string sequence,
        int minRepeats,
        int maxInterruption = MisaDefaultMaxInterruption) =>
        FindCompoundMicrosatellites(sequence, UniformMap(minRepeats), maxInterruption);

    /// <summary>
    /// Chains a given list of microsatellites into compound microsatellites with MISA's rule (<c>misa.pl</c> v1.0):
    /// the SSRs are ordered by start; consecutive SSRs i, i+1 belong to the same compound when the number of bases
    /// between them, <c>start(i+1) − end(i)</c> (0-based start, exclusive end), is ≤ <paramref name="maxInterruption"/>
    /// — adjacent SSRs (0 bases) and overlapping SSRs (negative) always join. The comparison uses the end of the
    /// previous SSR in the chain, not the maximal end so far, and the compound ends where its last SSR ends
    /// (MISA's <c>end</c> column). A compound containing an overlapping pair is type <c>c*</c>, otherwise <c>c</c>.
    /// </summary>
    /// <remarks>
    /// Use this overload to chain SSRs detected under another convention (e.g. MISA's own SSR list): given the
    /// same SSR list, the result equals MISA's <c>c</c>/<c>c*</c> rows (start, end, type and the notation
    /// <c>(AT)6ccgt(GA)7</c> / <c>(A)10(AT)6*</c>). The ordering by start is stable: SSRs with equal start keep their
    /// input order (MISA orders such ties by Perl hash order, so feeding MISA's SSRs in its own order reproduces its
    /// output). SSRs that are not part of a compound (MISA types p1…p6) are not returned.
    /// </remarks>
    /// <param name="sequence">The sequence the SSRs were found in (used for the interruption strings).</param>
    /// <param name="microsatellites">SSRs (positions within <paramref name="sequence"/>).</param>
    /// <param name="maxInterruption">Maximal number of interrupting bases (≥ 0; MISA default 100).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxInterruption"/> is negative, or an SSR lies
    /// outside the sequence or has an empty unit / non-positive copy number.</exception>
    public static IReadOnlyList<CompoundMicrosatelliteResult> AssembleCompoundMicrosatellites(
        string sequence,
        IEnumerable<MicrosatelliteResult> microsatellites,
        int maxInterruption = MisaDefaultMaxInterruption)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(microsatellites);
        ArgumentOutOfRangeException.ThrowIfNegative(maxInterruption);

        var list = microsatellites.ToList();
        foreach (var m in list)
        {
            if (string.IsNullOrEmpty(m.RepeatUnit) || m.RepeatCount < 1 || m.Position < 0
                || (long)m.Position + m.TotalLength > sequence.Length || m.TotalLength < 1)
            {
                throw new ArgumentOutOfRangeException(nameof(microsatellites), m,
                    "Every microsatellite must have a non-empty unit, at least one copy and lie inside the sequence.");
            }
        }

        return AssembleCompoundsCore(sequence, list, maxInterruption);
    }

    private static Dictionary<int, int> UniformMap(int minRepeats)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minRepeats, 2);
        var map = new Dictionary<int, int>();
        for (int p = 1; p <= 6; p++)
            map[p] = minRepeats;
        return map;
    }

    private static List<CompoundMicrosatelliteResult> AssembleCompoundsCore(
        string sequence, IEnumerable<MicrosatelliteResult> microsatellites, int maxInterruption)
    {
        // misa.pl: @order = sort { $start{$a} <=> $start{$b} } keys %start. The sort here is stable, so SSRs with
        // equal start keep their input order (FindMicrosatellites lists them by unit length).
        var ordered = microsatellites
            .OrderBy(m => m.Position)
            .ToList();

        var compounds = new List<CompoundMicrosatelliteResult>();
        int i = 0;
        while (i < ordered.Count)
        {
            // Single SSR (MISA type p1…p6): last one, or the next starts more than maxInterruption bases later.
            if (i + 1 >= ordered.Count || Gap(ordered[i], ordered[i + 1]) > maxInterruption)
            {
                i++;
                continue;
            }

            var components = new List<MicrosatelliteResult> { ordered[i] };
            var interruptions = new List<string>();
            var notation = new System.Text.StringBuilder();
            AppendComponent(notation, ordered[i]);
            bool overlapping = false;

            int j = i;
            while (j + 1 < ordered.Count && Gap(ordered[j], ordered[j + 1]) <= maxInterruption)
            {
                int gap = Gap(ordered[j], ordered[j + 1]);
                if (gap < 0)
                {
                    // misa.pl: "($motif)$repeats*" — overlapping SSRs, compound type c*.
                    overlapping = true;
                    interruptions.Add(string.Empty);
                    AppendComponent(notation, ordered[j + 1]);
                    notation.Append('*');
                }
                else
                {
                    // misa.pl: $interssr = lc substr($seq, end, start − end − 1) (1-based) = the gap bases.
                    string interruption = sequence.Substring(End(ordered[j]), gap).ToLowerInvariant();
                    interruptions.Add(interruption);
                    notation.Append(interruption);
                    AppendComponent(notation, ordered[j + 1]);
                }

                components.Add(ordered[j + 1]);
                j++;
            }

            compounds.Add(new CompoundMicrosatelliteResult(
                Start: ordered[i].Position,
                End: End(ordered[j]),
                Components: components,
                Interruptions: interruptions,
                IsOverlapping: overlapping,
                Notation: notation.ToString()));
            i = j + 1;
        }

        return compounds;

        static int End(MicrosatelliteResult m) => m.Position + m.TotalLength;
        static int Gap(MicrosatelliteResult a, MicrosatelliteResult b) => b.Position - End(a);
        static void AppendComponent(System.Text.StringBuilder sb, MicrosatelliteResult m) =>
            sb.Append('(').Append(m.RepeatUnit.ToUpperInvariant()).Append(')').Append(m.RepeatCount);
    }

    #endregion

    #region Canonical Motifs (MISA classes / Krait standard motifs)

    /// <summary>
    /// MISA repeat-type class of a motif "considering sequence complementary" (<c>misa.pl</c> <c>.statistics</c>,
    /// table "Frequency of classified repeat types"): <c>X/Y</c> where X and Y are the lexicographically smallest
    /// rotations of the motif and of its reverse complement, the smaller one first — e.g. AC, CA, GT, TG →
    /// <c>AC/GT</c>; A, T → <c>A/T</c>; AT → <c>AT/AT</c>.
    /// </summary>
    /// <param name="motif">Repeat unit of A/C/G/T (case-insensitive).</param>
    /// <exception cref="ArgumentException"><paramref name="motif"/> is null/empty or contains a non-ACGT symbol.</exception>
    public static string GetCanonicalMotifClass(string motif)
    {
        string m = NormalizeMotif(motif);
        string forward = MinimalRotation(m, string.CompareOrdinal);
        string reverse = MinimalRotation(DnaSequence.GetReverseComplementString(m), string.CompareOrdinal);
        return string.CompareOrdinal(forward, reverse) < 0 ? $"{forward}/{reverse}" : $"{reverse}/{forward}";
    }

    /// <summary>
    /// Krait standard motif (Du et al. 2018, lmdu/krait <c>motif.py</c> <c>StandardMotif.standard</c>): the smallest
    /// member of the motif's equivalence set under Krait's base order A &lt; T &lt; C &lt; G. Level 0 = the motif
    /// itself; 1 = its rotations ("similar motifs"); 2 = + rotations of the reverse complement; 3 = + rotations of
    /// the complement; 4 = + rotations of the reverse. Krait's GUI default is level 3; level 2 is the
    /// rotation + reverse-complement standardization (the MISA class grouping, e.g. AC/CA/GT/TG → AC; ACAT → ATAC).
    /// </summary>
    /// <param name="motif">Repeat unit of A/C/G/T (case-insensitive).</param>
    /// <param name="level">Standardization level 0–4 (default 2).</param>
    public static string GetStandardMotif(string motif, int level = 2)
    {
        string m = NormalizeMotif(motif);
        ArgumentOutOfRangeException.ThrowIfLessThan(level, 0);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(level, 4);
        if (level == 0)
            return m;

        string best = MinimalRotation(m, CompareKraitOrder);
        if (level >= 2)
            best = MinOf(best, MinimalRotation(DnaSequence.GetReverseComplementString(m), CompareKraitOrder));
        if (level >= 3)
        {
            string complement = string.Concat(m.Select(SequenceExtensions.GetComplementBase));
            best = MinOf(best, MinimalRotation(complement, CompareKraitOrder));
        }
        if (level >= 4)
        {
            var reversed = m.ToCharArray();
            Array.Reverse(reversed);
            best = MinOf(best, MinimalRotation(new string(reversed), CompareKraitOrder));
        }

        return best;

        static string MinOf(string current, string candidate) => CompareKraitOrder(candidate, current) < 0 ? candidate : current;
    }

    /// <summary>
    /// Counts microsatellites per MISA repeat-type class (<see cref="GetCanonicalMotifClass"/>): the MISA
    /// <c>.statistics</c> table "Frequency of classified repeat types (considering sequence complementary)", column
    /// <c>total</c>. Keys are ordered like MISA's rows (class length, then ordinal).
    /// </summary>
    public static IReadOnlyDictionary<string, int> GetCanonicalMotifFrequencies(
        IEnumerable<MicrosatelliteResult> microsatellites)
    {
        ArgumentNullException.ThrowIfNull(microsatellites);
        return CountByKey(microsatellites, GetCanonicalMotifClass);
    }

    /// <summary>
    /// Counts microsatellites per Krait standard motif (<see cref="GetStandardMotif"/> at <paramref name="level"/>),
    /// i.e. Krait's per-standard-motif SSR counts. Keys are ordered by length, then ordinal.
    /// </summary>
    public static IReadOnlyDictionary<string, int> GetStandardMotifFrequencies(
        IEnumerable<MicrosatelliteResult> microsatellites, int level = 2)
    {
        ArgumentNullException.ThrowIfNull(microsatellites);
        ArgumentOutOfRangeException.ThrowIfLessThan(level, 0);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(level, 4);
        return CountByKey(microsatellites, u => GetStandardMotif(u, level));
    }

    private static IReadOnlyDictionary<string, int> CountByKey(
        IEnumerable<MicrosatelliteResult> microsatellites, Func<string, string> key)
    {
        var counts = new SortedDictionary<string, int>(
            Comparer<string>.Create((a, b) => a.Length != b.Length ? a.Length.CompareTo(b.Length) : string.CompareOrdinal(a, b)));
        var cache = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var m in microsatellites)
        {
            if (!cache.TryGetValue(m.RepeatUnit ?? string.Empty, out var k))
            {
                k = key(m.RepeatUnit!);
                cache[m.RepeatUnit!] = k;
            }

            counts[k] = counts.TryGetValue(k, out int c) ? c + 1 : 1;
        }

        return new System.Collections.ObjectModel.ReadOnlyDictionary<string, int>(
            counts.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    private static string NormalizeMotif(string motif)
    {
        if (string.IsNullOrEmpty(motif))
            throw new ArgumentException("Motif must be a non-empty A/C/G/T string.", nameof(motif));
        foreach (char c in motif)
        {
            if (AcgtCode(c) < 0)
                throw new ArgumentException($"Motif must contain only A/C/G/T; found '{c}'.", nameof(motif));
        }

        return motif.ToUpperInvariant();
    }

    private static string MinimalRotation(string s, Comparison<string> compare)
    {
        string best = s;
        for (int r = 1; r < s.Length; r++)
        {
            string rotation = string.Concat(s.AsSpan(r), s.AsSpan(0, r));
            if (compare(rotation, best) < 0)
                best = rotation;
        }

        return best;
    }

    /// <summary>Krait <c>motif_to_number</c> order for equal-length motifs: A &lt; T &lt; C &lt; G.</summary>
    private static int CompareKraitOrder(string a, string b)
    {
        for (int k = 0; k < Math.Min(a.Length, b.Length); k++)
        {
            int d = KraitRank(a[k]) - KraitRank(b[k]);
            if (d != 0)
                return d;
        }

        return a.Length.CompareTo(b.Length);

        static int KraitRank(char c) => c switch { 'A' => 1, 'T' => 2, 'C' => 3, 'G' => 4, _ => 5 };
    }

    #endregion

    #region Approximate (Imperfect/Interrupted) Tandem Repeat Detection — TRF model

    // --- Tandem Repeats Finder (TRF; Benson 1999) ----------------------------------------------------
    // Benson G (1999) "Tandem repeats finder: a program to analyze DNA sequences", Nucleic Acids Res
    // 27(2):573-580, https://doi.org/10.1093/nar/27.2.573; TRF 4.10.0 README (github.com/
    // Benson-Genomics-Lab/TRF: parameters, "TRF Definitions", "How does Tandem Repeats Finder work?").
    //
    // TRF has a DETECTION component (k-tuple matches at a common distance d, tested against statistical
    // criteria) and an ANALYSIS component (wraparound dynamic programming (WDP) of the sequence against
    // tandem copies of a candidate pattern, majority-rule consensus, realignment against the consensus,
    // statistics "between adjacent copies"). This region implements both components as TRF 4.10.0 runs them:
    //   * ANALYSIS: full WDP for patterns <= 20 (SMALLDISTANCE) and, for larger patterns, the narrow-band WDP
    //     ("we limit WDP calculations to a narrow diagonal band ... for patterns larger than 20 characters ... the
    //     band radius is Δd_max. The band is periodically recentered around a run of matches"), consensus,
    //     realignment and TRF's table statistics, entropy and alignment rows;
    //   * DETECTION: k-tuple matches with Benson's tuple sizes (Table 1; PM 80 and 75), the sum-of-heads criterion
    //     (exact mean/variance, reproduces TRF's tables for all d), the apparent-size criterion (the distribution
    //     TRF estimates by simulation, computed exactly here; README example 56 reproduced), the random-walk
    //     distance range d ± ⌊2.3·√(PI·d)⌋ over the distances that have reached a criteria test, per-distance
    //     "already aligned" suppression, the best-period list for d > 250, the multiples test (3 of the 5 best
    //     periods), the minimum copy-number rule, TRF's MAXDISTANCE, -l, -r, -f, -m.
    // Measured parity with compiled TRF on 700 random sequences with embedded repeats: 99.8–100 % of TRF rows
    // identical over seven parameter sets; every remaining row is caused by one entry of TRF's Monte-Carlo
    // apparent-size table that differs from the exact value (substituting TRF's table gives 100 %;
    // docs/Evidence/REP-APPROX-001-Evidence.md §WP7). TRF is AGPL-3.0, this library is MIT: the TRF source was
    // read to understand behaviour only; the code is written from the published method, and TRF's tables were
    // used solely as an oracle for the derived ones.

    /// <summary>TRF match weight: "Match ... The recomended values for Match Mismatch and Delta are 2, 7, and 7" (TRF README).</summary>
    private const int TrfMatchWeight = 2;

    /// <summary>TRF mismatch weight −7 (TRF README; weights are "interpreted as negative numbers").</summary>
    private const int TrfMismatchWeight = -7;

    /// <summary>TRF indel weight (Delta) −7 per gap column (TRF README).</summary>
    private const int TrfIndelWeight = -7;

    /// <summary>Marks a dead WDP cell (TRF: a zero cell beyond the candidate is set to −1000 so that the
    /// local alignment cannot restart there).</summary>
    private const int TrfDeadCell = -1000;

    /// <summary>
    /// TRF Min_Distance_Window = 20: the smallest tandem-repeat span the detector wants to see; the k-tuple
    /// distance window and the backward WDP scan both extend at least this far (TRF 4.10.0 source, tr30dat.h).
    /// </summary>
    private const int TrfMinDistanceWindow = 20;

    /// <summary>
    /// Default minimum alignment score to report a tandem repeat. Benson (1999): "Only those repeats
    /// scoring at least 50 ... are reported"; the recommended TRF parameter set uses Minscore = 50.
    /// </summary>
    public const int DefaultApproximateMinScore = 50;

    /// <summary>
    /// Largest supported period. TRF README (4.10.0): "TRF will throw an error if a value of over 2000 is
    /// given for MaxPeriod ... very large TRs are outside the scope of the TRF statistical models."
    /// </summary>
    public const int MaxApproximatePeriod = 2000;

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats with the Tandem Repeats Finder model
    /// (Benson 1999; recommended parameters Match 2, Mismatch 7, Delta 7, PM 80, PI 10).
    /// </summary>
    /// <remarks>
    /// <para><b>Analysis (TRF-exact).</b> For every candidate (position i, distance d) the sequence is aligned
    /// by wraparound dynamic programming (local alignment against unlimited tandem copies of the pattern
    /// S[i−d+1..i]); a consensus is taken by majority rule from that alignment and the sequence is realigned
    /// against the consensus. Reported values follow the TRF table: indices, period = most common distance
    /// between matching characters of adjacent copies, copy number = aligned consensus columns / consensus
    /// size, consensus size, % matches and % indels "between adjacent copies overall" (not between the
    /// sequence and the consensus), alignment score, nucleotide composition and entropy. A repeat needs at
    /// least 1.9 copies (1.8 for large patterns) and score ≥ <paramref name="minScore"/>. Scoring: +2 for
    /// an identical A/C/G/T pair, −7 for any other pair (N and other symbols never match), −7 per gap.</para>
    /// <para><b>Detection (TRF criteria).</b> A candidate is examined when a k-tuple match at distance d ends at i, the
    /// heads counted in k-runs over the last max(d, 20) positions reach the sum-of-heads cut-off and the matches are
    /// spread over that window (apparent-size criterion) — for d &gt; 20 possibly summed over the random-walk range
    /// d ± ⌊2.3·√(PI·d)⌋ — and, for d &gt; 250, d is among the best periods of any earlier analysed region spanning it;
    /// d must be among the three best periods of the aligned region (period 1 needs ≥ 80 % of one base); overlapping
    /// reports are reduced with TRF's redundancy rule (≥ 90 % overlap, same period or a multiple scoring ≤ 1.1×).
    /// Patterns longer than 20 are aligned in TRF's narrow diagonal band. Measured on random sequences with embedded
    /// repeats, 99.8 % of TRF rows are identical (100 % at region level) at maxPeriod 500; the remaining rows trace
    /// to TRF's simulated (noisy) apparent-size table (docs/Evidence/REP-APPROX-001-Evidence.md). These overloads
    /// examine distances up to <paramref name="maxPeriod"/> only (TRF examines up to MAXDISTANCE ≥ 200 and filters
    /// by MaxPeriod afterwards; use the <see cref="TandemRepeatsFinderParameters"/> overload for that).</para>
    /// <para>Complexity: O(n · maxPeriod) for the k-tuple scan plus one wraparound DP over the aligned region per
    /// examined candidate: O(region · pattern) for patterns ≤ 20 (O(pattern) memory unless the candidate passes the
    /// copy and best-period tests, which store the traceback matrix) and O(region · band) for larger patterns, band
    /// = 2·min(2·max(6, Δd_max), ⌊pattern/3⌋) + 1. The apparent-size cut-offs are computed once per PM (≈ 0.1 s).</para>
    /// <para>Coordinates are 0-based (<c>Start</c>; TRF prints 1-based indices). Percentages are exact
    /// (TRF truncates them to integers). Output is ordered by start, then end, then period.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minPeriod">Minimum reported period (≥ 1; default 1).</param>
    /// <param name="maxPeriod">Maximum candidate distance and reported period (minPeriod..2000; default 6).</param>
    /// <param name="minScore">Minimum TRF alignment score to report (≥ 1; default 50, Benson 1999).</param>
    /// <returns>Approximate tandem repeats ordered by start position.</returns>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        DnaSequence sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateApproximateParameters(minPeriod, maxPeriod, minScore);
        return FindApproximateTandemRepeatsCore(sequence.Sequence, TrfModel.Legacy(maxPeriod, minScore), minPeriod);
    }

    /// <summary>
    /// Finds approximate (imperfect / interrupted) tandem repeats in a raw sequence string with the Tandem
    /// Repeats Finder model; case-insensitive, any non-A/C/G/T symbol never matches. See
    /// <see cref="FindApproximateTandemRepeats(DnaSequence,int,int,int)"/>. Null or empty input yields no repeats.
    /// </summary>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        string sequence,
        int minPeriod = 1,
        int maxPeriod = 6,
        int minScore = DefaultApproximateMinScore)
    {
        ValidateApproximateParameters(minPeriod, maxPeriod, minScore);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ApproximateTandemRepeatResult>();

        return FindApproximateTandemRepeatsCore(sequence.ToUpperInvariant(), TrfModel.Legacy(maxPeriod, minScore), minPeriod);
    }

    /// <summary>
    /// Finds approximate tandem repeats with an explicit Tandem Repeats Finder parameter set
    /// (<c>trf File Match Mismatch Delta PM PI Minscore MaxPeriod [-l n] [-r] [-f]</c>).
    /// </summary>
    /// <remarks>
    /// Same model as <see cref="FindApproximateTandemRepeats(string,int,int,int)"/>, with every TRF parameter
    /// explicit (<see cref="TandemRepeatsFinderParameters"/>): alignment weights, PM (selects TRF's tuple sizes and
    /// sum-of-heads cut-offs; 80 or 75), PI (random-walk distance range for d &gt; 20), Minscore, MaxPeriod, maximum
    /// TR length (<c>-l</c>), redundancy elimination (<c>-r</c> disables it) and flanking sequence (<c>-f</c>).
    /// As in TRF, candidate distances are examined up to MAXDISTANCE = max(200, min(max(MaxPeriod, 500),
    /// ⌊0.6·n⌋)) and MaxPeriod only filters the reported periods (the legacy overloads examine distances up to
    /// maxPeriod only). PI sets the random-walk range d ± ⌊2.3·√(PI·d)⌋ used by the detection criteria for d &gt; 20 and
    /// the narrow-band radius for patterns &gt; 20 (Benson 1999). Every result carries the final alignment rows
    /// (<see cref="ApproximateTandemRepeatResult.AlignedSequence"/> / <see cref="ApproximateTandemRepeatResult.AlignedConsensus"/>),
    /// TRF's entropy (<see cref="ApproximateTandemRepeatResult.EntropyTrf"/>) and, when
    /// <see cref="TandemRepeatsFinderParameters.FlankLength"/> &gt; 0, the flanking sequences.
    /// </remarks>
    /// <param name="sequence">Sequence (case-insensitive; any non-A/C/G/T symbol never matches). Null or empty → no repeats.</param>
    /// <param name="parameters">TRF parameter set (<see cref="TandemRepeatsFinderParameters.Recommended"/> = 2 7 7 80 10 50 500).</param>
    /// <param name="minPeriod">Library option: minimum reported period (≥ 1), applied after redundancy elimination.</param>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        string sequence,
        TandemRepeatsFinderParameters parameters,
        int minPeriod = 1)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        parameters.Validate();
        ArgumentOutOfRangeException.ThrowIfLessThan(minPeriod, 1);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ApproximateTandemRepeatResult>();

        return FindApproximateTandemRepeatsCore(
            sequence.ToUpperInvariant(), TrfModel.FromParameters(parameters, sequence.Length), minPeriod);
    }

    /// <summary>
    /// Finds approximate tandem repeats in a <see cref="DnaSequence"/> with an explicit TRF parameter set; see
    /// <see cref="FindApproximateTandemRepeats(string,TandemRepeatsFinderParameters,int)"/>.
    /// </summary>
    public static IEnumerable<ApproximateTandemRepeatResult> FindApproximateTandemRepeats(
        DnaSequence sequence,
        TandemRepeatsFinderParameters parameters,
        int minPeriod = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindApproximateTandemRepeats(sequence.Sequence, parameters, minPeriod);
    }

    /// <summary>
    /// TRF masked sequence (<c>trf ... -m</c>): a copy of <paramref name="sequence"/> with every position inside a
    /// reported tandem repeat replaced by <c>N</c> (TRF README: "every location that occurred in a tandem repeat
    /// changed to the letter 'N'"), or — with <paramref name="softMask"/> — lower-cased (soft masking). Positions
    /// outside repeats are returned as given (TRF itself upper-cases its whole output).
    /// </summary>
    /// <param name="sequence">Sequence to mask (null → <see cref="ArgumentNullException"/>; empty → empty).</param>
    /// <param name="parameters">TRF parameters; null = <see cref="TandemRepeatsFinderParameters.Recommended"/>.</param>
    /// <param name="softMask">Lower-case repeat positions instead of writing <c>N</c>.</param>
    public static string MaskApproximateTandemRepeats(
        string sequence,
        TandemRepeatsFinderParameters? parameters = null,
        bool softMask = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var repeats = FindApproximateTandemRepeats(sequence, parameters ?? TandemRepeatsFinderParameters.Recommended);
        return MaskApproximateTandemRepeats(sequence, repeats, softMask);
    }

    /// <summary>
    /// Masks the given tandem repeats in <paramref name="sequence"/> (TRF <c>-m</c> semantics: positions
    /// <c>Start .. Start + SpanLength − 1</c> of every repeat → <c>N</c>, or lower case with
    /// <paramref name="softMask"/>). Repeats reaching beyond the sequence are rejected.
    /// </summary>
    public static string MaskApproximateTandemRepeats(
        string sequence,
        IEnumerable<ApproximateTandemRepeatResult> repeats,
        bool softMask = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(repeats);
        var masked = sequence.ToCharArray();
        foreach (var repeat in repeats)
        {
            if (repeat.Start < 0 || repeat.SpanLength < 0 || repeat.Start + repeat.SpanLength > masked.Length)
                throw new ArgumentOutOfRangeException(nameof(repeats), "A repeat lies outside the sequence.");
            for (int p = repeat.Start; p < repeat.Start + repeat.SpanLength; p++)
                masked[p] = softMask ? char.ToLowerInvariant(masked[p]) : 'N';
        }

        return new string(masked);
    }

    private static void ValidateApproximateParameters(int minPeriod, int maxPeriod, int minScore)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minPeriod, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxPeriod, minPeriod);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(maxPeriod, MaxApproximatePeriod);
        ArgumentOutOfRangeException.ThrowIfLessThan(minScore, 1);
    }

    /// <summary>
    /// Resolved TRF run settings: alignment weights (mismatch / indel stored negative, as TRF uses them),
    /// PM / PI in percent, reporting thresholds, detection distance limit, WDP row cap, and output options.
    /// </summary>
    private sealed class TrfModel
    {
        public int Match { get; init; } = TrfMatchWeight;
        public int Mismatch { get; init; } = TrfMismatchWeight;
        public int Indel { get; init; } = TrfIndelWeight;
        public int Pm { get; init; } = TrfRecommendedPm;
        public int Pi { get; init; } = 10;
        public int MinScore { get; init; } = DefaultApproximateMinScore;
        public int MaxPeriod { get; init; } = MaxApproximatePeriod;
        public int MaxDistance { get; init; } = MaxApproximatePeriod;
        public int MaxWrapLength { get; init; } = TandemRepeatsFinderParameters.DefaultMaxRepeatLength;
        public bool EliminateRedundancy { get; init; } = true;
        public int FlankLength { get; init; }

        /// <summary>Recommended weights (2, 7, 7), PM 80, PI 10 — used by the Bernoulli-statistics analysis.</summary>
        public static readonly TrfModel Recommended = new();

        /// <summary>The legacy overloads: recommended weights, distances examined up to maxPeriod only.</summary>
        public static TrfModel Legacy(int maxPeriod, int minScore) =>
            new() { MinScore = minScore, MaxPeriod = maxPeriod, MaxDistance = maxPeriod };

        /// <summary>A full TRF parameter set; MAXDISTANCE as TRF derives it from MaxPeriod and the length.</summary>
        public static TrfModel FromParameters(TandemRepeatsFinderParameters p, int length) => new()
        {
            Match = p.MatchWeight,
            Mismatch = -p.MismatchPenalty,
            Indel = -p.IndelPenalty,
            Pm = p.MatchProbability,
            Pi = p.IndelProbability,
            MinScore = p.MinScore,
            MaxPeriod = p.MaxPeriod,
            MaxDistance = Math.Max(TrfMinMaxDistance, Math.Min(Math.Max(p.MaxPeriod, TrfDefaultMaxDistance), (int)(length * 0.6))),
            MaxWrapLength = p.MaxRepeatLength,
            EliminateRedundancy = p.EliminateRedundancy,
            FlankLength = p.FlankLength,
        };

        /// <summary>Random-walk distance radius ⌊2.3·√(PI·d)⌋ (Benson 1999), applied for d &gt; 20.</summary>
        public int DistanceRadius(int d) => d <= TrfSmallDistance ? 0 : RandomWalkRange(d);

        /// <summary>Δd_max = ⌊2.3·√(PI·d)⌋ (Benson 1999 random-walk range; TRF: Pindel = (float)PI/100).</summary>
        public int RandomWalkRange(int d) => (int)Math.Floor(2.3 * Math.Sqrt((double)((float)Pi / 100) * d));
    }

    /// <summary>TRF MAXDISTANCE lower bound (200) and the value used for MaxPeriod &lt; 500 (500).</summary>
    private const int TrfMinMaxDistance = 200;
    private const int TrfDefaultMaxDistance = 500;

    /// <summary>Largest pattern size aligned with the full WDP and without a random-walk distance range (TRF: 20).</summary>
    private const int TrfSmallDistance = 20;

    /// <summary>One column of a TRF alignment: sequence symbol (or '-'), pattern symbol (or '-'), 1-based
    /// sequence index, 0-based pattern index (for an inserted sequence symbol: the next pattern position).</summary>
    private readonly record struct TrfColumn(char Seq, char Pat, int SeqIndex, int PatIndex);

    /// <summary>A traced WDP alignment; <see cref="Columns"/> run from the RIGHTMOST column to the leftmost.</summary>
    private sealed record TrfAlignment(int Score, double CopyNumber, TrfColumn[] Columns)
    {
        public int First => Columns[^1].SeqIndex;
        public int Last => Columns[0].SeqIndex;
    }

    /// <summary>Extent of a WDP alignment without traceback (same optimum and same path as the traceback).</summary>
    private readonly record struct TrfExtent(int Score, int First, int Last, double CopyNumber);

    private static int TrfWeight(TrfModel m, char a, char b) =>
        a == b && AcgtCode(a) >= 0 ? m.Match : m.Mismatch;

    /// <summary>
    /// TRF k-tuple size for distance d. PM = 80: Benson 1999 Table 1 (k = 4 for d ≤ 29, 5 for 30..159, 7 for
    /// ≥ 160). PM = 75: TRF 4.10.0's tuple set for that PM (k = 3 for d ≤ 29, 4 for 30..43, 5 for 44..159, 7 for
    /// ≥ 160 — printed by the program as "tuple sizes 0,3,4,5,7 / tuple distances 0, 29, 43, 159"); each range
    /// starts where the exact sum-of-heads cut-off already reaches k + 1, the paper's rule for picking tuple
    /// ranges (first such d at PM = .75: k = 4 → 27, k = 5 → 41, k = 7 → 91).
    /// </summary>
    private static int TrfTupleSize(int d, int pm = TrfRecommendedPm)
    {
        if (pm == 75)
        {
            if (d <= 29)
                return 3;
            if (d <= 43)
                return 4;
            return d <= 159 ? 5 : 7;
        }

        if (d <= 29)
            return 4;
        return d <= 159 ? 5 : 7;
    }

    /// <summary>TRF's recommended matching probability PM = 80 (percent).</summary>
    private const int TrfRecommendedPm = 80;

    /// <summary>Smallest sum-of-heads cut-off in either TRF table (PM = 80 and PM = 75): 5 heads.</summary>
    private const int TrfMinSumOfHeads = 5;

    private static readonly int[] SumOfHeadsCache = new int[MaxApproximatePeriod + 1];
    private static readonly int[] SumOfHeadsCache75 = new int[MaxApproximatePeriod + 1];

    /// <summary>Sum-of-heads criterion at the recommended PM = 80 (see <see cref="TrfSumOfHeadsCriterion(int,int)"/>).</summary>
    internal static int TrfSumOfHeadsCriterion(int d) => TrfSumOfHeadsCriterion(d, TrfRecommendedPm);

    /// <summary>
    /// Sum-of-heads criterion for distance d (Benson 1999): R(d,k,PM) = total heads in head runs of length ≥ k
    /// in an iid Bernoulli(PM) sequence of length d; "the distribution of R is well approximated by the normal
    /// distribution and its exact mean and variance can be calculated"; the criterion is the largest x such
    /// that R ≥ x 95 % of the time. Mean and variance are computed exactly by a run-length Markov chain; the
    /// cut-off is ⌊μ − 1.65σ⌋, never below k + 1 ("the smallest pattern for tuple size k [has] a sum-of-heads
    /// criterion of at least k+1") and never below 5 (the floor of both TRF tables; binding only for k = 3 at
    /// PM = 75, d ≤ 20). Reproduces TRF 4.10.0's PM = 80 and PM = 75 tables for every d = 1..2000.
    /// </summary>
    /// <param name="d">Distance (pattern size) 1..2000.</param>
    /// <param name="pmPercent">Matching probability PM in percent: 80 or 75 (the values TRF has data for).</param>
    internal static int TrfSumOfHeadsCriterion(int d, int pmPercent)
    {
        int[] cache = pmPercent == 75 ? SumOfHeadsCache75 : SumOfHeadsCache;
        int cached = Volatile.Read(ref cache[d]);
        if (cached != 0)
            return cached;

        // One run of the chain per tuple size gives the mean and variance for every d of that tuple size.
        for (int d0 = 1; d0 <= MaxApproximatePeriod;)
        {
            int k = TrfTupleSize(d0, pmPercent);
            int d1 = d0;
            while (d1 < MaxApproximatePeriod && TrfTupleSize(d1 + 1, pmPercent) == k)
                d1++;
            FillSumOfHeadsGroup(cache, d0, d1, k, pmPercent / 100.0);
            d0 = d1 + 1;
        }

        return Volatile.Read(ref cache[d]);
    }

    /// <summary>Exact mean and variance of R after each of d = 1..d1 tosses (run-length Markov chain), cut-offs for d0..d1.</summary>
    private static void FillSumOfHeadsGroup(int[] cache, int d0, int d1, int k, double pm)
    {
        // State r = current head-run length (0..k-1) or k = inside a run already counted.
        var p = new double[k + 1];
        var m1 = new double[k + 1];
        var m2 = new double[k + 1];
        var np = new double[k + 1];
        var n1 = new double[k + 1];
        var n2 = new double[k + 1];
        p[0] = 1.0;
        for (int step = 1; step <= d1; step++)
        {
            Array.Clear(np);
            Array.Clear(n1);
            Array.Clear(n2);
            for (int r = 0; r <= k; r++)
            {
                if (p[r] == 0.0)
                    continue;
                // tails: run resets, nothing added
                np[0] += p[r] * (1 - pm);
                n1[0] += m1[r] * (1 - pm);
                n2[0] += m2[r] * (1 - pm);
                // heads: the k-th head of a run adds k, every later head adds 1
                int next = r < k ? r + 1 : k;
                int add = 1;
                if (r < k - 1)
                    add = 0;
                else if (r == k - 1)
                    add = k;
                np[next] += p[r] * pm;
                n1[next] += (m1[r] + add * p[r]) * pm;
                n2[next] += (m2[r] + 2.0 * add * m1[r] + (double)add * add * p[r]) * pm;
            }
            (p, np) = (np, p);
            (m1, n1) = (n1, m1);
            (m2, n2) = (n2, m2);

            if (step < d0)
                continue;
            double mean = m1.Sum();
            double sd = Math.Sqrt(Math.Max(0.0, m2.Sum() - mean * mean));
            Volatile.Write(ref cache[step], Math.Max(Math.Max(k + 1, TrfMinSumOfHeads), (int)(mean - 1.65 * sd)));
        }
    }

    private static readonly int[]?[] ApparentSizeTables = new int[]?[2];

    /// <summary>
    /// Apparent-size criterion y(d, k, PM) (Benson 1999; TRF README "Apparent Size Distribution"): S = the distance
    /// between the first and the last run of ≥ k heads in an iid Bernoulli(PM) sequence of length d, conditional on
    /// the sequence meeting the sum-of-heads criterion; y is the largest number such that S &gt; y 95 % of the time
    /// (README example: PM = .75, k = 5, d = 100 → 56, reproduced). TRF estimates this distribution by simulation;
    /// here it is computed exactly (see <see cref="ComputeApparentSizeGroup"/>). Sequence length max(d, 20), the
    /// distance window TRF tests.
    /// </summary>
    internal static int TrfApparentSize(int d, int pmPercent) => GetApparentSizeTable(pmPercent)[d];

    /// <summary>
    /// The apparent-size test as applied to a distance window of length L = max(d, 20) ending at i: the first k-tuple
    /// match must end at most L − y − 1 positions after the window's left end (i − L), i.e. first and last tuple
    /// lie more than y apart. Matches TRF 4.10.0's simulated <c>waitdata</c> tables within their simulation noise
    /// (docs/Evidence/REP-APPROX-001-Evidence.md).
    /// </summary>
    internal static int TrfApparentSizeOffset(int d, int pmPercent) =>
        Math.Max(d, TrfMinDistanceWindow) - TrfApparentSize(d, pmPercent) - 1;

    private static int[] GetApparentSizeTable(int pmPercent)
    {
        int slot = pmPercent == 75 ? 1 : 0;
        var table = Volatile.Read(ref ApparentSizeTables[slot]);
        if (table is not null)
            return table;

        table = new int[MaxApproximatePeriod + 1];
        for (int d0 = 1; d0 <= MaxApproximatePeriod;)
        {
            int k = TrfTupleSize(d0, pmPercent);
            int d1 = d0;
            while (d1 < MaxApproximatePeriod && TrfTupleSize(d1 + 1, pmPercent) == k)
                d1++;
            ComputeApparentSizeGroup(table, d0, d1, k, pmPercent);
            d0 = d1 + 1;
        }

        Interlocked.CompareExchange(ref ApparentSizeTables[slot], table, null);
        return ApparentSizeTables[slot]!;
    }

    /// <summary>
    /// Exact conditional distribution of the apparent size S for every d in [d0, d1] (one tuple size k). With f the
    /// position of the k-th head of the first k-run and E the last head of the last k-run, a sequence of length L
    /// decomposes into: the prefix up to f (first k-run completes at f; probability F(f), contributing k heads), the
    /// middle f+1..E (a run-length Markov chain started inside a counted run, ending inside one at step S = E − f,
    /// accumulating the remaining heads R′ of the sum of heads), and the tail after E (a tail, then no k-run;
    /// probability T(L − E)). Hence P(S = s, R ≥ x) = M(s, R′ ≥ x − k) · Σ_f F(f)·T(L − s − f), where the middle
    /// chain M is independent of L. y = max{y : P(S &gt; y | R ≥ x) ≥ 0.95}.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveOptimization)]
    private static void ComputeApparentSizeGroup(int[] table, int d0, int d1, int k, int pmPercent)
    {
        double p = pmPercent / 100.0, q = 1.0 - p;
        int lMax = Math.Max(d1, TrfMinDistanceWindow);

        // noRun[m] = P(no run of k heads in m tosses).
        var noRun = new double[lMax + 2];
        var runState = new double[k];
        var nextState = new double[k];
        runState[0] = 1.0;
        noRun[0] = 1.0;
        for (int toss = 1; toss <= lMax + 1; toss++)
        {
            double total = 0;
            foreach (double v in runState)
                total += v;
            nextState[0] = total * q;
            for (int r = 1; r < k; r++)
                nextState[r] = runState[r - 1] * p;
            (runState, nextState) = (nextState, runState);
            total = 0;
            foreach (double v in runState)
                total += v;
            noRun[toss] = total;
        }

        double pk = Math.Pow(p, k);
        var firstRun = new double[lMax + 1];
        var tail = new double[lMax + 1];
        for (int f = k; f <= lMax; f++)
            firstRun[f] = f == k ? pk : pk * q * noRun[f - k - 1];
        tail[0] = 1.0;
        for (int len = 1; len <= lMax; len++)
            tail[len] = q * noRun[len - 1];
        var outer = new double[lMax + 1];
        for (int len = 0; len <= lMax; len++)
        {
            double acc = 0;
            for (int f = k; f <= len; f++)
                acc += firstRun[f] * tail[len - f];
            outer[len] = acc;
        }

        // Required middle heads per d, and the cap of the tracked head count (larger counts are pooled).
        int count = d1 - d0 + 1;
        var need = new int[count];
        var length = new int[count];
        int cap = 0;
        for (int t = 0; t < count; t++)
        {
            need[t] = Math.Max(0, TrfSumOfHeadsCriterion(d0 + t, pmPercent) - k);
            length[t] = Math.Max(d0 + t, TrfMinDistanceWindow);
            cap = Math.Max(cap, need[t]);
        }

        var totals = new double[count];
        var lowCum = new double[count];
        var result = new int[count];
        var done = new bool[count];
        for (int pass = 0; pass < 2; pass++)
        {
            // Middle chain: state (run length 0..k, k = inside a counted run) × heads counted so far (0..cap),
            // flattened as r·(cap+1) + h. Heads beyond the reachable maximum stay zero and are skipped.
            int stride = cap + 1;
            var state = new double[(k + 1) * stride];
            var next = new double[(k + 1) * stride];
            var atLeast = new double[cap + 2];
            state[k * stride] = 1.0;
            for (int step = 0; step <= lMax; step++)
            {
                if (step > 0)
                {
                    int reach = Math.Min(cap, k * step);
                    Array.Clear(next);
                    for (int r = 0; r <= k; r++)
                    {
                        int from = r * stride;
                        // A head extends the run; the k-th head of a run adds k heads, every later head adds 1.
                        int to = k * stride, add = 1;
                        if (r < k - 1)
                            (to, add) = ((r + 1) * stride, 0);
                        else if (r == k - 1)
                            add = k;
                        for (int h = 0; h <= reach; h++)
                        {
                            double v = state[from + h];
                            if (v == 0.0)
                                continue;
                            next[h] += v * q;
                            next[to + Math.Min(cap, h + add)] += v * p;
                        }
                    }

                    (state, next) = (next, state);
                }

                atLeast[cap + 1] = 0.0;
                for (int h = cap; h >= 0; h--)
                    atLeast[h] = atLeast[h + 1] + state[k * stride + h];

                for (int t = 0; t < count; t++)
                {
                    if (step > length[t] || done[t])
                        continue;
                    double mass = atLeast[need[t]] * outer[length[t] - step];
                    if (pass == 0)
                    {
                        totals[t] += mass;
                        continue;
                    }

                    // y = largest y in [0, L−1] with P(S ≥ y+1) ≥ 0.95·total, i.e. P(S ≤ y) ≤ 0.05·total.
                    lowCum[t] += mass;
                    if (totals[t] - lowCum[t] < 0.95 * totals[t])
                    {
                        result[t] = step - 1;
                        done[t] = true;
                    }
                }
            }
        }

        for (int t = 0; t < count; t++)
            table[d0 + t] = done[t] ? Math.Min(result[t], length[t] - 1) : length[t] - 1;
    }

    private static IReadOnlyList<ApproximateTandemRepeatResult> FindApproximateTandemRepeatsCore(
        string sequence,
        TrfModel m,
        int minPeriod)
    {
        int n = sequence.Length;
        var s = new char[n + 1]; // 1-based, as in TRF
        sequence.CopyTo(0, s, 1, n);

        int maxDistance = m.MaxDistance;
        var found = new List<ApproximateTandemRepeatResult>();
        var runLength = new int[maxDistance + 1];
        var seenEnd = new int[maxDistance + 1];
        var windows = new TupleMatchWindow?[maxDistance + 1];
        var linked = new bool[maxDistance + 1];
        var bestPeriods = new Dictionary<(int First, int Last), int[]>();
        var bestPeriodList = new List<TrfBestPeriodEntry>();

        for (int i = 2; i <= n; i++)
        {
            bool acgt = AcgtCode(s[i]) >= 0;
            int dMax = Math.Min(maxDistance, i - 1);
            for (int d = 1; d <= dMax; d++)
            {
                if (!acgt || s[i] != s[i - d])
                {
                    runLength[d] = 0;
                    continue;
                }

                int k = TrfTupleSize(d, m.Pm);
                if (++runLength[d] < k)
                    continue;

                // A k-tuple match at distance d ends at i; record it in the distance window of the last
                // max(d, 20) positions (adjacent tuple matches form one entry of growing size).
                var window = windows[d] ??= new TupleMatchWindow(Math.Max(d, TrfMinDistanceWindow) + 1);
                window.Add(i, k, TrfWindowLeft(i, d));

                if (seenEnd[d] >= i)
                    continue;

                // A distance takes part in the range sums of its neighbours only once it has itself reached the
                // criteria test (TRF links a distance list into the active range when it is tested, and unlinks it
                // when a range scan finds it empty).
                linked[d] = true;
                if (!MeetsTrfCriteria(windows, linked, i, d, k, m))
                    continue;
                if (d > TrfBestPeriodListMinDistance && !IsAllowedByTrfBestPeriodList(bestPeriodList, i, d, m))
                    continue;

                var repeat = AnalyzeTrfCandidate(m, s, n, i, d, seenEnd, bestPeriods, bestPeriodList);
                if (repeat is not null)
                    found.Add(repeat.Value);
            }
        }

        // TRF: drop periods above MaxPeriod, sort by start (stable), eliminate redundancy (unless -r). The
        // minimum period is a library option applied afterwards, so it never resurrects a redundant multiple.
        var sorted = found.Where(r => r.Period <= m.MaxPeriod).OrderBy(r => r.Start).ToList();
        var reported = m.EliminateRedundancy ? RemoveTrfRedundancy(sorted) : sorted;

        return reported
            .Where(r => r.Period >= minPeriod)
            .OrderBy(r => r.Start)
            .ThenBy(r => r.Start + r.SpanLength)
            .ThenBy(r => r.Period)
            .ToList();
    }

    /// <summary>Left end of the tuple-match window of distance d at position i: the last max(d, 20) positions.</summary>
    private static int TrfWindowLeft(int i, int d) => i - Math.Max(d, TrfMinDistanceWindow) + 1;

    /// <summary>
    /// TRF candidate criteria for a k-tuple match at distance d ending at i (Benson 1999, "Statistical criteria";
    /// TRF README "Detection component"): the heads (matches in k-runs) in the window of the last max(d, 20)
    /// positions must reach the sum-of-heads cut-off R(d,k,PM) AND the matches must be spread over the window —
    /// apparent-size criterion: the first k-tuple match lies no later than <see cref="TrfApparentSizeOffset"/>
    /// positions after the window's left end. When d fails on its own and d &gt; 20 (random-walk range
    /// d ± ⌊2.3·√(PI·d)⌋, parameter-set API only), d still qualifies if no active distance in the range has more
    /// heads than d and the heads summed over d and the lower range — or over any window of the same width
    /// sliding up through the upper range — reach the cut-off while some summed distance carrying at least 35 %
    /// of min(cut-off, heads at d) meets the apparent-size test. Only distances that have themselves reached a
    /// criteria test take part ("active"); a scanned distance whose window has emptied leaves the active set.
    /// </summary>
    private static bool MeetsTrfCriteria(TupleMatchWindow?[] windows, bool[] linked, int i, int d, int k, TrfModel m)
    {
        var main = windows[d]!;
        int criterion = TrfSumOfHeadsCriterion(d, m.Pm);
        int maxFirstMatch = Math.Max(0, i - Math.Max(d, TrfMinDistanceWindow)) + TrfApparentSizeOffset(d, m.Pm);
        int mainHeads = main.Heads;
        int rangeMinHeads = (int)(0.35 * Math.Min(criterion, mainHeads));

        bool spreadOk = main.FirstMatch(k) <= maxFirstMatch;
        int spreadDistance = spreadOk ? d : 0;
        if (mainHeads >= criterion && spreadOk)
            return true;

        int radius = m.DistanceRadius(d);
        int low = Math.Max(d - radius, 1);
        int high = Math.Min(d + radius, m.MaxDistance);

        // Active, non-empty distance t (purged to the window at i); an emptied distance becomes inactive.
        int ActiveHeads(int t)
        {
            if (!linked[t])
                return 0;
            var w = windows[t]!;
            w.Purge(TrfWindowLeft(i, t));
            if (w.Heads == 0)
                linked[t] = false;
            return w.Heads;
        }

        // Lower range, scanned downwards: d must stay the best distance.
        int sum = mainHeads;
        int lowPointer = d;
        for (int t = d - 1; t >= low; t--)
        {
            int heads = ActiveHeads(t);
            if (heads == 0)
                continue;
            if (heads > mainHeads)
                return false;
            lowPointer = t;
            sum += heads;
            if (!spreadOk && heads >= rangeMinHeads && windows[t]!.FirstMatch(k) <= maxFirstMatch)
            {
                spreadOk = true;
                spreadDistance = t;
            }
        }

        // Upper range: d must stay the best distance.
        for (int t = d + 1; t <= high; t++)
        {
            if (ActiveHeads(t) > mainHeads)
                return false;
        }

        if (sum >= criterion && spreadOk)
            return true;

        // Slide a window of the lower range's width up through the upper range.
        int width = d - low + 1;
        for (int t = d + 1; t <= high; t++)
        {
            if (!linked[t])
                continue;
            int heads = windows[t]!.Heads;
            sum += heads;
            while (lowPointer < t - width + 1)
            {
                sum -= windows[lowPointer]!.Heads;
                if (spreadOk && lowPointer == spreadDistance)
                {
                    spreadOk = false;
                    spreadDistance = 0;
                }
                int next = lowPointer + 1;
                while (!linked[next])
                    next++;
                lowPointer = next;
            }

            if (heads >= rangeMinHeads && windows[t]!.FirstMatch(k) <= maxFirstMatch)
            {
                spreadOk = true;
                spreadDistance = t;
            }

            if (sum >= criterion && spreadOk)
                return true;
        }

        return false;
    }

    /// <summary>Distances above this use TRF's best-period list (TRF 4.10.0: d &gt; 250).</summary>
    private const int TrfBestPeriodListMinDistance = 250;

    /// <summary>A region analysed earlier and its five best periods (TRF best-period list).</summary>
    private sealed class TrfBestPeriodEntry(int low, int high, int[] best)
    {
        public int Low { get; } = low;
        public int High { get; set; } = high;
        public int[] Best { get; } = best;
    }

    /// <summary>
    /// Best-period list test for a candidate distance d &gt; 250 ending at i: when an earlier analysed region spans
    /// the candidate's two-copy extent (from i − 2d + 1 + the apparent-size offset up to i), d must be one of that
    /// region's five best periods (see <see cref="TrfBestPeriods"/>); with no spanning region the candidate passes.
    /// This keeps large multiples of an already analysed shorter period from being re-aligned.
    /// </summary>
    private static bool IsAllowedByTrfBestPeriodList(List<TrfBestPeriodEntry> list, int i, int d, TrfModel m)
    {
        int cutoff = i - 2 * m.MaxDistance;
        list.RemoveAll(e => e.High < cutoff);
        int need = i - 2 * d + 1 + TrfApparentSizeOffset(d, m.Pm);
        bool covered = false;
        foreach (var entry in list)
        {
            if (entry.Low <= need && entry.High >= i)
            {
                covered = true;
                if (Array.IndexOf(entry.Best, d) >= 0)
                    return true;
            }
        }

        return !covered;
    }

    /// <summary>
    /// Sliding window of k-tuple matches at one distance d (TRF distance list): runs of adjacent tuple
    /// matches, each stored as (end position, number of heads); runs whose end falls before the window's
    /// left end are dropped. <see cref="Heads"/> is the sum of heads in k-runs inside the window.
    /// </summary>
    private sealed class TupleMatchWindow(int capacity)
    {
        private readonly int[] _end = new int[capacity];
        private readonly int[] _size = new int[capacity];
        private int _head;
        private int _count;

        public int Heads { get; private set; }

        /// <summary>End of the first k-tuple match of the oldest run in the window (its end − its heads + k).</summary>
        public int FirstMatch(int tupleSize) => _end[_head] - _size[_head] + tupleSize;

        /// <summary>Drops runs that end before <paramref name="windowLeft"/>.</summary>
        public void Purge(int windowLeft)
        {
            while (_count > 0 && _end[_head] < windowLeft)
            {
                Heads -= _size[_head];
                _head = (_head + 1) % _end.Length;
                _count--;
            }
        }

        public void Add(int position, int tupleSize, int windowLeft)
        {
            Purge(windowLeft);

            int tail = (_head + _count - 1 + _end.Length) % _end.Length;
            if (_count > 0 && _end[tail] == position - 1)
            {
                _end[tail] = position;
                _size[tail]++;
                Heads++;
                return;
            }

            tail = (_head + _count) % _end.Length;
            _end[tail] = position;
            _size[tail] = tupleSize;
            _count++;
            Heads += tupleSize;
        }
    }

    /// <summary>
    /// TRF analysis of one candidate (position <paramref name="i"/>, distance <paramref name="d"/>): align the
    /// candidate pattern S[i−d+1..i], apply the copy-number and three-best-periods tests, build the consensus,
    /// realign against it and report the statistics. Updates <paramref name="seenEnd"/>[d] with the end of each
    /// alignment so that later matches at d inside the aligned region are not re-analysed (TRF).
    /// </summary>
    private static ApproximateTandemRepeatResult? AnalyzeTrfCandidate(
        TrfModel m, char[] s, int n, int i, int d, int[] seenEnd,
        Dictionary<(int First, int Last), int[]> bestPeriods, List<TrfBestPeriodEntry> bestPeriodList)
    {
        var pattern = new char[d];
        Array.Copy(s, i - d + 1, pattern, 0, d);

        int first, last;
        double copies;
        TrfAlignment? initial = null;
        if (d <= TrfSmallDistance)
        {
            // Extent-only pass first: the multiples test needs only the aligned region, so the traceback
            // matrix is built only for candidates that pass both tests (identical path either way).
            var extent = TrfWraparoundExtent(m, s, n, i, pattern);
            if (extent is null)
                return null;
            (first, last, copies) = (extent.Value.First, extent.Value.Last, extent.Value.CopyNumber);
        }
        else
        {
            initial = TrfBandAlign(m, s, n, i, pattern);
            if (initial is null)
                return null;
            (first, last, copies) = (initial.First, initial.Last, initial.CopyNumber);
        }

        MarkAligned(seenEnd, d, last);
        if (!MeetsTrfCopyNumber(copies, d, d))
            return null;

        int[]? best = null;
        if (d != 1)
        {
            best = GetTrfBestPeriods(s, first, last, m.MaxDistance, bestPeriods);
            bestPeriodList.Add(new TrfBestPeriodEntry(first, last, best));
        }

        if (!IsAmongTrfBestPeriods(s, first, last, d, best))
            return null;

        var alignment = initial ?? TrfWraparoundAlign(m, s, n, i, pattern);
        if (alignment is null)
            return null;
        char[] consensus = TrfConsensus(alignment.Columns, d);
        if (consensus.Length == 0)
            return null;

        var final = consensus.Length <= TrfSmallDistance
            ? TrfWraparoundAlign(m, s, n, i, consensus)
            : TrfBandAlign(m, s, n, i, consensus);
        if (final is null)
            return null;
        MarkAligned(seenEnd, d, final.Last);
        if (best is not null && bestPeriodList[^1].High > final.Last)
            bestPeriodList[^1].High = final.Last;
        // TRF quirk kept: the 50 < size <= 100 ramp uses the candidate distance d.
        if (!MeetsTrfCopyNumber(final.CopyNumber, consensus.Length, d) || final.Score < m.MinScore)
            return null;

        return TrfStatistics(m, s, n, final, consensus);
    }

    /// <summary>Records that sequence positions up to <paramref name="last"/> were aligned at distance d.</summary>
    private static void MarkAligned(int[] seenEnd, int d, int last) => seenEnd[d] = last;

    /// <summary>
    /// TRF minimum copy number: ≥ 1.9 copies for patterns ≤ 50, ramping from 1.9 down to 1.8 for 50..100,
    /// ≥ 1.8 above 100 (TRF 4.10.0; Benson 1999: "If at least two copies of the pattern are aligned with the
    /// sequence, the tandem repeat is reported").
    /// </summary>
    private static bool MeetsTrfCopyNumber(double copies, int size, int d)
    {
        if (size <= 50)
            return copies >= 1.9;
        if (size <= 100)
            return copies >= 1.9 - 0.002 * (d - 50);
        return copies >= 1.8;
    }

    /// <summary>Result of a WDP fill: best score, its cell, and (extent mode) the traced path's start and
    /// consumed pattern columns.</summary>
    private readonly record struct TrfFill(int Score, int RealRow, int Row, int Col, int First, int Consumed);

    /// <summary>
    /// Runs the TRF wraparound DP for <paramref name="pattern"/> around candidate end <paramref name="start"/>:
    /// a backward local scan locates the leftmost row reaching the best score, then a forward local alignment
    /// starting one pattern length before it yields the optimum (first strictly greatest cell, row-major).
    /// Each row is computed in two passes because the horizontal dependency wraps from the last pattern column
    /// to the first; zero cells beyond the candidate are killed so the alignment cannot restart there.
    /// With <paramref name="rows"/> the final forward rows are stored for traceback; without, every cell carries
    /// the start row and consumed pattern columns of the path the traceback would follow (same predecessor
    /// preference — diagonal, then vertical, then horizontal — on the final values), so the extent is known
    /// without the O(rows × pattern) matrix.
    /// </summary>
    private static TrfFill TrfWraparoundFill(TrfModel m, char[] s, int n, int start, char[] pattern, List<int[]>? rows)
    {
        int size = pattern.Length;
        var weights = TrfWeightRows(m, pattern);
        int indel = m.Indel;
        var up = new int[size];
        var diag = new int[size];
        var cur = new int[size];

        // Backward scan (sequence read right-to-left from the candidate end, pattern read in reverse).
        Array.Fill(up, indel);
        int maxScore = 0;
        int minRow = start;
        int killBelow = start - Math.Max(size, TrfMinDistanceWindow);
        int realRow = start + 1;
        bool endOfTrace = false;
        int backwardRows = 0;
        while (!endOfTrace && realRow > 1 && backwardRows < m.MaxWrapLength)
        {
            backwardRows++;
            realRow--;
            int[] w = weights[TrfSymbolClass(s[realRow])];
            int left = indel;
            for (int c = size - 1; c >= 0; c--)
            {
                diag[c] += w[c];
                left = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left)) + indel;
            }

            endOfTrace = true;
            for (int c = size - 1; c >= 0; c--)
            {
                int v = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left));
                left = up[c] = v + indel;
                if (realRow <= killBelow && v == 0)
                {
                    v = TrfDeadCell;
                    left = up[c] = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                cur[c] = v;
                if (v >= maxScore)
                {
                    maxScore = v;
                    minRow = realRow;
                }
            }

            for (int c = 0; c < size - 1; c++)
                diag[c] = cur[c + 1];
            diag[size - 1] = cur[0];
        }

        // Forward local alignment from one pattern length before the leftmost best row.
        Array.Fill(up, indel);
        Array.Clear(diag);
        bool trackExtent = rows is null;
        var prev = new int[size];
        var prevFirst = trackExtent ? new int[size] : [];
        var prevCols = trackExtent ? new int[size] : [];
        var curFirst = trackExtent ? new int[size] : [];
        var curCols = trackExtent ? new int[size] : [];
        var pending = trackExtent ? new bool[size] : [];
        rows?.Add(new int[size]);

        realRow = Math.Max(minRow - size - 1, 0);
        int row = 0;
        maxScore = 0;
        int bestReal = -1, bestRow = -1, bestCol = -1, bestFirst = 0, bestConsumed = 0;
        endOfTrace = false;
        while (!endOfTrace && realRow < n && row < m.MaxWrapLength)
        {
            row++;
            realRow++;
            int[] w = weights[TrfSymbolClass(s[realRow])];
            int[] values = rows is null ? cur : new int[size];
            int left = indel;
            for (int c = 0; c < size; c++)
            {
                diag[c] += w[c];
                left = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left)) + indel;
            }

            endOfTrace = true;
            int rowMax = 0, rowMaxCol = -1;
            for (int c = 0; c < size; c++)
            {
                int v = Math.Max(Math.Max(0, diag[c]), Math.Max(up[c], left));
                left = up[c] = v + indel;
                if (realRow > start && v == 0)
                {
                    v = TrfDeadCell;
                    left = up[c] = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                values[c] = v;
                if (v > maxScore && v > rowMax)
                {
                    rowMax = v;
                    rowMaxCol = c;
                }
            }

            if (trackExtent)
                TrackExtent(indel, values, prev, w, realRow, prevFirst, prevCols, curFirst, curCols, pending);

            if (rowMaxCol >= 0)
            {
                maxScore = rowMax;
                bestReal = realRow;
                bestRow = row;
                bestCol = rowMaxCol;
                if (trackExtent)
                {
                    bestFirst = curFirst[rowMaxCol];
                    bestConsumed = curCols[rowMaxCol];
                }
            }

            for (int c = size - 1; c > 0; c--)
                diag[c] = values[c - 1];
            diag[0] = values[size - 1];

            if (rows is null)
            {
                Array.Copy(values, prev, size);
                Array.Copy(curFirst, prevFirst, size);
                Array.Copy(curCols, prevCols, size);
            }
            else
            {
                rows.Add(values);
            }
        }

        return new TrfFill(maxScore, bestReal, bestRow, bestCol, bestFirst, bestConsumed);
    }

    /// <summary>Per-row extent bookkeeping for <see cref="TrfWraparoundFill"/> (start row and consumed pattern
    /// columns inherited along the traceback's preferred predecessor).</summary>
    private static void TrackExtent(
        int indel, int[] values, int[] prev, int[] w, int realRow,
        int[] prevFirst, int[] prevCols, int[] curFirst, int[] curCols, bool[] pending)
    {
        int size = values.Length;
        for (int c = 0; c < size; c++)
        {
            pending[c] = false;
            int v = values[c];
            if (v <= 0)
                continue;
            int jp = c == 0 ? size - 1 : c - 1;
            if (v == prev[jp] + w[c])
            {
                bool starts = prev[jp] <= 0;
                curFirst[c] = starts ? realRow : prevFirst[jp];
                curCols[c] = (starts ? 0 : prevCols[jp]) + 1;
            }
            else if (v == prev[c] + indel)
            {
                curFirst[c] = prevFirst[c];
                curCols[c] = prevCols[c];
            }
            else
            {
                pending[c] = true; // horizontal move: inherits from its left neighbour in this row
            }
        }

        // Horizontal chains may wrap from the last pattern column to the first: two sweeps resolve them.
        for (int sweep = 0; sweep < 2; sweep++)
        {
            for (int c = 0; c < size; c++)
            {
                int jp = c == 0 ? size - 1 : c - 1;
                if (pending[c] && !pending[jp] && values[jp] > 0)
                {
                    curFirst[c] = curFirst[jp];
                    curCols[c] = curCols[jp] + 1;
                    pending[c] = false;
                }
            }
        }
    }

    /// <summary>Weight of each pattern column against A, C, G, T and any other symbol (TRF: +2 for an
    /// identical A/C/G/T pair, −7 otherwise; N never matches).</summary>
    private static int[][] TrfWeightRows(TrfModel m, char[] pattern)
    {
        var weights = new int[5][];
        for (int symbol = 0; symbol < 5; symbol++)
        {
            weights[symbol] = new int[pattern.Length];
            for (int c = 0; c < pattern.Length; c++)
                weights[symbol][c] = symbol < 4 && TrfSymbolClass(pattern[c]) == symbol ? m.Match : m.Mismatch;
        }
        return weights;
    }

    /// <summary>A/C/G/T → 0..3, any other symbol → 4 (never matches).</summary>
    private static int TrfSymbolClass(char c)
    {
        int code = AcgtCode(c);
        return code < 0 ? 4 : code;
    }

    /// <summary>WDP optimum and its extent (first/last sequence index, copy number) without traceback.</summary>
    private static TrfExtent? TrfWraparoundExtent(TrfModel m, char[] s, int n, int start, char[] pattern)
    {
        var fill = TrfWraparoundFill(m, s, n, start, pattern, rows: null);
        if (fill.Score <= 0)
            return null;
        return new TrfExtent(fill.Score, fill.First, fill.RealRow, TrfCopyNumber(fill.Consumed, pattern.Length));
    }

    /// <summary>TRF copy number: full passes through the pattern plus the final partial pass, i.e. aligned
    /// pattern columns / pattern size (summed as TRF does, integer part then fraction).</summary>
    private static double TrfCopyNumber(int consumedColumns, int size) =>
        consumedColumns / size + (double)(consumedColumns % size) / size;

    /// <summary>
    /// Runs the TRF wraparound DP and traces the optimal local alignment back from its best cell
    /// (predecessor preference: match/mismatch, then a sequence symbol against a gap, then a pattern symbol
    /// against a gap). Columns are returned rightmost first.
    /// </summary>
    private static TrfAlignment? TrfWraparoundAlign(TrfModel m, char[] s, int n, int start, char[] pattern)
    {
        int size = pattern.Length;
        var rows = new List<int[]>();
        var fill = TrfWraparoundFill(m, s, n, start, pattern, rows);
        if (fill.Score <= 0)
            return null;

        var columns = new List<TrfColumn>();
        int i = fill.RealRow, r = fill.Row, j = fill.Col, consumed = 0;
        while (rows[r][j] > 0)
        {
            int v = rows[r][j];
            int jp = j == 0 ? size - 1 : j - 1;
            if (v == rows[r - 1][jp] + TrfWeight(m, s[i], pattern[j]))
            {
                columns.Add(new TrfColumn(s[i], pattern[j], i, j));
                consumed++;
                i--;
                r--;
                j = jp;
            }
            else if (v == rows[r - 1][j] + m.Indel)
            {
                columns.Add(new TrfColumn(s[i], '-', i, (j + 1) % size));
                i--;
                r--;
            }
            else if (v == rows[r][jp] + m.Indel)
            {
                columns.Add(new TrfColumn('-', pattern[j], i + 1, j));
                consumed++;
                j = jp;
            }
            else
            {
                throw new InvalidOperationException("Wraparound DP traceback is inconsistent.");
            }
        }

        return new TrfAlignment(fill.Score, TrfCopyNumber(consumed, size), columns.ToArray());
    }

    /// <summary>Smallest narrow-band radius (TRF 4.10.0: 6).</summary>
    private const int TrfMinBandRadius = 6;

    /// <summary>Consecutive diagonal row maxima that recentre the band (TRF 4.10.0: 3).</summary>
    private const int TrfBandRecenterRun = 3;

    /// <summary>C-style signed shift of a column difference: the residue in (−size/2, size/2].</summary>
    private static int BandShift(int difference, int size)
    {
        int k = ((difference % size) + size) % size;
        return size - k <= k ? -(size - k) : k;
    }

    private static int Wrap(int value, int size) => ((value % size) + size) % size;

    /// <summary>
    /// Narrow-band wraparound DP for patterns longer than 20 (Benson 1999 / TRF README "Narrow Band Alignment": "we
    /// limit WDP calculations to a narrow diagonal band in the alignment matrix for patterns larger than 20
    /// characters. In accordance with the random walk results, the band radius is Δd_max. The band is periodically
    /// recentered around a run of matches in the current best alignment"). Each DP row keeps 2w+1 cells around a
    /// band centre column that advances one pattern column per row and is moved onto the row's best cell once the
    /// row maxima have followed the diagonal through matches for 3 consecutive rows. As in the full WDP, a backward
    /// local scan from the candidate end (radius w = max(6, Δd_max)) finds the leftmost best cell; the forward local
    /// alignment then starts there, with its band anchored on that cell and widened to min(2w, ⌊size/3⌋) (TRF 4.04:
    /// "widened radius of narrowband alignment"). Optimum = last cell reaching the maximum; the traceback continues
    /// through zero-valued cells that are genuine continuations of the path (TRF 4.07b: "changed alignment to go
    /// further when score drops to 0").
    /// </summary>
    private static TrfAlignment? TrfBandAlign(TrfModel m, char[] s, int n, int start, char[] pattern)
    {
        int size = pattern.Length;
        var weights = TrfWeightRows(m, pattern);
        int indel = m.Indel;
        int backRadius = Math.Max(TrfMinBandRadius, m.RandomWalkRange(size));
        int foreRadius = Math.Min(2 * backRadius, size / 3);

        // ---- backward scan (sequence right to left from the candidate end, pattern read in reverse) ----
        int w = backRadius, width = 2 * w + 1;
        var diag = new int[width];
        var up = new int[width];
        var cur = new int[width];
        for (int j = 0; j < width; j++)
        {
            cur[j] = j <= w ? indel * (w - j) : TrfDeadCell;
            diag[j] = cur[j];
            up[j] = cur[j] + indel;
        }

        int prevCenter = 0, matchCol = -2, diagonalRun = 0;
        int maxScore = 0, minRealRow = start, minCenter = 0, minPosition = w;
        int realRow = start + 1, rows = 0;
        bool endOfTrace = false;
        while (!endOfTrace && realRow > 1 && rows < m.MaxWrapLength)
        {
            rows++;
            realRow--;
            int lastMatchCol = matchCol;
            int center = diagonalRun >= TrfBandRecenterRun ? Wrap(matchCol - 1, size) : Wrap(prevCenter - 1, size);
            int shift = BandShift(center - prevCenter, size);
            int[] wt = weights[TrfSymbolClass(s[realRow])];
            int left = TrfDeadCell, rowMax = -1;
            endOfTrace = true;
            for (int j = 2 * w; j >= 0; j--)
            {
                int col = Wrap(center - w + j, size);
                int di = j + shift + 1, ui = j + shift;
                bool hasDiag = di >= 0 && di <= 2 * w;
                int dv = hasDiag ? diag[di] + wt[col] : 0;
                int v = Math.Max(0, left);
                if (hasDiag)
                    v = Math.Max(v, dv);
                if (ui >= 0 && ui <= 2 * w)
                    v = Math.Max(v, up[ui]);
                left = v + indel;
                if (realRow <= start - size && v == 0)
                {
                    v = TrfDeadCell;
                    left = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                cur[j] = v;
                if (v >= maxScore)
                {
                    maxScore = v;
                    minRealRow = realRow;
                    minCenter = center;
                    minPosition = j;
                }

                if (v > rowMax)
                {
                    // A cell without a diagonal predecessor resets the match column but not the row maximum (TRF).
                    if (hasDiag)
                    {
                        rowMax = v;
                        matchCol = v == dv && wt[col] == m.Match ? col : -2;
                    }
                    else
                    {
                        matchCol = -2;
                    }
                }
            }

            for (int j = 0; j < width; j++)
            {
                diag[j] = cur[j];
                up[j] = cur[j] + indel;
            }

            diagonalRun = (matchCol - lastMatchCol + size) % size == size - 1 ? diagonalRun + 1 : 0;
            prevCenter = center;
        }

        // ---- forward local alignment from the leftmost best cell of the backward scan ----
        int zeroAt = Math.Max(0, minPosition - w + foreRadius);
        w = foreRadius;
        width = 2 * w + 1;
        var rowValues = new List<int[]>();
        var centers = new List<int>();
        var row0 = new int[width];
        for (int j = 0; j < width; j++)
            row0[j] = j < zeroAt ? TrfDeadCell : indel * (j - zeroAt);
        rowValues.Add(row0);
        centers.Add(Wrap(minCenter - 1, size));
        diag = (int[])row0.Clone();
        up = new int[width];
        for (int j = 0; j < width; j++)
            up[j] = row0[j] + indel;

        prevCenter = centers[0];
        matchCol = -2;
        diagonalRun = 0;
        maxScore = 0;
        int maxRealRow = -1, maxRow = -1, maxCol = -1;
        realRow = minRealRow - 1;
        int r = 0;
        endOfTrace = false;
        while (!endOfTrace && realRow < n && r < m.MaxWrapLength)
        {
            r++;
            realRow++;
            int lastMatchCol = matchCol;
            int center = diagonalRun >= TrfBandRecenterRun ? Wrap(matchCol + 1, size) : Wrap(prevCenter + 1, size);
            int shift = BandShift(center - prevCenter, size);
            int[] wt = weights[TrfSymbolClass(s[realRow])];
            var values = new int[width];
            int left = TrfDeadCell, rowMax = -1;
            endOfTrace = true;
            for (int j = 0; j < width; j++)
            {
                int col = Wrap(center - w + j, size);
                int di = j + shift - 1, ui = j + shift;
                bool hasDiag = di >= 0 && di <= 2 * w;
                int dv = hasDiag ? diag[di] + wt[col] : 0;
                int v = Math.Max(0, left);
                if (hasDiag)
                    v = Math.Max(v, dv);
                if (ui >= 0 && ui <= 2 * w)
                    v = Math.Max(v, up[ui]);
                left = v + indel;
                if (realRow >= start && v == 0)
                {
                    v = TrfDeadCell;
                    left = TrfDeadCell;
                }
                else
                {
                    endOfTrace = false;
                }

                values[j] = v;
                if (v >= maxScore)
                {
                    maxScore = v;
                    maxRealRow = realRow;
                    maxRow = r;
                    maxCol = col;
                }

                if (v > rowMax)
                {
                    if (hasDiag)
                    {
                        rowMax = v;
                        matchCol = v == dv && wt[col] == m.Match ? col : -2;
                    }
                    else
                    {
                        matchCol = -2;
                    }
                }
            }

            rowValues.Add(values);
            centers.Add(center);
            for (int j = 0; j < width; j++)
            {
                diag[j] = values[j];
                up[j] = values[j] + indel;
            }

            // TRF compares raw (C-style) residues here, so a non-match row (−2) can also extend the run.
            diagonalRun = (matchCol - lastMatchCol + size) % size == 1 ? diagonalRun + 1 : 0;
            prevCenter = center;
        }

        if (maxScore <= 0 || maxRow < 1)
            return null;

        return TrfBandTraceback(m, s, pattern, rowValues, centers, w, maxScore, maxRealRow, maxRow, maxCol);
    }

    /// <summary>
    /// Traceback of a narrow-band WDP from its optimum (predecessor preference as in the full WDP: match/mismatch,
    /// then sequence symbol against a gap, then pattern symbol against a gap). A zero cell ends the alignment unless
    /// it is reached from a predecessor on the path (a genuine zero, not a local-alignment restart); row 0 (the
    /// anchor row) always ends it. Columns are returned rightmost first.
    /// </summary>
    private static TrfAlignment? TrfBandTraceback(
        TrfModel m, char[] s, char[] pattern, List<int[]> rowValues, List<int> centers, int w,
        int score, int maxRealRow, int maxRow, int maxCol)
    {
        int size = pattern.Length;
        int r = maxRow, realRow = maxRealRow, c = maxCol;
        int i = w + BandShift(c - centers[r], size);
        int consumed = 0;
        var columns = new List<TrfColumn>();

        int Cell(int row, int index) => index >= 0 && index <= 2 * w ? rowValues[row][index] : int.MinValue / 2;

        while (r != 0)
        {
            int v = rowValues[r][i];
            int shift = BandShift(centers[r] - centers[r - 1], size);
            int upi = w + BandShift(c - centers[r - 1], size);

            // Moves available at band position i (TRF's band-edge cases).
            bool canDiag, canUp, canLeft = i > 0;
            if (shift >= 1)
            {
                canDiag = i <= 2 * w - shift + 1;
                canUp = i <= 2 * w - shift;
            }
            else
            {
                canDiag = i > -shift;
                canUp = i >= -shift;
            }

            int weight = TrfWeight(m, s[realRow], pattern[c]);
            bool fromDiag = canDiag && v == Cell(r - 1, upi - 1) + weight;
            bool fromUp = !fromDiag && canUp && v == Cell(r - 1, upi) + m.Indel;
            bool fromLeft = !fromDiag && !fromUp && canLeft && v == Cell(r, i - 1) + m.Indel;
            if (v < 0 || (v == 0 && !(fromDiag || fromUp || fromLeft)))
                break;

            if (fromDiag)
            {
                columns.Add(new TrfColumn(s[realRow], pattern[c], realRow, c));
                consumed++;
                c = Wrap(c - 1, size);
                realRow--;
                r--;
                i = upi - 1;
            }
            else if (fromUp)
            {
                columns.Add(new TrfColumn(s[realRow], '-', realRow, (c + 1) % size));
                realRow--;
                r--;
                i = upi;
            }
            else if (fromLeft)
            {
                columns.Add(new TrfColumn('-', pattern[c], realRow + 1, c));
                consumed++;
                c = Wrap(c - 1, size);
                i--;
            }
            else
            {
                break; // inconsistent cell (TRF reports a traceback error and stops there)
            }
        }

        return columns.Count == 0 ? null : new TrfAlignment(score, TrfCopyNumber(consumed, size), columns.ToArray());
    }

    /// <summary>
    /// Majority-rule consensus from an alignment against a pattern of <paramref name="patternLength"/>
    /// (Benson 1999: "we determine a consensus pattern by majority rule from the alignment of the copies with
    /// P"). Slot 2c+1 is pattern position c, slot 2c the insertion point before it. A position takes the most
    /// frequent aligned symbol, or is deleted when gaps are at least as frequent (ties: gap, A, C, G, T); an
    /// insertion point receives its most frequent inserted base when insertions occur there in at least half
    /// of the passes. Follows the TRF 4.10.0 (non-weighted) consensus rule.
    /// </summary>
    private static char[] TrfConsensus(TrfColumn[] columns, int patternLength)
    {
        int slots = 2 * patternLength + 1;
        var counts = new int[5, slots]; // A C G T gap
        var inserts = new int[slots];
        var passes = new int[slots];

        int last = -1;
        int k = 0;
        while (k < columns.Length)
        {
            int index = columns[k].PatIndex;
            if (index != last)
            {
                int symbol = ConsensusSymbol(columns[k].Seq);
                if (symbol >= 0)
                    counts[symbol, 2 * index + 1]++;
                if (last != -1)
                    passes[index == patternLength - 1 ? 0 : 2 * index + 2]++;
                last = index;
                k++;
            }
            else
            {
                inserts[2 * index]++;
                while (k < columns.Length && columns[k].PatIndex == last)
                {
                    int symbol = ConsensusSymbol(columns[k].Seq);
                    if (symbol >= 0)
                        counts[symbol, 2 * index]++;
                    k++;
                }
            }
        }

        const string Bases = "ACGT";
        var consensus = new List<char>(patternLength + 4);
        for (int slot = 0; slot < slots; slot++)
        {
            if (slot % 2 == 1)
            {
                int best = counts[4, slot];
                char chosen = '-';
                for (int b = 0; b < 4; b++)
                {
                    if (counts[b, slot] > best)
                    {
                        best = counts[b, slot];
                        chosen = Bases[b];
                    }
                }
                if (chosen != '-')
                    consensus.Add(chosen);
            }
            else if (passes[slot] != 0 && (float)inserts[slot] / passes[slot] >= 0.5f)
            {
                int best = counts[0, slot];
                char chosen = 'A';
                for (int b = 1; b < 4; b++)
                {
                    if (counts[b, slot] > best)
                    {
                        best = counts[b, slot];
                        chosen = Bases[b];
                    }
                }
                consensus.Add(chosen);
            }
        }

        return consensus.ToArray();
    }

    /// <summary>Consensus count slot: A/C/G/T → 0..3, gap → 4, other symbols are not counted (−1).</summary>
    private static int ConsensusSymbol(char c) => c == '-' ? 4 : AcgtCode(c);

    /// <summary>Match / mismatch / indel counts between adjacent copies and the TRF period (most common
    /// distance between matching characters).</summary>
    private readonly record struct TrfCopyComparison(int Matches, int Mismatches, int Indels, int Period);

    /// <summary>
    /// Compares ADJACENT copies through the consensus alignment (TRF: statistics refer to "the matches,
    /// mismatches and indels overall between adjacent copies in the sequence, not between the sequence and the
    /// consensus pattern"): two cursors one consensus period apart walk the alignment; aligned symbols of the
    /// two copies are a match or a mismatch, a symbol against a gap is an indel. The period is "the most
    /// common matching distance between corresponding characters in the alignment" (ties: smallest).
    /// </summary>
    private static TrfCopyComparison CompareAdjacentCopies(TrfColumn[] columns)
    {
        int length = columns.Length;
        int lp = 0;
        while (lp < length && columns[lp].Pat == '-')
            lp++;
        int rp = lp + 1;
        while (rp < length && columns[rp].PatIndex != columns[lp].PatIndex)
            rp++;
        while (rp < length && columns[rp].Pat == '-')
            rp++;
        if (rp >= length)
            return new TrfCopyComparison(0, 0, 0, 0);

        int matches = 0, mismatches = 0, indels = 0;
        var distances = new Dictionary<int, int>();
        void Match(int a, int b)
        {
            matches++;
            int distance = Math.Abs(columns[b].SeqIndex - columns[a].SeqIndex);
            distances[distance] = distances.GetValueOrDefault(distance) + 1;
        }

        while (rp < length && lp < rp)
        {
            bool leftGapPat = columns[lp].Pat == '-';
            bool rightGapPat = columns[rp].Pat == '-';
            if (!leftGapPat && !rightGapPat)
            {
                bool leftGapSeq = columns[lp].Seq == '-';
                bool rightGapSeq = columns[rp].Seq == '-';
                if (!leftGapSeq && !rightGapSeq)
                {
                    if (columns[lp].Seq == columns[rp].Seq) Match(lp, rp);
                    else mismatches++;
                }
                else if (leftGapSeq != rightGapSeq)
                {
                    indels++;
                }
                lp++;
                rp++;
            }
            else if (leftGapPat && rightGapPat)
            {
                if (columns[lp].Seq == columns[rp].Seq) Match(lp, rp);
                else mismatches++;
                lp++;
                rp++;
            }
            else if (leftGapPat)
            {
                indels++;
                lp++;
            }
            else
            {
                indels++;
                rp++;
            }
        }

        int period = 0, bestCount = 0;
        foreach (var (distance, count) in distances.OrderBy(kv => kv.Key))
        {
            if (count > bestCount)
            {
                bestCount = count;
                period = distance;
            }
        }

        return new TrfCopyComparison(matches, mismatches, indels, period);
    }

    /// <summary>Builds the reported TRF statistics from the final (consensus) alignment.</summary>
    private static ApproximateTandemRepeatResult TrfStatistics(TrfModel m, char[] s, int n, TrfAlignment alignment, char[] consensus)
    {
        var copies = CompareAdjacentCopies(alignment.Columns);
        int trials = copies.Matches + copies.Mismatches + copies.Indels;

        int first = alignment.First, last = alignment.Last, span = last - first + 1;
        var region = new string(s, first, span);
        int a = 0, c = 0, g = 0, t = 0;
        foreach (char ch in region)
        {
            switch (ch)
            {
                case 'A': a++; break;
                case 'C': c++; break;
                case 'G': g++; break;
                case 'T': t++; break;
            }
        }

        // TRF prints the consensus starting at the pattern position aligned with the first repeat base.
        int phase = alignment.Columns[^1].PatIndex;
        string rotated = new string(consensus, phase, consensus.Length - phase) + new string(consensus, 0, phase);

        return new ApproximateTandemRepeatResult(
            Start: first - 1,
            SpanLength: span,
            Period: copies.Period,
            ConsensusSize: consensus.Length,
            Consensus: rotated,
            CopyNumber: alignment.CopyNumber,
            PercentMatches: trials > 0 ? 100.0 * copies.Matches / trials : 0.0,
            PercentIndels: trials > 0 ? 100.0 * copies.Indels / trials : 0.0,
            AlignmentScore: alignment.Score)
        {
            PercentA = 100.0 * a / span,
            PercentC = 100.0 * c / span,
            PercentG = 100.0 * g / span,
            PercentT = 100.0 * t / span,
            Entropy = SequenceComplexity.CalculateShannonEntropy(region),
            EntropyTrf = TrfEntropy(a, c, g, t, span),
            AlignedSequence = AlignmentRow(alignment.Columns, sequenceRow: true),
            AlignedConsensus = AlignmentRow(alignment.Columns, sequenceRow: false),
            LeftFlank = m.FlankLength > 0 ? new string(s, Math.Max(1, first - m.FlankLength), first - Math.Max(1, first - m.FlankLength)) : null,
            RightFlank = m.FlankLength > 0 ? new string(s, last + 1, Math.Min(n, last + m.FlankLength) - last) : null,
        };
    }

    /// <summary>
    /// TRF's "Entropy (0-2)" column: −Σ p_b·log₂ p_b over A, C, G, T with p_b = count_b / region length, where the
    /// region length counts EVERY symbol of the repeat (TRF 4.10.0 get_statistics: the denominator is the number of
    /// non-gap characters of the aligned sequence, so an N lowers every p_b and the four p_b sum to less than 1).
    /// Terms are summed from the largest p_b down, with log₂ x computed as ln x / ln 2, as TRF does.
    /// </summary>
    private static double TrfEntropy(int a, int c, int g, int t, int length)
    {
        Span<double> p = [(double)a / length, (double)c / length, (double)g / length, (double)t / length];
        p.Sort();
        double entropy = 0.0;
        for (int b = 3; b >= 0; b--)
            entropy += p[b] == 0 ? 0 : p[b] * (Math.Log(p[b]) / Math.Log(2));
        return Math.Abs(entropy);
    }

    /// <summary>One row of the final alignment, left to right ('-' = gap): the repeat or the consensus copies.</summary>
    private static string AlignmentRow(TrfColumn[] columns, bool sequenceRow)
    {
        var row = new char[columns.Length];
        for (int k = 0; k < columns.Length; k++)
        {
            var column = columns[columns.Length - 1 - k];
            row[k] = sequenceRow ? column.Seq : column.Pat;
        }
        return new string(row);
    }

    /// <summary>
    /// TRF multiples test: the candidate distance must be one of the three best periods of the aligned region
    /// (Benson 1999 / TRF README: redundant reporting at multiples of the pattern size is limited "to, at most,
    /// three pattern sizes"). Best periods = the largest counts of distances between identical dinucleotides in
    /// the region after removing their least-squares linear trend (TRF 4.10.0 method; non-ACGT symbols count as
    /// A as in TRF). Period 1 instead requires ≥ 80 % of the region to be one base.
    /// </summary>
    private static bool IsAmongTrfBestPeriods(char[] s, int first, int last, int d, int[]? best)
    {
        if (d == 1)
        {
            int length = last - first + 1;
            var composition = new int[4];
            for (int p = first; p <= last; p++)
                composition[AcgtIndexOrA(s[p])]++;
            return composition.Max() * 100.0f / length >= 80.0f;
        }

        return Array.IndexOf(best!, d, 0, TrfPeriodsTested) >= 0;
    }

    /// <summary>Number of best periods kept per region (TRF: 5, for the best-period list).</summary>
    private const int TrfPeriodsKept = 5;

    /// <summary>Number of best periods the multiples test accepts (TRF: 3).</summary>
    private const int TrfPeriodsTested = 3;

    /// <summary>The best periods of a region, cached by region (different distances often align the same region).</summary>
    private static int[] GetTrfBestPeriods(
        char[] s, int first, int last, int maxPeriod, Dictionary<(int First, int Last), int[]> cache)
    {
        if (!cache.TryGetValue((first, last), out int[]? best))
        {
            best = TrfBestPeriods(s, first, last - first + 1, maxPeriod);
            cache[(first, last)] = best;
        }

        return best;
    }

    /// <summary>The five best periods of s[first..first+length−1], best first (see <see cref="IsAmongTrfBestPeriods"/>).</summary>
    private static int[] TrfBestPeriods(char[] s, int first, int length, int maxPeriod)
    {
        var best = new int[TrfPeriodsKept];
        int end = length - 2;
        if (end < 1)
            return best;

        var counts = new double[length];
        var history = new int[length];
        var heads = new int[16];
        Array.Fill(heads, -1);
        const int MaxCountedDistance = 3 * MaxApproximatePeriod;
        for (int p = 0; p <= end; p++)
        {
            int tuple = AcgtIndexOrA(s[first + p]) * 4 + AcgtIndexOrA(s[first + p + 1]);
            history[p] = heads[tuple];
            heads[tuple] = p;
            int distance = 0;
            for (int cur = p; history[cur] != -1 && distance < MaxCountedDistance; cur = history[cur])
            {
                distance = p - history[cur];
                counts[distance] += 1.0;
            }
        }

        double xy = 0, x = 0, y = 0, x2 = 0;
        for (int q = 1; q <= end; q++)
        {
            xy += q * counts[q];
            x += q;
            y += counts[q];
            x2 += (double)q * q;
        }
        double slope = (end * xy - x * y) / (end * x2 - x * x);
        for (int q = 1; q <= end; q++)
            counts[q] -= q * slope;

        int top = Math.Min(end, maxPeriod);
        for (int pick = 0; pick < best.Length; pick++)
        {
            int bestIndex = 0;
            double bestValue = 0.0;
            for (int q = 1; q <= top; q++)
            {
                if (counts[q] > bestValue)
                {
                    bestIndex = q;
                    bestValue = counts[q];
                }
            }
            best[pick] = bestIndex;
            counts[bestIndex] = 0.0;
        }

        return best;
    }

    /// <summary>A/C/G/T → 0..3 with any other symbol counted as A (TRF's zero-initialised index table).</summary>
    private static int AcgtIndexOrA(char c) => Math.Max(0, AcgtCode(c));

    /// <summary>
    /// TRF redundancy elimination over repeats sorted by start: of two repeats overlapping by ≥ 90 % of one of
    /// them, that one is dropped when it has the same period and no higher score, or a period that is a
    /// multiple of the other's and a score ≤ 1.1× (TRF README "Redundancy": the same repeat detected at
    /// several period sizes / "the same period size may be detected more than once").
    /// </summary>
    private static List<ApproximateTandemRepeatResult> RemoveTrfRedundancy(List<ApproximateTandemRepeatResult> repeats)
    {
        int i = 0;
        while (i < repeats.Count)
        {
            bool removedI = false;
            int j = i + 1;
            while (j < repeats.Count)
            {
                var a = repeats[i];
                var b = repeats[j];
                int overlap = Math.Min(a.Start + a.SpanLength, b.Start + b.SpanLength) - Math.Max(a.Start, b.Start);
                if (overlap <= 0)
                    break;
                if (!(overlap / (double)a.SpanLength < 0.9) && IsTrfRedundant(a, b))
                {
                    repeats.RemoveAt(i);
                    removedI = true;
                    break;
                }
                if (!(overlap / (double)b.SpanLength < 0.9) && IsTrfRedundant(b, a))
                {
                    repeats.RemoveAt(j);
                    continue;
                }
                j++;
            }

            if (!removedI)
                i++;
        }

        return repeats;
    }

    private static bool IsTrfRedundant(ApproximateTandemRepeatResult x, ApproximateTandemRepeatResult y) =>
        (x.Period > y.Period && y.Period > 0 && x.Period % y.Period == 0 && x.AlignmentScore <= 1.1 * y.AlignmentScore) ||
        (x.Period == y.Period && x.AlignmentScore <= y.AlignmentScore);

    #endregion

    #region TRF Bernoulli statistics (Benson 1999)

    // Benson (1999) / TRF README "Probabilistic Model of Tandem Repeats": "We model alignment of two tandem
    // copies of a pattern of length n by a sequence of n independent Bernoulli trials ... P(Heads), which we
    // also call PM or matching probability, represents the average percent identity between the copies ...
    // PI or indel probability specifies the average percentage of insertions and deletions between the
    // copies." Default "(PM = .80, PI = .10)"; "Probabilistic data is available for PM values of 80 and 75 and
    // PI values of 10 and 20."

    /// <summary>Benson (1999) default Bernoulli matching probability PM = 0.80.</summary>
    public const double TrfDefaultMatchProbability = 0.80;

    /// <summary>Benson (1999) default Bernoulli indel probability PI = 0.10.</summary>
    public const double TrfDefaultIndelProbability = 0.10;

    /// <summary>
    /// Estimates the TRF Bernoulli-model parameters of a tandem-repeat tract: PM (match probability) and PI
    /// (indel probability) between ADJACENT copies. The tract is analysed exactly as TRF analyses a detected
    /// repeat: wraparound-DP alignment against tandem copies of the candidate pattern (the last
    /// <paramref name="period"/> bases of the tract), majority-rule consensus, realignment against the
    /// consensus, then comparison of each copy with the next one through that alignment (TRF: statistics refer
    /// to "the matches, mismatches and indels overall between adjacent copies in the sequence, not between the
    /// sequence and the consensus pattern"). Each compared column is one Bernoulli trial: heads = match,
    /// tails = mismatch or indel. PM and PI therefore equal TRF's reported % matches / % indels (/100) for
    /// the same region.
    /// </summary>
    /// <remarks>
    /// The statistics cover the locally aligned part of the tract (flanks that do not align are ignored, as in
    /// TRF). A tract without two aligned copies yields zero trials and PM = PI = 0.
    /// <see cref="TandemRepeatBernoulliStatistics.ExpectedMatches"/> is PM × trials (= matches).
    /// </remarks>
    /// <param name="repeatTract">The tandem-repeat tract (≥ 2 × period symbols; case-insensitive).</param>
    /// <param name="period">Candidate period of the tract (1..2000).</param>
    /// <param name="expectedMatchProbability">PM the tract is compared with (default 0.80, Benson 1999).</param>
    public static TandemRepeatBernoulliStatistics ComputeBernoulliStatistics(
        string repeatTract,
        int period,
        double expectedMatchProbability = TrfDefaultMatchProbability)
    {
        ArgumentNullException.ThrowIfNull(repeatTract);
        ArgumentOutOfRangeException.ThrowIfLessThan(period, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(period, MaxApproximatePeriod);
        if (!(expectedMatchProbability >= 0.0 && expectedMatchProbability <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(expectedMatchProbability));
        if (repeatTract.Length < period * 2)
            throw new ArgumentException(
                "A tandem repeat needs at least two contiguous copies of the period.", nameof(repeatTract));

        int n = repeatTract.Length;
        var s = new char[n + 1];
        repeatTract.ToUpperInvariant().CopyTo(0, s, 1, n);
        var pattern = new char[period];
        Array.Copy(s, n - period + 1, pattern, 0, period);

        TrfCopyComparison copies = default;
        double copyNumber = 0.0;
        var initial = TrfWraparoundAlign(TrfModel.Recommended, s, n, n, pattern);
        if (initial is not null)
        {
            char[] consensus = TrfConsensus(initial.Columns, period);
            var final = consensus.Length > 0 ? TrfWraparoundAlign(TrfModel.Recommended, s, n, n, consensus) : null;
            if (final is not null)
            {
                copies = CompareAdjacentCopies(final.Columns);
                copyNumber = final.CopyNumber;
            }
        }

        int trials = copies.Matches + copies.Mismatches + copies.Indels;
        double pm = trials > 0 ? (double)copies.Matches / trials : 0.0;
        double pi = trials > 0 ? (double)copies.Indels / trials : 0.0;

        return new TandemRepeatBernoulliStatistics(
            Period: period,
            AdjacentCopyPairs: Math.Max(0, (int)Math.Ceiling(copyNumber) - 1),
            BernoulliTrials: trials,
            Matches: copies.Matches,
            Mismatches: copies.Mismatches,
            Indels: copies.Indels,
            MatchProbability: pm,
            IndelProbability: pi,
            PercentMatches: pm * 100.0,
            PercentIndels: pi * 100.0,
            ExpectedMatches: pm * trials,
            MeetsExpectedMatchProbability: trials > 0 && pm >= expectedMatchProbability);
    }

    #endregion

    #region Inverted Repeat Detection

    /// <summary>
    /// Finds inverted repeats: a left arm followed, after a loop, by a right arm equal to the reverse complement of
    /// the left arm — exact (perfect) stems by default, optionally with up to <c>maxMismatches</c> mismatched pairs,
    /// a maximum arm length and G·U wobble pairs. Such structures can form hairpin/stem-loop (single strand) or
    /// cruciform (duplex) structures.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Model (default: exact, maximal stems).</b> A stem is a triple (<c>LeftArmStart = i</c>,
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
    /// <b>Pairing.</b> By default only the unambiguous bases A, C, G, T pair, via the canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/> (A↔T, C↔G). As in EMBOSS <c>einverted</c>
    /// (which scores a match only for a/c/g/t), N and other IUPAC ambiguity codes never form a pair, so a run
    /// of N is not reported as a stem.
    /// </para>
    /// <para>
    /// <b>Options.</b> <c>maxMismatches &gt; 0</c> and a finite <c>maxArmLength</c> follow EMBOSS <c>palindrome</c>
    /// <c>-nummismatches</c> / <c>-maxpallen</c> exactly (see the parameters; with <c>minLoopLength = 0</c> the result
    /// set equals the EMBOSS 6.6.0 binary, cross-checked on thousands of random cases); <c>allowWobble</c> adds G·U
    /// pairs under the same maximal-stem rule. No gaps (bulges): for scored, gap-tolerant inverted repeats use
    /// <see cref="FindInvertedRepeatsScored(string, int, int, int, int, int)"/> (EMBOSS <c>einverted</c>).
    /// </para>
    /// <para>Coordinates are 0-based; results are ordered by <c>LeftArmStart</c>, then <c>RightArmStart</c>.
    /// Parameters are validated eagerly.</para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minArmLength">Minimum length of each arm (default: 4, must be ≥ 2).</param>
    /// <param name="maxLoopLength">Maximum loop length between arms (default: 50, must be ≥ <paramref name="minLoopLength"/>).</param>
    /// <param name="minLoopLength">Minimum loop length (default: 3, must be ≥ 0).</param>
    /// <param name="maxMismatches">
    /// Maximum number of mismatched pairs inside a stem (default 0 = exact stems; must be ≥ 0). EMBOSS
    /// <c>palindrome -nummismatches</c> semantics: every stem starts at a pairing outer pair and is extended inward
    /// pair by pair until the (<paramref name="maxMismatches"/>+1)-th mismatch (or until the loop would become shorter
    /// than <paramref name="minLoopLength"/>); mismatches at the inner end of the stem are trimmed, interior mismatches
    /// are kept and counted in <see cref="InvertedRepeatResult.Mismatches"/>; a stem lying inside another stem in both
    /// arms is dropped (<c>-overlap Y</c>).
    /// </param>
    /// <param name="maxArmLength">
    /// Maximum arm length (default <see cref="int.MaxValue"/> = unbounded; must be ≥ <paramref name="minArmLength"/>).
    /// EMBOSS <c>palindrome -maxpallen</c> semantics: stems are only started from outer pairs spanning at most
    /// <c>2·maxArmLength + maxLoopLength + 1</c> bases, and stems longer than <paramref name="maxArmLength"/> are not
    /// reported but still suppress the stems they contain (so a longer stem is <b>not</b> split into shorter pieces;
    /// only parts of it that do not lie inside a longer candidate can surface).
    /// </param>
    /// <param name="allowWobble">
    /// When true, G·T (G·U) wobble pairs count as pairs in addition to Watson–Crick pairs, using the canonical
    /// <see cref="RnaSecondaryStructure.CanPair(char, char)"/> (A·U, G·C, G·U — T is read as U; Crick 1966,
    /// Varani &amp; McClain 2000). Default false (Watson–Crick A·T/C·G only, U never pairs).
    /// </param>
    /// <returns>The maximal (non-nested) stems, ordered by <c>LeftArmStart</c>, then <c>RightArmStart</c>.</returns>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        DnaSequence sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3,
        int maxMismatches = 0,
        int maxArmLength = int.MaxValue,
        bool allowWobble = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength);

        return FindInvertedRepeatsDispatch(
            sequence.Sequence, minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength, allowWobble);
    }

    /// <summary>
    /// Finds maximal inverted repeats in a raw sequence string (case-insensitive).
    /// Same model, options and validation as
    /// <see cref="FindInvertedRepeats(DnaSequence, int, int, int, int, int, bool)"/>;
    /// <c>null</c> or empty input yields no results. Characters other than A/C/G/T never pair
    /// (with <paramref name="allowWobble"/>, U pairs as T).
    /// </summary>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        string sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3,
        int maxMismatches = 0,
        int maxArmLength = int.MaxValue,
        bool allowWobble = false)
    {
        ValidateInvertedRepeatParameters(minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<InvertedRepeatResult>();

        return FindInvertedRepeatsDispatch(
            sequence.ToUpperInvariant(), minArmLength, maxLoopLength, minLoopLength, maxMismatches, maxArmLength, allowWobble);
    }

    private static void ValidateInvertedRepeatParameters(
        int minArmLength, int maxLoopLength, int minLoopLength, int maxMismatches, int maxArmLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minArmLength, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLoopLength, minLoopLength);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxArmLength, minArmLength);
    }

    private static List<InvertedRepeatResult> FindInvertedRepeatsDispatch(
        string seq, int minArmLength, int maxLoopLength, int minLoopLength,
        int maxMismatches, int maxArmLength, bool allowWobble)
    {
        // An arm can never exceed n/2, so maxArmLength ≥ n/2 removes neither a start pair nor a stem.
        bool armBounded = maxArmLength < seq.Length / 2;
        if (maxMismatches == 0 && !armBounded)
            return FindInvertedRepeatsCore(seq, minArmLength, maxLoopLength, minLoopLength, allowWobble);

        int armLimit = armBounded ? maxArmLength : int.MaxValue;
        return FindInvertedRepeatsMismatchCore(
            seq, minArmLength, maxLoopLength, minLoopLength, maxMismatches, armLimit, allowWobble);
    }

    /// <summary>True when <paramref name="a"/> and <paramref name="b"/> form a Watson–Crick pair (ACGT only).</summary>
    private static bool IsWatsonCrickPair(char a, char b) =>
        AcgtCode(a) >= 0 && SequenceExtensions.GetComplementBase(a) == b;

    /// <summary>Stem pairing rule: Watson–Crick (ACGT), or the canonical RNA pair set incl. G·U when <paramref name="wobble"/>.</summary>
    private static bool IsStemPair(char a, char b, bool wobble) =>
        wobble ? RnaSecondaryStructure.CanPair(a, b) : IsWatsonCrickPair(a, b);

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
        int minLoopLength,
        bool wobble)
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
                if (!IsStemPair(seq[iIn], seq[rIn], wobble))
                    continue;
                // Inward-maximal: the next inner pair must fail, unless it would make the loop < minLoop.
                if (loop - 2 >= minLoopLength && IsStemPair(seq[iIn + 1], seq[rIn - 1], wobble))
                    continue;

                int arm = 1;
                while (iIn - arm >= 0 && rIn + arm < n && IsStemPair(seq[iIn - arm], seq[rIn + arm], wobble))
                    arm++;
                if (arm < minArmLength)
                    continue;

                int iOut = iIn - arm + 1;
                if (IsContainedInShiftedStem(seq, iOut, rIn, arm, loop - minLoopLength, wobble))
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
    private static bool IsContainedInShiftedStem(string seq, int left, int right, int arm, int maxShift, bool wobble)
    {
        int n = seq.Length;
        for (int m = 1; m <= maxShift; m++)
        {
            int len = arm + m;
            if (right + len <= n && IsExactStem(seq, left, right, len, wobble))
                return true;
            if (left - m >= 0 && IsExactStem(seq, left - m, right - m, len, wobble))
                return true;
        }
        return false;
    }

    /// <summary>True when s[i+k] pairs with s[j+len−1−k] for all k (checked from the innermost pair out).</summary>
    private static bool IsExactStem(string seq, int i, int j, int len, bool wobble)
    {
        for (int k = len - 1; k >= 0; k--)
        {
            if (!IsStemPair(seq[i + k], seq[j + len - 1 - k], wobble))
                return false;
        }
        return true;
    }

    /// <summary>
    /// Mismatch-tolerant and/or arm-bounded stems: the candidate set of EMBOSS <c>palindrome</c>
    /// (<c>-nummismatches k</c>, <c>-maxpallen</c>) followed by the <c>palindrome_AInB</c> nesting filter.
    /// All pairs of a stem lie on one anti-diagonal <c>D = left + right</c>; along a diagonal the pairs are indexed
    /// from the innermost admissible pair (loop ≥ <paramref name="minLoop"/>, index 0) outward. A candidate starts at a
    /// pairing outer index <c>u</c> and walks inward to just before the (k+1)-th mismatch (or to index 0); trailing
    /// mismatches are trimmed, so its inner end is the first pairing index above that stop. Kept when
    /// arm ≥ <paramref name="minArm"/>, loop ≤ <paramref name="maxLoopLength"/>, and (bounded arm) the outer span is at
    /// most <c>2·maxArm + maxLoopLength + 1</c>; a candidate lying inside another candidate in both arms is dropped
    /// (any container differs by at most <c>loop − minLoop</c> diagonals); candidates longer than
    /// <paramref name="maxArm"/> are then not reported (<c>palindrome_Print</c> filter).
    /// </summary>
    private static List<InvertedRepeatResult> FindInvertedRepeatsMismatchCore(
        string seq, int minArm, int maxLoopLength, int minLoop, int maxMismatches, int maxArm, bool wobble)
    {
        int n = seq.Length;
        var results = new List<InvertedRepeatResult>();
        if (n < 2 * minArm + minLoop)
            return results;

        long maxLoop = Math.Min(maxLoopLength, (long)n);
        long maxSpan = maxArm == int.MaxValue ? long.MaxValue : 2L * maxArm + maxLoopLength + 1;
        int diagonals = 2 * n - 1;

        // Candidates grouped by diagonal; within a diagonal stored with increasing outer-left position.
        var outerLeft = new List<int>();
        var innerLeft = new List<int>();
        var mismatchCount = new List<int>();
        var diagStart = new int[diagonals + 1];
        var mismatchAt = new List<int>();
        var nextPairAfter = new List<int>();
        var mismatchesBelow = new List<int>(); // per index u: mismatches at indices < u
        var diagOuter = new List<int>();
        var diagInner = new List<int>();
        var diagMism = new List<int>();

        for (int d = 0; d < diagonals; d++)
        {
            diagStart[d] = outerLeft.Count;
            long a0Long = (long)d - 1 - minLoop;
            if (a0Long < 0)
                continue;
            int a0 = (int)(a0Long / 2);
            int loop0 = d - 2 * a0 - 1;
            if (loop0 > maxLoop)
                continue;
            int innerLimit = (int)((maxLoop - loop0) / 2);
            long uMaxLong = Math.Min(a0, (long)n - 1 - d + a0);
            if (maxSpan != long.MaxValue)
            {
                long spanRoom = maxSpan - loop0 - 2;
                uMaxLong = spanRoom < 0 ? -1 : Math.Min(uMaxLong, spanRoom / 2);
            }
            int uMax = (int)uMaxLong;

            mismatchAt.Clear();
            nextPairAfter.Clear();
            mismatchesBelow.Clear();
            diagOuter.Clear();
            diagInner.Clear();
            diagMism.Clear();
            int pending = 0;
            int firstPair = -1;
            for (int u = 0; u <= uMax; u++)
            {
                mismatchesBelow.Add(mismatchAt.Count);
                if (!IsStemPair(seq[a0 - u], seq[d - a0 + u], wobble))
                {
                    mismatchAt.Add(u);
                    nextPairAfter.Add(-1);
                    // Every later start stops at or above this (k+1)-th mismatch → its inner end exceeds innerLimit.
                    if (mismatchAt.Count > maxMismatches && mismatchAt[^(maxMismatches + 1)] >= innerLimit)
                        break;
                    continue;
                }

                if (firstPair < 0)
                    firstPair = u;
                for (; pending < mismatchAt.Count; pending++)
                    nextPairAfter[pending] = u;

                // Inner end: first pairing index above the (k+1)-th mismatch met walking inward (trailing mismatches trimmed).
                int inner = mismatchAt.Count > maxMismatches
                    ? nextPairAfter[mismatchAt.Count - 1 - maxMismatches]
                    : firstPair;
                int mism = mismatchAt.Count - mismatchesBelow[inner];

                if (u - inner + 1 >= minArm && inner <= innerLimit)
                {
                    diagOuter.Add(a0 - u);
                    diagInner.Add(a0 - inner);
                    diagMism.Add(mism);
                }
            }

            for (int c = diagOuter.Count - 1; c >= 0; c--)
            {
                outerLeft.Add(diagOuter[c]);
                innerLeft.Add(diagInner[c]);
                mismatchCount.Add(diagMism[c]);
            }
        }
        diagStart[diagonals] = outerLeft.Count;

        // Prefix maximum of the inner-left end per diagonal (outer-left ascending) for containment queries, and the
        // longest candidate arm per diagonal / overall (a container shifted by m diagonals is ≥ m bases longer).
        var prefixMaxInner = new int[outerLeft.Count];
        var longestOnDiagonal = new int[diagonals];
        int longestArm = 0;
        for (int d = 0; d < diagonals; d++)
        {
            int best = int.MinValue;
            for (int c = diagStart[d]; c < diagStart[d + 1]; c++)
            {
                prefixMaxInner[c] = best = Math.Max(best, innerLeft[c]);
                longestOnDiagonal[d] = Math.Max(longestOnDiagonal[d], innerLeft[c] - outerLeft[c] + 1);
            }
            longestArm = Math.Max(longestArm, longestOnDiagonal[d]);
        }

        // True when diagonal dd holds a candidate with outer-left ≤ x and inner-left ≥ y.
        bool Covered(int dd, int x, int y)
        {
            if (longestOnDiagonal[dd] < y - x + 1)
                return false;
            int lo = diagStart[dd], hi = diagStart[dd + 1] - 1, found = -1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >>> 1;
                if (outerLeft[mid] <= x) { found = mid; lo = mid + 1; }
                else hi = mid - 1;
            }
            return found >= 0 && prefixMaxInner[found] >= y;
        }

        for (int d = 0; d < diagonals; d++)
        {
            for (int c = diagStart[d]; c < diagStart[d + 1]; c++)
            {
                int fs = outerLeft[c], fe = innerLeft[c];
                int arm = fe - fs + 1;
                if (arm > maxArm)
                    continue;
                int loop = d - 2 * fe - 1;
                int rightOuter = d - fs, rightInner = d - fe;

                bool contained = Covered(d, fs - 1, fe);
                int maxShift = Math.Min(loop - minLoop, longestArm - arm);
                for (int m = 1; !contained && m <= maxShift; m++)
                {
                    int up = d + m, down = d - m;
                    if (up < diagonals && Covered(up, Math.Min(fs, up - rightOuter), Math.Max(fe, up - rightInner)))
                        contained = true;
                    else if (down >= 0 && Covered(down, Math.Min(fs, down - rightOuter), Math.Max(fe, down - rightInner)))
                        contained = true;
                }
                if (contained)
                    continue;

                results.Add(new InvertedRepeatResult(
                    LeftArmStart: fs,
                    RightArmStart: rightInner,
                    ArmLength: arm,
                    LoopLength: loop,
                    LeftArm: seq.Substring(fs, arm),
                    RightArm: seq.Substring(rightInner, arm),
                    Loop: seq.Substring(fe + 1, loop),
                    CanFormHairpin: loop >= 3)
                { Mismatches = mismatchCount[c] });
            }
        }

        results.Sort(static (x, y) =>
        {
            int c = x.LeftArmStart.CompareTo(y.LeftArmStart);
            return c != 0 ? c : x.RightArmStart.CompareTo(y.RightArmStart);
        });
        return results;
    }

    #endregion

    #region Scored (einverted) Inverted Repeat Detection

    /// <summary>EMBOSS <c>einverted</c> default gap penalty (<c>-gap 12</c>).</summary>
    public const int EinvertedDefaultGapPenalty = 12;

    /// <summary>EMBOSS <c>einverted</c> default minimum score (<c>-threshold 50</c>).</summary>
    public const int EinvertedDefaultThreshold = 50;

    /// <summary>EMBOSS <c>einverted</c> default match score (<c>-match 3</c>).</summary>
    public const int EinvertedDefaultMatchScore = 3;

    /// <summary>EMBOSS <c>einverted</c> default mismatch score (<c>-mismatch -4</c>).</summary>
    public const int EinvertedDefaultMismatchScore = -4;

    /// <summary>
    /// EMBOSS 6.6 <c>einverted</c> default maximum repeat extent (<c>-maxrepeat 2000</c>, einverted.acd). Durbin's
    /// original compile-time value was 4000 (comment in <c>einverted.c</c>); pass it explicitly to reproduce that.
    /// </summary>
    public const int EinvertedDefaultMaxRepeatLength = 2000;

    /// <summary>Sentinel score of <c>einverted</c> (<c>rogue</c>) marking cells outside the sequence/window.</summary>
    private const int EinvertedRogue = 1_000_000;

    /// <summary>
    /// Finds imperfect (mismatch- and gap-tolerant) inverted repeats by score, reproducing EMBOSS <c>einverted</c>
    /// (Durbin &amp; Thierry-Mieg 1993, "Inverted repeats by dynamic programming").
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Scoring.</b> A local alignment of the sequence against its own reverse complement, grown outward from the
    /// loop: cell (<c>i</c>, <c>k</c>) pairs right-arm base <c>i</c> with left-arm base <c>i − 1 − k</c> and scores
    /// <c>H(i,k) = max(s(i, i−1−k) + max(0, H(i−1, k−2)), max(H(i, k−1), H(i−1, k−1)) − gap)</c>, where <c>s</c> is
    /// <paramref name="matchScore"/> for a Watson–Crick pair of A/C/G/T and <paramref name="mismatchScore"/>
    /// otherwise (N and other symbols never pair). Only cells with <c>k &lt; maxRepeatLength − 2</c> (outer span at
    /// most <paramref name="maxRepeatLength"/>) are scored.
    /// </para>
    /// <para>
    /// <b>Reporting.</b> Exactly the <c>einverted</c> procedure: for each right end <c>i</c> the first best cell with
    /// score ≥ <paramref name="threshold"/> is remembered against its left start (a later end with the same left start
    /// replaces it only with a higher score; left starts at or before the last reported right end are ignored); once
    /// the scan is <paramref name="maxRepeatLength"/> past a remembered left start, the best-scoring end in the window
    /// is traced back (same-row gap, previous-row gap, then diagonal) and reported, and the window is cleared.
    /// Results are returned in einverted's report order (not necessarily sorted by position). Degenerate parameters
    /// (threshold ≤ match score, zero gap penalty) can make einverted 6.6.0 abort (SIGFPE: its percentage divides by
    /// matches + mismatches): a trace-back that records no column is consumed without a result, and a repeat made
    /// only of gap columns is returned with <see cref="ScoredInvertedRepeatResult.PercentMatches"/> = 0 (the scan
    /// otherwise continues exactly as einverted's code).
    /// </para>
    /// <para>
    /// Cross-checked against the EMBOSS 6.6.0 <c>einverted</c> binary (coordinates, score, matches, mismatches,
    /// gaps and both alignment rows) on random sequences with planted imperfect inverted repeats. Output differences
    /// by design: coordinates are 0-based and inclusive; alignment rows use the upper-cased input characters (einverted
    /// prints lower case and shows non-ACGT symbols as '-').
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence (case-insensitive; symbols other than A/C/G/T never pair).</param>
    /// <param name="gapPenalty">Gap penalty (≥ 0, default 12).</param>
    /// <param name="threshold">Minimum reported score (≥ 0, default 50).</param>
    /// <param name="matchScore">Score of a Watson–Crick pair (≥ 0, default 3).</param>
    /// <param name="mismatchScore">Score of any other pair (≤ 0, default −4).</param>
    /// <param name="maxRepeatLength">
    /// Maximum extent from the start of the repeat to the end of its inverted copy (einverted <c>-maxrepeat</c>,
    /// ≥ 2, default 2000). Memory is O(min(maxRepeatLength, n)²).
    /// </param>
    /// <returns>Scored inverted repeats in einverted report order.</returns>
    public static IEnumerable<ScoredInvertedRepeatResult> FindInvertedRepeatsScored(
        string sequence,
        int gapPenalty = EinvertedDefaultGapPenalty,
        int threshold = EinvertedDefaultThreshold,
        int matchScore = EinvertedDefaultMatchScore,
        int mismatchScore = EinvertedDefaultMismatchScore,
        int maxRepeatLength = EinvertedDefaultMaxRepeatLength)
    {
        ValidateScoredInvertedRepeatParameters(gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ScoredInvertedRepeatResult>();
        return new EinvertedScanner(sequence.ToUpperInvariant(), gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength).Run();
    }

    /// <summary>
    /// Scored inverted repeats of a <see cref="DnaSequence"/>; see
    /// <see cref="FindInvertedRepeatsScored(string, int, int, int, int, int)"/>.
    /// </summary>
    public static IEnumerable<ScoredInvertedRepeatResult> FindInvertedRepeatsScored(
        DnaSequence sequence,
        int gapPenalty = EinvertedDefaultGapPenalty,
        int threshold = EinvertedDefaultThreshold,
        int matchScore = EinvertedDefaultMatchScore,
        int mismatchScore = EinvertedDefaultMismatchScore,
        int maxRepeatLength = EinvertedDefaultMaxRepeatLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindInvertedRepeatsScored(sequence.Sequence, gapPenalty, threshold, matchScore, mismatchScore, maxRepeatLength);
    }

    private static void ValidateScoredInvertedRepeatParameters(
        int gapPenalty, int threshold, int matchScore, int mismatchScore, int maxRepeatLength)
    {
        // Ranges of einverted.acd; scores must stay below the rogue sentinel (einverted's arrays assume it).
        ArgumentOutOfRangeException.ThrowIfNegative(gapPenalty);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(gapPenalty, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfNegative(threshold);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(threshold, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfNegative(matchScore);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(matchScore, EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(mismatchScore, 0);
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(mismatchScore, -EinvertedRogue);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxRepeatLength, 2);
    }

    /// <summary>
    /// The einverted scan: a ring of <c>W</c> DP rows (row <c>i</c> in slot <c>i mod W</c>), per-row best scores and a
    /// left-start → right-end table, both ring-indexed. With <c>W ≥ n + 2</c> no slot is ever reused and every row is
    /// scored in full, so <c>W</c> is capped at <c>n + 2</c> without changing the output.
    /// </summary>
    private sealed class EinvertedScanner
    {
        private readonly string _seq;
        private readonly int[] _code;
        private readonly int _n, _w, _gap, _threshold, _match, _mismatch;
        private readonly int[][] _rows;
        private readonly int[] _bestEnd;   // einverted "back": slot of a left start → right end (row) of its best repeat
        private readonly int[] _rowBest;   // einverted "localMax": best score of a row (0 = none)
        private readonly List<ScoredInvertedRepeatResult> _results = new();

        public EinvertedScanner(string seq, int gap, int threshold, int match, int mismatch, int maxRepeat)
        {
            _seq = seq;
            _n = seq.Length;
            _w = (int)Math.Min(maxRepeat, (long)_n + 2);
            _gap = gap;
            _threshold = threshold;
            _match = match;
            _mismatch = mismatch;
            _code = new int[_n];
            for (int p = 0; p < _n; p++)
            {
                int c = AcgtCode(seq[p]);
                _code[p] = c < 0 ? 4 : c;
            }
            _rows = new int[_w][];
            for (int r = 0; r < _w; r++)
                _rows[r] = new int[_w];
            _bestEnd = new int[_w];
            _rowBest = new int[_w];
        }

        /// <summary>Pair score of right base <paramref name="i"/> with left base <c>i − 1 − k</c>; rogue off the sequence.</summary>
        private int PairScore(int i, int k)
        {
            int partner = i - 1 - k;
            if (partner < 0)
                return EinvertedRogue;
            int ci = _code[i];
            return ci < 4 && _code[partner] == 3 - ci ? _match : _mismatch;
        }

        public List<ScoredInvertedRepeatResult> Run()
        {
            int w = _w;
            int lastReported = -1;
            for (int i = 0; i < _n + w; i++)
            {
                int slot = i % w;
                if (_bestEnd[slot] != 0)
                    lastReported = ReportWindow(i, slot);

                if (i >= _n)
                    continue;

                int[] prev = _rows[(i + w - 1) % w];
                int[] cur = _rows[slot];
                for (int k = 0; k < w - 2; k++)
                    cur[k] = PairScore(i, k);
                cur[w - 2] = cur[w - 1] = EinvertedRogue;
                if (i == 0)
                    Array.Clear(prev); // row "−1" of einverted is the zero-initialised matrix

                int best = _threshold - 1, bestK = 0;
                int c = cur[0], diag = -EinvertedRogue;
                // Always ends at the rogue column (k ≤ W − 2 or the first column past the sequence start).
                for (int k = 1; k < w; k++)
                {
                    int d = cur[k];
                    if (diag > 0)
                        d += diag;
                    diag = prev[k - 1];
                    if (diag > c)
                        c = diag;
                    c -= _gap;
                    if (d > c)
                        c = d;
                    cur[k] = c;
                    if (c >= EinvertedRogue)
                        break;
                    if (c > best)
                    {
                        best = c;
                        bestK = k;
                    }
                }

                if (bestK != 0)
                {
                    int leftStart = i - bestK - 1;
                    int startSlot = leftStart % w;
                    int previousEnd = _bestEnd[startSlot];
                    if (leftStart > lastReported && (previousEnd == 0 || _rowBest[previousEnd % w] < best))
                    {
                        _bestEnd[startSlot] = i;
                        _rowBest[slot] = best;
                    }
                }
                else
                {
                    _rowBest[slot] = 0;
                }
            }
            return _results;
        }

        /// <summary>Reports the best repeat of the window (i − W, end] and clears it; returns the new last-reported end.</summary>
        private int ReportWindow(int i, int slot)
        {
            int w = _w;
            int windowEnd = _bestEnd[slot];
            int bestRow = windowEnd;
            bool done = false;
            while (!done)
            {
                int bestScore = 0;
                bestRow = windowEnd;
                for (int j = windowEnd; j > i - w; j--)
                {
                    int score = _rowBest[j % w];
                    if (score > bestScore)
                    {
                        bestRow = j;
                        bestScore = score;
                    }
                }
                if (bestScore == 0)
                    break;
                done = TraceBack(bestScore, bestRow, i - w);
                if (!done)
                    _rowBest[bestRow % w] = 0;
            }

            for (int j = bestRow; j >= i - w; j--)
            {
                _bestEnd[j % w] = 0;
                _rowBest[j % w] = 0;
            }
            return bestRow;
        }

        private bool TraceBack(int score, int endRow, int minRow)
        {
            int w = _w;
            int[] row = _rows[endRow % w];
            int k = 0;
            while (k < w && row[k] != score)
                k++;

            var right = new List<int>();
            var left = new List<int>();
            int remaining = score, i = endRow, matches = 0, mismatches = 0, gaps = 0;
            while (remaining > 0 && k >= 1)
            {
                if (i < minRow)
                    return false;
                right.Add(i);
                left.Add(i - 1 - k);
                if (row[k - 1] == remaining + _gap)
                {
                    remaining += _gap;
                    gaps++;
                    k--;
                    continue;
                }
                row = _rows[(i - 1 + w) % w];
                if (row[k - 1] == remaining + _gap)
                {
                    remaining += _gap;
                    gaps++;
                    i--;
                    k--;
                    continue;
                }
                int s = PairScore(i, k);
                remaining -= s;
                if (s == _match)
                    matches++;
                else
                    mismatches++;
                i--;
                k -= 2;
            }

            // einverted divides by (matches + mismatches) and prints coordinates here: with no aligned column, or a
            // column before the sequence start (possible only after a ring slot was reused), it aborts or prints
            // garbage. Such a trace-back still counts as reported (scan state as in einverted) but yields no result.
            if (left.Count == 0 || left[0] < 0)
                return true;

            int len = left.Count;
            var leftRow = new char[len];
            var midRow = new char[len];
            var rightRow = new char[len];
            for (int t = 0; t < len; t++)
            {
                bool leftGap = t + 1 < len && left[t] == left[t + 1];
                bool rightGap = t + 1 < len && right[t] == right[t + 1];
                leftRow[t] = leftGap ? '-' : _seq[left[t]];
                rightRow[t] = rightGap ? '-' : _seq[right[t]];
                midRow[t] = !leftGap && !rightGap && _code[right[t]] + _code[left[t]] == 3 ? '|' : ' ';
            }

            _results.Add(new ScoredInvertedRepeatResult(
                LeftArmStart: left[0],
                LeftArmEnd: left[len - 1],
                RightArmStart: right[len - 1],
                RightArmEnd: right[0],
                Score: score,
                Matches: matches,
                Mismatches: mismatches,
                Gaps: gaps,
                LeftArmAlignment: new string(leftRow),
                MatchLine: new string(midRow),
                RightArmAlignment: new string(rightRow)));
            return true;
        }
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

    /// <summary>A/C/G/T → 0..3, other symbols → −1 (shared <see cref="SequenceComplexity.AcgtCode"/>; inputs here are upper-cased).</summary>
    private static int AcgtCode(char c) => SequenceComplexity.AcgtCode(c);

    /// <summary>A maximal pair of suffix starts (P, Q) in the encoded text with common prefix length exactly <c>Length</c>.</summary>
    private readonly record struct RawMaximalPair(int P, int Q, int Length);

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

        foreach (var pair in EnumerateForwardMaximalPairs(seq, minLength, maxLength))
        {
            int i = Math.Min(pair.P, pair.Q);
            int j = Math.Max(pair.P, pair.Q);
            long spacing = (long)j - i - pair.Length;
            if (spacing < minSpacing) continue;
            results.Add(new DirectRepeatResult(
                FirstPosition: i,
                SecondPosition: j,
                RepeatSequence: seq.Substring(i, pair.Length),
                Length: pair.Length,
                Spacing: (int)spacing));
        }

        results.Sort(static (a, b) =>
        {
            int c = a.FirstPosition.CompareTo(b.FirstPosition);
            return c != 0 ? c : a.SecondPosition.CompareTo(b.SecondPosition);
        });
        return results;
    }

    /// <summary>
    /// Every maximal exact repeated pair (Gusfield 1997 §7.12) of the (upper-cased) sequence with
    /// length in [minLength, maxLength]: A/C/G/T → 0..3, every other symbol a unique code (never matches).
    /// </summary>
    private static List<RawMaximalPair> EnumerateForwardMaximalPairs(string seq, int minLength, int maxLength)
    {
        int n = seq.Length;
        var symbols = new int[n];
        var leftClass = new int[n];
        for (int p = 0; p < n; p++)
        {
            int code = AcgtCode(seq[p]);
            symbols[p] = code >= 0 ? code : 4 + p;
            leftClass[p] = p > 0 && symbols[p - 1] < 4 ? symbols[p - 1] : UniqueLeftClass;
        }

        return EnumerateMaximalPairs(symbols, leftClass, strandOf: null, minLength, maxLength);
    }

    /// <summary>
    /// Bottom-up lcp-interval traversal of the suffix array (Gusfield 1997 §7.12.3; Abouelhoda, Kurtz &amp;
    /// Ohlebusch 2004). Each interval keeps its suffix positions in one linked list per class
    /// (left character 0..4, times the strand when <paramref name="strandOf"/> is given); merging a child
    /// into its parent at string depth ℓ emits every pair (p from the child, q already in the parent)
    /// whose left characters differ (or are undefined) — those pairs have LCP exactly ℓ (right-maximal)
    /// and are left-maximal. With <paramref name="strandOf"/>, only pairs on different strands are emitted.
    /// O(n log² n + z) for z visited pairs.
    /// </summary>
    private static List<RawMaximalPair> EnumerateMaximalPairs(
        int[] symbols,
        int[] leftClass,
        int[]? strandOf,
        int minLength,
        int maxLength)
    {
        var pairs = new List<RawMaximalPair>();
        int n = symbols.Length;
        if (n < 2)
            return pairs;

        int strands = strandOf is null ? 1 : 2;
        int classCount = LeftClassCount * strands;
        var classOf = new int[n];
        for (int p = 0; p < n; p++)
            classOf[p] = (strandOf is null ? 0 : strandOf[p] * LeftClassCount) + leftClass[p];

        int[] sa = SequenceComplexity.BuildSuffixArray(symbols);
        int[] lcp = SequenceComplexity.BuildLcpArray(symbols, sa);

        var next = new int[n];
        var stackLcp = new int[n + 1];
        var stackHead = new int[(n + 1) * classCount];
        var stackTail = new int[(n + 1) * classCount];
        var childHead = new int[classCount];
        var childTail = new int[classCount];
        int top = 0;
        stackLcp[0] = 0;
        Array.Fill(stackHead, -1, 0, classCount);

        for (int r = 1; r <= n; r++)
        {
            // Pending child: the leaf for suffix sa[r − 1].
            Array.Fill(childHead, -1);
            int leaf = sa[r - 1];
            next[leaf] = -1;
            childHead[classOf[leaf]] = leaf;
            childTail[classOf[leaf]] = leaf;

            int h = r < n ? lcp[r] : 0;
            while (stackLcp[top] > h)
            {
                MergeMaximalPairLists(top, classCount, strands, childHead, childTail, stackLcp, stackHead,
                    stackTail, next, minLength, maxLength, pairs);
                for (int c = 0; c < classCount; c++)
                {
                    childHead[c] = stackHead[top * classCount + c];
                    childTail[c] = stackTail[top * classCount + c];
                }
                top--;
            }

            if (stackLcp[top] < h)
            {
                top++;
                stackLcp[top] = h;
                Array.Fill(stackHead, -1, top * classCount, classCount);
            }

            MergeMaximalPairLists(top, classCount, strands, childHead, childTail, stackLcp, stackHead,
                stackTail, next, minLength, maxLength, pairs);
        }

        return pairs;
    }

    /// <summary>
    /// Emits the maximal pairs between a child interval's lists and the lists already accumulated in the
    /// interval at stack slot <paramref name="node"/>, then concatenates the child lists into the node.
    /// </summary>
    private static void MergeMaximalPairLists(
        int node,
        int classCount,
        int strands,
        int[] childHead,
        int[] childTail,
        int[] stackLcp,
        int[] stackHead,
        int[] stackTail,
        int[] next,
        int minLength,
        int maxLength,
        List<RawMaximalPair> pairs)
    {
        int length = stackLcp[node];
        int baseIdx = node * classCount;

        if (length >= minLength && length <= maxLength)
        {
            for (int a = 0; a < classCount; a++)
            {
                if (childHead[a] < 0) continue;
                int leftA = a % LeftClassCount;
                for (int b = 0; b < classCount; b++)
                {
                    if (strands == 2 && a / LeftClassCount == b / LeftClassCount) continue; // same strand
                    if (leftA == b % LeftClassCount && leftA != UniqueLeftClass) continue; // not left-maximal
                    int nodeListHead = stackHead[baseIdx + b];
                    if (nodeListHead < 0) continue;

                    for (int p = childHead[a]; p >= 0; p = next[p])
                    {
                        for (int q = nodeListHead; q >= 0; q = next[q])
                            pairs.Add(new RawMaximalPair(p, q, length));
                    }
                }
            }
        }

        for (int c = 0; c < classCount; c++)
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

    #region Reverse-Complement (Palindromic) Maximal Repeats

    /// <summary>
    /// Finds exact <b>reverse-complement maximal repeated pairs</b>: the second copy is the reverse
    /// complement of the first (MUMmer <c>repeat-match</c> without <c>-f</c>, the <c>r</c> lines; Vmatch
    /// <c>-p</c> "palindromic matches"; REPuter palindromic repeats, Kurtz &amp; Schleiermacher 1999).
    /// </summary>
    /// <remarks>
    /// <para>
    /// A triple (<c>i</c>, <c>k</c>, <c>L</c>) with 0-based starts <c>i ≤ k</c> is reported when
    /// <c>S[i..i+L) = revcomp(S[k..k+L))</c> (Watson–Crick complement via
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/>) and the pair is maximal: it can be extended
    /// neither outward (<c>i = 0</c>, <c>k + L = n</c>, or <c>S[i−1]</c> does not pair with <c>S[k+L]</c>) nor
    /// inward (<c>i + L = n</c>, <c>k = 0</c>, or <c>S[i+L]</c> does not pair with <c>S[k−1]</c>). This is the
    /// Vmatch Appendix A definition of a maximal palindromic exact match of a sequence with itself
    /// (<c>i ≤ j</c>; <c>i = k</c> is allowed and denotes a reverse palindrome such as <c>GAATTC</c> or the
    /// centre of a longer one).
    /// </para>
    /// <para>
    /// <b>Coordinates.</b> <see cref="ReverseComplementRepeatResult.FirstPosition"/> = <c>i</c> and
    /// <see cref="ReverseComplementRepeatResult.SecondPosition"/> = <c>k</c> are 0-based <b>forward-strand
    /// starts</b> of both copies (Vmatch <c>-p</c> prints exactly these). MUMmer <c>repeat-match</c> prints the
    /// line <c>Start1 = i + 1</c>, <c>Start2 = k + L</c> followed by <c>r</c>: its Start2 is the 1-based
    /// position of the <i>last</i> base of the second copy, i.e. where the copy starts when read on the
    /// reverse strand (<c>Data[a+t] = Complement(Data[b−t])</c> in repeat-match.cc <c>Verify_Match</c>).
    /// Convert with <c>k = Start2 − L</c> (0-based).
    /// </para>
    /// <para>
    /// Only A/C/G/T pair (case-insensitive); N, IUPAC codes, U and gaps never match (MUMmer <c>-n</c> /
    /// Vmatch wildcard convention), as in <see cref="FindDirectRepeats(string,int,int,int)"/>.
    /// Filters: <c>minLength ≤ L ≤ maxLength</c> (a longer maximal pair is not truncated) and
    /// <c>Spacing = k − i − L ≥ minSpacing</c>: the number of bases between the end of the first copy and
    /// the start of the second (the loop of a hairpin; negative when the copies overlap, e.g.
    /// <c>−L</c> for <c>i = k</c>). <c>minSpacing = int.MinValue</c> with <c>maxLength = int.MaxValue</c>
    /// returns the complete <c>repeat-match</c> reverse-complement set.
    /// </para>
    /// <para>
    /// Algorithm: the suffix array + LCP of <c>S · # · revcomp(S)</c> (shared helpers of
    /// <see cref="SequenceComplexity"/>) traversed bottom-up with per-(strand, left-character) lists
    /// (Gusfield 1997 §7.12.3, restricted to pairs across the two strands); each pair is found in both
    /// orientations and kept once (<c>i ≤ k</c>, as repeat-match does). O(n log² n + z).
    /// Unlike direct pairs, one (i, k) can carry several maximal lengths on different anti-diagonals
    /// (e.g. T₇ facing A₇); results are ordered by (FirstPosition, SecondPosition, Length).
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5, must be ≥ 2).</param>
    /// <param name="maxLength">Maximum repeat length (default: 50, must be ≥ <paramref name="minLength"/>).</param>
    /// <param name="minSpacing">Minimum number of bases between the copies (default: 1; negative admits overlap).</param>
    /// <returns>Maximal reverse-complement repeat pairs, sorted by (FirstPosition, SecondPosition, Length).</returns>
    public static IEnumerable<ReverseComplementRepeatResult> FindReverseComplementRepeats(
        DnaSequence sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        return FindReverseComplementRepeatsCore(sequence.Sequence, minLength, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds maximal exact reverse-complement repeat pairs in a raw sequence string (case-insensitive;
    /// non-ACGT symbols never match). <c>null</c> or empty input yields no results. See
    /// <see cref="FindReverseComplementRepeats(DnaSequence,int,int,int)"/> for the definition and coordinates.
    /// </summary>
    public static IEnumerable<ReverseComplementRepeatResult> FindReverseComplementRepeats(
        string sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ReverseComplementRepeatResult>();

        return FindReverseComplementRepeatsCore(sequence.ToUpperInvariant(), minLength, maxLength, minSpacing);
    }

    private static List<ReverseComplementRepeatResult> FindReverseComplementRepeatsCore(
        string seq,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        var results = new List<ReverseComplementRepeatResult>();
        int n = seq.Length;
        if (n < minLength)
            return results;

        foreach (var pair in EnumerateReverseComplementSeeds(seq, minLength, maxLength))
        {
            int i = pair.P;
            int k = n - pair.Q - pair.Length;
            if (k < i) continue; // mirror image of the pair (k, i), which is also enumerated

            long spacing = (long)k - i - pair.Length;
            if (spacing < minSpacing) continue;
            results.Add(new ReverseComplementRepeatResult(
                FirstPosition: i,
                SecondPosition: k,
                RepeatSequence: seq.Substring(i, pair.Length),
                SecondSequence: seq.Substring(k, pair.Length),
                Length: pair.Length,
                Spacing: (int)spacing));
        }

        results.Sort(static (a, b) =>
        {
            int c = a.FirstPosition.CompareTo(b.FirstPosition);
            if (c != 0) return c;
            c = a.SecondPosition.CompareTo(b.SecondPosition);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        });
        return results;
    }

    /// <summary>
    /// Every maximal exact match between S and R = revcomp(S) with length in [minLength, maxLength], as
    /// (<c>P</c> = start in S, <c>Q</c> = start in R, Length): the suffix array + LCP of <c>S · # · R</c>
    /// traversed with per-(strand, left-character) lists, pairs across the two strands only. Each
    /// palindromic pair appears in both orientations (its mirror is (n − Q − L, n − P − L)). A/C/G/T match;
    /// every other symbol and the separator get a unique code.
    /// </summary>
    private static List<RawMaximalPair> EnumerateReverseComplementSeeds(string seq, int minLength, int maxLength)
    {
        int n = seq.Length;
        int total = 2 * n + 1;
        var symbols = new int[total];
        var strandOf = new int[total];
        for (int p = 0; p < n; p++)
        {
            int code = AcgtCode(seq[p]);
            symbols[p] = code >= 0 ? code : 4 + p;

            int t = n + 1 + (n - 1 - p);
            int rc = code >= 0 ? AcgtCode(SequenceExtensions.GetComplementBase(seq[p])) : -1;
            symbols[t] = rc >= 0 ? rc : 4 + t;
            strandOf[t] = 1;
        }
        symbols[n] = 4 + n;
        strandOf[n] = 0;

        var leftClass = new int[total];
        for (int p = 0; p < total; p++)
            leftClass[p] = p > 0 && symbols[p - 1] < 4 ? symbols[p - 1] : UniqueLeftClass;

        var pairs = EnumerateMaximalPairs(symbols, leftClass, strandOf, minLength, maxLength);
        for (int x = 0; x < pairs.Count; x++)
        {
            var pair = pairs[x];
            int forward = pair.P < n ? pair.P : pair.Q;
            int reverse = pair.P < n ? pair.Q : pair.P;
            pairs[x] = new RawMaximalPair(forward, reverse - (n + 1), pair.Length);
        }
        return pairs;
    }

    #endregion

    #region Degenerate (k-mismatch) Direct Repeats

    /// <summary>
    /// Finds <b>maximal k-mismatch (degenerate) direct repeats</b> under the Hamming distance (REPuter,
    /// Kurtz &amp; Schleiermacher 1999; Kurtz et al. 2001 NAR 29:4633 "k-mismatch repeats"; Vmatch <c>-h k</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Definition (Vmatch manual, Appendix A "Basic Notions", Kurtz): two equal-length substrings
    /// <c>S[i..i+L)</c> and <c>S[j..j+L)</c>, <c>i &lt; j</c>, form a <i>k-mismatch repeat</i> when their Hamming
    /// distance is ≤ <paramref name="maxMismatches"/>; the repeat is <i>maximal</i> when it is not contained in
    /// another k-mismatch repeat — here, as for exact maximal pairs (Gusfield 1997 §7.12), containment is taken
    /// on the same alignment diagonal <c>j − i</c>: the pair cannot be extended one position to the left or to
    /// the right (sequence boundary, or the extension would exceed <paramref name="maxMismatches"/>). With
    /// <c>maxMismatches = 0</c> the result is exactly <see cref="FindDirectRepeats(string,int,int,int)"/>.
    /// A maximal window may end on a mismatch (then the next position is also a mismatch or a boundary).
    /// </para>
    /// <para>
    /// Only A/C/G/T match (case-insensitive); N, IUPAC codes, U and gaps count as mismatches (the Vmatch
    /// wildcard rule: "a degenerate match may contain a wildcard, but this always leads to a pair of
    /// mismatching characters"). Copies may overlap; <c>Spacing = j − i − L</c> is filtered by
    /// <paramref name="minSpacing"/> after maximality, and <c>minLength ≤ L ≤ maxLength</c>.
    /// </para>
    /// <para>
    /// Algorithm (REPuter seed-and-extend, complete): by the pigeonhole principle every such repeat of length
    /// ≥ <c>m</c> contains an exact maximal pair of length ≥ ⌊m/(k+1)⌋ (Vmatch <c>-seedlength</c> rule); seeds are
    /// enumerated with the suffix-array maximal-pair engine of <see cref="FindDirectRepeats(string,int,int,int)"/>
    /// and, per seed, every maximal window containing it is generated from the first k+1 mismatches to its
    /// left and right (all splits a + b = k, as in REPuter's maximum-error extension tables). Duplicates found
    /// from several seeds are reported once. Unlike Vmatch's default output (one E-value-best extension per
    /// seed, not a set definition), <b>all</b> maximal k-mismatch repeats are returned; with
    /// <paramref name="excludeContained"/> the output is identical to <c>vmatch -h k -allmax</c> (verified on
    /// 3 040 random cases, 1.09 M repeats). Cost O(n log² n + s·k + z) for s seeds; a small ⌊m/(k+1)⌋ makes s
    /// grow quadratically.
    /// Results are ordered by (FirstPosition, SecondPosition, Length).
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 10; must be ≥ 2 and &gt; <paramref name="maxMismatches"/>).</param>
    /// <param name="maxMismatches">Maximum Hamming distance k between the copies (default: 1; ≥ 0).</param>
    /// <param name="maxLength">Maximum repeat length (default: unbounded; ≥ <paramref name="minLength"/>).</param>
    /// <param name="minSpacing">Minimum number of bases between the copies (default: 1; negative admits overlap).</param>
    /// <param name="excludeContained">
    /// When <c>true</c>, also drop repeats contained (both copies as intervals) in a k-mismatch repeat on
    /// another diagonal — the literal Vmatch Appendix A maximality; output then equals Vmatch
    /// <c>vmatch -l minLength -h k -allmax</c> (k ≥ 1). Default <c>false</c>: per-diagonal maximality, which for
    /// k = 0 coincides with Gusfield maximal pairs / <see cref="FindDirectRepeats(string,int,int,int)"/> (Vmatch's
    /// exact-repeat output keeps such cross-diagonal-contained pairs too, e.g. (0,2,6) inside (0,1,7) in A×8).
    /// Containment is decided before the <paramref name="maxLength"/>/<paramref name="minSpacing"/> filters.
    /// </param>
    /// <returns>Maximal k-mismatch direct repeats, sorted by (FirstPosition, SecondPosition, Length).</returns>
    public static IEnumerable<ApproximateDirectRepeatResult> FindApproximateDirectRepeats(
        DnaSequence sequence,
        int minLength = 10,
        int maxMismatches = 1,
        int maxLength = int.MaxValue,
        int minSpacing = 1,
        bool excludeContained = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateApproximateDirectParameters(minLength, maxMismatches, maxLength);
        return FindApproximateDirectRepeatsCore(sequence.Sequence, minLength, maxMismatches, maxLength, minSpacing, excludeContained);
    }

    /// <summary>
    /// Finds maximal k-mismatch direct repeats in a raw sequence string (case-insensitive; non-ACGT symbols
    /// are mismatches). <c>null</c> or empty input yields no results. See
    /// <see cref="FindApproximateDirectRepeats(DnaSequence,int,int,int,int,bool)"/> for the definition.
    /// </summary>
    public static IEnumerable<ApproximateDirectRepeatResult> FindApproximateDirectRepeats(
        string sequence,
        int minLength = 10,
        int maxMismatches = 1,
        int maxLength = int.MaxValue,
        int minSpacing = 1,
        bool excludeContained = false)
    {
        ValidateApproximateDirectParameters(minLength, maxMismatches, maxLength);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<ApproximateDirectRepeatResult>();

        return FindApproximateDirectRepeatsCore(sequence.ToUpperInvariant(), minLength, maxMismatches, maxLength, minSpacing, excludeContained);
    }

    private static void ValidateApproximateDirectParameters(int minLength, int maxMismatches, int maxLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        // A window of length ≥ minLength must contain at least one matching position (and hence an exact seed).
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(maxMismatches, minLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);
    }

    /// <summary>A maximal k-mismatch window on diagonal <c>Diagonal = j − i</c>: first copy [Start, Start + Length).</summary>
    private readonly record struct MismatchWindow(int Start, int Diagonal, int Length, int Mismatches);

    private static List<ApproximateDirectRepeatResult> FindApproximateDirectRepeatsCore(
        string seq,
        int minLength,
        int k,
        int maxLength,
        int minSpacing,
        bool excludeContained)
    {
        var results = new List<ApproximateDirectRepeatResult>();
        int n = seq.Length;
        if (n <= minLength)
            return results;

        List<MismatchWindow> windows = FindMaximalMismatchWindows(seq, minLength, k);
        if (excludeContained)
            windows = RemoveCrossDiagonalContained(windows);

        foreach (var w in windows)
        {
            if (w.Length > maxLength) continue;
            long spacing = (long)w.Diagonal - w.Length;
            if (spacing < minSpacing) continue;
            int second = w.Start + w.Diagonal;
            results.Add(new ApproximateDirectRepeatResult(
                FirstPosition: w.Start,
                SecondPosition: second,
                Length: w.Length,
                Mismatches: w.Mismatches,
                Spacing: (int)spacing,
                FirstCopy: seq.Substring(w.Start, w.Length),
                SecondCopy: seq.Substring(second, w.Length)));
        }

        results.Sort(static (x, y) =>
        {
            int c = x.FirstPosition.CompareTo(y.FirstPosition);
            if (c != 0) return c;
            c = x.SecondPosition.CompareTo(y.SecondPosition);
            return c != 0 ? c : x.Length.CompareTo(y.Length);
        });
        return results;
    }

    /// <summary>
    /// Every per-diagonal maximal window with ≤ k mismatches and length ≥ minLength (no other filter), each once.
    /// </summary>
    private static List<MismatchWindow> FindMaximalMismatchWindows(string seq, int minLength, int k)
    {
        int[] codes = AcgtCodes(seq);
        int seedLength = Math.Max(1, minLength / (k + 1));
        return FindMaximalMismatchWindows(codes, codes, EnumerateForwardMaximalPairs(seq, seedLength, int.MaxValue), minLength, k);
    }

    /// <summary>A/C/G/T → 0..3, every other symbol → −1 (never matches), for an upper-cased sequence.</summary>
    private static int[] AcgtCodes(string seq)
    {
        var codes = new int[seq.Length];
        for (int p = 0; p < seq.Length; p++)
            codes[p] = AcgtCode(seq[p]);
        return codes;
    }

    /// <summary>
    /// Every per-diagonal maximal window with ≤ k mismatches and length ≥ minLength between texts
    /// <paramref name="u"/> and <paramref name="v"/> (codes; −1 never matches), grown from exact seeds
    /// (<c>P</c> in u, <c>Q</c> in v; for a self-comparison <c>u == v</c> the pair is taken with P &lt; Q).
    /// <c>Diagonal = Q − P</c> (negative allowed): the window [Start, Start + Length) of u faces
    /// [Start + Diagonal, …) of v. Each window once.
    /// </summary>
    private static List<MismatchWindow> FindMaximalMismatchWindows(
        int[] u, int[] v, IEnumerable<RawMaximalPair> seeds, int minLength, int k)
    {
        bool self = ReferenceEquals(u, v);
        var windows = new List<MismatchWindow>();
        var seen = new HashSet<(int, int, int)>();
        var left = new int[k + 2];
        var right = new int[k + 2];

        foreach (var seed in seeds)
        {
            int i = self ? Math.Min(seed.P, seed.Q) : seed.P;
            int d = (self ? Math.Max(seed.P, seed.Q) : seed.Q) - i;

            // left[a] = u-position of the a-th mismatch to the left of the seed (a = 1..k+1), right[b]
            // likewise to the right; scanning stops after k + 1 mismatches or at the diagonal boundary.
            int leftLimit = Math.Max(0, -d); // first u-index whose partner v-index is ≥ 0
            int leftCount = 0;
            for (int p = i - 1; p >= leftLimit && leftCount <= k; p--)
            {
                if (!IsAcgtMatch(u, v, p, p + d)) left[++leftCount] = p;
            }

            int rightCount = 0;
            int rightLimit = Math.Min(u.Length, v.Length - d); // u-index where either text would end
            for (int p = i + seed.Length; p < rightLimit && rightCount <= k; p++)
            {
                if (!IsAcgtMatch(u, v, p, p + d)) right[++rightCount] = p;
            }

            int maxA = Math.Min(k, leftCount);
            for (int a = 0; a <= maxA; a++)
            {
                int b = Math.Min(k - a, rightCount);
                // Maximal: each side is blocked by a mismatch that would exceed k, or by the diagonal boundary
                // (b < k − a already means the right side reached its boundary).
                if (a + b != k && a != leftCount) continue;

                int start = a < leftCount ? left[a + 1] + 1 : leftLimit;
                int end = b < rightCount ? right[b + 1] : rightLimit;
                int length = end - start;
                if (length < minLength) continue;
                if (seen.Add((start, d, length)))
                    windows.Add(new MismatchWindow(start, d, length, a + b));
            }
        }

        return windows;
    }

    /// <summary>
    /// Drops every window contained (both copies, as intervals) in a window on another diagonal — the literal
    /// Vmatch Appendix A maximality ("not contained in another k-mismatch match"), which Vmatch <c>-h k -allmax</c>
    /// applies. Windows on one diagonal are never nested (both maximal), so per diagonal they are sorted by start
    /// with increasing ends and a containing window on diagonal d′ is found by one binary search. A window of
    /// length L can only be contained in one of length ≥ L + |d − d′|, which bounds the diagonals probed.
    /// </summary>
    private static List<MismatchWindow> RemoveCrossDiagonalContained(List<MismatchWindow> windows)
    {
        if (windows.Count < 2)
            return windows;

        var byDiagonal = new Dictionary<int, (int[] Starts, int[] Ends)>();
        foreach (var group in windows.GroupBy(w => w.Diagonal))
        {
            var sorted = group.OrderBy(w => w.Start).ToArray();
            byDiagonal[group.Key] = (sorted.Select(w => w.Start).ToArray(), sorted.Select(w => w.Start + w.Length).ToArray());
        }

        int longest = windows.Max(w => w.Length);
        var kept = new List<MismatchWindow>(windows.Count);
        foreach (var w in windows)
        {
            bool contained = false;
            int reach = longest - w.Length;
            for (int delta = -reach; delta <= reach && !contained; delta++)
            {
                if (delta == 0 || !byDiagonal.TryGetValue(w.Diagonal + delta, out var other)) continue;
                // Window on d′ = d + delta must cover [Start + min(0, −delta), Start + Length + max(0, −delta)).
                long needStart = w.Start + Math.Min(0, -delta);
                long needEnd = (long)w.Start + w.Length + Math.Max(0, -delta);
                int idx = Array.BinarySearch(other.Starts, (int)Math.Max(needStart, int.MinValue));
                if (idx < 0) idx = ~idx - 1;
                contained = idx >= 0 && other.Ends[idx] >= needEnd;
            }
            if (!contained) kept.Add(w);
        }

        return kept;
    }

    private static bool IsAcgtMatch(int[] u, int[] v, int p, int q) => u[p] >= 0 && u[p] == v[q];

    #endregion

    #region Degenerate Repeats — Vmatch -h / -e, direct and palindromic (REPuter)

    /// <summary>
    /// Finds <b>maximal degenerate repeats</b> — k-mismatch (Hamming) or k-differences (unit-cost edit
    /// distance) repeats, direct or reverse-complement (palindromic) — with the exact Vmatch / REPuter
    /// semantics (Kurtz et al. 2001 NAR 29:4633; Vmatch manual Appendix A; <c>vmatch [-p] -l m (-h|-e) k -allmax</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Definitions</b> (Vmatch manual App. A "Basic Notions"). A match <c>(l, i, r, j)</c> pairs the left
    /// instance <c>u = S[i..i+l)</c> with the right instance <c>w = S[j..j+r)</c>: <i>direct</i> when
    /// <c>u ≈ w</c> (<c>i &lt; j</c>), <i>palindromic</i> when <c>u ≈ wcc(reverse(w))</c>, the reverse complement
    /// (<c>i ≤ j</c>). With <see cref="ApproximateRepeatDistance.Hamming"/> <c>≈</c> means <c>l = r</c> and
    /// <c>d_H(u, w') ≤ k</c> (a <i>k-mismatch match</i>); with <see cref="ApproximateRepeatDistance.Edit"/> it
    /// means <c>d_E(u, w') ≤ k</c>, the minimum number of mismatches, insertions and deletions (a
    /// <i>k-differences match</i>; <c>l ≠ r</c> allowed). A match is <i>contained</i> in <c>(l′, i′, r′, j′)</c>
    /// when <c>i′ ≤ i ≤ i+l ≤ i′+l′</c> and <c>j′ ≤ j ≤ j+r ≤ j′+r′</c>; it is <i>maximal</i> when it is not
    /// contained in another match of the same kind. Both instances must have length ≥ <paramref name="minLength"/>
    /// (Vmatch <c>-l</c>). Wildcards (N, IUPAC codes, U, gaps) always mismatch.
    /// </para>
    /// <para>
    /// <b>Vmatch conventions reproduced</b> (source <c>kurtz/extendED.c</c>, <c>mcontain.c</c>, Vmatch 2.3.1):
    /// (1) a direct k-differences match is admitted only when its right instance is not embedded in the left
    /// one (<c>i + l &lt; j + r</c>) and, when the instances overlap, the non-overlapping part
    /// <c>(j − i) + (j + r) − (i + l)</c> exceeds the distance (<c>acceptmatch</c>; this removes trivial
    /// self-alignments such as <c>S[0..n)</c> against <c>S[1..n)</c>); (2) a palindromic match and its mirror
    /// <c>(r, j, l, i)</c> are the same pair of strings, so maximality is decided over both orientations and the
    /// orientation with <c>i ≤ j</c> is reported (for <c>i = j</c>, <c>l ≠ r</c>, both). Maximality is decided
    /// before the <paramref name="maxLength"/> and <paramref name="minSpacing"/> filters.
    /// </para>
    /// <para>
    /// <b>Algorithm</b> (REPuter / Vmatch seed-and-extend, complete): every match with an instance of length
    /// ≥ m and ≤ k errors contains an exact match of length ≥ ⌊m/(k+1)⌋ (pigeonhole; Vmatch
    /// <c>-seedlength</c> rule). Seeds are the exact maximal pairs (direct: <see cref="FindDirectRepeats(string,int,int,int)"/>'s
    /// engine; palindromic: the S·#·revcomp(S) engine of <see cref="FindReverseComplementRepeats(string,int,int,int)"/>,
    /// both orientations). Hamming: per seed, the windows bounded by the (k+1)-th mismatches on each side
    /// (as in <see cref="FindApproximateDirectRepeats(string,int,int,int,int,bool)"/>). Edit: per seed, the
    /// greedy furthest-reaching fronts of Ukkonen 1985 / Myers 1986 (Vmatch <c>frontSEP.c</c>) to the left and
    /// right, ≤ k errors in total, every combination of a left and a right front point is a candidate whose
    /// distance is the smallest a + b producing it; candidates are de-duplicated, filtered by rule (1) and
    /// reduced to the maximal ones (a containing candidate is found by an output-sensitive segment-tree search).
    /// Vmatch additionally stops a left extension that crosses another exact match of length ≥ the seed length
    /// (<c>evalentrybackward</c> "seed … detected while scanning"), assuming it is found from that seed; for
    /// k-differences matches this loses some maximal matches (and then reports contained ones). This method
    /// returns the complete set: identical to Vmatch with that shortcut disabled, and to a brute-force
    /// enumeration of the definition (Evidence REP-DIRECT-001 §WP8). Results ordered by
    /// (FirstPosition, SecondPosition, FirstLength, SecondLength).
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum length of each instance (Vmatch <c>-l</c>; default 10, ≥ 2, &gt; <paramref name="maxDifferences"/>).</param>
    /// <param name="maxDifferences">Maximum distance k (Vmatch <c>-h k</c> / <c>-e k</c>; default 1, ≥ 1).</param>
    /// <param name="distance">Hamming (k-mismatch) or unit-cost edit distance (k-differences; default).</param>
    /// <param name="reverseComplement"><c>false</c> (default): direct repeats; <c>true</c>: palindromic (Vmatch <c>-p</c>).</param>
    /// <param name="maxLength">Maximum length of each instance (default: unbounded; ≥ <paramref name="minLength"/>).</param>
    /// <param name="minSpacing">
    /// Minimum <c>Spacing = SecondPosition − FirstPosition − FirstLength</c> (default 1: non-overlapping copies, as
    /// the other repeat finders; <c>int.MinValue</c> returns the complete Vmatch set).
    /// </param>
    /// <returns>Maximal degenerate repeats, sorted by (FirstPosition, SecondPosition, FirstLength, SecondLength).</returns>
    public static IEnumerable<DegenerateRepeatResult> FindDegenerateRepeats(
        DnaSequence sequence,
        int minLength = 10,
        int maxDifferences = 1,
        ApproximateRepeatDistance distance = ApproximateRepeatDistance.Edit,
        bool reverseComplement = false,
        int maxLength = int.MaxValue,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateDegenerateParameters(minLength, maxDifferences, distance, maxLength);
        return FindDegenerateRepeatsCore(sequence.Sequence, minLength, maxDifferences, distance, reverseComplement, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds maximal degenerate repeats in a raw sequence string (case-insensitive; non-ACGT symbols always
    /// mismatch). <c>null</c> or empty input yields no results. See
    /// <see cref="FindDegenerateRepeats(DnaSequence,int,int,ApproximateRepeatDistance,bool,int,int)"/>.
    /// </summary>
    public static IEnumerable<DegenerateRepeatResult> FindDegenerateRepeats(
        string sequence,
        int minLength = 10,
        int maxDifferences = 1,
        ApproximateRepeatDistance distance = ApproximateRepeatDistance.Edit,
        bool reverseComplement = false,
        int maxLength = int.MaxValue,
        int minSpacing = 1)
    {
        ValidateDegenerateParameters(minLength, maxDifferences, distance, maxLength);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<DegenerateRepeatResult>();

        return FindDegenerateRepeatsCore(sequence.ToUpperInvariant(), minLength, maxDifferences, distance, reverseComplement, maxLength, minSpacing);
    }

    private static void ValidateDegenerateParameters(int minLength, int maxDifferences, ApproximateRepeatDistance distance, int maxLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxDifferences, 1); // Vmatch: "-h k" / "-e k" with k > 0
        // An instance of length ≥ minLength with ≤ k errors must contain a matching position (an exact seed).
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(maxDifferences, minLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxLength, minLength);
        if (distance is not (ApproximateRepeatDistance.Hamming or ApproximateRepeatDistance.Edit))
            throw new ArgumentOutOfRangeException(nameof(distance), distance, "Unknown distance.");
    }

    /// <summary>A degenerate match between text u = S and text v (S, or R = revcomp(S)): u[U..U+UL) ≈ v[V..V+VL).</summary>
    private readonly record struct DegenerateMatch(int U, int ULength, int V, int VLength, int Distance)
    {
        public int UEnd => U + ULength;
        public int VEnd => V + VLength;
    }

    private static List<DegenerateRepeatResult> FindDegenerateRepeatsCore(
        string seq,
        int minLength,
        int k,
        ApproximateRepeatDistance distance,
        bool palindromic,
        int maxLength,
        int minSpacing)
    {
        var results = new List<DegenerateRepeatResult>();
        int n = seq.Length;
        if (n < minLength || (!palindromic && n <= minLength))
            return results;

        int[] u = AcgtCodes(seq);
        int[] v = u;
        int seedLength = Math.Max(1, minLength / (k + 1));
        List<RawMaximalPair> seeds;
        if (palindromic)
        {
            // R = revcomp(S): R[q] pairs with S[n − 1 − q].
            v = new int[n];
            for (int q = 0; q < n; q++)
            {
                int c = u[n - 1 - q];
                v[q] = c >= 0 ? 3 - c : -1; // A C G T = 0 1 2 3 → complement 3 − c
            }
            seeds = EnumerateReverseComplementSeeds(seq, seedLength, int.MaxValue);
        }
        else
        {
            seeds = EnumerateForwardMaximalPairs(seq, seedLength, int.MaxValue);
        }

        List<DegenerateMatch> matches;
        if (distance == ApproximateRepeatDistance.Hamming)
        {
            var windows = RemoveCrossDiagonalContained(FindMaximalMismatchWindows(u, v, seeds, minLength, k));
            matches = windows.ConvertAll(w => new DegenerateMatch(w.Start, w.Length, w.Start + w.Diagonal, w.Length, w.Mismatches));
        }
        else
        {
            matches = RemoveContainedMatches(FindDifferenceMatchCandidates(u, v, seeds, minLength, k, selfDirect: !palindromic));
        }

        foreach (var m in matches)
        {
            int i = m.U;
            int j = palindromic ? n - m.VEnd : m.V;
            if (palindromic && j < i) continue; // the mirror (r, j, l, i) is reported instead
            if (m.ULength > maxLength || m.VLength > maxLength) continue;
            long spacing = (long)j - i - m.ULength;
            if (spacing < minSpacing) continue;
            results.Add(new DegenerateRepeatResult(
                FirstPosition: i,
                FirstLength: m.ULength,
                SecondPosition: j,
                SecondLength: m.VLength,
                Distance: m.Distance,
                Spacing: (int)spacing,
                FirstCopy: seq.Substring(i, m.ULength),
                SecondCopy: seq.Substring(j, m.VLength),
                IsReverseComplement: palindromic));
        }

        results.Sort(static (x, y) =>
        {
            int c = x.FirstPosition.CompareTo(y.FirstPosition);
            if (c != 0) return c;
            c = x.SecondPosition.CompareTo(y.SecondPosition);
            if (c != 0) return c;
            c = x.FirstLength.CompareTo(y.FirstLength);
            return c != 0 ? c : x.SecondLength.CompareTo(y.SecondLength);
        });
        return results;
    }

    /// <summary>
    /// All candidate k-differences matches grown from the seeds: per seed, the furthest-reaching points of the
    /// left and right greedy fronts (≤ a and ≤ b errors, a + b ≤ k) are combined; candidates with an instance
    /// shorter than <paramref name="minLength"/> are dropped. For a self-comparison (<paramref name="selfDirect"/>)
    /// the instances are ordered by start and the Vmatch <c>acceptmatch</c> rule is applied with the smallest
    /// distance producing the candidate. Each candidate once, with its smallest distance; candidates contained
    /// in another candidate of the same seed are already dropped.
    /// </summary>
    private static List<DegenerateMatch> FindDifferenceMatchCandidates(
        int[] u, int[] v, List<RawMaximalPair> seeds, int minLength, int k, bool selfDirect)
    {
        int width = 2 * k + 1;
        var leftFront = new int[(k + 1) * width];
        var rightFront = new int[(k + 1) * width];
        var perSeed = new Dictionary<(int, int, int, int), int>();
        var seedMatches = new List<DegenerateMatch>();
        var best = new Dictionary<(int, int, int, int), int>();

        foreach (var seed in seeds)
        {
            int p = seed.P, q = seed.Q, m = seed.Length;
            if (selfDirect && p > q) (p, q) = (q, p);

            ComputeEditFronts(u, v, p, q, p, q, -1, k, leftFront);
            ComputeEditFronts(u, v, p + m, q + m, u.Length - (p + m), v.Length - (q + m), +1, k, rightFront);

            perSeed.Clear();
            for (int a = 0; a <= k; a++)
            {
                for (int dl = -a; dl <= a; dl++)
                {
                    int xl = leftFront[a * width + dl + k];
                    if (xl < 0) continue;
                    int us = p - xl, vs = q - xl - dl;
                    for (int b = 0; b <= k - a; b++)
                    {
                        for (int dr = -b; dr <= b; dr++)
                        {
                            int xr = rightFront[b * width + dr + k];
                            if (xr < 0) continue;
                            int ue = p + m + xr, ve = q + m + xr + dr;
                            if (ue - us < minLength || ve - vs < minLength) continue;

                            int s1 = us, e1 = ue, s2 = vs, e2 = ve;
                            if (selfDirect && s1 > s2) (s1, e1, s2, e2) = (s2, e2, s1, e1);
                            var key = (s1, e1, s2, e2);
                            int dist = a + b;
                            if (!perSeed.TryGetValue(key, out int old) || dist < old)
                                perSeed[key] = dist;
                        }
                    }
                }
            }

            seedMatches.Clear();
            foreach (var (key, dist) in perSeed)
            {
                if (selfDirect && !IsAcceptedSelfMatch(key.Item1, key.Item2, key.Item3, key.Item4, dist)) continue;
                seedMatches.Add(new DegenerateMatch(key.Item1, key.Item2 - key.Item1, key.Item3, key.Item4 - key.Item3, dist));
            }

            for (int x = 0; x < seedMatches.Count; x++)
            {
                var c = seedMatches[x];
                bool contained = false;
                for (int y = 0; y < seedMatches.Count && !contained; y++)
                    contained = y != x && Contains(seedMatches[y], c);
                if (contained) continue;
                var key = (c.U, c.UEnd, c.V, c.VEnd);
                if (!best.TryGetValue(key, out int old) || c.Distance < old)
                    best[key] = c.Distance;
            }
        }

        var result = new List<DegenerateMatch>(best.Count);
        foreach (var (key, dist) in best)
            result.Add(new DegenerateMatch(key.Item1, key.Item2 - key.Item1, key.Item3, key.Item4 - key.Item3, dist));
        return result;
    }

    /// <summary>
    /// Vmatch <c>acceptmatch</c> (kurtz/extendED.c) for a direct self-match with instances [s1, e1) and
    /// [s2, e2), s1 ≤ s2: rejects identical starts, a right instance embedded in the left one, and overlapping
    /// instances whose non-overlapping part <c>(s2 − s1) + e2 − e1</c> is at most the distance.
    /// </summary>
    private static bool IsAcceptedSelfMatch(int s1, int e1, int s2, int e2, int dist)
    {
        if (s1 >= s2) return false;
        if (e1 <= s2) return true;     // no overlap
        if (e1 >= e2) return false;    // embedded
        return (s2 - s1) + e2 - e1 > dist;
    }

    private static bool Contains(DegenerateMatch outer, DegenerateMatch inner) =>
        outer.U <= inner.U && inner.UEnd <= outer.UEnd && outer.V <= inner.V && inner.VEnd <= outer.VEnd;

    /// <summary>
    /// Greedy furthest-reaching fronts (Ukkonen 1985; Myers 1986; Vmatch <c>frontSEP.c</c>) of an alignment
    /// anchored at u-position <paramref name="pu"/> / v-position <paramref name="pv"/> and extended in direction
    /// <paramref name="dir"/> (+1 right, −1 left) over at most <paramref name="lu"/> / <paramref name="lv"/>
    /// symbols. <c>fronts[a·(2k+1) + δ + k]</c> = largest number x of u-symbols consumed by an alignment with
    /// a differences that ends on diagonal δ (v-symbols consumed = x + δ), or −1 when no such point exists.
    /// Unit costs; codes −1 never match.
    /// </summary>
    private static void ComputeEditFronts(int[] u, int[] v, int pu, int pv, int lu, int lv, int dir, int k, int[] fronts)
    {
        int width = 2 * k + 1;
        Array.Fill(fronts, -1);
        fronts[k] = SlideFront(u, v, pu, pv, lu, lv, dir, 0, 0);

        for (int a = 1; a <= k; a++)
        {
            int prev = (a - 1) * width + k, cur = a * width + k;
            for (int d = -a; d <= a; d++)
            {
                int x = -1;
                if (d > -a && d < a)
                    x = Math.Max(x, ValidFrontPoint(fronts[prev + d], 1, d, lu, lv));           // mismatch
                if (d + 1 <= a - 1)
                    x = Math.Max(x, ValidFrontPoint(fronts[prev + d + 1], 1, d, lu, lv));       // u-symbol deleted
                if (d - 1 >= -(a - 1))
                    x = Math.Max(x, ValidFrontPoint(fronts[prev + d - 1], 0, d, lu, lv));       // v-symbol inserted
                fronts[cur + d] = x < 0 ? -1 : SlideFront(u, v, pu, pv, lu, lv, dir, x, d);
            }
        }
    }

    private static int ValidFrontPoint(int previous, int step, int d, int lu, int lv)
    {
        if (previous < 0) return -1;
        int x = previous + step;
        return x <= lu && x + d <= lv ? x : -1;
    }

    private static int SlideFront(int[] u, int[] v, int pu, int pv, int lu, int lv, int dir, int x, int d)
    {
        if (dir > 0)
        {
            while (x < lu && x + d < lv && IsAcgtMatch(u, v, pu + x, pv + x + d)) x++;
        }
        else
        {
            while (x < lu && x + d < lv && IsAcgtMatch(u, v, pu - 1 - x, pv - 1 - x - d)) x++;
        }
        return x;
    }

    /// <summary>
    /// Keeps the matches not contained in another one (Vmatch App. A maximality). Matches sorted by U; a max
    /// segment tree over UEnd enumerates exactly the matches with U′ ≤ U and UEnd′ ≥ UEnd, which are then
    /// checked on the v-interval. Input has no duplicates.
    /// </summary>
    private static List<DegenerateMatch> RemoveContainedMatches(List<DegenerateMatch> matches)
    {
        int count = matches.Count;
        if (count < 2)
            return matches;

        var sorted = matches.ToArray();
        Array.Sort(sorted, static (x, y) =>
        {
            int c = x.U.CompareTo(y.U);
            return c != 0 ? c : y.UEnd.CompareTo(x.UEnd);
        });

        int size = 1;
        while (size < count) size <<= 1;
        var tree = new int[2 * size];
        Array.Fill(tree, int.MinValue);
        for (int x = 0; x < count; x++) tree[size + x] = sorted[x].UEnd;
        for (int x = size - 1; x >= 1; x--) tree[x] = Math.Max(tree[2 * x], tree[2 * x + 1]);

        var kept = new List<DegenerateMatch>(count);
        var stack = new Stack<(int Node, int Lo, int Hi)>();
        int upper = 0;
        for (int x = 0; x < count; x++)
        {
            var m = sorted[x];
            while (upper < count && sorted[upper].U <= m.U) upper++;

            bool contained = false;
            stack.Clear();
            stack.Push((1, 0, size - 1));
            while (stack.Count > 0 && !contained)
            {
                var (node, lo, hi) = stack.Pop();
                if (lo >= upper || tree[node] < m.UEnd) continue;
                if (lo == hi)
                {
                    contained = lo != x && sorted[lo].V <= m.V && m.VEnd <= sorted[lo].VEnd;
                    continue;
                }
                int mid = (lo + hi) >> 1;
                stack.Push((2 * node + 1, mid + 1, hi));
                stack.Push((2 * node, lo, mid));
            }

            if (!contained) kept.Add(m);
        }

        return kept;
    }

    #endregion

    #region Supermaximal Repeats

    /// <summary>
    /// Finds <b>supermaximal repeats</b> (Gusfield 1997 §7.12.1): maximal repeats that never occur as a
    /// substring of any other maximal repeat (also Vmatch <c>-supermax</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// A string α is a maximal repeat when it has a maximal pair (see <see cref="FindDirectRepeats(string,int,int,int)"/>);
    /// it is supermaximal when no other maximal repeat contains it. Gusfield Theorem 7.12.4: α is supermaximal
    /// iff its locus is an internal suffix-tree node whose children are all leaves and whose leaves have pairwise
    /// distinct left characters. On the suffix array this is an lcp-interval with only singleton children
    /// (a local maximum of the LCP array) whose suffixes have pairwise distinct left characters (Abouelhoda,
    /// Kurtz &amp; Ohlebusch 2004). A suffix at position 0 or preceded by a non-ACGT symbol has an undefined left
    /// character, distinct from every other.
    /// </para>
    /// <para>
    /// Only A/C/G/T match (case-insensitive); other symbols never match. Each supermaximal repeat is reported
    /// once with every (possibly overlapping) occurrence, positions ascending; results are ordered by first
    /// occurrence, then length. O(n log² n) (shared suffix-array + LCP helpers of <see cref="SequenceComplexity"/>).
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5, must be ≥ 1).</param>
    /// <returns>Supermaximal repeats with all their occurrences.</returns>
    public static IEnumerable<SupermaximalRepeatResult> FindSupermaximalRepeats(DnaSequence sequence, int minLength = 5)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        return FindSupermaximalRepeatsCore(sequence.Sequence, minLength);
    }

    /// <summary>
    /// Finds supermaximal repeats in a raw sequence string (case-insensitive; non-ACGT symbols never match).
    /// <c>null</c> or empty input yields no results. See <see cref="FindSupermaximalRepeats(DnaSequence,int)"/>.
    /// </summary>
    public static IEnumerable<SupermaximalRepeatResult> FindSupermaximalRepeats(string sequence, int minLength = 5)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<SupermaximalRepeatResult>();

        return FindSupermaximalRepeatsCore(sequence.ToUpperInvariant(), minLength);
    }

    private static List<SupermaximalRepeatResult> FindSupermaximalRepeatsCore(string seq, int minLength)
    {
        var results = new List<SupermaximalRepeatResult>();
        int n = seq.Length;
        if (n < 2)
            return results;

        var symbols = new int[n];
        for (int p = 0; p < n; p++)
        {
            int code = AcgtCode(seq[p]);
            symbols[p] = code >= 0 ? code : 4 + p;
        }

        int[] sa = SequenceComplexity.BuildSuffixArray(symbols);
        int[] lcp = SequenceComplexity.BuildLcpArray(symbols, sa);

        // Local maxima of the LCP array: runs lcp[lb+1..rb] = ℓ with lcp[lb] < ℓ and lcp[rb+1] < ℓ
        // (lcp[0] = 0, virtual lcp[n] = 0) are exactly the lcp-intervals whose children are all leaves.
        Span<bool> seenLeft = stackalloc bool[4];
        int r = 1;
        while (r < n)
        {
            int ell = lcp[r];
            if (ell == 0 || lcp[r - 1] >= ell)
            {
                r++;
                continue;
            }

            int rb = r;
            while (rb + 1 < n && lcp[rb + 1] == ell) rb++;
            bool localMax = rb + 1 >= n || lcp[rb + 1] < ell;
            int lb = r - 1;

            if (localMax && ell >= minLength)
            {
                seenLeft.Clear();
                bool leftDiverse = true;
                for (int x = lb; x <= rb && leftDiverse; x++)
                {
                    int p = sa[x];
                    if (p == 0 || symbols[p - 1] >= 4) continue; // undefined left character: always distinct
                    int c = symbols[p - 1];
                    if (seenLeft[c]) leftDiverse = false;
                    seenLeft[c] = true;
                }

                if (leftDiverse)
                {
                    var positions = new int[rb - lb + 1];
                    for (int x = lb; x <= rb; x++) positions[x - lb] = sa[x];
                    Array.Sort(positions);
                    results.Add(new SupermaximalRepeatResult(seq.Substring(positions[0], ell), ell, positions));
                }
            }

            r = rb + 1;
        }

        results.Sort(static (a, b) =>
        {
            int c = a.Positions[0].CompareTo(b.Positions[0]);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        });
        return results;
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
    /// most runs; ties go to the unit whose first run appears first in that order. For MISA's rotation + reverse-complement
    /// class table use <see cref="GetCanonicalMotifFrequencies"/>; per-unit-size thresholds: the
    /// <see cref="GetTandemRepeatSummary(DnaSequence,IReadOnlyDictionary{int,int})"/> overload.</description></item>
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

        return SummarizeMicrosatellites(sequence, FindMicrosatellites(sequence, 1, 6, minRepeats).ToList());
    }

    /// <summary>
    /// Gets a summary of all perfect microsatellites with a separate minimum copy number per unit length (MISA-style
    /// definition, e.g. <see cref="MisaDefaultMinRepeats"/> = <c>1-10 2-6 3-5 4-5 5-5 6-5</c>). Fields are defined as in
    /// <see cref="GetTandemRepeatSummary(DnaSequence,int)"/>; the SSR list is
    /// <see cref="FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},CancellationToken,IProgress{double})"/>
    /// with the given thresholds. Unit lengths absent from the map are not searched (their class count is 0).
    /// </summary>
    /// <param name="sequence">DNA sequence to summarize.</param>
    /// <param name="minRepeatsByUnitLength">Unit length (1–6) → minimum number of complete copies (≥ 2).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The map is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A unit length is outside 1–6 or a threshold is &lt; 2.</exception>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength, maxUnitLength: 6);

        return SummarizeMicrosatellites(
            sequence,
            FindMicrosatellitesCore(sequence.Sequence, thresholds, CancellationToken.None, null).ToList());
    }

    /// <summary>
    /// Summary statistics over the microsatellites found with a per-unit-length threshold map (unit lengths 1–6) and
    /// an explicit scan convention; with <see cref="MicrosatelliteScanMode.MisaRegex"/> the counts are those of
    /// misa.pl's SSR list (its <c>.statistics</c> totals per unit size).
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="minRepeatsByUnitLength">Unit length (1–6) → minimum copies (≥ 2).</param>
    /// <param name="scanMode">Scan convention.</param>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        DnaSequence sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        MicrosatelliteScanMode scanMode)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength, maxUnitLength: 6);
        ValidateScanMode(scanMode);

        return SummarizeMicrosatellites(
            sequence,
            FindMicrosatellitesCore(sequence.Sequence, thresholds, CancellationToken.None, null, scanMode).ToList());
    }

    /// <summary>
    /// Raw-string (N/IUPAC-tolerant) counterpart of <see cref="GetTandemRepeatSummary(DnaSequence,int)"/>: the summary
    /// aggregates exactly <see cref="FindMicrosatellites(string,int,int,int)"/> with unit lengths 1–6 (case-insensitive;
    /// only A/C/G/T form units, so N runs and IUPAC codes never belong to a microsatellite and interrupt runs —
    /// MISA <c>[acgt]</c>). <c>PercentageOfSequence</c> uses the full input length, non-ACGT symbols included, as
    /// denominator (misa.pl "Total size of examined sequences (bp)" = <c>length $seq</c>). <c>null</c> or empty input
    /// yields the empty summary (all counts 0, no longest repeat).
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="minRepeats">Minimum number of complete copies (≥ 2; default 3), applied to every unit length.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minRepeats"/> is less than 2.</exception>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        string sequence,
        int minRepeats = 3)
    {
        var found = FindMicrosatellites(sequence, 1, 6, minRepeats).ToList();
        return SummarizeMicrosatellites(sequence?.Length ?? 0, found);
    }

    /// <summary>
    /// Raw-string (N/IUPAC-tolerant) counterpart of
    /// <see cref="GetTandemRepeatSummary(DnaSequence,IReadOnlyDictionary{int,int})"/>; the SSR list is
    /// <see cref="FindMicrosatellites(string,IReadOnlyDictionary{int,int},CancellationToken,IProgress{double})"/>.
    /// Denominator and <c>null</c>/empty handling as in <see cref="GetTandemRepeatSummary(string,int)"/>.
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="minRepeatsByUnitLength">Unit length (1–6) → minimum number of complete copies (≥ 2).</param>
    /// <exception cref="ArgumentNullException">The map is null.</exception>
    /// <exception cref="ArgumentException">The map is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A unit length is outside 1–6 or a threshold is &lt; 2.</exception>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        string sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength)
        => GetTandemRepeatSummary(sequence, minRepeatsByUnitLength, MicrosatelliteScanMode.MaximalRuns);

    /// <summary>
    /// Raw-string (N/IUPAC-tolerant) counterpart of
    /// <see cref="GetTandemRepeatSummary(DnaSequence,IReadOnlyDictionary{int,int},MicrosatelliteScanMode)"/>; the SSR
    /// list is <see cref="FindMicrosatellites(string,IReadOnlyDictionary{int,int},MicrosatelliteScanMode,CancellationToken,IProgress{double})"/>.
    /// With <see cref="MicrosatelliteScanMode.MisaRegex"/> the total and the per-unit-size counts equal misa.pl's
    /// <c>.statistics</c> ("Total number of identified SSRs", "Distribution to different repeat type classes") on
    /// sequences containing N or IUPAC codes too. Denominator and <c>null</c>/empty handling as in
    /// <see cref="GetTandemRepeatSummary(string,int)"/>.
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="minRepeatsByUnitLength">Unit length (1–6) → minimum number of complete copies (≥ 2).</param>
    /// <param name="scanMode">Scan convention.</param>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        string sequence,
        IReadOnlyDictionary<int, int> minRepeatsByUnitLength,
        MicrosatelliteScanMode scanMode)
    {
        var thresholds = ThresholdsFromMap(minRepeatsByUnitLength, maxUnitLength: 6);
        ValidateScanMode(scanMode);

        if (string.IsNullOrEmpty(sequence))
            return SummarizeMicrosatellites(0, []);

        return SummarizeMicrosatellites(
            sequence.Length,
            FindMicrosatellitesCore(sequence.ToUpperInvariant(), thresholds, CancellationToken.None, null, scanMode).ToList());
    }

    private static TandemRepeatSummary SummarizeMicrosatellites(
        DnaSequence sequence, List<MicrosatelliteResult> microsatellites)
        => SummarizeMicrosatellites(sequence.Length, microsatellites);

    private static TandemRepeatSummary SummarizeMicrosatellites(
        int sequenceLength, List<MicrosatelliteResult> microsatellites)
    {
        var byType = microsatellites
            .GroupBy(m => m.RepeatType)
            .ToDictionary(g => g.Key, g => g.ToList());

        int totalBases = microsatellites.Sum(m => m.TotalLength);

        // PercentageOfSequence is the fraction of the sequence that is tandem-repeat, so it must use the
        // DISTINCT bases covered by any microsatellite (the union of their spans), not the raw sum of repeat
        // lengths: overlapping repeats (e.g. a homopolymer run also matched as a dinucleotide repeat) would
        // otherwise double-count bases and push the percentage above 100. Covered bases ≤ sequence length, so
        // the percentage is always in [0, 100]. TotalRepeatBases keeps the (possibly overlapping) repeat content.
        long coveredBases = CountCoveredBases(microsatellites, sequenceLength);
        double percentageOfSequence = sequenceLength > 0
            ? (double)coveredBases / sequenceLength * 100
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
    /// <see cref="FindInvertedRepeats(DnaSequence,int,int,int,int,int,bool)"/>, which reports only maximal stems (EMBOSS
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
/// Scan convention of the per-unit-length microsatellite search (see
/// <see cref="RepeatFinder.FindMicrosatellites(DnaSequence,IReadOnlyDictionary{int,int},MicrosatelliteScanMode,CancellationToken,IProgress{double})"/>).
/// </summary>
public enum MicrosatelliteScanMode
{
    /// <summary>Each maximal primitive ACGT run reported once at its left end (Kolpakov–Kucherov maximal repetitions; pytrf/Krait).</summary>
    MaximalRuns,

    /// <summary>misa.pl v1.0 regex scan: leftmost greedy match, resumed after each match; non-primitive matches consumed, then rejected.</summary>
    MisaRegex,
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
/// A compound microsatellite (MISA, Thiel et al. 2003): two or more microsatellites whose neighbours in start order
/// are separated by at most the maximal interruption. See
/// <see cref="RepeatFinder.AssembleCompoundMicrosatellites"/> for the chaining rule.
/// </summary>
/// <param name="Start">0-based start of the first component (MISA prints 1-based: add 1).</param>
/// <param name="End">Exclusive 0-based end of the LAST component in start order (= MISA's 1-based inclusive
/// <c>end</c>); when the last component is nested inside an earlier one this is smaller than the maximal end.</param>
/// <param name="Components">Component SSRs in start order (at least two).</param>
/// <param name="Interruptions">Bases between each consecutive component pair (lower case, as in MISA; empty when
/// adjacent or overlapping); one entry per join.</param>
/// <param name="IsOverlapping">True when some consecutive components overlap (MISA type <c>c*</c>).</param>
/// <param name="Notation">MISA SSR notation, e.g. <c>(AT)6ccgt(GA)7</c> or <c>(A)10(AT)6*</c>.</param>
public readonly record struct CompoundMicrosatelliteResult(
    int Start,
    int End,
    IReadOnlyList<MicrosatelliteResult> Components,
    IReadOnlyList<string> Interruptions,
    bool IsOverlapping,
    string Notation)
{
    /// <summary>MISA SSR type: <c>c*</c> when components overlap, otherwise <c>c</c>.</summary>
    public string MisaType => IsOverlapping ? "c*" : "c";

    /// <summary>MISA <c>size</c> column: <see cref="End"/> − <see cref="Start"/>.</summary>
    public int Length => End - Start;
}

/// <summary>
/// Result of approximate (imperfect / interrupted) tandem-repeat detection with the Tandem Repeats Finder
/// model (Benson 1999). Fields mirror the TRF table: <see cref="Start"/> (0-based; TRF prints 1-based) and
/// <see cref="SpanLength"/> give the indices; <see cref="Period"/> is the most common distance between matching
/// characters of adjacent copies (may differ from <see cref="ConsensusSize"/>); <see cref="CopyNumber"/> is the
/// number of copies aligned with the consensus; <see cref="PercentMatches"/> / <see cref="PercentIndels"/> are
/// between ADJACENT copies (exact values; TRF truncates to integers); <see cref="AlignmentScore"/> is the WDP
/// score; <see cref="Consensus"/> starts at the phase of the first repeat base.
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
    int AlignmentScore)
{
    /// <summary>Percentage of A in the repeat region (TRF "A" column, exact; denominator = region length).</summary>
    public double PercentA { get; init; }

    /// <summary>Percentage of C in the repeat region (TRF "C" column, exact).</summary>
    public double PercentC { get; init; }

    /// <summary>Percentage of G in the repeat region (TRF "G" column, exact).</summary>
    public double PercentG { get; init; }

    /// <summary>Percentage of T in the repeat region (TRF "T" column, exact).</summary>
    public double PercentT { get; init; }

    /// <summary>
    /// Shannon entropy of the region's A/C/G/T composition in bits, 0–2, normalised over the A/C/G/T symbols
    /// only. Equals TRF's "Entropy (0-2)" column unless the region contains N or other non-ACGT symbols; see
    /// <see cref="EntropyTrf"/> for the TRF value in every case.
    /// </summary>
    public double Entropy { get; init; }

    /// <summary>
    /// TRF's "Entropy (0-2)" exactly as TRF 4.10.0 computes it: −Σ p·log₂ p over A, C, G, T with each p taken over
    /// the WHOLE region length (N and other symbols included in the denominator, excluded from the sum). Identical
    /// to <see cref="Entropy"/> for pure A/C/G/T regions.
    /// </summary>
    public double EntropyTrf { get; init; }

    /// <summary>
    /// Sequence row of the final (consensus) alignment, left to right: the repeat's symbols with '-' where a
    /// consensus symbol is deleted (TRF alignment file, top line of each pair). Removing the gaps gives the region.
    /// </summary>
    public string? AlignedSequence { get; init; }

    /// <summary>
    /// Consensus row of the final alignment, column-aligned with <see cref="AlignedSequence"/>: the consensus copies
    /// (starting at the phase of the first repeat base) with '-' where a sequence symbol is inserted (TRF alignment
    /// file, bottom line of each pair).
    /// </summary>
    public string? AlignedConsensus { get; init; }

    /// <summary>
    /// Up to <see cref="TandemRepeatsFinderParameters.FlankLength"/> symbols immediately left of the repeat (TRF
    /// <c>-f</c> "Left flanking sequence"; upper case); empty at the sequence start; null when flanks were not requested.
    /// </summary>
    public string? LeftFlank { get; init; }

    /// <summary>Up to <see cref="TandemRepeatsFinderParameters.FlankLength"/> symbols immediately right of the repeat
    /// (TRF <c>-f</c> "Right flanking sequence"); null when flanks were not requested.</summary>
    public string? RightFlank { get; init; }
}

/// <summary>
/// Tandem Repeats Finder parameter set (<c>trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]</c>;
/// TRF 4.10.0 README "parameters"). Defaults are TRF's recommended command line <c>2 7 7 80 10 50 500</c>.
/// </summary>
public sealed record TandemRepeatsFinderParameters
{
    /// <summary>TRF default longest expected TR array (<c>-l 2</c> = 2 million bp).</summary>
    public const int DefaultMaxRepeatLength = 2_000_000;

    /// <summary>The recommended parameter set 2 7 7 80 10 50 500 (TRF README).</summary>
    public static TandemRepeatsFinderParameters Recommended { get; } = new();

    /// <summary>Match weight (TRF "Match"; positive; "a match weight of 2 has proven effective").</summary>
    public int MatchWeight { get; init; } = 2;

    /// <summary>Mismatch penalty, given positive and applied as −value (TRF "Mismatch"; 3–7 recommended, 7 default).</summary>
    public int MismatchPenalty { get; init; } = 7;

    /// <summary>Indel penalty per gap column, given positive (TRF "Delta"; 3–7 recommended, 7 default).</summary>
    public int IndelPenalty { get; init; } = 7;

    /// <summary>Matching probability PM in percent (TRF "PM"): 80 or 75 — "probabilistic data is available for PM
    /// values of 80 and 75". Selects the tuple sizes and the sum-of-heads cut-offs.</summary>
    public int MatchProbability { get; init; } = 80;

    /// <summary>Indel probability PI in percent (TRF "PI"; data for 10 and 20; any 1..100 accepted as by TRF).
    /// Sets the random-walk distance range Δd_max = ⌊2.3·√(PI/100·d)⌋ for d &gt; 20 and the narrow-band radius for
    /// patterns &gt; 20 (max(6, Δd_max) backward, min(2·max(6, Δd_max), ⌊size/3⌋) forward). TRF 4.10.0 exits when a
    /// band exceeds 150 cells (its MAXBANDWIDTH; e.g. PI 20 with patterns ≥ 1365); this implementation has no such limit.</summary>
    public int IndelProbability { get; init; } = 10;

    /// <summary>Minimum alignment score to report (TRF "Minscore"; ≥ 1; 50 recommended).</summary>
    public int MinScore { get; init; } = 50;

    /// <summary>Maximum reported period (TRF "MaxPeriod"; 1..2000; 500 recommended).</summary>
    public int MaxPeriod { get; init; } = 500;

    /// <summary>Longest tandem-repeat array expected, in bp (TRF <c>-l n</c> = n million; default 2 000 000): the
    /// wraparound alignment extends at most this many rows on either side of the candidate.</summary>
    public int MaxRepeatLength { get; init; } = DefaultMaxRepeatLength;

    /// <summary>TRF redundancy elimination (default on; TRF <c>-r</c> turns it off).</summary>
    public bool EliminateRedundancy { get; init; } = true;

    /// <summary>Flanking sequence length reported on each side (0 = none; TRF <c>-f</c> uses 500, <c>-ngs</c> 50).</summary>
    public int FlankLength { get; init; }

    /// <summary>Validates the set against TRF's accepted ranges; throws <see cref="ArgumentOutOfRangeException"/>.</summary>
    public void Validate()
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(MatchWeight, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(MismatchPenalty, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(IndelPenalty, 1);
        // TRF has tuple-size and sum-of-heads data only for PM = 80 and PM = 75.
        if (MatchProbability != 80)
            ArgumentOutOfRangeException.ThrowIfNotEqual(MatchProbability, 75);
        ArgumentOutOfRangeException.ThrowIfLessThan(IndelProbability, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(IndelProbability, 100);
        ArgumentOutOfRangeException.ThrowIfLessThan(MinScore, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(MaxPeriod, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(MaxPeriod, RepeatFinder.MaxApproximatePeriod);
        ArgumentOutOfRangeException.ThrowIfLessThan(MaxRepeatLength, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(FlankLength);
    }
}

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

    /// <summary>
    /// Number of mismatched (non-pairing) positions inside the stem; always 0 for exact stems
    /// (<c>maxMismatches = 0</c>). Mismatches never occur at either end of a stem.
    /// </summary>
    public int Mismatches { get; init; }
}

/// <summary>
/// A scored imperfect inverted repeat (EMBOSS <c>einverted</c> semantics). Coordinates are 0-based and inclusive:
/// <c>LeftArmStart ≤ LeftArmEnd &lt; RightArmStart ≤ RightArmEnd</c>. The three alignment rows are column-aligned:
/// <see cref="LeftArmAlignment"/> reads the left arm 5'→3' (outer → inner), <see cref="RightArmAlignment"/> reads the
/// right arm from its outer end inward (3'→5'), '-' marks a gap, and <see cref="MatchLine"/> has '|' for each
/// Watson–Crick pair. <see cref="Score"/>, <see cref="Matches"/>, <see cref="Mismatches"/> and <see cref="Gaps"/> are
/// einverted's figures: a final trace-back gap step is counted in <see cref="Gaps"/> but not drawn, and an innermost
/// loop-0 pair reached by the last diagonal step contributes to <see cref="Score"/> without being counted or drawn.
/// </summary>
public readonly record struct ScoredInvertedRepeatResult(
    int LeftArmStart,
    int LeftArmEnd,
    int RightArmStart,
    int RightArmEnd,
    int Score,
    int Matches,
    int Mismatches,
    int Gaps,
    string LeftArmAlignment,
    string MatchLine,
    string RightArmAlignment)
{
    /// <summary>Length of the left arm (bases).</summary>
    public int LeftArmLength => LeftArmEnd - LeftArmStart + 1;

    /// <summary>Length of the right arm (bases).</summary>
    public int RightArmLength => RightArmEnd - RightArmStart + 1;

    /// <summary>Bases between the arms.</summary>
    public int LoopLength => RightArmStart - LeftArmEnd - 1;

    /// <summary>Percentage of aligned (non-gap) columns that pair, 100·matches / (matches + mismatches); 0 when none.</summary>
    public double PercentMatches => Matches + Mismatches == 0 ? 0 : 100.0 * Matches / (Matches + Mismatches);
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
/// Maximal exact reverse-complement repeat pair (see
/// <see cref="RepeatFinder.FindReverseComplementRepeats(DnaSequence,int,int,int)"/>):
/// <c>RepeatSequence = S[FirstPosition..+Length)</c> is the reverse complement of
/// <c>SecondSequence = S[SecondPosition..+Length)</c>. Both positions are 0-based forward-strand starts,
/// <c>FirstPosition ≤ SecondPosition</c>; MUMmer <c>repeat-match</c> prints <c>FirstPosition + 1</c> and
/// <c>SecondPosition + Length</c> followed by <c>r</c>. <c>Spacing = SecondPosition − FirstPosition − Length</c>
/// (negative when the copies overlap).
/// </summary>
public readonly record struct ReverseComplementRepeatResult(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    string SecondSequence,
    int Length,
    int Spacing);

/// <summary>
/// Maximal k-mismatch (degenerate) direct repeat (see
/// <see cref="RepeatFinder.FindApproximateDirectRepeats(DnaSequence,int,int,int,int,bool)"/>): the copies
/// <c>FirstCopy = S[FirstPosition..+Length)</c> and <c>SecondCopy = S[SecondPosition..+Length)</c> differ at
/// <c>Mismatches</c> positions (Hamming distance; non-ACGT symbols count as mismatches);
/// <c>Spacing = SecondPosition − FirstPosition − Length</c>.
/// </summary>
public readonly record struct ApproximateDirectRepeatResult(
    int FirstPosition,
    int SecondPosition,
    int Length,
    int Mismatches,
    int Spacing,
    string FirstCopy,
    string SecondCopy);

/// <summary>Distance used by <see cref="RepeatFinder.FindDegenerateRepeats(DnaSequence,int,int,ApproximateRepeatDistance,bool,int,int)"/>.</summary>
public enum ApproximateRepeatDistance
{
    /// <summary>Hamming distance: equal-length instances, mismatches only (Vmatch <c>-h k</c>, "k-mismatch match").</summary>
    Hamming,

    /// <summary>Unit-cost edit (Levenshtein) distance: mismatches, insertions, deletions (Vmatch <c>-e k</c>, "k-differences match").</summary>
    Edit,
}

/// <summary>
/// Maximal degenerate repeat (see
/// <see cref="RepeatFinder.FindDegenerateRepeats(DnaSequence,int,int,ApproximateRepeatDistance,bool,int,int)"/>):
/// left instance <c>FirstCopy = S[FirstPosition..+FirstLength)</c>, right instance
/// <c>SecondCopy = S[SecondPosition..+SecondLength)</c> (forward strand, 0-based; Vmatch prints exactly
/// <c>FirstLength FirstPosition SecondLength SecondPosition Distance</c>). For a direct repeat
/// <c>d(FirstCopy, SecondCopy) = Distance</c>; for a reverse-complement (palindromic) repeat
/// <c>d(FirstCopy, revcomp(SecondCopy)) = Distance</c>. <c>Spacing = SecondPosition − FirstPosition − FirstLength</c>
/// (negative when the instances overlap).
/// </summary>
public readonly record struct DegenerateRepeatResult(
    int FirstPosition,
    int FirstLength,
    int SecondPosition,
    int SecondLength,
    int Distance,
    int Spacing,
    string FirstCopy,
    string SecondCopy,
    bool IsReverseComplement);

/// <summary>
/// Supermaximal repeat (see <see cref="RepeatFinder.FindSupermaximalRepeats(DnaSequence,int)"/>): the repeated
/// string, its length and every 0-based start position (ascending, ≥ 2 occurrences, possibly overlapping).
/// </summary>
public readonly record struct SupermaximalRepeatResult(
    string Sequence,
    int Length,
    IReadOnlyList<int> Positions);

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
