using Seqeron.Genomics.Infrastructure;
namespace Seqeron.Genomics.Chromosome;

/// <summary>
/// Provides genome assembly quality assessment algorithms.
/// Includes N50/L50, contiguity metrics, completeness analysis, and repeat detection.
/// </summary>
public static class GenomeAssemblyAnalyzer
{
    #region Records

    /// <summary>
    /// Assembly statistics summary.
    /// </summary>
    public readonly record struct AssemblyStatistics(
        int TotalSequences,
        long TotalLength,
        long TotalLengthNoGaps,
        int N50,
        int L50,
        int N90,
        int L90,
        int LargestContig,
        int SmallestContig,
        double MeanLength,
        double MedianLength,
        double GcContent,
        int TotalGaps,
        long TotalGapLength,
        double GapPercentage);

    /// <summary>
    /// Nx/Lx statistics for various thresholds.
    /// </summary>
    public readonly record struct NxStatistics(
        int Threshold,
        int Nx,
        int Lx,
        long CumulativeLength);

    /// <summary>
    /// Gap information in assembly.
    /// </summary>
    public readonly record struct GapInfo(
        string SequenceId,
        int Start,
        int End,
        int Length,
        string GapType);

    /// <summary>
    /// Repeat annotation.
    /// </summary>
    public readonly record struct RepeatAnnotation(
        string SequenceId,
        int Start,
        int End,
        string RepeatClass,
        string RepeatFamily,
        double DivergencePercent,
        char Strand);

    /// <summary>
    /// BUSCO-like completeness result: the BUSCO classes Complete (= Single-copy + Duplicated),
    /// Fragmented and Missing over n = <see cref="TotalGenes"/> markers.
    /// </summary>
    /// <remarks>
    /// Every percentage is relative to n, as BUSCO <c>hmmer.py</c> <c>_get_busco_percentages</c>
    /// (<c>d_percent = multi_copy / total_buscos · 100</c>) and compleasm summary; percentages are not
    /// rounded (BUSCO rounds to one decimal only in <see cref="BuscoNotation"/>).
    /// </remarks>
    /// <param name="TotalGenes">n, the number of markers.</param>
    /// <param name="Complete">C = S + D.</param>
    /// <param name="CompleteSingleCopy">S.</param>
    /// <param name="CompleteDuplicated">D.</param>
    /// <param name="Fragmented">F.</param>
    /// <param name="Missing">M = n − C − F.</param>
    /// <param name="CompletenessPercent">C% = C / n · 100.</param>
    /// <param name="DuplicationPercent">BUSCO D% = D / n · 100 (before B10 F43 it was D / C · 100).
    /// <see cref="BuscoNotation"/> recomputes C% and D% from the counts.</param>
    public readonly record struct CompletenessResult(
        int TotalGenes,
        int Complete,
        int CompleteSingleCopy,
        int CompleteDuplicated,
        int Fragmented,
        int Missing,
        double CompletenessPercent,
        double DuplicationPercent)
    {
        /// <summary>BUSCO S% = S / n · 100 (0 when n = 0).</summary>
        public double SingleCopyPercent => Percent(CompleteSingleCopy);

        /// <summary>BUSCO F% = F / n · 100 (0 when n = 0).</summary>
        public double FragmentedPercent => Percent(Fragmented);

        /// <summary>BUSCO M% = M / n · 100 (0 when n = 0).</summary>
        public double MissingPercent => Percent(Missing);

        /// <summary>
        /// BUSCO one-line summary <c>C:x%[S:x%,D:x%],F:x%,M:x%,n:N</c> (BUSCO <c>hmmer.py</c>
        /// <c>record_results</c>; each percentage <c>round(·, 1)</c>, i.e. half-to-even at one decimal,
        /// printed in Python's shortest form, e.g. <c>100.0</c>, <c>33.3</c>).
        /// </summary>
        public string BuscoNotation => string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"C:{Round1(Percent(Complete))}%[S:{Round1(SingleCopyPercent)}%,D:{Round1(Percent(CompleteDuplicated))}%],F:{Round1(FragmentedPercent)}%,M:{Round1(MissingPercent)}%,n:{TotalGenes}");

        // BUSCO evaluates (count / total_buscos) * 100 in this order; (count * 100) / n can differ in the
        // last bit (199/2000: 9.950000000000001 vs 9.95 → rounds to 10.0 vs 9.9).
        private double Percent(int count) => TotalGenes > 0 ? count / (double)TotalGenes * 100 : 0;

        // Python round(x, 1) rounds the exact binary value half-to-even (e.g. round(0.05, 1) = 0.1
        // because the double 0.05 is slightly above 0.05). The exact decimal expansion (.NET Core 3+
        // formats doubles exactly) is rounded in decimal arithmetic to reproduce that.
        private static string Round1(double percent)
        {
            var invariant = System.Globalization.CultureInfo.InvariantCulture;
            decimal exact = decimal.Parse(percent.ToString("F25", invariant), invariant);
            return Math.Round(exact, 1, MidpointRounding.ToEven).ToString("0.0", invariant);
        }
    }

    /// <summary>
    /// One local alignment of a marker gene to an assembly sequence (see
    /// <see cref="FindMarkerHits(string, IEnumerable{ValueTuple{string, string}})"/>).
    /// </summary>
    /// <param name="SequenceId">Assembly sequence (contig/scaffold) id.</param>
    /// <param name="Strand">'+' when the marker aligns as given, '-' when its reverse complement aligns.</param>
    /// <param name="Start">0-based inclusive start of the aligned region on the assembly sequence (forward strand).</param>
    /// <param name="End">0-based inclusive end of the aligned region on the assembly sequence.</param>
    /// <param name="MarkerStart">0-based inclusive start of the aligned part, in the marker's own (forward) coordinates.</param>
    /// <param name="MarkerEnd">0-based inclusive end of the aligned part, in marker coordinates.</param>
    /// <param name="Identity">Identical columns / alignment columns (gap columns included; BLAST / EMBOSS identity), in [0, 1].</param>
    /// <param name="Coverage">Aligned marker span / marker length = (MarkerEnd − MarkerStart + 1) / length
    /// (compleasm <c>Protein_mapped_rate = (Protein_End − Protein_Start) / Protein_length</c>), in (0, 1].</param>
    public readonly record struct MarkerHit(
        string SequenceId,
        char Strand,
        int Start,
        int End,
        int MarkerStart,
        int MarkerEnd,
        double Identity,
        double Coverage);

    /// <summary>
    /// k-mer spectrum evaluation of a read set (error rate, genome size) and, when an assembly is
    /// given, Merqury k-mer completeness and consensus quality (QV).
    /// See <see cref="EstimateCompletenessFromKmers(IEnumerable{ValueTuple{string, int}}, IEnumerable{ValueTuple{string, string}}?, int, bool)"/>.
    /// </summary>
    /// <param name="K">k-mer length (length of every spectrum k-mer).</param>
    /// <param name="SolidThreshold">Merqury <c>filt</c>: the first local minimum of the histogram
    /// (<c>kmerHistToPloidyDepth</c> <c>effectiveMinDepth</c>); solid k-mers have count &gt; filt, 0 when the
    /// histogram has no local minimum.</param>
    /// <param name="PeakCoverage">Coverage used for the genome size: the histogram mode above
    /// <see cref="SolidThreshold"/> (smallest such count on ties), or the caller's expected coverage.</param>
    /// <param name="TotalKmers">Σ count over the spectrum (k-mer instances).</param>
    /// <param name="ErrorKmers">Σ count of the k-mers with count ≤ <see cref="SolidThreshold"/>.</param>
    /// <param name="ErrorRate">Per-base read error rate 1 − (1 − ErrorKmers/TotalKmers)^(1/k) (GenomeScope
    /// <c>uniroot</c> of Σ_{i=1..k} C(k,i)·x^i·(1−x)^(k−i)·total = error, whose closed form this is).</param>
    /// <param name="EstimatedGenomeSize">(TotalKmers − ErrorKmers) / PeakCoverage, rounded half-to-even
    /// (GenomeScope <c>total_len</c> with the observed peak; R <c>round</c>).</param>
    /// <param name="SolidKmers">Distinct spectrum k-mers with count &gt; filt (Merqury <c>read.gt$filt.meryl</c>).</param>
    /// <param name="AssemblySolidKmers">Distinct assembly k-mers that are solid read k-mers (Merqury <c>meryl intersect</c>).</param>
    /// <param name="Completeness">Merqury k-mer completeness AssemblySolidKmers / SolidKmers as a fraction
    /// (Merqury prints ×100); NaN when no assembly was given.</param>
    /// <param name="AssemblyKmers">Σ count of assembly k-mers (Merqury TOTAL, <c>meryl statistics</c> "present").</param>
    /// <param name="AssemblyOnlyKmers">Σ count of assembly k-mers absent from the read spectrum (Merqury ASM_ONLY).</param>
    /// <param name="AssemblyErrorRate">Merqury E = 1 − (1 − ASM_ONLY/TOTAL)^(1/k); NaN without assembly k-mers.</param>
    /// <param name="QualityValue">Merqury QV = −10·log10(E); +∞ when E = 0, NaN without assembly k-mers.</param>
    public readonly record struct KmerCompletenessEstimate(
        int K,
        int SolidThreshold,
        int PeakCoverage,
        long TotalKmers,
        long ErrorKmers,
        double ErrorRate,
        long EstimatedGenomeSize,
        long SolidKmers,
        long AssemblySolidKmers,
        double Completeness,
        long AssemblyKmers,
        long AssemblyOnlyKmers,
        double AssemblyErrorRate,
        double QualityValue);

    /// <summary>
    /// Scaffold gap structure.
    /// </summary>
    public readonly record struct ScaffoldStructure(
        string ScaffoldId,
        IReadOnlyList<(string ContigId, int Start, int End)> Contigs,
        IReadOnlyList<GapInfo> Gaps,
        int TotalLength,
        int ContigLength,
        int GapLength);

    /// <summary>
    /// K-mer comparison of two assemblies (Mash / sourmash model, strand-collapsed exact k-mer sets).
    /// </summary>
    /// <param name="Assembly1Name">Label of assembly 1.</param>
    /// <param name="Assembly2Name">Label of assembly 2.</param>
    /// <param name="KmerSize">k used for both sets.</param>
    /// <param name="DistinctKmers1">|K1|: distinct canonical k-mers of assembly 1 (all contigs pooled).</param>
    /// <param name="DistinctKmers2">|K2|: distinct canonical k-mers of assembly 2.</param>
    /// <param name="SharedKmers">|K1 ∩ K2|.</param>
    /// <param name="Containment1">|K1 ∩ K2| / |K1| — fraction of assembly 1's k-mers present in assembly 2 (0 when K1 is empty).</param>
    /// <param name="Containment2">|K1 ∩ K2| / |K2|.</param>
    /// <param name="JaccardIndex">|K1 ∩ K2| / |K1 ∪ K2| (0 when both sets are empty).</param>
    /// <param name="MashDistance">D = −(1/k)·ln(2J/(1+J)) (Ondov et al. 2016) with Mash's boundary rules (J = 1 → 0, J = 0 → 1).</param>
    public readonly record struct AssemblyComparison(
        string Assembly1Name,
        string Assembly2Name,
        int KmerSize,
        int DistinctKmers1,
        int DistinctKmers2,
        int SharedKmers,
        double Containment1,
        double Containment2,
        double JaccardIndex,
        double MashDistance);

    #endregion

    #region Basic Statistics

    /// <summary>
    /// Calculates comprehensive assembly statistics.
    /// </summary>
    /// <remarks>
    /// Definitions (QUAST <c>basic_stats.py</c> / <c>N50.py</c> / <c>qutils.py</c>, ablab/quast master):
    /// <list type="bullet">
    /// <item>TotalLength = Σ sequence lengths (every character, N included); N50/L50, N90/L90 =
    /// <see cref="CalculateNx(IReadOnlyList{int},long,int)"/> at x = 50 / 90 (QUAST <c>N50_and_L50</c>).</item>
    /// <item>GcContent = pooled (G + C) / (A + C + G + T) over all sequences as a fraction in [0, 1]
    /// (QUAST <c>GC_content</c>: <c>total_GC_amount / total_contig_length</c> with the ACGT length
    /// <c>len − count('N')</c> after <c>correct_seq</c> upper-cases and turns IUPAC codes into N).
    /// N, IUPAC ambiguity codes and every other non-ACGT symbol are excluded from both counts. 0 when there
    /// is no A/C/G/T.</item>
    /// <item>Gaps = maximal runs of gap symbols (QUAST <c>split_by_ns</c> with the threshold 1); TotalGapLength
    /// = number of gap symbols; TotalLengthNoGaps = TotalLength − TotalGapLength; GapPercentage = 100 ·
    /// TotalGapLength / TotalLength. A gap symbol is N/n and, with <paramref name="nonAcgtAsGap"/> (the
    /// default, QUAST), every other symbol outside A/C/G/T (IUPAC ambiguity codes, U, '-', …); see
    /// <see cref="IsGapSymbol"/>.</item>
    /// <item>MeanLength = TotalLength / count; MedianLength = sample median of the lengths (mean of the
    /// two central values for an even count, <see cref="StatisticsHelper.Median"/>).</item>
    /// </list>
    /// An empty input returns all zeros (QUAST reports these fields as undefined).
    /// </remarks>
    /// <param name="sequences">Assembled sequences.</param>
    /// <param name="nonAcgtAsGap">true (default, QUAST <c>correct_seq</c>): every symbol that is not A/C/G/T
    /// (either case) is an N, i.e. counts as a gap and not as a base; false: only N/n is a gap and the
    /// canonical counter's base set (IUPAC codes and U are bases for the length but outside the A/C/G/T
    /// denominator) applies.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A sequence is null.</exception>
    public static AssemblyStatistics CalculateStatistics(
        IEnumerable<(string Id, string Sequence)> sequences,
        bool nonAcgtAsGap = true)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var seqList = sequences.Select(s => RequireSequence(s, nameof(sequences))).ToList();

        if (seqList.Count == 0)
        {
            return new AssemblyStatistics(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        var lengths = seqList.Select(s => s.Length).OrderByDescending(l => l).ToList();
        long totalLength = lengths.Sum(l => (long)l);

        // Calculate N50/L50
        var n50Stats = CalculateNx(lengths, totalLength, 50);
        var n90Stats = CalculateNx(lengths, totalLength, 90);

        // GC content: pooled G+C over A+C+G+T (QUAST GC_content) via the canonical Core counter.
        long gcCount = 0;
        long acgtCount = 0;
        int totalGaps = 0;
        long totalGapLength = 0;
        foreach (var sequence in seqList)
        {
            var (gc, valid) = sequence.AsSpan().CountGcAndValidNucleotides();
            gcCount += gc;

            // Gap statistics: maximal gap-symbol runs (QUAST split_by_ns, threshold 1).
            long sequenceGapLength = 0;
            foreach (var (_, runLength) in EnumerateGapRuns(sequence, nonAcgtAsGap))
            {
                totalGaps++;
                sequenceGapLength += runLength;
            }
            totalGapLength += sequenceGapLength;

            // QUAST GC_content denominator: len − count('N') after correct_seq (every non-ACGT is an N).
            acgtCount += nonAcgtAsGap ? sequence.Length - sequenceGapLength : valid;
        }
        double gcContent = acgtCount > 0 ? gcCount / (double)acgtCount : 0;
        long totalNoGaps = totalLength - totalGapLength;
        double gapPercentage = totalLength > 0 ? totalGapLength * 100.0 / totalLength : 0;

        return new AssemblyStatistics(
            seqList.Count,
            totalLength,
            totalNoGaps,
            n50Stats.Nx,
            n50Stats.Lx,
            n90Stats.Nx,
            n90Stats.Lx,
            lengths[0],
            lengths[^1],
            totalLength / (double)seqList.Count,
            StatisticsHelper.Median(lengths.Select(l => (double)l).ToArray()),
            gcContent,
            totalGaps,
            totalGapLength,
            gapPercentage);
    }

    // Nx/Lx threshold percentage for the N50/L50 statistics.
    // Source: Miller, Koren & Sutton (2010), Genomics 95(6):315-327, §1.2 —
    // "the smallest contig ... whose combined length represents at least 50% of the assembly".
    private const int N50ThresholdPercent = 50;

    /// <summary>
    /// Calculates Nx and Lx statistics from contig/scaffold lengths.
    /// </summary>
    /// <remarks>
    /// Nx is the length of the shortest sequence such that sequences of that length or
    /// longer cover at least x% of the assembly; Lx is the count of those sequences.
    /// Source: Miller, Koren &amp; Sutton (2010), Genomics 95(6):315-327, §1.2. Delegates to the
    /// single canonical core <see cref="SequenceAssembler.CalculateNx(IReadOnlyList{int},long,int)"/>
    /// (port of QUAST <c>quast_libs/N50.py</c> <c>NG50_and_LG50</c>, inclusive "at least x%").
    /// Passing an (estimated) genome size as <paramref name="totalLength"/> gives NGx/LGx; when the
    /// lengths never reach x% of it (QUAST <c>None</c>), as for an empty list, Nx = Lx = 0.
    /// </remarks>
    /// <param name="sortedLengths">Non-negative sequence lengths sorted in non-increasing order.</param>
    /// <param name="totalLength">Sum of all sequence lengths (the assembly size), or a genome size for NGx.</param>
    /// <param name="threshold">Percentage threshold x in [0, 100] (50 for N50/L50).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> outside [0, 100] or negative <paramref name="totalLength"/>.</exception>
    /// <exception cref="ArgumentException">A negative length, or lengths not in non-increasing order.</exception>
    public static NxStatistics CalculateNx(
        IReadOnlyList<int> sortedLengths,
        long totalLength,
        int threshold)
    {
        var (nx, lx, cumulative) = SequenceAssembler.CalculateNx(sortedLengths, totalLength, threshold);
        return new NxStatistics(threshold, nx, lx, cumulative);
    }

    /// <summary>
    /// Calculates Nx and Lx statistics directly from an unsorted set of contig lengths.
    /// </summary>
    /// <remarks>
    /// Convenience overload that sorts the lengths descending and computes the assembly
    /// total internally, then delegates to <see cref="CalculateNx(IReadOnlyList{int},long,int)"/>.
    /// Source: Miller, Koren &amp; Sutton (2010), Genomics 95(6):315-327, §1.2.
    /// </remarks>
    /// <param name="lengths">Contig/scaffold lengths (any order).</param>
    /// <param name="threshold">Percentage threshold x in [0, 100] (50 for N50/L50).</param>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> outside [0, 100].</exception>
    /// <exception cref="ArgumentException">A length is negative.</exception>
    public static NxStatistics CalculateNx(IEnumerable<int> lengths, int threshold)
    {
        ArgumentNullException.ThrowIfNull(lengths);
        var sorted = lengths.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);
        return CalculateNx(sorted, total, threshold);
    }

    /// <summary>
    /// Calculates the N50 statistic: the length of the shortest contig in the smallest
    /// set of longest contigs whose combined length is at least 50% of the assembly.
    /// </summary>
    /// <remarks>
    /// Source: Miller, Koren &amp; Sutton (2010), Genomics 95(6):315-327, §1.2 —
    /// "The contig N50 is the length of the smallest contig in the set that contains the
    /// fewest (largest) contigs whose combined length represents at least 50% of the
    /// assembly." Returns 0 for an empty input.
    /// </remarks>
    /// <param name="lengths">Contig/scaffold lengths (any order).</param>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> is null.</exception>
    /// <exception cref="ArgumentException">A length is negative.</exception>
    public static int CalculateN50(IEnumerable<int> lengths)
    {
        return CalculateNx(lengths, N50ThresholdPercent).Nx;
    }

    /// <summary>
    /// Default Nx curve thresholds: the deciles 10, 20, ..., 90 (the percentages of the Nx plot,
    /// QUAST <c>Nx</c> plot x = 0..100 in steps; N50 and N90 are among them).
    /// </summary>
    private static readonly int[] DefaultNxCurveThresholds = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };

    /// <summary>
    /// Calculates Nx/Lx for multiple thresholds (the Nx curve).
    /// </summary>
    /// <remarks>
    /// One <see cref="NxStatistics"/> per entry of <paramref name="thresholds"/>, in ascending
    /// threshold order (a repeated threshold yields a repeated entry); an empty list means the
    /// deciles 10, 20, ..., 90. Each point is the canonical
    /// <see cref="CalculateNx(IReadOnlyList{int},long,int)"/> (QUAST <c>NG50_and_LG50</c>) over the
    /// assembly size. The arguments are validated and the curve is computed eagerly, so an invalid
    /// input throws at the call (not on first enumeration).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> or <paramref name="thresholds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is outside [0, 100].</exception>
    /// <exception cref="ArgumentException">A length is negative.</exception>
    public static IEnumerable<NxStatistics> CalculateNxCurve(
        IEnumerable<int> lengths,
        params int[] thresholds)
    {
        ArgumentNullException.ThrowIfNull(lengths);
        ArgumentNullException.ThrowIfNull(thresholds);

        var sorted = lengths.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);
        if (thresholds.Length == 0)
            thresholds = DefaultNxCurveThresholds;

        var curve = new List<NxStatistics>(thresholds.Length);
        foreach (int t in thresholds.OrderBy(t => t))
            curve.Add(CalculateNx(sorted, total, t));
        return curve;
    }

    /// <summary>
    /// Calculates auN (area under the Nx curve): auN = Σ Lᵢ² / Σ Lⱼ.
    /// </summary>
    /// <remarks>
    /// Li H (2020), "auN: a new metric to measure assembly contiguity"
    /// (lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity); QUAST <c>au_metric</c>
    /// (<c>quast_libs/N50.py</c>). Returns 0 for an empty or all-zero list (QUAST <c>None</c>).
    /// Σ Lᵢ² and Σ Lⱼ are accumulated exactly in integers (128-bit / 64-bit), so no precision is lost
    /// for assemblies above 2^53 squared bases; the single division is performed on the exact
    /// quotient and remainder (error ≤ 1 ulp of the result).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> is null.</exception>
    /// <exception cref="ArgumentException">A length is negative.</exception>
    public static double CalculateAuN(IEnumerable<int> lengths)
        => CalculateAuNCore(lengths, genomeSize: null);

    /// <summary>
    /// Calculates auNG, the area under the NGx curve: auNG = Σ Lᵢ² / G, where G is the (estimated) genome
    /// size instead of the assembly length (QUAST <c>au_metric(numlist, reference_length)</c>, which divides
    /// by <c>reference_length</c> when it is given; Li 2020, "auN"/"auNG" in the blog post and the
    /// <c>calN50</c> <c>-L</c> option). Σ Lᵢ² is accumulated exactly in 128-bit integers and divided once
    /// by the integer <paramref name="genomeSize"/> (quotient + remainder), as <see cref="CalculateAuN(IEnumerable{int})"/>.
    /// An empty list returns 0 (QUAST <c>None</c>). No clamp is applied when the assembly is longer than G, as in QUAST.
    /// </summary>
    /// <param name="lengths">Sequence lengths (non-negative).</param>
    /// <param name="genomeSize">Reference / estimated genome size G in bases (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> is null.</exception>
    /// <exception cref="ArgumentException">A length is negative.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genomeSize"/> ≤ 0.</exception>
    public static double CalculateAuN(IEnumerable<int> lengths, long genomeSize)
    {
        if (genomeSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(genomeSize), genomeSize,
                "The genome size must be positive.");
        return CalculateAuNCore(lengths, genomeSize);
    }

    private static double CalculateAuNCore(IEnumerable<int> lengths, long? genomeSize)
    {
        ArgumentNullException.ThrowIfNull(lengths);

        long total = 0;
        UInt128 sumOfSquares = UInt128.Zero;
        foreach (int length in lengths)
        {
            if (length < 0)
                throw new ArgumentException("Sequence lengths must be non-negative.", nameof(lengths));
            total += length;
            sumOfSquares += (UInt128)((ulong)length * (ulong)length);
        }

        // QUAST au_metric: empty list -> None; denominator = reference_length when given, else Σ L (all-zero -> None).
        long denominator = genomeSize ?? total;
        if (denominator == 0)
            return 0;

        // Exact integer quotient + remainder, divided once: Σ L² / denominator = q + r / denominator.
        UInt128 quotient = sumOfSquares / (ulong)denominator;
        UInt128 remainder = sumOfSquares % (ulong)denominator;
        return (double)quotient + (double)remainder / denominator;
    }

    #endregion

    #region Gap Analysis

    /// <summary>
    /// Default minimum N-run length that separates two contigs of a scaffold: NCBI treats "every run of
    /// 10 or more Ns" as an assembly gap when assembly statistics are calculated (GenBank "Gapped Format
    /// for Genome Submissions", table2asn <c>-gaps-min 10</c>), and QUAST breaks scaffolds into contigs at
    /// runs of ≥ <c>qconfig.Ns_break_threshold = 10</c> Ns (<c>qutils.broke_scaffolds</c>). Shorter runs
    /// are ambiguous bases inside a contig (AGP v2.1: "there may be short runs of Ns due to ambiguous
    /// base calls").
    /// </summary>
    public const int DefaultScaffoldGapMinLength = 10;

    /// <summary>
    /// True for a gap symbol: <c>N</c>/<c>n</c> and, when <paramref name="nonAcgtAsGap"/>, every other
    /// symbol that is not A/C/G/T (either case). QUAST <c>qutils.correct_seq</c> upper-cases the sequence,
    /// turns the IUPAC ambiguity codes M K R Y W S V B H D into N and rejects any other symbol outside
    /// A/C/G/T/N (the contig file is skipped); seqtk <c>gap</c> treats every non-ACGT as a gap. Here a
    /// symbol QUAST would reject (U, '-', '*', …) is also an N, the lenient reading shared with seqtk.
    /// </summary>
    private static bool IsGapSymbol(char c, bool nonAcgtAsGap)
        => nonAcgtAsGap
            ? SequenceCase.ToUpperAscii(c) is not ('A' or 'C' or 'G' or 'T')
            : c is 'N' or 'n';

    /// <summary>
    /// Enumerates the maximal runs of gap symbols (<see cref="IsGapSymbol"/>) in
    /// <paramref name="sequence"/> as (0-based start, length). The single gap-run scanner behind
    /// <see cref="FindGaps"/>, <see cref="AnalyzeScaffolds"/>, <see cref="ExtractContigs"/> and the gap
    /// fields of <see cref="CalculateStatistics"/>.
    /// </summary>
    private static IEnumerable<(int Start, int Length)> EnumerateGapRuns(string sequence, bool nonAcgtAsGap)
    {
        int i = 0;
        while (i < sequence.Length)
        {
            if (IsGapSymbol(sequence[i], nonAcgtAsGap))
            {
                int start = i;
                while (i < sequence.Length && IsGapSymbol(sequence[i], nonAcgtAsGap))
                    i++;
                yield return (start, i - start);
            }
            else
            {
                i++;
            }
        }
    }

    /// <summary>
    /// Splits a scaffold into its contigs, port of QUAST <c>qutils.split_by_ns</c>: the scaffold is cut
    /// at every maximal N-run of length ≥ <paramref name="minGapLength"/>; shorter N-runs stay inside the
    /// contig. Returns the non-empty contigs as half-open [Start, End) and the separating gaps.
    /// </summary>
    private static (List<(int Start, int End)> Contigs, List<(int Start, int Length)> Gaps)
        SplitAtGaps(string sequence, int minGapLength, bool nonAcgtAsGap)
    {
        var contigs = new List<(int Start, int End)>();
        var gaps = new List<(int Start, int Length)>();
        int contigStart = 0;
        foreach (var run in EnumerateGapRuns(sequence, nonAcgtAsGap))
        {
            if (run.Length < minGapLength)
                continue;
            if (run.Start > contigStart)
                contigs.Add((contigStart, run.Start));
            gaps.Add(run);
            contigStart = run.Start + run.Length;
        }
        if (sequence.Length > contigStart)
            contigs.Add((contigStart, sequence.Length));
        return (contigs, gaps);
    }

    private static string RequireSequence((string Id, string Sequence) entry, string paramName)
        => entry.Sequence ?? throw new ArgumentException(
            $"Sequence '{entry.Id}' is null.", paramName);

    private static void ValidateMinGapLength(int minGapLength)
    {
        if (minGapLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minGapLength), minGapLength,
                "A gap is a run of at least one N; the minimum gap length must be ≥ 1.");
    }

    /// <summary>
    /// Identifies all gaps (maximal runs of gap symbols: <c>N</c>/<c>n</c> and, by default, every other
    /// non-ACGT symbol) in assembly sequences.
    /// </summary>
    /// <remarks>
    /// One <see cref="GapInfo"/> per maximal N-run of length ≥ <paramref name="minGapLength"/>, in
    /// sequence order. Coordinates are 0-based with an <b>inclusive</b> <c>End</c>
    /// (<c>Length = End − Start + 1</c>); the AGP v2.1 equivalent is <c>[Start + 1, End + 1]</c>
    /// (1-based inclusive) and the BED / <c>seqtk gap</c> equivalent is <c>[Start, End + 1)</c>.
    /// With <paramref name="nonAcgtAsGap"/> (default, QUAST <c>qutils.correct_seq</c> + <c>split_by_ns</c>;
    /// <c>seqtk gap</c> likewise treats every non-ACGT symbol as a gap) IUPAC ambiguity codes and any other
    /// symbol outside A/C/G/T are gap symbols like N; with false only <c>N</c>/<c>n</c> is a gap and IUPAC
    /// codes are bases (AGP v2.1, NCBI gapped format). With the default 1 the gaps are
    /// exactly those counted by <see cref="CalculateStatistics"/> (same <paramref name="nonAcgtAsGap"/>); pass
    /// <see cref="DefaultScaffoldGapMinLength"/> (10) for the NCBI assembly-gap convention.
    /// <see cref="GapInfo.GapType"/> is a Seqeron length class (&lt; 10 "Short", &lt; 100 "Medium",
    /// &lt; 1000 "Long", else "Scaffold"); it is <b>not</b> the AGP column 7b gap type, which records
    /// linkage evidence that a sequence does not carry.
    /// </remarks>
    /// <param name="sequences">Assembled sequences.</param>
    /// <param name="minGapLength">Minimum gap length (≥ 1).</param>
    /// <param name="nonAcgtAsGap">true (default): every non-ACGT symbol is a gap symbol (QUAST); false: N/n only.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<GapInfo> FindGaps(
        IEnumerable<(string Id, string Sequence)> sequences,
        int minGapLength = 1,
        bool nonAcgtAsGap = true)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ValidateMinGapLength(minGapLength);
        return FindGapsIterator(sequences, minGapLength, nonAcgtAsGap);
    }

    private static IEnumerable<GapInfo> FindGapsIterator(
        IEnumerable<(string Id, string Sequence)> sequences, int minGapLength, bool nonAcgtAsGap)
    {
        foreach (var entry in sequences)
        {
            string sequence = RequireSequence(entry, nameof(sequences));
            foreach (var (start, length) in EnumerateGapRuns(sequence, nonAcgtAsGap))
            {
                if (length >= minGapLength)
                    yield return new GapInfo(entry.Id, start, start + length - 1, length, ClassifyGap(length));
            }
        }
    }

    /// <summary>
    /// Classifies gap by length (Seqeron length bins; not an AGP gap type).
    /// </summary>
    private static string ClassifyGap(int length)
    {
        return length switch
        {
            < 10 => "Short",
            < 100 => "Medium",
            < 1000 => "Long",
            _ => "Scaffold"
        };
    }

    /// <summary>
    /// Analyzes gap distribution in assembly.
    /// </summary>
    /// <remarks>
    /// Count, arithmetic mean, sample median (mean of the two central values for an even count, numpy
    /// <c>median</c> / <see cref="StatisticsHelper.Median"/>, as <see cref="CalculateStatistics"/>
    /// MedianLength), maximum gap length and the number of gaps per <see cref="GapInfo.GapType"/>.
    /// An empty list returns all zeros and an empty map.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="gaps"/> is null.</exception>
    public static (int Count, double MeanLength, double MedianLength, int MaxLength, IReadOnlyDictionary<string, int> TypeCounts)
        AnalyzeGapDistribution(IEnumerable<GapInfo> gaps)
    {
        ArgumentNullException.ThrowIfNull(gaps);
        var gapList = gaps.ToList();

        if (gapList.Count == 0)
            return (0, 0, 0, 0, new Dictionary<string, int>());

        var lengths = gapList.Select(g => (double)g.Length).ToArray();
        var typeCounts = gapList.GroupBy(g => g.GapType)
                                .ToDictionary(g => g.Key, g => g.Count());

        return (
            gapList.Count,
            lengths.Average(),
            StatisticsHelper.Median(lengths),
            gapList.Max(g => g.Length),
            typeCounts);
    }

    #endregion

    #region Scaffold Analysis

    /// <summary>
    /// Analyzes scaffold structure (contigs and gaps).
    /// </summary>
    /// <remarks>
    /// A scaffold is split into contigs at every maximal N-run of length ≥ <paramref name="minGapLength"/>
    /// (QUAST <c>qutils.split_by_ns</c>; default <see cref="DefaultScaffoldGapMinLength"/> = 10, the NCBI
    /// assembly-gap convention). Shorter N-runs are ambiguous bases inside a contig (AGP v2.1), so
    /// <c>ContigLength + GapLength = TotalLength</c> always holds. Contigs are named
    /// <c>{id}_contig{n}</c> (n from 1) and, like the gaps, use 0-based inclusive <c>[Start, End]</c>
    /// coordinates; a leading or trailing gap opens no empty contig. A gap symbol is N/n and, with
    /// <paramref name="nonAcgtAsGap"/> (default, QUAST <c>correct_seq</c>), every other non-ACGT symbol (see
    /// <see cref="FindGaps"/>).
    /// </remarks>
    /// <param name="scaffolds">Scaffolds.</param>
    /// <param name="minGapLength">Minimum gap-run length that splits contigs (≥ 1).</param>
    /// <param name="nonAcgtAsGap">true (default): every non-ACGT symbol is a gap symbol (QUAST); false: N/n only.</param>
    /// <exception cref="ArgumentNullException"><paramref name="scaffolds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<ScaffoldStructure> AnalyzeScaffolds(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minGapLength = DefaultScaffoldGapMinLength,
        bool nonAcgtAsGap = true)
    {
        ArgumentNullException.ThrowIfNull(scaffolds);
        ValidateMinGapLength(minGapLength);
        return AnalyzeScaffoldsIterator(scaffolds, minGapLength, nonAcgtAsGap);
    }

    private static IEnumerable<ScaffoldStructure> AnalyzeScaffoldsIterator(
        IEnumerable<(string Id, string Sequence)> scaffolds, int minGapLength, bool nonAcgtAsGap)
    {
        foreach (var entry in scaffolds)
        {
            string id = entry.Id;
            string sequence = RequireSequence(entry, nameof(scaffolds));
            var (pieces, runs) = SplitAtGaps(sequence, minGapLength, nonAcgtAsGap);

            var contigs = new List<(string, int, int)>(pieces.Count);
            for (int c = 0; c < pieces.Count; c++)
                contigs.Add(($"{id}_contig{c + 1}", pieces[c].Start, pieces[c].End - 1));

            var gaps = runs
                .Select(r => new GapInfo(id, r.Start, r.Start + r.Length - 1, r.Length, ClassifyGap(r.Length)))
                .ToList();

            yield return new ScaffoldStructure(
                id, contigs, gaps,
                sequence.Length,
                pieces.Sum(p => p.End - p.Start),
                runs.Sum(r => r.Length));
        }
    }

    /// <summary>
    /// Extracts contigs from scaffolds.
    /// </summary>
    /// <remarks>
    /// Port of QUAST "broken scaffolds" (<c>qutils.split_by_ns</c> followed by the
    /// <c>--min-contig</c> length filter of <c>correct_fasta</c>): each scaffold is cut at every maximal
    /// N-run of length ≥ <paramref name="minGapLength"/> (default <see cref="DefaultScaffoldGapMinLength"/>
    /// = 10, QUAST <c>Ns_break_threshold</c> and the NCBI assembly-gap convention; pass 1 to cut at every
    /// N); shorter N-runs stay inside the contig. Contigs of length ≥ <paramref name="minContigLength"/>
    /// are returned in scaffold order, named <c>{id}_contig{n}</c> where n is the contig's 1-based
    /// position among <b>all</b> contigs of its scaffold, dropped short ones included (QUAST
    /// <c>split_by_ns</c> numbers every piece and <c>--min-contig</c> filters afterwards), so a
    /// name is the same as the one <see cref="AnalyzeScaffolds"/> gives that contig with the same
    /// <paramref name="minGapLength"/>; gaps in the numbering mark dropped contigs. Gap symbols are N/n and,
    /// with <paramref name="nonAcgtAsGap"/> (default, QUAST <c>correct_seq</c>), every other non-ACGT symbol
    /// (see <see cref="FindGaps"/>); the returned pieces are slices of the input, so case and the original
    /// symbols of short runs kept inside a contig are preserved (QUAST writes them upper-cased as N).
    /// </remarks>
    /// <param name="scaffolds">Scaffolds.</param>
    /// <param name="minContigLength">Minimum contig length (≥ 0).</param>
    /// <param name="minGapLength">Minimum gap-run length that splits contigs (≥ 1).</param>
    /// <param name="nonAcgtAsGap">true (default): every non-ACGT symbol is a gap symbol (QUAST); false: N/n only.</param>
    /// <exception cref="ArgumentNullException"><paramref name="scaffolds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minContigLength"/> &lt; 0 or
    /// <paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<(string Id, string Sequence)> ExtractContigs(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minContigLength = 200,
        int minGapLength = DefaultScaffoldGapMinLength,
        bool nonAcgtAsGap = true)
    {
        ArgumentNullException.ThrowIfNull(scaffolds);
        if (minContigLength < 0)
            throw new ArgumentOutOfRangeException(nameof(minContigLength), minContigLength,
                "The minimum contig length cannot be negative.");
        ValidateMinGapLength(minGapLength);
        return ExtractContigsIterator(scaffolds, minContigLength, minGapLength, nonAcgtAsGap);
    }

    private static IEnumerable<(string Id, string Sequence)> ExtractContigsIterator(
        IEnumerable<(string Id, string Sequence)> scaffolds, int minContigLength, int minGapLength, bool nonAcgtAsGap)
    {
        foreach (var entry in scaffolds)
        {
            string sequence = RequireSequence(entry, nameof(scaffolds));
            var pieces = SplitAtGaps(sequence, minGapLength, nonAcgtAsGap).Contigs;
            for (int c = 0; c < pieces.Count; c++)
            {
                var (start, end) = pieces[c];
                if (end - start >= minContigLength)
                    yield return ($"{entry.Id}_contig{c + 1}", sequence[start..end]);
            }
        }
    }

    #endregion

    #region Completeness Assessment

    /// <summary>
    /// Seed word length of the marker search: BLAST+ <c>megablast</c> (the default <c>blastn</c> task)
    /// <c>word_size</c> 28. A marker shorter than this is seeded by its whole length (exact occurrence).
    /// </summary>
    public const int MarkerSeedWordSize = 28;

    /// <summary>
    /// compleasm <c>--min_length_percent</c> default (0.6): a hit is valid when its aligned fraction of the
    /// marker reaches this value OR its identity reaches the identity threshold.
    /// </summary>
    public const double DefaultMarkerMinLengthFraction = 0.6;

    /// <summary>
    /// Assesses assembly completeness with marker genes (BUSCO-like): every marker is classified as
    /// Complete single-copy, Complete duplicated, Fragmented or Missing.
    /// </summary>
    /// <remarks>
    /// <para><b>Hit search</b> (<see cref="FindMarkerHits(string, IEnumerable{ValueTuple{string, string}})"/>):
    /// each marker and its reverse complement are seeded with exact ACGT words of
    /// <see cref="MarkerSeedWordSize"/> and every seed cluster is extended by Smith-Waterman-Gotoh local
    /// alignment with the BLAST+ blastn scoring (<see cref="SequenceAligner.BlastDna"/>). A hit reports
    /// identity (identical / alignment columns) and coverage (aligned marker span / marker length).</para>
    /// <para><b>Classification</b> = compleasm <c>record_1st_gene_label</c> (Huang &amp; Li 2023, length-rate
    /// mode, v0.2.6): valid hits are those with coverage ≥ <paramref name="minLengthFraction"/> OR identity
    /// ≥ <paramref name="identityThreshold"/> (compleasm <c>pass_tids</c> filter); no valid hit → Missing;
    /// complete hits have coverage ≥ <paramref name="coverageThreshold"/> (<c>min_complete</c>, default 0.9);
    /// none → Fragmented; one → Single; several → Duplicated when they lie on different sequences or form
    /// more than one cluster of overlapping regions on one sequence (<c>get_region_clusters</c>), else Single.</para>
    /// <para><b>Percentages</b> are BUSCO's: every class over n (<see cref="CompletenessResult"/>).</para>
    /// <para><b>Not BUSCO itself:</b> BUSCO/compleasm search lineage protein HMM profiles (HMMER, with
    /// lineage score and length cut-offs) through gene predictors or miniprot; this method aligns the caller's
    /// nucleotide marker sequences directly, so frameshifts/introns are not modelled and no profile score
    /// cut-off exists (the 28-bp exact seed is the significance filter).</para>
    /// </remarks>
    /// <param name="assembly">Assembly sequences (id, sequence); case-insensitive.</param>
    /// <param name="markerGenes">Marker genes (id, nucleotide sequence); each entry is one marker.</param>
    /// <param name="identityThreshold">compleasm <c>min_identity</c> analogue in [0, 1] (default 0.9; compleasm's protein default is 0.4).</param>
    /// <param name="coverageThreshold">compleasm <c>min_complete</c> in [0, 1] (default 0.9).</param>
    /// <param name="minLengthFraction">compleasm <c>min_length_percent</c> in [0, 1] (default 0.6).</param>
    /// <exception cref="ArgumentNullException"><paramref name="assembly"/> or <paramref name="markerGenes"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is outside [0, 1] or NaN.</exception>
    /// <exception cref="ArgumentException">A sequence is null or a marker is empty.</exception>
    public static CompletenessResult AssessCompleteness(
        IEnumerable<(string Id, string Sequence)> assembly,
        IEnumerable<(string GeneId, string Sequence)> markerGenes,
        double identityThreshold = 0.9,
        double coverageThreshold = 0.9,
        double minLengthFraction = DefaultMarkerMinLengthFraction)
    {
        ArgumentNullException.ThrowIfNull(assembly);
        ArgumentNullException.ThrowIfNull(markerGenes);
        ValidateFraction(identityThreshold, nameof(identityThreshold));
        ValidateFraction(coverageThreshold, nameof(coverageThreshold));
        ValidateFraction(minLengthFraction, nameof(minLengthFraction));

        var assemblySeqs = MaterializeAssembly(assembly, nameof(assembly));
        var markers = markerGenes
            .Select(g => RequireMarker(g.GeneId, g.Sequence, nameof(markerGenes)))
            .ToList();

        if (markers.Count == 0)
            return new CompletenessResult(0, 0, 0, 0, 0, 0, 0, 0);

        int singleCopy = 0, duplicated = 0, fragmented = 0, missing = 0;
        foreach (string marker in markers)
        {
            var hits = FindMarkerHitsCore(marker, assemblySeqs);
            switch (ClassifyMarker(hits, identityThreshold, coverageThreshold, minLengthFraction))
            {
                case MarkerClass.Single: singleCopy++; break;
                case MarkerClass.Duplicated: duplicated++; break;
                case MarkerClass.Fragmented: fragmented++; break;
                default: missing++; break;
            }
        }

        int n = markers.Count;
        int complete = singleCopy + duplicated;
        return new CompletenessResult(
            n,
            complete,
            singleCopy,
            duplicated,
            fragmented,
            missing,
            complete / (double)n * 100,
            duplicated / (double)n * 100);
    }

    /// <summary>
    /// Finds the local alignments of one marker gene to the assembly, on both strands.
    /// </summary>
    /// <remarks>
    /// Seeds: every ACGT-only word of length w = min(<see cref="MarkerSeedWordSize"/>, marker length) shared
    /// by the marker (or its reverse complement) and a sequence gives a diagonal d = sequence position −
    /// marker position. Sorted diagonals are split where consecutive values differ by more than half the
    /// marker length (two non-overlapping copies differ by at least the marker length), and each cluster
    /// [d_min, d_max] is aligned against the window [d_min − w, d_max + length + w) of the sequence with
    /// <see cref="SequenceAligner.LocalAlignAffine(string, string, ScoringMatrix?)"/> and
    /// <see cref="SequenceAligner.BlastDna"/> (+2/−3, gap 5 + 2·k). Identical hits from overlapping windows
    /// are reported once. Coordinates are 0-based inclusive; reverse-strand hits keep forward-strand
    /// sequence coordinates and forward marker coordinates.
    /// </remarks>
    /// <param name="markerSequence">Marker nucleotide sequence (non-empty; case-insensitive).</param>
    /// <param name="assembly">Assembly sequences (id, sequence).</param>
    /// <returns>Hits ordered by strand ('+' first), input sequence order, then start.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="markerSequence"/> or <paramref name="assembly"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="markerSequence"/> is empty or a sequence is null.</exception>
    public static IReadOnlyList<MarkerHit> FindMarkerHits(
        string markerSequence,
        IEnumerable<(string Id, string Sequence)> assembly)
    {
        ArgumentNullException.ThrowIfNull(markerSequence);
        ArgumentNullException.ThrowIfNull(assembly);
        string marker = RequireMarker("marker", markerSequence, nameof(markerSequence));
        return FindMarkerHitsCore(marker, MaterializeAssembly(assembly, nameof(assembly)));
    }

    private enum MarkerClass { Missing, Fragmented, Single, Duplicated }

    /// <summary>compleasm v0.2.6 <c>pass_tids</c> filter + <c>record_1st_gene_label(by_length=False)</c>.</summary>
    private static MarkerClass ClassifyMarker(
        List<MarkerHit> hits, double minIdentity, double minComplete, double minLengthFraction)
    {
        var valid = hits.Where(h => h.Coverage >= minLengthFraction || h.Identity >= minIdentity).ToList();
        if (valid.Count == 0)
            return MarkerClass.Missing;

        var complete = valid.Where(h => h.Coverage >= minComplete).ToList();
        if (complete.Count == 0)
            return MarkerClass.Fragmented;
        if (complete.Count == 1)
            return MarkerClass.Single;
        if (complete.Select(h => h.SequenceId).Distinct(StringComparer.Ordinal).Count() > 1)
            return MarkerClass.Duplicated;

        // get_region_clusters on PAF half-open (start, stop): overlap when last.start <= start <= last.stop.
        int clusters = 0;
        long clusterStop = long.MinValue;
        foreach (var hit in complete.OrderBy(h => h.Start))
        {
            long stop = hit.End + 1L;
            if (clusters > 0 && hit.Start <= clusterStop)
            {
                clusterStop = Math.Max(clusterStop, stop);
            }
            else
            {
                clusters++;
                clusterStop = stop;
            }
        }
        return clusters == 1 ? MarkerClass.Single : MarkerClass.Duplicated;
    }

    private static List<MarkerHit> FindMarkerHitsCore(string marker, List<(string Id, string Sequence)> assembly)
    {
        int length = marker.Length;
        int word = Math.Min(MarkerSeedWordSize, length);
        var hits = new List<MarkerHit>();

        foreach (char strand in "+-")
        {
            string query = strand == '+' ? marker : DnaSequence.GetReverseComplementString(marker);
            var seedIndex = BuildSeedIndex(query, word);
            if (seedIndex.Count == 0)
                continue;

            foreach (var (id, sequence) in assembly)
            {
                var diagonals = CollectSeedDiagonals(sequence, word, seedIndex);
                if (diagonals.Count == 0)
                    continue;
                diagonals.Sort();

                var seen = new HashSet<(int, int, int, int)>();
                int clusterStart = 0;
                for (int i = 1; i <= diagonals.Count; i++)
                {
                    if (i < diagonals.Count && diagonals[i] - diagonals[i - 1] <= length / 2)
                        continue;

                    int windowStart = Math.Max(0, diagonals[clusterStart] - word);
                    int windowEnd = (int)Math.Min(sequence.Length, (long)diagonals[i - 1] + length + word);
                    clusterStart = i;
                    if (windowEnd <= windowStart)
                        continue;

                    var alignment = SequenceAligner.LocalAlignAffine(
                        query, sequence.Substring(windowStart, windowEnd - windowStart), SequenceAligner.BlastDna);
                    if (alignment.Score <= 0)
                        continue;

                    int queryStart = alignment.StartPosition1;
                    int queryEnd = alignment.EndPosition1;
                    int start = windowStart + alignment.StartPosition2;
                    int end = windowStart + alignment.EndPosition2;
                    if (!seen.Add((queryStart, queryEnd, start, end)))
                        continue;

                    double identity = SequenceAligner.CalculateStatistics(alignment, SequenceAligner.BlastDna).Identity / 100.0;
                    int markerStart = strand == '+' ? queryStart : length - 1 - queryEnd;
                    int markerEnd = strand == '+' ? queryEnd : length - 1 - queryStart;
                    hits.Add(new MarkerHit(id, strand, start, end, markerStart, markerEnd, identity,
                        (queryEnd - queryStart + 1) / (double)length));
                }
            }
        }

        return hits;
    }

    /// <summary>2-bit code of an upper-case base (A=0, C=1, G=2, T=3), −1 for anything else.</summary>
    private static int BaseCode(char c) => c switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1,
    };

    /// <summary>
    /// Enumerates (position, 2-bit word) for every ACGT-only window of length <paramref name="word"/>
    /// (≤ 31) of an upper-case sequence; a non-ACGT base restarts the window (Jellyfish/BLAST seeding rule).
    /// This is a positional seed scan, not k-mer counting, so it does not route through KmerAnalyzer.
    /// </summary>
    private static IEnumerable<(int Position, ulong Word)> EnumerateSeedWords(string sequence, int word)
    {
        ulong mask = word == 32 ? ulong.MaxValue : (1UL << (2 * word)) - 1;
        ulong code = 0;
        int filled = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            int b = BaseCode(sequence[i]);
            if (b < 0)
            {
                filled = 0;
                code = 0;
                continue;
            }
            code = ((code << 2) | (uint)b) & mask;
            if (++filled >= word)
                yield return (i - word + 1, code);
        }
    }

    private static Dictionary<ulong, List<int>> BuildSeedIndex(string query, int word)
    {
        var index = new Dictionary<ulong, List<int>>();
        foreach (var (position, code) in EnumerateSeedWords(query, word))
        {
            if (!index.TryGetValue(code, out var positions))
                index[code] = positions = new List<int>();
            positions.Add(position);
        }
        return index;
    }

    private static List<int> CollectSeedDiagonals(string sequence, int word, Dictionary<ulong, List<int>> seedIndex)
    {
        var diagonals = new List<int>();
        foreach (var (position, code) in EnumerateSeedWords(sequence, word))
        {
            if (seedIndex.TryGetValue(code, out var queryPositions))
            {
                foreach (int q in queryPositions)
                    diagonals.Add(position - q);
            }
        }
        return diagonals;
    }

    private static List<(string Id, string Sequence)> MaterializeAssembly(
        IEnumerable<(string Id, string Sequence)> assembly, string paramName)
        => assembly.Select(s => (s.Id, SequenceCase.ToUpperAscii(RequireSequence(s, paramName)))).ToList();

    private static string RequireMarker(string id, string? sequence, string paramName)
    {
        if (sequence is null)
            throw new ArgumentException($"Marker '{id}' sequence is null.", paramName);
        if (sequence.Length == 0)
            throw new ArgumentException($"Marker '{id}' sequence is empty.", paramName);
        return SequenceCase.ToUpperAscii(sequence);
    }

    private static void ValidateFraction(double value, string paramName)
    {
        if (!(value >= 0 && value <= 1))
            throw new ArgumentOutOfRangeException(paramName, value, "The threshold must lie in [0, 1].");
    }

    /// <summary>
    /// Estimates read error rate, genome size and (spectrum only) k-mer statistics from a read k-mer count
    /// spectrum. Without an assembly there is no k-mer completeness: <c>Completeness</c> is NaN.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="EstimateCompletenessFromKmers(IEnumerable{ValueTuple{string, int}}, IEnumerable{ValueTuple{string, string}}?, int, bool)"/>
    /// with no assembly; see there for the definitions (Merqury solid threshold, GenomeScope error rate,
    /// k-mer genome size).
    /// </remarks>
    /// <param name="kmerSpectrum">Read k-mer counts (each k-mer once, all of the same length k, count ≥ 1).</param>
    /// <param name="expectedCoverage">Peak coverage to use; 0 detects the histogram peak.</param>
    public static (double Completeness, double ErrorRate, long EstimatedGenomeSize)
        EstimateCompletenessFromKmers(
            IEnumerable<(string Kmer, int Count)> kmerSpectrum,
            int expectedCoverage = 0)
    {
        var estimate = EstimateCompletenessFromKmers(kmerSpectrum, null, expectedCoverage, canonical: true);
        return (estimate.Completeness, estimate.ErrorRate, estimate.EstimatedGenomeSize);
    }

    /// <summary>
    /// k-mer evaluation of a read spectrum and (optionally) an assembly: read error rate, genome size,
    /// Merqury k-mer completeness and Merqury QV.
    /// </summary>
    /// <remarks>
    /// <para><b>Histogram</b>: h(c) = number of distinct spectrum k-mers with count c, listed for the counts
    /// present in ascending order (<c>meryl histogram</c>).</para>
    /// <para><b>Solid threshold</b> (Merqury <c>build/filt.sh</c> → <c>kmerHistToPloidyDepth.jar</c>
    /// <c>effectiveMinDepth</c>): walking the histogram lines with slope d = h_i − h_{i−1} (the first line's
    /// slope is h_0 − Double.MAX_VALUE), the first line where the slope turns from negative to positive gives
    /// filt = c_i − 1; filt = 0 when no such line exists. Solid k-mers have count &gt; filt.</para>
    /// <para><b>Read error rate</b>: E = Σ_{c ≤ filt} c·h(c), T = Σ c·h(c); error = 1 − (1 − E/T)^(1/k)
    /// (GenomeScope 1.0 <c>error_rate</c> root of Σ_{i=1..k} C(k,i)x^i(1−x)^(k−i)·T = E; Vurture et al. 2017).</para>
    /// <para><b>Genome size</b>: round((T − E) / peak) (GenomeScope <c>total_len = (total_kmers −
    /// total_error_kmers)/(2·kcov)</c>; the homozygous peak 2·kcov is taken as the observed histogram mode
    /// above filt, or <paramref name="expectedCoverage"/>). The GenomeScope negative-binomial mixture fit
    /// (heterozygosity, duplication) is not performed: for a heterozygous genome whose highest peak is the
    /// heterozygous one, pass the homozygous coverage as <paramref name="expectedCoverage"/>.</para>
    /// <para><b>Assembly</b> (Merqury <c>eval/spectra-cn.sh</c>, Rhie et al. 2020): the assembly k-mers are
    /// counted per sequence by <see cref="KmerAnalyzer.CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress{double}?)"/>
    /// (ACGT-only windows; canonical when <paramref name="canonical"/>, as meryl). Completeness =
    /// |distinct assembly k-mers ∩ solid read k-mers| / |solid read k-mers|; QV = −10·log10(1 − (1 −
    /// ASM_ONLY/TOTAL)^(1/k)), where TOTAL counts all assembly k-mer instances and ASM_ONLY those absent
    /// from the read spectrum (all read k-mers, not only solid ones, as <c>meryl difference asm reads</c>).</para>
    /// </remarks>
    /// <param name="readKmerSpectrum">Read k-mer counts: each k-mer once, all of length k, count ≥ 1. With
    /// <paramref name="canonical"/> and an assembly the k-mers must be ACGT and canonical (min of k-mer and
    /// reverse complement, as <c>meryl</c> / <c>jellyfish -C</c> output).</param>
    /// <param name="assembly">Assembly sequences, or null for a spectrum-only estimate.</param>
    /// <param name="expectedCoverage">Peak coverage to use; 0 detects the histogram peak.</param>
    /// <param name="canonical">Count assembly k-mers canonically (meryl default).</param>
    /// <exception cref="ArgumentNullException"><paramref name="readKmerSpectrum"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="expectedCoverage"/> is negative.</exception>
    /// <exception cref="ArgumentException">The spectrum is empty, has a null/empty or duplicate k-mer, mixed
    /// k-mer lengths or a count &lt; 1; with an assembly and <paramref name="canonical"/>, a k-mer is not ACGT or
    /// not canonical; an assembly sequence is null.</exception>
    public static KmerCompletenessEstimate EstimateCompletenessFromKmers(
        IEnumerable<(string Kmer, int Count)> readKmerSpectrum,
        IEnumerable<(string Id, string Sequence)>? assembly,
        int expectedCoverage = 0,
        bool canonical = true)
    {
        ArgumentNullException.ThrowIfNull(readKmerSpectrum);
        if (expectedCoverage < 0)
            throw new ArgumentOutOfRangeException(nameof(expectedCoverage), expectedCoverage,
                "Expected coverage cannot be negative (0 = detect the peak).");

        var spectrum = new Dictionary<string, int>(StringComparer.Ordinal);
        int k = -1;
        foreach (var (kmer, count) in readKmerSpectrum)
        {
            if (string.IsNullOrEmpty(kmer))
                throw new ArgumentException("A spectrum k-mer is null or empty.", nameof(readKmerSpectrum));
            if (k < 0)
                k = kmer.Length;
            else if (kmer.Length != k)
                throw new ArgumentException(
                    $"All spectrum k-mers must have the same length (found {k} and {kmer.Length}).", nameof(readKmerSpectrum));
            if (count < 1)
                throw new ArgumentException($"K-mer '{kmer}' has count {count}; counts must be ≥ 1.", nameof(readKmerSpectrum));
            if (!spectrum.TryAdd(SequenceCase.ToUpperAscii(kmer), count))
                throw new ArgumentException($"K-mer '{kmer}' is listed more than once.", nameof(readKmerSpectrum));
        }
        if (spectrum.Count == 0)
            throw new ArgumentException("The k-mer spectrum is empty.", nameof(readKmerSpectrum));

        var histogram = new SortedDictionary<int, long>();
        foreach (int count in spectrum.Values)
            histogram[count] = histogram.TryGetValue(count, out long h) ? h + 1 : 1;

        int filt = MerqurySolidThreshold(histogram);
        long total = 0, errorKmers = 0, solidKmers = 0;
        int peak = 0;
        long peakHeight = 0;
        foreach (var (count, distinct) in histogram)
        {
            total += count * distinct;
            if (count <= filt)
            {
                errorKmers += count * distinct;
                continue;
            }
            solidKmers += distinct;
            if (distinct > peakHeight)
            {
                peakHeight = distinct;
                peak = count;
            }
        }
        if (expectedCoverage > 0)
            peak = expectedCoverage;

        double errorRate = KmerErrorRate(errorKmers, total, k);
        long genomeSize = (long)Math.Round((total - errorKmers) / (double)peak, MidpointRounding.ToEven);

        if (assembly is null)
            return new KmerCompletenessEstimate(k, filt, peak, total, errorKmers, errorRate, genomeSize,
                solidKmers, 0, double.NaN, 0, 0, double.NaN, double.NaN);

        if (canonical)
        {
            foreach (string kmer in spectrum.Keys)
            {
                if (kmer.Any(c => BaseCode(c) < 0))
                    throw new ArgumentException(
                        $"K-mer '{kmer}' is not ACGT; canonical evaluation needs meryl-style ACGT k-mers.", nameof(readKmerSpectrum));
                if (string.CompareOrdinal(DnaSequence.GetReverseComplementString(kmer), kmer) < 0)
                    throw new ArgumentException(
                        $"K-mer '{kmer}' is not canonical; pass a canonical spectrum (meryl / jellyfish -C) or canonical = false.",
                        nameof(readKmerSpectrum));
            }
        }

        var options = new KmerCountingOptions(Canonical: canonical, AcgtOnly: true);
        var assemblyKmers = new Dictionary<string, long>(StringComparer.Ordinal);
        foreach (var entry in assembly)
        {
            string sequence = RequireSequence(entry, nameof(assembly));
            if (sequence.Length < k)
                continue;
            foreach (var (kmer, count) in KmerAnalyzer.CountKmers(sequence, k, options))
                assemblyKmers[kmer] = assemblyKmers.TryGetValue(kmer, out long c) ? c + count : count;
        }

        long assemblyTotal = 0, assemblyOnly = 0, assemblySolid = 0;
        foreach (var (kmer, count) in assemblyKmers)
        {
            assemblyTotal += count;
            if (!spectrum.TryGetValue(kmer, out int readCount))
                assemblyOnly += count;
            else if (readCount > filt)
                assemblySolid++;
        }

        double completeness = assemblySolid / (double)solidKmers;
        double assemblyErrorRate = assemblyTotal == 0 ? double.NaN : KmerErrorRate(assemblyOnly, assemblyTotal, k);
        double qv = double.IsNaN(assemblyErrorRate) ? double.NaN : -10.0 * Math.Log10(assemblyErrorRate);

        return new KmerCompletenessEstimate(k, filt, peak, total, errorKmers, errorRate, genomeSize,
            solidKmers, assemblySolid, completeness, assemblyTotal, assemblyOnly, assemblyErrorRate, qv);
    }

    /// <summary>
    /// Per-base error rate from the fraction of erroneous k-mers: 1 − (1 − errorKmers/totalKmers)^(1/k)
    /// (Merqury <c>qv.sh</c> ERROR; closed form of GenomeScope's <c>uniroot</c>).
    /// </summary>
    private static double KmerErrorRate(long errorKmers, long totalKmers, int k)
        => 1.0 - Math.Pow(1.0 - errorKmers / (double)totalKmers, 1.0 / k);

    /// <summary>
    /// Merqury solid-k-mer threshold: <c>kmerHistToPloidyDepth</c> <c>effectiveMinDepth</c> (the value
    /// <c>filt.sh</c> reads from line 2 of the ploidy table), transcribed from HistToPloidyDepth.java
    /// <c>getPloidyDepthTable</c>: first histogram line whose slope d turns positive after a negative
    /// previous slope (d·prevD &lt; 0 and d &gt; prevD) gives depth − 1.
    /// </summary>
    private static int MerqurySolidThreshold(SortedDictionary<int, long> histogram)
    {
        bool first = true;
        double prevCount = double.MaxValue;
        double prevD = 0;
        foreach (var (depth, distinct) in histogram)
        {
            double count = distinct;
            if (first)
            {
                prevD = count - prevCount;
                prevCount = count;
                first = false;
                continue;
            }
            double d = count - prevCount;
            if (d * prevD < 0 && d > prevD)
                return depth - 1;
            prevCount = count;
            prevD = d;
        }
        return 0;
    }

    #endregion

    #region Repeat Analysis

    /// <summary>
    /// Identifies repetitive regions with the k-mer-count masking model of NCBI WindowMasker
    /// (Morgulis et al. 2006, Bioinformatics 22:134; NCBI C++ Toolkit <c>src/algo/winmask</c>:
    /// <c>seq_masker.cpp</c> <c>CSeqMasker::DoMask</c>, <c>seq_masker_window.cpp</c>,
    /// <c>seq_masker_score_mean.cpp</c>, <c>win_mask_gen_counts.cpp</c>).
    /// </summary>
    /// <remarks>
    /// <para><b>Counting.</b> Every k-mer (unit) of all input sequences is counted over the whole input
    /// ("genome"). A window containing a non-ACGT symbol contributes no unit (WindowMasker <c>ambig</c> reset;
    /// <see cref="KmerAnalyzer.CountKmers(string,int,KmerCountingOptions,System.Threading.CancellationToken,System.IProgress{double}?)"/>
    /// with <c>AcgtOnly</c>, the Jellyfish window rule). With <paramref name="bothStrands"/> (default, as
    /// WindowMasker) a unit and its reverse complement share one count: count(u) = F(u) + F(rc(u)), which is
    /// WindowMasker's canonical counting (a self-complementary unit counts its forward occurrences twice, as
    /// <c>win_mask_gen_counts.cpp</c> increments both <c>unit</c> and <c>runit</c> when they are equal).</para>
    /// <para><b>Scoring.</b> A window of <paramref name="windowSize"/> bases (step 1; only all-ACGT windows) has
    /// W − k + 1 units; its score is the integer mean of their counts (<c>CSeqMaskerScoreMean</c>: <c>sum / num</c>).</para>
    /// <para><b>Masking</b> (the first pass of <c>DoMask</c>, copied branch by branch): a window with score ≥
    /// <paramref name="minCopies"/> (T_threshold) starts or extends the masked interval; a window with score below
    /// <paramref name="extendCopies"/> (T_extend) closes it unless it begins inside the interval; windows scoring in
    /// [T_extend, T_threshold) keep an open interval connected without enlarging it. Intervals are 0-based,
    /// inclusive, [Start, End] = first base of the first triggering window … last base of the last one.
    /// <c>Copies</c> is the largest unit count inside the interval.</para>
    /// <para>Deviations from WindowMasker, declared: T_threshold / T_extend are given as counts (WindowMasker derives
    /// them from percentiles of the genome-wide unit-count distribution); the optional merge pass, discontiguous
    /// units and the count-database score cut-offs (T_low, T_high) are not implemented.</para>
    /// </remarks>
    /// <param name="sequences">Named sequences (case-insensitive); counting is pooled over all of them.</param>
    /// <param name="kmerSize">Unit size k (≥ 1; default 15).</param>
    /// <param name="minCopies">T_threshold: window score that triggers masking (≥ 1; default 3).</param>
    /// <param name="windowSize">Window length in bases (≥ k; default 100).</param>
    /// <param name="bothStrands">Count units on both strands (WindowMasker default).</param>
    /// <param name="extendCopies">T_extend (1..minCopies); null = <paramref name="minCopies"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A size or threshold is out of range.</exception>
    public static IEnumerable<(string SequenceId, int Start, int End, int Copies)> FindRepetitiveRegions(
        IEnumerable<(string Id, string Sequence)> sequences,
        int kmerSize = 15,
        int minCopies = 3,
        int windowSize = 100,
        bool bothStrands = true,
        int? extendCopies = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(kmerSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(minCopies, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, kmerSize);
        int textend = extendCopies ?? minCopies;
        ArgumentOutOfRangeException.ThrowIfLessThan(textend, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(textend, minCopies);

        return FindRepetitiveRegionsCore(
            sequences.Select(s => (s.Id, Sequence: SequenceCase.ToUpperAscii(s.Sequence ?? string.Empty))).ToList(),
            kmerSize, minCopies, windowSize, bothStrands, textend);
    }

    private static IEnumerable<(string SequenceId, int Start, int End, int Copies)> FindRepetitiveRegionsCore(
        List<(string Id, string Sequence)> sequences,
        int k, int cutoff, int windowSize, bool bothStrands, int textend)
    {
        var forward = new Dictionary<string, int>(StringComparer.Ordinal);
        var acgt = new KmerCountingOptions(AcgtOnly: true);
        foreach (var (_, sequence) in sequences)
        {
            foreach (var (kmer, n) in KmerAnalyzer.CountKmers(sequence, k, acgt))
                forward[kmer] = forward.GetValueOrDefault(kmer) + n;
        }

        int numUnits = windowSize - k + 1;
        foreach (var (id, seq) in sequences)
        {
            int n = seq.Length;
            if (n < windowSize) continue;

            // unit[j] = count of the k-mer starting at j (-1 when it contains a non-ACGT symbol).
            int[] unit = UnitCounts(seq, k, forward, bothStrands);

            // Prefix sums of the counts and of the invalid-unit flags give each window's mean in O(1).
            var sum = new long[unit.Length + 1];
            var bad = new int[unit.Length + 1];
            for (int j = 0; j < unit.Length; j++)
            {
                sum[j + 1] = sum[j] + Math.Max(unit[j], 0);
                bad[j + 1] = bad[j] + (unit[j] < 0 ? 1 : 0);
            }

            int start = 0, end = 0, cend = 0;
            for (int s = 0; s + windowSize <= n; s++)
            {
                // A window is usable only if all of its units are valid (all W bases are A/C/G/T).
                if (bad[s + numUnits] - bad[s] != 0) continue;

                int e = s + windowSize - 1;
                long score = (sum[s + numUnits] - sum[s]) / numUnits;

                if (score < textend)
                {
                    if (end > start && s > cend)
                    {
                        yield return (id, start, end, MaxUnit(unit, start, end, k));
                        start = end = cend = 0;
                    }
                }
                else if (score < cutoff)
                {
                    if (end > start)
                    {
                        if (s > cend + 1)
                        {
                            yield return (id, start, end, MaxUnit(unit, start, end, k));
                            start = end = cend = 0;
                        }
                        else
                        {
                            cend = e;
                        }
                    }
                }
                else
                {
                    if (end > start)
                    {
                        if (s > cend + 1)
                        {
                            yield return (id, start, end, MaxUnit(unit, start, end, k));
                            start = s;
                        }
                    }
                    else
                    {
                        start = s;
                    }

                    cend = end = e;
                }
            }

            if (end > start)
                yield return (id, start, end, MaxUnit(unit, start, end, k));
        }
    }

    /// <summary>
    /// Count of every k-mer start of an upper-case sequence (-1 for a k-mer with a non-ACGT symbol), looked up
    /// by span (no per-position string); with <paramref name="bothStrands"/> the reverse-complement count is added.
    /// </summary>
    private static int[] UnitCounts(string seq, int k, Dictionary<string, int> forward, bool bothStrands)
    {
        var lookup = forward.GetAlternateLookup<ReadOnlySpan<char>>();
        var unit = new int[seq.Length - k + 1];
        char[] rc = new char[k];
        int j = 0;
        foreach (ReadOnlySpan<char> kmer in seq.AsSpan().EnumerateKmers(k))
        {
            if (kmer.IsValidDna())
            {
                int c = lookup.TryGetValue(kmer, out int f) ? f : 0;
                if (bothStrands)
                {
                    if (!kmer.TryGetReverseComplement(rc))
                        throw new InvalidOperationException("Reverse-complement buffer too small.");
                    c += lookup.TryGetValue(rc, out int r) ? r : 0;
                }

                unit[j] = c;
            }
            else
            {
                unit[j] = -1;
            }

            j++;
        }

        return unit;
    }

    private static int MaxUnit(int[] unit, int start, int end, int k)
    {
        int max = 0;
        for (int j = start; j <= end - k + 1; j++)
            max = Math.Max(max, unit[j]);
        return max;
    }

    /// <summary>
    /// Finds approximate (interrupted) tandem repeats with the Tandem Repeats Finder model
    /// (Benson 1999, Nucleic Acids Res 27:573): delegates to
    /// <see cref="RepeatFinder.FindApproximateTandemRepeats(string,TandemRepeatsFinderParameters,int)"/> with TRF's
    /// recommended parameters (Match 2, Mismatch 7, Delta 7, PM 80, PI 10).
    /// </summary>
    /// <remarks>
    /// <c>Unit</c> is TRF's consensus pattern; <c>Copies</c> TRF's copy number (aligned consensus columns /
    /// consensus size; fractional); <c>Purity</c> is TRF's "percent matches between adjacent copies" as a fraction
    /// (0–1); <c>Start</c>/<c>End</c> are 0-based inclusive. A repeat is reported when its alignment score is
    /// ≥ <paramref name="minScore"/> (TRF Minscore, 50 recommended), its period lies in
    /// [<paramref name="minUnitLength"/>, <paramref name="maxUnitLength"/>] (TRF MaxPeriod ≤ 2000) and its copy number is
    /// ≥ <paramref name="minCopies"/>. Non-ACGT symbols never match.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A bound is out of range (period 1..2000, copies and score ≥ 1).</exception>
    public static IEnumerable<(string SequenceId, int Start, int End, string Unit, double Copies, double Purity)>
        FindTandemRepeats(
            IEnumerable<(string Id, string Sequence)> sequences,
            int minUnitLength = 2,
            int maxUnitLength = 50,
            int minCopies = 3,
            int minScore = RepeatFinder.DefaultApproximateMinScore)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(minUnitLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxUnitLength, minUnitLength);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(maxUnitLength, RepeatFinder.MaxApproximatePeriod);
        ArgumentOutOfRangeException.ThrowIfLessThan(minCopies, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(minScore, 1);

        var parameters = new TandemRepeatsFinderParameters { MaxPeriod = maxUnitLength, MinScore = minScore };
        return FindTandemRepeatsCore(sequences.ToList(), parameters, minUnitLength, minCopies);
    }

    private static IEnumerable<(string SequenceId, int Start, int End, string Unit, double Copies, double Purity)>
        FindTandemRepeatsCore(
            List<(string Id, string Sequence)> sequences,
            TandemRepeatsFinderParameters parameters,
            int minUnitLength,
            int minCopies)
    {
        foreach (var (id, sequence) in sequences)
        {
            foreach (var r in RepeatFinder.FindApproximateTandemRepeats(sequence, parameters, minUnitLength))
            {
                if (r.CopyNumber < minCopies) continue;
                yield return (id, r.Start, r.Start + r.SpanLength - 1, r.Consensus, r.CopyNumber, r.PercentMatches / 100.0);
            }
        }
    }

    /// <summary>
    /// Calculates repeat content as RepeatMasker's <c>ProcessRepeats</c> table does: every base is counted once
    /// ("Define length of coverage based on possibly overlapping fragments"), so overlapping annotations on the
    /// same sequence are not double counted and the class lengths sum to the total masked length.
    /// </summary>
    /// <remarks>
    /// Annotations are swept per sequence in (Start, input order); an annotation is credited with the bases beyond
    /// the furthest end seen so far, to its own <c>RepeatClass</c> (a fully contained annotation adds nothing).
    /// <c>RepeatPercentage</c> = masked bases / <paramref name="genomeLength"/> · 100 (RepeatMasker
    /// "bases masked: n bp (p %)"). Coordinates are 0-based inclusive. An empty annotation list or a genome length of
    /// 0 yields all zeros.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="repeats"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genomeLength"/> &lt; 0, or the masked bases exceed it.</exception>
    /// <exception cref="ArgumentException">An annotation has Start &lt; 0 or End &lt; Start.</exception>
    public static (long TotalRepeatLength, double RepeatPercentage, IReadOnlyDictionary<string, long> RepeatClassLengths)
        CalculateRepeatContent(
            IEnumerable<RepeatAnnotation> repeats,
            long genomeLength)
    {
        ArgumentNullException.ThrowIfNull(repeats);
        ArgumentOutOfRangeException.ThrowIfLessThan(genomeLength, 0);
        var repeatList = repeats.ToList();
        foreach (var r in repeatList)
        {
            if (r.Start < 0 || r.End < r.Start)
                throw new ArgumentException(
                    $"Invalid repeat annotation interval [{r.Start}, {r.End}] on '{r.SequenceId}'.", nameof(repeats));
        }

        if (repeatList.Count == 0 || genomeLength == 0)
            return (0, 0, new Dictionary<string, long>());

        long total = 0;
        var classLengths = new Dictionary<string, long>();
        foreach (var perSequence in repeatList
                     .Select((r, i) => (r, i))
                     .GroupBy(x => x.r.SequenceId))
        {
            int furthestEnd = -1;
            foreach (var (r, _) in perSequence.OrderBy(x => x.r.Start).ThenBy(x => x.i))
            {
                int from = Math.Max(r.Start, furthestEnd + 1);
                long added = r.End >= from ? r.End - from + 1 : 0;
                if (added > 0)
                {
                    total += added;
                    classLengths[r.RepeatClass] = classLengths.GetValueOrDefault(r.RepeatClass) + added;
                }

                furthestEnd = Math.Max(furthestEnd, r.End);
            }
        }

        // Zero-length contributions still list their class (every annotated class appears, as before).
        foreach (var r in repeatList)
            classLengths.TryAdd(r.RepeatClass, 0);

        if (total > genomeLength)
            throw new ArgumentOutOfRangeException(nameof(genomeLength),
                $"The annotations cover {total} bp, more than the genome length {genomeLength}.");

        return (total, total * 100.0 / genomeLength, classLengths);
    }

    #endregion

    #region Assembly Comparison

    /// <summary>
    /// Compares two assemblies by their exact, strand-collapsed k-mer sets (Mash / sourmash model).
    /// </summary>
    /// <remarks>
    /// All contigs of an assembly are pooled into one set (a multi-FASTA is one genome for <c>mash sketch</c>);
    /// k-mers never span two contigs. K-mers are canonical, min(w, RC(w)), and any window containing a non-ACGT
    /// symbol is skipped — Mash <c>Sketch.cpp</c> <c>addMinHashes</c> rule, applied through
    /// <see cref="KmerAnalyzer.DistinctKmers(string, int, KmerCountingOptions)"/> (<c>Canonical = true</c>).
    /// <list type="bullet">
    /// <item>Jaccard J = |K1 ∩ K2| / |K1 ∪ K2| (Ondov et al. 2016, Genome Biol 17:132, eq. 1);</item>
    /// <item>Mash distance D = −(1/k)·ln(2J/(1+J)) (eq. 4) via <see cref="KmerAnalyzer.MashDistanceFromJaccard"/>;
    /// identical sets (also both empty) → 0, no shared k-mer → 1 (Mash <c>CommandDistance.cpp</c>);</item>
    /// <item>containment C(Ki, Kj) = |K1 ∩ K2| / |Ki| (Koslicki &amp; Zabeti 2019; sourmash <c>--containment</c>).</item>
    /// </list>
    /// This is a whole-assembly alignment-free similarity, not an alignment: it reports no aligned bases,
    /// identity percentage or breakpoints (those need nucmer/dnadiff alignments — see <see cref="FindSyntenicBlocks"/>
    /// for the MUM-cluster blocks).
    /// </remarks>
    /// <exception cref="ArgumentNullException">An assembly is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> is not positive.</exception>
    public static AssemblyComparison CompareAssemblies(
        IEnumerable<(string Id, string Sequence)> assembly1,
        IEnumerable<(string Id, string Sequence)> assembly2,
        string name1 = "Assembly1",
        string name2 = "Assembly2",
        int kmerSize = 21)
    {
        ArgumentNullException.ThrowIfNull(assembly1);
        ArgumentNullException.ThrowIfNull(assembly2);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(kmerSize);

        var kmers1 = PooledCanonicalKmers(assembly1, kmerSize);
        var kmers2 = PooledCanonicalKmers(assembly2, kmerSize);

        var (small, large) = kmers1.Count <= kmers2.Count ? (kmers1, kmers2) : (kmers2, kmers1);
        int shared = small.Count(large.Contains);
        int union = kmers1.Count + kmers2.Count - shared;

        double jaccard = union == 0 ? 0.0 : (double)shared / union;
        // Mash CommandDistance: shared == union (identical sets, incl. both empty) → 0.
        double mash = shared == union ? 0.0 : KmerAnalyzer.MashDistanceFromJaccard(jaccard, kmerSize);

        return new AssemblyComparison(
            name1, name2, kmerSize,
            kmers1.Count, kmers2.Count, shared,
            kmers1.Count > 0 ? (double)shared / kmers1.Count : 0.0,
            kmers2.Count > 0 ? (double)shared / kmers2.Count : 0.0,
            jaccard, mash);
    }

    /// <summary>Distinct canonical ACGT-only k-mers of every contig, pooled (canonical counter, no private loop).</summary>
    private static HashSet<string> PooledCanonicalKmers(IEnumerable<(string Id, string Sequence)> sequences, int k)
    {
        var options = new KmerCountingOptions(Canonical: true, AcgtOnly: true);
        var pooled = new HashSet<string>(StringComparer.Ordinal);
        foreach (var (_, sequence) in sequences)
            pooled.UnionWith(KmerAnalyzer.DistinctKmers(sequence ?? string.Empty, k, options));
        return pooled;
    }

    /// <summary>nucmer defaults (MUMmer 3/4 usage text): <c>-g</c> maxgap.</summary>
    private const int NucmerMaxGap = 90;

    /// <summary>nucmer default <c>--diagfactor</c> (<c>-d</c>): clustering diagonal-difference separation factor.</summary>
    private const double NucmerDiagonalFactor = 0.12;

    /// <summary>nucmer default <c>--diagdiff</c> (<c>-D</c>): fixed diagonal difference allowed between adjacent matches.</summary>
    private const int NucmerDiagonalDifference = 5;

    /// <summary>
    /// Finds syntenic (collinear) blocks between two assemblies as MUMmer's <c>nucmer</c> does before
    /// extension: maximal unique matches (MUMs) are anchors, and anchors are clustered by the <c>mgaps</c>
    /// algorithm (Delcher et al. 2002, NAR 30:2478, <i>Fast algorithms for large-scale genome alignment and
    /// comparison</i>; Marçais et al. 2018, PLoS Comput Biol 14:e1005944).
    /// </summary>
    /// <remarks>
    /// <para>Assembly 1 is the reference, assembly 2 the query (as <c>nucmer ref query</c>). For every pair
    /// (reference sequence, query sequence) and both query strands the MUMs of length ≥ <paramref name="kmerSize"/>
    /// are found with the library suffix tree (<c>FindMaximalUniqueMatches</c>; nucmer default
    /// <c>--mumreference</c> = unique in the reference, <paramref name="uniqueInBoth"/> = <c>--mum</c>); symbols
    /// other than A, C, G, T never match (case-insensitive). Anchors are then clustered with the <c>mgaps</c>
    /// routine of MUMmer 4 (<c>src/tigr/mgaps.cc</c>, <c>include/mummer/mgaps.hh</c>): repeat-internal anchors are
    /// filtered (<c>Filter_Matches</c>), anchors whose gap to the previous one is at most
    /// <paramref name="maxGap"/> and whose diagonals differ by at most max(<paramref name="diagonalDifference"/>,
    /// <paramref name="diagonalFactor"/>·gap) are linked (union-find), and from each linked component the
    /// highest-scoring colinear chain is extracted repeatedly (<c>Process_Cluster</c>: score = Σ length −
    /// overlap − diagonal drift). A chain whose total anchor length is at least <paramref name="minBlockSize"/>
    /// (nucmer <c>-c</c>) is reported.</para>
    /// <para>Coordinates are 0-based inclusive on the forward strand of each sequence; blocks are ordered by
    /// reference sequence (input order), Start1, then Start2. <c>IsInverted</c> is true when the block matches the
    /// reverse-complement strand of the query. Unlike nucmer no match extension (Smith–Waterman between anchors)
    /// is done: a block spans first anchor start to last anchor end.</para>
    /// </remarks>
    /// <param name="assembly1">Reference sequences.</param>
    /// <param name="assembly2">Query sequences.</param>
    /// <param name="minBlockSize">Minimum summed anchor length of a block (nucmer <c>-c</c>; mgaps <c>-l</c>).</param>
    /// <param name="kmerSize">Minimum MUM length (nucmer <c>-l</c> minmatch; default 21 as before, nucmer 20).</param>
    /// <param name="maxGap">Maximum query gap between linked anchors (nucmer <c>-g</c>, 90).</param>
    /// <param name="diagonalFactor">Diagonal separation factor (nucmer <c>-d</c>, 0.12).</param>
    /// <param name="diagonalDifference">Fixed diagonal difference (nucmer <c>-D</c>, 5).</param>
    /// <param name="uniqueInBoth">MUMs unique in both sequences (<c>--mum</c>) instead of the reference only.</param>
    public static IEnumerable<(string Seq1, int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        FindSyntenicBlocks(
            IEnumerable<(string Id, string Sequence)> assembly1,
            IEnumerable<(string Id, string Sequence)> assembly2,
            int minBlockSize = 1000,
            int kmerSize = 21,
            int maxGap = NucmerMaxGap,
            double diagonalFactor = NucmerDiagonalFactor,
            int diagonalDifference = NucmerDiagonalDifference,
            bool uniqueInBoth = false)
    {
        ArgumentNullException.ThrowIfNull(assembly1);
        ArgumentNullException.ThrowIfNull(assembly2);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(minBlockSize);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(kmerSize);
        ArgumentOutOfRangeException.ThrowIfNegative(maxGap);
        ArgumentOutOfRangeException.ThrowIfNegative(diagonalDifference);
        if (!(diagonalFactor >= 0))
            throw new ArgumentOutOfRangeException(nameof(diagonalFactor), diagonalFactor, "Must be non-negative.");

        var uniqueness = uniqueInBoth ? global::SuffixTree.MumUniqueness.Both : global::SuffixTree.MumUniqueness.Reference;
        var queries = assembly2.Select(q => (q.Id, Forward: NormalizeForMatching(q.Sequence, '$'))).ToList();
        return FindSyntenicBlocksCore(assembly1, queries, minBlockSize, kmerSize, maxGap, diagonalFactor,
            diagonalDifference, uniqueness);
    }

    private static IEnumerable<(string Seq1, int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        FindSyntenicBlocksCore(
            IEnumerable<(string Id, string Sequence)> assembly1,
            List<(string Id, string Forward)> queries,
            int minBlockSize, int minMatch, int maxGap, double diagonalFactor, int diagonalDifference,
            global::SuffixTree.MumUniqueness uniqueness)
    {
        var clusterer = new MgapsClusterer(diagonalDifference, maxGap, minBlockSize, diagonalFactor);

        // Reverse complement of the (already normalised ACGT/'$') query; '$' passes through the complement.
        var reverse = queries.Select(q => string.Create(q.Forward.Length, q.Forward, static (span, f) =>
        {
            f.AsSpan().TryGetReverseComplement(span);
        })).ToList();

        foreach (var (refId, refSeq) in assembly1)
        {
            string reference = NormalizeForMatching(refSeq, '#');
            if (reference.Length < minMatch)
                continue;

            var tree = global::SuffixTree.SuffixTree.Build(reference);
            var found = new List<(int Start1, int End1, int Q, int Start2, int End2, bool Inverted)>();

            for (int q = 0; q < queries.Count; q++)
            {
                int n = queries[q].Forward.Length;
                for (int strand = 0; strand < 2; strand++)
                {
                    string query = strand == 0 ? queries[q].Forward : reverse[q];
                    var mums = tree.FindMaximalUniqueMatches(query, minMatch, uniqueness);
                    if (mums.Count == 0)
                        continue;

                    foreach (var chain in clusterer.Cluster(mums.Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToList()))
                    {
                        int s1 = chain.Min(m => m.Start1), e1 = chain.Max(m => m.Start1 + m.Len) - 1;
                        int s2 = chain.Min(m => m.Start2), e2 = chain.Max(m => m.Start2 + m.Len) - 1;
                        // Reverse strand: positions are on RC(query); map back to the forward strand.
                        if (strand == 1)
                            (s2, e2) = (n - 1 - e2, n - 1 - s2);
                        found.Add((s1, e1, q, s2, e2, strand == 1));
                    }
                }
            }

            foreach (var b in found.OrderBy(f => f.Start1).ThenBy(f => f.Start2).ThenBy(f => f.Q).ThenBy(f => f.Inverted))
                yield return (refId, b.Start1, b.End1, queries[b.Q].Id, b.Start2, b.End2, b.Inverted);
        }
    }

    /// <summary>Upper-cases and replaces every non-ACGT symbol by <paramref name="mismatchSymbol"/> (never matched; distinct per side).</summary>
    private static string NormalizeForMatching(string? sequence, char mismatchSymbol)
    {
        if (string.IsNullOrEmpty(sequence))
            return string.Empty;
        return string.Create(sequence.Length, (sequence, mismatchSymbol), static (span, state) =>
        {
            for (int i = 0; i < span.Length; i++)
            {
                char c = SequenceCase.ToUpperAscii(state.sequence[i]);
                span[i] = c is 'A' or 'C' or 'G' or 'T' ? c : state.mismatchSymbol;
            }
        });
    }

    /// <summary>
    /// Literal port of MUMmer 4 <c>mgaps</c> (<c>ClusterMatches::Cluster_each</c>, <c>Filter_Matches</c>,
    /// <c>Process_Cluster</c>; <c>Use_Extents</c> = false, as <c>mgaps</c> without <c>-C</c>).
    /// Start1 = reference position, Start2 = query position.
    /// </summary>
    private sealed class MgapsClusterer
    {
        internal sealed class Match
        {
            public long Start1, Start2, Len;
            public long Score, From, Adj;
            public int ClusterId;
            public bool Good = true, Tentative;
        }

        private readonly int _fixedSeparation;
        private readonly long _maxSeparation;
        private readonly long _minOutputScore;
        private readonly double _separationFactor;

        public MgapsClusterer(int fixedSeparation, long maxSeparation, long minOutputScore, double separationFactor)
        {
            _fixedSeparation = fixedSeparation;
            _maxSeparation = maxSeparation;
            _minOutputScore = minOutputScore;
            _separationFactor = separationFactor;
        }

        /// <summary>Clusters (ref, query, length) matches; returns the chains (each ordered as mgaps prints it).</summary>
        public List<List<(int Start1, int Start2, int Len)>> Cluster(List<(int Ref, int Query, int Length)> matches)
        {
            var result = new List<List<(int, int, int)>>();
            int n = matches.Count;
            if (n == 0)
                return result;

            // 1-based, as the C++ (index 0 unused).
            var a = new Match[n + 1];
            for (int i = 0; i < n; i++)
                a[i + 1] = new Match { Start1 = matches[i].Ref, Start2 = matches[i].Query, Len = matches[i].Length };

            var uf = new int[n + 1];
            Array.Fill(uf, -1);

            Array.Sort(a, 1, n, Comparer<Match>.Create(ByStart2));
            n = FilterMatches(a, n);

            for (int i = 1; i < n; i++)
            {
                long iEnd = a[i].Start2 + a[i].Len;
                long iDiag = a[i].Start2 - a[i].Start1;
                for (int j = i + 1; j <= n; j++)
                {
                    long sep = a[j].Start2 - iEnd;
                    if (sep > _maxSeparation)
                        break;
                    long diagDiff = Math.Abs(a[j].Start2 - a[j].Start1 - iDiag);
                    if (diagDiff <= Math.Max(_fixedSeparation, (long)(_separationFactor * sep)))
                        Union(uf, Find(uf, i), Find(uf, j));
                }
            }

            for (int i = 1; i <= n; i++)
            {
                a[i].ClusterId = Find(uf, i);
                a[i].Good = false;
            }

            Array.Sort(a, 1, n, Comparer<Match>.Create((x, y) =>
                x.ClusterId != y.ClusterId ? x.ClusterId.CompareTo(y.ClusterId) : ByStart2(x, y)));

            for (int i = 1; i <= n;)
            {
                int j = i + 1;
                while (j <= n && a[i].ClusterId == a[j].ClusterId)
                    j++;
                ProcessCluster(new List<Match>(new ArraySegment<Match>(a, i, j - i)), result);
                i = j;
            }

            return result;
        }

        private static int ByStart2(Match x, Match y)
            => x.Start2 != y.Start2 ? x.Start2.CompareTo(y.Start2) : x.Start1.CompareTo(y.Start1);

        private static int Find(int[] uf, int a)
        {
            if (uf[a] < 0)
                return a;
            int i = a;
            while (uf[i] > 0)
                i = uf[i];
            for (int j = a, k; uf[j] != i; j = k)
            {
                k = uf[j];
                uf[j] = i;
            }
            return i;
        }

        private static void Union(int[] uf, int a, int b)
        {
            if (a == b)
                return;
            if (uf[a] < uf[b])
            {
                uf[a] += uf[b];
                uf[b] = a;
            }
            else
            {
                uf[b] += uf[a];
                uf[a] = b;
            }
        }

        /// <summary><c>Filter_Matches</c>: drops repeat-internal matches, merges overlapping same-diagonal matches; returns the new count.</summary>
        private static int FilterMatches(Match[] a, int n)
        {
            // a[1..n]; the C++ works on A+1 with 0-based indices.
            for (int i = 1; i < n; i++)
            {
                if (!a[i].Good)
                    continue;
                long iDiag = a[i].Start2 - a[i].Start1;
                long iEnd = a[i].Start2 + a[i].Len;

                for (int j = i + 1; j <= n && a[j].Start2 <= iEnd; j++)
                {
                    if (!a[j].Good)
                        continue;
                    long jDiag = a[j].Start2 - a[j].Start1;
                    if (iDiag == jDiag)
                    {
                        long jExtent = a[j].Len + a[j].Start2 - a[i].Start2;
                        if (jExtent > a[i].Len)
                        {
                            a[i].Len = jExtent;
                            iEnd = a[i].Start2 + jExtent;
                        }
                        a[j].Good = false;
                    }
                    else if (a[i].Start1 == a[j].Start1)
                    {
                        long olap = a[i].Start2 + a[i].Len - a[j].Start2;
                        if (ResolveTie(a[i], a[j], olap))
                            break;
                    }
                    else if (a[i].Start2 == a[j].Start2)
                    {
                        long olap = a[i].Start1 + a[i].Len - a[j].Start1;
                        if (ResolveTie(a[i], a[j], olap))
                            break;
                    }
                }
            }

            int w = 1;
            for (int i = 1; i <= n; i++)
                if (a[i].Good)
                    a[w++] = a[i];
            return w - 1;
        }

        /// <summary>The shared body of the two repeat-resolution branches; true = <c>break</c> out of the j loop.</summary>
        private static bool ResolveTie(Match i, Match j, long olap)
        {
            if (i.Len < j.Len)
            {
                if (olap >= i.Len / 2)
                {
                    i.Good = false;
                    return true;
                }
            }
            else if (j.Len < i.Len)
            {
                if (olap >= j.Len / 2)
                    j.Good = false;
            }
            else if (olap >= i.Len / 2)
            {
                j.Tentative = true;
                if (i.Tentative)
                {
                    i.Good = false;
                    return true;
                }
            }
            return false;
        }

        /// <summary><c>Process_Cluster</c>: repeatedly extract the best chain of the (sorted) component.</summary>
        private void ProcessCluster(List<Match> a, List<List<(int, int, int)>> output)
        {
            while (a.Count > 0)
            {
                int n = a.Count;
                for (int i = 0; i < n; i++)
                {
                    a[i].Score = a[i].Len;
                    a[i].Adj = 0;
                    a[i].From = -1;
                    for (int j = 0; j < i; j++)
                    {
                        long olap1 = a[j].Start1 + a[j].Len - a[i].Start1;
                        long olap2 = a[j].Start2 + a[j].Len - a[i].Start2;
                        long olap = Math.Max(Math.Max(0, olap1), olap2);
                        long pen = olap + Math.Abs(a[i].Start2 - a[i].Start1 - (a[j].Start2 - a[j].Start1));
                        if (a[j].Score + a[i].Len - pen > a[i].Score)
                        {
                            a[i].From = j;
                            a[i].Score = a[j].Score + a[i].Len - pen;
                            a[i].Adj = olap;
                        }
                    }
                }

                int best = 0;
                for (int i = 1; i < n; i++)
                    if (a[i].Score > a[best].Score)
                        best = i;

                long total = 0;
                for (int i = best; i >= 0; i = (int)a[i].From)
                {
                    a[i].Good = true;
                    total += a[i].Len;
                }

                if (total >= _minOutputScore)
                    output.Add(a.Where(m => m.Good).Select(m => ((int)m.Start1, (int)m.Start2, (int)m.Len)).ToList());

                a = a.Where(m => !m.Good).ToList();
            }
        }
    }

    #endregion

    #region Quality Assessment

    /// <summary>
    /// Sliding-window sequence-composition metrics of every sequence: GC content, N count and linguistic
    /// complexity per window (a window-level screen; QUAST reports GC as a whole-assembly value and a
    /// GC-content distribution over fixed windows).
    /// </summary>
    /// <remarks>
    /// Windows start at 0, <paramref name="stepSize"/>, 2·<paramref name="stepSize"/>, … for every start &lt;
    /// length; the last window(s) are truncated at the sequence end (Biopython <c>GC_skew</c> slicing
    /// <c>seq[i:i+window]</c>), so a sequence shorter than the window yields one truncated window and no base
    /// is left uncovered. <c>Position</c> is the 0-based start of the window (not its centre) and
    /// <c>WindowSize</c> its actual (possibly truncated) length. <c>GcContent</c> is the fraction in [0,1]
    /// from the canonical <see cref="GcSkewCalculator.CalculateWindowedGcContent(string,int,int,bool,bool)"/>
    /// (G+C over A+C+G+T+U, the QUAST <c>GC_content</c> rule of ASSEMBLY-STATS-001: N and IUPAC codes are
    /// not in the denominator; 0 for a window without an A/C/G/T/U base). <c>NCount</c> counts N/n.
    /// <c>Complexity</c> is the canonical <see cref="SequenceComplexity.CalculateLinguisticComplexity(string,int)"/>
    /// of the window (m = 10; Troyanskaya et al. 2002).
    /// </remarks>
    /// <param name="sequences">(id, sequence) pairs.</param>
    /// <param name="windowSize">Window length (≥ 1; default 1000).</param>
    /// <param name="stepSize">Distance between window starts (≥ 1); default max(1, windowSize / 2),
    /// i.e. 50 % overlap (a convention, not a published value).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is below 1 (windowSize 1 used to make the step 0 and never end).</exception>
    public static IEnumerable<(string SequenceId, int Position, int WindowSize, double GcContent, int NCount, double Complexity)>
        CalculateLocalQuality(
            IEnumerable<(string Id, string Sequence)> sequences,
            int windowSize = 1000,
            int? stepSize = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        int step = stepSize ?? Math.Max(1, windowSize / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(step, 1);
        return CalculateLocalQualityIterator(sequences, windowSize, step);
    }

    private static IEnumerable<(string SequenceId, int Position, int WindowSize, double GcContent, int NCount, double Complexity)>
        CalculateLocalQualityIterator(IEnumerable<(string Id, string Sequence)> sequences, int windowSize, int step)
    {
        foreach (var entry in sequences)
        {
            string id = entry.Id, sequence = RequireSequence(entry, nameof(sequences));
            foreach (var w in EnumerateCompositionWindows(sequence, windowSize, step))
            {
                string window = sequence.Substring(w.Start, w.Length);
                yield return (id, w.Start, w.Length, w.Gc, w.NCount,
                    SequenceComplexity.CalculateLinguisticComplexity(SequenceCase.ToUpperAscii(window)));
            }
        }
    }

    // The single window driver of the local-quality screens: GC from the canonical windowed GC profile
    // (includePartialWindow), N count and the number of A/C/G/T/U bases of each window.
    private static IEnumerable<(int Start, int Length, double Gc, int NCount, int ValidBases)>
        EnumerateCompositionWindows(string sequence, int windowSize, int step)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var p in GcSkewCalculator.CalculateWindowedGcContent(
            sequence, windowSize, step, fraction: true, includePartialWindow: true))
        {
            int length = p.WindowEnd - p.WindowStart + 1;
            var span = sequence.AsSpan(p.WindowStart, length);
            int n = 0;
            foreach (char c in span)
                if (c == 'N' || c == 'n')
                    n++;
            yield return (p.WindowStart, length, p.GcContent, n, span.CountGcAndValidNucleotides().ValidCount);
        }
    }

    /// <summary>
    /// A region flagged by <see cref="FindSuspiciousRegions"/>: inclusive 0-based [Start, End] and the
    /// measured value of every criterion that fired (0 for a criterion that did not).
    /// </summary>
    /// <param name="SequenceId">Sequence id.</param>
    /// <param name="Start">0-based first base.</param>
    /// <param name="End">0-based last base (inclusive).</param>
    /// <param name="Reason">"; "-joined criteria in the fixed order High N content, GC deviation, Low complexity.</param>
    /// <param name="MaxGcDeviation">Largest |window GC − pooled GC| of the GC-flagged windows (fraction).</param>
    /// <param name="MaxNFraction">Largest N fraction of the N-flagged windows.</param>
    /// <param name="MaxDustScore">Largest DUST score (Morgulis 2006) of the low-complexity intervals.</param>
    public readonly record struct SuspiciousRegion(
        string SequenceId, int Start, int End, string Reason,
        double MaxGcDeviation, double MaxNFraction, double MaxDustScore);

    /// <summary>
    /// Sequence-only screen for regions of an assembly that look anomalous. It is <b>not</b> a misassembly
    /// caller: QUAST/GAGE misassemblies need an alignment to a reference, and no published threshold defines
    /// "suspicious" for GC or N content, so the cut-offs are caller-set parameters and no composite score is
    /// reported — each region carries the measured values of the criteria that fired.
    /// </summary>
    /// <remarks>
    /// Criteria (union, then overlapping/adjacent regions of any criterion are merged into one):
    /// <list type="bullet">
    /// <item><b>High N content</b> — a window whose N fraction exceeds <paramref name="maxNFraction"/>
    /// (gaps are N runs in AGP/QUAST; use <see cref="FindGaps"/> for exact gap intervals).</item>
    /// <item><b>GC deviation</b> — a window whose GC fraction differs from the pooled assembly GC by more than
    /// <paramref name="gcDeviation"/> (GC-content outlier detection). Pooled GC = Σ(G+C)/Σ(A+C+G+T) over all
    /// sequences (the ASSEMBLY-STATS-001 / QUAST rule); windows are those of
    /// <see cref="CalculateLocalQuality"/> (step = windowSize/2), a truncated trailing window is tested only if
    /// the previous window does not already cover it, and a window without any A/C/G/T/U base has no GC and is
    /// not tested. No pooled GC (no valid base anywhere) → no GC test.</item>
    /// <item><b>Low complexity</b> — the perfect intervals of the symmetric DUST algorithm (Morgulis et al.
    /// 2006; canonical <see cref="SequenceComplexity.FindLowComplexityIntervals(string,int,double,int,DustEngine)"/>,
    /// window 64, score &gt; <paramref name="dustThreshold"/>, default 2.0 = dustmasker level 20), exact intervals.</item>
    /// </list>
    /// Regions are inclusive 0-based [Start, End] = union of the flagged windows/intervals.
    /// </remarks>
    /// <param name="sequences">(id, sequence) pairs.</param>
    /// <param name="gcDeviation">Maximum allowed |window GC − pooled GC| as a fraction (default 0.15; ≥ 0).</param>
    /// <param name="dustThreshold">DUST score threshold (default 2.0; ≥ 0).</param>
    /// <param name="maxNFraction">Maximum allowed N fraction per window (default 0.1; in [0,1]).</param>
    /// <param name="windowSize">GC/N window length (≥ 1; default 500).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A parameter is out of range (also a NaN).</exception>
    public static IEnumerable<SuspiciousRegion> FindSuspiciousRegions(
        IEnumerable<(string Id, string Sequence)> sequences,
        double gcDeviation = 0.15,
        double dustThreshold = 2.0,
        double maxNFraction = 0.1,
        int windowSize = 500)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (!(gcDeviation >= 0) || double.IsInfinity(gcDeviation))
            throw new ArgumentOutOfRangeException(nameof(gcDeviation), gcDeviation, "Must be a finite value ≥ 0.");
        if (!(maxNFraction >= 0 && maxNFraction <= 1))
            throw new ArgumentOutOfRangeException(nameof(maxNFraction), maxNFraction, "Must be in [0, 1].");
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        if (!(dustThreshold >= 0) || double.IsInfinity(dustThreshold))
            throw new ArgumentOutOfRangeException(nameof(dustThreshold), dustThreshold, "Must be a finite value ≥ 0.");
        return FindSuspiciousRegionsIterator(
            sequences.Select(e => (e.Id, Sequence: RequireSequence(e, nameof(sequences)))).ToList(), gcDeviation, dustThreshold, maxNFraction, windowSize);
    }

    private const int NRich = 1, GcOutlier = 2, LowComplexity = 4;

    private static IEnumerable<SuspiciousRegion> FindSuspiciousRegionsIterator(
        List<(string Id, string Sequence)> seqList, double gcDeviation, double dustThreshold,
        double maxNFraction, int windowSize)
    {
        // Pooled GC over the whole assembly: Σ(G+C) / Σ(A+C+G+T) (QUAST rule, ASSEMBLY-STATS-001).
        long totalGc = 0, totalValid = 0;
        foreach (var (_, sequence) in seqList)
        {
            var (gc, valid) = sequence.AsSpan().CountGcAndValidNucleotides();
            totalGc += gc;
            totalValid += valid;
        }
        double? globalGc = totalValid > 0 ? totalGc / (double)totalValid : null;
        int step = Math.Max(1, windowSize / 2);

        foreach (var (id, sequence) in seqList)
        {
            // (start, end inclusive, criterion, measured value)
            var hits = new List<(int Start, int End, int Flag, double Value)>();

            int previousStart = -1;
            foreach (var w in EnumerateCompositionWindows(sequence, windowSize, step))
            {
                bool covered = w.Length < windowSize && previousStart >= 0 && previousStart + windowSize >= sequence.Length;
                previousStart = w.Start;
                if (covered)
                    continue;

                int end = w.Start + w.Length - 1;
                double nFraction = w.NCount / (double)w.Length;
                if (nFraction > maxNFraction)
                    hits.Add((w.Start, end, NRich, nFraction));

                if (globalGc is { } g && w.ValidBases > 0 && Math.Abs(w.Gc - g) > gcDeviation)
                    hits.Add((w.Start, end, GcOutlier, Math.Abs(w.Gc - g)));
            }

            foreach (var (start, endExclusive) in SequenceComplexity.FindLowComplexityIntervals(
                sequence, threshold: dustThreshold))
                hits.Add((start, endExclusive - 1, LowComplexity,
                    SequenceComplexity.CalculateDustScore(sequence.Substring(start, endExclusive - start))));

            hits.Sort((x, y) => x.Start != y.Start ? x.Start.CompareTo(y.Start) : x.End.CompareTo(y.End));

            int i = 0;
            while (i < hits.Count)
            {
                int rStart = hits[i].Start, rEnd = hits[i].End, flags = 0;
                double maxGc = 0, maxN = 0, maxDust = 0;
                for (; i < hits.Count && hits[i].Start <= rEnd + 1; i++)
                {
                    var h = hits[i];
                    rEnd = Math.Max(rEnd, h.End);
                    flags |= h.Flag;
                    if (h.Flag == NRich) maxN = Math.Max(maxN, h.Value);
                    else if (h.Flag == GcOutlier) maxGc = Math.Max(maxGc, h.Value);
                    else maxDust = Math.Max(maxDust, h.Value);
                }

                var reasons = new List<string>();
                if ((flags & NRich) != 0) reasons.Add("High N content");
                if ((flags & GcOutlier) != 0) reasons.Add("GC deviation");
                if ((flags & LowComplexity) != 0) reasons.Add("Low complexity");
                yield return new SuspiciousRegion(id, rStart, rEnd, string.Join("; ", reasons), maxGc, maxN, maxDust);
            }
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Filters sequences by length.
    /// </summary>
    /// <remarks>
    /// Keeps, in input order, the sequences with <c>minLength ≤ length ≤ maxLength</c> (both bounds
    /// inclusive, as QUAST <c>--min-contig</c> <c>len(seq) &gt;= min_contig</c> and
    /// <c>seqkit seq -m/-M</c>); length counts every character, N included. A minimum above the
    /// maximum yields nothing.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<(string Id, string Sequence)> FilterByLength(
        IEnumerable<(string Id, string Sequence)> sequences,
        int minLength = 0,
        int maxLength = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return sequences.Where(s =>
        {
            int length = RequireSequence(s, nameof(sequences)).Length;
            return length >= minLength && length <= maxLength;
        });
    }

    /// <summary>
    /// Sorts sequences by length (descending by default).
    /// </summary>
    /// <remarks>
    /// Stable: sequences of equal length keep their input order (LINQ <c>OrderBy</c>/<c>OrderByDescending</c>).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<(string Id, string Sequence)> SortByLength(
        IEnumerable<(string Id, string Sequence)> sequences,
        bool descending = true)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var checkedSequences = sequences.Select(s => (s.Id, Sequence: RequireSequence(s, nameof(sequences))));
        return descending ?
            checkedSequences.OrderByDescending(s => s.Sequence.Length) :
            checkedSequences.OrderBy(s => s.Sequence.Length);
    }

    /// <summary>
    /// Calculates the sequence length histogram: bin "&lt;b" counts lengths in [previous bound, b)
    /// and "&gt;=last" counts lengths at or above the last bound.
    /// </summary>
    /// <param name="lengths">Sequence lengths.</param>
    /// <param name="bins">Strictly increasing upper bounds (empty → 100, 500, 1 000, …, 1 000 000).</param>
    /// <exception cref="ArgumentNullException"><paramref name="lengths"/> or <paramref name="bins"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="bins"/> is not strictly increasing
    /// (an unsorted bound list would put a length in the first larger bound, not its own bin).</exception>
    public static IReadOnlyDictionary<string, int> CalculateLengthDistribution(
        IEnumerable<int> lengths,
        params int[] bins)
    {
        ArgumentNullException.ThrowIfNull(lengths);
        ArgumentNullException.ThrowIfNull(bins);
        if (bins.Length == 0)
            bins = new[] { 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };
        for (int i = 1; i < bins.Length; i++)
        {
            if (bins[i] <= bins[i - 1])
                throw new ArgumentException("Bin upper bounds must be strictly increasing.", nameof(bins));
        }

        var distribution = bins.ToDictionary(b => $"<{b}", _ => 0);
        distribution[$">={bins[^1]}"] = 0;

        foreach (int length in lengths)
        {
            bool found = false;
            foreach (int bin in bins)
            {
                if (length < bin)
                {
                    distribution[$"<{bin}"]++;
                    found = true;
                    break;
                }
            }
            if (!found)
                distribution[$">={bins[^1]}"]++;
        }

        return distribution;
    }

    #endregion
}
