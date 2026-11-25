namespace Seqeron.Genomics.Chromosome;

/// <summary>
/// Provides genome assembly quality assessment algorithms.
/// Includes N50/L50, contiguity metrics, completeness analysis, and repeat detection.
/// </summary>
public static class GenomeAssemblyAnalyzer
{
    #region Records

    /// <summary>
    /// Represents an assembled sequence (contig or scaffold).
    /// </summary>
    public readonly record struct AssembledSequence(
        string Id,
        string Sequence,
        int Length,
        double GcContent,
        int GapCount,
        int TotalGapLength);

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
    /// Assembly comparison result.
    /// </summary>
    public readonly record struct AssemblyComparison(
        string Assembly1Name,
        string Assembly2Name,
        double AlignedFraction1,
        double AlignedFraction2,
        int Breakpoints,
        int Inversions,
        int Translocations,
        double SequenceIdentity);

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
    /// <c>len − count('N')</c> after <c>correct_seq</c> upper-cases and turns IUPAC codes into N);
    /// counted by the canonical <see cref="SequenceExtensions.CountGcAndValidNucleotides"/>, so
    /// N, IUPAC ambiguity codes and gap symbols are excluded from both counts (U, which QUAST rejects,
    /// counts as a valid base). 0 when there is no A/C/G/T.</item>
    /// <item>Gaps = maximal runs of N/n (QUAST <c>split_by_ns</c> with the threshold 1); TotalGapLength
    /// = number of N/n; TotalLengthNoGaps = TotalLength − TotalGapLength; GapPercentage = 100 ·
    /// TotalGapLength / TotalLength. IUPAC ambiguity codes are bases, not gaps.</item>
    /// <item>MeanLength = TotalLength / count; MedianLength = sample median of the lengths (mean of the
    /// two central values for an even count, <see cref="StatisticsHelper.Median"/>).</item>
    /// </list>
    /// An empty input returns all zeros (QUAST reports these fields as undefined).
    /// </remarks>
    public static AssemblyStatistics CalculateStatistics(
        IEnumerable<(string Id, string Sequence)> sequences)
    {
        var seqList = sequences.Select(s => ParseSequence(s.Id, s.Sequence)).ToList();

        if (seqList.Count == 0)
        {
            return new AssemblyStatistics(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        var lengths = seqList.Select(s => s.Length).OrderByDescending(l => l).ToList();
        long totalLength = lengths.Sum(l => (long)l);
        long totalNoGaps = seqList.Sum(s => (long)(s.Length - s.TotalGapLength));

        // Calculate N50/L50
        var n50Stats = CalculateNx(lengths, totalLength, 50);
        var n90Stats = CalculateNx(lengths, totalLength, 90);

        // GC content: pooled G+C over A+C+G+T (QUAST GC_content) via the canonical Core counter.
        long gcCount = 0;
        long acgtCount = 0;
        foreach (var s in seqList)
        {
            var (gc, valid) = s.Sequence.AsSpan().CountGcAndValidNucleotides();
            gcCount += gc;
            acgtCount += valid;
        }
        double gcContent = acgtCount > 0 ? gcCount / (double)acgtCount : 0;

        // Gap statistics
        int totalGaps = seqList.Sum(s => s.GapCount);
        long totalGapLength = seqList.Sum(s => (long)s.TotalGapLength);
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

    /// <summary>
    /// Parses sequence and extracts metrics.
    /// </summary>
    private static AssembledSequence ParseSequence(string id, string sequence)
    {
        int gapCount = 0;
        int totalGapLength = 0;
        foreach (var (_, length) in EnumerateNRuns(sequence))
        {
            gapCount++;
            totalGapLength += length;
        }

        double gc = sequence.AsSpan().CalculateGcFraction();

        return new AssembledSequence(id, sequence, sequence.Length, gc, gapCount, totalGapLength);
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
    public static NxStatistics CalculateNx(IEnumerable<int> lengths, int threshold)
    {
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
    public static int CalculateN50(IEnumerable<int> lengths)
    {
        return CalculateNx(lengths, N50ThresholdPercent).Nx;
    }

    /// <summary>
    /// Calculates Nx statistics for multiple thresholds.
    /// </summary>
    public static IEnumerable<NxStatistics> CalculateNxCurve(
        IEnumerable<int> lengths,
        params int[] thresholds)
    {
        var sorted = lengths.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);

        if (thresholds.Length == 0)
            thresholds = new[] { 10, 20, 30, 40, 50, 60, 70, 80, 90 };

        foreach (int t in thresholds.OrderBy(t => t))
        {
            yield return CalculateNx(sorted, total, t);
        }
    }

    /// <summary>
    /// Calculates auN (area under the Nx curve): auN = Σ Lᵢ² / Σ Lⱼ.
    /// </summary>
    /// <remarks>
    /// Li H (2020), "auN: a new metric to measure assembly contiguity"
    /// (lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity); QUAST <c>au_metric</c>
    /// (<c>quast_libs/N50.py</c>). Returns 0 for an empty or all-zero list (QUAST <c>None</c>).
    /// </remarks>
    public static double CalculateAuN(IEnumerable<int> lengths)
    {
        var sorted = lengths.OrderByDescending(l => l).ToList();
        if (sorted.Count == 0) return 0;

        long total = sorted.Sum(l => (long)l);
        if (total == 0) return 0;

        // auN is the sum of (length * length) / total
        double auN = sorted.Sum(l => (double)l * l) / total;
        return auN;
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
    /// Enumerates the maximal runs of <c>N</c>/<c>n</c> in <paramref name="sequence"/> as
    /// (0-based start, length). The single N-run scanner behind <see cref="FindGaps"/>,
    /// <see cref="AnalyzeScaffolds"/>, <see cref="ExtractContigs"/> and the gap fields of
    /// <see cref="CalculateStatistics"/>. Only N/n is a gap symbol; IUPAC ambiguity codes are bases.
    /// </summary>
    private static IEnumerable<(int Start, int Length)> EnumerateNRuns(string sequence)
    {
        int i = 0;
        while (i < sequence.Length)
        {
            if (sequence[i] is 'N' or 'n')
            {
                int start = i;
                while (i < sequence.Length && sequence[i] is 'N' or 'n')
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
        SplitAtGaps(string sequence, int minGapLength)
    {
        var contigs = new List<(int Start, int End)>();
        var gaps = new List<(int Start, int Length)>();
        int contigStart = 0;
        foreach (var run in EnumerateNRuns(sequence))
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
    /// Identifies all gaps (maximal runs of <c>N</c>/<c>n</c>) in assembly sequences.
    /// </summary>
    /// <remarks>
    /// One <see cref="GapInfo"/> per maximal N-run of length ≥ <paramref name="minGapLength"/>, in
    /// sequence order. Coordinates are 0-based with an <b>inclusive</b> <c>End</c>
    /// (<c>Length = End − Start + 1</c>); the AGP v2.1 equivalent is <c>[Start + 1, End + 1]</c>
    /// (1-based inclusive) and the BED / <c>seqtk gap</c> equivalent is <c>[Start, End + 1)</c>.
    /// Only <c>N</c>/<c>n</c> marks a gap: IUPAC ambiguity codes are bases (AGP v2.1, NCBI gapped
    /// format; QUAST <c>correct_seq</c> and <c>seqtk gap</c> instead treat them as N — documented
    /// deviation, consistent with <see cref="CalculateStatistics"/>). With the default 1 the gaps are
    /// exactly those counted by <see cref="CalculateStatistics"/>; pass
    /// <see cref="DefaultScaffoldGapMinLength"/> (10) for the NCBI assembly-gap convention.
    /// <see cref="GapInfo.GapType"/> is a Seqeron length class (&lt; 10 "Short", &lt; 100 "Medium",
    /// &lt; 1000 "Long", else "Scaffold"); it is <b>not</b> the AGP column 7b gap type, which records
    /// linkage evidence that a sequence does not carry.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<GapInfo> FindGaps(
        IEnumerable<(string Id, string Sequence)> sequences,
        int minGapLength = 1)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ValidateMinGapLength(minGapLength);
        return FindGapsIterator(sequences, minGapLength);
    }

    private static IEnumerable<GapInfo> FindGapsIterator(
        IEnumerable<(string Id, string Sequence)> sequences, int minGapLength)
    {
        foreach (var entry in sequences)
        {
            string sequence = RequireSequence(entry, nameof(sequences));
            foreach (var (start, length) in EnumerateNRuns(sequence))
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
    /// coordinates; a leading or trailing gap opens no empty contig.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="scaffolds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<ScaffoldStructure> AnalyzeScaffolds(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minGapLength = DefaultScaffoldGapMinLength)
    {
        ArgumentNullException.ThrowIfNull(scaffolds);
        ValidateMinGapLength(minGapLength);
        return AnalyzeScaffoldsIterator(scaffolds, minGapLength);
    }

    private static IEnumerable<ScaffoldStructure> AnalyzeScaffoldsIterator(
        IEnumerable<(string Id, string Sequence)> scaffolds, int minGapLength)
    {
        foreach (var entry in scaffolds)
        {
            string id = entry.Id;
            string sequence = RequireSequence(entry, nameof(scaffolds));
            var (pieces, runs) = SplitAtGaps(sequence, minGapLength);

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
    /// are returned in scaffold order, named <c>{id}_contig{n}</c> with n counting the returned contigs
    /// (QUAST also counts the dropped ones). The contigs are exactly those of
    /// <see cref="AnalyzeScaffolds"/> with the same <paramref name="minGapLength"/>; case is preserved.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="scaffolds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minContigLength"/> &lt; 0 or
    /// <paramref name="minGapLength"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">A sequence is null (raised during enumeration).</exception>
    public static IEnumerable<(string Id, string Sequence)> ExtractContigs(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minContigLength = 200,
        int minGapLength = DefaultScaffoldGapMinLength)
    {
        ArgumentNullException.ThrowIfNull(scaffolds);
        if (minContigLength < 0)
            throw new ArgumentOutOfRangeException(nameof(minContigLength), minContigLength,
                "The minimum contig length cannot be negative.");
        ValidateMinGapLength(minGapLength);
        return ExtractContigsIterator(scaffolds, minContigLength, minGapLength);
    }

    private static IEnumerable<(string Id, string Sequence)> ExtractContigsIterator(
        IEnumerable<(string Id, string Sequence)> scaffolds, int minContigLength, int minGapLength)
    {
        foreach (var entry in scaffolds)
        {
            string sequence = RequireSequence(entry, nameof(scaffolds));
            int contigNum = 1;
            foreach (var (start, end) in SplitAtGaps(sequence, minGapLength).Contigs)
            {
                if (end - start >= minContigLength)
                    yield return ($"{entry.Id}_contig{contigNum++}", sequence[start..end]);
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
        => assembly.Select(s => (s.Id, RequireSequence(s, paramName).ToUpperInvariant())).ToList();

    private static string RequireMarker(string id, string? sequence, string paramName)
    {
        if (sequence is null)
            throw new ArgumentException($"Marker '{id}' sequence is null.", paramName);
        if (sequence.Length == 0)
            throw new ArgumentException($"Marker '{id}' sequence is empty.", paramName);
        return sequence.ToUpperInvariant();
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
            if (!spectrum.TryAdd(kmer.ToUpperInvariant(), count))
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
    /// Identifies repetitive sequences using k-mer frequency.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, int Copies)> FindRepetitiveRegions(
        IEnumerable<(string Id, string Sequence)> sequences,
        int kmerSize = 15,
        int minCopies = 3,
        int windowSize = 100)
    {
        // Build k-mer frequency table
        var kmerCounts = new Dictionary<string, int>();

        foreach (var (_, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer) + 1;
                }
            }
        }

        // Find regions with high k-mer frequency
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            int regionStart = -1;
            int maxCopies = 0;

            for (int i = 0; i <= upper.Length - kmerSize; i += windowSize / 10)
            {
                // Check window
                int windowEnd = Math.Min(i + windowSize, upper.Length);
                int highCopyKmers = 0;
                int localMaxCopies = 0;

                for (int j = i; j <= windowEnd - kmerSize; j++)
                {
                    string kmer = upper.Substring(j, kmerSize);
                    if (!kmer.Contains('N') && kmerCounts.TryGetValue(kmer, out int count) && count >= minCopies)
                    {
                        highCopyKmers++;
                        localMaxCopies = Math.Max(localMaxCopies, count);
                    }
                }

                bool isRepetitive = highCopyKmers > windowSize / (kmerSize * 2);

                if (isRepetitive && regionStart < 0)
                {
                    regionStart = i;
                    maxCopies = localMaxCopies;
                }
                else if (!isRepetitive && regionStart >= 0)
                {
                    yield return (id, regionStart, i - 1, maxCopies);
                    regionStart = -1;
                    maxCopies = 0;
                }
                else if (isRepetitive)
                {
                    maxCopies = Math.Max(maxCopies, localMaxCopies);
                }
            }

            if (regionStart >= 0)
            {
                yield return (id, regionStart, upper.Length - 1, maxCopies);
            }
        }
    }

    /// <summary>
    /// Identifies tandem repeats.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, string Unit, int Copies, double Purity)>
        FindTandemRepeats(
            IEnumerable<(string Id, string Sequence)> sequences,
            int minUnitLength = 2,
            int maxUnitLength = 50,
            int minCopies = 3)
    {
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();

            // `i` normally steps by 1 but jumps past a called repeat (below); a while loop makes that
            // in-body skip explicit where a for loop would flag S127.
            int i = 0;
            while (i < upper.Length - minUnitLength * minCopies)
            {
                for (int unitLen = minUnitLength; unitLen <= maxUnitLength && i + unitLen * minCopies <= upper.Length; unitLen++)
                {
                    string unit = upper.Substring(i, unitLen);
                    if (unit.Contains('N')) continue;

                    // Count copies
                    int copies = 1;
                    int matches = 0;
                    int pos = i + unitLen;

                    while (pos + unitLen <= upper.Length)
                    {
                        int unitMatches = 0;
                        for (int j = 0; j < unitLen; j++)
                        {
                            if (upper[pos + j] == unit[j])
                                unitMatches++;
                        }

                        if (unitMatches >= unitLen * 0.8) // 80% match
                        {
                            copies++;
                            matches += unitMatches;
                            pos += unitLen;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (copies >= minCopies)
                    {
                        int totalLength = copies * unitLen;
                        double purity = (unitLen + matches) / (double)(totalLength);

                        yield return (id, i, i + totalLength - 1, unit, copies, purity);

                        // Skip past this repeat
                        i += totalLength - 1;
                        break;
                    }
                }

                i++;
            }
        }
    }

    /// <summary>
    /// Calculates repeat content statistics.
    /// </summary>
    public static (long TotalRepeatLength, double RepeatPercentage, IReadOnlyDictionary<string, long> RepeatClassLengths)
        CalculateRepeatContent(
            IEnumerable<RepeatAnnotation> repeats,
            long genomeLength)
    {
        var repeatList = repeats.ToList();

        if (repeatList.Count == 0 || genomeLength == 0)
            return (0, 0, new Dictionary<string, long>());

        long totalLength = repeatList.Sum(r => (long)(r.End - r.Start + 1));
        double percentage = totalLength * 100.0 / genomeLength;

        var classLengths = repeatList
            .GroupBy(r => r.RepeatClass)
            .ToDictionary(g => g.Key, g => g.Sum(r => (long)(r.End - r.Start + 1)));

        return (totalLength, percentage, classLengths);
    }

    #endregion

    #region Assembly Comparison

    /// <summary>
    /// Compares two assemblies using k-mer content.
    /// </summary>
    public static AssemblyComparison CompareAssemblies(
        IEnumerable<(string Id, string Sequence)> assembly1,
        IEnumerable<(string Id, string Sequence)> assembly2,
        string name1 = "Assembly1",
        string name2 = "Assembly2",
        int kmerSize = 21)
    {
        var kmers1 = ExtractKmers(assembly1, kmerSize);
        var kmers2 = ExtractKmers(assembly2, kmerSize);

        int shared = kmers1.Keys.Count(k => kmers2.ContainsKey(k));

        double aligned1 = kmers1.Count > 0 ? shared / (double)kmers1.Count : 0;
        double aligned2 = kmers2.Count > 0 ? shared / (double)kmers2.Count : 0;

        // Estimate identity from shared k-mers
        double identity = (aligned1 + aligned2) / 2;

        return new AssemblyComparison(
            name1, name2,
            aligned1, aligned2,
            0, 0, 0, identity);
    }

    /// <summary>
    /// Extracts k-mers from assembly.
    /// </summary>
    private static Dictionary<string, int> ExtractKmers(
        IEnumerable<(string Id, string Sequence)> sequences,
        int kmerSize)
    {
        var kmers = new Dictionary<string, int>();

        foreach (var (_, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    kmers[kmer] = kmers.GetValueOrDefault(kmer) + 1;
                }
            }
        }

        return kmers;
    }

    /// <summary>
    /// Finds syntenic blocks between assemblies.
    /// </summary>
    public static IEnumerable<(string Seq1, int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        FindSyntenicBlocks(
            IEnumerable<(string Id, string Sequence)> assembly1,
            IEnumerable<(string Id, string Sequence)> assembly2,
            int minBlockSize = 1000,
            int kmerSize = 21)
    {
        var asm1 = assembly1.ToList();
        var asm2 = assembly2.ToList();

        // Build k-mer index for assembly2
        var kmerIndex = new Dictionary<string, List<(string SeqId, int Position)>>();

        foreach (var (seqId, sequence) in asm2)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    if (!kmerIndex.ContainsKey(kmer))
                        kmerIndex[kmer] = new List<(string, int)>();
                    kmerIndex[kmer].Add((seqId, i));
                }
            }
        }

        // Find matching regions
        foreach (var (seq1Id, seq1) in asm1)
        {
            var upper1 = seq1.ToUpperInvariant();
            var anchors = new List<(int Pos1, string Seq2, int Pos2)>();

            for (int i = 0; i <= upper1.Length - kmerSize; i += kmerSize)
            {
                string kmer = upper1.Substring(i, kmerSize);
                if (!kmer.Contains('N') && kmerIndex.TryGetValue(kmer, out var hits))
                {
                    foreach (var (seq2Id, pos2) in hits)
                    {
                        anchors.Add((i, seq2Id, pos2));
                    }
                }
            }

            // Cluster anchors into blocks
            var blocks = ClusterAnchors(anchors, minBlockSize);

            foreach (var block in blocks)
            {
                yield return (seq1Id, block.Start1, block.End1,
                             block.Seq2, block.Start2, block.End2, block.IsInverted);
            }
        }
    }

    /// <summary>
    /// Clusters k-mer anchors into syntenic blocks.
    /// </summary>
    private static List<(int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        ClusterAnchors(
            List<(int Pos1, string Seq2, int Pos2)> anchors,
            int minBlockSize)
    {
        var result = new List<(int, int, string, int, int, bool)>();

        if (anchors.Count == 0)
            return result;

        // Group by target sequence
        var bySeq2 = anchors.GroupBy(a => a.Seq2);

        foreach (var group in bySeq2)
        {
            var sorted = group.OrderBy(a => a.Pos1).ToList();

            int blockStart1 = sorted[0].Pos1;
            int blockStart2 = sorted[0].Pos2;
            int lastPos1 = sorted[0].Pos1;
            int lastPos2 = sorted[0].Pos2;
            bool? isInverted = null;

            for (int i = 1; i < sorted.Count; i++)
            {
                int gap1 = sorted[i].Pos1 - lastPos1;
                int gap2 = sorted[i].Pos2 - lastPos2;

                bool currentInverted = gap2 < 0;

                if (isInverted == null)
                    isInverted = currentInverted;

                // Check if continuous
                if (gap1 > minBlockSize || Math.Abs(gap2) > minBlockSize || currentInverted != isInverted)
                {
                    // End current block
                    if (lastPos1 - blockStart1 >= minBlockSize)
                    {
                        result.Add((blockStart1, lastPos1, group.Key,
                                   Math.Min(blockStart2, lastPos2),
                                   Math.Max(blockStart2, lastPos2),
                                   isInverted ?? false));
                    }

                    blockStart1 = sorted[i].Pos1;
                    blockStart2 = sorted[i].Pos2;
                    isInverted = null;
                }

                lastPos1 = sorted[i].Pos1;
                lastPos2 = sorted[i].Pos2;
            }

            // Handle last block
            if (lastPos1 - blockStart1 >= minBlockSize)
            {
                result.Add((blockStart1, lastPos1, group.Key,
                           Math.Min(blockStart2, lastPos2),
                           Math.Max(blockStart2, lastPos2),
                           isInverted ?? false));
            }
        }

        return result;
    }

    #endregion

    #region Quality Assessment

    /// <summary>
    /// Calculates per-base quality metrics.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Position, int WindowSize, double GcContent, int NCount, double Complexity)>
        CalculateLocalQuality(
            IEnumerable<(string Id, string Sequence)> sequences,
            int windowSize = 1000)
    {
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();

            for (int i = 0; i < upper.Length; i += windowSize / 2)
            {
                int end = Math.Min(i + windowSize, upper.Length);
                string window = upper[i..end];

                int nCount = window.Count(c => c == 'N');
                double gcContent = window.CalculateGcFractionFast();
                double complexity = CalculateLinguisticComplexity(window);

                yield return (id, i, end - i, gcContent, nCount, complexity);
            }
        }
    }

    /// <summary>
    /// Calculates linguistic complexity (ratio of distinct k-mers to possible k-mers).
    /// </summary>
    private static double CalculateLinguisticComplexity(string sequence, int kmerSize = 4)
    {
        if (sequence.Length < kmerSize)
            return 0;

        var kmers = new HashSet<string>();
        int validKmers = 0;

        for (int i = 0; i <= sequence.Length - kmerSize; i++)
        {
            string kmer = sequence.Substring(i, kmerSize);
            if (!kmer.Contains('N'))
            {
                kmers.Add(kmer);
                validKmers++;
            }
        }

        if (validKmers == 0)
            return 0;

        // Maximum possible k-mers for DNA
        int maxKmers = (int)Math.Pow(4, kmerSize);
        int possibleKmers = Math.Min(validKmers, maxKmers);

        return kmers.Count / (double)possibleKmers;
    }

    /// <summary>
    /// Identifies potentially misassembled regions.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, string Reason, double Score)>
        FindSuspiciousRegions(
            IEnumerable<(string Id, string Sequence)> sequences,
            double gcDeviation = 0.15,
            double minComplexity = 0.3)
    {
        // First pass: calculate global GC
        var seqList = sequences.ToList();
        long totalGc = 0;
        long totalBases = 0;

        foreach (var (_, sequence) in seqList)
        {
            totalGc += sequence.Count(c => c == 'G' || c == 'g' || c == 'C' || c == 'c');
            totalBases += sequence.Count(c => c != 'N' && c != 'n');
        }

        double globalGc = totalBases > 0 ? totalGc / (double)totalBases : 0.5;

        // Second pass: find anomalies
        foreach (var (id, sequence) in seqList)
        {
            var quality = CalculateLocalQuality(
                new[] { (id, sequence) }, windowSize: 500).ToList();

            int regionStart = -1;
            string? reason = null;
            double maxScore = 0;

            foreach (var (_, pos, size, gc, nCount, complexity) in quality)
            {
                bool suspicious = false;
                double score = 0;
                string? currentReason = null;

                // High N content
                if (nCount > size * 0.1)
                {
                    suspicious = true;
                    score = Math.Max(score, nCount / (double)size);
                    currentReason = "High N content";
                }

                // GC deviation
                if (Math.Abs(gc - globalGc) > gcDeviation)
                {
                    suspicious = true;
                    score = Math.Max(score, Math.Abs(gc - globalGc) / gcDeviation);
                    currentReason = currentReason == null ? "GC deviation" : currentReason + "; GC deviation";
                }

                // Low complexity
                if (complexity < minComplexity && nCount < size * 0.5)
                {
                    suspicious = true;
                    score = Math.Max(score, 1 - complexity / minComplexity);
                    currentReason = currentReason == null ? "Low complexity" : currentReason + "; Low complexity";
                }

                if (suspicious && regionStart < 0)
                {
                    regionStart = pos;
                    reason = currentReason;
                    maxScore = score;
                }
                else if (!suspicious && regionStart >= 0)
                {
                    yield return (id, regionStart, pos - 1, reason!, maxScore);
                    regionStart = -1;
                    reason = null;
                    maxScore = 0;
                }
                else if (suspicious)
                {
                    maxScore = Math.Max(maxScore, score);
                    if (currentReason != null && !reason!.Contains(currentReason))
                        reason += "; " + currentReason;
                }
            }

            if (regionStart >= 0)
            {
                yield return (id, regionStart, sequence.Length - 1, reason!, maxScore);
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
    /// <exception cref="ArgumentException"><paramref name="bins"/> is not strictly increasing
    /// (an unsorted bound list would put a length in the first larger bound, not its own bin).</exception>
    public static IReadOnlyDictionary<string, int> CalculateLengthDistribution(
        IEnumerable<int> lengths,
        params int[] bins)
    {
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
