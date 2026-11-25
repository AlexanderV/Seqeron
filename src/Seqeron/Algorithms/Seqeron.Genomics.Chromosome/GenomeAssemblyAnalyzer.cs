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
    /// BUSCO-like completeness result.
    /// </summary>
    public readonly record struct CompletenessResult(
        int TotalGenes,
        int Complete,
        int CompleteSingleCopy,
        int CompleteDuplicated,
        int Fragmented,
        int Missing,
        double CompletenessPercent,
        double DuplicationPercent);

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
    /// Assesses assembly completeness using marker genes (BUSCO-like).
    /// </summary>
    public static CompletenessResult AssessCompleteness(
        IEnumerable<(string Id, string Sequence)> assembly,
        IEnumerable<(string GeneId, string Sequence)> markerGenes,
        double identityThreshold = 0.9,
        double coverageThreshold = 0.9)
    {
        var assemblySeqs = assembly.ToList();
        var markers = markerGenes.ToList();

        if (markers.Count == 0)
            return new CompletenessResult(0, 0, 0, 0, 0, 0, 0, 0);

        int complete = 0;
        int singleCopy = 0;
        int duplicated = 0;
        int fragmented = 0;
        int missing = 0;

        foreach (var (_, geneSeq) in markers)
        {
            var hits = FindMarkerHits(geneSeq, assemblySeqs, identityThreshold, coverageThreshold);

            if (hits.Count == 0)
            {
                missing++;
            }
            else if (hits.Any(h => h.Coverage >= coverageThreshold))
            {
                complete++;
                var completeHits = hits.Count(h => h.Coverage >= coverageThreshold);
                if (completeHits == 1)
                    singleCopy++;
                else
                    duplicated++;
            }
            else
            {
                fragmented++;
            }
        }

        double completeness = markers.Count > 0 ? complete * 100.0 / markers.Count : 0;
        double duplication = complete > 0 ? duplicated * 100.0 / complete : 0;

        return new CompletenessResult(
            markers.Count,
            complete,
            singleCopy,
            duplicated,
            fragmented,
            missing,
            completeness,
            duplication);
    }

    /// <summary>
    /// Finds marker gene hits in assembly.
    /// </summary>
    private static List<(string SeqId, double Identity, double Coverage)> FindMarkerHits(
        string markerSeq,
        List<(string Id, string Sequence)> assembly,
        double identityThreshold,
        double coverageThreshold)
    {
        var hits = new List<(string, double, double)>();

        // Simplified k-mer based search
        int kmerSize = Math.Min(31, markerSeq.Length / 3);
        if (kmerSize < 11) kmerSize = 11;

        var markerKmers = new HashSet<string>();
        for (int i = 0; i <= markerSeq.Length - kmerSize; i++)
        {
            markerKmers.Add(markerSeq.Substring(i, kmerSize).ToUpperInvariant());
        }

        foreach (var (seqId, sequence) in assembly)
        {
            int matchingKmers = 0;
            var seqUpper = sequence.ToUpperInvariant();

            for (int i = 0; i <= seqUpper.Length - kmerSize; i++)
            {
                if (markerKmers.Contains(seqUpper.Substring(i, kmerSize)))
                    matchingKmers++;
            }

            double coverage = markerKmers.Count > 0 ?
                matchingKmers / (double)markerKmers.Count : 0;

            if (coverage >= coverageThreshold * 0.5) // Relaxed for fragmented
            {
                // Estimate identity from k-mer matches
                double identity = Math.Min(1.0, coverage * 1.1);
                if (identity >= identityThreshold * 0.8)
                {
                    hits.Add((seqId, identity, coverage));
                }
            }
        }

        return hits;
    }

    /// <summary>
    /// Estimates genome completeness by k-mer spectrum analysis.
    /// </summary>
    public static (double Completeness, double ErrorRate, long EstimatedGenomeSize)
        EstimateCompletenessFromKmers(
            IEnumerable<(string Kmer, int Count)> kmerSpectrum,
            int expectedCoverage = 0)
    {
        var spectrum = kmerSpectrum.ToList();

        if (spectrum.Count == 0)
            return (0, 0, 0);

        // Find coverage peak (mode)
        var countDistribution = spectrum
            .Where(k => k.Count > 1) // Exclude singletons (likely errors)
            .GroupBy(k => k.Count)
            .Select(g => (Count: g.Key, Frequency: g.Count()))
            .OrderByDescending(g => g.Frequency)
            .ToList();

        if (countDistribution.Count == 0)
            return (0, 1.0, 0);

        int peakCoverage = expectedCoverage > 0 ?
            expectedCoverage :
            countDistribution[0].Count;

        // Calculate error rate from singleton ratio
        int singletons = spectrum.Count(k => k.Count == 1);
        double errorRate = singletons / (double)spectrum.Count;

        // Estimate genome size
        long totalKmers = spectrum.Sum(k => (long)k.Count);
        long estimatedSize = peakCoverage > 0 ? totalKmers / peakCoverage : 0;

        // Estimate completeness
        int solidKmers = spectrum.Count(k => k.Count >= peakCoverage / 2);
        double completeness = solidKmers / (double)Math.Max(1, spectrum.Count - singletons);

        return (completeness, errorRate, estimatedSize);
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
