namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region HRD score (ONCO-HRD-001)

    /// <summary>
    /// Myriad myChoice CDx / Telli et al. (2016) genomic-instability cutoff: a tumour is HRD-high when its
    /// combined HRD score is at or above this value. Source: Telli ML et al. (2016), Clin Cancer Res
    /// 22(15):3764–3773 — "HR deficiency, defined as HRD score ≥42 or BRCA1/2 mutation". Boundary inclusive.
    /// </summary>
    public const int HrdHighScoreThreshold = 42;

    /// <summary>HRD (homologous recombination deficiency) classification of a combined genomic-scar score.</summary>
    public enum HrdStatus
    {
        /// <summary>HRD score below the <see cref="HrdHighScoreThreshold"/> cutoff (HR-proficient signal).</summary>
        HrdNegative,

        /// <summary>HRD score at or above the <see cref="HrdHighScoreThreshold"/> cutoff (HR-deficient).</summary>
        HrdHigh
    }

    /// <summary>
    /// The three genomic-scar component counts that sum to the combined HRD score.
    /// </summary>
    /// <param name="Loh">
    /// HRD-LOH score: number of LOH regions longer than 15 Mb but shorter than a whole chromosome
    /// (Abkevich et al. 2012).
    /// </param>
    /// <param name="Tai">
    /// Telomeric allelic-imbalance score (NtAI): number of allelic-imbalance regions that extend to a
    /// sub-telomere but do not cross the centromere (Birkbak et al. 2012).
    /// </param>
    /// <param name="Lst">
    /// Large-scale state-transition score: number of chromosomal breaks between adjacent regions each
    /// ≥ 10 Mb after filtering regions &lt; 3 Mb (Popova et al. 2012).
    /// </param>
    public readonly record struct HrdComponents(int Loh, int Tai, int Lst);

    /// <summary>Result of a combined HRD determination from the three genomic-scar component counts.</summary>
    /// <param name="Components">The LOH / TAI / LST component counts that were summed.</param>
    /// <param name="Score">Combined HRD score = LOH + TAI + LST (unweighted sum).</param>
    /// <param name="Status">HRD-high when <paramref name="Score"/> ≥ <see cref="HrdHighScoreThreshold"/>, else HRD-negative.</param>
    public readonly record struct HrdResult(HrdComponents Components, int Score, HrdStatus Status);

    /// <summary>
    /// Computes the combined HRD score as the unweighted sum of the three genomic-scar component counts:
    /// score = <paramref name="loh"/> + <paramref name="tai"/> + <paramref name="lst"/>. Source: Telli ML
    /// et al. (2016), Clin Cancer Res 22(15):3764–3773 — the "combined homologous recombination deficiency
    /// (HRD) score, an unweighted sum of LOH, TAI, and LST scores". The components are non-negative event
    /// counts (LOH regions / telomeric allelic imbalances / large-scale state transitions), so each must be ≥ 0.
    /// </summary>
    /// <param name="loh">HRD-LOH component count (Abkevich et al. 2012); must be ≥ 0.</param>
    /// <param name="tai">Telomeric allelic-imbalance (NtAI) component count (Birkbak et al. 2012); must be ≥ 0.</param>
    /// <param name="lst">Large-scale state-transition component count (Popova et al. 2012); must be ≥ 0.</param>
    /// <returns>The combined HRD score (LOH + TAI + LST).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Any component is negative.</exception>
    public static int CalculateHRDScore(int loh, int tai, int lst)
    {
        if (loh < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(loh), loh, "HRD component counts must be ≥ 0.");
        }

        if (tai < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(tai), tai, "HRD component counts must be ≥ 0.");
        }

        if (lst < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(lst), lst, "HRD component counts must be ≥ 0.");
        }

        // Telli et al. (2016): HRD score = unweighted sum of LOH + TAI + LST. Sum in a wider type so
        // that extreme (near-Int32.MaxValue) component counts can NEVER wrap to a negative score
        // (INV-04: score ≥ 0). A sum that genuinely exceeds the int range is out of the documented
        // contract and is rejected rather than silently overflowing.
        long sum = (long)loh + tai + lst;
        if (sum > int.MaxValue)
        {
            throw new ArgumentOutOfRangeException(
                nameof(loh),
                sum,
                "Combined HRD score (loh + tai + lst) exceeds Int32.MaxValue.");
        }

        return (int)sum;
    }

    /// <summary>
    /// Classifies a combined HRD score as <see cref="HrdStatus.HrdHigh"/> when it is at or above the
    /// myChoice/Telli 2016 cutoff (≥ <see cref="HrdHighScoreThreshold"/> = 42; boundary inclusive),
    /// otherwise <see cref="HrdStatus.HrdNegative"/>. Source: Telli ML et al. (2016), Clin Cancer Res
    /// 22(15):3764–3773 — "HR deficiency, defined as HRD score ≥42".
    /// </summary>
    /// <param name="score">Combined HRD score (LOH + TAI + LST); must be ≥ 0.</param>
    /// <returns><see cref="HrdStatus.HrdHigh"/> if score ≥ 42, else <see cref="HrdStatus.HrdNegative"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is negative.</exception>
    public static HrdStatus ClassifyHRDStatus(int score)
    {
        if (score < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(score), score, "HRD score must be ≥ 0.");
        }

        return score >= HrdHighScoreThreshold ? HrdStatus.HrdHigh : HrdStatus.HrdNegative;
    }

    /// <summary>
    /// End-to-end HRD determination from the three genomic-scar component counts: sums them into the
    /// combined HRD score (<see cref="CalculateHRDScore(int,int,int)"/>) and classifies it against the
    /// myChoice/Telli 2016 cutoff (<see cref="ClassifyHRDStatus"/>). The three counts are produced upstream
    /// from segmented copy-number/allelic data per Abkevich et al. (2012) (LOH), Birkbak et al. (2012) (TAI),
    /// and Popova et al. (2012) (LST).
    /// </summary>
    /// <param name="components">The LOH / TAI / LST component counts.</param>
    /// <returns>The components, the combined HRD score, and the HRD status.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Any component count is negative.</exception>
    public static HrdResult DetectHRD(HrdComponents components)
    {
        int score = CalculateHRDScore(components.Loh, components.Tai, components.Lst);
        return new HrdResult(components, score, ClassifyHRDStatus(score));
    }

    /// <summary>
    /// End-to-end HRD determination that <b>derives the HRD-LOH component directly from allele-specific
    /// copy-number segments</b> (via <see cref="DetectLOH(IEnumerable{AlleleSpecificSegment})"/>, the
    /// Abkevich et al. 2012 / scarHRD <c>calc.hrd</c> rule: LOH regions &gt; 15 Mb that do not span a whole
    /// chromosome), then sums it with the caller-supplied telomeric-allelic-imbalance (<paramref name="tai"/>)
    /// and large-scale-state-transition (<paramref name="lst"/>) counts into the combined HRD score and
    /// classifies it against the myChoice/Telli 2016 cutoff (≥ <see cref="HrdHighScoreThreshold"/>).
    /// <para>
    /// This overload derives only the LOH component and accepts caller-supplied TAI and LST (e.g. when those
    /// were computed by an external pipeline). To derive all three components — LOH, TAI and LST — directly
    /// from the segments, use <see cref="DetectHRD(IEnumerable{AlleleSpecificSegment}, ReferenceGenome)"/>.
    /// </para>
    /// Source: Telli ML et al. (2016), Clin Cancer Res 22(15):3764–3773 ("unweighted sum of LOH, TAI, and
    /// LST scores"; "HRD score ≥42"); Abkevich et al. (2012) for the derived LOH component.
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments the HRD-LOH count is derived from. Must not be null.</param>
    /// <param name="tai">Caller-supplied telomeric-allelic-imbalance (NtAI) count (Birkbak et al. 2012); must be ≥ 0.</param>
    /// <param name="lst">Caller-supplied large-scale-state-transition count (Popova et al. 2012); must be ≥ 0.</param>
    /// <returns>The components (LOH derived, TAI/LST as supplied), the combined HRD score, and the HRD status.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tai"/> or <paramref name="lst"/> is negative.</exception>
    public static HrdResult DetectHRD(IEnumerable<AlleleSpecificSegment> segments, int tai, int lst)
    {
        ArgumentNullException.ThrowIfNull(segments);

        int loh = DetectLOH(segments).Score;
        return DetectHRD(new HrdComponents(loh, tai, lst));
    }

    #endregion


    #region Loss of heterozygosity (ONCO-LOH-001)

    /// <summary>
    /// Minimum LOH-region length (in base pairs) for a region to be counted toward the HRD-LOH score.
    /// The size comparison is strict (length must be &gt; this value). Source: Abkevich V et al. (2012),
    /// Br J Cancer 107(10):1776–1782 (PMID 23047548) — HRD-LOH counts LOH regions of "intermediate size"
    /// (&gt; 15 Mb and &lt; whole chromosome); reproduced in the scarHRD reference implementation as
    /// <c>sizelimitLOH = 15e6</c> with the filter <c>segLOH[,4]-segLOH[,3] &gt; sizelimit1</c>
    /// (https://github.com/sztup/scarHRD/blob/master/R/calc.hrd.R). 15 Mb = 15,000,000 bp.
    /// </summary>
    public const long HrdLohMinRegionLengthBp = 15_000_000L;

    /// <summary>
    /// An allele-specific copy-number segment: the unit of input for LOH detection, mirroring the scarHRD
    /// segmentation table (chromosome, start, end, major-allele CN, minor-allele CN). Source: scarHRD
    /// <c>scar_score.R</c> input columns (chromosome / start / end / A-allele CN / B-allele CN).
    /// </summary>
    /// <param name="Chromosome">Chromosome identifier (e.g. "1", "chrX"). Used to group segments.</param>
    /// <param name="Start">Segment start coordinate (bp). Must satisfy <see cref="End"/> &gt; <see cref="Start"/>.</param>
    /// <param name="End">Segment end coordinate (bp). Segment length = <see cref="End"/> − <see cref="Start"/>.</param>
    /// <param name="MajorCopyNumber">Major-allele copy number (A allele, the larger of the two). Must be ≥ 0.</param>
    /// <param name="MinorCopyNumber">Minor-allele copy number (B allele, the smaller of the two). Must be ≥ 0.</param>
    public readonly record struct AlleleSpecificSegment(
        string Chromosome,
        long Start,
        long End,
        int MajorCopyNumber,
        int MinorCopyNumber)
    {
        /// <summary>Segment length in base pairs, computed as End − Start (per scarHRD: <c>seg[,4]-seg[,3]</c>).</summary>
        public long Length => End - Start;
    }

    /// <summary>A copy-number segment that was counted as an HRD-LOH region.</summary>
    /// <param name="Chromosome">Chromosome of the region.</param>
    /// <param name="Start">Region start (bp).</param>
    /// <param name="End">Region end (bp).</param>
    /// <param name="Length">Region length (bp) = End − Start.</param>
    public readonly record struct LohRegion(string Chromosome, long Start, long End, long Length);

    /// <summary>Result of an HRD-LOH determination over a set of allele-specific segments.</summary>
    /// <param name="Regions">The qualifying LOH regions that were counted.</param>
    /// <param name="Score">HRD-LOH score = number of qualifying regions (= <c>Regions.Count</c>).</param>
    public readonly record struct LohResult(IReadOnlyList<LohRegion> Regions, int Score);

    /// <summary>
    /// Tests whether a single segment exhibits loss of heterozygosity: the minor-allele copy number is 0
    /// while the major-allele copy number is non-zero. Source: scarHRD <c>calc.hrd.R</c> —
    /// <c>segLOH &lt;- segSamp[segSamp[,nB] == 0 &amp; segSamp[,nA] != 0,]</c> (minor lost, major retained).
    /// A homozygous deletion (both alleles 0) is therefore NOT LOH.
    /// </summary>
    private static bool IsLohSegment(in AlleleSpecificSegment segment)
        => segment.MinorCopyNumber == 0 && segment.MajorCopyNumber != 0;

    /// <summary>
    /// Detects HRD-associated loss-of-heterozygosity regions from allele-specific copy-number segments and
    /// returns both the qualifying regions and the HRD-LOH score. The score is the number of LOH regions
    /// longer than 15 Mb (<see cref="HrdLohMinRegionLengthBp"/>, strict) that do not span a whole chromosome.
    /// Algorithm (Abkevich et al. 2012; scarHRD <c>calc.hrd</c>):
    /// <list type="number">
    /// <item><description>Group segments by chromosome.</description></item>
    /// <item><description>Mark a chromosome as whole-chromosome-LOH when every one of its segments is LOH
    /// (minor == 0); regions on such chromosomes are excluded (Abkevich: "&lt; whole chromosome").</description></item>
    /// <item><description>Merge adjacent same-LOH-state segments (oncoscanR <c>score_loh</c> merge step) so a
    /// long LOH region split into pieces counts once.</description></item>
    /// <item><description>Keep LOH regions whose length is strictly greater than 15 Mb, excluding
    /// whole-chromosome-LOH chromosomes; the count of survivors is the HRD-LOH score.</description></item>
    /// </list>
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments. Must not be null.</param>
    /// <returns>The qualifying LOH regions and the HRD-LOH score.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    public static LohResult DetectLOH(IEnumerable<AlleleSpecificSegment> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);

        var regions = new List<LohRegion>();

        foreach (var group in GroupValidatedByChromosome(segments))
        {
            // scarHRD: a chromosome whose every segment is LOH (minor == 0) is "global LOH" → excluded.
            bool wholeChromosomeLoh = group.All(static s => s.MinorCopyNumber == 0);
            if (wholeChromosomeLoh)
            {
                continue;
            }

            foreach (var merged in MergeAdjacentSameState(group))
            {
                // scarHRD: LOH region kept when minor == 0, major != 0, and length strictly > 15 Mb.
                if (IsLohSegment(merged) && merged.Length > HrdLohMinRegionLengthBp)
                {
                    regions.Add(new LohRegion(merged.Chromosome, merged.Start, merged.End, merged.Length));
                }
            }
        }

        return new LohResult(regions, regions.Count);
    }

    /// <summary>
    /// Computes the HRD-LOH score (number of qualifying LOH regions &gt; 15 Mb that do not span a whole
    /// chromosome) from allele-specific copy-number segments. Convenience wrapper over
    /// <see cref="DetectLOH(IEnumerable{AlleleSpecificSegment})"/>. Source: Abkevich et al. (2012) /
    /// scarHRD <c>calc.hrd</c>.
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments. Must not be null.</param>
    /// <returns>The HRD-LOH score (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    public static int CalculateHrdLohScore(IEnumerable<AlleleSpecificSegment> segments)
        => DetectLOH(segments).Score;

    /// <summary>
    /// Computes the length-weighted fraction of a single chromosome that lies under loss of heterozygosity:
    /// (total length of LOH segments on the chromosome) ÷ (total covered length on the chromosome). The
    /// result is in [0, 1] (Registry invariant). A segment is LOH per the Abkevich/scarHRD rule
    /// (minor == 0 &amp; major != 0). Unlike the HRD-LOH score, this fraction applies no 15 Mb size filter and
    /// no whole-chromosome exclusion — it is a raw per-chromosome LOH burden. If the chromosome has no
    /// covered length (absent from <paramref name="segments"/>), the fraction is 0.
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments. Must not be null.</param>
    /// <param name="chromosome">The chromosome identifier to score. Must not be null.</param>
    /// <returns>The LOH fraction of the chromosome in [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or <paramref name="chromosome"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    public static double CalculateLOHFraction(IEnumerable<AlleleSpecificSegment> segments, string chromosome)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ArgumentNullException.ThrowIfNull(chromosome);

        long totalLength = 0;
        long lohLength = 0;

        foreach (var segment in segments)
        {
            ValidateSegment(segment);
            if (!string.Equals(segment.Chromosome, chromosome, StringComparison.Ordinal))
            {
                continue;
            }

            totalLength += segment.Length;
            if (IsLohSegment(segment))
            {
                lohLength += segment.Length;
            }
        }

        if (totalLength == 0)
        {
            return 0.0;
        }

        return (double)lohLength / totalLength;
    }

    /// <summary>
    /// Groups validated segments by chromosome (each segment validated as it is read), preserving first-seen
    /// chromosome order. Validation here means LOH detection is order-independent in count (INV-6).
    /// </summary>
    private static IEnumerable<IReadOnlyList<AlleleSpecificSegment>> GroupValidatedByChromosome(
        IEnumerable<AlleleSpecificSegment> segments)
    {
        var byChromosome = new Dictionary<string, List<AlleleSpecificSegment>>(StringComparer.Ordinal);
        var order = new List<string>();

        foreach (var segment in segments)
        {
            ValidateSegment(segment);
            if (!byChromosome.TryGetValue(segment.Chromosome, out var list))
            {
                list = new List<AlleleSpecificSegment>();
                byChromosome[segment.Chromosome] = list;
                order.Add(segment.Chromosome);
            }

            list.Add(segment);
        }

        foreach (var chromosome in order)
        {
            yield return byChromosome[chromosome];
        }
    }

    /// <summary>
    /// Merges adjacent segments that share the same LOH state into single regions, after sorting by start.
    /// Two segments are merged when they share the same LOH/non-LOH state and are adjacent or overlapping
    /// (gap ≤ 1 bp). Source: oncoscanR <c>score_loh</c> — "merges overlapping or adjacent LOH segments
    /// (separated by 1bp)"; analogous to scarHRD's <c>shrink.seg.ai.wrapper</c> state merge.
    /// </summary>
    private static IReadOnlyList<AlleleSpecificSegment> MergeAdjacentSameState(
        IReadOnlyList<AlleleSpecificSegment> chromosomeSegments)
    {
        var sorted = chromosomeSegments.OrderBy(static s => s.Start).ThenBy(static s => s.End).ToList();
        var merged = new List<AlleleSpecificSegment>();

        foreach (var segment in sorted)
        {
            if (merged.Count == 0)
            {
                merged.Add(segment);
                continue;
            }

            var last = merged[^1];
            bool sameState = IsLohSegment(last) == IsLohSegment(segment);
            // Adjacent/overlapping = the next segment starts at most 1 bp after the previous one ends.
            const long MaxMergeGapBp = 1L;
            bool adjacent = segment.Start - last.End <= MaxMergeGapBp;

            if (sameState && adjacent)
            {
                long newEnd = Math.Max(last.End, segment.End);
                // Preserve LOH state: keep the major/minor of the LOH side when merging an LOH run.
                merged[^1] = last with { End = newEnd };
            }
            else
            {
                merged.Add(segment);
            }
        }

        return merged;
    }

    /// <summary>Validates a segment: positive length and non-negative copy numbers.</summary>
    private static void ValidateSegment(in AlleleSpecificSegment segment)
    {
        if (segment.End <= segment.Start)
        {
            throw new ArgumentException(
                $"Segment on '{segment.Chromosome}' must have End > Start (got Start={segment.Start}, End={segment.End}).",
                nameof(segment));
        }

        if (segment.MajorCopyNumber < 0 || segment.MinorCopyNumber < 0)
        {
            throw new ArgumentException(
                $"Segment on '{segment.Chromosome}' must have non-negative copy numbers " +
                $"(got Major={segment.MajorCopyNumber}, Minor={segment.MinorCopyNumber}).",
                nameof(segment));
        }
    }

    #endregion


    #region HRD-TAI and HRD-LST derivation from segments (ONCO-HRD-001)

    /// <summary>
    /// Minimum segment length (bp) for a segment to be retained before HRD-TAI assignment. Source: scarHRD
    /// <c>calc.ai_new</c> — <c>seg &lt;- seg[seg[,4]-seg[,3] &gt;= min.size,]</c> with default
    /// <c>min.size=1e6</c> (https://github.com/sztup/scarHRD/blob/master/R/calc.ai_new.R). 1 Mb.
    /// </summary>
    private const long HrdTaiMinSegmentLengthBp = 1_000_000L;

    /// <summary>
    /// Minimum segment length (bp) below which an arm segment is smoothed away before HRD-LST counting, and
    /// also the maximum gap between two large segments that still counts as one transition. Source: Popova
    /// et al. (2012), Cancer Res 72(21):5454 — LSTs are counted "after smoothing and filtering &lt;3 Mb"
    /// small-scale variation; scarHRD <c>calc.lst</c> uses <c>3e6</c> for both the smoothing filter
    /// (<c>(arm[,4]-arm[,3]) &lt; 3e6</c>) and the adjacency gap (<c>(seg[k,3]-seg[k-1,4]) &lt; 3e6</c>). 3 Mb.
    /// </summary>
    private const long HrdLstSmoothingLengthBp = 3_000_000L;

    /// <summary>
    /// Minimum length (bp) for a chromosome-arm segment to qualify as one side of an LST break. Source:
    /// Popova et al. (2012) — a break "between adjacent regions each of at least 10 megabases"; scarHRD
    /// <c>calc.lst</c> flags <c>(arm[,4]-arm[,3]) &gt;= 10e6</c>. 10 Mb.
    /// </summary>
    private const long HrdLstLargeSegmentLengthBp = 10_000_000L;

    /// <summary>
    /// Per-chromosome centromere boundaries: <c>Start</c> = centromere start (p-arm boundary, scarHRD
    /// <c>chrominfo[i,2]</c>); <c>End</c> = centromere end (q-arm boundary, scarHRD <c>chrominfo[i,3]</c>).
    /// </summary>
    private readonly record struct Centromere(long Start, long End);

    /// <summary>
    /// Centromere boundaries (chromosomes 1–22), GRCh38 / hg38. Embedded published reference data:
    /// the two UCSC cytoBand <c>acen</c> bands per chromosome — centromere start = p11 <c>acen</c>
    /// <c>chromStart</c>, centromere end = q11 <c>acen</c> <c>chromEnd</c>. Source: UCSC Genome Browser
    /// cytoBand track for hg38 (https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand,
    /// retrieved 2026-06-23), cross-verified against the NCBI GRC modeled-centromere table
    /// (https://www.ncbi.nlm.nih.gov/grc/human — CEN1 122,026,460–125,184,587 ≈ chr1 121.7–125.1 Mb;
    /// CEN21 10,864,561–12,915,808 ≈ chr21 10.9–13.0 Mb, agreeing to cytoband resolution). scarHRD's
    /// <c>chrominfo_grch38</c> derives from the same cytoBand <c>acen</c> regions. Indexed 0 = chr1 … 21 = chr22.
    /// </summary>
    private static readonly Centromere[] GRCh38Centromeres =
    {
        new(121_700_000L, 125_100_000L), // chr1
        new(91_800_000L, 96_000_000L),   // chr2
        new(87_800_000L, 94_000_000L),   // chr3
        new(48_200_000L, 51_800_000L),   // chr4
        new(46_100_000L, 51_400_000L),   // chr5
        new(58_500_000L, 62_600_000L),   // chr6
        new(58_100_000L, 62_100_000L),   // chr7
        new(43_200_000L, 47_200_000L),   // chr8
        new(42_200_000L, 45_500_000L),   // chr9
        new(38_000_000L, 41_600_000L),   // chr10
        new(51_000_000L, 55_800_000L),   // chr11
        new(33_200_000L, 37_800_000L),   // chr12
        new(16_500_000L, 18_900_000L),   // chr13
        new(16_100_000L, 18_200_000L),   // chr14
        new(17_500_000L, 20_500_000L),   // chr15
        new(35_300_000L, 38_400_000L),   // chr16
        new(22_700_000L, 27_400_000L),   // chr17
        new(15_400_000L, 21_500_000L),   // chr18
        new(24_200_000L, 28_100_000L),   // chr19
        new(25_700_000L, 30_400_000L),   // chr20
        new(10_900_000L, 13_000_000L),   // chr21
        new(13_700_000L, 17_400_000L),   // chr22
    };

    /// <summary>
    /// Centromere boundaries (chromosomes 1–22), GRCh37 / hg19. Embedded published reference data: the two
    /// UCSC cytoBand <c>acen</c> bands per chromosome (centromere start = p11 <c>acen</c> <c>chromStart</c>,
    /// end = q11 <c>acen</c> <c>chromEnd</c>). Source: UCSC Genome Browser cytoBand track for hg19
    /// (https://api.genome.ucsc.edu/getData/track?genome=hg19;track=cytoBand, retrieved 2026-06-23).
    /// Indexed 0 = chr1 … 21 = chr22.
    /// </summary>
    private static readonly Centromere[] GRCh37Centromeres =
    {
        new(121_500_000L, 128_900_000L), // chr1
        new(90_500_000L, 96_800_000L),   // chr2
        new(87_900_000L, 93_900_000L),   // chr3
        new(48_200_000L, 52_700_000L),   // chr4
        new(46_100_000L, 50_700_000L),   // chr5
        new(58_700_000L, 63_300_000L),   // chr6
        new(58_000_000L, 61_700_000L),   // chr7
        new(43_100_000L, 48_100_000L),   // chr8
        new(47_300_000L, 50_700_000L),   // chr9
        new(38_000_000L, 42_300_000L),   // chr10
        new(51_600_000L, 55_700_000L),   // chr11
        new(33_300_000L, 38_200_000L),   // chr12
        new(16_300_000L, 19_500_000L),   // chr13
        new(16_100_000L, 19_100_000L),   // chr14
        new(15_800_000L, 20_700_000L),   // chr15
        new(34_600_000L, 38_600_000L),   // chr16
        new(22_200_000L, 25_800_000L),   // chr17
        new(15_400_000L, 19_000_000L),   // chr18
        new(24_400_000L, 28_600_000L),   // chr19
        new(25_600_000L, 29_400_000L),   // chr20
        new(10_900_000L, 14_300_000L),   // chr21
        new(12_200_000L, 17_900_000L),   // chr22
    };

    /// <summary>
    /// Returns the embedded per-chromosome centromere boundary table (chromosomes 1–22) for a reference
    /// assembly, indexed 0 = chr1 … 21 = chr22. Source: UCSC cytoBand <c>acen</c> bands (see
    /// <see cref="GRCh38Centromeres"/> / <see cref="GRCh37Centromeres"/>).
    /// </summary>
    /// <param name="genome">The reference assembly.</param>
    /// <returns>The 22-element centromere table for the assembly.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    private static Centromere[] GetCentromeres(ReferenceGenome genome) => genome switch
    {
        ReferenceGenome.GRCh38 => GRCh38Centromeres,
        ReferenceGenome.GRCh37 => GRCh37Centromeres,
        _ => throw new ArgumentOutOfRangeException(nameof(genome), genome, "Unknown reference genome."),
    };

    /// <summary>
    /// Returns the centromere (p-arm end / q-arm start) for an autosome of <paramref name="genome"/>, or
    /// <c>false</c> for sex chromosomes / contigs / non-autosomes — which are excluded from TAI and LST,
    /// per scarHRD (LST skips chr23/24/X/Y; the centromere table is autosome-only).
    /// </summary>
    private static bool TryGetCentromere(string chromosome, ReferenceGenome genome, out Centromere centromere)
    {
        centromere = default;
        if (!TryGetAutosomeNumber(chromosome, out int number))
        {
            return false;
        }

        centromere = GetCentromeres(genome)[number - 1];
        return true;
    }

    /// <summary>
    /// A segment is in allelic imbalance when its two allele copy numbers differ (major ≠ minor). Source:
    /// Birkbak et al. (2012) — "Allelic imbalance was defined as any time the copy number of the two alleles
    /// were not equal, and at least one allele was present"; scarHRD <c>calc.ai_new</c> even-ploidy path
    /// <c>AI &lt;- c(0,2)[match(seg[,7]==seg[,8], ...)]</c> (col 7 = major, col 8 = minor). A balanced
    /// homozygous deletion (0,0) is not imbalance (both alleles absent), consistent with major == minor.
    /// </summary>
    private static bool IsAllelicImbalance(in AlleleSpecificSegment segment)
        => segment.MajorCopyNumber != segment.MinorCopyNumber;

    /// <summary>
    /// Computes the HRD-TAI (telomeric allelic-imbalance / NtAI) score from allele-specific copy-number
    /// segments: the number of allelic-imbalance regions that extend to a sub-telomere but do not cross the
    /// centromere. Source: Birkbak et al. (2012), Cancer Discov 2(4):366; scarHRD <c>calc.ai_new</c>
    /// (even-ploidy path). Algorithm, per chromosome (autosomes only; the centromere table is autosome-only):
    /// <list type="number">
    /// <item><description>Drop segments shorter than 1 Mb (<see cref="HrdTaiMinSegmentLengthBp"/>), then merge
    /// adjacent same-allele-state segments (scarHRD <c>shrink.seg.ai</c>).</description></item>
    /// <item><description>A chromosome with a single remaining segment that is imbalanced is whole-chromosome
    /// AI — NOT telomeric (scarHRD AI = 3).</description></item>
    /// <item><description>The FIRST segment counts as telomeric when it is imbalanced and its <b>end</b> is
    /// strictly before the centromere start (it touches the p-telomere and stays on the p-arm).</description></item>
    /// <item><description>The LAST segment counts as telomeric when it is imbalanced and its <b>start</b> is
    /// strictly after the centromere end (it touches the q-telomere and stays on the q-arm).</description></item>
    /// </list>
    /// The score is the number of telomeric-AI regions summed over autosomes.
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments. Must not be null.</param>
    /// <param name="genome">Reference assembly supplying the centromere table (default GRCh38).</param>
    /// <returns>The HRD-TAI score (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static int CalculateHrdTaiScore(
        IEnumerable<AlleleSpecificSegment> segments,
        ReferenceGenome genome = ReferenceGenome.GRCh38)
    {
        ArgumentNullException.ThrowIfNull(segments);

        int taiCount = 0;
        foreach (var group in GroupValidatedByChromosome(segments))
        {
            if (!TryGetCentromere(group[0].Chromosome, genome, out Centromere centromere))
            {
                continue; // non-autosome: excluded from TAI (scarHRD centromere table is autosome-only).
            }

            // scarHRD calc.ai_new: drop < 1 Mb, then merge adjacent equal-allele-state segments.
            var sized = new List<AlleleSpecificSegment>();
            foreach (var segment in group)
            {
                if (segment.Length >= HrdTaiMinSegmentLengthBp)
                {
                    sized.Add(segment);
                }
            }

            if (sized.Count == 0)
            {
                continue;
            }

            IReadOnlyList<AlleleSpecificSegment> merged = MergeAdjacentSameAlleleState(sized);

            // Single remaining segment with imbalance → whole-chromosome AI (scarHRD AI=3), not telomeric.
            if (merged.Count == 1)
            {
                continue;
            }

            // First segment: imbalanced AND end strictly before the centromere start → p-telomeric.
            var first = merged[0];
            if (IsAllelicImbalance(first) && first.End < centromere.Start)
            {
                taiCount++;
            }

            // Last segment: imbalanced AND start strictly after the centromere end → q-telomeric.
            var last = merged[^1];
            if (IsAllelicImbalance(last) && last.Start > centromere.End)
            {
                taiCount++;
            }
        }

        return taiCount;
    }

    /// <summary>
    /// Computes the HRD-LST (large-scale state-transition) score from allele-specific copy-number segments:
    /// the number of chromosomal breaks between two adjacent regions each ≥ 10 Mb, after smoothing away
    /// regions &lt; 3 Mb, counted per chromosome arm. Source: Popova et al. (2012), Cancer Res 72(21):5454;
    /// scarHRD <c>calc.lst</c> (<c>chr.arm='no'</c> path). Algorithm, per autosome (sex chromosomes excluded):
    /// <list type="number">
    /// <item><description>Split into the p-arm (segments with start ≤ centromere start, clamped to it) and the
    /// q-arm (segments with end ≥ centromere end, clamped to it); merge adjacent same-allele-state
    /// segments on each arm.</description></item>
    /// <item><description>Iteratively remove arm segments &lt; 3 Mb (<see cref="HrdLstSmoothingLengthBp"/>),
    /// re-merging after each removal, until none remain (scarHRD <c>while</c> smoothing loop).</description></item>
    /// <item><description>Flag each surviving arm segment as large when its length ≥ 10 Mb
    /// (<see cref="HrdLstLargeSegmentLengthBp"/>). Count one LST for each adjacent pair where BOTH are large
    /// AND the gap between them is &lt; 3 Mb.</description></item>
    /// </list>
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments. Must not be null.</param>
    /// <param name="genome">Reference assembly supplying the centromere table (default GRCh38).</param>
    /// <returns>The HRD-LST score (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static int CalculateHrdLstScore(
        IEnumerable<AlleleSpecificSegment> segments,
        ReferenceGenome genome = ReferenceGenome.GRCh38)
    {
        ArgumentNullException.ThrowIfNull(segments);

        int lstCount = 0;
        foreach (var group in GroupValidatedByChromosome(segments))
        {
            if (!TryGetCentromere(group[0].Chromosome, genome, out Centromere centromere))
            {
                continue; // sex chromosomes / non-autosomes excluded from LST (scarHRD).
            }

            // scarHRD calc.lst: chromosomes with < 2 segments cannot yield a transition.
            if (group.Count < 2)
            {
                continue;
            }

            // p-arm: segments whose start is at/left of the centromere start; clamp the last to the boundary.
            var pArm = new List<AlleleSpecificSegment>();
            foreach (var segment in group)
            {
                if (segment.Start <= centromere.Start)
                {
                    pArm.Add(segment);
                }
            }

            // q-arm: segments whose end is at/right of the centromere end; clamp the first to the boundary.
            var qArm = new List<AlleleSpecificSegment>();
            foreach (var segment in group)
            {
                if (segment.End >= centromere.End)
                {
                    qArm.Add(segment);
                }
            }

            lstCount += CountArmTransitions(pArm, clampLastEnd: centromere.Start, clampFirstStart: null);
            lstCount += CountArmTransitions(qArm, clampLastEnd: null, clampFirstStart: centromere.End);
        }

        return lstCount;
    }

    /// <summary>
    /// Counts LST breaks on a single chromosome arm: merge same-state, clamp to the centromere boundary,
    /// iteratively drop &lt; 3 Mb segments (re-merging), then count adjacent pairs that are both ≥ 10 Mb with a
    /// &lt; 3 Mb gap. Mirrors the per-arm block of scarHRD <c>calc.lst</c>.
    /// </summary>
    private static int CountArmTransitions(
        List<AlleleSpecificSegment> armSegments,
        long? clampLastEnd,
        long? clampFirstStart)
    {
        if (armSegments.Count == 0)
        {
            return 0;
        }

        var arm = new List<AlleleSpecificSegment>(MergeAdjacentSameAlleleState(armSegments));

        // scarHRD clamps the arm's inner edge to the centromere boundary (q.arm[1,3] / p.arm[last,4]).
        if (clampLastEnd is long pEnd)
        {
            arm[^1] = arm[^1] with { End = pEnd };
        }

        if (clampFirstStart is long qStart)
        {
            arm[0] = arm[0] with { Start = qStart };
        }

        // Iteratively remove < 3 Mb segments and re-merge until none remain (scarHRD while-loop).
        SmoothShortSegments(arm, clampLastEnd, clampFirstStart);

        if (arm.Count < 2)
        {
            return 0;
        }

        int breaks = 0;
        for (int k = 1; k < arm.Count; k++)
        {
            bool previousLarge = arm[k - 1].Length >= HrdLstLargeSegmentLengthBp;
            bool currentLarge = arm[k].Length >= HrdLstLargeSegmentLengthBp;
            long gap = arm[k].Start - arm[k - 1].End;
            if (previousLarge && currentLarge && gap < HrdLstSmoothingLengthBp)
            {
                breaks++;
            }
        }

        return breaks;
    }

    /// <summary>
    /// Iteratively removes arm segments shorter than 3 Mb and re-merges adjacent same-state segments after
    /// each removal, re-clamping the arm boundary, until no segment is shorter than 3 Mb. Mirrors the
    /// scarHRD <c>while(length(n.3mb) &gt; 0){ arm &lt;- arm[-n.3mb[1],]; arm &lt;- shrink.seg.ai(arm) }</c> loop.
    /// </summary>
    private static void SmoothShortSegments(
        List<AlleleSpecificSegment> arm,
        long? clampLastEnd,
        long? clampFirstStart)
    {
        while (true)
        {
            int shortIndex = -1;
            for (int i = 0; i < arm.Count; i++)
            {
                if (arm[i].Length < HrdLstSmoothingLengthBp)
                {
                    shortIndex = i;
                    break;
                }
            }

            if (shortIndex < 0)
            {
                return;
            }

            arm.RemoveAt(shortIndex);
            if (arm.Count == 0)
            {
                return;
            }

            var remerged = new List<AlleleSpecificSegment>(MergeAdjacentSameAlleleState(arm));
            if (clampLastEnd is long pEnd)
            {
                remerged[^1] = remerged[^1] with { End = pEnd };
            }

            if (clampFirstStart is long qStart)
            {
                remerged[0] = remerged[0] with { Start = qStart };
            }

            arm.Clear();
            arm.AddRange(remerged);
        }
    }

    /// <summary>
    /// Merges adjacent segments (after sorting by start) that share the same allele state — equal major AND
    /// minor copy numbers — into single regions. Source: scarHRD <c>shrink.seg.ai</c> —
    /// <c>new.chr[(j-1),7]==new.chr[j,7] &amp; new.chr[(j-1),8]==new.chr[j,8]</c> (cols 7,8 = A/B allele CN).
    /// Used by HRD-TAI and HRD-LST so a region split into equal-state pieces counts once.
    /// </summary>
    private static IReadOnlyList<AlleleSpecificSegment> MergeAdjacentSameAlleleState(
        IReadOnlyList<AlleleSpecificSegment> chromosomeSegments)
    {
        var sorted = chromosomeSegments.OrderBy(static s => s.Start).ThenBy(static s => s.End).ToList();
        var merged = new List<AlleleSpecificSegment>();

        foreach (var segment in sorted)
        {
            if (merged.Count == 0)
            {
                merged.Add(segment);
                continue;
            }

            var last = merged[^1];
            bool sameState = last.MajorCopyNumber == segment.MajorCopyNumber
                && last.MinorCopyNumber == segment.MinorCopyNumber;
            if (sameState)
            {
                merged[^1] = last with { End = Math.Max(last.End, segment.End) };
            }
            else
            {
                merged.Add(segment);
            }
        }

        return merged;
    }

    /// <summary>
    /// End-to-end HRD determination that derives <b>all three</b> genomic-scar components — HRD-LOH
    /// (<see cref="DetectLOH(IEnumerable{AlleleSpecificSegment})"/>, Abkevich et al. 2012), HRD-TAI
    /// (<see cref="CalculateHrdTaiScore"/>, Birkbak et al. 2012) and HRD-LST
    /// (<see cref="CalculateHrdLstScore"/>, Popova et al. 2012) — directly from allele-specific copy-number
    /// segments, then sums them into the combined HRD score and classifies it against the myChoice/Telli 2016
    /// cutoff (≥ <see cref="HrdHighScoreThreshold"/>). TAI and LST use the embedded per-build centromere
    /// table (UCSC cytoBand <c>acen</c>); both are computed on autosomes only.
    /// Source: scarHRD <c>scar_score</c> (<c>sum_HRD0 &lt;- res_lst + res_hrd + res_ai[1]</c>); Telli ML et al.
    /// (2016), Clin Cancer Res 22(15):3764–3773 ("unweighted sum of LOH, TAI, and LST scores"; "HRD score ≥42").
    /// </summary>
    /// <param name="segments">Allele-specific copy-number segments all three components are derived from. Must not be null.</param>
    /// <param name="genome">Reference assembly supplying the centromere table for TAI/LST (default GRCh38).</param>
    /// <returns>The derived components (LOH, TAI, LST), the combined HRD score, and the HRD status.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive length or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static HrdResult DetectHRD(
        IEnumerable<AlleleSpecificSegment> segments,
        ReferenceGenome genome = ReferenceGenome.GRCh38)
    {
        ArgumentNullException.ThrowIfNull(segments);

        // Materialize once so the three single-pass derivations all see the same segments.
        var materialized = segments as IReadOnlyList<AlleleSpecificSegment> ?? segments.ToList();

        int loh = DetectLOH(materialized).Score;
        int tai = CalculateHrdTaiScore(materialized, genome);
        int lst = CalculateHrdLstScore(materialized, genome);

        return DetectHRD(new HrdComponents(loh, tai, lst));
    }

    #endregion

}
