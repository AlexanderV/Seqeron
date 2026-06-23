using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

/// <summary>
/// Provides somatic mutation calling from tumor / matched-normal variant evidence.
/// <list type="bullet">
/// <item><description>
/// <b>CallSomaticMutations</b>: classifies each tumor variant as somatic or germline by comparing
/// the tumor allele frequency against the matched-normal allele frequency, per the somatic-state
/// definition S = {(f_t, f_n): f_t ≠ f_n} restricted to a homozygous-reference normal genotype
/// (Saunders et al. 2012, <i>Bioinformatics</i> 28(14):1811–1817; Kim et al. 2018, <i>Nat. Methods</i> 15:591–594).
/// </description></item>
/// <item><description>
/// <b>FilterGermlineVariants</b>: removes variants whose evidence in the matched normal indicates a
/// germline event, mirroring Mutect2's rule of skipping variants clearly present in the matched normal
/// (Benjamin et al. 2019, GATK Mutect2 documentation / bioRxiv 861054).
/// </description></item>
/// <item><description>
/// <b>CalculateSomaticScore</b>: a deterministic, monotone confidence score in [0, 1] derived from the
/// separation between tumor and normal allele frequencies (tumor presence and normal absence).
/// </description></item>
/// </list>
/// The classification is rule-based on observed allele fractions; full caller probability/LOD models
/// (Bayesian somatic likelihood) are out of scope and intentionally not implemented.
/// </summary>
public static class OncologyAnalyzer
{
    #region Constants

    /// <summary>
    /// Minimum tumor variant allele frequency (VAF) for a variant to be considered present in the tumor.
    /// Source: Yan et al. (2021), Scientific Reports 11:11640 — whole-exome sequencing has a mutation
    /// limit of detection at VAF of 5%; putative mutations called at ≤ 5% VAF are frequently sequencing
    /// errors. Default = 0.05 (5%).
    /// </summary>
    public const double DefaultTumorVafThreshold = 0.05;

    /// <summary>
    /// Maximum matched-normal VAF for a variant to be considered absent from the normal (i.e. the normal
    /// genotype is treated as homozygous reference, only sequencing noise present).
    /// Source: Saunders et al. (2012), Bioinformatics 28(14):1811–1817 — somatic calls are restricted to
    /// the homozygous-reference normal genotype P(S, G_n = ref/ref | D). The 0.01 (1%) noise ceiling
    /// reflects the baseline / sub-detection band used to declare a site germline-negative (Yan et al.
    /// 2021 report normal baseline thresholds in the sub-percent band). Default = 0.01 (1%).
    /// </summary>
    public const double DefaultNormalVafThreshold = 0.01;

    /// <summary>
    /// Standard-normal quantile z for a two-sided 95% confidence level used by the Wilson score interval.
    /// Source: Wilson E.B. (1927), via the Binomial proportion confidence interval specification, which
    /// states z₀.₀₅ = 1.96 for a 95% interval
    /// (https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval).
    /// </summary>
    private const double ZScore95 = 1.96;

    /// <summary>Default confidence level (0.95) for <see cref="CalculateVAFConfidenceInterval"/>.</summary>
    public const double DefaultVafConfidence = 0.95;

    /// <summary>
    /// Normal (germline) total copy number assumed by <see cref="AdjustVAFForPurity"/> — autosomal diploid.
    /// Source: Tarabichi et al. (2017), PMC5538405; CNAqc (Genome Biology 2024) — the VAF/purity relation
    /// fixes the normal contribution at n_tot,n = 2.
    /// </summary>
    private const double NormalDiploidCopyNumber = 2.0;

    /// <summary>
    /// Factor relating a clonal heterozygous somatic SNV VAF to tumor purity at a copy-neutral diploid locus:
    /// purity ρ = 2·VAF. Source: Antonello et al. (2024), CNAqc, Genome Biology 25:38 — the expected-VAF
    /// relation v = m·π / [2(1−π) + π·n_tot] with m = 1, n_tot = 2 gives v = π/2 (CNAqc: purity 60% ⇔ VAF 30%).
    /// </summary>
    private const double HeterozygousDiploidPurityFactor = 2.0;

    /// <summary>
    /// Mass in picograms of one haploid human genome equivalent. Source: Devonshire A.S. et al. (2014),
    /// <i>Anal. Bioanal. Chem.</i> (PMC4182654) — "one copy is … a single human haploid genome that is
    /// calculated as 3.3 pg." Corroborated by Alcaide et al. (2020), <i>Sci. Rep.</i> 10:12564
    /// ("1 ng of cfDNA roughly contains 303 haploid genome equivalents"; 1000 / 3.3 ≈ 303).
    /// </summary>
    private const double PicogramsPerHaploidGenome = 3.3;

    /// <summary>Picograms per nanogram (unit conversion).</summary>
    private const double PicogramsPerNanogram = 1000.0;

    /// <summary>
    /// Factor relating a clonal heterozygous somatic SNV VAF to ctDNA tumor fraction at a copy-neutral
    /// diploid locus: tumor fraction = 2·VAF. Same diploid-heterozygous identity as
    /// <see cref="HeterozygousDiploidPurityFactor"/> (CNAqc expected-VAF relation v = π/2; Antonello et al.
    /// 2024, <i>Genome Biology</i> 25:38). For ctDNA the tumour-derived fraction plays the role of purity.
    /// </summary>
    private const double TumorFractionFromVafFactor = 2.0;

    /// <summary>
    /// Default minimum detection probability for <see cref="IsCtDnaDetected"/> — the 95% sensitivity
    /// convention. Same 0.95 confidence level already used by <see cref="DefaultVafConfidence"/>.
    /// Source: Newman et al. (2014), <i>Nat. Med.</i> 20(5):548–554 report assay sensitivity/specificity at
    /// the conventional 95% operating point. ASSUMPTION: only the boolean detect flag depends on this value;
    /// the returned probability (<see cref="CtDnaDetectionProbability"/>) is source-exact.
    /// </summary>
    public const double DefaultCtDnaDetectionProbability = 0.95;

    /// <summary>
    /// Minimum expected number of mutant molecules (λ = n·d·k) for a variant to be considered detectable:
    /// at least one mutant molecule must be expected in the sequenced input. Source: the Poisson detection
    /// model (US Patent 11,085,084 B2; Avanzini et al. 2020, <i>Sci. Adv.</i> 6(50):eabc4308), under which
    /// detection in the low-burden regime is Poisson-limited by the expected mutant count λ.
    /// </summary>
    private const double MinExpectedMutantMolecules = 1.0;

    /// <summary>
    /// Default absolute expression z-score threshold above which a gene is an over/under-expression outlier.
    /// Source: cBioPortal FAQ (https://docs.cbioportal.org/user-guide/faq/) — "By default, samples with
    /// expression z-scores &gt;2 or &lt;-2 in any queried genes are considered altered." The rule is strict:
    /// z &gt; +2 ⇒ overexpressed, z &lt; −2 ⇒ underexpressed; |z| = 2 exactly is NOT an outlier. Default = 2.0.
    /// </summary>
    public const double DefaultExpressionOutlierThreshold = 2.0;

    #endregion

    #region Public types

    /// <summary>
    /// A single observed variant with the read evidence required for tumor/normal classification.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier.</param>
    /// <param name="Position">1-based reference position of the variant.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Alternate (observed) allele.</param>
    /// <param name="TumorAltReads">Number of reads supporting the alternate allele in the tumor.</param>
    /// <param name="TumorTotalReads">Total covering reads at this site in the tumor.</param>
    /// <param name="NormalAltReads">Number of reads supporting the alternate allele in the matched normal.</param>
    /// <param name="NormalTotalReads">Total covering reads at this site in the matched normal.</param>
    public readonly record struct VariantObservation(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        int TumorAltReads,
        int TumorTotalReads,
        int NormalAltReads,
        int NormalTotalReads);

    /// <summary>Classification outcome for a variant.</summary>
    public enum SomaticStatus
    {
        /// <summary>Present in tumor and absent in the matched normal: a somatic variant.</summary>
        Somatic,

        /// <summary>Present in both tumor and matched normal: a germline variant.</summary>
        Germline,

        /// <summary>Not present in the tumor above the detection threshold.</summary>
        NotDetected
    }

    /// <summary>Result of classifying one variant.</summary>
    /// <param name="Variant">The variant that was classified.</param>
    /// <param name="TumorVaf">Tumor allele frequency f_t = altReads / totalReads.</param>
    /// <param name="NormalVaf">Matched-normal allele frequency f_n (0 when no matched normal).</param>
    /// <param name="Status">Somatic / Germline / NotDetected classification.</param>
    /// <param name="SomaticScore">Confidence in [0, 1] that the variant is somatic.</param>
    public readonly record struct SomaticCall(
        VariantObservation Variant,
        double TumorVaf,
        double NormalVaf,
        SomaticStatus Status,
        double SomaticScore);

    /// <summary>
    /// A variant allele frequency point estimate with a Wilson score confidence interval.
    /// </summary>
    /// <param name="Vaf">Point estimate VAF = altReads / totalReads (the empirical allele fraction).</param>
    /// <param name="Lower">Lower bound of the Wilson score interval (≥ 0).</param>
    /// <param name="Upper">Upper bound of the Wilson score interval (≤ 1).</param>
    /// <param name="Confidence">Two-sided confidence level the interval was computed at (e.g. 0.95).</param>
    public readonly record struct VafConfidenceInterval(
        double Vaf,
        double Lower,
        double Upper,
        double Confidence);

    /// <summary>
    /// A clonal somatic mutation together with the allele-specific copy-number state required to estimate
    /// tumor purity by inverting the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot]
    /// (Antonello et al. 2024, <i>Genome Biology</i> 25:38).
    /// </summary>
    /// <param name="Vaf">Observed variant allele frequency v ∈ [0, 1] of a clonal somatic mutation.</param>
    /// <param name="Multiplicity">Mutation multiplicity m: number of tumour-genome copies carrying the mutant allele (≥ 1).</param>
    /// <param name="TumorTotalCopyNumber">Tumour total copy number n_tot = n_A + n_B at the locus (≥ 1).</param>
    public readonly record struct PurityVariant(double Vaf, int Multiplicity, int TumorTotalCopyNumber);

    #endregion

    #region CallSomaticMutations

    /// <summary>
    /// Classifies each tumor variant as somatic, germline, or not-detected by comparing the tumor and
    /// matched-normal allele frequencies. A variant is <see cref="SomaticStatus.Somatic"/> when its tumor
    /// VAF is at or above <paramref name="tumorVafThreshold"/> (present in tumor) and its matched-normal
    /// VAF is at or below <paramref name="normalVafThreshold"/> (absent from the normal, i.e. the normal is
    /// homozygous reference). This realizes the somatic state S = {(f_t, f_n): f_t ≠ f_n} restricted to a
    /// ref/ref normal genotype (Saunders et al. 2012; Kim et al. 2018).
    /// </summary>
    /// <param name="variants">Tumor variant observations with matched-normal read evidence.</param>
    /// <param name="tumorVafThreshold">Minimum tumor VAF for presence (default <see cref="DefaultTumorVafThreshold"/>).</param>
    /// <param name="normalVafThreshold">Maximum normal VAF for absence (default <see cref="DefaultNormalVafThreshold"/>).</param>
    /// <returns>One <see cref="SomaticCall"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is outside [0, 1].</exception>
    public static IReadOnlyList<SomaticCall> CallSomaticMutations(
        IEnumerable<VariantObservation> variants,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ValidateThreshold(tumorVafThreshold, nameof(tumorVafThreshold));
        ValidateThreshold(normalVafThreshold, nameof(normalVafThreshold));

        var calls = new List<SomaticCall>();
        foreach (var variant in variants)
        {
            calls.Add(Classify(variant, tumorVafThreshold, normalVafThreshold));
        }

        return calls;
    }

    /// <summary>
    /// Classifies a single variant. See <see cref="CallSomaticMutations(IEnumerable{VariantObservation}, double, double)"/>.
    /// </summary>
    public static SomaticCall Classify(
        VariantObservation variant,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ValidateThreshold(tumorVafThreshold, nameof(tumorVafThreshold));
        ValidateThreshold(normalVafThreshold, nameof(normalVafThreshold));

        double tumorVaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
        double normalVaf = CalculateVaf(variant.NormalAltReads, variant.NormalTotalReads);

        SomaticStatus status;
        if (tumorVaf < tumorVafThreshold)
        {
            // Not present in the tumor above the detection limit (Yan et al. 2021).
            status = SomaticStatus.NotDetected;
        }
        else if (normalVaf <= normalVafThreshold)
        {
            // Present in tumor, absent in normal (ref/ref) → somatic (Saunders et al. 2012).
            status = SomaticStatus.Somatic;
        }
        else
        {
            // Present in both tumor and normal → germline (Benjamin et al. 2019, Mutect2).
            status = SomaticStatus.Germline;
        }

        double score = status == SomaticStatus.Somatic
            ? CalculateSomaticScore(tumorVaf, normalVaf)
            : 0.0;

        return new SomaticCall(variant, tumorVaf, normalVaf, status, score);
    }

    #endregion

    #region FilterGermlineVariants

    /// <summary>
    /// Removes germline variants, returning only somatic calls. A variant is filtered out as germline
    /// when its matched-normal VAF exceeds <paramref name="normalVafThreshold"/>, mirroring Mutect2's rule
    /// of skipping variants clearly present in the matched normal (Benjamin et al. 2019). Variants not
    /// detected in the tumor (below <paramref name="tumorVafThreshold"/>) are also excluded.
    /// </summary>
    /// <param name="variants">Tumor variant observations with matched-normal read evidence.</param>
    /// <param name="tumorVafThreshold">Minimum tumor VAF for presence (default <see cref="DefaultTumorVafThreshold"/>).</param>
    /// <param name="normalVafThreshold">Maximum normal VAF for absence (default <see cref="DefaultNormalVafThreshold"/>).</param>
    /// <returns>The subset of variants classified as <see cref="SomaticStatus.Somatic"/>, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<SomaticCall> FilterGermlineVariants(
        IEnumerable<VariantObservation> variants,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);

        return CallSomaticMutations(variants, tumorVafThreshold, normalVafThreshold)
            .Where(c => c.Status == SomaticStatus.Somatic)
            .ToList();
    }

    #endregion

    #region CalculateSomaticScore

    /// <summary>
    /// Computes a deterministic somatic-confidence score in [0, 1] from a variant's read evidence. The
    /// score increases with tumor presence (f_t) and with the separation f_t − f_n between tumor and normal
    /// allele frequencies, reflecting the somatic criterion that f_t ≠ f_n with the normal at ref/ref
    /// (Saunders et al. 2012). When the normal VAF equals or exceeds the tumor VAF the score is 0
    /// (no somatic evidence).
    /// </summary>
    /// <param name="variant">The variant observation to score.</param>
    /// <returns>Somatic confidence in [0, 1].</returns>
    public static double CalculateSomaticScore(VariantObservation variant)
    {
        double tumorVaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
        double normalVaf = CalculateVaf(variant.NormalAltReads, variant.NormalTotalReads);
        return CalculateSomaticScore(tumorVaf, normalVaf);
    }

    /// <summary>
    /// Computes the somatic-confidence score from tumor and normal allele frequencies.
    /// score = f_t × (f_t − f_n) / f_t = f_t − f_n clamped to [0, 1], i.e. the allele-frequency separation
    /// scaled by tumor presence. Equivalent simplified form: max(0, f_t − f_n).
    /// </summary>
    private static double CalculateSomaticScore(double tumorVaf, double normalVaf)
    {
        // Separation between tumor and normal allele frequencies; 0 when the normal carries the allele
        // at or above the tumor level (no somatic signal). Bounded in [0, 1] since both VAFs are in [0, 1].
        double separation = tumorVaf - normalVaf;
        return separation > 0.0 ? separation : 0.0;
    }

    #endregion

    #region Variant Allele Frequency

    /// <summary>
    /// Computes the empirical variant allele frequency (VAF) at a locus as the fraction of covering reads
    /// that support the alternate allele: VAF = altReads / totalReads. This is the model-free allele
    /// fraction derived from per-allele read depths (GATK <c>AD</c> field: alt AD / Σ AD; samtools read
    /// counts), distinct from Mutect2's Bayesian <c>AF</c> estimate. A site with no coverage
    /// (totalReads == 0) returns 0, since an uncovered site provides no evidence of the allele.
    /// </summary>
    /// <param name="altReads">Reads supporting the alternate allele (≥ 0, ≤ <paramref name="totalReads"/>).</param>
    /// <param name="totalReads">Total covering reads at the locus (≥ 0).</param>
    /// <returns>VAF in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative, or altReads &gt; totalReads.</exception>
    public static double CalculateVAF(int altReads, int totalReads) => CalculateVaf(altReads, totalReads);

    /// <summary>
    /// Computes the empirical VAF and its Wilson score confidence interval for the underlying allele
    /// proportion. The Wilson score interval (Wilson 1927) for n trials with p̂ = altReads/totalReads is:
    /// <code>
    /// center = (p̂ + z²/(2n)) / (1 + z²/n)
    /// margin = (z / (1 + z²/n)) · √( p̂(1−p̂)/n + z²/(4n²) )
    /// interval = center ± margin
    /// </code>
    /// where z is the standard-normal quantile for the requested confidence (z = 1.96 for 95%). The
    /// interval is bounded within [0, 1] (no overshoot) and has non-zero width even at p̂ = 0 or 1,
    /// unlike the Wald interval. Source: Wilson E.B. (1927), JASA 22(158):209–212, via the Binomial
    /// proportion confidence interval specification.
    /// </summary>
    /// <param name="altReads">Reads supporting the alternate allele (≥ 0, ≤ <paramref name="totalReads"/>).</param>
    /// <param name="totalReads">Total covering reads at the locus (&gt; 0 for a defined interval).</param>
    /// <param name="confidence">Two-sided confidence level in (0, 1); default 0.95 (z = 1.96).</param>
    /// <returns>The VAF point estimate with its Wilson score interval.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A read count is invalid (see <see cref="CalculateVAF"/>), totalReads == 0, or confidence is outside (0, 1)
    /// at a level other than the supported 0.95.
    /// </exception>
    public static VafConfidenceInterval CalculateVAFConfidenceInterval(
        int altReads,
        int totalReads,
        double confidence = DefaultVafConfidence)
    {
        double vaf = CalculateVaf(altReads, totalReads);

        if (totalReads == 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalReads), "A confidence interval is undefined with zero coverage.");
        }

        double z = ZScoreFor(confidence);

        double n = totalReads;
        double pHat = vaf;
        double z2 = z * z;
        double denominator = 1.0 + z2 / n;
        double center = (pHat + z2 / (2.0 * n)) / denominator;
        double margin = (z / denominator) * Math.Sqrt(pHat * (1.0 - pHat) / n + z2 / (4.0 * n * n));

        // The Wilson interval is mathematically within [0, 1]; clamp guards against floating-point drift
        // at the exact boundaries p̂ = 0 (lower = 0) and p̂ = 1 (upper = 1).
        double lower = Math.Max(0.0, center - margin);
        double upper = Math.Min(1.0, center + margin);

        return new VafConfidenceInterval(vaf, lower, upper, confidence);
    }

    /// <summary>
    /// Adjusts an observed VAF for tumor purity and tumor-segment ploidy, recovering the per-tumour-copy
    /// mutant fraction (multiplicity × cancer cell fraction). Inverting the CNAqc expected-VAF relation
    /// v = (m·π) / (2(1−π) + π·n_tot) gives:
    /// <code>
    /// adjusted = vaf · (2(1−π) + π·ploidy) / π
    /// </code>
    /// where π = purity, ploidy = tumor total copy number n_tot, and the normal contribution is fixed at
    /// 2 (autosomal diploid). For a heterozygous somatic SNV in a diploid tumor (ploidy = 2) this reduces
    /// to adjusted = vaf / (π/2): e.g. observed VAF 0.4 at purity 0.8 ⇒ 1.0. Source: CNAqc
    /// (Genome Biology 2024); Tarabichi et al. (2017), PMC5538405.
    /// </summary>
    /// <param name="vaf">Observed VAF in [0, 1].</param>
    /// <param name="purity">Tumor purity π in (0, 1].</param>
    /// <param name="ploidy">Tumor total copy number n_tot at the locus (&gt; 0); 2 for a diploid segment.</param>
    /// <returns>The purity/ploidy-corrected mutant fraction (multiplicity × CCF).</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0, 1], purity ∉ (0, 1], or ploidy ≤ 0.</exception>
    public static double AdjustVAFForPurity(double vaf, double purity, double ploidy)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in the range [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in the range (0, 1]; correction divides by purity.");
        }

        if (double.IsNaN(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy (tumor total copy number) must be positive.");
        }

        // Weighted average total copies per cell: 2(1−π) from normal cells + π·ploidy from tumor cells.
        double averageCopiesPerCell = NormalDiploidCopyNumber * (1.0 - purity) + purity * ploidy;
        return vaf * averageCopiesPerCell / purity;
    }

    /// <summary>
    /// Estimates tumor purity ρ from the variant allele frequencies of clonal, heterozygous somatic SNVs
    /// at copy-neutral diploid loci. For such a variant (multiplicity m = 1, total copy number n_tot = 2),
    /// the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot] reduces to v = π/2, so the per-variant
    /// purity is ρ = 2·v (Antonello et al. 2024, <i>Genome Biology</i> 25:38; CNAqc reports purity 60% ⇔ VAF 30%).
    /// The per-variant estimates are aggregated by their median, which is robust to subclonal / outlier VAFs.
    /// </summary>
    /// <param name="variants">Observed clonal heterozygous somatic SNVs (each contributes ρ = 2·VAF).</param>
    /// <returns>Estimated tumor purity ρ ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (purity is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A variant has invalid read counts, or a VAF &gt; 0.5 (which would imply purity &gt; 1 under the diploid model).
    /// </exception>
    public static double EstimatePurityFromVAF(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var purities = new List<double>();
        foreach (VariantObservation variant in variants)
        {
            double vaf = CalculateVAF(variant.TumorAltReads, variant.TumorTotalReads);
            purities.Add(EstimatePurityFromVaf(vaf));
        }

        if (purities.Count == 0)
        {
            throw new ArgumentException("Cannot estimate purity from an empty variant set.", nameof(variants));
        }

        return Median(purities);
    }

    /// <summary>
    /// Estimates tumor purity ρ from a single clonal heterozygous somatic SNV VAF at a copy-neutral diploid
    /// locus using the closed form ρ = 2·v (the m = 1, n_tot = 2 special case of the CNAqc expected-VAF
    /// relation; Antonello et al. 2024, <i>Genome Biology</i> 25:38).
    /// </summary>
    /// <param name="vaf">Observed VAF v ∈ [0, 0.5] of the clonal heterozygous SNV.</param>
    /// <returns>Estimated tumor purity ρ = 2·v ∈ [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// vaf ∉ [0, 1], or vaf &gt; 0.5 (which would imply purity &gt; 1 under the diploid heterozygous model).
    /// </exception>
    public static double EstimatePurityFromVaf(double vaf)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in the range [0, 1].");
        }

        // ρ = 2·v; a heterozygous SNV in a diploid tumour cannot exceed VAF 0.5, so ρ > 1 is impossible input.
        double purity = HeterozygousDiploidPurityFactor * vaf;
        if (purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(vaf), vaf,
                "VAF > 0.5 implies purity > 1 under the heterozygous diploid model; the locus is not copy-neutral diploid heterozygous.");
        }

        return purity;
    }

    /// <summary>
    /// Estimates tumor purity ρ from clonal somatic mutations with known allele-specific copy-number state by
    /// inverting the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot] for π:
    /// <code>
    /// π = 2·v / [m + v·(2 − n_tot)]
    /// </code>
    /// where m is the mutation multiplicity and n_tot the tumour total copy number (Antonello et al. 2024,
    /// <i>Genome Biology</i> 25:38). The per-variant estimates are aggregated by their median.
    /// </summary>
    /// <param name="variants">Clonal somatic mutations with VAF, multiplicity m, and tumour total copy number n_tot.</param>
    /// <returns>Estimated tumor purity ρ ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (purity is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A variant has vaf ∉ [0, 1], multiplicity &lt; 1, n_tot &lt; 1, or yields a purity outside [0, 1].
    /// </exception>
    public static double EstimatePurity(IEnumerable<PurityVariant> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var purities = new List<double>();
        foreach (PurityVariant variant in variants)
        {
            purities.Add(EstimatePurityFromAlleleSpecificVaf(variant));
        }

        if (purities.Count == 0)
        {
            throw new ArgumentException("Cannot estimate purity from an empty variant set.", nameof(variants));
        }

        return Median(purities);
    }

    /// <summary>
    /// Inverts the CNAqc expected-VAF relation for a single allele-specific variant:
    /// π = 2·v / [m + v·(2 − n_tot)].
    /// </summary>
    private static double EstimatePurityFromAlleleSpecificVaf(in PurityVariant variant)
    {
        if (double.IsNaN(variant.Vaf) || variant.Vaf < 0.0 || variant.Vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf, "VAF must be in the range [0, 1].");
        }

        if (variant.Multiplicity < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Multiplicity, "Mutation multiplicity m must be at least 1.");
        }

        if (variant.TumorTotalCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.TumorTotalCopyNumber, "Tumour total copy number n_tot must be at least 1.");
        }

        // π = 2v / [m + v(2 − n_tot)], the algebraic inverse of v = mπ / [2(1−π) + π·n_tot].
        double denominator = variant.Multiplicity + variant.Vaf * (NormalDiploidCopyNumber - variant.TumorTotalCopyNumber);
        if (denominator <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf,
                "The (VAF, multiplicity, copy-number) combination does not correspond to a purity in [0, 1].");
        }

        double purity = NormalDiploidCopyNumber * variant.Vaf / denominator;
        if (purity < 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf,
                "The (VAF, multiplicity, copy-number) combination yields a purity outside [0, 1].");
        }

        return purity;
    }

    /// <summary>Median of a non-empty list of values (lower-mid average for even counts). Does not mutate the input.</summary>
    private static double Median(List<double> values)
    {
        double[] sorted = values.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        int mid = n / 2;
        return (n % 2 == 1) ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
    }

    #endregion

    #region Helpers

    /// <summary>
    /// Maps a two-sided confidence level to its standard-normal quantile z. Only the source-cited 95%
    /// level (z = 1.96, Wilson 1927) is supported; other levels would require additional cited z values.
    /// </summary>
    private static double ZScoreFor(double confidence)
    {
        if (double.IsNaN(confidence) || confidence <= 0.0 || confidence >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(confidence), confidence, "Confidence must be in the open interval (0, 1).");
        }

        if (Math.Abs(confidence - DefaultVafConfidence) < 1e-12)
        {
            return ZScore95;
        }

        throw new ArgumentOutOfRangeException(
            nameof(confidence), confidence,
            "Only the 0.95 confidence level (z = 1.96) is supported by the cited source.");
    }

    /// <summary>
    /// Computes a variant allele frequency f = altReads / totalReads. Returns 0 when there is no coverage
    /// (totalReads == 0) so an uncovered site is treated as the allele being absent.
    /// </summary>
    private static double CalculateVaf(int altReads, int totalReads)
    {
        if (altReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(altReads), "Alt read count cannot be negative.");
        }

        if (totalReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(totalReads), "Total read count cannot be negative.");
        }

        if (altReads > totalReads)
        {
            throw new ArgumentOutOfRangeException(nameof(altReads), "Alt read count cannot exceed total read count.");
        }

        return totalReads == 0 ? 0.0 : (double)altReads / totalReads;
    }

    private static void ValidateThreshold(double value, string paramName)
    {
        if (double.IsNaN(value) || value < 0.0 || value > 1.0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, "Threshold must be in the range [0, 1].");
        }
    }

    #endregion

    #region Driver Mutation Detection (20/20 rule)

    /// <summary>
    /// Fraction-of-mutations threshold of the Vogelstein 20/20 rule: a gene is classified as a driver
    /// only when more than this fraction of its mutations meet the oncogene or tumor-suppressor criterion.
    /// Source: Vogelstein B et al. (2013), Science 339(6127):1546–1558 — "more than 20%"; restated verbatim
    /// by Tokheim &amp; Karchin (2020), Bioinformatics 36(6):1712–1719 ("OGs have &gt;20% ... TSGs have &gt;20% ...").
    /// The comparison is strict (&gt;), so an exact 20% fraction is not sufficient. Value = 0.20.
    /// </summary>
    public const double DriverGeneFractionThreshold = 0.20;

    /// <summary>
    /// Minimum number of mutations at the same protein position for that position to count as
    /// <i>recurrent</i> (a hotspot). Source: Miller ML et al. (2017), Oncotarget 8(20):33321–33333 — a
    /// recurrent position requires "at least two mutations of the same class" at an identical location.
    /// Value = 2.
    /// </summary>
    public const int RecurrentPositionMinCount = 2;

    /// <summary>
    /// The functional consequence of a coding mutation, restricted to the categories the 20/20 rule
    /// distinguishes. Truncating categories are those listed by Schroeder MP et al. (2014),
    /// Bioinformatics 30(17):i549–i555 and Miller ML et al. (2017): nonsense (stop gain/loss), frameshift
    /// indels, and splice donor/acceptor mutations. Missense at recurrent positions drives the oncogene call.
    /// </summary>
    public enum MutationConsequence
    {
        /// <summary>Amino-acid-changing substitution (drives the oncogene criterion when recurrent).</summary>
        Missense,

        /// <summary>Premature/lost stop codon (nonsense) — truncating/inactivating.</summary>
        Nonsense,

        /// <summary>Insertion/deletion shifting the reading frame — truncating/inactivating.</summary>
        Frameshift,

        /// <summary>Mutation at a splice donor/acceptor site — truncating/inactivating.</summary>
        SpliceSite,

        /// <summary>Synonymous or other non-truncating, non-missense change (counts toward the denominator only).</summary>
        Other
    }

    /// <summary>The 20/20-rule role assigned to a gene from its mutation spectrum.</summary>
    public enum DriverGeneRole
    {
        /// <summary>&gt;20% of mutations are missense at recurrent positions (Vogelstein 2013).</summary>
        Oncogene,

        /// <summary>&gt;20% of mutations are truncating/inactivating (Vogelstein 2013).</summary>
        TumorSuppressor,

        /// <summary>Neither criterion exceeds 20% (or an exact tie): not classified as a driver gene.</summary>
        Ambiguous
    }

    /// <summary>
    /// A single coding mutation observed in a gene, reduced to the features the 20/20 rule needs.
    /// </summary>
    /// <param name="Gene">Gene symbol the mutation falls in.</param>
    /// <param name="ProteinPosition">1-based codon / amino-acid position used to detect recurrence.</param>
    /// <param name="Consequence">Functional consequence category.</param>
    public readonly record struct GeneMutation(string Gene, int ProteinPosition, MutationConsequence Consequence);

    /// <summary>
    /// The 20/20-rule classification of one gene, with the two criterion fractions that produced it.
    /// </summary>
    /// <param name="Gene">Gene symbol.</param>
    /// <param name="Role">Assigned <see cref="DriverGeneRole"/>.</param>
    /// <param name="TruncatingFraction">Fraction of mutations that are truncating/inactivating, in [0, 1].</param>
    /// <param name="RecurrentMissenseFraction">Fraction of mutations that are missense at a recurrent position, in [0, 1].</param>
    /// <param name="MutationCount">Total number of mutations considered (the denominator).</param>
    public readonly record struct DriverGeneClassification(
        string Gene,
        DriverGeneRole Role,
        double TruncatingFraction,
        double RecurrentMissenseFraction,
        int MutationCount);

    /// <summary>
    /// Classifies a single gene by the Vogelstein 20/20 rule from its observed coding mutations. The gene is
    /// an <see cref="DriverGeneRole.Oncogene"/> when more than 20% of its mutations are missense at recurrent
    /// positions (a position observed ≥ <see cref="RecurrentPositionMinCount"/> times), and a
    /// <see cref="DriverGeneRole.TumorSuppressor"/> when more than 20% of its mutations are truncating
    /// (nonsense, frameshift, or splice-site). If both criteria are met the dominant fraction decides; an exact
    /// tie or neither criterion met yields <see cref="DriverGeneRole.Ambiguous"/>.
    /// Source: Vogelstein et al. (2013); Tokheim &amp; Karchin (2020); Schroeder et al. (2014); Miller et al. (2017).
    /// </summary>
    /// <param name="mutations">All coding mutations observed in one gene.</param>
    /// <returns>The gene's 20/20-rule classification with its criterion fractions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static DriverGeneClassification ClassifyGene(IEnumerable<GeneMutation> mutations)
    {
        ArgumentNullException.ThrowIfNull(mutations);

        var list = mutations as IReadOnlyList<GeneMutation> ?? mutations.ToList();
        int total = list.Count;
        string gene = total > 0 ? list[0].Gene : string.Empty;

        if (total == 0)
        {
            return new DriverGeneClassification(gene, DriverGeneRole.Ambiguous, 0.0, 0.0, 0);
        }

        int truncating = list.Count(IsTruncating);
        int recurrentMissense = CountRecurrentMissense(list);

        double truncatingFraction = (double)truncating / total;
        double recurrentMissenseFraction = (double)recurrentMissense / total;

        bool isTsg = truncatingFraction > DriverGeneFractionThreshold;
        bool isOg = recurrentMissenseFraction > DriverGeneFractionThreshold;

        DriverGeneRole role;
        if (isTsg && isOg)
        {
            // Both criteria pass (atypical per Vogelstein 2013 — well-documented genes far surpass one
            // criterion). Resolve by the dominant signal; an exact tie is genuinely ambiguous.
            role = truncatingFraction > recurrentMissenseFraction ? DriverGeneRole.TumorSuppressor
                 : recurrentMissenseFraction > truncatingFraction ? DriverGeneRole.Oncogene
                 : DriverGeneRole.Ambiguous;
        }
        else if (isTsg)
        {
            role = DriverGeneRole.TumorSuppressor;
        }
        else if (isOg)
        {
            role = DriverGeneRole.Oncogene;
        }
        else
        {
            role = DriverGeneRole.Ambiguous;
        }

        return new DriverGeneClassification(gene, role, truncatingFraction, recurrentMissenseFraction, total);
    }

    /// <summary>
    /// Computes the 20/20-rule driver-signal score for a gene: the larger of its truncating fraction and its
    /// recurrent-missense fraction, in [0, 1]. This is the transparent, source-derived strength of the driver
    /// signal underlying <see cref="ClassifyGene"/>; it is NOT an external pathogenicity model (CADD/SIFT/
    /// PolyPhen), which are caller-supplied / not implemented. Source: Vogelstein et al. (2013) 20/20 rule.
    /// </summary>
    /// <param name="mutations">All coding mutations observed in one gene.</param>
    /// <returns>max(truncating fraction, recurrent-missense fraction), in [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static double ScoreDriverPotential(IEnumerable<GeneMutation> mutations)
    {
        var c = ClassifyGene(mutations);
        return Math.Max(c.TruncatingFraction, c.RecurrentMissenseFraction);
    }

    /// <summary>
    /// Tests whether a mutation's (gene, protein position) is present in a caller-supplied set of known
    /// cancer hotspot positions. The 20/20 rule treats recurrent positions as the activating signal of
    /// oncogenes (Miller et al. 2017); curated hotspot catalogs (COSMIC, Cancer Hotspots, OncoKB) are
    /// supplied by the caller rather than hardcoded, since they cannot be reproduced authoritatively here.
    /// </summary>
    /// <param name="mutation">The mutation to test.</param>
    /// <param name="knownHotspots">Set of known hotspots as (gene, protein position) pairs.</param>
    /// <returns><c>true</c> when (gene, position) is in <paramref name="knownHotspots"/>; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="knownHotspots"/> is null.</exception>
    public static bool MatchCancerHotspots(
        GeneMutation mutation,
        IReadOnlySet<(string Gene, int ProteinPosition)> knownHotspots)
    {
        ArgumentNullException.ThrowIfNull(knownHotspots);
        return knownHotspots.Contains((mutation.Gene, mutation.ProteinPosition));
    }

    /// <summary>
    /// Identifies driver mutations across a set of somatic coding mutations by applying the 20/20 rule
    /// per gene: a mutation is a driver if its gene classifies as an <see cref="DriverGeneRole.Oncogene"/>
    /// or <see cref="DriverGeneRole.TumorSuppressor"/>, OR if its (gene, position) is a known hotspot in
    /// <paramref name="knownHotspots"/>. The returned mutations are always a subset of the input, in input
    /// order (invariant: driver_mutations ⊆ somatic_mutations). Source: Vogelstein et al. (2013); Miller
    /// et al. (2017).
    /// </summary>
    /// <param name="mutations">Somatic coding mutations (one entry per observed mutation).</param>
    /// <param name="knownHotspots">Optional caller-supplied hotspot set; null is treated as empty.</param>
    /// <returns>The subset of input mutations that fall in a driver gene or at a known hotspot.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static IReadOnlyList<GeneMutation> IdentifyDriverMutations(
        IEnumerable<GeneMutation> mutations,
        IReadOnlySet<(string Gene, int ProteinPosition)>? knownHotspots = null)
    {
        ArgumentNullException.ThrowIfNull(mutations);

        var list = mutations as IReadOnlyList<GeneMutation> ?? mutations.ToList();
        var hotspots = knownHotspots ?? EmptyHotspots;

        // Classify each gene once from its full mutation spectrum.
        var driverGenes = new HashSet<string>();
        foreach (var byGene in list.GroupBy(m => m.Gene, StringComparer.Ordinal))
        {
            if (ClassifyGene(byGene).Role != DriverGeneRole.Ambiguous)
            {
                driverGenes.Add(byGene.Key);
            }
        }

        var drivers = new List<GeneMutation>();
        foreach (var mutation in list)
        {
            if (driverGenes.Contains(mutation.Gene) || hotspots.Contains((mutation.Gene, mutation.ProteinPosition)))
            {
                drivers.Add(mutation);
            }
        }

        return drivers;
    }

    private static readonly IReadOnlySet<(string Gene, int ProteinPosition)> EmptyHotspots =
        new HashSet<(string, int)>();

    private static bool IsTruncating(GeneMutation mutation) =>
        mutation.Consequence is MutationConsequence.Nonsense
            or MutationConsequence.Frameshift
            or MutationConsequence.SpliceSite;

    /// <summary>
    /// Counts mutations that are missense AND located at a recurrent position (a protein position carrying
    /// ≥ <see cref="RecurrentPositionMinCount"/> missense mutations), per Miller et al. (2017).
    /// </summary>
    private static int CountRecurrentMissense(IReadOnlyList<GeneMutation> mutations)
    {
        var missenseByPosition = new Dictionary<int, int>();
        foreach (var mutation in mutations)
        {
            if (mutation.Consequence == MutationConsequence.Missense)
            {
                missenseByPosition.TryGetValue(mutation.ProteinPosition, out int count);
                missenseByPosition[mutation.ProteinPosition] = count + 1;
            }
        }

        int recurrent = 0;
        foreach (var positionCount in missenseByPosition.Values)
        {
            if (positionCount >= RecurrentPositionMinCount)
            {
                recurrent += positionCount;
            }
        }

        return recurrent;
    }

    #endregion

    #region Sequencing Artifact Detection

    /// <summary>
    /// GIV (Global Imbalance Value) threshold above which a library is declared damaged for a given
    /// substitution type. Source: Chen L. et al. (2017), Science 355(6326):752–756, as summarized in
    /// Nature Methods (2017) 14:330 ("DNA variants or DNA damage?"): "A GIV score of 1 indicates there
    /// is no DNA damage and a GIV score above 1.5 is defined as damaged DNA." Value = 1.5.
    /// </summary>
    public const double DamagedGivThreshold = 1.5;

    /// <summary>
    /// GIV value of a perfectly balanced (undamaged) library: equal G&gt;T counts in read 1 and read 2.
    /// Source: Chen et al. (2017) / Nature Methods (2017) — GIV = 1 means no DNA damage. Value = 1.0.
    /// </summary>
    public const double UndamagedGivScore = 1.0;

    /// <summary>
    /// Minimum two-sided Fisher exact p-value used before Phred-scaling the strand-bias score, mirroring
    /// GATK's <c>FisherStrand.MIN_PVALUE</c>. Source: Broad Institute GATK, FisherStrand.java
    /// (<c>static final double MIN_PVALUE = 1E-320;</c>). Caps FS so a p-value of 0 does not produce
    /// an infinite Phred score.
    /// </summary>
    private const double MinFisherPValue = 1E-320;

    /// <summary>
    /// Classes of sequencing artifact distinguished by substitution type (and, for OxoG, read-orientation
    /// imbalance). The two artifact classes are disjoint by substitution: deamination is C:G&gt;T:A
    /// (C&gt;T / G&gt;A), oxidation is G:C&gt;T:A read as G&gt;T / C&gt;A
    /// (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    public enum ArtifactType
    {
        /// <summary>Not a recognized substitution-class artifact (a candidate true variant).</summary>
        None,

        /// <summary>
        /// FFPE cytosine-deamination artifact: C&gt;T or G&gt;A (collectively C:G&gt;T:A). Deaminated
        /// cytosine becomes uracil, which pairs with adenine. Source: Do &amp; Dobrovic (2015).
        /// </summary>
        FfpeDeamination,

        /// <summary>
        /// OxoG (8-oxoguanine) oxidative artifact: G&gt;T (read 1) or C&gt;A (read 2, reverse complement).
        /// Source: Chen et al. (2017).
        /// </summary>
        OxoG
    }

    /// <summary>
    /// One observed candidate variant together with the strand- and read-orientation read evidence needed
    /// for artifact classification. The strand counts feed the GATK FisherStrand 2×2 contingency table; the
    /// read-mate counts (<paramref name="AltReadsR1"/> / <paramref name="AltReadsR2"/>) feed the OxoG GIV
    /// imbalance. (The repository has no BAM reader; these counts are supplied directly rather than parsed
    /// from a BAM file — an API-shape decision that does not change the classification rules.)
    /// </summary>
    /// <param name="ReferenceAllele">Single reference base (A/C/G/T).</param>
    /// <param name="AlternateAllele">Single alternate base (A/C/G/T).</param>
    /// <param name="RefForward">Reference-supporting reads on the forward strand.</param>
    /// <param name="RefReverse">Reference-supporting reads on the reverse strand.</param>
    /// <param name="AltForward">Alternate-supporting reads on the forward strand.</param>
    /// <param name="AltReverse">Alternate-supporting reads on the reverse strand.</param>
    /// <param name="AltReadsR1">Alternate-supporting reads from read 1 of the pair (for GIV).</param>
    /// <param name="AltReadsR2">Alternate-supporting reads from read 2 of the pair (for GIV).</param>
    public readonly record struct ArtifactObservation(
        char ReferenceAllele,
        char AlternateAllele,
        int RefForward,
        int RefReverse,
        int AltForward,
        int AltReverse,
        int AltReadsR1,
        int AltReadsR2);

    /// <summary>
    /// Classification of a single candidate variant as a sequencing artifact (or not), with the supporting
    /// substitution class, GIV score and Phred-scaled strand-bias score.
    /// </summary>
    /// <param name="Type">Artifact class by substitution (and read orientation for OxoG).</param>
    /// <param name="GivScore">GIV imbalance for this variant's substitution (R1/R2); 1.0 means balanced.</param>
    /// <param name="StrandBiasPhred">Phred-scaled two-sided Fisher strand-bias score (FS); 0 means none.</param>
    /// <param name="IsArtifact">True when the variant is flagged as a likely artifact.</param>
    public readonly record struct ArtifactCall(
        ArtifactType Type,
        double GivScore,
        double StrandBiasPhred,
        bool IsArtifact);

    /// <summary>
    /// Computes the GIV (Global Imbalance Value) for one substitution type as the ratio of read-1 to read-2
    /// alternate-supporting read counts: GIV = r1Count / r2Count. For OxoG the canonical substitution is
    /// G&gt;T, which appears in excess in read 1 (its reverse complement C&gt;A appearing in read 2), so an
    /// elevated GIV signals oxidative damage. GIV = 1 means a balanced, undamaged library; GIV &gt; 1.5 is
    /// defined as damaged. When both counts are 0 there is no imbalance evidence and GIV = 1; when only the
    /// read-2 count is 0 the imbalance is maximal and GIV = <see cref="double.PositiveInfinity"/>.
    /// Source: Chen et al. (2017); Nature Methods (2017); Ettwiller Damage-estimator.
    /// </summary>
    /// <param name="r1Count">Alternate-supporting reads from read 1 (≥ 0).</param>
    /// <param name="r2Count">Alternate-supporting reads from read 2 (≥ 0).</param>
    /// <returns>The GIV ratio r1Count / r2Count (≥ 0); 1.0 when both are 0; +∞ when only r2Count is 0.</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative.</exception>
    public static double CalculateGivScore(int r1Count, int r2Count)
    {
        if (r1Count < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(r1Count), "Read-1 count cannot be negative.");
        }

        if (r2Count < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(r2Count), "Read-2 count cannot be negative.");
        }

        if (r2Count == 0)
        {
            // No read-2 support: balanced when read 1 is also empty (no imbalance evidence), otherwise a
            // maximal one-sided imbalance (Chen et al. 2017 — GIV is an R1/R2 ratio).
            return r1Count == 0 ? UndamagedGivScore : double.PositiveInfinity;
        }

        return (double)r1Count / r2Count;
    }

    /// <summary>
    /// Computes the GATK FisherStrand score FS: the Phred-scaled p-value of a two-sided Fisher exact test on
    /// the 2×2 strand contingency table [refForward, refReverse, altForward, altReverse], testing whether the
    /// reference and alternate alleles are distributed differently across forward/reverse strands (strand
    /// bias). FS = −10·log₁₀(max(p, MIN_PVALUE)); FS = 0 when there is no bias (p = 1) and grows as the
    /// alleles segregate by strand. Source: Broad Institute GATK FisherStrand / StrandBiasTest
    /// (table cell ordering ref-fwd, ref-rev, alt-fwd, alt-rev; FS = phredScaleErrorRate(p)).
    /// </summary>
    /// <param name="refForward">Reference reads on the forward strand (≥ 0).</param>
    /// <param name="refReverse">Reference reads on the reverse strand (≥ 0).</param>
    /// <param name="altForward">Alternate reads on the forward strand (≥ 0).</param>
    /// <param name="altReverse">Alternate reads on the reverse strand (≥ 0).</param>
    /// <returns>The Phred-scaled FisherStrand score FS (≥ 0).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative.</exception>
    public static double CalculateStrandBias(int refForward, int refReverse, int altForward, int altReverse)
    {
        RequireNonNegative(refForward, nameof(refForward));
        RequireNonNegative(refReverse, nameof(refReverse));
        RequireNonNegative(altForward, nameof(altForward));
        RequireNonNegative(altReverse, nameof(altReverse));

        double pValue = FisherExactTwoSided(refForward, refReverse, altForward, altReverse);
        double floored = Math.Max(pValue, MinFisherPValue);

        // Phred-scaled error rate: FS = -10 * log10(p) (GATK QualityUtils.phredScaleErrorRate).
        return -10.0 * Math.Log10(floored);
    }

    /// <summary>
    /// Classifies one candidate variant as an artifact by substitution class. C&gt;T / G&gt;A are FFPE
    /// cytosine-deamination artifacts (Do &amp; Dobrovic 2015); G&gt;T / C&gt;A are OxoG oxidative artifacts
    /// (Chen et al. 2017). Any other substitution is <see cref="ArtifactType.None"/>. The returned call also
    /// carries the GIV score (from the read-mate alt counts) and the Phred-scaled strand-bias FS, and flags
    /// the variant as an artifact when it is a deamination class, OR an OxoG class whose GIV exceeds the
    /// damaged threshold (1.5).
    /// </summary>
    /// <param name="observation">The candidate variant with its strand / read-mate read evidence.</param>
    /// <returns>The artifact classification for this variant.</returns>
    public static ArtifactCall ClassifyArtifact(ArtifactObservation observation)
    {
        char reference = char.ToUpperInvariant(observation.ReferenceAllele);
        char alternate = char.ToUpperInvariant(observation.AlternateAllele);

        ArtifactType type = ClassifySubstitution(reference, alternate);
        double giv = CalculateGivScore(observation.AltReadsR1, observation.AltReadsR2);
        double strandBias = CalculateStrandBias(
            observation.RefForward, observation.RefReverse, observation.AltForward, observation.AltReverse);

        // Deamination (C>T/G>A) is flagged by substitution class alone; OxoG (G>T/C>A) is confirmed by the
        // read-orientation imbalance (GIV > 1.5 = damaged, Chen et al. 2017).
        bool isArtifact = type switch
        {
            ArtifactType.FfpeDeamination => true,
            ArtifactType.OxoG => giv > DamagedGivThreshold,
            _ => false
        };

        return new ArtifactCall(type, giv, strandBias, isArtifact);
    }

    /// <summary>
    /// Detects OxoG (8-oxoguanine) artifacts among candidate variants: returns the calls for variants whose
    /// substitution is the OxoG class (G&gt;T / C&gt;A) AND whose GIV read-orientation imbalance exceeds the
    /// damaged threshold (GIV &gt; 1.5). Source: Chen et al. (2017).
    /// </summary>
    /// <param name="variants">Candidate variant observations.</param>
    /// <returns>The OxoG-artifact calls, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ArtifactCall> DetectOxoGArtifacts(IEnumerable<ArtifactObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var oxoG = new List<ArtifactCall>();
        foreach (var variant in variants)
        {
            ArtifactCall call = ClassifyArtifact(variant);
            if (call.Type == ArtifactType.OxoG && call.IsArtifact)
            {
                oxoG.Add(call);
            }
        }

        return oxoG;
    }

    /// <summary>
    /// Filters sequencing artifacts out of a candidate variant set, returning only the variants that are NOT
    /// flagged as artifacts (per <see cref="ClassifyArtifact"/>). The result is always a subset of the input,
    /// in input order. Source: composition of the FFPE-deamination and OxoG substitution-class rules
    /// (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    /// <param name="variants">Candidate variant observations.</param>
    /// <returns>The subset of <paramref name="variants"/> not classified as artifacts, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ArtifactObservation> FilterArtifacts(IEnumerable<ArtifactObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var kept = new List<ArtifactObservation>();
        foreach (var variant in variants)
        {
            if (!ClassifyArtifact(variant).IsArtifact)
            {
                kept.Add(variant);
            }
        }

        return kept;
    }

    /// <summary>
    /// Maps a single-base substitution to its artifact class. C&gt;T / G&gt;A = FFPE deamination;
    /// G&gt;T / C&gt;A = OxoG oxidation; everything else = none (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    private static ArtifactType ClassifySubstitution(char reference, char alternate)
    {
        return (reference, alternate) switch
        {
            ('C', 'T') => ArtifactType.FfpeDeamination,
            ('G', 'A') => ArtifactType.FfpeDeamination,
            ('G', 'T') => ArtifactType.OxoG,
            ('C', 'A') => ArtifactType.OxoG,
            _ => ArtifactType.None
        };
    }

    private static void RequireNonNegative(int value, string paramName)
    {
        if (value < 0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, "Read count cannot be negative.");
        }
    }

    /// <summary>
    /// Two-sided Fisher exact test p-value for the 2×2 table
    /// <code>[[a, b], [c, d]]</code> (a = refForward, b = refReverse, c = altForward, d = altReverse).
    /// Sums the hypergeometric probabilities of all tables with the same margins whose probability is
    /// ≤ that of the observed table (the conventional two-sided definition). Source: GATK FisherStrand uses
    /// <c>FisherExactTest.twoSidedPValue</c> on this 2×2 strand table.
    /// </summary>
    private static double FisherExactTwoSided(int a, int b, int c, int d)
    {
        int rowOne = a + b;
        int rowTwo = c + d;
        int colOne = a + c;
        int total = a + b + c + d;

        if (total == 0)
        {
            // An empty table provides no evidence of strand bias: p = 1.
            return 1.0;
        }

        double observedLogProb = HypergeometricLogProbability(a, rowOne, rowTwo, colOne, total);

        // The cell 'a' ranges over the values compatible with the fixed margins.
        int minA = Math.Max(0, colOne - rowTwo);
        int maxA = Math.Min(rowOne, colOne);

        double pValue = 0.0;
        for (int x = minA; x <= maxA; x++)
        {
            double logProb = HypergeometricLogProbability(x, rowOne, rowTwo, colOne, total);
            // Include tables at least as extreme (≤ observed probability), with a small tolerance for
            // floating-point comparison of equal-probability tables.
            if (logProb <= observedLogProb + 1e-7)
            {
                pValue += Math.Exp(logProb);
            }
        }

        return Math.Min(1.0, pValue);
    }

    /// <summary>
    /// log P(table) under the hypergeometric distribution for a 2×2 table with the given margins, where the
    /// top-left cell equals <paramref name="a"/>:
    /// P = C(rowOne, a)·C(rowTwo, colOne−a) / C(total, colOne).
    /// </summary>
    private static double HypergeometricLogProbability(int a, int rowOne, int rowTwo, int colOne, int total)
    {
        return LogChoose(rowOne, a)
             + LogChoose(rowTwo, colOne - a)
             - LogChoose(total, colOne);
    }

    /// <summary>log of the binomial coefficient C(n, k) via log-gamma (numerically stable for read counts).</summary>
    private static double LogChoose(int n, int k)
    {
        if (k < 0 || k > n)
        {
            return double.NegativeInfinity;
        }

        return LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);
    }

    /// <summary>log(n!) = logΓ(n+1).</summary>
    private static double LogFactorial(int n) => LogGamma(n + 1.0);

    /// <summary>
    /// Lanczos approximation of the natural log of the gamma function. Coefficients g = 7, n = 9 per the
    /// standard Lanczos series (Numerical Recipes / Lanczos 1964); accurate to ~1e-13 for the positive
    /// arguments used here, so the resulting binomial coefficients are exact to floating precision for the
    /// read-count magnitudes encountered in strand-bias tables.
    /// </summary>
    private static double LogGamma(double x)
    {
        // Lanczos coefficients (g = 7).
        double[] coefficients =
        {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };

        const double g = 7.0;
        x -= 1.0;
        double sum = coefficients[0];
        for (int i = 1; i < coefficients.Length; i++)
        {
            sum += coefficients[i] / (x + i);
        }

        double t = x + g + 0.5;
        return 0.5 * Math.Log(2.0 * Math.PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(sum);
    }

    #endregion

    #region Cancer Variant Annotation (AMP/ASCO/CAP 2017 tiers)

    /// <summary>
    /// Minor-allele-frequency (MAF) cutoff at or above which a variant is treated as a common
    /// polymorphism and classified Tier IV (benign / likely benign). Source: Li MM et al. (2017),
    /// J Mol Diagn 19(1):4–23 — "the work group recommends using 1% (0.01) as a primary cutoff" for
    /// eliminating polymorphic or benign variants (Population Databases section), and Table 7 (Tier IV)
    /// lists "MAF ≥ 1% in the general population" as the population-database criterion. Value = 0.01.
    /// </summary>
    public const double BenignPopulationMafThreshold = 0.01;

    /// <summary>
    /// Strength of clinical/experimental evidence supporting a variant as a biomarker, per the
    /// four evidence levels of Li MM et al. (2017), Table 3 / Figure 2. Levels A and B map to Tier I
    /// (strong clinical significance); Levels C and D map to Tier II (potential clinical significance).
    /// </summary>
    public enum ClinicalEvidenceLevel
    {
        /// <summary>No biomarker evidence level assigned (the variant is not a known biomarker).</summary>
        None,

        /// <summary>
        /// Level A: biomarkers that predict response/resistance to FDA-approved therapies for a specific
        /// tumor type, or are included in professional guidelines (therapeutic/diagnostic/prognostic).
        /// Maps to Tier I. Source: Li et al. (2017), Table 3.
        /// </summary>
        A,

        /// <summary>
        /// Level B: biomarkers based on well-powered studies with expert consensus. Maps to Tier I.
        /// Source: Li et al. (2017), Table 3.
        /// </summary>
        B,

        /// <summary>
        /// Level C: FDA-approved/guideline therapies for a different tumor type (off-label), clinical-trial
        /// inclusion criteria, or diagnostic/prognostic significance from multiple small studies. Maps to
        /// Tier II. Source: Li et al. (2017), Table 3.
        /// </summary>
        C,

        /// <summary>
        /// Level D: plausible therapeutic significance from preclinical studies, or diagnostic/prognostic
        /// support from small studies / case reports without consensus. Maps to Tier II.
        /// Source: Li et al. (2017), Table 3.
        /// </summary>
        D
    }

    /// <summary>
    /// AMP/ASCO/CAP 2017 four-tier clinical-significance classification of a somatic sequence variant.
    /// Source: Li MM et al. (2017), J Mol Diagn 19(1):4–23, Figure 2.
    /// </summary>
    public enum VariantTier
    {
        /// <summary>Tier I: variants of strong clinical significance (Level A or B evidence).</summary>
        TierI_StrongClinicalSignificance,

        /// <summary>Tier II: variants of potential clinical significance (Level C or D evidence).</summary>
        TierII_PotentialClinicalSignificance,

        /// <summary>
        /// Tier III: variants of unknown clinical significance — not common in population databases and
        /// with no convincing published evidence of cancer association.
        /// </summary>
        TierIII_UnknownClinicalSignificance,

        /// <summary>
        /// Tier IV: benign or likely benign variants — observed at a significant allele frequency
        /// (MAF ≥ 1%) in population databases, or with no evidence of cancer association.
        /// </summary>
        TierIV_BenignOrLikelyBenign
    }

    /// <summary>
    /// Caller-supplied evidence for one somatic variant, reduced to the features the AMP/ASCO/CAP 2017
    /// tiering rule consumes (Li et al. 2017, Figure 2 / Tables 4–7). The guideline classifies variants
    /// from external knowledge (professional guidelines, population databases, somatic databases,
    /// literature); this library does not reproduce those curated resources — the relevant facts are
    /// supplied by the caller, who has performed the database lookups.
    /// </summary>
    /// <param name="Gene">Gene symbol the variant falls in.</param>
    /// <param name="ProteinChange">HGVS protein change (e.g. p.V600E); informational.</param>
    /// <param name="EvidenceLevel">Strongest assigned clinical evidence level (A–D, or None).</param>
    /// <param name="PopulationMaf">
    /// Minor allele frequency in a population database (e.g. gnomAD/ExAC/1000 Genomes), in [0, 1].
    /// A value ≥ <see cref="BenignPopulationMafThreshold"/> indicates a common polymorphism.
    /// </param>
    /// <param name="HasCancerAssociation">
    /// True when there is published evidence associating the variant with cancer (somatic database
    /// presence, functional/population study). Distinguishes Tier III from Tier IV when MAF is low.
    /// </param>
    public readonly record struct CancerVariantAnnotationInput(
        string Gene,
        string ProteinChange,
        ClinicalEvidenceLevel EvidenceLevel,
        double PopulationMaf,
        bool HasCancerAssociation);

    /// <summary>The tier classification of one variant, with the input evidence that produced it.</summary>
    /// <param name="Variant">The variant evidence that was classified.</param>
    /// <param name="Tier">Assigned AMP/ASCO/CAP 2017 tier.</param>
    public readonly record struct CancerVariantAnnotation(
        CancerVariantAnnotationInput Variant,
        VariantTier Tier);

    /// <summary>
    /// Classifies a single somatic variant into the AMP/ASCO/CAP 2017 four-tier system from caller-supplied
    /// evidence, applying the decision criteria of Li MM et al. (2017), Figure 2 in priority order:
    /// <list type="number">
    /// <item><description>Level A or B evidence ⇒ <see cref="VariantTier.TierI_StrongClinicalSignificance"/>.</description></item>
    /// <item><description>Level C or D evidence ⇒ <see cref="VariantTier.TierII_PotentialClinicalSignificance"/>.</description></item>
    /// <item><description>Otherwise, MAF ≥ 1% (common polymorphism) OR no cancer association ⇒
    /// <see cref="VariantTier.TierIV_BenignOrLikelyBenign"/> (Table 7).</description></item>
    /// <item><description>Otherwise (rare, no clinical evidence, but a cancer association exists) ⇒
    /// <see cref="VariantTier.TierIII_UnknownClinicalSignificance"/> (Table 6).</description></item>
    /// </list>
    /// Clinical evidence (Tier I/II) is evaluated before the benign-frequency rule because a Level A/B
    /// biomarker remains strongly significant even if it also appears in population databases; Table 4
    /// (Tier I) and Table 5 (Tier II) note such variants are "absent or extremely low MAF" but the
    /// guideline assigns them by evidence level, not frequency.
    /// </summary>
    /// <param name="variant">Caller-supplied evidence for the variant.</param>
    /// <returns>The variant's tier classification.</returns>
    /// <exception cref="ArgumentOutOfRangeException">PopulationMaf is NaN or outside [0, 1].</exception>
    public static VariantTier ClassifyVariantTier(CancerVariantAnnotationInput variant)
    {
        if (double.IsNaN(variant.PopulationMaf) || variant.PopulationMaf < 0.0 || variant.PopulationMaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.PopulationMaf, "Population MAF must be in the range [0, 1].");
        }

        // Tier I — strong clinical significance: Level A or B evidence (Li et al. 2017, Figure 2).
        if (variant.EvidenceLevel is ClinicalEvidenceLevel.A or ClinicalEvidenceLevel.B)
        {
            return VariantTier.TierI_StrongClinicalSignificance;
        }

        // Tier II — potential clinical significance: Level C or D evidence (Li et al. 2017, Figure 2).
        if (variant.EvidenceLevel is ClinicalEvidenceLevel.C or ClinicalEvidenceLevel.D)
        {
            return VariantTier.TierII_PotentialClinicalSignificance;
        }

        // No clinical evidence level. Tier IV — benign/likely benign: observed at a significant allele
        // frequency (MAF ≥ 1%, Table 7), OR no published evidence of cancer association (Figure 2).
        if (variant.PopulationMaf >= BenignPopulationMafThreshold || !variant.HasCancerAssociation)
        {
            return VariantTier.TierIV_BenignOrLikelyBenign;
        }

        // Tier III — unknown clinical significance: rare (low MAF), no clinical evidence, but a cancer
        // association exists so it cannot be called benign (Li et al. 2017, Table 6).
        return VariantTier.TierIII_UnknownClinicalSignificance;
    }

    /// <summary>
    /// Annotates a set of somatic variants with their AMP/ASCO/CAP 2017 clinical-significance tiers by
    /// applying <see cref="ClassifyVariantTier"/> to each variant. The output preserves input order and
    /// has one entry per input variant. Source: Li MM et al. (2017), J Mol Diagn 19(1):4–23.
    /// </summary>
    /// <param name="variants">Caller-supplied variant evidence records.</param>
    /// <returns>One <see cref="CancerVariantAnnotation"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A variant's PopulationMaf is outside [0, 1].</exception>
    public static IReadOnlyList<CancerVariantAnnotation> AnnotateCancerVariants(
        IEnumerable<CancerVariantAnnotationInput> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var annotations = new List<CancerVariantAnnotation>();
        foreach (var variant in variants)
        {
            annotations.Add(new CancerVariantAnnotation(variant, ClassifyVariantTier(variant)));
        }

        return annotations;
    }

    /// <summary>
    /// Looks up a variant's COSMIC (Catalogue Of Somatic Mutations In Cancer) annotation in a
    /// caller-supplied catalog keyed by (gene, protein change). COSMIC is a large, expert-curated
    /// somatic-mutation database (Tate JG et al. 2019, Nucleic Acids Res 47:D941–D947) that cannot be
    /// reproduced or hardcoded here; the caller passes the relevant records (e.g. a COSMIC export),
    /// and this method performs the exact-match lookup the AMP/ASCO/CAP workflow uses to flag a variant
    /// as present in a somatic database (Li et al. 2017, Tables 4–6, "Somatic database: COSMIC...").
    /// </summary>
    /// <param name="variant">The variant to look up.</param>
    /// <param name="cosmicCatalog">
    /// Caller-supplied COSMIC records keyed by (gene, protein change); e.g. COSMIC identifier strings.
    /// </param>
    /// <returns>
    /// The catalog value (e.g. a COSMIC ID) for the variant's (gene, protein change), or <c>null</c>
    /// when the variant is not present in the supplied catalog.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="cosmicCatalog"/> is null.</exception>
    public static string? GetCOSMICAnnotation(
        CancerVariantAnnotationInput variant,
        IReadOnlyDictionary<(string Gene, string ProteinChange), string> cosmicCatalog)
    {
        ArgumentNullException.ThrowIfNull(cosmicCatalog);
        return cosmicCatalog.TryGetValue((variant.Gene, variant.ProteinChange), out string? id) ? id : null;
    }

    #endregion

    #region Tumor Mutational Burden (TMB)

    /// <summary>
    /// FDA-approved tumor-mutational-burden-high (TMB-H) cutoff for pembrolizumab, in mutations per
    /// megabase. A tumor is TMB-High when its TMB is at or above this value (the threshold is inclusive).
    /// Source: Marcus L et al. (2021), FDA Approval Summary: Pembrolizumab for the Treatment of Tumor
    /// Mutational Burden–High Solid Tumors, Clin Cancer Res 27(17):4685–4689 — pembrolizumab approved
    /// (June 16, 2020) for solid tumors with "TMB ≥10 mutations/megabase (mut/Mb)", companion diagnostic
    /// FoundationOne CDx. Value = 10.0 mut/Mb.
    /// </summary>
    public const double TmbHighThreshold = 10.0;

    /// <summary>
    /// TMB-High classification result. Only the FDA/F1CDx TMB-High cutoff (≥ 10 mut/Mb) is source-backed,
    /// so the classification is two-tier (High vs Low); intermediate research cut-points (e.g. tumor-type
    /// specific 6/20 boundaries) are not implemented because no authoritative source defines them.
    /// Source: Marcus et al. (2021).
    /// </summary>
    public enum TmbStatus
    {
        /// <summary>TMB below the FDA TMB-High cutoff (TMB &lt; 10 mut/Mb).</summary>
        Low,

        /// <summary>TMB at or above the FDA TMB-High cutoff (TMB ≥ 10 mut/Mb).</summary>
        High
    }

    /// <summary>
    /// Computes the tumor mutational burden (TMB) as the number of somatic coding mutations per megabase of
    /// sequenced region: TMB = <paramref name="mutationCount"/> / <paramref name="targetRegionMb"/>, in
    /// mutations/megabase (mut/Mb). Source: Chalmers ZR et al. (2017), Genome Medicine 9:34 — "TMB was
    /// defined as the number of somatic, coding, base substitution, and indel mutations per megabase of
    /// genome examined" (the FoundationOne 315-gene panel denominator is 1.1 Mb of coding genome). The
    /// caller supplies the already-filtered somatic mutation count (germline / known-driver filtering is the
    /// upstream somatic caller's responsibility, per Chalmers 2017's pre-count filtering and ONCO-SOMATIC-001).
    /// </summary>
    /// <param name="mutationCount">Number of counted somatic mutations (≥ 0).</param>
    /// <param name="targetRegionMb">Size of the sequenced coding region in megabases (&gt; 0, finite).</param>
    /// <returns>TMB in mutations per megabase (≥ 0).</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="mutationCount"/> is negative, or <paramref name="targetRegionMb"/> is not a finite
    /// value greater than 0 (TMB is undefined when the megabase denominator is 0).
    /// </exception>
    public static double CalculateTMB(int mutationCount, double targetRegionMb)
    {
        if (mutationCount < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutationCount), mutationCount, "Mutation count cannot be negative.");
        }

        if (double.IsNaN(targetRegionMb) || double.IsInfinity(targetRegionMb) || targetRegionMb <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(targetRegionMb), targetRegionMb,
                "Target region size (Mb) must be a finite value greater than 0; TMB is undefined at 0 Mb.");
        }

        // TMB = mutations / megabase (Chalmers et al. 2017).
        return mutationCount / targetRegionMb;
    }

    /// <summary>
    /// Computes TMB from a set of classified somatic calls (e.g. from <see cref="CallSomaticMutations"/>):
    /// counts the calls with <see cref="SomaticStatus.Somatic"/> status and divides by
    /// <paramref name="targetRegionMb"/>. Germline and not-detected calls are excluded, matching the TMB
    /// definition of counting only somatic mutations (Chalmers et al. 2017; Friends of Cancer Research TMB
    /// Harmonization, Merino et al. 2020). Thin wrapper over
    /// <see cref="CalculateTMB(int, double)"/>.
    /// </summary>
    /// <param name="calls">Classified somatic calls; only <see cref="SomaticStatus.Somatic"/> are counted.</param>
    /// <param name="targetRegionMb">Size of the sequenced coding region in megabases (&gt; 0, finite).</param>
    /// <returns>TMB in mutations per megabase (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="calls"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="targetRegionMb"/> is not finite and &gt; 0.</exception>
    public static double CalculateTMB(IEnumerable<SomaticCall> calls, double targetRegionMb)
    {
        ArgumentNullException.ThrowIfNull(calls);

        int somaticCount = calls.Count(c => c.Status == SomaticStatus.Somatic);
        return CalculateTMB(somaticCount, targetRegionMb);
    }

    /// <summary>
    /// Classifies a TMB value (mut/Mb) as <see cref="TmbStatus.High"/> when it is at or above the FDA
    /// TMB-High cutoff (≥ <see cref="TmbHighThreshold"/> = 10 mut/Mb; the boundary is inclusive), otherwise
    /// <see cref="TmbStatus.Low"/>. Source: Marcus et al. (2021), FDA Approval Summary — pembrolizumab for
    /// "TMB ≥10 mut/Mb" solid tumors (companion diagnostic FoundationOne CDx).
    /// </summary>
    /// <param name="tmb">Tumor mutational burden in mutations per megabase (≥ 0, finite).</param>
    /// <returns><see cref="TmbStatus.High"/> if tmb ≥ 10, else <see cref="TmbStatus.Low"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tmb"/> is negative or not finite.</exception>
    public static TmbStatus ClassifyTMB(double tmb)
    {
        if (double.IsNaN(tmb) || double.IsInfinity(tmb) || tmb < 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(tmb), tmb, "TMB must be a finite value ≥ 0.");
        }

        return tmb >= TmbHighThreshold ? TmbStatus.High : TmbStatus.Low;
    }

    #endregion

    #region Microsatellite Instability (MSI)

    /// <summary>
    /// MSI score cutoff (as a fraction in [0,1]) at or above which a sample is microsatellite-instability-high
    /// (MSI-H). Source: niu-lab/msisensor2 README — "The recommended msi score cutoff value is 20%
    /// (msi high: msi score &gt;= 20%)", where the msi score is "number of msi sites / all valid sites". The
    /// boundary is inclusive. Value = 0.20 (20%).
    /// </summary>
    public const double MsiHighScoreThreshold = 0.20;

    /// <summary>
    /// Minimum number of unstable markers (out of the validated 5-marker reference panel) for a tumor to be
    /// classified MSI-H under the NCI/Bethesda categorical criteria. Source: Boland CR et al. (1998), Cancer
    /// Res 58(22):5248–5257 — a tumor is MSI-H "if two or more of the five markers show instability". Value = 2.
    /// </summary>
    public const int BethesdaMsiHighMarkerCount = 2;

    /// <summary>
    /// Number of unstable markers for MSI-L under the NCI/Bethesda criteria — "only one of the five markers
    /// shows instability". Source: Boland et al. (1998). Value = 1.
    /// </summary>
    public const int BethesdaMsiLowMarkerCount = 1;

    /// <summary>
    /// Microsatellite-instability status of a sample.
    /// <list type="bullet">
    /// <item><description><b>MSS</b> — microsatellite stable.</description></item>
    /// <item><description><b>MSI_Low</b> — low-frequency MSI (Bethesda: exactly 1 of 5 markers unstable).</description></item>
    /// <item><description><b>MSI_High</b> — high-frequency MSI (Bethesda ≥2/5; MSIsensor2 score ≥20%).</description></item>
    /// </list>
    /// Sources: Boland et al. (1998) (categorical MSS/MSI-L/MSI-H); niu-lab/msisensor2 (continuous score ≥20% → MSI-H).
    /// </summary>
    public enum MsiStatus
    {
        /// <summary>Microsatellite stable (no instability above the calling threshold).</summary>
        MSS,

        /// <summary>Low-frequency microsatellite instability (Bethesda: exactly one marker unstable).</summary>
        MSI_Low,

        /// <summary>High-frequency microsatellite instability (Bethesda ≥2/5 markers; MSIsensor2 score ≥20%).</summary>
        MSI_High
    }

    /// <summary>Result of an end-to-end MSI determination over a set of evaluated microsatellite loci.</summary>
    /// <param name="UnstableLoci">Number of loci called unstable (somatic indel in tumor vs normal).</param>
    /// <param name="TotalLoci">Total number of valid evaluated microsatellite loci.</param>
    /// <param name="Score">MSI score = <paramref name="UnstableLoci"/> / <paramref name="TotalLoci"/>, in [0,1].</param>
    /// <param name="Status">MSI-H (score ≥ 20%) or MSS, per the MSIsensor2 continuous-score cutoff.</param>
    public readonly record struct MsiResult(int UnstableLoci, int TotalLoci, double Score, MsiStatus Status);

    /// <summary>
    /// Computes the MSI score as the fraction of unstable microsatellite loci among all valid evaluated loci:
    /// score = <paramref name="unstableLoci"/> / <paramref name="totalLoci"/>. Source: niu-lab/msisensor2 —
    /// "the msi score (number of msi sites / all valid sites)"; Niu et al. (2014), Bioinformatics 30(7):1015 —
    /// the MSI score is the percentage of microsatellite sites with a somatic indel. Returned as a fraction in
    /// [0,1] (multiply by 100 for a percentage).
    /// </summary>
    /// <param name="unstableLoci">Number of loci called unstable (0 ≤ unstableLoci ≤ totalLoci).</param>
    /// <param name="totalLoci">Total number of valid evaluated loci (&gt; 0).</param>
    /// <returns>MSI score in [0,1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="totalLoci"/> ≤ 0 (score undefined with no valid loci), <paramref name="unstableLoci"/>
    /// is negative, or <paramref name="unstableLoci"/> &gt; <paramref name="totalLoci"/>.
    /// </exception>
    public static double CalculateMSIScore(int unstableLoci, int totalLoci)
    {
        if (totalLoci <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalLoci), totalLoci,
                "Total valid loci must be > 0; the MSI score is undefined with no evaluable loci.");
        }

        if (unstableLoci < 0 || unstableLoci > totalLoci)
        {
            throw new ArgumentOutOfRangeException(
                nameof(unstableLoci), unstableLoci,
                "Unstable loci must satisfy 0 ≤ unstableLoci ≤ totalLoci.");
        }

        // MSI score = unstable loci / valid loci (MSIsensor2; Niu et al. 2014).
        return (double)unstableLoci / totalLoci;
    }

    /// <summary>
    /// Classifies a continuous MSI score (fraction in [0,1]) as <see cref="MsiStatus.MSI_High"/> when it is at
    /// or above the MSIsensor2 cutoff (≥ <see cref="MsiHighScoreThreshold"/> = 20%; boundary inclusive),
    /// otherwise <see cref="MsiStatus.MSS"/>. Source: niu-lab/msisensor2 README — "msi high: msi score &gt;= 20%".
    /// MSIsensor2 defines only a binary MSI-H cutoff on the continuous score, so no MSI-L band is applied here;
    /// MSI-L is a marker-count concept handled by <see cref="ClassifyBethesdaPanel"/>.
    /// </summary>
    /// <param name="score">MSI score as a fraction in [0,1].</param>
    /// <returns><see cref="MsiStatus.MSI_High"/> if score ≥ 0.20, else <see cref="MsiStatus.MSS"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is not a finite value in [0,1].</exception>
    public static MsiStatus ClassifyMSIStatus(double score)
    {
        if (double.IsNaN(score) || double.IsInfinity(score) || score < 0.0 || score > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(score), score, "MSI score must be a finite value in [0,1].");
        }

        return score >= MsiHighScoreThreshold ? MsiStatus.MSI_High : MsiStatus.MSS;
    }

    /// <summary>
    /// Classifies a sample under the NCI/Bethesda categorical criteria from the number of unstable markers in a
    /// fixed reference panel: <see cref="MsiStatus.MSI_High"/> when ≥ 2 markers are unstable,
    /// <see cref="MsiStatus.MSI_Low"/> when exactly 1 is unstable, and <see cref="MsiStatus.MSS"/> when none is.
    /// Source: Boland CR et al. (1998), Cancer Res 58(22):5248–5257 — MSI-H "if two or more of the five markers
    /// show instability", MSI-L "if only one of the five markers shows instability", MSS otherwise.
    /// </summary>
    /// <param name="unstableMarkers">Number of unstable markers (0 ≤ unstableMarkers ≤ totalMarkers).</param>
    /// <param name="totalMarkers">Total markers evaluated in the panel (&gt; 0; classically 5).</param>
    /// <returns>The Bethesda MSI status.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="totalMarkers"/> ≤ 0, <paramref name="unstableMarkers"/> negative, or
    /// <paramref name="unstableMarkers"/> &gt; <paramref name="totalMarkers"/>.
    /// </exception>
    public static MsiStatus ClassifyBethesdaPanel(int unstableMarkers, int totalMarkers)
    {
        if (totalMarkers <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalMarkers), totalMarkers, "Total markers must be > 0.");
        }

        if (unstableMarkers < 0 || unstableMarkers > totalMarkers)
        {
            throw new ArgumentOutOfRangeException(
                nameof(unstableMarkers), unstableMarkers,
                "Unstable markers must satisfy 0 ≤ unstableMarkers ≤ totalMarkers.");
        }

        // Boland et al. (1998): ≥2 unstable → MSI-H; exactly 1 → MSI-L; 0 → MSS.
        if (unstableMarkers >= BethesdaMsiHighMarkerCount)
        {
            return MsiStatus.MSI_High;
        }

        return unstableMarkers == BethesdaMsiLowMarkerCount ? MsiStatus.MSI_Low : MsiStatus.MSS;
    }

    /// <summary>
    /// End-to-end MSI determination from per-locus stability calls: counts the unstable loci, computes the MSI
    /// score (<see cref="CalculateMSIScore"/>), and classifies it with the MSIsensor2 continuous-score cutoff
    /// (<see cref="ClassifyMSIStatus"/>). Each element of <paramref name="locusUnstableFlags"/> is one valid
    /// evaluated microsatellite locus: <c>true</c> = unstable (somatic indel), <c>false</c> = stable. This
    /// matches the MSIsensor pipeline where per-locus instability (chi-square tumor-vs-normal length
    /// distributions) is determined upstream and the sample score is the fraction of unstable loci
    /// (Niu et al. 2014; niu-lab/msisensor2).
    /// </summary>
    /// <param name="locusUnstableFlags">Per-locus stability flags (one entry per valid evaluated locus).</param>
    /// <returns>The unstable/total counts, the MSI score, and the MSI status.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="locusUnstableFlags"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The sequence is empty (no valid loci to evaluate).</exception>
    public static MsiResult DetectMSI(IEnumerable<bool> locusUnstableFlags)
    {
        ArgumentNullException.ThrowIfNull(locusUnstableFlags);

        int total = 0;
        int unstable = 0;
        foreach (bool isUnstable in locusUnstableFlags)
        {
            total++;
            if (isUnstable)
            {
                unstable++;
            }
        }

        if (total == 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(locusUnstableFlags), "At least one valid locus is required; the MSI score is undefined for an empty set.");
        }

        double score = CalculateMSIScore(unstable, total);
        return new MsiResult(unstable, total, score, ClassifyMSIStatus(score));
    }

    #endregion

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
    /// TAI and LST remain caller-supplied here because their faithful derivation (scarHRD
    /// <c>calc.ai_new</c> / <c>calc.lst</c>) requires the exact per-build centromere/telomere
    /// <c>chrominfo</c> coordinate table that scarHRD ships as binary package data and that could not be
    /// retrieved as a citable reference table; their telomeric-vs-interstitial classification (TAI) and
    /// p/q-arm split (LST) are sensitive to those coordinates, so deriving them from an unverified table
    /// would not reproduce scarHRD. Only the LOH component, whose definition needs no centromere table,
    /// is derived here. See <c>docs/Validation/LIMITATIONS.md</c> §2 (ONCO-HRD-001).
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

    #region SBS-96 Trinucleotide Context Catalog

    /// <summary>
    /// Number of single-base-substitution channels in the SBS-96 classification: six pyrimidine substitution
    /// subtypes × four 5' bases × four 3' bases. Source: Alexandrov et al. (2013), Nature 500:415-421 — "96
    /// possible mutation types (6 types of substitution * 4 types of 5' base * 4 types of 3' base)"; COSMIC
    /// SBS96 (https://cancer.sanger.ac.uk/signatures/sbs/sbs96/). Value = 96.
    /// </summary>
    public const int Sbs96ChannelCount = 96;

    /// <summary>
    /// The four 5'/3' flanking bases of an SBS-96 trinucleotide context, in canonical alphabetical order
    /// (A, C, G, T). Source: COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019, BMC Genomics
    /// 20:685) — each substitution has "sixteen possible trinucleotide (4 types of 5' base * 4 types of 3' base)".
    /// </summary>
    private static readonly char[] ContextBases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// The six single-base substitutions of the SBS-6 / SBS-96 classification, each referred to by the
    /// pyrimidine (C or T) of the mutated Watson-Crick base pair, in canonical order. Source: COSMIC SBS96 —
    /// "C>A, C>G, C>T, T>A, T>C, and T>G"; Alexandrov et al. (2013); Bergstrom et al. (2019).
    /// </summary>
    private static readonly (char Ref, char Alt)[] PyrimidineSubstitutions =
    {
        ('C', 'A'), ('C', 'G'), ('C', 'T'), ('T', 'A'), ('T', 'C'), ('T', 'G')
    };

    /// <summary>
    /// Classifies a single-base substitution into its SBS-96 trinucleotide-context channel, folded onto the
    /// pyrimidine strand. The mutated base sits in the centre of the trinucleotide; the channel is rendered as
    /// <c>5'[REF&gt;ALT]3'</c> (e.g. <c>A[C&gt;A]A</c>). Each substitution is referred to by the pyrimidine of
    /// the mutated Watson-Crick base pair: when the reference (mutated) base is a purine (A or G), the
    /// trinucleotide context and the substitution are reverse-complemented onto the pyrimidine strand before
    /// counting (e.g. a G&gt;T at 5'-T G A-3' folds to <c>T[C&gt;A]A</c>). Source: Alexandrov et al. (2013),
    /// Nature 500:415-421; COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019) — "using the purine
    /// base of the Watson-Crick base-pair for classifying mutation types will require taking the reverse
    /// complement sequence". Complement map A↔T, C↔G (Watson-Crick base pairing).
    /// </summary>
    /// <param name="fivePrime">Reference base immediately 5' of the mutated base (A/C/G/T, case-insensitive).</param>
    /// <param name="referenceBase">The mutated (reference) base (A/C/G/T, case-insensitive).</param>
    /// <param name="alternateBase">The substituted (alternate) base (A/C/G/T, case-insensitive, ≠ reference).</param>
    /// <param name="threePrime">Reference base immediately 3' of the mutated base (A/C/G/T, case-insensitive).</param>
    /// <returns>The SBS-96 channel label in the form <c>5'[REF&gt;ALT]3'</c> with a pyrimidine reference base.</returns>
    /// <exception cref="ArgumentException">A base is not A/C/G/T, or the reference and alternate bases are equal.</exception>
    public static string ClassifySbsContext(char fivePrime, char referenceBase, char alternateBase, char threePrime)
    {
        char five = NormalizeBase(fivePrime, nameof(fivePrime));
        char reference = NormalizeBase(referenceBase, nameof(referenceBase));
        char alternate = NormalizeBase(alternateBase, nameof(alternateBase));
        char three = NormalizeBase(threePrime, nameof(threePrime));

        if (reference == alternate)
        {
            throw new ArgumentException(
                $"A substitution requires reference ≠ alternate (got '{reference}' = '{alternate}').",
                nameof(alternateBase));
        }

        // Fold purine-reference substitutions onto the pyrimidine strand by reverse-complementing the
        // trinucleotide context AND the substitution (SigProfiler / COSMIC). For a pyrimidine reference
        // (C or T) the mutation is already on the pyrimidine strand and is kept as-is.
        if (reference is 'A' or 'G')
        {
            char foldedFive = Complement(three);   // 3' neighbour becomes the 5' neighbour after reversal
            char foldedThree = Complement(five);   // 5' neighbour becomes the 3' neighbour after reversal
            reference = Complement(reference);
            alternate = Complement(alternate);
            five = foldedFive;
            three = foldedThree;
        }

        return $"{five}[{reference}>{alternate}]{three}";
    }

    /// <summary>
    /// Enumerates all 96 canonical SBS-96 channel labels (<c>5'[REF&gt;ALT]3'</c>), in deterministic
    /// substitution-major order: the six pyrimidine substitutions (C&gt;A, C&gt;G, C&gt;T, T&gt;A, T&gt;C,
    /// T&gt;G), then 5' base (A,C,G,T), then 3' base (A,C,G,T). Source: COSMIC SBS96; Alexandrov et al. (2013)
    /// — 6 × 4 × 4 = 96. The ordering is a presentation convention and does not affect per-variant classification.
    /// </summary>
    /// <returns>The 96 distinct channel labels.</returns>
    public static IReadOnlyList<string> EnumerateSbs96Channels()
    {
        var channels = new List<string>(Sbs96ChannelCount);
        foreach (var (reference, alternate) in PyrimidineSubstitutions)
        {
            foreach (char five in ContextBases)
            {
                foreach (char three in ContextBases)
                {
                    channels.Add($"{five}[{reference}>{alternate}]{three}");
                }
            }
        }

        return channels;
    }

    /// <summary>
    /// Builds the SBS-96 mutational catalog (the 96-channel spectrum) from a collection of single-base
    /// substitutions, each given as its 5' base, reference (mutated) base, alternate base, and 3' base.
    /// Every variant is classified via <see cref="ClassifySbsContext"/> (with pyrimidine-strand folding) and
    /// tallied into its channel. All 96 channels are present in the result, including those with a zero count,
    /// so the spectrum has a fixed shape. The sum of the counts equals the number of input variants (each
    /// classifiable variant contributes exactly one count — the catalog is a partition). Source: Alexandrov
    /// et al. (2013); COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019).
    /// </summary>
    /// <param name="variants">
    /// SBS variants as (FivePrime, Reference, Alternate, ThreePrime) tuples; each must be a valid single-base
    /// substitution (A/C/G/T bases, reference ≠ alternate).
    /// </param>
    /// <returns>
    /// A dictionary keyed by all 96 channel labels mapping to the count of variants in each channel.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException">Any variant has an invalid base or reference == alternate.</exception>
    public static IReadOnlyDictionary<string, int> Build96ContextCatalog(
        IEnumerable<(char FivePrime, char Reference, char Alternate, char ThreePrime)> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        // Initialise all 96 channels to zero so the spectrum always has the full, fixed shape.
        var catalog = new Dictionary<string, int>(Sbs96ChannelCount, StringComparer.Ordinal);
        foreach (string channel in EnumerateSbs96Channels())
        {
            catalog[channel] = 0;
        }

        foreach (var (fivePrime, reference, alternate, threePrime) in variants)
        {
            string channel = ClassifySbsContext(fivePrime, reference, alternate, threePrime);
            catalog[channel]++;
        }

        return catalog;
    }

    /// <summary>
    /// Returns the Watson-Crick complement of a DNA base (A↔T, C↔G). Source: complementary base pairing,
    /// adenine pairs with thymine and cytosine pairs with guanine.
    /// </summary>
    private static char Complement(char baseChar) => baseChar switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"'{baseChar}' is not a DNA base (A/C/G/T).", nameof(baseChar))
    };

    /// <summary>
    /// Validates and upper-cases a single DNA base, rejecting anything that is not A/C/G/T.
    /// </summary>
    private static char NormalizeBase(char baseChar, string paramName)
    {
        char upper = char.ToUpperInvariant(baseChar);
        if (upper is not ('A' or 'C' or 'G' or 'T'))
        {
            throw new ArgumentException($"'{baseChar}' is not a DNA base (A/C/G/T).", paramName);
        }

        return upper;
    }

    #endregion

    #region Signature Fitting / Refitting (ONCO-SIG-002)

    /// <summary>
    /// Convergence tolerance ε for the Lawson-Hanson active-set NNLS main loop: the iteration stops when the
    /// largest gradient component over the inactive set R is ≤ ε. Source: Lawson C.L. &amp; Hanson R.J. (1974),
    /// <i>Solving Least Squares Problems</i>, Ch. 23 — the active-set algorithm terminates when
    /// max(w_R) ≤ ε (https://en.wikipedia.org/wiki/Non-negative_least_squares). A small positive value
    /// (1e-12) makes the stop effectively exact for the well-conditioned signature matrices used here.
    /// </summary>
    private const double NnlsTolerance = 1e-12;

    /// <summary>
    /// Maximum number of outer (active-set growth) iterations for the NNLS solver. The Lawson-Hanson method
    /// adds at most one index per outer iteration and is guaranteed to terminate in a finite number of steps;
    /// this cap (a small multiple of the number of signatures) is a safety bound against floating-point
    /// non-termination. Source: Lawson &amp; Hanson (1974), Ch. 23 (finite termination of the active-set method).
    /// </summary>
    private const int NnlsMaxOuterIterationsPerSignature = 30;

    /// <summary>
    /// The result of fitting (refitting) an observed mutational catalog to a set of reference signatures.
    /// </summary>
    /// <param name="Exposures">
    /// The fitted non-negative contribution (weight) of each reference signature, in signature order — the
    /// NNLS solution x of min‖S·x − d‖² with x ≥ 0 (Blokzijl et al. 2018; Lawson &amp; Hanson 1974).
    /// </param>
    /// <param name="NormalizedExposures">
    /// <paramref name="Exposures"/> divided by their sum, so they sum to 1 when the total is positive (all
    /// zero otherwise) — the proportion of mutations attributed to each signature (Rosenthal et al. 2016).
    /// </param>
    /// <param name="Reconstruction">The reconstructed catalog S·x (Rosenthal et al. 2016, R = S·W).</param>
    /// <param name="ReconstructionCosineSimilarity">
    /// Cosine similarity between the observed catalog d and the reconstruction S·x — the reconstruction
    /// quality measure (Blokzijl et al. 2018); 1.0 means the catalog is exactly representable.
    /// </param>
    public readonly record struct SignatureFitResult(
        IReadOnlyList<double> Exposures,
        IReadOnlyList<double> NormalizedExposures,
        IReadOnlyList<double> Reconstruction,
        double ReconstructionCosineSimilarity);

    /// <summary>
    /// Computes the cosine similarity between two equal-length non-negative vectors:
    /// <code>sim(A,B) = Σᵢ AᵢBᵢ / ( √(Σᵢ Aᵢ²) · √(Σᵢ Bᵢ²) )</code>
    /// i.e. the dot product divided by the product of the Euclidean norms — the cosine of the angle between
    /// the two vectors. The value lies in [0, 1] for non-negative inputs (0 = independent / orthogonal,
    /// 1 = identical direction) and is invariant to positive scaling of either vector. Source: Blokzijl et al.
    /// (2018), <i>Genome Medicine</i> 10:33 (§ "Mutational profile similarity"); Pan &amp; Wang (2020), iMutSig.
    /// When either vector has zero Euclidean norm the cosine is undefined (division by zero); this method
    /// returns 0.0 for that degenerate case (no shared direction).
    /// </summary>
    /// <param name="a">First vector (e.g. a mutational profile / 96-channel catalog).</param>
    /// <param name="b">Second vector of the same length as <paramref name="a"/>.</param>
    /// <returns>Cosine similarity in [0, 1]; 0.0 when either vector is all-zero.</returns>
    /// <exception cref="ArgumentNullException">Either argument is null.</exception>
    /// <exception cref="ArgumentException">The two vectors have different lengths, or are empty.</exception>
    public static double CosineSimilarity(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        if (a.Count == 0 || b.Count == 0)
        {
            throw new ArgumentException("Cosine similarity is undefined for empty vectors.", nameof(a));
        }

        if (a.Count != b.Count)
        {
            throw new ArgumentException(
                $"Vectors must have the same length (got {a.Count} and {b.Count}).", nameof(b));
        }

        double dot = 0.0;
        double normASquared = 0.0;
        double normBSquared = 0.0;
        for (int i = 0; i < a.Count; i++)
        {
            dot += a[i] * b[i];
            normASquared += a[i] * a[i];
            normBSquared += b[i] * b[i];
        }

        if (normASquared == 0.0 || normBSquared == 0.0)
        {
            // Undefined (zero-norm vector has no direction); treated as no shared direction.
            return 0.0;
        }

        return dot / (Math.Sqrt(normASquared) * Math.Sqrt(normBSquared));
    }

    /// <summary>
    /// Reconstructs a mutational catalog from reference signatures and their exposures: the matrix-vector
    /// product S·x, where each entry is Σⱼ signatures[j][channel] · exposures[j]. Source: Rosenthal et al.
    /// (2016), <i>Genome Biology</i> 17:31 — the reconstructed profile is S·W.
    /// </summary>
    /// <param name="signatures">
    /// Reference signatures as a list of equal-length channel vectors (one vector per signature; element
    /// <c>signatures[j][k]</c> is the weight of signature j in channel k).
    /// </param>
    /// <param name="exposures">The per-signature exposure (weight) of each signature, same count as signatures.</param>
    /// <returns>The reconstructed catalog vector of length equal to the signature channel count.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">No signatures, ragged signatures, or count mismatch with exposures.</exception>
    public static IReadOnlyList<double> ReconstructCatalog(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> exposures)
    {
        int channelCount = ValidateSignatures(signatures);
        ArgumentNullException.ThrowIfNull(exposures);

        if (exposures.Count != signatures.Count)
        {
            throw new ArgumentException(
                $"Exposure count ({exposures.Count}) must equal signature count ({signatures.Count}).",
                nameof(exposures));
        }

        var reconstruction = new double[channelCount];
        for (int j = 0; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j];
            double weight = exposures[j];
            for (int k = 0; k < channelCount; k++)
            {
                reconstruction[k] += signature[k] * weight;
            }
        }

        return reconstruction;
    }

    /// <summary>
    /// Fits (refits) an observed mutational catalog to a set of caller-supplied reference signatures by solving
    /// the non-negative least-squares problem
    /// <code>minₓ ‖ S·x − d ‖₂²,  subject to x ≥ 0</code>
    /// where the columns of S are the reference signatures, d is the observed catalog, and x is the fitted
    /// per-signature exposure (contribution) vector. The decomposition is the standard signature-refitting
    /// model: the catalog is projected onto the non-negative cone spanned by the reference signatures
    /// (Blokzijl et al. 2018). Reference signature profiles are <b>not</b> hardcoded — they are supplied by the
    /// caller (e.g. COSMIC SBS profiles); this method only performs the fit. The NNLS problem is solved with
    /// the Lawson-Hanson active-set algorithm (Lawson &amp; Hanson 1974). The result also exposes the
    /// proportion-normalised exposures (Rosenthal et al. 2016), the reconstruction S·x, and the cosine
    /// similarity between d and S·x as a reconstruction-quality measure (Blokzijl et al. 2018).
    /// </summary>
    /// <param name="catalog">The observed mutational catalog d (e.g. a 96-channel SBS spectrum), non-negative.</param>
    /// <param name="signatures">Reference signatures as equal-length channel vectors (one per signature).</param>
    /// <returns>The fit result: exposures, normalised exposures, reconstruction, and reconstruction cosine.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">
    /// No signatures, ragged signatures, or the catalog length differs from the signature channel count.
    /// </exception>
    public static SignatureFitResult FitSignatures(
        IReadOnlyList<double> catalog,
        IReadOnlyList<IReadOnlyList<double>> signatures)
    {
        ArgumentNullException.ThrowIfNull(catalog);
        int channelCount = ValidateSignatures(signatures);

        if (catalog.Count != channelCount)
        {
            throw new ArgumentException(
                $"Catalog length ({catalog.Count}) must equal the signature channel count ({channelCount}).",
                nameof(catalog));
        }

        double[] exposures = SolveNonNegativeLeastSquares(signatures, catalog, channelCount);

        IReadOnlyList<double> reconstruction = ReconstructCatalog(signatures, exposures);
        double[] normalized = NormalizeExposures(exposures);
        double reconstructionCosine = CosineSimilarity(catalog, reconstruction);

        return new SignatureFitResult(exposures, normalized, reconstruction, reconstructionCosine);
    }

    /// <summary>
    /// Normalises exposures into proportions that sum to 1 (each divided by the total), or all zeros when the
    /// total is zero. Source: Rosenthal et al. (2016) — "the weights W are normalized between 0 and 1".
    /// </summary>
    private static double[] NormalizeExposures(double[] exposures)
    {
        double sum = 0.0;
        foreach (double e in exposures)
        {
            sum += e;
        }

        var normalized = new double[exposures.Length];
        if (sum > 0.0)
        {
            for (int i = 0; i < exposures.Length; i++)
            {
                normalized[i] = exposures[i] / sum;
            }
        }

        return normalized;
    }

    /// <summary>
    /// Solves minₓ ‖ S·x − d ‖₂² subject to x ≥ 0 with the Lawson-Hanson active-set algorithm.
    /// Source: Lawson C.L. &amp; Hanson R.J. (1974), <i>Solving Least Squares Problems</i>, Ch. 23
    /// (https://en.wikipedia.org/wiki/Non-negative_least_squares). Index set P holds the passive (free,
    /// possibly non-zero) variables; R holds the active (clamped-to-zero) variables; the gradient
    /// w = Sᵀ(d − S·x) selects the next variable to free.
    /// </summary>
    /// <param name="signatures">Signature matrix S (column j = signatures[j]).</param>
    /// <param name="catalog">Observed vector d.</param>
    /// <param name="channelCount">Number of channels (rows of S, length of d).</param>
    /// <returns>The NNLS solution x (length = number of signatures).</returns>
    private static double[] SolveNonNegativeLeastSquares(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> catalog,
        int channelCount)
    {
        int n = signatures.Count;
        var x = new double[n];
        bool[] passive = new bool[n]; // true => index in P, false => in R

        int maxOuter = n * NnlsMaxOuterIterationsPerSignature;
        int outer = 0;

        while (outer++ < maxOuter)
        {
            // w = Sᵀ(d − S·x); only inactive (R) components matter for selection.
            double[] residual = ComputeResidual(signatures, x, catalog, channelCount);
            double[] gradient = ComputeGradient(signatures, residual);

            int j = -1;
            double maxGradient = NnlsTolerance;
            for (int i = 0; i < n; i++)
            {
                if (!passive[i] && gradient[i] > maxGradient)
                {
                    maxGradient = gradient[i];
                    j = i;
                }
            }

            if (j < 0)
            {
                // R empty or max(w_R) ≤ ε — KKT conditions satisfied.
                break;
            }

            passive[j] = true;

            // Inner loop: solve the unconstrained LS on P; if any becomes ≤ 0, take the bounded step.
            int innerGuard = 0;
            while (innerGuard++ <= n)
            {
                double[] s = SolveLeastSquaresOnPassiveSet(signatures, catalog, passive, channelCount, n);

                double minPassive = double.PositiveInfinity;
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && s[i] < minPassive)
                    {
                        minPassive = s[i];
                    }
                }

                if (minPassive > 0.0)
                {
                    x = s;
                    break;
                }

                // α = min over i in P with s_i ≤ 0 of x_i / (x_i − s_i).
                double alpha = double.PositiveInfinity;
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && s[i] <= 0.0)
                    {
                        double denom = x[i] - s[i];
                        if (denom != 0.0)
                        {
                            double candidate = x[i] / denom;
                            if (candidate < alpha)
                            {
                                alpha = candidate;
                            }
                        }
                    }
                }

                if (double.IsPositiveInfinity(alpha))
                {
                    // Numerical safeguard: no feasible step (should not occur for valid inputs).
                    x = s;
                    break;
                }

                for (int i = 0; i < n; i++)
                {
                    x[i] += alpha * (s[i] - x[i]);
                }

                // Move indices with x_i ≤ 0 from P back to R.
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && x[i] <= 0.0)
                    {
                        x[i] = 0.0;
                        passive[i] = false;
                    }
                }
            }
        }

        return x;
    }

    /// <summary>Computes the residual d − S·x.</summary>
    private static double[] ComputeResidual(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        double[] x,
        IReadOnlyList<double> catalog,
        int channelCount)
    {
        var residual = new double[channelCount];
        for (int k = 0; k < channelCount; k++)
        {
            residual[k] = catalog[k];
        }

        for (int j = 0; j < signatures.Count; j++)
        {
            double weight = x[j];
            if (weight == 0.0)
            {
                continue;
            }

            IReadOnlyList<double> signature = signatures[j];
            for (int k = 0; k < channelCount; k++)
            {
                residual[k] -= signature[k] * weight;
            }
        }

        return residual;
    }

    /// <summary>Computes the gradient w = Sᵀ·residual.</summary>
    private static double[] ComputeGradient(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        double[] residual)
    {
        var gradient = new double[signatures.Count];
        for (int j = 0; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j];
            double sum = 0.0;
            for (int k = 0; k < residual.Length; k++)
            {
                sum += signature[k] * residual[k];
            }

            gradient[j] = sum;
        }

        return gradient;
    }

    /// <summary>
    /// Solves the unconstrained least-squares problem restricted to the passive set P:
    /// s_P = ((S_P)ᵀ S_P)⁻¹ (S_P)ᵀ d, with s_R = 0, via the normal equations solved by Gaussian elimination.
    /// Source: Lawson &amp; Hanson (1974), Ch. 23.
    /// </summary>
    private static double[] SolveLeastSquaresOnPassiveSet(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> catalog,
        bool[] passive,
        int channelCount,
        int n)
    {
        var passiveIndices = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (passive[i])
            {
                passiveIndices.Add(i);
            }
        }

        int p = passiveIndices.Count;
        var s = new double[n];
        if (p == 0)
        {
            return s;
        }

        // Normal equations: (S_Pᵀ S_P) z = S_Pᵀ d.
        var ata = new double[p, p];
        var atb = new double[p];
        for (int a = 0; a < p; a++)
        {
            IReadOnlyList<double> sigA = signatures[passiveIndices[a]];
            double rhs = 0.0;
            for (int k = 0; k < channelCount; k++)
            {
                rhs += sigA[k] * catalog[k];
            }

            atb[a] = rhs;

            for (int b = 0; b < p; b++)
            {
                IReadOnlyList<double> sigB = signatures[passiveIndices[b]];
                double dot = 0.0;
                for (int k = 0; k < channelCount; k++)
                {
                    dot += sigA[k] * sigB[k];
                }

                ata[a, b] = dot;
            }
        }

        double[] z = SolveLinearSystem(ata, atb, p);
        for (int a = 0; a < p; a++)
        {
            s[passiveIndices[a]] = z[a];
        }

        return s;
    }

    /// <summary>
    /// Solves the dense linear system M·z = rhs (M is p×p, symmetric positive semi-definite here) by Gaussian
    /// elimination with partial pivoting. Standard direct method (CLRS, §28; Numerical Recipes §2.1).
    /// </summary>
    private static double[] SolveLinearSystem(double[,] matrix, double[] rhs, int p)
    {
        // Work on copies so the inputs are not mutated.
        var m = new double[p, p];
        var b = new double[p];
        for (int i = 0; i < p; i++)
        {
            b[i] = rhs[i];
            for (int k = 0; k < p; k++)
            {
                m[i, k] = matrix[i, k];
            }
        }

        for (int col = 0; col < p; col++)
        {
            // Partial pivot: largest magnitude in this column at or below the diagonal.
            int pivot = col;
            double best = Math.Abs(m[col, col]);
            for (int row = col + 1; row < p; row++)
            {
                double magnitude = Math.Abs(m[row, col]);
                if (magnitude > best)
                {
                    best = magnitude;
                    pivot = row;
                }
            }

            if (pivot != col)
            {
                for (int k = 0; k < p; k++)
                {
                    (m[col, k], m[pivot, k]) = (m[pivot, k], m[col, k]);
                }

                (b[col], b[pivot]) = (b[pivot], b[col]);
            }

            double diagonal = m[col, col];
            if (diagonal == 0.0)
            {
                // Singular column (collinear signatures); leave this component at 0.
                continue;
            }

            for (int row = col + 1; row < p; row++)
            {
                double factor = m[row, col] / diagonal;
                if (factor == 0.0)
                {
                    continue;
                }

                for (int k = col; k < p; k++)
                {
                    m[row, k] -= factor * m[col, k];
                }

                b[row] -= factor * b[col];
            }
        }

        // Back-substitution.
        var z = new double[p];
        for (int row = p - 1; row >= 0; row--)
        {
            double sum = b[row];
            for (int k = row + 1; k < p; k++)
            {
                sum -= m[row, k] * z[k];
            }

            double diagonal = m[row, row];
            z[row] = diagonal == 0.0 ? 0.0 : sum / diagonal;
        }

        return z;
    }

    /// <summary>
    /// Validates a signature matrix (non-null, at least one signature, each signature non-null and equal-length,
    /// non-empty) and returns the common channel count.
    /// </summary>
    private static int ValidateSignatures(IReadOnlyList<IReadOnlyList<double>> signatures)
    {
        ArgumentNullException.ThrowIfNull(signatures);

        if (signatures.Count == 0)
        {
            throw new ArgumentException("At least one reference signature is required.", nameof(signatures));
        }

        IReadOnlyList<double> first = signatures[0]
            ?? throw new ArgumentException("Signature vectors cannot be null.", nameof(signatures));
        int channelCount = first.Count;
        if (channelCount == 0)
        {
            throw new ArgumentException("Signature vectors cannot be empty.", nameof(signatures));
        }

        for (int j = 1; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j]
                ?? throw new ArgumentException("Signature vectors cannot be null.", nameof(signatures));
            if (signature.Count != channelCount)
            {
                throw new ArgumentException(
                    $"All signatures must have the same length (signature 0 has {channelCount}, " +
                    $"signature {j} has {signature.Count}).",
                    nameof(signatures));
            }
        }

        return channelCount;
    }

    #endregion

    #region De-novo Signature Extraction via NMF (ONCO-SIG-002)

    /// <summary>
    /// Maximum number of multiplicative-update iterations for the de-novo NMF signature extraction. Source:
    /// Lee &amp; Seung (2001), <i>Algorithms for Non-negative Matrix Factorization</i>, NIPS 13 — the
    /// multiplicative updates are "applied iteratively until W and H converge"
    /// (https://en.wikipedia.org/wiki/Non-negative_matrix_factorization). NMF is non-convex, so a finite cap is
    /// a safety bound; this default mirrors the iteration budgets used by SigProfiler-style extractors.
    /// </summary>
    public const int DefaultNmfMaxIterations = 10_000;

    /// <summary>
    /// Default relative-improvement convergence tolerance for the NMF objective: iteration stops when the
    /// per-iteration decrease of the Frobenius residual ‖V − WH‖²_F, relative to the previous value, drops below
    /// this threshold. Source: Lee &amp; Seung (2001), Theorem 1 — the objective is monotonically non-increasing,
    /// so a small relative-change stop is a valid convergence test. Value = 1e-10.
    /// </summary>
    public const double DefaultNmfTolerance = 1e-10;

    /// <summary>
    /// Default RNG seed for the non-negative random initialisation of the NMF factors. Fixed so that, for a
    /// given count matrix V, rank k, iteration budget and tolerance, the extracted signatures and exposures are
    /// reproducible (NMF is non-convex and initialisation-dependent). Mirrors the fixed-seed convention used by
    /// <see cref="DefaultBootstrapSeed"/>. Value = 42.
    /// </summary>
    public const int DefaultNmfSeed = 42;

    /// <summary>
    /// A small additive floor on the multiplicative-update denominators (and on the random initialisation) to
    /// avoid 0/0 when a row or column of a factor collapses to zero. Source: standard regularisation of the
    /// Lee &amp; Seung multiplicative updates whose denominators (WᵀWH, WHHᵀ) can vanish
    /// (https://en.wikipedia.org/wiki/Non-negative_matrix_factorization). Value = 1e-12, far below any
    /// meaningful mutation count.
    /// </summary>
    private const double NmfEpsilon = 1e-12;

    /// <summary>
    /// The result of de-novo NMF signature extraction from a mutation-count matrix V ≈ W·H at a caller-specified
    /// rank k (Lee &amp; Seung 2001; Alexandrov et al. 2013).
    /// </summary>
    /// <param name="Signatures">
    /// The extracted signatures W as a list of k channel vectors (one vector of length <c>ChannelCount</c> per
    /// signature). Each signature is L1-normalised so its channel weights sum to 1 — a probability distribution
    /// across the mutation channels, per the COSMIC / SigProfiler convention (Alexandrov et al. 2020).
    /// </param>
    /// <param name="Exposures">
    /// The exposures H as a k × samples matrix (<c>Exposures[j][s]</c> = activity of signature j in sample s).
    /// Non-negative; absorbs the per-signature scale removed by L1-normalising <paramref name="Signatures"/>.
    /// </param>
    /// <param name="FinalResidual">
    /// The final Frobenius reconstruction residual ‖V − W·H‖²_F (squared) after the last iteration.
    /// </param>
    /// <param name="Iterations">The number of multiplicative-update iterations actually performed.</param>
    /// <param name="ObjectiveHistory">
    /// The Frobenius residual ‖V − W·H‖²_F recorded after each iteration (length = <see cref="Iterations"/>),
    /// which is monotonically non-increasing (Lee &amp; Seung 2001, Theorem 1).
    /// </param>
    public readonly record struct SignatureExtractionResult(
        IReadOnlyList<IReadOnlyList<double>> Signatures,
        IReadOnlyList<IReadOnlyList<double>> Exposures,
        double FinalResidual,
        int Iterations,
        IReadOnlyList<double> ObjectiveHistory);

    /// <summary>
    /// De-novo mutational-signature extraction by Non-negative Matrix Factorization (NMF). Given a non-negative
    /// mutation-count matrix V (channels × samples) and a caller-specified rank k, factorises
    /// <code>V ≈ W·H,  W ≥ 0  (channels × k),  H ≥ 0  (k × samples)</code>
    /// where the columns of W are the de-novo signatures and H holds their per-sample exposures. Unlike
    /// <see cref="FitSignatures"/> (which refits exposures against caller-supplied reference signatures), the
    /// signatures here are <b>discovered from the data</b> — no reference profiles are required. The factors are
    /// found with the Lee &amp; Seung (2001) multiplicative update rules for the squared Euclidean (Frobenius)
    /// objective ‖V − W·H‖²_F:
    /// <code>
    /// H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)
    /// W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)
    /// </code>
    /// iterated until the relative decrease of the objective falls below <paramref name="tolerance"/> or
    /// <paramref name="maxIterations"/> is reached. Each extracted signature column of W is then L1-normalised to
    /// sum to 1 (a probability distribution over the channels), with the removed scale absorbed into H, per the
    /// COSMIC / SigProfiler convention (Alexandrov et al. 2013, 2020). NMF is non-convex, so the factorisation is
    /// a local optimum dependent on the (seeded, deterministic) non-negative random initialisation.
    /// </summary>
    /// <param name="countMatrix">
    /// The mutation-count matrix V as a list of rows, one per channel (e.g. 96 SBS channels); each row is a
    /// vector over the samples (<c>countMatrix[channel][sample]</c>). All entries must be finite and ≥ 0. Every
    /// row must have the same length (the sample count), and there must be at least one sample.
    /// </param>
    /// <param name="rank">The number of signatures k to extract; 1 ≤ k ≤ channel count.</param>
    /// <param name="maxIterations">Maximum multiplicative-update iterations (&gt; 0).</param>
    /// <param name="tolerance">Relative-improvement convergence tolerance (≥ 0).</param>
    /// <param name="seed">RNG seed for the non-negative random initialisation (reproducibility).</param>
    /// <returns>The extracted signatures W (L1-normalised), exposures H, residual, iteration count and history.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="countMatrix"/> or any row is null.</exception>
    /// <exception cref="ArgumentException">
    /// Empty matrix, zero samples, ragged rows, a negative or non-finite entry, rank &lt; 1, rank &gt; channel
    /// count, maxIterations ≤ 0, or tolerance &lt; 0.
    /// </exception>
    public static SignatureExtractionResult ExtractSignatures(
        IReadOnlyList<IReadOnlyList<double>> countMatrix,
        int rank,
        int maxIterations = DefaultNmfMaxIterations,
        double tolerance = DefaultNmfTolerance,
        int seed = DefaultNmfSeed)
    {
        double[][] v = ValidateCountMatrix(countMatrix, out int channelCount, out int sampleCount);

        if (rank < 1)
        {
            throw new ArgumentException($"Rank k must be ≥ 1 (got {rank}).", nameof(rank));
        }

        if (rank > channelCount)
        {
            throw new ArgumentException(
                $"Rank k ({rank}) cannot exceed the channel count ({channelCount}).", nameof(rank));
        }

        if (maxIterations <= 0)
        {
            throw new ArgumentException($"maxIterations must be > 0 (got {maxIterations}).", nameof(maxIterations));
        }

        if (tolerance < 0)
        {
            throw new ArgumentException($"tolerance must be ≥ 0 (got {tolerance}).", nameof(tolerance));
        }

        // Non-negative random initialisation (Lee & Seung do not prescribe one; uniform (0,1] is standard).
        var rng = new Random(seed);
        double[][] w = InitializeNonNegativeFactor(rng, channelCount, rank);   // channels × k
        double[][] h = InitializeNonNegativeFactor(rng, rank, sampleCount);    // k × samples

        var history = new List<double>(maxIterations);
        double previousObjective = double.PositiveInfinity;
        int iterations = 0;

        for (int iter = 0; iter < maxIterations; iter++)
        {
            // H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)   — Lee & Seung (2001), Theorem 1.
            UpdateH(w, h, v, channelCount, rank, sampleCount);
            // W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)   — Lee & Seung (2001), Theorem 1.
            UpdateW(w, h, v, channelCount, rank, sampleCount);

            iterations = iter + 1;
            double objective = FrobeniusResidualSquared(w, h, v, channelCount, rank, sampleCount);
            history.Add(objective);

            // Relative-improvement stop. The objective is monotonically non-increasing (Theorem 1), so
            // previousObjective ≥ objective; the decrease is non-negative.
            double decrease = previousObjective - objective;
            double denominator = previousObjective > NmfEpsilon ? previousObjective : 1.0;
            if (!double.IsInfinity(previousObjective) && decrease / denominator < tolerance)
            {
                previousObjective = objective;
                break;
            }

            previousObjective = objective;
        }

        // L1-normalise each signature column of W so its channel weights sum to 1, absorbing the scale into the
        // corresponding row of H (Alexandrov et al. 2013/2020; COSMIC SBS — signatures are probability
        // distributions over the channels). This fixes the NMF scaling ambiguity without changing W·H.
        NormalizeSignatureColumns(w, h, channelCount, rank, sampleCount);

        double finalResidual = FrobeniusResidualSquared(w, h, v, channelCount, rank, sampleCount);

        IReadOnlyList<IReadOnlyList<double>> signatures = TransposeColumnsToSignatures(w, channelCount, rank);
        IReadOnlyList<IReadOnlyList<double>> exposures = RowsToReadOnly(h);

        return new SignatureExtractionResult(signatures, exposures, finalResidual, iterations, history);
    }

    /// <summary>
    /// Validates the count matrix V (non-null, non-empty, rectangular, finite, non-negative) and returns it as a
    /// dense jagged array of rows, with the channel and sample counts.
    /// </summary>
    private static double[][] ValidateCountMatrix(
        IReadOnlyList<IReadOnlyList<double>> countMatrix, out int channelCount, out int sampleCount)
    {
        ArgumentNullException.ThrowIfNull(countMatrix);

        channelCount = countMatrix.Count;
        if (channelCount == 0)
        {
            throw new ArgumentException("The count matrix must have at least one channel (row).", nameof(countMatrix));
        }

        IReadOnlyList<double> firstRow = countMatrix[0]
            ?? throw new ArgumentException("Count-matrix rows cannot be null.", nameof(countMatrix));
        sampleCount = firstRow.Count;
        if (sampleCount == 0)
        {
            throw new ArgumentException("The count matrix must have at least one sample (column).", nameof(countMatrix));
        }

        var v = new double[channelCount][];
        for (int c = 0; c < channelCount; c++)
        {
            IReadOnlyList<double> row = countMatrix[c]
                ?? throw new ArgumentException("Count-matrix rows cannot be null.", nameof(countMatrix));
            if (row.Count != sampleCount)
            {
                throw new ArgumentException(
                    $"All channel rows must have the same sample count (row 0 has {sampleCount}, " +
                    $"row {c} has {row.Count}).",
                    nameof(countMatrix));
            }

            var dense = new double[sampleCount];
            for (int s = 0; s < sampleCount; s++)
            {
                double value = row[s];
                if (double.IsNaN(value) || double.IsInfinity(value) || value < 0)
                {
                    throw new ArgumentException(
                        $"Count-matrix entries must be finite and ≥ 0 (row {c}, sample {s} = {value}).",
                        nameof(countMatrix));
                }

                dense[s] = value;
            }

            v[c] = dense;
        }

        return v;
    }

    /// <summary>
    /// Builds a non-negative factor of the given shape with entries drawn uniformly from (0, 1], floored by
    /// <see cref="NmfEpsilon"/> so no entry is exactly zero (a zero row/column cannot recover under the
    /// multiplicative updates). Source: standard non-negative random initialisation for Lee &amp; Seung NMF.
    /// </summary>
    private static double[][] InitializeNonNegativeFactor(Random rng, int rows, int cols)
    {
        var factor = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            var rowValues = new double[cols];
            for (int j = 0; j < cols; j++)
            {
                rowValues[j] = rng.NextDouble() + NmfEpsilon;
            }

            factor[i] = rowValues;
        }

        return factor;
    }

    /// <summary>
    /// Multiplicative update of H: <c>H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)</c>. Source: Lee &amp; Seung (2001), Theorem 1
    /// (Euclidean objective). The denominator is floored by <see cref="NmfEpsilon"/> to avoid 0/0.
    /// </summary>
    private static void UpdateH(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        // numerator = Wᵀ V  (k × samples); denominator = Wᵀ (W H)  (k × samples).
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        for (int a = 0; a < k; a++)
        {
            for (int s = 0; s < samples; s++)
            {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int c = 0; c < channels; c++)
                {
                    numerator += w[c][a] * v[c][s];
                    denominator += w[c][a] * wh[c][s];
                }

                if (denominator < NmfEpsilon)
                {
                    denominator = NmfEpsilon;
                }

                h[a][s] *= numerator / denominator;
            }
        }
    }

    /// <summary>
    /// Multiplicative update of W: <c>W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)</c>. Source: Lee &amp; Seung (2001), Theorem 1
    /// (Euclidean objective). The denominator is floored by <see cref="NmfEpsilon"/> to avoid 0/0.
    /// </summary>
    private static void UpdateW(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        // numerator = V Hᵀ  (channels × k); denominator = (W H) Hᵀ  (channels × k).
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        for (int c = 0; c < channels; c++)
        {
            for (int a = 0; a < k; a++)
            {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int s = 0; s < samples; s++)
                {
                    numerator += v[c][s] * h[a][s];
                    denominator += wh[c][s] * h[a][s];
                }

                if (denominator < NmfEpsilon)
                {
                    denominator = NmfEpsilon;
                }

                w[c][a] *= numerator / denominator;
            }
        }
    }

    /// <summary>
    /// Computes the dense product W·H (channels × samples) of the current factors.
    /// </summary>
    private static double[][] MultiplyWh(double[][] w, double[][] h, int channels, int k, int samples)
    {
        var wh = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            var row = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                double sum = 0.0;
                for (int a = 0; a < k; a++)
                {
                    sum += w[c][a] * h[a][s];
                }

                row[s] = sum;
            }

            wh[c] = row;
        }

        return wh;
    }

    /// <summary>
    /// The squared Frobenius reconstruction residual ‖V − W·H‖²_F = Σ (V − W·H)². Source: Lee &amp; Seung (2001),
    /// the Euclidean objective F(W,H) = ‖V − WH‖²_F.
    /// </summary>
    private static double FrobeniusResidualSquared(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        double sum = 0.0;
        for (int c = 0; c < channels; c++)
        {
            for (int s = 0; s < samples; s++)
            {
                double reconstructed = 0.0;
                for (int a = 0; a < k; a++)
                {
                    reconstructed += w[c][a] * h[a][s];
                }

                double diff = v[c][s] - reconstructed;
                sum += diff * diff;
            }
        }

        return sum;
    }

    /// <summary>
    /// L1-normalises each signature column of W so its channel weights sum to 1, absorbing the removed scale into
    /// the matching row of H so that W·H is unchanged. Source: Alexandrov et al. (2013/2020); COSMIC SBS — each
    /// signature is a probability distribution over the channels. A column that sums to zero is left as zeros.
    /// </summary>
    private static void NormalizeSignatureColumns(double[][] w, double[][] h, int channels, int k, int samples)
    {
        for (int a = 0; a < k; a++)
        {
            double columnSum = 0.0;
            for (int c = 0; c < channels; c++)
            {
                columnSum += w[c][a];
            }

            if (columnSum <= NmfEpsilon)
            {
                continue;
            }

            for (int c = 0; c < channels; c++)
            {
                w[c][a] /= columnSum;
            }

            // Absorb the scale into H so W·H is invariant: H[a][·] *= columnSum.
            for (int s = 0; s < samples; s++)
            {
                h[a][s] *= columnSum;
            }
        }
    }

    /// <summary>
    /// Converts the channels × k factor W into k signature vectors (column a of W → signature a of length
    /// <paramref name="channels"/>), the per-signature channel-vector layout used elsewhere in this class.
    /// </summary>
    private static IReadOnlyList<IReadOnlyList<double>> TransposeColumnsToSignatures(double[][] w, int channels, int k)
    {
        var signatures = new List<IReadOnlyList<double>>(k);
        for (int a = 0; a < k; a++)
        {
            var signature = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                signature[c] = w[c][a];
            }

            signatures.Add(signature);
        }

        return signatures;
    }

    /// <summary>
    /// Wraps the rows of H as a read-only list of read-only rows.
    /// </summary>
    private static IReadOnlyList<IReadOnlyList<double>> RowsToReadOnly(double[][] h)
    {
        var rows = new List<IReadOnlyList<double>>(h.Length);
        foreach (double[] row in h)
        {
            rows.Add(row);
        }

        return rows;
    }

    #endregion

    #region Signature Exposure Bootstrap Confidence Intervals (ONCO-SIG-003)

    /// <summary>
    /// Default number of bootstrap replicates for exposure confidence-interval estimation. Source:
    /// Senkin S. (2021), MSA — confidence intervals are derived from "1000 bootstrap variations"
    /// (https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/); Wang S. et al., sigminer
    /// <c>sig_fit_bootstrap</c> recommends "Bootstrap replicates &gt;= 100"
    /// (https://shixiangwang.github.io/sigminer-doc/sigfit.html). Default = 1000.
    /// </summary>
    public const int DefaultBootstrapReplicates = 1000;

    /// <summary>
    /// Default RNG seed for the multinomial resampling step. Fixed so that, for a given catalog,
    /// signatures, replicate count and confidence level, the returned intervals are reproducible
    /// (deterministic test/clinical re-runs). Mirrors the fixed-seed convention used by
    /// <see cref="Seqeron.Genomics.Phylogenetics"/> bootstrap. Value = 42.
    /// </summary>
    public const int DefaultBootstrapSeed = 42;

    /// <summary>
    /// Default two-sided confidence level for the percentile bootstrap interval. Source: Efron B. (1979),
    /// percentile method — a 95% interval is the [2.5%, 97.5%] percentiles of the bootstrap distribution
    /// (Senkin 2021, MSA: "95% confidence intervals … taking [2.5%, 97.5%] percentiles of the resulting
    /// bootstrap activities", https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/). Default = 0.95.
    /// </summary>
    public const double DefaultBootstrapConfidence = 0.95;

    /// <summary>
    /// Selects how each bootstrap replicate of the mutational catalog is resampled in
    /// <see cref="BootstrapExposures"/>. Both schemes are described by Senkin (2021), MSA
    /// (<i>BMC Bioinformatics</i> 22:540, https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/), which notes that
    /// "the conditional distribution of a vector of independent Poisson variables is equivalent to multinomial
    /// distribution" — the two differ only in whether the total mutational burden N is held fixed.
    /// </summary>
    public enum BootstrapResampling
    {
        /// <summary>
        /// Fixed-N multinomial resample: the observed total N = Σ catalog mutations is redistributed across
        /// channels with probabilities pₖ = catalogₖ / N, so every replicate has exactly N mutations.
        /// This is the sigminer <c>sig_fit_bootstrap</c> scheme,
        /// <c>sample(K, total, replace=TRUE, prob=catalog/sum(catalog))</c>
        /// (https://github.com/ShixiangWang/sigminer/blob/master/R/sig_fit_bootstrap.R). Default — preserves the
        /// historical behaviour of <see cref="BootstrapExposures"/> byte-for-byte.
        /// </summary>
        Multinomial = 0,

        /// <summary>
        /// Poisson resample (Senkin 2021, MSA): each channel count is drawn independently as
        /// Poisson(observedₖ), so "for any given mutation category … the distribution of bootstrapped mutation
        /// counts follows a Poisson distribution" and "the total mutational burden is no longer fixed". This is
        /// the Poisson noise variant of the MSA parametric bootstrap, where the variance of each channel equals
        /// its mean (the observed count).
        /// </summary>
        Poisson = 1,
    }

    /// <summary>
    /// A bootstrap confidence interval for one signature's exposure (activity), produced by
    /// <see cref="BootstrapExposures"/>.
    /// </summary>
    /// <param name="PointEstimate">
    /// The exposure for this signature from the NNLS fit of the <i>observed</i> (un-resampled) catalog —
    /// the point estimate the interval is centred on (Senkin 2021; Huang et al. 2018).
    /// </param>
    /// <param name="Mean">The mean of this signature's exposure across the bootstrap replicates.</param>
    /// <param name="Lower">
    /// Lower bound of the percentile interval — the (½(1−c))·100-th percentile of the replicate exposures
    /// (2.5th percentile for c = 0.95). Source: Efron (1979) percentile method; Senkin (2021).
    /// </param>
    /// <param name="Upper">
    /// Upper bound — the (1−½(1−c))·100-th percentile of the replicate exposures (97.5th for c = 0.95).
    /// </param>
    /// <param name="Confidence">The two-sided confidence level c the interval was computed at.</param>
    public readonly record struct ExposureConfidenceInterval(
        double PointEstimate,
        double Mean,
        double Lower,
        double Upper,
        double Confidence);

    /// <summary>
    /// Estimates per-signature exposure confidence intervals by the parametric bootstrap: the observed integer
    /// mutational catalog is repeatedly resampled, each resampled catalog is refit to the reference signatures
    /// by NNLS (<see cref="FitSignatures(IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>), and a
    /// two-sided percentile confidence interval is taken per signature from the resulting bootstrap exposure
    /// distribution. The point estimate is the NNLS exposure of the un-resampled observed catalog. The
    /// <paramref name="resampling"/> parameter selects the resampling scheme:
    /// <list type="bullet">
    /// <item><see cref="BootstrapResampling.Multinomial"/> (default) — each replicate is a draw of
    /// N = Σ catalog mutations from the multinomial distribution with per-channel probabilities
    /// pₖ = catalogₖ / N (fixed total N).</item>
    /// <item><see cref="BootstrapResampling.Poisson"/> — each channel count is drawn independently as
    /// Poisson(observedₖ); the total burden is no longer fixed (Senkin 2021, MSA Poisson variant).</item>
    /// </list>
    /// <para>
    /// Sources: Huang X., Wojtowicz D., Przytycka T.M. (2018), <i>Bioinformatics</i> 34(2):330–337 — bootstrap
    /// resampling of the mutation-count vector to assess decomposition confidence; Senkin S. (2021), MSA,
    /// <i>BMC Bioinformatics</i> 22:540 — "mutations are accumulated following Poisson distributions for each
    /// mutation class", "drawing counts from independent binomial distributions, so that the total mutational
    /// burden is no longer fixed … for any given mutation category … the distribution of bootstrapped mutation
    /// counts follows a Poisson distribution", "the conditional distribution of a vector of independent Poisson
    /// variables is equivalent to multinomial distribution", "95% confidence intervals … taking [2.5%, 97.5%]
    /// percentiles"; Wang S. et al., sigminer <c>sig_fit_bootstrap</c> — resample via
    /// <c>sample(K, total, replace=TRUE, prob=catalog/sum(catalog))</c> (a multinomial draw) then refit;
    /// Efron B. (1979) percentile interval. Reference signature profiles are caller-supplied (not fabricated).
    /// </para>
    /// </summary>
    /// <param name="catalog">
    /// The observed mutational catalog as non-negative integer per-channel mutation counts (e.g. a 96-channel
    /// SBS spectrum). Counts (not proportions) are required because each resample needs the per-channel totals.
    /// </param>
    /// <param name="signatures">Reference signatures as equal-length channel vectors (one per signature).</param>
    /// <param name="replicates">Number of bootstrap replicates (≥ 1; default <see cref="DefaultBootstrapReplicates"/>).</param>
    /// <param name="confidence">Two-sided confidence level in (0, 1); default <see cref="DefaultBootstrapConfidence"/>.</param>
    /// <param name="seed">RNG seed for the resampling; fixed value makes results reproducible.</param>
    /// <param name="resampling">
    /// The bootstrap resampling scheme; <see cref="BootstrapResampling.Multinomial"/> (the historical default,
    /// fixed N) or <see cref="BootstrapResampling.Poisson"/> (Senkin 2021 Poisson variant).
    /// </param>
    /// <returns>One <see cref="ExposureConfidenceInterval"/> per signature, in signature order.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">
    /// No signatures, ragged signatures, catalog length ≠ channel count, or a negative catalog count.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// replicates &lt; 1, confidence outside (0, 1), or an unrecognised <paramref name="resampling"/> value.
    /// </exception>
    public static IReadOnlyList<ExposureConfidenceInterval> BootstrapExposures(
        IReadOnlyList<int> catalog,
        IReadOnlyList<IReadOnlyList<double>> signatures,
        int replicates = DefaultBootstrapReplicates,
        double confidence = DefaultBootstrapConfidence,
        int seed = DefaultBootstrapSeed,
        BootstrapResampling resampling = BootstrapResampling.Multinomial)
    {
        ArgumentNullException.ThrowIfNull(catalog);
        int channelCount = ValidateSignatures(signatures);

        if (catalog.Count != channelCount)
        {
            throw new ArgumentException(
                $"Catalog length ({catalog.Count}) must equal the signature channel count ({channelCount}).",
                nameof(catalog));
        }

        if (replicates < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(replicates), replicates, "At least one bootstrap replicate is required.");
        }

        if (double.IsNaN(confidence) || confidence <= 0.0 || confidence >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(confidence), confidence, "Confidence must be in the open interval (0, 1).");
        }

        if (resampling != BootstrapResampling.Multinomial && resampling != BootstrapResampling.Poisson)
        {
            throw new ArgumentOutOfRangeException(
                nameof(resampling), resampling, "Unrecognised bootstrap resampling scheme.");
        }

        // Observed counts and total N = Σ catalog (the multinomial sample size).
        long total = 0;
        var observed = new double[channelCount];
        for (int k = 0; k < channelCount; k++)
        {
            int count = catalog[k];
            if (count < 0)
            {
                throw new ArgumentException("Catalog counts cannot be negative.", nameof(catalog));
            }

            observed[k] = count;
            total += count;
        }

        int signatureCount = signatures.Count;

        // Point estimate: NNLS exposures of the observed (un-resampled) catalog (Senkin 2021; Huang 2018).
        IReadOnlyList<double> pointEstimate = FitSignatures(observed, signatures).Exposures;

        // Per-signature bootstrap exposure distributions.
        var replicateExposures = new double[signatureCount][];
        for (int j = 0; j < signatureCount; j++)
        {
            replicateExposures[j] = new double[replicates];
        }

        var random = new Random(seed);
        var resampled = new double[channelCount];
        for (int rep = 0; rep < replicates; rep++)
        {
            if (resampling == BootstrapResampling.Poisson)
            {
                PoissonResample(observed, random, resampled);
            }
            else
            {
                MultinomialResample(observed, total, random, resampled);
            }

            IReadOnlyList<double> exposures = FitSignatures(resampled, signatures).Exposures;
            for (int j = 0; j < signatureCount; j++)
            {
                replicateExposures[j][rep] = exposures[j];
            }
        }

        double lowerProbability = (1.0 - confidence) / 2.0;
        double upperProbability = 1.0 - lowerProbability;

        var intervals = new ExposureConfidenceInterval[signatureCount];
        for (int j = 0; j < signatureCount; j++)
        {
            double[] distribution = replicateExposures[j];
            double mean = Mean(distribution);
            double lower = Percentile(distribution, lowerProbability);
            double upper = Percentile(distribution, upperProbability);
            intervals[j] = new ExposureConfidenceInterval(
                pointEstimate[j], mean, lower, upper, confidence);
        }

        return intervals;
    }

    /// <summary>
    /// Draws one multinomial resample of the catalog: N = <paramref name="total"/> mutations are assigned to
    /// channels with probabilities pₖ = observedₖ / N, written into <paramref name="destination"/>. Implements
    /// the standard sequential conditional-binomial construction of a multinomial draw: channel k receives
    /// Binomial(remaining, pₖ / Σ_{i≥k} pᵢ) of the mutations not yet assigned. Source: Senkin (2021), MSA
    /// (multinomial/Poisson resampling); sigminer <c>sig_fit_bootstrap</c>
    /// (<c>sample(..., replace=TRUE, prob=catalog/sum(catalog))</c>). When N = 0 the resample is all zeros.
    /// </summary>
    private static void MultinomialResample(double[] observed, long total, Random random, double[] destination)
    {
        int channelCount = observed.Length;
        Array.Clear(destination, 0, channelCount);

        if (total == 0)
        {
            return;
        }

        long remaining = total;
        double remainingProbabilityMass = total; // Σ observed = N (un-normalized probability mass).

        for (int k = 0; k < channelCount && remaining > 0; k++)
        {
            double weight = observed[k];
            if (weight > 0.0)
            {
                // Conditional probability of falling in channel k among the not-yet-assigned channels.
                double p = remainingProbabilityMass > 0.0 ? weight / remainingProbabilityMass : 0.0;
                if (p >= 1.0)
                {
                    // Last channel with any mass: it takes all remaining draws.
                    destination[k] = remaining;
                    remaining = 0;
                }
                else
                {
                    long drawn = SampleBinomial(remaining, p, random);
                    destination[k] = drawn;
                    remaining -= drawn;
                }
            }

            remainingProbabilityMass -= weight;
        }
    }

    /// <summary>
    /// Draws one Poisson resample of the catalog: each channel k is independently assigned
    /// Poisson(observedₖ) mutations, written into <paramref name="destination"/>. Implements the Senkin (2021)
    /// MSA Poisson-noise variant — "for any given mutation category … the distribution of bootstrapped mutation
    /// counts follows a Poisson distribution" with mean equal to the observed count (variance = mean) — so the
    /// total mutational burden is not fixed across replicates. Source: Senkin S. (2021), MSA,
    /// <i>BMC Bioinformatics</i> 22:540 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/). A channel whose
    /// observed count is 0 has mean 0 and therefore always resamples to 0.
    /// </summary>
    private static void PoissonResample(double[] observed, Random random, double[] destination)
    {
        int channelCount = observed.Length;
        for (int k = 0; k < channelCount; k++)
        {
            destination[k] = SamplePoisson(observed[k], random);
        }
    }

    /// <summary>
    /// Samples a Poisson(lambda) variate by Knuth's multiplication-of-uniforms algorithm: draw uniforms until
    /// their running product falls below e^(−lambda); the number of draws minus one is the variate. Source:
    /// Knuth D.E., <i>The Art of Computer Programming</i>, Vol. 2 (Seminumerical Algorithms), §3.4.1 — the
    /// standard exact algorithm for generating Poisson deviates. lambda ≤ 0 returns 0 (degenerate mean-0 case).
    /// </summary>
    private static long SamplePoisson(double lambda, Random random)
    {
        if (lambda <= 0.0)
        {
            return 0;
        }

        double limit = Math.Exp(-lambda);
        long count = 0;
        double product = random.NextDouble();
        while (product > limit)
        {
            count++;
            product *= random.NextDouble();
        }

        return count;
    }

    /// <summary>
    /// Samples a Binomial(n, p) variate by summing n independent Bernoulli(p) trials. Exact for the
    /// per-channel sample sizes encountered in mutational catalogs. Source: standard definition of the
    /// Binomial distribution as the number of successes in n independent Bernoulli(p) trials.
    /// </summary>
    private static long SampleBinomial(long n, double p, Random random)
    {
        if (p <= 0.0)
        {
            return 0;
        }

        if (p >= 1.0)
        {
            return n;
        }

        long successes = 0;
        for (long i = 0; i < n; i++)
        {
            if (random.NextDouble() < p)
            {
                successes++;
            }
        }

        return successes;
    }

    /// <summary>Arithmetic mean of a non-empty array.</summary>
    private static double Mean(double[] values)
    {
        double sum = 0.0;
        foreach (double value in values)
        {
            sum += value;
        }

        return sum / values.Length;
    }

    /// <summary>
    /// Computes the empirical percentile (quantile) of <paramref name="values"/> at probability
    /// <paramref name="probability"/> ∈ [0, 1] using linear interpolation between order statistics on the
    /// 0-based rank h = probability·(n − 1) (the "linear interpolation of the modes for the order statistics"
    /// /R type-7 / NumPy default convention): result = x₍⌊h⌋₎ + (h − ⌊h⌋)·(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎) over the sorted
    /// values. Source: Hyndman R.J. &amp; Fan Y. (1996), <i>The American Statistician</i> 50(4):361–365
    /// (sample-quantile type 7); used to realise the Efron (1979) percentile bootstrap interval.
    /// </summary>
    private static double Percentile(double[] values, double probability)
    {
        int n = values.Length;
        if (n == 1)
        {
            return values[0];
        }

        var sorted = (double[])values.Clone();
        Array.Sort(sorted);

        double rank = probability * (n - 1);
        int lowerIndex = (int)Math.Floor(rank);
        int upperIndex = (int)Math.Ceiling(rank);
        double fraction = rank - lowerIndex;

        if (lowerIndex == upperIndex)
        {
            return sorted[lowerIndex];
        }

        return sorted[lowerIndex] + fraction * (sorted[upperIndex] - sorted[lowerIndex]);
    }

    #endregion

    #region Mutational Process Classification (ONCO-SIG-004)

    /// <summary>
    /// Minimum normalized relative contribution for a single mutational signature to be reported as
    /// present/active. A signature whose contribution falls below this fraction is excluded (set to zero).
    /// Source: Rosenthal R. et al. (2016), deconstructSigs, <i>Genome Biology</i> 17:31 — "the weights W
    /// are normalized between 0 and 1 and any signature with Wᵢ &lt; 6% is excluded"
    /// (https://doi.org/10.1186/s13059-016-0893-4); reference implementation <c>whichSignatures.R</c>
    /// declares <c>signature.cutoff = 0.06</c> and applies <c>weights[weights &lt; signature.cutoff] &lt;- 0</c>.
    /// The comparison is strict less-than, so a contribution of exactly 0.06 is retained. Value = 0.06.
    /// </summary>
    public const double DefaultSignatureContributionCutoff = 0.06;

    /// <summary>
    /// A recognized mutational process (mutagenic aetiology) inferred from active COSMIC SBS signatures.
    /// Aetiology assignments are from the COSMIC SBS catalogue (https://cancer.sanger.ac.uk/signatures/sbs/;
    /// Alexandrov et al. 2020, <i>Nature</i> 578:94–101).
    /// </summary>
    public enum MutationalProcess
    {
        /// <summary>Signature label not mapped to any recognized COSMIC aetiology.</summary>
        Unknown = 0,

        /// <summary>
        /// Aging / clock-like mutagenesis. COSMIC SBS1 ("Spontaneous deamination of 5-methylcytosine
        /// (clock-like signature)") and SBS5 ("Unknown (clock-like signature)").
        /// </summary>
        Aging,

        /// <summary>
        /// APOBEC cytidine-deaminase activity. COSMIC SBS2 and SBS13 ("Activity of APOBEC family of
        /// cytidine deaminases").
        /// </summary>
        Apobec,

        /// <summary>Tobacco smoking. COSMIC SBS4 ("Tobacco smoking").</summary>
        TobaccoSmoking,

        /// <summary>Ultraviolet light exposure. COSMIC SBS7a/7b/7c/7d ("Ultraviolet light exposure").</summary>
        UltravioletLight,

        /// <summary>
        /// Defective DNA mismatch repair. COSMIC SBS6, SBS15, SBS26 ("Defective DNA mismatch repair") and
        /// SBS20 ("Concurrent POLD1 mutations and defective DNA mismatch repair").
        /// </summary>
        MismatchRepairDeficiency,
    }

    /// <summary>
    /// COSMIC SBS signature label → mutational process map, taken verbatim from the COSMIC SBS catalogue
    /// proposed-aetiology strings (https://cancer.sanger.ac.uk/signatures/sbs/; Alexandrov et al. 2020).
    /// Only the five canonical processes named for ONCO-SIG-004 are mapped: Aging (SBS1, SBS5),
    /// APOBEC (SBS2, SBS13), Tobacco smoking (SBS4), UV (SBS7a–d), MMR deficiency (SBS6, SBS15, SBS20, SBS26).
    /// Labels are matched case-insensitively. Signatures outside this map resolve to
    /// <see cref="MutationalProcess.Unknown"/>.
    /// </summary>
    private static readonly IReadOnlyDictionary<string, MutationalProcess> SbsToProcess =
        new Dictionary<string, MutationalProcess>(StringComparer.OrdinalIgnoreCase)
        {
            ["SBS1"] = MutationalProcess.Aging,
            ["SBS5"] = MutationalProcess.Aging,
            ["SBS2"] = MutationalProcess.Apobec,
            ["SBS13"] = MutationalProcess.Apobec,
            ["SBS4"] = MutationalProcess.TobaccoSmoking,
            ["SBS7a"] = MutationalProcess.UltravioletLight,
            ["SBS7b"] = MutationalProcess.UltravioletLight,
            ["SBS7c"] = MutationalProcess.UltravioletLight,
            ["SBS7d"] = MutationalProcess.UltravioletLight,
            ["SBS6"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS15"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS20"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS26"] = MutationalProcess.MismatchRepairDeficiency,
        };

    /// <summary>
    /// Resolves a COSMIC SBS signature label (e.g. <c>"SBS2"</c>, <c>"SBS7a"</c>) to its mutational process
    /// per the COSMIC proposed-aetiology assignments. Matching is case-insensitive. Unmapped or
    /// unknown-aetiology labels return <see cref="MutationalProcess.Unknown"/>. Source: COSMIC SBS catalogue
    /// (https://cancer.sanger.ac.uk/signatures/sbs/; Alexandrov et al. 2020, <i>Nature</i> 578:94–101).
    /// </summary>
    /// <param name="signatureLabel">A COSMIC SBS signature label.</param>
    /// <returns>The mapped <see cref="MutationalProcess"/>, or <see cref="MutationalProcess.Unknown"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="signatureLabel"/> is null.</exception>
    public static MutationalProcess GetMutationalProcess(string signatureLabel)
    {
        ArgumentNullException.ThrowIfNull(signatureLabel);
        return SbsToProcess.TryGetValue(signatureLabel, out MutationalProcess process)
            ? process
            : MutationalProcess.Unknown;
    }

    /// <summary>
    /// The aggregated activity of one mutational process within a tumour sample, produced by
    /// <see cref="ClassifyMutationalProcess"/>.
    /// </summary>
    /// <param name="Process">The mutational process (aetiology).</param>
    /// <param name="Contribution">
    /// The summed normalized relative contribution of all of this process's signatures that survived the
    /// per-signature cutoff (a fraction in [0, 1]). Source: additive deconstructSigs weights (Rosenthal 2016).
    /// </param>
    public readonly record struct ProcessActivity(MutationalProcess Process, double Contribution);

    /// <summary>
    /// The result of classifying which mutational processes are active in a tumour sample.
    /// </summary>
    /// <param name="ActiveProcesses">
    /// The active processes (those with at least one surviving signature), in descending contribution order
    /// then by process for deterministic ordering. Empty when no signature survives the cutoff.
    /// </param>
    /// <param name="DominantProcess">
    /// The active process with the largest aggregated contribution, or <see cref="MutationalProcess.Unknown"/>
    /// when no process is active.
    /// </param>
    public readonly record struct MutationalProcessClassification(
        IReadOnlyList<ProcessActivity> ActiveProcesses,
        MutationalProcess DominantProcess);

    /// <summary>
    /// Classifies the mutational processes active in a tumour from per-signature exposures (activities).
    /// <para>
    /// The exposures are first converted to normalized relative contributions (each divided by their sum,
    /// so they sum to 1) per deconstructSigs ("the weights W are normalized between 0 and 1"). A signature
    /// is then declared present/active only when its normalized contribution is at least
    /// <paramref name="contributionCutoff"/> (default 6%); contributions below the cutoff are dropped, so
    /// the surviving contributions can sum to less than 1 (Rosenthal et al. 2016, deconstructSigs:
    /// "any signature with Wᵢ &lt; 6% is excluded"; <c>weights[weights &lt; signature.cutoff] &lt;- 0</c>).
    /// Each surviving signature is mapped to its mutational process via the COSMIC SBS aetiology map
    /// (<see cref="GetMutationalProcess"/>), and the per-signature contributions are summed per process.
    /// The dominant process is the active process with the largest aggregated contribution.
    /// </para>
    /// <para>
    /// Reference signature <i>profiles</i> are caller-supplied elsewhere (e.g. COSMIC SBS matrices used by
    /// <see cref="FitSignatures(IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>); this method
    /// consumes only the resulting exposures and their COSMIC labels.
    /// </para>
    /// Sources: Rosenthal R. et al. (2016), <i>Genome Biology</i> 17:31
    /// (https://doi.org/10.1186/s13059-016-0893-4); COSMIC SBS catalogue
    /// (https://cancer.sanger.ac.uk/signatures/sbs/); Alexandrov L.B. et al. (2020), <i>Nature</i> 578:94–101.
    /// </summary>
    /// <param name="exposures">
    /// Per-signature exposures as (COSMIC SBS label, non-negative activity) pairs. Activities are arbitrary
    /// non-negative magnitudes (mutation counts or proportions); they are normalized internally.
    /// </param>
    /// <param name="contributionCutoff">
    /// Minimum normalized relative contribution (strict lower bound is a contribution that is &lt; this value)
    /// for a signature to be active; default <see cref="DefaultSignatureContributionCutoff"/> (0.06).
    /// </param>
    /// <returns>The active processes (descending contribution) and the dominant process.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="exposures"/> (or a label) is null.</exception>
    /// <exception cref="ArgumentException">An exposure is negative.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="contributionCutoff"/> is outside [0, 1).</exception>
    public static MutationalProcessClassification ClassifyMutationalProcess(
        IReadOnlyList<(string Signature, double Exposure)> exposures,
        double contributionCutoff = DefaultSignatureContributionCutoff)
    {
        ArgumentNullException.ThrowIfNull(exposures);

        if (double.IsNaN(contributionCutoff) || contributionCutoff < 0.0 || contributionCutoff >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(contributionCutoff), contributionCutoff,
                "The contribution cutoff must be in the half-open interval [0, 1).");
        }

        double total = 0.0;
        for (int i = 0; i < exposures.Count; i++)
        {
            (string signature, double exposure) = exposures[i];
            ArgumentNullException.ThrowIfNull(signature);
            if (double.IsNaN(exposure) || exposure < 0.0)
            {
                throw new ArgumentException("Exposures must be non-negative.", nameof(exposures));
            }

            total += exposure;
        }

        var empty = new MutationalProcessClassification(
            Array.Empty<ProcessActivity>(), MutationalProcess.Unknown);

        // Σ exposure = 0 ⇒ normalization undefined ⇒ no active processes (INV-5).
        if (total <= 0.0)
        {
            return empty;
        }

        // Per-signature normalized contribution, cutoff, then aggregate surviving contributions by process.
        var byProcess = new Dictionary<MutationalProcess, double>();
        foreach ((string signature, double exposure) in exposures)
        {
            double contribution = exposure / total;
            if (contribution < contributionCutoff)
            {
                continue; // deconstructSigs: weights[weights < signature.cutoff] <- 0
            }

            MutationalProcess process = GetMutationalProcess(signature);
            if (process == MutationalProcess.Unknown)
            {
                continue; // unmapped-aetiology signatures contribute to no recognized process (COSMIC)
            }

            byProcess.TryGetValue(process, out double accumulated);
            byProcess[process] = accumulated + contribution;
        }

        if (byProcess.Count == 0)
        {
            return empty;
        }

        // Deterministic order: descending contribution, then by process enum for ties.
        var active = byProcess
            .Select(kv => new ProcessActivity(kv.Key, kv.Value))
            .OrderByDescending(p => p.Contribution)
            .ThenBy(p => p.Process)
            .ToArray();

        return new MutationalProcessClassification(active, active[0].Process);
    }

    #endregion

    #region Fusion Gene Detection (ONCO-FUSION-001)

    /// <summary>
    /// Default minimum number of junction-spanning (split) reads required to support a fusion.
    /// Source: STAR-Fusion source (Haas et al. 2017), default <c>my $MIN_JUNCTION_READS = 1;</c>
    /// https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinJunctionReads = 1;

    /// <summary>
    /// Default minimum total fusion support = (junction reads + spanning/discordant fragments), applied
    /// when at least one junction read is present.
    /// Source: STAR-Fusion source, default <c>my $MIN_SUM_FRAGS = 2;</c> ("requires at least one junction
    /// read"). https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinSumFrags = 2;

    /// <summary>
    /// Default minimum number of spanning (discordant) fragments required when there are NO junction reads.
    /// Source: STAR-Fusion source, default <c>my $MIN_SPANNING_FRAGS_ONLY = 5;</c>
    /// https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinSpanningFragsOnly = 5;

    /// <summary>Number of nucleotides per codon (reading frame is read in triplets).
    /// Source: Wikipedia "Reading frame" (Badger &amp; Olsen 1999; Lodish 6th ed.):
    /// a reading frame reads nucleotides "as a sequence of triplets". https://en.wikipedia.org/wiki/Reading_frame</summary>
    private const int CodonLength = 3;

    /// <summary>
    /// A candidate gene-fusion breakpoint with its per-class supporting-read counts.
    /// Read classes follow the Arriba output schema (Uhrig et al. 2021): split reads anchored in each
    /// partner and discordant (spanning) mate-pairs.
    /// </summary>
    /// <param name="Gene5Prime">5' (upstream) fusion partner gene symbol.</param>
    /// <param name="Gene3Prime">3' (downstream) fusion partner gene symbol.</param>
    /// <param name="SplitReads5Prime">Split (junction) reads whose longer segment (anchor) maps to the 5' gene (Arriba split_reads1).</param>
    /// <param name="SplitReads3Prime">Split (junction) reads whose longer segment (anchor) maps to the 3' gene (Arriba split_reads2).</param>
    /// <param name="DiscordantMates">Discordant (spanning/bridge) mate-pairs supporting the fusion (Arriba discordant_mates).</param>
    /// <param name="FivePrimeCodingBases">Coding bases the 5' partner contributes upstream of the breakpoint (for the reading-frame test); -1 if unknown.</param>
    /// <param name="ThreePrimeStartPhase">Coding-start phase (0, 1, or 2) of the 3' partner at the breakpoint; -1 if unknown.</param>
    public readonly record struct FusionCandidate(
        string Gene5Prime,
        string Gene3Prime,
        int SplitReads5Prime,
        int SplitReads3Prime,
        int DiscordantMates,
        int FivePrimeCodingBases = -1,
        int ThreePrimeStartPhase = -1);

    /// <summary>
    /// Minimum-support thresholds for fusion calling. Defaults mirror STAR-Fusion (Haas et al.).
    /// </summary>
    /// <param name="MinJunctionReads">Minimum junction (split) reads required (STAR-Fusion min_junction_reads).</param>
    /// <param name="MinSumFrags">Minimum total support when ≥1 junction read present (STAR-Fusion min_sum_frags).</param>
    /// <param name="MinSpanningFragsOnly">Minimum discordant fragments required when there are no junction reads (STAR-Fusion min_spanning_frags_only).</param>
    public readonly record struct FusionDetectionThresholds(
        int MinJunctionReads = DefaultMinJunctionReads,
        int MinSumFrags = DefaultMinSumFrags,
        int MinSpanningFragsOnly = DefaultMinSpanningFragsOnly)
    {
        // A record struct's *implicit* parameterless constructor (e.g. default(T) / new()) bypasses the
        // positional defaults above and zero-initialises every field, which would silently disable the
        // STAR-Fusion thresholds. Declaring an explicit parameterless constructor restores the defaults.
        public FusionDetectionThresholds()
            : this(DefaultMinJunctionReads, DefaultMinSumFrags, DefaultMinSpanningFragsOnly)
        {
        }
    }

    /// <summary>Reading-frame status of the 3' partner across the fusion junction.</summary>
    public enum FusionReadingFrame
    {
        /// <summary>Reading frame is preserved across the junction (3' partner stays in phase).</summary>
        InFrame,

        /// <summary>Reading frame is shifted across the junction (3' partner translated out of phase).</summary>
        OutOfFrame,

        /// <summary>Coding-phase information was not supplied, so frame could not be determined.</summary>
        Unknown
    }

    /// <summary>A detected (passing) gene fusion with its supporting evidence.</summary>
    /// <param name="Gene5Prime">5' fusion partner.</param>
    /// <param name="Gene3Prime">3' fusion partner.</param>
    /// <param name="JunctionReads">Total junction (split) reads = SplitReads5Prime + SplitReads3Prime.</param>
    /// <param name="DiscordantMates">Discordant (spanning) mate-pairs.</param>
    /// <param name="TotalSupport">Total supporting reads = SplitReads5Prime + SplitReads3Prime + DiscordantMates (Arriba).</param>
    /// <param name="ReadingFrame">In-frame / out-of-frame / unknown status of the junction.</param>
    public readonly record struct FusionCall(
        string Gene5Prime,
        string Gene3Prime,
        int JunctionReads,
        int DiscordantMates,
        int TotalSupport,
        FusionReadingFrame ReadingFrame);

    /// <summary>
    /// Total supporting reads for a candidate = split_reads1 + split_reads2 + discordant_mates.
    /// Source: Arriba output spec — "The total number of supporting reads can be obtained by summing up the
    /// reads given in the columns split_reads1, split_reads2, discordant_mates".
    /// https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public static int ComputeTotalSupport(FusionCandidate candidate)
        => candidate.SplitReads5Prime + candidate.SplitReads3Prime + candidate.DiscordantMates;

    /// <summary>
    /// Determines whether the 3' partner is fused in-frame using codon phase.
    /// A fusion is in-frame iff the coding bases the 5' partner contributes upstream of the breakpoint,
    /// taken modulo 3 relative to the 3' partner's coding-start phase, keep the downstream codons in phase:
    /// <c>(fivePrimeCodingBases - threePrimeStartPhase) mod 3 == 0</c>.
    /// Source: Genomics England exon-phase rule ("if one exon finishes after the second letter of a triplet
    /// (end phase 2), the next one should start with the third letter") +
    /// reading frames read in triplets / modulo 3 (Wikipedia "Reading frame", Badger &amp; Olsen 1999).
    /// https://www.genomicsengland.co.uk/blog/gene-fusion-reporting , https://en.wikipedia.org/wiki/Reading_frame
    /// </summary>
    /// <param name="fivePrimeCodingBases">Coding bases contributed by the 5' partner upstream of the breakpoint (≥ 0).</param>
    /// <param name="threePrimeStartPhase">Coding-start phase of the 3' partner at the breakpoint (0, 1, or 2).</param>
    /// <returns><see langword="true"/> if the junction preserves the reading frame.</returns>
    /// <exception cref="ArgumentOutOfRangeException">A negative base count, or a phase outside {0,1,2}.</exception>
    public static bool IsInFrame(int fivePrimeCodingBases, int threePrimeStartPhase)
    {
        if (fivePrimeCodingBases < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(fivePrimeCodingBases),
                "Coding-base count cannot be negative.");
        }

        if (threePrimeStartPhase < 0 || threePrimeStartPhase >= CodonLength)
        {
            throw new ArgumentOutOfRangeException(nameof(threePrimeStartPhase),
                "Coding-start phase must be 0, 1, or 2.");
        }

        // Modulo arithmetic over a non-negative dividend; result is in [0, CodonLength).
        return (fivePrimeCodingBases - threePrimeStartPhase) % CodonLength == 0;
    }

    /// <summary>
    /// Detects candidate gene fusions from breakpoint supporting-read counts using the STAR-Fusion
    /// minimum-support rule, and reports each passing fusion with its total support and reading-frame status.
    /// </summary>
    /// <remarks>
    /// A candidate is reported as a fusion iff:
    /// <list type="bullet">
    /// <item><description>gene5p ≠ gene3p (a gene is not fused with itself — Registry invariant); and</description></item>
    /// <item><description>when junction (split) reads ≥ 1: junctionReads ≥ MinJunctionReads AND totalSupport ≥ MinSumFrags
    /// (STAR-Fusion min_junction_reads / min_sum_frags); OR</description></item>
    /// <item><description>when there are no junction reads: discordantMates ≥ MinSpanningFragsOnly
    /// (STAR-Fusion min_spanning_frags_only).</description></item>
    /// </list>
    /// Results are returned ordered by descending total support (STAR-Fusion scores fusions by the
    /// abundance of supporting reads), with the gene pair as a deterministic tie-breaker.
    /// </remarks>
    /// <param name="candidates">Candidate breakpoints with per-class supporting-read counts.</param>
    /// <param name="thresholds">Optional support thresholds; defaults to STAR-Fusion defaults.</param>
    /// <returns>Detected fusions ordered by descending total support.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="candidates"/> is null.</exception>
    /// <exception cref="ArgumentException">A candidate has a negative supporting-read count.</exception>
    public static IReadOnlyList<FusionCall> DetectFusions(
        IEnumerable<FusionCandidate> candidates,
        FusionDetectionThresholds? thresholds = null)
    {
        ArgumentNullException.ThrowIfNull(candidates);

        FusionDetectionThresholds t = thresholds ?? new FusionDetectionThresholds();
        var calls = new List<FusionCall>();

        foreach (FusionCandidate candidate in candidates)
        {
            if (candidate.SplitReads5Prime < 0 || candidate.SplitReads3Prime < 0 || candidate.DiscordantMates < 0)
            {
                throw new ArgumentException(
                    "Supporting-read counts cannot be negative.", nameof(candidates));
            }

            // INV-1: a gene is not fused with itself (case-insensitive symbol comparison).
            if (string.Equals(candidate.Gene5Prime, candidate.Gene3Prime, StringComparison.OrdinalIgnoreCase))
            {
                continue;
            }

            int junctionReads = candidate.SplitReads5Prime + candidate.SplitReads3Prime;
            int totalSupport = ComputeTotalSupport(candidate);

            bool passes = junctionReads >= t.MinJunctionReads
                ? totalSupport >= t.MinSumFrags
                // No junction reads: fall back to the spanning-fragments-only rule.
                : candidate.DiscordantMates >= t.MinSpanningFragsOnly;

            if (!passes)
            {
                continue;
            }

            FusionReadingFrame frame = ResolveReadingFrame(candidate);
            calls.Add(new FusionCall(
                candidate.Gene5Prime,
                candidate.Gene3Prime,
                junctionReads,
                candidate.DiscordantMates,
                totalSupport,
                frame));
        }

        // INV-4: order by descending total support (abundance of supporting reads), deterministic ties.
        return calls
            .OrderByDescending(c => c.TotalSupport)
            .ThenBy(c => c.Gene5Prime, StringComparer.Ordinal)
            .ThenBy(c => c.Gene3Prime, StringComparer.Ordinal)
            .ToArray();
    }

    private static FusionReadingFrame ResolveReadingFrame(FusionCandidate candidate)
    {
        if (candidate.FivePrimeCodingBases < 0 || candidate.ThreePrimeStartPhase < 0
            || candidate.ThreePrimeStartPhase >= CodonLength)
        {
            return FusionReadingFrame.Unknown;
        }

        return IsInFrame(candidate.FivePrimeCodingBases, candidate.ThreePrimeStartPhase)
            ? FusionReadingFrame.InFrame
            : FusionReadingFrame.OutOfFrame;
    }

    #endregion

    #region Known Fusion Database Lookup (ONCO-FUSION-002)

    /// <summary>
    /// HGNC fusion-designation separator: a double colon between the 5' and 3' partner symbols.
    /// Source: Bruford et al. (2021), HGNC recommendations for the designation of gene fusions —
    /// "HGNC recommends that a new separator—a double colon (::)—be used in describing gene fusions,
    /// e.g., BCR::ABL1." https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </summary>
    public const string FusionDesignationSeparator = "::";

    /// <summary>
    /// The result of looking a fusion up against a known-fusion set.
    /// </summary>
    /// <param name="Designation">The HGNC <c>5'::3'</c> designation of the queried fusion (e.g. <c>BCR::ABL1</c>).</param>
    /// <param name="IsKnown"><see langword="true"/> if the directional designation was present in the supplied set.</param>
    /// <param name="Annotation">The caller-supplied annotation for the matched designation, or <see langword="null"/> if not known.</param>
    public readonly record struct KnownFusionMatch(string Designation, bool IsKnown, string? Annotation);

    /// <summary>
    /// Formats the HGNC designation of a gene fusion as <c>gene5p::gene3p</c>.
    /// The 5' partner is always written first, before the double colon, irrespective of chromosomal
    /// location or gene orientation; the designation is therefore directional (A::B ≠ B::A).
    /// Source: Bruford et al. (2021), HGNC recommendations for the designation of gene fusions —
    /// "the 5′ partner gene should always be listed first in the description of a fusion gene, i.e.,
    /// before the double colon" and "a double colon (::) … e.g., BCR::ABL1".
    /// https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </summary>
    /// <param name="gene5p">5' (upstream) partner gene symbol; must be non-empty.</param>
    /// <param name="gene3p">3' (downstream) partner gene symbol; must be non-empty.</param>
    /// <returns>The designation string <c>gene5p + "::" + gene3p</c>.</returns>
    /// <exception cref="ArgumentException">Either symbol is null, empty, or whitespace.</exception>
    public static string GetFusionAnnotation(string gene5p, string gene3p)
    {
        if (string.IsNullOrWhiteSpace(gene5p))
        {
            throw new ArgumentException("5' partner gene symbol must be non-empty.", nameof(gene5p));
        }

        if (string.IsNullOrWhiteSpace(gene3p))
        {
            throw new ArgumentException("3' partner gene symbol must be non-empty.", nameof(gene3p));
        }

        return gene5p + FusionDesignationSeparator + gene3p;
    }

    /// <summary>
    /// Looks a detected fusion up against a caller-supplied set of known fusions, keyed by their
    /// HGNC <c>5'::3'</c> designation.
    /// </summary>
    /// <remarks>
    /// The lookup is <b>directional</b>: the key is built with the 5' partner first (per Bruford et al. 2021),
    /// so a reciprocal fusion (partners swapped) is a different designation and does NOT match.
    /// Symbol comparison is case-insensitive (ordinal-ignore-case); the known-fusion set membership and the
    /// annotation text are entirely caller-supplied — this library bundles no curated fusion database
    /// (Mitelman / COSMIC / ChimerDB content is the caller's responsibility).
    /// Source (designation format and directional keying): Bruford et al. (2021),
    /// https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </remarks>
    /// <param name="fusion">The fusion to look up (its 5'/3' partners define the key).</param>
    /// <param name="knownFusions">
    /// Caller-supplied map from <c>5'::3'</c> designation to its annotation. For case-insensitive matching,
    /// supply a dictionary built with <see cref="StringComparer.OrdinalIgnoreCase"/> (the method also probes
    /// case-insensitively when the dictionary is not already case-insensitive).
    /// </param>
    /// <returns>A <see cref="KnownFusionMatch"/> reporting the designation and, if present, its annotation.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="knownFusions"/> is null.</exception>
    /// <exception cref="ArgumentException">A fusion partner symbol is null, empty, or whitespace.</exception>
    public static KnownFusionMatch MatchKnownFusions(
        FusionCall fusion,
        IReadOnlyDictionary<string, string> knownFusions)
    {
        ArgumentNullException.ThrowIfNull(knownFusions);

        string designation = GetFusionAnnotation(fusion.Gene5Prime, fusion.Gene3Prime);

        // Directional key (5'::3'). Try the supplied dictionary's own comparer first; if it is not
        // case-insensitive, fall back to an explicit case-insensitive scan so that e.g. "eml4::alk"
        // matches a stored "EML4::ALK" (HGNC symbols are case-defined, but inputs vary in case).
        if (knownFusions.TryGetValue(designation, out string? annotation))
        {
            return new KnownFusionMatch(designation, IsKnown: true, annotation);
        }

        foreach (KeyValuePair<string, string> entry in knownFusions)
        {
            if (string.Equals(entry.Key, designation, StringComparison.OrdinalIgnoreCase))
            {
                return new KnownFusionMatch(designation, IsKnown: true, entry.Value);
            }
        }

        return new KnownFusionMatch(designation, IsKnown: false, Annotation: null);
    }

    #endregion

    #region Fusion Breakpoint Analysis (ONCO-FUSION-003)

    /// <summary>
    /// Location category of a fusion breakpoint within a partner transcript, mirroring the Arriba
    /// <c>site1</c>/<c>site2</c> output column. Source: Arriba output spec — "Possible values are: 5'UTR,
    /// 3'UTR, UTR, CDS, exon, intron, and intergenic". https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public enum BreakpointSite
    {
        /// <summary>Coding sequence (the only site at which two reading frames are joined).</summary>
        Cds,

        /// <summary>5' untranslated region.</summary>
        FivePrimeUtr,

        /// <summary>3' untranslated region.</summary>
        ThreePrimeUtr,

        /// <summary>Untranslated region (strand/UTR side not resolved).</summary>
        Utr,

        /// <summary>Exon (non-coding part).</summary>
        Exon,

        /// <summary>Intron.</summary>
        Intron,

        /// <summary>Intergenic region.</summary>
        Intergenic
    }

    /// <summary>
    /// Reading-frame consequence reported by <see cref="AnalyzeBreakpoint"/>, mirroring the Arriba
    /// <c>reading_frame</c> column. Source: Arriba output spec — reading_frame ∈ {in-frame, out-of-frame,
    /// stop-codon, .}; the dot ("not predicted") is used when the peptide cannot be predicted because a
    /// breakpoint is not in coding context. https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public enum BreakpointFrameStatus
    {
        /// <summary>Both breakpoints are in CDS and the junction preserves the reading frame (in-frame).</summary>
        InFrame,

        /// <summary>Both breakpoints are in CDS but the junction shifts the reading frame (out-of-frame).</summary>
        OutOfFrame,

        /// <summary>The chimeric ORF reaches a stop codon at/after the junction (Arriba <c>stop-codon</c>).</summary>
        StopCodon,

        /// <summary>Frame cannot be called because a breakpoint is not in coding context (Arriba <c>.</c>).</summary>
        NotPredicted
    }

    /// <summary>
    /// A fusion breakpoint described by its per-partner site categories and the coding-frame quantities
    /// needed to call the junction reading frame. The 5'/3' partner symbols carry over from a
    /// <see cref="FusionCall"/>; <see cref="FivePrimeCodingBases"/> and <see cref="ThreePrimeStartPhase"/>
    /// are the same coding-frame inputs used by <see cref="IsInFrame"/> (ONCO-FUSION-001). Read class
    /// schema follows Arriba (Uhrig et al. 2021). https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    /// <param name="Gene5Prime">5' (upstream) fusion partner gene symbol.</param>
    /// <param name="Gene3Prime">3' (downstream) fusion partner gene symbol.</param>
    /// <param name="Site5Prime">Site category of the breakpoint in the 5' partner.</param>
    /// <param name="Site3Prime">Site category of the breakpoint in the 3' partner.</param>
    /// <param name="FivePrimeCodingBases">Coding bases the 5' partner contributes upstream of the breakpoint (≥ 0).</param>
    /// <param name="ThreePrimeStartPhase">Coding-start phase (0, 1, or 2) of the 3' partner at the breakpoint.</param>
    public readonly record struct FusionBreakpoint(
        string Gene5Prime,
        string Gene3Prime,
        BreakpointSite Site5Prime,
        BreakpointSite Site3Prime,
        int FivePrimeCodingBases,
        int ThreePrimeStartPhase);

    /// <summary>The breakpoint-analysis result for a fusion junction.</summary>
    /// <param name="Gene5Prime">5' partner (carried through unchanged).</param>
    /// <param name="Gene3Prime">3' partner (carried through unchanged).</param>
    /// <param name="Site5Prime">Site category of the 5' breakpoint.</param>
    /// <param name="Site3Prime">Site category of the 3' breakpoint.</param>
    /// <param name="BreakpointInCoding">True iff both breakpoints lie in CDS (a coding-to-coding junction).</param>
    /// <param name="FrameStatus">Reading-frame consequence of the junction.</param>
    public readonly record struct BreakpointAnalysis(
        string Gene5Prime,
        string Gene3Prime,
        BreakpointSite Site5Prime,
        BreakpointSite Site3Prime,
        bool BreakpointInCoding,
        BreakpointFrameStatus FrameStatus);

    /// <summary>The predicted protein product of a fusion junction.</summary>
    /// <param name="ChimericCds">The chimeric coding sequence = 5' CDS prefix ++ 3' CDS suffix (uppercased DNA).</param>
    /// <param name="Peptide">The translated fusion peptide, truncated at the first stop codon.</param>
    /// <param name="Effect">Reading-frame effect of the junction (in-frame / out-of-frame).</param>
    /// <param name="HasPrematureStop">True iff a stop codon was reached before the end of the chimeric ORF.</param>
    public readonly record struct FusionProteinPrediction(
        string ChimericCds,
        string Peptide,
        BreakpointFrameStatus Effect,
        bool HasPrematureStop);

    /// <summary>
    /// Analyzes a fusion breakpoint: reports the per-partner site categories and calls the junction
    /// reading frame. A reading-frame call (<see cref="BreakpointFrameStatus.InFrame"/> /
    /// <see cref="BreakpointFrameStatus.OutOfFrame"/>) is made ONLY when both breakpoints lie in CDS; if
    /// either breakpoint is in a UTR, intron, exon (non-coding part), or intergenic region the frame is
    /// <see cref="BreakpointFrameStatus.NotPredicted"/> (Arriba <c>reading_frame = .</c>). The in-frame
    /// test reuses the codon-phase rule of <see cref="IsInFrame"/>:
    /// <c>(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0</c>.
    /// Source: Arriba output spec (site / reading_frame), AGFusion frame rule (Murphy &amp; Elemento 2016).
    /// https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    /// <param name="fusion">The breakpoint to analyze.</param>
    /// <returns>The site categories and reading-frame consequence of the junction.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Both breakpoints are CDS but <see cref="FusionBreakpoint.FivePrimeCodingBases"/> is negative or
    /// <see cref="FusionBreakpoint.ThreePrimeStartPhase"/> is outside {0, 1, 2}.
    /// </exception>
    public static BreakpointAnalysis AnalyzeBreakpoint(FusionBreakpoint fusion)
    {
        bool bothCoding = fusion.Site5Prime == BreakpointSite.Cds
                       && fusion.Site3Prime == BreakpointSite.Cds;

        // A frame call is only defined for a coding-to-coding junction (Arriba reading_frame = '.' otherwise).
        BreakpointFrameStatus frame = bothCoding
            ? (IsInFrame(fusion.FivePrimeCodingBases, fusion.ThreePrimeStartPhase)
                ? BreakpointFrameStatus.InFrame
                : BreakpointFrameStatus.OutOfFrame)
            : BreakpointFrameStatus.NotPredicted;

        return new BreakpointAnalysis(
            fusion.Gene5Prime,
            fusion.Gene3Prime,
            fusion.Site5Prime,
            fusion.Site3Prime,
            bothCoding,
            frame);
    }

    /// <summary>
    /// Predicts the protein product of a fusion from the two partners' coding sequences. The chimeric CDS
    /// is the 5' partner's CDS taken up to the breakpoint (a prefix) concatenated with the 3' partner's CDS
    /// taken from the breakpoint onward (a suffix); it is then translated with the standard genetic code
    /// (NCBI Table 1) and truncated at the first stop codon.
    /// Source (verbatim from AGFusion model.py, Murphy &amp; Elemento 2016):
    /// <c>cds_5prime = transcript1.coding_sequence[0:junction5]</c>,
    /// <c>cds_3prime = transcript2.coding_sequence[junction3:]</c>,
    /// <c>seq = cds_5prime + cds_3prime</c>,
    /// <c>protein_seq = cds.seq.translate(); protein_seq = protein_seq[0:protein_seq.find("*")]</c>.
    /// When the junction is out-of-frame the chimeric CDS is first trimmed to a whole number of codons
    /// (<c>cds[0:3*(len//3)]</c>) so the 3' partner is read in its (shifted) frame.
    /// https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py
    /// </summary>
    /// <param name="fusion">The breakpoint (its site categories and codon-frame quantities drive the effect call).</param>
    /// <param name="transcripts">
    /// The two partner CDS sequences as (fivePrimeCds, threePrimeCds): the 5' partner's full coding sequence
    /// and the 3' partner's full coding sequence (DNA, A/C/G/T). The breakpoint offsets are taken from
    /// <paramref name="fusion"/>: the 5' prefix length is <see cref="FusionBreakpoint.FivePrimeCodingBases"/>
    /// and the 3' suffix starts at <see cref="FusionBreakpoint.ThreePrimeStartPhase"/>.
    /// </param>
    /// <returns>The chimeric CDS, the translated (first-stop-truncated) peptide, the frame effect, and a premature-stop flag.</returns>
    /// <exception cref="ArgumentNullException">A CDS sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// An offset is out of range for its CDS, or the 3' start phase is outside {0, 1, 2}.
    /// </exception>
    public static FusionProteinPrediction PredictFusionProtein(
        FusionBreakpoint fusion,
        (string FivePrimeCds, string ThreePrimeCds) transcripts)
    {
        ArgumentNullException.ThrowIfNull(transcripts.FivePrimeCds);
        ArgumentNullException.ThrowIfNull(transcripts.ThreePrimeCds);

        string fivePrimeCds = transcripts.FivePrimeCds.ToUpperInvariant();
        string threePrimeCds = transcripts.ThreePrimeCds.ToUpperInvariant();

        int junction5 = fusion.FivePrimeCodingBases;
        int junction3 = fusion.ThreePrimeStartPhase;

        if (junction5 < 0 || junction5 > fivePrimeCds.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(fusion),
                "FivePrimeCodingBases (5' prefix length) is out of range for the 5' CDS.");
        }

        if (junction3 < 0 || junction3 > threePrimeCds.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(fusion),
                "ThreePrimeStartPhase (3' suffix start) is out of range for the 3' CDS.");
        }

        // AGFusion: cds_5prime = transcript1.coding_sequence[0:junction5]; cds_3prime = transcript2[junction3:].
        string cds5 = fivePrimeCds.Substring(0, junction5);
        string cds3 = threePrimeCds.Substring(junction3);
        string chimericCds = cds5 + cds3;

        // Frame effect by the AGFusion / IsInFrame codon-phase rule (only meaningful junction3 in {0,1,2}).
        bool inFrame = junction3 < CodonLength
            && IsInFrame(junction5, junction3);
        BreakpointFrameStatus effect = inFrame
            ? BreakpointFrameStatus.InFrame
            : BreakpointFrameStatus.OutOfFrame;

        // AGFusion: an out-of-frame CDS is trimmed to whole codons before translation; an in-frame CDS is
        // translated as-is (a trailing partial codon, if any, is not translatable and is dropped).
        int translatableLength = chimericCds.Length - (chimericCds.Length % CodonLength);
        var peptide = new System.Text.StringBuilder(translatableLength / CodonLength);
        bool hasStop = false;

        for (int i = 0; i < translatableLength; i += CodonLength)
        {
            string codon = chimericCds.Substring(i, CodonLength);
            char aminoAcid = GeneticCode.Standard.Translate(codon);
            if (aminoAcid == StopCodonSymbol)
            {
                // AGFusion truncates the peptide at the first stop codon (protein_seq[0:find("*")]).
                hasStop = true;
                break;
            }

            peptide.Append(aminoAcid);
        }

        return new FusionProteinPrediction(chimericCds, peptide.ToString(), effect, hasStop);
    }

    /// <summary>
    /// Stop-codon marker returned by <see cref="GeneticCode.Translate"/> for UAA/UAG/UGA (NCBI Table 1).
    /// Source: AGFusion truncation at the first <c>'*'</c>; <see cref="GeneticCode"/> emits '*' for stops.
    /// </summary>
    private const char StopCodonSymbol = '*';

    #endregion

    #region Copy-Number Alteration Classification (ONCO-CNA-001)

    /// <summary>
    /// Reference (germline) ploidy used as the log2 anchor: an autosomal diploid genome has copy number 2.
    /// Source: CNVkit <c>cnvlib/call.py</c> — <c>_log2_ratio_to_absolute_pure</c> uses
    /// <c>ncopies = ref_copies * 2**log2_ratio</c> with <c>ref_copies = ploidy = 2</c> for autosomes.
    /// </summary>
    public const double DiploidReferencePloidy = 2.0;

    /// <summary>
    /// Default CNVkit hard-threshold cutoffs for calling integer copy number from a log2 ratio, in
    /// ascending order. The four cutoffs partition the log2 axis into the five copy-number states
    /// 0 / 1 / 2 / 3 / 4+. Source: CNVkit <c>cnvlib/call.py</c> <c>do_call</c> default
    /// <c>thresholds = (-1.1, -0.25, 0.2, 0.7)</c>; the <c>absolute_threshold</c> docstring states the
    /// cutoffs verbatim as DEL(0) &lt; −1.1, LOSS(1) &lt; −0.25, GAIN(3) ≥ +0.2, AMP(4) ≥ +0.7
    /// (tumor-sample heuristic, safe for purity ≥ 30%).
    /// </summary>
    public static readonly IReadOnlyList<double> DefaultCopyNumberThresholds =
        new[] { -1.1, -0.25, 0.2, 0.7 };

    /// <summary>
    /// Number of hard-threshold cutoffs required to define the five copy-number states. Four cutoffs
    /// partition the log2 axis into states 0/1/2/3/4+. Source: CNVkit <c>absolute_threshold</c>.
    /// </summary>
    private const int CopyNumberThresholdCount = 4;

    /// <summary>
    /// Integer copy number that marks the start of the amplification class (CN ≥ 4). Source: CNVkit
    /// <c>absolute_threshold</c> docstring — "AMP(4) ≥ +0.7"; values above the last threshold are called
    /// <c>ceil(2·2^log2)</c>, which is ≥ 4.
    /// </summary>
    private const int AmplificationCopyNumber = 4;

    /// <summary>
    /// A discrete copy-number alteration (CNA) state assigned to a genomic region from its log2 copy ratio.
    /// The five states correspond to CNVkit integer copy-number calls 0 / 1 / 2 / 3 / ≥4 for a diploid
    /// reference. Source: CNVkit <c>cnvlib/call.py</c> <c>absolute_threshold</c>; GISTIC2.0 amplitude
    /// semantics (Mermel et al. 2011).
    /// </summary>
    public enum CopyNumberState
    {
        /// <summary>Deep (homozygous) deletion: integer copy number 0 (log2 ≤ −1.1).</summary>
        DeepDeletion,

        /// <summary>Single-copy loss: integer copy number 1 (−1.1 &lt; log2 ≤ −0.25).</summary>
        Loss,

        /// <summary>Copy-number neutral (diploid): integer copy number 2 (−0.25 &lt; log2 ≤ 0.2).</summary>
        Neutral,

        /// <summary>Single-copy gain: integer copy number 3 (0.2 &lt; log2 ≤ 0.7).</summary>
        Gain,

        /// <summary>Amplification: integer copy number ≥ 4 (log2 &gt; 0.7).</summary>
        Amplification
    }

    /// <summary>
    /// A copy-number call for one region: the input log2 ratio, the continuous and integer absolute copy
    /// numbers, and the discrete CNA state.
    /// </summary>
    /// <param name="Log2Ratio">Input log2 copy ratio log2(tumor_depth / normal_depth).</param>
    /// <param name="AbsoluteCopyNumber">Continuous absolute copy number n = ploidy·2^log2.</param>
    /// <param name="IntegerCopyNumber">Hard-threshold integer copy number (CNVkit <c>absolute_threshold</c>).</param>
    /// <param name="State">Discrete CNA classification.</param>
    public readonly record struct CopyNumberCall(
        double Log2Ratio,
        double AbsoluteCopyNumber,
        int IntegerCopyNumber,
        CopyNumberState State);

    /// <summary>
    /// Converts a log2 copy ratio to a continuous absolute copy number for a pure sample:
    /// <c>n = ploidy · 2^log2</c>. For an autosomal diploid reference (ploidy = 2) this is
    /// <c>n = 2 · 2^log2</c>, so log2 = 0 ⇒ 2 copies, log2 = 1 ⇒ 4 copies, log2 = −1 ⇒ 1 copy.
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>_log2_ratio_to_absolute_pure</c>:
    /// <c>ncopies = ref_copies * 2**log2_ratio</c>.
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio (may be any finite value; NaN propagates to NaN).</param>
    /// <param name="ploidy">Reference (germline) ploidy; 2 for an autosomal diploid genome.</param>
    /// <returns>Continuous absolute copy number n = ploidy·2^log2 (≥ 0 for finite input).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static double Log2RatioToCopyNumber(double log2Ratio, double ploidy = DiploidReferencePloidy)
    {
        if (double.IsNaN(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy must be positive.");
        }

        return ploidy * Math.Pow(2.0, log2Ratio);
    }

    /// <summary>
    /// Calls an integer copy number from a log2 ratio using CNVkit's hard-threshold method. The copy number
    /// is the index of the first ascending threshold the log2 value is less than or equal to (counting up
    /// from 0); if the log2 value exceeds every threshold, the copy number is <c>ceil(ploidy · 2^log2)</c>.
    /// A NaN log2 ratio is a no-call and returns the neutral reference copy number (rounded ploidy).
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>absolute_threshold</c> — "Integer values are assigned for
    /// log2 ratio values less than each given threshold value in sequence, counting up from zero. Above the
    /// last threshold value, integer copy numbers are called assuming full purity, diploidy, and rounding up."
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio; NaN is a no-call (neutral).</param>
    /// <param name="thresholds">
    /// Exactly four strictly ascending cutoffs partitioning the log2 axis into states 0/1/2/3/4+; when null,
    /// <see cref="DefaultCopyNumberThresholds"/> (−1.1, −0.25, 0.2, 0.7) is used.
    /// </param>
    /// <param name="ploidy">Reference ploidy used both for the neutral no-call and the amplification ceiling.</param>
    /// <returns>The integer copy number (≥ 0).</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static int CallCopyNumber(
        double log2Ratio,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        var cutoffs = ValidateThresholds(thresholds);
        if (double.IsNaN(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy must be positive.");
        }

        if (double.IsNaN(log2Ratio))
        {
            // No-call: CNVkit replaces a NaN log2 with the neutral reference copy number.
            return (int)Math.Round(ploidy, MidpointRounding.AwayFromZero);
        }

        // CN = index of the first cutoff the log2 value is <= (inclusive boundary), counting from 0.
        for (int cn = 0; cn < cutoffs.Count; cn++)
        {
            if (log2Ratio <= cutoffs[cn])
            {
                return cn;
            }
        }

        // Above the last cutoff: round up the absolute copy number (CNVkit ceil), yielding CN ≥ 4.
        return (int)Math.Ceiling(Log2RatioToCopyNumber(log2Ratio, ploidy));
    }

    /// <summary>
    /// Classifies a single region's log2 copy ratio into a <see cref="CopyNumberCall"/> carrying the
    /// continuous absolute copy number, the hard-threshold integer copy number, and the discrete
    /// <see cref="CopyNumberState"/>. The state is derived from the integer copy number: 0 → DeepDeletion,
    /// 1 → Loss, 2 → Neutral, 3 → Gain, ≥4 → Amplification. Source: CNVkit <c>absolute_threshold</c>
    /// (DEL(0)/LOSS(1)/neutral(2)/GAIN(3)/AMP(4)); GISTIC2.0 amplitude semantics (Mermel et al. 2011).
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio; NaN is a no-call (Neutral, CN = rounded ploidy).</param>
    /// <param name="thresholds">Four ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy (default diploid).</param>
    /// <returns>The copy-number call with absolute CN, integer CN, and CNA state.</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static CopyNumberCall ClassifyCopyNumber(
        double log2Ratio,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        int integerCopyNumber = CallCopyNumber(log2Ratio, thresholds, ploidy);
        double absolute = double.IsNaN(log2Ratio) ? ploidy : Log2RatioToCopyNumber(log2Ratio, ploidy);
        CopyNumberState state = StateFromCopyNumber(integerCopyNumber);

        return new CopyNumberCall(log2Ratio, absolute, integerCopyNumber, state);
    }

    /// <summary>
    /// Classifies a sequence of per-region log2 copy ratios, returning one <see cref="CopyNumberCall"/> per
    /// input value in input order (length and order preserving). Thin per-element wrapper over
    /// <see cref="ClassifyCopyNumber(double, IReadOnlyList{double}?, double)"/>.
    /// </summary>
    /// <param name="log2Ratios">Per-region log2 copy ratios.</param>
    /// <param name="thresholds">Four ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy (default diploid).</param>
    /// <returns>One call per input log2 ratio, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="log2Ratios"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static IReadOnlyList<CopyNumberCall> ClassifyCopyNumbers(
        IEnumerable<double> log2Ratios,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        ArgumentNullException.ThrowIfNull(log2Ratios);
        var cutoffs = ValidateThresholds(thresholds);

        var calls = new List<CopyNumberCall>();
        foreach (double log2Ratio in log2Ratios)
        {
            calls.Add(ClassifyCopyNumber(log2Ratio, cutoffs, ploidy));
        }

        return calls;
    }

    /// <summary>Maps an integer copy number to its CNA state per CNVkit (0/1/2/3/≥4).</summary>
    private static CopyNumberState StateFromCopyNumber(int copyNumber)
    {
        // CN ≥ 4 is the amplification class (CNVkit AMP(4) ≥ +0.7).
        if (copyNumber >= AmplificationCopyNumber)
        {
            return CopyNumberState.Amplification;
        }

        return copyNumber switch
        {
            0 => CopyNumberState.DeepDeletion,
            1 => CopyNumberState.Loss,
            2 => CopyNumberState.Neutral,
            _ => CopyNumberState.Gain // copyNumber == 3
        };
    }

    /// <summary>
    /// Validates and returns the threshold cutoffs: exactly four strictly ascending values, or the default
    /// when null. Four cutoffs are required to define the five copy-number states (CNVkit
    /// <c>absolute_threshold</c>); a non-ascending list would not partition the log2 axis.
    /// </summary>
    private static IReadOnlyList<double> ValidateThresholds(IReadOnlyList<double>? thresholds)
    {
        if (thresholds is null)
        {
            return DefaultCopyNumberThresholds;
        }

        if (thresholds.Count != CopyNumberThresholdCount)
        {
            throw new ArgumentException(
                $"Exactly {CopyNumberThresholdCount} thresholds are required to define the five copy-number " +
                $"states (got {thresholds.Count}).",
                nameof(thresholds));
        }

        for (int i = 0; i < thresholds.Count; i++)
        {
            if (double.IsNaN(thresholds[i]))
            {
                throw new ArgumentException("Thresholds must not contain NaN.", nameof(thresholds));
            }

            if (i > 0 && thresholds[i] <= thresholds[i - 1])
            {
                throw new ArgumentException(
                    "Thresholds must be in strictly ascending order.", nameof(thresholds));
            }
        }

        return thresholds;
    }

    #endregion

    #region Focal Amplification Detection (ONCO-CNA-002)

    /// <summary>
    /// Default fraction-of-chromosome-arm cutoff separating focal from broad (arm-level) copy-number
    /// events. A segment whose length is strictly less than this fraction of its chromosome arm is focal;
    /// a segment occupying this fraction or more of the arm is arm-level. Source: Mermel et al. (2011)
    /// GISTIC2.0 — focal SCNAs have "length &lt; 98% of a chromosome arm"; events "occupying more than 98%
    /// of a chromosome arm" are arm-level. GISTIC2 parameter <c>broad_len_cutoff</c> default 0.98.
    /// </summary>
    public const double DefaultBroadLengthCutoff = 0.98;

    /// <summary>
    /// Default log2-ratio amplitude above which a copy-number gain is called an amplification. Source:
    /// GISTIC2 parameter <c>t_amp</c> default 0.1 — "Regions with a copy number gain above this positive
    /// value are considered amplified." A single-copy gain is log2(3/2) = 0.585 (CNVkit), well above 0.1.
    /// </summary>
    public const double DefaultAmplificationLog2Threshold = 0.1;

    /// <summary>
    /// Thresholds controlling focal-amplification detection: the amplitude cutoff (GISTIC2 <c>t_amp</c>)
    /// and the focal/broad length cutoff as a fraction of chromosome arm (GISTIC2 <c>broad_len_cutoff</c>).
    /// </summary>
    /// <param name="AmplificationLog2Threshold">log2 gain must strictly exceed this to be amplified (GISTIC2 <c>t_amp</c>, default 0.1).</param>
    /// <param name="BroadLengthCutoff">segment length ÷ arm length must be strictly below this to be focal (GISTIC2 <c>broad_len_cutoff</c>, default 0.98).</param>
    public readonly record struct FocalAmplificationThresholds(
        double AmplificationLog2Threshold,
        double BroadLengthCutoff)
    {
        /// <summary>GISTIC2 default thresholds: <c>t_amp</c> = 0.1, <c>broad_len_cutoff</c> = 0.98.</summary>
        public static FocalAmplificationThresholds Default { get; } =
            new(DefaultAmplificationLog2Threshold, DefaultBroadLengthCutoff);
    }

    /// <summary>
    /// A segmented copy-number region with the chromosome-arm context needed to apply the GISTIC2 length
    /// rule. The arm label (chromosome + arm letter, e.g. "17q") is matched against oncogene locations;
    /// the arm length lets the algorithm compute the segment-length / arm-length fraction.
    /// </summary>
    /// <param name="Arm">Chromosome-arm label, chromosome number followed by p/q (e.g. "17q", "8q", "7p").</param>
    /// <param name="Start">Segment start coordinate (bp); must satisfy <see cref="End"/> &gt; <see cref="Start"/>.</param>
    /// <param name="End">Segment end coordinate (bp).</param>
    /// <param name="ArmLength">Total length of the chromosome arm in bp; must be positive.</param>
    /// <param name="Log2Ratio">Segment mean log2 copy ratio.</param>
    public readonly record struct CopyNumberArmSegment(
        string Arm,
        long Start,
        long End,
        long ArmLength,
        double Log2Ratio)
    {
        /// <summary>Segment length in base pairs (End − Start).</summary>
        public long Length => End - Start;

        /// <summary>Segment length as a fraction of the chromosome arm (Length ÷ ArmLength).</summary>
        public double ArmFraction => (double)Length / ArmLength;
    }

    /// <summary>
    /// Tests whether a segment is a focal amplification: it is amplified (log2 strictly above the amplitude
    /// threshold) AND focal (length strictly below the broad-length cutoff fraction of its arm). Source:
    /// Mermel et al. (2011) length rule + GISTIC2 <c>t_amp</c>/<c>broad_len_cutoff</c>.
    /// </summary>
    /// <param name="segment">The arm-anchored copy-number segment.</param>
    /// <param name="thresholds">Amplitude and length cutoffs.</param>
    /// <returns><c>true</c> when the segment is an amplified, focal-length event.</returns>
    /// <exception cref="ArgumentException"><paramref name="segment"/> has non-positive arm length or End ≤ Start.</exception>
    public static bool IsFocalAmplification(
        in CopyNumberArmSegment segment,
        FocalAmplificationThresholds thresholds)
    {
        ValidateArmSegment(segment);

        bool amplified = segment.Log2Ratio > thresholds.AmplificationLog2Threshold;
        bool focal = segment.ArmFraction < thresholds.BroadLengthCutoff;
        return amplified && focal;
    }

    /// <summary>
    /// Detects focal amplifications among arm-anchored copy-number segments. A segment is reported when it
    /// is amplified (log2 &gt; <c>t_amp</c>) and focal (length &lt; <c>broad_len_cutoff</c> × arm length).
    /// The result is a subset of the input in input order (length- and order-preserving filter). Source:
    /// Mermel et al. (2011) GISTIC2.0 length-based focal/arm-level split; GISTIC2 <c>t_amp</c>/<c>broad_len_cutoff</c>.
    /// </summary>
    /// <param name="segments">Arm-anchored copy-number segments. Must not be null.</param>
    /// <param name="thresholds">Amplitude and length cutoffs; null uses <see cref="FocalAmplificationThresholds.Default"/> (GISTIC2 defaults).</param>
    /// <returns>The focal amplifications, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive arm length or End ≤ Start.</exception>
    public static IReadOnlyList<CopyNumberArmSegment> DetectFocalAmplifications(
        IEnumerable<CopyNumberArmSegment> segments,
        FocalAmplificationThresholds? thresholds = null)
    {
        ArgumentNullException.ThrowIfNull(segments);
        FocalAmplificationThresholds cutoffs = thresholds ?? FocalAmplificationThresholds.Default;

        var result = new List<CopyNumberArmSegment>();
        foreach (CopyNumberArmSegment segment in segments)
        {
            if (IsFocalAmplification(segment, cutoffs))
            {
                result.Add(segment);
            }
        }

        return result;
    }

    /// <summary>
    /// Maps focal-amplification segments to the recurrently amplified oncogenes resident on their
    /// chromosome arms. Each oncogene is reported once if any focal amplification falls on its arm. The
    /// panel and arms are: ERBB2 (17q), MYC (8q), EGFR (7p), CCND1 (11q), MDM2 (12q), CDK4 (12q). Source:
    /// NCBI Gene cytogenetic locations — ERBB2 17q12, MYC 8q24.21, EGFR 7p11.2, CCND1 11q13.3, MDM2 12q15,
    /// CDK4 12q14.1.
    /// </summary>
    /// <param name="amplifications">Focal amplifications (typically the output of <see cref="DetectFocalAmplifications"/>).</param>
    /// <returns>Distinct oncogene symbols whose arm carries a focal amplification, in panel order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="amplifications"/> is null.</exception>
    public static IReadOnlyList<string> IdentifyAmplifiedOncogenes(
        IEnumerable<CopyNumberArmSegment> amplifications)
    {
        ArgumentNullException.ThrowIfNull(amplifications);

        var amplifiedArms = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (CopyNumberArmSegment segment in amplifications)
        {
            if (!string.IsNullOrEmpty(segment.Arm))
            {
                amplifiedArms.Add(segment.Arm);
            }
        }

        var genes = new List<string>();
        foreach ((string gene, string arm) in OncogeneArms)
        {
            if (amplifiedArms.Contains(arm))
            {
                genes.Add(gene);
            }
        }

        return genes;
    }

    /// <summary>
    /// Recurrently amplified oncogenes and their chromosome arms (chromosome + arm letter), from NCBI Gene
    /// cytogenetic locations. Order is the registry panel order. Source: NCBI Gene — ERBB2 17q12 (Gene ID
    /// 2064), MYC 8q24.21 (4609), EGFR 7p11.2 (1956), CCND1 11q13.3 (595), MDM2 12q15 (4193), CDK4 12q14.1 (1019).
    /// </summary>
    private static readonly IReadOnlyList<(string Gene, string Arm)> OncogeneArms = new[]
    {
        ("ERBB2", "17q"),
        ("MYC", "8q"),
        ("EGFR", "7p"),
        ("CCND1", "11q"),
        ("MDM2", "12q"),
        ("CDK4", "12q"),
    };

    /// <summary>Validates an arm segment: positive arm length and End &gt; Start.</summary>
    private static void ValidateArmSegment(in CopyNumberArmSegment segment)
    {
        if (segment.ArmLength <= 0)
        {
            throw new ArgumentException(
                $"Segment on '{segment.Arm}' must have a positive arm length (got {segment.ArmLength}).",
                nameof(segment));
        }

        if (segment.End <= segment.Start)
        {
            throw new ArgumentException(
                $"Segment on '{segment.Arm}' must have End > Start (got Start={segment.Start}, End={segment.End}).",
                nameof(segment));
        }
    }

    #endregion

    #region Homozygous Deletion Detection (ONCO-CNA-003)

    /// <summary>
    /// Integer copy number of a homozygous (deep) deletion: a region with zero copies of both alleles, i.e.
    /// total/absolute copy number 0. Source: Cheng et al. (2017) Nat Commun 8:1221 — homozygous deletions are
    /// "regions having zero copies of both alleles in the tumour cells"; cBioPortal discrete-CNA scale — "−2"
    /// (Deep Deletion) is "a deep loss, possibly a homozygous deletion" (the deepest discrete loss), mapping to
    /// the integer copy-number 0 (CNVkit <c>absolute_threshold</c> DEL(0), shared with ONCO-CNA-001).
    /// </summary>
    private const int HomozygousDeletionCopyNumber = 0;

    /// <summary>
    /// Tests whether an arm-anchored segment is a homozygous (deep) deletion: its hard-threshold integer copy
    /// number is 0 (DeepDeletion). A single-copy loss (integer CN 1, cBioPortal "−1" shallow / heterozygous) is
    /// NOT a homozygous deletion. Source: Cheng et al. (2017) (total CN 0 = both alleles lost); cBioPortal
    /// (−2 = Deep Deletion); CNVkit <c>absolute_threshold</c> integer-CN calling (via <see cref="CallCopyNumber"/>).
    /// </summary>
    /// <param name="segment">The arm-anchored copy-number segment.</param>
    /// <param name="thresholds">
    /// Exactly four strictly ascending log2 cutoffs partitioning states 0/1/2/3/4+; null uses
    /// <see cref="DefaultCopyNumberThresholds"/> (CNVkit −1.1, −0.25, 0.2, 0.7).
    /// </param>
    /// <param name="ploidy">Reference (germline) ploidy; 2 for an autosomal diploid genome.</param>
    /// <returns><c>true</c> when the segment's integer copy number is 0.</returns>
    /// <exception cref="ArgumentException"><paramref name="segment"/> has non-positive arm length or End ≤ Start; or invalid thresholds.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static bool IsHomozygousDeletion(
        in CopyNumberArmSegment segment,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        ValidateArmSegment(segment);
        return CallCopyNumber(segment.Log2Ratio, thresholds, ploidy) == HomozygousDeletionCopyNumber;
    }

    /// <summary>
    /// Detects homozygous (deep) deletions among arm-anchored copy-number segments. A segment is reported when
    /// its hard-threshold integer copy number is 0 — total copy number 0, i.e. both alleles lost — which is the
    /// cBioPortal "−2" Deep Deletion / DeepDeletion state. Single-copy (heterozygous) losses, neutral, gain and
    /// amplification segments are excluded. The result is a subset of the input in input order (order-preserving
    /// filter). Source: Cheng et al. (2017) Nat Commun 8:1221 (homozygous = zero copies of both alleles);
    /// cBioPortal discrete-CNA scale; CNVkit <c>absolute_threshold</c> integer-CN calling.
    /// </summary>
    /// <param name="segments">Arm-anchored copy-number segments. Must not be null.</param>
    /// <param name="thresholds">Four strictly ascending log2 cutoffs; null uses CNVkit defaults (−1.1, −0.25, 0.2, 0.7).</param>
    /// <param name="ploidy">Reference (germline) ploidy; 2 for an autosomal diploid genome.</param>
    /// <returns>The homozygous-deletion segments, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive arm length or End ≤ Start; or invalid thresholds.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not positive.</exception>
    public static IReadOnlyList<CopyNumberArmSegment> DetectHomozygousDeletions(
        IEnumerable<CopyNumberArmSegment> segments,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        ArgumentNullException.ThrowIfNull(segments);

        var result = new List<CopyNumberArmSegment>();
        foreach (CopyNumberArmSegment segment in segments)
        {
            if (IsHomozygousDeletion(segment, thresholds, ploidy))
            {
                result.Add(segment);
            }
        }

        return result;
    }

    /// <summary>
    /// Maps homozygous-deletion segments to the recurrently deleted tumour suppressors resident on their
    /// chromosome arms. Each gene is reported once if any homozygous deletion falls on its arm. The panel and
    /// arms are: TP53 (17p), RB1 (13q), CDKN2A (9p), PTEN (10q), BRCA1 (17q), BRCA2 (13q). Source: NCBI Gene
    /// cytogenetic locations — TP53 17p13.1, RB1 13q14.2, CDKN2A 9p21.3, PTEN 10q23.31, BRCA1 17q21.31,
    /// BRCA2 13q13.1; tumour-suppressor role of recurrent homozygous deletions per Cheng et al. (2017).
    /// </summary>
    /// <param name="deletions">Homozygous deletions (typically the output of <see cref="DetectHomozygousDeletions"/>).</param>
    /// <returns>Distinct tumour-suppressor symbols whose arm carries a homozygous deletion, in panel order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="deletions"/> is null.</exception>
    public static IReadOnlyList<string> IdentifyDeletedTumorSuppressors(
        IEnumerable<CopyNumberArmSegment> deletions)
    {
        ArgumentNullException.ThrowIfNull(deletions);

        var deletedArms = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (CopyNumberArmSegment segment in deletions)
        {
            if (!string.IsNullOrEmpty(segment.Arm))
            {
                deletedArms.Add(segment.Arm);
            }
        }

        var genes = new List<string>();
        foreach ((string gene, string arm) in TumorSuppressorArms)
        {
            if (deletedArms.Contains(arm))
            {
                genes.Add(gene);
            }
        }

        return genes;
    }

    /// <summary>
    /// Recurrently deleted tumour suppressors and their chromosome arms (chromosome + arm letter), from NCBI
    /// Gene cytogenetic locations. Order is the registry panel order. Source: NCBI Gene — TP53 17p13.1 (Gene ID
    /// 7157), RB1 13q14.2 (5925), CDKN2A 9p21.3 (1029), PTEN 10q23.31 (5728), BRCA1 17q21.31 (672), BRCA2
    /// 13q13.1 (675).
    /// </summary>
    private static readonly IReadOnlyList<(string Gene, string Arm)> TumorSuppressorArms = new[]
    {
        ("TP53", "17p"),
        ("RB1", "13q"),
        ("CDKN2A", "9p"),
        ("PTEN", "10q"),
        ("BRCA1", "17q"),
        ("BRCA2", "13q"),
    };

    #endregion

    #region Tumor Ploidy Estimation (ONCO-PLOIDY-001)

    /// <summary>
    /// Minimum major-allele copy number for a segment to count as "elevated" toward whole-genome doubling.
    /// Source: facets-suite <c>is_genome_doubled</c> (<c>segs$mcn &gt;= 2</c>; PMID 30013179, Bielski et al.
    /// 2018) — WGD is assessed on the major copy number, where <c>mcn = tcn - lcn</c>.
    /// </summary>
    private const int WholeGenomeDoublingMajorCopyNumber = 2;

    /// <summary>
    /// Genome fraction (by length) at major copy number ≥ 2 above which a tumour is called whole-genome doubled.
    /// Source: facets-suite <c>is_genome_doubled(..., treshold = 0.5)</c> (PMID 30013179, Bielski et al. 2018):
    /// <c>wgd = frac_elevated_mcn &gt; treshold</c> — strictly greater than half of the genome.
    /// </summary>
    private const double WholeGenomeDoublingFractionThreshold = 0.5;

    /// <summary>
    /// Reference human genome assembly whose chromosome-size table is used as the denominator of the
    /// whole-genome-doubling genome fraction. Source: facets-suite <c>is_genome_doubled(segs, chrom_info, ...)</c>
    /// is parameterised by a <c>genome</c> build (<c>'hg19' | 'hg18' | 'hg38'</c>), each supplying its own
    /// chromosome-size object; the denominator <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c> is the
    /// reference assembly's autosomal length, NOT the interrogated-segment length.
    /// </summary>
    public enum ReferenceGenome
    {
        /// <summary>GRCh38 / hg38 (the current human reference assembly).</summary>
        GRCh38,

        /// <summary>GRCh37 / hg19 (the legacy human reference assembly).</summary>
        GRCh37,
    }

    /// <summary>
    /// Autosomal chromosome lengths (chromosomes 1–22, base pairs) of GRCh38 / hg38, indexed by chromosome
    /// number (entry 0 = chr1 … entry 21 = chr22). Embedded published reference data. Source: UCSC
    /// <c>hg38.chrom.sizes</c> (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes,
    /// retrieved 2026-06-22), cross-verified against the Ensembl REST assembly endpoint for GRCh38.p14
    /// (https://rest.ensembl.org/info/assembly/homo_sapiens — chr1 248,956,422; chr21 46,709,983; chr22
    /// 50,818,468; chrX 156,040,895). Only autosomes are used for the WGD denominator (facets-suite restricts to
    /// <c>chrom %in% 1:22</c>).
    /// </summary>
    private static readonly long[] GRCh38AutosomeLengths =
    {
        248_956_422L, // chr1
        242_193_529L, // chr2
        198_295_559L, // chr3
        190_214_555L, // chr4
        181_538_259L, // chr5
        170_805_979L, // chr6
        159_345_973L, // chr7
        145_138_636L, // chr8
        138_394_717L, // chr9
        133_797_422L, // chr10
        135_086_622L, // chr11
        133_275_309L, // chr12
        114_364_328L, // chr13
        107_043_718L, // chr14
        101_991_189L, // chr15
        90_338_345L,  // chr16
        83_257_441L,  // chr17
        80_373_285L,  // chr18
        58_617_616L,  // chr19
        64_444_167L,  // chr20
        46_709_983L,  // chr21
        50_818_468L,  // chr22
    };

    /// <summary>
    /// Autosomal chromosome lengths (chromosomes 1–22, base pairs) of GRCh37 / hg19, indexed by chromosome
    /// number (entry 0 = chr1 … entry 21 = chr22). Embedded published reference data. Source: UCSC
    /// <c>hg19.chrom.sizes</c> (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes,
    /// retrieved 2026-06-22). Only autosomes are used for the WGD denominator (facets-suite restricts to
    /// <c>chrom %in% 1:22</c>).
    /// </summary>
    private static readonly long[] GRCh37AutosomeLengths =
    {
        249_250_621L, // chr1
        243_199_373L, // chr2
        198_022_430L, // chr3
        191_154_276L, // chr4
        180_915_260L, // chr5
        171_115_067L, // chr6
        159_138_663L, // chr7
        146_364_022L, // chr8
        141_213_431L, // chr9
        135_534_747L, // chr10
        135_006_516L, // chr11
        133_851_895L, // chr12
        115_169_878L, // chr13
        107_349_540L, // chr14
        102_531_392L, // chr15
        90_354_753L,  // chr16
        81_195_210L,  // chr17
        78_077_248L,  // chr18
        59_128_983L,  // chr19
        63_025_520L,  // chr20
        48_129_895L,  // chr21
        51_304_566L,  // chr22
    };

    /// <summary>Number of autosomes in the human genome (chromosomes 1–22). Trivial structural constant.</summary>
    private const int AutosomeCount = 22;

    /// <summary>
    /// Returns the embedded autosomal chromosome-length table (chromosomes 1–22, base pairs) for a reference
    /// assembly, indexed 0 = chr1 … 21 = chr22. Source: UCSC <c>*.chrom.sizes</c> (see
    /// <see cref="GRCh38AutosomeLengths"/> / <see cref="GRCh37AutosomeLengths"/>).
    /// </summary>
    /// <param name="genome">The reference assembly.</param>
    /// <returns>The 22-element autosome length table for the assembly.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static IReadOnlyList<long> GetAutosomeLengths(ReferenceGenome genome) => genome switch
    {
        ReferenceGenome.GRCh38 => GRCh38AutosomeLengths,
        ReferenceGenome.GRCh37 => GRCh37AutosomeLengths,
        _ => throw new ArgumentOutOfRangeException(nameof(genome), genome, "Unknown reference genome."),
    };

    /// <summary>
    /// Total autosomal genome length (Σ of chromosome-1–22 lengths, base pairs) of a reference assembly — the
    /// denominator of the whole-genome-doubling genome fraction. Source: facets-suite
    /// <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c>; sizes from UCSC <c>*.chrom.sizes</c>.
    /// GRCh38 = 2,875,001,522 bp; GRCh37 = 2,881,033,286 bp.
    /// </summary>
    /// <param name="genome">The reference assembly.</param>
    /// <returns>The summed autosomal length in base pairs.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static long GetAutosomalGenomeLength(ReferenceGenome genome)
    {
        IReadOnlyList<long> lengths = GetAutosomeLengths(genome);
        long sum = 0L;
        for (int i = 0; i < lengths.Count; i++)
        {
            sum += lengths[i];
        }

        return sum;
    }

    /// <summary>
    /// Parses a chromosome identifier to its autosome number (1–22), accepting both bare ("7") and "chr"-prefixed
    /// ("chr7") forms. Returns <c>false</c> for sex chromosomes, mitochondria, contigs, or anything outside 1–22,
    /// which the WGD fraction excludes (facets-suite <c>chrom %in% 1:22</c>).
    /// </summary>
    /// <param name="chromosome">The chromosome identifier from a segment.</param>
    /// <param name="number">The parsed autosome number (1–22) when the method returns <c>true</c>.</param>
    /// <returns><c>true</c> when the identifier denotes an autosome (1–22).</returns>
    private static bool TryGetAutosomeNumber(string? chromosome, out int number)
    {
        number = 0;
        if (string.IsNullOrEmpty(chromosome))
        {
            return false;
        }

        ReadOnlySpan<char> name = chromosome.AsSpan();
        if (name.Length > 3 &&
            (name[0] is 'c' or 'C') && (name[1] is 'h' or 'H') && (name[2] is 'r' or 'R'))
        {
            name = name[3..];
        }

        return int.TryParse(name, out number) && number is >= 1 and <= AutosomeCount;
    }

    /// <summary>
    /// Estimates the average tumour ploidy ψ as the segment-length-weighted mean of per-segment total copy
    /// number: ψ = Σ(CN_i · L_i) / Σ(L_i), where CN_i = MajorCopyNumber + MinorCopyNumber and L_i = End − Start.
    /// Source: Patchwork (Genome Biology) — "The average ploidy, PloidyTum, is the average total copy number of
    /// all genomic segments weighted by segment length"; the originating allele-specific method is ASCAT
    /// (Van Loo et al., PNAS 2010, 10.1073/pnas.1009843107), which reports a final tumour ploidy on the n-scale
    /// (2n = diploid). A pure-diploid (all 1:1) genome has ψ = 2.0; ">2.7n" marks aneuploidy (Van Loo et al.).
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments (the <see cref="AlleleSpecificSegment"/> shared with ONCO-LOH-001 /
    /// ONCO-HRD-001). Per-segment total copy number is Major + Minor; length is End − Start. Must not be null,
    /// must be non-empty, and every segment must have End &gt; Start and non-negative copy numbers.
    /// </param>
    /// <returns>The length-weighted average ploidy ψ (&gt; 0 for any genome with at least one positive copy number).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="segments"/> is empty (ploidy is undefined for an empty genome), or a segment has
    /// End ≤ Start or a negative copy number.
    /// </exception>
    public static double EstimatePloidy(IEnumerable<AlleleSpecificSegment> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);

        double weightedCopyNumberSum = 0.0;
        long totalLength = 0L;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            long length = segment.Length;
            int totalCopyNumber = segment.MajorCopyNumber + segment.MinorCopyNumber;
            weightedCopyNumberSum += (double)totalCopyNumber * length;
            totalLength += length;
        }

        if (totalLength == 0L)
        {
            throw new ArgumentException(
                "Cannot estimate ploidy from an empty segment set (the length-weighted mean is undefined).",
                nameof(segments));
        }

        // ψ = Σ(CN_i · L_i) / Σ(L_i) — Patchwork length-weighted mean of total copy number.
        return weightedCopyNumberSum / totalLength;
    }

    /// <summary>
    /// Determines whether a tumour genome has undergone whole-genome doubling (WGD), computing the genome
    /// fraction against a <b>reference chromosome-size table</b> (the authoritative autosomal genome length),
    /// exactly as facets-suite does. WGD is called when the fraction of the <i>reference autosomal genome</i>
    /// (chromosomes 1–22) covered by segments with major-allele copy number ≥ 2 is strictly greater than 0.5.
    /// Source: facets-suite <c>is_genome_doubled(segs, chrom_info, treshold = 0.5)</c> (PMID 30013179, Bielski
    /// et al. 2018, Nat Genet 50:1189–1195):
    /// <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c>;
    /// <c>frac_elevated_mcn = sum(length where mcn ≥ 2 &amp; chrom %in% 1:22) / autosomal_genome</c>;
    /// <c>wgd = frac_elevated_mcn &gt; treshold</c>, with <c>mcn = tcn − lcn</c> (major-allele copy number).
    /// Because the denominator is the true genome length (not the sum of supplied segments), segments that do not
    /// tile the genome no longer bias the fraction; only autosomal (chr1–22) segments contribute to the numerator
    /// (sex chromosomes / contigs are ignored). The test uses the major (not total) copy number, so a balanced
    /// diploid genome (all 1:1, total CN 2, major CN 1) is NOT doubled, whereas a 2:0 LOH or 2:2 genome IS.
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments (<see cref="AlleleSpecificSegment"/>). Only segments on autosomes
    /// (chromosomes 1–22, "chr"-prefixed or bare) contribute to the elevated-major-CN numerator. Must not be
    /// null, and every segment must have End &gt; Start and non-negative copy numbers.
    /// </param>
    /// <param name="genome">
    /// Reference assembly whose autosomal chromosome-size table is the fraction denominator
    /// (default <see cref="ReferenceGenome.GRCh38"/>).
    /// </param>
    /// <returns><c>true</c> when more than half the reference autosomal genome has major copy number ≥ 2.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has End ≤ Start or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static bool DetectWholeGenomeDoubling(
        IEnumerable<AlleleSpecificSegment> segments,
        ReferenceGenome genome = ReferenceGenome.GRCh38)
    {
        ArgumentNullException.ThrowIfNull(segments);

        long autosomalGenomeLength = GetAutosomalGenomeLength(genome);

        long elevatedLength = 0L;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            // facets-suite: numerator restricted to autosomes (chrom %in% 1:22). Non-autosomal segments are
            // ignored (sex chromosomes / contigs do not contribute to the autosomal WGD fraction).
            if (!TryGetAutosomeNumber(segment.Chromosome, out _))
            {
                continue;
            }

            // mcn = major-allele copy number; elevated when major CN ≥ 2 (facets-suite segs$mcn >= 2).
            if (segment.MajorCopyNumber >= WholeGenomeDoublingMajorCopyNumber)
            {
                elevatedLength += segment.Length;
            }
        }

        // wgd = frac_elevated_mcn > 0.5 (strict), denominator = reference autosomal genome length.
        double fractionElevatedMajorCn = (double)elevatedLength / autosomalGenomeLength;
        return fractionElevatedMajorCn > WholeGenomeDoublingFractionThreshold;
    }

    /// <summary>
    /// Determines whole-genome doubling using the <b>supplied segments' total length</b> as the genome-fraction
    /// denominator (the legacy behaviour), rather than a reference chromosome-size table. This is correct only
    /// when the supplied segments tile the interrogated (autosomal) genome; otherwise prefer the reference-table
    /// overload <see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment}, ReferenceGenome)"/>.
    /// WGD is called when Σ(length where major CN ≥ 2) ÷ Σ(all supplied segment length) is strictly greater than
    /// 0.5. Source: facets-suite <c>is_genome_doubled</c> rule (PMID 30013179) applied with the interrogated
    /// segments as the denominator; <c>mcn = tcn − lcn</c>.
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments. The fraction denominator is the total length of <b>all</b> supplied
    /// segments (the interrogated genome), regardless of chromosome. Must not be null, must be non-empty, and
    /// every segment must have End &gt; Start and non-negative copy numbers.
    /// </param>
    /// <returns><c>true</c> when more than half the supplied genome (by length) has major copy number ≥ 2.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="segments"/> is empty (the fraction is undefined), or a segment has End ≤ Start or a
    /// negative copy number.
    /// </exception>
    public static bool DetectWholeGenomeDoublingFromSuppliedLength(IEnumerable<AlleleSpecificSegment> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);

        long elevatedLength = 0L;
        long totalLength = 0L;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            long length = segment.Length;
            totalLength += length;
            // mcn = major-allele copy number; elevated when major CN ≥ 2 (facets-suite segs$mcn >= 2).
            if (segment.MajorCopyNumber >= WholeGenomeDoublingMajorCopyNumber)
            {
                elevatedLength += length;
            }
        }

        if (totalLength == 0L)
        {
            throw new ArgumentException(
                "Cannot assess whole-genome doubling from an empty segment set (the genome fraction is undefined).",
                nameof(segments));
        }

        // wgd = frac_elevated_mcn > 0.5 (strict) — facets-suite is_genome_doubled, supplied-length denominator.
        double fractionElevatedMajorCn = (double)elevatedLength / totalLength;
        return fractionElevatedMajorCn > WholeGenomeDoublingFractionThreshold;
    }

    #endregion

    #region Clonal vs Subclonal Classification (ONCO-CLONAL-001)

    /// <summary>
    /// Cancer-cell-fraction (CCF) threshold separating clonal from subclonal mutations. A mutation is called
    /// clonal when its CCF exceeds this value (with sufficient posterior probability); subclonal otherwise.
    /// Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — "We classified a mutation as clonal if the
    /// CCF harboring it was &gt;0.95 with probability &gt; 0.5, and subclonal otherwise."
    /// </summary>
    public const double ClonalCcfThreshold = 0.95;

    /// <summary>
    /// Minimum posterior probability that the CCF exceeds <see cref="ClonalCcfThreshold"/> required to call a
    /// mutation clonal. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — clonal if P(CCF &gt; 0.95) &gt; 0.5.
    /// </summary>
    public const double ClonalProbabilityThreshold = 0.5;

    /// <summary>
    /// Number of CCF grid points used to evaluate the posterior over the CCF c ∈ [0.01, 1]. Source: Landau et al.
    /// (2013), <i>Cell</i> 152(4):714–726, Extended Experimental Procedures — "calculating these values over a
    /// regular grid of 100 c values and normalizing by dividing them by their sum".
    /// </summary>
    private const int CcfGridPointCount = 100;

    /// <summary>
    /// Lower bound of the CCF grid. Source: Landau et al. (2013) — "with c ∈ [0.01,1]" (a mutation is present in
    /// at least one cancer cell, so the CCF cannot be exactly zero).
    /// </summary>
    private const double CcfGridLowerBound = 0.01;

    /// <summary>Upper bound of the CCF grid (a mutation present in every cancer cell). Source: Landau et al. (2013), c ∈ [0.01,1].</summary>
    private const double CcfGridUpperBound = 1.0;

    /// <summary>Clonal/subclonal status of one somatic mutation.</summary>
    public enum ClonalityStatus
    {
        /// <summary>Present in (essentially) all cancer cells: CCF &gt; 0.95 with posterior probability &gt; 0.5.</summary>
        Clonal,

        /// <summary>Present in only a subpopulation of cancer cells: it does not meet the clonal criterion.</summary>
        Subclonal,
    }

    /// <summary>
    /// Read evidence for one somatic mutation at a locus with known purity-corrected copy-number state, used to
    /// infer its cancer cell fraction. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — the posterior
    /// over CCF is built from <c>a</c> alternate reads out of <c>N</c> total reads at a locus of absolute somatic
    /// copy number <c>q</c>, with mutation multiplicity <c>M</c> (the number of tumour-genome copies carrying the
    /// mutant allele; DeCiFering / Satas et al. 2021, <i>Cell Systems</i> 12(10):1004–1018, Eq. 1).
    /// </summary>
    /// <param name="AltReads">Alternate-allele supporting reads <c>a</c> (≥ 0, ≤ <paramref name="TotalReads"/>).</param>
    /// <param name="TotalReads">Total reads at the locus <c>N</c> (≥ 1).</param>
    /// <param name="LocalCopyNumber">Absolute tumour total copy number at the locus <c>q</c> (≥ 1).</param>
    /// <param name="Multiplicity">Mutation multiplicity <c>M</c> — tumour-genome copies carrying the mutation (≥ 1, ≤ <paramref name="LocalCopyNumber"/>).</param>
    public readonly record struct ClonalityVariant(int AltReads, int TotalReads, int LocalCopyNumber, int Multiplicity)
    {
        /// <summary>Creates a variant at multiplicity 1 (one mutated copy), the default for a heterozygous SNV.</summary>
        public ClonalityVariant(int altReads, int totalReads, int localCopyNumber)
            : this(altReads, totalReads, localCopyNumber, 1)
        {
        }
    }

    /// <summary>The clonality classification of one mutation, with the CCF point estimate and posterior that produced it.</summary>
    /// <param name="Variant">The variant that was classified.</param>
    /// <param name="Ccf">Posterior-mean cancer cell fraction estimate (grid expectation), in [0.01, 1].</param>
    /// <param name="ProbabilityClonal">Posterior probability that the CCF exceeds <see cref="ClonalCcfThreshold"/>.</param>
    /// <param name="Status">Clonal / Subclonal classification.</param>
    public readonly record struct ClonalityCall(ClonalityVariant Variant, double Ccf, double ProbabilityClonal, ClonalityStatus Status);

    /// <summary>Summary of a clonal/subclonal classification over a set of variants.</summary>
    /// <param name="Calls">Per-variant classifications, in input order.</param>
    /// <param name="ClonalCount">Number of variants classified <see cref="ClonalityStatus.Clonal"/>.</param>
    /// <param name="SubclonalCount">Number of variants classified <see cref="ClonalityStatus.Subclonal"/>.</param>
    /// <param name="ClonalFraction">Fraction of variants that are clonal (ClonalCount / total); 0 for an empty set.</param>
    public readonly record struct ClonalityResult(
        IReadOnlyList<ClonalityCall> Calls,
        int ClonalCount,
        int SubclonalCount,
        double ClonalFraction);

    /// <summary>
    /// Classifies each somatic mutation as clonal or subclonal from its read evidence and the tumour purity.
    /// For each variant the posterior over the cancer cell fraction c is built on a regular grid of
    /// <see cref="CcfGridPointCount"/> points c ∈ [<see cref="CcfGridLowerBound"/>, 1] as
    /// <c>P(c) ∝ Binomial(a | N, f(c))</c> with a uniform prior, where the expected alternate-allele fraction is
    /// <c>f(c) = ρ·M·c / (2(1−ρ) + ρ·q)</c> — purity ρ, multiplicity M, local copy number q, normal diploid
    /// contribution 2(1−ρ). A mutation is called <see cref="ClonalityStatus.Clonal"/> when the posterior
    /// probability that c &gt; <see cref="ClonalCcfThreshold"/> exceeds <see cref="ClonalProbabilityThreshold"/>,
    /// and <see cref="ClonalityStatus.Subclonal"/> otherwise. Source: Landau et al. (2013), <i>Cell</i>
    /// 152(4):714–726 (Extended Experimental Procedures); multiplicity generalisation per Satas et al. (2021),
    /// <i>Cell Systems</i> 12(10):1004–1018 (Eq. 1).
    /// </summary>
    /// <param name="variants">Per-variant read evidence and copy-number state.</param>
    /// <param name="purity">Tumor purity ρ ∈ (0, 1].</param>
    /// <returns>Per-variant calls plus clonal/subclonal counts and the clonal fraction.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">purity ∉ (0, 1].</exception>
    /// <exception cref="ArgumentException">A variant has invalid read counts, copy number, or multiplicity.</exception>
    public static ClonalityResult ClassifyClonality(IEnumerable<ClonalityVariant> variants, double purity)
    {
        ArgumentNullException.ThrowIfNull(variants);
        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in the range (0, 1]; the CCF model divides by purity.");
        }

        var calls = new List<ClonalityCall>();
        int clonalCount = 0;
        foreach (ClonalityVariant variant in variants)
        {
            ClonalityCall call = ClassifyOne(variant, purity);
            calls.Add(call);
            if (call.Status == ClonalityStatus.Clonal)
            {
                clonalCount++;
            }
        }

        int total = calls.Count;
        int subclonalCount = total - clonalCount;
        // ClonalFraction is the clonal share; undefined for an empty set so reported as 0.
        double clonalFraction = total == 0 ? 0.0 : (double)clonalCount / total;
        return new ClonalityResult(calls, clonalCount, subclonalCount, clonalFraction);
    }

    /// <summary>
    /// Classifies a set of already-estimated cancer cell fractions (CCF point estimates, e.g. from
    /// <c>EstimateCCF</c>) as clonal vs subclonal and returns the indices of the clonal mutations. A CCF is
    /// clonal when it exceeds <see cref="ClonalCcfThreshold"/> (CCF &gt; 0.95), reflecting a mutation present in
    /// (essentially) all cancer cells. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — "classified a
    /// mutation as clonal if the CCF harboring it was &gt;0.95 … and subclonal otherwise".
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions, each in [0, 1].</param>
    /// <returns>The 0-based indices of the clonal CCF values, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">A CCF value is NaN or outside [0, 1].</exception>
    public static IReadOnlyList<int> IdentifyClonalMutations(IEnumerable<double> ccfValues)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);

        var clonalIndices = new List<int>();
        int index = 0;
        foreach (double ccf in ccfValues)
        {
            if (double.IsNaN(ccf) || ccf < 0.0 || ccf > 1.0)
            {
                throw new ArgumentException(
                    $"Cancer cell fraction must be in [0, 1]; got {ccf} at index {index}.", nameof(ccfValues));
            }

            if (ccf > ClonalCcfThreshold)
            {
                clonalIndices.Add(index);
            }

            index++;
        }

        return clonalIndices;
    }

    /// <summary>
    /// Builds the posterior over CCF for one variant and classifies it (Landau et al. 2013 grid model).
    /// </summary>
    private static ClonalityCall ClassifyOne(ClonalityVariant variant, double purity)
    {
        ValidateClonalityVariant(variant);

        // Expected alt-allele fraction f(c) = ρ·M·c / (2(1−ρ) + ρ·q): mutant copies M·c per cell scaled by
        // purity over the total DNA (normal 2(1−ρ) + tumour ρ·q). Landau 2013 (M=1) generalised by DeCiFering Eq.1.
        double denominator = NormalDiploidCopyNumber * (1.0 - purity) + purity * variant.LocalCopyNumber;
        double alleleFractionPerUnitCcf = purity * variant.Multiplicity / denominator;

        // Posterior P(c) ∝ Binomial(a | N, f(c)) on a uniform grid c ∈ [0.01, 1], uniform prior; normalise by sum.
        // The binomial coefficient C(N, a) is constant in c, so it cancels in normalisation and is omitted.
        double step = (CcfGridUpperBound - CcfGridLowerBound) / (CcfGridPointCount - 1);
        Span<double> weights = stackalloc double[CcfGridPointCount];
        double weightSum = 0.0;
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double c = CcfGridLowerBound + step * i;
            double f = Math.Min(1.0, alleleFractionPerUnitCcf * c);
            double likelihood = BinomialLikelihoodKernel(variant.AltReads, variant.TotalReads, f);
            weights[i] = likelihood;
            weightSum += likelihood;
        }

        double ccfMean = 0.0;
        double probabilityClonal = 0.0;
        // Guard against an all-zero posterior (e.g. f≈0 with a>0): fall back to a flat posterior over the grid.
        bool degenerate = weightSum <= 0.0 || double.IsNaN(weightSum);
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double c = CcfGridLowerBound + step * i;
            double posterior = degenerate ? 1.0 / CcfGridPointCount : weights[i] / weightSum;
            ccfMean += c * posterior;
            if (c > ClonalCcfThreshold)
            {
                probabilityClonal += posterior;
            }
        }

        ClonalityStatus status = probabilityClonal > ClonalProbabilityThreshold
            ? ClonalityStatus.Clonal
            : ClonalityStatus.Subclonal;

        // The reported posterior summaries are bounded by their documented invariants: the CCF point estimate is a
        // grid expectation in [0.01, 1] (INV-03) and the clonal probability is normalised posterior mass in [0, 1]
        // (INV-04). Summing the normalised grid weights can overshoot the bound by one ulp (e.g. 1.0000000000000002),
        // so clamp the reported values to their invariant range. The unclamped probabilityClonal already determined
        // status above, and clamping a near-1 value cannot change the > 0.5 decision.
        double reportedCcf = Math.Clamp(ccfMean, CcfGridLowerBound, CcfGridUpperBound);
        double reportedProbabilityClonal = Math.Clamp(probabilityClonal, 0.0, 1.0);
        return new ClonalityCall(variant, reportedCcf, reportedProbabilityClonal, status);
    }

    /// <summary>
    /// Binomial likelihood kernel L(a | N, p) = p^a · (1−p)^(N−a), without the constant C(N, a) factor (it cancels
    /// under grid normalisation). Computed via log-space to avoid underflow for large N.
    /// </summary>
    private static double BinomialLikelihoodKernel(int altReads, int totalReads, double p)
    {
        int refReads = totalReads - altReads;
        if (p <= 0.0)
        {
            // p = 0 explains zero alternate reads exactly, nothing else.
            return altReads == 0 ? 1.0 : 0.0;
        }

        if (p >= 1.0)
        {
            // p = 1 explains all-alternate reads exactly, nothing else.
            return refReads == 0 ? 1.0 : 0.0;
        }

        double logLikelihood = altReads * Math.Log(p) + refReads * Math.Log(1.0 - p);
        return Math.Exp(logLikelihood);
    }

    /// <summary>Validates read counts, local copy number, and multiplicity of a clonality variant.</summary>
    private static void ValidateClonalityVariant(ClonalityVariant variant)
    {
        if (variant.TotalReads < 1)
        {
            throw new ArgumentException(
                $"Total reads must be at least 1; got {variant.TotalReads}.", nameof(variant));
        }

        if (variant.AltReads < 0 || variant.AltReads > variant.TotalReads)
        {
            throw new ArgumentException(
                $"Alternate reads must be in [0, {variant.TotalReads}]; got {variant.AltReads}.", nameof(variant));
        }

        if (variant.LocalCopyNumber < 1)
        {
            throw new ArgumentException(
                $"Local copy number must be at least 1; got {variant.LocalCopyNumber}.", nameof(variant));
        }

        if (variant.Multiplicity < 1 || variant.Multiplicity > variant.LocalCopyNumber)
        {
            throw new ArgumentException(
                $"Multiplicity must be in [1, {variant.LocalCopyNumber}]; got {variant.Multiplicity}.", nameof(variant));
        }
    }

    #endregion

    #region EstimateCcf

    /// <summary>
    /// Normal locus copy number (diploid) contributing the 2(1−ρ) term in the CCF denominator.
    /// Source: McGranahan et al. (2016), <i>Science</i> 351(6280):1463–1469 — n_mut = VAF·(1/p)·[p·CN_t + CN_n·(1−p)]
    /// with CN_n = 2; Tarabichi et al. (2021), <i>Nat. Methods</i> 18:144–155 (Box 1).
    /// </summary>
    private const double NormalLocusCopyNumber = 2.0;

    /// <summary>Upper bound on a reported cancer cell fraction (a mutation in all cancer cells has CCF = 1).</summary>
    private const double MaxCancerCellFraction = 1.0;

    /// <summary>
    /// A single cancer-cell-fraction point estimate for one somatic mutation.
    /// </summary>
    /// <param name="Ccf">Reported cancer cell fraction, capped to [0, 1] (the registry invariant; a mutation
    /// present in all cancer cells has CCF = 1, per McGranahan et al. 2016).</param>
    /// <param name="RawCcf">Uncapped formula value VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m); may exceed 1 under sampling
    /// noise (CNAqc reports e.g. 1.06).</param>
    public readonly record struct CcfEstimate(double Ccf, double RawCcf);

    /// <summary>
    /// Result of clustering cancer cell fractions into clones/subclones by deterministic 1D k-means.
    /// </summary>
    /// <param name="Centroids">Cluster centroids (means) sorted in ascending order, one per cluster.</param>
    /// <param name="Assignments">For each input CCF (input order), the 0-based index of its assigned cluster
    /// in <paramref name="Centroids"/>.</param>
    /// <param name="ClonalClusterIndex">Index (into <paramref name="Centroids"/>) of the clonal cluster — the
    /// cluster with the highest centroid (Tarabichi et al. 2021: "the cluster with the highest CP … deemed clonal").</param>
    public readonly record struct CcfClustering(
        IReadOnlyList<double> Centroids,
        IReadOnlyList<int> Assignments,
        int ClonalClusterIndex);

    /// <summary>
    /// Estimates the cancer cell fraction (CCF) of a somatic mutation from its variant allele fraction, the tumor
    /// purity, the local tumor copy number, and the mutation multiplicity, using the standard point estimate
    /// <c>CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m)</c> — VAF the variant allele fraction, ρ the purity, N_T the local
    /// tumor copy number, the normal contributing 2(1−ρ), and m the integer mutation multiplicity (number of
    /// mutated copies per cancer cell). Source: McGranahan et al. (2016), <i>Science</i> 351(6280):1463–1469
    /// (n_mut = VAF·(1/p)·[p·CN_t + 2(1−p)], CCF = n_mut/m); Tarabichi et al. (2021), <i>Nat. Methods</i>
    /// 18:144–155 (Box 1); Zheng et al. (2022), <i>Bioinformatics</i> 38(15):3677–3683 (VAF = m·CCF·p/(c·p+2(1−p))).
    /// The reported <see cref="CcfEstimate.Ccf"/> is capped to [0, 1] to honour the 0 ≤ CCF ≤ 1 invariant while the
    /// uncapped value is exposed as <see cref="CcfEstimate.RawCcf"/>.
    /// </summary>
    /// <param name="vaf">Variant allele fraction (mutant read fraction), in [0, 1].</param>
    /// <param name="purity">Tumor purity ρ ∈ (0, 1].</param>
    /// <param name="tumorCopyNumber">Local tumor copy number N_T (≥ 1).</param>
    /// <param name="multiplicity">Mutation multiplicity m (number of mutated copies per cancer cell), in
    /// [1, <paramref name="tumorCopyNumber"/>].</param>
    /// <returns>The capped and raw cancer cell fraction.</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0,1], purity ∉ (0,1], or tumorCopyNumber &lt; 1.</exception>
    /// <exception cref="ArgumentException">multiplicity ∉ [1, tumorCopyNumber].</exception>
    public static CcfEstimate EstimateCcf(double vaf, double purity, int tumorCopyNumber, int multiplicity)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in (0, 1]; the CCF formula divides by purity.");
        }

        if (tumorCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tumorCopyNumber), tumorCopyNumber, "Tumor copy number must be at least 1.");
        }

        if (multiplicity < 1 || multiplicity > tumorCopyNumber)
        {
            throw new ArgumentException(
                $"Multiplicity must be in [1, {tumorCopyNumber}]; got {multiplicity}.", nameof(multiplicity));
        }

        // CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m): total DNA per cell = tumour ρ·N_T + normal 2(1−ρ); dividing the
        // observed mutant fraction by ρ·m / totalDna recovers the fraction of cancer cells carrying the mutation.
        double totalDnaPerCell = purity * tumorCopyNumber + NormalLocusCopyNumber * (1.0 - purity);
        double rawCcf = vaf * totalDnaPerCell / (purity * multiplicity);
        double cappedCcf = Math.Min(MaxCancerCellFraction, rawCcf);
        return new CcfEstimate(cappedCcf, rawCcf);
    }

    /// <summary>
    /// Clusters cancer cell fractions into <paramref name="clusterCount"/> clones/subclones using a deterministic
    /// one-dimensional Lloyd's k-means (Lloyd 1982, <i>IEEE Trans. Inf. Theory</i> 28(2):129–137): each value is
    /// assigned to the nearest centroid (least squared distance) and each centroid is recomputed as the mean of
    /// its members, iterating to convergence. Determinism is achieved without any RNG by seeding the centroids at
    /// evenly spaced quantiles of the sorted input. The clonal cluster is the one with the highest centroid
    /// (Tarabichi et al. 2021, <i>Nat. Methods</i> 18:144–155: "the cluster with the highest CP can be deemed clonal").
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions to cluster (each finite).</param>
    /// <param name="clusterCount">Number of clusters k, in [1, count of values].</param>
    /// <returns>Ascending-sorted centroids, per-value cluster assignments (input order), and the clonal cluster index.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, or a value is NaN/infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count].</exception>
    public static CcfClustering ClusterCcfValues(IReadOnlyList<double> ccfValues, int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);

        int n = ccfValues.Count;
        if (n == 0)
        {
            throw new ArgumentException("At least one CCF value is required.", nameof(ccfValues));
        }

        for (int i = 0; i < n; i++)
        {
            if (double.IsNaN(ccfValues[i]) || double.IsInfinity(ccfValues[i]))
            {
                throw new ArgumentException($"CCF value must be finite; got {ccfValues[i]} at index {i}.", nameof(ccfValues));
            }
        }

        if (clusterCount < 1 || clusterCount > n)
        {
            throw new ArgumentOutOfRangeException(
                nameof(clusterCount), clusterCount, $"Cluster count must be in [1, {n}].");
        }

        // Sort values (carrying original indices) so seeding and assignment are deterministic and order-independent.
        int[] order = Enumerable.Range(0, n).OrderBy(i => ccfValues[i]).ToArray();
        double[] sorted = order.Select(i => ccfValues[i]).ToArray();

        // Deterministic seeding: place centroid j at the value at quantile (j + 0.5)/k of the sorted data.
        double[] centroids = new double[clusterCount];
        for (int j = 0; j < clusterCount; j++)
        {
            int idx = (int)((j + 0.5) / clusterCount * n);
            if (idx >= n)
            {
                idx = n - 1;
            }

            centroids[j] = sorted[idx];
        }

        int[] sortedAssignments = new int[n];
        // Lloyd iterations: assignment step then update step; converges when no assignment changes.
        // Bounded by n iterations (each iteration strictly reduces WCSS until a fixed point on a finite set).
        for (int iteration = 0; iteration <= n; iteration++)
        {
            bool changed = AssignToNearestCentroid(sorted, centroids, sortedAssignments);
            RecomputeCentroids(sorted, sortedAssignments, centroids);
            if (!changed && iteration > 0)
            {
                break;
            }
        }

        // Map cluster labels to ascending-centroid order so the result is canonical and assignments are stable.
        int[] centroidRank = Enumerable.Range(0, clusterCount).OrderBy(j => centroids[j]).ToArray();
        int[] rankOfCluster = new int[clusterCount];
        double[] sortedCentroids = new double[clusterCount];
        for (int rank = 0; rank < clusterCount; rank++)
        {
            rankOfCluster[centroidRank[rank]] = rank;
            sortedCentroids[rank] = centroids[centroidRank[rank]];
        }

        int[] assignments = new int[n];
        for (int s = 0; s < n; s++)
        {
            assignments[order[s]] = rankOfCluster[sortedAssignments[s]];
        }

        // Clonal cluster = highest centroid; centroids are ascending so it is the last index.
        int clonalClusterIndex = clusterCount - 1;
        return new CcfClustering(sortedCentroids, assignments, clonalClusterIndex);
    }

    /// <summary>Assignment step: assigns each value to the nearest centroid; returns whether any assignment changed.</summary>
    private static bool AssignToNearestCentroid(double[] values, double[] centroids, int[] assignments)
    {
        bool changed = false;
        for (int i = 0; i < values.Length; i++)
        {
            int best = 0;
            double bestDistance = double.PositiveInfinity;
            for (int j = 0; j < centroids.Length; j++)
            {
                double diff = values[i] - centroids[j];
                double distance = diff * diff;
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    best = j;
                }
            }

            if (assignments[i] != best)
            {
                assignments[i] = best;
                changed = true;
            }
        }

        return changed;
    }

    /// <summary>Update step: recomputes each centroid as the mean of its assigned values (empty clusters keep their centroid).</summary>
    private static void RecomputeCentroids(double[] values, int[] assignments, double[] centroids)
    {
        Span<double> sums = stackalloc double[centroids.Length];
        Span<int> counts = stackalloc int[centroids.Length];
        for (int i = 0; i < values.Length; i++)
        {
            sums[assignments[i]] += values[i];
            counts[assignments[i]]++;
        }

        for (int j = 0; j < centroids.Length; j++)
        {
            if (counts[j] > 0)
            {
                centroids[j] = sums[j] / counts[j];
            }
        }
    }

    #endregion

    #region GenerateNeoantigenPeptides

    /// <summary>
    /// Shortest MHC class I peptide length (8-mer) over which candidate neoantigen windows are enumerated.
    /// Source: Jurtz et al. (2017) NetMHCpan-4.0, <i>J. Immunol.</i> 199(9):3360–3368 — class I predictions
    /// are made for peptides of length 8–14; the NetMHCpan-4.1 web service offers 8/9/10/11/12/13/14-mer
    /// peptide options (https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/). pVACtools restricts the
    /// canonical class I neoantigen search to lengths 8–11 (Hundal et al. 2020, <i>Cancer Immunol. Res.</i>
    /// 8(3):409–420 — "8–11-mer for Class I MHC").
    /// </summary>
    public const int MhcClassIMinPeptideLength = 8;

    /// <summary>
    /// Longest MHC class I peptide length (11-mer) over which candidate neoantigen windows are enumerated.
    /// Source: Hundal et al. (2020), <i>Cancer Immunol. Res.</i> 8(3):409–420 — pVACtools predicts the
    /// strongest MHC-binding peptides "(8–11-mer for Class I MHC …)".
    /// </summary>
    public const int MhcClassIMaxPeptideLength = 11;

    /// <summary>
    /// A candidate neoantigen peptide: a fixed-length window of the mutant protein that spans the mutated
    /// residue, paired with the wild-type peptide occupying the same coordinates (the agretope). Mutant and
    /// wild-type peptides differ only at the substituted residue (Hundal et al. 2020, <i>Cancer Immunol.
    /// Res.</i> 8(3):409–420; agretopicity is the differential binding of these two — Wells et al. 2020,
    /// <i>Cell</i> 183(3):818–834).
    /// </summary>
    /// <param name="Length">Peptide length k (number of residues), in [<see cref="MhcClassIMinPeptideLength"/>,
    /// <see cref="MhcClassIMaxPeptideLength"/>].</param>
    /// <param name="StartPosition">1-based position in the protein of the first residue of the window.</param>
    /// <param name="MutantPeptide">The k-mer taken from the mutant protein (carries the substituted residue).</param>
    /// <param name="WildTypePeptide">The k-mer at the same coordinates in the wild-type protein (the agretope).</param>
    /// <param name="MutationOffset">0-based offset of the mutated residue within the peptide window.</param>
    public readonly record struct NeoantigenPeptide(
        int Length,
        int StartPosition,
        string MutantPeptide,
        string WildTypePeptide,
        int MutationOffset);

    /// <summary>
    /// Generates the candidate MHC class I neoantigen peptide windows arising from a single somatic missense
    /// (amino-acid substitution) mutation. For each peptide length k the method enumerates every length-k
    /// window of the mutant protein that <b>spans</b> the mutated residue, and pairs it with the wild-type
    /// window at the same coordinates (the agretope). This is the windowing step of neoantigen prediction as
    /// defined by pVACtools (Hundal et al. 2020) and the 21-mer ± 10-flank construction of ProGeo-neo (Li et
    /// al. 2020, <i>BMC Med. Genomics</i> 13:52): a 21-mer with 10 residues flanking the substitution on each
    /// side contains exactly the 8–11-mer windows that overlap the mutation.
    /// <para>
    /// Binding affinity / IC50 is NOT computed here — that requires a trained MHC-binding model (e.g.
    /// NetMHCpan) and is caller-supplied or out of scope (ONCO-MHC-001).
    /// </para>
    /// </summary>
    /// <param name="wildTypeProtein">The wild-type (reference) protein sequence (one-letter amino-acid codes).</param>
    /// <param name="mutantResidue">The substituted (mutant) amino acid (one-letter code).</param>
    /// <param name="mutationPosition">1-based position of the substituted residue within the protein.</param>
    /// <param name="minLength">Minimum peptide length (default <see cref="MhcClassIMinPeptideLength"/>).</param>
    /// <param name="maxLength">Maximum peptide length (default <see cref="MhcClassIMaxPeptideLength"/>).</param>
    /// <returns>
    /// All candidate peptides, ordered by length ascending then by start position ascending. Empty only if no
    /// window of any requested length fits within the protein bounds while spanning the mutation.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="wildTypeProtein"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// The protein is empty; <paramref name="mutantResidue"/> is not a single character; the wild-type residue
    /// at the mutation position already equals the mutant residue (not a substitution); or the length range is
    /// invalid (min &gt; max, or min &lt; 1).
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="mutationPosition"/> is outside [1, protein length].
    /// </exception>
    public static IReadOnlyList<NeoantigenPeptide> GenerateNeoantigenPeptides(
        string wildTypeProtein,
        char mutantResidue,
        int mutationPosition,
        int minLength = MhcClassIMinPeptideLength,
        int maxLength = MhcClassIMaxPeptideLength)
    {
        ArgumentNullException.ThrowIfNull(wildTypeProtein);

        if (wildTypeProtein.Length == 0)
        {
            throw new ArgumentException("Protein sequence must be non-empty.", nameof(wildTypeProtein));
        }

        if (minLength < 1)
        {
            throw new ArgumentException($"Minimum peptide length must be at least 1; got {minLength}.", nameof(minLength));
        }

        if (maxLength < minLength)
        {
            throw new ArgumentException(
                $"Maximum peptide length ({maxLength}) must be ≥ minimum ({minLength}).", nameof(maxLength));
        }

        if (mutationPosition < 1 || mutationPosition > wildTypeProtein.Length)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutationPosition),
                mutationPosition,
                $"Mutation position must be in [1, {wildTypeProtein.Length}].");
        }

        int mutationIndex = mutationPosition - 1; // 0-based index of the substituted residue.
        char wildTypeResidue = wildTypeProtein[mutationIndex];
        if (wildTypeResidue == mutantResidue)
        {
            throw new ArgumentException(
                $"Mutant residue '{mutantResidue}' equals the wild-type residue at position {mutationPosition}; " +
                "a missense mutation requires a different amino acid.", nameof(mutantResidue));
        }

        // The mutant protein differs from the wild type only at the substituted residue.
        char[] mutantChars = wildTypeProtein.ToCharArray();
        mutantChars[mutationIndex] = mutantResidue;
        string mutantProtein = new(mutantChars);

        var peptides = new List<NeoantigenPeptide>();
        int proteinLength = wildTypeProtein.Length;

        // For each requested length k, a length-k window starts at 0-based index s and covers [s, s+k-1].
        // It spans the mutation iff s ≤ mutationIndex ≤ s+k-1, i.e. s ∈ [mutationIndex-k+1, mutationIndex],
        // additionally clamped to the protein bounds s ∈ [0, proteinLength-k]. This is exactly the set of
        // 8–11-mers contained in the 21-mer ±10-flank window of pVACtools / ProGeo-neo.
        for (int k = minLength; k <= maxLength; k++)
        {
            if (k > proteinLength)
            {
                continue; // No window of this length fits in the protein.
            }

            int firstStart = Math.Max(0, mutationIndex - k + 1);
            int lastStart = Math.Min(mutationIndex, proteinLength - k);
            for (int start = firstStart; start <= lastStart; start++)
            {
                string mutantPeptide = mutantProtein.Substring(start, k);
                string wildTypePeptide = wildTypeProtein.Substring(start, k);
                int offset = mutationIndex - start; // 0-based offset of the mutation within the window.
                peptides.Add(new NeoantigenPeptide(k, start + 1, mutantPeptide, wildTypePeptide, offset));
            }
        }

        return peptides;
    }

    #endregion

    #region ClassifyMhcBinding

    /// <summary>
    /// MHC molecule class for peptide-binding classification. Class I and class II have different accepted
    /// peptide-length ranges and different %Rank cutoffs (Reynisson et al. 2020, <i>Nucleic Acids Res.</i>
    /// 48(W1):W449–W454).
    /// </summary>
    public enum MhcClass
    {
        /// <summary>MHC class I (HLA-A/B/C). Canonical neoantigen peptide length 8–11.</summary>
        ClassI,

        /// <summary>MHC class II (HLA-DR/DQ/DP). Peptide length 13–25.</summary>
        ClassII
    }

    /// <summary>
    /// Binding-strength category assigned to a peptide–MHC pair from a caller-supplied predicted affinity
    /// (IC50) or %Rank. Categories follow the IEDB / NetMHCpan strong-/weak-binder convention.
    /// </summary>
    public enum BindingStrength
    {
        /// <summary>Strong binder (IC50 &lt; 50 nM, or class I %Rank &lt; 0.5% / class II %Rank &lt; 2%).</summary>
        Strong,

        /// <summary>Weak (intermediate) binder (IC50 &lt; 500 nM, or class I %Rank &lt; 2% / class II %Rank &lt; 10%).</summary>
        Weak,

        /// <summary>Not a binder (above the weak-binder cutoff).</summary>
        NonBinder
    }

    /// <summary>
    /// IC50 (nM) below which a peptide–MHC pair is a strong (high-affinity) binder. Source: Sette et al.
    /// (1994), <i>J. Immunol.</i> 153(12):5586–5592 — "an affinity threshold of approximately 500 nM
    /// (preferably 50 nM or less) apparently determines the capacity" to elicit a CTL response; the IEDB
    /// states "Peptides with IC50 values &lt;50 nM are considered high affinity". Strict inequality.
    /// </summary>
    public const double StrongBinderIc50Nm = 50.0;

    /// <summary>
    /// IC50 (nM) below which a peptide–MHC pair is at least a weak (intermediate-affinity) binder. Source:
    /// IEDB — "&lt;500 nM intermediate affinity"; Sette et al. (1994) (≈500 nM threshold); corroborated by
    /// Roomp, Antes &amp; Lengauer (2010), <i>BMC Bioinformatics</i> 11:90 (500 nM binder demarcation). Strict
    /// inequality.
    /// </summary>
    public const double WeakBinderIc50Nm = 500.0;

    /// <summary>
    /// Class I %Rank below which a peptide is a strong binder. Source: Reynisson et al. (2020) — "by default,
    /// %Rank &lt; 0.5% and %Rank &lt; 2% thresholds are considered for detecting SBs and WBs for class I".
    /// Strict inequality.
    /// </summary>
    public const double ClassIStrongBinderRankPercent = 0.5;

    /// <summary>
    /// Class I %Rank below which a peptide is at least a weak binder. Source: Reynisson et al. (2020) (class I
    /// WB &lt; 2%). Strict inequality.
    /// </summary>
    public const double ClassIWeakBinderRankPercent = 2.0;

    /// <summary>
    /// Class II %Rank below which a peptide is a strong binder. Source: Reynisson et al. (2020) — "%Rank &lt; 2%
    /// and %Rank &lt; 10%, for SBs and WBs for class II". Strict inequality.
    /// </summary>
    public const double ClassIIStrongBinderRankPercent = 2.0;

    /// <summary>
    /// Class II %Rank below which a peptide is at least a weak binder. Source: Reynisson et al. (2020) (class
    /// II WB &lt; 10%). Strict inequality.
    /// </summary>
    public const double ClassIIWeakBinderRankPercent = 10.0;

    /// <summary>
    /// Minimum accepted peptide length for MHC class II. Source: IEDB MHC class II tool description —
    /// "Peptides binding to MHC class II molecules ... typically range between 13 and 25 amino acids long".
    /// </summary>
    public const int MhcClassIIMinPeptideLength = 13;

    /// <summary>
    /// Maximum accepted peptide length for MHC class II. Source: IEDB MHC class II tool description (13–25).
    /// </summary>
    public const int MhcClassIIMaxPeptideLength = 25;

    /// <summary>
    /// Largest finite %Rank value (a %Rank is a percentile and must lie in [0, 100]). Source: Reynisson et
    /// al. (2020) — %Rank is "the top X% scores from random natural peptides".
    /// </summary>
    private const double MaxRankPercent = 100.0;

    /// <summary>
    /// Classifies a caller-supplied predicted peptide–MHC binding affinity (IC50 in nanomolar) into
    /// <see cref="BindingStrength.Strong"/> (IC50 &lt; 50 nM), <see cref="BindingStrength.Weak"/>
    /// (IC50 &lt; 500 nM), or <see cref="BindingStrength.NonBinder"/> (IC50 ≥ 500 nM). The cutoffs are the
    /// standard IEDB / NetMHCpan convention (Sette et al. 1994; IEDB) and the boundaries are strict
    /// inequalities, so 50 nM is weak (not strong) and 500 nM is a non-binder (not weak).
    /// <para>
    /// This method does NOT predict the IC50 — that requires a trained MHC-binding model (e.g. NetMHCpan) and
    /// is caller-supplied / out of scope (ONCO-MHC-001). It only classifies a supplied value.
    /// </para>
    /// </summary>
    /// <param name="ic50Nm">Predicted half-maximal inhibitory concentration in nM (must be &gt; 0).</param>
    /// <returns>The binding-strength category.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="ic50Nm"/> is not a finite value greater than 0 (IC50 is a positive concentration).
    /// </exception>
    public static BindingStrength ClassifyBindingAffinity(double ic50Nm)
    {
        if (double.IsNaN(ic50Nm) || double.IsInfinity(ic50Nm) || ic50Nm <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(ic50Nm), ic50Nm, "IC50 must be a finite concentration greater than 0 nM.");
        }

        if (ic50Nm < StrongBinderIc50Nm)
        {
            return BindingStrength.Strong;
        }

        return ic50Nm < WeakBinderIc50Nm ? BindingStrength.Weak : BindingStrength.NonBinder;
    }

    /// <summary>
    /// Classifies a caller-supplied predicted %Rank into <see cref="BindingStrength"/> using the NetMHCpan-4.1
    /// default cutoffs (Reynisson et al. 2020): class I — strong &lt; 0.5%, weak &lt; 2%; class II — strong
    /// &lt; 2%, weak &lt; 10%. The boundaries are strict inequalities (a value exactly at a cutoff falls into
    /// the weaker category).
    /// <para>
    /// This method does NOT predict the %Rank — the trained model is caller-supplied / out of scope.
    /// </para>
    /// </summary>
    /// <param name="percentRank">Predicted %Rank as a percentile in [0, 100].</param>
    /// <param name="mhcClass">MHC class selecting the cutoff set.</param>
    /// <returns>The binding-strength category.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="percentRank"/> is NaN or outside [0, 100] (a %Rank is a percentile).
    /// </exception>
    public static BindingStrength ClassifyBindingRank(double percentRank, MhcClass mhcClass)
    {
        if (double.IsNaN(percentRank) || percentRank < 0.0 || percentRank > MaxRankPercent)
        {
            throw new ArgumentOutOfRangeException(
                nameof(percentRank), percentRank, "%Rank must be a percentile in [0, 100].");
        }

        double strongCutoff = mhcClass == MhcClass.ClassI
            ? ClassIStrongBinderRankPercent
            : ClassIIStrongBinderRankPercent;
        double weakCutoff = mhcClass == MhcClass.ClassI
            ? ClassIWeakBinderRankPercent
            : ClassIIWeakBinderRankPercent;

        if (percentRank < strongCutoff)
        {
            return BindingStrength.Strong;
        }

        return percentRank < weakCutoff ? BindingStrength.Weak : BindingStrength.NonBinder;
    }

    /// <summary>
    /// Determines whether <paramref name="length"/> is a valid presented-peptide length for the given MHC
    /// class: class I 8–11 (canonical neoantigen range; Reynisson et al. 2020 gives 8–14 with default 8–11,
    /// matching <see cref="MhcClassIMinPeptideLength"/>/<see cref="MhcClassIMaxPeptideLength"/>), class II
    /// 13–25 (IEDB MHC class II tool description). Both bounds are inclusive.
    /// </summary>
    /// <param name="length">Peptide length (residue count).</param>
    /// <param name="mhcClass">MHC class selecting the accepted length range.</param>
    /// <returns><see langword="true"/> iff <paramref name="length"/> is within the class's accepted range.</returns>
    public static bool IsValidPeptideLength(int length, MhcClass mhcClass)
    {
        return mhcClass == MhcClass.ClassI
            ? length >= MhcClassIMinPeptideLength && length <= MhcClassIMaxPeptideLength
            : length >= MhcClassIIMinPeptideLength && length <= MhcClassIIMaxPeptideLength;
    }

    /// <summary>
    /// Classifies a candidate peptide–MHC pair end-to-end: a peptide whose length is not valid for the MHC
    /// class is not a presentable candidate and is classified <see cref="BindingStrength.NonBinder"/>
    /// regardless of affinity; otherwise the supplied IC50 is classified by
    /// <see cref="ClassifyBindingAffinity(double)"/>. Thin convenience wrapper over the length gate
    /// (<see cref="IsValidPeptideLength(int, MhcClass)"/>) and the affinity classifier.
    /// </summary>
    /// <param name="peptideLength">Peptide length (residue count).</param>
    /// <param name="ic50Nm">Caller-supplied predicted IC50 in nM (must be &gt; 0).</param>
    /// <param name="mhcClass">MHC class.</param>
    /// <returns>The binding-strength category, or <see cref="BindingStrength.NonBinder"/> for invalid length.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ic50Nm"/> is not finite and &gt; 0.</exception>
    public static BindingStrength ClassifyMhcBinding(int peptideLength, double ic50Nm, MhcClass mhcClass)
    {
        if (!IsValidPeptideLength(peptideLength, mhcClass))
        {
            return BindingStrength.NonBinder;
        }

        return ClassifyBindingAffinity(ic50Nm);
    }

    #endregion

    #region ctDNA analysis (ONCO-CTDNA-001)

    /// <summary>
    /// Computes the probability of detecting at least one circulating-tumor-DNA (ctDNA) molecule in a
    /// plasma cfDNA sample under the Poisson sampling model. With <paramref name="genomeEquivalents"/> n
    /// sequenced haploid genome equivalents, mutant allele fraction <paramref name="mutantAlleleFraction"/> d,
    /// and <paramref name="reporterCount"/> k independent tumour reporters, the expected number of mutant
    /// molecules is λ = n·d·k and the detection probability is
    /// <code>
    /// p = 1 − e^(−n·d·k)
    /// </code>
    /// Source: US Patent 11,085,084 B2 ("the probability of observing a single tumor reporter in cfDNA
    /// follows a Poisson distribution with mean λ = n × d … x = 1 − e^(−nd) … for k independent reporters
    /// p = 1 − e^(−ndk)"); Avanzini et al. (2020), <i>Science Advances</i> 6(50):eabc4308.
    /// </summary>
    /// <param name="genomeEquivalents">Number of sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1] (the ctDNA detection-limit fraction).</param>
    /// <param name="reporterCount">Number of independent tumour reporters k (≥ 1). Default 1.</param>
    /// <returns>Detection probability p = 1 − e^(−n·d·k) ∈ [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="genomeEquivalents"/> &lt; 0, <paramref name="mutantAlleleFraction"/> ∉ [0, 1] or not finite,
    /// or <paramref name="reporterCount"/> &lt; 1.
    /// </exception>
    public static double CtDnaDetectionProbability(
        int genomeEquivalents, double mutantAlleleFraction, int reporterCount = 1)
    {
        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        if (double.IsNaN(mutantAlleleFraction) || mutantAlleleFraction < 0.0 || mutantAlleleFraction > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutantAlleleFraction), mutantAlleleFraction, "Mutant allele fraction (d) must be in [0, 1].");
        }

        if (reporterCount < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(reporterCount), reporterCount, "Reporter count (k) must be at least 1.");
        }

        // λ = n·d·k; p = 1 − e^(−λ). λ = 0 ⇒ p = 0; large λ ⇒ p → 1 (never exceeds 1).
        double lambda = (double)genomeEquivalents * mutantAlleleFraction * reporterCount;
        return 1.0 - Math.Exp(-lambda); // p = 1 − e^(−λ).
    }

    /// <summary>
    /// Expected number of mutant molecules λ = n·d·k for a ctDNA sample. Source: US Patent 11,085,084 B2 /
    /// Avanzini et al. (2020) — Poisson mean λ = n × d (× k reporters); corroborated by the worked example of
    /// Pessoa et al. (2023): n = 15,000 genome equivalents at d = 0.001 (0.1% VAF) ⇒ λ = 15 mutant molecules.
    /// </summary>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1].</param>
    /// <param name="reporterCount">Independent tumour reporters k (≥ 1). Default 1.</param>
    /// <returns>Expected mutant-molecule count λ = n·d·k.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Same domain limits as <see cref="CtDnaDetectionProbability"/>.</exception>
    public static double ExpectedMutantMolecules(
        int genomeEquivalents, double mutantAlleleFraction, int reporterCount = 1)
    {
        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        if (double.IsNaN(mutantAlleleFraction) || mutantAlleleFraction < 0.0 || mutantAlleleFraction > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutantAlleleFraction), mutantAlleleFraction, "Mutant allele fraction (d) must be in [0, 1].");
        }

        if (reporterCount < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(reporterCount), reporterCount, "Reporter count (k) must be at least 1.");
        }

        return (double)genomeEquivalents * mutantAlleleFraction * reporterCount;
    }

    /// <summary>
    /// Decides whether ctDNA is detectable: at least one mutant molecule must be expected (λ = n·d·k ≥ 1)
    /// AND the Poisson detection probability must reach <paramref name="minDetectionProbability"/>. Source:
    /// the Poisson detection model (US Patent 11,085,084 B2; Avanzini et al. 2020) gives the probability;
    /// Newman et al. (2014), <i>Nat. Med.</i> 20(5):548–554 establish a validated detection range
    /// (0.025%–10% allele fraction) at the conventional 95% operating point. The λ ≥ 1 requirement enforces
    /// the physical limit that no mutant molecule can be observed when fewer than one is expected.
    /// </summary>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1].</param>
    /// <param name="reporterCount">Independent tumour reporters k (≥ 1). Default 1.</param>
    /// <param name="minDetectionProbability">Minimum p to call detected; default <see cref="DefaultCtDnaDetectionProbability"/> (0.95).</param>
    /// <returns><c>true</c> if λ ≥ 1 and p ≥ <paramref name="minDetectionProbability"/>; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Domain limits of <see cref="CtDnaDetectionProbability"/>, or <paramref name="minDetectionProbability"/> ∉ (0, 1].
    /// </exception>
    public static bool IsCtDnaDetected(
        int genomeEquivalents,
        double mutantAlleleFraction,
        int reporterCount = 1,
        double minDetectionProbability = DefaultCtDnaDetectionProbability)
    {
        if (double.IsNaN(minDetectionProbability) || minDetectionProbability <= 0.0 || minDetectionProbability > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minDetectionProbability), minDetectionProbability,
                "Minimum detection probability must be in (0, 1].");
        }

        double lambda = ExpectedMutantMolecules(genomeEquivalents, mutantAlleleFraction, reporterCount);
        if (lambda < MinExpectedMutantMolecules)
        {
            return false;
        }

        double probability = 1.0 - Math.Exp(-lambda);
        return probability >= minDetectionProbability;
    }

    /// <summary>
    /// Estimates the ctDNA tumour fraction from clonal heterozygous somatic SNVs observed in plasma at
    /// copy-neutral diploid loci. For such a variant the expected VAF is half the tumour-derived fraction
    /// (v = TF/2), so tumour fraction = 2 · (mean VAF). The mean is taken over the supplied variants'
    /// plasma VAFs (alt / total reads). Source: Antonello et al. (2024), CNAqc, <i>Genome Biology</i> 25:38
    /// (m = 1, n_tot = 2 special case of v = m·π / [2(1−π) + π·n_tot] gives v = π/2). The result is clamped
    /// to [0, 1] since a fraction cannot exceed 1.
    /// </summary>
    /// <param name="variants">Clonal heterozygous somatic SNVs observed in plasma (each VAF ≤ 0.5).</param>
    /// <returns>Estimated ctDNA tumour fraction ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (tumour fraction undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">A variant has invalid read counts or a VAF &gt; 0.5.</exception>
    public static double CalculateTumorFraction(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        double vafSum = 0.0;
        int count = 0;
        foreach (VariantObservation variant in variants)
        {
            double vaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
            if (vaf > 0.5)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(variants), vaf,
                    "A clonal heterozygous SNV cannot have VAF > 0.5 under the diploid model; locus is not copy-neutral diploid heterozygous.");
            }

            vafSum += vaf;
            count++;
        }

        if (count == 0)
        {
            throw new ArgumentException("Cannot estimate tumour fraction from an empty variant set.", nameof(variants));
        }

        double meanVaf = vafSum / count;
        double tumorFraction = TumorFractionFromVafFactor * meanVaf;
        return Math.Min(tumorFraction, 1.0); // a fraction cannot exceed 1.
    }

    /// <summary>
    /// Computes the mean plasma variant allele fraction across ctDNA reporters: the arithmetic mean of each
    /// variant's VAF = alt reads / total reads. Source: Newman et al. (2014), <i>Nat. Med.</i> 20(5):548–554 —
    /// ctDNA level is summarized as the fraction of mutant molecules across SNV/indel reporters.
    /// </summary>
    /// <param name="variants">Observed ctDNA reporters with plasma read counts.</param>
    /// <returns>Mean variant allele fraction ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (mean undefined).</exception>
    public static double CalculateMeanVaf(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        double vafSum = 0.0;
        int count = 0;
        foreach (VariantObservation variant in variants)
        {
            vafSum += CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
            count++;
        }

        if (count == 0)
        {
            throw new ArgumentException("Cannot compute mean VAF from an empty variant set.", nameof(variants));
        }

        return vafSum / count;
    }

    /// <summary>
    /// Converts a cfDNA input mass in nanograms to the number of haploid genome equivalents, using the
    /// standard conversion 3.3 pg per haploid genome (Devonshire et al. 2014, PMC4182654), i.e.
    /// GE = ng · 1000 / 3.3 (≈ 303 GE per ng; Alcaide et al. 2020, <i>Sci. Rep.</i> 10:12564). This gives the
    /// sampling depth n used in the Poisson detection model <see cref="CtDnaDetectionProbability"/>.
    /// </summary>
    /// <param name="cfDnaNanograms">cfDNA input mass in nanograms (≥ 0).</param>
    /// <returns>Haploid genome equivalents (a continuous count; not rounded).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="cfDnaNanograms"/> is negative or not finite.</exception>
    public static double HaploidGenomeEquivalents(double cfDnaNanograms)
    {
        if (double.IsNaN(cfDnaNanograms) || double.IsInfinity(cfDnaNanograms) || cfDnaNanograms < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(cfDnaNanograms), cfDnaNanograms, "cfDNA mass (ng) must be a finite non-negative value.");
        }

        return cfDnaNanograms * PicogramsPerNanogram / PicogramsPerHaploidGenome;
    }

    #endregion

    #region Minimal/Molecular Residual Disease (ONCO-MRD-001)

    /// <summary>
    /// Minimum number of tracked patient-specific variants that must be detected in plasma for a
    /// tumour-informed MRD assay to call the sample ctDNA-positive (MRD-positive). Source: Reinert et al.
    /// (2019), <i>JAMA Oncology</i> 5(8):1124–1131, as quoted in the tumour-informed ctDNA MRD review
    /// (PMC9265001, Table 1): "Plasma samples with at least two tumor-specific SNVs are defined as
    /// ctDNA-positive." (Signatera personalised 16-plex assay.)
    /// </summary>
    public const int DefaultMrdPositivityThreshold = 2;

    /// <summary>
    /// Minimum number of alternate (mutant) supporting reads at a tracked locus for that variant to be
    /// counted as detected in plasma. A variant with no supporting reads contributes no MRD signal.
    /// Source: Wan et al. (2020), <i>Sci. Transl. Med.</i> 12(548):eaaz8084 — per-locus signal must exceed
    /// background; the publishable universal read-count cutoff is error-model specific, so the minimal
    /// presence rule (≥ 1 alt read) is used by default and is configurable. The panel-level ≥2 calling rule
    /// (<see cref="DefaultMrdPositivityThreshold"/>) is independent of this value.
    /// </summary>
    public const int DefaultMrdMinSupportingReads = 1;

    /// <summary>
    /// A single patient-specific tumour marker tracked in plasma for MRD, with its plasma read evidence.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier of the tracked locus.</param>
    /// <param name="Position">1-based reference position of the tracked variant.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Tumour-specific alternate allele tracked at this locus.</param>
    /// <param name="PlasmaAltReads">Alternate (mutant) supporting reads observed in plasma at this locus (≥ 0).</param>
    /// <param name="PlasmaTotalReads">Total covering reads in plasma at this locus (≥ 0).</param>
    public readonly record struct TumorMarker(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        int PlasmaAltReads,
        int PlasmaTotalReads);

    /// <summary>Binary MRD call for a plasma timepoint.</summary>
    public enum MrdStatus
    {
        /// <summary>Fewer than the positivity threshold of tracked variants detected: no residual disease signal.</summary>
        Negative,

        /// <summary>At least the positivity threshold of tracked variants detected: molecular residual disease present.</summary>
        Positive
    }

    /// <summary>
    /// Result of a tumour-informed MRD assessment of one plasma sample against a patient-specific marker panel.
    /// </summary>
    /// <param name="Status">MRD-positive when <see cref="DetectedVariantCount"/> ≥ positivity threshold; else negative.</param>
    /// <param name="DetectedVariantCount">Number of tracked variants detected (alt reads ≥ minSupportingReads).</param>
    /// <param name="TrackedVariantCount">Total number of tracked variants in the panel.</param>
    /// <param name="IntegratedMutantAlleleFraction">
    /// Depth-weighted (read-pooled) mean plasma VAF across tracked loci: Σ alt / Σ total (INVAR IMAF). 0 when no reads.
    /// </param>
    /// <param name="DetectionProbability">
    /// Panel-level Poisson probability of detecting ≥ 1 ctDNA molecule, p = 1 − e^(−n·f·m), at the observed IMAF.
    /// </param>
    public readonly record struct MrdResult(
        MrdStatus Status,
        int DetectedVariantCount,
        int TrackedVariantCount,
        double IntegratedMutantAlleleFraction,
        double DetectionProbability);

    /// <summary>Per-timepoint MRD call in a longitudinal series.</summary>
    /// <param name="TimepointIndex">0-based index of the timepoint in input order.</param>
    /// <param name="Result">The MRD assessment for that timepoint's plasma sample.</param>
    public readonly record struct MrdTimepoint(int TimepointIndex, MrdResult Result);

    /// <summary>Longitudinal MRD tracking across ordered plasma timepoints.</summary>
    /// <param name="Timepoints">Per-timepoint MRD results, preserving input order.</param>
    /// <param name="FirstPositiveIndex">0-based index of the earliest MRD-positive timepoint, or −1 if none.</param>
    public readonly record struct MrdLongitudinalResult(
        IReadOnlyList<MrdTimepoint> Timepoints,
        int FirstPositiveIndex);

    /// <summary>
    /// Reports whether a tracked tumour marker is detected in plasma: its alternate (mutant) supporting-read
    /// count is at or above <paramref name="minSupportingReads"/>. Source: Wan et al. (2020) — a locus
    /// contributes ctDNA signal only when supported by mutant reads above background.
    /// </summary>
    /// <param name="marker">The tracked marker with its plasma read evidence.</param>
    /// <param name="minSupportingReads">Minimum alt reads to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <returns><c>true</c> if <c>PlasmaAltReads ≥ minSupportingReads</c>; otherwise <c>false</c>.</returns>
    public static bool IsVariantDetected(TumorMarker marker, int minSupportingReads = DefaultMrdMinSupportingReads)
    {
        if (minSupportingReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minSupportingReads), minSupportingReads, "Minimum supporting reads must be at least 1.");
        }

        return marker.PlasmaAltReads >= minSupportingReads;
    }

    /// <summary>
    /// Tumour-informed minimal/molecular residual disease (MRD) detection. Given a panel of patient-specific
    /// somatic markers (selected from the tumour and tracked in plasma), counts how many are detected and
    /// calls the sample MRD-positive when the number of detected variants reaches
    /// <paramref name="positivityThreshold"/> (default 2 of up to 16). Source: Reinert et al. (2019),
    /// <i>JAMA Oncology</i> 5(8):1124–1131 / Signatera analytical-validation white paper — "Plasma samples
    /// with at least two tumor-specific SNVs are defined as ctDNA-positive." The integrated mutant allele
    /// fraction (IMAF) is the depth-weighted mean plasma VAF across loci (Wan et al. 2020), and the
    /// panel-level Poisson detection probability p = 1 − e^(−n·f·m) reuses
    /// <see cref="CtDnaDetectionProbability"/> with m = number of tracked markers (Signatera white paper,
    /// Figure 2; Avanzini et al. 2020).
    /// </summary>
    /// <param name="tumorMarkers">Patient-specific tracked markers with plasma read evidence (non-empty).</param>
    /// <param name="positivityThreshold">Minimum detected variants to call positive (default <see cref="DefaultMrdPositivityThreshold"/>).</param>
    /// <param name="minSupportingReads">Minimum alt reads for a marker to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n for the panel Poisson p (≥ 0; default 0 ⇒ p = 0).</param>
    /// <returns>The MRD call with detected count, tracked count, IMAF, and panel detection probability.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="tumorMarkers"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="tumorMarkers"/> is empty (no panel to interrogate).</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold or <paramref name="genomeEquivalents"/> is out of range.</exception>
    public static MrdResult DetectMRD(
        IEnumerable<TumorMarker> tumorMarkers,
        int positivityThreshold = DefaultMrdPositivityThreshold,
        int minSupportingReads = DefaultMrdMinSupportingReads,
        int genomeEquivalents = 0)
    {
        ArgumentNullException.ThrowIfNull(tumorMarkers);

        if (positivityThreshold < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(positivityThreshold), positivityThreshold, "Positivity threshold must be at least 1.");
        }

        if (minSupportingReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minSupportingReads), minSupportingReads, "Minimum supporting reads must be at least 1.");
        }

        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        int trackedCount = 0;
        int detectedCount = 0;
        long altReadSum = 0;
        long totalReadSum = 0;
        foreach (TumorMarker marker in tumorMarkers)
        {
            trackedCount++;
            altReadSum += Math.Max(0, marker.PlasmaAltReads);
            totalReadSum += Math.Max(0, marker.PlasmaTotalReads);
            if (IsVariantDetected(marker, minSupportingReads))
            {
                detectedCount++;
            }
        }

        if (trackedCount == 0)
        {
            throw new ArgumentException(
                "Cannot assess MRD against an empty marker panel.", nameof(tumorMarkers));
        }

        // INVAR integrated mutant allele fraction: depth-weighted (read-pooled) plasma VAF across loci.
        double imaf = totalReadSum == 0 ? 0.0 : (double)altReadSum / totalReadSum;

        // Panel-level Poisson p = 1 - e^(-n*f*m) with f = IMAF, m = tracked markers (Signatera Fig 2).
        double detectionProbability = CtDnaDetectionProbability(genomeEquivalents, imaf, trackedCount);

        MrdStatus status = detectedCount >= positivityThreshold ? MrdStatus.Positive : MrdStatus.Negative;

        return new MrdResult(status, detectedCount, trackedCount, imaf, detectionProbability);
    }

    /// <summary>
    /// Longitudinal MRD tracking: applies <see cref="DetectMRD"/> to each timepoint's plasma marker panel,
    /// preserving input order, and reports the earliest MRD-positive timepoint (serial surveillance of
    /// residual disease, Reinert et al. 2019). Each timepoint is an independent plasma marker panel.
    /// </summary>
    /// <param name="timepoints">Ordered plasma marker panels, one per timepoint.</param>
    /// <param name="positivityThreshold">Minimum detected variants to call positive (default <see cref="DefaultMrdPositivityThreshold"/>).</param>
    /// <param name="minSupportingReads">Minimum alt reads for a marker to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n for the panel Poisson p (≥ 0; default 0).</param>
    /// <returns>Per-timepoint MRD results and the 0-based index of the first positive timepoint (−1 if none).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="timepoints"/> is null, or a timepoint panel is null.</exception>
    /// <exception cref="ArgumentException">A timepoint panel is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold or <paramref name="genomeEquivalents"/> is out of range.</exception>
    public static MrdLongitudinalResult TrackVariantsOverTime(
        IEnumerable<IEnumerable<TumorMarker>> timepoints,
        int positivityThreshold = DefaultMrdPositivityThreshold,
        int minSupportingReads = DefaultMrdMinSupportingReads,
        int genomeEquivalents = 0)
    {
        ArgumentNullException.ThrowIfNull(timepoints);

        var results = new List<MrdTimepoint>();
        int firstPositiveIndex = -1;
        int index = 0;
        foreach (IEnumerable<TumorMarker> panel in timepoints)
        {
            ArgumentNullException.ThrowIfNull(panel);
            MrdResult result = DetectMRD(panel, positivityThreshold, minSupportingReads, genomeEquivalents);
            if (firstPositiveIndex < 0 && result.Status == MrdStatus.Positive)
            {
                firstPositiveIndex = index;
            }

            results.Add(new MrdTimepoint(index, result));
            index++;
        }

        return new MrdLongitudinalResult(results, firstPositiveIndex);
    }

    // ---------------------------------------------------------------------------------------------
    // INVAR-style background-subtracted, tumour-AF-weighted ctDNA signal estimation (ONCO-MRD-001).
    //
    // Faithfully reproduces the core, caller-reproducible part of the INVAR pipeline (Wan et al. 2020,
    // Sci. Transl. Med. 12(548):eaaz8084; reference implementation INVAR2, nrlab-CRUK/INVAR2,
    // R/4_detection/generalisedLikelihoodRatioTest.R and R/shared/detectionFunctions.R):
    //  (a) per-locus / per-context BACKGROUND error subtraction (caller-supplied background AF per locus);
    //  (b) tumour-allele-fraction-weighted aggregation (IMAFv2) and a generalised-likelihood-ratio (GLRT)
    //      detection statistic whose mixture model weights each locus by its tumour AF vs background.
    // The fragment-length (size) weighting, outlier suppression and locus-noise filtering of the full
    // INVAR pipeline are NOT reproduced here; the caller supplies an already-cleaned background model.
    // ---------------------------------------------------------------------------------------------

    /// <summary>
    /// Initial value of the per-sample ctDNA fraction <c>p</c> for the INVAR EM maximum-likelihood
    /// estimator. Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R) — <c>initial_p = 0.01</c>.
    /// </summary>
    private const double InvarEmInitialP = 0.01;

    /// <summary>
    /// Number of expectation-maximisation iterations used to estimate the ctDNA fraction <c>p</c>.
    /// Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R) — <c>iterations = 200</c>.
    /// </summary>
    private const int InvarEmIterations = 200;

    /// <summary>
    /// A tracked patient-specific locus for INVAR-style ctDNA signal estimation: its plasma read evidence,
    /// the tumour allele fraction of the variant, and a caller-supplied background (non-reference) error rate.
    /// </summary>
    /// <param name="PlasmaAltReads">Mutant (alternate) supporting reads observed in plasma at this locus (≥ 0).</param>
    /// <param name="PlasmaTotalReads">Total covering reads in plasma at this locus (≥ 0).</param>
    /// <param name="TumorAlleleFraction">
    /// Tumour allele fraction of the tracked variant, <c>AF</c> in INVAR (0 &lt; AF ≤ 1). Loci with higher
    /// tumour AF carry more ctDNA signal and are weighted more strongly in the likelihood (Wan et al. 2020).
    /// </param>
    /// <param name="BackgroundErrorRate">
    /// Caller-supplied per-locus / per-trinucleotide-context background (non-reference) error rate <c>e</c>
    /// in INVAR (0 ≤ e &lt; 1), estimated from control plasma samples at the same loci. Subtracted from the
    /// observed plasma signal and used as the null read-error rate in the likelihood model.
    /// </param>
    public readonly record struct InvarLocus(
        int PlasmaAltReads,
        int PlasmaTotalReads,
        double TumorAlleleFraction,
        double BackgroundErrorRate);

    /// <summary>
    /// Result of an INVAR-style background-subtracted, tumour-AF-weighted ctDNA signal estimate for one sample.
    /// </summary>
    /// <param name="IntegratedMutantAlleleFractionV2">
    /// Background-subtracted, depth-weighted aggregate tumour fraction (INVAR2 IMAFv2): the depth-weighted mean
    /// over loci of <c>max(0, locusVAF − backgroundRate)</c>. ≥ 0; equals 0 when no locus exceeds background.
    /// </param>
    /// <param name="EstimatedTumorFraction">
    /// Maximum-likelihood ctDNA fraction <c>p̂</c> from the INVAR EM estimator under the AF-weighted mixture model.
    /// </param>
    /// <param name="LikelihoodRatio">
    /// Generalised-likelihood-ratio detection statistic: <c>logL(p̂) − logL(p = 0)</c> (per-locus mean scaled,
    /// as in INVAR2). Larger ⇒ stronger ctDNA evidence; ≈ 0 for a pure-background sample.
    /// </param>
    /// <param name="Detected">
    /// <c>true</c> when <see cref="LikelihoodRatio"/> reaches <paramref name="detectionThreshold"/> AND at least
    /// one mutant read is present; otherwise <c>false</c>.
    /// </param>
    /// <param name="LocusCount">Number of informative loci used (tumour AF &gt; 0).</param>
    public readonly record struct InvarSignalResult(
        double IntegratedMutantAlleleFractionV2,
        double EstimatedTumorFraction,
        double LikelihoodRatio,
        bool Detected,
        int LocusCount);

    /// <summary>
    /// Computes the INVAR2 background-subtracted, depth-weighted integrated mutant allele fraction (IMAFv2):
    /// the depth-weighted mean over loci of <c>max(0, plasmaVAF − backgroundRate)</c>. Source: INVAR2
    /// <c>calculateIMAFv2</c> (R/4_detection/generalisedLikelihoodRatioTest.R) — per-context
    /// <c>MEAN_AF.BS = pmax(0, MEAN_AF − BACKGROUND_AF)</c> then
    /// <c>weighted.mean(MEAN_AF.BS, TOTAL_DP)</c>.
    /// </summary>
    /// <param name="loci">The tracked informative loci (non-null). Loci with 0 total reads contribute 0 weight.</param>
    /// <returns>The background-subtracted depth-weighted tumour fraction (≥ 0); 0 when no covering reads.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    public static double IntegratedMutantAlleleFractionV2(IEnumerable<InvarLocus> loci)
    {
        ArgumentNullException.ThrowIfNull(loci);

        double weightedSum = 0.0;
        long totalDepth = 0;
        foreach (InvarLocus locus in loci)
        {
            int total = Math.Max(0, locus.PlasmaTotalReads);
            if (total == 0)
            {
                continue;
            }

            double vaf = (double)Math.Max(0, locus.PlasmaAltReads) / total;

            // Per-locus/per-context background subtraction: signal above background only (INVAR2 pmax(0, .)).
            double backgroundSubtracted = Math.Max(0.0, vaf - locus.BackgroundErrorRate);

            // Depth-weighted aggregation (INVAR2 weighted.mean(., TOTAL_DP)).
            weightedSum += backgroundSubtracted * total;
            totalDepth += total;
        }

        return totalDepth == 0 ? 0.0 : weightedSum / totalDepth;
    }

    /// <summary>
    /// INVAR-style estimate of residual ctDNA signal from a panel of tracked loci, each carrying plasma read
    /// evidence, a tumour allele fraction, and a caller-supplied per-locus background error rate. Performs
    /// (a) per-locus background subtraction, (b) tumour-AF-weighted aggregation (IMAFv2), (c) maximum-likelihood
    /// estimation of the ctDNA fraction <c>p̂</c> and (d) a generalised-likelihood-ratio detection statistic,
    /// faithfully following INVAR2 (Wan et al. 2020; nrlab-CRUK/INVAR2 detectionFunctions.R, no-size variant).
    ///
    /// <para>Mixture model (per read at a locus with tumour AF <c>AF</c>, background <c>e</c>, ctDNA fraction
    /// <c>p</c>): the probability a read is mutant is <c>q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)</c>; loci with
    /// higher <c>AF</c> and lower <c>e</c> contribute more signal, so the estimate is signal-to-noise weighted.</para>
    /// </summary>
    /// <param name="loci">Tracked informative loci with plasma evidence, tumour AF and background rate (non-empty).</param>
    /// <param name="detectionThreshold">
    /// Minimum generalised-likelihood-ratio statistic to call the sample ctDNA-positive (≥ 0; default 0 ⇒
    /// any positive evidence with a mutant read is detected). Larger values trade sensitivity for specificity.
    /// </param>
    /// <returns>The INVAR signal estimate: IMAFv2, ML ctDNA fraction, likelihood-ratio statistic and detection call.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="loci"/> has no informative locus (tumour AF &gt; 0).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="detectionThreshold"/> is negative, or a locus has an out-of-range tumour AF
    /// (must be 0 &lt; AF ≤ 1 to be informative) or background rate (must be 0 ≤ e &lt; 1).
    /// </exception>
    public static InvarSignalResult EstimateInvarSignal(
        IEnumerable<InvarLocus> loci,
        double detectionThreshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(loci);

        if (detectionThreshold < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(detectionThreshold), detectionThreshold, "Detection threshold cannot be negative.");
        }

        var materialised = loci as IReadOnlyCollection<InvarLocus> ?? loci.ToList();

        // IMAFv2 is computed over all covered loci (background-subtracted, depth-weighted).
        double imafV2 = IntegratedMutantAlleleFractionV2(materialised);

        // Build per-informative-locus vectors for the likelihood model (INVAR uses one row per molecule,
        // i.e. depth = total covering reads, M = mutant reads). Only loci with tumour AF > 0 are informative.
        var altReads = new List<double>();
        var totalReads = new List<double>();
        var tumorAf = new List<double>();
        var background = new List<double>();
        bool anyMutantRead = false;
        foreach (InvarLocus locus in materialised)
        {
            if (locus.TumorAlleleFraction < 0.0 || locus.TumorAlleleFraction > 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(loci), locus.TumorAlleleFraction, "Tumour allele fraction must be in [0, 1].");
            }

            if (locus.BackgroundErrorRate < 0.0 || locus.BackgroundErrorRate >= 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(loci), locus.BackgroundErrorRate, "Background error rate must be in [0, 1).");
            }

            // INVAR keeps only loci with tumour AF > 0 (filter(TUMOUR_AF > 0)).
            if (locus.TumorAlleleFraction <= 0.0)
            {
                continue;
            }

            int total = Math.Max(0, locus.PlasmaTotalReads);
            if (total == 0)
            {
                continue;
            }

            int alt = Math.Clamp(locus.PlasmaAltReads, 0, total);

            // INVAR guards a zero background by flooring it to one expected error in the locus depth
            // (BACKGROUND_AF = ifelse(BACKGROUND_AF > 0, BACKGROUND_AF, 1 / BACKGROUND_DP)), so log(e) is finite.
            double e = locus.BackgroundErrorRate > 0.0 ? locus.BackgroundErrorRate : 1.0 / total;

            altReads.Add(alt);
            totalReads.Add(total);
            tumorAf.Add(locus.TumorAlleleFraction);
            background.Add(e);
            if (alt > 0)
            {
                anyMutantRead = true;
            }
        }

        if (altReads.Count == 0)
        {
            throw new ArgumentException(
                "No informative locus (tumour AF > 0 with covering reads) to estimate INVAR signal.", nameof(loci));
        }

        double pMle = EstimateCtDnaFractionEm(altReads, totalReads, tumorAf, background);
        double nullLogLik = InvarLogLikelihood(altReads, totalReads, tumorAf, background, 0.0);
        double altLogLik = InvarLogLikelihood(altReads, totalReads, tumorAf, background, pMle);
        double likelihoodRatio = altLogLik - nullLogLik;

        bool detected = anyMutantRead && likelihoodRatio >= detectionThreshold;

        return new InvarSignalResult(imafV2, pMle, likelihoodRatio, detected, altReads.Count);
    }

    /// <summary>
    /// EM maximum-likelihood estimate of the ctDNA fraction <c>p</c> under the INVAR mixture model.
    /// Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R):
    /// <c>g = AF·(1−e) + (1−AF)·e</c>;
    /// E-step <c>Z0 = (1−g)·p / ((1−g)·p + (1−e)·(1−p))</c>, <c>Z1 = g·p / (g·p + e·(1−p))</c>;
    /// M-step <c>p = Σ(M·Z1 + (R−M)·Z0) / ΣR</c>.
    /// </summary>
    private static double EstimateCtDnaFractionEm(
        IReadOnlyList<double> m,
        IReadOnlyList<double> r,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e)
    {
        double p = InvarEmInitialP;
        for (int iter = 0; iter < InvarEmIterations; iter++)
        {
            double numerator = 0.0;
            double denominator = 0.0;
            for (int i = 0; i < m.Count; i++)
            {
                double g = (af[i] * (1.0 - e[i])) + ((1.0 - af[i]) * e[i]);
                double z0 = ((1.0 - g) * p) / (((1.0 - g) * p) + ((1.0 - e[i]) * (1.0 - p)));
                double z1 = (g * p) / ((g * p) + (e[i] * (1.0 - p)));
                numerator += (m[i] * z1) + ((r[i] - m[i]) * z0);
                denominator += r[i];
            }

            p = denominator == 0.0 ? 0.0 : numerator / denominator;
        }

        return p;
    }

    /// <summary>
    /// Per-locus-mean log-likelihood of the sample given a ctDNA fraction <c>p</c> under the INVAR mixture model.
    /// Source: INVAR2 <c>calc_log_likelihood</c> (R/shared/detectionFunctions.R):
    /// <c>q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)</c>;
    /// <c>logL = Σ[ lchoose(R, M) + M·log(q) + (R−M)·log(1−q) ] / length(R)</c>.
    /// </summary>
    private static double InvarLogLikelihood(
        IReadOnlyList<double> m,
        IReadOnlyList<double> r,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e,
        double p)
    {
        double sum = 0.0;
        for (int i = 0; i < m.Count; i++)
        {
            double q = (af[i] * (1.0 - e[i]) * p) + ((1.0 - af[i]) * e[i] * p) + (e[i] * (1.0 - p));

            // Clamp q strictly inside (0, 1) so the logs stay finite at the boundaries (p = 0, e = 0 ⇒ q = 0).
            q = Math.Clamp(q, double.Epsilon, 1.0 - double.Epsilon);

            double lchoose = LogChoose(r[i], m[i]);
            sum += lchoose + (m[i] * Math.Log(q)) + ((r[i] - m[i]) * Math.Log(1.0 - q));
        }

        return sum / m.Count;
    }

    /// <summary>
    /// Natural log of the binomial coefficient C(n, k) via log-gamma, matching R's <c>lchoose(n, k)</c>:
    /// <c>lgamma(n+1) − lgamma(k+1) − lgamma(n−k+1)</c>.
    /// </summary>
    private static double LogChoose(double n, double k)
    {
        if (k < 0.0 || k > n)
        {
            return double.NegativeInfinity;
        }

        return LogGamma(n + 1.0) - LogGamma(k + 1.0) - LogGamma(n - k + 1.0);
    }

    #endregion

    #region Clonal Hematopoiesis Filtering (ONCO-CHIP-001)

    /// <summary>
    /// Minimum variant allele fraction (VAF) for a driver-gene somatic mutation in blood to meet the
    /// clonal-hematopoiesis-of-indeterminate-potential (CHIP) definition. Source: Steensma et al. (2015),
    /// <i>Blood</i> 126(1):9–16 — "the mutant allele fraction must be ≥2% in the peripheral blood"
    /// (threshold is inclusive: ≥ 0.02).
    /// </summary>
    public const double ChipVafThreshold = 0.02;

    /// <summary>
    /// Default canonical CHIP driver-gene panel (HGNC symbols, upper-case): genes recurrently mutated in
    /// clonal hematopoiesis. Source: Steensma et al. (2015) Fig. 2A and Genovese et al. (2014),
    /// <i>NEJM</i> 371(26):2477–2487 ("Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately
    /// high numbers of somatic mutations"; JAK2 V617F and SF3B1 K700E recurrent; SRSF2/TP53 are established
    /// CHIP drivers). This is a labelled canonical set, NOT an invented value — callers may override it via
    /// the <c>chipGenes</c> parameter (the algorithm is gene-panel-agnostic, per Razavi et al. 2019).
    /// </summary>
    public static readonly IReadOnlyCollection<string> DefaultChipGenes = new[]
    {
        "DNMT3A", "TET2", "ASXL1", "TP53", "JAK2", "SF3B1", "SRSF2", "PPM1D"
    };

    /// <summary>
    /// A variant observed in plasma cell-free DNA (cfDNA) for CHIP analysis, carrying its locus, gene
    /// symbol, and plasma VAF, plus optional matched white-blood-cell (WBC) alt-read evidence used by
    /// <see cref="FilterCHIP(IEnumerable{ChipVariant}, IEnumerable{ChipVariant}, IReadOnlyCollection{string}?, double, int)"/>.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier of the locus.</param>
    /// <param name="Position">1-based reference position.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Alternate (mutant) allele.</param>
    /// <param name="Gene">HGNC gene symbol the variant falls in (case-insensitive on comparison).</param>
    /// <param name="Vaf">Plasma variant allele fraction in [0, 1].</param>
    /// <param name="AltReads">Alternate (mutant) supporting reads at this locus (≥ 0); used as WBC evidence.</param>
    public readonly record struct ChipVariant(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        string Gene,
        double Vaf,
        int AltReads = 0);

    /// <summary>
    /// Reports whether a gene symbol belongs to the CHIP driver-gene panel (case-insensitive). Source:
    /// Steensma et al. (2015) / Genovese et al. (2014) canonical driver genes; the panel is caller-supplied
    /// when <paramref name="chipGenes"/> is provided, otherwise <see cref="DefaultChipGenes"/>.
    /// </summary>
    /// <param name="gene">Gene symbol to test (null/empty ⇒ <c>false</c>).</param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <returns><c>true</c> when <paramref name="gene"/> is in the panel; otherwise <c>false</c>.</returns>
    public static bool IsCanonicalChipGene(string? gene, IReadOnlyCollection<string>? chipGenes = null)
    {
        if (string.IsNullOrEmpty(gene))
        {
            return false;
        }

        IReadOnlyCollection<string> panel = chipGenes ?? DefaultChipGenes;
        foreach (string g in panel)
        {
            if (string.Equals(g, gene, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Identifies candidate clonal-hematopoiesis (CHIP) variants by the gene + VAF heuristic: a variant is
    /// flagged CHIP when its gene is in the CHIP driver panel AND its plasma VAF is at or above
    /// <paramref name="minVaf"/> (default <see cref="ChipVafThreshold"/> = 0.02). Source: Steensma et al.
    /// (2015) — a somatic mutation in a gene recurrently mutated in hematologic malignancies at VAF ≥ 2%.
    /// This is a candidate flag; the definitive tumour-vs-CH origin test is matched-WBC subtraction
    /// (<see cref="FilterCHIP(IEnumerable{ChipVariant}, IEnumerable{ChipVariant}, IReadOnlyCollection{string}?, double, int)"/>,
    /// Razavi et al. 2019).
    /// </summary>
    /// <param name="variants">cfDNA variants to screen (non-null).</param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <param name="minVaf">Minimum VAF to meet the CHIP definition (default <see cref="ChipVafThreshold"/>); must be in (0, 1].</param>
    /// <returns>The subset of <paramref name="variants"/> flagged as candidate CHIP, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minVaf"/> ∉ (0, 1].</exception>
    public static IReadOnlyList<ChipVariant> IdentifyCHIPVariants(
        IEnumerable<ChipVariant> variants,
        IReadOnlyCollection<string>? chipGenes = null,
        double minVaf = ChipVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);

        if (!(minVaf > 0.0) || minVaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minVaf), minVaf, "Minimum CHIP VAF must be in the interval (0, 1].");
        }

        var result = new List<ChipVariant>();
        foreach (ChipVariant variant in variants)
        {
            if (IsCanonicalChipGene(variant.Gene, chipGenes) && variant.Vaf >= minVaf)
            {
                result.Add(variant);
            }
        }

        return result;
    }

    /// <summary>
    /// Removes clonal-hematopoiesis (CHIP) confounder variants from a plasma cfDNA call set so that the
    /// retained variants are candidate tumour-derived variants. A cfDNA variant is removed when EITHER
    /// (a) it is also detected in the matched white-blood-cell (WBC) sample at the same locus — the
    /// definitive matched-WBC origin test (Razavi et al. 2019, <i>Nat Med</i> 25:1928–1937: matched
    /// cfDNA–WBC sequencing assigns variant origin tumour-vs-CH) — OR (b) it meets the gene + VAF CHIP
    /// heuristic (<see cref="IdentifyCHIPVariants"/>, Steensma et al. 2015). Rule (a) applies regardless of
    /// gene; rule (b) is the fallback when no matched-WBC evidence exists. Output is a subset of the input
    /// in input order.
    /// </summary>
    /// <param name="variants">Plasma cfDNA variants to filter (non-null).</param>
    /// <param name="whiteBloodCellVariants">
    /// Matched WBC variants; a cfDNA variant sharing a locus (chromosome, position, ref, alt) with a WBC
    /// variant carrying ≥ <paramref name="minWbcAltReads"/> alt reads is treated as WBC/CH-derived.
    /// </param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <param name="minVaf">Minimum VAF for the gene+VAF CHIP heuristic (default <see cref="ChipVafThreshold"/>); in (0, 1].</param>
    /// <param name="minWbcAltReads">Minimum alt reads in matched WBC to count the locus as present (default 1; Wan et al. 2020).</param>
    /// <returns>cfDNA variants retained as candidate tumour-derived, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> or <paramref name="whiteBloodCellVariants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minVaf"/> ∉ (0, 1], or <paramref name="minWbcAltReads"/> &lt; 1.</exception>
    public static IReadOnlyList<ChipVariant> FilterCHIP(
        IEnumerable<ChipVariant> variants,
        IEnumerable<ChipVariant> whiteBloodCellVariants,
        IReadOnlyCollection<string>? chipGenes = null,
        double minVaf = ChipVafThreshold,
        int minWbcAltReads = DefaultMrdMinSupportingReads)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ArgumentNullException.ThrowIfNull(whiteBloodCellVariants);

        if (!(minVaf > 0.0) || minVaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minVaf), minVaf, "Minimum CHIP VAF must be in the interval (0, 1].");
        }

        if (minWbcAltReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minWbcAltReads), minWbcAltReads, "Minimum WBC alt reads must be at least 1.");
        }

        // Build the set of loci present in the matched WBC (alt reads >= cutoff) for O(1) lookup.
        var wbcLoci = new HashSet<(string, int, string, string)>();
        foreach (ChipVariant wbc in whiteBloodCellVariants)
        {
            if (wbc.AltReads >= minWbcAltReads)
            {
                wbcLoci.Add(LocusKey(wbc));
            }
        }

        var retained = new List<ChipVariant>();
        foreach (ChipVariant variant in variants)
        {
            bool inMatchedWbc = wbcLoci.Contains(LocusKey(variant));
            bool meetsChipHeuristic = IsCanonicalChipGene(variant.Gene, chipGenes) && variant.Vaf >= minVaf;
            if (!inMatchedWbc && !meetsChipHeuristic)
            {
                retained.Add(variant);
            }
        }

        return retained;
    }

    /// <summary>Locus identity key (chromosome, 1-based position, ref allele, alt allele) for matched-WBC subtraction.</summary>
    private static (string, int, string, string) LocusKey(ChipVariant v) =>
        (v.Chromosome, v.Position, v.ReferenceAllele, v.AlternateAllele);

    #endregion

    #region Tumor Phylogeny Reconstruction (ONCO-PHYLO-001)

    /// <summary>
    /// Cancer cell fraction (CCF) of the synthetic root (normal / germline) clone in every sample: it is, by
    /// definition, present in 100% of cells. Source: Popic et al. (2015), <i>Genome Biology</i> 16:91 — the
    /// lineage tree is a spanning tree rooted at the population that contains all observed clones.
    /// </summary>
    private const double RootCcf = 1.0;

    /// <summary>
    /// Default noise margin ε for the lineage-precedence (Eq. 2) and sum-rule (Eq. 5) inequalities. The cited
    /// sources relax both rules by a configurable ε (Popic et al. 2015 ϵ; Zheng et al. 2022 ε₁=0.1, ε₂=0.2). Because
    /// this unit consumes already-clustered CCF point estimates (clustering and its noise model are ONCO-CCF-001),
    /// the default is the strict ε = 0; callers may pass a positive tolerance to reproduce the source defaults.
    /// </summary>
    public const double DefaultPhylogenyTolerance = 0.0;

    /// <summary>
    /// One CCF cluster (a candidate clone/subclone) used as input to phylogeny reconstruction. Each cluster carries
    /// its cancer cell fraction in each sequenced sample. CCF clustering itself is out of scope (ONCO-CCF-001).
    /// </summary>
    /// <param name="Id">Caller-assigned cluster identifier (e.g. a subclone label). Used for deterministic tie-breaking.</param>
    /// <param name="CcfPerSample">Cancer cell fraction in each sample, each value in [0, 1]; all clusters must share length.</param>
    public readonly record struct CcfCluster(int Id, IReadOnlyList<double> CcfPerSample);

    /// <summary>A single parent → child edge of the reconstructed clonal tree.</summary>
    /// <param name="ParentId">Id of the ancestral cluster, or <see cref="ClonalPhylogeny.RootId"/> for the normal root.</param>
    /// <param name="ChildId">Id of the descendant cluster.</param>
    public readonly record struct ClonalEdge(int ParentId, int ChildId);

    /// <summary>
    /// The reconstructed rooted clonal tree: a synthetic normal root (<see cref="RootId"/>, CCF = 1 in every sample)
    /// plus one node per input CCF cluster, connected by parent→child <see cref="Edges"/>.
    /// </summary>
    /// <param name="RootId">Identifier of the synthetic normal/germline root node.</param>
    /// <param name="Clusters">Input clusters keyed by id, in input order.</param>
    /// <param name="Edges">Tree edges (parent→child), one per non-root cluster.</param>
    /// <param name="SampleCount">Number of samples per cluster.</param>
    public readonly record struct ClonalPhylogeny(
        int RootId,
        IReadOnlyList<CcfCluster> Clusters,
        IReadOnlyList<ClonalEdge> Edges,
        int SampleCount)
    {
        /// <summary>Returns the parent id of <paramref name="clusterId"/>, or null if it is the root or absent.</summary>
        public int? ParentOf(int clusterId)
        {
            foreach (ClonalEdge e in Edges)
            {
                if (e.ChildId == clusterId)
                {
                    return e.ParentId;
                }
            }

            return null;
        }

        /// <summary>Returns the ids of the direct children of <paramref name="clusterId"/>, in input order.</summary>
        public IReadOnlyList<int> ChildrenOf(int clusterId)
        {
            var children = new List<int>();
            foreach (ClonalEdge e in Edges)
            {
                if (e.ParentId == clusterId)
                {
                    children.Add(e.ChildId);
                }
            }

            return children;
        }
    }

    /// <summary>
    /// Reconstructs a rooted clonal (tumor) phylogeny from per-sample CCF clusters, applying the two lineage
    /// constraints from the multi-sample perfect-phylogeny model:
    /// <list type="number">
    /// <item><description><b>Lineage precedence (ancestor ≥ descendant), Eq. 2:</b> an edge u→v is admissible only if,
    /// for every sample i, <c>u.CCF[i] ≥ v.CCF[i] − ε</c> and (presence) <c>u.CCF[i] = 0 ⇒ v.CCF[i] = 0</c>. Source:
    /// Popic et al. (2015), <i>Genome Biology</i> 16:91, Eq. 2; Zheng et al. (2022) PICTograph, <i>Bioinformatics</i>
    /// 38(15):3677–3683 — "the CCF of any mutation cannot exceed the CCF of its ancestor".</description></item>
    /// <item><description><b>Sum rule, Eq. 5:</b> for every node u and every sample i, the children CCFs may not exceed
    /// the parent: <c>Σ_children v.CCF[i] ≤ u.CCF[i] + ε</c>. Source: Popic et al. (2015) Eq. 5; Zheng et al. (2022) —
    /// "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants".</description></item>
    /// </list>
    /// The constraints leave a set of valid trees; to return a single deterministic tree this method attaches each
    /// cluster (processed in descending order of total CCF) to its <i>deepest valid ancestor</i> — the admissible
    /// parent with the smallest total CCF whose remaining per-sample sum-rule budget still admits the child — with
    /// ties broken by ascending cluster id (Evidence Assumption 1).
    /// </summary>
    /// <param name="clusters">CCF clusters to place; each cluster's <see cref="CcfCluster.CcfPerSample"/> must have the same length.</param>
    /// <param name="tolerance">Noise margin ε for both inequalities; default <see cref="DefaultPhylogenyTolerance"/> (0).</param>
    /// <returns>The reconstructed <see cref="ClonalPhylogeny"/> rooted at a synthetic normal node.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="clusters"/> or any cluster's CCF list is null.</exception>
    /// <exception cref="ArgumentException">CCF lists differ in length, are empty, or contain NaN / out-of-[0,1] values; or two clusters share an id.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tolerance"/> is negative or NaN.</exception>
    public static ClonalPhylogeny ReconstructPhylogeny(
        IReadOnlyList<CcfCluster> clusters,
        double tolerance = DefaultPhylogenyTolerance)
    {
        ArgumentNullException.ThrowIfNull(clusters);
        if (double.IsNaN(tolerance) || tolerance < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tolerance), tolerance, "Phylogeny tolerance ε must be a non-negative number.");
        }

        int rootId = RootIdFor(clusters);
        if (clusters.Count == 0)
        {
            // Empty cohort: tree is the root alone, no clusters, no edges. One synthetic sample of CCF 1.
            return new ClonalPhylogeny(rootId, Array.Empty<CcfCluster>(), Array.Empty<ClonalEdge>(), 1);
        }

        int sampleCount = ValidateAndGetSampleCount(clusters);

        // Synthetic root: present in 100% of cells in every sample, so it is a valid ancestor of every cluster and
        // its sum-rule budget is the full RootCcf. The root participates in placement as node id = rootId.
        double[] rootCcf = new double[sampleCount];
        Array.Fill(rootCcf, RootCcf);

        // Per-node remaining sum-rule budget = node CCF minus the sum of CCFs of children already attached (per sample).
        var remainingBudget = new Dictionary<int, double[]>(clusters.Count + 1);
        remainingBudget[rootId] = (double[])rootCcf.Clone();
        var ccfById = new Dictionary<int, double[]>(clusters.Count + 1) { [rootId] = rootCcf };
        foreach (CcfCluster c in clusters)
        {
            double[] ccf = c.CcfPerSample.ToArray();
            ccfById[c.Id] = ccf;
            remainingBudget[c.Id] = (double[])ccf.Clone();
        }

        // Process clusters most-clonal-first (descending total CCF) so that ancestors are placed before descendants;
        // ties broken by ascending id for determinism (Evidence Assumption 1).
        var ordered = clusters
            .OrderByDescending(c => TotalCcf(ccfById[c.Id]))
            .ThenBy(c => c.Id)
            .ToList();

        var edges = new List<ClonalEdge>(clusters.Count);
        foreach (CcfCluster child in ordered)
        {
            double[] childCcf = ccfById[child.Id];

            // Candidate parents: the root plus every already-placed cluster. Choose the deepest valid ancestor =
            // the candidate with the smallest total CCF that (a) satisfies lineage precedence and (b) still has
            // per-sample budget for this child. Ties broken by ascending id.
            int bestParent = rootId;
            double bestParentTotal = double.PositiveInfinity;
            bool found = false;
            foreach (int candidateId in EnumerateCandidates(rootId, edges))
            {
                if (candidateId == child.Id)
                {
                    continue;
                }

                double[] parentCcf = ccfById[candidateId];
                if (!SatisfiesLineagePrecedence(parentCcf, childCcf, tolerance))
                {
                    continue;
                }

                if (!FitsSumRule(remainingBudget[candidateId], childCcf, tolerance))
                {
                    continue;
                }

                double candidateTotal = TotalCcf(parentCcf);
                if (!found
                    || candidateTotal < bestParentTotal
                    || (candidateTotal == bestParentTotal && candidateId < bestParent))
                {
                    bestParent = candidateId;
                    bestParentTotal = candidateTotal;
                    found = true;
                }
            }

            // The root always satisfies lineage precedence (CCF=1 ≥ any child) and, per the sum rule, can only be
            // exhausted if children already consumed its full budget; in that degenerate case we still attach to the
            // root because every cluster must have a parent (Popic 2015: spanning tree). 'found' is therefore the
            // generic path; the explicit root fallback preserves the spanning-tree invariant.
            int parentId = found ? bestParent : rootId;
            edges.Add(new ClonalEdge(parentId, child.Id));

            // Debit the chosen parent's per-sample budget by the child's CCF.
            double[] budget = remainingBudget[parentId];
            for (int i = 0; i < sampleCount; i++)
            {
                budget[i] -= childCcf[i];
            }
        }

        var clusterList = clusters.ToArray();
        return new ClonalPhylogeny(rootId, clusterList, edges, sampleCount);
    }

    /// <summary>
    /// Returns the ids of clusters on the <b>trunk</b> of the phylogeny — the clonal mutations shared by every
    /// tumor cell. The trunk is the path from the root down to (but excluding) the first branch point: a maximal
    /// chain of single-child nodes starting at the root's unique child. Source: Popic et al. (2015) — the trunk
    /// holds the mutations of the common predecessor present across all samples.
    /// </summary>
    /// <param name="phylogeny">A phylogeny produced by <see cref="ReconstructPhylogeny"/>.</param>
    /// <returns>Trunk cluster ids, ordered from the root downward; empty if the tree has no clusters.</returns>
    public static IReadOnlyList<int> IdentifyTrunkMutations(ClonalPhylogeny phylogeny)
    {
        var trunk = new List<int>();
        if (phylogeny.Clusters is null || phylogeny.Clusters.Count == 0)
        {
            return trunk;
        }

        // Walk down from the root while each node has exactly one child (no branching yet) and the root itself has
        // exactly one child. The first node with ≠1 children is the branch point; nodes below it are subclonal.
        int current = phylogeny.RootId;
        while (true)
        {
            IReadOnlyList<int> children = phylogeny.ChildrenOf(current);
            if (children.Count != 1)
            {
                break;
            }

            int only = children[0];
            trunk.Add(only);
            current = only;
        }

        return trunk;
    }

    /// <summary>
    /// Returns the ids of <b>branch</b> (subclonal) clusters — every input cluster that is not on the trunk.
    /// Source: Popic et al. (2015) — mutations off the common-predecessor trunk are subclonal lineage branches.
    /// </summary>
    /// <param name="phylogeny">A phylogeny produced by <see cref="ReconstructPhylogeny"/>.</param>
    /// <returns>Branch cluster ids in input order.</returns>
    public static IReadOnlyList<int> IdentifyBranchMutations(ClonalPhylogeny phylogeny)
    {
        var branches = new List<int>();
        if (phylogeny.Clusters is null || phylogeny.Clusters.Count == 0)
        {
            return branches;
        }

        var trunk = new HashSet<int>(IdentifyTrunkMutations(phylogeny));
        foreach (CcfCluster c in phylogeny.Clusters)
        {
            if (!trunk.Contains(c.Id))
            {
                branches.Add(c.Id);
            }
        }

        return branches;
    }

    /// <summary>
    /// Lineage precedence (Popic 2015 Eq. 2): for every sample i, <c>parent ≥ child − ε</c> and a parent absent in
    /// a sample (CCF = 0) cannot have a child present there (presence pattern, constraint (1)).
    /// </summary>
    private static bool SatisfiesLineagePrecedence(double[] parentCcf, double[] childCcf, double tolerance)
    {
        for (int i = 0; i < parentCcf.Length; i++)
        {
            if (parentCcf[i] < childCcf[i] - tolerance)
            {
                return false;
            }

            // Presence: if the parent is absent in sample i, the child must also be absent there.
            if (parentCcf[i] <= 0.0 && childCcf[i] > 0.0)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>
    /// Sum rule (Popic 2015 Eq. 5): the child's CCF fits a parent only if, in every sample, the parent's remaining
    /// budget (parent CCF minus already-attached children) is ≥ child CCF − ε.
    /// </summary>
    private static bool FitsSumRule(double[] remainingParentBudget, double[] childCcf, double tolerance)
    {
        for (int i = 0; i < remainingParentBudget.Length; i++)
        {
            if (remainingParentBudget[i] < childCcf[i] - tolerance)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>Sum of a cluster's CCF over all samples (proxy for clonal dominance / processing order).</summary>
    private static double TotalCcf(double[] ccf)
    {
        double total = 0.0;
        foreach (double v in ccf)
        {
            total += v;
        }

        return total;
    }

    /// <summary>Candidate parents = the root plus every cluster already attached to the tree (in attachment order).</summary>
    private static IEnumerable<int> EnumerateCandidates(int rootId, List<ClonalEdge> edges)
    {
        yield return rootId;
        foreach (ClonalEdge e in edges)
        {
            yield return e.ChildId;
        }
    }

    /// <summary>Chooses a synthetic root id distinct from every cluster id (one less than the minimum, or -1).</summary>
    private static int RootIdFor(IReadOnlyList<CcfCluster> clusters)
    {
        int minId = int.MaxValue;
        foreach (CcfCluster c in clusters)
        {
            if (c.Id < minId)
            {
                minId = c.Id;
            }
        }

        // -1 is conventional for "no real cluster"; if a caller used -1, step below the minimum to stay unique.
        return clusters.Count == 0 ? -1 : Math.Min(-1, minId - 1);
    }

    /// <summary>Validates every cluster's CCF list and returns the common sample count.</summary>
    private static int ValidateAndGetSampleCount(IReadOnlyList<CcfCluster> clusters)
    {
        int sampleCount = -1;
        var seenIds = new HashSet<int>(clusters.Count);
        foreach (CcfCluster c in clusters)
        {
            if (c.CcfPerSample is null)
            {
                throw new ArgumentNullException(nameof(clusters), $"Cluster {c.Id} has a null CCF list.");
            }

            if (c.CcfPerSample.Count == 0)
            {
                throw new ArgumentException($"Cluster {c.Id} has an empty CCF list; at least one sample is required.", nameof(clusters));
            }

            if (sampleCount < 0)
            {
                sampleCount = c.CcfPerSample.Count;
            }
            else if (c.CcfPerSample.Count != sampleCount)
            {
                throw new ArgumentException(
                    $"All clusters must have the same number of samples; cluster {c.Id} has {c.CcfPerSample.Count}, expected {sampleCount}.",
                    nameof(clusters));
            }

            if (!seenIds.Add(c.Id))
            {
                throw new ArgumentException($"Duplicate cluster id {c.Id}; cluster ids must be unique.", nameof(clusters));
            }

            foreach (double v in c.CcfPerSample)
            {
                if (double.IsNaN(v) || v < 0.0 || v > 1.0)
                {
                    throw new ArgumentException(
                        $"Cancer cell fraction must be in [0, 1]; cluster {c.Id} has {v}.", nameof(clusters));
                }
            }
        }

        return sampleCount;
    }

    #endregion

    #region Tumor Heterogeneity Analysis (ONCO-HETERO-001)

    /// <summary>
    /// MAD consistency-scaling constant 1.4826 = 1/Φ⁻¹(3/4): scales the raw median absolute deviation so that for
    /// a normally distributed variable the expected MAD equals its standard deviation. Source: Mroz &amp; Rocco
    /// (2013), <i>Oral Oncology</i> 49(3):211–215 / Mroz et al. (2015), <i>PLOS Medicine</i> 12(2):e1001786 —
    /// "The median [absolute deviation] is then multiplied by a factor of 1.4826, so that the expected MAD of a
    /// normally distributed variable is equal to its SD"; maftools <c>mathScore.R</c> uses the same 1.4826.
    /// </summary>
    private const double MadConsistencyConstant = 1.4826;

    /// <summary>
    /// Percentage-scaling factor in the MATH score. Source: Mroz &amp; Rocco (2013) / Mroz et al. (2015):
    /// "MATH = 100 × MAD/median"; maftools <c>mathScore.R</c>: <c>pat.math = pat.mad * 1.4826 / median(vaf)</c>
    /// with <c>pat.mad = median(abs.med.dev) * 100</c>.
    /// </summary>
    private const double MathPercentScale = 100.0;

    /// <summary>
    /// Result of a tumour intratumour-heterogeneity (ITH) analysis over a set of somatic mutations.
    /// </summary>
    /// <param name="MathScore">Mutant-Allele Tumour Heterogeneity (MATH) score = 100·1.4826·MAD(VAF)/median(VAF),
    /// computed over the mutant-allele (variant) fractions (Mroz &amp; Rocco 2013).</param>
    /// <param name="ShannonDiversity">Shannon diversity index H = −Σ pᵢ·ln(pᵢ) over the clone fractions pᵢ
    /// (fraction of mutations assigned to each CCF cluster), using the natural logarithm (Shannon 1948).</param>
    /// <param name="SubcloneCount">Number of distinct clones/subclones = number of non-empty CCF clusters.</param>
    /// <param name="SubclonalFraction">Fraction of mutations whose CCF is below the clonal threshold (CCF &lt; 0.95,
    /// Landau et al. 2013) and are therefore subclonal.</param>
    public readonly record struct HeterogeneityResult(
        double MathScore,
        double ShannonDiversity,
        int SubcloneCount,
        double SubclonalFraction);

    /// <summary>
    /// Computes the intratumour-heterogeneity (ITH) score as the Mutant-Allele Tumour Heterogeneity (MATH) score
    /// over a distribution of mutant-allele (variant) fractions:
    /// <c>MATH = 100 · 1.4826 · median(|fᵢ − median(f)|) / median(f)</c>, the ratio of the width (scaled median
    /// absolute deviation) to the centre (median) of the VAF distribution, expressed as a percentage. A wider,
    /// more dispersed VAF distribution (more genetic heterogeneity) gives a higher MATH. Source: Mroz &amp; Rocco
    /// (2013), <i>Oral Oncology</i> 49(3):211–215; Mroz et al. (2015), <i>PLOS Medicine</i> 12(2):e1001786
    /// ("MATH = 100 × MAD/median"); maftools <c>mathScore.R</c>.
    /// </summary>
    /// <param name="ccfDistribution">Mutant-allele fractions (VAFs) of the tumour's somatic mutations, each in
    /// [0, 1]; the median must be strictly positive (MATH divides by it).</param>
    /// <returns>The MATH score (≥ 0; 0 when every value equals the median).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfDistribution"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, a value is non-finite or outside [0, 1], or the
    /// median is 0 (the MATH ratio is undefined).</exception>
    public static double CalculateITH(IReadOnlyList<double> ccfDistribution)
    {
        ArgumentNullException.ThrowIfNull(ccfDistribution);

        int n = ccfDistribution.Count;
        if (n == 0)
        {
            throw new ArgumentException("At least one allele fraction is required.", nameof(ccfDistribution));
        }

        double[] values = new double[n];
        for (int i = 0; i < n; i++)
        {
            double v = ccfDistribution[i];
            if (double.IsNaN(v) || double.IsInfinity(v) || v < 0.0 || v > 1.0)
            {
                throw new ArgumentException(
                    $"Allele fraction must be a finite value in [0, 1]; got {v} at index {i}.", nameof(ccfDistribution));
            }

            values[i] = v;
        }

        double median = Median(values);
        if (median == 0.0)
        {
            throw new ArgumentException(
                "The median allele fraction is 0; the MATH score divides by the median and is undefined.",
                nameof(ccfDistribution));
        }

        // Raw MAD = median of absolute deviations from the median; scale by 1.4826 for normal consistency.
        double[] absDeviations = new double[n];
        for (int i = 0; i < n; i++)
        {
            absDeviations[i] = Math.Abs(values[i] - median);
        }

        double rawMad = Median(absDeviations);
        double scaledMad = MadConsistencyConstant * rawMad;
        return MathPercentScale * scaledMad / median;
    }

    /// <summary>
    /// Counts the number of distinct clones/subclones in a tumour as the number of non-empty CCF clusters produced
    /// by <see cref="ClusterCcfValues"/> (ONCO-CCF-001). Each cluster centroid is a clonal population; the count is
    /// the tumour's clonal richness. Source: Mroz et al. (2015) treat genetic heterogeneity as the number/spread of
    /// subpopulations; Liu et al. (2017), <i>BMC Genomics</i> 18:457 (PMC5468233) define richness as "the number of
    /// clones present" when computing Shannon-based ITH scores.
    /// </summary>
    /// <param name="ccfClusters">A CCF clustering (its <see cref="CcfClustering.Assignments"/> determine which
    /// clusters actually contain at least one mutation).</param>
    /// <returns>The number of clusters that contain at least one assigned mutation (≥ 1).</returns>
    /// <exception cref="ArgumentException"><paramref name="ccfClusters"/> has no centroids or no assignments.</exception>
    public static int InferSubclones(CcfClustering ccfClusters)
    {
        IReadOnlyList<int> assignments = ccfClusters.Assignments;
        if (assignments is null || assignments.Count == 0 || ccfClusters.Centroids is null || ccfClusters.Centroids.Count == 0)
        {
            throw new ArgumentException("The CCF clustering must contain at least one cluster and one assignment.", nameof(ccfClusters));
        }

        var occupied = new HashSet<int>();
        foreach (int label in assignments)
        {
            occupied.Add(label);
        }

        return occupied.Count;
    }

    /// <summary>
    /// Performs a tumour intratumour-heterogeneity (ITH) analysis from per-mutation variant allele fractions and
    /// cancer cell fractions, returning four standard ITH metrics: the MATH score over the VAFs (Mroz &amp; Rocco
    /// 2013), the Shannon diversity index H = −Σ pᵢ·ln(pᵢ) over the clone fractions (Shannon 1948), the number of
    /// subclones (CCF clusters), and the fraction of subclonal mutations (CCF &lt; 0.95, Landau et al. 2013). CCF
    /// values are clustered with <see cref="ClusterCcfValues"/> (ONCO-CCF-001) into <paramref name="clusterCount"/>
    /// clones; the clone fractions pᵢ are the proportions of mutations assigned to each cluster.
    /// </summary>
    /// <param name="variantAlleleFractions">Per-mutation VAFs in [0, 1] (used for the MATH score).</param>
    /// <param name="ccfValues">Per-mutation cancer cell fractions in [0, 1], aligned with
    /// <paramref name="variantAlleleFractions"/> (same count and order).</param>
    /// <param name="clusterCount">Number of clones k to cluster the CCF values into, in [1, count].</param>
    /// <returns>The MATH score, Shannon diversity, subclone count, and subclonal fraction.</returns>
    /// <exception cref="ArgumentNullException">either list is null.</exception>
    /// <exception cref="ArgumentException">the lists are empty or of different length, or values are out of range.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count].</exception>
    public static HeterogeneityResult AnalyzeHeterogeneity(
        IReadOnlyList<double> variantAlleleFractions,
        IReadOnlyList<double> ccfValues,
        int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(variantAlleleFractions);
        ArgumentNullException.ThrowIfNull(ccfValues);

        if (variantAlleleFractions.Count == 0)
        {
            throw new ArgumentException("At least one mutation is required.", nameof(variantAlleleFractions));
        }

        if (variantAlleleFractions.Count != ccfValues.Count)
        {
            throw new ArgumentException(
                $"VAF and CCF lists must have the same length; got {variantAlleleFractions.Count} and {ccfValues.Count}.",
                nameof(ccfValues));
        }

        double math = CalculateITH(variantAlleleFractions);

        // Cluster CCF values into clones (ONCO-CCF-001); validates ccfValues and clusterCount.
        CcfClustering clustering = ClusterCcfValues(ccfValues, clusterCount);
        int subcloneCount = InferSubclones(clustering);

        // Clone fractions pᵢ = proportion of mutations in each occupied cluster; Shannon H = −Σ pᵢ·ln pᵢ.
        int n = ccfValues.Count;
        var clusterSizes = new Dictionary<int, int>();
        foreach (int label in clustering.Assignments)
        {
            clusterSizes[label] = clusterSizes.TryGetValue(label, out int existing) ? existing + 1 : 1;
        }

        double shannon = 0.0;
        foreach (int size in clusterSizes.Values)
        {
            double p = (double)size / n;
            shannon -= p * Math.Log(p);
        }

        // Subclonal mutations: CCF strictly below the clonal threshold (Landau et al. 2013, reused from ONCO-CLONAL-001).
        int subclonal = 0;
        for (int i = 0; i < n; i++)
        {
            if (ccfValues[i] < ClonalCcfThreshold)
            {
                subclonal++;
            }
        }

        double subclonalFraction = (double)subclonal / n;
        return new HeterogeneityResult(math, shannon, subcloneCount, subclonalFraction);
    }

    /// <summary>
    /// Returns the median of the supplied values. For an even count the median is the arithmetic mean of the two
    /// central order statistics; for an odd count it is the central value (standard definition; matches R's
    /// <c>median</c> used by maftools <c>mathScore.R</c>). Does not mutate the input.
    /// </summary>
    private static double Median(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        int mid = n / 2;
        return (n % 2 == 1) ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2.0;
    }

    #endregion

    #region HLA nomenclature parsing and allele-specific HLA LOH (ONCO-HLA-001)

    /// <summary>
    /// Minimum number of colon-separated numeric fields in a valid HLA allele name. Source: WHO HLA
    /// Nomenclature (hla.alleles.org "Naming Alleles"); Marsh et al. (2010), Tissue Antigens 75(4):291–455 —
    /// "All alleles receive at least a four digit name, which corresponds to the first two sets of digits".
    /// </summary>
    private const int HlaMinFieldCount = 2;

    /// <summary>
    /// Maximum number of colon-separated numeric fields in a valid HLA allele name (type : protein :
    /// synonymous coding : non-coding). Source: WHO HLA Nomenclature — "up to four sets of digits separated
    /// by colons" (hla.alleles.org "Naming Alleles"; Marsh et al. 2010).
    /// </summary>
    private const int HlaMaxFieldCount = 4;

    /// <summary>
    /// Allele copy-number threshold below which an HLA allele is classified as lost. Source (verbatim):
    /// McGranahan et al. (2017), Cell 171(6):1259–1271 (PMC5720478) — "A copy number &lt; 0.5, is classified
    /// as subject to loss, and thereby indicative of LOH." The comparison is strict (CN must be &lt; 0.5).
    /// </summary>
    public const double HlaLohCopyNumberThreshold = 0.5;

    /// <summary>
    /// Allelic-imbalance significance threshold required before HLA LOH may be called. Source (verbatim):
    /// McGranahan et al. (2017), Cell 171(6):1259–1271 — "To avoid over-calling LOH, we also calculate a p
    /// value relating to allelic imbalance for each HLA gene. Allelic imbalance is determined if p &lt; 0.01
    /// using the paired Student's t-Test." The comparison is strict (p must be &lt; 0.01). Corroborated by the
    /// paired <c>t.test(..., paired=TRUE)</c> in mskcc/lohhla <c>LOHHLAscript.R</c>.
    /// </summary>
    public const double HlaLohAllelicImbalancePValueThreshold = 0.01;

    /// <summary>HLA expression-status suffix per WHO HLA Nomenclature ("Naming Alleles", hla.alleles.org).</summary>
    public enum HlaExpressionSuffix
    {
        /// <summary>No suffix — normally expressed allele.</summary>
        None,

        /// <summary><c>N</c> — Null allele (not expressed).</summary>
        Null,

        /// <summary><c>L</c> — Low cell-surface expression.</summary>
        Low,

        /// <summary><c>S</c> — Secreted molecule only (not on the cell surface).</summary>
        Secreted,

        /// <summary><c>C</c> — present in the Cytoplasm but not on the cell surface.</summary>
        Cytoplasm,

        /// <summary><c>A</c> — Aberrant expression (uncertain whether expressed).</summary>
        Aberrant,

        /// <summary><c>Q</c> — Questionable expression.</summary>
        Questionable,
    }

    /// <summary>
    /// Which of the two homologous HLA alleles at a locus was lost (had copy number below the loss threshold).
    /// </summary>
    public enum HlaLostAllele
    {
        /// <summary>No allele was lost (locus retained, no allele-specific LOH).</summary>
        None,

        /// <summary>The first allele (allele 1) is the lost allele.</summary>
        Allele1,

        /// <summary>The second allele (allele 2) is the lost allele.</summary>
        Allele2,

        /// <summary>Both alleles fell below the loss threshold (homozygous loss; not allele-specific LOH).</summary>
        Both,
    }

    /// <summary>
    /// A parsed HLA allele name per the WHO HLA Nomenclature. Format: <c>HLA-&lt;Gene&gt;*F1:F2[:F3[:F4]][suffix]</c>.
    /// Source: hla.alleles.org "Naming Alleles"; Marsh et al. (2010), Tissue Antigens 75(4):291–455.
    /// </summary>
    /// <param name="Gene">Gene name (e.g. "A", "B", "C", "DRB1"), upper-cased.</param>
    /// <param name="Fields">The 2–4 numeric fields, each as its original digit string (e.g. ["02","01"]).
    /// Field 1 = type/allele group; Field 2 = specific HLA protein; Field 3 = synonymous coding substitutions;
    /// Field 4 = non-coding differences.</param>
    /// <param name="Suffix">Optional expression-status suffix (None when absent).</param>
    public readonly record struct HlaAllele(string Gene, IReadOnlyList<string> Fields, HlaExpressionSuffix Suffix)
    {
        /// <summary>The allele-group (first) field — the serological "type" digits.</summary>
        public string AlleleGroup => Fields[0];

        /// <summary>The specific-HLA-protein (second) field.</summary>
        public string Protein => Fields[1];

        /// <summary>The canonical normalized allele name, e.g. <c>HLA-A*02:01:01:02L</c>.</summary>
        public string Name
        {
            get
            {
                string body = "HLA-" + Gene + "*" + string.Join(":", Fields);
                return Suffix == HlaExpressionSuffix.None ? body : body + SuffixLetter(Suffix);
            }
        }
    }

    /// <summary>
    /// Caller-supplied allele-specific copy-number evidence at one HLA gene, as produced by an HLA copy-number
    /// caller (e.g. LOHHLA): the estimated copy number of each of the two homologous alleles plus the
    /// allelic-imbalance p value. Source: McGranahan et al. (2017) — LOHHLA reports per-allele copy number
    /// (<c>HLA_type1copyNum</c>, <c>HLA_type2copyNum</c>) and a paired-t-test allelic-imbalance p value.
    /// </summary>
    /// <param name="Allele1">Allele 1 name (informational).</param>
    /// <param name="Allele1CopyNumber">Estimated copy number of allele 1 (≥ 0).</param>
    /// <param name="Allele2">Allele 2 name (informational).</param>
    /// <param name="Allele2CopyNumber">Estimated copy number of allele 2 (≥ 0).</param>
    /// <param name="AllelicImbalancePValue">Paired Student's t-test p value for allelic imbalance, in [0, 1].</param>
    public readonly record struct HlaAlleleCopyNumber(
        string Allele1,
        double Allele1CopyNumber,
        string Allele2,
        double Allele2CopyNumber,
        double AllelicImbalancePValue);

    /// <summary>Result of an allele-specific HLA LOH determination (LOHHLA classification).</summary>
    /// <param name="IsLoh">True iff allele-specific HLA LOH was called (one allele lost with significant imbalance).</param>
    /// <param name="LostAllele">Which allele was lost (<see cref="HlaLostAllele.None"/> when no LOH).</param>
    /// <param name="AllelicImbalanceSignificant">True iff the allelic-imbalance p value is below the threshold.</param>
    public readonly record struct HlaLohResult(bool IsLoh, HlaLostAllele LostAllele, bool AllelicImbalanceSignificant);

    private static char SuffixLetter(HlaExpressionSuffix suffix) => suffix switch
    {
        HlaExpressionSuffix.Null => 'N',
        HlaExpressionSuffix.Low => 'L',
        HlaExpressionSuffix.Secreted => 'S',
        HlaExpressionSuffix.Cytoplasm => 'C',
        HlaExpressionSuffix.Aberrant => 'A',
        HlaExpressionSuffix.Questionable => 'Q',
        _ => '\0',
    };

    private static bool TryMapSuffix(char letter, out HlaExpressionSuffix suffix)
    {
        // WHO HLA Nomenclature expression-status suffixes (hla.alleles.org "Naming Alleles").
        switch (letter)
        {
            case 'N': suffix = HlaExpressionSuffix.Null; return true;
            case 'L': suffix = HlaExpressionSuffix.Low; return true;
            case 'S': suffix = HlaExpressionSuffix.Secreted; return true;
            case 'C': suffix = HlaExpressionSuffix.Cytoplasm; return true;
            case 'A': suffix = HlaExpressionSuffix.Aberrant; return true;
            case 'Q': suffix = HlaExpressionSuffix.Questionable; return true;
            default: suffix = HlaExpressionSuffix.None; return false;
        }
    }

    /// <summary>
    /// Parses and validates an HLA allele name per the WHO HLA Nomenclature
    /// (<c>HLA-&lt;Gene&gt;*F1:F2[:F3[:F4]][suffix]</c>). The input is trimmed and the gene name upper-cased;
    /// the <c>HLA-</c> prefix is mandatory, the gene is separated from the fields by <c>*</c>, fields are
    /// colon-separated digit groups, and an optional trailing expression-status letter (N/L/S/C/A/Q) may
    /// follow the last field. Source: hla.alleles.org "Naming Alleles"; Marsh et al. (2010), Tissue Antigens
    /// 75(4):291–455.
    /// </summary>
    /// <param name="alleleName">The HLA allele name to parse.</param>
    /// <returns>The parsed <see cref="HlaAllele"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alleleName"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="alleleName"/> is empty or whitespace.</exception>
    /// <exception cref="FormatException">The name does not conform to WHO HLA nomenclature (missing
    /// <c>HLA-</c> prefix or <c>*</c>; fewer than 2 or more than 4 fields; a non-numeric field; an invalid
    /// trailing suffix).</exception>
    public static HlaAllele ParseHlaAllele(string alleleName)
    {
        ArgumentNullException.ThrowIfNull(alleleName);
        if (string.IsNullOrWhiteSpace(alleleName))
        {
            throw new ArgumentException("HLA allele name must not be empty or whitespace.", nameof(alleleName));
        }

        string text = alleleName.Trim();

        // The "HLA-" prefix is mandatory (case-insensitive); the gene is separated from the fields by '*'.
        const string Prefix = "HLA-";
        if (!text.StartsWith(Prefix, StringComparison.OrdinalIgnoreCase))
        {
            throw new FormatException($"HLA allele name must start with '{Prefix}': '{alleleName}'.");
        }

        string afterPrefix = text.Substring(Prefix.Length);
        int star = afterPrefix.IndexOf('*');
        if (star <= 0 || star == afterPrefix.Length - 1)
        {
            throw new FormatException($"HLA allele name must have a gene and field block separated by '*': '{alleleName}'.");
        }

        string gene = afterPrefix.Substring(0, star).ToUpperInvariant();
        string fieldBlock = afterPrefix.Substring(star + 1);

        // Strip an optional single trailing expression-status suffix (only when the last char is a letter).
        var suffix = HlaExpressionSuffix.None;
        char last = fieldBlock[fieldBlock.Length - 1];
        if (char.IsLetter(last))
        {
            if (!TryMapSuffix(char.ToUpperInvariant(last), out suffix))
            {
                throw new FormatException($"Invalid HLA expression-status suffix '{last}' (allowed: N, L, S, C, A, Q): '{alleleName}'.");
            }

            fieldBlock = fieldBlock.Substring(0, fieldBlock.Length - 1);
        }

        string[] fields = fieldBlock.Split(':');
        if (fields.Length < HlaMinFieldCount || fields.Length > HlaMaxFieldCount)
        {
            throw new FormatException(
                $"HLA allele name must have between {HlaMinFieldCount} and {HlaMaxFieldCount} colon-separated fields, found {fields.Length}: '{alleleName}'.");
        }

        foreach (string field in fields)
        {
            // WHO HLA Nomenclature fields are ASCII decimal-digit groups ('0'–'9'); char.IsDigit also accepts
            // non-ASCII Unicode decimal digits (e.g. fullwidth '０', Arabic-Indic '٢'), which are not valid
            // nomenclature and must be rejected as malformed.
            if (field.Length == 0 || !field.All(static c => c is >= '0' and <= '9'))
            {
                throw new FormatException($"HLA allele field '{field}' must be a non-empty digit group: '{alleleName}'.");
            }
        }

        return new HlaAllele(gene, fields, suffix);
    }

    /// <summary>
    /// Non-throwing variant of <see cref="ParseHlaAllele(string)"/>. Returns true and the parsed allele on
    /// success; false and the default value when the name is null/empty or does not conform to WHO HLA
    /// nomenclature.
    /// </summary>
    /// <param name="alleleName">The HLA allele name to parse.</param>
    /// <param name="allele">The parsed allele on success; default otherwise.</param>
    /// <returns>True if parsing succeeded; otherwise false.</returns>
    public static bool TryParseHlaAllele(string? alleleName, out HlaAllele allele)
    {
        if (alleleName is null)
        {
            allele = default;
            return false;
        }

        try
        {
            allele = ParseHlaAllele(alleleName);
            return true;
        }
        catch (Exception ex) when (ex is FormatException or ArgumentException)
        {
            allele = default;
            return false;
        }
    }

    /// <summary>
    /// Classifies allele-specific HLA loss of heterozygosity from caller-supplied per-allele copy number and
    /// the allelic-imbalance p value, per the LOHHLA decision rule. HLA LOH is called iff <b>exactly one</b>
    /// of the two homologous alleles has copy number strictly below
    /// <see cref="HlaLohCopyNumberThreshold"/> (0.5) <b>and</b> the allelic-imbalance p value is strictly
    /// below <see cref="HlaLohAllelicImbalancePValueThreshold"/> (0.01). If both alleles fall below the loss
    /// threshold the locus is reported as homozygous loss (<see cref="HlaLostAllele.Both"/>), which is not
    /// allele-specific LOH. Source: McGranahan et al. (2017), Cell 171(6):1259–1271 (PMC5720478) — "A copy
    /// number &lt; 0.5, is classified as subject to loss, and thereby indicative of LOH" and "Allelic
    /// imbalance is determined if p &lt; 0.01 using the paired Student's t-Test".
    /// </summary>
    /// <param name="alleleCopyNumber">Caller-supplied per-allele copy number and allelic-imbalance p value.</param>
    /// <returns>The HLA LOH classification.</returns>
    /// <exception cref="ArgumentException">A copy number is negative, or the p value is outside [0, 1].</exception>
    public static HlaLohResult DetectHlaLoh(HlaAlleleCopyNumber alleleCopyNumber)
    {
        // Copy numbers must be finite and non-negative. NaN must be rejected explicitly: `NaN < 0` is false, so
        // a NaN copy number would otherwise slip past the bound and leak into the loss-threshold comparison.
        if (!double.IsFinite(alleleCopyNumber.Allele1CopyNumber) || !double.IsFinite(alleleCopyNumber.Allele2CopyNumber)
            || alleleCopyNumber.Allele1CopyNumber < 0 || alleleCopyNumber.Allele2CopyNumber < 0)
        {
            throw new ArgumentException("HLA allele copy numbers must be finite and non-negative.", nameof(alleleCopyNumber));
        }

        // The p value must be finite and in [0, 1]. NaN must be rejected explicitly: both `NaN < 0.0` and
        // `NaN > 1.0` are false, so a NaN p value would otherwise pass this check and silently be treated as
        // not-significant in the `< 0.01` test rather than reported as malformed input.
        if (!double.IsFinite(alleleCopyNumber.AllelicImbalancePValue)
            || alleleCopyNumber.AllelicImbalancePValue is < 0.0 or > 1.0)
        {
            throw new ArgumentException("Allelic-imbalance p value must be in [0, 1].", nameof(alleleCopyNumber));
        }

        // LOHHLA over-calling guard: significant allelic imbalance (paired t-test p < 0.01) is required.
        bool imbalanceSignificant = alleleCopyNumber.AllelicImbalancePValue < HlaLohAllelicImbalancePValueThreshold;

        bool allele1Lost = alleleCopyNumber.Allele1CopyNumber < HlaLohCopyNumberThreshold;
        bool allele2Lost = alleleCopyNumber.Allele2CopyNumber < HlaLohCopyNumberThreshold;

        if (allele1Lost && allele2Lost)
        {
            // Both homologs below 0.5 → homozygous loss, not allele-specific LOH (Evidence assumption).
            return new HlaLohResult(false, HlaLostAllele.Both, imbalanceSignificant);
        }

        if (imbalanceSignificant && allele1Lost)
        {
            return new HlaLohResult(true, HlaLostAllele.Allele1, true);
        }

        if (imbalanceSignificant && allele2Lost)
        {
            return new HlaLohResult(true, HlaLostAllele.Allele2, true);
        }

        return new HlaLohResult(false, HlaLostAllele.None, imbalanceSignificant);
    }

    #endregion

    #region Clinical Actionability (OncoKB Therapeutic Levels of Evidence)

    /// <summary>
    /// OncoKB therapeutic level of evidence assigned to a biomarker–drug association, indicating how
    /// strongly an alteration is predictive of sensitivity (Levels 1, 2, 3A, 3B, 4) or resistance
    /// (Levels R1, R2) to a therapy. Source: Chakravarty D et al. (2017), JCO Precis Oncol 2017:1–16;
    /// definitions verbatim from the OncoKB Therapeutic Levels of Evidence (V2) document.
    /// </summary>
    public enum OncoKbLevel
    {
        /// <summary>No leveled therapeutic association (variant is not actionable on this axis).</summary>
        None,

        /// <summary>
        /// Level R2 — Investigational Resistance: "Compelling clinical evidence supports the biomarker as
        /// being predictive of resistance to a drug." Lowest in the combined order. Source: OncoKB Levels V2.
        /// </summary>
        R2,

        /// <summary>
        /// Level 4 — Hypothetical: "Compelling biological evidence supports the biomarker as being
        /// predictive of response to a drug." Source: OncoKB Levels V2.
        /// </summary>
        Level4,

        /// <summary>
        /// Level 3B — Investigational: "Standard care or investigational biomarker predictive of response
        /// to an FDA-approved or investigational drug in another indication." Source: OncoKB Levels V2.
        /// </summary>
        Level3B,

        /// <summary>
        /// Level 3A — Investigational: "Compelling clinical evidence supports the biomarker as being
        /// predictive of response to a drug in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level3A,

        /// <summary>
        /// Level 2 — Standard Care: "Standard care biomarker recommended by the NCCN or other professional
        /// guidelines predictive of response to an FDA-approved drug in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level2,

        /// <summary>
        /// Level 1 — Standard Care: "FDA-recognized biomarker predictive of response to an FDA-approved drug
        /// in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level1,

        /// <summary>
        /// Level R1 — Standard Care Resistance: "Standard care biomarker predictive of resistance to an
        /// FDA-approved drug in this indication." Highest in the combined order. Source: OncoKB Levels V2.
        /// </summary>
        R1
    }

    /// <summary>
    /// Combined actionability ranking of the OncoKB levels, highest first. The integer order of this array
    /// encodes the OncoKB HIGHEST_LEVEL precedence: R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2. Source:
    /// oncokb-annotator README, column HIGHEST_LEVEL — "Order: LEVEL_R1 &gt; LEVEL_1 &gt; LEVEL_2 &gt;
    /// LEVEL_3A &gt; LEVEL_3B &gt; LEVEL_4 &gt; LEVEL_R2".
    /// </summary>
    private static readonly OncoKbLevel[] CombinedRankingHighestFirst =
    {
        OncoKbLevel.R1,
        OncoKbLevel.Level1,
        OncoKbLevel.Level2,
        OncoKbLevel.Level3A,
        OncoKbLevel.Level3B,
        OncoKbLevel.Level4,
        OncoKbLevel.R2
    };

    /// <summary>Levels that denote sensitivity (response) to a therapy. Source: OncoKB Levels V2 (1/2/3A/3B/4).</summary>
    private static readonly HashSet<OncoKbLevel> SensitivityLevels = new()
    {
        OncoKbLevel.Level1, OncoKbLevel.Level2, OncoKbLevel.Level3A, OncoKbLevel.Level3B, OncoKbLevel.Level4
    };

    /// <summary>Levels that denote resistance to a therapy. Source: OncoKB Levels V2 (R1/R2).</summary>
    private static readonly HashSet<OncoKbLevel> ResistanceLevels = new()
    {
        OncoKbLevel.R1, OncoKbLevel.R2
    };

    /// <summary>
    /// Levels categorized as "standard care" by OncoKB (as opposed to investigational/hypothetical):
    /// Levels 1, 2 (sensitivity) and R1 (resistance). Source: OncoKB Curation SOP v3 — "The highest levels
    /// of evidence, Levels 1 and 2, refer to the standard implications... Level R1 refers to the standard
    /// implications for resistance"; Levels 3A/3B/4/R2 are investigational/hypothetical.
    /// </summary>
    private static readonly HashSet<OncoKbLevel> StandardCareLevels = new()
    {
        OncoKbLevel.Level1, OncoKbLevel.Level2, OncoKbLevel.R1
    };

    /// <summary>
    /// One caller-supplied biomarker–drug therapeutic association from a precision-oncology knowledgebase
    /// (e.g. an OncoKB export). The library does not embed the OncoKB curated content (3,000+ alterations
    /// across 418 genes, Chakravarty 2017); the caller performs the lookup and supplies the relevant rows.
    /// </summary>
    /// <param name="Drug">Therapy name the association refers to.</param>
    /// <param name="Level">OncoKB therapeutic level of evidence for this drug–variant association.</param>
    public readonly record struct TherapyAssociation(string Drug, OncoKbLevel Level);

    /// <summary>
    /// Caller-supplied evidence for one variant's clinical actionability: the gene/protein change (for
    /// reporting) and the set of leveled drug associations curated for it. Mirrors the framework boundary of
    /// <see cref="CancerVariantAnnotationInput"/> — actionability comes from a caller-supplied knowledgebase.
    /// </summary>
    /// <param name="Gene">Gene symbol the variant falls in.</param>
    /// <param name="ProteinChange">HGVS protein change (e.g. p.V600E); informational.</param>
    /// <param name="Associations">Leveled drug associations from the knowledgebase (may be empty, never null).</param>
    public readonly record struct VariantActionabilityInput
    {
        /// <summary>Gene symbol the variant falls in.</summary>
        public string Gene { get; }

        /// <summary>HGVS protein change (e.g. p.V600E); informational.</summary>
        public string ProteinChange { get; }

        /// <summary>Leveled drug associations from the caller-supplied knowledgebase.</summary>
        public IReadOnlyList<TherapyAssociation> Associations { get; }

        /// <summary>Creates an actionability input. <paramref name="associations"/> must not be null.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="associations"/> is null.</exception>
        public VariantActionabilityInput(
            string gene, string proteinChange, IReadOnlyList<TherapyAssociation> associations)
        {
            ArgumentNullException.ThrowIfNull(associations);
            Gene = gene;
            ProteinChange = proteinChange;
            Associations = associations;
        }
    }

    /// <summary>
    /// Per-variant clinical actionability assessment under the OncoKB therapeutic levels system.
    /// </summary>
    /// <param name="Variant">The variant evidence that was assessed.</param>
    /// <param name="HighestSensitiveLevel">Highest sensitivity level (1 &gt; 2 &gt; 3A &gt; 3B &gt; 4), or None.</param>
    /// <param name="HighestResistanceLevel">Highest resistance level (R1 &gt; R2), or None.</param>
    /// <param name="HighestCombinedLevel">Highest level over both axes (R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2), or None.</param>
    public readonly record struct ActionabilityAssessment(
        VariantActionabilityInput Variant,
        OncoKbLevel HighestSensitiveLevel,
        OncoKbLevel HighestResistanceLevel,
        OncoKbLevel HighestCombinedLevel)
    {
        /// <summary>True when the variant has at least one leveled therapeutic association.</summary>
        public bool IsActionable => HighestCombinedLevel != OncoKbLevel.None;
    }

    /// <summary>
    /// Compares two OncoKB levels by the combined actionability order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2.
    /// Returns a positive number when <paramref name="a"/> is more actionable than <paramref name="b"/>,
    /// negative when less, zero when equal. <see cref="OncoKbLevel.None"/> ranks below every leveled value.
    /// Source: oncokb-annotator README HIGHEST_LEVEL order.
    /// </summary>
    /// <param name="a">First level.</param>
    /// <param name="b">Second level.</param>
    /// <returns>Sign indicates which level is higher in the combined order.</returns>
    public static int CompareLevels(OncoKbLevel a, OncoKbLevel b)
        // The enum is declared in ascending actionability (None lowest, R1 highest), so the underlying
        // integer value already encodes the combined order; comparing values is the documented precedence.
        => ((int)a).CompareTo((int)b);

    /// <summary>
    /// Returns the highest OncoKB level over a set of levels using the combined order R1 &gt; 1 &gt; 2 &gt;
    /// 3A &gt; 3B &gt; 4 &gt; R2, restricted to the levels in <paramref name="allowed"/> (used to compute the
    /// sensitivity-only and resistance-only maxima). Returns <see cref="OncoKbLevel.None"/> when no allowed
    /// level is present. Source: oncokb-annotator README HIGHEST_*_LEVEL orders.
    /// </summary>
    private static OncoKbLevel HighestLevel(
        IReadOnlyList<TherapyAssociation> associations, HashSet<OncoKbLevel>? allowed)
    {
        OncoKbLevel best = OncoKbLevel.None;
        foreach (var association in associations)
        {
            if (allowed is not null && !allowed.Contains(association.Level))
            {
                continue;
            }

            if (CompareLevels(association.Level, best) > 0)
            {
                best = association.Level;
            }
        }

        return best;
    }

    /// <summary>
    /// Classifies a single variant's clinical actionability to the highest OncoKB therapeutic level over all
    /// its caller-supplied drug associations, under the combined order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4
    /// &gt; R2. Returns <see cref="OncoKbLevel.None"/> when the variant has no leveled association (not
    /// actionable). Source: Chakravarty D et al. (2017); oncokb-annotator README HIGHEST_LEVEL.
    /// </summary>
    /// <param name="variant">Caller-supplied variant actionability evidence.</param>
    /// <returns>The highest combined OncoKB level, or <see cref="OncoKbLevel.None"/>.</returns>
    public static OncoKbLevel ClassifyActionabilityLevel(VariantActionabilityInput variant)
    {
        ArgumentNullException.ThrowIfNull(variant.Associations);
        return HighestLevel(variant.Associations, allowed: null);
    }

    /// <summary>
    /// Assesses the clinical actionability of each variant under the OncoKB therapeutic levels system,
    /// computing the highest sensitivity level (1 &gt; 2 &gt; 3A &gt; 3B &gt; 4), highest resistance level
    /// (R1 &gt; R2), and highest combined level (R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2) from the
    /// caller-supplied drug associations. Output preserves input order, one entry per variant. Source:
    /// Chakravarty D et al. (2017); oncokb-annotator README HIGHEST_LEVEL / HIGHEST_SENSITIVE_LEVEL /
    /// HIGHEST_RESISTANCE_LEVEL.
    /// </summary>
    /// <param name="variants">Caller-supplied variant actionability evidence records.</param>
    /// <returns>One <see cref="ActionabilityAssessment"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ActionabilityAssessment> AssessActionability(
        IEnumerable<VariantActionabilityInput> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var assessments = new List<ActionabilityAssessment>();
        foreach (var variant in variants)
        {
            ArgumentNullException.ThrowIfNull(variant.Associations);

            OncoKbLevel sensitive = HighestLevel(variant.Associations, SensitivityLevels);
            OncoKbLevel resistance = HighestLevel(variant.Associations, ResistanceLevels);
            OncoKbLevel combined = HighestLevel(variant.Associations, allowed: null);

            assessments.Add(new ActionabilityAssessment(variant, sensitive, resistance, combined));
        }

        return assessments;
    }

    /// <summary>
    /// Returns the caller-supplied therapy associations for a variant ordered by descending OncoKB level
    /// (most actionable first) under the combined order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2. Thin
    /// presentation wrapper over the knowledgebase rows; returns an empty list (never null) when there are no
    /// associations. Source: oncokb-annotator README HIGHEST_LEVEL order.
    /// </summary>
    /// <param name="variant">Caller-supplied variant actionability evidence.</param>
    /// <returns>The associations ordered most-actionable first.</returns>
    public static IReadOnlyList<TherapyAssociation> GetTherapyRecommendations(VariantActionabilityInput variant)
    {
        ArgumentNullException.ThrowIfNull(variant.Associations);

        var ordered = new List<TherapyAssociation>(variant.Associations);
        // Descending by combined order: higher enum value = more actionable, so negate the comparison.
        ordered.Sort((x, y) => CompareLevels(y.Level, x.Level));
        return ordered;
    }

    /// <summary>
    /// True when the level is one of OncoKB's "standard care" levels (1, 2, R1) as opposed to the
    /// investigational/hypothetical levels (3A, 3B, 4, R2). Source: OncoKB Curation SOP v3.
    /// </summary>
    /// <param name="level">An OncoKB therapeutic level.</param>
    /// <returns>True for standard-care levels (1, 2, R1).</returns>
    public static bool IsStandardCare(OncoKbLevel level) => StandardCareLevels.Contains(level);

    #endregion

    #region Complex Somatic Rearrangement (Chromothripsis) Classification

    /// <summary>
    /// Minimum number of oscillating copy-number changes (per-segment CN state transitions) required by
    /// the first-pass chromothripsis screen. Source: Magrangeas et al. (2011), Blood 118(3):675–678 — the
    /// lowest first-pass operational cutoff of "10, 20, or 50 oscillating copy number changes" cited by the
    /// Korbel &amp; Campbell (2013) framework (Cell 152:1226–1236; review PMC3861665). Default = 10.
    /// </summary>
    public const int MinOscillatingCopyNumberChanges = 10;

    /// <summary>
    /// Maximum number of distinct copy-number states permitted for a chromothripsis call. Korbel &amp;
    /// Campbell (2013) define the hallmark profile as oscillation between (canonically) two — and at most
    /// two-or-three — copy-number states, in contrast with progressive amplification (many ascending
    /// states). Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236 (criterion B). Default = 3.
    /// </summary>
    public const int MaxChromothripsisCopyNumberStates = 3;

    /// <summary>
    /// Minimum number of clustered intrachromosomal structural variants for an event to be eligible for a
    /// chromothripsis call. Source: Cortés-Ciriano et al. (2020), Nat. Genet. 52:331–341 — focal events
    /// "comprising fewer than six SVs" are excluded. Default = 6.
    /// </summary>
    public const int MinChromothripsisSvBurden = 6;

    /// <summary>
    /// Number of adjacent oscillating copy-number segments at or above which a chromothripsis call is
    /// high-confidence. Source: Cortés-Ciriano et al. (2020) — high-confidence calls "display oscillations
    /// between two states in at least seven adjacent segments". Default = 7.
    /// </summary>
    public const int HighConfidenceOscillatingSegments = 7;

    /// <summary>
    /// Minimum number of adjacent oscillating copy-number segments for a low-confidence chromothripsis
    /// signal. Source: Cortés-Ciriano et al. (2020) — low-confidence calls "involve between four and six
    /// segments". Default = 4 (the band [4, 6] is low-confidence; ≥ 7 is high-confidence).
    /// </summary>
    public const int LowConfidenceOscillatingSegments = 4;

    /// <summary>
    /// Coefficient of variation of inter-breakpoint distances expected under the random-breakpoint null.
    /// Source: Korbel &amp; Campbell (2013) — "the null hypothesis of random breakpoints predicts that the
    /// distance between breakpoints should be distributed exponentially"; the exponential distribution has
    /// CV = 1, so over-dispersion toward many short gaps with few long gaps (clustering) gives CV &gt; 1.
    /// </summary>
    private const double ExponentialNullCoefficientOfVariation = 1.0;

    /// <summary>
    /// Classification of a chromosome's somatic structural-rearrangement profile.
    /// </summary>
    public enum ComplexRearrangementType
    {
        /// <summary>Does not meet the chromothripsis hallmark criteria.</summary>
        NotComplex,

        /// <summary>
        /// Chromothripsis: clustered breakpoints with oscillation between ≤ 3 (canonically 2) copy-number
        /// states and sufficient oscillation/SV burden (Korbel &amp; Campbell 2013; Cortés-Ciriano 2020).
        /// </summary>
        Chromothripsis
    }

    /// <summary>
    /// Confidence tier for a chromothripsis signal, from the number of adjacent oscillating segments.
    /// </summary>
    public enum ChromothripsisConfidence
    {
        /// <summary>Fewer than four adjacent oscillating segments — no chromothripsis signal.</summary>
        None,

        /// <summary>Four to six adjacent oscillating segments — low-confidence (Cortés-Ciriano 2020).</summary>
        Low,

        /// <summary>Seven or more adjacent oscillating segments — high-confidence (Cortés-Ciriano 2020).</summary>
        High
    }

    /// <summary>
    /// Input for <see cref="ClassifyComplexRearrangement"/>: the per-segment copy-number states along one
    /// chromosomal region and the number of clustered intrachromosomal structural variants supporting it.
    /// </summary>
    /// <param name="SegmentCopyNumbers">Per-segment integer copy numbers in genomic order along the region.</param>
    /// <param name="StructuralVariantCount">Number of clustered intrachromosomal SVs in the region.</param>
    public readonly record struct ComplexRearrangementInput(
        IReadOnlyList<int> SegmentCopyNumbers,
        int StructuralVariantCount);

    /// <summary>
    /// Result of complex-rearrangement classification.
    /// </summary>
    /// <param name="Type">The classification (chromothripsis or not).</param>
    /// <param name="Confidence">The confidence tier derived from adjacent oscillating-segment count.</param>
    /// <param name="OscillationCount">Number of per-segment copy-number state transitions.</param>
    /// <param name="OscillatingSegmentCount">Number of segments participating in the oscillation.</param>
    /// <param name="DistinctStateCount">Number of distinct copy-number states in the profile.</param>
    /// <param name="StructuralVariantCount">Clustered intrachromosomal SV burden of the region.</param>
    public readonly record struct ComplexRearrangementResult(
        ComplexRearrangementType Type,
        ChromothripsisConfidence Confidence,
        int OscillationCount,
        int OscillatingSegmentCount,
        int DistinctStateCount,
        int StructuralVariantCount);

    /// <summary>
    /// Summary of a breakpoint-clustering test against the random-breakpoint exponential null.
    /// </summary>
    /// <param name="BreakpointCount">Number of breakpoints provided.</param>
    /// <param name="MeanGap">Mean inter-breakpoint distance.</param>
    /// <param name="CoefficientOfVariation">Standard deviation / mean of inter-breakpoint distances.</param>
    /// <param name="IsClustered">True when CV &gt; 1 (over-dispersed relative to the exponential null).</param>
    public readonly record struct BreakpointClusteringResult(
        int BreakpointCount,
        double MeanGap,
        double CoefficientOfVariation,
        bool IsClustered);

    /// <summary>
    /// Counts the number of oscillating copy-number changes along a region: the number of adjacent segments
    /// whose copy-number state differs from the immediately preceding segment. This is the "oscillating
    /// copy number changes" quantity used by the first-pass chromothripsis screen.
    /// Source: Magrangeas et al. (2011), Blood 118(3):675–678; Korbel &amp; Campbell (2013), Cell 152:1226–1236.
    /// For <c>n</c> segments the count is in [0, n−1]; fewer than two segments yields 0.
    /// </summary>
    /// <param name="segmentCopyNumbers">Per-segment integer copy numbers in genomic order.</param>
    /// <returns>The number of adjacent copy-number state transitions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segmentCopyNumbers"/> is null.</exception>
    public static int CountCopyNumberStateOscillations(IReadOnlyList<int> segmentCopyNumbers)
    {
        ArgumentNullException.ThrowIfNull(segmentCopyNumbers);

        int transitions = 0;
        for (int i = 1; i < segmentCopyNumbers.Count; i++)
        {
            if (segmentCopyNumbers[i] != segmentCopyNumbers[i - 1])
            {
                transitions++;
            }
        }

        return transitions;
    }

    /// <summary>
    /// Tests a set of genomic breakpoint positions for clustering against the random-breakpoint null.
    /// Under the null of uniformly random breakpoints the inter-breakpoint distances are exponentially
    /// distributed, which has a coefficient of variation (CV = sd/mean) of 1; over-dispersion toward many
    /// short gaps with a few long gaps (a tight cluster plus outliers) gives CV &gt; 1, which flags
    /// clustering. Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236 (criterion A, exponential null).
    /// </summary>
    /// <param name="breakpointPositions">Genomic breakpoint coordinates (any order); sorted internally.</param>
    /// <returns>The clustering summary; with fewer than three breakpoints clustering cannot be assessed
    /// (fewer than two gaps), so <see cref="BreakpointClusteringResult.IsClustered"/> is false.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="breakpointPositions"/> is null.</exception>
    public static BreakpointClusteringResult TestBreakpointClustering(IReadOnlyList<long> breakpointPositions)
    {
        ArgumentNullException.ThrowIfNull(breakpointPositions);

        // Need at least two gaps (three breakpoints) to define a CV; otherwise clustering is undefined.
        if (breakpointPositions.Count < 3)
        {
            return new BreakpointClusteringResult(breakpointPositions.Count, 0.0, 0.0, false);
        }

        var sorted = breakpointPositions.OrderBy(p => p).ToArray();
        int gapCount = sorted.Length - 1;
        var gaps = new double[gapCount];
        for (int i = 0; i < gapCount; i++)
        {
            gaps[i] = sorted[i + 1] - sorted[i];
        }

        double mean = gaps.Average();
        if (mean <= 0.0)
        {
            // All breakpoints coincide: degenerate, treat as not assessable.
            return new BreakpointClusteringResult(breakpointPositions.Count, 0.0, 0.0, false);
        }

        double variance = gaps.Sum(g => (g - mean) * (g - mean)) / gapCount;
        double cv = Math.Sqrt(variance) / mean;
        bool clustered = cv > ExponentialNullCoefficientOfVariation;

        return new BreakpointClusteringResult(breakpointPositions.Count, mean, cv, clustered);
    }

    /// <summary>
    /// Classifies a chromosomal region's somatic structural-rearrangement profile as chromothripsis or not,
    /// applying the Korbel &amp; Campbell (2013) hallmark criteria together with the Cortés-Ciriano (2020)
    /// operational thresholds. A region is called <see cref="ComplexRearrangementType.Chromothripsis"/> when
    /// ALL of the following hold: (i) the copy-number profile oscillates between at most
    /// <see cref="MaxChromothripsisCopyNumberStates"/> (canonically 2) distinct states — the two-state
    /// hallmark, excluding progressive amplification; (ii) it has at least
    /// <see cref="MinOscillatingCopyNumberChanges"/> oscillating copy-number changes (first-pass screen);
    /// and (iii) the clustered intrachromosomal SV burden is at least <see cref="MinChromothripsisSvBurden"/>.
    /// The confidence tier is derived independently from the number of adjacent oscillating segments
    /// (≥ <see cref="HighConfidenceOscillatingSegments"/> → High; [<see cref="LowConfidenceOscillatingSegments"/>, 6] → Low; else None).
    /// Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236; Cortés-Ciriano et al. (2020), Nat. Genet. 52:331–341;
    /// Magrangeas et al. (2011), Blood 118(3):675–678.
    /// </summary>
    /// <param name="input">Per-segment copy numbers and the clustered SV burden of the region.</param>
    /// <returns>The classification result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="input"/>.SegmentCopyNumbers is null.</exception>
    public static ComplexRearrangementResult ClassifyComplexRearrangement(ComplexRearrangementInput input)
    {
        ArgumentNullException.ThrowIfNull(input.SegmentCopyNumbers);

        var states = input.SegmentCopyNumbers;
        int oscillations = CountCopyNumberStateOscillations(states);

        // Segments participating in an oscillation: a run of k transitions spans k+1 segments.
        int oscillatingSegments = oscillations > 0 ? oscillations + 1 : 0;

        int distinctStates = states.Count == 0 ? 0 : states.Distinct().Count();

        // Confidence tier from adjacent oscillating-segment count (Cortés-Ciriano 2020).
        ChromothripsisConfidence confidence;
        if (oscillatingSegments >= HighConfidenceOscillatingSegments)
        {
            confidence = ChromothripsisConfidence.High;
        }
        else if (oscillatingSegments >= LowConfidenceOscillatingSegments)
        {
            confidence = ChromothripsisConfidence.Low;
        }
        else
        {
            confidence = ChromothripsisConfidence.None;
        }

        // Chromothripsis hallmark gate: two-state oscillation + ≥10 oscillations + ≥6 clustered SVs.
        bool isChromothripsis =
            distinctStates >= 2 &&
            distinctStates <= MaxChromothripsisCopyNumberStates &&
            oscillations >= MinOscillatingCopyNumberChanges &&
            input.StructuralVariantCount >= MinChromothripsisSvBurden;

        var type = isChromothripsis
            ? ComplexRearrangementType.Chromothripsis
            : ComplexRearrangementType.NotComplex;

        return new ComplexRearrangementResult(
            type,
            confidence,
            oscillations,
            oscillatingSegments,
            distinctStates,
            input.StructuralVariantCount);
    }

    #endregion

    #region Tumor Gene Expression Outlier / Signature Score (ONCO-EXPR-001)

    /// <summary>
    /// Direction of an expression outlier relative to the reference cohort.
    /// </summary>
    public enum ExpressionDirection
    {
        /// <summary>z &gt; +threshold — expression is elevated versus the reference cohort (overexpressed).</summary>
        Over,

        /// <summary>z &lt; −threshold — expression is reduced versus the reference cohort (underexpressed).</summary>
        Under,
    }

    /// <summary>
    /// A single gene whose sample expression is an outlier relative to its reference cohort.
    /// </summary>
    /// <param name="Gene">Gene identifier.</param>
    /// <param name="ZScore">The expression z-score z = (value − μ)/σ of the gene in this sample.</param>
    /// <param name="Direction">Whether the gene is over- or under-expressed (sign of the z-score).</param>
    public readonly record struct ExpressionOutlier(string Gene, double ZScore, ExpressionDirection Direction);

    /// <summary>
    /// Computes the expression z-score of a single sample value relative to a reference cohort:
    /// z = (value − μ) / σ, where μ is the arithmetic mean and σ is the <b>sample</b> standard deviation
    /// (divisor n − 1) of the cohort.
    /// </summary>
    /// <remarks>
    /// Source: cBioPortal mRNA z-score normalization specification
    /// (https://docs.cbioportal.org/z-score-normalization-script/) — z = "(r - mu)/sigma where r is the raw
    /// expression value, and mu and sigma are the mean and standard deviation". The reference implementation
    /// <c>NormalizeExpressionLevels.java</c> (cbioportal-core) computes σ with divisor (n − 1)
    /// (<c>std = std/(double)(v.length-1)</c>), i.e. the sample standard deviation, and aborts when σ = 0.
    /// </remarks>
    /// <param name="value">The sample expression value (raw <c>r</c>) on a normalization scale.</param>
    /// <param name="referenceCohort">Reference expression values for this gene; at least two values required.</param>
    /// <returns>The z-score of <paramref name="value"/> relative to the cohort.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="referenceCohort"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// The cohort has fewer than two values (sample standard deviation is undefined), or the cohort has a
    /// standard deviation of 0 (no defined z-score; mirrors the reference implementation's fatal error).
    /// </exception>
    public static double CalculateExpressionZScore(double value, IReadOnlyList<double> referenceCohort)
    {
        ArgumentNullException.ThrowIfNull(referenceCohort);

        int n = referenceCohort.Count;

        // Sample standard deviation (divisor n − 1) is undefined for fewer than two observations.
        if (n < 2)
        {
            throw new ArgumentException(
                "Reference cohort must contain at least two values to compute a sample standard deviation.",
                nameof(referenceCohort));
        }

        // Zero-spread (all values identical) is the documented σ = 0 error case
        // (NormalizeExpressionLevels.java fatal error). Detect it robustly by the
        // cohort's RANGE rather than from the computed SD: with a constant cohort,
        // accumulating the mean incurs a ≤1-ULP rounding error, which leaves a tiny
        // non-zero Σ(rᵢ − μ)² and an SD on the order of 1e-17 — the `sd == 0` test
        // alone would then be bypassed and emit a spurious ~1e20 z-score. The
        // max == min check is exact for a no-spread cohort and immune to that drift.
        double first = referenceCohort[0];
        bool anyDifferent = false;
        for (int i = 1; i < n; i++)
        {
            if (referenceCohort[i] != first)
            {
                anyDifferent = true;
                break;
            }
        }

        if (!anyDifferent)
        {
            throw new ArgumentException(
                "Cannot compute a z-score relative to a reference cohort with a standard deviation of 0.",
                nameof(referenceCohort));
        }

        double mean = 0.0;
        for (int i = 0; i < n; i++)
        {
            mean += referenceCohort[i];
        }

        mean /= n;

        double sumSquaredDeviations = 0.0;
        for (int i = 0; i < n; i++)
        {
            double d = referenceCohort[i] - mean;
            sumSquaredDeviations += d * d;
        }

        // Sample standard deviation: divisor (n − 1) per NormalizeExpressionLevels.java std().
        double sd = Math.Sqrt(sumSquaredDeviations / (n - 1));

        if (sd == 0.0)
        {
            throw new ArgumentException(
                "Cannot compute a z-score relative to a reference cohort with a standard deviation of 0.",
                nameof(referenceCohort));
        }

        return (value - mean) / sd;
    }

    /// <summary>
    /// Identifies genes whose sample expression is an outlier relative to per-gene reference cohorts, using
    /// the z-score rule z &gt; +threshold (overexpressed) or z &lt; −threshold (underexpressed).
    /// </summary>
    /// <remarks>
    /// Source: cBioPortal FAQ (https://docs.cbioportal.org/user-guide/faq/) — "samples with expression
    /// z-scores &gt;2 or &lt;-2 in any queried genes are considered altered." The default threshold is 2.0
    /// (<see cref="DefaultExpressionOutlierThreshold"/>); the comparison is strict, so |z| = threshold is not
    /// reported as an outlier.
    /// </remarks>
    /// <param name="sampleExpression">The sample's expression value per gene.</param>
    /// <param name="referenceCohorts">Per-gene reference cohort values; must contain every gene in the sample.</param>
    /// <param name="threshold">Absolute z-score threshold (must be positive). Default 2.0.</param>
    /// <returns>The outlier genes, in the iteration order of <paramref name="sampleExpression"/>.</returns>
    /// <exception cref="ArgumentNullException">Either dictionary is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> is not positive.</exception>
    /// <exception cref="ArgumentException">A sampled gene has no reference cohort, or a cohort is degenerate.</exception>
    public static IReadOnlyList<ExpressionOutlier> IdentifyOutlierGenes(
        IReadOnlyDictionary<string, double> sampleExpression,
        IReadOnlyDictionary<string, IReadOnlyList<double>> referenceCohorts,
        double threshold = DefaultExpressionOutlierThreshold)
    {
        ArgumentNullException.ThrowIfNull(sampleExpression);
        ArgumentNullException.ThrowIfNull(referenceCohorts);

        if (threshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(threshold), threshold, "Outlier threshold must be positive.");
        }

        var outliers = new List<ExpressionOutlier>();

        foreach (KeyValuePair<string, double> entry in sampleExpression)
        {
            if (!referenceCohorts.TryGetValue(entry.Key, out IReadOnlyList<double>? cohort))
            {
                throw new ArgumentException(
                    $"No reference cohort supplied for gene '{entry.Key}'.", nameof(referenceCohorts));
            }

            double z = CalculateExpressionZScore(entry.Value, cohort);

            // Strict thresholds: z > +t ⇒ over, z < −t ⇒ under; |z| = t is not an outlier (cBioPortal FAQ).
            if (z > threshold)
            {
                outliers.Add(new ExpressionOutlier(entry.Key, z, ExpressionDirection.Over));
            }
            else if (z < -threshold)
            {
                outliers.Add(new ExpressionOutlier(entry.Key, z, ExpressionDirection.Under));
            }
        }

        return outliers;
    }

    /// <summary>
    /// Computes the combined z-score (gene-signature / pathway activity score) over a set of member-gene
    /// z-scores: a = (Σᵢ zᵢ) / √k, where k is the number of member genes.
    /// </summary>
    /// <remarks>
    /// Source: Lee E. et al. (2008), "Inferring Pathway Activity toward Precise Disease Classification",
    /// PLoS Comput Biol 4(11):e1000217 (https://doi.org/10.1371/journal.pcbi.1000217) — the per-gene
    /// z-scores of a gene set are "averaged into a combined z-score … the square root of the number of member
    /// genes is used in the denominator to stabilize the variance of the mean." Corroborated by the GSVA
    /// "combined z-score" method. Member z-scores are caller-supplied (e.g. from
    /// <see cref="CalculateExpressionZScore"/>); the signature gene set is caller-defined.
    /// </remarks>
    /// <param name="memberZScores">The per-gene z-scores of the signature's member genes (k ≥ 1).</param>
    /// <returns>The combined z-score activity a = (Σ z) / √k.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="memberZScores"/> is null.</exception>
    /// <exception cref="ArgumentException">The signature is empty (k = 0; the score is undefined).</exception>
    public static double CalculateSignatureScore(IReadOnlyList<double> memberZScores)
    {
        ArgumentNullException.ThrowIfNull(memberZScores);

        int k = memberZScores.Count;
        if (k == 0)
        {
            throw new ArgumentException(
                "Signature must contain at least one member gene z-score.", nameof(memberZScores));
        }

        double sum = 0.0;
        for (int i = 0; i < k; i++)
        {
            sum += memberZScores[i];
        }

        return sum / Math.Sqrt(k);
    }

    #endregion
}
