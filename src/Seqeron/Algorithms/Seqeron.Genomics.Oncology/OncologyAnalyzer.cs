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
public static partial class OncologyAnalyzer
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

}
