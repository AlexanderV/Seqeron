using System;
using System.Collections.Generic;
using System.Linq;

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
}
