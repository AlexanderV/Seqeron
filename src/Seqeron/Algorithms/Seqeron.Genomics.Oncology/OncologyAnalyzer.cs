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

    #region Helpers

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
}
