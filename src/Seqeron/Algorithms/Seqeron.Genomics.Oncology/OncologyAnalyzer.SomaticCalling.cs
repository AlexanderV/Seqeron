using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
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

}
