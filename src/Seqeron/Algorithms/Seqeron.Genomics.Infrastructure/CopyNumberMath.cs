namespace Seqeron.Genomics.Infrastructure
{
    /// <summary>
    /// Canonical log2-ratio ↔ absolute copy-number conversions shared by the copy-number analyzers
    /// (oncology CNA calling, structural-variant CNV calling, chromosome-scale copy number).
    /// Ported from CNVkit <c>cnvlib/call.py</c> (<c>_log2_ratio_to_absolute_pure</c>,
    /// <c>_log2_ratio_to_absolute</c>, <c>log2_ratios</c>); see https://github.com/etal/cnvkit.
    /// </summary>
    /// <remarks>
    /// The sex-chromosome bookkeeping of CNVkit (<c>get_as_dframe_and_set_reference_and_expect_copies</c>) is left to
    /// the caller: pass the chromosome's <i>reference</i> copies (its ploidy in the CNVkit reference, e.g. 1 for chrX
    /// of a male reference, ploidy/2 for chrY) and its <i>expected</i> copies (its neutral ploidy in the sample, e.g.
    /// 2 for chrX of a female sample, 0 for chrY of a female sample). For autosomes both equal the sample ploidy.
    /// </remarks>
    public static class CopyNumberMath
    {
        /// <summary>
        /// Lower bound applied to <c>absolute / ploidy</c> before taking the logarithm in
        /// <see cref="AbsoluteToLog2Ratio"/>, so a zero copy number maps to log2(1e-3) ≈ −9.966 instead of −∞.
        /// Source: CNVkit <c>cnvlib/call.py</c> <c>log2_ratios(..., min_abs_val=1e-3)</c>.
        /// </summary>
        public const double DefaultMinAbsoluteRatio = 1e-3;

        /// <summary>
        /// Converts a log2 copy ratio to a continuous absolute copy number for a pure sample:
        /// <c>n = r · 2^v</c>. Source: CNVkit <c>_log2_ratio_to_absolute_pure</c>
        /// (<c>ncopies = ref_copies * 2**log2_ratio</c>).
        /// </summary>
        /// <param name="log2Ratio">log2 copy ratio v; NaN propagates to NaN, −∞ gives 0, +∞ gives +∞.</param>
        /// <param name="referenceCopies">Reference copies r (ploidy of the chromosome in the reference); finite, ≥ 0.</param>
        /// <returns>Continuous (unrounded) absolute copy number.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="referenceCopies"/> is negative, NaN or infinite.</exception>
        public static double Log2RatioToAbsolute(double log2Ratio, double referenceCopies)
        {
            ValidateCopies(referenceCopies, nameof(referenceCopies));
            return referenceCopies * Math.Pow(2.0, log2Ratio);
        }

        /// <summary>
        /// Converts a log2 copy ratio to a continuous absolute tumour copy number, correcting for normal-cell
        /// contamination: for purity <c>p &lt; 1</c>,
        /// <c>n = max(0, (r · 2^v − x · (1 − p)) / p)</c> (NaN stays NaN); for <c>p = 1</c> the pure path
        /// <c>n = r · 2^v</c> (<see cref="Log2RatioToAbsolute(double, double)"/>) is used.
        /// Source: CNVkit <c>cnvlib/call.py</c> <c>_log2_ratio_to_absolute</c> — derived from
        /// <c>2^v = p·n/r + (1−p)·x/r</c>; the clamp at 0 is CNVkit issue #503 (a segment deleted below the
        /// normal-contamination floor would otherwise extrapolate to a negative count).
        /// </summary>
        /// <param name="log2Ratio">Observed log2 copy ratio v.</param>
        /// <param name="referenceCopies">Reference copies r; finite, ≥ 0.</param>
        /// <param name="expectedCopies">Expected (neutral) copies x of the chromosome in the sample; finite, ≥ 0.</param>
        /// <param name="purity">Tumour purity p ∈ (0, 1] (CNVkit <c>call --purity</c> accepts exactly this range).</param>
        /// <returns>Continuous (unrounded) absolute tumour copy number, ≥ 0 unless NaN.</returns>
        /// <exception cref="ArgumentOutOfRangeException">A copy count is negative/non-finite, or purity ∉ (0, 1].</exception>
        public static double Log2RatioToAbsolute(double log2Ratio, double referenceCopies, double expectedCopies, double purity)
        {
            ValidateCopies(referenceCopies, nameof(referenceCopies));
            ValidateCopies(expectedCopies, nameof(expectedCopies));
            ValidatePurity(purity);

            if (purity < 1.0)
            {
                double copies = (referenceCopies * Math.Pow(2.0, log2Ratio) - expectedCopies * (1.0 - purity)) / purity;
                return double.IsNaN(copies) ? copies : Math.Max(0.0, copies);
            }

            return Log2RatioToAbsolute(log2Ratio, referenceCopies);
        }

        /// <summary>
        /// Converts an absolute copy number back to a log2 ratio relative to the sample ploidy:
        /// <c>log2(max(n / ploidy, minAbsoluteRatio))</c> (NaN stays NaN). This is the autosomal form of CNVkit
        /// <c>log2_ratios</c> (without <c>round_to_int</c>); CNVkit additionally adds +1 on chrY (and on chrX for a
        /// haploid-X reference), which callers handle themselves.
        /// </summary>
        /// <param name="absoluteCopies">Absolute copy number n.</param>
        /// <param name="ploidy">Sample ploidy; finite, &gt; 0.</param>
        /// <param name="minAbsoluteRatio">Floor for <c>n / ploidy</c>; default <see cref="DefaultMinAbsoluteRatio"/>.</param>
        /// <returns>The log2 ratio.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> or <paramref name="minAbsoluteRatio"/> is not finite positive.</exception>
        public static double AbsoluteToLog2Ratio(double absoluteCopies, double ploidy, double minAbsoluteRatio = DefaultMinAbsoluteRatio)
        {
            if (!double.IsFinite(ploidy) || ploidy <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy must be a finite positive number.");
            }

            if (!double.IsFinite(minAbsoluteRatio) || minAbsoluteRatio <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(minAbsoluteRatio), minAbsoluteRatio, "The ratio floor must be a finite positive number.");
            }

            // numpy np.maximum propagates NaN; Math.Max does too.
            return Math.Log2(Math.Max(absoluteCopies / ploidy, minAbsoluteRatio));
        }

        private static void ValidateCopies(double copies, string name)
        {
            if (!double.IsFinite(copies) || copies < 0.0)
            {
                throw new ArgumentOutOfRangeException(name, copies, "Copy count must be a finite non-negative number.");
            }
        }

        private static void ValidatePurity(double purity)
        {
            if (!(purity > 0.0 && purity <= 1.0))
            {
                throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
            }
        }
    }
}
