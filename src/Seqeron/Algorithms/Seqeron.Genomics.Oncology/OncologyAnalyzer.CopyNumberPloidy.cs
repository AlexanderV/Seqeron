namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
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
    /// cutoffs verbatim as DEL(0) ≤ −1.1, LOSS(1) ≤ −0.25, GAIN(3) &gt; +0.2, AMP(4) &gt; +0.7
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
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number.</exception>
    public static double Log2RatioToCopyNumber(double log2Ratio, double ploidy = DiploidReferencePloidy)
    {
        CopyNumberMath.ValidatePloidy(ploidy);
        return CopyNumberMath.Log2RatioToAbsolute(log2Ratio, ploidy);
    }

    /// <summary>
    /// Converts a log2 copy ratio to a continuous absolute tumour copy number corrected for tumour purity
    /// (normal-cell contamination). For <c>purity &lt; 1</c>:
    /// <c>n = max(0, (ploidy · 2^log2 − ploidy · (1 − purity)) / purity)</c>; for <c>purity = 1</c> this is exactly
    /// <see cref="Log2RatioToCopyNumber(double, double)"/>. Autosomal form (reference = expected copies = ploidy).
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>_log2_ratio_to_absolute</c> (via
    /// <see cref="CopyNumberMath.Log2RatioToAbsolute(double, double, double, double)"/>).
    /// </summary>
    /// <param name="log2Ratio">Observed log2 copy ratio; NaN propagates to NaN.</param>
    /// <param name="ploidy">Reference (germline) ploidy; finite positive.</param>
    /// <param name="purity">Tumour purity ∈ (0, 1] (CNVkit <c>call --purity</c> range).</param>
    /// <returns>Purity-corrected continuous absolute copy number (≥ 0 unless NaN).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not finite positive, or <paramref name="purity"/> ∉ (0, 1].</exception>
    public static double Log2RatioToCopyNumber(double log2Ratio, double ploidy, double purity)
    {
        CopyNumberMath.ValidatePloidy(ploidy);
        return CopyNumberMath.Log2RatioToAbsolute(log2Ratio, ploidy, ploidy, purity);
    }

    /// <summary>
    /// Rescales an observed log2 ratio for tumour purity as CNVkit <c>do_call</c> does before thresholding:
    /// absolute = <c>_log2_ratio_to_absolute(v, ploidy, ploidy, purity)</c>, then
    /// <c>log2 = log2(max(absolute / ploidy, 1e-3))</c> (<c>log2_ratios</c>, autosomal). Purity 1 (CNVkit skips
    /// rescaling when <c>purity &lt; 1.0</c> is false) and NaN (CNVkit 0.9.14 propagates NaN to the threshold
    /// no-call) return the input unchanged.
    /// </summary>
    private static double RescaleLog2ForPurity(double log2Ratio, double ploidy, double purity)
    {
        CopyNumberMath.ValidatePurity(purity);
        if (purity >= 1.0 || double.IsNaN(log2Ratio))
        {
            return log2Ratio;
        }

        double absolute = CopyNumberMath.Log2RatioToAbsolute(log2Ratio, ploidy, ploidy, purity);
        return CopyNumberMath.AbsoluteToLog2Ratio(absolute, ploidy);
    }

    /// <summary>
    /// Calls an integer copy number from a log2 ratio using CNVkit's hard-threshold method. The copy number
    /// is the index of the first ascending threshold the log2 value is less than or equal to (counting up
    /// from 0); if the log2 value exceeds every threshold, the copy number is <c>ceil(ploidy · 2^log2)</c>.
    /// A NaN log2 ratio is a no-call and returns the neutral reference copy number (ploidy rounded half to
    /// even, as numpy <c>ndarray.round()</c> in CNVkit <c>do_call</c>).
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>absolute_threshold</c> — "Integer values are assigned for
    /// log2 ratio values up to each given threshold value in sequence, counting up from zero. Above the
    /// last threshold value, integer copy numbers are called assuming full purity, rounding up from the
    /// reference copy number."
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio; NaN is a no-call (neutral).</param>
    /// <param name="thresholds">
    /// Exactly four strictly ascending cutoffs partitioning the log2 axis into states 0/1/2/3/4+; when null,
    /// <see cref="DefaultCopyNumberThresholds"/> (−1.1, −0.25, 0.2, 0.7) is used.
    /// </param>
    /// <param name="ploidy">Reference ploidy used both for the neutral no-call and the amplification ceiling.</param>
    /// <returns>The integer copy number (≥ 0).</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number.</exception>
    /// <remarks>
    /// When <c>ceil(ploidy·2^log2)</c> exceeds <see cref="int.MaxValue"/> (log2 ≳ 30 for diploid, or +∞; CNVkit
    /// raises <c>OverflowError</c> for +∞) the returned copy number saturates at <see cref="int.MaxValue"/>
    /// (Amplification), so the result is always ≥ 0.
    /// </remarks>
    public static int CallCopyNumber(
        double log2Ratio,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        double copyNumber = CallCopyNumberUnbounded(log2Ratio, ValidateThresholds(thresholds), ploidy);

        // Above Int32 (log2 ≳ 30 for diploid, or +∞) the CNVkit value ceil(ploidy·2^log2) is not representable
        // by this int-valued API (CNVkit itself raises OverflowError for +∞). The call is saturated explicitly
        // at Int32.MaxValue: the state stays Amplification (CN ≥ 4) and CN ≥ 0 (INV-3) — never a wrapped value.
        return copyNumber >= int.MaxValue ? int.MaxValue : (int)copyNumber;
    }

    /// <summary>
    /// Purity-aware CNVkit threshold call (<c>do_call(method="threshold", purity=…)</c>): when
    /// <paramref name="purity"/> &lt; 1 the log2 ratio is first rescaled for normal-cell contamination
    /// (<c>absolute_clonal</c> → <c>log2_ratios</c>: <c>log2(max(n/ploidy, 1e-3))</c> with
    /// <c>n = max(0, (ploidy·2^v − ploidy·(1−p))/p)</c>), then the hard thresholds of
    /// <see cref="CallCopyNumber(double, IReadOnlyList{double}?, double)"/> are applied to the rescaled value.
    /// <paramref name="purity"/> = 1 is identical to the purity-less overload. Autosomal form.
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>do_call</c>, <c>absolute_dataframe</c>, <c>_log2_ratio_to_absolute</c>,
    /// <c>log2_ratios</c>, <c>absolute_threshold</c>.
    /// </summary>
    /// <param name="log2Ratio">Observed log2 copy ratio; NaN is a no-call (neutral).</param>
    /// <param name="thresholds">Four strictly ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy (default diploid in the purity-less overload).</param>
    /// <param name="purity">Tumour purity ∈ (0, 1].</param>
    /// <returns>The integer copy number (≥ 0), saturated at <see cref="int.MaxValue"/>.</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not finite positive, or <paramref name="purity"/> ∉ (0, 1].</exception>
    public static int CallCopyNumber(
        double log2Ratio,
        IReadOnlyList<double>? thresholds,
        double ploidy,
        double purity)
    {
        var cutoffs = ValidateThresholds(thresholds);
        CopyNumberMath.ValidatePloidy(ploidy);
        return CallCopyNumber(RescaleLog2ForPurity(log2Ratio, ploidy, purity), cutoffs, ploidy);
    }

    /// <summary>
    /// Single implementation of the CNVkit <c>absolute_threshold</c> hard-threshold call, returning the integer
    /// copy number as an integral <see cref="double"/> WITHOUT the Int32 range check (a value above the last
    /// cutoff is <c>ceil(ploidy·2^log2)</c>, which may be arbitrarily large or +∞). Shared by
    /// <see cref="CallCopyNumber"/> (which saturates values beyond Int32) and predicates that only need to know
    /// which state a value falls in, e.g. <see cref="IsHomozygousDeletion"/> (CN == 0) — a huge / +∞ log2 is
    /// a (non-representable) amplification, never CN 0.
    /// Source: CNVkit <c>cnvlib/call.py</c> <c>absolute_threshold</c> + <c>do_call</c> (<c>absolutes.round()</c>).
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio; NaN is a no-call (neutral).</param>
    /// <param name="cutoffs">Already-validated four strictly ascending cutoffs.</param>
    /// <param name="ploidy">Reference ploidy (validated here).</param>
    /// <returns>The integer copy number (≥ 0) as a double; may exceed Int32 or be +∞.</returns>
    private static double CallCopyNumberUnbounded(double log2Ratio, IReadOnlyList<double> cutoffs, double ploidy)
    {
        CopyNumberMath.ValidatePloidy(ploidy);

        if (double.IsNaN(log2Ratio))
        {
            // No-call: CNVkit absolute_threshold stores the neutral reference copy number (ref_copies = ploidy);
            // do_call then takes cn = absolutes.round(), and numpy rounds half to even (round(2.5) = 2).
            return Math.Round(ploidy, MidpointRounding.ToEven);
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
        return Math.Ceiling(Log2RatioToCopyNumber(log2Ratio, ploidy));
    }

    /// <summary>
    /// Classifies a single region's log2 copy ratio into a <see cref="CopyNumberCall"/> carrying the
    /// continuous absolute copy number, the hard-threshold integer copy number, and the discrete
    /// <see cref="CopyNumberState"/>. The state is derived from the integer copy number: 0 → DeepDeletion,
    /// 1 → Loss, 2 → Neutral, 3 → Gain, ≥4 → Amplification. Source: CNVkit <c>absolute_threshold</c>
    /// (DEL(0)/LOSS(1)/neutral(2)/GAIN(3)/AMP(4)); GISTIC2.0 amplitude semantics (Mermel et al. 2011).
    /// </summary>
    /// <param name="log2Ratio">log2 copy ratio; NaN is a no-call (Neutral, CN = ploidy rounded half to even).</param>
    /// <param name="thresholds">Four ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy (default diploid).</param>
    /// <returns>The copy-number call with absolute CN, integer CN, and CNA state.</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number.</exception>
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
    /// Purity-aware <see cref="ClassifyCopyNumber(double, IReadOnlyList{double}?, double)"/>: the integer copy number
    /// is <see cref="CallCopyNumber(double, IReadOnlyList{double}?, double, double)"/> (thresholds applied to the
    /// purity-rescaled log2, CNVkit <c>do_call</c>) and <see cref="CopyNumberCall.AbsoluteCopyNumber"/> is the
    /// purity-corrected absolute copy number (<see cref="Log2RatioToCopyNumber(double, double, double)"/>; ploidy for a
    /// NaN no-call). <see cref="CopyNumberCall.Log2Ratio"/> keeps the observed (input) value.
    /// <paramref name="purity"/> = 1 is identical to the purity-less overload.
    /// </summary>
    /// <param name="log2Ratio">Observed log2 copy ratio; NaN is a no-call.</param>
    /// <param name="thresholds">Four ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy.</param>
    /// <param name="purity">Tumour purity ∈ (0, 1].</param>
    /// <returns>The copy-number call.</returns>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not finite positive, or <paramref name="purity"/> ∉ (0, 1].</exception>
    public static CopyNumberCall ClassifyCopyNumber(
        double log2Ratio,
        IReadOnlyList<double>? thresholds,
        double ploidy,
        double purity)
    {
        int integerCopyNumber = CallCopyNumber(log2Ratio, thresholds, ploidy, purity);
        double absolute = double.IsNaN(log2Ratio) ? ploidy : Log2RatioToCopyNumber(log2Ratio, ploidy, purity);
        return new CopyNumberCall(log2Ratio, absolute, integerCopyNumber, StateFromCopyNumber(integerCopyNumber));
    }

    /// <summary>
    /// Purity-aware <see cref="ClassifyCopyNumbers(IEnumerable{double}, IReadOnlyList{double}?, double)"/>: one
    /// <see cref="ClassifyCopyNumber(double, IReadOnlyList{double}?, double, double)"/> call per input, in input order.
    /// </summary>
    /// <param name="log2Ratios">Per-region observed log2 copy ratios.</param>
    /// <param name="thresholds">Four ascending cutoffs; null uses <see cref="DefaultCopyNumberThresholds"/>.</param>
    /// <param name="ploidy">Reference ploidy.</param>
    /// <param name="purity">Tumour purity ∈ (0, 1].</param>
    /// <returns>One call per input log2 ratio, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="log2Ratios"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="thresholds"/> is not four strictly ascending values.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not finite positive, or <paramref name="purity"/> ∉ (0, 1].</exception>
    public static IReadOnlyList<CopyNumberCall> ClassifyCopyNumbers(
        IEnumerable<double> log2Ratios,
        IReadOnlyList<double>? thresholds,
        double ploidy,
        double purity)
    {
        ArgumentNullException.ThrowIfNull(log2Ratios);
        var cutoffs = ValidateThresholds(thresholds);
        CopyNumberMath.ValidatePloidy(ploidy);
        CopyNumberMath.ValidatePurity(purity);

        var calls = new List<CopyNumberCall>();
        foreach (double log2Ratio in log2Ratios)
        {
            calls.Add(ClassifyCopyNumber(log2Ratio, cutoffs, ploidy, purity));
        }

        return calls;
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
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number.</exception>
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

        /// <summary>
        /// Optional number of markers (probes / SNPs) in the segment. Supply it together with
        /// <see cref="ArmMarkerCount"/> to measure the arm fraction in marker units, as GISTIC2 does by default.
        /// <c>null</c> (default) = absent; the arm fraction is then measured in bp.
        /// </summary>
        public int? MarkerCount { get; init; }

        /// <summary>
        /// Optional total number of markers on the segment's chromosome arm (GISTIC2 <c>band.snp_length</c>: markers
        /// whose position lies in the arm's cytoband span). Supply it together with <see cref="MarkerCount"/>.
        /// </summary>
        public int? ArmMarkerCount { get; init; }

        /// <summary>
        /// Segment length as a fraction of the chromosome arm. In marker units
        /// (<see cref="MarkerCount"/> ÷ <see cref="ArmMarkerCount"/>) when both marker counts are supplied — the GISTIC2
        /// default (<c>make_sample_B.m</c> → <c>normalize_by_arm_length(D,B,cyto,1,…)</c>, <c>norm_type = 1</c> "by number
        /// of snps": <c>(en − st + 1) ./ armlengths_by_snp</c>); otherwise in bp (Length ÷ ArmLength).
        /// </summary>
        public double ArmFraction => MarkerCount is int markers && ArmMarkerCount is int armMarkers
            ? (double)markers / armMarkers
            : (double)Length / ArmLength;
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
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="thresholds"/> is NaN or outside the GISTIC2 ranges
    /// (<c>t_amp</c> ∈ [0, ∞), <c>broad_len_cutoff</c> ∈ [0, 2]).</exception>
    public static bool IsFocalAmplification(
        in CopyNumberArmSegment segment,
        FocalAmplificationThresholds thresholds)
    {
        ValidateFocalThresholds(thresholds);
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
    /// <remarks>
    /// This is the per-event focal filter of the GISTIC2 reference implementation
    /// (<c>snputil/reconstruct_genomes.m</c>, <c>broad_or_focal = 'focal'</c>: event arm-fraction &lt;
    /// <c>broad_len_cutoff</c> AND amplitude vs <c>t_amp</c>). The amplitude test is strict (&gt; <c>t_amp</c>), following
    /// the GISTIC2 documentation ("gain above this positive value") and GISTIC2 <c>gene_calls.m</c>;
    /// <c>reconstruct_genomes.m</c> uses &gt;=, which differs only at exact equality.
    /// <para><b>Per-segment events.</b> GISTIC2 applies this filter to SCNA <i>events</i> produced by its ziggurat
    /// deconstruction (amplitude measured relative to the underlying level; broad levels learned from the cohort). This
    /// overload treats each input segment as one event with amplitude = its <see cref="CopyNumberArmSegment.Log2Ratio"/>,
    /// so raw segments of an arm-level gain interrupted by a focal peak (e.g. 0.5 | 1.5 | 0.5) report the flanks as focal.
    /// For GISTIC2 events (one broad 0.5 + one focal +1.0) use <see cref="DeconstructZiggurat"/> followed by
    /// <see cref="DetectFocalAmplificationEvents"/> (B24 F50–F52).</para>
    /// Segment coordinates are half-open (<c>Length = End − Start</c>). <para><b>Arm-fraction units.</b> GISTIC2 measures
    /// the arm fraction in markers by default (<c>make_sample_B.m</c> → <c>normalize_by_arm_length(D,B,cyto,1,2)</c>,
    /// <c>norm_type = 1</c>: segment markers ÷ markers in the arm's cytoband span). When a segment carries
    /// <see cref="CopyNumberArmSegment.MarkerCount"/> and <see cref="CopyNumberArmSegment.ArmMarkerCount"/>, its
    /// <see cref="CopyNumberArmSegment.ArmFraction"/> is in marker units (e.g. 39/40 = 0.975 &lt; 0.98 is focal even if the
    /// segment covers 99% of the arm in bp); without marker counts it is in bp (Length ÷ ArmLength). Centromere-spanning
    /// events (GISTIC2 <c>ref_length = 2</c>: p-fraction + q-fraction) are outside the arm-anchored model; they are handled
    /// by <see cref="DeconstructZiggurat"/>.</para>
    /// </remarks>
    /// <param name="segments">Arm-anchored copy-number segments. Must not be null.</param>
    /// <param name="thresholds">Amplitude and length cutoffs; null uses <see cref="FocalAmplificationThresholds.Default"/> (GISTIC2 defaults).</param>
    /// <returns>The focal amplifications, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has non-positive arm length or End ≤ Start.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="thresholds"/> is NaN or outside the GISTIC2 ranges
    /// (<c>t_amp</c> ∈ [0, ∞), <c>broad_len_cutoff</c> ∈ [0, 2]).</exception>
    public static IReadOnlyList<CopyNumberArmSegment> DetectFocalAmplifications(
        IEnumerable<CopyNumberArmSegment> segments,
        FocalAmplificationThresholds? thresholds = null)
    {
        ArgumentNullException.ThrowIfNull(segments);
        FocalAmplificationThresholds cutoffs = thresholds ?? FocalAmplificationThresholds.Default;
        ValidateFocalThresholds(cutoffs);

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
    /// Arm-level mapping (see the <see cref="CopyNumberRegion"/> overload for GISTIC2 locus-overlap mapping).
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
        return GenesOnAffectedArms(amplifications, OncogeneArms);
    }

    /// <summary>
    /// Shared arm → gene-panel mapping of <see cref="IdentifyAmplifiedOncogenes(IEnumerable{CopyNumberArmSegment})"/> (ONCO-CNA-002) and
    /// <see cref="IdentifyDeletedTumorSuppressors(IEnumerable{CopyNumberArmSegment})"/> (ONCO-CNA-003): collects the distinct (case-insensitive,
    /// non-empty) arm labels of the affected segments and returns, in panel order, every panel gene resident on one
    /// of them (each gene at most once).
    /// </summary>
    private static List<string> GenesOnAffectedArms(
        IEnumerable<CopyNumberArmSegment> affectedSegments,
        IReadOnlyList<(string Gene, string Arm)> panel)
    {
        var affectedArms = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (CopyNumberArmSegment segment in affectedSegments)
        {
            if (!string.IsNullOrEmpty(segment.Arm))
            {
                affectedArms.Add(segment.Arm);
            }
        }

        var genes = new List<string>();
        foreach ((string gene, string arm) in panel)
        {
            if (affectedArms.Contains(arm))
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

    /// <summary>
    /// A copy-number region (e.g. an amplified or deleted segment, or a GISTIC2 peak) in chromosome coordinates,
    /// 1-based and closed (<c>Start..End</c> inclusive, as in GISTIC2 seg files and <c>genes_at</c>).
    /// </summary>
    /// <param name="Chromosome">Chromosome name ("17", "chr17", "X"); a leading "chr" and letter case are ignored.</param>
    /// <param name="Start">First base of the region (1-based, inclusive).</param>
    /// <param name="End">Last base of the region (inclusive); must satisfy <see cref="End"/> ≥ <see cref="Start"/>.</param>
    public readonly record struct CopyNumberRegion(string Chromosome, long Start, long End);

    /// <summary>
    /// A gene locus in chromosome coordinates, 1-based and closed (<c>Start..End</c> inclusive), as in a GISTIC2
    /// reference-gene table (<c>rg.symb</c>, <c>rg.chrn</c>, <c>rg.start</c>, <c>rg.end</c>).
    /// </summary>
    /// <param name="Symbol">Gene symbol (non-empty).</param>
    /// <param name="Chromosome">Chromosome name ("17", "chr17", "X"); a leading "chr" and letter case are ignored.</param>
    /// <param name="Start">First base of the gene (1-based, inclusive).</param>
    /// <param name="End">Last base of the gene (inclusive); must satisfy <see cref="End"/> ≥ <see cref="Start"/>.</param>
    public readonly record struct GeneLocus(string Symbol, string Chromosome, long Start, long End);

    /// <summary>
    /// GRCh38 loci of the default oncogene panel (same genes and order as the arm-level panel of
    /// <see cref="IdentifyAmplifiedOncogenes(IEnumerable{CopyNumberArmSegment})"/>). Source: the GISTIC2 hg38 reference
    /// gene table (broadinstitute/gistic2 <c>refgenes/Gencode.v22.170324/gencode_genes.tsv</c>, GENCODE v22 gene records,
    /// 1-based closed): ERBB2 chr17:39,687,914–39,730,426; MYC chr8:127,735,434–127,741,434; EGFR
    /// chr7:55,019,021–55,256,620; CCND1 chr11:69,641,087–69,654,474; MDM2 chr12:68,808,172–68,850,686; CDK4
    /// chr12:57,747,727–57,756,013.
    /// </summary>
    public static IReadOnlyList<GeneLocus> DefaultOncogeneLoci { get; } = new GeneLocus[]
    {
        new("ERBB2", "17", 39_687_914, 39_730_426),
        new("MYC", "8", 127_735_434, 127_741_434),
        new("EGFR", "7", 55_019_021, 55_256_620),
        new("CCND1", "11", 69_641_087, 69_654_474),
        new("MDM2", "12", 68_808_172, 68_850_686),
        new("CDK4", "12", 57_747_727, 57_756_013),
    };

    /// <summary>
    /// Locus-overlap gene mapping (GISTIC2 rule): reports, in panel order and each at most once, every panel gene whose
    /// locus overlaps an amplified region. Overlap is GISTIC2 <c>genes_at.m</c> with <c>partial_hits = 1</c> (the
    /// <c>genetables.m</c> default): same chromosome AND <c>gene.start ≤ region.end</c> AND <c>gene.end ≥ region.start</c>
    /// (closed intervals — touching at one base is an overlap). Unlike the arm-level
    /// <see cref="IdentifyAmplifiedOncogenes(IEnumerable{CopyNumberArmSegment})"/>, a focal amplification elsewhere on
    /// 17q does NOT report ERBB2. GISTIC2's "[closest gene]" fallback for gene-less peaks and its widening of peak
    /// boundaries to the flanking markers (<c>genomic_location(…,1)</c>) are not applied: regions are used as given.
    /// </summary>
    /// <param name="amplifiedRegions">Amplified regions (1-based closed coordinates). Must not be null.</param>
    /// <param name="genePanel">Gene loci to test; null uses <see cref="DefaultOncogeneLoci"/> (GRCh38).</param>
    /// <returns>Distinct symbols of panel genes overlapping any region, in panel order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="amplifiedRegions"/> is null.</exception>
    /// <exception cref="ArgumentException">A region or locus has an empty chromosome / symbol or End &lt; Start.</exception>
    public static IReadOnlyList<string> IdentifyAmplifiedOncogenes(
        IEnumerable<CopyNumberRegion> amplifiedRegions,
        IReadOnlyList<GeneLocus>? genePanel = null)
    {
        ArgumentNullException.ThrowIfNull(amplifiedRegions);
        return GenesOverlappingRegions(amplifiedRegions, genePanel ?? DefaultOncogeneLoci);
    }

    /// <summary>
    /// Shared locus-overlap mapping of the <see cref="CopyNumberRegion"/> overloads of
    /// <see cref="IdentifyAmplifiedOncogenes(IEnumerable{CopyNumberRegion}, IReadOnlyList{GeneLocus}?)"/> and
    /// <see cref="IdentifyDeletedTumorSuppressors(IEnumerable{CopyNumberRegion}, IReadOnlyList{GeneLocus}?)"/>: GISTIC2
    /// <c>genes_at.m</c> (<c>partial_hits = 1</c>) closed-interval overlap, panel order, each gene at most once.
    /// </summary>
    private static List<string> GenesOverlappingRegions(
        IEnumerable<CopyNumberRegion> regions,
        IReadOnlyList<GeneLocus> panel)
    {
        foreach (GeneLocus locus in panel)
        {
            if (string.IsNullOrEmpty(locus.Symbol))
            {
                throw new ArgumentException("Gene locus symbols must be non-empty.", nameof(panel));
            }

            ValidateClosedInterval(locus.Chromosome, locus.Start, locus.End, nameof(panel));
        }

        var regionList = new List<(string Chromosome, long Start, long End)>();
        foreach (CopyNumberRegion region in regions)
        {
            ValidateClosedInterval(region.Chromosome, region.Start, region.End, nameof(regions));
            regionList.Add((NormalizeChromosomeName(region.Chromosome), region.Start, region.End));
        }

        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var genes = new List<string>();
        foreach (GeneLocus locus in panel)
        {
            string chromosome = NormalizeChromosomeName(locus.Chromosome);
            foreach ((string regionChromosome, long start, long end) in regionList)
            {
                if (string.Equals(chromosome, regionChromosome, StringComparison.OrdinalIgnoreCase)
                    && locus.Start <= end && locus.End >= start)
                {
                    if (seen.Add(locus.Symbol))
                    {
                        genes.Add(locus.Symbol);
                    }

                    break;
                }
            }
        }

        return genes;
    }

    /// <summary>Validates a 1-based closed interval: non-empty chromosome and End ≥ Start (GISTIC2 <c>genes_at</c> errors on End &lt; Start).</summary>
    private static void ValidateClosedInterval(string chromosome, long start, long end, string paramName)
    {
        if (string.IsNullOrWhiteSpace(chromosome))
        {
            throw new ArgumentException("Chromosome names must be non-empty.", paramName);
        }

        if (end < start)
        {
            throw new ArgumentException(
                $"Interval on '{chromosome}' must have End ≥ Start (got Start={start}, End={end}).", paramName);
        }
    }

    /// <summary>Chromosome key for matching: trims whitespace and a leading "chr" (any case), so "chr17" ≡ "17".</summary>
    private static string NormalizeChromosomeName(string chromosome)
    {
        string name = chromosome.Trim();
        return name.StartsWith("chr", StringComparison.OrdinalIgnoreCase) ? name[3..] : name;
    }

    /// <summary>
    /// Validates focal-amplification thresholds against the ranges enforced by the GISTIC2 reference
    /// implementation (<c>gp_gistic2_from_seg.m</c>: <c>-ta</c> ∈ [0, Inf], <c>-brlen</c> ∈ [0, 2], non-numeric
    /// values rejected). A NaN threshold would otherwise make every comparison false and silently report nothing;
    /// a negative <c>t_amp</c> would call copy-number losses "amplified".
    /// </summary>
    private static void ValidateFocalThresholds(FocalAmplificationThresholds thresholds)
    {
        double tAmp = thresholds.AmplificationLog2Threshold;
        if (double.IsNaN(tAmp) || tAmp < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(thresholds), tAmp,
                "AmplificationLog2Threshold (GISTIC2 t_amp) must be a number in [0, +Infinity).");
        }

        double cutoff = thresholds.BroadLengthCutoff;
        if (double.IsNaN(cutoff) || cutoff < 0 || cutoff > 2)
        {
            throw new ArgumentOutOfRangeException(
                nameof(thresholds), cutoff,
                "BroadLengthCutoff (GISTIC2 broad_len_cutoff) must be a fraction of chromosome arm in [0, 2].");
        }
    }

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

        // Marker-unit arm fraction (GISTIC2 norm_type = 1): both counts or neither; a GISTIC2 segment spans ≥ 1 marker
        // and an arm-anchored segment cannot hold more markers than its arm.
        if (segment.MarkerCount.HasValue != segment.ArmMarkerCount.HasValue)
        {
            throw new ArgumentException(
                $"Segment on '{segment.Arm}' must supply both MarkerCount and ArmMarkerCount, or neither.",
                nameof(segment));
        }

        if (segment.MarkerCount is int markers && segment.ArmMarkerCount is int armMarkers
            && (markers < 1 || armMarkers < markers))
        {
            throw new ArgumentException(
                $"Segment on '{segment.Arm}' must have 1 ≤ MarkerCount ≤ ArmMarkerCount " +
                $"(got MarkerCount={markers}, ArmMarkerCount={armMarkers}).",
                nameof(segment));
        }
    }

    #endregion

    #region GISTIC2 ziggurat deconstruction (ONCO-CNA-002, B24 F50–F52)

    /// <summary>
    /// Marker layout of one chromosome for the GISTIC2 ziggurat deconstruction. The chromosome's markers are numbered
    /// 1..(<see cref="PArmMarkerCount"/> + <see cref="QArmMarkerCount"/>) in genomic order; the first
    /// <see cref="PArmMarkerCount"/> lie on the p arm and the rest on the q arm (GISTIC2 <c>normalize_by_arm_length</c>:
    /// <c>band.snp_length</c> = markers in the arm's cytoband span, a marker is on q when its position is ≥ the q-arm start).
    /// An arm without markers (e.g. the p arm of an acrocentric chromosome) has count 0.
    /// </summary>
    /// <param name="Chromosome">Chromosome label (e.g. "1", "17", "X"); unique within a layout.</param>
    /// <param name="PArmMarkerCount">Number of markers on the p arm (≥ 0).</param>
    /// <param name="QArmMarkerCount">Number of markers on the q arm (≥ 0); P + Q ≥ 1.</param>
    public readonly record struct ZigguratChromosome(string Chromosome, int PArmMarkerCount, int QArmMarkerCount);

    /// <summary>
    /// A copy-number segment of one sample in marker coordinates (1-based, closed, within its chromosome; see
    /// <see cref="ZigguratChromosome"/>). A sample's segments must tile every chromosome of the layout exactly
    /// (GISTIC2 <c>D.dat</c> is a full marker × sample matrix).
    /// </summary>
    /// <param name="Chromosome">Chromosome label matching a <see cref="ZigguratChromosome.Chromosome"/>.</param>
    /// <param name="StartMarker">First marker (1-based within the chromosome).</param>
    /// <param name="EndMarker">Last marker (inclusive), ≥ <paramref name="StartMarker"/>.</param>
    /// <param name="Value">Segment value: log2 ratio (GISTIC2 default input) or copy number − 2.</param>
    public readonly record struct ZigguratSegment(string Chromosome, int StartMarker, int EndMarker, double Value);

    /// <summary>One row of a GISTIC2 B / Z / Q array (columns 1–10 of <c>Qs.m</c>; B uses 1–5 and fraction).</summary>
    internal struct GisticZiggRow
    {
        /// <summary>Column 1 — chromosome index (1-based position in the layout).</summary>
        public int Chromosome;
        /// <summary>Column 2 — start marker (global, 1-based).</summary>
        public int Start;
        /// <summary>Column 3 — end marker (global, 1-based, inclusive).</summary>
        public int End;
        /// <summary>Column 4 — amplitude (B: segment value; Z/Q: event amplitude).</summary>
        public double Amplitude;
        /// <summary>Column 5 — sample index (1-based, column of <c>D.dat</c>).</summary>
        public int Sample;
        /// <summary>Column 6 — starting copy-number level.</summary>
        public double StartLevel;
        /// <summary>Column 7 — ending copy-number level.</summary>
        public double EndLevel;
        /// <summary>Column 8 (B column 6) — length as chromosome-arm fraction (marker units, p + q when spanning).</summary>
        public double Fraction;
        /// <summary>Column 9 — deconstruction (table) score.</summary>
        public double Score;
        /// <summary>Column 10 — broad level of the arm.</summary>
        public double ArmLevel;
    }

    /// <summary>
    /// GISTIC2 marker layout: chromosomes in layout order are GISTIC2 chromosome numbers 1..C, markers are numbered
    /// globally 1..M (chromosome-major), as in <c>D.chrn</c> / <c>D.pos</c>.
    /// </summary>
    internal sealed class GisticMarkerLayout
    {
        private readonly int[] _offset;
        private readonly int[] _p;
        private readonly int[] _q;
        private readonly string[] _names;
        private readonly Dictionary<string, int> _index;

        internal GisticMarkerLayout(IReadOnlyList<ZigguratChromosome> chromosomes)
        {
            ArgumentNullException.ThrowIfNull(chromosomes);
            if (chromosomes.Count == 0)
            {
                throw new ArgumentException("The marker layout must contain at least one chromosome.", nameof(chromosomes));
            }

            int n = chromosomes.Count;
            _offset = new int[n];
            _p = new int[n];
            _q = new int[n];
            _names = new string[n];
            _index = new Dictionary<string, int>(StringComparer.Ordinal);
            long total = 0;
            for (int i = 0; i < n; i++)
            {
                ZigguratChromosome c = chromosomes[i];
                if (string.IsNullOrWhiteSpace(c.Chromosome))
                {
                    throw new ArgumentException("Chromosome labels must be non-empty.", nameof(chromosomes));
                }

                if (c.PArmMarkerCount < 0 || c.QArmMarkerCount < 0 || (long)c.PArmMarkerCount + c.QArmMarkerCount < 1)
                {
                    throw new ArgumentException(
                        $"Chromosome '{c.Chromosome}' must have non-negative arm marker counts with at least one marker " +
                        $"(got p={c.PArmMarkerCount}, q={c.QArmMarkerCount}).", nameof(chromosomes));
                }

                if (!_index.TryAdd(c.Chromosome, i + 1))
                {
                    throw new ArgumentException($"Chromosome '{c.Chromosome}' appears more than once in the layout.", nameof(chromosomes));
                }

                _offset[i] = (int)total;
                _p[i] = c.PArmMarkerCount;
                _q[i] = c.QArmMarkerCount;
                _names[i] = c.Chromosome;
                total += (long)c.PArmMarkerCount + c.QArmMarkerCount;
                if (total > int.MaxValue)
                {
                    throw new ArgumentException("The layout holds more than Int32.MaxValue markers.", nameof(chromosomes));
                }
            }

            MarkerCount = (int)total;
        }

        /// <summary>Number of chromosomes C.</summary>
        internal int ChromosomeCount => _names.Length;

        /// <summary>Total number of markers M.</summary>
        internal int MarkerCount { get; }

        /// <summary>Label of chromosome <paramref name="chr"/> (1-based).</summary>
        internal string Name(int chr) => _names[chr - 1];

        /// <summary>Global index of the marker before the chromosome's first marker.</summary>
        internal int Offset(int chr) => _offset[chr - 1];

        /// <summary>p-arm marker count (<c>armlengths_by_snp(2·chr − 1)</c>).</summary>
        internal int PCount(int chr) => _p[chr - 1];

        /// <summary>q-arm marker count (<c>armlengths_by_snp(2·chr)</c>).</summary>
        internal int QCount(int chr) => _q[chr - 1];

        /// <summary>Global index of the chromosome's last marker (<c>chrnEnd</c>).</summary>
        internal int ChromosomeEnd(int chr) => _offset[chr - 1] + _p[chr - 1] + _q[chr - 1];

        /// <summary>Global index of the first q-arm marker (<c>armstart_by_snp(2·chr)</c>).</summary>
        internal int QStart(int chr) => _offset[chr - 1] + _p[chr - 1] + 1;

        /// <summary>1-based chromosome number of a label, or 0 when absent.</summary>
        internal int IndexOf(string label) => label is not null && _index.TryGetValue(label, out int i) ? i : 0;
    }

    /// <summary>
    /// GISTIC2 <c>normalize_by_arm_length(D,Q,cyto,1,2)</c> (norm_type 1 = marker units, ref_length 2 = p + q sum) for one
    /// segment: a segment starting on q is divided by the q-arm marker count, one ending on p by the p-arm count, and a
    /// centromere-spanning segment is the sum of its p part ÷ p markers and its q part ÷ q markers.
    /// </summary>
    internal static double GisticArmFraction(GisticMarkerLayout layout, int chr, int start, int end)
    {
        int qStart = layout.QStart(chr);
        if (start >= qStart)
        {
            return (double)(end - start + 1) / layout.QCount(chr);
        }

        if (end >= qStart)
        {
            // spans_cent, ref_length 2: Q(:,3) = armstart_by_snp − 1 on p; QQ(:,2) = armstart_by_snp on q; summed.
            double fract = (double)((qStart - 1) - start + 1) / layout.PCount(chr);
            return fract + (double)(end - qStart + 1) / layout.QCount(chr);
        }

        return (double)(end - start + 1) / layout.PCount(chr);
    }

    /// <summary>
    /// GISTIC2 <c>make_sample_B.m</c>: converts one sample's segments into the B array
    /// <c>[chrn st en amp sample fract]</c> — breakpoints where the (transformed) value changes
    /// (<c>find(diff(D.dat(:,idx)) ~= 0)</c>) plus chromosome ends, value = <c>D.dat(bpt)</c>, arm fraction from
    /// <see cref="GisticArmFraction"/>. Validates that the segments tile every layout chromosome exactly.
    /// </summary>
    /// <param name="layout">Marker layout.</param>
    /// <param name="segments">The sample's segments.</param>
    /// <param name="sample">1-based sample index (B column 5).</param>
    /// <param name="transform">Value transform applied per segment before breakpoint detection (cap / log→CN); null = identity.</param>
    internal static List<GisticZiggRow> GisticMakeSampleB(
        GisticMarkerLayout layout,
        IReadOnlyList<ZigguratSegment> segments,
        int sample,
        Func<double, double>? transform = null)
    {
        ArgumentNullException.ThrowIfNull(segments);
        int c = layout.ChromosomeCount;
        var byChromosome = new List<ZigguratSegment>[c];
        for (int i = 0; i < c; i++)
        {
            byChromosome[i] = new List<ZigguratSegment>();
        }

        foreach (ZigguratSegment s in segments)
        {
            int chr = layout.IndexOf(s.Chromosome);
            if (chr == 0)
            {
                throw new ArgumentException(
                    $"Sample {sample - 1}: chromosome '{s.Chromosome}' is not in the marker layout.", nameof(segments));
            }

            if (double.IsNaN(s.Value))
            {
                throw new ArgumentException(
                    $"Sample {sample - 1}: segment {s.Chromosome}:{s.StartMarker}-{s.EndMarker} has a NaN value.", nameof(segments));
            }

            byChromosome[chr - 1].Add(s);
        }

        var b = new List<GisticZiggRow>();
        for (int chr = 1; chr <= c; chr++)
        {
            List<ZigguratSegment> list = byChromosome[chr - 1];
            list.Sort((x, y) => x.StartMarker.CompareTo(y.StartMarker));
            int markers = layout.PCount(chr) + layout.QCount(chr);
            int expectedStart = 1;
            int offset = layout.Offset(chr);
            int runStart = 0;
            double runValue = 0;
            bool open = false;
            foreach (ZigguratSegment s in list)
            {
                if (s.StartMarker != expectedStart || s.EndMarker < s.StartMarker || s.EndMarker > markers)
                {
                    throw new ArgumentException(
                        $"Sample {sample - 1}: segments on chromosome '{layout.Name(chr)}' must tile markers 1..{markers} " +
                        $"without gaps or overlaps (segment {s.StartMarker}-{s.EndMarker}, expected start {expectedStart}).",
                        nameof(segments));
                }

                double value = transform is null ? s.Value : transform(s.Value);
                if (!double.IsFinite(value))
                {
                    throw new ArgumentException(
                        $"Sample {sample - 1}: segment {s.Chromosome}:{s.StartMarker}-{s.EndMarker} has a non-finite " +
                        "copy-number value (set ZigguratOptions.Cap to bound infinite input).", nameof(segments));
                }

                if (open && value != runValue)
                {
                    AddBRow(b, layout, chr, offset + runStart, offset + s.StartMarker - 1, runValue, sample);
                    runStart = s.StartMarker;
                }
                else if (!open)
                {
                    runStart = s.StartMarker;
                    open = true;
                }

                runValue = value;
                expectedStart = s.EndMarker + 1;
            }

            if (expectedStart != markers + 1)
            {
                throw new ArgumentException(
                    $"Sample {sample - 1}: segments on chromosome '{layout.Name(chr)}' must tile markers 1..{markers} " +
                    $"(covered up to {expectedStart - 1}).", nameof(segments));
            }

            AddBRow(b, layout, chr, offset + runStart, offset + markers, runValue, sample);
        }

        return b;
    }

    private static void AddBRow(List<GisticZiggRow> b, GisticMarkerLayout layout, int chr, int start, int end, double value, int sample) =>
        b.Add(new GisticZiggRow
        {
            Chromosome = chr,
            Start = start,
            End = end,
            Amplitude = value,
            Sample = sample,
            Fraction = GisticArmFraction(layout, chr, start, end),
        });

    /// <summary>
    /// GISTIC2 <c>merge_adj_segs.m</c>: merges adjacent rows with equal amplitude (end extended, fractions added).
    /// </summary>
    private static void GisticMergeAdjacentSegments(List<GisticZiggRow> m)
    {
        int i = 0;
        while (i < m.Count - 1)
        {
            if (m[i].Amplitude == m[i + 1].Amplitude)
            {
                GisticZiggRow row = m[i];
                row.End = m[i + 1].End;
                row.Fraction = row.Fraction + m[i + 1].Fraction;
                m[i] = row;
                m.RemoveAt(i + 1);
            }
            else
            {
                i++;
            }
        }
    }

    /// <summary>
    /// GISTIC2 <c>atomic_zigg_deconstruction.m</c>: repeatedly takes the (first) highest segment, records the step from
    /// its higher neighbour (left on ties) as an event <c>[chrn st en amp sample cn_st cn_en fract]</c>, lowers it to that
    /// neighbour and merges equal neighbours; a non-zero residual level is recorded as a final event from 0.
    /// </summary>
    internal static List<GisticZiggRow> GisticAtomicZigguratDeconstruction(List<GisticZiggRow> source)
    {
        var z = new List<GisticZiggRow>();
        if (source.Count == 0)
        {
            return z;
        }

        var bt = new List<GisticZiggRow>(source);
        int sample = bt[0].Sample;
        while (bt.Count > 1 && MaxAmplitude(bt, out int mi) > 0)
        {
            int adj;
            if (mi == 0)
            {
                adj = 1;
            }
            else if (mi == bt.Count - 1)
            {
                adj = mi - 1;
            }
            else
            {
                adj = bt[mi - 1].Amplitude >= bt[mi + 1].Amplitude ? mi - 1 : mi + 1;
            }

            double diff = bt[mi].Amplitude - bt[adj].Amplitude;
            if (diff > 0)
            {
                z.Add(new GisticZiggRow
                {
                    Chromosome = bt[mi].Chromosome,
                    Start = bt[mi].Start,
                    End = bt[mi].End,
                    Fraction = bt[mi].Fraction,
                    StartLevel = bt[adj].Amplitude,
                    EndLevel = bt[mi].Amplitude,
                });
            }

            int kk = Math.Min(adj, mi);
            GisticZiggRow lowered = bt[mi];
            lowered.Amplitude = bt[adj].Amplitude;
            bt[mi] = lowered;
            while (kk < bt.Count - 1 && bt[kk].Amplitude == bt[kk + 1].Amplitude)
            {
                GisticZiggRow row = bt[kk];
                row.End = bt[kk + 1].End;
                row.Fraction = row.Fraction + bt[kk + 1].Fraction;
                bt[kk] = row;
                bt.RemoveAt(kk + 1);
            }
        }

        if (bt[0].Amplitude != 0)
        {
            z.Add(new GisticZiggRow
            {
                Chromosome = bt[0].Chromosome,
                Start = bt[0].Start,
                End = bt[0].End,
                Fraction = bt[0].Fraction,
                StartLevel = 0,
                EndLevel = bt[0].Amplitude,
            });
        }

        for (int i = 0; i < z.Count; i++)
        {
            GisticZiggRow row = z[i];
            row.Amplitude = row.EndLevel - row.StartLevel;
            row.Sample = sample;
            z[i] = row;
        }

        return z;
    }

    /// <summary>MATLAB <c>[m, i] = max(x)</c> over amplitudes: first index of the maximum (no NaN reaches here).</summary>
    private static double MaxAmplitude(List<GisticZiggRow> rows, out int index)
    {
        index = 0;
        double max = rows[0].Amplitude;
        for (int i = 1; i < rows.Count; i++)
        {
            if (rows[i].Amplitude > max)
            {
                max = rows[i].Amplitude;
                index = i;
            }
        }

        return max;
    }

    /// <summary>GISTIC2 <c>prepare_B.m</c>: positive part (Ba) and negated negative part (Bd), each merged.</summary>
    private static (List<GisticZiggRow> Amp, List<GisticZiggRow> Del) GisticPrepareB(List<GisticZiggRow> b)
    {
        var ba = new List<GisticZiggRow>(b.Count);
        var bd = new List<GisticZiggRow>(b.Count);
        foreach (GisticZiggRow row in b)
        {
            GisticZiggRow a = row;
            a.Amplitude = row.Amplitude * (row.Amplitude > 0 ? 1.0 : 0.0); // B.*(B>0): keeps −0 for negatives
            ba.Add(a);
            GisticZiggRow d = row;
            d.Amplitude = -1 * (row.Amplitude * (row.Amplitude < 0 ? 1.0 : 0.0));
            bd.Add(d);
        }

        GisticMergeAdjacentSegments(ba);
        GisticMergeAdjacentSegments(bd);
        return (ba, bd);
    }

    /// <summary>GISTIC2 <c>add_broad_levels_to_zigg.m</c>: shifts cn_st / cn_en by the arm's broad level.</summary>
    private static List<GisticZiggRow> GisticAddBroadLevel(List<GisticZiggRow> z, double broadLevel)
    {
        for (int i = 0; i < z.Count; i++)
        {
            GisticZiggRow row = z[i];
            row.StartLevel += broadLevel;
            row.EndLevel += broadLevel;
            z[i] = row;
        }

        return z;
    }

    /// <summary>
    /// GISTIC2 <c>deconstruct_chr.m</c>: splits a chromosome's B rows at the breakpoint row (<c>en == chr_bpt</c>),
    /// subtracts the p / q broad levels, and deconstructs the positive and negative parts of each arm. Note the
    /// reference quirk kept as is: a breakpoint on the first row (with more rows following) treats all rows as q.
    /// </summary>
    internal static (List<GisticZiggRow> Amp, List<GisticZiggRow> Del) GisticDeconstructChromosome(
        List<GisticZiggRow> b, int chrBreakpoint, double pLevel, double qLevel)
    {
        int bptRow = b.FindIndex(r => r.End == chrBreakpoint);
        if (bptRow < 0)
        {
            throw new InvalidOperationException("Chromosome breakpoint must correspond to segment breakpoint!");
        }

        List<GisticZiggRow> bp = new();
        List<GisticZiggRow> bq = new();
        if (bptRow == b.Count - 1)
        {
            bp = Shift(b, 0, b.Count, pLevel);
        }
        else if (bptRow == 0)
        {
            bq = Shift(b, 0, b.Count, qLevel);
        }
        else
        {
            bp = Shift(b, 0, bptRow + 1, pLevel);
            bq = Shift(b, bptRow + 1, b.Count, qLevel);
        }

        (List<GisticZiggRow> bpa, List<GisticZiggRow> bpd) = GisticPrepareB(bp);
        (List<GisticZiggRow> bqa, List<GisticZiggRow> bqd) = GisticPrepareB(bq);
        List<GisticZiggRow> zap = GisticAddBroadLevel(GisticAtomicZigguratDeconstruction(bpa), pLevel);
        List<GisticZiggRow> zdp = GisticAddBroadLevel(GisticAtomicZigguratDeconstruction(bpd), pLevel);
        List<GisticZiggRow> zaq = GisticAddBroadLevel(GisticAtomicZigguratDeconstruction(bqa), qLevel);
        List<GisticZiggRow> zdq = GisticAddBroadLevel(GisticAtomicZigguratDeconstruction(bqd), qLevel);
        zap.AddRange(zaq);
        zdp.AddRange(zdq);
        return (zap, zdp);

        static List<GisticZiggRow> Shift(List<GisticZiggRow> rows, int from, int to, double level)
        {
            var result = new List<GisticZiggRow>(to - from);
            for (int i = from; i < to; i++)
            {
                GisticZiggRow row = rows[i];
                row.Amplitude -= level;
                result.Add(row);
            }

            return result;
        }
    }

    /// <summary>
    /// GISTIC2 <c>deconstruct_sample.m</c>: ziggurat deconstruction of one sample's B array against given broad levels
    /// (<paramref name="broadLevels"/>[2·ch − 2] = p level, [2·ch − 1] = q level) and per-chromosome breakpoints
    /// (global end marker of the last p segment). Returns the amplification events (ZA) and the deletion events (ZD,
    /// amplitudes positive as returned by the reference; <c>perform_deconstruction</c> negates them).
    /// </summary>
    internal static (List<GisticZiggRow> Amp, List<GisticZiggRow> Del) GisticDeconstructSample(
        List<GisticZiggRow> b, IReadOnlyList<double> broadLevels, IReadOnlyList<int> chromosomeBreakpoints)
    {
        var za = new List<GisticZiggRow>();
        var zd = new List<GisticZiggRow>();
        int maxChr = 0;
        foreach (GisticZiggRow row in b)
        {
            maxChr = Math.Max(maxChr, row.Chromosome);
        }

        for (int ch = 1; ch <= maxChr; ch++)
        {
            List<GisticZiggRow> bt = b.FindAll(r => r.Chromosome == ch);
            (List<GisticZiggRow> a, List<GisticZiggRow> d) = GisticDeconstructChromosome(
                bt, chromosomeBreakpoints[ch - 1], broadLevels[2 * ch - 2], broadLevels[2 * ch - 1]);
            za.AddRange(a);
            zd.AddRange(d);
        }

        return (za, zd);
    }

    /// <summary>
    /// GISTIC2 event length × amplitude log-likelihood table (<c>generate_2d_hists.m</c> output): amplitude bin edges
    /// <c>xamp = −2:0.08:2</c>, length bin edges <c>ylen = 0:0.04:2</c> (arm fractions) and <c>log_hd</c> (51 × 51).
    /// </summary>
    internal sealed class GisticLengthAmplitudeTable
    {
        internal GisticLengthAmplitudeTable(double[] amplitudeEdges, double[] lengthEdges, double[,] logDensity)
        {
            AmplitudeEdges = amplitudeEdges;
            LengthEdges = lengthEdges;
            LogDensity = logDensity;
        }

        /// <summary><c>xamp</c> — amplitude bin edges.</summary>
        internal double[] AmplitudeEdges { get; }

        /// <summary><c>ylen</c> — length (arm fraction) bin edges.</summary>
        internal double[] LengthEdges { get; }

        /// <summary><c>log_hd(xidx, yidx)</c> (0-based here).</summary>
        internal double[,] LogDensity { get; }
    }

    /// <summary>GISTIC2 <c>generate_2d_hists</c> pseudocount argument used by <c>perform_deconstruction.m</c> (<c>.01</c>).</summary>
    private const double GisticHistogramPseudocount = 0.01;

    /// <summary>
    /// GNU Octave colon range <c>start:step:limit</c> for the GISTIC2 bin edges: element i = start + i·step (verified
    /// equal to Octave's <c>-2:.08:2</c> and <c>0:.04:2</c> for all 51 elements).
    /// </summary>
    private static double[] GisticRange(double start, double step, int count)
    {
        var r = new double[count];
        for (int i = 0; i < count; i++)
        {
            r[i] = start + i * step;
        }

        return r;
    }

    /// <summary>
    /// GISTIC2 <c>hist2d.m</c> bin index (1-based) of <paramref name="x"/> over edges <paramref name="r"/>: 1 when
    /// x ≤ r(1), length(r) when x ≥ r(end), otherwise i with r(i) ≤ x &lt; r(i+1).
    /// </summary>
    private static int GisticHistBin(double x, double[] r)
    {
        int bin = 0;
        if (x <= r[0])
        {
            bin = 1;
        }

        if (x >= r[^1])
        {
            bin = r.Length;
        }

        for (int i = 0; i < r.Length - 1; i++)
        {
            if (x >= r[i] && x < r[i + 1])
            {
                bin = i + 1;
            }
        }

        return bin;
    }

    /// <summary>
    /// GISTIC2 <c>generate_2d_hists(QA,QD,[],[],.01,1)</c>: 2-D histogram (<c>hist2d</c> / <c>crosstab_full</c>) of event
    /// amplitude (column 4) × arm fraction (column 8) over QA ∪ QD, plus <c>pseudocount/100·sum(sum(hd))</c> per bin,
    /// normalised by <c>sum(sum(hd1))</c> (column sums first, as MATLAB/Octave <c>sum</c>), then <c>log</c>. The
    /// plotting branch (<c>do_plot</c>) has no effect on the table.
    /// </summary>
    internal static GisticLengthAmplitudeTable GisticGenerate2dHistogram(
        IReadOnlyList<GisticZiggRow> amplifications, IReadOnlyList<GisticZiggRow> deletions)
    {
        double[] xamp = GisticRange(-2, 0.08, 51);
        double[] ylen = GisticRange(0, 0.04, 51);
        int nx = xamp.Length;
        int ny = ylen.Length;
        var hd = new double[nx, ny];
        foreach (GisticZiggRow row in amplifications.Concat(deletions))
        {
            int xi = GisticHistBin(row.Amplitude, xamp);
            int yi = GisticHistBin(row.Fraction, ylen);
            hd[xi - 1, yi - 1] += 1;
        }

        double total = SumColumnsThenTotal(hd);
        double pseudo = GisticHistogramPseudocount / 100 * total;
        var hd1 = new double[nx, ny];
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                hd1[i, j] = hd[i, j] + pseudo;
            }
        }

        double total1 = SumColumnsThenTotal(hd1);
        var logHd = new double[nx, ny];
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                logHd[i, j] = Math.Log(hd1[i, j] / total1);
            }
        }

        return new GisticLengthAmplitudeTable(xamp, ylen, logHd);

        static double SumColumnsThenTotal(double[,] m)
        {
            double total = 0;
            for (int j = 0; j < m.GetLength(1); j++)
            {
                double column = 0;
                for (int i = 0; i < m.GetLength(0); i++)
                {
                    column += m[i, j];
                }

                total += column;
            }

            return total;
        }
    }

    /// <summary>
    /// GISTIC2 <c>score_ziggs_by_table.m</c>: <c>log_hd(xidx, yidx)</c> with xidx = last amplitude edge strictly below
    /// the event amplitude (else 1) and yidx = last length edge strictly below the arm fraction (else 1).
    /// </summary>
    internal static double GisticScoreByTable(in GisticZiggRow row, GisticLengthAmplitudeTable table)
    {
        return table.LogDensity[LastBelow(table.AmplitudeEdges, row.Amplitude), LastBelow(table.LengthEdges, row.Fraction)];

        static int LastBelow(double[] edges, double value)
        {
            for (int i = edges.Length - 1; i >= 0; i--)
            {
                if (edges[i] < value)
                {
                    return i;
                }
            }

            return 0;
        }
    }

    /// <summary>
    /// GISTIC2 <c>ziggurat_on_extremes.m</c>: deconstructs the segments above the minimum level (QAS, levels shifted
    /// back by +min) and below the maximum level (QDS, inverted back: levels −(l − max), amplitude negated), each scored.
    /// </summary>
    private static (List<GisticZiggRow> Amp, List<GisticZiggRow> Del) GisticZigguratOnExtremes(
        List<GisticZiggRow> b, double minLevel, double maxLevel, GisticLengthAmplitudeTable table)
    {
        var ba = new List<GisticZiggRow>(b.Count);
        var bd = new List<GisticZiggRow>(b.Count);
        foreach (GisticZiggRow row in b)
        {
            GisticZiggRow a = row;
            a.Amplitude = row.Amplitude - minLevel;
            ba.Add(a);
            GisticZiggRow d = row;
            d.Amplitude = -1 * (row.Amplitude - maxLevel);
            bd.Add(d);
        }

        GisticMergeAdjacentSegments(ba);
        GisticMergeAdjacentSegments(bd);
        List<GisticZiggRow> qas = GisticAtomicZigguratDeconstruction(ba);
        List<GisticZiggRow> qds = GisticAtomicZigguratDeconstruction(bd);
        for (int i = 0; i < qas.Count; i++)
        {
            GisticZiggRow row = qas[i];
            row.StartLevel += minLevel;
            row.EndLevel += minLevel;
            row.Score = GisticScoreByTable(row, table);
            qas[i] = row;
        }

        for (int i = 0; i < qds.Count; i++)
        {
            GisticZiggRow row = qds[i];
            row.StartLevel = -1 * (row.StartLevel - maxLevel);
            row.EndLevel = -1 * (row.EndLevel - maxLevel);
            row.Amplitude = -1 * row.Amplitude;
            row.Score = GisticScoreByTable(row, table);
            qds[i] = row;
        }

        return (qas, qds);
    }

    /// <summary>
    /// GISTIC2 <c>iterative_ziggurat.m</c>: amplifications must start at least at <paramref name="level"/>, deletions at
    /// most at it; amplitudes recomputed (<c>cn_en − cn_st</c>), events with the wrong sign dropped, changed events rescored.
    /// </summary>
    private static List<GisticZiggRow> GisticIterativeZiggurat(
        List<GisticZiggRow> qas, List<GisticZiggRow> qds, double level, GisticLengthAmplitudeTable table)
    {
        var q = new List<GisticZiggRow>(qas.Count + qds.Count);
        foreach (GisticZiggRow source in qas)
        {
            GisticZiggRow row = source;
            row.StartLevel = source.StartLevel >= level ? source.StartLevel : level; // Octave max(x, y): x >= y ? x : y
            bool changed = source.StartLevel - row.StartLevel != 0;
            row.Amplitude = row.EndLevel - row.StartLevel;
            if (row.Amplitude > 0)
            {
                if (changed)
                {
                    row.Score = GisticScoreByTable(row, table);
                }

                q.Add(row);
            }
        }

        foreach (GisticZiggRow source in qds)
        {
            GisticZiggRow row = source;
            row.StartLevel = source.StartLevel <= level ? source.StartLevel : level; // Octave min(x, y): x <= y ? x : y
            bool changed = source.StartLevel - row.StartLevel != 0;
            row.Amplitude = row.EndLevel - row.StartLevel;
            if (row.Amplitude < 0)
            {
                if (changed)
                {
                    row.Score = GisticScoreByTable(row, table);
                }

                q.Add(row);
            }
        }

        return q;
    }

    /// <summary>Result of <see cref="GisticFindMaxBroadLevel"/> (<c>[broad_level max_Q max_score num_levels]</c>).</summary>
    internal readonly record struct GisticBroadLevelChoice(
        double BroadLevel, List<GisticZiggRow> Events, double Score, int LevelCount);

    /// <summary>
    /// GISTIC2 <c>find_max_broad_level_by_table.m</c>: for each distinct segment value of an arm (ascending), deconstruct
    /// the arm against it (<c>ziggurat_on_extremes</c> + <c>iterative_ziggurat</c>), score the events plus the broad event
    /// (<c>[chrn min(st) max(en) level sample 0 level arm_fract]</c>), and keep the level with the highest total score
    /// (first on ties). One distinct value: that level, a single broad row (score column 0) scored by table. Empty: level
    /// 0, no events, score 0.
    /// </summary>
    internal static GisticBroadLevelChoice GisticFindMaxBroadLevel(
        List<GisticZiggRow> b, GisticLengthAmplitudeTable table, double armFraction)
    {
        // unique(B(:,4)) as in GNU Octave unique.m: stable sort, then the LAST element of each run of equal values is
        // kept (only observable for −0 vs +0).
        var unique = new List<double>(b.Count);
        foreach (double v in b.Select(r => r.Amplitude).OrderBy(v => v))
        {
            if (unique.Count > 0 && unique[^1] == v)
            {
                unique[^1] = v;
            }
            else
            {
                unique.Add(v);
            }
        }

        if (unique.Count > 1)
        {
            (List<GisticZiggRow> qas, List<GisticZiggRow> qds) = GisticZigguratOnExtremes(b, unique[0], unique[^1], table);
            int start = int.MaxValue;
            int end = int.MinValue;
            foreach (GisticZiggRow row in b)
            {
                start = Math.Min(start, row.Start);
                end = Math.Max(end, row.End);
            }

            List<GisticZiggRow>? best = null;
            double bestScore = double.NaN;
            int bestIndex = -1;
            for (int k = 0; k < unique.Count; k++)
            {
                List<GisticZiggRow> q = GisticIterativeZiggurat(qas, qds, unique[k], table);
                var broad = new GisticZiggRow
                {
                    Chromosome = b[0].Chromosome,
                    Start = start,
                    End = end,
                    Amplitude = unique[k],
                    Sample = b[0].Sample,
                    StartLevel = 0,
                    EndLevel = unique[k],
                    Fraction = armFraction,
                };
                broad.Score = GisticScoreByTable(broad, table);
                double sum = 0;
                foreach (GisticZiggRow row in q)
                {
                    sum += row.Score;
                }

                double score = sum + broad.Score;
                if (broad.Amplitude != 0)
                {
                    q.Add(broad);
                }

                // MATLAB [m, i] = max(score): first maximum, NaN ignored (index 1 when all NaN).
                if (bestIndex < 0 || (!double.IsNaN(score) && (double.IsNaN(bestScore) || score > bestScore)))
                {
                    best = q;
                    bestScore = score;
                    bestIndex = k;
                }
            }

            return new GisticBroadLevelChoice(unique[bestIndex], best!, bestScore, unique.Count);
        }

        if (unique.Count == 1)
        {
            var single = new GisticZiggRow
            {
                Chromosome = b[0].Chromosome,
                Start = b[0].Start,
                End = b[0].End,
                Amplitude = unique[0],
                Sample = b[0].Sample,
                StartLevel = 0,
                EndLevel = unique[0],
                Fraction = armFraction,
            };
            double score = GisticScoreByTable(single, table);
            return new GisticBroadLevelChoice(unique[0], new List<GisticZiggRow> { single }, score, 1);
        }

        return new GisticBroadLevelChoice(0, new List<GisticZiggRow>(), 0, 0);
    }

    /// <summary>Kind of a GISTIC2 ziggurat event (the four <c>Qs</c> fields of <c>perform_ziggurat_deconstruction.m</c>).</summary>
    public enum ZigguratEventType
    {
        /// <summary><c>Qs.amp</c> — amplification starting at or above copy-number level 0.</summary>
        Amplification,

        /// <summary><c>Qs.del</c> — deletion starting at or below level 0.</summary>
        Deletion,

        /// <summary><c>Qs.aod</c> — amplification ending at or below level 0 (gain over a deletion).</summary>
        AmplificationOverDeletion,

        /// <summary><c>Qs.doa</c> — deletion ending at or above level 0 (loss over an amplification).</summary>
        DeletionOverAmplification,
    }

    /// <summary>
    /// One GISTIC2 SCNA event from the ziggurat deconstruction (a row of <c>Qs.amp/del/aod/doa</c>; <c>Qs.m</c> columns).
    /// </summary>
    /// <param name="Type">Event kind (which <c>Qs</c> field).</param>
    /// <param name="Sample">0-based index into the input sample list (GISTIC2 column 5 − 1, mapped through <c>Qs.sdesc</c> when the noise filter removed samples).</param>
    /// <param name="Chromosome">Chromosome label (column 1).</param>
    /// <param name="StartMarker">First marker, 1-based within the chromosome (column 2).</param>
    /// <param name="EndMarker">Last marker, inclusive (column 3).</param>
    /// <param name="Amplitude">Event amplitude, positive for all four kinds (column 12, in copy number − 2 units).</param>
    /// <param name="StartLevel">Starting copy-number level (column 6).</param>
    /// <param name="EndLevel">Ending copy-number level (column 7).</param>
    /// <param name="ArmFraction">Length as chromosome-arm fraction in markers; p + q for centromere-spanning events (column 8).</param>
    /// <param name="Score">Deconstruction (table log-likelihood) score (column 9; 0 for single-segment chromosomes and single-level arms).</param>
    /// <param name="ArmLevel">Broad level of the arm the event lies on (column 10).</param>
    public readonly record struct ZigguratEvent(
        ZigguratEventType Type,
        int Sample,
        string Chromosome,
        int StartMarker,
        int EndMarker,
        double Amplitude,
        double StartLevel,
        double EndLevel,
        double ArmFraction,
        double Score,
        double ArmLevel)
    {
        /// <summary>
        /// GISTIC2 <c>reconstruct_genomes.m</c> length test: focal when <see cref="ArmFraction"/> &lt;
        /// <paramref name="broadLengthCutoff"/> (default 0.98), broad otherwise.
        /// </summary>
        public bool IsFocal(double broadLengthCutoff = DefaultBroadLengthCutoff) => ArmFraction < broadLengthCutoff;
    }

    /// <summary>
    /// Options of <see cref="DeconstructZiggurat"/>, with the GISTIC2 defaults (<c>gistic2_param_defaults.m</c>,
    /// <c>run_focal_gistic.m</c> → <c>perform_ziggurat_deconstruction(D,cyto,ziggs,cap)</c> → <c>perform_deconstruction(…,1)</c>).
    /// </summary>
    public sealed record ZigguratOptions
    {
        /// <summary>GISTIC2 defaults: log2 input, cap 1.5, one iteration.</summary>
        public static ZigguratOptions Default { get; } = new();

        /// <summary>
        /// Segment values are log2 ratios (GISTIC2 <c>D.islog</c>, default true) and are converted to copy number − 2 by
        /// <c>2^(x+1) − 2</c>; false = values already in copy number − 2 units.
        /// </summary>
        public bool InputIsLog2 { get; init; } = true;

        /// <summary>
        /// Symmetric cap (GISTIC2 <c>cap</c>, default 1.5) applied before conversion: log2 values are clipped to
        /// [−cap, cap]; copy-number input to [2^(1−cap) − 2, 2^(1+cap) − 2]. null = no cap. Must be &gt; 0.
        /// </summary>
        public double? Cap { get; init; } = 1.5;

        /// <summary>Number of broad-level re-estimation iterations (<c>niters</c>, GISTIC2 uses 1); ≥ 1.</summary>
        public int Iterations { get; init; } = 1;

        /// <summary>
        /// GISTIC2 noise filter (<c>remove_noisy_samples.m</c>, called from <c>perform_ziggurat_deconstruction.m</c>
        /// l.102): samples with more than this many segments are removed before deconstruction (kept when
        /// <c>count ≤ max</c>). Default <see cref="DefaultMaxSegmentsPerSample"/> = 2500, the GISTIC2 pipeline default
        /// (<c>gistic2_param_defaults.m</c> l.136, <c>gp_gistic2_from_seg.m</c> <c>-maxseg</c> l.224; the 500 fallback in
        /// <c>perform_ziggurat_deconstruction.m</c> l.60 applies only when <c>ziggs</c> lacks the field, which the
        /// pipeline always sets). null = no filtering. Must be ≥ 0.
        /// </summary>
        public int? MaxSegmentsPerSample { get; init; } = DefaultMaxSegmentsPerSample;
    }

    /// <summary>GISTIC2 default <c>max_segs_per_sample</c> (<c>gistic2_param_defaults.m</c> l.136): 2500.</summary>
    public const int DefaultMaxSegmentsPerSample = 2500;

    /// <summary>
    /// Result of <see cref="DeconstructZiggurat"/>: the GISTIC2 SCNA events of all samples, in <c>Qs</c> order
    /// (all <c>amp</c>, then <c>del</c>, <c>aod</c>, <c>doa</c>; within a field chromosome-major, then sample, then p before q).
    /// </summary>
    /// <param name="Events">The deconstructed events.</param>
    public sealed record ZigguratDeconstruction(IReadOnlyList<ZigguratEvent> Events)
    {
        /// <summary>
        /// Per input sample, the GISTIC2 segment count used by the noise filter (<c>remove_noisy_samples.m</c>
        /// <c>segment_count</c>): maximal runs of equal capped value along the genome-ordered marker column
        /// (chromosomes in layout order, equal values across a chromosome boundary forming one run) — the
        /// <c>getbpt_counts</c> of the annealed SegArray after <c>cap_vals</c>. Empty when the filter is disabled.
        /// </summary>
        public IReadOnlyList<int> SegmentCounts { get; init; } = Array.Empty<int>();

        /// <summary>
        /// 0-based indices (into the input sample list) of the samples removed by the noise filter
        /// (<see cref="ZigguratOptions.MaxSegmentsPerSample"/>), ascending; empty when none was removed.
        /// </summary>
        public IReadOnlyList<int> RemovedSamples { get; init; } = Array.Empty<int>();

        /// <summary>Events of one sample (0-based index), in <see cref="Events"/> order.</summary>
        public IReadOnlyList<ZigguratEvent> ForSample(int sample) => Events.Where(e => e.Sample == sample).ToList();
    }

    /// <summary>
    /// GISTIC2 ziggurat deconstruction of segmented copy-number data into broad and focal SCNA events, each measured
    /// relative to the level it sits on. Port of <c>perform_ziggurat_deconstruction.m</c> (cap, log2 → copy number − 2)
    /// and <c>perform_deconstruction.m</c>: initial deconstruction against level 0 (<c>make_sample_B</c>,
    /// <c>deconstruct_sample</c>), cohort length × amplitude log-likelihood table (<c>generate_2d_hists</c>), then per
    /// sample and chromosome the p/q breakpoint maximising the table score of the arm deconstructions
    /// (<c>find_max_broad_level_by_table</c>) minus a BIC-style penalty, final event split around level 0
    /// (<c>make_final_Qs</c>), repeated <see cref="ZigguratOptions.Iterations"/> times with the table rebuilt.
    /// </summary>
    /// <remarks>
    /// <para><b>Penalty as coded in GISTIC2.</b> The comment reads "(2k − 1)·ln n where k = #broad levels"; the code uses
    /// n = number of segments on the chromosome and gives 3·ln n when both arm parts are non-empty, but ln(1) = 0 when the
    /// breakpoint is the chromosome's last segment (q part empty: <c>len_bpts = 1</c>). This port follows the code.</para>
    /// <para><b>Cohort model.</b> The table is learned from all supplied samples, so a sample's events depend on the
    /// cohort; n = 1 is allowed (<c>perform_deconstruction</c> runs on any sample set; only the outer GISTIC2 wrapper
    /// requires ≥ 2 samples after noise filtering; here one remaining sample is accepted, consistently with n = 1 input).</para>
    /// <para><b>Noise filter.</b> <c>remove_noisy_samples.m</c> runs after the cap: samples whose segment count
    /// (<see cref="ZigguratDeconstruction.SegmentCounts"/>) exceeds <see cref="ZigguratOptions.MaxSegmentsPerSample"/>
    /// (default 2500) are dropped and listed in <see cref="ZigguratDeconstruction.RemovedSamples"/>; the cohort table is
    /// learned from the kept samples only. The count follows the default SegArray branch (<c>use_segarray = 1</c>:
    /// <c>getbpt_counts</c> = number of segments); GISTIC2's uncompressed branch counts <c>diff ~= 0</c> breakpoints
    /// (one less) — Octave-checked, not followed. All samples removed ⇒ <see cref="ArgumentException"/>
    /// (GISTIC2 <c>all_data_removed</c>). <see cref="ZigguratEvent.Sample"/> always indexes the input list (GISTIC2 maps
    /// <c>Qs</c> column 5 through the filtered <c>Qs.sdesc</c>).</para>
    /// <para>When no sample has any non-zero segment the result is empty (every chromosome is one zero segment; GISTIC2 returns
    /// empty event arrays too — Octave-checked).</para>
    /// <para><b>Units.</b> Amplitudes and levels are in copy number − 2 units (GISTIC2 <c>Qs</c> columns 4/6/7/12); the
    /// arm fraction is in markers (<c>norm_type = 1</c>), p + q for centromere-spanning events.</para>
    /// </remarks>
    /// <param name="chromosomes">Marker layout (GISTIC2 chromosome order = list order; every chromosome present in every sample).</param>
    /// <param name="samples">Per-sample segments tiling every layout chromosome (<see cref="ZigguratSegment"/>).</param>
    /// <param name="options">Options; null = <see cref="ZigguratOptions.Default"/>.</param>
    /// <returns>The events of all samples.</returns>
    /// <exception cref="ArgumentNullException">A required argument or a sample list is null.</exception>
    /// <exception cref="ArgumentException">Invalid layout, no samples, segments that do not tile a chromosome, NaN values,
    /// a non-finite converted value (e.g. +∞ log2 with no cap), or every sample removed by the noise filter.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Cap ≤ 0 or NaN, Iterations &lt; 1, or MaxSegmentsPerSample &lt; 0.</exception>
    public static ZigguratDeconstruction DeconstructZiggurat(
        IReadOnlyList<ZigguratChromosome> chromosomes,
        IReadOnlyList<IReadOnlyList<ZigguratSegment>> samples,
        ZigguratOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(chromosomes);
        ArgumentNullException.ThrowIfNull(samples);
        ZigguratOptions opts = options ?? ZigguratOptions.Default;
        if (opts.Cap is double capValue && (double.IsNaN(capValue) || capValue <= 0))
        {
            throw new ArgumentOutOfRangeException(nameof(options), capValue, "ZigguratOptions.Cap must be a positive number (or null).");
        }

        if (opts.Iterations < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(options), opts.Iterations, "ZigguratOptions.Iterations must be at least 1.");
        }

        if (opts.MaxSegmentsPerSample is int maxSegs && maxSegs < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(options), maxSegs, "ZigguratOptions.MaxSegmentsPerSample must be non-negative (or null).");
        }

        var layout = new GisticMarkerLayout(chromosomes);
        if (samples.Count == 0)
        {
            throw new ArgumentException("At least one sample is required.", nameof(samples));
        }

        Func<double, double> transform = GisticInputTransform(opts);
        var all = new List<GisticZiggRow>[samples.Count];
        for (int j = 0; j < samples.Count; j++)
        {
            IReadOnlyList<ZigguratSegment> segments = samples[j]
                ?? throw new ArgumentNullException(nameof(samples), $"Sample {j} is null.");
            all[j] = GisticMakeSampleB(layout, segments, j + 1, transform);
        }

        // remove_noisy_samples.m: keepers = num_bpts <= max_segs_per_sample (after the cap, perform_ziggurat_deconstruction l.102).
        int[] segmentCounts = Array.Empty<int>();
        var kept = new List<int>(samples.Count);
        var removed = new List<int>();
        if (opts.MaxSegmentsPerSample is int maxSegments)
        {
            Func<double, double> cap = GisticCap(opts);
            segmentCounts = new int[samples.Count];
            for (int j = 0; j < samples.Count; j++)
            {
                segmentCounts[j] = GisticSegmentCount(layout, samples[j], cap);
                (segmentCounts[j] <= maxSegments ? kept : removed).Add(j);
            }

            if (kept.Count == 0)
            {
                throw new ArgumentException(
                    $"All samples were removed by noise filtering (more than {maxSegments} segments each; GISTIC2 all_data_removed).",
                    nameof(samples));
            }
        }
        else
        {
            kept.AddRange(Enumerable.Range(0, samples.Count));
        }

        var b = new List<GisticZiggRow>[kept.Count];
        for (int k = 0; k < kept.Count; k++)
        {
            b[k] = all[kept[k]];
            if (removed.Count > 0)
            {
                for (int r = 0; r < b[k].Count; r++)
                {
                    GisticZiggRow row = b[k][r];
                    row.Sample = k + 1;
                    b[k][r] = row;
                }
            }
        }

        List<GisticZiggRow> q = GisticPerformDeconstruction(layout, b, opts.Iterations, out bool any);
        if (!any)
        {
            return new ZigguratDeconstruction(Array.Empty<ZigguratEvent>())
            {
                SegmentCounts = segmentCounts,
                RemovedSamples = removed,
            };
        }

        var (qa, qd, qaod, qdoa) = GisticMakeFinalQs(q);
        var events = new List<ZigguratEvent>(qa.Count + qd.Count + qaod.Count + qdoa.Count);
        AddEvents(qa, ZigguratEventType.Amplification, 1);
        AddEvents(qd, ZigguratEventType.Deletion, -1);
        AddEvents(qaod, ZigguratEventType.AmplificationOverDeletion, 1);
        AddEvents(qdoa, ZigguratEventType.DeletionOverAmplification, -1);
        return new ZigguratDeconstruction(events) { SegmentCounts = segmentCounts, RemovedSamples = removed };

        void AddEvents(List<GisticZiggRow> rows, ZigguratEventType type, int sign)
        {
            foreach (GisticZiggRow r in rows)
            {
                int offset = layout.Offset(r.Chromosome);
                events.Add(new ZigguratEvent(
                    type, kept[r.Sample - 1], layout.Name(r.Chromosome), r.Start - offset, r.End - offset,
                    sign * r.Amplitude, r.StartLevel, r.EndLevel, r.Fraction, r.Score, r.ArmLevel));
            }
        }
    }

    /// <summary>
    /// Focal amplifications among deconstructed GISTIC2 events — <c>reconstruct_genomes.m</c> with
    /// <c>broad_or_focal = 'focal'</c> on the <c>amp</c> and <c>aod</c> fields: <c>Q(:,8) &lt; broad_len_cutoff</c> and
    /// <c>Q(:,12) ≥ t_amp</c> (≥, as in <c>reconstruct_genomes.m</c>; the amplitude is in the deconstruction's copy
    /// number − 2 units, as GISTIC2 compares it). For 1p 0.5 | 1.5 | 0.5 this yields the +1.0 step only; the 0.5
    /// underlying level is one broad event (whole-arm fraction ≥ 0.98).
    /// </summary>
    /// <param name="deconstruction">Result of <see cref="DeconstructZiggurat"/>.</param>
    /// <param name="thresholds">Amplitude and length cutoffs; null = GISTIC2 defaults (t_amp 0.1, broad_len_cutoff 0.98).</param>
    /// <returns>Focal amplification events (<c>amp</c> events first, then <c>aod</c>), in event order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="deconstruction"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="thresholds"/> is NaN or outside the GISTIC2 ranges.</exception>
    public static IReadOnlyList<ZigguratEvent> DetectFocalAmplificationEvents(
        ZigguratDeconstruction deconstruction,
        FocalAmplificationThresholds? thresholds = null)
    {
        ArgumentNullException.ThrowIfNull(deconstruction);
        FocalAmplificationThresholds cutoffs = thresholds ?? FocalAmplificationThresholds.Default;
        ValidateFocalThresholds(cutoffs);
        return deconstruction.Events
            .Where(e => e.Type is ZigguratEventType.Amplification or ZigguratEventType.AmplificationOverDeletion
                        && e.ArmFraction < cutoffs.BroadLengthCutoff
                        && e.Amplitude >= cutoffs.AmplificationLog2Threshold)
            .ToList();
    }

    /// <summary>
    /// <c>perform_ziggurat_deconstruction.m</c> value preparation: cap (<c>D.dat(D.dat &gt; maxcap) = maxcap</c>, then the
    /// minimum), then log2 → copy number − 2 (<c>2.^(D.dat+1)−2</c>); for copy-number input the cap is
    /// <c>2^(1±cap) − 2</c>.
    /// </summary>
    private static Func<double, double> GisticInputTransform(ZigguratOptions options)
    {
        Func<double, double> cap = GisticCap(options);
        bool log = options.InputIsLog2;
        return value =>
        {
            value = cap(value);
            return log ? Math.Pow(2, value + 1) - 2 : value;
        };
    }

    /// <summary>
    /// <c>perform_ziggurat_deconstruction.m</c> cap only (input units): <c>D.dat(D.dat &gt; maxcap) = maxcap</c>, then
    /// the minimum; identity when <see cref="ZigguratOptions.Cap"/> is null.
    /// </summary>
    private static Func<double, double> GisticCap(ZigguratOptions options)
    {
        double maxCap = double.PositiveInfinity;
        double minCap = double.NegativeInfinity;
        if (options.Cap is double cap)
        {
            maxCap = cap;
            minCap = -cap;
            if (!options.InputIsLog2)
            {
                maxCap = Math.Pow(2, 1 + maxCap) - 2;
                minCap = Math.Pow(2, 1 + minCap) - 2;
            }
        }

        return value =>
        {
            if (value > maxCap)
            {
                value = maxCap;
            }

            if (value < minCap)
            {
                value = minCap;
            }

            return value;
        };
    }

    /// <summary>
    /// GISTIC2 <c>remove_noisy_samples.m</c> segment count of one (already validated) sample: <c>cap_vals</c> anneals the
    /// SegArray column (<c>anneal.m</c>: adjacent equal values joined, <c>diff(vals) ~= 0</c>), then
    /// <c>getbpt_counts</c> = number of stored segments — i.e. 1 + the number of value changes between consecutive
    /// segments of the genome-ordered column (layout chromosome order, chromosome boundaries not kept).
    /// </summary>
    private static int GisticSegmentCount(GisticMarkerLayout layout, IReadOnlyList<ZigguratSegment> segments, Func<double, double> cap)
    {
        var ordered = segments
            .OrderBy(s => layout.IndexOf(s.Chromosome))
            .ThenBy(s => s.StartMarker)
            .ToList();
        int count = 1;
        for (int i = 1; i < ordered.Count; i++)
        {
            if (cap(ordered[i].Value) != cap(ordered[i - 1].Value))
            {
                count++;
            }
        }

        return count;
    }

    /// <summary>
    /// GISTIC2 <c>perform_deconstruction.m</c> (after <c>make_sample_B</c>): returns the final Q array (before
    /// <c>make_final_Qs</c>); <paramref name="anyEvent"/> is false when the initial deconstruction has no events.
    /// </summary>
    internal static List<GisticZiggRow> GisticPerformDeconstruction(
        GisticMarkerLayout layout, IReadOnlyList<List<GisticZiggRow>> b, int iterations, out bool anyEvent)
    {
        int nChr = layout.ChromosomeCount;
        var chrnEnd = new int[nChr];
        for (int ch = 1; ch <= nChr; ch++)
        {
            chrnEnd[ch - 1] = layout.ChromosomeEnd(ch);
        }

        var zeroLevels = new double[2 * nChr];
        var qa = new List<GisticZiggRow>();
        var qd = new List<GisticZiggRow>();
        foreach (List<GisticZiggRow> sample in b)
        {
            (List<GisticZiggRow> za, List<GisticZiggRow> zd) = GisticDeconstructSample(sample, zeroLevels, chrnEnd);
            qa.AddRange(za);
            foreach (GisticZiggRow row in zd)
            {
                GisticZiggRow negated = row;
                negated.Amplitude = -1 * row.Amplitude;
                qd.Add(negated);
            }
        }

        anyEvent = qa.Count + qd.Count > 0;
        var q = new List<GisticZiggRow>();
        if (!anyEvent)
        {
            return q;
        }

        GisticLengthAmplitudeTable table = GisticGenerate2dHistogram(qa, qd);
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            var perChromosome = new List<GisticZiggRow>[nChr, b.Count];
            for (int j = 0; j < b.Count; j++)
            {
                for (int ch = 1; ch <= nChr; ch++)
                {
                    perChromosome[ch - 1, j] = GisticDeconstructChromosomeByTable(b[j].FindAll(r => r.Chromosome == ch), table);
                }
            }

            q = new List<GisticZiggRow>();
            for (int ch = 0; ch < nChr; ch++)
            {
                for (int j = 0; j < b.Count; j++)
                {
                    q.AddRange(perChromosome[ch, j]);
                }
            }

            if (iteration < iterations - 1)
            {
                var (finalA, finalD, _, _) = GisticMakeFinalQs(q);
                table = GisticGenerate2dHistogram(finalA, finalD);
            }
        }

        return q;
    }

    /// <summary>
    /// The per-chromosome body of the <c>perform_deconstruction.m</c> iteration: a single segment is one broad event
    /// (<c>[Bt(1:5) 0 amp sum(fract) 0 amp]</c>); otherwise every segment end i is tried as the p/q breakpoint, each arm
    /// part gets its table-optimal broad level, and the breakpoint with the highest
    /// <c>p_score + q_score − penalty</c> wins (first on ties), penalty = <c>ln(len_bpts)</c> when an arm part has no
    /// levels, else <c>3·ln(len_bpts)</c>, with <c>len_bpts = 1</c> when the q part is empty and the segment count otherwise.
    /// </summary>
    private static List<GisticZiggRow> GisticDeconstructChromosomeByTable(List<GisticZiggRow> bt, GisticLengthAmplitudeTable table)
    {
        int n = bt.Count;
        if (n == 1)
        {
            GisticZiggRow only = bt[0];
            return new List<GisticZiggRow>
            {
                new()
                {
                    Chromosome = only.Chromosome,
                    Start = only.Start,
                    End = only.End,
                    Amplitude = only.Amplitude,
                    Sample = only.Sample,
                    StartLevel = 0,
                    EndLevel = only.Amplitude,
                    Fraction = 0 + only.Fraction,
                    Score = 0,
                    ArmLevel = only.Amplitude,
                },
            };
        }

        List<GisticZiggRow>? best = null;
        double bestScore = double.NaN;
        for (int i = 1; i <= n; i++)
        {
            double pFract = 0;
            for (int k = 0; k < i; k++)
            {
                pFract += bt[k].Fraction;
            }

            double qFract = 0;
            for (int k = i; k < n; k++)
            {
                qFract += bt[k].Fraction;
            }

            GisticBroadLevelChoice p = GisticFindMaxBroadLevel(bt.GetRange(0, i), table, pFract);
            GisticBroadLevelChoice qChoice = GisticFindMaxBroadLevel(bt.GetRange(i, n - i), table, qFract);
            double ziggScore = p.Score + qChoice.Score;
            int lenBpts = qChoice.LevelCount == 0 ? 1 : n;
            double penalty = p.LevelCount == 0 || qChoice.LevelCount == 0
                ? Math.Log(lenBpts)
                : 3 * Math.Log(lenBpts);
            double score = ziggScore - penalty;

            if (best is null || (!double.IsNaN(score) && (double.IsNaN(bestScore) || score > bestScore)))
            {
                best = WithArmLevel(p.Events, p.BroadLevel);
                best.AddRange(WithArmLevel(qChoice.Events, qChoice.BroadLevel));
                bestScore = score;
            }
        }

        return best!;

        static List<GisticZiggRow> WithArmLevel(List<GisticZiggRow> rows, double level)
        {
            var result = new List<GisticZiggRow>(rows.Count);
            foreach (GisticZiggRow row in rows)
            {
                GisticZiggRow r = row;
                r.ArmLevel = level;
                result.Add(r);
            }

            return result;
        }
    }

    /// <summary>
    /// GISTIC2 <c>make_final_Qs.m</c>: QA = amp &gt; 0 ∧ cn_st ≥ 0; QD = amp &lt; 0 ∧ cn_st ≤ 0; QAOD = amp &gt; 0 ∧
    /// cn_en ≤ 0; QDOA = amp &lt; 0 ∧ cn_en ≥ 0; events crossing level 0 (cn_st·cn_en &lt; 0) are split at 0 and the
    /// parts appended (to QA/QAOD resp. QD/QDOA).
    /// </summary>
    internal static (List<GisticZiggRow> Amp, List<GisticZiggRow> Del, List<GisticZiggRow> AmpOverDel, List<GisticZiggRow> DelOverAmp)
        GisticMakeFinalQs(List<GisticZiggRow> q)
    {
        var qa = q.FindAll(r => r.Amplitude > 0 && r.StartLevel >= 0);
        var qd = q.FindAll(r => r.Amplitude < 0 && r.StartLevel <= 0);
        var qaod = q.FindAll(r => r.Amplitude > 0 && r.EndLevel <= 0);
        var qdoa = q.FindAll(r => r.Amplitude < 0 && r.EndLevel >= 0);
        var newA = new List<GisticZiggRow>();
        var newAod = new List<GisticZiggRow>();
        var newD = new List<GisticZiggRow>();
        var newDoa = new List<GisticZiggRow>();
        foreach (GisticZiggRow cur in q)
        {
            if (cur.StartLevel * cur.EndLevel < 0 && (cur.Amplitude > 0 || cur.Amplitude < 0))
            {
                GisticZiggRow over = cur;   // [cur(1:6) 0 cur(8:10)]
                over.EndLevel = 0;
                over.Amplitude = over.EndLevel - over.StartLevel;
                GisticZiggRow pure = cur;   // [cur(1:5) 0 cur(7:10)]
                pure.StartLevel = 0;
                pure.Amplitude = pure.EndLevel - pure.StartLevel;
                if (cur.Amplitude > 0)
                {
                    newA.Add(pure);
                    newAod.Add(over);
                }
                else
                {
                    newD.Add(pure);
                    newDoa.Add(over);
                }
            }
        }

        qa.AddRange(newA);
        qd.AddRange(newD);
        qaod.AddRange(newAod);
        qdoa.AddRange(newDoa);
        return (qa, qd, qaod, qdoa);
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
    /// <remarks>
    /// Never throws for an extreme log2: a huge finite or +∞ log2 lies above the last cutoff, so its CNVkit copy
    /// number is <c>ceil(ploidy·2^log2)</c> ≥ 4 (an amplification) and the predicate is <c>false</c>, even when
    /// that copy number is not representable as Int32 (where <see cref="CallCopyNumber"/> saturates). −∞ is ≤ every
    /// cutoff ⇒ CN 0 ⇒ <c>true</c>. NaN is a neutral no-call (CN = ploidy rounded half to even) ⇒ <c>false</c>.
    /// </remarks>
    /// <exception cref="ArgumentException"><paramref name="segment"/> has non-positive arm length or End ≤ Start; or invalid thresholds.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number.</exception>
    public static bool IsHomozygousDeletion(
        in CopyNumberArmSegment segment,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        ValidateArmSegment(segment);

        // Only "is the integer CN 0?" is needed, so use the unbounded CNVkit call: a huge / +∞ log2 is above the
        // last cutoff (CN = ceil(ploidy·2^log2) ≥ 4), i.e. an amplification — never a homozygous deletion.
        return CallCopyNumberUnbounded(segment.Log2Ratio, ValidateThresholds(thresholds), ploidy)
            == HomozygousDeletionCopyNumber;
    }

    /// <summary>
    /// Purity-aware <see cref="IsHomozygousDeletion(in CopyNumberArmSegment, IReadOnlyList{double}?, double)"/>: the
    /// segment's log2 ratio is rescaled for tumour purity exactly as CNVkit <c>do_call</c> does
    /// (<c>_log2_ratio_to_absolute</c> → <c>log2_ratios</c>) before the hard-threshold call; the segment is a
    /// homozygous deletion when that integer copy number is 0. With normal contamination an observed log2 of −1.0
    /// (CN 1 when pure) is CN 0 at purity 0.7 (absolute 0.571 → rescaled log2 −1.807 ≤ −1.1).
    /// <paramref name="purity"/> = 1 is identical to the purity-less overload.
    /// </summary>
    /// <param name="segment">The arm-anchored copy-number segment.</param>
    /// <param name="thresholds">Four strictly ascending log2 cutoffs; null uses CNVkit defaults.</param>
    /// <param name="ploidy">Reference (germline) ploidy.</param>
    /// <param name="purity">Tumour purity ∈ (0, 1].</param>
    /// <returns><c>true</c> when the purity-corrected integer copy number is 0.</returns>
    /// <exception cref="ArgumentException"><paramref name="segment"/> has non-positive arm length or End ≤ Start; or invalid thresholds.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not finite positive, or <paramref name="purity"/> ∉ (0, 1].</exception>
    public static bool IsHomozygousDeletion(
        in CopyNumberArmSegment segment,
        IReadOnlyList<double>? thresholds,
        double ploidy,
        double purity)
    {
        ValidateArmSegment(segment);
        var cutoffs = ValidateThresholds(thresholds);
        CopyNumberMath.ValidatePloidy(ploidy);
        return CallCopyNumberUnbounded(RescaleLog2ForPurity(segment.Log2Ratio, ploidy, purity), cutoffs, ploidy)
            == HomozygousDeletionCopyNumber;
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
    /// <exception cref="ArgumentException">A segment has non-positive arm length or End ≤ Start; or invalid thresholds (checked even for empty input).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> is not a finite positive number (checked even for empty input).</exception>
    public static IReadOnlyList<CopyNumberArmSegment> DetectHomozygousDeletions(
        IEnumerable<CopyNumberArmSegment> segments,
        IReadOnlyList<double>? thresholds = null,
        double ploidy = DiploidReferencePloidy)
    {
        ArgumentNullException.ThrowIfNull(segments);

        // Validate the calling parameters eagerly (once), so malformed thresholds / ploidy are rejected even for
        // an empty segment list — consistent with DetectFocalAmplifications and ClassifyCopyNumbers.
        var cutoffs = ValidateThresholds(thresholds);
        CopyNumberMath.ValidatePloidy(ploidy);

        var result = new List<CopyNumberArmSegment>();
        foreach (CopyNumberArmSegment segment in segments)
        {
            if (IsHomozygousDeletion(segment, cutoffs, ploidy))
            {
                result.Add(segment);
            }
        }

        return result;
    }

    /// <summary>
    /// Purity-aware <see cref="DetectHomozygousDeletions(IEnumerable{CopyNumberArmSegment}, IReadOnlyList{double}?, double)"/>:
    /// reports, in input order, the segments for which
    /// <see cref="IsHomozygousDeletion(in CopyNumberArmSegment, IReadOnlyList{double}?, double, double)"/> holds
    /// (CNVkit <c>do_call</c> purity rescaling, then threshold CN 0). <paramref name="purity"/> = 1 is identical to
    /// the purity-less overload.
    /// </summary>
    /// <param name="segments">Arm-anchored copy-number segments. Must not be null.</param>
    /// <param name="thresholds">Four strictly ascending log2 cutoffs; null uses CNVkit defaults.</param>
    /// <param name="ploidy">Reference (germline) ploidy.</param>
    /// <param name="purity">Tumour purity ∈ (0, 1].</param>
    /// <returns>The homozygous-deletion segments, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment is malformed; or invalid thresholds (checked even for empty input).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ploidy"/> or <paramref name="purity"/> is invalid (checked even for empty input).</exception>
    public static IReadOnlyList<CopyNumberArmSegment> DetectHomozygousDeletions(
        IEnumerable<CopyNumberArmSegment> segments,
        IReadOnlyList<double>? thresholds,
        double ploidy,
        double purity)
    {
        ArgumentNullException.ThrowIfNull(segments);
        var cutoffs = ValidateThresholds(thresholds);
        CopyNumberMath.ValidatePloidy(ploidy);
        CopyNumberMath.ValidatePurity(purity);

        var result = new List<CopyNumberArmSegment>();
        foreach (CopyNumberArmSegment segment in segments)
        {
            if (IsHomozygousDeletion(segment, cutoffs, ploidy, purity))
            {
                result.Add(segment);
            }
        }

        return result;
    }

    /// <summary>
    /// Arm-level mapping (see the <see cref="CopyNumberRegion"/> overload for GISTIC2 locus-overlap mapping).
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
        return GenesOnAffectedArms(deletions, TumorSuppressorArms);
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

    /// <summary>
    /// GRCh38 loci of the default tumour-suppressor panel (same genes and order as the arm-level panel of
    /// <see cref="IdentifyDeletedTumorSuppressors(IEnumerable{CopyNumberArmSegment})"/>). Source: the GISTIC2 hg38
    /// reference gene table (broadinstitute/gistic2 <c>refgenes/Gencode.v22.170324/gencode_genes.tsv</c>, GENCODE v22,
    /// 1-based closed): TP53 chr17:7,661,779–7,687,550; RB1 chr13:48,303,751–48,481,986; CDKN2A chr9:21,967,753–21,995,301;
    /// PTEN chr10:87,863,113–87,971,930; BRCA1 chr17:43,044,295–43,125,483; BRCA2 chr13:32,315,474–32,400,266.
    /// </summary>
    public static IReadOnlyList<GeneLocus> DefaultTumorSuppressorLoci { get; } = new GeneLocus[]
    {
        new("TP53", "17", 7_661_779, 7_687_550),
        new("RB1", "13", 48_303_751, 48_481_986),
        new("CDKN2A", "9", 21_967_753, 21_995_301),
        new("PTEN", "10", 87_863_113, 87_971_930),
        new("BRCA1", "17", 43_044_295, 43_125_483),
        new("BRCA2", "13", 32_315_474, 32_400_266),
    };

    /// <summary>
    /// Locus-overlap tumour-suppressor mapping (GISTIC2 rule): reports, in panel order and each at most once, every
    /// panel gene whose locus overlaps a deleted region — GISTIC2 <c>genes_at.m</c> with <c>partial_hits = 1</c>: same
    /// chromosome AND <c>gene.start ≤ region.end</c> AND <c>gene.end ≥ region.start</c> (closed intervals). Unlike the
    /// arm-level <see cref="IdentifyDeletedTumorSuppressors(IEnumerable{CopyNumberArmSegment})"/>, a deletion elsewhere on
    /// 17p does NOT report TP53. Regions are used as given (no closest-gene fallback, no marker widening).
    /// </summary>
    /// <param name="deletedRegions">Deleted regions (1-based closed coordinates). Must not be null.</param>
    /// <param name="genePanel">Gene loci to test; null uses <see cref="DefaultTumorSuppressorLoci"/> (GRCh38).</param>
    /// <returns>Distinct symbols of panel genes overlapping any region, in panel order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="deletedRegions"/> is null.</exception>
    /// <exception cref="ArgumentException">A region or locus has an empty chromosome / symbol or End &lt; Start.</exception>
    public static IReadOnlyList<string> IdentifyDeletedTumorSuppressors(
        IEnumerable<CopyNumberRegion> deletedRegions,
        IReadOnlyList<GeneLocus>? genePanel = null)
    {
        ArgumentNullException.ThrowIfNull(deletedRegions);
        return GenesOverlappingRegions(deletedRegions, genePanel ?? DefaultTumorSuppressorLoci);
    }

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
    /// Reference human genome assembly (coordinate system of the input segments). Its autosomal chromosome-size
    /// table is the denominator of the explicit reference-assembly WGD fraction
    /// (<see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment}, ReferenceGenome)"/>). Note: in
    /// facets-suite the <c>genome</c> build (<c>'hg19' | 'hg18' | 'hg38'</c>) supplies centromere positions only;
    /// its WGD denominator comes from <c>get_sample_genome</c> (interrogated span) — the default
    /// <see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment})"/>.
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
    /// denominator of the reference-genome WGD fraction (facets-suite's
    /// <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c> evaluated on the reference chromosome sizes
    /// instead of the sample's interrogated spans); sizes from UCSC <c>*.chrom.sizes</c>.
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
    /// Chromosome arm geometry of a reference assembly from the UCSC cytoBand table, split as GISTIC2 splits arms
    /// (<c>normalize_by_arm_length.m</c> l.40-53: an arm is every band whose name starts with chromosome + <c>p</c> /
    /// <c>q</c>, so the p-arm acen band (<c>p11.1</c>/<c>p11</c>) belongs to p and the q-arm acen band to q). Half-open,
    /// 0-based UCSC coordinates: p-arm = [0, <see cref="PArmEnd"/>), q-arm = [<see cref="PArmEnd"/>, <see cref="Length"/>);
    /// GISTIC2 arm length = <c>band.end − band.start + 1</c> with <c>band.start = cyto.start + 1</c>, i.e. the same values.
    /// </summary>
    /// <param name="Chromosome">Chromosome label without "chr" ("1"–"22", "X", "Y").</param>
    /// <param name="CentromereStart">Start of the p-arm acen band (centromere start).</param>
    /// <param name="PArmEnd">End of the p-arm acen band = start of the q-arm acen band (GISTIC2 p/q boundary).</param>
    /// <param name="CentromereEnd">End of the q-arm acen band (centromere end).</param>
    /// <param name="Length">Chromosome length (end of the last band; equals UCSC <c>chrom.sizes</c>).</param>
    public readonly record struct ChromosomeArms(string Chromosome, long CentromereStart, long PArmEnd, long CentromereEnd, long Length)
    {
        /// <summary>GISTIC2 p-arm length in bp (= <see cref="PArmEnd"/>).</summary>
        public long PArmLength => PArmEnd;

        /// <summary>GISTIC2 q-arm length in bp (<see cref="Length"/> − <see cref="PArmEnd"/>).</summary>
        public long QArmLength => Length - PArmEnd;
    }

    /// <summary>
    /// GRCh38 / hg38 arm geometry (chr1–22, X, Y) from UCSC <c>hg38/database/cytoBand.txt</c> as shipped in GISTIC2
    /// (broadinstitute/gistic2 26c590bd <c>refgenes/hg38.UCSC.add_mir.160920/cytoBand.txt</c>, identical to
    /// <c>refgenes/Gencode.v22.170324/cytoBand.txt</c> and to the <c>cyto</c> struct of
    /// <c>support/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat</c>); acen rows, e.g. chr1 l.34-35
    /// (<c>121700000 123400000 p11.1</c> / <c>123400000 125100000 q11</c>).
    /// </summary>
    private static readonly ChromosomeArms[] GRCh38Arms =
    {
        new("1", 121_700_000L, 123_400_000L, 125_100_000L, 248_956_422L),
        new("2", 91_800_000L, 93_900_000L, 96_000_000L, 242_193_529L),
        new("3", 87_800_000L, 90_900_000L, 94_000_000L, 198_295_559L),
        new("4", 48_200_000L, 50_000_000L, 51_800_000L, 190_214_555L),
        new("5", 46_100_000L, 48_800_000L, 51_400_000L, 181_538_259L),
        new("6", 58_500_000L, 59_800_000L, 62_600_000L, 170_805_979L),
        new("7", 58_100_000L, 60_100_000L, 62_100_000L, 159_345_973L),
        new("8", 43_200_000L, 45_200_000L, 47_200_000L, 145_138_636L),
        new("9", 42_200_000L, 43_000_000L, 45_500_000L, 138_394_717L),
        new("10", 38_000_000L, 39_800_000L, 41_600_000L, 133_797_422L),
        new("11", 51_000_000L, 53_400_000L, 55_800_000L, 135_086_622L),
        new("12", 33_200_000L, 35_500_000L, 37_800_000L, 133_275_309L),
        new("13", 16_500_000L, 17_700_000L, 18_900_000L, 114_364_328L),
        new("14", 16_100_000L, 17_200_000L, 18_200_000L, 107_043_718L),
        new("15", 17_500_000L, 19_000_000L, 20_500_000L, 101_991_189L),
        new("16", 35_300_000L, 36_800_000L, 38_400_000L, 90_338_345L),
        new("17", 22_700_000L, 25_100_000L, 27_400_000L, 83_257_441L),
        new("18", 15_400_000L, 18_500_000L, 21_500_000L, 80_373_285L),
        new("19", 24_200_000L, 26_200_000L, 28_100_000L, 58_617_616L),
        new("20", 25_700_000L, 28_100_000L, 30_400_000L, 64_444_167L),
        new("21", 10_900_000L, 12_000_000L, 13_000_000L, 46_709_983L),
        new("22", 13_700_000L, 15_000_000L, 17_400_000L, 50_818_468L),
        new("X", 58_100_000L, 61_000_000L, 63_800_000L, 156_040_895L),
        new("Y", 10_300_000L, 10_400_000L, 10_600_000L, 57_227_415L),
    };

    /// <summary>
    /// GRCh37 / hg19 arm geometry (chr1–22, X, Y) from the UCSC hg19 cytoBand table (862 bands) bundled in GISTIC2 as
    /// the <c>cyto</c> struct of <c>support/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat</c> (gistic2 26c590bd;
    /// read with GNU Octave 8.4 <c>load</c>), e.g. chr1 acen <c>121500000 125000000 1p11.1</c> /
    /// <c>125000000 128900000 1q11</c>.
    /// </summary>
    private static readonly ChromosomeArms[] GRCh37Arms =
    {
        new("1", 121_500_000L, 125_000_000L, 128_900_000L, 249_250_621L),
        new("2", 90_500_000L, 93_300_000L, 96_800_000L, 243_199_373L),
        new("3", 87_900_000L, 91_000_000L, 93_900_000L, 198_022_430L),
        new("4", 48_200_000L, 50_400_000L, 52_700_000L, 191_154_276L),
        new("5", 46_100_000L, 48_400_000L, 50_700_000L, 180_915_260L),
        new("6", 58_700_000L, 61_000_000L, 63_300_000L, 171_115_067L),
        new("7", 58_000_000L, 59_900_000L, 61_700_000L, 159_138_663L),
        new("8", 43_100_000L, 45_600_000L, 48_100_000L, 146_364_022L),
        new("9", 47_300_000L, 49_000_000L, 50_700_000L, 141_213_431L),
        new("10", 38_000_000L, 40_200_000L, 42_300_000L, 135_534_747L),
        new("11", 51_600_000L, 53_700_000L, 55_700_000L, 135_006_516L),
        new("12", 33_300_000L, 35_800_000L, 38_200_000L, 133_851_895L),
        new("13", 16_300_000L, 17_900_000L, 19_500_000L, 115_169_878L),
        new("14", 16_100_000L, 17_600_000L, 19_100_000L, 107_349_540L),
        new("15", 15_800_000L, 19_000_000L, 20_700_000L, 102_531_392L),
        new("16", 34_600_000L, 36_600_000L, 38_600_000L, 90_354_753L),
        new("17", 22_200_000L, 24_000_000L, 25_800_000L, 81_195_210L),
        new("18", 15_400_000L, 17_200_000L, 19_000_000L, 78_077_248L),
        new("19", 24_400_000L, 26_500_000L, 28_600_000L, 59_128_983L),
        new("20", 25_600_000L, 27_500_000L, 29_400_000L, 63_025_520L),
        new("21", 10_900_000L, 13_200_000L, 14_300_000L, 48_129_895L),
        new("22", 12_200_000L, 14_700_000L, 17_900_000L, 51_304_566L),
        new("X", 58_100_000L, 60_600_000L, 63_000_000L, 155_270_560L),
        new("Y", 11_600_000L, 12_500_000L, 13_400_000L, 59_373_566L),
    };

    /// <summary>
    /// Chromosome arm lengths of a reference assembly (24 entries: chr1–22, X, Y in that order), split at the p/q
    /// acen band boundary as GISTIC2 does (<see cref="ChromosomeArms"/>). GRCh38: Σ p 1,030,800,000 bp, Σ q
    /// 2,057,469,832 bp; GRCh37: Σ p 1,040,600,000 bp, Σ q 2,055,077,412 bp. Chromosome lengths 1–22 equal
    /// <see cref="GetAutosomeLengths"/>.
    /// </summary>
    /// <param name="genome">The reference assembly.</param>
    /// <returns>The arm table of the assembly.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static IReadOnlyList<ChromosomeArms> GetChromosomeArmLengths(ReferenceGenome genome) => genome switch
    {
        ReferenceGenome.GRCh38 => GRCh38Arms,
        ReferenceGenome.GRCh37 => GRCh37Arms,
        _ => throw new ArgumentOutOfRangeException(nameof(genome), genome, "Unknown reference genome."),
    };

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

        return int.TryParse(StripChrPrefix(chromosome.AsSpan()), out number) && number is >= 1 and <= AutosomeCount;
    }

    /// <summary>
    /// Removes a leading case-insensitive "chr" from a chromosome label ("chr7" → "7", "ChrX" → "X"); a label that is
    /// exactly "chr" (or shorter) is returned unchanged. Shared by the WGD autosome filter
    /// (<see cref="TryGetAutosomeNumber"/>) and the sex-chromosome filter (<see cref="IsSexChromosomeLabel"/>).
    /// </summary>
    private static ReadOnlySpan<char> StripChrPrefix(ReadOnlySpan<char> name) =>
        name.Length > 3 && name.StartsWith("chr", StringComparison.OrdinalIgnoreCase) ? name[3..] : name;

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
        double totalLength = 0.0; // double (as ASCAT/facets-suite sum numerics): Σ L must not overflow Int64
        int segmentCount = 0;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            long length = segment.Length;
            int totalCopyNumber = segment.MajorCopyNumber + segment.MinorCopyNumber;
            weightedCopyNumberSum += (double)totalCopyNumber * length;
            totalLength += length;
            segmentCount++;
        }

        if (segmentCount == 0)
        {
            throw new ArgumentException(
                "Cannot estimate ploidy from an empty segment set (the length-weighted mean is undefined).",
                nameof(segments));
        }

        // ψ = Σ(CN_i · L_i) / Σ(L_i) — Patchwork length-weighted mean of total copy number.
        return weightedCopyNumberSum / totalLength;
    }

    /// <summary>
    /// Determines whether a tumour genome has undergone whole-genome doubling (WGD) exactly as the facets-suite
    /// reference implementation of the Bielski et al. rule does (Bielski et al. 2018, Nat Genet 50:1189–1195,
    /// PMID 30013179): more than half of the <i>interrogated</i> autosomal genome has major-allele copy number ≥ 2.
    /// Source: facets-suite <c>R/copy-number-scores.R</c> (master) — <c>calculate_fraction_cna</c> calls
    /// <c>is_genome_doubled(segs, get_sample_genome(segs, genome), treshold = 0.5)</c>; the denominator is
    /// <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c> with <c>size = max(end) − min(start)</c> of each
    /// chromosome's own segments; <c>frac_elevated_mcn = sum(length[mcn ≥ 2 &amp; chrom %in% 1:22]) / autosomal_genome</c>;
    /// <c>wgd = frac_elevated_mcn &gt; 0.5</c>; <c>mcn = tcn − lcn</c> = max(Major, Minor).
    /// <para>This is the default (review 2026-09, F32: the reference's denominator); it is the same computation as
    /// <see cref="DetectWholeGenomeDoublingFromSuppliedLength"/>. The reference-assembly denominator (UCSC
    /// <c>*.chrom.sizes</c>, robust to partial-genome inputs but not the facets-suite call) is the explicit overload
    /// <see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment}, ReferenceGenome)"/>.</para>
    /// The test uses the major (not total) copy number, so a balanced diploid genome (all 1:1, total CN 2,
    /// major CN 1) is NOT doubled, whereas a 2:0 LOH or 2:2 genome IS.
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments. Must not be null, must contain at least one autosomal (chr1–22,
    /// bare or "chr"-prefixed) segment, and every segment must have End &gt; Start and non-negative copy numbers.
    /// </param>
    /// <returns><c>true</c> when more than half of the interrogated autosomal genome has major copy number ≥ 2.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="segments"/> contains no autosomal segment (facets-suite 0/0 → <c>NA</c>), or a segment has
    /// End ≤ Start or a negative copy number.
    /// </exception>
    public static bool DetectWholeGenomeDoubling(IEnumerable<AlleleSpecificSegment> segments)
        => DetectWholeGenomeDoublingFromSuppliedLength(segments);

    /// <summary>
    /// Whole-genome-doubling call by the Bielski et al. major-CN ≥ 2 / &gt; 50 % rule, with the <b>reference
    /// assembly's</b> autosomal length (Σ chr1–22 of the UCSC <c>*.chrom.sizes</c> table) as the denominator
    /// instead of the sample's interrogated span. Numerator, threshold and major-CN rule follow facets-suite
    /// <c>is_genome_doubled</c>: <c>frac_elevated_mcn = sum(length where mcn ≥ 2 &amp; chrom %in% 1:22) / G</c>;
    /// <c>wgd = frac_elevated_mcn &gt; 0.5</c>, <c>mcn = tcn − lcn</c> = max(Major, Minor).
    /// <para><b>Not the facets-suite call:</b> facets-suite divides by <c>get_sample_genome</c> spans (the default
    /// <see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment})"/>). This explicit option judges
    /// partial-genome inputs (a few chromosomes, or segments that do not reach the telomeres) against the whole
    /// autosomal genome; for a segmentation spanning every autosome end-to-end the two agree.</para>
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments (<see cref="AlleleSpecificSegment"/>). Only segments on autosomes
    /// (chromosomes 1–22, "chr"-prefixed or bare) contribute to the elevated-major-CN numerator. Must not be
    /// null, and every segment must have End &gt; Start and non-negative copy numbers.
    /// </param>
    /// <param name="genome">Reference assembly whose autosomal chromosome-size table is the fraction denominator.</param>
    /// <returns><c>true</c> when more than half the reference autosomal genome has major copy number ≥ 2.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment has End ≤ Start or a negative copy number.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is not a defined value.</exception>
    public static bool DetectWholeGenomeDoubling(
        IEnumerable<AlleleSpecificSegment> segments,
        ReferenceGenome genome)
    {
        ArgumentNullException.ThrowIfNull(segments);

        long autosomalGenomeLength = GetAutosomalGenomeLength(genome);

        // as.numeric(...) sums in facets-suite: accumulate in double so Σ length cannot overflow Int64.
        double elevatedLength = 0.0;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            // facets-suite: numerator restricted to autosomes (chrom %in% 1:22). Non-autosomal segments are
            // ignored (sex chromosomes / contigs do not contribute to the autosomal WGD fraction).
            if (!TryGetAutosomeNumber(segment.Chromosome, out _))
            {
                continue;
            }

            // mcn = tcn − lcn = max(Major, Minor); elevated when mcn ≥ 2 (facets-suite segs$mcn >= 2).
            if (IsElevatedMajorCopyNumber(segment))
            {
                elevatedLength += segment.Length;
            }
        }

        // wgd = frac_elevated_mcn > 0.5 (strict), denominator = reference autosomal genome length.
        double fractionElevatedMajorCn = elevatedLength / autosomalGenomeLength;
        return fractionElevatedMajorCn > WholeGenomeDoublingFractionThreshold;
    }

    /// <summary>
    /// Determines whole-genome doubling exactly as the facets-suite reference implementation does (the same call as
    /// the default <see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment})"/>), with the
    /// genome-fraction denominator taken from the <b>supplied (interrogated) segments</b> rather than a reference
    /// chromosome-size table. Source: facets-suite <c>R/copy-number-scores.R</c> (master; PMID 30013179, Bielski
    /// et al. 2018): <c>calculate_fraction_cna</c> passes <c>sample_chrom_info = get_sample_genome(segs, genome)</c>
    /// — per chromosome <c>size = max(end) − min(start)</c> of that chromosome's segments — to
    /// <c>is_genome_doubled</c>, which computes
    /// <c>autosomal_genome = sum(chrom_info$size[chr %in% 1:22])</c>,
    /// <c>frac_elevated_mcn = sum(length[mcn ≥ 2 &amp; chrom %in% 1:22]) / autosomal_genome</c> and
    /// <c>wgd = frac_elevated_mcn &gt; 0.5</c>, with <c>length = end − start</c> and <c>mcn = tcn − lcn</c>
    /// (<c>parse_segs</c>). Hence: (i) only autosomal (chr1–22) segments enter the numerator and the denominator;
    /// (ii) the denominator is the sum of each autosome's interrogated span, so unsegmented gaps between the first
    /// and last segment of a chromosome count toward the genome but not toward the elevated length.
    /// </summary>
    /// <param name="segments">
    /// Allele-specific copy-number segments. Must not be null, must contain at least one autosomal (chr1–22,
    /// bare or "chr"-prefixed) segment, and every segment must have End &gt; Start and non-negative copy numbers.
    /// The major copy number is <c>tcn − lcn</c> with <c>lcn</c> the lesser allele copy number, i.e.
    /// max(Major, Minor), so the call does not depend on the allele labelling.
    /// </param>
    /// <returns><c>true</c> when more than half of the interrogated autosomal genome has major copy number ≥ 2.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="segments"/> contains no autosomal segment (facets-suite divides 0 by 0 and returns
    /// <c>NA</c>: the fraction is undefined), or a segment has End ≤ Start or a negative copy number.
    /// </exception>
    public static bool DetectWholeGenomeDoublingFromSuppliedLength(IEnumerable<AlleleSpecificSegment> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);

        // get_sample_genome: per-autosome interrogated span [min(start), max(end)).
        var spanStart = new long[AutosomeCount];
        var spanEnd = new long[AutosomeCount];
        var seen = new bool[AutosomeCount];

        // as.numeric(...) sums in facets-suite: accumulate in double so Σ length cannot overflow Int64.
        double elevatedLength = 0.0;
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            if (!TryGetAutosomeNumber(segment.Chromosome, out int autosome))
            {
                continue; // chrom %in% 1:22 — sex chromosomes / contigs excluded from numerator and denominator
            }

            int k = autosome - 1;
            if (!seen[k])
            {
                seen[k] = true;
                spanStart[k] = segment.Start;
                spanEnd[k] = segment.End;
            }
            else
            {
                spanStart[k] = Math.Min(spanStart[k], segment.Start);
                spanEnd[k] = Math.Max(spanEnd[k], segment.End);
            }

            if (IsElevatedMajorCopyNumber(segment))
            {
                elevatedLength += segment.Length;
            }
        }

        double interrogatedAutosomalGenome = 0.0;
        for (int k = 0; k < AutosomeCount; k++)
        {
            if (seen[k])
            {
                interrogatedAutosomalGenome += (double)spanEnd[k] - spanStart[k];
            }
        }

        if (interrogatedAutosomalGenome <= 0.0)
        {
            throw new ArgumentException(
                "Cannot assess whole-genome doubling without autosomal (chr1–22) segments " +
                "(the interrogated autosomal genome is empty, so the fraction is undefined).",
                nameof(segments));
        }

        // wgd = frac_elevated_mcn > 0.5 (strict) — facets-suite is_genome_doubled with get_sample_genome sizes.
        return elevatedLength / interrogatedAutosomalGenome > WholeGenomeDoublingFractionThreshold;
    }

    /// <summary>
    /// facets-suite elevated-major-CN test: <c>mcn = tcn − lcn ≥ 2</c> (<c>parse_segs</c>, <c>is_genome_doubled</c>),
    /// where <c>lcn</c> is the lesser (minor) allele copy number, so <c>mcn = max(Major, Minor)</c> irrespective of
    /// how the two allele copy numbers are labelled in the input record.
    /// </summary>
    private static bool IsElevatedMajorCopyNumber(in AlleleSpecificSegment segment)
        => Math.Max(segment.MajorCopyNumber, segment.MinorCopyNumber) >= WholeGenomeDoublingMajorCopyNumber;

    /// <summary>
    /// Estimates the tumour ploidy as ASCAT reports it: the <b>probe-count</b>-weighted mean of per-segment total copy
    /// number, ψ = Σ(CN_i · n_i) / Σ n_i, with CN_i = Major + Minor and n_i the number of probes (SNPs / loci) in
    /// segment i. Source: ASCAT <c>R/ascat.runAscat.R</c> (VanLoo-lab/ascat, master): the search-time ploidy
    /// <c>ploidy = sum((nA+nB) * s[, "length"]) / sum(s[, "length"])</c> (runASCAT, l. 283), where
    /// <c>s[, "length"]</c> is the probe count of each <c>make_segments</c> segment, and the reported
    /// <c>ploidy = mean(nA+nB, na.rm=TRUE)</c> over per-probe copy numbers (ascat.runAscat, l. 98) — a per-probe mean
    /// is exactly the probe-count-weighted mean of the segment values. This differs from
    /// <see cref="EstimatePloidy(IEnumerable{AlleleSpecificSegment})"/> (Patchwork, bp-weighted) unless probe density
    /// is uniform. The caller chooses which probes the counts represent (ASCAT's reported ploidy counts every non-NA
    /// logR probe on every chromosome; the search-time value counts autosomal heterozygous probes).
    /// </summary>
    /// <param name="segments">Allele-specific integer copy-number segments; every segment must have End &gt; Start and
    /// non-negative copy numbers.</param>
    /// <param name="probeCounts">Number of probes in each segment, aligned by position with
    /// <paramref name="segments"/>; every count must be ≥ 1 (an ASCAT segment contains at least one probe).</param>
    /// <returns>The probe-count-weighted mean total copy number.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The inputs are empty or of different lengths, a segment is invalid, or a
    /// probe count is &lt; 1.</exception>
    public static double EstimatePloidy(IEnumerable<AlleleSpecificSegment> segments, IEnumerable<int> probeCounts)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ArgumentNullException.ThrowIfNull(probeCounts);

        double weightedCopyNumberSum = 0.0;
        double totalProbes = 0.0;
        int segmentCount = 0;
        using IEnumerator<int> counts = probeCounts.GetEnumerator();
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            if (!counts.MoveNext())
            {
                throw new ArgumentException("Each segment needs exactly one probe count (fewer counts than segments).",
                    nameof(probeCounts));
            }

            int probes = counts.Current;
            if (probes < 1)
            {
                throw new ArgumentException(
                    $"Every segment must contain at least one probe (got {probes}).", nameof(probeCounts));
            }

            // ASCAT: sum((nA+nB) * length) / sum(length), length = probe count.
            weightedCopyNumberSum += ((double)segment.MajorCopyNumber + segment.MinorCopyNumber) * probes;
            totalProbes += probes;
            segmentCount++;
        }

        if (counts.MoveNext())
        {
            throw new ArgumentException("Each segment needs exactly one probe count (more counts than segments).",
                nameof(probeCounts));
        }

        if (segmentCount == 0)
        {
            throw new ArgumentException(
                "Cannot estimate ploidy from an empty segment set (the probe-weighted mean is undefined).",
                nameof(segments));
        }

        return weightedCopyNumberSum / totalProbes;
    }

    /// <summary>
    /// Whole-genome-doubling status as defined by ASCAT <c>ascat.metrics</c> (column <c>WGD</c>), derived from the
    /// size-weighted mode of the autosomal major-allele copy number.
    /// </summary>
    public enum AscatWgdStatus
    {
        /// <summary>ASCAT <c>NA</c>: the mode of the major allele is 0 (or outside 1–5), so WGD is undefined.</summary>
        NotAvailable,

        /// <summary>ASCAT <c>0</c>: mode of the major allele = 1 (no WGD).</summary>
        NoWgd,

        /// <summary>ASCAT <c>1</c>: mode of the major allele = 2 (one WGD).</summary>
        Wgd,

        /// <summary>ASCAT <c>"1+"</c>: mode of the major allele ∈ {3, 4, 5} (at least one WGD).</summary>
        WgdPlus,
    }

    /// <summary>
    /// The WGD-related genome metrics of ASCAT <c>ascat.metrics</c>, computed on autosomes only.
    /// </summary>
    /// <param name="ModeMinorAllele">ASCAT <c>mode_minA</c>: size-weighted mode of the minor allele copy number (capped at 5).</param>
    /// <param name="ModeMajorAllele">ASCAT <c>mode_majA</c>: size-weighted mode of the major allele copy number (capped at 5).</param>
    /// <param name="WgdStatus">ASCAT <c>WGD</c>: NA / 0 / 1 / "1+" from <paramref name="ModeMajorAllele"/>.</param>
    /// <param name="GenomicInstability">ASCAT <c>GI</c>: fraction of the autosomal genome not at the baseline
    /// state (1:1 without WGD, 2:2 with WGD or "1+"), rounded to 4 decimals; <c>null</c> when WGD is NA.</param>
    /// <param name="LossOfHeterozygosity">ASCAT <c>LOH</c>: fraction of the autosomal genome with minor allele 0,
    /// rounded to 4 decimals.</param>
    public readonly record struct AscatGenomeMetrics(
        int ModeMinorAllele,
        int ModeMajorAllele,
        AscatWgdStatus WgdStatus,
        double? GenomicInstability,
        double LossOfHeterozygosity)
    {
        /// <summary>The ASCAT <c>WGD</c> label: "NA", "0", "1" or "1+".</summary>
        public string WgdLabel => WgdStatus switch
        {
            AscatWgdStatus.NoWgd => "0",
            AscatWgdStatus.Wgd => "1",
            AscatWgdStatus.WgdPlus => "1+",
            _ => "NA",
        };
    }

    /// <summary>Cap applied to allele copy numbers before taking the mode (ascat.metrics <c>modeAllele</c>: <c>y[y&gt;5]=5</c>).</summary>
    private const int AscatModeAlleleCap = 5;

    /// <summary>Decimal places of the rounded ASCAT GI / LOH metrics (<c>round(…, 4)</c>).</summary>
    private const int AscatMetricDecimals = 4;

    /// <summary>
    /// Computes ASCAT's whole-genome-doubling status and genomic-instability (GI) score — the
    /// <c>mode_minA</c>, <c>mode_majA</c>, <c>WGD</c>, <c>GI</c> and <c>LOH</c> columns of <c>ascat.metrics</c>.
    /// Source: ASCAT <c>R/ascat.metrics.R</c> (VanLoo-lab/ascat, master), ported verbatim:
    /// <list type="bullet">
    /// <item>Only autosomes: <c>profile[chr %in% setdiff(chrs, sexchromosomes)]</c> with ASCAT's default
    /// <c>sexchromosomes = c("X","Y")</c> (an optional "chr" prefix is ignored).</item>
    /// <item><c>modeAllele</c>: per-segment weight <c>(endpos − startpos)/1e6</c> (no +1); allele copy number
    /// <c>round</c>ed and capped at 5; weights summed per value (<c>tapply</c>, groups in ascending value order);
    /// the mode is the first maximum after the stable decreasing <c>order</c>, so an exact tie resolves to the
    /// <b>smaller</b> copy number.</item>
    /// <item>WGD: <c>mode_majA</c> = 0 → NA; 1 → "0"; 2 → "1"; 3–5 → "1+".</item>
    /// <item><c>computeGIscore</c>: <c>round(1 − Σ size[nMajor = b ∧ nMinor = b] / Σ size, 4)</c>, size =
    /// <c>endpos − startpos + 1</c>, baseline b = 1 for WGD "0" and b = 2 for WGD "1" and "1+"; NA when WGD is NA.</item>
    /// <item>LOH: <c>round(Σ size[nMinor = 0] / Σ size, 4)</c> over the same autosomal profile.</item>
    /// </list>
    /// Major/minor are taken as max/min of the two allele copy numbers (ASCAT's nMajor ≥ nMinor), so the result does
    /// not depend on the allele labelling. Per-value weights are summed in input order in double precision (R's
    /// <c>sum</c> accumulates in extended precision; this can matter only for ties broken by the last ulp).
    /// </summary>
    /// <param name="segments">Allele-specific integer copy-number segments (ASCAT <c>segments</c>). Must contain at
    /// least one autosomal segment; every segment must have End &gt; Start and non-negative copy numbers.</param>
    /// <returns>The ASCAT WGD / GI genome metrics.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">A segment is invalid, or there is no autosomal segment (ASCAT's mode is
    /// then empty and the metrics undefined).</exception>
    public static AscatGenomeMetrics ComputeAscatGenomeMetrics(IEnumerable<AlleleSpecificSegment> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);

        // tapply groups: index = capped allele copy number 0..5.
        var majorWeight = new double[AscatModeAlleleCap + 1];
        var minorWeight = new double[AscatModeAlleleCap + 1];
        var majorSeen = new bool[AscatModeAlleleCap + 1];
        var minorSeen = new bool[AscatModeAlleleCap + 1];
        var profile = new List<(int Major, int Minor, double Size)>();
        foreach (AlleleSpecificSegment segment in segments)
        {
            ValidateSegment(segment);
            if (segment.Chromosome is null)
            {
                throw new ArgumentException("A segment has a null chromosome label.", nameof(segments));
            }

            if (IsSexChromosomeLabel(segment.Chromosome))
            {
                continue; // setdiff(chrs, sexchromosomes)
            }

            int major = Math.Max(segment.MajorCopyNumber, segment.MinorCopyNumber);
            int minor = Math.Min(segment.MajorCopyNumber, segment.MinorCopyNumber);
            double modeWeight = ((double)segment.End - segment.Start) / 1e6; // (endpos − startpos)/1e6
            int majorKey = Math.Min(major, AscatModeAlleleCap);
            int minorKey = Math.Min(minor, AscatModeAlleleCap);
            majorWeight[majorKey] += modeWeight;
            majorSeen[majorKey] = true;
            minorWeight[minorKey] += modeWeight;
            minorSeen[minorKey] = true;
            profile.Add((major, minor, (double)segment.End - segment.Start + 1.0)); // size = endpos − startpos + 1
        }

        if (profile.Count == 0)
        {
            throw new ArgumentException(
                "ASCAT genome metrics are computed on autosomes only (X/Y excluded); no autosomal segment was supplied.",
                nameof(segments));
        }

        int modeMajor = AscatModeAllele(majorWeight, majorSeen);
        int modeMinor = AscatModeAllele(minorWeight, minorSeen);

        double totalSize = 0.0, lohSize = 0.0;
        foreach ((int _, int minor, double size) in profile)
        {
            totalSize += size;
            if (minor == 0)
            {
                lohSize += size;
            }
        }

        double loh = Math.Round(lohSize / totalSize, AscatMetricDecimals, MidpointRounding.ToEven);

        AscatWgdStatus status = modeMajor switch
        {
            1 => AscatWgdStatus.NoWgd,
            2 => AscatWgdStatus.Wgd,
            >= 3 and <= AscatModeAlleleCap => AscatWgdStatus.WgdPlus,
            _ => AscatWgdStatus.NotAvailable, // mode_majA == 0
        };

        double? gi = null;
        if (status != AscatWgdStatus.NotAvailable)
        {
            // computeGIscore: baseline 1 (WGD 0) or 2 (WGD 1; "1+" is scored with WGD = 1).
            int baseline = status == AscatWgdStatus.NoWgd ? 1 : 2;
            double baselineSize = 0.0;
            foreach ((int major, int minor, double size) in profile)
            {
                if (major == baseline && minor == baseline)
                {
                    baselineSize += size;
                }
            }

            gi = Math.Round(1.0 - baselineSize / totalSize, AscatMetricDecimals, MidpointRounding.ToEven);
        }

        return new AscatGenomeMetrics(modeMinor, modeMajor, status, gi, loh);
    }

    /// <summary>
    /// ascat.metrics <c>modeAllele</c> selection: the copy-number value with the largest summed weight; groups are
    /// visited in ascending value order (tapply) and only a strictly larger weight replaces the current best
    /// (stable <c>order(decreasing = TRUE)</c> + <c>which.max</c>), so ties go to the smaller value.
    /// </summary>
    private static int AscatModeAllele(double[] weights, bool[] seen)
    {
        int best = -1;
        for (int value = 0; value < weights.Length; value++)
        {
            if (seen[value] && (best < 0 || weights[value] > weights[best]))
            {
                best = value;
            }
        }

        return best;
    }

    #endregion

    /// <summary>True when <paramref name="chromosome"/> is X or Y (an optional "chr" prefix is ignored, case-insensitive): the sex chromosomes are excluded from the ploidy estimate.</summary>
    private static bool IsSexChromosomeLabel(string chromosome)
    {
        ReadOnlySpan<char> name = StripChrPrefix(chromosome.AsSpan().Trim());
        return name.Equals("X", StringComparison.OrdinalIgnoreCase) || name.Equals("Y", StringComparison.OrdinalIgnoreCase);
    }
}
