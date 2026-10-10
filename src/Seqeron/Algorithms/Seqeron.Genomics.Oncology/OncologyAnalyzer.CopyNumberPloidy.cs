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
        ValidatePloidy(ploidy);
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
        ValidatePloidy(ploidy);
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
        ValidatePurity(purity);
        if (purity >= 1.0 || double.IsNaN(log2Ratio))
        {
            return log2Ratio;
        }

        double absolute = CopyNumberMath.Log2RatioToAbsolute(log2Ratio, ploidy, ploidy, purity);
        return CopyNumberMath.AbsoluteToLog2Ratio(absolute, ploidy);
    }

    /// <summary>Validates tumour purity ∈ (0, 1] (CNVkit <c>commands.py</c> <c>purity_value</c>).</summary>
    private static void ValidatePurity(double purity)
    {
        if (!(purity > 0.0 && purity <= 1.0))
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        }
    }

    /// <summary>
    /// Validates the reference ploidy: it must be a finite positive number, because
    /// <c>n = ploidy · 2^log2</c> (CNVkit <c>_log2_ratio_to_absolute_pure</c>) is meaningless otherwise.
    /// </summary>
    private static void ValidatePloidy(double ploidy)
    {
        if (!double.IsFinite(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy must be a finite positive number.");
        }
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
        ValidatePloidy(ploidy);
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
        ValidatePloidy(ploidy);

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
        ValidatePloidy(ploidy);
        ValidatePurity(purity);

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
    /// <para><b>Not implemented — ziggurat deconstruction.</b> GISTIC2 applies this filter to SCNA <i>events</i>
    /// produced by its ziggurat deconstruction (amplitude measured relative to the underlying level; broad levels
    /// estimated from the whole cohort). Here each input segment is treated as one event with amplitude = its
    /// <see cref="CopyNumberArmSegment.Log2Ratio"/>. Consequence: when raw segments are supplied, an arm-level gain
    /// interrupted by a focal peak (e.g. 0.5 | 1.5 | 0.5) yields flank segments that are individually &lt; 98% of the arm
    /// and are reported as focal, whereas GISTIC2 would call one broad event (0.5) plus one focal event (+1.0).
    /// Supply deconstructed events (or whole-arm-merged segments) to match GISTIC2.</para>
    /// Segment coordinates are half-open (<c>Length = End − Start</c>). <para><b>Arm-fraction units.</b> GISTIC2 measures
    /// the arm fraction in markers by default (<c>make_sample_B.m</c> → <c>normalize_by_arm_length(D,B,cyto,1,2)</c>,
    /// <c>norm_type = 1</c>: segment markers ÷ markers in the arm's cytoband span). When a segment carries
    /// <see cref="CopyNumberArmSegment.MarkerCount"/> and <see cref="CopyNumberArmSegment.ArmMarkerCount"/>, its
    /// <see cref="CopyNumberArmSegment.ArmFraction"/> is in marker units (e.g. 39/40 = 0.975 &lt; 0.98 is focal even if the
    /// segment covers 99% of the arm in bp); without marker counts it is in bp (Length ÷ ArmLength). Centromere-spanning
    /// events (GISTIC2 <c>ref_length = 2</c>: p-fraction + q-fraction) are outside the arm-anchored model.</para>
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
        ValidatePloidy(ploidy);
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
        ValidatePloidy(ploidy);

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
        ValidatePloidy(ploidy);
        ValidatePurity(purity);

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
    /// (<see cref="TryGetAutosomeNumber"/>) and the ASCAT sex-chromosome filter (<see cref="IsAscatSexChromosome"/>).
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

            if (IsAscatSexChromosome(segment.Chromosome))
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


    #region Upstream allele-specific derivation: segmentation, purity/ploidy fit, multiplicity (ONCO-ASCAT-001)

    /// <summary>
    /// Platform/technology parameter γ in the ASCAT logR model. For massively parallel sequencing data
    /// (WGS/WES/TS) γ = 1; the SNP-array default 0.55 does not apply. Source: ASCAT README / Van Loo lab
    /// (VanLoo-lab/ascat): "For massively parallel sequencing data, gamma should always be set to 1."
    /// </summary>
    public const double AscatSequencingGamma = 1.0;

    /// <summary>
    /// Worst-case squared distance of a value to the nearest integer, (1/2)² = 0.25; the per-segment term in
    /// the ASCAT theoretical-maximum-distance used to normalise goodness of fit to a percentage. Source:
    /// ascat.runAscat.R — <c>TheoretMaxdist = sum(rep(0.25, n) * length * ...)</c>.
    /// </summary>
    private const double AscatWorstCaseIntegerDistance = 0.25;

    /// <summary>
    /// Down-weight applied to balanced (BAF = 0.5) segments in the ASCAT goodness-of-fit, because such
    /// segments carry little allele-specific information. Source: ascat.runAscat.R —
    /// <c>ifelse(b == 0.5, 0.05, 1)</c>.
    /// </summary>
    private const double AscatBalancedSegmentWeight = 0.05;

    /// <summary>BAF value of a perfectly balanced (1:1) heterozygous segment; the down-weight pivot in the GoF.</summary>
    private const double BalancedBaf = 0.5;

    /// <summary>
    /// A single per-locus allele-specific measurement at a germline-heterozygous SNP: the log-R ratio (total
    /// signal, "r") and the B-allele frequency (allelic contrast, "b"). These are the two ASCAT input tracks
    /// (Van Loo et al. 2010, PNAS) and are <b>observed measurements</b> supplied by the caller — they are the
    /// raw data, not a derived quantity.
    /// </summary>
    /// <param name="Chromosome">Contig label (used to group loci into per-chromosome segments).</param>
    /// <param name="Position">0-based genomic coordinate of the SNP.</param>
    /// <param name="LogR">Log-R ratio r (log2 total-signal ratio vs the reference baseline).</param>
    /// <param name="BAF">B-allele frequency b ∈ [0, 1] at the germline-heterozygous SNP.</param>
    public readonly record struct AlleleSpecificLocus(string Chromosome, long Position, double LogR, double BAF);

    /// <summary>
    /// A genomic segment summarised by its mean logR and mean BAF, produced by allele-specific segmentation
    /// of per-locus <see cref="AlleleSpecificLocus"/> data. A single fitted logR value and a BAF value are
    /// obtained per segment (Van Loo et al. 2010, ASCAT).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Start">0-based start coordinate (first locus position in the segment).</param>
    /// <param name="End">End coordinate (last locus position in the segment; End ≥ Start).</param>
    /// <param name="MeanLogR">Length-unweighted mean logR over the segment's loci.</param>
    /// <param name="MeanBAF">Mean "folded" BAF (distance from 0.5, re-centred) over the segment's loci.</param>
    /// <param name="LocusCount">Number of loci summarised by the segment.</param>
    public readonly record struct AlleleSpecificSegmentSummary(
        string Chromosome,
        long Start,
        long End,
        double MeanLogR,
        double MeanBAF,
        int LocusCount)
    {
        /// <summary>Segment length in base pairs (End − Start).</summary>
        public long Length => End - Start;
    }

    /// <summary>
    /// Result of the joint ASCAT purity/ploidy fit: the recovered purity ρ, the ASCAT output ploidy, the goodness of
    /// fit, the allele-specific integer copy-number segments, and the model ploidy parameter ψ.
    /// </summary>
    /// <param name="Purity">Recovered tumour purity ρ (aberrant cell fraction) ∈ (0, 1] (ASCAT <c>purity</c>).</param>
    /// <param name="Ploidy">ASCAT output <c>ploidy</c>: the mean integer total copy number (major + minor) over the
    /// heterozygous probes, i.e. the <see cref="AlleleSpecificSegmentSummary.LocusCount"/>-weighted mean of the emitted
    /// segments' total copy number. It equals ψ for an exact integer-copy-number genome. From the germline-aware path
    /// (<see cref="FitPurityPloidyFromAspcf"/>, B24 F37) it is ASCAT's <c>mean(nA + nB)</c> over <b>all</b> probes,
    /// heterozygous and homozygous.</param>
    /// <param name="GoodnessOfFit">Percentage goodness of fit (1 − distance/TheoretMaxdist)·100, in (−∞, 100].</param>
    /// <param name="Segments">The allele-specific integer copy-number segments (major/minor CN) implied by (ρ, ψ).</param>
    public readonly record struct PurityPloidyFit(
        double Purity,
        double Ploidy,
        double GoodnessOfFit,
        IReadOnlyList<AlleleSpecificSegment> Segments)
    {
        /// <summary>The selected ploidy parameter ψ of the ASCAT model (ASCAT <c>psi</c>, a grid value).</summary>
        public double Psi { get; init; }

        /// <summary>ASCAT <c>nonaberrant</c> flag: ≤ 3 % of the probes are allelically imbalanced and no imbalanced
        /// segment exceeds 0.5 % of the probes (MINABB / MINABBREGION).</summary>
        public bool IsNonAberrant { get; init; }
    }

    /// <summary>ASCAT sample sex (<c>gender</c> argument of <c>ascat.loadData</c> / <c>runASCAT</c>).</summary>
    public enum AscatGender
    {
        /// <summary><c>"XX"</c> (ASCAT default when <c>gender = NULL</c>): X and Y are emitted with the diploid model.</summary>
        XX,

        /// <summary><c>"XY"</c>: X (non-PAR) and Y are haploid — <c>nA = (ρ − 1 + (2(1 − ρ) + ρψ)·2^(r/γ))/ρ</c>, <c>nB = 0</c>.</summary>
        XY,
    }

    /// <summary>
    /// Sex-chromosome model of the ASCAT copy-number output (runASCAT, ascat.runAscat.R, VanLoo-lab/ascat). It does
    /// <b>not</b> change the purity/ploidy fit — ASCAT always fits on autosomes only
    /// (<c>autoprobes = !(SNPposhet[,1] %in% sexchromosomes)</c>) — only the integer segments emitted on X and Y and
    /// hence the reported ploidy:
    /// <list type="bullet">
    /// <item><see cref="AscatGender.XX"/> (<see cref="Female"/>, the default): <c>haploidchrs</c> is empty, so X and Y
    /// use the diploid equations (Y is not special-cased: the <c>nullprobes</c> branch is reached only for non-diploid
    /// probes). <see cref="XNonPar"/> is ignored, as in ASCAT (<c>!is.null(X_nonPAR) &amp;&amp; gender == "XY"</c>).</item>
    /// <item><see cref="AscatGender.XY"/>: <c>haploidchrs = c("X", "Y")</c>; a haploid segment gets
    /// <c>nAraw = (rho − 1 + ((1 − rho)·2 + rho·psi)·2^(logR/gamma))/rho</c>, <c>nBraw = 0</c> (normal cells carry one
    /// copy), followed by the usual negative-value correction and R rounding. Without <see cref="XNonPar"/>
    /// (ASCAT <c>genomeVersion = NULL</c>) the whole of X is haploid. With it, ASCAT's
    /// <c>diploidprobes_fixnonPAR</c> applies: an X segment is haploid only when its overlap with the closed interval
    /// [Start, End] exceeds 50 % of the segment's closed span [first probe position, last probe position]
    /// (IRanges widths), otherwise it is diploid (pseudo-autosomal). Y is always haploid.</item>
    /// </list>
    /// Chromosome labels are matched case-insensitively with an optional "chr" prefix (as for the autosomal filter).
    /// The interval is compared with the segment coordinates exactly as given (ASCAT compares it with <c>SNPpos</c>,
    /// which is 1-based; the presets are ASCAT's 1-based constants).
    /// </summary>
    public sealed record AscatSexModel
    {
        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "hg19"</c> (GRCh37): <c>c(2699521, 154931043)</c>.</summary>
        public static (long Start, long End) XNonParHg19 { get; } = (2_699_521, 154_931_043);

        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "hg38"</c> (GRCh38): <c>c(2781480, 155701382)</c>.</summary>
        public static (long Start, long End) XNonParHg38 { get; } = (2_781_480, 155_701_382);

        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "CHM13"</c> (T2T-CHM13): <c>c(2394411, 153925834)</c>.</summary>
        public static (long Start, long End) XNonParChm13 { get; } = (2_394_411, 153_925_834);

        /// <summary>ASCAT default (<c>gender = "XX"</c>): X and Y diploid. Bit-identical to the gender-less entry points.</summary>
        public static AscatSexModel Female { get; } = new(AscatGender.XX, null);

        /// <summary><c>gender = "XY"</c>, <c>X_nonPAR = NULL</c> (ASCAT <c>genomeVersion = NULL</c>): all of X and Y haploid.</summary>
        public static AscatSexModel Male { get; } = new(AscatGender.XY, null);

        /// <summary>The sample sex.</summary>
        public AscatGender Gender { get; }

        /// <summary>Closed X non-PAR interval (ASCAT <c>X_nonPAR</c>); <c>null</c> = whole X haploid in a male.</summary>
        public (long Start, long End)? XNonPar { get; }

        /// <summary>Creates a sex model.</summary>
        /// <param name="gender">The sample sex.</param>
        /// <param name="xNonPar">Optional closed X non-PAR interval (Start ≤ End); used only for <see cref="AscatGender.XY"/>.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="gender"/> is undefined or the interval has Start &gt; End.</exception>
        public AscatSexModel(AscatGender gender, (long Start, long End)? xNonPar = null)
        {
            if (gender is not (AscatGender.XX or AscatGender.XY))
            {
                throw new ArgumentOutOfRangeException(nameof(gender), gender, "Gender must be XX or XY.");
            }

            if (xNonPar is { } interval && interval.Start > interval.End)
            {
                throw new ArgumentOutOfRangeException(nameof(xNonPar), xNonPar, "The X non-PAR interval needs Start ≤ End.");
            }

            Gender = gender;
            XNonPar = xNonPar;
        }

        /// <summary>A male (<c>"XY"</c>) model with ASCAT's X non-PAR interval of <paramref name="genome"/>
        /// (<c>ascat.loadData(genomeVersion = "hg19" | "hg38")</c>).</summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is undefined.</exception>
        public static AscatSexModel MaleWithXNonPar(ReferenceGenome genome) => genome switch
        {
            ReferenceGenome.GRCh37 => new AscatSexModel(AscatGender.XY, XNonParHg19),
            ReferenceGenome.GRCh38 => new AscatSexModel(AscatGender.XY, XNonParHg38),
            _ => throw new ArgumentOutOfRangeException(nameof(genome), genome, "Unknown reference genome."),
        };

        /// <summary>
        /// True when a segment is haploid under this model (runASCAT <c>!diploidprobes</c> at the segment start; ASCAT's
        /// <c>nullprobes</c> branch is unreachable with the default <c>sexchromosomes = c("X", "Y")</c>).
        /// </summary>
        internal bool IsHaploid(string chromosome, long start, long end)
        {
            if (Gender != AscatGender.XY)
            {
                return false;
            }

            ReadOnlySpan<char> name = StripChrPrefix(chromosome.AsSpan().Trim());
            if (name.Equals("Y", StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }

            if (!name.Equals("X", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }

            if (XNonPar is not { } nonPar)
            {
                return true;
            }

            // diploidprobes_fixnonPAR: findOverlaps(nonPAR, SEGMENTS) and width(pintersect)/width(segment) > 0.5.
            double overlap = (double)Math.Min(end, nonPar.End) - Math.Max(start, nonPar.Start) + 1.0;
            return overlap > 0.0 && overlap / ((double)end - start + 1.0) > 0.5;
        }
    }

    /// <summary>
    /// Segments per-locus allele-specific signal (logR, BAF) into contiguous regions, producing one
    /// (mean logR, mirrored BAF) summary per segment, by running the published ASCAT allele-specific segmentation
    /// (<c>ascat.aspcf</c>, VanLoo-lab/ascat ascat.aspcf.R; Nilsen et al. 2012, <i>BMC Genomics</i> 13:591; Ross et al.
    /// 2021, <i>Bioinformatics</i> 37:1909) with ASCAT's default penalty
    /// (<see cref="AspcfDefaultPenalty"/> = 70). It is exactly
    /// <c><see cref="SegmentAlleleSpecificAspcf"/>(loci, <see cref="AspcfDefaultPenalty"/>)</c>.
    /// <para>
    /// Compatibility (B24 F35): this name formerly ran an unsourced greedy left-to-right mean-shift heuristic with
    /// caller-chosen absolute thresholds. That heuristic has no published reference and was replaced; the threshold
    /// parameters below have no ASPCF meaning (ASPCF's cost is MAD-standardised and penalised, with the fixed minimum
    /// segment length kmin = 6) and are <b>ignored</b> — they are kept, and still range-checked, only so existing call
    /// sites keep compiling and keep their argument-validation behaviour. Pass a penalty to
    /// <see cref="SegmentAlleleSpecificAspcf"/> to tune the segmentation.
    /// </para>
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome. LogR must be
    /// finite and BAF in [0, 1] (as <see cref="SegmentAlleleSpecificAspcf"/>).</param>
    /// <param name="logRChangeThreshold">Ignored (former greedy-heuristic logR threshold; no ASPCF equivalent). Must be &gt; 0.</param>
    /// <param name="bafChangeThreshold">Ignored (former greedy-heuristic BAF threshold; no ASPCF equivalent). Must be &gt; 0.</param>
    /// <param name="minLociPerSegment">Ignored (ASPCF uses ASCAT's fixed kmin = 6). Must be ≥ 1.</param>
    /// <returns>The ASPCF segment summaries (mean raw logR, ASPCF mirrored BAF) in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException">a locus has a null chromosome label, a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">a threshold ≤ 0 or NaN, or minLociPerSegment &lt; 1.</exception>
    public static IReadOnlyList<AlleleSpecificSegmentSummary> SegmentAlleleSpecific(
        IEnumerable<AlleleSpecificLocus> loci,
        double logRChangeThreshold,
        double bafChangeThreshold = 0.1,
        int minLociPerSegment = 1)
    {
        ArgumentNullException.ThrowIfNull(loci);

        if (double.IsNaN(logRChangeThreshold) || logRChangeThreshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(logRChangeThreshold), logRChangeThreshold, "The logR change threshold must be positive.");
        }

        if (double.IsNaN(bafChangeThreshold) || bafChangeThreshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(bafChangeThreshold), bafChangeThreshold, "The BAF change threshold must be positive.");
        }

        if (minLociPerSegment < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minLociPerSegment), minLociPerSegment, "At least one locus per segment is required.");
        }

        return SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty);
    }

    /// <summary>
    /// Raw (real-valued) ASCAT allele-specific copy numbers (nA, nB) for one segment given (r, b, ρ, ψ, γ).
    /// Source: ascat.runAscat.R (VanLoo-lab/ascat), verbatim:
    /// <code>
    /// nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
    /// nB = (rho-1 +  b   *2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
    /// </code>
    /// With ASCAT's segmented BAF convention b ≤ 0.5 (<c>Tumor_BAF_segmented = 1 − bafPCFed</c>) nA is the major
    /// and nB the minor allele.
    /// </summary>
    private static (double NA, double NB) AscatRawCopyNumbers(double r, double b, double rho, double psi, double gamma)
    {
        double scaledTotal = Math.Pow(2.0, r / gamma) * MixtureCopiesPerCell(rho, psi);
        double nA = (rho - 1.0 - (b - 1.0) * scaledTotal) / rho;
        double nB = (rho - 1.0 + b * scaledTotal) / rho;
        return (nA, nB);
    }

    /// <summary>
    /// Raw ASCAT copy number of a haploid (male X non-PAR / Y) segment, ascat.runAscat.R <c>seg_raw</c>, verbatim:
    /// <c>nAraw = (rho-1 + ((1-rho)*2+rho*psi)*2^(logR/gamma))/rho</c> (the normal cells contribute one copy).
    /// </summary>
    private static double AscatRawHaploidCopyNumber(double r, double rho, double psi, double gamma) =>
        (rho - 1.0 + Math.Pow(2.0, r / gamma) * MixtureCopiesPerCell(rho, psi)) / rho;

    /// <summary>
    /// Average number of copies of a locus per cell in a tumour sample of purity ρ whose tumour cells carry
    /// <paramref name="tumorCopies"/> copies and whose normal cells are diploid: <c>2(1 − ρ) + ρ·n</c>. This mixture
    /// denominator is shared by the ASCAT logR model (<c>(1 − rho)·2 + rho·psi</c>, ascat.runAscat.R), Battenberg
    /// (<c>psi = rho·psit + 2(1 − rho)</c>) and the Landau/CNAqc expected-VAF model (<c>f = ρ·M·c / (2(1 − ρ) + ρ·q)</c>,
    /// <see cref="ClassifyClonality"/>). IEEE addition and multiplication are commutative, so every caller is
    /// bit-identical to its former inline form.
    /// </summary>
    private static double MixtureCopiesPerCell(double purity, double tumorCopies) =>
        NormalDiploidCopyNumber * (1.0 - purity) + purity * tumorCopies;

    // ---- ASCAT runASCAT solution-selection constants (ascat.runAscat.R, verbatim) ----

    /// <summary>ASCAT <c>MINRHO = 0.2</c>: minimum aberrant-cell fraction of an accepted optimum.</summary>
    private const double AscatMinRho = 0.2;

    /// <summary>ASCAT <c>MINGOODNESSOFFIT = 80</c> (%): minimum goodness of fit of an accepted optimum.</summary>
    private const double AscatMinGoodnessOfFit = 80.0;

    /// <summary>ASCAT <c>MINPERCZERO = 0.02</c>: minimum fraction of allele copies rounded to 0 (pass 1).</summary>
    private const double AscatMinPercentZero = 0.02;

    /// <summary>ASCAT <c>MINPERCZEROABB = 0.1</c>: minimum zero-allele fraction over aberrant segments (passes 2–3).</summary>
    private const double AscatMinPercentZeroAberrant = 0.1;

    /// <summary>ASCAT <c>MINPERCODDEVEN = 0.05</c>: minimum fraction of odd/even allele pairs (pass 3).</summary>
    private const double AscatMinPercentOddEven = 0.05;

    /// <summary>ASCAT <c>MINPLOIDYSTRICT = 1.7</c>: strict lower ploidy bound of the fallback passes 2 and 4.</summary>
    private const double AscatMinPloidyStrict = 1.7;

    /// <summary>ASCAT <c>MAXPLOIDYSTRICT = 2.3</c>: strict upper ploidy bound of the fallback passes 2 and 4.</summary>
    private const double AscatMaxPloidyStrict = 2.3;

    /// <summary>ASCAT <c>MINABB = 0.03</c>: a sample with ≤ 3 % aberrant (BAF ≠ 0.5) probes may be non-aberrant.</summary>
    private const double AscatMinAberrantFraction = 0.03;

    /// <summary>ASCAT <c>MINABBREGION = 0.005</c>: … and no aberrant segment larger than 0.5 % of the probes.</summary>
    private const double AscatMinAberrantRegionFraction = 0.005;

    /// <summary>ASCAT local-minimum window half-width: <c>seld = d[(i-3):(i+3), (j-3):(j+3)]</c> (a 7 × 7 window).</summary>
    private const int AscatLocalMinimumHalfWindow = 3;

    /// <summary>ASCAT pass-3 mask value for grid columns with ρ &gt; 1: <c>d[, cold] = 1E20</c>.</summary>
    private const double AscatMaskedDistance = 1e20;

    /// <summary>ASCAT <c>limitround = 0.5</c>: odd-total evidence threshold for balanced (BAF = 0.5) segments.</summary>
    private const double AscatLimitRound = 0.5;

    /// <summary>
    /// One autosomal ASCAT fitting segment (<c>make_segments</c> row): segmented logR r, segmented BAF b in ASCAT's
    /// ≤ 0.5 orientation, and <c>length</c> = the number of germline-heterozygous probes of the segment.
    /// </summary>
    private readonly record struct AscatFitSegment(double R, double B, double Length);

    /// <summary>
    /// True when <paramref name="chromosome"/> is one of ASCAT's default <c>sexchromosomes = c("X", "Y")</c>
    /// (an optional "chr" prefix is ignored, case-insensitive). ASCAT excludes these probes from the purity/ploidy
    /// fit: <c>autoprobes = !(SNPposhet[,1] %in% sexchromosomes)</c>.
    /// </summary>
    private static bool IsAscatSexChromosome(string chromosome)
    {
        ReadOnlySpan<char> name = StripChrPrefix(chromosome.AsSpan().Trim());
        return name.Equals("X", StringComparison.OrdinalIgnoreCase) || name.Equals("Y", StringComparison.OrdinalIgnoreCase);
    }

    /// <summary>
    /// Input contract of the ASCAT / ASPCF / Battenberg ports for one (logR, BAF) pair — per locus or per segment
    /// summary: a finite logR and a BAF in [0, 1].
    /// </summary>
    private static bool IsValidAlleleSignal(double logR, double baf) =>
        double.IsFinite(logR) && !double.IsNaN(baf) && baf >= 0.0 && baf <= 1.0;

    /// <summary>Mirrored BAF max(b, 1 − b) ∈ [0.5, 1] (ascat.aspcf <c>ifelse(b &gt; 0.5, b, 1 − b)</c>; Battenberg <c>l</c>).</summary>
    private static double MirrorBaf(double baf) => baf > BalancedBaf ? baf : 1.0 - baf;

    /// <summary>ASCAT segmented-BAF orientation: <c>Tumor_BAF_segmented = 1 − mirroredBAF</c> ∈ [0, 0.5].</summary>
    private static double ToAscatBaf(double baf) => baf > BalancedBaf ? 1.0 - baf : baf;

    /// <summary>R <c>.Machine$double.eps</c> = 2⁻⁵² (used by <c>seq.default</c>).</summary>
    private const double RMachineEpsilon = 2.220446049250313e-16;

    /// <summary>R <c>round()</c> (IEC 60559 round-half-to-even), as used throughout ascat.runAscat.R.</summary>
    private static double RRound(double x) => Math.Round(x, MidpointRounding.ToEven);

    /// <summary>R floored modulus <c>x %% 2</c> for an integral double (0 or 1, also for negative x).</summary>
    private static double RMod2(double x)
    {
        double m = x % 2.0;
        return m < 0.0 ? m + 2.0 : m;
    }

    /// <summary>
    /// R <c>seq(from, to, by)</c> for by &gt; 0 (<c>seq.default</c>): <c>n = as.integer((to−from)/by + 1e−10)</c>,
    /// <c>x = from + (0:n)·by</c>, clamped by <c>pmin(x, to)</c>; <c>from</c> alone when the span is negligible.
    /// </summary>
    private static double[] RSeq(double from, double to, double by)
    {
        double del = to - from;
        double scale = Math.Max(Math.Abs(to), Math.Abs(from));
        if (del == 0.0 || (scale > 0.0 && Math.Abs(del) / scale < 100.0 * RMachineEpsilon))
        {
            return new[] { from };
        }

        int n = (int)(del / by + 1e-10);
        var values = new double[n + 1];
        for (int k = 0; k <= n; k++)
        {
            values[k] = Math.Min(from + k * by, to);
        }

        return values;
    }

    /// <summary>
    /// The value R obtains from <c>as.numeric(rownames(d)[i])</c>: the grid value printed with 15 significant digits
    /// and parsed back (ASCAT reads the candidate ψ/ρ from the distance-matrix dimnames).
    /// </summary>
    private static double RDimnameValue(double value) =>
        double.Parse(value.ToString("G15", System.Globalization.CultureInfo.InvariantCulture),
            System.Globalization.CultureInfo.InvariantCulture);

    /// <summary>
    /// ASCAT <c>create_distance_matrix</c> cell: nA/nB at (ρ, ψ), the minor allele chosen genome-wide
    /// (<c>if (sum(nA) &lt; sum(nB)) nMinor = nA else nMinor = nB</c>), and
    /// <c>d = sum(abs(nMinor − pmax(round(nMinor), 0))^2 · length · ifelse(b == 0.5, 0.05, 1), na.rm = TRUE)</c>.
    /// </summary>
    private static double AscatDistance(AscatFitSegment[] segs, double rho, double psi, double gamma)
    {
        double sumA = 0.0, sumB = 0.0;
        var nA = new double[segs.Length];
        var nB = new double[segs.Length];
        for (int i = 0; i < segs.Length; i++)
        {
            (nA[i], nB[i]) = AscatRawCopyNumbers(segs[i].R, segs[i].B, rho, psi, gamma);
            if (!double.IsNaN(nA[i])) sumA += nA[i];
            if (!double.IsNaN(nB[i])) sumB += nB[i];
        }

        double[] nMinor = sumA < sumB ? nA : nB;
        double d = 0.0;
        for (int i = 0; i < segs.Length; i++)
        {
            double dev = Math.Abs(nMinor[i] - Math.Max(RRound(nMinor[i]), 0.0));
            double term = dev * dev * segs[i].Length * (segs[i].B == BalancedBaf ? AscatBalancedSegmentWeight : 1.0);
            if (!double.IsNaN(term))
            {
                d += term; // na.rm = TRUE
            }
        }

        return d;
    }

    /// <summary>
    /// Per-candidate statistics of one ASCAT local minimum (ascat.runAscat.R): the recomputed ploidy
    /// <c>sum((nA+nB)·length)/sum(length)</c>, <c>percentzero</c>, <c>perczeroAbb</c> (NaN → 0) and <c>percOddEven</c>.
    /// </summary>
    private static (double Ploidy, double PercentZero, double PercentZeroAberrant, double PercentOddEven) AscatCandidateStatistics(
        AscatFitSegment[] segs, double rho, double psi, double gamma)
    {
        double totalLength = 0.0, ploidyNumerator = 0.0, zero = 0.0, zeroAbb = 0.0, abbLength = 0.0, oddEven = 0.0;
        foreach (AscatFitSegment s in segs)
        {
            (double nA, double nB) = AscatRawCopyNumbers(s.R, s.B, rho, psi, gamma);
            double rA = RRound(nA), rB = RRound(nB);
            double aberrant = s.B == BalancedBaf ? 0.0 : 1.0;
            totalLength += s.Length;
            ploidyNumerator += (nA + nB) * s.Length;
            double zeros = (rA == 0.0 ? 1.0 : 0.0) + (rB == 0.0 ? 1.0 : 0.0);
            zero += zeros * s.Length;
            zeroAbb += zeros * s.Length * aberrant;
            abbLength += s.Length * aberrant;
            double modA = RMod2(rA), modB = RMod2(rB);
            if ((modA == 0.0 && modB == 1.0) || (modA == 1.0 && modB == 0.0))
            {
                oddEven += s.Length;
            }
        }

        double percentZeroAberrant = zeroAbb / abbLength;
        if (double.IsNaN(percentZeroAberrant))
        {
            percentZeroAberrant = 0.0; // "the next can happen if BAF is a flat line at 0.5"
        }

        return (ploidyNumerator / totalLength, zero / totalLength, percentZeroAberrant, oddEven / totalLength);
    }

    /// <summary>
    /// Integer allele-specific copy number of one segment at the chosen (ρ, ψ), ported verbatim from the
    /// <c>seg_raw</c> construction of ascat.runAscat.R: raw nA/nB, the negative-value correction
    /// (<c>nA+nB &lt; 0 ⇒ 0,0</c>; a negative allele is folded into the other), R half-to-even rounding, and the
    /// balanced-segment odd-total rule (<c>limitround = 0.5</c>: for BAF = 0.5, if nA+nB exceeds the rounded total by
    /// more than 0.5 nA is raised by one; if it falls short by more than 0.5 nB is lowered by one).
    /// </summary>
    /// <param name="r">Segmented logR.</param>
    /// <param name="bAscat">Segmented BAF in ASCAT's ≤ 0.5 orientation (0 when the segment has no heterozygous probe).</param>
    /// <param name="rho">Purity ρ.</param>
    /// <param name="psi">Ploidy parameter ψ.</param>
    /// <param name="gamma">Platform parameter γ.</param>
    /// <param name="haploid">Male X (non-PAR) / Y segment (<see cref="AscatSexModel"/>): runASCAT's non-diploid branch
    /// <c>nAraw = (rho − 1 + ((1 − rho)·2 + rho·psi)·2^(logR/gamma))/rho</c>, <c>nBraw = 0</c>.</param>
    private static (double Major, double Minor) AscatRoundSegment(
        double r, double bAscat, double rho, double psi, double gamma, bool haploid = false)
    {
        (double nAraw, double nBraw) = haploid
            ? (AscatRawHaploidCopyNumber(r, rho, psi, gamma), 0.0)
            : AscatRawCopyNumbers(r, bAscat, rho, psi, gamma);
        if (nAraw + nBraw < 0.0)
        {
            nAraw = 0.0;
            nBraw = 0.0;
        }
        else if (nAraw < 0.0)
        {
            nBraw = nAraw + nBraw;
            nAraw = 0.0;
        }
        else if (nBraw < 0.0)
        {
            nAraw = nAraw + nBraw;
            nBraw = 0.0;
        }

        double rA = RRound(nAraw), rB = RRound(nBraw);
        double nA = rA, nB = rB;
        if (bAscat == BalancedBaf)
        {
            if (nAraw + nBraw > rA + rB + AscatLimitRound)
            {
                nA = rA + 1.0;
            }
            else if (nAraw + nBraw < rA + rB - AscatLimitRound)
            {
                nB = rB - 1.0;
            }
        }

        return (nA, nB);
    }

    /// <summary>Converts an ASCAT integer copy number (≥ 0) to <see cref="int"/>, saturating at Int32.MaxValue (NaN → 0).</summary>
    private static int AscatCopyNumberToInt(double value)
    {
        if (double.IsNaN(value) || value <= 0.0)
        {
            return 0;
        }

        return value >= int.MaxValue ? int.MaxValue : (int)value;
    }

    /// <summary>
    /// Builds the fit result at the selected (ρ, ψ): the integer major/minor segments (ASCAT <c>seg_raw</c> nMajor/nMinor,
    /// one per input summary, all chromosomes), the ASCAT output ploidy — the mean integer total copy number over the
    /// heterozygous probes (<c>ploidy = mean(nA + nB)</c>, here weighted by <c>LocusCount</c>) — and the given GoF.
    /// </summary>
    private static PurityPloidyFit BuildAscatFit(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments, double rho, double psi, double gamma, double goodnessOfFit,
        bool nonAberrant, AscatSexModel sexModel)
    {
        var result = new List<AlleleSpecificSegment>(segments.Count);
        var probeCounts = new int[segments.Count];
        for (int i = 0; i < segments.Count; i++)
        {
            AlleleSpecificSegmentSummary s = segments[i];
            (double major, double minor) = AscatRoundSegment(s.MeanLogR, ToAscatBaf(s.MeanBAF), rho, psi, gamma,
                sexModel.IsHaploid(s.Chromosome, s.Start, s.End));
            int majorInt = AscatCopyNumberToInt(major);
            int minorInt = AscatCopyNumberToInt(minor);
            // Segments with End == Start (single-position) get a 1 bp span so AlleleSpecificSegment.Length > 0.
            long end = s.End > s.Start ? s.End : s.Start + 1;
            result.Add(new AlleleSpecificSegment(s.Chromosome, s.Start, end, majorInt, minorInt));
            probeCounts[i] = s.LocusCount;
        }

        // ASCAT ploidy = mean(nA + nB) over probes = probe-count-weighted mean (canonical overload).
        return new PurityPloidyFit(rho, EstimatePloidy(result, probeCounts), goodnessOfFit, result)
        {
            Psi = psi,
            IsNonAberrant = nonAberrant,
        };
    }

    /// <summary>Validates the segment summaries consumed by the ASCAT fit and returns the autosomal fitting segments.</summary>
    private static AscatFitSegment[] PrepareAscatSegments(IReadOnlyList<AlleleSpecificSegmentSummary> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);
        if (segments.Count == 0)
        {
            throw new ArgumentException("At least one segment is required to fit purity and ploidy.", nameof(segments));
        }

        var autosomal = new List<AscatFitSegment>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            if (s.Chromosome is null)
            {
                throw new ArgumentException("A segment has a null chromosome label.", nameof(segments));
            }

            if (!IsValidAlleleSignal(s.MeanLogR, s.MeanBAF))
            {
                throw new ArgumentException(
                    "Every segment needs a finite mean logR and a mean BAF in [0, 1].", nameof(segments));
            }

            if (s.LocusCount < 1)
            {
                throw new ArgumentException(
                    "Every segment must summarise at least one heterozygous locus (ASCAT weights segments by probe count).",
                    nameof(segments));
            }

            if (!IsAscatSexChromosome(s.Chromosome))
            {
                autosomal.Add(new AscatFitSegment(s.MeanLogR, ToAscatBaf(s.MeanBAF), s.LocusCount));
            }
        }

        if (autosomal.Count == 0)
        {
            throw new ArgumentException(
                "ASCAT fits purity and ploidy on autosomal segments only (sex chromosomes X/Y are excluded); none was supplied.",
                nameof(segments));
        }

        return autosomal.ToArray();
    }

    /// <summary>
    /// Jointly estimates tumour purity ρ and ploidy ψ from segment-level (logR, BAF) summaries with the ASCAT
    /// algorithm (Van Loo et al. 2010, <i>PNAS</i> 107:16910), ported from <c>runASCAT</c> in ascat.runAscat.R
    /// (VanLoo-lab/ascat):
    /// <list type="number">
    /// <item><b>Distance matrix</b> (<c>create_distance_matrix</c>): for ψ ∈ seq(ploidyMin − 0.5, ploidyMax + 0.5,
    /// ploidyStep) × ρ ∈ seq(purityMin, purityMax, purityStep), d(ψ, ρ) = Σ (nMinor − max(round(nMinor), 0))² ·
    /// length · w_b over the <b>autosomal</b> segments, where length = number of heterozygous probes
    /// (<see cref="AlleleSpecificSegmentSummary.LocusCount"/>, ASCAT <c>make_segments</c>) and w_b = 0.05 for BAF = 0.5.</item>
    /// <item><b>Local minima</b>: a grid cell is a candidate when it is the strict minimum of its 7 × 7 neighbourhood.</item>
    /// <item><b>Filter cascade</b> (first pass that yields a candidate wins): (1) non-aberrant sample excluded,
    /// ploidyMin &lt; ploidy &lt; ploidyMax, ρ ≥ 0.2, GoF &gt; 80 %, percentzero &gt; 0.02; (2) 1.7 &lt; ploidy &lt; 2.3 and
    /// perczeroAbb &gt; 0.1; (3) grid columns with ρ &gt; 1 masked (so ρ = 1 can be a border optimum), percentzero /
    /// perczeroAbb / percOddEven alternatives; (4) 1.7 &lt; ploidy &lt; 2.3 only. Here ploidy = Σ(nA+nB)·length/Σlength.</item>
    /// <item><b>Selection</b>: the candidate with the smallest distance (the last one in grid order on an exact tie);
    /// ρ &gt; 1 is reported as 1. GoF = (1 − d/TheoretMaxdist)·100 with TheoretMaxdist = Σ 0.25·length·w_b.</item>
    /// <item><b>Integer segments</b> (ASCAT <c>seg_raw</c>): negative-value correction, R half-to-even rounding and the
    /// balanced odd-total rule, see <see cref="AscatRoundSegment"/>.</item>
    /// </list>
    /// Sex-chromosome (X/Y) segments are excluded from the fit (as in ASCAT, for either sex) and emitted with the diploid
    /// model (ASCAT gender "XX"); pass an <see cref="AscatSexModel"/> (overload, B24 F38) for a male sample.
    /// </summary>
    /// <param name="segments">Segment summaries (from <see cref="SegmentAlleleSpecificAspcf"/> or a caller's segmenter). Non-empty;
    /// each needs a finite mean logR, a mean BAF in [0, 1] and LocusCount ≥ 1; at least one autosomal segment.</param>
    /// <param name="purityMin">Lower bound of the purity grid, in (0, 1] (ASCAT <c>min_purity</c> = 0.1).</param>
    /// <param name="purityMax">Upper bound of the purity grid, finite and ≥ purityMin (ASCAT <c>max_purity</c> = 1.05;
    /// grid points above 1 let ρ = 1 be an interior optimum and are reported as ρ = 1).</param>
    /// <param name="purityStep">Purity grid step (&gt; 0; ASCAT 0.01).</param>
    /// <param name="ploidyMin">ASCAT <c>min_ploidy</c> (&gt; 0, default 1.5): lower ploidy filter; the ψ grid starts at ploidyMin − 0.5.</param>
    /// <param name="ploidyMax">ASCAT <c>max_ploidy</c> (≥ ploidyMin, default 5.5): upper ploidy filter; the ψ grid ends at ploidyMax + 0.5.</param>
    /// <param name="ploidyStep">Ploidy grid step (&gt; 0; ASCAT 0.05).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>The recovered ρ, the ASCAT output ploidy, the percentage GoF, the integer segments, and ψ
    /// (<see cref="PurityPloidyFit.Psi"/>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum ("ASCAT could not find an optimal
    /// ploidy and purity value"); use <see cref="TryFitPurityPloidy"/> to test without an exception.</exception>
    public static PurityPloidyFit FitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => FitPurityPloidy(segments, AscatSexModel.Female, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// ASCAT purity/ploidy fit (see <see cref="FitPurityPloidy(IReadOnlyList{AlleleSpecificSegmentSummary}, double, double, double, double, double, double, double)"/>)
    /// with an explicit sex-chromosome model (runASCAT <c>gender</c> / <c>X_nonPAR</c>, B24 F38). The fit (ρ, ψ, GoF,
    /// non-aberrant flag) is unchanged — sex chromosomes never enter it; only the X/Y integer segments and therefore the
    /// reported ploidy depend on <paramref name="sexModel"/> (see <see cref="AscatSexModel"/>). Each summary is one
    /// runASCAT segment; its [Start, End] is the span tested against the X non-PAR interval.
    /// </summary>
    /// <param name="segments">Segment summaries (as the gender-less overload).</param>
    /// <param name="sexModel">Sex-chromosome model; <see cref="AscatSexModel.Female"/> reproduces the gender-less overload.</param>
    /// <param name="purityMin">See the gender-less overload.</param>
    /// <param name="purityMax">See the gender-less overload.</param>
    /// <param name="purityStep">See the gender-less overload.</param>
    /// <param name="ploidyMin">See the gender-less overload.</param>
    /// <param name="ploidyMax">See the gender-less overload.</param>
    /// <param name="ploidyStep">See the gender-less overload.</param>
    /// <param name="gamma">See the gender-less overload.</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum.</exception>
    public static PurityPloidyFit FitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        if (!TryFitPurityPloidy(segments, sexModel, out PurityPloidyFit fit, purityMin, purityMax, purityStep,
                ploidyMin, ploidyMax, ploidyStep, gamma))
        {
            throw new InvalidOperationException(
                "ASCAT could not find an optimal ploidy and purity value: no local minimum of the distance matrix passes " +
                "the ASCAT solution filters (ascat.runAscat.R).");
        }

        return fit;
    }

    /// <summary>
    /// ASCAT purity/ploidy fit (see <see cref="FitPurityPloidy"/>) that reports failure instead of throwing:
    /// returns <c>false</c> (and a default <paramref name="fit"/>) when no local minimum of the distance matrix passes
    /// the ASCAT filters — ASCAT's <c>rho = NA</c> outcome. Argument errors still throw.
    /// </summary>
    public static bool TryFitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => TryFitPurityPloidy(segments, AscatSexModel.Female, out fit, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidy(IReadOnlyList{AlleleSpecificSegmentSummary}, AscatSexModel, double, double, double, double, double, double, double)"/>
    /// reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing (B24 F38).
    /// </summary>
    public static bool TryFitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);
        if (!TryFindAscatOptimum(s, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma,
                out double rhoOpt, out double psiOpt, out double goodnessOfFit, out bool nonAberrant))
        {
            fit = default;
            return false;
        }

        fit = BuildAscatFit(segments, rhoOpt, psiOpt, gamma, goodnessOfFit, nonAberrant, sexModel);
        return true;
    }

    /// <summary>
    /// The runASCAT grid search (ascat.runAscat.R) on validated autosomal fitting segments: distance matrix, strict
    /// 7 × 7 local minima, the four-pass filter cascade and the optimum selection (see <see cref="FitPurityPloidy"/>).
    /// Returns <c>false</c> when no candidate passes (ASCAT <c>rho = NA</c>); ρ &gt; 1 is reported as 1.
    /// </summary>
    private static bool TryFindAscatOptimum(
        AscatFitSegment[] s, double purityMin, double purityMax, double purityStep,
        double ploidyMin, double ploidyMax, double ploidyStep, double gamma,
        out double rhoOpt, out double psiOpt, out double goodnessOfFitOpt, out bool nonAberrant)
    {
        double[] psiPos = RSeq(ploidyMin - 0.5, ploidyMax + 0.5, ploidyStep);
        double[] rhoPos = RSeq(purityMin, purityMax, purityStep);
        int rows = psiPos.Length, cols = rhoPos.Length;
        var d = new double[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                d[i, j] = AscatDistance(s, rhoPos[j], psiPos[i], gamma);
            }
        }

        (double theoreticalMaxDistance, nonAberrant) = AscatSampleSummary(s);
        bool sampleNonAberrant = nonAberrant;

        var candidates = new List<(double M, int I, int J, double GoodnessOfFit)>();
        bool strictPloidyWindowReachable = ploidyMin < AscatMaxPloidyStrict && ploidyMax > AscatMinPloidyStrict;

        // Pass 1: all filters.
        CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
            !sampleNonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
            && st.GoodnessOfFit > AscatMinGoodnessOfFit && st.PercentZero > AscatMinPercentZero);

        // Pass 2: drop percentzero (allow non-aberrant solutions) with strict ploidy borders.
        if (candidates.Count == 0 && strictPloidyWindowReachable)
        {
            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                st.Ploidy > AscatMinPloidyStrict && st.Ploidy < AscatMaxPloidyStrict && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit && st.PercentZeroAberrant > AscatMinPercentZeroAberrant);
        }

        // Pass 3: allow 100 % aberrant cells — mask the rho > 1 columns so rho = 1 can be a (border) minimum.
        if (candidates.Count == 0)
        {
            for (int j = 0; j < cols; j++)
            {
                if (RDimnameValue(rhoPos[j]) > 1.0)
                {
                    for (int i = 0; i < rows; i++)
                    {
                        d[i, j] = AscatMaskedDistance;
                    }
                }
            }

            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                !sampleNonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit
                && (st.PercentZeroAberrant > AscatMinPercentZeroAberrant || st.PercentZero > AscatMinPercentZero
                    || st.PercentOddEven > AscatMinPercentOddEven));
        }

        // Pass 4: drop the percentzero filters, strict ploidy borders.
        if (candidates.Count == 0 && strictPloidyWindowReachable)
        {
            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                st.Ploidy > AscatMinPloidyStrict && st.Ploidy < AscatMaxPloidyStrict && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit);
        }

        if (candidates.Count == 0)
        {
            rhoOpt = psiOpt = goodnessOfFitOpt = double.NaN;
            return false;
        }

        // optlim = sort(localmin)[1]; the loop keeps the LAST optimum whose distance equals optlim.
        double optimum = double.PositiveInfinity;
        foreach (var c in candidates)
        {
            optimum = Math.Min(optimum, c.M);
        }

        (double M, int I, int J, double GoodnessOfFit) best = default;
        foreach (var c in candidates)
        {
            if (c.M == optimum)
            {
                best = c;
            }
        }

        psiOpt = RDimnameValue(psiPos[best.I]);
        rhoOpt = Math.Min(1.0, RDimnameValue(rhoPos[best.J])); // if (rho_opt1 > 1) rho_opt1 = 1
        goodnessOfFitOpt = best.GoodnessOfFit;
        return true;
    }

    /// <summary>
    /// Sample-level ASCAT quantities shared by the grid fit (<see cref="TryFitPurityPloidy"/>) and the manual fit
    /// (<see cref="EvaluatePurityPloidy"/>), ascat.runAscat.R: <c>TheoretMaxdist = sum(rep(0.25, n) · length ·
    /// ifelse(b == 0.5, 0.05, 1))</c> (the GoF normaliser) and the <c>nonaberrant</c> flag (MINABB: aberrant probes ≤ 3 %,
    /// MINABBREGION: no aberrant segment above 0.5 % of the probes).
    /// </summary>
    private static (double TheoreticalMaxDistance, bool NonAberrant) AscatSampleSummary(AscatFitSegment[] s)
    {
        double theoreticalMaxDistance = 0.0, totalLength = 0.0, aberrantLength = 0.0, maxAberrantSegment = 0.0;
        foreach (AscatFitSegment seg in s)
        {
            bool balanced = seg.B == BalancedBaf;
            theoreticalMaxDistance += AscatWorstCaseIntegerDistance * seg.Length * (balanced ? AscatBalancedSegmentWeight : 1.0);
            totalLength += seg.Length;
            if (!balanced)
            {
                aberrantLength += seg.Length;
                maxAberrantSegment = Math.Max(maxAberrantSegment, seg.Length);
            }
        }

        bool nonAberrant = aberrantLength / totalLength <= AscatMinAberrantFraction
                           && maxAberrantSegment / totalLength <= AscatMinAberrantRegionFraction;
        return (theoreticalMaxDistance, nonAberrant);
    }

    /// <summary>Statistics of one ASCAT local-minimum candidate, as tested by the runASCAT filter passes.</summary>
    private readonly record struct AscatCandidate(
        double Rho, double Ploidy, double GoodnessOfFit, double PercentZero, double PercentZeroAberrant, double PercentOddEven);

    /// <summary>
    /// One runASCAT filter pass: scans the interior of the distance matrix (i ∈ 4..nrow−3, j ∈ 4..ncol−3 in R's 1-based
    /// indices) for cells that are the strict minimum of their 7 × 7 window (<c>seld[4,4] = max(seld); min(seld) &gt; m</c>)
    /// and appends those whose statistics pass <paramref name="accept"/>.
    /// </summary>
    private static void CollectAscatOptima(
        double[,] d, double[] psiPos, double[] rhoPos, AscatFitSegment[] segs, double gamma, double theoreticalMaxDistance,
        List<(double M, int I, int J, double GoodnessOfFit)> candidates, Func<AscatCandidate, bool> accept)
    {
        int rows = d.GetLength(0), cols = d.GetLength(1);
        const int h = AscatLocalMinimumHalfWindow;
        for (int i = h; i < rows - h; i++)
        {
            for (int j = h; j < cols - h; j++)
            {
                double m = d[i, j];
                double windowMax = double.NegativeInfinity;
                for (int di = -h; di <= h; di++)
                {
                    for (int dj = -h; dj <= h; dj++)
                    {
                        windowMax = Math.Max(windowMax, d[i + di, j + dj]);
                    }
                }

                bool strictMinimum = true;
                for (int di = -h; di <= h && strictMinimum; di++)
                {
                    for (int dj = -h; dj <= h; dj++)
                    {
                        double v = di == 0 && dj == 0 ? windowMax : d[i + di, j + dj];
                        if (!(v > m))
                        {
                            strictMinimum = false;
                            break;
                        }
                    }
                }

                if (!strictMinimum)
                {
                    continue;
                }

                double psi = RDimnameValue(psiPos[i]);
                double rho = RDimnameValue(rhoPos[j]);
                var stats = AscatCandidateStatistics(segs, rho, psi, gamma);
                double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
                var candidate = new AscatCandidate(
                    rho, stats.Ploidy, goodnessOfFit, stats.PercentZero, stats.PercentZeroAberrant, stats.PercentOddEven);
                if (accept(candidate))
                {
                    candidates.Add((m, i, j, goodnessOfFit));
                }
            }
        }
    }

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (<c>rho_manual</c> / <c>psi_manual</c> in ascat.runAscat.R): skips the
    /// grid search and returns the goodness of fit — d computed as in <c>create_distance_matrix</c> over the autosomal
    /// segments, GoF = (1 − d/TheoretMaxdist)·100 — together with the integer segments and ASCAT ploidy at (ρ, ψ).
    /// </summary>
    /// <param name="segments">Segment summaries (same requirements as <see cref="FitPurityPloidy"/>).</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
        => EvaluatePurityPloidy(segments, AscatSexModel.Female, purity, ploidy, gamma);

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (as the gender-less overload) and an explicit sex-chromosome model
    /// for the emitted X/Y segments (<see cref="AscatSexModel"/>, B24 F38).
    /// </summary>
    /// <param name="segments">Segment summaries.</param>
    /// <param name="sexModel">Sex-chromosome model.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);
        double m = AscatDistance(s, purity, ploidy, gamma);
        double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
        return BuildAscatFit(segments, purity, ploidy, gamma, goodnessOfFit, nonAberrant, sexModel);
    }

    /// <summary>
    /// ASCAT purity/ploidy fit (runASCAT, ascat.runAscat.R) on a germline-aware ASPCF segmentation
    /// (<see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/>), i.e.
    /// with germline-homozygous probes, exactly as runASCAT treats them:
    /// <list type="bullet">
    /// <item><b>Fit</b>: the distance matrix, filters, goodness of fit and non-aberrant flag use the heterozygous
    /// autosomal probes only (<c>r = lrrsegmented[names(bafsegmented)]</c>), grouped by <c>make_segments</c> into runs of
    /// identical (segmented logR, segmented BAF) with length = number of heterozygous probes. Homozygous probes enter
    /// only through the logR levels; the grid search is the one of <see cref="FitPurityPloidy"/>.</item>
    /// <item><b>Segments</b> (<c>seg_raw</c>): one per <see cref="AspcfSegmentation.Segments"/> entry (runs of equal
    /// segmented logR within a chromosome), using the BAF of the run's first heterozygous probe; a run without
    /// heterozygous probes (BAF NA) uses ASCAT's <c>bafke = 0</c>, which after the negative-value correction puts the
    /// whole total copy number on the major allele (only nA + nB is meaningful there).</item>
    /// <item><b>Ploidy</b>: ASCAT's reported <c>ploidy = mean(nA + nB)</c> over <b>all</b> probes (heterozygous and
    /// homozygous, sex chromosomes included), i.e. the <see cref="AspcfSegment.LocusCount"/>-weighted mean total copy
    /// number of the segments.</item>
    /// </list>
    /// Parameters and the no-optimum behaviour are those of <see cref="FitPurityPloidy"/>.
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="purityMin">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="purityMax">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="purityStep">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyMin">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyMax">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyStep">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="gamma">See <see cref="FitPurityPloidy"/>.</param>
    /// <returns>ρ, the all-probe ASCAT ploidy, GoF, the <c>seg_raw</c> integer segments and ψ.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum (use
    /// <see cref="TryFitPurityPloidyFromAspcf"/>).</exception>
    public static PurityPloidyFit FitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => FitPurityPloidyFromAspcf(segmentation, AscatSexModel.Female, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf(AspcfSegmentation, double, double, double, double, double, double, double)"/>
    /// with an explicit sex-chromosome model (runASCAT <c>gender</c> / <c>X_nonPAR</c>, B24 F38): the fit is unchanged
    /// (sex chromosomes never enter it); the X/Y <c>seg_raw</c> segments — and so the all-probe ploidy — follow
    /// <paramref name="sexModel"/> (see <see cref="AscatSexModel"/>). For the X non-PAR rule a segment's span is its
    /// first/last probe position (<see cref="AspcfSegment.Start"/>/<see cref="AspcfSegment.End"/>, i.e. the run of
    /// equal segmented logR that <c>diploidprobes_fixnonPAR</c> builds with <c>rle</c>).
    /// <para>Not ported: for a male with <c>X_nonPAR</c>, <c>ascat.aspcf</c> additionally re-labels the germline
    /// genotypes of non-PAR X probes (all homozygous, then a random autosome-matched fraction heterozygous); here the
    /// caller's <see cref="AspcfSegmentation.GermlineHeterozygous"/> flags are used as given.</para>
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="sexModel">Sex-chromosome model; <see cref="AscatSexModel.Female"/> reproduces the gender-less overload.</param>
    /// <param name="purityMin">See the gender-less overload.</param>
    /// <param name="purityMax">See the gender-less overload.</param>
    /// <param name="purityStep">See the gender-less overload.</param>
    /// <param name="ploidyMin">See the gender-less overload.</param>
    /// <param name="ploidyMax">See the gender-less overload.</param>
    /// <param name="ploidyStep">See the gender-less overload.</param>
    /// <param name="gamma">See the gender-less overload.</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum.</exception>
    public static PurityPloidyFit FitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        if (!TryFitPurityPloidyFromAspcf(segmentation, sexModel, out PurityPloidyFit fit, purityMin, purityMax, purityStep,
                ploidyMin, ploidyMax, ploidyStep, gamma))
        {
            throw new InvalidOperationException(
                "ASCAT could not find an optimal ploidy and purity value: no local minimum of the distance matrix passes " +
                "the ASCAT solution filters (ascat.runAscat.R).");
        }

        return fit;
    }

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf"/> reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing.
    /// </summary>
    public static bool TryFitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => TryFitPurityPloidyFromAspcf(segmentation, AscatSexModel.Female, out fit, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf(AspcfSegmentation, AscatSexModel, double, double, double, double, double, double, double)"/>
    /// reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing (B24 F38).
    /// </summary>
    public static bool TryFitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = AscatMakeSegments(segmentation);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);
        if (!TryFindAscatOptimum(s, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma,
                out double rhoOpt, out double psiOpt, out double goodnessOfFit, out bool nonAberrant))
        {
            fit = default;
            return false;
        }

        fit = BuildAscatFitFromAspcf(segmentation, rhoOpt, psiOpt, gamma, goodnessOfFit, nonAberrant, sexModel);
        return true;
    }

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (<c>rho_manual</c>/<c>psi_manual</c>) on a germline-aware ASPCF
    /// segmentation: GoF over the heterozygous autosomal <c>make_segments</c> runs, <c>seg_raw</c> segments and the
    /// all-probe ploidy as in <see cref="FitPurityPloidyFromAspcf"/>.
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
        => EvaluatePurityPloidyFromAspcf(segmentation, AscatSexModel.Female, purity, ploidy, gamma);

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy on a germline-aware ASPCF segmentation (as the gender-less overload)
    /// with an explicit sex-chromosome model for the X/Y <c>seg_raw</c> segments (<see cref="AscatSexModel"/>, B24 F38).
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="sexModel">Sex-chromosome model.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = AscatMakeSegments(segmentation);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);
        double m = AscatDistance(s, purity, ploidy, gamma);
        double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
        return BuildAscatFitFromAspcf(segmentation, purity, ploidy, gamma, goodnessOfFit, nonAberrant, sexModel);
    }

    /// <summary>
    /// runASCAT fitting segments from a germline-aware segmentation: the heterozygous autosomal probes' (segmented logR,
    /// 1 − mirrored segmented BAF) pairs, grouped by ASCAT <c>make_segments</c> into runs of identical pairs (the runs
    /// ignore chromosome boundaries, as in R); length = number of probes in the run.
    /// </summary>
    private static AscatFitSegment[] AscatMakeSegments(AspcfSegmentation segmentation)
    {
        ArgumentNullException.ThrowIfNull(segmentation);
        var result = new List<AscatFitSegment>();
        double previousR = 1e10, previousB = -1.0;
        int count = 0;
        for (int i = 0; i < segmentation.Loci.Count; i++)
        {
            if (!segmentation.GermlineHeterozygous[i] || IsAscatSexChromosome(segmentation.Loci[i].Chromosome))
            {
                continue;
            }

            double r = segmentation.SegmentedLogR[i];
            double b = 1.0 - segmentation.SegmentedBaf[i]; // Tumor_BAF_segmented = 1 - bafPCFed
            if (b != previousB || r != previousR)
            {
                if (count > 0)
                {
                    result.Add(new AscatFitSegment(previousR, previousB, count));
                }

                count = 0;
            }

            count++;
            previousR = r;
            previousB = b;
        }

        if (count == 0)
        {
            throw new ArgumentException(
                "ASCAT fits purity and ploidy on heterozygous autosomal probes; the segmentation has none.",
                nameof(segmentation));
        }

        result.Add(new AscatFitSegment(previousR, previousB, count));
        return result.ToArray();
    }

    /// <summary>
    /// runASCAT output at the selected (ρ, ψ) for a germline-aware segmentation: <c>seg_raw</c> integer segments (one per
    /// logR segment; <c>bafke</c> = first heterozygous BAF, 0 when none) and <c>ploidy = mean(nA + nB)</c> over all probes.
    /// </summary>
    private static PurityPloidyFit BuildAscatFitFromAspcf(
        AspcfSegmentation segmentation, double rho, double psi, double gamma, double goodnessOfFit, bool nonAberrant,
        AscatSexModel sexModel)
    {
        var result = new List<AlleleSpecificSegment>(segmentation.Segments.Count);
        var probeCounts = new int[segmentation.Segments.Count];
        for (int i = 0; i < segmentation.Segments.Count; i++)
        {
            AspcfSegment s = segmentation.Segments[i];
            double bafke = s.HasBaf ? 1.0 - s.MeanBAF : 0.0; // "if (is.na(bafke)) bafke = 0"
            (double major, double minor) = AscatRoundSegment(s.MeanLogR, bafke, rho, psi, gamma,
                sexModel.IsHaploid(s.Chromosome, s.Start, s.End));
            long end = s.End > s.Start ? s.End : s.Start + 1;
            result.Add(new AlleleSpecificSegment(
                s.Chromosome, s.Start, end, AscatCopyNumberToInt(major), AscatCopyNumberToInt(minor)));
            probeCounts[i] = s.LocusCount;
        }

        return new PurityPloidyFit(rho, EstimatePloidy(result, probeCounts), goodnessOfFit, result)
        {
            Psi = psi,
            IsNonAberrant = nonAberrant,
        };
    }

    /// <summary>
    /// Validates a fixed ASCAT/Battenberg model point (ρ, ψ, γ) — shared by <see cref="EvaluatePurityPloidy"/> and
    /// <see cref="FitSubclonalCopyNumber"/>: ρ ∈ (0, 1], ψ &gt; 0 and γ &gt; 0, all finite (the ASCAT equations divide by ρ
    /// and by γ).
    /// </summary>
    private static void ValidateAscatModelParameters(double purity, double ploidy, double gamma)
    {
        if (!double.IsFinite(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity ρ must be in (0, 1].");
        }

        if (!double.IsFinite(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy ψ must be positive and finite.");
        }

        if (!double.IsFinite(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive and finite.");
        }
    }

    private static void ValidateGrid(
        double purityMin, double purityMax, double purityStep,
        double ploidyMin, double ploidyMax, double ploidyStep, double gamma)
    {
        if (!double.IsFinite(purityMin) || purityMin <= 0.0 || purityMin > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMin), purityMin, "purityMin must be in (0, 1].");
        }

        if (!double.IsFinite(purityMax) || purityMax < purityMin)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMax), purityMax, "purityMax must be finite and ≥ purityMin.");
        }

        if (!double.IsFinite(purityStep) || purityStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityStep), purityStep, "purityStep must be positive.");
        }

        if (!double.IsFinite(ploidyMin) || ploidyMin <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMin), ploidyMin, "ploidyMin must be positive.");
        }

        if (!double.IsFinite(ploidyMax) || ploidyMax < ploidyMin)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMax), ploidyMax, "ploidyMax must be finite and ≥ ploidyMin.");
        }

        if (!double.IsFinite(ploidyStep) || ploidyStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyStep), ploidyStep, "ploidyStep must be positive.");
        }

        if (!double.IsFinite(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive.");
        }

        // Bound the grid so a pathological step cannot allocate an unbounded distance matrix.
        double cells = ((ploidyMax - ploidyMin + 1.0) / ploidyStep + 1.0) * ((purityMax - purityMin) / purityStep + 1.0);
        if (cells > AscatMaxGridCells)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purityStep), purityStep, $"The (ρ, ψ) grid would have {cells:G3} cells (limit {AscatMaxGridCells:G3}).");
        }
    }

    /// <summary>Upper bound on the number of (ρ, ψ) grid cells evaluated (ASCAT default grid: 101 × 96 = 9 696).</summary>
    private const double AscatMaxGridCells = 4_000_000;

    /// <summary>
    /// Derives the integer mutation multiplicity m (number of mutated copies per cancer cell) of a somatic
    /// variant from its VAF, the tumour purity ρ, and the local total / major copy number, so that
    /// <see cref="EstimateCcf"/> can be driven without a caller-supplied multiplicity. The expected number of
    /// mutated copies for a clonal mutation is n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)] (McGranahan et al. 2016,
    /// <i>Science</i> 351:1463; equivalently the inversion of the PICTograph model VAF = m·CCF·ρ /
    /// (N_T·ρ + 2(1−ρ)) at CCF = 1, Zheng et al. 2022, <i>Bioinformatics</i> 38:3677). The result is rounded to
    /// the nearest integer with ties to even and floored at 1, exactly as facets-suite
    /// <c>expected_mutant_copies</c> (mskcc/facets-suite <c>R/ccf-annotate-maf.R</c>, "Based on PMID 28270531":
    /// <c>mu &lt; 1 → 1</c>, then R <c>round</c> = IEC 60559 half-to-even, so n_mut = 2.5 → 2 and 1.5 → 2). In
    /// addition (not in facets-suite) the result is capped at majorCopyNumber (a variant present on at least one
    /// copy cannot exceed the major-allele copy number).
    /// </summary>
    /// <param name="vaf">Observed variant allele fraction ∈ [0, 1].</param>
    /// <param name="purity">Tumour purity ρ ∈ (0, 1].</param>
    /// <param name="totalCopyNumber">Local tumour total copy number N_T (≥ 1).</param>
    /// <param name="majorCopyNumber">Local major-allele copy number, the upper bound on multiplicity (in [1, N_T]).</param>
    /// <returns>The integer mutation multiplicity m ∈ [1, majorCopyNumber].</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0,1], purity ∉ (0,1], totalCopyNumber &lt; 1, or majorCopyNumber ∉ [1, totalCopyNumber].</exception>
    public static int DeriveMultiplicity(double vaf, double purity, int totalCopyNumber, int majorCopyNumber)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        }

        if (totalCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(totalCopyNumber), totalCopyNumber, "Total copy number must be ≥ 1.");
        }

        if (majorCopyNumber < 1 || majorCopyNumber > totalCopyNumber)
        {
            throw new ArgumentOutOfRangeException(
                nameof(majorCopyNumber), majorCopyNumber, $"Major copy number must be in [1, {totalCopyNumber}].");
        }

        // n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)] — McGranahan 2016 observed mutation copy number (CCF=1 ⇒ m = n_mut).
        // This is exactly the canonical CNAqc purity/copy-number VAF correction (AdjustVAFForPurity, ONCO-VAF-001).
        double rawMultiplicity = AdjustVAFForPurity(vaf, purity, totalCopyNumber);
        // facets-suite expected_mutant_copies: R round() is IEC 60559 half-to-even (2.5 → 2), and mu < 1 → 1.
        int rounded = (int)Math.Round(rawMultiplicity, MidpointRounding.ToEven);
        // Clamp to [1, major CN]: the floor at 1 is facets-suite's mu < 1 → 1 (round(mu) ≤ 1 there); the cap at the
        // major-allele copy number is a documented extra (facets-suite leaves expected_alt_copies uncapped).
        return Math.Clamp(rounded, 1, majorCopyNumber);
    }

    /// <summary>
    /// Default ASPCF penalty used when a caller does not supply one: ASCAT's <c>ascat.aspcf(..., penalty = 70)</c>
    /// (VanLoo-lab/ascat, ascat.aspcf.R; Ross et al. 2021, <i>Bioinformatics</i> 37:1909). The penalty is charged per
    /// breakpoint on the <b>standardised</b> joint cost (each track's squared error divided by its MAD-based variance),
    /// so it is scale-free. The copynumber package's single-track default is γ = 40 (Nilsen et al. 2012).
    /// </summary>
    public const double AspcfDefaultPenalty = 70.0;

    /// <summary>ASCAT/copynumber minimum ASPCF segment length <c>kmin = 6</c> (<c>fastAspcf(..., 6, ...)</c>).</summary>
    private const int AspcfMinSegmentLength = 6;

    /// <summary>ASCAT <c>madWins(x, 2.5, 25)</c>: winsorisation at τ = 2.5 MAD-standard deviations.</summary>
    private const double AspcfWinsorTau = 2.5;

    /// <summary>ASCAT <c>madWins</c>/<c>getMad</c> running-median half-window k = 25 (filter width 2k + 1 = 51).</summary>
    private const int AspcfMedianHalfWindow = 25;

    /// <summary>ASCAT <c>fastAspcf</c> window size <c>w = 1000</c> and overlap <c>d = 100</c>.</summary>
    private const int AspcfWindowSize = 1000;

    /// <summary>ASCAT <c>fastAspcf</c> window overlap <c>d = 100</c>.</summary>
    private const int AspcfWindowOverlap = 100;

    /// <summary>ASCAT: re-run segmentation with the next larger penalty while ≥ 800 distinct logR levels remain.</summary>
    private const int AspcfMaxSegmentLevels = 800;

    /// <summary>ASCAT penalty ladder <c>segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140))</c>, kept where ≥ penalty.</summary>
    private static readonly double[] AspcfPenaltyLadder = { 35.0, 50.0, 70.0, 100.0, 140.0 };

    /// <summary>
    /// Allele-Specific Piecewise Constant Fitting (ASPCF) as run by ASCAT (<c>ascat.aspcf</c>, VanLoo-lab/ascat
    /// ascat.aspcf.R; the bivariate PCF of the copynumber package, Nilsen et al. 2012, <i>BMC Genomics</i> 13:591;
    /// Ross et al. 2021, <i>Bioinformatics</i> 37:1909), ported per chromosome on the heterozygous loci:
    /// <list type="number">
    /// <item>logR and mirrored BAF are MAD-winsorised (<c>madWins(x, 2.5, 25)</c>: running median of width 51 with
    /// Tukey end rule, residuals clipped at ±2.5·MAD).</item>
    /// <item>Chromosomes with fewer than 6 loci form one segment (mean winsorised mirrored BAF).</item>
    /// <item>Otherwise <c>fastAspcf</c>: overlapping windows (1000 loci, overlap 100); in each window the joint
    /// penalised cost <c>Σ_segments [SSE_logR/sd₁² + SSE_BAF/sd₂²] + γ·(#breakpoints)</c> is minimised exactly by the
    /// PCF dynamic program with minimum segment length kmin = 6 (<c>aspcfpart</c>), sd₁ / sd₂ being the MAD of the
    /// residuals from a running median (<c>getMad</c>) of logR and of the flipped BAF min(b, 1 − b); a window whose
    /// sd is 0 or undefined (e.g. noise-free data) contributes no breakpoint.</item>
    /// <item>Each segment's BAF is 0.5 + mean|b − 0.5|, shrunk to exactly 0.5 when
    /// <c>sqrt(sd₂² + μ²) &lt; 2·sd₂</c> (balanced); its logR is the mean of the raw (unwinsorised) logR.</item>
    /// <item>While ≥ 800 distinct segment logR levels remain, segmentation is repeated with the next larger penalty
    /// of the ASCAT ladder (35, 50, 70, 100, 140).</item>
    /// </list>
    /// Breakpoints never cross a contig boundary (loci are grouped into contiguous same-chromosome runs, input order).
    /// Every supplied locus is treated as a germline-heterozygous SNP; for mixed heterozygous / homozygous input (logR
    /// averaging over homozygous probes, homozygous-stretch resegmentation) use the germline-aware overload
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/> (B24 F36).
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome. LogR must be
    /// finite and BAF in [0, 1].</param>
    /// <param name="penalty">ASCAT/copynumber penalty γ &gt; 0 per breakpoint on the standardised cost
    /// (see <see cref="AspcfDefaultPenalty"/>).</param>
    /// <returns>The segment summaries (mean logR, ASPCF mirrored BAF) in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException">a locus has a null chromosome label, a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static IReadOnlyList<AlleleSpecificSegmentSummary> SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(loci);
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        var ordered = new List<AlleleSpecificLocus>();
        foreach (AlleleSpecificLocus locus in loci)
        {
            if (locus.Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            if (!IsValidAlleleSignal(locus.LogR, locus.BAF))
            {
                throw new ArgumentException("Every locus needs a finite logR and a BAF in [0, 1].", nameof(loci));
            }

            ordered.Add(locus);
        }

        // Contiguous same-chromosome runs (ASCATobj$chr), in input order.
        var runs = new List<(int Lo, int Hi)>();
        int start = 0;
        while (start < ordered.Count)
        {
            int end = start;
            while (end + 1 < ordered.Count && ordered[end + 1].Chromosome == ordered[start].Chromosome)
            {
                end++;
            }

            runs.Add((start, end));
            start = end + 1;
        }

        // segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140)); segmentlengths[segmentlengths >= penalty]
        var ladder = new List<double> { penalty };
        foreach (double p in AspcfPenaltyLadder)
        {
            if (p > penalty)
            {
                ladder.Add(p);
            }
        }

        List<AlleleSpecificSegmentSummary> result = new();
        foreach (double segmentPenalty in ladder)
        {
            result = new List<AlleleSpecificSegmentSummary>();
            foreach ((int lo, int hi) in runs)
            {
                SegmentChromosomeAspcf(ordered, lo, hi, segmentPenalty, result);
            }

            var levels = new HashSet<double>();
            foreach (AlleleSpecificSegmentSummary s in result)
            {
                levels.Add(s.MeanLogR);
            }

            if (levels.Count < AspcfMaxSegmentLevels)
            {
                break;
            }
        }

        return result;
    }

    /// <summary>
    /// ASCAT per-chromosome ASPCF (the loop body of <c>ascat.aspcf</c>) on one run of loci [lo, hi]; appends the
    /// segment summaries to <paramref name="output"/>.
    /// </summary>
    private static void SegmentChromosomeAspcf(
        List<AlleleSpecificLocus> loci, int lo, int hi, double penalty,
        List<AlleleSpecificSegmentSummary> output)
    {
        int n = hi - lo + 1;
        var lr = new double[n];
        var baf = new double[n];
        var mirrored = new double[n];
        for (int i = 0; i < n; i++)
        {
            lr[i] = loci[lo + i].LogR;
            baf[i] = loci[lo + i].BAF;
            mirrored[i] = MirrorBaf(baf[i]); // ifelse(bafsel > 0.5, bafsel, 1 - bafsel)
        }

        double[] lrWins = MadWinsorize(lr, AspcfWinsorTau, AspcfMedianHalfWindow);
        double[] bafWinsMirrored = MadWinsorize(mirrored, AspcfWinsorTau, AspcfMedianHalfWindow);
        var bafWins = new double[n];
        for (int i = 0; i < n; i++)
        {
            bafWins[i] = baf[i] > BalancedBaf ? bafWinsMirrored[i] : 1.0 - bafWinsMirrored[i];
        }

        if (n < AspcfMinSegmentLength)
        {
            // logRASPCF = mean(logRaveraged); bafASPCF = mean(bafselwinsmirrored); level re-adapted to mean(lr).
            output.Add(BuildAspcfSummary(loci, lo, 0, n, lr, Mean(bafWinsMirrored, 0, n)));
            return;
        }

        (int[] breakpoints, double[] segmentBaf) = FastAspcf(lrWins, bafWins, AspcfMinSegmentLength, penalty);
        for (int s = 0; s + 1 < breakpoints.Length; s++)
        {
            // Segment s covers the 0-based loci [breakpoints[s], breakpoints[s+1]); logR level = mean raw logR.
            output.Add(BuildAspcfSummary(loci, lo, breakpoints[s], breakpoints[s + 1], lr, segmentBaf[s]));
        }
    }

    /// <summary>Summary of loci [from, to) of the run starting at <paramref name="lo"/>: mean raw logR and the ASPCF BAF.</summary>
    private static AlleleSpecificSegmentSummary BuildAspcfSummary(
        List<AlleleSpecificLocus> loci, int lo, int from, int to, double[] rawLogR, double segmentBaf) =>
        new(
            Chromosome: loci[lo + from].Chromosome,
            Start: loci[lo + from].Position,
            End: loci[lo + to - 1].Position,
            MeanLogR: Mean(rawLogR, from, to),
            MeanBAF: segmentBaf,
            LocusCount: to - from);

    /// <summary>Arithmetic mean of x[from..to) (R <c>mean</c>).</summary>
    private static double Mean(double[] x, int from, int to)
    {
        double sum = 0.0;
        for (int i = from; i < to; i++)
        {
            sum += x[i];
        }

        return sum / (to - from);
    }

    /// <summary>
    /// ASCAT <c>fastAspcf(logR, allB, kmin, gamma)</c>: windowed exact bivariate PCF. Returns the segment boundaries
    /// (0-based, <c>[b₀ = 0, …, b_S = N]</c>) and each segment's BAF level <c>0.5 + μ</c> (μ shrunk to 0 when
    /// <c>sqrt(sd₂² + μ²) &lt; 2·sd₂</c>).
    /// </summary>
    private static (int[] Breakpoints, double[] SegmentBaf) FastAspcf(double[] logR, double[] allB, int kmin, double gamma)
    {
        int bigN = logR.Length;
        int startw = -AspcfWindowOverlap;
        int stopw = AspcfWindowSize - AspcfWindowOverlap;
        int validWindows = 0;
        double var2 = 0.0;
        var breakpts = new List<int> { 0 };
        while (true)
        {
            int from = Math.Max(1, startw);
            int to = Math.Min(stopw, bigN);
            int len = to - from + 1;
            var logRPart = new double[len];
            var allBFlip = new double[len];
            for (int i = 0; i < len; i++)
            {
                logRPart[i] = logR[from - 1 + i];
                double b = allB[from - 1 + i];
                allBFlip[i] = ToAscatBaf(b); // min(b, 1 − b)
            }

            double sd1 = GetMad(logRPart, AspcfMedianHalfWindow);
            double sd2 = GetMad(allBFlip, AspcfMedianHalfWindow);
            if (!double.IsNaN(sd1) && !double.IsNaN(sd2) && sd1 != 0.0 && sd2 != 0.0)
            {
                List<int> part = AspcfPart(logRPart, allBFlip, startw, stopw, AspcfWindowOverlap, sd1, sd2, bigN, kmin, gamma);
                int last = breakpts[^1];
                foreach (int bp in part)
                {
                    if (bp > last)
                    {
                        breakpts.Add(bp); // breakptspart[larger]
                    }
                }

                var2 += sd2 * sd2;
                validWindows++;
            }

            if (stopw < bigN + AspcfWindowOverlap)
            {
                startw = Math.Min(stopw - 2 * AspcfWindowOverlap + 1, bigN - 2 * AspcfWindowOverlap);
                stopw = startw + AspcfWindowSize;
            }
            else
            {
                break;
            }
        }

        // breakpts <- unique(c(breakpts, N))
        var unique = new List<int>();
        foreach (int bp in breakpts.Append(bigN))
        {
            if (!unique.Contains(bp))
            {
                unique.Add(bp);
            }
        }

        if (validWindows == 0)
        {
            validWindows = 1; // "just in case the sd-test never passes"
        }

        double sd2Pooled = Math.Sqrt(var2 / validWindows);
        var segmentBaf = new double[unique.Count - 1];
        for (int s = 0; s + 1 < unique.Count; s++)
        {
            int first = unique[s], lastExclusive = unique[s + 1];
            double mu = 0.0;
            for (int i = first; i < lastExclusive; i++)
            {
                mu += Math.Abs(allB[i] - BalancedBaf);
            }

            mu = lastExclusive > first ? mu / (lastExclusive - first) : 0.0;
            if (Math.Sqrt(sd2Pooled * sd2Pooled + mu * mu) < 2.0 * sd2Pooled)
            {
                mu = 0.0;
            }

            segmentBaf[s] = mu + BalancedBaf;
        }

        return (unique.ToArray(), segmentBaf);
    }

    /// <summary>
    /// ASCAT/copynumber <c>aspcfpart</c>: exact bivariate PCF (minimum segment length <paramref name="kmin"/>) on one
    /// window, cost = SSE₁/sd₁² + SSE₂/sd₂² per segment plus γ per breakpoint. Returns the global breakpoints (1-based
    /// "last index of a segment", as in R) that fall inside the window's use-range [a + d, b − d].
    /// </summary>
    private static List<int> AspcfPart(
        double[] y1, double[] y2, int a, int b, int d, double sd1, double sd2, int bigN, int kmin, double gamma)
    {
        int from = Math.Max(1, a);
        int useFrom = Math.Max(1, a + d);
        int useTo = Math.Min(bigN, b - d);
        int n = y1.Length;
        var result = new List<int>();
        if (n < 2 * kmin)
        {
            // R returns breakpts <- 0 unshifted ("Check that vectors are long enough to run algorithm").
            result.Add(0);
            return result;
        }

        double s1Sq = sd1 * sd1, s2Sq = sd2 * sd2;

        // 1-based arrays of length n + 1 (index 0 unused), mirroring the R code.
        double initSum1 = 0.0, initKvad1 = 0.0, initSum2 = 0.0, initKvad2 = 0.0;
        for (int i = 1; i <= kmin; i++)
        {
            initSum1 += y1[i - 1];
            initKvad1 += y1[i - 1] * y1[i - 1];
            initSum2 += y2[i - 1];
            initKvad2 += y2[i - 1] * y2[i - 1];
        }

        double initAve1 = initSum1 / kmin, initAve2 = initSum2 / kmin;
        var bestCost = new double[n + 1];
        var bestSplit = new int[n + 1];
        bestCost[kmin] = (initKvad1 - initSum1 * initAve1) / s1Sq + (initKvad2 - initSum2 * initAve2) / s2Sq;

        var sum1 = new double[n + 1];
        var sum2 = new double[n + 1];
        var kvad1 = new double[n + 1];
        var kvad2 = new double[n + 1];
        var aver1 = new double[n + 1];
        var aver2 = new double[n + 1];
        var cost = new double[n + 1];
        int kminP1 = kmin + 1;

        for (int k = kminP1; k <= 2 * kmin - 1; k++)
        {
            double yk1 = y1[k - 1], yk2 = y2[k - 1];
            for (int t = kminP1; t <= k; t++)
            {
                sum1[t] += yk1;
                aver1[t] = sum1[t] / (k - t + 1);
                kvad1[t] += yk1 * yk1;
                sum2[t] += yk2;
                aver2[t] = sum2[t] / (k - t + 1);
                kvad2[t] += yk2 * yk2;
            }

            double bestAver1 = (initSum1 + sum1[kminP1]) / k;
            double bestAver2 = (initSum2 + sum2[kminP1]) / k;
            double cost1 = ((initKvad1 + kvad1[kminP1]) - k * bestAver1 * bestAver1) / s1Sq;
            double cost2 = ((initKvad2 + kvad2[kminP1]) - k * bestAver2 * bestAver2) / s2Sq;
            bestCost[k] = cost1 + cost2;
        }

        for (int m = 2 * kmin; m <= n; m++)
        {
            int nMkminP1 = m - kmin + 1;
            double ym1 = y1[m - 1], ym2 = y2[m - 1];
            for (int t = kminP1; t <= m; t++)
            {
                sum1[t] += ym1;
                aver1[t] = sum1[t] / (m - t + 1);
                kvad1[t] += ym1 * ym1;
                sum2[t] += ym2;
                aver2[t] = sum2[t] / (m - t + 1);
                kvad2[t] += ym2 * ym2;
            }

            // Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2; Pos <- which.min(...) + kmin
            int pos = -1;
            double minCost = double.PositiveInfinity;
            for (int t = kminP1; t <= nMkminP1; t++)
            {
                double c1 = (kvad1[t] - sum1[t] * aver1[t]) / s1Sq;
                double c2 = (kvad2[t] - sum2[t] * aver2[t]) / s2Sq;
                cost[t] = bestCost[t - 1] + c1 + c2;
                if (pos < 0 || cost[t] < minCost)
                {
                    minCost = cost[t];
                    pos = t;
                }
            }

            double best = cost[pos] + gamma;
            double totAver1 = (sum1[kminP1] + initSum1) / m;
            double totCost1 = ((kvad1[kminP1] + initKvad1) - m * totAver1 * totAver1) / s1Sq;
            double totAver2 = (sum2[kminP1] + initSum2) / m;
            double totCost2 = ((kvad2[kminP1] + initKvad2) - m * totAver2 * totAver2) / s2Sq;
            double totCost = totCost1 + totCost2;
            if (totCost < best)
            {
                pos = 1;
                best = totCost;
            }

            bestCost[m] = best;
            bestSplit[m] = pos - 1;
        }

        // Trace back: breakpts <- c(bestSplit[n], breakpts) until n == 0 (the leading 0 is included).
        var trace = new List<int> { n };
        int cursor = n;
        while (cursor > 0)
        {
            cursor = bestSplit[cursor];
            trace.Insert(0, cursor);
        }

        foreach (int bp in trace)
        {
            int global = bp + from - 1;
            if (global >= useFrom && global <= useTo)
            {
                result.Add(global);
            }
        }

        return result;
    }

    /// <summary>
    /// ASCAT <c>getMad(x, k)</c>: zeros removed (likely imputed), then <c>mad(x − medianFilter(x, k))</c>. NaN when no
    /// value remains (R <c>mad(numeric(0))</c> = NA).
    /// </summary>
    private static double GetMad(double[] x, int k)
    {
        var nonZero = new List<double>(x.Length);
        foreach (double v in x)
        {
            if (v != 0.0)
            {
                nonZero.Add(v);
            }
        }

        if (nonZero.Count == 0)
        {
            return double.NaN;
        }

        double[] values = nonZero.ToArray();
        double[] runMedian = MedianFilter(values, k);
        var dif = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            dif[i] = values[i] - runMedian[i];
        }

        return RMad(dif);
    }

    /// <summary>ASCAT <c>madWins(x, tau, k)</c>: <c>xhat = medianFilter(x, k)</c>; residuals clipped at ±tau·mad(x − xhat).</summary>
    private static double[] MadWinsorize(double[] x, double tau, int k)
    {
        if (x.Length == 0)
        {
            return Array.Empty<double>();
        }

        double[] xhat = MedianFilter(x, k);
        var d = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            d[i] = x[i] - xhat[i];
        }

        double z = tau * RMad(d);
        var win = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            win[i] = xhat[i] + Math.Clamp(d[i], -z, z); // psi(d, z)
        }

        return win;
    }

    /// <summary>
    /// R <c>mad(x)</c> = 1.4826 · median(|x − median(x)|) (constant <see cref="MadConsistencyConstant"/>, raw MAD from
    /// <see cref="RawMedianAbsoluteDeviation"/>, both shared with the MATH score).
    /// </summary>
    private static double RMad(double[] x) =>
        MadConsistencyConstant * RawMedianAbsoluteDeviation(x, StatisticsHelper.Median(x));

    /// <summary>
    /// Unscaled median absolute deviation about <paramref name="center"/>: median(|xᵢ − center|) (canonical
    /// <see cref="StatisticsHelper.Median"/>). Shared by R <c>mad</c> (<see cref="RMad"/>, ASCAT ASPCF) and the maftools
    /// MATH score (<see cref="CalculateITH"/>), which apply the 1.4826 factor in different operation orders.
    /// </summary>
    private static double RawMedianAbsoluteDeviation(double[] x, double center)
    {
        var dev = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            dev[i] = Math.Abs(x[i] - center);
        }

        return StatisticsHelper.Median(dev);
    }

    /// <summary>
    /// ASCAT <c>medianFilter(x, k)</c>: <c>runmed(x, 2k + 1, endrule = "median")</c>, the width reduced to n (odd n) or
    /// n − 1 (even n) when it exceeds the series length.
    /// </summary>
    private static double[] MedianFilter(double[] x, int k)
    {
        int n = x.Length;
        int width = 2 * k + 1;
        if (width > n)
        {
            if (n == 0)
            {
                width = 1;
            }
            else
            {
                width = n % 2 == 0 ? n - 1 : n;
            }
        }

        return RunningMedian(x, width);
    }

    /// <summary>
    /// R <c>runmed(x, k, endrule = "median")</c> for odd k ≤ n: centred running medians of width k in the interior,
    /// the first/last k%/%2 values kept, then <c>smoothEnds(res, k)</c> (Tukey's end-point rule).
    /// </summary>
    private static double[] RunningMedian(double[] x, int k)
    {
        int n = x.Length;
        var res = (double[])x.Clone();
        int half = k / 2;
        if (half < 1)
        {
            return res;
        }

        var window = new double[k];
        for (int i = half; i < n - half; i++)
        {
            Array.Copy(x, i - half, window, 0, k);
            res[i] = StatisticsHelper.Median(window);
        }

        return SmoothEnds(res, k);
    }

    /// <summary>R <c>smoothEnds(y, k)</c> (stats package), verbatim.</summary>
    private static double[] SmoothEnds(double[] y, int k)
    {
        int half = k / 2;
        if (half < 1)
        {
            return y;
        }

        int n = y.Length;
        var sm = (double[])y.Clone();
        // 0-based: R index i ↔ C# index i − 1.
        if (half >= 2)
        {
            sm[1] = Med3(y[0], y[1], y[2]);
            sm[n - 2] = Med3(y[n - 1], y[n - 2], y[n - 3]);
            for (int i = 3; i <= half; i++)
            {
                int j = 2 * i - 1;
                sm[i - 1] = MedianOfRange(y, 0, j);
                sm[n - i] = MedianOfRange(y, n - j, j);
            }
        }

        sm[0] = Med3(y[0], sm[1], sm[1] - 2.0 * (sm[2] - sm[1]));
        sm[n - 1] = Med3(y[n - 1], sm[n - 2], sm[n - 2] - 2.0 * (sm[n - 3] - sm[n - 2]));
        return sm;
    }

    /// <summary>Median of y[start .. start + count) (odd count, R <c>med.odd</c>).</summary>
    private static double MedianOfRange(double[] y, int start, int count)
    {
        var tmp = new double[count];
        Array.Copy(y, start, tmp, 0, count);
        return StatisticsHelper.Median(tmp);
    }

    /// <summary>R <c>smoothEnds</c> helper <c>med3(a, b, c)</c>, verbatim.</summary>
    private static double Med3(double a, double b, double c)
    {
        double m = b;
        if (a < b)
        {
            if (c < b)
            {
                m = a >= c ? a : c;
            }
        }
        else
        {
            if (c > b)
            {
                m = a <= c ? a : c;
            }
        }

        return m;
    }

    // ---- ASCAT ascat.aspcf with germline-homozygous probes (B24 F36) ----

    /// <summary>ASCAT <c>predictGermlineHomozygousStretches</c>: hard-coded p-value 0.001 of a homozygous run.</summary>
    private const double AscatHomozygousStretchPValue = 0.001;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: PCF context of ±100 probes around the stretch.</summary>
    private const int AscatHomStretchContext = 100;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: replacement margin of ±5 probes around the stretch.</summary>
    private const int AscatHomStretchMargin = 5;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: a probe is replaced when |Δ level| &gt; 0.3 …</summary>
    private const double AscatHomStretchLevelDifference = 0.3;

    /// <summary>… and only when more than 5 probes of the stretch differ (<c>sum(dif &gt; 0.3) &gt; 5</c>).</summary>
    private const int AscatHomStretchMinDifferingProbes = 5;

    /// <summary>
    /// One ASCAT logR segment of a germline-aware ASPCF segmentation (<see cref="AspcfSegmentation"/>): a maximal run of
    /// loci with the same segmented logR inside one chromosome — exactly the segments runASCAT builds from
    /// <c>rle(lrrsegmented)</c> ∪ chromosome ends (ascat.runAscat.R, <c>tlrstart</c>/<c>tlrend</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Start">Position of the first locus of the segment.</param>
    /// <param name="End">Position of the last locus of the segment.</param>
    /// <param name="MeanLogR">Segmented logR level (mean raw logR of the segment's loci, heterozygous and homozygous).</param>
    /// <param name="MeanBAF">Mirrored (≥ 0.5) segmented BAF of the segment's <b>first</b> germline-heterozygous locus
    /// (runASCAT <c>bafke = bafsegmented[bafpos][1]</c>); <see cref="double.NaN"/> when the segment has no heterozygous
    /// locus (a germline-homozygous stretch: ASCAT then sets <c>bafke = 0</c> and only nA + nB is meaningful).</param>
    /// <param name="LocusCount">Number of loci (all germline genotypes) in the segment.</param>
    /// <param name="HeterozygousLocusCount">Number of germline-heterozygous loci in the segment.</param>
    public readonly record struct AspcfSegment(
        string Chromosome,
        long Start,
        long End,
        double MeanLogR,
        double MeanBAF,
        int LocusCount,
        int HeterozygousLocusCount)
    {
        /// <summary>True when the segment has a segmented BAF (at least one germline-heterozygous locus).</summary>
        public bool HasBaf => !double.IsNaN(MeanBAF);
    }

    /// <summary>
    /// Result of ASCAT allele-specific segmentation with germline genotypes
    /// (<see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/>): the
    /// per-locus ASCAT <c>Tumor_LogR_segmented</c> (every locus) and <c>Tumor_BAF_segmented</c> (heterozygous loci only),
    /// plus the logR segments runASCAT builds from them. Consumed by <see cref="FitPurityPloidyFromAspcf"/>.
    /// </summary>
    public sealed class AspcfSegmentation
    {
        internal AspcfSegmentation(
            AlleleSpecificLocus[] loci, bool[] heterozygous, double[] segmentedLogR, double[] segmentedBaf,
            AspcfSegment[] segments)
        {
            Loci = Array.AsReadOnly(loci);
            GermlineHeterozygous = Array.AsReadOnly(heterozygous);
            SegmentedLogR = Array.AsReadOnly(segmentedLogR);
            SegmentedBaf = Array.AsReadOnly(segmentedBaf);
            Segments = Array.AsReadOnly(segments);
        }

        /// <summary>The input loci, in input order.</summary>
        public IReadOnlyList<AlleleSpecificLocus> Loci { get; }

        /// <summary>The germline genotype of each locus (true = heterozygous), in input order.</summary>
        public IReadOnlyList<bool> GermlineHeterozygous { get; }

        /// <summary>ASCAT <c>Tumor_LogR_segmented</c>: the segmented logR level of every locus.</summary>
        public IReadOnlyList<double> SegmentedLogR { get; }

        /// <summary>Mirrored (≥ 0.5) segmented BAF of each heterozygous locus (ASCAT <c>Tumor_BAF_segmented = 1 − value</c>);
        /// <see cref="double.NaN"/> at germline-homozygous loci (ASCAT segments BAF on heterozygous probes only).</summary>
        public IReadOnlyList<double> SegmentedBaf { get; }

        /// <summary>The logR segments (runs of equal segmented logR within a chromosome), in input order.</summary>
        public IReadOnlyList<AspcfSegment> Segments { get; }
    }

    /// <summary>
    /// ASCAT allele-specific segmentation with germline genotypes — a port of the complete <c>ascat.aspcf</c>
    /// (VanLoo-lab/ascat ascat.aspcf.R; Ross et al. 2021, <i>Bioinformatics</i> 37:1909) for mixed heterozygous /
    /// homozygous input, as ASCAT runs it with <c>ascat.gg$germlinegenotypes</c> (TRUE = homozygous):
    /// <list type="number">
    /// <item>Germline-homozygous stretches (<c>predictGermlineHomozygousStretches</c>): with h = the fraction of homozygous
    /// loci, a run of ≥ <c>ceiling(log(0.001, h))</c> consecutive homozygous loci is a stretch.</item>
    /// <item>Per chromosome: logR of <b>all</b> loci is MAD-winsorised; BAF of the heterozygous loci only. Each
    /// heterozygous locus gets the mean winsorised logR of the loci between the midpoints to its heterozygous
    /// neighbours; ASPCF (<c>fastAspcf</c>, kmin 6; &lt; 6 heterozygous loci ⇒ one segment) runs on these averages and
    /// the heterozygous BAF. Between consecutive heterozygous loci of different levels the logR breakpoint is placed
    /// in the homozygous gap at the minimum absolute deviation (ascat.aspcf "find best breakpoint", verbatim
    /// including its <c>1:0</c> indexing at bp = 0); levels are re-estimated as the mean raw logR. A chromosome without
    /// heterozygous loci becomes one logR segment (mean raw logR) with no BAF.</item>
    /// <item>Each homozygous stretch is re-segmented: exact PCF (<c>exactPcf</c>, kmin 6, penalty floor(γ/4)) of the
    /// winsorised raw logR over the stretch ± 100 loci; the stretch ± 5 loci takes the PCF level where it differs by
    /// more than 0.3, if more than 5 loci differ.</item>
    /// <item>Genome-wide: zero levels are filled from their neighbours (<c>fillNA(zeroIsNA = TRUE)</c>) and every run of
    /// equal levels is re-estimated as its mean raw logR; the penalty ladder (35, 50, 70, 100, 140) is climbed while
    /// ≥ 800 distinct levels remain.</item>
    /// </list>
    /// With every locus heterozygous the logR levels and BAF are those of
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, double)"/> (which is unchanged), apart from
    /// two R genome-wide steps that overload omits — a level exactly 0 is replaced by <c>fillNA</c>, and equal adjacent
    /// levels on consecutive chromosomes are re-averaged together.
    /// Loci of one chromosome must be contiguous in the input (each contiguous same-label run is one ASCAT <c>chr</c>
    /// part); positions are not re-sorted.
    /// </summary>
    /// <param name="loci">Per-locus measurements in genome order. LogR must be finite; BAF must be in [0, 1] at
    /// heterozygous loci and is ignored at homozygous loci (it may be NaN there, e.g. a copy-number-only probe).</param>
    /// <param name="germlineHeterozygous">Germline genotype per locus: true = heterozygous (ASCAT
    /// <c>germlinegenotypes == FALSE</c>), false = homozygous. Same length as <paramref name="loci"/>.</param>
    /// <param name="penalty">ASPCF penalty (ASCAT <c>penalty</c>, default <see cref="AspcfDefaultPenalty"/> = 70).</param>
    /// <returns>The per-locus segmented logR/BAF and the runASCAT logR segments.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Lengths differ, a chromosome label is null, a logR is non-finite, or a
    /// heterozygous locus has a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static AspcfSegmentation SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        IReadOnlyList<bool> germlineHeterozygous,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(loci);
        ArgumentNullException.ThrowIfNull(germlineHeterozygous);
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        AlleleSpecificLocus[] probes = loci.ToArray();
        if (probes.Length != germlineHeterozygous.Count)
        {
            throw new ArgumentException(
                "Exactly one germline genotype is required per locus.", nameof(germlineHeterozygous));
        }

        var het = new bool[probes.Length];
        for (int i = 0; i < probes.Length; i++)
        {
            het[i] = germlineHeterozygous[i];
            if (probes[i].Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            if (!double.IsFinite(probes[i].LogR))
            {
                throw new ArgumentException("Every locus needs a finite logR.", nameof(loci));
            }

            if (het[i] && !IsValidAlleleSignal(probes[i].LogR, probes[i].BAF))
            {
                throw new ArgumentException("Every heterozygous locus needs a BAF in [0, 1].", nameof(loci));
            }
        }

        int n = probes.Length;
        var logR = new double[n];
        for (int i = 0; i < n; i++)
        {
            logR[i] = probes[i].LogR;
        }

        List<(int Lo, int Hi)> runs = ContiguousChromosomeRuns(probes);
        List<(int Run, int Start, int End)> stretches = PredictGermlineHomozygousStretches(runs, het);

        var ladder = new List<double> { penalty };
        foreach (double p in AspcfPenaltyLadder)
        {
            if (p > penalty)
            {
                ladder.Add(p);
            }
        }

        double[] logRPcfed = Array.Empty<double>();
        double[] bafPcfed = new double[n];
        foreach (double segmentPenalty in ladder)
        {
            logRPcfed = new double[n];
            Array.Fill(bafPcfed, double.NaN);
            for (int r = 0; r < runs.Count; r++)
            {
                SegmentChromosomeAspcfGermline(logR, probes, het, runs[r].Lo, runs[r].Hi, segmentPenalty, logRPcfed, bafPcfed);
                ResegmentHomozygousStretches(logR, runs[r].Lo, runs[r].Hi, stretches, r, segmentPenalty, logRPcfed);
            }

            AscatFillNa(logRPcfed);
            ReadaptLevels(logRPcfed, logR, 0, n);

            if (logRPcfed.Distinct().Count() < AspcfMaxSegmentLevels)
            {
                break;
            }
        }

        return new AspcfSegmentation(probes, het, logRPcfed, bafPcfed, BuildAspcfSegments(probes, het, runs, logRPcfed, bafPcfed));
    }

    /// <summary>Contiguous same-chromosome runs of the loci (ASCAT <c>chr</c> parts), in input order.</summary>
    private static List<(int Lo, int Hi)> ContiguousChromosomeRuns(AlleleSpecificLocus[] loci)
    {
        var runs = new List<(int Lo, int Hi)>();
        int start = 0;
        while (start < loci.Length)
        {
            int end = start;
            while (end + 1 < loci.Length && loci[end + 1].Chromosome == loci[start].Chromosome)
            {
                end++;
            }

            runs.Add((start, end));
            start = end + 1;
        }

        return runs;
    }

    /// <summary>
    /// ASCAT <c>predictGermlineHomozygousStretches</c>: runs of ≥ <c>ceiling(log(0.001, h))</c> consecutive homozygous loci
    /// per chromosome, h = genome-wide homozygous fraction (R <c>log(x, base) = log(x)/log(base)</c>). Returns 0-based
    /// global (first, last) indices.
    /// </summary>
    private static List<(int Run, int Start, int End)> PredictGermlineHomozygousStretches(
        List<(int Lo, int Hi)> runs, bool[] het)
    {
        var stretches = new List<(int Run, int Start, int End)>();
        if (het.Length == 0)
        {
            return stretches;
        }

        int homCount = 0;
        foreach (bool h in het)
        {
            if (!h)
            {
                homCount++;
            }
        }

        double fractionHom = (double)homCount / het.Length;
        double threshold = Math.Ceiling(Math.Log(AscatHomozygousStretchPValue) / Math.Log(fractionHom));
        for (int r = 0; r < runs.Count; r++)
        {
            int count = 0, first = -1;
            for (int i = runs[r].Lo; i <= runs[r].Hi; i++)
            {
                if (!het[i])
                {
                    if (count == 0)
                    {
                        first = i;
                    }

                    count++;
                }
                else
                {
                    // An empty run gives R an NA row (min of an empty vector), which ascat.aspcf skips.
                    if (count > 0 && count >= threshold)
                    {
                        stretches.Add((r, first, i - 1));
                    }

                    count = 0;
                }
            }

            if (count > 0 && count >= threshold)
            {
                stretches.Add((r, first, runs[r].Hi));
            }
        }

        return stretches;
    }

    /// <summary>
    /// The per-chromosome body of <c>ascat.aspcf</c> with germline genotypes, on loci [lo, hi]: writes the segmented logR
    /// of every locus to <paramref name="logRPcfed"/> and the mirrored segmented BAF of each heterozygous locus to
    /// <paramref name="bafPcfed"/>.
    /// </summary>
    private static void SegmentChromosomeAspcfGermline(
        double[] logR, AlleleSpecificLocus[] loci, bool[] het, int lo, int hi, double penalty,
        double[] logRPcfed, double[] bafPcfed)
    {
        int len = hi - lo + 1;
        var lr = new double[len];
        Array.Copy(logR, lo, lr, 0, len);
        double[] lrWins = MadWinsorize(lr, AspcfWinsorTau, AspcfMedianHalfWindow);

        var indices = new List<int>(); // 1-based local indices of the heterozygous loci (R "indices")
        for (int i = 0; i < len; i++)
        {
            if (het[lo + i])
            {
                indices.Add(i + 1);
            }
        }

        int h = indices.Count;
        if (h == 0)
        {
            // No heterozygous probe: a single logR segment, no BAF.
            double level = Mean(lr, 0, len);
            Array.Fill(logRPcfed, level, lo, len);
            return;
        }

        var bafSel = new double[h];
        var mirrored = new double[h];
        for (int k = 0; k < h; k++)
        {
            bafSel[k] = loci[lo + indices[k] - 1].BAF;
            mirrored[k] = MirrorBaf(bafSel[k]);
        }

        double[] bafWinsMirrored = MadWinsorize(mirrored, AspcfWinsorTau, AspcfMedianHalfWindow);
        var bafWins = new double[h];
        for (int k = 0; k < h; k++)
        {
            bafWins[k] = bafSel[k] > BalancedBaf ? bafWinsMirrored[k] : 1.0 - bafWinsMirrored[k];
        }

        // averageIndices = c(1, (indices[-h] + indices[-1])/2, length(lr) + 0.01); start = ceiling, end = floor(· − 0.01).
        var logRAveraged = new double[h];
        for (int k = 0; k < h; k++)
        {
            int startIndex, endIndex;
            if (h == 1)
            {
                startIndex = 1;
                endIndex = len;
            }
            else
            {
                double lower = k == 0 ? 1.0 : (indices[k - 1] + indices[k]) / 2.0;
                double upper = k == h - 1 ? len + 0.01 : (indices[k] + indices[k + 1]) / 2.0;
                startIndex = (int)Math.Ceiling(lower);
                endIndex = (int)Math.Floor(upper - 0.01);
            }

            logRAveraged[k] = Mean(lrWins, startIndex - 1, endIndex);
        }

        var levelPerHet = new double[h];
        var bafPerHet = new double[h];
        if (h < AspcfMinSegmentLength)
        {
            Array.Fill(levelPerHet, Mean(logRAveraged, 0, h));
            Array.Fill(bafPerHet, Mean(bafWinsMirrored, 0, h));
        }
        else
        {
            (int[] breakpoints, double[] segmentBaf) = FastAspcf(logRAveraged, bafWins, AspcfMinSegmentLength, penalty);
            for (int s = 0; s + 1 < breakpoints.Length; s++)
            {
                double level = Mean(logRAveraged, breakpoints[s], breakpoints[s + 1]); // fastAspcf yhat1
                for (int k = breakpoints[s]; k < breakpoints[s + 1]; k++)
                {
                    levelPerHet[k] = level;
                    bafPerHet[k] = segmentBaf[s];
                }
            }
        }

        // logRc: extend the heterozygous levels over the homozygous gaps ("find best breakpoint").
        var logRc = new double[len];
        for (int p = 0; p < h; p++)
        {
            int at = indices[p]; // 1-based
            if (p == 0)
            {
                Array.Fill(logRc, levelPerHet[0], 0, at);
            }

            if (p == h - 1)
            {
                Array.Fill(logRc, levelPerHet[p], at, len - at);
            }
            else if (levelPerHet[p] == levelPerHet[p + 1])
            {
                Array.Fill(logRc, levelPerHet[p], at, indices[p + 1] - at);
            }
            else
            {
                int total = indices[p + 1] - at;
                int breakpoint = AscatBestGapBreakpoint(lr, at, total, levelPerHet[p], levelPerHet[p + 1]);
                Array.Fill(logRc, levelPerHet[p], at, breakpoint);
                Array.Fill(logRc, levelPerHet[p + 1], at + breakpoint, total - breakpoint);
            }
        }

        // 2nd step: adapt levels (rle(logRc) runs → mean raw logR).
        int runStart = 0;
        for (int i = 1; i <= len; i++)
        {
            if (i == len || logRc[i] != logRc[runStart])
            {
                double level = Mean(lr, runStart, i);
                Array.Fill(logRPcfed, level, lo + runStart, i - runStart);
                runStart = i;
            }
        }

        for (int k = 0; k < h; k++)
        {
            bafPcfed[lo + indices[k] - 1] = bafPerHet[k];
        }
    }

    /// <summary>
    /// ascat.aspcf "find best breakpoint" between heterozygous loci at 1-based local positions <paramref name="at"/> and
    /// at + <paramref name="total"/>: for bp ∈ 0..total−1, d(bp) = Σ|lr[(1:bp) + at] − left| + Σ|lr[((bp+1):total) + at] −
    /// right|; returns <c>which.min(d) − 1</c>. R's <c>1:0 = c(1, 0)</c> makes bp = 0 compare lr[at + 1] and lr[at] with
    /// the left level — kept verbatim.
    /// </summary>
    private static int AscatBestGapBreakpoint(double[] lr, int at, int total, double left, double right)
    {
        int best = 0;
        double bestDistance = double.NaN;
        for (int bp = 0; bp < total; bp++)
        {
            double leftSum = 0.0;
            if (bp == 0)
            {
                leftSum += Math.Abs(lr[at] - left);     // lr[1 + at] (R 1-based) = C# lr[at]
                leftSum += Math.Abs(lr[at - 1] - left); // lr[0 + at]
            }
            else
            {
                for (int q = 1; q <= bp; q++)
                {
                    leftSum += Math.Abs(lr[q + at - 1] - left);
                }
            }

            double rightSum = 0.0;
            for (int q = bp + 1; q <= total; q++)
            {
                rightSum += Math.Abs(lr[q + at - 1] - right);
            }

            double distance = leftSum + rightSum;
            if (bp == 0 || distance < bestDistance)
            {
                bestDistance = distance;
                best = bp;
            }
        }

        return best;
    }

    /// <summary>
    /// ascat.aspcf "correct wrong segments in germline homozygous stretches" for the stretches of chromosome run
    /// <paramref name="run"/> [lo, hi]: exact PCF (kmin 6, penalty floor(γ/4)) of the winsorised raw logR over the stretch
    /// ± 100 loci; the stretch ± 5 loci takes the PCF level where it differs by &gt; 0.3, if &gt; 5 loci differ.
    /// </summary>
    private static void ResegmentHomozygousStretches(
        double[] logR, int lo, int hi, List<(int Run, int Start, int End)> stretches, int run, double penalty,
        double[] logRPcfed)
    {
        double pcfPenalty = Math.Floor(penalty / 4.0);
        foreach ((int stretchRun, int start, int end) in stretches)
        {
            if (stretchRun != run)
            {
                continue;
            }

            int start2 = Math.Max(start - AscatHomStretchContext, lo);
            int end2 = Math.Min(end + AscatHomStretchContext, hi);
            int start3 = Math.Max(start - AscatHomStretchMargin, lo);
            int end3 = Math.Min(end + AscatHomStretchMargin, hi);

            var window = new double[end2 - start2 + 1];
            Array.Copy(logR, start2, window, 0, window.Length);
            double[] pcfed = ExactPcf(MadWinsorize(window, AspcfWinsorTau, AspcfMedianHalfWindow), AspcfMinSegmentLength, pcfPenalty);

            int count = end3 - start3 + 1;
            var dif = new double[count];
            int differing = 0;
            for (int i = 0; i < count; i++)
            {
                dif[i] = Math.Abs(pcfed[start3 - start2 + i] - logRPcfed[start3 + i]);
                if (dif[i] > AscatHomStretchLevelDifference)
                {
                    differing++;
                }
            }

            if (differing > AscatHomStretchMinDifferingProbes)
            {
                for (int i = 0; i < count; i++)
                {
                    if (dif[i] > AscatHomStretchLevelDifference)
                    {
                        logRPcfed[start3 + i] = pcfed[start3 - start2 + i];
                    }
                }
            }
        }
    }

    /// <summary>
    /// ASCAT <c>exactPcf(y, kmin, gamma)</c> (ascat.aspcf.R): exact univariate piecewise-constant fit by Potts filtering,
    /// minimum segment length kmin, penalty γ per discontinuity on the raw (unstandardised) squared error.
    /// </summary>
    private static double[] ExactPcf(double[] y, int kmin, double gamma)
    {
        int n = y.Length;
        var yhat = new double[n];
        if (n < 2 * kmin)
        {
            Array.Fill(yhat, Mean(y, 0, n));
            return yhat;
        }

        // 1-based arrays (index 0 unused), mirroring the R code.
        double initSum = 0.0, initKvad = 0.0;
        for (int i = 1; i <= kmin; i++)
        {
            initSum += y[i - 1];
            initKvad += y[i - 1] * y[i - 1];
        }

        double initAve = initSum / kmin;
        var bestCost = new double[n + 1];
        var bestSplit = new int[n + 1];
        var bestAver = new double[n + 1];
        var sum = new double[n + 1];
        var kvad = new double[n + 1];
        var aver = new double[n + 1];
        var cost = new double[n + 1];
        bestCost[kmin] = initKvad - initSum * initAve;
        bestAver[kmin] = initAve;
        int kminP1 = kmin + 1;
        for (int k = kminP1; k <= 2 * kmin - 1; k++)
        {
            double yk = y[k - 1];
            for (int t = kminP1; t <= k; t++)
            {
                sum[t] += yk;
                aver[t] = sum[t] / (k - t + 1);
                kvad[t] += yk * yk;
            }

            bestAver[k] = (initSum + sum[kminP1]) / k;
            bestCost[k] = (initKvad + kvad[kminP1]) - k * (bestAver[k] * bestAver[k]);
        }

        for (int m = 2 * kmin; m <= n; m++)
        {
            double yn = y[m - 1];
            double yn2 = yn * yn;
            for (int t = kminP1; t <= m; t++)
            {
                sum[t] += yn;
                aver[t] = sum[t] / (m - t + 1);
                kvad[t] += yn2;
            }

            int nMkminP1 = m - kmin + 1;
            int pos = -1;
            for (int t = kminP1; t <= nMkminP1; t++)
            {
                cost[t] = bestCost[t - 1] + kvad[t] - sum[t] * aver[t] + gamma;
                if (pos < 0 || cost[t] < cost[pos])
                {
                    pos = t; // which.min: first minimum
                }
            }

            double best = cost[pos];
            double bestLevel = aver[pos];
            double totAver = (sum[kminP1] + initSum) / m;
            double totCost = (kvad[kminP1] + initKvad) - m * totAver * totAver;
            if (totCost < best)
            {
                pos = 1;
                best = totCost;
                bestLevel = totAver;
            }

            bestCost[m] = best;
            bestAver[m] = bestLevel;
            bestSplit[m] = pos - 1;
        }

        int cursor = n;
        while (cursor > 0)
        {
            for (int i = bestSplit[cursor]; i < cursor; i++)
            {
                yhat[i] = bestAver[cursor];
            }

            cursor = bestSplit[cursor];
        }

        return yhat;
    }

    /// <summary>
    /// ASCAT <c>fillNA(vec, zeroIsNA = TRUE)</c> (ascat.aspcf.R), verbatim: levels exactly 0 become missing; a missing run
    /// at the start/end takes the next/previous value; an interior run takes the previous value up to
    /// <c>start + ceiling(N/2)</c> and the next value after it (for N = 2, 3 R's descending <c>(end+1):end</c> sets only
    /// the last element to the next value).
    /// </summary>
    private static void AscatFillNa(double[] vec)
    {
        int n = vec.Length;
        var missing = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (vec[i] == 0.0 || double.IsNaN(vec[i]))
            {
                vec[i] = double.NaN;
                missing.Add(i);
            }
        }

        if (missing.Count == 0)
        {
            return;
        }

        var starts = new List<int> { missing[0] };
        var ends = new List<int>();
        for (int k = 1; k < missing.Count; k++)
        {
            if (missing[k] - missing[k - 1] > 1)
            {
                ends.Add(missing[k - 1]);
                starts.Add(missing[k]);
            }
        }

        ends.Add(missing[^1]);
        double At(int index) => index < n ? vec[index] : double.NaN; // R: out-of-range index → NA

        int startAt = 0;
        if (starts[0] == 0)
        {
            double fill = At(ends[0] + 1);
            Array.Fill(vec, fill, 0, ends[0] + 1);
            startAt = 1;
        }

        if (startAt >= starts.Count)
        {
            return;
        }

        int endAt = starts.Count - 1;
        if (double.IsNaN(vec[n - 1]))
        {
            Array.Fill(vec, vec[starts[endAt] - 1], starts[endAt], ends[endAt] - starts[endAt] + 1);
            endAt--;
        }

        for (int k = startAt; k <= endAt; k++)
        {
            int start = starts[k], end = ends[k];
            int count = 1 + end - start;
            if (count == 1)
            {
                vec[start] = vec[start - 1];
                continue;
            }

            int midpoint = start + (count + 1) / 2; // start + ceiling(N/2)
            double previous = vec[start - 1], next = vec[end + 1];
            Array.Fill(vec, previous, start, midpoint - start + 1);
            if (midpoint < end)
            {
                Array.Fill(vec, next, midpoint + 1, end - midpoint);
            }
            else
            {
                vec[end] = next; // vec[(end+1):end] = vec[end+1]
            }
        }
    }

    /// <summary>
    /// ascat.aspcf "adapt levels again": every run of equal values of <paramref name="levels"/> in [from, to) — across
    /// chromosome boundaries, as R's genome-wide <c>rle</c> — is set to the mean raw logR of the run.
    /// </summary>
    private static void ReadaptLevels(double[] levels, double[] rawLogR, int from, int to)
    {
        int runStart = from;
        for (int i = from + 1; i <= to; i++)
        {
            if (i == to || !levels[i].Equals(levels[runStart]))
            {
                Array.Fill(levels, Mean(rawLogR, runStart, i), runStart, i - runStart);
                runStart = i;
            }
        }
    }

    /// <summary>
    /// runASCAT segments: runs of equal segmented logR split at chromosome ends (<c>union(tlrend, tlrend.chr)</c>),
    /// each carrying the mirrored BAF of its first heterozygous locus (NaN if none).
    /// </summary>
    private static AspcfSegment[] BuildAspcfSegments(
        AlleleSpecificLocus[] loci, bool[] het, List<(int Lo, int Hi)> runs, double[] logRPcfed, double[] bafPcfed)
    {
        var segments = new List<AspcfSegment>();
        foreach ((int lo, int hi) in runs)
        {
            int start = lo;
            for (int i = lo + 1; i <= hi + 1; i++)
            {
                if (i == hi + 1 || !logRPcfed[i].Equals(logRPcfed[start]))
                {
                    double baf = double.NaN;
                    int hetCount = 0;
                    for (int k = start; k < i; k++)
                    {
                        if (het[k])
                        {
                            if (hetCount == 0)
                            {
                                baf = bafPcfed[k];
                            }

                            hetCount++;
                        }
                    }

                    segments.Add(new AspcfSegment(
                        loci[start].Chromosome, loci[start].Position, loci[i - 1].Position, logRPcfed[start], baf,
                        i - start, hetCount));
                    start = i;
                }
            }
        }

        return segments.ToArray();
    }

    /// <summary>
    /// One integer allele-specific copy-number state of a (possibly sub-clonal) segment, present in a given
    /// fraction of tumour cells. Mirrors the Battenberg output (Nik-Zainal et al. 2012, <i>Cell</i> 149:994:
    /// <c>nMaj1_A, nMin1_A, frac1_A</c>): a state is a (major, minor) integer pair plus the cellular fraction.
    /// </summary>
    /// <param name="MajorCopyNumber">Major-allele integer copy number of this state (≥ 0).</param>
    /// <param name="MinorCopyNumber">Minor-allele integer copy number of this state (≥ 0).</param>
    /// <param name="CellFraction">Fraction of tumour cells carrying this state (Battenberg <c>frac1_A</c>/<c>frac2_A</c>; 1 for a
    /// clonal state). Reported unclamped, as Battenberg does, so it can leave [0, 1] when the segment BAF lies beyond
    /// the nearest edge's corners.</param>
    public readonly record struct SubclonalCopyNumberState(
        int MajorCopyNumber,
        int MinorCopyNumber,
        double CellFraction)
    {
        /// <summary>Total integer copy number of this state (major + minor).</summary>
        public int TotalCopyNumber => MajorCopyNumber + MinorCopyNumber;
    }

    /// <summary>
    /// Sub-clonal copy-number fit of one segment under the Battenberg two-population model (Nik-Zainal et al. 2012,
    /// <i>Cell</i> 149:994): a segment is either <b>clonal</b> (one integer state in all tumour cells) or
    /// <b>sub-clonal</b> (a mixture of two adjacent integer states, fractions summing to 1).
    /// </summary>
    /// <param name="Segment">The segment these states describe.</param>
    /// <param name="PrimaryState">State 1 (Battenberg <c>nMaj1_A, nMin1_A, frac1_A</c>).</param>
    /// <param name="SecondaryState">State 2 (Battenberg <c>frac2</c>), or <c>null</c> for a clonal segment.</param>
    /// <param name="IsSubclonal">True when two states were needed (the observed CN was not (near-)integer).</param>
    public readonly record struct SubclonalSegmentFit(
        AlleleSpecificSegmentSummary Segment,
        SubclonalCopyNumberState PrimaryState,
        SubclonalCopyNumberState? SecondaryState,
        bool IsSubclonal);

    /// <summary>
    /// Maximum distance of an allele-specific copy number from the nearest integer below which the pre-2026-09
    /// implementation called a segment clonal. <b>No longer used</b> by <see cref="FitSubclonalCopyNumber"/>, which now
    /// applies Battenberg's own BAF-space tolerance <see cref="BattenbergMaxBafDistance"/>; retained for API compatibility.
    /// </summary>
    [Obsolete("Unused since FIN-B24/F14: Battenberg determine_copynumber uses maxdist; kept for binary compatibility.")]
    public const double SubclonalIntegerTolerance = 0.05;

    /// <summary>
    /// Battenberg <c>maxdist = 0.01</c> (<c>callSubclones</c> default): a segment whose BAF lies within 0.01 of the BAF of
    /// the closest clonal corner of its nearest edge is called clonal (<c>if (abs(l − test.level) &lt; maxdist) pval = 1</c>).
    /// </summary>
    public const double BattenbergMaxBafDistance = 0.01;

    /// <summary>Battenberg <c>cn_upper_limit = 1000</c>: major copy number used when the minor allele is negative at BAF = 1.</summary>
    private const double BattenbergCopyNumberUpperLimit = 1000.0;

    /// <summary>Battenberg: a negative minor copy number is raised to <c>nMinor = 0.01</c> (major raised along the BAF line).</summary>
    private const double BattenbergMinimumMinorCopyNumber = 0.01;

    /// <summary>
    /// Fits each segment's allele-specific copy number to one integer state (clonal) or a mixture of two integer states
    /// with a sub-clonal cellular fraction, porting Battenberg's <c>determine_copynumber</c> (Wedge-lab/battenberg,
    /// R/fitcopynumber.R; <c>orderEdges</c> in R/orderEdges.R; Nik-Zainal et al. 2012, <i>Cell</i> 149:994):
    /// <list type="number">
    /// <item>l = max(BAF, 1 − BAF); ψ_all = ρ·ψ + 2(1 − ρ);
    /// nMajor = (ρ − 1 + l·ψ_all·2^(logR/γ))/ρ, nMinor = (ρ − 1 + (1 − l)·ψ_all·2^(logR/γ))/ρ (the ASCAT equations).
    /// A negative nMinor is raised to 0.01 with nMajor moved along the BAF line (or set to 1000 when l = 1).</item>
    /// <item>The four corners (⌊nMaj⌋/⌈nMaj⌉ × ⌊nMin⌋/⌈nMin⌉) and their BAF levels (1 − ρ + ρ·nMaj)/(2 − 2ρ + ρ·(nMaj + nMin))
    /// select the nearest edge of the copy-number square (<c>orderEdges</c> option 1; the two states differ in one
    /// allele), using the total-copy-number (logR) priority <c>ntot &lt; x + y + 1</c>.</item>
    /// <item>Clonal when the closest corner of that edge has |l − level| &lt; <see cref="BattenbergMaxBafDistance"/>
    /// (reported as that corner with fraction 1); otherwise sub-clonal with state 1 at fraction
    /// τ = (1 − ρ + ρ·M₂ − 2l(1 − ρ) − lρ(m₂ + M₂)) / (lρ(m₁ + M₁) − lρ(m₂ + M₂) − ρM₁ + ρM₂) and state 2 at 1 − τ, the
    /// fraction that reproduces the segment BAF as a mixture of the two states (Battenberg reports τ unclamped).</item>
    /// </list>
    /// A segment summary carries a single BAF value, i.e. its SNP BAFs have zero spread; Battenberg's per-SNP
    /// t-test then returns <c>pval = 0</c> (<c>if (is.na(sd(BAFke)) || sd(BAFke) == 0) pval = 0</c>), so only the
    /// maxdist rule can make a segment clonal — exactly the path ported here (with per-SNP BAFs, Battenberg's t-test is
    /// run by <see cref="FitSubclonalCopyNumberWithSnpTest"/>; the alternative solutions B–F, SDfrac and the seeded
    /// bootstrap confidence intervals by <see cref="FitSubclonalCopyNumberWithBootstrap"/>).
    /// </summary>
    /// <param name="segments">Segment summaries (e.g. from <see cref="SegmentAlleleSpecificAspcf"/>). Non-null; finite logR,
    /// BAF in [0, 1].</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (the ASCAT/Battenberg <c>psit</c>, &gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>Per-segment clonal/sub-clonal copy-number fits in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">a segment has a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0,1], ψ ≤ 0, or γ ≤ 0.</exception>
    public static IReadOnlyList<SubclonalSegmentFit> FitSubclonalCopyNumber(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        double psiAll = MixtureCopiesPerCell(purity, ploidy); // psi = rho*psit + 2*(1-rho)
        var fits = new List<SubclonalSegmentFit>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            BattenbergSegmentGeometry g = BattenbergSegmentEdge(s, purity, psiAll, gamma, nameof(segments));
            // Constant (summary) BAF ⇒ Battenberg pval = 0; only maxdist can make the segment clonal.
            fits.Add(BattenbergSegmentFit(s, g, purity, g.ClosestDistance < BattenbergMaxBafDistance));
        }

        return fits;
    }

    /// <summary>
    /// Battenberg's significance level <c>siglevel = 0.05</c> (<c>callSubclones</c> default): a segment is sub-clonal when
    /// its SNP-BAF t-test p-value is ≤ siglevel (<c>if (pval[i] &lt;= siglevel)</c>).
    /// </summary>
    public const double BattenbergSignificanceLevel = 0.05;

    /// <summary>
    /// A segment summary together with the phased BAFs of its heterozygous SNPs (Battenberg <c>BAFphased</c> — the raw
    /// SNP BAF flipped to the segment's side, <c>ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>; used as given, not mirrored
    /// again), the input of Battenberg's per-segment clonality t-test.
    /// </summary>
    /// <param name="Segment">The segment (its <see cref="AlleleSpecificSegmentSummary.MeanBAF"/> is Battenberg's
    /// segment level <c>BAFseg</c>, its <see cref="AlleleSpecificSegmentSummary.MeanLogR"/> the segment logR).</param>
    /// <param name="PhasedSnpBafs">The segment's phased SNP BAFs (Battenberg <c>BAFke</c>), each in [0, 1]; may be empty.</param>
    public readonly record struct SubclonalSegmentSnpBafs(
        AlleleSpecificSegmentSummary Segment,
        IReadOnlyList<double> PhasedSnpBafs);

    /// <summary>
    /// A Battenberg segment fit with the clonality-test p-value Battenberg reports (<c>subcloneres$pval</c>).
    /// </summary>
    /// <param name="Fit">The clonal / sub-clonal copy-number fit.</param>
    /// <param name="PValue">Battenberg <c>pval</c>: 1 when the segment BAF is within maxdist of the closest corner level,
    /// 0 when the SNP BAFs have no spread (fewer than two SNPs, or constant), else the two-sided one-sample t-test p-value
    /// of the SNP BAFs against that corner level.</param>
    public readonly record struct SubclonalSegmentTestedFit(SubclonalSegmentFit Fit, double PValue);

    /// <summary>
    /// <see cref="FitSubclonalCopyNumber"/> with Battenberg's per-SNP clonality test (Wedge-lab/battenberg
    /// R/fitcopynumber.R, <c>determine_copynumber</c>; B24 F39). After the nearest edge and its closest corner level
    /// <c>test.level</c> are found exactly as in <see cref="FitSubclonalCopyNumber"/>:
    /// <code>
    /// if (is.na(sd(BAFke)) || sd(BAFke) == 0) pval = 0
    /// else pval = t.test(BAFke, alternative = "two.sided", mu = test.level)$p.value
    /// if (abs(l − test.level) &lt; maxdist) pval = 1
    /// sub-clonal ⇔ pval ≤ siglevel
    /// </code>
    /// so a segment whose BAF is beyond maxdist but whose SNP BAFs are too noisy to reject the clonal level is called
    /// clonal (and a segment with ≥ 2 distinct SNP BAFs that rejects it is sub-clonal). The t-test is
    /// <see cref="StatisticsHelper.OneSampleTTestPValue"/>; where R's <c>t.test</c> would stop ("data are essentially
    /// constant", Battenberg would abort) the p-value is taken as 0, as for exactly constant BAFs. With every
    /// segment's SNP BAFs constant this reproduces <see cref="FitSubclonalCopyNumber"/> exactly.
    /// </summary>
    /// <param name="segments">Segments with their phased SNP BAFs. Non-null; finite logR, BAF in [0, 1], SNP BAFs in [0, 1].</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <param name="significanceLevel">Battenberg <c>siglevel</c> ∈ [0, 1] (default <see cref="BattenbergSignificanceLevel"/>).</param>
    /// <param name="maxBafDistance">Battenberg <c>maxdist</c> ≥ 0 (default <see cref="BattenbergMaxBafDistance"/>).</param>
    /// <returns>Per-segment fits with Battenberg's p-value, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or a SNP-BAF list is null.</exception>
    /// <exception cref="ArgumentException">a segment has a non-finite logR, a BAF outside [0, 1], or a SNP BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0,1], ψ ≤ 0, γ ≤ 0, siglevel ∉ [0, 1] or maxdist &lt; 0 / not finite.</exception>
    public static IReadOnlyList<SubclonalSegmentTestedFit> FitSubclonalCopyNumberWithSnpTest(
        IReadOnlyList<SubclonalSegmentSnpBafs> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma,
        double significanceLevel = BattenbergSignificanceLevel,
        double maxBafDistance = BattenbergMaxBafDistance)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ValidateAscatModelParameters(purity, ploidy, gamma);
        if (!(significanceLevel >= 0.0 && significanceLevel <= 1.0))
        {
            throw new ArgumentOutOfRangeException(nameof(significanceLevel), significanceLevel, "siglevel must lie in [0, 1].");
        }

        if (!(maxBafDistance >= 0.0) || double.IsPositiveInfinity(maxBafDistance))
        {
            throw new ArgumentOutOfRangeException(nameof(maxBafDistance), maxBafDistance, "maxdist must be finite and ≥ 0.");
        }

        double psiAll = MixtureCopiesPerCell(purity, ploidy);
        var fits = new List<SubclonalSegmentTestedFit>(segments.Count);
        foreach (SubclonalSegmentSnpBafs item in segments)
        {
            if (item.PhasedSnpBafs is null)
            {
                throw new ArgumentNullException(nameof(segments), "Every segment needs a (possibly empty) phased SNP-BAF list.");
            }

            foreach (double b in item.PhasedSnpBafs)
            {
                if (!(b >= 0.0 && b <= 1.0))
                {
                    throw new ArgumentException("Every phased SNP BAF must lie in [0, 1].", nameof(segments));
                }
            }

            BattenbergSegmentGeometry g = BattenbergSegmentEdge(item.Segment, purity, psiAll, gamma, nameof(segments));
            double pValue = StatisticsHelper.OneSampleTTestPValue(item.PhasedSnpBafs, g.ClosestTestLevel);
            if (double.IsNaN(pValue))
            {
                pValue = 0.0; // sd NA (< 2 SNPs) or 0 (constant BAF): Battenberg pval = 0
            }

            if (g.ClosestDistance < maxBafDistance)
            {
                pValue = 1.0;
            }

            fits.Add(new SubclonalSegmentTestedFit(
                BattenbergSegmentFit(item.Segment, g, purity, isClonal: !(pValue <= significanceLevel)),
                pValue));
        }

        return fits;
    }

    /// <summary>Battenberg per-segment geometry: mirrored BAF l, nearest edge, and the closest corner of that edge.</summary>
    private readonly record struct BattenbergSegmentGeometry(
        double L, double Maj1, double Min1, double Maj2, double Min2, bool FirstClosest, double ClosestDistance)
    {
        public double ClosestTestLevel => FirstClosest ? Level1 : Level2;

        public double Level1 { get; init; }

        public double Level2 { get; init; }

        /// <summary>Battenberg <c>levels[2]</c>/<c>levels[3]</c> (corner BAFs of the square, 0/0 ⇒ 0.5), x, y and ntot.</summary>
        public double SquareLevel2 { get; init; }

        public double SquareLevel3 { get; init; }

        public double X { get; init; }

        public double Y { get; init; }

        public double Ntot { get; init; }
    }

    /// <summary>The <c>determine_copynumber</c> body up to <c>whichclosestlevel.test</c> (shared by both clonality rules).</summary>
    private static BattenbergSegmentGeometry BattenbergSegmentEdge(
        AlleleSpecificSegmentSummary s, double rho, double psiAll, double gamma, string paramName)
    {
        if (!IsValidAlleleSignal(s.MeanLogR, s.MeanBAF))
        {
            throw new ArgumentException("Every segment needs a finite mean logR and a mean BAF in [0, 1].", paramName);
        }

        double l = MirrorBaf(s.MeanBAF); // max(BAF, 1 − BAF)
        double scaled = psiAll * Math.Pow(2.0, s.MeanLogR / gamma);
        double nMajor = (rho - 1.0 + l * scaled) / rho;
        double nMinor = (rho - 1.0 + (1.0 - l) * scaled) / rho;

        // Increase nMajor and nMinor together, to avoid impossible combinations (negative sub-clonal fractions).
        if (nMinor < 0.0)
        {
            nMajor = l == 1.0
                ? BattenbergCopyNumberUpperLimit
                : nMajor + l * (BattenbergMinimumMinorCopyNumber - nMinor) / (1.0 - l);
            nMinor = BattenbergMinimumMinorCopyNumber;
        }

        double x = Math.Floor(nMinor), y = Math.Floor(nMajor);
        double ntot = nMajor + nMinor;
        // Corners, sorted in the order of ascending BAF: (⌊M⌋,⌈m⌉), (⌈M⌉,⌈m⌉), (⌊M⌋,⌊m⌋), (⌈M⌉,⌊m⌋).
        double level2 = BattenbergCornerLevel(Math.Ceiling(nMajor), Math.Ceiling(nMinor), rho, zeroCornerIsBalanced: true);
        double level3 = BattenbergCornerLevel(Math.Floor(nMajor), Math.Floor(nMinor), rho, zeroCornerIsBalanced: true);
        (double maj1, double min1, double maj2, double min2) = BattenbergNearestEdge(level2, level3, l, ntot, x, y);

        // Clonality test on the corners of the nearest edge (test.levels carry no 0/0 correction in Battenberg).
        double testLevel1 = BattenbergCornerLevel(maj1, min1, rho, zeroCornerIsBalanced: false);
        double testLevel2 = BattenbergCornerLevel(maj2, min2, rho, zeroCornerIsBalanced: false);
        double dist1 = double.IsNaN(testLevel1) ? double.PositiveInfinity : Math.Abs(testLevel1 - l);
        double dist2 = double.IsNaN(testLevel2) ? double.PositiveInfinity : Math.Abs(testLevel2 - l);
        bool firstClosest = dist1 <= dist2; // which.min: first index on a tie
        return new BattenbergSegmentGeometry(l, maj1, min1, maj2, min2, firstClosest, firstClosest ? dist1 : dist2)
        {
            Level1 = testLevel1,
            Level2 = testLevel2,
            SquareLevel2 = level2,
            SquareLevel3 = level3,
            X = x,
            Y = y,
            Ntot = ntot,
        };
    }

    /// <summary>
    /// Battenberg output row: the closest corner with fraction 1 when clonal, else state 1 at τ and state 2 at 1 − τ.
    /// </summary>
    private static SubclonalSegmentFit BattenbergSegmentFit(
        AlleleSpecificSegmentSummary s, BattenbergSegmentGeometry g, double rho, bool isClonal)
    {
        double l = g.L, maj1 = g.Maj1, min1 = g.Min1, maj2 = g.Maj2, min2 = g.Min2;
        if (isClonal)
        {
            var clonal = g.FirstClosest
                ? new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), 1.0)
                : new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0);
            return new SubclonalSegmentFit(s, clonal, SecondaryState: null, IsSubclonal: false);
        }

        double tau = BattenbergTau(l, rho, maj1, min1, maj2, min2);
        return new SubclonalSegmentFit(
            s,
            new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), tau),
            new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0 - tau),
            IsSubclonal: true);
    }

    /// <summary>
    /// Battenberg corner BAF level <c>(1 − ρ + ρ·nMaj)/(2 − 2ρ + ρ·(nMaj + nMin))</c>; for the corner (0, 0) at ρ = 1
    /// (0/0) the square's <c>levels</c> are set to 0.5 while the edge <c>test.levels</c> stay NaN.
    /// </summary>
    private static double BattenbergCornerLevel(double nMaj, double nMin, double rho, bool zeroCornerIsBalanced)
    {
        if (zeroCornerIsBalanced && nMaj == 0.0 && nMin == 0.0)
        {
            return BalancedBaf;
        }

        return (1.0 - rho + rho * nMaj) / (2.0 - 2.0 * rho + rho * (nMaj + nMin));
    }

    /// <summary>
    /// Battenberg <c>orderEdges(levels, l, ntot, x, y)</c>, first option (the nearest edge): case 1/2a when
    /// l &gt; levels[3], case 2c when l &gt; levels[2], else case 2b; within each case the logR criterion
    /// <c>ntot &lt; x + y + 1</c> picks the lower- or higher-total edge. Returns (nMaj1, nMin1, nMaj2, nMin2).
    /// </summary>
    private static (double Maj1, double Min1, double Maj2, double Min2) BattenbergNearestEdge(
        double level2, double level3, double l, double ntot, double x, double y)
    {
        bool lowerTotal = ntot < x + y + 1.0;
        if (l > level3)
        {
            // case 1 or 2a
            return lowerTotal ? (y, x, y + 1.0, x) : (y + 1.0, x, y + 1.0, x + 1.0);
        }

        if (l > level2)
        {
            // case 2c
            return lowerTotal ? (y, x, y, x + 1.0) : (y + 1.0, x, y + 1.0, x + 1.0);
        }

        // case 2b
        return lowerTotal ? (y, x, y, x + 1.0) : (y, x + 1.0, y + 1.0, x + 1.0);
    }

    /// <summary>
    /// One of Battenberg's six candidate sub-clonal solutions A–F of a segment (<c>callSubclones</c> extended output
    /// columns <c>nMaj1_X, nMin1_X, frac1_X, nMaj2_X, nMin2_X, frac2_X, SDfrac_X, SDfrac_X_BS, frac1_X_0.025, frac1_X_0.975</c>).
    /// A solution whose edge would need a negative copy number is <c>NA</c> in Battenberg: its copy numbers are
    /// <c>null</c> and its numbers NaN (such solutions are listed last).
    /// </summary>
    /// <param name="MajorCopyNumber1">nMaj1 (state 1 major), or null (NA).</param>
    /// <param name="MinorCopyNumber1">nMin1, or null.</param>
    /// <param name="Fraction1">frac1 = τ (unclamped), NaN when NA.</param>
    /// <param name="MajorCopyNumber2">nMaj2 (state 2 major), or null.</param>
    /// <param name="MinorCopyNumber2">nMin2, or null.</param>
    /// <param name="Fraction2">frac2 = 1 − τ.</param>
    /// <param name="FractionSd">SDfrac: |τ(l + s) − τ|/2 + |τ(l − s) − τ|/2 with s = sd(BAFke)/√n (NaN for &lt; 2 SNPs).</param>
    /// <param name="FractionBootstrapSd">SDfrac_BS: sd of the bootstrap τ values (NaN when any is NaN).</param>
    /// <param name="Fraction1Lower">frac1_0.025: the 25th smallest bootstrap τ (NaN when fewer than 25 non-NaN values).</param>
    /// <param name="Fraction1Upper">frac1_0.975: the 975th smallest bootstrap τ (NaN when fewer than 975).</param>
    public readonly record struct BattenbergSubclonalSolution(
        int? MajorCopyNumber1,
        int? MinorCopyNumber1,
        double Fraction1,
        int? MajorCopyNumber2,
        int? MinorCopyNumber2,
        double Fraction2,
        double FractionSd,
        double FractionBootstrapSd,
        double Fraction1Lower,
        double Fraction1Upper);

    /// <summary>
    /// A Battenberg <c>determine_copynumber</c> output row with its alternative solutions and confidence intervals.
    /// </summary>
    /// <param name="Fit">Solution A as a clonal / sub-clonal fit (identical to <see cref="FitSubclonalCopyNumberWithSnpTest"/>).</param>
    /// <param name="PValue">Battenberg <c>pval</c>.</param>
    /// <param name="Baf">Battenberg <c>BAF</c> column: the mirrored segment level l = max(BAFseg, 1 − BAFseg).</param>
    /// <param name="TotalCopyNumber">Battenberg <c>ntot</c> = nMajor + nMinor (raw, after the negative-minor correction).</param>
    /// <param name="Solutions">Solutions A–F (six entries) for a sub-clonal segment; empty for a clonal one (Battenberg NA).</param>
    public sealed record BattenbergSegmentCall(
        SubclonalSegmentFit Fit,
        double PValue,
        double Baf,
        double TotalCopyNumber,
        IReadOnlyList<BattenbergSubclonalSolution> Solutions);

    /// <summary>Battenberg <c>noperms = 1000</c> (<c>callSubclones</c> default).</summary>
    public const int BattenbergDefaultPermutations = 1000;

    /// <summary>
    /// <see cref="FitSubclonalCopyNumberWithSnpTest"/> plus the rest of Battenberg's <c>determine_copynumber</c> row
    /// (Wedge-lab/battenberg R/fitcopynumber.R; B24 F41): for a sub-clonal segment all six <c>orderEdges</c> solutions
    /// A–F (rows with a negative copy number become NA and move to the end), each with τ, the delta-method
    /// <c>SDfrac</c> and the bootstrap <c>SDfrac_BS</c> / <c>frac1_0.025</c> / <c>frac1_0.975</c>:
    /// <code>
    /// for (option in 1:6) for (j in 1:noperms) {
    ///   permBAFs = sample(BAFke, length(BAFke), replace = T); permMeanBAF = mean(permBAFs)
    ///   permFraction[j] = τ(permMeanBAF; option) }
    /// SDfrac_BS = sd(permFraction); frac1_0.025 = sort(permFraction)[25]; frac1_0.975 = sort(permFraction)[975]
    /// </code>
    /// Note the bounds are fixed order statistics, not <c>quantile()</c>: they are the 2.5 %/97.5 % points only for
    /// noperms = 1000 and NA for fewer resamples. Resampling uses a port of R's default generator — Mersenne-Twister
    /// seeded by <c>set.seed(seed)</c> (initial scrambling + LCG fill) with R ≥ 3.6 "Rejection" <c>sample</c>
    /// (<c>R_unif_index</c>/<c>rbits</c>) — consuming draws in Battenberg's order (segments in input order, options
    /// A–F, NA options included), so the resamples are R's: <c>set.seed(seed); determine_copynumber(…)</c> is reproduced
    /// (the bootstrap columns agree to ~1e−15; R's long-double <c>mean</c>/<c>sd</c> differ only in the last bits).
    /// <c>callSubclones</c> calls <c>set.seed</c> once and then <c>determine_copynumber</c> twice (before/after
    /// <c>merge_segments</c>), so its final table matches this method only for the RNG state at the second call.
    /// </summary>
    /// <param name="segments">Segments with their phased SNP BAFs (e.g. from <see cref="BuildBattenbergSegments"/>).</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (&gt; 0).</param>
    /// <param name="seed">R <c>set.seed</c> value (Battenberg <c>seed</c>; any 32-bit integer).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <param name="significanceLevel">Battenberg <c>siglevel</c> ∈ [0, 1].</param>
    /// <param name="maxBafDistance">Battenberg <c>maxdist</c> ≥ 0.</param>
    /// <param name="permutations">Battenberg <c>noperms</c> ≥ 1 (default 1000).</param>
    /// <returns>One call per segment, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or a SNP-BAF list is null.</exception>
    /// <exception cref="ArgumentException">invalid segment or SNP BAF (as <see cref="FitSubclonalCopyNumberWithSnpTest"/>).</exception>
    /// <exception cref="ArgumentOutOfRangeException">invalid ρ, ψ, γ, siglevel, maxdist, or permutations &lt; 1.</exception>
    public static IReadOnlyList<BattenbergSegmentCall> FitSubclonalCopyNumberWithBootstrap(
        IReadOnlyList<SubclonalSegmentSnpBafs> segments,
        double purity,
        double ploidy,
        int seed,
        double gamma = AscatSequencingGamma,
        double significanceLevel = BattenbergSignificanceLevel,
        double maxBafDistance = BattenbergMaxBafDistance,
        int permutations = BattenbergDefaultPermutations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(permutations, 1);
        IReadOnlyList<SubclonalSegmentTestedFit> tested = FitSubclonalCopyNumberWithSnpTest(
            segments, purity, ploidy, gamma, significanceLevel, maxBafDistance);
        double psiAll = MixtureCopiesPerCell(purity, ploidy);
        var rng = new RMersenneTwister(seed);
        var calls = new List<BattenbergSegmentCall>(tested.Count);
        for (int i = 0; i < tested.Count; i++)
        {
            SubclonalSegmentTestedFit t = tested[i];
            BattenbergSegmentGeometry g = BattenbergSegmentEdge(segments[i].Segment, purity, psiAll, gamma, nameof(segments));
            IReadOnlyList<BattenbergSubclonalSolution> solutions = t.Fit.IsSubclonal
                ? BattenbergSolutions(g, purity, segments[i].PhasedSnpBafs, permutations, rng)
                : Array.Empty<BattenbergSubclonalSolution>();
            calls.Add(new BattenbergSegmentCall(t.Fit, t.PValue, g.L, g.Ntot, solutions));
        }

        return calls;
    }

    /// <summary>The sub-clonal branch of <c>determine_copynumber</c>: six solutions with SDfrac and bootstrap columns.</summary>
    private static BattenbergSubclonalSolution[] BattenbergSolutions(
        BattenbergSegmentGeometry g, double rho, IReadOnlyList<double> bafke, int noperms, RMersenneTwister rng)
    {
        (double Maj1, double Min1, double Maj2, double Min2)[] edges = BattenbergAllEdges(g.SquareLevel2, g.SquareLevel3, g.L, g.Ntot, g.X, g.Y);
        double l = g.L;
        int n = bafke.Count;
        double sdl = n < 2 ? double.NaN : RSampleSd(bafke) / Math.Sqrt(n);
        var result = new BattenbergSubclonalSolution[edges.Length];
        var perm = new double[n];
        var permFraction = new double[noperms];
        for (int o = 0; o < edges.Length; o++)
        {
            (double maj1, double min1, double maj2, double min2) = edges[o];
            double tau = BattenbergTau(l, rho, maj1, min1, maj2, min2);
            double sdtau = Math.Abs(BattenbergTau(l + sdl, rho, maj1, min1, maj2, min2) - tau) / 2
                           + Math.Abs(BattenbergTau(l - sdl, rho, maj1, min1, maj2, min2) - tau) / 2;
            for (int j = 0; j < noperms; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    perm[k] = bafke[rng.UnifIndex(n)];
                }

                double permMean = n == 0 ? double.NaN : RMean(perm);
                permFraction[j] = BattenbergTau(permMean, rho, maj1, min1, maj2, min2);
            }

            double sdBoot = RSampleSd(permFraction);
            var ordered = permFraction.Where(v => !double.IsNaN(v)).ToArray(); // sort() drops NA/NaN
            Array.Sort(ordered);
            bool na = double.IsNaN(maj1);
            result[o] = new BattenbergSubclonalSolution(
                na ? null : AscatCopyNumberToInt(maj1),
                na ? null : AscatCopyNumberToInt(min1),
                tau,
                na ? null : AscatCopyNumberToInt(maj2),
                na ? null : AscatCopyNumberToInt(min2),
                1.0 - tau,
                sdtau,
                sdBoot,
                ordered.Length >= 25 ? ordered[24] : double.NaN,
                ordered.Length >= 975 ? ordered[974] : double.NaN);
        }

        return result;
    }

    /// <summary>
    /// Battenberg sub-clonal fraction of state 1 at BAF <paramref name="l"/>:
    /// τ = (1 − ρ + ρM₂ − 2l(1 − ρ) − lρ(m₂ + M₂)) / (lρ(m₁ + M₁) − lρ(m₂ + M₂) − ρM₁ + ρM₂) (R operation order).
    /// </summary>
    private static double BattenbergTau(double l, double rho, double maj1, double min1, double maj2, double min2) =>
        (1.0 - rho + rho * maj2 - 2.0 * l * (1.0 - rho) - l * rho * (min2 + maj2))
        / (l * rho * (min1 + maj1) - l * rho * (min2 + maj2) - rho * maj1 + rho * maj2);

    /// <summary>R <c>sd(x)</c>: √(Σ(x − x̄)²/(n − 1)) with R's refined mean; NaN for n &lt; 2 or any NaN.</summary>
    private static double RSampleSd(IReadOnlyList<double> x)
    {
        int n = x.Count;
        if (n < 2)
        {
            return double.NaN;
        }

        double mean = RMean(x);
        double ss = 0.0;
        for (int i = 0; i < n; i++)
        {
            double d = x[i] - mean;
            ss += d * d;
        }

        return Math.Sqrt(ss / (n - 1));
    }

    /// <summary>
    /// Battenberg <c>orderEdges(levels, l, ntot, x, y)</c>, all six options (R/orderEdges.R), negative-copy-number rows
    /// set to NaN and moved to the end (<c>determine_copynumber</c>: <c>rbind(all.edges[-na.indices,], all.edges[na.indices,])</c>).
    /// </summary>
    private static (double Maj1, double Min1, double Maj2, double Min2)[] BattenbergAllEdges(
        double level2, double level3, double l, double ntot, double x, double y)
    {
        double[] maj1, min1, maj2, min2;
        bool lower = ntot < x + y + 1.0;
        if (l > level3)
        {
            // case 1 or 2a
            if (lower)
            {
                maj1 = new[] { y, y - 1, y, y + 1, y + 1, y + 1 };
                min1 = new[] { x, x, x, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 2, y + 1, y + 1, y + 1 };
                min2 = new[] { x, x, x, x + 1, x + 1, x + 2 };
            }
            else
            {
                maj1 = new[] { y + 1, y + 1, y + 1, y, y - 1, y };
                min1 = new[] { x, x - 1, x, x, x, x };
                maj2 = new[] { y + 1, y + 1, y + 1, y + 1, y + 1, y + 2 };
                min2 = new[] { x + 1, x + 1, x + 2, x, x, x };
            }
        }
        else if (l > level2)
        {
            // case 2c
            if (lower)
            {
                maj1 = new[] { y, y, y, y + 1, y + 1, y + 1 };
                min1 = new[] { x, x - 1, x, x, x - 1, x };
                maj2 = new[] { y, y, y, y + 1, y + 1, y + 1 };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 2 };
            }
            else
            {
                maj1 = new[] { y + 1, y + 1, y + 1, y, y, y };
                min1 = new[] { x, x - 1, x, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 1, y, y, y };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 2 };
            }
        }
        else
        {
            // case 2b
            if (lower)
            {
                maj1 = new[] { y, y, y, y, y - 1, y };
                min1 = new[] { x, x - 1, x, x + 1, x + 1, x + 1 };
                maj2 = new[] { y, y, y, y + 1, y + 1, y + 2 };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 1 };
            }
            else
            {
                maj1 = new[] { y, y - 1, y, y, y, y };
                min1 = new[] { x + 1, x + 1, x + 1, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 2, y, y, y };
                min2 = new[] { x + 1, x + 1, x + 1, x + 1, x + 1, x + 2 };
            }
        }

        var valid = new List<(double, double, double, double)>(6);
        var invalid = new List<(double, double, double, double)>(6);
        for (int i = 0; i < 6; i++)
        {
            if (maj1[i] < 0 || min1[i] < 0 || maj2[i] < 0 || min2[i] < 0)
            {
                invalid.Add((double.NaN, double.NaN, double.NaN, double.NaN));
            }
            else
            {
                valid.Add((maj1[i], min1[i], maj2[i], min2[i]));
            }
        }

        valid.AddRange(invalid);
        return valid.ToArray();
    }

    /// <summary>
    /// R's default RNG (RNG.c): Mersenne-Twister (<c>MT_genrand</c>, unsigned 32-bit state, output ·2.3283064365386963e−10,
    /// <c>fixup</c> into (0, 1)), seeded as <c>set.seed(seed)</c> does (<c>Randomize</c>: 50 LCG scrambles
    /// <c>seed = 69069·seed + 1</c>, then 625 LCG values fill <c>dummy[0..624]</c>, <c>dummy[0] = mti = 624</c>), with R ≥ 3.6
    /// <c>sample.kind = "Rejection"</c> integer draws (<c>R_unif_index</c>: <c>rbits(ceil(log2(n)))</c> until &lt; n).
    /// E.g. <c>set.seed(42); sample(5, 10, TRUE)</c> = 1 5 1 1 2 4 2 2 1 4.
    /// </summary>
    private sealed class RMersenneTwister
    {
        private const int N = 624;
        private const int M = 397;
        private const uint MatrixA = 0x9908b0df;
        private const uint UpperMask = 0x80000000;
        private const uint LowerMask = 0x7fffffff;
        private const double TwoTo32Inverse = 2.3283064365386963e-10;
        private const double I2To32M1 = 2.328306437080797e-10;
        private readonly uint[] _mt = new uint[N];
        private int _mti;

        public RMersenneTwister(int seed)
        {
            uint s = unchecked((uint)seed);
            for (int j = 0; j < 50; j++)
            {
                s = unchecked((69069u * s) + 1u);
            }

            s = unchecked((69069u * s) + 1u); // dummy[0] (mti), overwritten by FixupSeeds
            for (int j = 0; j < N; j++)
            {
                s = unchecked((69069u * s) + 1u);
                _mt[j] = s;
            }

            _mti = N;
        }

        /// <summary>R <c>unif_rand()</c>.</summary>
        public double UnifRand()
        {
            double v = Genrand();
            if (v <= 0.0)
            {
                return 0.5 * I2To32M1;
            }

            return 1.0 - v <= 0.0 ? 1.0 - (0.5 * I2To32M1) : v;
        }

        /// <summary>R <c>R_unif_index(n)</c> (rejection sampling); 0-based index in [0, n).</summary>
        public int UnifIndex(int n)
        {
            if (n <= 0)
            {
                return 0;
            }

            int bits = (int)Math.Ceiling(Math.Log2(n));
            double dv;
            do
            {
                long v = 0;
                for (int k = 0; k <= bits; k += 16)
                {
                    int v1 = (int)Math.Floor(UnifRand() * 65536);
                    v = (65536 * v) + v1;
                }

                dv = v & ((1L << bits) - 1);
            }
            while (n <= dv);

            return (int)dv;
        }

        private double Genrand()
        {
            uint y;
            if (_mti >= N)
            {
                int kk;
                for (kk = 0; kk < N - M; kk++)
                {
                    y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                    _mt[kk] = _mt[kk + M] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
                }

                for (; kk < N - 1; kk++)
                {
                    y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                    _mt[kk] = _mt[kk + (M - N)] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
                }

                y = (_mt[N - 1] & UpperMask) | (_mt[0] & LowerMask);
                _mt[N - 1] = _mt[M - 1] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
                _mti = 0;
            }

            y = _mt[_mti++];
            y ^= y >> 11;
            y ^= (y << 7) & 0x9d2c5680u;
            y ^= (y << 15) & 0xefc60000u;
            y ^= y >> 18;
            return y * TwoTo32Inverse;
        }
    }

    #endregion

    #region Battenberg phased-BAF segmentation (ONCO-ASCAT-001, B24 F40)

    /// <summary>
    /// One germline-heterozygous SNP with its tumour BAF, expressed on a <b>caller-supplied</b> haplotype phase (Battenberg
    /// <c>combine.baf.files</c> output, columns <c>Chromosome, Position, BAF</c>; the BAF is that of haplotype 1 after
    /// IMPUTE2/Beagle5 phasing). <see cref="double.NaN"/> marks a missing BAF (dropped, as Battenberg drops <c>NA</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">SNP position (the SNPs of a chromosome are taken in input order, as in Battenberg).</param>
    /// <param name="Baf">Haplotype-phased BAF in [0, 1], or NaN.</param>
    public readonly record struct PhasedBafSnp(string Chromosome, long Position, double Baf);

    /// <summary>A prior breakpoint (e.g. a structural-variant junction): Battenberg <c>prior_breakpoints_file</c> row.</summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">Breakpoint position.</param>
    public readonly record struct BattenbergPriorBreakpoint(string Chromosome, long Position);

    /// <summary>Battenberg <c>calc_seg_baf_option</c>: how the BAF of a segment is recalculated after segmentation.</summary>
    public enum BattenbergSegmentBafOption
    {
        /// <summary>1 — median of the segment's phased BAFs (<c>adjustSegmValues</c>).</summary>
        Median = 1,

        /// <summary>2 — the PCF segment mean, unchanged.</summary>
        Mean = 2,

        /// <summary>3 (Battenberg default) — the median, unless it is exactly 0 or 1, then the mean.</summary>
        MedianUnlessExtreme = 3,
    }

    /// <summary>
    /// Parameters of Battenberg <c>segment.baf.phased</c> (Wedge-lab/battenberg R/segmentation.R); defaults are
    /// Battenberg's (<c>gamma = 10, phasegamma = 3, kmin = 3, phasekmin = 3, no_segmentation = F,
    /// calc_seg_baf_option = 3</c>, no prior breakpoints).
    /// </summary>
    public sealed record BattenbergPhasedSegmentationOptions
    {
        /// <summary>Battenberg defaults.</summary>
        public static BattenbergPhasedSegmentationOptions Default { get; } = new();

        /// <summary><c>gamma</c>: PCF penalty factor of the second (copy-number) segmentation (multiplied by the MAD sd).</summary>
        public double Gamma { get; init; } = 10.0;

        /// <summary><c>phasegamma</c>: PCF penalty factor of the first (phase-correcting) segmentation.</summary>
        public double PhaseGamma { get; init; } = 3.0;

        /// <summary><c>kmin</c>: minimum SNPs per segment of the second segmentation, in [1, <see cref="BattenbergMaxKmin"/>].</summary>
        public int Kmin { get; init; } = 3;

        /// <summary><c>phasekmin</c>: minimum SNPs per segment of the first segmentation, in [1, <see cref="BattenbergMaxKmin"/>].</summary>
        public int PhaseKmin { get; init; } = 3;

        /// <summary><c>no_segmentation</c>: switch the haplotype blocks but take the mean phased BAF as the segment BAF.</summary>
        public bool NoSegmentation { get; init; }

        /// <summary><c>calc_seg_baf_option</c> (default <see cref="BattenbergSegmentBafOption.MedianUnlessExtreme"/>).</summary>
        public BattenbergSegmentBafOption SegmentBafOption { get; init; } = BattenbergSegmentBafOption.MedianUnlessExtreme;

        /// <summary>Prior breakpoints (<c>prior_breakpoints_file</c>), in file order; empty = none.</summary>
        public IReadOnlyList<BattenbergPriorBreakpoint> PriorBreakpoints { get; init; } = Array.Empty<BattenbergPriorBreakpoint>();
    }

    /// <summary>
    /// Largest <c>kmin</c>/<c>phasekmin</c> accepted: Battenberg's <c>filterMarkS4</c> indexes <c>3·kmin + 6</c> SNPs,
    /// so larger values leave R's index ranges for the minimum PCF input of 50 SNPs (R then misbehaves or stops).
    /// </summary>
    public const int BattenbergMaxKmin = 14;

    /// <summary>
    /// One output row of Battenberg <c>segment.baf.phased</c> (<c>Chromosome, Position, BAF, BAFphased, BAFseg</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">SNP position.</param>
    /// <param name="Baf">The input phased BAF.</param>
    /// <param name="BafPhased">Haplotype-block-corrected BAF <c>ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>.</param>
    /// <param name="BafSegment">The segment BAF <c>BAFseg</c> (per <see cref="BattenbergSegmentBafOption"/>).</param>
    public readonly record struct PhasedBafSegmentedSnp(
        string Chromosome, long Position, double Baf, double BafPhased, double BafSegment);

    /// <summary>A raw logR probe (Battenberg <c>LogRvals</c>: chromosome, position, logR). NaN/±∞ logR are ignored.</summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">Probe position.</param>
    /// <param name="LogR">Raw logR.</param>
    public readonly record struct LogRProbe(string Chromosome, long Position, double LogR);

    /// <summary>Battenberg <c>bkps_to_presegment_breakpoints</c>: <c>maxsnpdist = 3000000</c> (a gap ≥ this splits).</summary>
    private const long BattenbergMaxSnpDistance = 3_000_000;

    /// <summary>Battenberg <c>run_pcf</c>: presegments with fewer SNPs are not segmented (their mean BAF is used).</summary>
    private const int BattenbergMinPcfSnps = 50;

    /// <summary>Battenberg <c>run_pcf</c>: MAD sd floor 0.09 ("binomial around 0.5 at depth 30").</summary>
    private const double BattenbergMinBafSd = 0.09;

    /// <summary>
    /// Segments caller-phased SNP BAFs exactly as Battenberg <c>segment.baf.phased</c> (Wedge-lab/battenberg
    /// R/segmentation.R, default path; B24 F40), the step between haplotype phasing and <c>callSubclones</c>:
    /// <list type="number">
    /// <item>Per chromosome (first-appearance order), SNPs with a missing BAF are dropped and the chromosome is
    /// pre-segmented at prior breakpoints and at every gap ≥ 3 Mb between consecutive SNPs
    /// (<c>bkps_to_presegment_breakpoints</c>/<c>addin_bigholes</c>, including R's index quirks).</item>
    /// <item>Per presegment (<c>run_pcf</c>): sd = <c>getMad(min(BAF, 1 − BAF), k = 25)</c> (NA → 0, floored at 0.09);
    /// a first PCF <c>selectFastPcf(BAF, phasekmin, phasegamma·sd)</c> finds the switched haplotype blocks,
    /// <c>BAFphased = ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>; a second PCF <c>selectFastPcf(BAFphased, kmin,
    /// gamma·sd)</c> gives <c>BAFseg</c> (presegments of &lt; 50 SNPs, or <c>no_segmentation</c>, use the mean); the
    /// segment BAF is then recalculated per <c>calc_seg_baf_option</c>.</item>
    /// </list>
    /// <c>selectFastPcf</c> is the copynumber-package fast PCF that Battenberg ships in R/fastPCF.R (<c>filterMarkS4</c>
    /// candidate breakpoints + <c>PottsCompact</c>, and the 5000-SNP windowed <c>runPcfSubset</c> above 15 000 SNPs) —
    /// a different routine from ASCAT's <c>exactPcf</c>/<c>fastAspcf</c>, ported here once; <c>getMad</c>/<c>runmed</c>
    /// are the shared ASCAT helpers (identical code). Haplotype phasing itself (IMPUTE2/Beagle5 against the 1000 Genomes
    /// reference panel) is <b>not</b> run: the caller supplies phased BAFs. Agreement with Battenberg (R 4.3.3) is to
    /// ≲ 1e−12; R accumulates <c>cumsum</c>/<c>mean</c> in 80-bit long double, which can only matter at exact cost ties.
    /// </summary>
    /// <param name="snps">Phased SNP BAFs. Non-null; non-null chromosome; BAF in [0, 1] or NaN.</param>
    /// <param name="options">Battenberg parameters (default <see cref="BattenbergPhasedSegmentationOptions.Default"/>).</param>
    /// <returns>Battenberg's per-SNP rows (chromosomes in first-appearance order, presegments in order).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="snps"/>, a chromosome, or the breakpoint list is null.</exception>
    /// <exception cref="ArgumentException">a BAF outside [0, 1] (and not NaN).</exception>
    /// <exception cref="ArgumentOutOfRangeException">gamma/phasegamma not finite or &lt; 0; kmin/phasekmin outside [1, 14];
    /// an undefined <see cref="BattenbergSegmentBafOption"/>.</exception>
    public static IReadOnlyList<PhasedBafSegmentedSnp> SegmentPhasedBaf(
        IReadOnlyList<PhasedBafSnp> snps,
        BattenbergPhasedSegmentationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(snps);
        options ??= BattenbergPhasedSegmentationOptions.Default;
        ValidateBattenbergSegmentationOptions(options);

        var chromosomes = new List<string>();
        var byChromosome = new Dictionary<string, (List<long> Pos, List<double> Baf)>(StringComparer.Ordinal);
        foreach (PhasedBafSnp snp in snps)
        {
            if (snp.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(snps), "Every SNP needs a chromosome.");
            }

            if (!(double.IsNaN(snp.Baf) || (snp.Baf >= 0.0 && snp.Baf <= 1.0)))
            {
                throw new ArgumentException("Every phased BAF must lie in [0, 1] (or be NaN = missing).", nameof(snps));
            }

            if (!byChromosome.TryGetValue(snp.Chromosome, out var lists))
            {
                lists = (new List<long>(), new List<double>());
                byChromosome.Add(snp.Chromosome, lists);
                chromosomes.Add(snp.Chromosome);
            }

            if (!double.IsNaN(snp.Baf))
            {
                lists.Pos.Add(snp.Position);
                lists.Baf.Add(snp.Baf);
            }
        }

        var output = new List<PhasedBafSegmentedSnp>(snps.Count);
        foreach (string chr in chromosomes)
        {
            (List<long> posList, List<double> bafList) = byChromosome[chr];
            long[] pos = posList.ToArray();
            double[] baf = bafList.ToArray();
            var bkps = new List<long>();
            foreach (BattenbergPriorBreakpoint b in options.PriorBreakpoints)
            {
                if (b.Chromosome == chr)
                {
                    bkps.Add(b.Position);
                }
            }

            foreach ((long? start, long? end) in BattenbergPresegments(pos, bkps))
            {
                BattenbergRunPcf(chr, pos, baf, start, end, options, output);
            }
        }

        return output;
    }

    /// <summary>
    /// Groups <see cref="SegmentPhasedBaf"/> rows into the segments Battenberg <c>determine_copynumber</c> fits — a new
    /// segment starts wherever <c>BAFseg</c> or the chromosome changes (<c>switchpoints</c>) — with each segment's logR
    /// = mean of the finite logR probes of the same chromosome within [first SNP position, last SNP position]
    /// (<c>mean(LogRvals[LogRpos &gt;= startpos &amp; LogRpos &lt;= endpos &amp; !is.infinite(.), 3], na.rm = T)</c>, 0 when
    /// none). The result feeds <see cref="FitSubclonalCopyNumberWithSnpTest"/>: segment BAF = <c>BAFseg</c> (the level
    /// <c>l</c> is mirrored there), SNP BAFs = <c>BAFphased</c> (<c>BAFke</c>), Start/End = min/max SNP position.
    /// </summary>
    /// <param name="segmentedSnps">Rows of <see cref="SegmentPhasedBaf"/> (BAFs in [0, 1]). Non-null.</param>
    /// <param name="logR">Raw logR probes. Non-null (may be empty; NaN/±∞ ignored).</param>
    /// <returns>One entry per Battenberg segment, in row order.</returns>
    /// <exception cref="ArgumentNullException">an argument or a chromosome is null.</exception>
    /// <exception cref="ArgumentException">a BAF or phased BAF outside [0, 1].</exception>
    public static IReadOnlyList<SubclonalSegmentSnpBafs> BuildBattenbergSegments(
        IReadOnlyList<PhasedBafSegmentedSnp> segmentedSnps,
        IReadOnlyList<LogRProbe> logR)
    {
        ArgumentNullException.ThrowIfNull(segmentedSnps);
        ArgumentNullException.ThrowIfNull(logR);
        var probes = new Dictionary<string, List<(long Pos, double LogR)>>(StringComparer.Ordinal);
        foreach (LogRProbe p in logR)
        {
            if (p.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(logR), "Every logR probe needs a chromosome.");
            }

            if (!double.IsFinite(p.LogR))
            {
                continue; // !is.infinite(.) and na.rm = T
            }

            if (!probes.TryGetValue(p.Chromosome, out var list))
            {
                list = new List<(long, double)>();
                probes.Add(p.Chromosome, list);
            }

            list.Add((p.Position, p.LogR));
        }

        var result = new List<SubclonalSegmentSnpBafs>();
        int n = segmentedSnps.Count;
        int from = 0;
        for (int i = 0; i < n; i++)
        {
            PhasedBafSegmentedSnp row = segmentedSnps[i];
            if (row.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(segmentedSnps), "Every row needs a chromosome.");
            }

            if (!(row.BafSegment >= 0.0 && row.BafSegment <= 1.0) || !(row.BafPhased >= 0.0 && row.BafPhased <= 1.0))
            {
                throw new ArgumentException("Every BAFseg and BAFphased must lie in [0, 1].", nameof(segmentedSnps));
            }

            bool last = i == n - 1
                        || segmentedSnps[i + 1].BafSegment != row.BafSegment
                        || segmentedSnps[i + 1].Chromosome != row.Chromosome;
            if (!last)
            {
                continue;
            }

            long start = long.MaxValue, end = long.MinValue;
            var phased = new double[i - from + 1];
            for (int j = from; j <= i; j++)
            {
                start = Math.Min(start, segmentedSnps[j].Position);
                end = Math.Max(end, segmentedSnps[j].Position);
                phased[j - from] = segmentedSnps[j].BafPhased;
            }

            var inRange = new List<double>();
            if (probes.TryGetValue(row.Chromosome, out var chrProbes))
            {
                foreach ((long pos, double value) in chrProbes)
                {
                    if (pos >= start && pos <= end)
                    {
                        inRange.Add(value);
                    }
                }
            }

            double segLogR = inRange.Count == 0 ? 0.0 : RMean(inRange); // is.na(LogR) ⇒ 0
            result.Add(new SubclonalSegmentSnpBafs(
                new AlleleSpecificSegmentSummary(row.Chromosome, start, end, segLogR, row.BafSegment, phased.Length),
                phased));
            from = i + 1;
        }

        return result;
    }

    private static void ValidateBattenbergSegmentationOptions(BattenbergPhasedSegmentationOptions options)
    {
        BattenbergPhasedSegmentationOptions o = options;
        if (!(o.Gamma >= 0.0) || double.IsPositiveInfinity(o.Gamma))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.Gamma, "gamma must be finite and ≥ 0.");
        }

        if (!(o.PhaseGamma >= 0.0) || double.IsPositiveInfinity(o.PhaseGamma))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.PhaseGamma, "phasegamma must be finite and ≥ 0.");
        }

        if (o.Kmin < 1 || o.Kmin > BattenbergMaxKmin)
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.Kmin, "kmin must lie in [1, 14].");
        }

        if (o.PhaseKmin < 1 || o.PhaseKmin > BattenbergMaxKmin)
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.PhaseKmin, "phasekmin must lie in [1, 14].");
        }

        if (!Enum.IsDefined(o.SegmentBafOption))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.SegmentBafOption, "Unknown calc_seg_baf_option.");
        }

        ArgumentNullException.ThrowIfNull(o.PriorBreakpoints);
        foreach (BattenbergPriorBreakpoint b in o.PriorBreakpoints)
        {
            if (b.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(options), "Every prior breakpoint needs a chromosome.");
            }
        }
    }

    /// <summary>
    /// Battenberg <c>bkps_to_presegment_breakpoints(…, addin_bigholes = T)</c> on one chromosome; R <c>NA</c> = null
    /// (a presegment with an <c>NA</c> bound selects no SNP).
    /// </summary>
    private static List<(long? Start, long? End)> BattenbergPresegments(long[] pos, List<long> bkps)
    {
        int n = pos.Length;
        long? At(int oneBased) => oneBased >= 1 && oneBased <= n ? pos[oneBased - 1] : null; // R x[i], NA beyond
        var result = new List<(long?, long?)>();

        // addin_bigholes over SNP positions sel (1-based indices into pos)
        long? AddBigHoles(List<int> sel, long? startpos)
        {
            for (int k = 0; k + 1 < sel.Count; k++)
            {
                if (pos[sel[k + 1] - 1] - pos[sel[k] - 1] >= BattenbergMaxSnpDistance)
                {
                    result.Add((startpos, pos[sel[k] - 1]));
                    startpos = pos[sel[k + 1] - 1];
                }
            }

            return startpos;
        }

        if (bkps.Count > 0)
        {
            long? startpos;
            int startFromSv;
            long? first = At(1);
            if (first.HasValue && first.Value < bkps[0])
            {
                startpos = first;
                startFromSv = 1;
            }
            else
            {
                // R: NA < bkp is NA ⇒ if() would stop; an empty chromosome never reaches here with data.
                startpos = bkps[0];
                startFromSv = 2;
            }

            foreach (int idx in RColonRange(startFromSv, bkps.Count))
            {
                long? sv = idx >= 1 && idx <= bkps.Count ? bkps[idx - 1] : null;
                var selected = new List<int>();
                if (startpos.HasValue && sv.HasValue)
                {
                    for (int j = 1; j <= n; j++)
                    {
                        if (pos[j - 1] >= startpos.Value && pos[j - 1] <= sv.Value)
                        {
                            selected.Add(j);
                        }
                    }
                }

                if (selected.Count > 0)
                {
                    startpos = AddBigHoles(selected, startpos);
                    int endIndex = selected[^1];
                    result.Add((startpos, pos[endIndex - 1]));
                    startpos = At(endIndex + 1);
                }
            }

            if (n > 0 && pos[n - 1] > bkps[^1])
            {
                result.Add((startpos, pos[n - 1]));
            }
        }
        else
        {
            long? startpos = At(1);
            var all = new List<int>(n);
            for (int j = 1; j <= n; j++)
            {
                all.Add(j);
            }

            startpos = AddBigHoles(all, startpos);
            result.Add((startpos, At(n)));
        }

        return result;
    }

    /// <summary>R <c>a:b</c> for integers (descending when a &gt; b).</summary>
    private static IEnumerable<int> RColonRange(int a, int b)
    {
        int step = a <= b ? 1 : -1;
        int count = Math.Abs(b - a) + 1;
        for (int k = 0; k < count; k++)
        {
            yield return a + k * step;
        }
    }

    /// <summary>Battenberg <c>run_pcf</c> on one presegment; appends its rows to <paramref name="output"/>.</summary>
    private static void BattenbergRunPcf(
        string chr, long[] pos, double[] bafAll, long? start, long? end,
        BattenbergPhasedSegmentationOptions o, List<PhasedBafSegmentedSnp> output)
    {
        var idx = new List<int>();
        if (start.HasValue && end.HasValue)
        {
            for (int j = 0; j < pos.Length; j++)
            {
                if (pos[j] >= start.Value && pos[j] <= end.Value)
                {
                    idx.Add(j);
                }
            }
        }

        int n = idx.Count;
        if (n == 0)
        {
            return;
        }

        var baf = new double[n];
        var folded = new double[n];
        for (int k = 0; k < n; k++)
        {
            baf[k] = bafAll[idx[k]];
            folded[k] = baf[k] < BalancedBaf ? baf[k] : 1.0 - baf[k];
        }

        double sdev = GetMad(folded, 25);
        if (double.IsNaN(sdev))
        {
            sdev = 0.0;
        }

        if (sdev < BattenbergMinBafSd)
        {
            sdev = BattenbergMinBafSd;
        }

        double[] bafSegm = n < BattenbergMinPcfSnps
            ? Filled(n, RMean(baf))
            : BattenbergSelectFastPcf(baf, o.PhaseKmin, o.PhaseGamma * sdev);

        var phased = new double[n];
        for (int k = 0; k < n; k++)
        {
            phased[k] = bafSegm[k] > BalancedBaf ? baf[k] : 1.0 - baf[k];
        }

        double[] phSeg = n < BattenbergMinPcfSnps || o.NoSegmentation
            ? Filled(n, RMean(phased))
            : BattenbergSelectFastPcf(phased, o.Kmin, o.Gamma * sdev);

        if (o.SegmentBafOption != BattenbergSegmentBafOption.Mean)
        {
            double[] median = BattenbergAdjustSegmValues(phased, phSeg);
            for (int k = 0; k < n; k++)
            {
                bool keepMean = o.SegmentBafOption == BattenbergSegmentBafOption.MedianUnlessExtreme
                                && (median[k] == 0.0 || median[k] == 1.0);
                phSeg[k] = keepMean ? phSeg[k] : median[k];
            }
        }

        for (int k = 0; k < n; k++)
        {
            output.Add(new PhasedBafSegmentedSnp(chr, pos[idx[k]], baf[k], phased[k], phSeg[k]));
        }
    }

    /// <summary>Battenberg <c>adjustSegmValues</c>: each run of equal <c>BAFseg</c> (R <c>rle</c>) → median of its BAFphased.</summary>
    private static double[] BattenbergAdjustSegmValues(double[] phased, double[] seg)
    {
        int n = seg.Length;
        var result = new double[n];
        int from = 0;
        for (int i = 0; i < n; i++)
        {
            if (i < n - 1 && seg[i + 1] == seg[i])
            {
                continue;
            }

            double median = StatisticsHelper.Median(new ArraySegment<double>(phased, from, i - from + 1));
            for (int k = from; k <= i; k++)
            {
                result[k] = median;
            }

            from = i + 1;
        }

        return result;
    }

    private static double[] Filled(int n, double value)
    {
        var a = new double[n];
        Array.Fill(a, value);
        return a;
    }

    /// <summary>R <c>mean(x)</c>: sum/n refined by the mean residual (R summary.c; R accumulates in long double).</summary>
    private static double RMean(IReadOnlyList<double> x)
    {
        int n = x.Count;
        double s = 0.0;
        for (int i = 0; i < n; i++)
        {
            s += x[i];
        }

        s /= n;
        if (double.IsFinite(s))
        {
            double t = 0.0;
            for (int i = 0; i < n; i++)
            {
                t += x[i] - s;
            }

            s += t / n;
        }

        return s;
    }

    /// <summary>
    /// Battenberg/copynumber <c>selectFastPcf(x, kmin, gamma, yest = T)$yhat</c> (R/fastPCF.R): <c>runFastPcf</c> with
    /// filter fractions (0.15, 0.15) below 1000 values, (0.12, 0.05) below 15 000, else <c>runPcfSubset</c>.
    /// </summary>
    private static double[] BattenbergSelectFastPcf(double[] x, int kmin, double gamma)
    {
        int n = x.Length;
        if (n < 1000)
        {
            return BattenbergRunFastPcf(x, kmin, gamma, 0.15, 0.15);
        }

        return n < 15000
            ? BattenbergRunFastPcf(x, kmin, gamma, 0.12, 0.05)
            : BattenbergRunPcfSubset(x, kmin, gamma, 0.12, 0.05);
    }

    private static double[] BattenbergRunFastPcf(double[] x, int kmin, double gamma, double frac1, double frac2)
    {
        bool[] mark = BattenbergFilterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9);
        mark[^1] = true;
        (int[] nr, double[] sum, double[] sq) = BattenbergCompact(x, x.Length, mark);
        return BattenbergPottsCompact(kmin, gamma, nr, sum, sq);
    }

    /// <summary>copynumber <c>runPcfSubset</c>: Potts marks refined over 5000-value windows advancing by 4000.</summary>
    private static double[] BattenbergRunPcfSubset(double[] x, int kmin, double gamma, double frac1, double frac2)
    {
        const int subSize = 5000;
        int antGen = x.Length;
        bool[] mark = BattenbergFilterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9);
        var markInit = new bool[subSize];
        Array.Copy(mark, markInit, subSize - 1);
        markInit[subSize - 1] = true;
        (int[] nr, double[] sum, double[] sq) = BattenbergCompact(x, subSize, markInit);
        var mark2 = new bool[antGen];
        Array.Copy(BattenbergMarkWithPotts(kmin, gamma, nr, sum, sq, subSize), mark2, subSize);
        mark2[(4 * subSize / 5) - 1] = true;
        int start = (4 * subSize / 5) + 1; // 1-based
        while (start + subSize < antGen)
        {
            int slutt = start + subSize - 1;
            var markSub = new bool[slutt];
            Array.Copy(mark2, markSub, start - 1);
            Array.Copy(mark, start - 1, markSub, start - 1, slutt - start + 1);
            markSub[slutt - 1] = true;
            (nr, sum, sq) = BattenbergCompact(x, slutt, markSub);
            Array.Copy(BattenbergMarkWithPotts(kmin, gamma, nr, sum, sq, slutt), mark2, slutt);
            start += 4 * subSize / 5;
            mark2[start - 2] = true;
        }

        var finalMark = new bool[antGen];
        Array.Copy(mark2, finalMark, start - 1);
        Array.Copy(mark, start - 1, finalMark, start - 1, antGen - start + 1);
        (nr, sum, sq) = BattenbergCompact(x, antGen, finalMark);
        return BattenbergPottsCompact(kmin, gamma, nr, sum, sq);
    }

    /// <summary>
    /// copynumber <c>compact(y[1:length], mark)</c>: counts, sums and sums of squares between marked positions, from the
    /// running <c>cumsum</c>s (differences of the cumulative sums at the marks, as R does).
    /// </summary>
    private static (int[] Nr, double[] Sum, double[] Sq) BattenbergCompact(double[] y, int length, bool[] mark)
    {
        var nr = new List<int>();
        var sum = new List<double>();
        var sq = new List<double>();
        double cy = 0.0, cy2 = 0.0, lowCy = 0.0, lowCy2 = 0.0;
        int low = 0;
        for (int i = 0; i < length; i++)
        {
            cy += y[i];
            cy2 += y[i] * y[i];
            if (mark[i])
            {
                nr.Add(i + 1 - low);
                sum.Add(cy - lowCy);
                sq.Add(cy2 - lowCy2);
                low = i + 1;
                lowCy = cy;
                lowCy2 = cy2;
            }
        }

        return (nr.ToArray(), sum.ToArray(), sq.ToArray());
    }

    /// <summary>
    /// copynumber <c>PottsCompact(kmin, gamma, nr, res, sq, yest = T)$yhat</c>: exact Potts filtering on the compacted
    /// array (1-based arrays below mirror the R code line by line), then <c>findEst</c>.
    /// </summary>
    private static double[] BattenbergPottsCompact(int kmin, double gamma, int[] nr, double[] res, double[] sq)
    {
        int bigN = nr.Length;
        long total = 0;
        double totalSum = 0.0;
        for (int i = 0; i < bigN; i++)
        {
            total += nr[i];
            totalSum += res[i];
        }

        if (total < 2 * kmin)
        {
            return Filled((int)total, totalSum / total); // R returns the scalar estimate
        }

        int[] bestSplit = BattenbergPottsSplits(kmin, gamma, nr, res, sq, markSub: null);
        return BattenbergFindEst(bestSplit, bigN, nr, res);
    }

    /// <summary>
    /// The Potts recursion shared by copynumber <c>PottsCompact</c> and <c>markWithPotts</c>; returns 1-based
    /// <c>bestSplit</c> and, when <paramref name="markSub"/> is given, sets <c>markSub[Pos − 1]</c> as markWithPotts does.
    /// </summary>
    private static int[] BattenbergPottsSplits(int kmin, double gamma, int[] nr, double[] res, double[] sq, bool[]? markSub)
    {
        int bigN = nr.Length;
        var ant = new double[bigN + 1];
        var sum = new double[bigN + 1];
        var kvad = new double[bigN + 1];
        var cost = new double[bigN + 1];
        var bestCost = new double[bigN + 1];
        var bestSplit = new int[bigN + 1];
        double initAnt = nr[0], initSum = res[0], initKvad = sq[0];
        double initAve = initSum / initAnt;
        bestCost[1] = initKvad - initSum * initAve;
        int k = 2;
        double cum = bigN > 1 ? nr[0] + nr[1] : nr[0];
        while (cum < 2 * kmin && k < bigN)
        {
            for (int j = 2; j <= k; j++)
            {
                ant[j] += nr[k - 1];
                sum[j] += res[k - 1];
                kvad[j] += sq[k - 1];
            }

            double s2 = initSum + sum[2];
            bestCost[k] = (initKvad + kvad[2]) - s2 * s2 / (initAnt + ant[2]);
            k++;
            cum += nr[k - 1];
        }

        for (int n = k; n <= bigN; n++)
        {
            for (int j = 2; j <= n; j++)
            {
                ant[j] += nr[n - 1];
                sum[j] += res[n - 1];
                kvad[j] += sq[n - 1];
            }

            int limit = n;
            while (limit > 2 && ant[limit] < kmin)
            {
                limit--;
            }

            int pos = -1;
            double best = double.NaN;
            for (int j = 2; j <= limit; j++)
            {
                cost[j] = bestCost[j - 1] + kvad[j] - sum[j] * sum[j] / ant[j];
                if (!double.IsNaN(cost[j]) && (pos < 0 || cost[j] < best))
                {
                    pos = j;
                    best = cost[j];
                }
            }

            double c = best + gamma;
            double st = sum[2] + initSum;
            double totCost = (kvad[2] + initKvad) - st * st / (ant[2] + initAnt);
            if (totCost < c)
            {
                pos = 1;
                c = totCost;
            }

            bestCost[n] = c;
            bestSplit[n] = pos - 1;
            if (markSub is not null && pos - 1 >= 1)
            {
                markSub[pos - 2] = true; // R markSub[Pos − 1] (1-based)
            }
        }

        return bestSplit;
    }

    /// <summary>copynumber <c>findEst(bestSplit, N, Nr, Sum, yest = T)$yhat</c>.</summary>
    private static double[] BattenbergFindEst(int[] bestSplit, int bigN, int[] nr, double[] sum)
    {
        var lengths = new List<int>();
        int n = bigN;
        while (n > 0)
        {
            lengths.Add(n - bestSplit[n]);
            n = bestSplit[n];
        }

        lengths.Reverse();
        var yhat = new List<double>();
        int start = 0; // 0-based into the compact arrays
        foreach (int len in lengths)
        {
            int lengdeOrig = 0;
            double s = 0.0;
            for (int i = start; i < start + len; i++)
            {
                lengdeOrig += nr[i];
                s += sum[i];
            }

            double verdi = s / lengdeOrig;
            for (int i = 0; i < lengdeOrig; i++)
            {
                yhat.Add(verdi);
            }

            start += len;
        }

        return yhat.ToArray();
    }

    /// <summary>copynumber <c>markWithPotts</c> + <c>findMarks</c>: Potts split points mapped back to original indices.</summary>
    private static bool[] BattenbergMarkWithPotts(int kmin, double gamma, int[] nr, double[] res, double[] sq, int subSize)
    {
        var markSub = new bool[nr.Length];
        BattenbergPottsSplits(kmin, gamma, nr, res, sq, markSub);
        var mark = new bool[subSize];
        int orig = 0;
        for (int i = 0; i < markSub.Length; i++)
        {
            orig += nr[i];
            if (markSub[i])
            {
                mark[orig - 1] = true; // findMarks: mark[startOrig − 1] at the end of compact cell i
            }
        }

        return mark;
    }

    /// <summary>
    /// copynumber <c>filterMarkS4(x, kmin, L, L2, frac1, frac2, frac3, thres)</c>: candidate breakpoints from two
    /// high-pass filters (widths 6L, 6L2) and a kmin-segment filter, with type-7 quantile limits
    /// (<see cref="StatisticsHelper.SampleQuantileType7"/>). Arrays are 1-based (index 0 unused) as in R.
    /// </summary>
    private static bool[] BattenbergFilterMarkS4(
        double[] x, int kmin, int bigL, int bigL2, double frac1, double frac2, double frac3, double thres)
    {
        int n = x.Length;
        var xc = new double[n + 1]; // R xc = c(0, cumsum(x)); xc[j] (1-based R) = xc[j − 1] here
        for (int i = 0; i < n; i++)
        {
            xc[i + 1] = xc[i] + x[i];
        }

        double Xc(int rIndex) => xc[rIndex - 1];

        // cost1 (length n, 1-based)
        var cost1 = new double[n + 1];
        for (int t = 1; t <= n - 6 * bigL + 1; t++)
        {
            cost1[3 * bigL - 1 + t] = Math.Abs(
                4 * Xc(t + 3 * bigL) - Xc(t) - Xc(t + bigL) - Xc(t + 5 * bigL) - Xc(t + 6 * bigL));
        }

        double[] test = SevenMax(cost1, n);
        var cost1B = new List<double>();
        for (int j = 1; j <= n; j++)
        {
            if (cost1[j] >= thres * test[j])
            {
                cost1B.Add(cost1[j]);
            }
        }

        double frac1B = Math.Min(0.8, frac1 * n / cost1B.Count);
        double limit = StatisticsHelper.SampleQuantileType7(cost1B, 1 - frac1B);
        var mark = new bool[n + 1];
        for (int j = 1; j <= n; j++)
        {
            mark[j] = cost1[j] > limit && cost1[j] > 0.9 * test[j];
        }

        int m2 = n - 6 * bigL2 + 1;
        var cost2 = new double[m2];
        for (int t = 1; t <= m2; t++)
        {
            cost2[t - 1] = Math.Abs(
                4 * Xc(t + 3 * bigL2) - Xc(t) - Xc(t + bigL2) - Xc(t + 5 * bigL2) - Xc(t + 6 * bigL2));
        }

        double limit2 = StatisticsHelper.SampleQuantileType7(cost2, 1 - frac2);
        var mark2 = new bool[n + 1];
        for (int t = 1; t <= m2; t++)
        {
            mark2[3 * bigL2 - 1 + t] = cost2[t - 1] > limit2;
        }

        if (3 * bigL > kmin)
        {
            SetRange(mark, kmin, 3 * bigL - 1, true);
            SetRange(mark, n - 3 * bigL + 1, n - kmin, true);
        }
        else
        {
            SetRange(mark, kmin, kmin, true);
            SetRange(mark, n - kmin, n - kmin, true);
        }

        if (kmin > 1)
        {
            int m = n - 3 * kmin + 1;
            var shortAb = new double[m + 1];
            for (int t = 1; t <= m; t++)
            {
                shortAb[t] = Math.Abs(3 * (Xc(t + 2 * kmin) - Xc(t + kmin)) - (Xc(t + 3 * kmin) - Xc(t)));
            }

            double[] test2 = SevenMax(shortAb, m);
            var cost1C = new List<double>();
            for (int t = 1; t <= m; t++)
            {
                if (shortAb[t] >= thres * test2[t])
                {
                    cost1C.Add(shortAb[t]);
                }
            }

            double frac1C = Math.Min(0.8, frac3 * m / cost1C.Count);
            double limit3 = StatisticsHelper.SampleQuantileType7(cost1C, 1 - frac1C);
            for (int t = 1; t <= m; t++)
            {
                if (shortAb[t] > limit3 && shortAb[t] > thres * test2[t])
                {
                    mark[kmin - 1 + t] = true;     // markH2
                    mark[2 * kmin - 1 + t] = true; // markH3
                }
            }
        }

        for (int j = 1; j <= n; j++)
        {
            mark[j] |= mark2[j];
        }

        if (3 * bigL > kmin)
        {
            SetRange(mark, 1, kmin - 1, false);
            SetRange(mark, kmin, 3 * bigL - 1, true);
            SetRange(mark, n - 3 * bigL + 1, n - kmin, true);
            SetRange(mark, n - kmin + 1, n - 1, false);
            mark[n] = true;
        }
        else
        {
            SetRange(mark, 1, kmin - 1, false);
            SetRange(mark, n - kmin + 1, n - 1, false);
            mark[n] = true;
            mark[kmin] = true;
            mark[n - kmin] = true;
        }

        var result = new bool[n];
        Array.Copy(mark, 1, result, 0, n);
        return result;
    }

    /// <summary>
    /// filterMarkS4's <c>c(rep(0, 3), pmax(v[i..i+6]), rep(0, 3))</c> over a 1-based vector v of length len.
    /// </summary>
    private static double[] SevenMax(double[] v, int len)
    {
        var test = new double[len + 1];
        for (int i = 1; i <= len - 6; i++)
        {
            double mx = v[i];
            for (int d = 1; d <= 6; d++)
            {
                mx = Math.Max(mx, v[i + d]);
            }

            test[3 + i] = mx;
        }

        return test;
    }

    /// <summary>R <c>mark[a:b] &lt;- value</c> on a 1-based array (descending when a &gt; b; index 0 ignored).</summary>
    private static void SetRange(bool[] mark, int a, int b, bool value)
    {
        foreach (int i in RColonRange(a, b))
        {
            if (i >= 1)
            {
                mark[i] = value;
            }
        }
    }

    #endregion

    #region CNAqc peak-based purity QC (ONCO-PURITY-001)

    /// <summary>
    /// A somatic mutation for CNAqc's peak-based purity QC: its VAF and the clonal allele-specific copy-number state
    /// (Major:minor) of the segment it maps to (CNAqc <c>karyotype = "Major:minor"</c>).
    /// </summary>
    /// <param name="Vaf">Variant allele frequency NV/DP in [0, 1].</param>
    /// <param name="MajorCopyNumber">Major allele copy number (≥ 0).</param>
    /// <param name="MinorCopyNumber">Minor allele copy number (≥ 0).</param>
    public readonly record struct PurityPeakMutation(double Vaf, int MajorCopyNumber, int MinorCopyNumber);

    /// <summary>How detected VAF peaks are assigned to the expected clonal peaks of a karyotype.</summary>
    public enum PurityPeakMatchingStrategy
    {
        /// <summary>
        /// Each expected peak takes the nearest non-discarded data peak (first on ties) — CNAqc 1.1.5
        /// <c>analyze_peaks_common</c>, the only strategy its <c>analyze_peaks</c> runs.
        /// </summary>
        Closest,

        /// <summary>
        /// Expected peaks (descending) are paired with the highest-VAF non-discarded data peaks (descending), padding
        /// with the rightmost peak when there are fewer data peaks — CNAqc's legacy <c>peak_detector</c>
        /// (<c>matching_strategy = "rightmost"</c>; still documented by <c>analyze_peaks</c> but no longer forwarded to
        /// <c>analyze_peaks_common</c> in 1.1.5).
        /// </summary>
        Rightmost,
    }

    /// <summary>Origin of a data peak.</summary>
    public enum PurityPeakSource
    {
        /// <summary>Gaussian-KDE maximum found by <c>peakPick</c> (CNAqc <c>from = "KDE"</c>).</summary>
        Kde,

        /// <summary>A caller-supplied mixture-component mean (e.g. BMix Binomial <c>B.params</c>; CNAqc <c>from = "BMix"</c>).</summary>
        Mixture,
    }

    /// <summary>
    /// Parameters of <see cref="AnalyzePurityPeaks"/>; defaults are CNAqc 1.1.5 <c>analyze_peaks</c> defaults.
    /// </summary>
    public sealed record PurityPeakOptions
    {
        /// <summary>CNAqc default simple clonal karyotypes <c>c('1:0', '1:1', '2:0', '2:1', '2:2')</c> as (Major, minor).</summary>
        public static IReadOnlyList<(int Major, int Minor)> SimpleClonalKaryotypes { get; } =
            new[] { (1, 0), (1, 1), (2, 0), (2, 1), (2, 2) };

        /// <summary>Default options (CNAqc defaults).</summary>
        public static PurityPeakOptions Default { get; } = new();

        /// <summary>Karyotypes to QC (<c>karyotypes</c>); must be a subset of <see cref="SimpleClonalKaryotypes"/>.</summary>
        public IReadOnlyList<(int Major, int Minor)> Karyotypes { get; init; } = SimpleClonalKaryotypes;

        /// <summary>Minimum share n_k / N of all mutations (<c>min_karyotype_size</c>, default 0), in [0, 1).</summary>
        public double MinKaryotypeSize { get; init; }

        /// <summary>Minimum mutations per karyotype (<c>min_absolute_karyotype_mutations</c>, default 100; n ≥ this).</summary>
        public int MinAbsoluteKaryotypeMutations { get; init; } = 100;

        /// <summary>Purity error ε (<c>purity_error</c>, default 0.05), in (0, 1); sets the VAF bands δ.</summary>
        public double PurityError { get; init; } = 0.05;

        /// <summary>VAF tolerance around a data peak for band overlap (<c>VAF_tolerance</c>, default 0.015), ≥ 0.</summary>
        public double VafTolerance { get; init; } = 0.015;

        /// <summary>KDE bandwidth multiplier (<c>kernel_adjust</c>, default 1), &gt; 0.</summary>
        public double KernelAdjust { get; init; } = 1.0;

        /// <summary>Peak-matching rule (default <see cref="PurityPeakMatchingStrategy.Closest"/>).</summary>
        public PurityPeakMatchingStrategy MatchingStrategy { get; init; } = PurityPeakMatchingStrategy.Closest;

        /// <summary>Only mutations with VAF &gt; this are analysed (<c>min_VAF</c>, default 0).</summary>
        public double MinVaf { get; init; }

        /// <summary>
        /// Optional mixture-component means per karyotype (CNAqc adds BMix Binomial-mixture peaks,
        /// <c>bmixfit(K.Binomials = 1:4)</c>, to the KDE peaks). Each mean is snapped to the nearest KDE grid point.
        /// BMix is stochastic (k-means starts, jittered initial means) and is not ported; null = KDE peaks only.
        /// </summary>
        public IReadOnlyDictionary<(int Major, int Minor), IReadOnlyList<double>>? MixturePeaks { get; init; }

        /// <summary>true reproduces R ≤ 4.3 <c>density</c> values (<c>old.coords = TRUE</c>); default R ≥ 4.4.</summary>
        public bool LegacyDensityCoordinates { get; init; }
    }

    /// <summary>A VAF peak detected in a karyotype's data (CNAqc <c>xy_peaks</c> row).</summary>
    /// <param name="X">Peak VAF (KDE peaks rounded to 2 decimals; mixture peaks at the snapped KDE grid point).</param>
    /// <param name="Y">Density at the peak (KDE peaks rounded to 2 decimals).</param>
    /// <param name="CountsPerBin">Mutations in the 0.01-wide VAF histogram bin <c>round(100·X)</c> (null if outside 1..100).</param>
    /// <param name="Discarded">KDE peak with Y ≤ max(Y)/20 (never used for matching); mixture peaks are never discarded.</param>
    /// <param name="Source">KDE or mixture.</param>
    public readonly record struct PurityDataPeak(double X, double Y, int? CountsPerBin, bool Discarded, PurityPeakSource Source);

    /// <summary>One expected clonal peak and its matched data peak (CNAqc <c>peaks_analysis$matches</c> row).</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="Multiplicity">Mutation multiplicity m (1 or Major).</param>
    /// <param name="ExpectedPeak">m·π / (2(1−π) + π·(Major + minor)).</param>
    /// <param name="DeltaVaf">VAF band half-width δ = 2·m·ε / (2 + π·(ploidy − 2))².</param>
    /// <param name="MatchedPeak">The data peak assigned to this expectation.</param>
    /// <param name="OffsetVaf">ExpectedPeak − MatchedPeak.X.</param>
    /// <param name="Offset">Purity-space offset 2·m·OffsetVaf / (m + X·(2 − ploidy))².</param>
    /// <param name="Weight">Karyotype weight n_k / Σ n over analysed karyotypes.</param>
    /// <param name="Matched">[X ± VafTolerance] overlaps [ExpectedPeak ± δ].</param>
    public readonly record struct PurityPeakMatch(
        int MajorCopyNumber,
        int MinorCopyNumber,
        int Multiplicity,
        double ExpectedPeak,
        double DeltaVaf,
        PurityDataPeak MatchedPeak,
        double OffsetVaf,
        double Offset,
        double Weight,
        bool Matched);

    /// <summary>Per-karyotype result of <see cref="AnalyzePurityPeaks"/>.</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="MutationCount">Mutations pooled for this karyotype (VAF &gt; MinVaf).</param>
    /// <param name="Weight">n_k / Σ n over analysed karyotypes.</param>
    /// <param name="Score">Σ Weight·Offset over this karyotype's expected peaks.</param>
    /// <param name="Pass">QC of the expected peak whose matched data peak holds the most mutations (CountsPerBin).</param>
    /// <param name="Density">The Gaussian KDE of the karyotype's VAFs.</param>
    /// <param name="Peaks">All data peaks (KDE, then mixture).</param>
    /// <param name="Matches">Expected-peak matches (Closest: m = 1 then Major; Rightmost: descending expected VAF).</param>
    public sealed record PurityPeakKaryotype(
        int MajorCopyNumber,
        int MinorCopyNumber,
        int MutationCount,
        double Weight,
        double Score,
        bool Pass,
        KernelDensityEstimate Density,
        IReadOnlyList<PurityDataPeak> Peaks,
        IReadOnlyList<PurityPeakMatch> Matches);

    /// <summary>Result of <see cref="AnalyzePurityPeaks"/> (CNAqc <c>x$peaks_analysis</c> for simple clonal CNAs).</summary>
    /// <param name="Purity">The purity being QC'd.</param>
    /// <param name="Score">Σ Weight·Offset over all matches (CNAqc λ, printed as "Purity correction"); NaN if nothing analysed.</param>
    /// <param name="Pass">Sample QC: the PASS/FAIL class with the larger summed match weight (ties → FAIL); null if no karyotype passed the filters.</param>
    /// <param name="Karyotypes">Analysed karyotypes in (Major, minor) order.</param>
    /// <param name="Matches">All matches, karyotype by karyotype.</param>
    public sealed record PurityPeakAnalysis(
        double Purity,
        double Score,
        bool? Pass,
        IReadOnlyList<PurityPeakKaryotype> Karyotypes,
        IReadOnlyList<PurityPeakMatch> Matches);

    /// <summary>
    /// CNAqc peak-based QC of a tumour purity estimate for simple clonal karyotypes (CNAqc 1.1.5
    /// <c>analyze_peaks</c> → <c>analyze_peaks_common</c>; Antonello et al. 2024, <i>Genome Biology</i> 25:38).
    /// For each karyotype K = Major:minor in <see cref="PurityPeakOptions.Karyotypes"/> with n_K ≥ MinAbsoluteKaryotypeMutations
    /// and n_K / N ≥ MinKaryotypeSize (N = all mutations with VAF &gt; MinVaf, any karyotype):
    /// <list type="number">
    /// <item>pool its VAFs; Gaussian KDE (R <c>density</c>, bw.nrd0 × KernelAdjust); peaks = union of
    /// <c>peakPick::peakpick(neighlim = 1..5)</c> maxima, (x, y) rounded to 2 decimals, distinct x, clamped to [0, 1];
    /// a KDE peak is discarded if y ≤ max(y)/20; optional mixture peaks are appended;</item>
    /// <item>expected peaks for m ∈ {1, Major}: v_m = m·π / (2(1−π) + π·(Major + minor)), bands
    /// δ_m = 2·m·ε / (2 + π·(ploidy − 2))²;</item>
    /// <item>match each expected peak to a non-discarded data peak (<see cref="PurityPeakMatchingStrategy"/>);
    /// offset = 2·m·(v_m − x) / (m + x·(2 − ploidy))² (purity units), weight = n_K / Σ n_analysed;
    /// matched ⇔ [x ± VafTolerance] ∩ [v_m ± δ_m] ≠ ∅;</item>
    /// <item>karyotype PASS ⇔ the expected peak whose matched data peak has the largest histogram count is matched;
    /// sample PASS ⇔ PASS rows carry more total weight than FAIL rows; score λ = Σ weight·offset.</item>
    /// </list>
    /// CNAqc does not threshold λ against ε (ε only sets the bands) and proposes no corrected purity: <c>print</c>
    /// reports λ as "Purity correction". <c>p_binsize_peaks</c> is accepted by CNAqc 1.1.5 but unused. KDE peak
    /// detection is deterministic and ported exactly; CNAqc's BMix mixture peaks (stochastic) can be supplied via
    /// <see cref="PurityPeakOptions.MixturePeaks"/>. Bootstrap (<c>n_bootstrap</c> &gt; 1), complex and subclonal
    /// karyotypes (<c>analyze_peaks_general</c>/<c>_subclonal</c>) are not ported.
    /// </summary>
    /// <param name="mutations">Mutations with VAF and karyotype (all karyotypes; non-simple ones count towards N).</param>
    /// <param name="purity">The purity π ∈ (0, 1] to QC.</param>
    /// <param name="options">CNAqc parameters (null = defaults).</param>
    /// <returns>Per-karyotype peaks, matches, scores and the QC verdict.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid purity, option, VAF or copy number.</exception>
    /// <exception cref="ArgumentException">A requested karyotype is not a simple clonal karyotype.</exception>
    /// <exception cref="InvalidOperationException">An analysed karyotype yields no usable VAF peak (CNAqc errors too).</exception>
    public static PurityPeakAnalysis AnalyzePurityPeaks(
        IEnumerable<PurityPeakMutation> mutations,
        double purity,
        PurityPeakOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(mutations);
        options ??= PurityPeakOptions.Default;
        ValidatePurityPeakArguments(purity, options);

        var byKaryotype = new SortedDictionary<(int Major, int Minor), List<double>>();
        int total = 0;
        foreach (PurityPeakMutation mutation in mutations)
        {
            if (double.IsNaN(mutation.Vaf) || mutation.Vaf < 0.0 || mutation.Vaf > 1.0)
                throw new ArgumentOutOfRangeException(nameof(mutations), mutation.Vaf, "VAF must be in [0, 1].");
            if (mutation.MajorCopyNumber < 0 || mutation.MinorCopyNumber < 0)
                throw new ArgumentOutOfRangeException(nameof(mutations), "Allele copy numbers must be non-negative.");
            if (!(mutation.Vaf > options.MinVaf)) continue; // CNAqc: filter(VAF > min_VAF)

            var key = (mutation.MajorCopyNumber, mutation.MinorCopyNumber);
            if (!byKaryotype.TryGetValue(key, out List<double>? vafs))
            {
                vafs = new List<double>();
                byKaryotype[key] = vafs;
            }

            vafs.Add(mutation.Vaf);
            total++;
        }

        var analysed = new List<(int Major, int Minor)>();
        int analysedTotal = 0;
        foreach (KeyValuePair<(int Major, int Minor), List<double>> entry in byKaryotype)
        {
            int n = entry.Value.Count;
            if (options.Karyotypes.Contains(entry.Key)
                && n >= options.MinAbsoluteKaryotypeMutations
                && (double)n / total >= options.MinKaryotypeSize)
            {
                analysed.Add(entry.Key);
                analysedTotal += n;
            }
        }

        if (analysed.Count == 0)
        {
            return new PurityPeakAnalysis(purity, double.NaN, null, Array.Empty<PurityPeakKaryotype>(), Array.Empty<PurityPeakMatch>());
        }

        var karyotypes = new List<PurityPeakKaryotype>(analysed.Count);
        var allMatches = new List<PurityPeakMatch>();
        foreach ((int major, int minor) in analysed)
        {
            List<double> vafs = byKaryotype[(major, minor)];
            double weight = (double)vafs.Count / analysedTotal;
            IReadOnlyList<double>? mixturePeaks = null;
            if (options.MixturePeaks is not null && options.MixturePeaks.TryGetValue((major, minor), out IReadOnlyList<double>? m))
                mixturePeaks = m;

            PurityPeakKaryotype result = AnalyzePurityPeaksKaryotype(major, minor, vafs, weight, purity, options, mixturePeaks);
            karyotypes.Add(result);
            allMatches.AddRange(result.Matches);
        }

        double score = 0.0, passWeight = 0.0, failWeight = 0.0;
        foreach (PurityPeakKaryotype k in karyotypes)
        {
            foreach (PurityPeakMatch match in k.Matches)
            {
                score += match.Weight * match.Offset;
                if (k.Pass) passWeight += match.Weight; else failWeight += match.Weight;
            }
        }

        // dplyr group_by(QC) orders "FAIL" < "PASS"; arrange(desc(prop)) is stable, so a tie resolves to FAIL.
        return new PurityPeakAnalysis(purity, score, passWeight > failWeight, karyotypes, allMatches);
    }

    private static void ValidatePurityPeakArguments(double purity, PurityPeakOptions options)
    {
        if (!(purity > 0.0 && purity <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        if (!(options.PurityError > 0.0 && options.PurityError < 1.0))
            throw new ArgumentOutOfRangeException(nameof(options), options.PurityError, "PurityError must be in (0, 1).");
        if (!(options.MinKaryotypeSize >= 0.0 && options.MinKaryotypeSize < 1.0))
            throw new ArgumentOutOfRangeException(nameof(options), options.MinKaryotypeSize, "MinKaryotypeSize must be in [0, 1).");
        if (options.MinAbsoluteKaryotypeMutations < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.MinAbsoluteKaryotypeMutations, "MinAbsoluteKaryotypeMutations must be ≥ 0.");
        if (!(options.VafTolerance >= 0.0) || double.IsPositiveInfinity(options.VafTolerance))
            throw new ArgumentOutOfRangeException(nameof(options), options.VafTolerance, "VafTolerance must be finite and ≥ 0.");
        if (!(options.KernelAdjust > 0.0) || double.IsPositiveInfinity(options.KernelAdjust))
            throw new ArgumentOutOfRangeException(nameof(options), options.KernelAdjust, "KernelAdjust must be finite and > 0.");
        if (double.IsNaN(options.MinVaf))
            throw new ArgumentOutOfRangeException(nameof(options), options.MinVaf, "MinVaf must not be NaN.");
        if (options.Karyotypes is null)
            throw new ArgumentException("Karyotypes must not be null.", nameof(options));
        foreach ((int Major, int Minor) k in options.Karyotypes)
        {
            if (!PurityPeakOptions.SimpleClonalKaryotypes.Contains(k))
                throw new ArgumentException(
                    $"Karyotype {k.Major}:{k.Minor} is not a simple clonal karyotype (1:0, 1:1, 2:0, 2:1, 2:2).", nameof(options));
        }
    }

    private static PurityPeakKaryotype AnalyzePurityPeaksKaryotype(
        int major, int minor, List<double> vafs, double weight, double purity, PurityPeakOptions options,
        IReadOnlyList<double>? mixturePeaks)
    {
        int ploidy = major + minor;
        KernelDensityEstimate density = StatisticsHelper.GaussianKernelDensity(
            vafs, options.KernelAdjust, legacyCoordinates: options.LegacyDensityCoordinates);
        int[] histogram = VafHistogram(vafs);
        List<PurityDataPeak> peaks = DetectKdePeaks(density, histogram);
        if (peaks.Count == 0)
        {
            // CNAqc simple_peak_detector: `if (indexes[1] == 0)` on an empty peak set is an R error.
            throw new InvalidOperationException($"Cannot find KDE peaks for karyotype {major}:{minor}.");
        }

        if (mixturePeaks is not null)
        {
            foreach (double p in mixturePeaks)
            {
                // CNAqc mixture_peak_detector: w_den = which.min(abs(den$x − p)); x, y from the grid (unrounded).
                int w = 0;
                for (int i = 1; i < density.X.Count; i++)
                {
                    if (Math.Abs(density.X[i] - p) < Math.Abs(density.X[w] - p)) w = i;
                }

                double x = density.X[w];
                peaks.Add(new PurityDataPeak(x, density.Y[w], HistogramCount(histogram, Math.Round(x * 100, MidpointRounding.ToEven)), false, PurityPeakSource.Mixture));
            }
        }

        var candidates = peaks.Where(p => !p.Discarded).ToList();
        if (candidates.Count == 0)
            throw new InvalidOperationException($"No non-discarded VAF peak for karyotype {major}:{minor}.");

        // expected_vaf_peak: multiplicities unique(c(1, Major)).
        var multiplicities = major == 1 ? new[] { 1 } : new[] { 1, major };
        var expected = multiplicities
            .Select(m => (Multiplicity: m, Peak: m * purity / MixtureCopiesPerCell(purity, ploidy)))
            .ToList();

        var assignment = new List<((int Multiplicity, double Peak) Expectation, PurityDataPeak Peak)>();
        if (options.MatchingStrategy == PurityPeakMatchingStrategy.Closest)
        {
            foreach (var e in expected)
            {
                PurityDataPeak best = candidates[0];
                foreach (PurityDataPeak c in candidates)
                {
                    if (Math.Abs(c.X - e.Peak) < Math.Abs(best.X - e.Peak)) best = c;
                }

                assignment.Add((e, best));
            }
        }
        else
        {
            var expectedDesc = expected.OrderByDescending(e => e.Peak).ToList(); // stable
            var peaksDesc = candidates.OrderByDescending(p => p.X).ToList();      // stable
            for (int i = 0; i < expectedDesc.Count; i++)
            {
                assignment.Add((expectedDesc[i], i < peaksDesc.Count ? peaksDesc[i] : peaksDesc[0]));
            }
        }

        var matches = new List<PurityPeakMatch>(assignment.Count);
        double karyotypeScore = 0.0;
        foreach (((int m, double peak), PurityDataPeak data) in assignment)
        {
            double band = 2 * m * options.PurityError;
            double spread = 2 + (purity * (ploidy - 2));
            double deltaVaf = band / (spread * spread);
            double offsetVaf = peak - data.X;
            double denom = m + (data.X * (2 - ploidy));
            double offset = 2 * m * offsetVaf / (denom * denom);
            bool matched = Math.Max(data.X - options.VafTolerance, peak - deltaVaf)
                           <= Math.Min(data.X + options.VafTolerance, peak + deltaVaf);
            matches.Add(new PurityPeakMatch(major, minor, m, peak, deltaVaf, data, offsetVaf, offset, weight, matched));
            karyotypeScore += weight * offset;
        }

        // QC per karyotype: arrange(desc(counts_per_bin)) (stable, NA last) → first row's `matched`.
        PurityPeakMatch lead = matches[0];
        foreach (PurityPeakMatch match in matches)
        {
            if ((match.MatchedPeak.CountsPerBin ?? int.MinValue) > (lead.MatchedPeak.CountsPerBin ?? int.MinValue)) lead = match;
        }

        return new PurityPeakKaryotype(major, minor, vafs.Count, weight, karyotypeScore, lead.Matched, density, peaks, matches);
    }

    // CNAqc simple_peak_detector: union of peakpick(neighlim = 1..5) on the density's y; arrange(x); round(x, 2),
    // round(y, 2); distinct(x) keeping the first; x ∈ (1, 1.01) → 1, x ∈ (−0.01, 0) → 0; keep 0 ≤ x ≤ 1; counts from
    // hist(breaks = seq(0, 1, 0.01)) at round(100·x) (first index 0 → 1); discarded ⇔ y ≤ max(y)·(1/20).
    private static List<PurityDataPeak> DetectKdePeaks(KernelDensityEstimate density, int[] histogram)
    {
        var picked = new SortedSet<int>();
        for (int neighlim = 1; neighlim <= 5; neighlim++)
        {
            bool[] mask = StatisticsHelper.PeakPick(density.Y, neighlim);
            for (int i = 0; i < mask.Length; i++)
            {
                if (mask[i]) picked.Add(i);
            }
        }

        var rows = new List<(double X, double Y)>();
        var seen = new HashSet<double>();
        foreach (int i in picked) // grid x increases with the index, so this is arrange(x)
        {
            double x = RRound(density.X[i], 2);
            if (!seen.Add(x)) continue;
            double y = RRound(density.Y[i], 2);
            if (x > 1 && x < 1.01) x = 1;
            else if (x < 0 && x > -0.01) x = 0;
            if (x <= 1 && x >= 0) rows.Add((x, y));
        }

        var peaks = new List<PurityDataPeak>(rows.Count);
        if (rows.Count == 0) return peaks;
        double maxY = rows.Max(r => r.Y);
        for (int k = 0; k < rows.Count; k++)
        {
            double index = Math.Round(rows[k].X * 100, MidpointRounding.ToEven);
            if (k == 0 && index == 0) index = 1;
            peaks.Add(new PurityDataPeak(rows[k].X, rows[k].Y, HistogramCount(histogram, index),
                rows[k].Y <= maxY * (1.0 / 20), PurityPeakSource.Kde));
        }

        return peaks;
    }

    // R hist(VAF, breaks = seq(0, 1, 0.01), plot = FALSE)$counts: right-closed bins, include.lowest, breaks fuzzed by
    // 1e-7·median(diff(breaks)) (first break down, the others up) before R's C_BinCount bisection.
    private static int[] VafHistogram(List<double> vafs)
    {
        const int Bins = 100;
        var breaks = new double[Bins + 1];
        for (int k = 0; k <= Bins; k++) breaks[k] = Math.Min(k * 0.01, 1.0); // seq(0, 1, 0.01) = from + (0:n)·by, pmin(to)
        var widths = new double[Bins];
        for (int k = 0; k < Bins; k++) widths[k] = breaks[k + 1] - breaks[k];
        double diddle = 1e-7 * StatisticsHelper.Median(widths);
        var fuzzy = new double[Bins + 1];
        fuzzy[0] = breaks[0] - diddle;
        for (int k = 1; k <= Bins; k++) fuzzy[k] = breaks[k] + diddle;

        var counts = new int[Bins];
        foreach (double v in vafs)
        {
            if (!(fuzzy[0] <= v && v <= fuzzy[Bins])) continue;
            int lo = 0, hi = Bins;
            while (hi - lo >= 2)
            {
                int mid = (hi + lo) / 2;
                if (v > fuzzy[mid]) lo = mid; else hi = mid;
            }

            counts[lo]++;
        }

        return counts;
    }

    // R 1-based hst[index]; indices outside 1..100 give NA (null).
    private static int? HistogramCount(int[] histogram, double index) =>
        index >= 1 && index <= histogram.Length ? histogram[(int)index - 1] : null;

    // R ≥ 4.0 round(x, digits) (src/nmath/fround.c): pick the closer of floor/ceil(x·10^d)/10^d, even on ties.
    private static double RRound(double x, int digits)
    {
        if (double.IsNaN(x) || double.IsInfinity(x) || x == 0.0) return x;
        double sgn = 1.0;
        if (x < 0)
        {
            sgn = -1.0;
            x = -x;
        }

        double l10x = 0.30102999566398119521 * (0.5 + Math.Floor(Math.Log2(x))); // M_LOG10_2·(0.5 + logb(x))
        if (l10x + digits > 15) return sgn * x; // DBL_DIG
        double pow10 = Math.Pow(10, digits);
        double x10 = x * pow10, i10 = Math.Floor(x10);
        double xd = i10 / pow10, xu = Math.Ceiling(x10) / pow10;
        double du = xu - x, dd = x - xd;
        return sgn * ((du < dd || (i10 % 2 == 1 && du == dd)) ? xu : xd);
    }

    #endregion

}
