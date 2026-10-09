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
        return ploidy * Math.Pow(2.0, log2Ratio);
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
    /// Segment coordinates are half-open (<c>Length = End − Start</c>).
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
    /// Shared arm → gene-panel mapping of <see cref="IdentifyAmplifiedOncogenes"/> (ONCO-CNA-002) and
    /// <see cref="IdentifyDeletedTumorSuppressors"/> (ONCO-CNA-003): collects the distinct (case-insensitive,
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
    /// table is the denominator of the reference-genome WGD fraction
    /// (<see cref="DetectWholeGenomeDoubling(IEnumerable{AlleleSpecificSegment}, ReferenceGenome)"/>). Note: in
    /// facets-suite the <c>genome</c> build (<c>'hg19' | 'hg18' | 'hg38'</c>) supplies centromere positions only;
    /// its WGD denominator comes from <c>get_sample_genome</c> (interrogated span), see
    /// <see cref="DetectWholeGenomeDoublingFromSuppliedLength"/>.
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
    /// Determines whether a tumour genome has undergone whole-genome doubling (WGD) by the Bielski et al. rule —
    /// more than half of the <i>autosomal genome</i> (chromosomes 1–22) has major-allele copy number ≥ 2
    /// (Bielski et al. 2018, Nat Genet 50:1189–1195, PMID 30013179) — taking "autosomal genome" to be the
    /// <b>reference assembly's</b> autosomal length (Σ chr1–22 of the UCSC <c>*.chrom.sizes</c> table).
    /// Numerator, threshold and major-CN rule follow facets-suite <c>is_genome_doubled</c>:
    /// <c>frac_elevated_mcn = sum(length where mcn ≥ 2 &amp; chrom %in% 1:22) / autosomal_genome</c>;
    /// <c>wgd = frac_elevated_mcn &gt; 0.5</c>, <c>mcn = tcn − lcn</c> = max(Major, Minor).
    /// <para><b>Denominator deviation from the facets-suite code:</b> facets-suite passes
    /// <c>get_sample_genome(segs, genome)</c> as <c>chrom_info</c>, whose per-chromosome <c>size</c> is the
    /// interrogated span <c>max(end) − min(start)</c> of the sample's own segments — not the reference chromosome
    /// length. That exact rule is <see cref="DetectWholeGenomeDoublingFromSuppliedLength"/>. This overload instead
    /// divides by the fixed reference autosomal length, so partial-genome inputs (a few chromosomes, or segments
    /// that do not reach the telomeres) are judged against the whole autosomal genome; for a segmentation spanning
    /// every autosome end-to-end the two agree.</para>
    /// The test uses the major (not total) copy number, so a balanced diploid genome (all 1:1, total CN 2,
    /// major CN 1) is NOT doubled, whereas a 2:0 LOH or 2:2 genome IS.
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
    /// Determines whole-genome doubling exactly as the facets-suite reference implementation does, with the
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
    /// segments' total copy number. It equals ψ for an exact integer-copy-number genome.</param>
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

    /// <summary>
    /// Segments per-locus allele-specific signal (logR, BAF) into contiguous regions, producing one
    /// (mean logR, mean BAF) summary per segment, with a deterministic <b>greedy joint mean-shift heuristic</b> on the
    /// logR and the (mirrored) BAF tracks. <b>This is not ASPCF or CBS</b> and has no published reference
    /// implementation: it is a single left-to-right pass with caller-chosen absolute thresholds, no noise
    /// standardisation and no global (penalised least-squares) optimisation, so its breakpoints differ from ASCAT's.
    /// For the published allele-specific segmentation (ASCAT <c>ascat.aspcf</c>; Nilsen et al. 2012, <i>BMC Genomics</i>
    /// 13:591; Ross et al. 2021) use <see cref="SegmentAlleleSpecificAspcf"/>. Rule: a new
    /// segment starts when the next locus's logR deviates from the running segment mean by more than
    /// <paramref name="logRChangeThreshold"/>, OR its mirrored BAF deviates by more than
    /// <paramref name="bafChangeThreshold"/>, or when the chromosome changes. Segmenting on BAF as well as logR is
    /// essential: a copy-neutral LOH region (e.g. 2:0) has the same logR as a balanced 1:1 region but a very
    /// different BAF, so a logR-only scan would wrongly merge them. The BAF is "folded" to its distance from 0.5
    /// and re-centred (b' = 0.5 + |b − 0.5|) before averaging so that the two symmetric heterozygous BAF clusters
    /// (b and 1 − b) do not cancel — the standard mirrored-BAF summary used by allele-specific callers.
    /// Unlike <see cref="SegmentAlleleSpecificAspcf"/>, locus values are not validated (a BAF outside [0, 1] yields a
    /// folded mean above 1, which <see cref="FitPurityPloidy"/> then rejects).
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome.</param>
    /// <param name="logRChangeThreshold">logR mean-shift threshold that starts a new segment. Must be &gt; 0.</param>
    /// <param name="bafChangeThreshold">Mirrored-BAF mean-shift threshold that starts a new segment. Must be &gt; 0.</param>
    /// <param name="minLociPerSegment">Minimum loci a running segment must have before a change can split it. Must be ≥ 1.</param>
    /// <returns>The segment summaries in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a threshold ≤ 0 or minLociPerSegment &lt; 1.</exception>
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

        var result = new List<AlleleSpecificSegmentSummary>();
        var current = new List<AlleleSpecificLocus>();
        double runningLogRSum = 0.0;
        double runningFoldedBafSum = 0.0;

        foreach (AlleleSpecificLocus locus in loci)
        {
            if (locus.Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            double foldedBaf = FoldBafAboutHalf(locus.BAF);
            bool chromosomeChanged = current.Count > 0 && current[^1].Chromosome != locus.Chromosome;
            bool meanShift = false;
            if (!chromosomeChanged && current.Count >= minLociPerSegment)
            {
                double currentLogRMean = runningLogRSum / current.Count;
                double currentBafMean = runningFoldedBafSum / current.Count;
                // ASPCF/CBS joint mean-shift: split on a logR change OR a (mirrored) BAF change.
                meanShift = Math.Abs(locus.LogR - currentLogRMean) > logRChangeThreshold
                            || Math.Abs(foldedBaf - currentBafMean) > bafChangeThreshold;
            }

            if ((chromosomeChanged || meanShift) && current.Count > 0)
            {
                result.Add(BuildSegmentSummary(current, runningLogRSum, runningFoldedBafSum));
                current = new List<AlleleSpecificLocus>();
                runningLogRSum = 0.0;
                runningFoldedBafSum = 0.0;
            }

            current.Add(locus);
            runningLogRSum += locus.LogR;
            runningFoldedBafSum += foldedBaf;
        }

        if (current.Count > 0)
        {
            result.Add(BuildSegmentSummary(current, runningLogRSum, runningFoldedBafSum));
        }

        return result;
    }

    /// <summary>Greedy-segmenter BAF fold: b' = 0.5 + |b − 0.5|, so the two symmetric het clusters (b, 1−b) reinforce
    /// instead of cancel when averaged.</summary>
    private static double FoldBafAboutHalf(double baf) => BalancedBaf + Math.Abs(baf - BalancedBaf);

    /// <summary>
    /// Builds a (mean logR, mirrored-mean BAF) summary from a non-empty run of same-chromosome loci, given the run's
    /// logR sum and folded-BAF sum (accumulated in locus order by the caller).
    /// </summary>
    private static AlleleSpecificSegmentSummary BuildSegmentSummary(
        List<AlleleSpecificLocus> loci, double logRSum, double foldedBafSum)
    {
        return new AlleleSpecificSegmentSummary(
            Chromosome: loci[0].Chromosome,
            Start: loci[0].Position,
            End: loci[^1].Position,
            MeanLogR: logRSum / loci.Count,
            MeanBAF: foldedBafSum / loci.Count,
            LocusCount: loci.Count);
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
    private static (double Major, double Minor) AscatRoundSegment(double r, double bAscat, double rho, double psi, double gamma)
    {
        (double nAraw, double nBraw) = AscatRawCopyNumbers(r, bAscat, rho, psi, gamma);
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
        bool nonAberrant)
    {
        var result = new List<AlleleSpecificSegment>(segments.Count);
        double cnSum = 0.0, probeSum = 0.0;
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            (double major, double minor) = AscatRoundSegment(s.MeanLogR, ToAscatBaf(s.MeanBAF), rho, psi, gamma);
            int majorInt = AscatCopyNumberToInt(major);
            int minorInt = AscatCopyNumberToInt(minor);
            // Segments with End == Start (single-position) get a 1 bp span so AlleleSpecificSegment.Length > 0.
            long end = s.End > s.Start ? s.End : s.Start + 1;
            result.Add(new AlleleSpecificSegment(s.Chromosome, s.Start, end, majorInt, minorInt));
            cnSum += ((double)majorInt + minorInt) * s.LocusCount;
            probeSum += s.LocusCount;
        }

        return new PurityPloidyFit(rho, cnSum / probeSum, goodnessOfFit, result)
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
    /// Sex-chromosome (X/Y) segments are excluded from the fit and emitted with the diploid model (ASCAT gender "XX").
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
    {
        if (!TryFitPurityPloidy(segments, out PurityPloidyFit fit, purityMin, purityMax, purityStep,
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
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

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

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);

        var candidates = new List<(double M, int I, int J, double GoodnessOfFit)>();
        bool strictPloidyWindowReachable = ploidyMin < AscatMaxPloidyStrict && ploidyMax > AscatMinPloidyStrict;

        // Pass 1: all filters.
        CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
            !nonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
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
                !nonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
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
            fit = default;
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

        double psiOpt = RDimnameValue(psiPos[best.I]);
        double rhoOpt = Math.Min(1.0, RDimnameValue(rhoPos[best.J])); // if (rho_opt1 > 1) rho_opt1 = 1
        fit = BuildAscatFit(segments, rhoOpt, psiOpt, gamma, best.GoodnessOfFit, nonAberrant);
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
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);
        double m = AscatDistance(s, purity, ploidy, gamma);
        double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
        return BuildAscatFit(segments, purity, ploidy, gamma, goodnessOfFit, nonAberrant);
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
    /// Not ported (no inputs for them): germline-homozygous-stretch resegmentation and the averaging of homozygous
    /// probes' logR (every supplied locus is treated as a germline-heterozygous SNP).
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
    /// maxdist rule can make a segment clonal — exactly the path ported here. The bootstrap confidence intervals and the
    /// alternative solutions B–F of Battenberg's output are not produced.
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

        double rho = purity;
        double psiAll = MixtureCopiesPerCell(rho, ploidy); // psi = rho*psit + 2*(1-rho)
        var fits = new List<SubclonalSegmentFit>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            if (!IsValidAlleleSignal(s.MeanLogR, s.MeanBAF))
            {
                throw new ArgumentException("Every segment needs a finite mean logR and a mean BAF in [0, 1].", nameof(segments));
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
            double closestDistance = firstClosest ? dist1 : dist2;

            if (closestDistance < BattenbergMaxBafDistance)
            {
                var clonal = firstClosest
                    ? new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), 1.0)
                    : new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0);
                fits.Add(new SubclonalSegmentFit(s, clonal, SecondaryState: null, IsSubclonal: false));
                continue;
            }

            double tau = (1.0 - rho + rho * maj2 - 2.0 * l * (1.0 - rho) - l * rho * (min2 + maj2))
                         / (l * rho * (min1 + maj1) - l * rho * (min2 + maj2) - rho * maj1 + rho * maj2);
            fits.Add(new SubclonalSegmentFit(
                s,
                new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), tau),
                new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0 - tau),
                IsSubclonal: true));
        }

        return fits;
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

    #endregion

}
