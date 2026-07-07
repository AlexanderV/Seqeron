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
    /// Result of the joint ASCAT purity/ploidy fit: the recovered purity ρ and ploidy ψ, the goodness of fit,
    /// and the allele-specific integer copy-number segments those parameters imply.
    /// </summary>
    /// <param name="Purity">Recovered tumour purity ρ (aberrant cell fraction) ∈ (0, 1].</param>
    /// <param name="Ploidy">Recovered tumour ploidy ψ (length-weighted mean total copy number).</param>
    /// <param name="GoodnessOfFit">Percentage goodness of fit (1 − distance/TheoretMaxdist)·100, in (−∞, 100].</param>
    /// <param name="Segments">The allele-specific integer copy-number segments (major/minor CN) implied by (ρ, ψ).</param>
    public readonly record struct PurityPloidyFit(
        double Purity,
        double Ploidy,
        double GoodnessOfFit,
        IReadOnlyList<AlleleSpecificSegment> Segments);

    /// <summary>
    /// Segments per-locus allele-specific signal (logR, BAF) into contiguous regions, producing one
    /// (mean logR, mean BAF) summary per segment. Implements a deterministic <b>joint</b> mean-shift changepoint
    /// scan on both the logR and the (mirrored) BAF tracks — the allele-specific segmentation step that precedes
    /// the ASCAT model (ASPCF; Nilsen et al. 2012, <i>BMC Genomics</i> 13:591; CBS, Olshen et al. 2004): a new
    /// segment starts when the next locus's logR deviates from the running segment mean by more than
    /// <paramref name="logRChangeThreshold"/>, OR its mirrored BAF deviates by more than
    /// <paramref name="bafChangeThreshold"/>, or when the chromosome changes. Segmenting on BAF as well as logR is
    /// essential: a copy-neutral LOH region (e.g. 2:0) has the same logR as a balanced 1:1 region but a very
    /// different BAF, so a logR-only scan would wrongly merge them. The BAF is "folded" to its distance from 0.5
    /// and re-centred (b' = 0.5 + |b − 0.5|) before averaging so that the two symmetric heterozygous BAF clusters
    /// (b and 1 − b) do not cancel — the standard mirrored-BAF summary used by allele-specific callers.
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

            double foldedBaf = BalancedBaf + Math.Abs(locus.BAF - BalancedBaf);
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
                result.Add(BuildSegmentSummary(current));
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
            result.Add(BuildSegmentSummary(current));
        }

        return result;
    }

    /// <summary>Builds a (mean logR, mirrored-mean BAF) summary from a non-empty run of same-chromosome loci.</summary>
    private static AlleleSpecificSegmentSummary BuildSegmentSummary(List<AlleleSpecificLocus> loci)
    {
        double logRSum = 0.0;
        double foldedBafSum = 0.0;
        foreach (AlleleSpecificLocus locus in loci)
        {
            logRSum += locus.LogR;
            // Mirror BAF about 0.5 so the two symmetric het clusters (b, 1−b) reinforce instead of cancel.
            foldedBafSum += BalancedBaf + Math.Abs(locus.BAF - BalancedBaf);
        }

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
    /// </summary>
    private static (double NA, double NB) AscatRawCopyNumbers(double r, double b, double rho, double psi, double gamma)
    {
        double scaledTotal = Math.Pow(2.0, r / gamma) * ((1.0 - rho) * NormalDiploidCopyNumber + rho * psi);
        double nA = (rho - 1.0 - (b - 1.0) * scaledTotal) / rho;
        double nB = (rho - 1.0 + b * scaledTotal) / rho;
        return (nA, nB);
    }

    /// <summary>
    /// Jointly estimates tumour purity ρ and ploidy ψ from segment-level (logR, BAF) summaries by grid search,
    /// mapping each segment to allele-specific copy numbers (nA, nB) with the ASCAT equations and minimising the
    /// segment-length-weighted squared distance of the minor allele to the nearest non-negative integer (the
    /// ASCAT "sunrise" goodness of fit). Source: Van Loo et al. (2010), <i>PNAS</i> 107:16910 (grid over ploidy ×
    /// aberrant-cell-fraction, "copy number calls as close as possible to nonnegative whole numbers"); equations
    /// and objective ported verbatim from ascat.runAscat.R. The returned segments carry the rounded, clamped
    /// integer major/minor copy numbers at the optimal (ρ, ψ), ready for the downstream ploidy / LOH / CCF code.
    /// </summary>
    /// <param name="segments">Segment summaries (from <see cref="SegmentAlleleSpecific"/> or a caller's segmenter). Non-empty.</param>
    /// <param name="purityMin">Lower bound of the purity grid, in (0, 1].</param>
    /// <param name="purityMax">Upper bound of the purity grid, in (0, 1] and ≥ purityMin.</param>
    /// <param name="purityStep">Purity grid step (&gt; 0).</param>
    /// <param name="ploidyMin">Lower bound of the ploidy grid (&gt; 0).</param>
    /// <param name="ploidyMax">Upper bound of the ploidy grid (≥ ploidyMin).</param>
    /// <param name="ploidyStep">Ploidy grid step (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>The recovered (ρ, ψ), the percentage goodness of fit, and the implied integer copy-number segments.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a grid bound or step is out of range.</exception>
    public static PurityPloidyFit FitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purityMin = 0.05,
        double purityMax = 1.0,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.0,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        ArgumentNullException.ThrowIfNull(segments);
        if (segments.Count == 0)
        {
            throw new ArgumentException("At least one segment is required to fit purity and ploidy.", nameof(segments));
        }

        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

        // Per-segment GoF weight = segment length (≥ 1) × balanced down-weight, exactly as ascat.runAscat.R.
        double[] weights = new double[segments.Count];
        double theoreticalMaxDistance = 0.0;
        for (int i = 0; i < segments.Count; i++)
        {
            AlleleSpecificSegmentSummary s = segments[i];
            long length = s.Length;
            // A single-locus or zero-span segment still contributes; use LocusCount as a positive weight floor.
            double baseWeight = length > 0 ? length : Math.Max(1, s.LocusCount);
            double balancedWeight = Math.Abs(s.MeanBAF - BalancedBaf) < 1e-9 ? AscatBalancedSegmentWeight : 1.0;
            weights[i] = baseWeight * balancedWeight;
            theoreticalMaxDistance += AscatWorstCaseIntegerDistance * weights[i];
        }

        double bestSelectionDistance = double.PositiveInfinity;
        double bestMinorDistance = double.PositiveInfinity;
        double bestPurity = purityMin;
        double bestPloidy = ploidyMin;

        for (double rho = purityMin; rho <= purityMax + 1e-12; rho += purityStep)
        {
            for (double psi = ploidyMin; psi <= ploidyMax + 1e-12; psi += ploidyStep)
            {
                double minorDistance = 0.0;     // ASCAT GoF objective (minor allele only).
                double selectionDistance = 0.0; // selection objective (both alleles → integers) to break 2n/4n ties.
                bool feasible = true;
                for (int i = 0; i < segments.Count; i++)
                {
                    AlleleSpecificSegmentSummary s = segments[i];
                    (double nA, double nB) = AscatRawCopyNumbers(s.MeanLogR, s.MeanBAF, rho, psi, gamma);
                    double minor = Math.Min(nA, nB);
                    double major = Math.Max(nA, nB);
                    // Physical feasibility: copy numbers cannot be meaningfully negative beyond rounding noise.
                    if (minor < -0.5 || major < -0.5)
                    {
                        feasible = false;
                        break;
                    }

                    double minorInt = Math.Max(0.0, Math.Round(minor, MidpointRounding.AwayFromZero));
                    double majorInt = Math.Max(0.0, Math.Round(major, MidpointRounding.AwayFromZero));
                    double minorDev = minor - minorInt;
                    double majorDev = major - majorInt;
                    // ascat.runAscat.R: d = sum( |nMinor - round(nMinor)|^2 * length * balancedWeight ).
                    minorDistance += minorDev * minorDev * weights[i];
                    // ASCAT rounds BOTH alleles to integers; including the major-allele deviation in the
                    // selection objective disambiguates the 2n vs 4n (doubled) solutions that share a minor fit.
                    selectionDistance += (minorDev * minorDev + majorDev * majorDev) * weights[i];
                }

                if (!feasible)
                {
                    continue;
                }

                // Prefer the lower selection distance; on a (near-)exact tie prefer the lower ploidy ψ, the ASCAT
                // parsimony convention (Van Loo 2010 selects the non-doubled solution when both fit equally well).
                bool strictlyBetter = selectionDistance < bestSelectionDistance - 1e-12;
                bool tieLowerPloidy = Math.Abs(selectionDistance - bestSelectionDistance) <= 1e-12 && psi < bestPloidy - 1e-12;
                if (strictlyBetter || tieLowerPloidy)
                {
                    bestSelectionDistance = selectionDistance;
                    bestMinorDistance = minorDistance;
                    bestPurity = rho;
                    bestPloidy = psi;
                }
            }
        }

        var bestSegments = new List<AlleleSpecificSegment>(segments.Count);
        for (int i = 0; i < segments.Count; i++)
        {
            AlleleSpecificSegmentSummary s = segments[i];
            (double nA, double nB) = AscatRawCopyNumbers(s.MeanLogR, s.MeanBAF, bestPurity, bestPloidy, gamma);
            int rMajor = (int)Math.Max(0.0, Math.Round(Math.Max(nA, nB), MidpointRounding.AwayFromZero));
            int rMinor = (int)Math.Max(0.0, Math.Round(Math.Min(nA, nB), MidpointRounding.AwayFromZero));
            // Segments with End == Start (single-position) get a 1 bp span so AlleleSpecificSegment.Length > 0.
            long end = s.End > s.Start ? s.End : s.Start + 1;
            bestSegments.Add(new AlleleSpecificSegment(s.Chromosome, s.Start, end, rMajor, rMinor));
        }

        // goodnessOfFit = (1 - distance/TheoretMaxdist) * 100, per ascat.runAscat.R (minor-allele distance).
        double goodnessOfFit = theoreticalMaxDistance > 0.0
            ? (1.0 - bestMinorDistance / theoreticalMaxDistance) * 100.0
            : 100.0;

        // Snap the reported (ρ, ψ) back to their grid bounds: floating-point accumulation in the
        // `rho += purityStep` / `psi += ploidyStep` walk can drift the top grid point a few ULPs past
        // its maximum, which would otherwise report a purity > 1 (no cell population can exceed 100 %).
        double reportedPurity = Math.Clamp(bestPurity, purityMin, purityMax);
        double reportedPloidy = Math.Clamp(bestPloidy, ploidyMin, ploidyMax);
        return new PurityPloidyFit(reportedPurity, reportedPloidy, goodnessOfFit, bestSegments);
    }

    private static void ValidateGrid(
        double purityMin, double purityMax, double purityStep,
        double ploidyMin, double ploidyMax, double ploidyStep, double gamma)
    {
        if (double.IsNaN(purityMin) || purityMin <= 0.0 || purityMin > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMin), purityMin, "purityMin must be in (0, 1].");
        }

        if (double.IsNaN(purityMax) || purityMax <= 0.0 || purityMax > 1.0 || purityMax < purityMin)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMax), purityMax, "purityMax must be in (0, 1] and ≥ purityMin.");
        }

        if (double.IsNaN(purityStep) || purityStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityStep), purityStep, "purityStep must be positive.");
        }

        if (double.IsNaN(ploidyMin) || ploidyMin <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMin), ploidyMin, "ploidyMin must be positive.");
        }

        if (double.IsNaN(ploidyMax) || ploidyMax < ploidyMin)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMax), ploidyMax, "ploidyMax must be ≥ ploidyMin.");
        }

        if (double.IsNaN(ploidyStep) || ploidyStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyStep), ploidyStep, "ploidyStep must be positive.");
        }

        if (double.IsNaN(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive.");
        }
    }

    /// <summary>
    /// Derives the integer mutation multiplicity m (number of mutated copies per cancer cell) of a somatic
    /// variant from its VAF, the tumour purity ρ, and the local total / major copy number, so that
    /// <see cref="EstimateCcf"/> can be driven without a caller-supplied multiplicity. The expected number of
    /// mutated copies for a clonal mutation is n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)] (McGranahan et al. 2016,
    /// <i>Science</i> 351:1463; equivalently the inversion of the PICTograph model VAF = m·CCF·ρ /
    /// (N_T·ρ + 2(1−ρ)) at CCF = 1, Zheng et al. 2022, <i>Bioinformatics</i> 38:3677). The result is rounded to
    /// the nearest integer and clamped to [1, majorCopyNumber] (a variant present on at least one copy cannot
    /// exceed the major-allele copy number).
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
        double totalDnaPerCell = purity * totalCopyNumber + NormalDiploidCopyNumber * (1.0 - purity);
        double rawMultiplicity = vaf * totalDnaPerCell / purity;
        int rounded = (int)Math.Round(rawMultiplicity, MidpointRounding.AwayFromZero);
        // Clamp to [1, major CN]: an observed variant sits on ≥ 1 copy and ≤ the major-allele copy number.
        return Math.Clamp(rounded, 1, majorCopyNumber);
    }

    /// <summary>
    /// Default ASPCF penalty γ used when a caller does not supply one. The copynumber package default is γ = 40
    /// (Nilsen et al. 2012, <i>BMC Genomics</i> 13:591 — "A fairly conservative penalty of γ = 40 is the default
    /// in the copynumber package"); ASCAT later raised its internal default to 70 (Ross et al. 2021,
    /// <i>Bioinformatics</i> 37:1909). Because the repository ASPCF API segments caller-supplied logR/BAF tracks on
    /// the caller's own scale, γ is exposed as a parameter; this constant only documents the published default.
    /// </summary>
    public const double AspcfDefaultPenalty = 40.0;

    /// <summary>
    /// Allele-Specific Piecewise Constant Fitting (ASPCF): the penalised-least-squares changepoint segmentation
    /// that ASCAT uses, jointly segmenting the logR and (mirrored) BAF tracks on a single common breakpoint set.
    /// Source: Nilsen et al. (2012), <i>BMC Genomics</i> 13:591 — minimise
    /// <c>L(S | y, γ) = Σ_{I∈S} Σ_{j∈I} (y_j − ȳ_I)² + γ·|S|</c> with the dynamic-program recurrence
    /// <c>e_k = min_{j≤k} ( d_{jk} + e_{j−1} + γ )</c>, <c>e_0 = 0</c>, where <c>d_{jk}</c> is the within-segment
    /// SSE of loci j..k; extended to the allele-specific joint cost
    /// <c>L(S | y₁,y₂, γ) = L(S | y₁,γ) + L(S | y₂,γ)</c> (Nilsen 2012; Ross et al. 2021, <i>Bioinformatics</i>
    /// 37:1909): a single segmentation with common breakpoints but a separate per-track segment mean, so the
    /// per-segment data cost is <c>(logR-SSE + mirroredBAF-SSE)</c> and γ is charged once per segment. This returns
    /// the GLOBAL optimum of the penalised cost (unlike the greedy <see cref="SegmentAlleleSpecific"/> mean-shift).
    /// BAF is mirrored to its distance from 0.5 and re-centred (<c>b' = 0.5 + |b − 0.5|</c>) so the two symmetric
    /// het clusters collapse to one track (Ross 2021 — "mirroring BAFs to obtain a single track in regions of
    /// allelic imbalance"); without this a copy-neutral LOH (2:0) and a balanced 1:1 region — equal logR — would be
    /// merged. The DP runs per chromosome (breakpoints never cross a contig boundary). Time O(n²) per chromosome.
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome.</param>
    /// <param name="penalty">Penalty γ &gt; 0 charged per segment (see <see cref="AspcfDefaultPenalty"/>).</param>
    /// <returns>The segment summaries (mean logR, mirrored-mean BAF) in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException">a locus has a null chromosome label.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0 or NaN.</exception>
    public static IReadOnlyList<AlleleSpecificSegmentSummary> SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(loci);
        if (double.IsNaN(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive.");
        }

        // Materialise and group into contiguous same-chromosome runs (input order preserved within each).
        var ordered = new List<AlleleSpecificLocus>();
        foreach (AlleleSpecificLocus locus in loci)
        {
            if (locus.Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            ordered.Add(locus);
        }

        var result = new List<AlleleSpecificSegmentSummary>();
        int start = 0;
        while (start < ordered.Count)
        {
            int end = start;
            while (end + 1 < ordered.Count && ordered[end + 1].Chromosome == ordered[start].Chromosome)
            {
                end++;
            }

            SegmentChromosomeAspcf(ordered, start, end, penalty, result);
            start = end + 1;
        }

        return result;
    }

    /// <summary>
    /// Runs the PCF dynamic program on one chromosome's run of loci (inclusive indices [lo, hi]) and appends the
    /// optimal segments to <paramref name="output"/>. Implements Nilsen et al. (2012) eq. for the joint cost.
    /// </summary>
    private static void SegmentChromosomeAspcf(
        List<AlleleSpecificLocus> loci, int lo, int hi, double penalty,
        List<AlleleSpecificSegmentSummary> output)
    {
        int n = hi - lo + 1;

        // Prefix sums for O(1) within-segment SSE: SSE(a..b) = Σx² − (Σx)²/m, for both tracks (Nilsen 2012 L′).
        double[] logRPrefix = new double[n + 1];
        double[] logRSqPrefix = new double[n + 1];
        double[] bafPrefix = new double[n + 1];
        double[] bafSqPrefix = new double[n + 1];
        for (int i = 0; i < n; i++)
        {
            double r = loci[lo + i].LogR;
            // Mirror BAF about 0.5 → single allelic-imbalance track (Ross 2021).
            double b = BalancedBaf + Math.Abs(loci[lo + i].BAF - BalancedBaf);
            logRPrefix[i + 1] = logRPrefix[i] + r;
            logRSqPrefix[i + 1] = logRSqPrefix[i] + r * r;
            bafPrefix[i + 1] = bafPrefix[i] + b;
            bafSqPrefix[i + 1] = bafSqPrefix[i] + b * b;
        }

        // e[k] = min penalised cost of segmenting the first k loci; back[k] = start index of the last segment.
        double[] e = new double[n + 1];
        int[] back = new int[n + 1];
        e[0] = 0.0;
        for (int k = 1; k <= n; k++)
        {
            double best = double.PositiveInfinity;
            int bestStart = 0;
            for (int j = 1; j <= k; j++)
            {
                // d_{jk} = within-segment SSE of loci (j..k) on both tracks (joint cost = sum of the two SSEs).
                double cost = e[j - 1] + AspcfSegmentSse(logRPrefix, logRSqPrefix, j - 1, k)
                              + AspcfSegmentSse(bafPrefix, bafSqPrefix, j - 1, k)
                              + penalty;
                if (cost < best - 1e-12)
                {
                    best = cost;
                    bestStart = j - 1;
                }
            }

            e[k] = best;
            back[k] = bestStart;
        }

        // Backtrack the optimal segmentation, then emit in genomic (left-to-right) order.
        var bounds = new List<(int Start, int End)>();
        int cursor = n;
        while (cursor > 0)
        {
            int segStart = back[cursor];
            bounds.Add((segStart, cursor)); // half-open [segStart, cursor) over the local 0-based run.
            cursor = segStart;
        }

        bounds.Reverse();
        foreach ((int segStart, int segEnd) in bounds)
        {
            var run = new List<AlleleSpecificLocus>(segEnd - segStart);
            for (int i = segStart; i < segEnd; i++)
            {
                run.Add(loci[lo + i]);
            }

            output.Add(BuildSegmentSummary(run));
        }
    }

    /// <summary>Within-segment SSE for the half-open prefix range (a, b]: Σx² − (Σx)²/m (m = b − a). 0 if m ≤ 1.</summary>
    private static double AspcfSegmentSse(double[] prefix, double[] sqPrefix, int a, int b)
    {
        int m = b - a;
        if (m <= 1)
        {
            return 0.0; // a single point has zero within-segment variance.
        }

        double sum = prefix[b] - prefix[a];
        double sumSq = sqPrefix[b] - sqPrefix[a];
        double sse = sumSq - (sum * sum) / m;
        return sse > 0.0 ? sse : 0.0; // guard tiny negative round-off.
    }

    /// <summary>
    /// One integer allele-specific copy-number state of a (possibly sub-clonal) segment, present in a given
    /// fraction of tumour cells. Mirrors the Battenberg output (Nik-Zainal et al. 2012, <i>Cell</i> 149:994:
    /// <c>nMaj1_A, nMin1_A, frac1_A</c>): a state is a (major, minor) integer pair plus the cellular fraction.
    /// </summary>
    /// <param name="MajorCopyNumber">Major-allele integer copy number of this state (≥ 0).</param>
    /// <param name="MinorCopyNumber">Minor-allele integer copy number of this state (≥ 0).</param>
    /// <param name="CellFraction">Fraction of tumour cells carrying this state, ∈ [0, 1].</param>
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
    /// <param name="PrimaryState">State 1 (Battenberg <c>frac1</c>) — the higher-fraction state.</param>
    /// <param name="SecondaryState">State 2 (Battenberg <c>frac2</c>), or <c>null</c> for a clonal segment.</param>
    /// <param name="IsSubclonal">True when two states were needed (the observed CN was not (near-)integer).</param>
    public readonly record struct SubclonalSegmentFit(
        AlleleSpecificSegmentSummary Segment,
        SubclonalCopyNumberState PrimaryState,
        SubclonalCopyNumberState? SecondaryState,
        bool IsSubclonal);

    /// <summary>
    /// Maximum distance of an allele-specific copy number from the nearest integer below which the segment is
    /// called <b>clonal</b> (a single integer state). Beyond it the segment is modelled as a sub-clonal mixture of
    /// the two bracketing integers. 0.05 mirrors ASCAT's "as close as possible to nonnegative whole numbers"
    /// integer-snapping tolerance (Van Loo et al. 2010) used by Battenberg to decide clonal vs sub-clonal.
    /// </summary>
    public const double SubclonalIntegerTolerance = 0.05;

    /// <summary>
    /// Fits each segment's allele-specific copy number to one integer state (clonal) or a mixture of two adjacent
    /// integer states with a sub-clonal cellular fraction (sub-clonal), implementing the Battenberg two-population
    /// model (Nik-Zainal et al. 2012, <i>Cell</i> 149:994; Wedge-lab/battenberg): "if there are two states it
    /// represents subclonal copy number … two populations of cells, each with a different state … which together
    /// give the total copy number for that segment and a fraction of tumour cells that carry each allele." The
    /// real-valued ASCAT allele-specific copy numbers (nA, nB) are computed for the segment at the fitted (ρ, ψ)
    /// via the ASCAT equations (Van Loo et al. 2010); a value that is within <see cref="SubclonalIntegerTolerance"/>
    /// of an integer collapses to a single (clonal) state, otherwise it is decomposed as
    /// <c>n_obs = f·⌈n_obs⌉ + (1 − f)·⌊n_obs⌋</c> with <c>f = n_obs − ⌊n_obs⌋ ∈ [0,1]</c> (the unique two-state
    /// mixture reproducing n_obs). The two alleles are decomposed jointly with a single shared fraction f estimated
    /// as the mean of the per-allele fractions, so the two states are (⌈nA⌉, ⌈nB⌉) at fraction f and
    /// (⌊nA⌋, ⌊nB⌋) at fraction 1 − f, matching the Battenberg (nMaj1/nMin1/frac1, nMaj2/nMin2/frac2) layout.
    /// </summary>
    /// <param name="segments">Segment summaries (e.g. from <see cref="SegmentAlleleSpecificAspcf"/>). Non-null.</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>Per-segment clonal/sub-clonal copy-number fits in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0,1], ψ ≤ 0, or γ ≤ 0.</exception>
    public static IReadOnlyList<SubclonalSegmentFit> FitSubclonalCopyNumber(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        ArgumentNullException.ThrowIfNull(segments);
        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity ρ must be in (0, 1].");
        }

        if (double.IsNaN(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy ψ must be positive.");
        }

        if (double.IsNaN(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive.");
        }

        var fits = new List<SubclonalSegmentFit>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            (double nA, double nB) = AscatRawCopyNumbers(s.MeanLogR, s.MeanBAF, purity, ploidy, gamma);
            double major = Math.Max(0.0, Math.Max(nA, nB));
            double minor = Math.Max(0.0, Math.Min(nA, nB));

            double majorFrac = major - Math.Floor(major); // distance above the lower bracketing integer.
            double minorFrac = minor - Math.Floor(minor);

            // Clonal when BOTH alleles snap to integers within tolerance; else a two-state mixture is required.
            bool majorClonal = majorFrac <= SubclonalIntegerTolerance || majorFrac >= 1.0 - SubclonalIntegerTolerance;
            bool minorClonal = minorFrac <= SubclonalIntegerTolerance || minorFrac >= 1.0 - SubclonalIntegerTolerance;

            if (majorClonal && minorClonal)
            {
                int majInt = (int)Math.Max(0.0, Math.Round(major, MidpointRounding.AwayFromZero));
                int minInt = (int)Math.Max(0.0, Math.Round(minor, MidpointRounding.AwayFromZero));
                fits.Add(new SubclonalSegmentFit(
                    s,
                    new SubclonalCopyNumberState(majInt, minInt, 1.0),
                    SecondaryState: null,
                    IsSubclonal: false));
                continue;
            }

            // Two-state mixture (Battenberg single shared fraction): each allele is a convex combination of its two
            // bracketing integers, both alleles sharing ONE fraction f. The pairing of the alleles' ceil/floor
            // integers into the two cell populations is ambiguous, so both pairings are tried and the one with the
            // smaller least-squares residual is kept — the unique two-state decomposition reproducing (major, minor).
            int majCeil = (int)Math.Ceiling(major);
            int majFloor = (int)Math.Floor(major);
            int minCeil = (int)Math.Ceiling(minor);
            int minFloor = (int)Math.Floor(minor);

            // Pairing P1 (co-monotone): state_hi = (majCeil, minCeil), state_lo = (majFloor, minFloor).
            (double fP1, double resP1) = SolveSharedFraction(major, minor, majCeil, minCeil, majFloor, minFloor);
            // Pairing P2 (anti-monotone): state_hi = (majCeil, minFloor), state_lo = (majFloor, minCeil).
            (double fP2, double resP2) = SolveSharedFraction(major, minor, majCeil, minFloor, majFloor, minCeil);

            SubclonalCopyNumberState hiState, loState; // hi = the "ceiling-on-major" state (fraction f).
            double f;
            if (resP1 <= resP2)
            {
                f = fP1;
                hiState = new SubclonalCopyNumberState(majCeil, minCeil, f);
                loState = new SubclonalCopyNumberState(majFloor, minFloor, 1.0 - f);
            }
            else
            {
                f = fP2;
                hiState = new SubclonalCopyNumberState(majCeil, minFloor, f);
                loState = new SubclonalCopyNumberState(majFloor, minCeil, 1.0 - f);
            }

            // Battenberg frac1 ≥ frac2: the higher-fraction state is the primary (state 1).
            (SubclonalCopyNumberState primary, SubclonalCopyNumberState secondary) =
                f >= 0.5 ? (hiState, loState) : (loState, hiState);

            fits.Add(new SubclonalSegmentFit(s, primary, secondary, IsSubclonal: true));
        }

        return fits;
    }

    /// <summary>
    /// Solves for the single shared cellular fraction f that best reproduces both observed alleles as
    /// <c>major = f·aHi + (1−f)·aLo</c> and <c>minor = f·bHi + (1−f)·bLo</c> by least squares, returning f
    /// (clamped to [0,1]) and the residual sum of squares of the fit. Two integer states share one f per the
    /// Battenberg single-fraction segment model (Nik-Zainal et al. 2012).
    /// </summary>
    private static (double F, double Residual) SolveSharedFraction(
        double major, double minor, int aHi, int bHi, int aLo, int bLo)
    {
        // For each allele: observed = aLo + f·(aHi − aLo)  ⇒ stack the two equations and solve LS for f.
        double da = aHi - aLo;
        double db = bHi - bLo;
        double denom = da * da + db * db;
        double f;
        if (denom < 1e-12)
        {
            f = 0.0; // both states identical for both alleles → degenerate; fraction is irrelevant.
        }
        else
        {
            f = (da * (major - aLo) + db * (minor - bLo)) / denom;
            f = Math.Clamp(f, 0.0, 1.0);
        }

        double majFit = aLo + f * da;
        double minFit = bLo + f * db;
        double residual = (major - majFit) * (major - majFit) + (minor - minFit) * (minor - minFit);
        return (f, residual);
    }

    #endregion

}
