using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs hybridization probes for various applications (FISH, microarray, Northern blot, etc.).
/// </summary>
public static class ProbeDesigner
{
    #region Records

    /// <summary>
    /// Probe design parameters.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Tm.</b> <see cref="MinTm"/>/<see cref="MaxTm"/> are compared with the probe Tm computed by
    /// Primer3's <c>seqtm</c> (<see cref="PrimerDesigner.CalculateMeltingTemperaturePrimer3"/>: SantaLucia 1998
    /// nearest-neighbour Tm with the SantaLucia salt correction for ≤ 36 nt, Primer3 <c>long_seq_tm</c>
    /// 81.5 + 16.6·log10([Mon]_eq) + 0.41·%GC − 600/N for longer probes) at the hybridization conditions
    /// <see cref="DnaConcentrationNanomolar"/>, <see cref="MonovalentMillimolar"/>, <see cref="DivalentMillimolar"/>,
    /// <see cref="DntpMillimolar"/>. The defaults are Primer3's hybridization-probe (internal-oligo) conditions
    /// PRIMER_INTERNAL_DNA_CONC = 50 nM, PRIMER_INTERNAL_SALT_MONOVALENT = 50 mM, PRIMER_INTERNAL_SALT_DIVALENT = 0,
    /// PRIMER_INTERNAL_DNTP_CONC = 0 (<c>libprimer3.cc</c> <c>pr_set_default_global_args</c>). The nearest-neighbour
    /// length limit is <see cref="MaxNearestNeighborLength"/> (Primer3 <c>nn_max_len</c>, default 36). The presets in
    /// <see cref="Defaults"/> state their own conditions; each preset's Tm window is reachable on this scale for its
    /// length and G+C window (audit round 3, A3-10).
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> With <see cref="StructureScreen"/> = <see cref="ProbeStructureScreen.Thermodynamic"/>
    /// (default) a probe of ≤ <see cref="ThermodynamicScreenMaxLength"/> (default 60) nt made of A/C/G/T is screened as Primer3 screens a hybridization probe
    /// (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1): ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin
    /// Tm (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>) must not exceed <see cref="MaxStructureTm"/>
    /// (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH = 47 °C). The self-dimer limit replaces
    /// <see cref="MaxSelfComplementarity"/>; the hairpin limit is applied when <see cref="AvoidSecondaryStructure"/>.
    /// Probes longer than <see cref="ThermodynamicScreenMaxLength"/> (thal.c <c>THAL_MAX_ALIGN</c>, default 60; a
    /// larger value is an opt-in that runs the unchanged ntthal recursions on longer probes), probes with non-ACGT
    /// bases and <see cref="ProbeStructureScreen.Heuristic"/> use the fallback screens: the self-dimer criterion is Primer3's
    /// alignment-mode internal-oligo screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0, <c>oligo_compl</c>): the dpal
    /// <c>self_any</c> (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>) and <c>self_end</c>
    /// (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>) scores, which have no length limit, must not
    /// exceed <see cref="MaxSelfAny"/> / <see cref="MaxSelfEnd"/> (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00);
    /// the hairpin criterion stays the sequence-only inverted-repeat stem screen (≥ 4 bp, loop 3, ≥ 80 % matched) —
    /// with the default THAL_MAX_ALIGN = 60 a &gt; 60-nt probe gets no thermodynamic hairpin; raise
    /// <see cref="ThermodynamicScreenMaxLength"/> for the ntthal hairpin of longer probes (audit round 3, A3-9).
    /// <see cref="MaxSelfComplementarity"/> (position-wise fold-back fraction limit) is
    /// no longer used by any screen; it is kept for source compatibility.
    /// </para>
    /// </remarks>
    public readonly record struct ProbeParameters(
        int MinLength,
        int MaxLength,
        double MinTm,
        double MaxTm,
        double MinGc,
        double MaxGc,
        int MaxHomopolymer,
        bool AvoidSecondaryStructure,
        double MaxSelfComplementarity)
    {
        /// <summary>Probe (oligo) concentration in nM for the Tm (Primer3 PRIMER_INTERNAL_DNA_CONC, default 50).</summary>
        public double DnaConcentrationNanomolar { get; init; } = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar;

        /// <summary>Monovalent cation concentration in mM (Primer3 PRIMER_INTERNAL_SALT_MONOVALENT, default 50).</summary>
        public double MonovalentMillimolar { get; init; } = PrimerDesigner.Primer3InternalMonovalentMillimolar;

        /// <summary>Mg²⁺ concentration in mM (Primer3 PRIMER_INTERNAL_SALT_DIVALENT, default 0).</summary>
        public double DivalentMillimolar { get; init; } = PrimerDesigner.Primer3InternalDivalentMillimolar;

        /// <summary>dNTP concentration in mM (Primer3 PRIMER_INTERNAL_DNTP_CONC, default 0).</summary>
        public double DntpMillimolar { get; init; } = PrimerDesigner.Primer3InternalDntpMillimolar;

        /// <summary>
        /// Longest probe (nt) whose Tm is the SantaLucia 1998 nearest-neighbour Tm; longer probes get Primer3's
        /// <c>long_seq_tm</c> (Primer3 <c>seqtm</c> <c>nn_max_len</c>; default <see cref="PrimerDesigner.Primer3MaxNnTmLength"/>
        /// = 36, Primer3's MAX_NN_TM_LENGTH; ≥ 0).
        /// </summary>
        public int MaxNearestNeighborLength { get; init; } = PrimerDesigner.Primer3MaxNnTmLength;

        /// <summary>Self-structure screen (default <see cref="ProbeStructureScreen.Thermodynamic"/>).</summary>
        public ProbeStructureScreen StructureScreen { get; init; } = ProbeStructureScreen.Thermodynamic;

        /// <summary>
        /// Maximum ntthal self-dimer / 3′ self-dimer / hairpin Tm (°C) of the thermodynamic screen
        /// (Primer3 PRIMER_INTERNAL_MAX_SELF_ANY_TH = _SELF_END_TH = _HAIRPIN_TH = 47 °C).
        /// </summary>
        public double MaxStructureTm { get; init; } = PrimerDesigner.Primer3MaxStructureTm;

        /// <summary>
        /// Longest probe (nt) screened by the <see cref="ProbeStructureScreen.Thermodynamic"/> ntthal screen — thal.h
        /// <c>THAL_MAX_ALIGN</c>. Default <see cref="PrimerDesigner.NtthalMaxAlignLength"/> = 60 (Primer3 / primer3-py;
        /// unchanged behaviour). THAL_MAX_ALIGN is only a compile-time guard of thal.c, so a larger value (opt-in,
        /// ≤ <see cref="PrimerDesigner.NtthalMaxSequenceLength"/>) screens longer ACGT probes with the unchanged ntthal
        /// self-dimer / 3′ self-dimer / hairpin recursions (= ntthal compiled with <c>-DTHAL_MAX_ALIGN=…</c>) instead of
        /// the fallback screens. Cost grows as O(n²·30²) per probe and screen — keep it to the probe lengths in use.
        /// Also the THAL_MAX_ALIGN of the <see cref="AssessCrossHybridization"/> / <see cref="ValidateProbe"/> non-target site
        /// duplex Tm (audit round 3, A3-27). Outside 60–10 000 → <see cref="ArgumentOutOfRangeException"/>.
        /// </summary>
        public int ThermodynamicScreenMaxLength { get; init; } = PrimerDesigner.NtthalMaxAlignLength;

        /// <summary>
        /// Maximum Primer3 alignment-mode <c>self_any</c> of the fallback self-dimer screen
        /// (PRIMER_INTERNAL_MAX_SELF_ANY, Primer3 default 12.00; flagged when strictly greater).
        /// </summary>
        public double MaxSelfAny { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>
        /// Maximum Primer3 alignment-mode <c>self_end</c> of the fallback self-dimer screen
        /// (PRIMER_INTERNAL_MAX_SELF_END, Primer3 default 12.00; flagged when strictly greater).
        /// </summary>
        public double MaxSelfEnd { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>
        /// Ranking of <see cref="DesignProbes(string, ProbeParameters?, int)"/> (default
        /// <see cref="ProbeRanking.AdditiveScore"/>, the library heuristic). <see cref="ProbeRanking.Primer3Penalty"/>
        /// ranks by Primer3's internal-oligo objective <c>p_obj_fn</c> (PRIMER_INTERNAL_n_PENALTY) around
        /// <see cref="OptTm"/> / <see cref="OptLength"/> — the sourced ranking.
        /// </summary>
        public ProbeRanking Ranking { get; init; } = ProbeRanking.AdditiveScore;

        /// <summary>
        /// Optimum probe Tm (°C) of the <see cref="ProbeRanking.Primer3Penalty"/> ranking (Primer3 PRIMER_INTERNAL_OPT_TM,
        /// default 60 °C, <c>pr_set_default_global_args</c> <c>o_args.opt_tm</c>). With that ranking it must lie in
        /// [<see cref="MinTm"/>, <see cref="MaxTm"/>] (Primer3 <c>_pr_data_control</c>: "Optimum internal oligo Tm lower
        /// than minimum or higher than maximum"); unused by the additive score.
        /// </summary>
        public double OptTm { get; init; } = PrimerDesigner.DefaultPrimer3Optima.OptTm;

        /// <summary>
        /// Optimum probe length (nt) of the <see cref="ProbeRanking.Primer3Penalty"/> ranking (Primer3 PRIMER_INTERNAL_OPT_SIZE,
        /// default 20, <c>o_args.opt_size</c>). With that ranking it must lie in [<see cref="MinLength"/>, <see cref="MaxLength"/>]
        /// (Primer3 <c>_pr_data_control</c>: "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE &gt; MAX_SIZE" / "&lt; MIN_SIZE");
        /// unused by the additive score.
        /// </summary>
        public int OptLength { get; init; } = PrimerDesigner.DefaultPrimer3Optima.OptSize;
    }

    /// <summary>Ranking order of <see cref="DesignProbes(string, ProbeParameters?, int)"/> (see <see cref="ProbeParameters.Ranking"/>).</summary>
    public enum ProbeRanking
    {
        /// <summary>
        /// Library heuristic (default, unchanged): 1.0 minus fixed penalties (GC 0.3, Tm 0.3, homopolymer 0.2,
        /// self-complementarity 0.2, secondary structure 0.15, repeats 0.1, terminal G/C 0.02 each), score descending.
        /// The penalty values have no published source.
        /// </summary>
        AdditiveScore,

        /// <summary>
        /// Primer3's internal-oligo objective <c>p_obj_fn</c> (<c>OT_INTL</c> branch, Primer3 default internal-oligo
        /// weights: PRIMER_INTERNAL_WT_TM_GT/_LT = PRIMER_INTERNAL_WT_SIZE_GT/_LT = 1, every other weight 0), i.e.
        /// |Tm − <see cref="ProbeParameters.OptTm"/>| + |length − <see cref="ProbeParameters.OptLength"/>| with the probe Tm
        /// at the parameters' conditions (<see cref="PrimerDesigner.CalculatePrimer3Penalty"/>), ordered as Primer3's
        /// <c>primer_rec_comp</c>: penalty ascending, then start descending, then length ascending. Sourced ranking
        /// (Primer3 <c>libprimer3.cc</c>; = primer3-py PRIMER_INTERNAL_n_PENALTY).
        /// </summary>
        Primer3Penalty
    }

    /// <summary>Self-structure screen used by the probe designers (see <see cref="ProbeParameters"/>).</summary>
    public enum ProbeStructureScreen
    {
        /// <summary>Primer3 thermodynamic screen (ntthal self-dimer, 3′ self-dimer, hairpin Tm) for ACGT probes of ≤ <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (default 60) nt.</summary>
        Thermodynamic,

        /// <summary>The fallback screens for every probe: Primer3 alignment-mode self_any / self_end (dpal, limits
        /// <see cref="ProbeParameters.MaxSelfAny"/> / <see cref="ProbeParameters.MaxSelfEnd"/>) and the sequence-only
        /// inverted-repeat hairpin stem.</summary>
        Heuristic
    }

    /// <summary>
    /// Default probe parameters for different applications.
    /// </summary>
    /// <remarks>
    /// Every preset's Tm window is reachable on its own Tm scale and conditions for its length and G+C window
    /// (audit round 3, A3-10; primer3-py 2.3.1 <c>calc_tm</c> witnesses in the PROBE-DESIGN-001 tests). Primer3
    /// <c>long_seq_tm</c> (the Tm of every probe &gt; <see cref="ProbeParameters.MaxNearestNeighborLength"/>) at the
    /// 50 mM Primer3 probe conditions spans 71.25–85.35 °C over the FISH window (200–500 nt, G+C 35–65 %),
    /// 70.30–82.50 °C over the Northern window (100–300 nt, G+C 40–60 %) and 70.32–85.35 °C over the Southern window
    /// (150–500 nt, G+C 35–65 %); their Tm windows are library conventions (no published Tm window for these long
    /// probes was found — set them for the protocol's hybridization buffer). The qPCR window (68–70 °C, Primer Express)
    /// is reachable for 20–30-nt probes with G+C 30–80 % (nearest-neighbour Tm at the Primer3 probe conditions).
    /// </remarks>
    public static class Defaults
    {
        /// <summary>OligoArray 2.0 microarray Tm conditions: [Na⁺] = 1 M, oligo 1 µM (Rouillard, Zuker &amp; Gulari 2003).</summary>
        public const double OligoArrayMonovalentMillimolar = 1000.0;

        /// <summary>OligoArray 2.0 oligo concentration for the Tm, 1 µM = 1000 nM (Rouillard, Zuker &amp; Gulari 2003).</summary>
        public const double OligoArrayDnaConcentrationNanomolar = 1000.0;

        /// <summary>
        /// Long-oligo microarray probe preset (50–60 nt: Kane et al. 2000 50-mers, Agilent 60-mers; G+C 40–60 %) with
        /// OligoArray 2.0's Tm (Rouillard, Zuker &amp; Gulari 2003, NAR 31:3057): nearest-neighbour Tm over the whole
        /// oligo at [Na⁺] = 1 M and 1 µM oligo, Tm window 82–90 °C (the paper's microarray design setting). Here the
        /// nearest-neighbour Tm is Primer3 <c>oligotm</c> (SantaLucia 1998 unified parameters, C/4 for a
        /// non-self-complementary oligo) with <see cref="ProbeParameters.MaxNearestNeighborLength"/> = 60, no Mg²⁺/dNTP.
        /// The previous window 75–85 °C at the 50 mM Primer3 conditions was unreachable for 50–60-mers with G+C ≤ 60 %
        /// (<c>long_seq_tm</c> maximum 74.50 °C). The structure screen runs at the same conditions.
        /// </summary>
        public static ProbeParameters Microarray => new(
            MinLength: 50, MaxLength: 60,
            MinTm: 82, MaxTm: 90,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3)
        {
            MonovalentMillimolar = OligoArrayMonovalentMillimolar,
            DnaConcentrationNanomolar = OligoArrayDnaConcentrationNanomolar,
            DivalentMillimolar = 0,
            DntpMillimolar = 0,
            MaxNearestNeighborLength = 60
        };

        /// <summary>FISH probe preset (200–500 nt, G+C 35–65 %, Tm 70–90 °C at the Primer3 probe conditions; library
        /// convention, reachable: <c>long_seq_tm</c> 71.25–85.35 °C over the window).</summary>
        public static ProbeParameters FISH => new(
            MinLength: 200, MaxLength: 500,
            MinTm: 70, MaxTm: 90,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 8,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);

        /// <summary>Northern-blot probe preset (100–300 nt, G+C 40–60 %, Tm 65–80 °C at the Primer3 probe conditions;
        /// library convention, reachable: <c>long_seq_tm</c> 70.30–82.50 °C over the window).</summary>
        public static ProbeParameters NorthernBlot => new(
            MinLength: 100, MaxLength: 300,
            MinTm: 65, MaxTm: 80,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 6,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.35);

        /// <summary>
        /// Hydrolysis (TaqMan) qPCR probe preset: Tm 68–70 °C and G+C 30–80 % (Applied Biosystems
        /// Primer Express probe guidelines: "Tm should be 68 to 70 °C", "Keep G-C content in the 30-80%
        /// range"; the Primer Express Tm scale itself is proprietary — here the Tm is Primer3 <c>seqtm</c>
        /// at <see cref="ProbeParameters.MonovalentMillimolar"/> etc.); length 20–30 nt within the
        /// 18–30 nt hydrolysis-probe range (IDT). Use <see cref="EvaluateTaqManProbe"/> for the
        /// chemistry-specific rules (no 5′ G, more C than G, no GGGG, probe Tm ≥ primer Tm + 10 °C).
        /// </summary>
        public static ProbeParameters qPCR => new(
            MinLength: 20, MaxLength: 30,
            MinTm: 68, MaxTm: 70,
            MinGc: 0.30, MaxGc: 0.80,
            MaxHomopolymer: 4,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.25);

        /// <summary>Southern-blot probe preset (150–500 nt, G+C 35–65 %, Tm 65–75 °C at the Primer3 probe conditions;
        /// library convention, reachable: <c>long_seq_tm</c> 70.32–85.35 °C over the window).</summary>
        public static ProbeParameters SouthernBlot => new(
            MinLength: 150, MaxLength: 500,
            MinTm: 65, MaxTm: 75,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 7,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);
    }

    /// <summary>
    /// Designed probe.
    /// </summary>
    public readonly record struct Probe(
        string Sequence,
        int Start,
        int End,
        double Tm,
        double GcContent,
        double Score,
        ProbeType Type,
        IReadOnlyList<string> Warnings)
    {
        /// <summary>
        /// Primer3 internal-oligo penalty (<c>p_obj_fn</c>, PRIMER_INTERNAL_n_PENALTY) of the probe when
        /// <see cref="ProbeParameters.Ranking"/> = <see cref="ProbeRanking.Primer3Penalty"/> and its Tm is computable;
        /// otherwise <c>null</c>. Lower is better. <see cref="Score"/> stays the additive library score.
        /// </summary>
        public double? Primer3Penalty { get; init; }
    }

    /// <summary>
    /// Probe set for tiling.
    /// </summary>
    public readonly record struct TilingProbeSet(
        IReadOnlyList<Probe> Probes,
        int Coverage,
        double MeanTm,
        double TmRange);

    /// <summary>
    /// Probe validation result (<see cref="ValidateProbe"/>).
    /// </summary>
    /// <param name="IsValid">True when no issue was recorded (at most one reference site meeting the Kane et al. (2000)
    /// criteria, no self-structure flag, no cross-hybridizing non-target).</param>
    /// <param name="SpecificityScore">Library-defined uniqueness score 1/N over the N ungapped candidate binding
    /// sites (0 when there is none); a library convention, not a published metric, and not used by
    /// <paramref name="IsValid"/> — see <see cref="ValidateProbe"/>.</param>
    /// <param name="OffTargetHits">Total ungapped k-mismatch hits across the references (the intended site
    /// included).</param>
    /// <param name="SelfComplementarity">Position-wise fold-back fraction (fraction of positions i with
    /// s[i] = revcomp(s)[i]); a library metric, reported for compatibility and not used by any screen (the
    /// self-dimer criterion is ntthal or, in the fallback, Primer3 alignment-mode <see cref="SelfAny"/> /
    /// <see cref="SelfEnd"/>).</param>
    /// <param name="HasSecondaryStructure">Hairpin flag: ntthal hairpin Tm &gt; MaxStructureTm (thermodynamic
    /// screen) or the inverted-repeat stem screen (fallback).</param>
    /// <param name="Issues">Recorded validation issues.</param>
    public readonly record struct ProbeValidation(
        bool IsValid,
        double SpecificityScore,
        int OffTargetHits,
        double SelfComplementarity,
        bool HasSecondaryStructure,
        IReadOnlyList<string> Issues)
    {
        /// <summary>True when the Primer3 thermodynamic self-structure screen was applied (A/C/G/T, ≤ ThermodynamicScreenMaxLength nt
        /// probe with <see cref="ProbeStructureScreen.Thermodynamic"/>); false for the fallback screen (Primer3
        /// alignment-mode self_any / self_end + inverted-repeat hairpin stem).</summary>
        public bool ThermodynamicScreen { get; init; }

        /// <summary>ntthal self-dimer (ANY) Tm in °C (Primer3 SELF_ANY_TH); null when the fallback screen was used.</summary>
        public double? SelfDimerTm { get; init; }

        /// <summary>ntthal 3′ self-dimer (END1) Tm in °C (Primer3 SELF_END_TH); null when the fallback screen was used.</summary>
        public double? SelfEndDimerTm { get; init; }

        /// <summary>ntthal hairpin Tm in °C (Primer3 HAIRPIN_TH); null when the fallback screen was used.</summary>
        public double? HairpinTm { get; init; }

        /// <summary>Primer3 alignment-mode internal-oligo <c>self_any</c> (PRIMER_INTERNAL_SELF_ANY with
        /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0; <see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>),
        /// reported for every non-empty probe (any length); the self-dimer criterion of the fallback screen.</summary>
        public double? SelfAny { get; init; }

        /// <summary>Primer3 alignment-mode internal-oligo <c>self_end</c>
        /// (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>), reported for every non-empty probe.</summary>
        public double? SelfEnd { get; init; }

        /// <summary>Kane et al. (2000) cross-hybridization assessment of every supplied non-target sequence/strand
        /// (<see cref="AssessCrossHybridization"/>); empty when no non-target sequences were supplied.</summary>
        public IReadOnlyList<CrossHybridizationAssessment> CrossHybridization { get; init; }
            = Array.Empty<CrossHybridizationAssessment>();

        /// <summary>Number of the <see cref="OffTargetHits"/> ungapped sites that meet a Kane et al. (2000)
        /// cross-hybridization criterion (identity (L − mismatches) / L &gt; the identity threshold, default 0.75, or a
        /// run of identical positions longer than the contiguous threshold, default 15 nt); more than one such site
        /// records the off-target issue (audit round 3, A3-12).</summary>
        public int CrossHybridizingHits { get; init; }
    }

    /// <summary>
    /// Kane et al. (2000) cross-hybridization assessment of a probe against one strand of one non-target
    /// sequence (<see cref="AssessCrossHybridization"/>).
    /// </summary>
    /// <param name="NonTargetIndex">Zero-based index of the non-target sequence (in the supplied order).</param>
    /// <param name="ReverseComplementStrand">False: the non-target as given; true: its reverse complement.</param>
    /// <param name="Identity">Overall identity = identical columns of the best local (Smith–Waterman–Gotoh)
    /// alignment of the probe with this strand ÷ probe length, in [0, 1] (Kane: "similar over the 50 base target").</param>
    /// <param name="IdenticalColumns">Identical aligned columns of that alignment.</param>
    /// <param name="AlignmentScore">Raw score of that alignment (0 when nothing aligns).</param>
    /// <param name="LongestContiguousMatch">Length of the longest stretch identical between the probe and this
    /// strand (longest common substring).</param>
    /// <param name="ExceedsIdentityThreshold">Identity strictly above the identity threshold (Kane: &gt; 75 %).</param>
    /// <param name="ExceedsContiguousThreshold">Contiguous match strictly longer than the threshold (Kane: &gt; 15 nt).</param>
    public readonly record struct CrossHybridizationAssessment(
        int NonTargetIndex,
        bool ReverseComplementStrand,
        double Identity,
        int IdenticalColumns,
        int AlignmentScore,
        int LongestContiguousMatch,
        bool ExceedsIdentityThreshold,
        bool ExceedsContiguousThreshold)
    {
        /// <summary>True when a Kane criterion or the optional duplex-Tm threshold is met (the probe may
        /// cross-hybridize with this strand).</summary>
        public bool CrossHybridizes => ExceedsIdentityThreshold || ExceedsContiguousThreshold || ExceedsDuplexTmThreshold;

        /// <summary>True when a duplex-Tm threshold was given and <see cref="DuplexTm"/> exceeds it (OligoArray 2.0:
        /// a cross-hybridization with Tm above the user's specificity threshold makes the probe non-specific).</summary>
        public bool ExceedsDuplexTmThreshold { get; init; }

        /// <summary>0-based inclusive start of the aligned site in the assessed strand (−1 when nothing aligns).</summary>
        public int SiteStart { get; init; }

        /// <summary>0-based inclusive end of the aligned site in the assessed strand (−1 when nothing aligns).</summary>
        public int SiteEnd { get; init; }

        /// <summary>
        /// Thermodynamic stability of the probe on this off-target site: ntthal duplex (THAL_ANY) Tm in °C of the
        /// probe with the strand complementary to the aligned site (primer3-py <c>calc_heterodimer(probe,
        /// revcomp(site))</c>; 0 when no duplex forms) at the assessment conditions — the duplex-Tm cross-hybridization
        /// check of OligoArray (Rouillard et al. 2003). Null when both the probe and the site are longer than the
        /// conditions' <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (thal.c <c>THAL_MAX_ALIGN</c>, default
        /// 60: thal_check_errors refuses only when both strands are longer — primer3-py raises), either is longer than
        /// 10 000 nt (THAL_MAX_SEQ), probe or site contain a non-ACGT base, or nothing aligns.
        /// </summary>
        public double? DuplexTm { get; init; }
    }

    /// <summary>
    /// A single gapped (Smith-Waterman) hit of a probe against a reference sequence.
    /// Unlike the ungapped Hamming scan used by <see cref="ValidateProbe"/>, a hit may span an
    /// insertion or deletion (see <see cref="HasGaps"/>), so off-target sites reachable only through
    /// an indel are detected (Altschul et al. 1990; Smith &amp; Waterman 1981).
    /// </summary>
    /// <param name="ReferenceIndex">Zero-based index of the reference sequence (in the supplied order) the hit was found in.</param>
    /// <param name="Start">Zero-based start position of the aligned region within the reference.</param>
    /// <param name="End">Zero-based end position (inclusive) of the aligned region within the reference.</param>
    /// <param name="Identity">
    /// Fraction of identical aligned columns over the probe length
    /// (identical columns / probe length), in [0, 1]. Computed exactly as in the library's
    /// reused ANI identity convention (Goris et al. 2007; pyani <c>ani_pid = ani_alnids / qlen</c>).
    /// </param>
    /// <param name="Coverage">Fraction of the probe covered by ungapped aligned columns (ungapped columns / probe length), in [0, 1].</param>
    /// <param name="HasGaps">True when the alignment contains at least one gap (insertion or deletion) — i.e. the hit is reachable only with an indel.</param>
    /// <param name="AlignedProbe">The probe side of the local alignment (may contain '-').</param>
    /// <param name="AlignedReference">The reference side of the local alignment (may contain '-').</param>
    public readonly record struct GappedProbeHit(
        int ReferenceIndex,
        int Start,
        int End,
        double Identity,
        double Coverage,
        bool HasGaps,
        string AlignedProbe,
        string AlignedReference);

    /// <summary>
    /// Result of the opt-in gapped (Smith-Waterman) off-target scan, separating the intended
    /// on-target match from genuine off-target hits.
    /// </summary>
    /// <param name="OnTargetHits">
    /// Perfect, ungapped, full-coverage exact matches (identity = 1.0, coverage = 1.0, no gaps).
    /// The first such hit is the probe's intended on-target binding site; any additional perfect
    /// exact matches are reported here too (a probe matching several identical sites is itself a
    /// specificity concern) but the first one is excluded from <see cref="OffTargetHits"/>.
    /// </param>
    /// <param name="OffTargetHits">
    /// Genuine off-target hits: every hit at or above the identity threshold that is NOT the single
    /// intended on-target site. This includes imperfect (mismatched) hits, indel-containing hits, and
    /// any extra perfect repeats. This is the corrected count that no longer pools the on-target match
    /// with off-targets (cf. <see cref="ProbeValidation.OffTargetHits"/>).
    /// </param>
    /// <param name="MinIdentity">The identity threshold used to call a hit (see <see cref="ScanOffTargetsGapped"/>).</param>
    public readonly record struct GappedSpecificityResult(
        IReadOnlyList<GappedProbeHit> OnTargetHits,
        IReadOnlyList<GappedProbeHit> OffTargetHits,
        double MinIdentity)
    {
        /// <summary>Number of genuine off-target hits (excludes the intended on-target match).</summary>
        public int OffTargetCount => OffTargetHits.Count;

        /// <summary>True when exactly one on-target site and no off-target hits were found.</summary>
        public bool IsSpecific => OnTargetHits.Count == 1 && OffTargetHits.Count == 0;
    }

    /// <summary>
    /// Karlin–Altschul statistics of an off-target alignment hit: the raw alignment score's
    /// statistical significance expressed as a bit score and an expectation (E) value.
    /// </summary>
    /// <remarks>
    /// Karlin &amp; Altschul (1990, PNAS 87:2264); Altschul et al. (1990, J Mol Biol 215:403).
    /// </remarks>
    /// <param name="RawScore">The raw alignment score S (sum of substitution/gap scores).</param>
    /// <param name="Lambda">
    /// The Karlin–Altschul scale parameter λ — the unique positive root of
    /// Σ_{i,j} p_i p_j e^{λ s_ij} = 1 for the scoring scheme and base frequencies.
    /// </param>
    /// <param name="K">The Karlin–Altschul search-space scale parameter K (computed for the scoring scheme unless supplied by the caller).</param>
    /// <param name="BitScore">The normalized bit score S' = (λS − ln K) / ln 2.</param>
    /// <param name="EValue">
    /// The expected number of distinct alignments scoring ≥ S by chance: E = K·m·n·e^{−λS} = m·n·2^{−S'}.
    /// </param>
    /// <param name="QueryLength">The query (probe) length m used in the search space.</param>
    /// <param name="DatabaseLength">The database (reference) length n used in the search space.</param>
    public readonly record struct KarlinAltschulStatistics(
        double RawScore,
        double Lambda,
        double K,
        double BitScore,
        double EValue,
        int QueryLength,
        long DatabaseLength);

    /// <summary>
    /// Probe type.
    /// </summary>
    public enum ProbeType
    {
        Standard,
        Tiling,
        Antisense,
        LNA, // Locked Nucleic Acid
        MolecularBeacon
    }

    /// <summary>
    /// Result of evaluating a candidate hydrolysis (TaqMan) probe against the
    /// Applied Biosystems / Thermo Fisher TaqMan probe-design guidelines.
    /// Each boolean records whether one published rule is satisfied; <see cref="PassesAll"/>
    /// is the conjunction of every rule.
    /// </summary>
    /// <remarks>
    /// Sources (retrieved 2026-06-24):
    /// PREMIER Biosoft, "TaqMan probe design tips" (http://www.premierbiosoft.com/tech_notes/TaqMan.html);
    /// Thermo Fisher / Applied Biosystems "Designing a TaqMan Gene Expression Assay" and
    /// "TaqMan MGB Probe and Primer Sets" application notes.
    /// </remarks>
    public readonly record struct TaqManProbeEvaluation(
        string Sequence,
        bool NoGuanineAt5Prime,
        bool MoreCytosineThanGuanine,
        bool NoRunOfFourOrMoreG,
        bool GcContentInRange,
        bool LengthInRange,
        bool ProbeTmAbovePrimer,
        double Tm,
        double GcContent,
        int CytosineCount,
        int GuanineCount,
        bool PassesAll,
        IReadOnlyList<string> Violations);

    #endregion

    #region TaqMan (hydrolysis-probe) design — opt-in

    // --- Published TaqMan probe-design thresholds (Applied Biosystems / Thermo Fisher; PREMIER Biosoft) ---

    // Probe length range. PREMIER Biosoft: "TaqMan probes consist of a 18-22 bp oligonucleotide probe".
    // (IDT / Thermo state 18-30; we use the tighter ABI/PREMIER 18-22 default, configurable.)
    private const int TaqManMinLength = 18;
    private const int TaqManMaxLength = 22;

    // G+C content range. PREMIER Biosoft: "The G+C content should ideally be 30-80%".
    private const double TaqManMinGc = 0.30;
    private const double TaqManMaxGc = 0.80;

    // No runs of identical nucleotides, "especially four or more consecutive Gs" (PREMIER Biosoft / ABI).
    private const int TaqManMaxGuanineRun = 4;

    // Probe Tm should be ~10 °C higher than the primer Tm so the probe binds before Taq extends.
    // PREMIER Biosoft / Thermo Fisher: "TaqMan probe Tm should be 10 °C higher than the Primer Tm".
    private const double TaqManProbeTmDeltaAbovePrimer = 10.0;

    /// <summary>
    /// Evaluates a single candidate probe sequence against the published TaqMan
    /// (5'-nuclease hydrolysis probe) design guidelines. This is an opt-in chemistry-specific
    /// check; the generic <see cref="DesignProbes(string, ProbeParameters?, int)"/> designer is unchanged.
    /// </summary>
    /// <param name="probeSequence">The candidate probe sequence (5'→3').</param>
    /// <param name="primerTm">
    /// Melting temperature (°C) of the amplification primers. The probe Tm must be at least
    /// <c>primerTm + 10 °C</c> (Applied Biosystems / Thermo Fisher). Pass <see langword="null"/>
    /// to skip the probe-Tm-vs-primer gate (it is then reported as satisfied).
    /// </param>
    /// <param name="minLength">Minimum probe length (default 18, per ABI/PREMIER Biosoft).</param>
    /// <param name="maxLength">Maximum probe length (default 22, per ABI/PREMIER Biosoft).</param>
    /// <param name="dnaConcentrationNanomolar">Probe concentration (nM) for the probe Tm (default 50, Primer3 PRIMER_INTERNAL_DNA_CONC).</param>
    /// <param name="monovalentMillimolar">Monovalent cations (mM) for the probe Tm (default 50, PRIMER_INTERNAL_SALT_MONOVALENT).</param>
    /// <param name="divalentMillimolar">Mg²⁺ (mM) for the probe Tm (default 0, PRIMER_INTERNAL_SALT_DIVALENT).</param>
    /// <param name="dntpMillimolar">dNTP (mM) for the probe Tm (default 0, PRIMER_INTERNAL_DNTP_CONC).</param>
    /// <returns>A <see cref="TaqManProbeEvaluation"/> recording each rule outcome. <see cref="TaqManProbeEvaluation.Tm"/>
    /// is the Primer3 <c>seqtm</c> probe Tm (SantaLucia 1998 nearest-neighbour + SantaLucia salt correction,
    /// <see cref="PrimerDesigner.CalculateMeltingTemperaturePrimer3"/>) at the given conditions — compare it with a
    /// primer Tm computed on the same scale (e.g. <see cref="PrimerDesigner.CalculateMeltingTemperaturePrimer3"/>);
    /// it is <c>NaN</c> when the probe has fewer than 2 bases or a non-ACGT base, and the Tm gate then fails
    /// whenever a primer Tm is supplied.</returns>
    public static TaqManProbeEvaluation EvaluateTaqManProbe(
        string probeSequence,
        double? primerTm = null,
        int minLength = TaqManMinLength,
        int maxLength = TaqManMaxLength,
        double dnaConcentrationNanomolar = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar,
        double monovalentMillimolar = PrimerDesigner.Primer3InternalMonovalentMillimolar,
        double divalentMillimolar = PrimerDesigner.Primer3InternalDivalentMillimolar,
        double dntpMillimolar = PrimerDesigner.Primer3InternalDntpMillimolar)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);

        string seq = probeSequence.ToUpperInvariant();
        var violations = new List<string>();

        // Rule 1: no G at the 5' end (a 5' G adjacent to the reporter dye quenches
        // reporter fluorescence even after cleavage).
        bool noGuanineAt5Prime = seq.Length > 0 && seq[0] != 'G';
        if (!noGuanineAt5Prime)
            violations.Add("5' end is a guanine (quenches the reporter dye even after cleavage)");

        // Rule 2: more Cs than Gs in the probe sequence.
        int cCount = seq.Count(c => c == 'C');
        int gCount = seq.Count(c => c == 'G');
        bool moreCThanG = cCount > gCount;
        if (!moreCThanG)
            violations.Add($"not more C than G (C={cCount}, G={gCount})");

        // Rule 3: no run of four or more consecutive Gs.
        int maxGRun = GetMaxGuanineRunLength(seq);
        bool noRunOf4G = maxGRun < TaqManMaxGuanineRun;
        if (!noRunOf4G)
            violations.Add($"run of {maxGRun} consecutive Gs (>= {TaqManMaxGuanineRun})");

        // Rule 4: G+C content within 30-80%.
        double gc = seq.CalculateGcFractionFast();
        bool gcInRange = gc >= TaqManMinGc && gc <= TaqManMaxGc;
        if (!gcInRange)
            violations.Add($"G+C content {gc:P0} outside {TaqManMinGc:P0}-{TaqManMaxGc:P0}");

        // Rule 5: probe length within range.
        bool lengthInRange = seq.Length >= minLength && seq.Length <= maxLength;
        if (!lengthInRange)
            violations.Add($"length {seq.Length} outside {minLength}-{maxLength} nt");

        // Rule 6: probe Tm at least ~10 °C above the primer Tm (when a primer Tm is supplied).
        double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(
            seq, dnaConcentrationNanomolar, monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        bool probeTmAbovePrimer = !primerTm.HasValue || tm >= primerTm.Value + TaqManProbeTmDeltaAbovePrimer;
        if (!probeTmAbovePrimer)
            violations.Add(
                $"probe Tm {tm:F1} °C is not >= primer Tm {primerTm!.Value:F1} + {TaqManProbeTmDeltaAbovePrimer:F0} °C");

        bool passesAll = noGuanineAt5Prime && moreCThanG && noRunOf4G
            && gcInRange && lengthInRange && probeTmAbovePrimer;

        return new TaqManProbeEvaluation(
            seq,
            noGuanineAt5Prime,
            moreCThanG,
            noRunOf4G,
            gcInRange,
            lengthInRange,
            probeTmAbovePrimer,
            tm,
            gc,
            cCount,
            gCount,
            passesAll,
            violations);
    }

    /// <summary>
    /// Chooses the better TaqMan probe strand. Per Applied Biosystems / Thermo Fisher,
    /// the probe should be designed on the strand with <b>more Cs than Gs</b>; if a guanine
    /// occurs at the 5' end of one strand, the complement (antisense) strand should be used.
    /// Returns the sense (given) strand or its reverse complement, whichever better satisfies the rules.
    /// </summary>
    /// <param name="senseStrand">The candidate probe sequence on the sense strand (5'→3').</param>
    /// <param name="primerTm">Optional primer Tm for the probe-Tm gate (see <see cref="EvaluateTaqManProbe"/>).</param>
    /// <returns>
    /// A tuple of the chosen probe sequence (5'→3'), whether it is the reverse-complement
    /// (antisense) strand, and the evaluation of the chosen strand.
    /// </returns>
    public static (string Probe, bool IsReverseComplement, TaqManProbeEvaluation Evaluation) SelectTaqManStrand(
        string senseStrand,
        double? primerTm = null)
    {
        ArgumentNullException.ThrowIfNull(senseStrand);

        string sense = senseStrand.ToUpperInvariant();
        string antisense = DnaSequence.GetReverseComplementString(sense);

        var senseEval = EvaluateTaqManProbe(sense, primerTm);
        var antisenseEval = EvaluateTaqManProbe(antisense, primerTm);

        // Prefer a strand that passes all rules; otherwise prefer the one satisfying the two
        // hard reporter-dye rules (no 5'-G, then more C than G); finally fall back to the sense strand.
        if (senseEval.PassesAll && !antisenseEval.PassesAll)
            return (sense, false, senseEval);
        if (antisenseEval.PassesAll && !senseEval.PassesAll)
            return (antisense, true, antisenseEval);

        int senseRank = RankTaqManStrand(senseEval);
        int antisenseRank = RankTaqManStrand(antisenseEval);

        return antisenseRank > senseRank
            ? (antisense, true, antisenseEval)
            : (sense, false, senseEval);
    }

    // Ranks a strand: the no-5'-G rule then the more-C-than-G rule are the chemistry-critical
    // reporter-dye constraints, weighted above the remaining quality rules.
    private static int RankTaqManStrand(TaqManProbeEvaluation e)
    {
        int rank = 0;
        if (e.NoGuanineAt5Prime) rank += 4;
        if (e.MoreCytosineThanGuanine) rank += 2;
        if (e.NoRunOfFourOrMoreG) rank += 1;
        if (e.GcContentInRange) rank += 1;
        if (e.ProbeTmAbovePrimer) rank += 1;
        return rank;
    }

    private static int GetMaxGuanineRunLength(string sequence)
    {
        int maxRun = 0;
        int currentRun = 0;

        foreach (char c in sequence)
        {
            if (c == 'G')
            {
                currentRun++;
                maxRun = Math.Max(maxRun, currentRun);
            }
            else
            {
                currentRun = 0;
            }
        }

        return maxRun;
    }

    #endregion

    #region MGB (minor-groove binder) design rules — opt-in (qualitative; quantitative ΔTm is a residual)

    // --- Citable 3'-MGB probe-design rules (Kutyavin et al. 2000, Nucleic Acids Res 28(2):655-661) ---
    //
    // Kutyavin et al. (2000) established that conjugating a minor-groove binder (MGB) to the 3' end
    // of a DNA probe greatly stabilises the duplex, so MGB probes are designed SHORTER than
    // unmodified probes: "for MGB probes this length variation is narrowed to a range of 12-20mers"
    // (a 12mer MGB has ~the same Tm as a 27mer unmodified probe). The MGB is attached at the 3' end
    // ("3'-MGB-ODNs are easier to prepare … MGB-modified solid supports and automated DNA synthesis
    // can be used"). The QUANTITATIVE MGB ΔTm is empirical with no published closed-form model
    // (the stabilisation "varies by sequence"; A+T-rich MGB sites gain more than G+C-rich ones), so
    // only these qualitative DESIGN rules are implemented here; the quantitative MGB ΔTm is left as
    // an honest residual (see docs/Validation/LIMITATIONS.md).
    // Source (retrieved 2026-06-24): Kutyavin IV et al. (2000) Nucleic Acids Res 28(2):655-661,
    //   https://doi.org/10.1093/nar/28.2.655.

    // Kutyavin (2000): MGB-probe length range narrowed to 12-20mers.
    private const int MgbMinLength = 12;
    private const int MgbMaxLength = 20;

    /// <summary>
    /// Result of checking a candidate probe against the citable 3'-MGB (minor-groove binder)
    /// design rules of Kutyavin et al. (2000). These are <b>qualitative</b> design-rule checks
    /// only; the quantitative MGB ΔTm is empirical (no published formula) and is not computed.
    /// </summary>
    /// <param name="Sequence">The (upper-cased) probe sequence checked.</param>
    /// <param name="LengthInMgbRange">True when the probe length is within the MGB 12–20mer window.</param>
    /// <param name="Length">The probe length in nucleotides.</param>
    /// <param name="MgbAttachmentEnd">The recommended MGB attachment end (always 3', per Kutyavin 2000).</param>
    /// <param name="Guidance">Human-readable design guidance / any rule violations.</param>
    public readonly record struct MgbProbeDesign(
        string Sequence,
        bool LengthInMgbRange,
        int Length,
        string MgbAttachmentEnd,
        IReadOnlyList<string> Guidance);

    /// <summary>
    /// Evaluates a candidate probe against the citable 3'-MGB (minor-groove binder) probe-design
    /// rules of Kutyavin et al. (2000): the MGB is attached at the <b>3' end</b>, and MGB probes are
    /// designed <b>shorter</b> (12–20mer) than unmodified probes. This is an opt-in, qualitative
    /// design-rule check; the generic <see cref="DesignProbes(string, ProbeParameters?, int)"/>
    /// designer and all defaults are unchanged. The <b>quantitative</b> MGB ΔTm is empirical (no
    /// published closed-form model in Kutyavin 2000) and is deliberately NOT computed.
    /// </summary>
    /// <param name="probeSequence">The candidate probe sequence (5'→3').</param>
    /// <returns>An <see cref="MgbProbeDesign"/> recording the length-window outcome and 3'-MGB guidance.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="probeSequence"/> is null.</exception>
    public static MgbProbeDesign EvaluateMgbProbeDesign(string probeSequence)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);

        // Qualitative MGB rules only — the quantitative MGB ΔTm is not computed (Kutyavin 2000 is
        // empirical, no closed form). Strict mode throws; Moderate/Permissive allow this narrower result.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("PROBE-DESIGN-001");

        string seq = probeSequence.ToUpperInvariant();
        var guidance = new List<string>();

        bool lengthInRange = seq.Length >= MgbMinLength && seq.Length <= MgbMaxLength;
        if (!lengthInRange)
            guidance.Add(
                $"length {seq.Length} outside the MGB {MgbMinLength}-{MgbMaxLength}mer window (Kutyavin 2000)");

        // 3'-MGB placement guidance (always applies).
        guidance.Add("attach the minor-groove binder at the 3' end (Kutyavin 2000)");

        return new MgbProbeDesign(
            seq,
            lengthInRange,
            seq.Length,
            MgbAttachmentEnd: "3'",
            guidance);
    }

    #endregion

    #region Probe Design

    /// <summary>
    /// Designs probes for a target sequence: every window of <see cref="ProbeParameters.MinLength"/>..<see cref="ProbeParameters.MaxLength"/>
    /// (GC within ±0.1 of the window) is scored 1.0 minus penalties — GC outside range 0.3, Tm (Primer3 <c>seqtm</c> at the
    /// parameters' conditions) outside range 0.3, homopolymer run above the maximum 0.2, self-complementarity 0.2,
    /// secondary structure 0.15 (only with <see cref="ProbeParameters.AvoidSecondaryStructure"/>), di-/trinucleotide
    /// microsatellite ≥ 4 copies 0.1, G/C at each terminus 0.02 (see <see cref="ProbeParameters"/> for the self-structure
    /// screens). Probes with score &gt; 0 are returned by score descending, ties by (length, start) ascending. This
    /// additive score (the default, <see cref="ProbeRanking.AdditiveScore"/>) is a library ranking heuristic: its penalty
    /// values have no published source. With <see cref="ProbeParameters.Ranking"/> = <see cref="ProbeRanking.Primer3Penalty"/>
    /// the same candidates (score &gt; 0) are instead ranked by Primer3's internal-oligo objective <c>p_obj_fn</c>
    /// (PRIMER_INTERNAL_n_PENALTY = |Tm − <see cref="ProbeParameters.OptTm"/>| + |length − <see cref="ProbeParameters.OptLength"/>|
    /// with Primer3's default internal-oligo weights, reported in <see cref="Probe.Primer3Penalty"/>) in Primer3's
    /// <c>primer_rec_comp</c> order (penalty ascending, start descending, length ascending; a probe whose Tm is not
    /// computable ranks last) — the sourced ranking. <see cref="DesignProbesPrimer3"/> is Primer3's complete picker
    /// (Primer3's acceptance limits as well as its ranking).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">With the <see cref="ProbeRanking.Primer3Penalty"/> ranking:
    /// <see cref="ProbeParameters.OptLength"/> outside [MinLength, MaxLength] or <see cref="ProbeParameters.OptTm"/> outside
    /// [MinTm, MaxTm] (Primer3 <c>_pr_data_control</c>); an undefined <see cref="ProbeParameters.Ranking"/> value;
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> outside 60–10 000.</exception>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 10)
    {
        var param = parameters ?? Defaults.Microarray;
        ValidateRanking(param, nameof(parameters));
        ValidateThermodynamicScreenMaxLength(param, nameof(parameters));
        return DesignProbesIterator(targetSequence, param, maxProbes);
    }

    private static IEnumerable<Probe> DesignProbesIterator(string targetSequence, ProbeParameters param, int maxProbes)
    {
        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();

        // Use optimized evaluation with prefix sums for O(1) GC lookup
        var candidates = DesignProbesOptimized(targetSequence, param, maxProbes);

        foreach (var probe in candidates)
        {
            yield return probe;
        }
    }

    /// <summary>
    /// Designs probes with a genome-wide specificity check through a suffix tree: the candidates of
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/> are walked in their ranking order and each one's
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> value (1 / occurrences in the index)
    /// is applied — with <paramref name="requireUnique"/> a probe occurring more than once is dropped (the remaining
    /// probes keep their score and order), otherwise <see cref="Probe.Score"/> is multiplied by the specificity and the
    /// probes are re-ranked on the scaled score (stable: equal scores keep the documented tie order; the
    /// <see cref="ProbeRanking.Primer3Penalty"/> order is unaffected because the specificity does not enter the penalty).
    /// <b>Every</b> candidate is considered, lazily: the walk stops as soon as <paramref name="maxProbes"/> probes are
    /// produced, so the common case costs no more than the base design plus O(m) per inspected candidate.
    /// O(n × m) for probe generation + O(m) per specificity check.
    /// </summary>
    /// <param name="targetSequence">Target sequence to design probes for.</param>
    /// <param name="genomeIndex">Pre-built suffix tree index for the genome (enables O(m) specificity lookup).</param>
    /// <param name="parameters">Probe design parameters.</param>
    /// <param name="maxProbes">Maximum number of probes to return (none when ≤ 0).</param>
    /// <param name="requireUnique">If true, only return probes unique in the genome.</param>
    /// <param name="bothStrands">Count the probe's reverse-complement occurrences in the index too (passed to
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/>: a probe whose reverse complement
    /// occurs in a double-stranded genome binds the other strand there; a reverse-palindromic probe is counted once).
    /// Default false: the indexed strand only (the behaviour before audit round 4, B07 F60) — pass true when the index
    /// holds one strand of a double-stranded genome.</param>
    /// <remarks>Candidates come in the order of <see cref="ProbeParameters.Ranking"/> (see
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/>), after the specificity scaling when
    /// <paramref name="requireUnique"/> is false. Before audit round 4 (B07 F59) only the top
    /// <paramref name="maxProbes"/> × 5 candidates were inspected and the scaled scores were not re-ranked.</remarks>
    /// <exception cref="ArgumentOutOfRangeException">Invalid Primer3-ranking optima (see
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/>).</exception>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters? parameters = null,
        int maxProbes = 10,
        bool requireUnique = true,
        bool bothStrands = false)
    {
        var param = parameters ?? Defaults.Microarray;
        ValidateRanking(param, nameof(parameters));
        ValidateThermodynamicScreenMaxLength(param, nameof(parameters));
        return DesignProbesIterator(targetSequence, genomeIndex, param, maxProbes, requireUnique, bothStrands);
    }

    private static IEnumerable<Probe> DesignProbesIterator(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters param,
        int maxProbes,
        bool requireUnique,
        bool bothStrands)
    {
        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength || maxProbes <= 0)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();

        // Every candidate is considered, lazily in ranking order, until maxProbes survive (audit round 4, A4-1, F59):
        // requireUnique drops a probe with specificity < 1 (its score is unchanged, so the order is the base ranking);
        // otherwise the score is scaled by the specificity (never raised) and EnumerateRankedProbes re-ranks on it.
        Func<Probe, Probe?> applySpecificity = requireUnique
            ? probe => CheckSpecificity(probe.Sequence, genomeIndex, bothStrands) < 1.0 ? null : probe
            : probe => probe with { Score = probe.Score * CheckSpecificity(probe.Sequence, genomeIndex, bothStrands) };

        int returned = 0;
        foreach (var probe in EnumerateRankedProbes(targetSequence, param, applySpecificity))
        {
            yield return probe;
            if (++returned >= maxProbes)
                yield break;
        }
    }

    // Cheap per-candidate properties (everything except the self-structure screens). The ranking score
    // is rebuilt from these in the original penalty order, so lazily adding the structure penalties gives
    // bit-identical scores to an eager evaluation.
    private readonly record struct ProbeBase(
        string Sequence,
        int Start,
        int Order,
        double Gc,
        double Tm,
        bool TmComputable,
        bool GcOutside,
        bool TmOutside,
        int Homopolymer,
        bool HasRepeats,
        double PositionPenalty,
        double BaseScore,
        double? Primer3Penalty);

    // Primer3 _pr_data_control checks of the optima used by the Primer3Penalty ranking.
    private static void ValidateRanking(ProbeParameters param, string paramName)
    {
        if (!Enum.IsDefined(param.Ranking))
            throw new ArgumentOutOfRangeException(paramName, $"Undefined probe ranking {param.Ranking}.");
        if (param.Ranking != ProbeRanking.Primer3Penalty)
            return;
        if (param.OptLength > param.MaxLength)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE (Primer3 _pr_data_control): OptLength must not exceed MaxLength.");
        if (param.OptLength < param.MinLength)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE < MIN_SIZE (Primer3 _pr_data_control): OptLength must not be below MinLength.");
        if (!(param.OptTm >= param.MinTm && param.OptTm <= param.MaxTm))
            throw new ArgumentOutOfRangeException(paramName,
                "Optimum internal oligo Tm lower than minimum or higher than maximum (Primer3 _pr_data_control): OptTm must lie in [MinTm, MaxTm].");
    }

    // Primer3 p_obj_fn, internal-oligo branch (OT_INTL) with Primer3's default internal-oligo weights (Tm and size
    // weights 1, every other weight 0) around the parameters' optima; null when the Tm is not computable.
    private static double? ComputePrimer3ProbePenalty(double tm, int length, double gcFraction, ProbeParameters param) =>
        double.IsNaN(tm)
            ? null
            : PrimerDesigner.CalculatePrimer3Penalty(
                new Primer3PenaltyInputs(tm, length, 100.0 * gcFraction),
                PrimerDesigner.DefaultPrimer3Weights,
                new Primer3Optima(param.OptTm, param.OptLength, null));

    // Primer3 primer_rec_comp: penalty (quality) ascending, then start descending, then length ascending.
    private static int ComparePrimer3Rank(double penaltyA, int startA, int lengthA, double penaltyB, int startB, int lengthB)
    {
        int c = penaltyA.CompareTo(penaltyB);
        if (c != 0) return c;
        c = startB.CompareTo(startA);
        return c != 0 ? c : lengthA.CompareTo(lengthB);
    }

    private static ProbeBase EvaluateProbeBase(string sequence, int start, int order, ProbeParameters param, double gc)
    {
        double tmRaw = CalculateProbeTm(sequence, param);
        bool tmComputable = !double.IsNaN(tmRaw);
        double tm = tmComputable ? tmRaw : 0.0;

        bool gcOutside = gc < param.MinGc || gc > param.MaxGc;
        bool tmOutside = !tmComputable || tm < param.MinTm || tm > param.MaxTm;
        int homopolymer = PrimerDesigner.FindLongestHomopolymer(sequence);
        bool repeats = HasSimpleRepeats(sequence);

        double positionPenalty = 0;
        if (sequence.StartsWith('G') || sequence.StartsWith('C'))
            positionPenalty += 0.02;
        if (sequence.EndsWith('G') || sequence.EndsWith('C'))
            positionPenalty += 0.02;

        double score = 1.0;
        if (gcOutside) score -= 0.3;
        if (tmOutside) score -= 0.3;
        if (homopolymer > param.MaxHomopolymer) score -= 0.2;
        if (repeats) score -= 0.1;
        score -= positionPenalty;

        double? primer3Penalty = param.Ranking == ProbeRanking.Primer3Penalty
            ? ComputePrimer3ProbePenalty(tmRaw, sequence.Length, gc, param)
            : null;

        return new ProbeBase(sequence, start, order, gc, tm, tmComputable, gcOutside, tmOutside,
            homopolymer, repeats, positionPenalty, score, primer3Penalty);
    }

    // Adds the self-structure screens to a base evaluation (penalties in the original order:
    // GC, Tm, homopolymer, self-complementarity, secondary structure, repeats, terminal G/C).
    private static Probe? FinishProbe(
        ProbeBase b, ProbeParameters param,
        Dictionary<string, (bool, string?, bool, string?)>? structureCache = null)
    {
        var warnings = new List<string>();
        double score = 1.0;

        if (b.GcOutside)
        {
            score -= 0.3;
            warnings.Add($"GC content {b.Gc:P0} outside range");
        }

        if (b.TmOutside)
        {
            score -= 0.3;
            warnings.Add(b.TmComputable
                ? $"Tm {b.Tm:F1}°C outside range"
                : "Tm not computable: probe contains a non-ACGT base");
        }

        if (b.Homopolymer > param.MaxHomopolymer)
        {
            score -= 0.2;
            warnings.Add($"Homopolymer run of {b.Homopolymer}");
        }

        (bool, string?, bool, string?) screen;
        if (structureCache is null || !structureCache.TryGetValue(b.Sequence, out screen))
        {
            screen = EvaluateSelfStructure(b.Sequence, param);
            structureCache?.Add(b.Sequence, screen);
        }
        var (selfComp, selfCompWarning, structure, structureWarning) = screen;
        if (selfComp)
        {
            score -= 0.2;
            warnings.Add(selfCompWarning!);
        }

        if (param.AvoidSecondaryStructure && structure)
        {
            score -= 0.15;
            warnings.Add(structureWarning!);
        }

        if (b.HasRepeats)
        {
            score -= 0.1;
            warnings.Add("Contains simple repeats");
        }

        score -= b.PositionPenalty;

        if (score <= 0)
            return null;

        return new Probe(
            b.Sequence,
            b.Start,
            b.Start + b.Sequence.Length - 1,
            b.Tm,
            b.Gc,
            Math.Max(0, score),
            ProbeType.Standard,
            warnings)
        {
            Primer3Penalty = b.Primer3Penalty,
        };
    }

    // Self-structure screens (see ProbeParameters remarks). Returns (self-complementarity flag + warning,
    // secondary-structure flag + warning); the caller applies the secondary-structure flag only when
    // AvoidSecondaryStructure is set.
    private static (bool SelfComp, string? SelfCompWarning, bool Structure, string? StructureWarning)
        EvaluateSelfStructure(string sequence, ProbeParameters param)
    {
        if (ComputeThermodynamicSelfStructure(sequence, param) is { } v)
        {
            {
                double max = param.MaxStructureTm;
                bool selfComp = v.SelfAnyTh > max || v.SelfEndTh > max;
                bool hairpin = v.HairpinTh > max;
                return (
                    selfComp,
                    selfComp ? $"Self-dimer Tm {Math.Max(v.SelfAnyTh, v.SelfEndTh):F1}°C exceeds {max:0.##}°C (ntthal)" : null,
                    hairpin,
                    hairpin ? $"Hairpin Tm {v.HairpinTh:F1}°C exceeds {max:0.##}°C (ntthal)" : null);
            }
        }

        var (selfDimer, selfDimerWarning) = AlignmentSelfDimerScreen(sequence, param);
        // The stem flag is used only with AvoidSecondaryStructure (FinishProbe); skip the O(n³) scan otherwise.
        bool stem = param.AvoidSecondaryStructure && HasSecondaryStructurePotential(sequence);
        return (
            selfDimer,
            selfDimerWarning,
            stem,
            stem ? "Potential secondary structure" : null);
    }

    // Fallback self-dimer criterion: Primer3 alignment-mode internal-oligo oligo_compl (dpal self_any, self_end;
    // PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END, default 12.00), no length limit. The decision only needs the
    // threshold comparison (early-exit dpal), so the warning names the exceeded limit; ValidateProbe reports the
    // exact values.
    private static (bool Flag, string? Warning) AlignmentSelfDimerScreen(string sequence, ProbeParameters param)
    {
        if (!PrimerDesigner.ExceedsPrimer3SelfComplementarity(sequence, param.MaxSelfAny, param.MaxSelfEnd, out bool any))
            return (false, null);
        return (true, any
            ? string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: Primer3 self_any exceeds {param.MaxSelfAny:0.00}")
            : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: Primer3 self_end exceeds {param.MaxSelfEnd:0.00}"));
    }

    private static (bool Flag, string? Warning) AlignmentSelfDimerScreen(double selfAny, double selfEnd, ProbeParameters param)
    {
        if (selfAny > param.MaxSelfAny)
            return (true, string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"Self-complementarity: Primer3 self_any {selfAny:0.00} exceeds {param.MaxSelfAny:0.00}"));
        if (selfEnd > param.MaxSelfEnd)
            return (true, string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"Self-complementarity: Primer3 self_end {selfEnd:0.00} exceeds {param.MaxSelfEnd:0.00}"));
        return (false, null);
    }

    // Primer3 thermodynamic self-structure of a probe (ntthal ANY / END1 self-dimer and hairpin Tm at the
    // parameters' conditions), or null when the fallback screen applies (Heuristic screen, longer than
    // ThermodynamicScreenMaxLength = THAL_MAX_ALIGN (default 60), or a non-ACGT base).
    private static PrimerDesigner.Primer3OligoStructure? ComputeThermodynamicSelfStructure(
        string sequence, ProbeParameters param)
    {
        if (param.StructureScreen != ProbeStructureScreen.Thermodynamic || sequence.Length > param.ThermodynamicScreenMaxLength)
            return null;

        return PrimerDesigner.CalculatePrimer3OligoStructure(
            sequence, param.MonovalentMillimolar, param.DivalentMillimolar, param.DntpMillimolar,
            param.DnaConcentrationNanomolar, param.ThermodynamicScreenMaxLength);
    }

    // ProbeParameters.ThermodynamicScreenMaxLength: THAL_MAX_ALIGN override range 60 (Primer3) … 10 000 (THAL_MAX_SEQ).
    private static void ValidateThermodynamicScreenMaxLength(ProbeParameters param, string paramName)
    {
        if (param.ThermodynamicScreenMaxLength < PrimerDesigner.NtthalMaxAlignLength
            || param.ThermodynamicScreenMaxLength > PrimerDesigner.NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(paramName,
                "ThermodynamicScreenMaxLength (THAL_MAX_ALIGN) must be in 60..10000.");
    }

    /// <summary>
    /// The top <paramref name="maxProbes"/> of <see cref="EnumerateRankedProbes"/>: every window of the configured
    /// length range (prefix-sum GC) is ranked; the self-structure screens only lower a score, so they are evaluated
    /// lazily and the scan stops once no remaining candidate can enter the top <paramref name="maxProbes"/>; the
    /// result equals an exhaustive evaluation: probes with score &gt; 0, ordered by score descending, ties by
    /// (length, start) ascending (or by the Primer3 penalty with <see cref="ProbeRanking.Primer3Penalty"/>).
    /// </summary>
    private static List<Probe> DesignProbesOptimized(
        string targetSequence,
        ProbeParameters param,
        int maxProbes) =>
        maxProbes <= 0
            ? new List<Probe>()
            : EnumerateRankedProbes(targetSequence, param).Take(maxProbes).ToList();

    /// <summary>
    /// Lazily yields every candidate (score &gt; 0 after the self-structure screens, then <paramref name="adjust"/>;
    /// a null from <paramref name="adjust"/> drops the probe) in the ranking order of <see cref="ProbeParameters.Ranking"/>:
    /// <see cref="ProbeRanking.AdditiveScore"/> — (adjusted) score descending, ties by enumeration order (length, start)
    /// ascending; <see cref="ProbeRanking.Primer3Penalty"/> — Primer3 <c>primer_rec_comp</c> on the penalty (which
    /// <paramref name="adjust"/> does not change). <paramref name="adjust"/> must never raise a score: the additive order is
    /// then produced by a lazy merge — candidates are finished in descending base-score order (an upper bound of their final
    /// score) and a finished probe is released once it ranks before the bound of every unfinished candidate — so a prefix
    /// of the output costs only the candidates it needs, and the full output equals an exhaustive evaluation + stable sort.
    /// </summary>
    private static IEnumerable<Probe> EnumerateRankedProbes(
        string targetSequence,
        ProbeParameters param,
        Func<Probe, Probe?>? adjust = null)
    {
        var bases = EvaluateCandidateBases(targetSequence, param);
        var structureCache = new Dictionary<string, (bool, string?, bool, string?)>(StringComparer.Ordinal);

        if (param.Ranking == ProbeRanking.Primer3Penalty)
        {
            // The penalty does not depend on the self-structure screens: sort by primer_rec_comp (a non-computable
            // penalty ranks last) and finish lazily.
            bases.Sort((a, b) => ComparePrimer3Rank(
                a.Primer3Penalty ?? double.PositiveInfinity, a.Start, a.Sequence.Length,
                b.Primer3Penalty ?? double.PositiveInfinity, b.Start, b.Sequence.Length));
            foreach (var b in bases)
            {
                if (FinishProbe(b, param, structureCache) is { } p && (adjust is null ? p : adjust(p)) is { } q)
                    yield return q;
            }
            yield break;
        }

        // Rank key: score descending, then enumeration order ascending (= stable sort of the eager scan).
        static int CompareKeys((double Score, int Order) a, (double Score, int Order) b)
        {
            int c = b.Score.CompareTo(a.Score);
            return c != 0 ? c : a.Order.CompareTo(b.Order);
        }

        var pending = new PriorityQueue<Probe, (double Score, int Order)>(
            Comparer<(double Score, int Order)>.Create(CompareKeys));
        // Stable: equal base scores keep the enumeration order, so the bound (BaseScore, Order) is non-decreasing in rank.
        foreach (var b in bases.OrderByDescending(x => x.BaseScore))
        {
            // Structure screens and adjust never raise a score: a pending probe that ranks before this candidate's best
            // possible key ranks before every later candidate too.
            while (pending.TryPeek(out var head, out var key) && CompareKeys(key, (b.BaseScore, b.Order)) < 0)
            {
                pending.Dequeue();
                yield return head;
            }

            if (FinishProbe(b, param, structureCache) is { } p && (adjust is null ? p : adjust(p)) is { } q)
                pending.Enqueue(q, (q.Score, b.Order));
        }

        while (pending.TryDequeue(out var rest, out _))
            yield return rest;
    }

    // Every window of the configured length range (prefix-sum GC) with base score > 0, in enumeration order
    // (length ascending, then start ascending).
    private static List<ProbeBase> EvaluateCandidateBases(string targetSequence, ProbeParameters param)
    {
        int n = targetSequence.Length;

        // gcPrefixSum[i] / validPrefixSum[i] = count of G/C / of valid nucleotides (A/C/G/T/U) in sequence[0..i-1],
        // classified by the canonical SequenceExtensions.CountGcAndValidNucleotides, so a window's GC fraction equals
        // CalculateGcFractionFast(window) (= Primer3 gc_and_n_content: G+C over the non-N bases) in O(1).
        int[] gcPrefixSum = new int[n + 1];
        int[] validPrefixSum = new int[n + 1];
        for (int i = 0; i < n; i++)
        {
            var (gcBase, validBase) = targetSequence.AsSpan(i, 1).CountGcAndValidNucleotides();
            gcPrefixSum[i + 1] = gcPrefixSum[i] + gcBase;
            validPrefixSum[i + 1] = validPrefixSum[i] + validBase;
        }

        var bases = new List<ProbeBase>();
        int order = 0;
        for (int length = param.MinLength; length <= param.MaxLength && length <= n; length++)
        {
            for (int start = 0; start <= n - length; start++)
            {
                int gcCount = gcPrefixSum[start + length] - gcPrefixSum[start];
                int validCount = validPrefixSum[start + length] - validPrefixSum[start];
                double gc = validCount == 0 ? 0 : (double)gcCount / validCount;

                // Early rejection far outside the GC window.
                if (gc < param.MinGc - 0.1 || gc > param.MaxGc + 0.1)
                    continue;

                var b = EvaluateProbeBase(targetSequence.Substring(start, length), start, order++, param, gc);
                if (b.BaseScore > 0)
                    bases.Add(b);
            }
        }
        return bases;
    }

    /// <summary>
    /// Evaluates a potential probe sequence (all screens, eager).
    /// </summary>
    private static Probe? EvaluateProbe(string sequence, int start, ProbeParameters param)
    {
        double gc = sequence.CalculateGcFractionFast();
        return FinishProbe(EvaluateProbeBase(sequence, start, 0, param, gc), param);
    }

    /// <summary>
    /// Designs tiling probes to cover entire sequence.
    /// </summary>
    /// <remarks>
    /// Windows of <paramref name="probeLength"/> start every <c>probeLength − overlap</c> bases; each is scored
    /// like <see cref="DesignProbes(string, ProbeParameters?, int)"/> (Tm = Primer3 <c>seqtm</c> at the
    /// parameters' conditions, thermodynamic self-structure screen for ≤ ThermodynamicScreenMaxLength nt, default 60). A window whose score is
    /// ≤ 0 is still emitted (score 0.3, warning "Suboptimal probe…") so coverage is preserved. Windows keep their
    /// tiling order; with <see cref="ProbeParameters.Ranking"/> = <see cref="ProbeRanking.Primer3Penalty"/> each scored
    /// window also carries <see cref="Probe.Primer3Penalty"/>.
    /// </remarks>
    public static TilingProbeSet DesignTilingProbes(
        string targetSequence,
        int probeLength = 60,
        int overlap = 20,
        ProbeParameters? parameters = null)
    {
        var param = parameters ?? Defaults.Microarray with
        {
            MinLength = probeLength,
            MaxLength = probeLength
        };
        ValidateRanking(param, nameof(parameters));
        ValidateThermodynamicScreenMaxLength(param, nameof(parameters));

        targetSequence = targetSequence.ToUpperInvariant();
        var probes = new List<Probe>();
        int step = probeLength - overlap;

        for (int start = 0; start <= targetSequence.Length - probeLength; start += step)
        {
            string probeSeq = targetSequence.Substring(start, probeLength);
            var probe = EvaluateProbe(probeSeq, start, param);

            if (probe.HasValue)
            {
                probes.Add(probe.Value with { Type = ProbeType.Tiling });
            }
            else
            {
                // Add with warnings for coverage
                double tm = CalculateProbeTm(probeSeq, param);
                double gc = probeSeq.CalculateGcFractionFast();
                probes.Add(new Probe(
                    probeSeq, start, start + probeLength - 1,
                    double.IsNaN(tm) ? 0.0 : tm, gc, 0.3, ProbeType.Tiling,
                    new List<string> { "Suboptimal probe, included for coverage" }));
            }
        }

        // Calculate coverage
        int covered = 0;
        var coveredPositions = new bool[targetSequence.Length];
        foreach (var probe in probes)
        {
            for (int i = probe.Start; i <= probe.End && i < targetSequence.Length; i++)
            {
                if (!coveredPositions[i])
                {
                    coveredPositions[i] = true;
                    covered++;
                }
            }
        }

        double meanTm = probes.Average(p => p.Tm);
        double tmRange = probes.Max(p => p.Tm) - probes.Min(p => p.Tm);

        return new TilingProbeSet(probes, covered, meanTm, tmRange);
    }

    /// <summary>
    /// Designs antisense probe for RNA detection.
    /// </summary>
    public static IEnumerable<Probe> DesignAntisenseProbes(
        string mRnaSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 5)
    {
        // Get reverse complement for antisense probes
        string antisense = DnaSequence.GetReverseComplementString(mRnaSequence);

        foreach (var probe in DesignProbes(antisense, parameters, maxProbes))
        {
            yield return probe with { Type = ProbeType.Antisense };
        }
    }

    /// <summary>
    /// Settings of the Primer3 hybridization-probe (internal-oligo) picker
    /// (<see cref="DesignProbesPrimer3"/>). Defaults are Primer3's <c>PRIMER_INTERNAL_*</c> defaults
    /// (<c>libprimer3.cc</c> <c>pr_set_default_global_args</c>, <c>o_args</c>).
    /// </summary>
    /// <param name="MinSize">PRIMER_INTERNAL_MIN_SIZE (18).</param>
    /// <param name="OptSize">PRIMER_INTERNAL_OPT_SIZE (20).</param>
    /// <param name="MaxSize">PRIMER_INTERNAL_MAX_SIZE (27; at most 36, Primer3 MAX_PRIMER_LENGTH).</param>
    /// <param name="MinTm">PRIMER_INTERNAL_MIN_TM (57 °C).</param>
    /// <param name="OptTm">PRIMER_INTERNAL_OPT_TM (60 °C).</param>
    /// <param name="MaxTm">PRIMER_INTERNAL_MAX_TM (63 °C).</param>
    /// <param name="MinGcPercent">PRIMER_INTERNAL_MIN_GC (20 %).</param>
    /// <param name="MaxGcPercent">PRIMER_INTERNAL_MAX_GC (80 %).</param>
    /// <param name="MaxPolyX">PRIMER_INTERNAL_MAX_POLY_X (5).</param>
    /// <param name="MaxSelfAnyTh">PRIMER_INTERNAL_MAX_SELF_ANY_TH (47 °C).</param>
    /// <param name="MaxSelfEndTh">PRIMER_INTERNAL_MAX_SELF_END_TH (47 °C).</param>
    /// <param name="MaxHairpinTh">PRIMER_INTERNAL_MAX_HAIRPIN_TH (47 °C).</param>
    /// <param name="MonovalentMillimolar">PRIMER_INTERNAL_SALT_MONOVALENT (50 mM).</param>
    /// <param name="DivalentMillimolar">PRIMER_INTERNAL_SALT_DIVALENT (0 mM).</param>
    /// <param name="DntpMillimolar">PRIMER_INTERNAL_DNTP_CONC (0 mM).</param>
    /// <param name="DnaConcentrationNanomolar">PRIMER_INTERNAL_DNA_CONC (50 nM).</param>
    public sealed record Primer3ProbeSettings(
        int MinSize = 18,
        int OptSize = 20,
        int MaxSize = 27,
        double MinTm = 57.0,
        double OptTm = 60.0,
        double MaxTm = 63.0,
        double MinGcPercent = 20.0,
        double MaxGcPercent = 80.0,
        int MaxPolyX = 5,
        double MaxSelfAnyTh = PrimerDesigner.Primer3MaxStructureTm,
        double MaxSelfEndTh = PrimerDesigner.Primer3MaxStructureTm,
        double MaxHairpinTh = PrimerDesigner.Primer3MaxStructureTm,
        double MonovalentMillimolar = PrimerDesigner.Primer3InternalMonovalentMillimolar,
        double DivalentMillimolar = PrimerDesigner.Primer3InternalDivalentMillimolar,
        double DntpMillimolar = PrimerDesigner.Primer3InternalDntpMillimolar,
        double DnaConcentrationNanomolar = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar)
    {
        /// <summary>
        /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: true (Primer3 default) = ntthal self-any / self-end / hairpin Tm limits;
        /// false = Primer3 alignment mode, dpal <c>self_any</c> / <c>self_end</c>
        /// (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/> /
        /// <see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>) limited by <see cref="MaxSelfAny"/> /
        /// <see cref="MaxSelfEnd"/>, no hairpin check. For an internal oligo of a primer pair the mode follows the
        /// primer screen (<see cref="PrimerStructureScreen.Primer3Alignment"/> ⇒ false).
        /// </summary>
        public bool ThermodynamicOligoAlignment { get; init; } = true;

        /// <summary>PRIMER_INTERNAL_MAX_SELF_ANY (alignment mode; Primer3 default 12.00, must be in [0, 32767]).</summary>
        public double MaxSelfAny { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>PRIMER_INTERNAL_MAX_SELF_END (alignment mode; Primer3 default 12.00, must be in [0, 32767]).</summary>
        public double MaxSelfEnd { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>
        /// PRIMER_INTERNAL_OPT_GC_PERCENT: the GC optimum of the internal-oligo GC penalty terms
        /// (<see cref="WeightGcPercentGt"/>/<see cref="WeightGcPercentLt"/>, <c>p_obj_fn</c> OT_INTL); null = undefined, as in
        /// Primer3's code (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>). Inert while both GC weights are 0
        /// (Primer3's default); a non-zero GC weight without it is rejected with Primer3's <c>_pr_data_control</c> error ("Hyb
        /// probe GC content is part of objective function while optimum gc_content is not defined", <see cref="ArgumentException"/>).
        /// </summary>
        public double? OptGcPercent { get; init; }

        /// <summary>PRIMER_INTERNAL_WT_GC_PERCENT_GT (Primer3 default 0): weight × (GC% − <see cref="OptGcPercent"/>) when GC% is above the optimum.</summary>
        public double WeightGcPercentGt { get; init; }

        /// <summary>PRIMER_INTERNAL_WT_GC_PERCENT_LT (Primer3 default 0): weight × (<see cref="OptGcPercent"/> − GC%) when GC% is below the optimum.</summary>
        public double WeightGcPercentLt { get; init; }

        /// <summary>
        /// PRIMER_INTERNAL_MISHYB_LIBRARY (primer3-py <c>mishyb_lib</c>): when set and non-empty, every oligo gets Primer3's
        /// library mishybridization score (<see cref="PrimerDesigner.CalculateLibraryMishyb"/>,
        /// <see cref="Primer3Probe.LibraryMishyb"/>) and is rejected when any entry's weighted score exceeds
        /// <see cref="MaxLibraryMishyb"/> (Primer3 OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ, a "five-prime problem" that also ends
        /// the 5′ extension of that 3′ end). Null = no library (Primer3's default). Not part of the JSON form (the MCP tools
        /// take the library as a separate name → sequence argument).
        /// </summary>
        [System.Text.Json.Serialization.JsonIgnore]
        public PrimerMisprimingLibrary? MishybLibrary { get; init; }

        /// <summary>
        /// PRIMER_INTERNAL_MAX_LIBRARY_MISHYB (Primer3 default <see cref="PrimerDesigner.Primer3InternalMaxLibraryMishyb"/>,
        /// 12.00): maximum weighted library score of one oligo; Primer3 compares with the value truncated to a C
        /// <c>short</c>. In alignment mode (<see cref="ThermodynamicOligoAlignment"/> = false) it must not exceed 32767.
        /// </summary>
        public double MaxLibraryMishyb { get; init; } = PrimerDesigner.Primer3InternalMaxLibraryMishyb;

        /// <summary>
        /// PRIMER_INTERNAL_WT_LIBRARY_MISHYB (Primer3 default 0): penalty weight × <see cref="Primer3Probe.LibraryMishyb"/>
        /// (<c>p_obj_fn</c> OT_INTL <c>repeat_sim</c> term). Non-zero requires <see cref="MishybLibrary"/>.
        /// </summary>
        public double WeightLibraryMishyb { get; init; }

        /// <summary>
        /// PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS for the mishyb library (Primer3 default 0 = false: an IUPAC code in an
        /// entry never aligns). Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.LibraryAmbiguityCodesConsensus"/> is used instead.
        /// </summary>
        public bool LibraryAmbiguityCodesConsensus { get; init; }

        // The mishyb library in use: null when absent or empty (Primer3 seq_lib_num_seq == 0).
        internal PrimerMisprimingLibrary? ActiveMishybLibrary => MishybLibrary is { Count: > 0 } lib ? lib : null;

        /// <summary>
        /// PRIMER_ANNEALING_TEMP (°C; Primer3 default <see cref="PrimerDesigner.Primer3DefaultAnnealingTemperature"/> = −10
        /// = off; at most 100). When &gt; 0 every oligo gets its fraction bound at this temperature
        /// (<see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/> at the settings' conditions,
        /// <see cref="Primer3Probe.Bound"/>) and is rejected outside [<see cref="MinBound"/>, <see cref="MaxBound"/>].
        /// Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.AnnealingTemperature"/> is used instead.
        /// </summary>
        public double AnnealingTemperature { get; init; } = PrimerDesigner.Primer3DefaultAnnealingTemperature;

        /// <summary>PRIMER_INTERNAL_MIN_BOUND (% bound; default −10), checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
        public double MinBound { get; init; } = PrimerDesigner.Primer3MinBound;

        /// <summary>PRIMER_INTERNAL_MAX_BOUND (% bound; default 110), checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
        public double MaxBound { get; init; } = PrimerDesigner.Primer3MaxBound;

        /// <summary>PRIMER_INTERNAL_OPT_BOUND (% bound; default 97; must lie in [<see cref="MinBound"/>, <see cref="MaxBound"/>]).</summary>
        public double OptBound { get; init; } = PrimerDesigner.Primer3OptBound;

        /// <summary>
        /// PRIMER_INTERNAL_WT_BOUND_GT (default 0): weight × (bound − <see cref="OptBound"/>) above the optimum. Unlike the
        /// primer terms, Primer3's internal-oligo <c>p_obj_fn</c> applies the bound terms whatever
        /// <see cref="AnnealingTemperature"/> is; without an annealing temperature (or for an oligo with no bound value)
        /// the bound is Primer3's OLIGOTM_ERROR = −999999.9999, so only <see cref="WeightBoundLt"/> then contributes —
        /// weight × (OptBound + 999999.9999), reproduced as Primer3 computes it.
        /// </summary>
        public double WeightBoundGt { get; init; }

        /// <summary>PRIMER_INTERNAL_WT_BOUND_LT (default 0): weight × (<see cref="OptBound"/> − bound) below the optimum
        /// (see <see cref="WeightBoundGt"/> for the no-annealing-temperature case).</summary>
        public double WeightBoundLt { get; init; }

        /// <summary>
        /// PRIMER_INTERNAL_MIN_QUALITY (default 0): with SEQUENCE_QUALITY an oligo whose minimum base quality
        /// (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>, <see cref="Primer3Probe.MinSequenceQuality"/>) is lower is
        /// rejected (Primer3 <c>sequence_quality_is_ok</c>, a "five-prime problem"). A non-zero value requires quality data and must
        /// lie in [<see cref="QualityRangeMin"/>, <see cref="QualityRangeMax"/>]. Internal oligos have no end-quality limit.
        /// </summary>
        public int MinQuality { get; init; }

        /// <summary>
        /// PRIMER_LOWERCASE_MASKING (default false): an oligo whose 3′ base (its rightmost template base) is a lower-case
        /// a/c/g/t of the template as given is rejected (Primer3 <c>is_lowercase_masked</c> on <c>trimmed_orig_seq</c>);
        /// lower-case bases elsewhere are accepted. Used by <see cref="DesignProbesPrimer3"/>; for PRIMER_PICK_INTERNAL_OLIGO the
        /// pair search's <see cref="PrimerPairOptions.LowercaseMasking"/> applies.
        /// </summary>
        public bool LowercaseMasking { get; init; }

        /// <summary>PRIMER_INTERNAL_WT_SEQ_QUAL (default 0): × (<see cref="QualityRangeMax"/> − the oligo's minimum base quality),
        /// the last term of Primer3's internal-oligo <c>p_obj_fn</c>. Non-zero requires SEQUENCE_QUALITY.</summary>
        public double WeightSequenceQuality { get; init; }

        /// <summary>PRIMER_INTERNAL_WT_END_QUAL (default 0): parsed by Primer3 2.3.1 but never read by <c>p_obj_fn</c> — accepted,
        /// <b>no effect</b> (reproduced).</summary>
        public double WeightEndQuality { get; init; }

        /// <summary>PRIMER_QUALITY_RANGE_MIN (default 0). Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.QualityRangeMin"/> is used instead.</summary>
        public int QualityRangeMin { get; init; } = PrimerDesigner.Primer3QualityRangeMin;

        /// <summary>PRIMER_QUALITY_RANGE_MAX (default 100). Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.QualityRangeMax"/> is used instead.</summary>
        public int QualityRangeMax { get; init; } = PrimerDesigner.Primer3QualityRangeMax;
    }

    /// <summary>
    /// A hybridization probe picked by <see cref="DesignProbesPrimer3"/> — the Primer3
    /// <c>PRIMER_INTERNAL_n_*</c> output values.
    /// </summary>
    /// <param name="Sequence">Probe sequence (5′→3′, upper case; the template strand).</param>
    /// <param name="Start">0-based start in the template (PRIMER_INTERNAL_n = start,length).</param>
    /// <param name="Length">Probe length.</param>
    /// <param name="Tm">PRIMER_INTERNAL_n_TM (Primer3 <c>seqtm</c>, °C).</param>
    /// <param name="GcPercent">PRIMER_INTERNAL_n_GC_PERCENT.</param>
    /// <param name="SelfAnyTh">PRIMER_INTERNAL_n_SELF_ANY_TH (ntthal ANY self-dimer Tm, °C).</param>
    /// <param name="SelfEndTh">PRIMER_INTERNAL_n_SELF_END_TH (ntthal END1 self-dimer Tm, °C).</param>
    /// <param name="HairpinTh">PRIMER_INTERNAL_n_HAIRPIN_TH (ntthal hairpin Tm, °C).</param>
    /// <param name="Penalty">PRIMER_INTERNAL_n_PENALTY (Primer3 <c>p_obj_fn</c>, internal-oligo branch).</param>
    public readonly record struct Primer3Probe(
        string Sequence,
        int Start,
        int Length,
        double Tm,
        double GcPercent,
        // NaN (alignment mode: no ntthal value) is written as the JSON literal "NaN" instead of failing serialization.
        [property: System.Text.Json.Serialization.JsonNumberHandling(System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals)] double SelfAnyTh,
        [property: System.Text.Json.Serialization.JsonNumberHandling(System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals)] double SelfEndTh,
        [property: System.Text.Json.Serialization.JsonNumberHandling(System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals)] double HairpinTh,
        double Penalty)
    {
        /// <summary>PRIMER_INTERNAL_n_SELF_ANY (Primer3 alignment mode, dpal score); null in thermodynamic mode, where
        /// <see cref="SelfAnyTh"/>/<see cref="SelfEndTh"/>/<see cref="HairpinTh"/> are set (they are NaN in alignment mode).</summary>
        public double? SelfAny { get; init; }

        /// <summary>PRIMER_INTERNAL_n_SELF_END (Primer3 alignment mode, dpal score); null in thermodynamic mode.</summary>
        public double? SelfEnd { get; init; }

        /// <summary>
        /// PRIMER_INTERNAL_n_LIBRARY_MISHYB score (primer3-py key PRIMER_INTERNAL_n_LIBRARY_MISPRIMING;
        /// <see cref="PrimerDesigner.CalculateLibraryMishyb"/>) when <see cref="Primer3ProbeSettings.MishybLibrary"/> is
        /// set, otherwise null.
        /// </summary>
        public double? LibraryMishyb { get; init; }

        /// <summary>The library entry named in PRIMER_INTERNAL_n_LIBRARY_MISHYB (with <see cref="LibraryMishyb"/>).</summary>
        public string? LibraryMishybName { get; init; }

        /// <summary>PRIMER_INTERNAL_n_BOUND: the fraction (%) bound at <see cref="Primer3ProbeSettings.AnnealingTemperature"/>
        /// when that is &gt; 0 (<see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/>), otherwise <c>null</c>.</summary>
        public double? Bound { get; init; }

        /// <summary>PRIMER_INTERNAL_n_MIN_SEQ_QUALITY: the oligo's minimum base quality when SEQUENCE_QUALITY is given
        /// (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>), otherwise <c>null</c>.</summary>
        public int? MinSequenceQuality { get; init; }
    }

    // Primer3 MAX_PRIMER_LENGTH (oligo length limit of the picker and of seqtm's nearest-neighbour branch).
    private const int Primer3MaxOligoLength = 36;

    /// <summary>
    /// Picks hybridization probes exactly as Primer3 does for <c>PRIMER_TASK=pick_hyb_probe_only</c>
    /// (Rozen &amp; Skaletsky 2000; Untergasser et al. 2012; <c>libprimer3.cc</c> <c>make_internal_oligo_list</c>
    /// → <c>pick_primer_range</c> → <c>calc_and_check_oligo_features</c> → <c>p_obj_fn</c>, default
    /// thermodynamic mode): every window of <see cref="Primer3ProbeSettings.MinSize"/>..<see cref="Primer3ProbeSettings.MaxSize"/>
    /// bases made only of A/C/G/T (PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0) is accepted when its G+C % is within
    /// [MinGcPercent, MaxGcPercent], its longest mononucleotide run ≤ MaxPolyX, its Tm (Primer3 <c>seqtm</c> with
    /// MAX_NN_TM_LENGTH = 36,
    /// <see cref="PrimerDesigner.CalculateMeltingTemperaturePrimer3"/>, at the settings' conditions) within
    /// [MinTm, MaxTm], and its ntthal self-dimer, 3′ self-dimer and hairpin Tm
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>) ≤ MaxSelfAnyTh / MaxSelfEndTh / MaxHairpinTh.
    /// As in Primer3, windows are enumerated per 3′ end with growing length, and an N, a poly-X run or a too-stable
    /// self-dimer (Primer3's "five-prime problems") stops the extension for that 3′ end, so a longer window that
    /// would pass is not considered. Accepted probes are ranked by the Primer3 penalty |Tm − OptTm| + |length − OptSize|
    /// (<see cref="PrimerDesigner.CalculatePrimer3Penalty"/> with Primer3's default internal-oligo weights) and
    /// ordered as Primer3's <c>primer_rec_comp</c>: penalty ascending, then start descending, then length
    /// ascending. With <see cref="Primer3ProbeSettings.ThermodynamicOligoAlignment"/> = false (Primer3
    /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0) the structure limits are Primer3's alignment-mode ones instead: dpal
    /// <c>self_any</c> &gt; <see cref="Primer3ProbeSettings.MaxSelfAny"/> (a five-prime problem) or <c>self_end</c> &gt;
    /// <see cref="Primer3ProbeSettings.MaxSelfEnd"/> (both 12.00 by default) reject the window, no hairpin check, and the
    /// values are reported in <see cref="Primer3Probe.SelfAny"/>/<see cref="Primer3Probe.SelfEnd"/>.
    /// Verified against primer3-py 2.3.1 <c>design_primers</c> (PRIMER_INTERNAL_n_* values, both modes).
    /// </summary>
    /// <param name="template">Template sequence (case-insensitive); probes are picked on this strand.</param>
    /// <param name="settings">Picker settings (default: Primer3 defaults).</param>
    /// <param name="numReturn">PRIMER_NUM_RETURN (default 5; at least 1, as Primer3's <c>_pr_data_control</c> requires).</param>
    /// <param name="sequenceQuality">SEQUENCE_QUALITY: one integer quality per template base (null or empty = none). With it
    /// every oligo gets <see cref="Primer3Probe.MinSequenceQuality"/>, is checked against
    /// <see cref="Primer3ProbeSettings.MinQuality"/> and weighted by <see cref="Primer3ProbeSettings.WeightSequenceQuality"/>;
    /// its length must equal the template length and its values lie in [<see cref="Primer3ProbeSettings.QualityRangeMin"/>,
    /// <see cref="Primer3ProbeSettings.QualityRangeMax"/>] (Primer3 <c>_pr_data_control</c>).</param>
    /// <returns>Up to <paramref name="numReturn"/> probes, best first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid sizes (MinSize &lt; 1, MaxSize &lt; MinSize,
    /// MaxSize &gt; 36), an invalid mishyb limit, or <paramref name="numReturn"/> &lt; 1 (Primer3 "PRIMER_NUM_RETURN &lt; 1").</exception>
    /// <exception cref="ArgumentException">A GC weight without <see cref="Primer3ProbeSettings.OptGcPercent"/> or a mishyb weight
    /// without a mishyb library (Primer3 <c>_pr_data_control</c>).</exception>
    public static IReadOnlyList<Primer3Probe> DesignProbesPrimer3(
        string template,
        Primer3ProbeSettings? settings = null,
        int numReturn = 5,
        IReadOnlyList<int>? sequenceQuality = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var s = settings ?? new Primer3ProbeSettings();
        ValidatePrimer3ProbeSettings(s, nameof(settings));
        if (numReturn < 1)
            throw new ArgumentOutOfRangeException(nameof(numReturn), "PRIMER_NUM_RETURN < 1 (Primer3 _pr_data_control).");
        PrimerDesigner.ValidatePrimer3Quality(sequenceQuality, template.Length, 0, s.MinQuality, s.QualityRangeMin,
            s.QualityRangeMax, 0.0, s.WeightSequenceQuality, nameof(sequenceQuality));

        var accepted = EnumeratePrimer3InternalOligos(template.ToUpperInvariant(), 0, template.Length, s, screenStructure: true,
            sequenceQuality is { Count: > 0 } ? sequenceQuality : null, s.LowercaseMasking ? template : null);

        // primer_rec_comp: quality ascending, then start descending, then length ascending.
        accepted.Sort((a, b) => ComparePrimer3Rank(a.Penalty, a.Start, a.Length, b.Penalty, b.Start, b.Length));

        return accepted.Count > numReturn ? accepted.GetRange(0, numReturn) : accepted;
    }

    // _pr_data_control: PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT need PRIMER_INTERNAL_OPT_GC_PERCENT.
    internal static void ValidatePrimer3ProbeGcOptimum(Primer3ProbeSettings s, string paramName)
    {
        if ((s.WeightGcPercentGt != 0 || s.WeightGcPercentLt != 0) && s.OptGcPercent is null)
            throw new ArgumentException(PrimerDesigner.ProbeGcOptimumUndefinedMessage, paramName);
    }

    // Primer3 settings check shared by the hybridization-probe picker and PrimerDesigner's
    // PRIMER_PICK_INTERNAL_OLIGO path.
    internal static void ValidatePrimer3ProbeSettings(Primer3ProbeSettings s, string paramName)
    {
        if (s.MinSize < 1 || s.MaxSize < s.MinSize || s.MaxSize > Primer3MaxOligoLength)
            throw new ArgumentOutOfRangeException(paramName,
                $"Probe sizes must satisfy 1 ≤ MinSize ≤ MaxSize ≤ {Primer3MaxOligoLength}.");
        if (!(s.MaxSelfAny >= 0 && s.MaxSelfAny <= short.MaxValue && s.MaxSelfEnd >= 0 && s.MaxSelfEnd <= short.MaxValue))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for internal oligo complementarity restrictions (Primer3: 0 ≤ PRIMER_INTERNAL_MAX_SELF_ANY/_END ≤ 32767).");
        PrimerDesigner.ValidatePrimer3Conditions(
            s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar, s.DnaConcentrationNanomolar, paramName);
        if (s.OptGcPercent is { } opt && !double.IsFinite(opt))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_INTERNAL_OPT_GC_PERCENT must be finite.");
        ValidatePrimer3ProbeGcOptimum(s, paramName);
        // _pr_data_control: PRIMER_INTERNAL_MAX_LIBRARY_MISHYB > SHRT_MAX (alignment mode); a mishyb weight without a library.
        if (double.IsNaN(s.MaxLibraryMishyb) || Math.Abs(s.MaxLibraryMishyb) >= int.MaxValue
            || (s.MaxLibraryMishyb > short.MaxValue && !s.ThermodynamicOligoAlignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_INTERNAL_MAX_LIBRARY_MISHYB.");
        if (s.WeightLibraryMishyb != 0 && s.ActiveMishybLibrary is null)
            throw new ArgumentException(
                "Internal oligo mispriming score is part of objective function while mishyb library is not defined (Primer3 _pr_data_control).",
                paramName);
        // _pr_data_control: PRIMER_INTERNAL_OPT_BOUND within [MIN, MAX]; PRIMER_ANNEALING_TEMP ≤ 100.
        if (double.IsNaN(s.WeightBoundGt) || double.IsNaN(s.WeightBoundLt))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_INTERNAL_WT_BOUND_GT/_LT must not be NaN.");
        PrimerDesigner.ValidatePrimer3Bound(s.AnnealingTemperature, s.MinBound, s.MaxBound, s.OptBound,
            internalOligo: true, paramName);
    }

    /// <summary>
    /// Primer3 internal-oligo candidate list (<c>make_internal_oligo_list</c> → <c>pick_primer_range</c>
    /// → <c>calc_and_check_oligo_features</c>, <c>OT_INTL</c>, thermodynamic mode) over
    /// <c>seq[regionStart, regionEnd)</c>, in Primer3's enumeration order (3′ end descending, then length
    /// ascending), with <see cref="Primer3Probe.Penalty"/> = <c>p_obj_fn</c> (internal-oligo weights).
    /// With <paramref name="screenStructure"/> = true (PRIMER_TASK=pick_hyb_probe_only, a <c>primer_list</c>
    /// output) the ntthal self-any / self-end / hairpin limits are applied while enumerating, exactly as
    /// Primer3 does for list output (a too-stable self-dimer is a "five-prime problem" ending the extension);
    /// with false (internal oligo for a primer pair, <c>primer_pairs</c> output) Primer3 postpones them to
    /// <c>choose_internal_oligo</c>: the structure fields are then <see cref="double.NaN"/> and the caller
    /// screens the oligo it picks (<see cref="PassesPrimer3ProbeStructure"/>).
    /// </summary>
    internal static List<Primer3Probe> EnumeratePrimer3InternalOligos(
        string seq, int regionStart, int regionEnd, Primer3ProbeSettings s, bool screenStructure,
        IReadOnlyList<int>? quality = null, string? lowercase = null)
    {
        // Primer3 o_args.weights defaults, with PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT and PRIMER_INTERNAL_OPT_GC_PERCENT.
        // PRIMER_INTERNAL_WT_LIBRARY_MISHYB is the o_args repeat_sim weight.
        var weights = PrimerDesigner.DefaultPrimer3Weights with
        {
            GcGt = s.WeightGcPercentGt,
            GcLt = s.WeightGcPercentLt,
            LibraryMispriming = s.WeightLibraryMishyb,
            BoundGt = s.WeightBoundGt,
            BoundLt = s.WeightBoundLt,
            SequenceQuality = s.WeightSequenceQuality,
        };
        // PRIMER_INTERNAL_MISHYB_LIBRARY: scored while enumerating for list output (three_conditions) or when weighted
        // (calc_and_check_oligo_features), otherwise postponed to choose_internal_oligo.
        var library = s.ActiveMishybLibrary;
        bool scoreLibrary = library is not null && (screenStructure || s.WeightLibraryMishyb != 0);
        var optima = new Primer3Optima(s.OptTm, s.OptSize, s.OptGcPercent)
        {
            OptBound = s.OptBound,
        };
        double monovalentEq = PrimerDesigner.Primer3MonovalentEquivalent(s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar);
        var accepted = new List<Primer3Probe>();

        // pick_primer_range: for every 3' end, oligos of increasing length (5' extensions). A failure that no
        // 5' extension can cure — too many Ns, poly-X, self-any (Primer3 five_prime_problem bits) — ends the
        // extension loop for that 3' end; checks run in calc_and_check_oligo_features order and stop at the
        // first failure (Ns, GC, poly-X, Tm, self-any, self-end, hairpin).
        for (int end = regionEnd - 1; end >= regionStart + s.MinSize - 1; end--)
        {
            // PRIMER_LOWERCASE_MASKING (is_lowercase_masked, before every other check): a lower-case 3' base rejects every
            // oligo ending here (lowercase = the case-preserving template, null when masking is off).
            if (lowercase is not null && PrimerDesigner.IsLowercaseMaskedBase(lowercase[end]))
                continue;
            for (int len = s.MinSize; len <= s.MaxSize; len++)
            {
                int start = end - len + 1;
                if (start < regionStart)
                    break;

                string oligo = seq.Substring(start, len);
                int gc = 0, ns = 0;
                bool other = false;
                foreach (char c in oligo)
                {
                    if (c is 'G' or 'C') gc++;
                    else if (c == 'N') ns++;
                    else if (c is not ('A' or 'T')) other = true;
                }
                if (ns > 0)
                    break; // OP_TOO_MANY_NS (PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0): five-prime problem

                double gcPercent = 100.0 * gc / len;
                if (gcPercent < s.MinGcPercent || gcPercent > s.MaxGcPercent)
                    continue;
                // sequence_quality_is_ok (OT_INTL: the minimum quality only), after the GC check: OP_LOW_SEQUENCE_QUALITY
                // is a five-prime problem (every 5' extension keeps the low base).
                int? minQuality = null;
                if (quality is not null)
                {
                    minQuality = PrimerDesigner.CalculateSequenceQualityPrimer3(quality, start, len, true, s.QualityRangeMax).Min;
                    if (minQuality < s.MinQuality)
                        break;
                }
                if (PrimerDesigner.FindLongestHomopolymer(oligo) > s.MaxPolyX)
                    break; // OP_HIGH_POLY_X: five-prime problem
                if (other)
                    continue; // non-ACGT, non-N symbol: no Tm (seqtm error → low Tm)

                var (tm, bound) = PrimerDesigner.Primer3SeqTm(oligo, s.DnaConcentrationNanomolar, monovalentEq, s.AnnealingTemperature);
                if (tm < s.MinTm || tm > s.MaxTm)
                    continue;
                // Fraction bound (PRIMER_ANNEALING_TEMP > 0; not a five-prime problem).
                if (s.AnnealingTemperature > 0.0 && (bound < s.MinBound || bound > s.MaxBound))
                    continue;

                double selfAny = double.NaN, selfEnd = double.NaN, hairpin = double.NaN;
                double? alnAny = null, alnEnd = null;
                if (screenStructure && !s.ThermodynamicOligoAlignment)
                {
                    // oligo_compl (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0): self_any, then self_end; no hairpin.
                    alnAny = PrimerDesigner.CalculatePrimerSelfAnyComplementarity(oligo);
                    if (alnAny > s.MaxSelfAny)
                        break; // OP_HIGH_SELF_ANY: five-prime problem
                    alnEnd = PrimerDesigner.CalculatePrimerSelfEndComplementarity(oligo);
                    if (alnEnd > s.MaxSelfEnd)
                        continue;
                }
                else if (screenStructure)
                {
                    var st = ComputePrimer3ProbeStructure(oligo, s);
                    if (st.SelfAnyTh > s.MaxSelfAnyTh)
                        break; // OP_HIGH_SELF_ANY: five-prime problem
                    if (st.SelfEndTh > s.MaxSelfEndTh || st.HairpinTh > s.MaxHairpinTh)
                        continue;
                    (selfAny, selfEnd, hairpin) = (st.SelfAnyTh, st.SelfEndTh, st.HairpinTh);
                }

                // oligo_repeat_library_mispriming (OT_INTL), after the structure checks.
                PrimerDesigner.LibraryMispriming? lib = null;
                if (scoreLibrary)
                {
                    lib = ComputePrimer3ProbeLibrary(oligo, s, library!);
                    if (lib.Exceeds)
                        break; // OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ: five-prime problem
                }

                // p_obj_fn OT_INTL: no end_stability term (end_oligodg is computed for primers only).
                double penalty = PrimerDesigner.CalculatePrimer3Penalty(
                    new Primer3PenaltyInputs(tm, len, gcPercent)
                    {
                        LibraryMispriming = lib?.MaxScore ?? 0.0,
                        Bound = bound,
                        SequenceQuality = minQuality,
                        QualityRangeMax = s.QualityRangeMax,
                    }, weights, optima);
                accepted.Add(new Primer3Probe(oligo, start, len, tm, gcPercent, selfAny, selfEnd, hairpin, penalty)
                {
                    SelfAny = alnAny,
                    SelfEnd = alnEnd,
                    LibraryMishyb = lib?.MaxScore,
                    LibraryMishybName = lib?.Name,
                    Bound = s.AnnealingTemperature > 0.0 && bound != PrimerDesigner.Primer3OligoTmError ? bound : null,
                    MinSequenceQuality = minQuality,
                });
            }
        }
        return accepted;
    }

    private static PrimerDesigner.LibraryMispriming ComputePrimer3ProbeLibrary(
        string oligo, Primer3ProbeSettings s, PrimerMisprimingLibrary library) =>
        PrimerDesigner.ComputeLibraryMispriming(oligo, true, library, s.LibraryAmbiguityCodesConsensus,
            s.MaxLibraryMishyb, isInternal: true);

    private static PrimerDesigner.Primer3OligoStructure ComputePrimer3ProbeStructure(string oligo, Primer3ProbeSettings s) =>
        PrimerDesigner.CalculatePrimer3OligoStructure(
            oligo, s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar, s.DnaConcentrationNanomolar)!.Value;

    /// <summary>
    /// Primer3 <c>choose_internal_oligo</c> postponed checks of one internal oligo (<c>oligo_compl_thermod</c>:
    /// self-any and self-end Tm, then <c>oligo_hairpin</c>; alignment mode: <c>oligo_compl</c> dpal self_any /
    /// self_end; then, unless already scored while enumerating, <c>oligo_repeat_library_mispriming</c> against the
    /// mishyb library) at the settings' conditions; returns the oligo with its structure (and library) values filled
    /// in, or <c>null</c> when a limit is exceeded.
    /// </summary>
    internal static Primer3Probe? PassesPrimer3ProbeStructure(Primer3Probe probe, Primer3ProbeSettings s)
    {
        if (!s.ThermodynamicOligoAlignment)
        {
            // oligo_compl with o_args (alignment mode): self_any, then self_end.
            double any = PrimerDesigner.CalculatePrimerSelfAnyComplementarity(probe.Sequence);
            if (any > s.MaxSelfAny)
                return null;
            double end = PrimerDesigner.CalculatePrimerSelfEndComplementarity(probe.Sequence);
            if (end > s.MaxSelfEnd)
                return null;
            probe = probe with { SelfAny = any, SelfEnd = end };
        }
        else
        {
            var st = ComputePrimer3ProbeStructure(probe.Sequence, s);
            if (st.SelfAnyTh > s.MaxSelfAnyTh || st.SelfEndTh > s.MaxSelfEndTh || st.HairpinTh > s.MaxHairpinTh)
                return null;
            probe = probe with { SelfAnyTh = st.SelfAnyTh, SelfEndTh = st.SelfEndTh, HairpinTh = st.HairpinTh };
        }
        if (s.ActiveMishybLibrary is { } library && probe.LibraryMishyb is null)
        {
            var lib = ComputePrimer3ProbeLibrary(probe.Sequence, s, library);
            if (lib.Exceeds)
                return null;
            probe = probe with { LibraryMishyb = lib.MaxScore, LibraryMishybName = lib.Name };
        }
        return probe;
    }

    // Tyagi & Kramer (1996) / Marras et al. design rule: the probe–target hybrid Tm and the stem
    // (hairpin) Tm should both be 7–10 °C above the detection (PCR annealing) temperature.
    private const double BeaconMinTmAboveDetection = 7.0;
    private const double BeaconMaxTmAboveDetection = 10.0;

    /// <summary>
    /// Designs a molecular beacon (Tyagi &amp; Kramer 1996): a target-complementary loop of
    /// <paramref name="probeLength"/> bases flanked by complementary GC-rich arms
    /// (5′ arm = ⌊stem/2⌋ G + remaining C, 3′ arm = its reverse complement).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The loop is the window with the best score: −0.2 each for GC outside 40–60 %, loop Tm outside the
    /// Tm window and a homopolymer run &gt; 4. The loop (probe–target) Tm is Primer3 <c>seqtm</c> at the Primer3
    /// hybridization-probe conditions (50 nM, 50 mM monovalent, no Mg²⁺/dNTP; <see cref="ProbeParameters"/>).
    /// With <paramref name="detectionTemperatureCelsius"/> = T the Tm window is the published
    /// [T + 7, T + 10] °C (molecular-beacon design rules of Tyagi &amp; Kramer / Marras et al.: the probe Tm and
    /// the stem Tm 7–10 °C above the detection temperature) and the beacon's stem-loop (hairpin) Tm is
    /// checked against T + 7 °C; without it the library's default window 55–65 °C is used.
    /// </para>
    /// <para>
    /// Warnings carry the stem/loop sizes and, for a beacon of ≤ <c>maxAlignLength</c> (default 60) nt, the ntthal hairpin Tm of the whole
    /// beacon (the stem melting temperature; <see cref="PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double, double, double)"/>
    /// at the same conditions). <see cref="Probe.Tm"/> is the loop (probe–target) Tm.
    /// </para>
    /// </remarks>
    /// <param name="targetSequence">Target sequence (the loop is a window of it).</param>
    /// <param name="probeLength">Loop length (Tyagi &amp; Kramer: 15–30 nt; default 25).</param>
    /// <param name="stemLength">Arm length in bp (Tyagi &amp; Kramer: 5–7 bp; default 5).</param>
    /// <param name="detectionTemperatureCelsius">Optional detection (annealing) temperature T in °C.</param>
    /// <param name="maxAlignLength">ntthal <c>THAL_MAX_ALIGN</c> of the stem-loop hairpin Tm
    /// (<see cref="PrimerDesigner.NtthalMaxAlignLength"/> = 60, Primer3's compile-time default, …
    /// <see cref="PrimerDesigner.NtthalMaxSequenceLength"/> = 10 000): beacons up to this length get the ntthal hairpin
    /// Tm (opt-in for beacons longer than 60 nt; = thal.c compiled with <c>-DTHAL_MAX_ALIGN=…</c>; audit round 3, A3-27).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    public static Probe? DesignMolecularBeacon(
        string targetSequence,
        int probeLength = 25,
        int stemLength = 5,
        double? detectionTemperatureCelsius = null,
        int maxAlignLength = PrimerDesigner.NtthalMaxAlignLength)
    {
        if (maxAlignLength < PrimerDesigner.NtthalMaxAlignLength || maxAlignLength > PrimerDesigner.NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(nameof(maxAlignLength), maxAlignLength,
                "THAL_MAX_ALIGN must be in 60..10000.");
        if (targetSequence.Length < probeLength)
            return null;

        targetSequence = targetSequence.ToUpperInvariant();
        var conditions = Primer3ProbeConditions; // only the (default Primer3 internal-oligo) conditions are used

        double minLoopTm = detectionTemperatureCelsius is { } t0 ? t0 + BeaconMinTmAboveDetection : 55;
        double maxLoopTm = detectionTemperatureCelsius is { } t1 ? t1 + BeaconMaxTmAboveDetection : 65;

        // Find best region in target
        double bestScore = 0;
        string? bestLoop = null;
        int bestStart = 0;
        double bestTm = 0;

        int loopLength = probeLength;
        for (int start = 0; start <= targetSequence.Length - loopLength; start++)
        {
            string loop = targetSequence.Substring(start, loopLength);
            double gc = loop.CalculateGcFractionFast();
            double tmRaw = CalculateProbeTm(loop, conditions);
            double tm = double.IsNaN(tmRaw) ? 0.0 : tmRaw;

            double score = 1.0;
            if (gc < 0.40 || gc > 0.60) score -= 0.2;
            if (double.IsNaN(tmRaw) || tm < minLoopTm || tm > maxLoopTm) score -= 0.2;
            if (PrimerDesigner.FindLongestHomopolymer(loop) > 4) score -= 0.2;

            if (score > bestScore)
            {
                bestScore = score;
                bestLoop = loop;
                bestStart = start;
                bestTm = tm;
            }
        }

        if (bestLoop == null)
            return null;

        // Add stem sequences (GC-rich for stability)
        string stem5 = new string('G', stemLength / 2) + new string('C', stemLength - stemLength / 2);
        string stem3 = DnaSequence.GetReverseComplementString(stem5);

        string beaconSequence = stem5 + bestLoop + stem3;
        var warnings = new List<string> { $"Stem: {stemLength}bp, Loop: {loopLength}bp" };

        if (beaconSequence.Length <= maxAlignLength)
        {
            var hp = PrimerDesigner.CalculateHairpinThermodynamicsNtthal(
                beaconSequence,
                conditions.MonovalentMillimolar / 1000.0,
                conditions.DivalentMillimolar / 1000.0,
                conditions.DntpMillimolar / 1000.0,
                PrimerDesigner.NtthalDefaultTemperatureCelsius,
                PrimerDesigner.NtthalDefaultMaxLoop,
                maxAlignLength);
            if (hp is { } h)
            {
                warnings.Add($"Stem-loop (hairpin) Tm {h.TmCelsius:F1}°C (ntthal)");
                if (detectionTemperatureCelsius is { } t && h.TmCelsius < t + BeaconMinTmAboveDetection)
                    warnings.Add(
                        $"Stem-loop Tm {h.TmCelsius:F1}°C is less than {BeaconMinTmAboveDetection:0} °C above the detection temperature {t:0.#}°C");
            }
        }

        if (detectionTemperatureCelsius is { } td && bestTm < td + BeaconMinTmAboveDetection)
            warnings.Add(
                $"Probe Tm {bestTm:F1}°C is less than {BeaconMinTmAboveDetection:0} °C above the detection temperature {td:0.#}°C");

        return new Probe(
            beaconSequence,
            bestStart,
            bestStart + loopLength - 1,
            bestTm,
            beaconSequence.CalculateGcFractionFast(),
            bestScore,
            ProbeType.MolecularBeacon,
            warnings);
    }

    #endregion

    #region Probe Validation

    /// <summary>
    /// Validates a probe: ungapped k-mismatch hit count over the references, self-structure screen, and
    /// (optionally) the Kane et al. (2000) cross-hybridization criteria against known non-target sequences.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Hits.</b> Every reference is scanned with the canonical ungapped k-mismatch (Hamming) matcher
    /// <see cref="ApproximateMatcher.FindWithMismatches(string, string, int)"/> (case-insensitive, overlapping; the
    /// strand given only by default, and with <paramref name="bothStrands"/> also for the probe's reverse complement —
    /// the probe's sites on the other strand of a double-stranded reference, as blastn <c>-strand both</c> and
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> with <c>bothStrands</c>; a reference
    /// position hit in both orientations — a reverse-palindromic probe, or a reverse-palindromic site — is one site,
    /// counted once, as <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> counts a
    /// palindromic probe once); <see cref="ProbeValidation.OffTargetHits"/> is the total and includes the intended site.
    /// Each hit is judged by the Kane et al. (2000) criteria on its ungapped diagonal — identity (L − mismatches) / L
    /// over the probe length &gt; <paramref name="maxNonTargetIdentity"/> or a run of identical positions longer than
    /// <paramref name="maxContiguousMatch"/> (<see cref="ProbeValidation.CrossHybridizingHits"/>); more than one such
    /// site records an issue (the probe can bind somewhere besides its intended site). With the default 3 mismatches
    /// every hit of a probe ≥ 13 nt meets the identity criterion ((L − 3) / L &gt; 0.75), so the decision equals the
    /// former "more than one hit" rule there. <see cref="ProbeValidation.SpecificityScore"/> = 1/N for N ≥ 1 hits
    /// (0 for none) is a library-defined uniqueness score (the share of the probe's N candidate binding sites
    /// taken by one site), not a published specificity metric (no published 1/N score was found — Kane 2000,
    /// OligoArray 2.0, Li &amp; Stormo 2001 and Primer3's PRIMER_INTERNAL_MAX_LIBRARY_MISHYB judge a site by
    /// identity / free energy / alignment score, not by a hit count); it is reported only and does not enter
    /// <see cref="ProbeValidation.IsValid"/>.
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> The same screen as <see cref="DesignProbes(string, ProbeParameters?, int)"/>: with
    /// <see cref="ProbeStructureScreen.Thermodynamic"/> (default) an A/C/G/T probe of ≤ <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (default 60) nt is screened as Primer3
    /// screens a hybridization probe — ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>, primer3-py <c>calc_homodimer</c> /
    /// <c>calc_end_stability</c> / <c>calc_hairpin</c> parity) at the <paramref name="conditions"/> salt/oligo
    /// concentrations must not exceed <see cref="ProbeParameters.MaxStructureTm"/> (PRIMER_INTERNAL_MAX_SELF_ANY_TH
    /// = _SELF_END_TH = _HAIRPIN_TH = 47 °C). Probes longer than <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/>
    /// (thal.c THAL_MAX_ALIGN, default 60; opt-in larger values run ntthal on longer probes), non-ACGT probes and
    /// <see cref="ProbeStructureScreen.Heuristic"/> use the fallback screens: Primer3 alignment-mode internal-oligo
    /// self-complementarity (dpal <c>self_any</c> / <c>self_end</c> &gt; <see cref="ProbeParameters.MaxSelfAny"/> /
    /// <see cref="ProbeParameters.MaxSelfEnd"/>, PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00, no length limit)
    /// and the sequence-only inverted-repeat hairpin stem. The alignment-mode
    /// values are reported for every probe in <see cref="ProbeValidation.SelfAny"/> / <see cref="ProbeValidation.SelfEnd"/>;
    /// <see cref="ProbeValidation.SelfComplementarity"/> (fold-back fraction) is a library metric only.
    /// </para>
    /// <para>
    /// <b>Cross-hybridization.</b> When <paramref name="nonTargetSequences"/> is given, every non-target (both
    /// strands) is assessed with <see cref="AssessCrossHybridization"/> (Kane et al. 2000: overall identity &gt; 75 %
    /// or a contiguous identical stretch &gt; 15 nt → the probe may cross-hybridize); each cross-hybridizing
    /// non-target strand records an issue.
    /// </para>
    /// <para><see cref="ProbeValidation.IsValid"/> is true when no issue was recorded.</para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence to validate (case-insensitive). Null throws; empty → invalid result.</param>
    /// <param name="referenceSequences">Reference sequences scanned for ungapped hits (target included).</param>
    /// <param name="maxMismatches">Search radius of the ungapped site scan (≥ 0; default 3 — a library convention kept for
    /// compatibility, no published hybridization source: the figure is the CRISPR guide mismatch tolerance). Whether a
    /// found site is a problem is decided by the Kane criteria; to find every ungapped site that meets the Kane identity
    /// criterion pass ⌈L/4⌉ − 1 (the largest d with (L − d) / L &gt; 0.75).</param>
    /// <param name="selfComplementarityThreshold">Former fold-back-fraction limit of the fallback screen; kept for
    /// source compatibility and copied into <see cref="ProbeParameters.MaxSelfComplementarity"/>, but no longer used by
    /// any screen (the fallback self-dimer limits are <see cref="ProbeParameters.MaxSelfAny"/> / <see cref="ProbeParameters.MaxSelfEnd"/>)
    /// (default 0.3, the Microarray preset's <see cref="ProbeParameters.MaxSelfComplementarity"/>).</param>
    /// <param name="conditions">Hybridization conditions and structure-screen settings (default: Primer3 probe
    /// conditions 50 nM / 50 mM / 0 Mg²⁺ / 0 dNTP, thermodynamic screen, 47 °C — the settings of
    /// <see cref="Defaults.Microarray"/> at those conditions); its <see cref="ProbeParameters.MaxSelfComplementarity"/> is replaced by
    /// <paramref name="selfComplementarityThreshold"/>. Used: the salt / dNTP / oligo concentrations (ntthal screen and
    /// non-target duplex Tm; Primer3 <c>_pr_data_control</c> legality — monovalent and oligo &gt; 0, Mg²⁺ and dNTP ≥ 0,
    /// else <see cref="ArgumentOutOfRangeException"/>), <see cref="ProbeParameters.StructureScreen"/>,
    /// <see cref="ProbeParameters.MaxStructureTm"/>, <see cref="ProbeParameters.MaxSelfAny"/> / <see cref="ProbeParameters.MaxSelfEnd"/>,
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (THAL_MAX_ALIGN of the screen and of the non-target duplex
    /// Tm; 60–10 000, else <see cref="ArgumentOutOfRangeException"/>).</param>
    /// <param name="nonTargetSequences">Optional known non-target sequences for the Kane assessment.</param>
    /// <param name="maxNonTargetIdentity">Kane identity threshold in [0, 1] (default 0.75; flagged when strictly above),
    /// for the reference sites and the non-targets; outside [0, 1] or NaN → <see cref="ArgumentOutOfRangeException"/>.</param>
    /// <param name="maxContiguousMatch">Kane contiguous-identity threshold in nt (≥ 0, default 15; flagged when strictly
    /// longer), for the reference sites and the non-targets; negative → <see cref="ArgumentOutOfRangeException"/>.</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style off-target duplex-Tm threshold (°C; see
    /// <see cref="AssessCrossHybridization"/>); null = Kane criteria only.</param>
    /// <param name="bothStrands">Also scan the references for the probe's reverse complement (default false: the
    /// references' given strand only, the behaviour before audit round 4, B07 F60). A site found in both orientations
    /// at the same reference position is counted once and is cross-hybridizing when either orientation meets a Kane
    /// criterion. The non-target assessment is two-stranded regardless (<see cref="AssessCrossHybridization"/>).</param>
    public static ProbeValidation ValidateProbe(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        int maxMismatches = 3,
        double selfComplementarityThreshold = 0.3,
        ProbeParameters? conditions = null,
        IEnumerable<string>? nonTargetSequences = null,
        double maxNonTargetIdentity = KaneMaxIdentity,
        int maxContiguousMatch = KaneMaxContiguousMatch,
        double? maxDuplexTm = null,
        bool bothStrands = false)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);
        if (double.IsNaN(maxNonTargetIdentity) || maxNonTargetIdentity < 0 || maxNonTargetIdentity > 1)
            throw new ArgumentOutOfRangeException(nameof(maxNonTargetIdentity), "Identity threshold must be in [0, 1].");
        if (maxContiguousMatch < 0)
            throw new ArgumentOutOfRangeException(nameof(maxContiguousMatch), "Contiguous-match threshold cannot be negative.");
        if (conditions is { } stated)
        {
            PrimerDesigner.ValidatePrimer3Conditions(stated.MonovalentMillimolar, stated.DivalentMillimolar,
                stated.DntpMillimolar, stated.DnaConcentrationNanomolar, nameof(conditions));
            ValidateThermodynamicScreenMaxLength(stated, nameof(conditions));
        }

        probeSequence = probeSequence.ToUpperInvariant();
        var issues = new List<string>();

        // Empty probe is a degenerate input — cannot hybridize specifically
        if (probeSequence.Length == 0)
        {
            return new ProbeValidation(
                IsValid: false,
                SpecificityScore: 0.0,
                OffTargetHits: 0,
                SelfComplementarity: 0.0,
                HasSecondaryStructure: false,
                Issues: new List<string> { "Empty probe sequence" });
        }

        int offTargetHits = 0;
        int crossHybridizingHits = 0;

        // Canonical ungapped k-mismatch (Hamming) scan; case-insensitive, overlapping hits included. Each hit is a
        // candidate binding site; it counts as a cross-hybridizing site when it meets a Kane et al. (2000) criterion.
        // With bothStrands the probe's reverse complement is scanned too (the probe's sites on the other strand of a
        // double-stranded reference, as blastn -strand both / CheckSpecificity(bothStrands)); a reference position hit
        // in both orientations (a reverse-palindromic probe, or a site that is itself a reverse palindrome) is one
        // site, cross-hybridizing when either orientation meets a Kane criterion (audit round 4, A4-2, F60).
        string[] patterns = [probeSequence];
        if (bothStrands)
        {
            string reverseComplement = DnaSequence.GetReverseComplementString(probeSequence);
            if (!string.Equals(reverseComplement, probeSequence, StringComparison.Ordinal))
                patterns = [probeSequence, reverseComplement];
        }
        foreach (var reference in referenceSequences)
        {
            var sites = new Dictionary<int, bool>();
            foreach (string pattern in patterns)
            {
                foreach (var hit in ApproximateMatcher.FindWithMismatches(reference, pattern, maxMismatches))
                {
                    bool kane = UngappedSiteMeetsKaneCriteria(hit.MismatchPositions, probeSequence.Length, maxNonTargetIdentity, maxContiguousMatch);
                    sites[hit.Position] = sites.TryGetValue(hit.Position, out bool seen) ? seen || kane : kane;
                }
            }
            offTargetHits += sites.Count;
            crossHybridizingHits += sites.Values.Count(kane => kane);
        }

        // More than one site meeting the Kane criteria = the probe can bind somewhere besides its intended site.
        if (crossHybridizingHits > 1)
        {
            issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"{crossHybridizingHits} potential off-target sites (Kane 2000: identity > {maxNonTargetIdentity * 100:0.##}% or > {maxContiguousMatch} contiguous identical nt)"));
        }

        // Self-structure: the DesignProbes screen (Primer3 ntthal for ACGT probes ≤ ThermodynamicScreenMaxLength; otherwise Primer3
        // alignment-mode self_any / self_end + the inverted-repeat hairpin stem).
        var param = (conditions ?? Primer3ProbeConditions) with { MaxSelfComplementarity = selfComplementarityThreshold };
        double selfComp = CalculateSelfComplementarity(probeSequence);
        double alnSelfAny = PrimerDesigner.CalculatePrimerSelfAnyComplementarity(probeSequence);
        double alnSelfEnd = PrimerDesigner.CalculatePrimerSelfEndComplementarity(probeSequence);
        var thermo = ComputeThermodynamicSelfStructure(probeSequence, param);
        bool selfCompIssue;
        bool hasStructure;
        if (thermo is { } t)
        {
            double max = param.MaxStructureTm;
            selfCompIssue = t.SelfAnyTh > max || t.SelfEndTh > max;
            hasStructure = t.HairpinTh > max;
            if (selfCompIssue)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: ntthal self-dimer Tm {Math.Max(t.SelfAnyTh, t.SelfEndTh):F1}°C exceeds {max:0.##}°C"));
            if (hasStructure)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Potential secondary structure formation: ntthal hairpin Tm {t.HairpinTh:F1}°C exceeds {max:0.##}°C"));
        }
        else
        {
            (selfCompIssue, string? selfCompWarning) = AlignmentSelfDimerScreen(alnSelfAny, alnSelfEnd, param);
            hasStructure = HasSecondaryStructurePotential(probeSequence);
            if (selfCompIssue)
                issues.Add(selfCompWarning!);
            if (hasStructure)
                issues.Add("Potential secondary structure formation");
        }

        // Kane et al. (2000) cross-hybridization criteria against known non-targets (optional).
        IReadOnlyList<CrossHybridizationAssessment> cross = Array.Empty<CrossHybridizationAssessment>();
        if (nonTargetSequences is not null)
        {
            cross = AssessCrossHybridization(probeSequence, nonTargetSequences, maxNonTargetIdentity, maxContiguousMatch,
                conditions: param, maxDuplexTm: maxDuplexTm);
            foreach (var c in cross.Where(c => c.CrossHybridizes))
            {
                issues.Add($"Cross-hybridization risk with non-target {c.NonTargetIndex}"
                    + (c.ReverseComplementStrand ? " (reverse complement)" : "")
                    + string.Create(System.Globalization.CultureInfo.InvariantCulture, $": identity {c.Identity * 100:F0}%, longest contiguous match {c.LongestContiguousMatch} nt (Kane 2000)")
                    + (c.ExceedsDuplexTmThreshold
                        ? string.Create(System.Globalization.CultureInfo.InvariantCulture, $", site duplex Tm {c.DuplexTm:F1}°C > {maxDuplexTm:0.##}°C")
                        : ""));
            }
        }

        // Library uniqueness score: 0 hits → 0, N ≥ 1 hits → 1/N.
        double specificity = offTargetHits == 0 ? 0.0 : 1.0 / offTargetHits;

        return new ProbeValidation(
            issues.Count == 0,
            specificity,
            offTargetHits,
            selfComp,
            hasStructure,
            issues)
        {
            ThermodynamicScreen = thermo is not null,
            SelfDimerTm = thermo?.SelfAnyTh,
            SelfEndDimerTm = thermo?.SelfEndTh,
            HairpinTm = thermo?.HairpinTh,
            SelfAny = alnSelfAny,
            SelfEnd = alnSelfEnd,
            CrossHybridization = cross,
            CrossHybridizingHits = crossHybridizingHits,
        };
    }

    // Kane et al. (2000) criteria applied to an ungapped k-mismatch site of a probe of length L with mismatches at the
    // given (ascending, probe-relative) positions: identity (L − d) / L over the probe length > maxIdentity, or the
    // longest run of identical positions > maxContiguousMatch (the same measures AssessCrossHybridization takes from the
    // best local alignment / longest common substring, here on the ungapped diagonal of the hit).
    private static bool UngappedSiteMeetsKaneCriteria(
        IReadOnlyList<int> mismatchPositions, int probeLength, double maxIdentity, int maxContiguousMatch)
    {
        double identity = (double)(probeLength - mismatchPositions.Count) / probeLength;
        if (identity > maxIdentity)
            return true;
        int longest = 0, previous = -1;
        foreach (int position in mismatchPositions)
        {
            longest = Math.Max(longest, position - previous - 1);
            previous = position;
        }
        longest = Math.Max(longest, probeLength - previous - 1);
        return longest > maxContiguousMatch;
    }

    /// <summary>
    /// Exact-match uniqueness of a probe in a suffix-tree-indexed genome: N = number of occurrences of the probe
    /// (and, with <paramref name="bothStrands"/>, of its reverse complement — the probe's binding sites on the
    /// other strand of a double-stranded genome; a reverse-palindromic probe is counted once); returns 0 for
    /// N = 0, else 1/N (the library uniqueness score of <see cref="ValidateProbe"/>).
    /// </summary>
    /// <param name="probeSequence">Probe sequence (case-insensitive; matched against the index as upper case).</param>
    /// <param name="genomeIndex">Suffix tree of the genome (built on upper-case text).</param>
    /// <param name="bothStrands">Also count reverse-complement occurrences (default false: indexed strand only).</param>
    public static double CheckSpecificity(
        string probeSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        bool bothStrands = false)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(genomeIndex);
        probeSequence = probeSequence.ToUpperInvariant();

        // Check if probe sequence exists in genome
        int hitCount = genomeIndex.CountOccurrences(probeSequence);
        if (bothStrands)
        {
            string rc = DnaSequence.GetReverseComplementString(probeSequence);
            if (!string.Equals(rc, probeSequence, StringComparison.Ordinal))
                hitCount += genomeIndex.CountOccurrences(rc);
        }

        if (hitCount == 0)
            return 0; // Probe doesn't match target

        return 1.0 / hitCount;
    }

    // --- Kane et al. (2000) cross-hybridization criteria ---
    // Kane MD et al. (2000) Nucleic Acids Res 28(22):4552-4557 (abstract): non-target transcripts ">75% similar
    // over the 50 base target may show cross-hybridization", and a non-target region must not include a stretch
    // of identical sequence ">15 contiguous bases"; summarised by later probe-design pipelines as "a probe is
    // likely to cross-hybridize with a nontarget if overall sequence identity is > 75% or if there is a
    // contiguous match > 15 bp".
    private const double KaneMaxIdentity = 0.75;
    private const int KaneMaxContiguousMatch = 15;

    // Non-target strands longer than this are aligned in overlapping chunks (memory O(probe × chunk)).
    private const int CrossHybridizationChunkLength = 4096;

    /// <summary>
    /// Kane et al. (2000) cross-hybridization criteria of a probe against known non-target sequences: for each
    /// non-target strand the probe may cross-hybridize when (a) its overall identity — identical columns of the
    /// best local alignment ÷ probe length — is strictly above <paramref name="maxIdentity"/> (Kane: &gt; 75 %
    /// similar over the probe), or (b) the longest stretch identical between probe and non-target is strictly
    /// longer than <paramref name="maxContiguousMatch"/> (Kane: &gt; 15 contiguous bases).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The local alignment is the canonical Smith–Waterman–Gotoh aligner
    /// <see cref="SequenceAligner.LocalAlignAffine(string, string, ScoringMatrix?)"/> with BLAST+ blastn scoring
    /// (<see cref="SequenceAligner.BlastDna"/>: +2/−3, gap existence 5, extension 2) by default — the reported
    /// optimal alignment is one of those Biopython's local <c>PairwiseAligner</c> enumerates. Non-target strands
    /// longer than 4096 nt are aligned in chunks overlapping by the longest reference span a positive-scoring
    /// local alignment can have (m + ⌈m·match/|extend|⌉), so the best score equals the whole-strand score; among
    /// equal-scoring chunks the first is reported. The longest contiguous match is the longest common substring
    /// from the canonical suffix tree (<c>SuffixTree.LongestCommonSubstringInfo</c>). Comparison is
    /// case-insensitive and literal (an N matches only N).
    /// </para>
    /// <para>
    /// With <paramref name="bothStrands"/> (default) each non-target is assessed as given and as its reverse
    /// complement, since a double-stranded non-target (genomic DNA, ds cDNA) offers the probe both strands;
    /// pass false for single-stranded non-targets given in the probe's sense.
    /// </para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence (non-empty).</param>
    /// <param name="nonTargetSequences">Known non-target sequences (null entries are treated as empty).</param>
    /// <param name="maxIdentity">Identity threshold in [0, 1] (default 0.75).</param>
    /// <param name="maxContiguousMatch">Contiguous-identity threshold in nt (≥ 0; default 15).</param>
    /// <param name="bothStrands">Assess the reverse complement of each non-target too (default true).</param>
    /// <param name="scoring">Local-alignment scoring (default <see cref="SequenceAligner.BlastDna"/>, affine gaps).</param>
    /// <param name="conditions">Hybridization conditions for the site duplex Tm
    /// (<see cref="CrossHybridizationAssessment.DuplexTm"/>; default Primer3 probe conditions 50 nM / 50 mM / 0 / 0).
    /// Its <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> is the ntthal <c>THAL_MAX_ALIGN</c> of the
    /// duplex (default 60; opt-in up to 10 000 for longer probes and sites — audit round 3, A3-27).</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style specificity threshold (°C): a strand whose site duplex Tm
    /// is strictly above it is flagged (<see cref="CrossHybridizationAssessment.ExceedsDuplexTmThreshold"/>); null
    /// (default) applies only the Kane criteria. The threshold is assay-specific (Rouillard et al. 2003: user-set).</param>
    /// <returns>One assessment per non-target strand: index order, forward strand before reverse complement.</returns>
    /// <exception cref="ArgumentNullException">A null probe or non-target collection.</exception>
    /// <exception cref="ArgumentException">An empty probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxIdentity"/> outside [0, 1], a negative
    /// <paramref name="maxContiguousMatch"/>, or a <paramref name="conditions"/>
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> outside 60–10 000.</exception>
    public static IReadOnlyList<CrossHybridizationAssessment> AssessCrossHybridization(
        string probeSequence,
        IEnumerable<string> nonTargetSequences,
        double maxIdentity = KaneMaxIdentity,
        int maxContiguousMatch = KaneMaxContiguousMatch,
        bool bothStrands = true,
        ScoringMatrix? scoring = null,
        ProbeParameters? conditions = null,
        double? maxDuplexTm = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(nonTargetSequences);
        if (probeSequence.Length == 0)
            throw new ArgumentException("Probe sequence cannot be empty.", nameof(probeSequence));
        if (double.IsNaN(maxIdentity) || maxIdentity < 0 || maxIdentity > 1)
            throw new ArgumentOutOfRangeException(nameof(maxIdentity), "Identity threshold must be in [0, 1].");
        if (maxContiguousMatch < 0)
            throw new ArgumentOutOfRangeException(nameof(maxContiguousMatch), "Contiguous-match threshold cannot be negative.");

        if (conditions is { } stated)
            ValidateThermodynamicScreenMaxLength(stated, nameof(conditions));

        var matrix = scoring ?? SequenceAligner.BlastDna;
        var cond = conditions ?? Primer3ProbeConditions;
        string probe = probeSequence.ToUpperInvariant();
        // thal.c thal_check_errors: at least one strand ≤ THAL_MAX_ALIGN (ThermodynamicScreenMaxLength, default 60)
        // and neither longer than THAL_MAX_SEQ (10 000).
        int maxAlign = cond.ThermodynamicScreenMaxLength;
        bool duplexComputable = probe.Length <= PrimerDesigner.NtthalMaxSequenceLength
            && probe.All(c => c is 'A' or 'C' or 'G' or 'T');
        var probeTree = global::SuffixTree.SuffixTree.Build(probe);
        var result = new List<CrossHybridizationAssessment>();

        int index = 0;
        foreach (var nonTarget in nonTargetSequences)
        {
            string forward = (nonTarget ?? string.Empty).ToUpperInvariant();
            result.Add(AssessStrand(index, false, forward));
            if (bothStrands)
                result.Add(AssessStrand(index, true, DnaSequence.GetReverseComplementString(forward)));
            index++;
        }

        return result;

        CrossHybridizationAssessment AssessStrand(int idx, bool reverse, string strand)
        {
            var (score, identical, siteStart, siteEnd) = BestLocalAlignment(probe, strand, matrix);
            int contiguous = strand.Length == 0 ? 0 : probeTree.LongestCommonSubstringInfo(strand).Substring.Length;
            double identity = (double)identical / probe.Length;
            double? duplexTm = null;
            if (duplexComputable && siteStart >= 0)
            {
                string site = strand.Substring(siteStart, siteEnd - siteStart + 1);
                if (Math.Min(probe.Length, site.Length) <= maxAlign
                    && site.Length <= PrimerDesigner.NtthalMaxSequenceLength
                    && site.All(c => c is 'A' or 'C' or 'G' or 'T'))
                {
                    // The probe hybridizes to the strand complementary to the site (primer3-py calc_heterodimer).
                    var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal(
                        probe, DnaSequence.GetReverseComplementString(site), PrimerDesigner.NtthalAlignmentMode.Any,
                        cond.MonovalentMillimolar / 1000.0, cond.DivalentMillimolar / 1000.0, cond.DntpMillimolar / 1000.0,
                        cond.DnaConcentrationNanomolar * 1e-9, PrimerDesigner.NtthalDefaultTemperatureCelsius,
                        PrimerDesigner.NtthalDefaultMaxLoop, maxAlign);
                    duplexTm = d?.TmCelsius ?? 0.0;
                }
            }

            return new CrossHybridizationAssessment(
                idx, reverse, identity, identical, score, contiguous,
                identity > maxIdentity, contiguous > maxContiguousMatch)
            {
                SiteStart = siteStart,
                SiteEnd = siteEnd,
                DuplexTm = duplexTm,
                ExceedsDuplexTmThreshold = maxDuplexTm is double limit && duplexTm is double tm && tm > limit,
            };
        }
    }

    // Best local alignment (score, identical columns) of the probe with a strand via the canonical affine
    // Smith–Waterman–Gotoh aligner, chunking long strands so that every positive-scoring alignment lies wholly
    // inside one chunk (its reference span is at most m + ⌈m·match/|extend|⌉ residues).
    private static (int Score, int Identical, int Start, int End) BestLocalAlignment(
        string probe, string strand, ScoringMatrix matrix)
    {
        if (strand.Length == 0)
            return (0, 0, -1, -1);

        int m = probe.Length;
        int overlap = matrix.GapExtend < 0 && matrix.Match > 0
            ? m + (int)Math.Ceiling((double)m * matrix.Match / -matrix.GapExtend) + 1
            : int.MaxValue;
        int chunk = Math.Max(CrossHybridizationChunkLength, 2 * Math.Min(overlap, int.MaxValue / 4));
        if (overlap == int.MaxValue || strand.Length <= chunk)
            return Summarize(SequenceAligner.LocalAlignAffine(probe, strand, matrix), 0);

        (int Score, int Identical, int Start, int End) best = (-1, 0, -1, -1);
        int step = chunk - overlap;
        bool last = false;
        int start = 0;
        while (!last)
        {
            int len = Math.Min(chunk, strand.Length - start);
            last = start + len >= strand.Length;
            var r = Summarize(SequenceAligner.LocalAlignAffine(probe, strand.Substring(start, len), matrix), start);
            if (r.Score > best.Score)
                best = r;

            start += step;
        }

        return best;

        static (int Score, int Identical, int Start, int End) Summarize(AlignmentResult aln, int offset)
        {
            int identical = 0;
            for (int k = 0; k < aln.AlignedSequence1.Length; k++)
            {
                char c1 = aln.AlignedSequence1[k];
                if (c1 != AlignmentGapChar && c1 == aln.AlignedSequence2[k])
                    identical++;
            }

            return aln.Score > 0
                ? (aln.Score, identical, offset + aln.StartPosition2, offset + aln.EndPosition2)
                : (0, 0, -1, -1);
        }
    }

    // --- Gapped off-target scan thresholds (sourced) ---

    // Default minimum identity (identical aligned columns / probe length) to call an off-target.
    // Kane et al. (2000, Nucleic Acids Res. 28(22):4552-4557): "for a given oligonucleotide probe
    // any 'non-target' transcripts (cDNAs) >75% similar over the [...] target may show
    // cross-hybridization." 0.75 is therefore the empirically-grounded similarity threshold above
    // which cross-hybridization (an off-target) becomes a concern; callers may override it.
    private const double DefaultOffTargetMinIdentity = 0.75;

    // Extra reference length scanned per window beyond the probe length, to allow a few indels in
    // the local alignment (a gap shifts the reference frame). Two extra bases lets the Smith-Waterman
    // local alignment absorb short insertions/deletions while keeping each window O(probeLen) wide.
    private const int GappedScanGapAllowance = 2;

    // Gap character emitted by SequenceAligner in its aligned-output strings.
    private const char AlignmentGapChar = '-';

    /// <summary>
    /// Opt-in <b>gapped</b> off-target scan using the library's validated Smith-Waterman local
    /// aligner (<see cref="SequenceAligner.LocalAlign(string, string, ScoringMatrix?)"/>). Unlike the
    /// ungapped Hamming-distance scan in <see cref="ValidateProbe"/> (which only tolerates
    /// substitutions in a fixed-length window), this finds off-target sites reachable through
    /// insertions or deletions — the "BLAST-grade" improvement (gapped local alignment handles
    /// indels the ungapped scan misses; Altschul et al. 1990; Smith &amp; Waterman 1981).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The default ungapped <see cref="ValidateProbe"/> behaviour is unchanged; this is a separate,
    /// additive entry point. It also corrects the on/off-target pooling of
    /// <see cref="ProbeValidation.OffTargetHits"/>: the single intended on-target site (the perfect,
    /// ungapped, full-coverage exact match) is reported separately and is NOT counted as an off-target.
    /// </para>
    /// <para>
    /// This is an exhaustive sliding Smith-Waterman scan (O(g · n · m) over reference length g and
    /// probe length n·m), not a seeded BLAST index over a whole genome; for genome-scale search a
    /// k-mer/seed index would be required.
    /// </para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence to scan (5'→3'). Null throws; empty yields no hits.</param>
    /// <param name="referenceSequences">Reference sequences to scan for hits. Null throws.</param>
    /// <param name="minIdentity">
    /// Minimum alignment identity (identical aligned columns / probe length) to report a hit (hits with
    /// identity ≥ <paramref name="minIdentity"/> are reported). Default 0.75 from Kane et al. (2000): non-targets
    /// &gt; 75 % similar over the probe may cross-hybridize. This scan reports sites; the Kane decision rule itself
    /// (identity strictly &gt; 75 % or a contiguous identical stretch &gt; 15 nt, per non-target, both strands) is
    /// <see cref="AssessCrossHybridization"/>.
    /// </param>
    /// <param name="scoring">
    /// Scoring matrix for the local alignment. Defaults to <see cref="SequenceAligner.BlastDna"/>
    /// (+2/-3, gap -2) — the same BLAST-style DNA scoring reused for the library's gapped ANI alignment.
    /// </param>
    /// <returns>A <see cref="GappedSpecificityResult"/> separating on-target from off-target hits.</returns>
    public static GappedSpecificityResult ScanOffTargetsGapped(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        double minIdentity = DefaultOffTargetMinIdentity,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);

        var onTarget = new List<GappedProbeHit>();
        var offTarget = new List<GappedProbeHit>();

        string probe = probeSequence.ToUpperInvariant();
        if (probe.Length == 0)
            return new GappedSpecificityResult(onTarget, offTarget, minIdentity);

        var matrix = scoring ?? SequenceAligner.BlastDna;
        bool onTargetClaimed = false;

        int refIndex = 0;
        foreach (var reference in referenceSequences)
        {
            string text = (reference ?? string.Empty).ToUpperInvariant();
            foreach (var hit in ScanReferenceGapped(probe, text, refIndex, minIdentity, matrix))
            {
                // The intended on-target is the (first) perfect ungapped full-coverage exact match.
                bool isPerfectExact = !hit.HasGaps
                    && hit.Identity >= 1.0
                    && hit.Coverage >= 1.0;

                if (isPerfectExact && !onTargetClaimed)
                {
                    onTargetClaimed = true;
                    onTarget.Add(hit);
                }
                else if (isPerfectExact)
                {
                    // An additional perfect repeat: report it as an on-target-class match for
                    // visibility, but it is still a genuine extra binding site → counts as off-target.
                    onTarget.Add(hit);
                    offTarget.Add(hit);
                }
                else
                {
                    offTarget.Add(hit);
                }
            }

            refIndex++;
        }

        return new GappedSpecificityResult(onTarget, offTarget, minIdentity);
    }

    /// <summary>
    /// Slides a Smith-Waterman local alignment of <paramref name="probe"/> across
    /// <paramref name="text"/> and returns one best non-overlapping hit per distinct site whose
    /// identity is at least <paramref name="minIdentity"/>.
    /// </summary>
    private static List<GappedProbeHit> ScanReferenceGapped(
        string probe, string text, int refIndex, double minIdentity, ScoringMatrix matrix)
    {
        int probeLen = probe.Length;
        int windowLen = probeLen + GappedScanGapAllowance;
        var raw = new List<GappedProbeHit>();

        for (int start = 0; start < text.Length; start++)
        {
            int span = Math.Min(windowLen, text.Length - start);
            if (span <= 0)
                break;

            string window = text.Substring(start, span);
            AlignmentResult aln = SequenceAligner.LocalAlign(probe, window, matrix);

            // SequenceAligner.LocalAlign aligns (sequence1 = probe, sequence2 = window):
            // AlignedSequence1 is the probe side, AlignedSequence2 the reference side.
            string ap = aln.AlignedSequence1;
            string ar = aln.AlignedSequence2;
            if (ap.Length == 0)
                continue;

            int identical = 0;
            int ungapped = 0;
            bool hasGap = false;
            for (int k = 0; k < ap.Length; k++)
            {
                char c1 = ap[k];
                char c2 = ar[k];
                if (c1 == AlignmentGapChar || c2 == AlignmentGapChar)
                {
                    hasGap = true;
                    continue;
                }

                ungapped++;
                if (c1 == c2)
                    identical++;
            }

            double identity = (double)identical / probeLen;
            if (identity < minIdentity)
                continue;

            double coverage = (double)ungapped / probeLen;
            int absStart = start + aln.StartPosition2;
            int absEnd = start + aln.EndPosition2;

            raw.Add(new GappedProbeHit(
                refIndex, absStart, absEnd, identity, coverage, hasGap, ap, ar));
        }

        // Greedy best-per-site selection: take highest-identity (then highest-coverage, then
        // leftmost) hits first, accepting a hit only if its reference span does not overlap an
        // already-accepted one. Overlapping windows that re-detect the same site collapse to one hit.
        raw.Sort((a, b) =>
        {
            int byIdentity = b.Identity.CompareTo(a.Identity);
            if (byIdentity != 0) return byIdentity;
            int byCoverage = b.Coverage.CompareTo(a.Coverage);
            if (byCoverage != 0) return byCoverage;
            return a.Start.CompareTo(b.Start);
        });

        var accepted = new List<GappedProbeHit>();
        foreach (var hit in raw)
        {
            bool overlaps = accepted.Any(a => !(hit.End < a.Start || hit.Start > a.End));
            if (!overlaps)
                accepted.Add(hit);
        }

        accepted.Sort((a, b) => a.Start.CompareTo(b.Start));
        return accepted;
    }

    // --- Karlin–Altschul statistics (opt-in) ---
    //
    // Karlin & Altschul (1990, PNAS 87:2264) / Altschul et al. (1990, J Mol Biol 215:403):
    //   E = K·m·n·e^{−λS}                       (expected HSPs with score ≥ S by chance)
    //   S' = (λS − ln K) / ln 2                  (normalized "bit" score)
    //   E = m·n·2^{−S'}                          (E in terms of the bit score)
    //   λ is the unique positive root of  Σ_{i,j} p_i p_j e^{λ s_ij} = 1.
    // Verbatim formulas retrieved 2026-06-24 from the NCBI BLAST course
    //   "The Statistics of Sequence Similarity Scores" (Altschul),
    //   https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html, and from Durand,
    //   "BLAST (Karlin–Altschul) Statistics" (CMU 03-711, citing Karlin & Altschul 1990 and
    //   Altschul et al. 1990), http://www.cs.cmu.edu/~durand/03-711/2011/Lectures/Blast-informationContent-2011.pdf.
    // The theory requires a scoring scheme whose expected per-pair score is negative and that
    // has at least one positive score (Altschul, ibid.); otherwise λ is undefined.

    // Uniform nucleotide background frequency p_i = 0.25 for the four bases A,C,G,T
    // (the standard assumption when computing the +1/−3 nucleotide λ ≈ 1.374, NCBI blastn).
    private const double UniformBaseFrequency = 0.25;

    // Bisection bounds/iterations for solving Σ p_i p_j e^{λ s_ij} = 1. The function is strictly
    // increasing in λ once it crosses 1 (the expected score is negative, so it starts < 1 and a
    // positive score makes it diverge), so a simple bisection converges to the unique positive root.
    private const double LambdaSearchUpperBound = 100.0;
    private const int LambdaBisectionIterations = 200;

    /// <summary>
    /// Computes the Karlin–Altschul scale parameter λ for a match/mismatch nucleotide scoring
    /// scheme under uniform (0.25) base frequencies, by solving the defining equation
    /// Σ_{i,j} p_i p_j e^{λ s_ij} = 1 numerically (bisection on the unique positive root).
    /// </summary>
    /// <remarks>
    /// With four equiprobable bases and a simple match/mismatch matrix, the 16 ordered pairs split
    /// into 4 matches (probability 4·0.25² = 0.25) and 12 mismatches (probability 0.75), so the
    /// equation reduces to 0.25·e^{λ·match} + 0.75·e^{λ·mismatch} = 1. For the BLAST +1/−3 scheme
    /// this yields λ ≈ 1.374 (NCBI blastn). Karlin &amp; Altschul (1990); Altschul et al. (1990).
    /// </remarks>
    /// <param name="match">Match score (must be &gt; 0 — the required positive score).</param>
    /// <param name="mismatch">Mismatch score (must be &lt; 0).</param>
    /// <param name="baseFrequency">
    /// Per-base background frequency (default 0.25, uniform — the standard nucleotide composition NCBI blastn uses).
    /// The model takes p(match) = 4·p² and p(mismatch) = 1 − 4·p², which is the composition-exact value only at
    /// p = 0.25; for a non-uniform composition use
    /// <see cref="ComputeUngappedKarlinParameters(int, int, IReadOnlyList{double}, KarlinKMethod)"/>. Must lie in (0, 0.5).
    /// </param>
    /// <returns>The positive λ solving the Karlin–Altschul equation.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when the scheme cannot define λ: the match score is not positive, the mismatch score is
    /// not negative, or the expected per-pair score is not negative (the theory's preconditions); or when
    /// <paramref name="baseFrequency"/> is outside (0, 0.5) (p(match) = 4·p² would not be a probability in (0, 1)).
    /// </exception>
    public static double ComputeLambdaNucleotide(
        int match,
        int mismatch,
        double baseFrequency = UniformBaseFrequency)
    {
        // Karlin–Altschul preconditions: at least one positive score, and negative expected score.
        ValidateMatchMismatch(match, mismatch);

        // p(match) = 4 · p² (the four identical ordered pairs); p(mismatch) = 1 − p(match).
        double pMatch = MatchProbability(baseFrequency);
        double pMismatch = 1.0 - pMatch;

        // Expected per-pair score must be negative for the theory to hold.
        double expectedScore = pMatch * match + pMismatch * mismatch;
        if (expectedScore >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the expected per-pair score must be negative.");

        // f(λ) = p(match)·e^{λ·match} + p(mismatch)·e^{λ·mismatch} − 1: f(0) = 0, f'(0) = expectedScore < 0,
        // and the positive match score drives f → ∞, so f crosses 0 exactly once at the positive root.
        return SolveLambda(pMatch, match, pMismatch, mismatch);
    }

    /// <summary>
    /// Computes the Karlin–Altschul bit score and E-value for an off-target alignment hit's raw
    /// score, given the search-space dimensions and the scoring scheme.
    /// </summary>
    /// <remarks>
    /// <para>
    /// E = K·m·n·e^{−λS}, S' = (λS − ln K) / ln 2, E = m·n·2^{−S'} (Karlin &amp; Altschul 1990;
    /// Altschul et al. 1990). λ is computed from <paramref name="scoring"/> by
    /// <see cref="ComputeLambdaNucleotide"/>; K is, unless the caller supplies it, computed for the same
    /// scoring scheme and background by the Karlin–Altschul (1990) lattice formula as NCBI BLAST+ computes it
    /// (<see cref="ComputeUngappedKarlinParameters(int, int, double, KarlinKMethod)"/>; +1/−3 → K = 0.7106, +2/−3 → K = 0.4081,
    /// the values blastn 2.12 reports for ungapped searches).
    /// </para>
    /// <para>
    /// These are the <b>ungapped</b> statistics on the raw search space m·n. For a gapped (affine) alignment score
    /// such as <see cref="CrossHybridizationAssessment.AlignmentScore"/>, and for BLAST's edge-effect length
    /// correction, use <see cref="ComputeBlastnStatistics(int, int, long, int, ScoringMatrix?, bool, KarlinKMethod)"/>.
    /// </para>
    /// </remarks>
    /// <param name="rawScore">The raw alignment score S of the hit.</param>
    /// <param name="queryLength">Query (probe) length m (&gt; 0).</param>
    /// <param name="databaseLength">Database (reference) length n (&gt; 0).</param>
    /// <param name="scoring">
    /// The scoring scheme. Its <see cref="ScoringMatrix.Match"/>/<see cref="ScoringMatrix.Mismatch"/>
    /// determine λ. Defaults to <see cref="SequenceAligner.BlastDna"/> (+2/−3). Pass a +1/−3 matrix to
    /// reproduce the published λ ≈ 1.374.
    /// </param>
    /// <param name="k">
    /// The Karlin–Altschul K parameter; null (default) computes it for <paramref name="scoring"/> and
    /// <paramref name="baseFrequency"/> (<see cref="ComputeUngappedKarlinParameters(int, int, double, KarlinKMethod)"/>).
    /// </param>
    /// <param name="baseFrequency">Per-base background frequency for λ (default 0.25, uniform).</param>
    /// <param name="kMethod">
    /// How the computed K treats schemes whose scores share a divisor &gt; 1 (default
    /// <see cref="KarlinKMethod.ReducedLattice"/>, scale-invariant; <see cref="KarlinKMethod.NcbiBlast"/> reproduces
    /// blastn's printed K, e.g. 1.17 for +4/−6). Ignored when <paramref name="k"/> is given.
    /// </param>
    /// <returns>The <see cref="KarlinAltschulStatistics"/> for the hit.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="scoring"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown for non-positive lengths or K, or a scheme for which λ is undefined.</exception>
    public static KarlinAltschulStatistics ComputeKarlinAltschul(
        double rawScore,
        int queryLength,
        long databaseLength,
        ScoringMatrix? scoring = null,
        double? k = null,
        double baseFrequency = UniformBaseFrequency,
        KarlinKMethod kMethod = KarlinKMethod.ReducedLattice)
    {
        var matrix = scoring ?? SequenceAligner.BlastDna;
        ArgumentNullException.ThrowIfNull(matrix);
        if (queryLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(queryLength), "Query length m must be positive.");
        if (databaseLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseLength), "Database length n must be positive.");
        if (k is double given && !(given > 0))
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        double lambda = ComputeLambdaNucleotide(matrix.Match, matrix.Mismatch, baseFrequency);
        double kValue = k ?? ComputeUngappedKarlinParameters(matrix.Match, matrix.Mismatch, baseFrequency, kMethod).K;

        // S' = (λS − ln K) / ln 2  (Altschul et al. 1990).
        double bitScore = (lambda * rawScore - Math.Log(kValue)) / Math.Log(2.0);

        // E = K·m·n·e^{−λS}  (Karlin & Altschul 1990).
        double eValue = kValue * queryLength * databaseLength * Math.Exp(-lambda * rawScore);

        return new KarlinAltschulStatistics(
            rawScore, lambda, kValue, bitScore, eValue, queryLength, databaseLength);
    }

    // --- NCBI BLAST+ Karlin–Altschul machinery (blastn statistics) ---
    //
    // Ported from NCBI C++ Toolkit BLAST+ core (retrieved 2026-10-01 from
    // https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/blast/core/):
    //   blast_stat.c  — BlastScoreFreqCalc, Blast_KarlinBlkUngappedCalc (λ, BlastKarlinLtoH, BlastKarlinLHtoK),
    //                   blastn_values_* tables + s_GetNuclValuesArray / s_SplitArrayOf8 / s_AdjustGapParametersByGcd,
    //                   Blast_KarlinBlkNuclGappedCalc, s_GetUngappedBeta, Blast_GetNuclAlphaBeta,
    //                   BLAST_KarlinStoE_simple, BLAST_ComputeLengthAdjustment;
    //   blast_setup.c — BLAST_CalcEffLengths (effective search space);
    //   blast_hits.c  — Blast_HSPListGetEvalues (even-score round-down), Blast_HSPListGetBitScores;
    //   ncbi_math.c   — BLAST_Expm1.
    // Cross-checked against NCBI blastn 2.12.0+ (Lambda/K/H footers, effective search spaces and E-values).

    /// <summary>
    /// Karlin–Altschul parameters of a nucleotide scoring system, as NCBI BLAST+ uses them.
    /// </summary>
    /// <param name="Lambda">Scale parameter λ.</param>
    /// <param name="K">Search-space scale parameter K.</param>
    /// <param name="H">Relative entropy H (nats per aligned pair).</param>
    /// <param name="Alpha">
    /// Edge-effect parameter α (BLAST+ <c>Blast_GetNuclAlphaBeta</c>; ungapped: λ/H). The length adjustment uses α/λ.
    /// </param>
    /// <param name="Beta">Edge-effect parameter β (BLAST+ <c>Blast_GetNuclAlphaBeta</c>).</param>
    /// <param name="RoundDown">
    /// True when BLAST+ rounds odd gapped scores down to the next even score before computing the E-value
    /// (tables for 2/−3, 2/−5, 2/−7, 3/−4: "these parameters can only be applied to even scores").
    /// </param>
    /// <param name="Gapped">True for gapped (affine) parameters, false for ungapped ones.</param>
    public readonly record struct KarlinAltschulParameters(
        double Lambda,
        double K,
        double H,
        double Alpha,
        double Beta,
        bool RoundDown,
        bool Gapped);

    /// <summary>
    /// How the ungapped Karlin–Altschul K is evaluated when the scores share a common divisor δ &gt; 1
    /// (e.g. +4/−6, δ = 2). Schemes with δ = 1 give identical K under both choices.
    /// </summary>
    public enum KarlinKMethod
    {
        /// <summary>
        /// Karlin &amp; Altschul (1990, PNAS 87:2264, appendix) on the δ-reduced lattice: K is invariant under
        /// multiplying all scores by an integer (+4/−6 → K = 0.4081 = K(+2/−3)).
        /// </summary>
        ReducedLattice = 0,

        /// <summary>
        /// Exactly as NCBI BLAST+ <c>BlastKarlinLHtoK</c> (blast_stat.c) computes it: low/high/λ are reduced by δ but
        /// the series reads the score probabilities at the reduced offsets from the unreduced lowest score
        /// (<c>probArrayStartLow[j]</c>, j = 0…range/δ), so K is not scale-invariant. Reproduces the K blastn 2.12.0+
        /// prints: +4/−6 → 1.17, +4/−10 → 1.06, +6/−4 → 1.63, +6/−10 → 1.03, +8/−10 → 1.07, +10/−8 → 1.31.
        /// The closed forms (reduced lowest score −1 or highest +1, e.g. +2/−2, +2/−4, +4/−2) agree with
        /// <see cref="ReducedLattice"/>.
        /// </summary>
        NcbiBlast = 1,
    }

    /// <summary>
    /// NCBI BLAST+ blastn statistics of an alignment score: bit score and E-value over the effective
    /// (edge-corrected) search space.
    /// </summary>
    /// <param name="RawScore">Raw alignment score S.</param>
    /// <param name="EValueScore">The score the E-value uses (S rounded down to even when <see cref="KarlinAltschulParameters.RoundDown"/>).</param>
    /// <param name="Parameters">The Karlin–Altschul parameters used.</param>
    /// <param name="BitScore">S' = (λS − ln K)/ln 2 on the raw score (BLAST+ <c>Blast_HSPListGetBitScores</c>).</param>
    /// <param name="EValue">E = K·(m − ℓ)·(n − N·ℓ)·e^{−λS} (BLAST+ <c>BLAST_KarlinStoE_simple</c>).</param>
    /// <param name="LengthAdjustment">The edge-effect length adjustment ℓ (BLAST+ <c>BLAST_ComputeLengthAdjustment</c>).</param>
    /// <param name="EffectiveSearchSpace">(m − ℓ)·max(1, n − N·ℓ).</param>
    /// <param name="QueryLength">Query (probe) length m.</param>
    /// <param name="DatabaseLength">Total database (subject) length n.</param>
    /// <param name="DatabaseSequenceCount">Number of database sequences N.</param>
    public readonly record struct BlastnStatistics(
        int RawScore,
        int EValueScore,
        KarlinAltschulParameters Parameters,
        double BitScore,
        double EValue,
        int LengthAdjustment,
        double EffectiveSearchSpace,
        int QueryLength,
        long DatabaseLength,
        int DatabaseSequenceCount);

    // BlastKarlinLHtoK: BLAST_KARLIN_K_SUMLIMIT_DEFAULT and BLAST_KARLIN_K_ITER_MAX (blast_stat.c).
    private const double KarlinKSumLimit = 0.0001;
    private const int KarlinKIterMax = 100;

    // BLAST_ComputeLengthAdjustment: kMaxIterations.
    private const int LengthAdjustmentMaxIterations = 20;

    // blastn_values_* (blast_stat.c): rows {gap open, gap extend, λ, K, H, α, β, θ}; a leading {0, 0} row is the
    // non-affine (greedy megablast) entry split off by s_SplitArrayOf8. Keyed by (reward, penalty) after division
    // by their gcd; GapOpenMax/GapExtendMax start the "infinite" gap-cost domain where ungapped values apply.
    private sealed record BlastnValueTable(double[][] Rows, int GapOpenMax, int GapExtendMax, bool RoundDown);

    private static readonly Dictionary<(int Reward, int Penalty), BlastnValueTable> BlastnValueTables = new()
    {
        [(1, -5)] = new([[0, 0, 1.39, 0.747, 1.38, 1.00, 0, 100], [3, 3, 1.39, 0.747, 1.38, 1.00, 0, 100]], 3, 3, false),
        [(1, -4)] = new([
            [0, 0, 1.383, 0.738, 1.36, 1.02, 0, 100], [1, 2, 1.36, 0.67, 1.2, 1.1, 0, 98],
            [0, 2, 1.26, 0.43, 0.90, 1.4, -1, 91], [2, 1, 1.35, 0.61, 1.1, 1.2, -1, 98],
            [1, 1, 1.22, 0.35, 0.72, 1.7, -3, 88]], 2, 2, false),
        [(2, -7)] = new([
            [0, 0, 0.69, 0.73, 1.34, 0.515, 0, 100], [2, 4, 0.68, 0.67, 1.2, 0.55, 0, 99],
            [0, 4, 0.63, 0.43, 0.90, 0.7, -1, 91], [4, 2, 0.675, 0.62, 1.1, 0.6, -1, 98],
            [2, 2, 0.61, 0.35, 0.72, 1.7, -3, 88]], 4, 4, true),
        [(1, -3)] = new([
            [0, 0, 1.374, 0.711, 1.31, 1.05, 0, 100], [2, 2, 1.37, 0.70, 1.2, 1.1, 0, 99],
            [1, 2, 1.35, 0.64, 1.1, 1.2, -1, 98], [0, 2, 1.25, 0.42, 0.83, 1.5, -2, 91],
            [2, 1, 1.34, 0.60, 1.1, 1.2, -1, 97], [1, 1, 1.21, 0.34, 0.71, 1.7, -2, 88]], 2, 2, false),
        [(2, -5)] = new([
            [0, 0, 0.675, 0.65, 1.1, 0.6, -1, 99], [2, 4, 0.67, 0.59, 1.1, 0.6, -1, 98],
            [0, 4, 0.62, 0.39, 0.78, 0.8, -2, 91], [4, 2, 0.67, 0.61, 1.0, 0.65, -2, 98],
            [2, 2, 0.56, 0.32, 0.59, 0.95, -4, 82]], 4, 4, true),
        [(1, -2)] = new([
            [0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96], [2, 2, 1.33, 0.62, 1.1, 1.2, 0, 99],
            [1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97], [0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89],
            [3, 1, 1.32, 0.57, 1.0, 1.3, -1, 99], [2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96],
            [1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85]], 2, 2, false),
        [(2, -3)] = new([
            [0, 0, 0.55, 0.21, 0.46, 1.2, -5, 87], [4, 4, 0.63, 0.42, 0.84, 0.75, -2, 99],
            [2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97], [0, 4, 0.55, 0.21, 0.46, 1.2, -5, 87],
            [3, 3, 0.615, 0.37, 0.68, 0.9, -3, 97], [6, 2, 0.63, 0.42, 0.84, 0.75, -2, 99],
            [5, 2, 0.625, 0.41, 0.78, 0.8, -2, 99], [4, 2, 0.61, 0.35, 0.68, 0.9, -3, 96],
            [2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81]], 6, 4, true),
        [(3, -4)] = new([
            [6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95], [5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92],
            [4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86], [6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88],
            [5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81], [4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69]], 6, 3, true),
        [(1, -1)] = new([
            [3, 2, 1.09, 0.31, 0.55, 2.0, -2, 99], [2, 2, 1.07, 0.27, 0.49, 2.2, -3, 97],
            [1, 2, 1.02, 0.21, 0.36, 2.8, -6, 92], [0, 2, 0.80, 0.064, 0.17, 4.8, -16, 72],
            [4, 1, 1.08, 0.28, 0.54, 2.0, -2, 98], [3, 1, 1.06, 0.25, 0.46, 2.3, -4, 96],
            [2, 1, 0.99, 0.17, 0.30, 3.3, -10, 90]], 4, 2, false),
        [(3, -2)] = new([[5, 5, 0.208, 0.030, 0.072, 2.9, -47, 77]], 5, 5, false),
        [(4, -5)] = new([
            [0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74], [6, 5, 0.28, 0.21, 0.47, 0.6, -7, 93],
            [5, 5, 0.27, 0.17, 0.39, 0.7, -9, 90], [4, 5, 0.25, 0.10, 0.31, 0.8, -10, 83],
            [3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76]], 12, 8, false),
        [(5, -4)] = new([[10, 6, 0.163, 0.068, 0.16, 1.0, -19, 85], [8, 6, 0.146, 0.039, 0.11, 1.3, -29, 76]], 25, 10, false),
    };

    /// <summary>
    /// Ungapped Karlin–Altschul parameters λ, K and H of a match/mismatch nucleotide scoring scheme under the
    /// standard uniform composition (or a given per-base frequency), computed as NCBI BLAST+ computes them
    /// (<c>Blast_KarlinBlkUngappedCalc</c>): λ is the positive root of Σ p_s e^{λs} = 1, H = λ·Σ s·p_s·e^{λs}
    /// (<c>BlastKarlinLtoH</c>), and K follows Karlin &amp; Altschul (1990, PNAS 87:2264, eq. and appendix) on the
    /// score lattice of span δ = gcd of the scores (<c>BlastKarlinLHtoK</c>: closed forms when the lattice's lowest
    /// score is −1 or its highest +1, otherwise the convergent series over gapless-alignment score distributions).
    /// </summary>
    /// <remarks>
    /// Reproduces the Lambda/K/H that NCBI blastn 2.12.0+ prints for ungapped searches: +1/−3 → 1.374 / 0.711 / 1.31,
    /// +2/−3 → 0.634 / 0.408 / 0.912, +1/−2 → 1.33 / 0.621 / 1.12. By default (<see cref="KarlinKMethod.ReducedLattice"/>)
    /// K is computed on the reduced lattice, so it is invariant under multiplying all scores by an integer (λ scales by
    /// its inverse) as the theory requires; BLAST+ indexes the series' probabilities by the reduced offset from the
    /// unreduced lowest score when δ &gt; 1 (blastn 2.12.0+ prints K = 1.17 for +4/−6 although +2/−3 gives 0.408) —
    /// <see cref="KarlinKMethod.NcbiBlast"/> reproduces that value exactly. <see cref="KarlinAltschulParameters.Alpha"/>
    /// = λ/H and <see cref="KarlinAltschulParameters.Beta"/> = BLAST+ <c>s_GetUngappedBeta</c> (−2 for +1/−1 and
    /// +2/−3, else 0).
    /// </remarks>
    /// <param name="match">Match score (&gt; 0).</param>
    /// <param name="mismatch">Mismatch score (&lt; 0).</param>
    /// <param name="baseFrequency">Per-base frequency, p(match) = 4·p² (default 0.25; must lie in (0, 0.5)).</param>
    /// <param name="kMethod">
    /// K for schemes whose scores share a divisor δ &gt; 1: <see cref="KarlinKMethod.ReducedLattice"/> (default,
    /// scale-invariant per Karlin &amp; Altschul 1990) or <see cref="KarlinKMethod.NcbiBlast"/> (BLAST+ indexing,
    /// as blastn prints it: +4/−6 → 1.1666856431064105). Identical for δ = 1.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException">λ undefined (no positive score, non-negative mismatch or
    /// non-negative expected score), a base frequency outside (0, 0.5), or an undefined <paramref name="kMethod"/>.</exception>
    public static KarlinAltschulParameters ComputeUngappedKarlinParameters(
        int match, int mismatch, double baseFrequency = UniformBaseFrequency,
        KarlinKMethod kMethod = KarlinKMethod.ReducedLattice)
    {
        ValidateMatchMismatch(match, mismatch);
        return UngappedKarlinParameters(match, mismatch, MatchProbability(baseFrequency), kMethod);
    }

    /// <summary>
    /// Ungapped Karlin–Altschul parameters for a match/mismatch scheme under an arbitrary base composition
    /// (Σ_{i,j} p_i p_j e^{λ s_ij} = 1 with p(match) = Σ p_i²), computed as in
    /// <see cref="ComputeUngappedKarlinParameters(int, int, double, KarlinKMethod)"/>.
    /// </summary>
    /// <param name="match">Match score (&gt; 0).</param>
    /// <param name="mismatch">Mismatch score (&lt; 0).</param>
    /// <param name="baseFrequencies">Frequencies of A, C, G, T (four finite non-negative values with a positive sum;
    /// normalized to sum 1 as BLAST+ <c>BlastScoreFreqCalc</c> normalizes the score probabilities).</param>
    /// <param name="kMethod">K for schemes whose scores share a divisor &gt; 1 (default
    /// <see cref="KarlinKMethod.ReducedLattice"/>; <see cref="KarlinKMethod.NcbiBlast"/> = BLAST+ indexing).</param>
    /// <exception cref="ArgumentNullException"><paramref name="baseFrequencies"/> is null.</exception>
    /// <exception cref="ArgumentException">Not exactly four frequencies, a negative/non-finite one, or a zero sum.</exception>
    /// <exception cref="ArgumentOutOfRangeException">λ undefined for the scheme and composition.</exception>
    public static KarlinAltschulParameters ComputeUngappedKarlinParameters(
        int match, int mismatch, IReadOnlyList<double> baseFrequencies,
        KarlinKMethod kMethod = KarlinKMethod.ReducedLattice)
    {
        ArgumentNullException.ThrowIfNull(baseFrequencies);
        if (baseFrequencies.Count != 4)
            throw new ArgumentException("Exactly four base frequencies (A, C, G, T) are required.", nameof(baseFrequencies));
        double sum = 0, sumSquares = 0;
        foreach (double f in baseFrequencies)
        {
            if (!double.IsFinite(f) || f < 0)
                throw new ArgumentException("Base frequencies must be finite and non-negative.", nameof(baseFrequencies));
            sum += f;
            sumSquares += f * f;
        }

        if (!(sum > 0))
            throw new ArgumentException("Base frequencies must have a positive sum.", nameof(baseFrequencies));

        ValidateMatchMismatch(match, mismatch);
        double pMatch = sumSquares / (sum * sum);
        if (!(pMatch < 1.0))
            throw new ArgumentOutOfRangeException(nameof(baseFrequencies),
                "Karlin–Altschul λ is undefined: a single-base composition has no mismatch (expected score not negative).");
        return UngappedKarlinParameters(match, mismatch, pMatch, kMethod);
    }

    /// <summary>
    /// Gapped Karlin–Altschul parameters NCBI BLAST+ uses for a blastn reward/penalty/gap-cost combination
    /// (<c>Blast_KarlinBlkNuclGappedCalc</c> + <c>Blast_GetNuclAlphaBeta</c>): the simulated values of the
    /// <c>blastn_values_*</c> tables (gap costs scaled and λ, α divided by gcd(reward, penalty)); gap costs at
    /// or beyond the table's "infinite" domain use the ungapped parameters; gap costs 0/0 select BLAST's
    /// non-affine (greedy megablast) row.
    /// </summary>
    /// <remarks>
    /// Examples (blast_stat.c; printed by blastn 2.12.0+): +2/−3 gap 5/2 (blastn task default) → λ 0.625, K 0.41,
    /// H 0.78, α 0.8, β −2, even-score round-down; +1/−3 gap 2/2 → 1.37 / 0.70 / 1.2; +1/−2 gap 2/2 → 1.33 / 0.62 / 1.1.
    /// </remarks>
    /// <param name="reward">Match reward (&gt; 0).</param>
    /// <param name="penalty">Mismatch penalty (&lt; 0).</param>
    /// <param name="gapOpen">Gap existence cost (≥ 0; a gap of length k costs gapOpen + k·gapExtend).</param>
    /// <param name="gapExtend">Gap extension cost (≥ 0).</param>
    /// <param name="kMethod">
    /// K of the ungapped block copied in the infinite gap-cost domain for a scheme whose scores share a divisor &gt; 1:
    /// <see cref="KarlinKMethod.NcbiBlast"/> (default — what BLAST+ copies; blastn +4/−6 gap 12/8 prints "Gapped …
    /// 0.317 1.17 0.912") or <see cref="KarlinKMethod.ReducedLattice"/> (scale-invariant). The tabulated rows are unaffected.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException">Non-positive reward, non-negative penalty, negative gap costs or
    /// an undefined <paramref name="kMethod"/>.</exception>
    /// <exception cref="ArgumentException">A reward/penalty or gap-cost combination BLAST+ does not support.</exception>
    public static KarlinAltschulParameters GetBlastnGappedKarlinParameters(
        int reward, int penalty, int gapOpen, int gapExtend, KarlinKMethod kMethod = KarlinKMethod.NcbiBlast)
    {
        ValidateMatchMismatch(reward, penalty);
        if (gapOpen < 0)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap existence cost cannot be negative.");
        if (gapExtend < 0)
            throw new ArgumentOutOfRangeException(nameof(gapExtend), "Gap extension cost cannot be negative.");

        int divisor = Gcd(reward, -penalty);
        if (!BlastnValueTables.TryGetValue((reward / divisor, penalty / divisor), out var table))
            throw new ArgumentException(
                $"Substitution scores {reward} and {penalty} are not supported by NCBI BLAST+ blastn statistics.",
                nameof(reward));

        // s_SplitArrayOf8: a leading {0, 0} row is the non-affine entry; the remaining rows are the affine ones.
        bool split = table.Rows[0][0] == 0 && table.Rows[0][1] == 0;
        double[]? linear = split ? table.Rows[0] : null;
        IEnumerable<double[]> normal = split ? table.Rows.Skip(1) : table.Rows;

        double[]? row = null;
        if (gapOpen == 0 && gapExtend == 0 && linear is not null)
        {
            row = linear;
        }
        else
        {
            // s_AdjustGapParametersByGcd: the table's gap costs are multiplied by the divisor.
            row = normal.FirstOrDefault(r => (int)r[0] * divisor == gapOpen && (int)r[1] * divisor == gapExtend);
        }

        if (row is not null)
        {
            // s_AdjustGapParametersByGcd: λ and α are divided by the divisor.
            return new KarlinAltschulParameters(
                row[2] / divisor, row[3], row[4], row[5] / divisor, row[6], table.RoundDown, Gapped: true);
        }

        if (gapOpen >= table.GapOpenMax * divisor && gapExtend >= table.GapExtendMax * divisor)
        {
            // Infinite gap-cost domain: Blast_KarlinBlkCopy(kbp, kbp_ungap); α/β fall back to the ungapped values.
            var ungapped = UngappedKarlinParameters(reward, penalty, MatchProbability(UniformBaseFrequency), kMethod);
            return ungapped with { RoundDown = table.RoundDown, Gapped = true };
        }

        throw new ArgumentException(
            $"Gap existence and extension values {gapOpen} and {gapExtend} are not supported for substitution scores " +
            $"{reward} and {penalty}; supported: " +
            string.Join(", ", normal.Select(r => $"{(int)r[0] * divisor}/{(int)r[1] * divisor}")) +
            $", or any values at least {table.GapOpenMax * divisor}/{table.GapExtendMax * divisor}.",
            nameof(gapOpen));
    }

    /// <summary>
    /// BLAST+ edge-effect length adjustment ℓ (<c>BLAST_ComputeLengthAdjustment</c>): the integer approximation to the
    /// fixed point of ℓ = β + (α/λ)·(ln K + ln((m − ℓ)(n − N·ℓ))), kept small enough that K(m − ℓ)(n − N·ℓ) &gt; max(m, n).
    /// The effective search space is then (m − ℓ)(n − N·ℓ) (Altschul &amp; Gish 1996; Altschul et al. 2001).
    /// </summary>
    /// <param name="k">K (&gt; 0).</param>
    /// <param name="alphaOverLambda">α/λ (ungapped: 1/H).</param>
    /// <param name="beta">β.</param>
    /// <param name="queryLength">Query length m (&gt; 0).</param>
    /// <param name="databaseLength">Database length n (&gt; 0).</param>
    /// <param name="databaseSequenceCount">Number of database sequences N (≥ 1).</param>
    /// <returns>The length adjustment ℓ ≥ 0.</returns>
    /// <exception cref="ArgumentOutOfRangeException">A non-positive K, length or sequence count, or a non-finite α/λ or β.</exception>
    public static int ComputeLengthAdjustment(
        double k, double alphaOverLambda, double beta, int queryLength, long databaseLength, int databaseSequenceCount = 1)
    {
        if (!(k > 0) || !double.IsFinite(k))
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive and finite.");
        if (!double.IsFinite(alphaOverLambda))
            throw new ArgumentOutOfRangeException(nameof(alphaOverLambda), "α/λ must be finite.");
        if (!double.IsFinite(beta))
            throw new ArgumentOutOfRangeException(nameof(beta), "β must be finite.");
        if (queryLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(queryLength), "Query length m must be positive.");
        if (databaseLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseLength), "Database length n must be positive.");
        if (databaseSequenceCount <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseSequenceCount), "Database sequence count N must be positive.");

        double m = queryLength;
        double n = databaseLength;
        double bigN = databaseSequenceCount;
        double logK = Math.Log(k);

        // ell_max: largest ℓ with K(m − ℓ)(n − Nℓ) > max(m, n) (quadratic formula 2c/(−b + √(b² − 4ac))).
        double a = bigN;
        double mb = m * bigN + n;
        double c = n * m - Math.Max(m, n) / k;
        if (c < 0)
            return 0;
        double ellMax = 2 * c / (mb + Math.Sqrt(mb * mb - 4 * a * c));

        double ellMin = 0, ellNext = 0;
        bool converged = false;
        for (int i = 1; i <= LengthAdjustmentMaxIterations; i++)
        {
            double ell = ellNext;
            double ss = (m - ell) * (n - bigN * ell);
            double ellBar = alphaOverLambda * (logK + Math.Log(ss)) + beta;
            if (ellBar >= ell)
            {
                ellMin = ell;
                if (ellBar - ellMin <= 1.0)
                {
                    converged = true;
                    break;
                }

                if (ellMin == ellMax)
                    break;
            }
            else
            {
                ellMax = ell;
            }

            if (ellMin <= ellBar && ellBar <= ellMax)
                ellNext = ellBar;          // ell_bar is in range: accept it
            else if (i == 1)
                ellNext = ellMax;
            else
                ellNext = (ellMin + ellMax) / 2;
        }

        int adjustment = (int)ellMin;
        if (converged)
        {
            // floor(ell_min) is taken as floor(ell_fixed) unless ceil(ell_min) is still below the fixed point.
            double ell = Math.Ceiling(ellMin);
            if (ell <= ellMax)
            {
                double ss = (m - ell) * (n - bigN * ell);
                if (alphaOverLambda * (logK + Math.Log(ss)) + beta >= ell)
                    adjustment = (int)ell;
            }
        }

        return adjustment;
    }

    /// <summary>
    /// NCBI BLAST+ blastn statistics of an alignment score: the Karlin–Altschul parameters of the scoring scheme
    /// (gapped: <see cref="GetBlastnGappedKarlinParameters"/>; ungapped: <see cref="ComputeUngappedKarlinParameters(int, int, double, KarlinKMethod)"/>),
    /// the edge-effect length adjustment (<see cref="ComputeLengthAdjustment"/>), the effective search space
    /// (m − ℓ)·max(1, n − N·ℓ) (<c>BLAST_CalcEffLengths</c>), the bit score S' = (λS − ln K)/ln 2 and the E-value
    /// E = K·(m − ℓ)(n − N·ℓ)·e^{−λS}, with odd gapped scores rounded down to even for the E-value where the BLAST+
    /// table requires it (<c>Blast_HSPListGetEvalues</c>).
    /// </summary>
    /// <remarks>
    /// Reproduces NCBI blastn 2.12.0+: e.g. a 40-nt query against a 3079-nt subject with +2/−3, gap 5/2 →
    /// ℓ = 11, effective search space 88972, raw score 80 → 73.4 bits, E = 7e-18; a raw score of 15 is evaluated as 14.
    /// The scoring follows <see cref="ScoringMatrix"/>'s convention (gap of length k = GapOpen + k·GapExtend),
    /// so BLAST's gap existence/extension costs are −GapOpen/−GapExtend.
    /// </remarks>
    /// <param name="rawScore">Raw alignment score S (e.g. <see cref="CrossHybridizationAssessment.AlignmentScore"/>).</param>
    /// <param name="queryLength">Query (probe) length m (&gt; 0).</param>
    /// <param name="databaseLength">Total database length n (&gt; 0).</param>
    /// <param name="databaseSequenceCount">Number of database sequences N (≥ 1; default 1, a single subject).</param>
    /// <param name="scoring">Scoring scheme (default <see cref="SequenceAligner.BlastDna"/>: +2/−3, gap 5/2).</param>
    /// <param name="gapped">Gapped statistics (default true, as blastn); false for ungapped HSP scores.</param>
    /// <param name="kMethod">
    /// Ungapped K for a scheme whose scores share a divisor &gt; 1 (ungapped statistics, or gapped in the infinite
    /// gap-cost domain): <see cref="KarlinKMethod.NcbiBlast"/> (default — blastn's own value; +4/−6 ungapped prints
    /// K 1.17, effective search space 573996 for a 200-nt query vs a 3100-nt subject) or
    /// <see cref="KarlinKMethod.ReducedLattice"/> (scale-invariant K, = the reduced scheme's). No effect for δ = 1.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException">Non-positive lengths/count, or positive gap scores.</exception>
    /// <exception cref="ArgumentException">A scoring scheme BLAST+ has no gapped statistics for.</exception>
    public static BlastnStatistics ComputeBlastnStatistics(
        int rawScore,
        int queryLength,
        long databaseLength,
        int databaseSequenceCount = 1,
        ScoringMatrix? scoring = null,
        bool gapped = true,
        KarlinKMethod kMethod = KarlinKMethod.NcbiBlast)
    {
        var matrix = scoring ?? SequenceAligner.BlastDna;
        if (queryLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(queryLength), "Query length m must be positive.");
        if (databaseLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseLength), "Database length n must be positive.");
        if (databaseSequenceCount <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseSequenceCount), "Database sequence count N must be positive.");
        if (matrix.GapOpen > 0 || matrix.GapExtend > 0)
            throw new ArgumentOutOfRangeException(nameof(scoring), "Gap scores must be non-positive (penalties).");

        KarlinAltschulParameters p = gapped
            ? GetBlastnGappedKarlinParameters(matrix.Match, matrix.Mismatch, -matrix.GapOpen, -matrix.GapExtend, kMethod)
            : ComputeUngappedKarlinParameters(matrix.Match, matrix.Mismatch, UniformBaseFrequency, kMethod);

        int adjustment = ComputeLengthAdjustment(
            p.K, p.Alpha / p.Lambda, p.Beta, queryLength, databaseLength, databaseSequenceCount);

        // BLAST_CalcEffLengths: effective_db_length = n − N·ℓ (at least 1); search space × (m − ℓ).
        double effectiveDb = Math.Max(1.0, databaseLength - (double)databaseSequenceCount * adjustment);
        double searchSpace = effectiveDb * (queryLength - adjustment);

        // Blast_HSPListGetEvalues: score &= ~1 for gapped round-down tables (E-value only).
        int eScore = p.Gapped && p.RoundDown ? rawScore & ~1 : rawScore;
        double logK = Math.Log(p.K);
        double eValue = searchSpace * Math.Exp(-p.Lambda * eScore + logK);
        double bitScore = (rawScore * p.Lambda - logK) / Math.Log(2.0);

        return new BlastnStatistics(rawScore, eScore, p, bitScore, eValue, adjustment, searchSpace,
            queryLength, databaseLength, databaseSequenceCount);
    }

    /// <summary>
    /// BLAST+ blastn statistics of the best local alignment of a probe with one subject sequence (bl2seq-style
    /// search space: m = probe length, n = subject length, N = 1). The alignment is the canonical affine
    /// Smith–Waterman–Gotoh <see cref="SequenceAligner.LocalAlignAffine(string, string, ScoringMatrix?)"/> — the same
    /// optimal local score <see cref="AssessCrossHybridization"/> reports — and its score is evaluated by
    /// <see cref="ComputeBlastnStatistics(int, int, long, int, ScoringMatrix?, bool, KarlinKMethod)"/> with gapped parameters.
    /// Only the given strand is aligned; pass the reverse complement for the other strand.
    /// </summary>
    /// <param name="probeSequence">Probe (query) sequence (non-empty).</param>
    /// <param name="subjectSequence">Subject (off-target) sequence (non-empty).</param>
    /// <param name="scoring">Scoring scheme (default <see cref="SequenceAligner.BlastDna"/>).</param>
    /// <exception cref="ArgumentNullException">A null sequence.</exception>
    /// <exception cref="ArgumentException">An empty sequence or an unsupported scoring scheme.</exception>
    public static BlastnStatistics ComputeBlastnStatistics(
        string probeSequence, string subjectSequence, ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(subjectSequence);
        if (probeSequence.Length == 0)
            throw new ArgumentException("Probe sequence cannot be empty.", nameof(probeSequence));
        if (subjectSequence.Length == 0)
            throw new ArgumentException("Subject sequence cannot be empty.", nameof(subjectSequence));

        var matrix = scoring ?? SequenceAligner.BlastDna;
        var (score, _, _, _) = BestLocalAlignment(
            probeSequence.ToUpperInvariant(), subjectSequence.ToUpperInvariant(), matrix);
        return ComputeBlastnStatistics(score, probeSequence.Length, subjectSequence.Length, 1, matrix, gapped: true);
    }

    private static void ValidateMatchMismatch(int match, int mismatch)
    {
        if (match <= 0)
            throw new ArgumentOutOfRangeException(nameof(match),
                "Karlin–Altschul λ is undefined: the scoring scheme must have at least one positive score.");
        if (mismatch >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the mismatch score must be negative.");
    }

    // p(match) = 4·p² for four equiprobable bases of frequency p (p must keep it a probability in (0, 1)).
    private static double MatchProbability(double baseFrequency)
    {
        if (double.IsNaN(baseFrequency) || baseFrequency <= 0 || baseFrequency >= 0.5)
            throw new ArgumentOutOfRangeException(nameof(baseFrequency), "Base frequency must lie in (0, 0.5).");
        return 4.0 * baseFrequency * baseFrequency;
    }

    // Blast_KarlinBlkUngappedCalc for the two-score distribution {match: pMatch, mismatch: 1 − pMatch}.
    private static KarlinAltschulParameters UngappedKarlinParameters(
        int match, int mismatch, double pMatch, KarlinKMethod kMethod = KarlinKMethod.ReducedLattice)
    {
        if (kMethod is not (KarlinKMethod.ReducedLattice or KarlinKMethod.NcbiBlast))
            throw new ArgumentOutOfRangeException(nameof(kMethod), kMethod, "Unknown K method.");

        double pMismatch = 1.0 - pMatch;
        double expected = pMatch * match + pMismatch * mismatch;
        if (expected >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the expected per-pair score must be negative.");

        double lambda = SolveLambda(pMatch, match, pMismatch, mismatch);

        // BlastKarlinLtoH: H = λ Σ s·p_s·e^{λs}.
        double h = lambda * (pMatch * match * Math.Exp(lambda * match) + pMismatch * mismatch * Math.Exp(lambda * mismatch));

        // BlastKarlinLHtoK on the lattice reduced by δ = gcd(match, −mismatch): scores low = mismatch/δ … high = match/δ,
        // λ·δ, mean score / δ (H is scale-invariant). The closed forms read the true probabilities of the lowest and
        // highest scores; the series reads seriesProb[j], j = 0 … high − low:
        //  - ReducedLattice: the reduced lattice (seriesProb[j] = probability of reduced score low + j);
        //  - NcbiBlast: BLAST+'s probArrayStartLow, the UNREDUCED array (seriesProb[j] = probability of the unreduced
        //    score mismatch + j), which BlastKarlinLHtoK indexes with the reduced offsets.
        int delta = Gcd(match, -mismatch);
        int low = mismatch / delta, high = match / delta;
        int seriesSpan = kMethod == KarlinKMethod.NcbiBlast ? match - mismatch : high - low;
        var seriesProb = new double[seriesSpan + 1];
        seriesProb[0] = pMismatch;
        seriesProb[seriesSpan] = pMatch;
        double k = KarlinLHtoK(seriesProb, pMismatch, pMatch, low, high, lambda * delta, h, expected / delta);

        double beta = (match == 1 && mismatch == -1) || (match == 2 && mismatch == -3) ? -2.0 : 0.0;
        return new KarlinAltschulParameters(lambda, k, h, lambda / h, beta, RoundDown: false, Gapped: false);
    }

    // Unique positive root of pM·e^{λ·match} + pMm·e^{λ·mismatch} = 1 by bisection to double resolution.
    private static double SolveLambda(double pMatch, int match, double pMismatch, int mismatch)
    {
        double lo = 0.0;
        double hi = LambdaSearchUpperBound;
        for (int i = 0; i < LambdaBisectionIterations; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (pMatch * Math.Exp(mid * match) + pMismatch * Math.Exp(mid * mismatch) - 1.0 > 0.0)
                hi = mid;
            else
                lo = mid;
        }

        return 0.5 * (lo + hi);
    }

    // BlastKarlinLHtoK (blast_stat.c) after the gcd reduction: low/high/λ/mean are the reduced values, pLow/pHigh the
    // probabilities of the lowest/highest score (closed forms, sprob[low·δ]/sprob[high·δ]) and prob[j], j = 0 … high − low,
    // the score probabilities the series reads (probArrayStartLow[j]).
    private static double KarlinLHtoK(
        double[] prob, double pLow, double pHigh, int low, int high, double lambda, double h, double scoreAverage)
    {
        int range = high - low;
        double firstTermClosedForm = h / lambda;
        double expMinusLambda = Math.Exp(-lambda);

        if (low == -1 && high == 1)
            return (pLow - pHigh) * (pLow - pHigh) / pLow;

        if (low == -1 || high == 1)
        {
            if (high != 1)
                firstTermClosedForm = scoreAverage * scoreAverage / firstTermClosedForm;
            return firstTermClosedForm * (1.0 - expMinusLambda);
        }

        // P(i, j): probability of total score i over a gapless alignment of j pairs (shifted so index 0 ↔ j·low).
        var p = new double[KarlinKIterMax * range + 1];
        double outerSum = 0.0, innerSum = 1.0;
        int lowScore = 0, highScore = 0;
        p[0] = 1.0;
        for (int iter = 0; iter < KarlinKIterMax && innerSum > KarlinKSumLimit;)
        {
            int first = range, last = range;
            lowScore += low;
            highScore += high;
            for (int idx = highScore - lowScore; idx >= 0; idx--)
            {
                double sum = 0.0;
                for (int i1 = idx - first, j = first; i1 >= idx - last; i1--, j++)
                    sum += p[i1] * prob[j];
                p[idx] = sum;
                if (first > 0)
                    first--;
                if (idx <= range)
                    last--;
            }

            // Horner's rule: Σ_{i<0} P(i)·e^{λi} + Σ_{i≥0} P(i).
            int ptr = 0;
            innerSum = p[ptr];
            int score = lowScore + 1;
            for (; score < 0; score++)
                innerSum = p[++ptr] + innerSum * expMinusLambda;
            innerSum *= expMinusLambda;
            for (; score <= highScore; score++)
                innerSum += p[++ptr];

            innerSum /= ++iter;
            outerSum += innerSum;
        }

        return -Math.Exp(-2.0 * outerSum) / (firstTermClosedForm * BlastExpm1(-lambda));
    }

    // BLAST_Expm1 (ncbi_math.c).
    private static double BlastExpm1(double x)
    {
        double absx = Math.Abs(x);
        if (absx > .33)
            return Math.Exp(x) - 1.0;
        if (absx < 1e-16)
            return x;
        return x * (1.0 + x * (1.0 / 2 + x * (1.0 / 6 + x * (1.0 / 24 + x * (1.0 / 120 + x * (1.0 / 720
            + x * (1.0 / 5040 + x * (1.0 / 40320 + x * (1.0 / 362880 + x * (1.0 / 3628800 + x * (1.0 / 39916800
            + x * (1.0 / 479001600 + x / 6227020800.0))))))))))));
    }

    private static int Gcd(int a, int b)
    {
        a = Math.Abs(a);
        b = Math.Abs(b);
        while (b != 0)
            (a, b) = (b, a % b);
        return a;
    }

    #endregion

    #region Oligo Analysis

    /// <summary>
    /// Analyzes oligonucleotide properties: Tm (Primer3 <c>seqtm</c> at the Primer3 hybridization-probe
    /// conditions 50 nM / 50 mM monovalent / no Mg²⁺ / no dNTP — SantaLucia 1998 nearest-neighbour for
    /// ≤ 36 nt, <c>long_seq_tm</c> above; <c>NaN</c> for fewer than 2 bases or any non-ACGT base, e.g. RNA),
    /// G+C fraction, single-stranded molecular weight (<see cref="CalculateMolecularWeight(string)"/>) and the
    /// mononucleotide-sum ε260 (<see cref="CalculateExtinctionCoefficient(string)"/>).
    /// </summary>
    public static (double Tm, double GcContent, double MolecularWeight, double ExtinctionCoefficient)
        AnalyzeOligo(string sequence)
    {
        sequence = sequence.ToUpperInvariant();

        double tm = CalculateProbeTm(sequence, Primer3ProbeConditions);
        double gc = sequence.CalculateGcFractionFast();
        double mw = CalculateMolecularWeight(sequence);
        double extinction = CalculateExtinctionCoefficient(sequence);

        return (tm, gc, mw, extinction);
    }

    /// <summary>
    /// Average molecular weight (Da) of a single-stranded linear oligonucleotide, delegating to the
    /// canonical <c>SequenceStatistics.CalculateNucleotideMolecularWeight</c> (Biopython
    /// <c>Bio.SeqUtils.molecular_weight(seq, seq_type)</c>: sum of nucleoside-monophosphate masses −
    /// (n − 1) × water 18.01528; DNA dAMP 331.2218, dCMP 307.1971, dGMP 347.2212, dTMP 322.2085; RNA AMP 347.2212,
    /// CMP 323.1965, GMP 363.2206, UMP 324.1813; unknown symbols skipped). The strand type is inferred:
    /// RNA when the sequence contains U and no T, otherwise DNA. Use
    /// <see cref="CalculateMolecularWeight(string, bool)"/> to state it explicitly.
    /// </summary>
    public static double CalculateMolecularWeight(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;
        string s = sequence.ToUpperInvariant();
        bool isRna = s.Contains('U') && !s.Contains('T');
        return CalculateMolecularWeight(s, isDna: !isRna);
    }

    /// <summary>
    /// Average molecular weight (Da) of a single-stranded linear DNA (<paramref name="isDna"/> = true) or RNA
    /// oligonucleotide — the canonical <c>SequenceStatistics.CalculateNucleotideMolecularWeight</c>
    /// (= Biopython <c>molecular_weight(seq, "DNA" | "RNA")</c>).
    /// </summary>
    public static double CalculateMolecularWeight(string sequence, bool isDna) =>
        Seqeron.Genomics.Analysis.SequenceStatistics.CalculateNucleotideMolecularWeight(sequence, isDna);

    /// <summary>
    /// Molar extinction coefficient at 260 nm (M⁻¹·cm⁻¹) as the <b>sum of mononucleotide</b> values
    /// (A 15400, C 7400, G 11500, T 8700, U 9900; any other symbol 10000). This ignores base-stacking
    /// hypochromicity; the nearest-neighbour value is <see cref="CalculateExtinctionCoefficientNearestNeighbor"/>.
    /// </summary>
    public static double CalculateExtinctionCoefficient(string sequence)
    {
        double coefficient = 0;
        sequence = sequence.ToUpperInvariant();

        // Individual nucleotide contributions
        foreach (char c in sequence)
        {
            coefficient += c switch
            {
                'A' => 15400,
                'C' => 7400,
                'G' => 11500,
                'T' => 8700,
                'U' => 9900,
                _ => 10000
            };
        }

        return coefficient;
    }

    // Nearest-neighbour ε260 parameters (M⁻¹·cm⁻¹) of Cantor, Warshaw & Shapiro (1970, Biopolymers 9:1059)
    // for DNA and Warshaw & Tinoco (1966) for RNA, as tabulated by OligoCalc (Kibbe 2007) and oligo vendors.
    // Row = 5′ base, column = 3′ base, order A, C, G, T/U.
    private static readonly double[,] DnaNearestNeighborEpsilon =
    {
        { 27400, 21200, 25000, 22800 },
        { 21200, 14600, 18000, 15200 },
        { 25200, 17600, 21600, 20000 },
        { 23400, 16200, 19000, 16800 },
    };

    private static readonly double[] DnaMononucleotideEpsilon = { 15400, 7400, 11500, 8700 };

    private static readonly double[,] RnaNearestNeighborEpsilon =
    {
        { 27400, 21000, 25000, 24000 },
        { 21000, 14200, 17800, 16200 },
        { 25200, 17400, 21600, 21200 },
        { 24600, 17200, 20000, 19600 },
    };

    private static readonly double[] RnaMononucleotideEpsilon = { 15400, 7200, 11500, 9900 };

    /// <summary>
    /// Nearest-neighbour molar extinction coefficient at 260 nm of a single-stranded oligonucleotide
    /// (Cantor, Warshaw &amp; Shapiro 1970; Warshaw &amp; Tinoco 1966 for RNA):
    /// ε = Σ<sub>i=1..n−1</sub> ε(N<sub>i</sub>N<sub>i+1</sub>) − Σ<sub>i=2..n−1</sub> ε(N<sub>i</sub>); a single base is its
    /// mononucleotide value.
    /// </summary>
    /// <param name="sequence">Oligo (case-insensitive). DNA: A/C/G/T; RNA: A/C/G/U.</param>
    /// <param name="isDna">True for the DNA table, false for the RNA table.</param>
    /// <returns>ε260 in M⁻¹·cm⁻¹; 0 for an empty sequence; <c>NaN</c> when a base outside the table's alphabet occurs.</returns>
    public static double CalculateExtinctionCoefficientNearestNeighbor(string sequence, bool isDna = true)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        var nn = isDna ? DnaNearestNeighborEpsilon : RnaNearestNeighborEpsilon;
        var mono = isDna ? DnaMononucleotideEpsilon : RnaMononucleotideEpsilon;
        char fourth = isDna ? 'T' : 'U';

        int[] idx = new int[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            int k = c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => c == fourth ? 3 : -1 };
            if (k < 0)
                return double.NaN;
            idx[i] = k;
        }

        if (idx.Length == 1)
            return mono[idx[0]];

        double epsilon = 0;
        for (int i = 0; i < idx.Length - 1; i++)
            epsilon += nn[idx[i], idx[i + 1]];
        for (int i = 1; i < idx.Length - 1; i++)
            epsilon -= mono[idx[i]];
        return epsilon;
    }

    /// <summary>
    /// Calculates concentration from absorbance.
    /// </summary>
    public static double CalculateConcentration(
        double absorbance260,
        double extinctionCoefficient,
        double pathLength = 1.0)
    {
        // Beer-Lambert law: A = εcl
        // c = A / (ε * l)
        return absorbance260 / (extinctionCoefficient * pathLength) * 1e6; // µM
    }

    #endregion

    #region Helper Methods

    // Probe Tm = Primer3 seqtm at the parameters' hybridization conditions (NaN when not computable).
    private static double CalculateProbeTm(string sequence, ProbeParameters param) =>
        PrimerDesigner.CalculateMeltingTemperaturePrimer3(
            sequence,
            param.DnaConcentrationNanomolar,
            param.MonovalentMillimolar,
            param.DivalentMillimolar,
            param.DntpMillimolar,
            param.MaxNearestNeighborLength);

    // Primer3 hybridization-probe (internal-oligo) conditions 50 nM / 50 mM / 0 Mg²⁺ / 0 dNTP, nearest-neighbour Tm up to
    // 36 nt, with the remaining settings of the Microarray preset (thermodynamic screen, 47 °C): the default conditions of
    // ValidateProbe, AssessCrossHybridization, DesignMolecularBeacon and AnalyzeOligo (the Microarray preset's own
    // conditions are OligoArray's 1 M / 1 µM since A3-10).
    private static ProbeParameters Primer3ProbeConditions => Defaults.Microarray with
    {
        DnaConcentrationNanomolar = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar,
        MonovalentMillimolar = PrimerDesigner.Primer3InternalMonovalentMillimolar,
        DivalentMillimolar = PrimerDesigner.Primer3InternalDivalentMillimolar,
        DntpMillimolar = PrimerDesigner.Primer3InternalDntpMillimolar,
        MaxNearestNeighborLength = PrimerDesigner.Primer3MaxNnTmLength
    };

    private static double CalculateSelfComplementarity(string sequence)
    {
        string revComp = DnaSequence.GetReverseComplementString(sequence);
        int matches = 0;

        for (int i = 0; i < sequence.Length; i++)
        {
            if (sequence[i] == revComp[i])
                matches++;
        }

        return matches / (double)sequence.Length;
    }

    private static bool HasSecondaryStructurePotential(string sequence)
    {
        // Check for inverted repeats that could form hairpins
        int halfLen = sequence.Length / 2;

        for (int stemLen = 4; stemLen <= halfLen; stemLen++)
        {
            for (int i = 0; i <= sequence.Length - stemLen * 2 - 3; i++)
            {
                string left = sequence.Substring(i, stemLen);
                string right = sequence.Substring(i + stemLen + 3, stemLen);
                string rightRC = DnaSequence.GetReverseComplementString(right);

                int matches = 0;
                for (int j = 0; j < stemLen; j++)
                {
                    if (left[j] == rightRC[j])
                        matches++;
                }

                if (matches >= stemLen * 0.8)
                    return true;
            }
        }

        return false;
    }

    // Simple-sequence repeat: a di- or trinucleotide microsatellite of ≥ 4 complete copies, found by the
    // canonical MISA-convention scanner (primitive units only, A/C/G/T only).
    private static bool HasSimpleRepeats(string sequence) =>
        Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(sequence, 2, 3, 4).Any();

    #endregion
}
