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
    /// PRIMER_INTERNAL_DNTP_CONC = 0 (<c>libprimer3.cc</c> <c>pr_set_default_global_args</c>).
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> With <see cref="StructureScreen"/> = <see cref="ProbeStructureScreen.Thermodynamic"/>
    /// (default) a probe of ≤ 60 nt made of A/C/G/T is screened as Primer3 screens a hybridization probe
    /// (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1): ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin
    /// Tm (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>) must not exceed <see cref="MaxStructureTm"/>
    /// (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH = 47 °C). The self-dimer limit replaces
    /// <see cref="MaxSelfComplementarity"/>; the hairpin limit is applied when <see cref="AvoidSecondaryStructure"/>.
    /// Probes longer than 60 nt (thal.c <c>THAL_MAX_ALIGN</c>), probes with non-ACGT bases and
    /// <see cref="ProbeStructureScreen.Heuristic"/> use the sequence-only screens: the position-wise
    /// fold-back fraction (fraction of positions i with s[i] = revcomp(s)[i], compared with
    /// <see cref="MaxSelfComplementarity"/>) and an inverted-repeat stem (≥ 4 bp, loop 3, ≥ 80 % matched).
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

        /// <summary>Self-structure screen (default <see cref="ProbeStructureScreen.Thermodynamic"/>).</summary>
        public ProbeStructureScreen StructureScreen { get; init; } = ProbeStructureScreen.Thermodynamic;

        /// <summary>
        /// Maximum ntthal self-dimer / 3′ self-dimer / hairpin Tm (°C) of the thermodynamic screen
        /// (Primer3 PRIMER_INTERNAL_MAX_SELF_ANY_TH = _SELF_END_TH = _HAIRPIN_TH = 47 °C).
        /// </summary>
        public double MaxStructureTm { get; init; } = PrimerDesigner.Primer3MaxStructureTm;
    }

    /// <summary>Self-structure screen used by the probe designers (see <see cref="ProbeParameters"/>).</summary>
    public enum ProbeStructureScreen
    {
        /// <summary>Primer3 thermodynamic screen (ntthal self-dimer, 3′ self-dimer, hairpin Tm) for ≤ 60-nt ACGT probes.</summary>
        Thermodynamic,

        /// <summary>Sequence-only screens (fold-back fraction and inverted-repeat stem) for every probe.</summary>
        Heuristic
    }

    /// <summary>
    /// Default probe parameters for different applications.
    /// </summary>
    public static class Defaults
    {
        public static ProbeParameters Microarray => new(
            MinLength: 50, MaxLength: 60,
            MinTm: 75, MaxTm: 85,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3);

        public static ProbeParameters FISH => new(
            MinLength: 200, MaxLength: 500,
            MinTm: 70, MaxTm: 90,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 8,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);

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
        IReadOnlyList<string> Warnings);

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
    /// <param name="IsValid">True when no issue was recorded (no multiple hits, no self-structure flag, no
    /// cross-hybridizing non-target).</param>
    /// <param name="SpecificityScore">Library-defined uniqueness score 1/N over the N ungapped candidate binding
    /// sites (0 when there is none); not a published metric — see <see cref="ValidateProbe"/>.</param>
    /// <param name="OffTargetHits">Total ungapped k-mismatch hits across the references (the intended site
    /// included).</param>
    /// <param name="SelfComplementarity">Position-wise fold-back fraction (fraction of positions i with
    /// s[i] = revcomp(s)[i]); reported always, used as the self-complementarity criterion only by the
    /// sequence-only fallback screen.</param>
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
        /// <summary>True when the Primer3 thermodynamic self-structure screen was applied (≤ 60-nt A/C/G/T
        /// probe with <see cref="ProbeStructureScreen.Thermodynamic"/>); false for the sequence-only fallback.</summary>
        public bool ThermodynamicScreen { get; init; }

        /// <summary>ntthal self-dimer (ANY) Tm in °C (Primer3 SELF_ANY_TH); null when the fallback screen was used.</summary>
        public double? SelfDimerTm { get; init; }

        /// <summary>ntthal 3′ self-dimer (END1) Tm in °C (Primer3 SELF_END_TH); null when the fallback screen was used.</summary>
        public double? SelfEndDimerTm { get; init; }

        /// <summary>ntthal hairpin Tm in °C (Primer3 HAIRPIN_TH); null when the fallback screen was used.</summary>
        public double? HairpinTm { get; init; }

        /// <summary>Kane et al. (2000) cross-hybridization assessment of every supplied non-target sequence/strand
        /// (<see cref="AssessCrossHybridization"/>); empty when no non-target sequences were supplied.</summary>
        public IReadOnlyList<CrossHybridizationAssessment> CrossHybridization { get; init; }
            = Array.Empty<CrossHybridizationAssessment>();
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
        /// check of OligoArray (Rouillard et al. 2003). Null when the probe is longer than 60 nt (thal.c
        /// THAL_MAX_ALIGN), probe or site contain a non-ACGT base, or nothing aligns.
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
    /// <param name="K">The Karlin–Altschul search-space scale parameter K (supplied by the caller).</param>
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
        double gc = CalculateGcContent(seq);
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
    /// additive score is a library ranking heuristic; <see cref="DesignProbesPrimer3"/> is Primer3's published picker.
    /// </summary>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 10)
    {
        var param = parameters ?? Defaults.Microarray;

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
    /// Designs probes with genome-wide specificity check using suffix tree.
    /// O(n × m) for probe generation + O(m) per specificity check.
    /// </summary>
    /// <param name="targetSequence">Target sequence to design probes for.</param>
    /// <param name="genomeIndex">Pre-built suffix tree index for the genome (enables O(m) specificity lookup).</param>
    /// <param name="parameters">Probe design parameters.</param>
    /// <param name="maxProbes">Maximum number of probes to return.</param>
    /// <param name="requireUnique">If true, only return probes unique in the genome.</param>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters? parameters = null,
        int maxProbes = 10,
        bool requireUnique = true)
    {
        var param = parameters ?? Defaults.Microarray;

        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();

        // Get candidates using optimized method
        var candidates = DesignProbesOptimized(targetSequence, param, maxProbes * 5); // Get more candidates for filtering

        int returned = 0;
        foreach (var probe in candidates)
        {
            if (returned >= maxProbes)
                yield break;

            // Fast O(m) specificity check using suffix tree
            double specificity = CheckSpecificity(probe.Sequence, genomeIndex);

            if (requireUnique && specificity < 1.0)
                continue; // Skip non-unique probes

            // Boost score based on specificity
            var adjustedProbe = probe with
            {
                Score = probe.Score * specificity
            };

            returned++;
            yield return adjustedProbe;
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
        double BaseScore);

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

        return new ProbeBase(sequence, start, order, gc, tm, tmComputable, gcOutside, tmOutside,
            homopolymer, repeats, positionPenalty, score);
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
            warnings);
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

        double fraction = CalculateSelfComplementarity(sequence);
        bool high = fraction > param.MaxSelfComplementarity;
        bool stem = HasSecondaryStructurePotential(sequence);
        return (
            high,
            high ? $"High self-complementarity {fraction:P0}" : null,
            stem,
            stem ? "Potential secondary structure" : null);
    }

    // Primer3 thermodynamic self-structure of a probe (ntthal ANY / END1 self-dimer and hairpin Tm at the
    // parameters' conditions), or null when the sequence-only fallback applies (Heuristic screen, > 60 nt,
    // or a non-ACGT base).
    private static PrimerDesigner.Primer3OligoStructure? ComputeThermodynamicSelfStructure(
        string sequence, ProbeParameters param)
    {
        if (param.StructureScreen != ProbeStructureScreen.Thermodynamic || sequence.Length > NtthalMaxLength)
            return null;

        return PrimerDesigner.CalculatePrimer3OligoStructure(
            sequence, param.MonovalentMillimolar, param.DivalentMillimolar, param.DntpMillimolar,
            param.DnaConcentrationNanomolar);
    }

    /// <summary>
    /// Ranks every window of the configured length range (prefix-sum GC). The self-structure screens
    /// only lower a score, so they are evaluated lazily in descending base-score order and the scan stops
    /// once no remaining candidate can enter the top <paramref name="maxProbes"/> (branch and bound); the
    /// result equals an exhaustive evaluation: probes with score &gt; 0, ordered by score descending, ties
    /// by (length, start) ascending.
    /// </summary>
    private static List<Probe> DesignProbesOptimized(
        string targetSequence,
        ProbeParameters param,
        int maxProbes)
    {
        int n = targetSequence.Length;
        if (maxProbes <= 0)
            return new List<Probe>();

        // gcPrefixSum[i] = count of G/C in sequence[0..i-1]
        int[] gcPrefixSum = new int[n + 1];
        for (int i = 0; i < n; i++)
        {
            char c = targetSequence[i];
            gcPrefixSum[i + 1] = gcPrefixSum[i] + (c == 'G' || c == 'C' ? 1 : 0);
        }

        var bases = new List<ProbeBase>();
        int order = 0;
        for (int length = param.MinLength; length <= param.MaxLength && length <= n; length++)
        {
            for (int start = 0; start <= n - length; start++)
            {
                int gcCount = gcPrefixSum[start + length] - gcPrefixSum[start];
                double gc = (double)gcCount / length;

                // Early rejection far outside the GC window.
                if (gc < param.MinGc - 0.1 || gc > param.MaxGc + 0.1)
                    continue;

                var b = EvaluateProbeBase(targetSequence.Substring(start, length), start, order++, param, gc);
                if (b.BaseScore > 0)
                    bases.Add(b);
            }
        }

        // Rank key: score descending, then enumeration order ascending (= stable sort of the eager scan).
        static int CompareKeys((double Score, int Order) a, (double Score, int Order) b)
        {
            int c = b.Score.CompareTo(a.Score);
            return c != 0 ? c : a.Order.CompareTo(b.Order);
        }

        var best = new SortedSet<(double Score, int Order)>(Comparer<(double Score, int Order)>.Create(CompareKeys));
        var probesByOrder = new Dictionary<int, Probe>();
        var structureCache = new Dictionary<string, (bool, string?, bool, string?)>(StringComparer.Ordinal);
        foreach (var b in bases.OrderByDescending(x => x.BaseScore))
        {
            // Structure screens never raise a score: once a candidate's best possible key ranks after the
            // current k-th best, so does every later candidate (base score non-increasing, ties by order).
            if (best.Count == maxProbes && CompareKeys((b.BaseScore, b.Order), best.Max) > 0)
                break;

            if (FinishProbe(b, param, structureCache) is not { } p)
                continue;

            best.Add((p.Score, b.Order));
            probesByOrder[b.Order] = p;
            if (best.Count > maxProbes)
            {
                var worst = best.Max;
                best.Remove(worst);
                probesByOrder.Remove(worst.Order);
            }
        }

        return best.Select(k => probesByOrder[k.Order]).ToList();
    }

    /// <summary>
    /// Evaluates a potential probe sequence (all screens, eager).
    /// </summary>
    private static Probe? EvaluateProbe(string sequence, int start, ProbeParameters param)
    {
        double gc = CalculateGcContent(sequence);
        return FinishProbe(EvaluateProbeBase(sequence, start, 0, param, gc), param);
    }

    /// <summary>
    /// Designs tiling probes to cover entire sequence.
    /// </summary>
    /// <remarks>
    /// Windows of <paramref name="probeLength"/> start every <c>probeLength − overlap</c> bases; each is scored
    /// like <see cref="DesignProbes(string, ProbeParameters?, int)"/> (Tm = Primer3 <c>seqtm</c> at the
    /// parameters' conditions, thermodynamic self-structure screen for ≤ 60 nt). A window whose score is
    /// ≤ 0 is still emitted (score 0.3, warning "Suboptimal probe…") so coverage is preserved.
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
                double gc = CalculateGcContent(probeSeq);
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
        double DnaConcentrationNanomolar = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar);

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
        double SelfAnyTh,
        double SelfEndTh,
        double HairpinTh,
        double Penalty);

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
    /// ascending. Verified against primer3-py 2.3.1 <c>design_primers</c> (PRIMER_INTERNAL_n_* values).
    /// </summary>
    /// <param name="template">Template sequence (case-insensitive); probes are picked on this strand.</param>
    /// <param name="settings">Picker settings (default: Primer3 defaults).</param>
    /// <param name="numReturn">PRIMER_NUM_RETURN (default 5).</param>
    /// <returns>Up to <paramref name="numReturn"/> probes, best first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid sizes (MinSize &lt; 1, MaxSize &lt; MinSize,
    /// MaxSize &gt; 36) or a negative <paramref name="numReturn"/>.</exception>
    public static IReadOnlyList<Primer3Probe> DesignProbesPrimer3(
        string template,
        Primer3ProbeSettings? settings = null,
        int numReturn = 5)
    {
        ArgumentNullException.ThrowIfNull(template);
        var s = settings ?? new Primer3ProbeSettings();
        if (s.MinSize < 1 || s.MaxSize < s.MinSize || s.MaxSize > Primer3MaxOligoLength)
            throw new ArgumentOutOfRangeException(nameof(settings),
                $"Probe sizes must satisfy 1 ≤ MinSize ≤ MaxSize ≤ {Primer3MaxOligoLength}.");
        if (numReturn < 0)
            throw new ArgumentOutOfRangeException(nameof(numReturn), "Must be ≥ 0.");

        string seq = template.ToUpperInvariant();
        int n = seq.Length;
        var weights = PrimerDesigner.DefaultPrimer3Weights; // = Primer3 o_args.weights defaults
        var optima = new Primer3Optima(s.OptTm, s.OptSize, PrimerDesigner.DefaultPrimer3Optima.OptGcPercent);
        var accepted = new List<Primer3Probe>();

        // pick_primer_range: for every 3' end, oligos of increasing length (5' extensions). A failure that no
        // 5' extension can cure — too many Ns, poly-X, self-any (Primer3 five_prime_problem bits) — ends the
        // extension loop for that 3' end; checks run in calc_and_check_oligo_features order and stop at the
        // first failure (Ns, GC, poly-X, Tm, self-any, self-end, hairpin).
        for (int end = n - 1; end >= s.MinSize - 1; end--)
        {
            for (int len = s.MinSize; len <= s.MaxSize; len++)
            {
                int start = end - len + 1;
                if (start < 0)
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
                if (PrimerDesigner.FindLongestHomopolymer(oligo) > s.MaxPolyX)
                    break; // OP_HIGH_POLY_X: five-prime problem
                if (other)
                    continue; // non-ACGT, non-N symbol: no Tm (seqtm error → low Tm)

                double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(
                    oligo, s.DnaConcentrationNanomolar, s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar);
                if (tm < s.MinTm || tm > s.MaxTm)
                    continue;

                var st = PrimerDesigner.CalculatePrimer3OligoStructure(
                    oligo, s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar, s.DnaConcentrationNanomolar)!.Value;
                if (st.SelfAnyTh > s.MaxSelfAnyTh)
                    break; // OP_HIGH_SELF_ANY: five-prime problem
                if (st.SelfEndTh > s.MaxSelfEndTh || st.HairpinTh > s.MaxHairpinTh)
                    continue;

                double penalty = PrimerDesigner.CalculatePrimer3Penalty(
                    new Primer3PenaltyInputs(tm, len, gcPercent), weights, optima);
                accepted.Add(new Primer3Probe(oligo, start, len, tm, gcPercent,
                    st.SelfAnyTh, st.SelfEndTh, st.HairpinTh, penalty));
            }
        }

        // primer_rec_comp: quality ascending, then start descending, then length ascending.
        accepted.Sort((a, b) =>
        {
            int c = a.Penalty.CompareTo(b.Penalty);
            if (c != 0) return c;
            c = b.Start.CompareTo(a.Start);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        });

        return accepted.Count > numReturn ? accepted.GetRange(0, numReturn) : accepted;
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
    /// Warnings carry the stem/loop sizes and, for a beacon of ≤ 60 nt, the ntthal hairpin Tm of the whole
    /// beacon (the stem melting temperature; <see cref="PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double, double, double)"/>
    /// at the same conditions). <see cref="Probe.Tm"/> is the loop (probe–target) Tm.
    /// </para>
    /// </remarks>
    /// <param name="targetSequence">Target sequence (the loop is a window of it).</param>
    /// <param name="probeLength">Loop length (Tyagi &amp; Kramer: 15–30 nt; default 25).</param>
    /// <param name="stemLength">Arm length in bp (Tyagi &amp; Kramer: 5–7 bp; default 5).</param>
    /// <param name="detectionTemperatureCelsius">Optional detection (annealing) temperature T in °C.</param>
    public static Probe? DesignMolecularBeacon(
        string targetSequence,
        int probeLength = 25,
        int stemLength = 5,
        double? detectionTemperatureCelsius = null)
    {
        if (targetSequence.Length < probeLength)
            return null;

        targetSequence = targetSequence.ToUpperInvariant();
        var conditions = Defaults.Microarray; // only the (default Primer3 internal-oligo) conditions are used

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
            double gc = CalculateGcContent(loop);
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

        if (beaconSequence.Length <= NtthalMaxLength)
        {
            var hp = PrimerDesigner.CalculateHairpinThermodynamicsNtthal(
                beaconSequence,
                conditions.MonovalentMillimolar / 1000.0,
                conditions.DivalentMillimolar / 1000.0,
                conditions.DntpMillimolar / 1000.0);
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
            CalculateGcContent(beaconSequence),
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
    /// <see cref="ApproximateMatcher.FindWithMismatches(string, string, int)"/> (case-insensitive, overlapping, the
    /// strand given only); <see cref="ProbeValidation.OffTargetHits"/> is the total and includes the intended site.
    /// More than one hit records an issue. <see cref="ProbeValidation.SpecificityScore"/> = 1/N for N ≥ 1 hits
    /// (0 for none) is a library-defined uniqueness score (the share of the probe's N candidate binding sites
    /// taken by one site), not a published specificity metric; the sourced cross-hybridization decision is the
    /// Kane assessment below.
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> The same screen as <see cref="DesignProbes(string, ProbeParameters?, int)"/>: with
    /// <see cref="ProbeStructureScreen.Thermodynamic"/> (default) a ≤ 60-nt A/C/G/T probe is screened as Primer3
    /// screens a hybridization probe — ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>, primer3-py <c>calc_homodimer</c> /
    /// <c>calc_end_stability</c> / <c>calc_hairpin</c> parity) at the <paramref name="conditions"/> salt/oligo
    /// concentrations must not exceed <see cref="ProbeParameters.MaxStructureTm"/> (PRIMER_INTERNAL_MAX_SELF_ANY_TH
    /// = _SELF_END_TH = _HAIRPIN_TH = 47 °C). Longer probes (thal.c THAL_MAX_ALIGN = 60), non-ACGT probes and
    /// <see cref="ProbeStructureScreen.Heuristic"/> use the sequence-only screens: fold-back fraction &gt;
    /// <paramref name="selfComplementarityThreshold"/> and the inverted-repeat stem.
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
    /// <param name="maxMismatches">Maximum mismatches of the ungapped site-counting scan (≥ 0; default 3, a screening
    /// tolerance kept for compatibility — the sourced hybridization cross-reactivity decision is the Kane assessment).</param>
    /// <param name="selfComplementarityThreshold">Fold-back-fraction limit of the sequence-only fallback screen
    /// (default 0.3, the Microarray preset's <see cref="ProbeParameters.MaxSelfComplementarity"/>).</param>
    /// <param name="conditions">Hybridization conditions and structure-screen settings (default
    /// <see cref="Defaults.Microarray"/>: Primer3 probe conditions 50 nM / 50 mM / 0 Mg²⁺ / 0 dNTP, thermodynamic
    /// screen, 47 °C); its <see cref="ProbeParameters.MaxSelfComplementarity"/> is replaced by
    /// <paramref name="selfComplementarityThreshold"/>.</param>
    /// <param name="nonTargetSequences">Optional known non-target sequences for the Kane assessment.</param>
    /// <param name="maxNonTargetIdentity">Kane identity threshold (default 0.75; flagged when strictly above).</param>
    /// <param name="maxContiguousMatch">Kane contiguous-identity threshold in nt (default 15; flagged when strictly longer).</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style off-target duplex-Tm threshold (°C; see
    /// <see cref="AssessCrossHybridization"/>); null = Kane criteria only.</param>
    public static ProbeValidation ValidateProbe(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        int maxMismatches = 3,
        double selfComplementarityThreshold = 0.3,
        ProbeParameters? conditions = null,
        IEnumerable<string>? nonTargetSequences = null,
        double maxNonTargetIdentity = KaneMaxIdentity,
        int maxContiguousMatch = KaneMaxContiguousMatch,
        double? maxDuplexTm = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);

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

        // Check off-target hits via approximate matching
        foreach (var reference in referenceSequences)
        {
            // Canonical ungapped k-mismatch (Hamming) scan; case-insensitive, overlapping hits included.
            offTargetHits += ApproximateMatcher.FindWithMismatches(reference, probeSequence, maxMismatches).Count();
        }

        if (offTargetHits > 1)
        {
            issues.Add($"{offTargetHits} potential off-target sites");
        }

        // Self-structure: the DesignProbes screen (Primer3 ntthal for ≤ 60-nt ACGT probes, sequence-only otherwise).
        var param = (conditions ?? Defaults.Microarray) with { MaxSelfComplementarity = selfComplementarityThreshold };
        double selfComp = CalculateSelfComplementarity(probeSequence);
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
            selfCompIssue = selfComp > selfComplementarityThreshold;
            hasStructure = HasSecondaryStructurePotential(probeSequence);
            if (selfCompIssue)
                issues.Add($"Self-complementarity: {selfComp:P0}");
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
            CrossHybridization = cross,
        };
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
    /// (<see cref="CrossHybridizationAssessment.DuplexTm"/>; default Primer3 probe conditions 50 nM / 50 mM / 0 / 0).</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style specificity threshold (°C): a strand whose site duplex Tm
    /// is strictly above it is flagged (<see cref="CrossHybridizationAssessment.ExceedsDuplexTmThreshold"/>); null
    /// (default) applies only the Kane criteria. The threshold is assay-specific (Rouillard et al. 2003: user-set).</param>
    /// <returns>One assessment per non-target strand: index order, forward strand before reverse complement.</returns>
    /// <exception cref="ArgumentNullException">A null probe or non-target collection.</exception>
    /// <exception cref="ArgumentException">An empty probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxIdentity"/> outside [0, 1] or a negative
    /// <paramref name="maxContiguousMatch"/>.</exception>
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

        var matrix = scoring ?? SequenceAligner.BlastDna;
        var cond = conditions ?? Defaults.Microarray;
        string probe = probeSequence.ToUpperInvariant();
        bool duplexComputable = probe.Length <= NtthalMaxLength && probe.All(c => c is 'A' or 'C' or 'G' or 'T');
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
                if (site.All(c => c is 'A' or 'C' or 'G' or 'T'))
                {
                    // The probe hybridizes to the strand complementary to the site (primer3-py calc_heterodimer).
                    var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal(
                        probe, DnaSequence.GetReverseComplementString(site), PrimerDesigner.NtthalAlignmentMode.Any,
                        cond.MonovalentMillimolar / 1000.0, cond.DivalentMillimolar / 1000.0, cond.DntpMillimolar / 1000.0,
                        cond.DnaConcentrationNanomolar * 1e-9);
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
    /// Per-base background frequency (default 0.25, uniform). The four bases are assumed equiprobable
    /// at this value; the four base frequencies must sum to 1, i.e. <paramref name="baseFrequency"/> = 0.25.
    /// </param>
    /// <returns>The positive λ solving the Karlin–Altschul equation.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when the scheme cannot define λ: the match score is not positive, the mismatch score is
    /// not negative, or the expected per-pair score is not negative (the theory's preconditions).
    /// </exception>
    public static double ComputeLambdaNucleotide(
        int match,
        int mismatch,
        double baseFrequency = UniformBaseFrequency)
    {
        // Karlin–Altschul preconditions: at least one positive score, and negative expected score.
        if (match <= 0)
            throw new ArgumentOutOfRangeException(nameof(match),
                "Karlin–Altschul λ is undefined: the scoring scheme must have at least one positive score.");
        if (mismatch >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the mismatch score must be negative.");

        // p(match) = 4 · p² (the four identical ordered pairs); p(mismatch) = 1 − p(match).
        double pMatch = 4.0 * baseFrequency * baseFrequency;
        double pMismatch = 1.0 - pMatch;

        // Expected per-pair score must be negative for the theory to hold.
        double expectedScore = pMatch * match + pMismatch * mismatch;
        if (expectedScore >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the expected per-pair score must be negative.");

        // f(λ) = p(match)·e^{λ·match} + p(mismatch)·e^{λ·mismatch} − 1.
        // f(0) = 0; f'(0) = expectedScore < 0 so f dips below 0 for small λ>0, then a positive
        // match score drives e^{λ·match} → ∞, so f crosses 0 exactly once at the positive root.
        static double F(double lambda, double pM, int m, double pMm, int mm)
            => pM * Math.Exp(lambda * m) + pMm * Math.Exp(lambda * mm) - 1.0;

        double lo = 0.0;
        double hi = LambdaSearchUpperBound;
        for (int i = 0; i < LambdaBisectionIterations; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (F(mid, pMatch, match, pMismatch, mismatch) > 0.0)
                hi = mid;
            else
                lo = mid;
        }

        return 0.5 * (lo + hi);
    }

    /// <summary>
    /// Computes the Karlin–Altschul bit score and E-value for an off-target alignment hit's raw
    /// score, given the search-space dimensions and the scoring scheme.
    /// </summary>
    /// <remarks>
    /// <para>
    /// E = K·m·n·e^{−λS}, S' = (λS − ln K) / ln 2, E = m·n·2^{−S'} (Karlin &amp; Altschul 1990;
    /// Altschul et al. 1990). λ is computed from <paramref name="scoring"/> by
    /// <see cref="ComputeLambdaNucleotide"/>; K is supplied by the caller (the closed form requires
    /// the score-lattice machinery of Karlin–Altschul; for the BLAST +1/−3 nucleotide scheme the
    /// published value is K ≈ 0.711, NCBI blastn).
    /// </para>
    /// <para>This is additive and opt-in; <see cref="ScanOffTargetsGapped"/> and its defaults are unchanged.</para>
    /// </remarks>
    /// <param name="rawScore">The raw alignment score S of the hit.</param>
    /// <param name="queryLength">Query (probe) length m (&gt; 0).</param>
    /// <param name="databaseLength">Database (reference) length n (&gt; 0).</param>
    /// <param name="scoring">
    /// The scoring scheme. Its <see cref="ScoringMatrix.Match"/>/<see cref="ScoringMatrix.Mismatch"/>
    /// determine λ. Defaults to <see cref="SequenceAligner.BlastDna"/> (+2/−3). Pass a +1/−3 matrix to
    /// reproduce the published λ ≈ 1.374.
    /// </param>
    /// <param name="k">The Karlin–Altschul K parameter (default 0.711, the published nucleotide value).</param>
    /// <param name="baseFrequency">Per-base background frequency for λ (default 0.25, uniform).</param>
    /// <returns>The <see cref="KarlinAltschulStatistics"/> for the hit.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="scoring"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown for non-positive lengths or K, or a scheme for which λ is undefined.</exception>
    public static KarlinAltschulStatistics ComputeKarlinAltschul(
        double rawScore,
        int queryLength,
        long databaseLength,
        ScoringMatrix? scoring = null,
        double k = DefaultNucleotideK,
        double baseFrequency = UniformBaseFrequency)
    {
        var matrix = scoring ?? SequenceAligner.BlastDna;
        ArgumentNullException.ThrowIfNull(matrix);
        if (queryLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(queryLength), "Query length m must be positive.");
        if (databaseLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseLength), "Database length n must be positive.");
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        double lambda = ComputeLambdaNucleotide(matrix.Match, matrix.Mismatch, baseFrequency);

        // S' = (λS − ln K) / ln 2  (Altschul et al. 1990).
        double bitScore = (lambda * rawScore - Math.Log(k)) / Math.Log(2.0);

        // E = K·m·n·e^{−λS}  (Karlin & Altschul 1990).
        double eValue = k * queryLength * databaseLength * Math.Exp(-lambda * rawScore);

        return new KarlinAltschulStatistics(
            rawScore, lambda, k, bitScore, eValue, queryLength, databaseLength);
    }

    // Published Karlin–Altschul K for the BLAST +1/−3 nucleotide scheme (NCBI blastn reports
    // Lambda ≈ 1.37, K ≈ 0.711 for match=1/mismatch=−3). K's full closed form needs the
    // score-probability lattice/geometric-spacing machinery of Karlin & Altschul (1990); it is
    // therefore exposed as a caller parameter, defaulted to this published value.
    private const double DefaultNucleotideK = 0.711;

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

        double tm = CalculateProbeTm(sequence, Defaults.Microarray);
        double gc = CalculateGcContent(sequence);
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

    // thal.c THAL_MAX_ALIGN: ntthal cannot align two strands both longer than 60 nt (self-dimer, hairpin).
    private const int NtthalMaxLength = 60;

    private static double CalculateGcContent(string sequence) =>
        sequence.Length > 0 ? sequence.CalculateGcFractionFast() : 0;

    // Probe Tm = Primer3 seqtm at the parameters' hybridization conditions (NaN when not computable).
    private static double CalculateProbeTm(string sequence, ProbeParameters param) =>
        PrimerDesigner.CalculateMeltingTemperaturePrimer3(
            sequence,
            param.DnaConcentrationNanomolar,
            param.MonovalentMillimolar,
            param.DivalentMillimolar,
            param.DntpMillimolar);

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
