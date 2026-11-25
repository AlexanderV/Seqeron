using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs hybridization probes for various applications (FISH, microarray, Northern blot, etc.).
/// </summary>
public static partial class ProbeDesigner
{

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
        /// stem-loop screen <see cref="PrimerDesigner.HasHairpinPotential"/> (stem ≥ 4 bp, loop ≥ 3 nt).</summary>
        Heuristic
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
    /// Probe set for tiling (<see cref="DesignTilingProbes"/>).
    /// </summary>
    /// <param name="Probes">Tiling windows in start order (the last one may be anchored at the target end).</param>
    /// <param name="Coverage">Number of target positions covered by at least one window (= target length unless the
    /// overlap is negative).</param>
    /// <param name="MeanTm">Mean probe Tm.</param>
    /// <param name="TmRange">max(Tm) − min(Tm).</param>
    public readonly record struct TilingProbeSet(
        IReadOnlyList<Probe> Probes,
        int Coverage,
        double MeanTm,
        double TmRange);

    /// <summary>
    /// Probe validation result (<see cref="ValidateProbe"/>).
    /// </summary>
    /// <param name="IsValid">True when no issue was recorded (at most one reference site meeting the Kane et al. (2000)
    /// criteria, no sourced self-structure flag — ntthal self-dimer / 3′ self-dimer / hairpin Tm or, in the fallback,
    /// Primer3 alignment-mode self_any / self_end — and no cross-hybridizing non-target). The fallback sequence-only
    /// stem-loop flag is reported in <see cref="HasSecondaryStructure"/> / <see cref="Warnings"/> and does not enter it
    /// (audit round 5, A5-5).</param>
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
    /// screen; an issue, so <see cref="IsValid"/> is false) or the sequence-only stem-loop screen
    /// <see cref="PrimerDesigner.HasHairpinPotential"/> (fallback: an exactly complementary stem of ≥ 4 bp closing a loop
    /// of ≥ 3 nt; a library-convention presence test, recorded as a <see cref="Warnings"/> entry and not deciding
    /// <see cref="IsValid"/>).</param>
    /// <param name="Issues">Recorded validation issues (sourced criteria; any issue makes <see cref="IsValid"/> false).</param>
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
        /// alignment-mode self_any / self_end + the <see cref="PrimerDesigner.HasHairpinPotential"/> stem-loop screen).</summary>
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

        /// <summary>Reported findings that do not decide <see cref="IsValid"/>: the fallback sequence-only stem-loop flag
        /// ("Potential secondary structure formation", <see cref="PrimerDesigner.HasHairpinPotential"/> — stem ≥ 4 bp is a
        /// library convention without a published source, and almost every random probe of &gt; 60 nt has such a stem),
        /// set for probes screened by the fallback (longer than ThermodynamicScreenMaxLength, non-ACGT, or
        /// <see cref="ProbeStructureScreen.Heuristic"/>). Until audit round 5 (A5-5) this entry was an issue and made
        /// the probe invalid; raise <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> for a sourced hairpin
        /// decision (ntthal hairpin Tm) on long probes.</summary>
        public IReadOnlyList<string> Warnings { get; init; } = Array.Empty<string>();
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

    // Non-target strands longer than this are aligned in overlapping chunks (memory O(probe × chunk)).
    private const int CrossHybridizationChunkLength = 4096;

    // Best local alignment (score, identical columns) of the probe with a strand via the canonical affine
    // Smith–Waterman–Gotoh aligner, chunking long strands so that every positive-scoring alignment lies wholly
    // inside one chunk (its reference span is at most m + ⌈m·match/|extend|⌉ residues).
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static (int Score, int Identical, int Start, int End) BestLocalAlignment(
        string probe, string strand, ScoringMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(probe);
        ArgumentNullException.ThrowIfNull(strand);
        ArgumentNullException.ThrowIfNull(matrix);
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
            // Identical (non-gap) columns: the canonical SequenceAligner.CalculateStatistics Matches.
            int identical = SequenceAligner.CalculateStatistics(aln).Matches;

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

            // Column counts from the canonical SequenceAligner.CalculateStatistics: identical = Matches,
            // ungapped = Matches + Mismatches, gapped = Gaps > 0.
            var stats = SequenceAligner.CalculateStatistics(aln);
            int identical = stats.Matches;
            int ungapped = stats.Matches + stats.Mismatches;
            bool hasGap = stats.Gaps > 0;

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
    /// not negative, or the expected per-pair score is not negative (the theory's preconditions); a score outside
    /// BLAST+'s score range (<c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive): a match score ≥ 32767 (BLAST_SCORE_MAX) or a mismatch score ≤ −32768 (BLAST_SCORE_MIN) (ParamName <c>match</c> / <c>mismatch</c>); or when
    /// <paramref name="baseFrequency"/> is outside (0, 0.5) (p(match) = 4·p² would not be a probability in (0, 1)).
    /// </exception>
    public static double ComputeLambdaNucleotide(
        int match,
        int mismatch,
        double baseFrequency = UniformBaseFrequency)
        => LambdaNucleotide(match, mismatch, baseFrequency, nameof(match), nameof(mismatch));

    // ComputeLambdaNucleotide with the caller's parameter names for the scheme (A9-1: each public method reports its own).
    private static double LambdaNucleotide(
        int match, int mismatch, double baseFrequency, string matchParam, string mismatchParam)
    {
        // Karlin–Altschul preconditions: at least one positive score, and negative expected score.
        ValidateMatchMismatch(match, mismatch, matchParam, mismatchParam);

        // p(match) = 4 · p² (the four identical ordered pairs); p(mismatch) = 1 − p(match).
        double pMatch = MatchProbability(baseFrequency);
        double pMismatch = 1.0 - pMatch;

        // Expected per-pair score must be negative for the theory to hold.
        double expectedScore = pMatch * match + pMismatch * mismatch;
        if (expectedScore >= 0)
            throw new ArgumentOutOfRangeException(mismatchParam,
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
    /// determine λ. Null (the default) means <see cref="SequenceAligner.BlastDna"/> (+2/−3); it is not an error.
    /// Pass a +1/−3 matrix to reproduce the published λ ≈ 1.374.
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
    /// <exception cref="ArgumentOutOfRangeException">Thrown for non-positive lengths or K, a base frequency outside (0, 0.5),
    /// an undefined <paramref name="kMethod"/>, or (ParamName <c>scoring</c>) a scheme for which λ is undefined (non-positive
    /// match, non-negative mismatch or non-negative expected score) or with a score outside
    /// BLAST+'s score range (<c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive): a match score ≥ 32767 (BLAST_SCORE_MAX) or a mismatch score ≤ −32768 (BLAST_SCORE_MIN).</exception>
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

        double lambda = LambdaNucleotide(matrix.Match, matrix.Mismatch, baseFrequency, nameof(scoring), nameof(scoring));
        double kValue = k ?? UngappedKarlinParameters(
            matrix.Match, matrix.Mismatch, MatchProbability(baseFrequency), kMethod, nameof(scoring)).K;

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
    /// In ungapped <see cref="BlastnStatistics"/> of a reward/penalty scheme without a BLAST+ <c>blastn_values_*</c>
    /// table, α = β = 0 (BLAST+ then applies no length adjustment).
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
    /// non-negative expected score), a score at or beyond BLAST+'s BLAST_SCORE_MIN / BLAST_SCORE_MAX (INT2 limits, exclusive), a base
    /// frequency outside (0, 0.5), or an undefined <paramref name="kMethod"/>.</exception>
    public static KarlinAltschulParameters ComputeUngappedKarlinParameters(
        int match, int mismatch, double baseFrequency = UniformBaseFrequency,
        KarlinKMethod kMethod = KarlinKMethod.ReducedLattice)
    {
        ValidateMatchMismatch(match, mismatch, nameof(match), nameof(mismatch));
        return UngappedKarlinParameters(match, mismatch, MatchProbability(baseFrequency), kMethod, nameof(mismatch));
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
    /// <exception cref="ArgumentOutOfRangeException">λ undefined for the scheme and composition, a score outside
    /// BLAST+'s score range (<c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive): a match score ≥ 32767 (BLAST_SCORE_MAX) or a mismatch score ≤ −32768 (BLAST_SCORE_MIN) (ParamName <c>match</c> / <c>mismatch</c>), or an undefined <paramref name="kMethod"/>.</exception>
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

        ValidateMatchMismatch(match, mismatch, nameof(match), nameof(mismatch));
        double pMatch = sumSquares / (sum * sum);
        if (!(pMatch < 1.0))
            throw new ArgumentOutOfRangeException(nameof(baseFrequencies),
                "Karlin–Altschul λ is undefined: a single-base composition has no mismatch (expected score not negative).");
        return UngappedKarlinParameters(match, mismatch, pMatch, kMethod, nameof(mismatch));
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
    /// <exception cref="ArgumentOutOfRangeException">Non-positive reward, non-negative penalty, a reward ≥ 32767 or penalty
    /// ≤ −32768 (BLAST+'s score range, <c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive; ParamName
    /// <c>reward</c> / <c>penalty</c>), negative gap costs or an undefined <paramref name="kMethod"/>.</exception>
    /// <exception cref="ArgumentException">A reward/penalty or gap-cost combination BLAST+ does not support.</exception>
    public static KarlinAltschulParameters GetBlastnGappedKarlinParameters(
        int reward, int penalty, int gapOpen, int gapExtend, KarlinKMethod kMethod = KarlinKMethod.NcbiBlast)
    {
        ValidateMatchMismatch(reward, penalty, nameof(reward), nameof(penalty));
        if (gapOpen < 0)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap existence cost cannot be negative.");
        if (gapExtend < 0)
            throw new ArgumentOutOfRangeException(nameof(gapExtend), "Gap extension cost cannot be negative.");
        return BlastnGappedKarlinParameters(
            reward, penalty, gapOpen, gapExtend, kMethod, nameof(reward), nameof(penalty), nameof(gapOpen));
    }

    // GetBlastnGappedKarlinParameters after the argument checks; rewardParam/penaltyParam/gapParam are the public
    // caller's parameter names for an unsupported scheme or gap-cost combination (A9-1).
    private static KarlinAltschulParameters BlastnGappedKarlinParameters(
        int reward, int penalty, int gapOpen, int gapExtend, KarlinKMethod kMethod,
        string rewardParam, string penaltyParam, string gapParam)
    {

        int divisor = Gcd(reward, -penalty);
        if (!BlastnValueTables.TryGetValue((reward / divisor, penalty / divisor), out var table))
            throw new ArgumentException(
                $"Substitution scores {reward} and {penalty} are not supported by NCBI BLAST+ blastn statistics.",
                rewardParam);

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
            var ungapped = UngappedKarlinParameters(
                reward, penalty, MatchProbability(UniformBaseFrequency), kMethod, penaltyParam);
            return ungapped with { RoundDown = table.RoundDown, Gapped = true };
        }

        throw new ArgumentException(
            $"Gap existence and extension values {gapOpen} and {gapExtend} are not supported for substitution scores " +
            $"{reward} and {penalty}; supported: " +
            string.Join(", ", normal.Select(r => $"{(int)r[0] * divisor}/{(int)r[1] * divisor}")) +
            $", or any values at least {table.GapOpenMax * divisor}/{table.GapExtendMax * divisor}.",
            gapParam);
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
    /// Ungapped statistics of a reward/penalty scheme whose gcd-reduced form has no BLAST+ <c>blastn_values_*</c> table
    /// (e.g. +3/−5, +1/−6, +2/−1, +6/−10) use α = β = 0 as BLAST+ does (<c>Blast_GetNuclAlphaBeta</c> fails and
    /// <c>BLAST_CalcEffLengths</c> keeps its zero initial values), so ℓ = 0 and the search space is m·n: blastn 2.12.0+
    /// -ungapped, 73-nt query vs 446-nt subject → 32558 for all of them, while a tabulated scheme (+2/−3) gives 28470.
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
    /// <exception cref="ArgumentOutOfRangeException">Non-positive lengths/count, positive gap scores, an undefined
    /// <paramref name="kMethod"/>, or (ParamName <c>scoring</c>) a non-positive match, non-negative mismatch, a
    /// non-negative expected score (ungapped) or a score outside BLAST+'s score range (<c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive): a match score ≥ 32767 (BLAST_SCORE_MAX) or a mismatch score ≤ −32768 (BLAST_SCORE_MIN).</exception>
    /// <exception cref="ArgumentException">A scoring scheme BLAST+ has no gapped statistics for (ParamName <c>scoring</c>).</exception>
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

        ValidateMatchMismatch(matrix.Match, matrix.Mismatch, nameof(scoring), nameof(scoring));
        KarlinAltschulParameters p = gapped
            ? BlastnGappedKarlinParameters(matrix.Match, matrix.Mismatch, -matrix.GapOpen, -matrix.GapExtend, kMethod,
                nameof(scoring), nameof(scoring), nameof(scoring))
            : UngappedKarlinParameters(
                matrix.Match, matrix.Mismatch, MatchProbability(UniformBaseFrequency), kMethod, nameof(scoring));

        // Blast_GetNuclAlphaBeta returns s_GetNuclValuesArray's error status without setting α/β when the reduced
        // reward/penalty has no blastn_values_* table; BLAST_CalcEffLengths ignores the status and keeps its
        // initial α = β = 0, so ℓ = 0 (ungapped only — a gapped search with such a scheme is rejected above).
        if (!gapped && !HasBlastnValueTable(matrix.Match, matrix.Mismatch))
            p = p with { Alpha = 0, Beta = 0 };

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
    /// <exception cref="ArgumentException">An empty sequence or an unsupported scoring scheme (ParamName <c>scoring</c>).</exception>
    /// <exception cref="ArgumentOutOfRangeException">(ParamName <c>scoring</c>) a non-positive match, non-negative mismatch,
    /// positive gap scores, or a score outside BLAST+'s score range (<c>BlastScoreBlkMaxScoreSet</c>, blast_stat.h INT2 limits, exclusive): a match score ≥ 32767 (BLAST_SCORE_MAX) or a mismatch score ≤ −32768 (BLAST_SCORE_MIN).</exception>
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
        ValidateMatchMismatch(matrix.Match, matrix.Mismatch, nameof(scoring), nameof(scoring));
        var (score, _, _, _) = BestLocalAlignment(
            probeSequence.ToUpperInvariant(), subjectSequence.ToUpperInvariant(), matrix);
        return ComputeBlastnStatistics(score, probeSequence.Length, subjectSequence.Length, 1, matrix, gapped: true);
    }

    // s_GetNuclValuesArray: is there a blastn_values_* table for the gcd-reduced reward/penalty?
    private static bool HasBlastnValueTable(int reward, int penalty)
    {
        int divisor = Gcd(reward, -penalty);
        return BlastnValueTables.ContainsKey((reward / divisor, penalty / divisor));
    }

    // matchParam/mismatchParam: the public caller's parameter names (match/mismatch, reward/penalty, or scoring for a
    // ScoringMatrix argument), so each public method throws with its own ParamName (audit round 9, A9-1).
    private static void ValidateMatchMismatch(int match, int mismatch, string matchParam, string mismatchParam)
    {
        if (match <= 0)
            throw new ArgumentOutOfRangeException(matchParam,
                "Karlin–Altschul λ is undefined: the scoring scheme must have at least one positive score.");
        if (mismatch >= 0)
            throw new ArgumentOutOfRangeException(mismatchParam,
                "Karlin–Altschul λ is undefined: the mismatch score must be negative.");
        // BLAST+ blast_stat.h: BLAST_SCORE_MIN = INT2_MIN, BLAST_SCORE_MAX = INT2_MAX (one-letter comparison scores).
        // BlastScoreBlkMaxScoreSet skips scores <= BLAST_SCORE_MIN or >= BLAST_SCORE_MAX, so the limits are exclusive.
        if (match >= short.MaxValue)
            throw new ArgumentOutOfRangeException(matchParam, "Match score must be below BLAST_SCORE_MAX (32767).");
        if (mismatch <= short.MinValue)
            throw new ArgumentOutOfRangeException(mismatchParam, "Mismatch score must be above BLAST_SCORE_MIN (−32768).");
    }

    // p(match) = 4·p² for four equiprobable bases of frequency p (p must keep it a probability in (0, 1)).
    private static double MatchProbability(double baseFrequency)
    {
        if (double.IsNaN(baseFrequency) || baseFrequency <= 0 || baseFrequency >= 0.5)
            throw new ArgumentOutOfRangeException(nameof(baseFrequency), "Base frequency must lie in (0, 0.5).");
        return 4.0 * baseFrequency * baseFrequency;
    }

    // Blast_KarlinBlkUngappedCalc for the two-score distribution {match: pMatch, mismatch: 1 − pMatch}.
    // mismatchParam: the public caller's parameter name reported when the expected score is not negative (A9-1).
    private static KarlinAltschulParameters UngappedKarlinParameters(
        int match, int mismatch, double pMatch, KarlinKMethod kMethod, string mismatchParam)
    {
        if (kMethod is not (KarlinKMethod.ReducedLattice or KarlinKMethod.NcbiBlast))
            throw new ArgumentOutOfRangeException(nameof(kMethod), kMethod, "Unknown K method.");

        double pMismatch = 1.0 - pMatch;
        double expected = pMatch * match + pMismatch * mismatch;
        if (expected >= 0)
            throw new ArgumentOutOfRangeException(mismatchParam,
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
        // long: |int.MinValue| (e.g. -penalty for penalty = int.MinValue) does not fit in int. The result divides a
        // validated positive int, so it fits.
        long x = Math.Abs((long)a), y = Math.Abs((long)b);
        while (y != 0)
            (x, y) = (y, x % y);
        return (int)x;
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
    /// <returns>ε260 in M⁻¹·cm⁻¹; 0 for a null or empty sequence (as <see cref="CalculateExtinctionCoefficientNearestNeighbor"/>
    /// and <see cref="CalculateMolecularWeight(string)"/>).</returns>
    public static double CalculateExtinctionCoefficient(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;
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
    /// Oligonucleotide concentration from absorbance by the Beer–Lambert law, A = ε·c·l ⇒ c = A / (ε·l).
    /// </summary>
    /// <param name="absorbance260">Absorbance at 260 nm (finite; blank-subtracted readings may be ≤ 0).</param>
    /// <param name="extinctionCoefficient">Molar extinction coefficient ε₂₆₀ in L·mol⁻¹·cm⁻¹ (finite, &gt; 0),
    /// e.g. from <see cref="CalculateExtinctionCoefficient"/>.</param>
    /// <param name="pathLength">Optical path length l in cm (finite, &gt; 0; default 1).</param>
    /// <returns>Concentration in µM (mol/L × 10⁶).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A non-finite argument, or ε or l ≤ 0 (c undefined).</exception>
    public static double CalculateConcentration(
        double absorbance260,
        double extinctionCoefficient,
        double pathLength = 1.0)
    {
        if (!double.IsFinite(absorbance260))
            throw new ArgumentOutOfRangeException(nameof(absorbance260), "Absorbance must be finite.");
        if (!double.IsFinite(extinctionCoefficient) || extinctionCoefficient <= 0)
            throw new ArgumentOutOfRangeException(nameof(extinctionCoefficient), "Extinction coefficient must be finite and positive.");
        if (!double.IsFinite(pathLength) || pathLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(pathLength), "Path length must be finite and positive.");
        return absorbance260 / (extinctionCoefficient * pathLength) * 1e6;
    }

}
