// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

/// <summary>
/// Parameters for primer design. <c>MaxDinucleotideRepeats</c> limits the unsourced library screen
/// <see cref="PrimerDesigner.FindLongestDinucleotideRepeat"/> (no Primer3 counterpart; <see cref="int.MaxValue"/> disables
/// it, as in <see cref="ThermodynamicPrimerDesigner.Primer3DefaultParameters"/>). Primer3's 3′-end checks are <see cref="GcClamp"/> (PRIMER_GC_CLAMP),
/// <see cref="MaxEndGc"/> (PRIMER_MAX_END_GC) and <see cref="MaxEndStability"/> (PRIMER_MAX_END_STABILITY).
/// Two positional members predate them and are kept only for source compatibility (deprecated):
/// <c>Avoid3PrimeGC</c> is a library rule, not Primer3's — despite its name it <i>requires</i> at least one G/C
/// among the two 3′-most bases (use <see cref="GcClamp"/> = 1 for Primer3's "3′-most base is G/C", or
/// <see cref="MaxEndGc"/> to limit 3′ G/C); <c>Check3PrimeStability</c> no longer has any effect — it gated
/// ΔG(3′ pentamer) &lt; −9 kcal/mol, which no ACGT pentamer reaches (minimum −6.86, GCGCG), so it never rejected a
/// primer; Primer3's end-stability limit is <see cref="MaxEndStability"/>, always applied.
/// </summary>
public readonly record struct PrimerParameters(
    int MinLength,
    int MaxLength,
    int OptimalLength,
    double MinGcContent,
    double MaxGcContent,
    double MinTm,
    double MaxTm,
    double OptimalTm,
    int MaxHomopolymer,
    int MaxDinucleotideRepeats,
    bool Avoid3PrimeGC,
    bool Check3PrimeStability,
    PrimerStructureScreen StructureScreen = PrimerStructureScreen.Primer3Thermodynamic,
    double MaxStructureTm = PrimerDesigner.Primer3MaxStructureTm)
{
    /// <summary>
    /// Limit (°C) for every Primer3 thermodynamic structure value under
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> — Primer3's PRIMER_MAX_SELF_ANY_TH,
    /// PRIMER_MAX_SELF_END_TH, PRIMER_MAX_HAIRPIN_TH, PRIMER_PAIR_MAX_COMPL_ANY_TH and
    /// PRIMER_PAIR_MAX_COMPL_END_TH, all 47 °C by default. A value strictly greater fails.
    /// The value <c>0</c> (e.g. from <c>default(PrimerParameters)</c>) also means 47 °C.
    /// </summary>
    public double EffectiveMaxStructureTm => MaxStructureTm > 0 ? MaxStructureTm : PrimerDesigner.Primer3MaxStructureTm;

    /// <summary>
    /// PRIMER_MAX_SELF_ANY under <see cref="PrimerStructureScreen.Primer3Alignment"/>: maximum Primer3 alignment-mode
    /// <c>self_any</c> (<see cref="ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>); null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxSelfAny"/> (8.00). A value strictly greater fails. Must be in [0, 32767].
    /// </summary>
    public double? MaxSelfAny { get; init; }

    /// <summary>
    /// PRIMER_MAX_SELF_END under <see cref="PrimerStructureScreen.Primer3Alignment"/>: maximum Primer3 alignment-mode
    /// <c>self_end</c> (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>); null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxSelfEnd"/> (3.00). As in Primer3's <c>characterize_pair</c>, it also bounds the
    /// pair's reverse-orientation 3′ complementarity <c>align(right, revcomp(left), DPAL_GLOBAL_END)</c>.
    /// </summary>
    public double? MaxSelfEnd { get; init; }

    /// <summary>
    /// Primer3 per-primer penalty weights (PRIMER_WT_*) used for <see cref="PrimerCandidate.Penalty"/>; null =
    /// <see cref="ThermodynamicPrimerDesigner.DefaultPrimer3Weights"/>. The secondary-structure mode of the weights follows
    /// <see cref="PrimerParameters.StructureScreen"/> (Primer3 PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT):
    /// <see cref="PrimerStructureScreen.Primer3Alignment"/> uses PRIMER_WT_SELF_ANY / _SELF_END
    /// (<see cref="Primer3PenaltyWeights.SelfAny"/>/<see cref="Primer3PenaltyWeights.SelfEnd"/>) × the dpal scores,
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> PRIMER_WT_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH
    /// with the ntthal Tm values; the <see cref="PrimerStructureScreen.Heuristic"/> screen has no Primer3 structure
    /// values, so its structure terms are 0. When a structure weight is non-zero the values are computed for every
    /// candidate (as Primer3's <c>calc_and_check_oligo_features</c> does). PRIMER_WT_END_STABILITY
    /// (<see cref="Primer3PenaltyWeights.EndStability"/>) multiplies Primer3's <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="PrimerDesigner.Calculate3PrimeStability"/> (positive kcal/mol), as
    /// <c>p_obj_fn</c> does for left/right primers (internal oligos have no such term).
    /// </summary>
    public Primer3PenaltyWeights? PenaltyWeights { get; init; }

    /// <summary>
    /// PRIMER_SALT_MONOVALENT: monovalent cation concentration (mM, &gt; 0) of the primer reaction; null = Primer3's
    /// default <see cref="PrimerDesigner.Primer3MonovalentMillimolar"/> (50 mM). Primer3 uses the primer conditions
    /// (<c>p_args</c>) for the primer Tm (<c>seqtm</c>), the ntthal self-dimer / 3′ self-dimer / hairpin and pair
    /// complementarity values (<c>create_thal_arg_holder(p_args)</c>) and the product Tm (<c>long_seq_tm</c>).
    /// The probe-side counterpart is <see cref="ThermodynamicProbeDesigner.ProbeParameters.MonovalentMillimolar"/>; the internal oligo
    /// of a pair uses <see cref="ThermodynamicProbeDesigner.Primer3ProbeSettings.MonovalentMillimolar"/> (PRIMER_INTERNAL_SALT_MONOVALENT).
    /// </summary>
    public double? MonovalentMillimolar { get; init; }

    /// <summary>PRIMER_SALT_DIVALENT: Mg²⁺ concentration (mM, ≥ 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DivalentMillimolar"/> (1.5 mM).</summary>
    public double? DivalentMillimolar { get; init; }

    /// <summary>PRIMER_DNTP_CONC: dNTP concentration (mM, ≥ 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DntpMillimolar"/> (0.6 mM).</summary>
    public double? DntpMillimolar { get; init; }

    /// <summary>PRIMER_DNA_CONC: primer (oligo) concentration (nM, &gt; 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DnaConcentrationNanomolar"/> (50 nM).</summary>
    public double? DnaConcentrationNanomolar { get; init; }

    /// <summary>
    /// PRIMER_OPT_GC_PERCENT: the GC optimum of the PRIMER_WT_GC_PERCENT_GT/_LT penalty terms
    /// (<see cref="Primer3PenaltyWeights.GcGt"/>/<see cref="Primer3PenaltyWeights.GcLt"/>); null =
    /// undefined, as in Primer3's code (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>; the manual's 50 is not
    /// applied). It is inert while both GC weights are 0 (Primer3's default); a non-zero GC weight without it is rejected with
    /// Primer3's <c>_pr_data_control</c> error ("Primer GC content is part of objective function while optimum gc_content is
    /// not defined", <see cref="ArgumentException"/>).
    /// </summary>
    public double? OptimalGcPercent { get; init; }

    /// <summary>
    /// PRIMER_MAX_END_STABILITY (kcal/mol, ≥ 0): a primer fails when its Primer3 <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="PrimerDesigner.Calculate3PrimeStability"/> (the positive −ΔG°37 of the
    /// 3′ pentamer) is strictly greater; null = Primer3's default <see cref="PrimerDesigner.Primer3MaxEndStability"/>
    /// (100, never exceeded). Primer3 applies it to left/right primers only (no internal-oligo counterpart).
    /// </summary>
    public double? MaxEndStability { get; init; }

    /// <summary>
    /// PRIMER_GC_CLAMP: number of consecutive G/C required at the 3′ end of a primer (0 = none, Primer3's default);
    /// must not exceed <see cref="MinLength"/> (Primer3 <c>_pr_data_control</c>). A negative value requires nothing,
    /// as in Primer3. Left/right primers only.
    /// </summary>
    public int GcClamp { get; init; }

    /// <summary>
    /// PRIMER_MAX_END_GC (0–5): maximum number of G/C among the five 3′-most bases of a primer; null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxEndGc"/> (5, not checked). Left/right primers only.
    /// </summary>
    public int? MaxEndGc { get; init; }

    /// <summary>
    /// PRIMER_MISPRIMING_LIBRARY (primer3-py <c>misprime_lib</c>): when set and non-empty, every primer gets Primer3's
    /// library mispriming score (<see cref="ThermodynamicPrimerDesigner.CalculateLibraryMispriming"/>,
    /// <see cref="PrimerCandidate.LibraryMispriming"/>) and fails when any entry's weighted score exceeds
    /// <see cref="MaxLibraryMispriming"/>; pairs are limited by <see cref="PrimerPairOptions.MaxLibraryMispriming"/>.
    /// As in Primer3's <c>calc_and_check_oligo_features</c> / <c>characterize_pair</c>, the pair search scores primers
    /// when they are picked only if PRIMER_WT_LIBRARY_MISPRIMING (<see cref="Primer3PenaltyWeights.LibraryMispriming"/>)
    /// is non-zero, otherwise lazily when a pair is characterized (same results). Null = no library (Primer3's default).
    /// Not part of the JSON form of the parameters (the MCP tools take the library as a separate name → sequence argument).
    /// </summary>
    [System.Text.Json.Serialization.JsonIgnore]
    public PrimerMisprimingLibrary? MisprimingLibrary { get; init; }

    /// <summary>
    /// PRIMER_MAX_LIBRARY_MISPRIMING: maximum weighted library score of one primer; null = Primer3's default
    /// <see cref="ThermodynamicPrimerDesigner.Primer3MaxLibraryMispriming"/> (12.00). Primer3 compares with the value truncated to a
    /// C <c>short</c> (12.9 acts as 12). Under <see cref="PrimerStructureScreen.Primer3Alignment"/> it must not exceed
    /// 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double? MaxLibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS: false (null, Primer3's default 0 — <c>pr_set_default_global_args_2</c>, used by
    /// primer3-py; the "version 1" defaults had 1) = an IUPAC code in a library entry never aligns and N scores −0.25;
    /// true (1) = an IUPAC code matches every base it represents (so runs of N match any primer).
    /// </summary>
    public bool? LibraryAmbiguityCodesConsensus { get; init; }

    /// <summary>
    /// PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT: false (Primer3's default 0) = template mispriming is scored with dpal
    /// (<see cref="MaxTemplateMispriming"/>, PRIMER_WT_TEMPLATE_MISPRIMING); true (1) = with the ntthal THAL_END1 Tm
    /// (<see cref="MaxTemplateMisprimingTh"/>, PRIMER_WT_TEMPLATE_MISPRIMING_TH; template ≤ 10000 nt). Independent of
    /// <see cref="StructureScreen"/> (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT). See
    /// <see cref="ThermodynamicPrimerDesigner.CalculateTemplateMispriming"/>.
    /// </summary>
    public bool ThermodynamicTemplateAlignment { get; init; }

    /// <summary>
    /// PRIMER_MAX_TEMPLATE_MISPRIMING (alignment mode): in <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> /
    /// <see cref="ThermodynamicPrimerDesigner.DesignPrimerPairs"/> a primer whose template mispriming score
    /// (<see cref="ThermodynamicPrimerDesigner.CalculateTemplateMispriming"/>, <see cref="PrimerCandidate.TemplateMispriming"/>) is
    /// strictly greater fails. Null = Primer3's default <see cref="PrimerDesigner.Primer3UndefinedTemplateMispriming"/>
    /// (−100); a negative value is not checked. Must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double? MaxTemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_MAX_TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C): as <see cref="MaxTemplateMispriming"/> for the ntthal
    /// template-mispriming Tm; null = −100 (not checked).
    /// </summary>
    public double? MaxTemplateMisprimingTh { get; init; }

    /// <summary>
    /// PRIMER_ANNEALING_TEMP (°C; null = Primer3's default −10 = off). When &gt; 0 Primer3 computes each primer's
    /// fraction bound at this temperature (<see cref="ThermodynamicPrimerDesigner.CalculateFractionBoundPrimer3"/>, reported in
    /// <see cref="PrimerCandidate.Bound"/>), rejects primers outside [<see cref="MinBound"/>, <see cref="MaxBound"/>] and
    /// adds the PRIMER_WT_BOUND_GT / _LT terms (<see cref="Primer3PenaltyWeights.BoundGt"/> /
    /// <see cref="Primer3PenaltyWeights.BoundLt"/>) around <see cref="OptBound"/>. Must not exceed 100 °C. In a primer-pair
    /// design it is also the internal oligo's annealing temperature (one global Primer3 setting).
    /// </summary>
    public double? AnnealingTemperature { get; init; }

    /// <summary>PRIMER_MIN_BOUND (% bound; null = Primer3's −10); checked only when <see cref="AnnealingTemperature"/> &gt; 0.
    /// A primer longer than 36 bases has no bound value (Primer3 OLIGOTM_ERROR = −999999.9999) and fails it.</summary>
    public double? MinBound { get; init; }

    /// <summary>PRIMER_MAX_BOUND (% bound; null = Primer3's 110); checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
    public double? MaxBound { get; init; }

    /// <summary>PRIMER_OPT_BOUND (% bound; null = Primer3's 97): the optimum of the PRIMER_WT_BOUND_GT/_LT terms; must lie
    /// in [<see cref="MinBound"/>, <see cref="MaxBound"/>] (Primer3 <c>_pr_data_control</c>).</summary>
    public double? OptBound { get; init; }

    /// <summary>
    /// PRIMER_MIN_QUALITY (default 0): with SEQUENCE_QUALITY (<see cref="PrimerPairOptions.SequenceQuality"/>) a primer
    /// whose minimum base quality (<see cref="ThermodynamicPrimerDesigner.CalculateSequenceQualityPrimer3"/>,
    /// <see cref="PrimerCandidate.MinSequenceQuality"/>) is lower fails (Primer3 <c>sequence_quality_is_ok</c>). A non-zero
    /// value requires quality data and must lie in [<see cref="QualityRangeMin"/>, <see cref="QualityRangeMax"/>]
    /// (Primer3 <c>_pr_data_control</c>; <see cref="ThermodynamicPrimerDesigner.EvaluatePrimer"/> has no quality data, so it rejects it).
    /// </summary>
    public int MinQuality { get; init; }

    /// <summary>
    /// PRIMER_MIN_END_QUALITY (default 0): with SEQUENCE_QUALITY a primer whose minimum quality over its five 3′-most bases
    /// is lower fails (checked after <see cref="MinQuality"/>; left/right primers only). Not range-checked by Primer3.
    /// </summary>
    public int MinEndQuality { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MIN (null = Primer3's 0): every SEQUENCE_QUALITY value must be ≥ it. One global
    /// Primer3 setting: it also applies to the internal oligo of a pair.</summary>
    public int? QualityRangeMin { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MAX (null = Primer3's 100): every SEQUENCE_QUALITY value must be ≤ it; the PRIMER_WT_SEQ_QUAL
    /// term is weight × (QualityRangeMax − min quality). One global Primer3 setting (also used for the internal oligo).</summary>
    public int? QualityRangeMax { get; init; }

    /// <summary>Effective PRIMER_QUALITY_RANGE_MIN.</summary>
    public int EffectiveQualityRangeMin => QualityRangeMin ?? ThermodynamicPrimerDesigner.Primer3QualityRangeMin;

    /// <summary>Effective PRIMER_QUALITY_RANGE_MAX.</summary>
    public int EffectiveQualityRangeMax => QualityRangeMax ?? ThermodynamicPrimerDesigner.Primer3QualityRangeMax;

    /// <summary>PRIMER_ANNEALING_TEMP in effect (Primer3 default −10 = off).</summary>
    public double EffectiveAnnealingTemperature => AnnealingTemperature ?? ThermodynamicPrimerDesigner.Primer3DefaultAnnealingTemperature;

    /// <summary>PRIMER_MIN_BOUND in effect (default −10).</summary>
    public double EffectiveMinBound => MinBound ?? ThermodynamicPrimerDesigner.Primer3MinBound;

    /// <summary>PRIMER_MAX_BOUND in effect (default 110).</summary>
    public double EffectiveMaxBound => MaxBound ?? ThermodynamicPrimerDesigner.Primer3MaxBound;

    /// <summary>PRIMER_OPT_BOUND in effect (default 97).</summary>
    public double EffectiveOptBound => OptBound ?? ThermodynamicPrimerDesigner.Primer3OptBound;

    /// <summary>Effective PRIMER_MAX_TEMPLATE_MISPRIMING.</summary>
    public double EffectiveMaxTemplateMispriming => MaxTemplateMispriming ?? PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>Effective PRIMER_MAX_TEMPLATE_MISPRIMING_TH.</summary>
    public double EffectiveMaxTemplateMisprimingTh => MaxTemplateMisprimingTh ?? PrimerDesigner.Primer3UndefinedTemplateMispriming;

    // The per-primer limit and weight of the active template alignment mode.
    public double ActiveMaxTemplateMispriming =>
        ThermodynamicTemplateAlignment ? EffectiveMaxTemplateMisprimingTh : EffectiveMaxTemplateMispriming;

    public double ActiveTemplateMisprimingWeight =>
        ThermodynamicTemplateAlignment
            ? PenaltyWeights?.TemplateMisprimingTh ?? 0.0
            : PenaltyWeights?.TemplateMispriming ?? 0.0;

    /// <summary>Effective PRIMER_MAX_LIBRARY_MISPRIMING.</summary>
    public double EffectiveMaxLibraryMispriming => MaxLibraryMispriming ?? ThermodynamicPrimerDesigner.Primer3MaxLibraryMispriming;

    /// <summary>Effective PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS.</summary>
    public bool EffectiveLibraryAmbiguityCodesConsensus => LibraryAmbiguityCodesConsensus ?? false;

    // The library in use: null when absent or empty (Primer3 seq_lib_num_seq == 0).
    public PrimerMisprimingLibrary? ActiveLibrary => MisprimingLibrary is { Count: > 0 } lib ? lib : null;

    /// <summary>Effective PRIMER_MAX_END_STABILITY (kcal/mol).</summary>
    public double EffectiveMaxEndStability => MaxEndStability ?? PrimerDesigner.Primer3MaxEndStability;

    /// <summary>Effective PRIMER_MAX_END_GC.</summary>
    public int EffectiveMaxEndGc => MaxEndGc ?? PrimerDesigner.Primer3MaxEndGc;

    /// <summary>Effective PRIMER_SALT_MONOVALENT (mM).</summary>
    public double EffectiveMonovalentMillimolar => MonovalentMillimolar ?? PrimerDesigner.Primer3MonovalentMillimolar;

    /// <summary>Effective PRIMER_SALT_DIVALENT (mM).</summary>
    public double EffectiveDivalentMillimolar => DivalentMillimolar ?? PrimerDesigner.Primer3DivalentMillimolar;

    /// <summary>Effective PRIMER_DNTP_CONC (mM).</summary>
    public double EffectiveDntpMillimolar => DntpMillimolar ?? PrimerDesigner.Primer3DntpMillimolar;

    /// <summary>Effective PRIMER_DNA_CONC (nM).</summary>
    public double EffectiveDnaConcentrationNanomolar => DnaConcentrationNanomolar ?? PrimerDesigner.Primer3DnaConcentrationNanomolar;

    /// <summary>Effective PRIMER_MAX_SELF_ANY (<see cref="MaxSelfAny"/> or 8.00).</summary>
    public double EffectiveMaxSelfAny => MaxSelfAny ?? PrimerDesigner.Primer3MaxSelfAny;

    /// <summary>Effective PRIMER_MAX_SELF_END (<see cref="MaxSelfEnd"/> or 3.00).</summary>
    public double EffectiveMaxSelfEnd => MaxSelfEnd ?? PrimerDesigner.Primer3MaxSelfEnd;

    // Primer3 _pr_data_control: PRIMER_WT_GC_PERCENT_GT/_LT need PRIMER_OPT_GC_PERCENT.
    public void ValidateGcOptimum(string paramName)
    {
        if ((PenaltyWeights is { } w && (w.GcGt != 0 || w.GcLt != 0)) && OptimalGcPercent is null)
            throw new ArgumentException(PrimerDesigner.PrimerGcOptimumUndefinedMessage, paramName);
    }

    // Primer3 _pr_data_control checks of the primer conditions.
    public void ValidateConditions(string paramName)
    {
        PrimerDesigner.ValidatePrimer3Conditions(EffectiveMonovalentMillimolar, EffectiveDivalentMillimolar,
            EffectiveDntpMillimolar, EffectiveDnaConcentrationNanomolar, paramName);
        if (OptimalGcPercent is { } optGc && !double.IsFinite(optGc))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_OPT_GC_PERCENT must be finite.");
        // _pr_data_control: PRIMER_MAX_END_GC must be between 0 to 5; PRIMER_MAX_END_STABILITY must be non-negative;
        // PRIMER_GC_CLAMP > PRIMER_MIN_SIZE.
        if (EffectiveMaxEndGc is < 0 or > PrimerDesigner.Primer3MaxEndGc)
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_MAX_END_GC must be between 0 to 5.");
        if (!(EffectiveMaxEndStability >= 0))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_MAX_END_STABILITY must be non-negative.");
        if (GcClamp > MinLength)
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_GC_CLAMP > PRIMER_MIN_SIZE.");
        // _pr_data_control: PRIMER_MAX_LIBRARY_MISPRIMING > SHRT_MAX (alignment mode); a library weight without a library.
        double maxLib = EffectiveMaxLibraryMispriming;
        if (double.IsNaN(maxLib) || Math.Abs(maxLib) >= int.MaxValue
            || (maxLib > short.MaxValue && StructureScreen == PrimerStructureScreen.Primer3Alignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_MAX_LIBRARY_MISPRIMING.");
        if ((PenaltyWeights?.LibraryMispriming ?? 0) != 0 && ActiveLibrary is null)
            throw new ArgumentException("Mispriming score is part of objective function, but mispriming library is not defined (Primer3 _pr_data_control).", paramName);
        // _pr_data_control: PRIMER_MAX_TEMPLATE_MISPRIMING > SHRT_MAX in alignment mode.
        if (double.IsNaN(EffectiveMaxTemplateMispriming) || double.IsNaN(EffectiveMaxTemplateMisprimingTh)
            || (EffectiveMaxTemplateMispriming > short.MaxValue && !ThermodynamicTemplateAlignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_MAX_TEMPLATE_MISPRIMING.");
        // A negative template weight makes Primer3 skip the score (_pr_need_template_mispriming: weight > 0) but still
        // add the term (p_obj_fn: weight ≠ 0), which fails its PR_ASSERT.
        if (!((PenaltyWeights?.TemplateMispriming ?? 0) >= 0) || !((PenaltyWeights?.TemplateMisprimingTh ?? 0) >= 0))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_WT_TEMPLATE_MISPRIMING[_TH] must be ≥ 0.");
        // _pr_data_control: PRIMER_OPT_BOUND within [PRIMER_MIN_BOUND, PRIMER_MAX_BOUND]; PRIMER_ANNEALING_TEMP ≤ 100.
        ThermodynamicPrimerDesigner.ValidatePrimer3Bound(EffectiveAnnealingTemperature, EffectiveMinBound, EffectiveMaxBound,
            EffectiveOptBound, internalOligo: false, paramName);
    }
}

/// <summary>
/// Pair-level options of <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> / <see cref="ThermodynamicPrimerDesigner.DesignPrimerPairs"/>,
/// named after the Primer3 tags they implement. Defaults are Primer3's except
/// <see cref="MaxTmDifference"/> (library default 5 °C; Primer3 100 °C — use <see cref="Primer3Defaults"/>).
/// </summary>
public sealed record PrimerPairOptions
{
    /// <summary>Library defaults (Primer3 defaults with PRIMER_PAIR_MAX_DIFF_TM = 5 °C).</summary>
    public static PrimerPairOptions Default { get; } = new();

    /// <summary>Primer3's own defaults (PRIMER_PAIR_MAX_DIFF_TM = 100 °C).</summary>
    public static PrimerPairOptions Primer3Defaults { get; } = new() { MaxTmDifference = PrimerDesigner.Primer3MaxPairTmDifference };

    /// <summary>
    /// PRIMER_PRODUCT_SIZE_RANGE: product-size ranges in order of preference (default 100–300 bp). A later
    /// range is used only when no further pair fits the earlier ones (Primer3 <c>choose_pair_or_triple</c>).
    /// </summary>
    public IReadOnlyList<ProductSizeRange> ProductSizeRanges { get; init; } = [new ProductSizeRange(100, 300)];

    /// <summary>PRIMER_PAIR_MAX_DIFF_TM: maximum |Tm_left − Tm_right| (°C); library default <see cref="PrimerDesigner.MaxPairTmDifference"/> = 5.</summary>
    public double MaxTmDifference { get; init; } = PrimerDesigner.MaxPairTmDifference;

    /// <summary>PRIMER_NUM_RETURN: number of pairs returned by <see cref="ThermodynamicPrimerDesigner.DesignPrimerPairs"/> (default 5, ≥ 1).</summary>
    public int NumReturn { get; init; } = 5;

    /// <summary>SEQUENCE_INCLUDED_REGION (start, length): primers and products must lie inside it (default: whole template).</summary>
    public (int Start, int Length)? IncludedRegion { get; init; }

    /// <summary>PRIMER_PRODUCT_OPT_SIZE (bp); required when a product-size weight is non-zero.</summary>
    public int? ProductOptSize { get; init; }

    /// <summary>PRIMER_PRODUCT_OPT_TM (°C); required when a product-Tm weight is non-zero.</summary>
    public double? ProductOptTm { get; init; }

    /// <summary>PRIMER_PRODUCT_MIN_TM (°C): pairs whose product Tm is lower fail (default: no limit).</summary>
    public double? ProductMinTm { get; init; }

    /// <summary>PRIMER_PRODUCT_MAX_TM (°C): pairs whose product Tm is higher fail (default: no limit).</summary>
    public double? ProductMaxTm { get; init; }

    /// <summary>
    /// PRIMER_PAIR_MAX_COMPL_ANY: maximum Primer3 alignment-mode pair <c>compl_any</c> under
    /// <see cref="PrimerStructureScreen.Primer3Alignment"/> (default 8.00; must be in [0, 32767]).
    /// </summary>
    public double MaxComplAny { get; init; } = PrimerDesigner.Primer3MaxPairComplAny;

    /// <summary>
    /// PRIMER_PAIR_MAX_COMPL_END: maximum Primer3 alignment-mode pair <c>compl_end</c> (left vs. right orientation)
    /// under <see cref="PrimerStructureScreen.Primer3Alignment"/> (default 3.00; must be in [0, 32767]).
    /// </summary>
    public double MaxComplEnd { get; init; } = PrimerDesigner.Primer3MaxPairComplEnd;

    /// <summary>PRIMER_PAIR_WT_* pair objective weights (default Primer3's).</summary>
    public Primer3PairWeights Weights { get; init; } = new();

    /// <summary>
    /// PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (default 24.00): with a <see cref="PrimerParameters.MisprimingLibrary"/>, a pair
    /// fails when its library mispriming score (Primer3 <c>pair_repeat_sim</c>: the maximum over library entries of the
    /// integer part of left score + right score, <see cref="PrimerPairResult.LibraryMispriming"/>) is strictly greater.
    /// Under <see cref="PrimerStructureScreen.Primer3Alignment"/> it must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double MaxLibraryMispriming { get; init; } = ThermodynamicPrimerDesigner.Primer3PairMaxLibraryMispriming;

    /// <summary>
    /// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING (alignment mode; default −100 = not checked): a pair whose template mispriming
    /// score (<see cref="PrimerPairResult.TemplateMispriming"/>) is strictly greater fails; a negative value is not checked.
    /// Must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double MaxTemplateMispriming { get; init; } = PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>
    /// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C; default −100). As in Primer3's
    /// <c>characterize_pair</c> the pair score is computed when this is ≥ 0 or PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH &gt; 0,
    /// and a pair fails when the limit is non-zero and the score exceeds it — so with the default −100 and a non-zero
    /// pair weight every pair fails, and 0 means no limit (Primer3 behaviour, reproduced).
    /// </summary>
    public double MaxTemplateMisprimingTh { get; init; } = PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>
    /// PRIMER_INSIDE_PENALTY (default <see cref="ThermodynamicPrimerDesigner.Primer3DefaultInsidePenalty"/> = −1). With the default
    /// inside (−1) and outside (0) penalties primers never overlap the target. Any other combination makes Primer3 score
    /// every primer's 3′ end against the (single) target instead (<see cref="ThermodynamicPrimerDesigner.CalculatePositionPenaltyPrimer3"/>):
    /// a primer may then extend into the target as long as its 3′ end does not pass the target's far end, a 3′ end
    /// inside the target costs <c>InsidePenalty</c> per base, one outside it <see cref="OutsidePenalty"/> per base of
    /// distance (× PRIMER_WT_POS_PENALTY, <see cref="Primer3PenaltyWeights.PositionPenalty"/>), and a pair must still
    /// span the target (left 3′ end before the right primer's 3′ end). Primer3 applies the value as given, so the default
    /// −1 makes inside positions negative when only <see cref="OutsidePenalty"/> is changed.
    /// </summary>
    public double InsidePenalty { get; init; } = ThermodynamicPrimerDesigner.Primer3DefaultInsidePenalty;

    /// <summary>PRIMER_OUTSIDE_PENALTY (default 0): per-base penalty of a 3′ end outside the target when the position
    /// penalties are not Primer3's defaults (see <see cref="InsidePenalty"/>).</summary>
    public double OutsidePenalty { get; init; } = ThermodynamicPrimerDesigner.Primer3DefaultOutsidePenalty;

    // Primer3 _PR_DEFAULT_POSITION_PENALTIES.
    public bool DefaultPositionPenalties => ThermodynamicPrimerDesigner.IsDefaultPositionPenalties(InsidePenalty, OutsidePenalty);

    /// <summary>
    /// PRIMER_PICK_INTERNAL_OLIGO: pick a hybridization (internal) oligo for every pair — the lowest-penalty
    /// oligo of the Primer3 internal-oligo list (<see cref="ThermodynamicProbeDesigner.DesignProbesPrimer3"/> rules) lying
    /// strictly between the primers; a pair without one fails (Primer3 <c>choose_internal_oligo</c>).
    /// </summary>
    public bool PickInternalOligo { get; init; }

    /// <summary>PRIMER_INTERNAL_* settings of the internal oligo (default Primer3's).</summary>
    public ThermodynamicProbeDesigner.Primer3ProbeSettings? InternalOligo { get; init; }

    /// <summary>
    /// SEQUENCE_QUALITY: one integer quality per template base (null or empty = none, Primer3's default). With it every
    /// primer (and internal oligo) gets its minimum base quality (<see cref="ThermodynamicPrimerDesigner.CalculateSequenceQualityPrimer3"/>,
    /// <see cref="PrimerCandidate.MinSequenceQuality"/>), checked against <see cref="PrimerParameters.MinQuality"/> /
    /// <see cref="PrimerParameters.MinEndQuality"/> (<see cref="ThermodynamicProbeDesigner.Primer3ProbeSettings.MinQuality"/>) and
    /// weighted by PRIMER_WT_SEQ_QUAL (<see cref="Primer3PenaltyWeights.SequenceQuality"/>). Its length must equal the
    /// template length and every value must lie in [<see cref="PrimerParameters.QualityRangeMin"/>,
    /// <see cref="PrimerParameters.QualityRangeMax"/>] (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public IReadOnlyList<int>? SequenceQuality { get; init; }

    /// <summary>
    /// PRIMER_MASK_TEMPLATE: mask the included region with Primer3's k-mer masker (<see cref="ThermodynamicPrimerDesigner.MaskTemplatePrimer3"/>,
    /// <see cref="MaskKmerLists"/>, <see cref="MaskFailureRate"/>, <see cref="MaskFivePrimeDirection"/>,
    /// <see cref="MaskThreePrimeDirection"/>): a left primer whose 3′ base is masked on the forward copy, or a right primer whose
    /// 3′ base (its leftmost template base) is masked on the reverse copy, is rejected (Primer3 <c>is_lowercase_masked</c>), and
    /// every primer gets its predicted failure rate (<see cref="ThermodynamicPrimerDesigner.CalculateMaskFailureRatePrimer3"/>,
    /// <see cref="PrimerCandidate.MaskFailureRate"/>) for PRIMER_WT_MASK_FAILURE_RATE
    /// (<see cref="Primer3PenaltyWeights.MaskFailureRate"/>). Requires <see cref="MaskKmerLists"/> (primer3-py: "masking template
    /// chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified"). Internal oligos are not masked (Primer3 checks them on the
    /// unmasked, case-preserving template — <see cref="LowercaseMasking"/>, which masking implies).
    /// </summary>
    public bool MaskTemplate { get; init; }

    /// <summary>
    /// PRIMER_LOWERCASE_MASKING (default false): with a case-preserving template (the <c>string</c> overloads of
    /// <see cref="ThermodynamicPrimerDesigner.DesignPrimers(string, int, int, PrimerParameters?, PrimerPairOptions?)"/> /
    /// <see cref="ThermodynamicPrimerDesigner.DesignPrimerPairs(string, int, int, PrimerParameters?, PrimerPairOptions?)"/>) a left primer
    /// or internal oligo whose 3′ base (rightmost template base), or a right primer whose 3′ base (leftmost template base), is a
    /// lower-case a/c/g/t is rejected (Primer3 <c>calc_and_check_oligo_features</c> → <c>is_lowercase_masked</c> on
    /// <c>trimmed_orig_seq</c>); lower-case bases elsewhere in the oligo are accepted. <see cref="MaskTemplate"/> implies it
    /// (primer3-py sets <c>lowercase_masking = mask_template</c>; primers are then checked on the masked copies, which keep the
    /// input's lower case). No effect with a <see cref="DnaSequence"/> template (upper case).
    /// </summary>
    public bool LowercaseMasking { get; init; }

    /// <summary>The masker's 11-mer / 16-mer genome count lists (PRIMER_MASK_KMERLIST_PATH / _PREFIX; see
    /// <see cref="PrimerMaskingKmerLists.FromGenomeTester4Files"/>). Not part of the JSON form.</summary>
    [System.Text.Json.Serialization.JsonIgnore]
    public PrimerMaskingKmerLists? MaskKmerLists { get; init; }

    /// <summary>PRIMER_MASK_FAILURE_RATE (default 0.1): a template position is masked when the failure rate of a primer ending
    /// there exceeds it (0 = no masking by score).</summary>
    public double MaskFailureRate { get; init; } = ThermodynamicPrimerDesigner.Primer3MaskFailureRate;

    /// <summary>PRIMER_MASK_5P_DIRECTION (default 1, ≥ 0): bases masked from a high-failure 3′ end towards the primer's 5′ end
    /// (including the 3′ base).</summary>
    public int MaskFivePrimeDirection { get; init; } = ThermodynamicPrimerDesigner.Primer3MaskFivePrimeDirection;

    /// <summary>PRIMER_MASK_3P_DIRECTION (default 0, 0–4984): bases masked beyond a high-failure 3′ end.</summary>
    public int MaskThreePrimeDirection { get; init; } = ThermodynamicPrimerDesigner.Primer3MaskThreePrimeDirection;

    /// <summary>
    /// PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE (≥ −1, default −1): once <see cref="ThermodynamicPrimerDesigner.DesignPrimerPairs"/> has
    /// selected a pair, no later pair may use a left primer whose 3′ end is fewer than this many bases from the selected
    /// left primer's 3′ end; 0 excludes only the identical left primer; −1 allows reuse (Primer3
    /// <c>choose_pair_or_triple</c>, <c>left_oligo_in_pair_overlaps_used_oligo</c>).
    /// </summary>
    public int MinLeftThreePrimeDistance { get; init; } = -1;

    /// <summary>
    /// PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE (≥ −1, default −1): as <see cref="MinLeftThreePrimeDistance"/> for right
    /// primers (3′ end = the leftmost top-strand base, <see cref="PrimerCandidate.Position"/>).
    /// </summary>
    public int MinRightThreePrimeDistance { get; init; } = -1;

    /// <summary>
    /// PRIMER_MIN_THREE_PRIME_DISTANCE: Primer3's shorthand that sets <see cref="MinLeftThreePrimeDistance"/> and
    /// <see cref="MinRightThreePrimeDistance"/> to the same value. Primer3's shorthand also sets
    /// PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE, which <c>choose_pair_or_triple</c> reads only with a non-empty
    /// SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST; junction lists are outside this API's scope, so — exactly as in Primer3
    /// without that list — internal oligos are never excluded by distance (primer3-py 2.3.1 parity, audit round 5, A5-3).
    /// Reading it returns the common value, or null when the left and right distances differ.
    /// </summary>
    public int? MinThreePrimeDistance
    {
        get => MinLeftThreePrimeDistance == MinRightThreePrimeDistance ? MinLeftThreePrimeDistance : null;
        init
        {
            MinLeftThreePrimeDistance = value ?? -1;
            MinRightThreePrimeDistance = value ?? -1;
        }
    }
}

/// <summary>
/// Measured properties of a single primer used as input to the Primer3 penalty
/// objective (<see cref="ThermodynamicPrimerDesigner.CalculatePrimer3Penalty"/>). Units mirror
/// Primer3's <c>p_obj_fn</c>: Tm in °C, length in bases, GC in percent (0–100),
/// self/3'-complementarity as local-alignment scores, and the count of N bases.
/// </summary>
/// <param name="Tm">Primer melting temperature in °C.</param>
/// <param name="Length">Primer length in bases.</param>
/// <param name="GcPercent">GC content as a percentage in [0, 100].</param>
/// <param name="SelfAny">Self-complementarity: the local-alignment score PRIMER_SELF_ANY when
/// <see cref="Primer3PenaltyWeights.ThermodynamicOligoAlignment"/> is false, or the self-dimer Tm in °C
/// (PRIMER_SELF_ANY_TH) when it is true — Primer3 stores both in the same field.</param>
/// <param name="SelfEnd">3'-self-complementarity: PRIMER_SELF_END score (alignment mode) or
/// PRIMER_SELF_END_TH Tm in °C (thermodynamic mode).</param>
/// <param name="NumNs">Number of ambiguous N bases in the primer.</param>
/// <param name="HairpinTh">Hairpin Tm in °C (PRIMER_HAIRPIN_TH); used only in thermodynamic mode.</param>
/// <param name="EndStability">3'-end stability ΔG magnitude in kcal/mol (PRIMER_END_STABILITY).</param>
public readonly record struct Primer3PenaltyInputs(
    double Tm,
    int Length,
    double GcPercent,
    double SelfAny = 0.0,
    double SelfEnd = 0.0,
    int NumNs = 0,
    double HairpinTh = 0.0,
    double EndStability = 0.0)
{
    /// <summary>
    /// Library mispriming score (Primer3 <c>repeat_sim.score[repeat_sim.max]</c>, PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING;
    /// <see cref="ThermodynamicPrimerDesigner.CalculateLibraryMispriming"/>); 0 without a mispriming library.
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// Template mispriming score (Primer3 <c>oligo_max_template_mispriming</c>, PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING[_TH];
    /// <see cref="ThermodynamicPrimerDesigner.CalculateTemplateMispriming"/>); 0 when not computed.
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// Fraction bound in % (Primer3 <c>h->bound</c>, <see cref="ThermodynamicPrimerDesigner.CalculateFractionBoundPrimer3"/>), or
    /// <c>null</c> when the bound terms do not apply: for a left/right primer Primer3 adds the PRIMER_WT_BOUND_GT/_LT
    /// terms only when PRIMER_ANNEALING_TEMP &gt; 0 (pass null otherwise). The internal-oligo branch of <c>p_obj_fn</c>
    /// has no such gate: an internal oligo without a bound value carries Primer3's OLIGOTM_ERROR (−999999.9999).
    /// </summary>
    public double? Bound { get; init; }

    /// <summary>
    /// Position penalty (Primer3 <c>h->position_penalty</c>, <see cref="ThermodynamicPrimerDesigner.CalculatePositionPenaltyPrimer3"/>);
    /// 0 under Primer3's default inside / outside penalties.
    /// </summary>
    public double PositionPenalty { get; init; }

    /// <summary>
    /// Minimum base quality of the oligo (Primer3 <c>h->seq_quality</c>, <see cref="ThermodynamicPrimerDesigner.CalculateSequenceQualityPrimer3"/>),
    /// or <c>null</c> without SEQUENCE_QUALITY (Primer3 then sets it to <see cref="QualityRangeMax"/>, so the
    /// PRIMER_WT_SEQ_QUAL term is 0).
    /// </summary>
    public int? SequenceQuality { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MAX of the PRIMER_WT_SEQ_QUAL term (default 100).</summary>
    public int QualityRangeMax { get; init; } = ThermodynamicPrimerDesigner.Primer3QualityRangeMax;

    /// <summary>Predicted PCR failure rate of the primer (Primer3 <c>h->failure_rate</c>,
    /// <see cref="ThermodynamicPrimerDesigner.CalculateMaskFailureRatePrimer3"/>); 0 without template masking.</summary>
    public double MaskFailureRate { get; init; }
}

/// <summary>
/// Weights for the Primer3 per-primer penalty objective (the <c>PRIMER_WT_*</c>
/// parameters). Tm, GC and size each have separate "greater-than" (_gt) and
/// "less-than" (_lt) weights, applied one-sidedly relative to the optimum.
/// Defaults are <see cref="ThermodynamicPrimerDesigner.DefaultPrimer3Weights"/>.
/// </summary>
/// <param name="TmGt">PRIMER_WT_TM_GT.</param>
/// <param name="TmLt">PRIMER_WT_TM_LT.</param>
/// <param name="SizeGt">PRIMER_WT_SIZE_GT.</param>
/// <param name="SizeLt">PRIMER_WT_SIZE_LT.</param>
/// <param name="GcGt">PRIMER_WT_GC_PERCENT_GT.</param>
/// <param name="GcLt">PRIMER_WT_GC_PERCENT_LT.</param>
/// <param name="SelfAny">PRIMER_WT_SELF_ANY (alignment mode only).</param>
/// <param name="SelfEnd">PRIMER_WT_SELF_END (alignment mode only).</param>
/// <param name="NumNs">PRIMER_WT_NUM_NS.</param>
/// <param name="SelfAnyTh">PRIMER_WT_SELF_ANY_TH (thermodynamic mode only).</param>
/// <param name="SelfEndTh">PRIMER_WT_SELF_END_TH (thermodynamic mode only).</param>
/// <param name="HairpinTh">PRIMER_WT_HAIRPIN_TH (thermodynamic mode only).</param>
/// <param name="EndStability">PRIMER_WT_END_STABILITY.</param>
/// <param name="ThermodynamicOligoAlignment">PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: false = alignment-score
/// secondary-structure terms (SelfAny/SelfEnd weights), true = thermodynamic terms (…Th weights). Defaults to
/// false for an explicitly constructed weight set; <see cref="ThermodynamicPrimerDesigner.DefaultPrimer3Weights"/> uses
/// Primer3's default (true).</param>
public readonly record struct Primer3PenaltyWeights(
    double TmGt,
    double TmLt,
    double SizeGt,
    double SizeLt,
    double GcGt,
    double GcLt,
    double SelfAny,
    double SelfEnd,
    double NumNs,
    double SelfAnyTh = 0.0,
    double SelfEndTh = 0.0,
    double HairpinTh = 0.0,
    double EndStability = 0.0,
    bool ThermodynamicOligoAlignment = false)
{
    /// <summary>
    /// PRIMER_WT_LIBRARY_MISPRIMING (Primer3 <c>weights.repeat_sim</c>, default 0): × the library mispriming score
    /// (<see cref="Primer3PenaltyInputs.LibraryMispriming"/>). Non-zero requires
    /// <see cref="PrimerParameters.MisprimingLibrary"/> in <see cref="ThermodynamicPrimerDesigner.EvaluatePrimer"/> /
    /// <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_WT_TEMPLATE_MISPRIMING (Primer3 <c>weights.template_mispriming</c>, default 0): × the template mispriming
    /// score (<see cref="Primer3PenaltyInputs.TemplateMispriming"/>) when <see cref="ThermodynamicTemplateAlignment"/> is false.
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_WT_TEMPLATE_MISPRIMING_TH (<c>weights.template_mispriming_th</c>, default 0): thermodynamic template term
    /// (when <see cref="ThermodynamicTemplateAlignment"/> is true) with the 5 °C <c>temp_cutoff</c> rule of the other
    /// thermodynamic terms: s ≥ Tm − 5 → w·(s − (Tm − 6)), else w / (Tm − 4 − s).
    /// </summary>
    public double TemplateMisprimingTh { get; init; }

    /// <summary>
    /// PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT for the template terms (default false = 0, Primer3's default);
    /// <see cref="ThermodynamicPrimerDesigner.EvaluatePrimer"/> / <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> take it from
    /// <see cref="PrimerParameters.ThermodynamicTemplateAlignment"/>.
    /// </summary>
    public bool ThermodynamicTemplateAlignment { get; init; }

    /// <summary>PRIMER_WT_BOUND_GT (Primer3 <c>weights.bound_gt</c>, default 0): × (bound − <see cref="Primer3Optima.OptBound"/>)
    /// when the fraction bound (<see cref="Primer3PenaltyInputs.Bound"/>) is above the optimum.</summary>
    public double BoundGt { get; init; }

    /// <summary>PRIMER_WT_BOUND_LT (<c>weights.bound_lt</c>, default 0): × (<see cref="Primer3Optima.OptBound"/> − bound)
    /// when the fraction bound is below the optimum.</summary>
    public double BoundLt { get; init; }

    /// <summary>PRIMER_WT_POS_PENALTY (<c>weights.pos_penalty</c>, Primer3 default 1): × the position penalty
    /// (<see cref="Primer3PenaltyInputs.PositionPenalty"/>).</summary>
    public double PositionPenalty { get; init; } = ThermodynamicPrimerDesigner.Primer3WeightPositionPenalty;

    /// <summary>
    /// PRIMER_WT_SEQ_QUAL / PRIMER_INTERNAL_WT_SEQ_QUAL (<c>weights.seq_quality</c>, default 0): × (PRIMER_QUALITY_RANGE_MAX −
    /// the oligo's minimum base quality) (<see cref="Primer3PenaltyInputs.SequenceQuality"/>,
    /// <see cref="Primer3PenaltyInputs.QualityRangeMax"/>). Non-zero requires SEQUENCE_QUALITY (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double SequenceQuality { get; init; }

    /// <summary>
    /// PRIMER_WT_END_QUAL / PRIMER_INTERNAL_WT_END_QUAL (<c>weights.end_quality</c>, default 0). Primer3 2.3.1 parses and stores
    /// this weight but its <c>p_obj_fn</c> never reads it, so it has <b>no effect</b> on the penalty — reproduced: it is
    /// accepted and ignored.
    /// </summary>
    public double EndQuality { get; init; }

    /// <summary>
    /// PRIMER_WT_MASK_FAILURE_RATE (<c>weights.failure_rate</c>, default 0): × the primer's predicted failure rate
    /// (<see cref="Primer3PenaltyInputs.MaskFailureRate"/>, <see cref="ThermodynamicPrimerDesigner.CalculateMaskFailureRatePrimer3"/>), which
    /// Primer3 computes only with PRIMER_MASK_TEMPLATE (<see cref="PrimerPairOptions.MaskTemplate"/>; otherwise 0). Left/right
    /// primers only; added after the size terms as in <c>p_obj_fn</c>.
    /// </summary>
    public double MaskFailureRate { get; init; }
}

/// <summary>
/// Optimal parameter values for the Primer3 penalty objective
/// (<c>PRIMER_OPT_TM</c>, <c>PRIMER_OPT_SIZE</c>, <c>PRIMER_OPT_GC_PERCENT</c>).
/// Defaults are <see cref="ThermodynamicPrimerDesigner.DefaultPrimer3Optima"/>.
/// </summary>
public readonly record struct Primer3Optima(
    double OptTm,
    int OptSize,
    double? OptGcPercent)
{
    /// <summary>PRIMER_OPT_BOUND / PRIMER_INTERNAL_OPT_BOUND (% bound; Primer3 default
    /// <see cref="ThermodynamicPrimerDesigner.Primer3OptBound"/> = 97): the optimum of the bound terms.</summary>
    public double OptBound { get; init; } = ThermodynamicPrimerDesigner.Primer3OptBound;
}
