// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

public static partial class ThermodynamicProbeDesigner
{
    // Primer3 MAX_PRIMER_LENGTH (oligo length limit of the picker and of seqtm's nearest-neighbour branch).
    private const int Primer3MaxOligoLength = 36;

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
    /// <see cref="ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3"/>) at the given conditions — compare it with a
    /// primer Tm computed on the same scale (e.g. <see cref="ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3"/>);
    /// it is <c>NaN</c> when the probe has fewer than 2 bases or a non-ACGT base, and the Tm gate then fails
    /// whenever a primer Tm is supplied.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="probeSequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">DNA or monovalent concentration not &gt; 0, divalent or dNTP concentration
    /// negative, or any of them non-finite (Primer3 <c>_pr_data_control</c>).</exception>
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
        // Primer3 _pr_data_control condition limits, each reported under its own parameter (heavy tier, HEAVY-2: an
        // infinite salt surfaced as ThermoConstants' "naConcentration").
        if (!(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(nameof(dnaConcentrationNanomolar), dnaConcentrationNanomolar, "PRIMER_INTERNAL_DNA_CONC must be > 0 and finite.");
        if (!(monovalentMillimolar > 0) || double.IsInfinity(monovalentMillimolar))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar, "PRIMER_INTERNAL_SALT_MONOVALENT must be > 0 and finite.");
        if (!(divalentMillimolar >= 0) || double.IsInfinity(divalentMillimolar))
            throw new ArgumentOutOfRangeException(nameof(divalentMillimolar), divalentMillimolar, "PRIMER_INTERNAL_SALT_DIVALENT must be ≥ 0 and finite.");
        if (!(dntpMillimolar >= 0) || double.IsInfinity(dntpMillimolar))
            throw new ArgumentOutOfRangeException(nameof(dntpMillimolar), dntpMillimolar, "PRIMER_INTERNAL_DNTP_CONC must be ≥ 0 and finite.");

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
        double tm = ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(
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
        /// (<see cref="ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity"/> /
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
        /// library mishybridization score (<see cref="ThermodynamicPrimerDesigner.CalculateLibraryMishyb"/>,
        /// <see cref="Primer3Probe.LibraryMishyb"/>) and is rejected when any entry's weighted score exceeds
        /// <see cref="MaxLibraryMishyb"/> (Primer3 OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ, a "five-prime problem" that also ends
        /// the 5′ extension of that 3′ end). Null = no library (Primer3's default). Not part of the JSON form (the MCP tools
        /// take the library as a separate name → sequence argument).
        /// </summary>
        [System.Text.Json.Serialization.JsonIgnore]
        public PrimerMisprimingLibrary? MishybLibrary { get; init; }

        /// <summary>
        /// PRIMER_INTERNAL_MAX_LIBRARY_MISHYB (Primer3 default <see cref="ThermodynamicPrimerDesigner.Primer3InternalMaxLibraryMishyb"/>,
        /// 12.00): maximum weighted library score of one oligo; Primer3 compares with the value truncated to a C
        /// <c>short</c>. In alignment mode (<see cref="ThermodynamicOligoAlignment"/> = false) it must not exceed 32767.
        /// </summary>
        public double MaxLibraryMishyb { get; init; } = ThermodynamicPrimerDesigner.Primer3InternalMaxLibraryMishyb;

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
        /// PRIMER_ANNEALING_TEMP (°C; Primer3 default <see cref="ThermodynamicPrimerDesigner.Primer3DefaultAnnealingTemperature"/> = −10
        /// = off; at most 100). When &gt; 0 every oligo gets its fraction bound at this temperature
        /// (<see cref="ThermodynamicPrimerDesigner.CalculateFractionBoundPrimer3"/> at the settings' conditions,
        /// <see cref="Primer3Probe.Bound"/>) and is rejected outside [<see cref="MinBound"/>, <see cref="MaxBound"/>].
        /// Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.AnnealingTemperature"/> is used instead.
        /// </summary>
        public double AnnealingTemperature { get; init; } = ThermodynamicPrimerDesigner.Primer3DefaultAnnealingTemperature;

        /// <summary>PRIMER_INTERNAL_MIN_BOUND (% bound; default −10), checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
        public double MinBound { get; init; } = ThermodynamicPrimerDesigner.Primer3MinBound;

        /// <summary>PRIMER_INTERNAL_MAX_BOUND (% bound; default 110), checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
        public double MaxBound { get; init; } = ThermodynamicPrimerDesigner.Primer3MaxBound;

        /// <summary>PRIMER_INTERNAL_OPT_BOUND (% bound; default 97; must lie in [<see cref="MinBound"/>, <see cref="MaxBound"/>]).</summary>
        public double OptBound { get; init; } = ThermodynamicPrimerDesigner.Primer3OptBound;

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
        /// (<see cref="ThermodynamicPrimerDesigner.CalculateSequenceQualityPrimer3"/>, <see cref="Primer3Probe.MinSequenceQuality"/>) is lower is
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
        public int QualityRangeMin { get; init; } = ThermodynamicPrimerDesigner.Primer3QualityRangeMin;

        /// <summary>PRIMER_QUALITY_RANGE_MAX (default 100). Primer3 has one global setting: for the internal oligo of a primer pair
        /// <see cref="PrimerParameters.QualityRangeMax"/> is used instead.</summary>
        public int QualityRangeMax { get; init; } = ThermodynamicPrimerDesigner.Primer3QualityRangeMax;
    }
}
