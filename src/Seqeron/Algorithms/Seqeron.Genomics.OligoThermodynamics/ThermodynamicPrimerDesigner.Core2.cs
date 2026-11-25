// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

public static partial class ThermodynamicPrimerDesigner
{
    /// <summary>
    /// Default primer design parameters (library conventions — <c>MaxDinucleotideRepeats</c> = 4 is an unsourced
    /// library screen, see <see cref="FindLongestDinucleotideRepeat"/>; the Primer3 3′-end checks PRIMER_GC_CLAMP,
    /// PRIMER_MAX_END_GC and PRIMER_MAX_END_STABILITY are at their Primer3 defaults, i.e. inactive;
    /// <c>Check3PrimeStability</c> is deprecated and has no effect, see <see cref="PrimerParameters"/>).
    /// </summary>
    public static readonly PrimerParameters DefaultParameters = new(
        MinLength: 18,
        MaxLength: 25,
        OptimalLength: 20,
        MinGcContent: 40,
        MaxGcContent: 60,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 4,
        MaxDinucleotideRepeats: 4,
        Avoid3PrimeGC: false,
        Check3PrimeStability: true
    );

    /// <summary>
    /// Primer3's default per-primer settings (<c>libprimer3.cc</c> <c>pr_set_default_global_args_1/_2</c>,
    /// primer3-py 2.3.1 <c>design_primers</c> defaults): PRIMER_MIN/OPT/MAX_SIZE = 18/20/27,
    /// PRIMER_MIN/MAX_GC = 20/80 %, PRIMER_MIN/OPT/MAX_TM = 57/60/63 °C, PRIMER_MAX_POLY_X = 5, no
    /// dinucleotide-repeat limit (Primer3 has none), PRIMER_GC_CLAMP = 0, PRIMER_MAX_END_GC = 5 and
    /// PRIMER_MAX_END_STABILITY = 100 (<see cref="PrimerParameters.GcClamp"/>, <see cref="PrimerParameters.MaxEndGc"/>,
    /// <see cref="PrimerParameters.MaxEndStability"/> left at their Primer3 defaults; 100 cannot be exceeded by an ACGT
    /// 3′ pentamer), thermodynamic structure screen with every limit 47 °C. With <see cref="PrimerPairOptions.Primer3Defaults"/> this makes
    /// <see cref="DesignPrimerPairs"/> reproduce <c>primer3.design_primers</c> run with only
    /// SEQUENCE_TEMPLATE / SEQUENCE_TARGET set. (<see cref="DefaultParameters"/> keeps the library's
    /// stricter conventions: 18–25 nt, 40–60 % GC, poly-X ≤ 4, dinucleotide repeat ≤ 4.)
    /// </summary>
    public static readonly PrimerParameters Primer3DefaultParameters = new(
        MinLength: 18,
        MaxLength: 27,
        OptimalLength: 20,
        MinGcContent: 20,
        MaxGcContent: 80,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 5,
        MaxDinucleotideRepeats: int.MaxValue,
        Avoid3PrimeGC: false,
        Check3PrimeStability: false
    );

    /// <summary>
    /// Primer3 <c>seqtm</c> (<see cref="CalculateMeltingTemperaturePrimer3(string, double, double, double, double)"/>)
    /// with the nearest-neighbour length limit as an argument: oligos of at most
    /// <paramref name="maxNearestNeighborLength"/> nt get the SantaLucia 1998 nearest-neighbour Tm, longer ones
    /// Primer3's <c>long_seq_tm</c> (<c>oligotm.c</c> <c>seqtm</c>: <c>if (len &gt; nn_max_len) long_seq_tm … else
    /// oligotm</c>; primer3-py <c>calc_tm(…, max_nn_length)</c>). Primer3's own primer/probe picker uses
    /// <see cref="Primer3MaxNnTmLength"/> = 36 (the five-argument overload); primer3-py <c>calc_tm</c> defaults to
    /// <c>max_nn_length=60</c>, i.e. <paramref name="maxNearestNeighborLength"/> = 60 reproduces it; e.g. OligoArray-style microarray probes of 50–60 nt need the
    /// nearest-neighbour Tm over the whole oligo (limit ≥ 60).
    /// </summary>
    /// <param name="primer">Oligo sequence (case-insensitive).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <param name="maxNearestNeighborLength">Longest oligo (nt) that gets the nearest-neighbour Tm (Primer3
    /// <c>nn_max_len</c>, ≥ 0).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when the sequence is null/shorter than 2 bases or contains a
    /// non-ACGT character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">For a non-positive DNA concentration, negative ion
    /// concentrations, a zero total monovalent-equivalent concentration or a negative
    /// <paramref name="maxNearestNeighborLength"/>.</exception>
    public static double CalculateMeltingTemperaturePrimer3(
        string primer,
        double dnaConcentrationNanomolar,
        double monovalentMillimolar,
        double divalentMillimolar,
        double dntpMillimolar,
        int maxNearestNeighborLength)
    {
        if (maxNearestNeighborLength < 0)
            throw new ArgumentOutOfRangeException(nameof(maxNearestNeighborLength), maxNearestNeighborLength,
                "The nearest-neighbour length limit must be ≥ 0 nt.");
        if (!(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(nameof(dnaConcentrationNanomolar), dnaConcentrationNanomolar,
                "DNA concentration must be a positive finite value in nM.");
        if (!(monovalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Monovalent cation concentration must be ≥ 0 mM.");
        if (!(divalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(divalentMillimolar), divalentMillimolar,
                "Divalent cation concentration must be ≥ 0 mM.");
        if (!(dntpMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMillimolar), dntpMillimolar,
                "dNTP concentration must be ≥ 0 mM.");

        double monovalentEq = Primer3MonovalentEquivalent(monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        if (!(monovalentEq > 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Total monovalent-equivalent cation concentration must be > 0 mM.");

        return Primer3SeqTm(primer, dnaConcentrationNanomolar, monovalentEq, Primer3DefaultAnnealingTemperature,
            maxNearestNeighborLength).Tm;
    }

    // Primer3 seqtm() / oligotm() (SantaLucia 1998 Tm, SantaLucia salt correction, nn_max_len = maxNnLength, Primer3's
    // MAX_NN_TM_LENGTH = 36 by default) with the
    // fraction bound at the annealing temperature (oligotm.c, PRIMER_ANNEALING_TEMP): Tm is NaN when Primer3 reports
    // OLIGOTM_ERROR (fewer than 2 bases, a non-ACGT base); Bound is Primer3OligoTmError unless annealingTemperature > 0
    // and the nearest-neighbour branch applies (long_seq_tm leaves bound = OLIGOTM_ERROR).
    // Input: validated conditions, monovalentEq = [Mon] + 120·√([Mg²⁺] − [dNTP]) in mM.
    public static (double Tm, double Bound) Primer3SeqTm(
        string? primer, double dnaConcentrationNanomolar, double monovalentEq, double annealingTemperature,
        int maxNnLength = Primer3MaxNnTmLength)
    {
        if (string.IsNullOrEmpty(primer) || primer.Length < 2)
            return (double.NaN, Primer3OligoTmError);

        string seq = primer.ToUpperInvariant();
        int n = seq.Length;
        int gc = 0;
        foreach (char c in seq)
        {
            if (c is 'G' or 'C') gc++;
            else if (c is not ('A' or 'T')) return (double.NaN, Primer3OligoTmError);
        }

        if (n > maxNnLength)
            return (ThermoConstants.CalculateSaltAdjustedTm((double)gc / n, n, monovalentEq / 1000.0), Primer3OligoTmError);

        // oligotm(): integer accumulation, then ΔH = dh·(−100) cal/mol, ΔS = ds·(−0.1) cal/(K·mol).
        bool symmetric = OligoSymmetry.IsSymmetric(seq);
        int dh = 0, ds = symmetric ? Primer3SymmetryDs : 0;
        foreach (char end in new[] { seq[0], seq[^1] })
        {
            if (end is 'A' or 'T') { dh += Primer3TerminalAtDh; ds += Primer3TerminalAtDs; }
            else { dh += Primer3TerminalGcDh; ds += Primer3TerminalGcDs; }
        }
        for (int i = 0; i < n - 1; i++)
        {
            int x = BaseIndex(seq[i]), y = BaseIndex(seq[i + 1]);
            dh += Primer3SantaLucia1998Dh[x, y];
            ds += Primer3SantaLucia1998Ds[x, y];
        }

        double deltaH = dh * -100.0;
        double deltaS = ds * -0.1;
        deltaS += SantaLuciaEntropySaltCoefficient * (n - 1) * Math.Log(monovalentEq / 1000.0);
        // Equation A (self-complementary, C_T/1) or Equation B (C_T/4) of oligotm.c.
        double strandDivisor = symmetric ? 1000000000.0 : 4000000000.0;
        double tm = deltaH / (deltaS + Primer3GasConstant * Math.Log(dnaConcentrationNanomolar / strandDivisor)) - KelvinOffset;
        double bound = Primer3OligoTmError;
        if (annealingTemperature > 0.0)
        {
            // oligotm.c (santalucia salt correction): ddG = ΔH − (Ta + 273.15)·ΔS; Ka = exp(−ddG / (R·(Ta + 273.15)));
            // bound = 100 / (1 + √(1 / ((C/x)·Ka))), x = 1e9 (Equation A) or 4e9 (Equation B).
            double ddG = deltaH - (annealingTemperature + KelvinOffset) * deltaS;
            double ka = Math.Exp(-ddG / (Primer3GasConstant * (annealingTemperature + KelvinOffset)));
            bound = (1 / (1 + Math.Sqrt(1 / ((dnaConcentrationNanomolar / strandDivisor) * ka)))) * 100;
        }
        return (tm, bound);

        static int BaseIndex(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };
    }

    // dpal.c _dpal_long_nopath_maxgap1_local (score only, max gap 1, cells floored at 0). With a finite
    // stopAbove the scan returns as soon as the running optimum exceeds it (the value is then a lower bound).
    // matrix: a dpal substitution matrix (MatrixSsm); null = DpalSsm (non-ACGT scored as N).
    private static int DpalLocalFast(string x, string y, double stopAbove = double.PositiveInfinity, int[]? matrix = null)
    {
        int Ssm(char a, char b) => matrix is null ? DpalSsm(a, b) : MatrixSsm(matrix, a, b);
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[ylen];
        var s1 = new int[ylen];
        var s2 = new int[ylen];
        int smax = 0, score, a;

        // Row 0.
        for (int j = 0; j < ylen; j++)
        {
            score = Ssm(x[0], y[j]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s0[j] = score;
        }
        // Row 1 (for |X| = 1 the C code reads X[1] = NUL, whose ssm row is INT_MIN: every cell is 0).
        if (xlen == 1)
            return smax;
        score = Ssm(x[1], y[0]);
        if (score < 0) score = 0;
        else if (score > smax) smax = score;
        s1[0] = score;
        for (int j = 1; j < ylen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += Ssm(x[1], y[j]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s1[j] = score;
        }

        for (int i = 2; i < xlen; i++)
        {
            score = Ssm(x[i], y[0]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s2[0] = score;
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += Ssm(x[i], y[1]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s2[1] = score;
            for (int j = 2; j < ylen; j++)
            {
                score = s0[j - 1];
                if ((a = s1[j - 2]) > score) score = a;
                score += gap;
                if ((a = s1[j - 1]) > score) score = a;
                score += Ssm(x[i], y[j]);
                if (score < 0) score = 0;
                else if (score > smax) smax = score;
                s2[j] = score;
            }
            if (smax > stopAbove)
                return smax;
            (s0, s1, s2) = (s1, s2, s0);
        }
        return smax;
    }

    // ---- Primer3 weighted penalty objective (PRIMER-TM-001) -------------------
    // Reproduces the per-primer objective function `p_obj_fn` (left/right primer
    // branch) from Primer3's reference source `libprimer3.cc`, with the documented
    // default weights and optima. Lower penalty = better primer, exactly as Primer3.
    // Source: Primer3 source libprimer3.cc p_obj_fn / pr_set_default_global_args_2;
    //         Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE";
    //         Untergasser et al. (2012) NAR 40(15):e115; Koressaar & Remm (2007).

    /// <summary>
    /// Primer3 default per-primer objective weights, taken verbatim from
    /// <c>pr_set_default_global_args_2</c> in Primer3's <c>libprimer3.cc</c>:
    /// length/Tm weights = 1; GC, self-complementarity and N weights = 0.
    /// Source: Primer3 source (branch main), function pr_set_default_global_args_2.
    /// </summary>
    public static readonly Primer3PenaltyWeights DefaultPrimer3Weights = new(
        TmGt: 1.0,           // PRIMER_WT_TM_GT  (libprimer3.cc: weights.temp_gt = 1)
        TmLt: 1.0,           // PRIMER_WT_TM_LT  (weights.temp_lt = 1)
        SizeGt: 1.0,         // PRIMER_WT_SIZE_GT (weights.length_gt = 1)
        SizeLt: 1.0,         // PRIMER_WT_SIZE_LT (weights.length_lt = 1)
        GcGt: 0.0,           // PRIMER_WT_GC_PERCENT_GT (weights.gc_content_gt = 0)
        GcLt: 0.0,           // PRIMER_WT_GC_PERCENT_LT (weights.gc_content_lt = 0)
        SelfAny: 0.0,        // PRIMER_WT_SELF_ANY (weights.compl_any = 0)
        SelfEnd: 0.0,        // PRIMER_WT_SELF_END (weights.compl_end = 0)
        NumNs: 0.0,          // PRIMER_WT_NUM_NS  (weights.num_ns = 0)
        SelfAnyTh: 0.0,      // PRIMER_WT_SELF_ANY_TH (weights.compl_any_th = 0)
        SelfEndTh: 0.0,      // PRIMER_WT_SELF_END_TH (weights.compl_end_th = 0)
        HairpinTh: 0.0,      // PRIMER_WT_HAIRPIN_TH  (weights.hairpin_th = 0)
        EndStability: 0.0,   // PRIMER_WT_END_STABILITY (weights.end_stability = 0)
        ThermodynamicOligoAlignment: true // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1 (pr_set_default_global_args_2)
    );

    /// <summary>
    /// Primer3 default per-primer optima. OPT_TM = 60 °C and OPT_SIZE = 20 bases are
    /// from <c>libprimer3.c</c> (opt_tm = 60.0, opt_size = 20); OPT_GC_PERCENT is undefined (null) as in Primer3's code
    /// (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>; the manual's 50 is not applied): a non-zero GC weight
    /// then needs an explicit optimum (Primer3 <c>_pr_data_control</c>).
    /// <see cref="EvaluatePrimer"/> / <see cref="DesignPrimers"/> take the GC optimum from
    /// <see cref="PrimerParameters.OptimalGcPercent"/>.
    /// </summary>
    public static readonly Primer3Optima DefaultPrimer3Optima = new(
        OptTm: 60.0,         // PRIMER_OPT_TM (libprimer3.cc: opt_tm = 60.0) °C
        OptSize: 20,         // PRIMER_OPT_SIZE (libprimer3.cc: opt_size = 20) bases
        OptGcPercent: null   // PRIMER_OPT_GC_PERCENT (libprimer3.c: PR_UNDEFINED_INT_OPT) %
    );

    /// <summary>
    /// Computes the Primer3 per-primer penalty (objective function value) for a single
    /// primer, faithfully reproducing the left/right-primer branch of Primer3's
    /// <c>p_obj_fn</c>. The penalty is the weighted sum of one-sided deviations of Tm,
    /// length and GC% from their optima, plus the weighted secondary-structure terms
    /// (alignment-score mode: <c>compl_any</c>/<c>compl_end</c>; thermodynamic mode, Primer3's
    /// default: <c>compl_any_th</c>/<c>compl_end_th</c>/<c>hairpin_th</c> with the fixed 5 °C
    /// <c>temp_cutoff</c>), the number-of-Ns term and the 3'-end-stability term. Each term is
    /// added only when its weight is non-zero (and, for Tm/GC/size, the deviation has the
    /// matching sign), so the result is always ≥ 0; <b>lower is better</b>, exactly as Primer3
    /// sorts candidates. Cross-checked against primer3-py 2.3.1 <c>design_primers</c>
    /// PRIMER_LEFT/RIGHT_n_PENALTY in both alignment modes.
    /// The library mispriming term <c>repeat_sim</c> is PRIMER_WT_LIBRARY_MISPRIMING
    /// (<see cref="Primer3PenaltyWeights.LibraryMispriming"/>) × <see cref="Primer3PenaltyInputs.LibraryMispriming"/>.
    /// The fraction-bound terms PRIMER_WT_BOUND_GT / _LT (<see cref="Primer3PenaltyWeights.BoundGt"/> /
    /// <see cref="Primer3PenaltyWeights.BoundLt"/> around <see cref="Primer3Optima.OptBound"/>) apply when
    /// <see cref="Primer3PenaltyInputs.Bound"/> is set (Primer3: PRIMER_ANNEALING_TEMP &gt; 0); the position term is
    /// PRIMER_WT_POS_PENALTY × <see cref="Primer3PenaltyInputs.PositionPenalty"/> (PRIMER_INSIDE/OUTSIDE_PENALTY).
    /// The sequence-quality term is PRIMER_WT_SEQ_QUAL (<see cref="Primer3PenaltyWeights.SequenceQuality"/>) ×
    /// (<see cref="Primer3PenaltyInputs.QualityRangeMax"/> − <see cref="Primer3PenaltyInputs.SequenceQuality"/>).
    /// The masker term is PRIMER_WT_MASK_FAILURE_RATE (<see cref="Primer3PenaltyWeights.MaskFailureRate"/>) ×
    /// <see cref="Primer3PenaltyInputs.MaskFailureRate"/>. <para>The template mispriming terms (PRIMER_WT_TEMPLATE_MISPRIMING / _TH,
    /// <see cref="Primer3PenaltyWeights.TemplateMispriming"/> / <see cref="Primer3PenaltyWeights.TemplateMisprimingTh"/>)
    /// use <see cref="Primer3PenaltyInputs.TemplateMispriming"/>.</para>
    /// </summary>
    /// <param name="inputs">Measured primer properties (Tm in °C, length in bases, GC in
    /// percent 0–100, self/3' local-alignment scores, count of N bases).</param>
    /// <param name="weights">Objective weights; defaults to <see cref="DefaultPrimer3Weights"/>.</param>
    /// <param name="optima">Parameter optima; defaults to <see cref="DefaultPrimer3Optima"/>.</param>
    /// <returns>The Primer3 objective-function value (penalty); 0 means every term is at its optimum.</returns>
    public static double CalculatePrimer3Penalty(
        Primer3PenaltyInputs inputs,
        Primer3PenaltyWeights? weights = null,
        Primer3Optima? optima = null)
    {
        var w = weights ?? DefaultPrimer3Weights;
        var o = optima ?? DefaultPrimer3Optima;

        double sum = 0.0;

        // Tm term: one-sided, separate _gt / _lt weights (p_obj_fn temp_gt / temp_lt).
        if (w.TmGt != 0 && inputs.Tm > o.OptTm)
            sum += w.TmGt * (inputs.Tm - o.OptTm);
        if (w.TmLt != 0 && inputs.Tm < o.OptTm)
            sum += w.TmLt * (o.OptTm - inputs.Tm);

        // Fraction-bound terms (bound_gt / bound_lt around PRIMER_OPT_BOUND); for primers only when
        // PRIMER_ANNEALING_TEMP > 0 (inputs.Bound is then set).
        if (inputs.Bound is { } bound)
        {
            if (w.BoundGt != 0 && bound > o.OptBound)
                sum += w.BoundGt * (bound - o.OptBound);
            if (w.BoundLt != 0 && bound < o.OptBound)
                sum += w.BoundLt * (o.OptBound - bound);
        }

        // GC% term (gc_content is a percentage 0–100 in libprimer3.c); a GC weight needs a defined optimum
        // (_pr_data_control: "Primer GC content is part of objective function while optimum gc_content is not defined").
        if (w.GcGt != 0 || w.GcLt != 0)
        {
            double optGc = o.OptGcPercent ?? throw new ArgumentException(PrimerGcOptimumUndefinedMessage, nameof(optima));
            if (w.GcGt != 0 && inputs.GcPercent > optGc)
                sum += w.GcGt * (inputs.GcPercent - optGc);
            if (w.GcLt != 0 && inputs.GcPercent < optGc)
                sum += w.GcLt * (optGc - inputs.GcPercent);
        }

        // Length/size term (p_obj_fn length_lt / length_gt).
        if (w.SizeLt != 0 && inputs.Length < o.OptSize)
            sum += w.SizeLt * (o.OptSize - inputs.Length);
        if (w.SizeGt != 0 && inputs.Length > o.OptSize)
            sum += w.SizeGt * (inputs.Length - o.OptSize);

        // Masker failure-rate term (failure_rate; left/right primers, PRIMER_MASK_TEMPLATE).
        if (w.MaskFailureRate != 0)
            sum += w.MaskFailureRate * inputs.MaskFailureRate;

        // Secondary-structure terms: p_obj_fn switches on thermodynamic_oligo_alignment.
        if (!w.ThermodynamicOligoAlignment)
        {
            // Mode 0: local-alignment scores, linear weights (compl_any, compl_end).
            if (w.SelfAny != 0)
                sum += w.SelfAny * inputs.SelfAny;
            if (w.SelfEnd != 0)
                sum += w.SelfEnd * inputs.SelfEnd;
        }
        else
        {
            // Mode 1 (Primer3 default): SelfAny / SelfEnd / HairpinTh are structure Tm values (°C)
            // (compl_any_th, compl_end_th, hairpin_th).
            sum += ThermodynamicStructurePenalty(w.SelfAnyTh, inputs.Tm, inputs.SelfAny);
            sum += ThermodynamicStructurePenalty(w.SelfEndTh, inputs.Tm, inputs.SelfEnd);
            sum += ThermodynamicStructurePenalty(w.HairpinTh, inputs.Tm, inputs.HairpinTh);
        }

        // Number-of-Ns term (num_ns).
        if (w.NumNs != 0)
            sum += w.NumNs * inputs.NumNs;

        // Library mispriming term (repeat_sim): weight · repeat_sim.score[repeat_sim.max].
        if (w.LibraryMispriming != 0)
            sum += w.LibraryMispriming * inputs.LibraryMispriming;

        // Position term (pos_penalty): weight · position penalty relative to the target.
        if (w.PositionPenalty != 0)
            sum += w.PositionPenalty * inputs.PositionPenalty;

        // 3'-end stability term (end_stability): weight · ΔG magnitude (kcal/mol, as Primer3 reports it).
        if (w.EndStability != 0)
            sum += w.EndStability * inputs.EndStability;

        // Sequence quality term (seq_quality): weight · (PRIMER_QUALITY_RANGE_MAX − h->seq_quality); 0 without quality
        // data (seq_quality = range max). Internal oligos: the last term of the OT_INTL branch (the end-stability,
        // position and template inputs are 0 there). PRIMER_WT_END_QUAL is never read by p_obj_fn.
        if (w.SequenceQuality != 0)
            sum += w.SequenceQuality * (inputs.QualityRangeMax - (inputs.SequenceQuality ?? inputs.QualityRangeMax));

        // Template mispriming terms (after seq_quality): alignment mode linear
        // (weights.template_mispriming), thermodynamic mode the temp_cutoff rule (weights.template_mispriming_th).
        if (!w.ThermodynamicTemplateAlignment)
        {
            if (w.TemplateMispriming != 0)
                sum += w.TemplateMispriming * inputs.TemplateMispriming;
        }
        else
            sum += ThermodynamicStructurePenalty(w.TemplateMisprimingTh, inputs.Tm, inputs.TemplateMispriming);

        return sum;
    }

    // SantaLucia (1998) Table 2 unified NN parameters exactly as tabulated in Primer3 oligotm.c
    // (SantaLucia_1998_dH in units of −100 cal/mol, SantaLucia_1998_dS in units of −0.1 cal/(K·mol);
    // index order A, C, G, T). Primer3 sums these as integers, so primers with the same
    // nearest-neighbour multiset get bit-identical Tm values (this matters for tie-breaking in
    // DesignPrimers). Differs from the SantaLucia & Hicks (2004) set (ThermoConstants NnParameterSet.SantaLuciaHicks2004) in
    // AA/TT (−7.9/−22.2 vs −7.6/−21.3).
    private static readonly int[,] Primer3SantaLucia1998Dh =
    {
        { 79, 84, 78, 72 },
        { 85, 80, 106, 78 },
        { 82, 98, 80, 84 },
        { 72, 82, 85, 79 }
    };

    private static readonly int[,] Primer3SantaLucia1998Ds =
    {
        { 222, 224, 210, 204 },
        { 227, 199, 272, 210 },
        { 222, 244, 199, 224 },
        { 213, 222, 227, 222 }
    };

    // SantaLucia (1998) Table 2 initiation per duplex end, in the same integer units
    // (oligotm.c "Terminal penalty"): terminal A·T ΔH +2.3 / ΔS +4.1 → (−23, −41);
    // terminal G·C ΔH +0.1 / ΔS −2.8 → (−1, +28); symmetry ΔS −1.4 → +14.
    private const int Primer3TerminalAtDh = -23, Primer3TerminalAtDs = -41;

    private const int Primer3TerminalGcDh = -1, Primer3TerminalGcDs = 28;

    private const int Primer3SymmetryDs = 14;

    // oligotm.c constants: R = 1.987 cal/(K·mol); divalent_to_monovalent factor 120 (mM units).
    private const double Primer3GasConstant = 1.987;

    /// <summary>
    /// Primer3's default primer melting temperature (°C) — Primer3 <c>seqtm</c> with
    /// PRIMER_TM_FORMULA = SantaLucia 1998 and PRIMER_SALT_CORRECTIONS = SantaLucia 1998:
    /// <para>
    /// [Mon]_eq = [Mon] + 120·√([Mg²⁺] − [dNTP]) (mM; von Ahsen 2001, 0 when [Mg²⁺] ≤ [dNTP]);
    /// ΔS = ΔS°(1 M) + 0.368·(N − 1)·ln([Mon]_eq / 1000);
    /// Tm = ΔH / (ΔS + R·ln(C/4)) − 273.15, with C/1 for a self-complementary oligo,
    /// R = 1.987, C = DNA concentration in mol/L;
    /// for N &gt; <see cref="Primer3MaxNnTmLength"/>: Tm = 81.5 + 16.6·log10([Mon]_eq/1000) + 41·GC/N − 600/N.
    /// </para>
    /// Defaults are Primer3's PRIMER_DNA_CONC = 50 nM, PRIMER_SALT_MONOVALENT = 50 mM,
    /// PRIMER_SALT_DIVALENT = 1.5 mM, PRIMER_DNTP_CONC = 0.6 mM. Equals primer3-py
    /// <c>calc_tm(…, max_nn_length=36)</c> — the MAX_NN_TM_LENGTH libprimer3 uses when picking primers and probes
    /// (e.g. AGCTAGCTAGCTAGCTAGCT → 58.101 °C). primer3-py's <c>calc_tm</c> itself defaults to
    /// <c>max_nn_length=60</c>, so for 37–60-nt oligos it differs from this method (39-mer
    /// CTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTAC: 73.049 °C here and with <c>max_nn_length=36</c>, 70.425 °C with the
    /// <c>calc_tm</c> default); use
    /// <see cref="CalculateMeltingTemperaturePrimer3(string, double, double, double, double, int)"/> with
    /// <c>maxNearestNeighborLength</c> = 60 to reproduce <c>calc_tm</c>'s default.
    /// </summary>
    /// <param name="primer">Primer sequence (case-insensitive).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when the sequence is null/shorter than 2 bases or contains a
    /// non-ACGT character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">For a non-positive DNA concentration, negative ion
    /// concentrations, or a zero total monovalent-equivalent concentration.</exception>
    public static double CalculateMeltingTemperaturePrimer3(
        string primer,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar) =>
        CalculateMeltingTemperaturePrimer3(primer, dnaConcentrationNanomolar, monovalentMillimolar,
            divalentMillimolar, dntpMillimolar, Primer3MaxNnTmLength);

    // Primer3 oligo_compl decision (self_any > maxSelfAny, else self_end > maxSelfEnd) without computing the full
    // scores: the dpal recurrences stop as soon as the running optimum exceeds the limit (the decision is identical
    // to comparing the full scores, since the optimum never decreases). Used by the probe fallback screen, where
    // probes of several hundred nt make the full O(n²) scores the dominant cost.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static bool ExceedsPrimer3SelfComplementarity(string oligo, double maxSelfAny, double maxSelfEnd,
        out bool selfAnyExceeded)
    {
        selfAnyExceeded = false;
        if (string.IsNullOrEmpty(oligo))
            return false;
        string x = oligo.ToUpperInvariant();
        string y = ReverseComplementPrimer3(x);
        if (y.Length < 3 ? y.Length > maxSelfAny : DpalLocalFast(x, y, maxSelfAny * 100.0) / 100.0 > maxSelfAny)
            return selfAnyExceeded = true;
        return DpalGlobalEndScore(x, y) > maxSelfEnd;
    }

    // Primer3 align() with DPAL_LOCAL: a second sequence shorter than 3 "maxes out" the score to its length
    // (not divided by 100); otherwise dpal score / 100, floored at 0.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static double DpalLocalScore(string x, string y)
    {
        if (y.Length < 3)
            return y.Length;
        int smax = DpalLocalFast(x, y);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    /// <summary>
    /// Primer3 alignment-mode self-complementarity <c>self_any</c> (<c>libprimer3.cc</c> <c>oligo_compl</c>:
    /// <c>align(oligo, revcomp(oligo), DPAL_LOCAL)</c>; Primer3 <c>PRIMER_LEFT/RIGHT/INTERNAL_n_SELF_ANY</c> with
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; default limits <c>PRIMER_MAX_SELF_ANY = 8.00</c>,
    /// <c>PRIMER_INTERNAL_MAX_SELF_ANY = 12.00</c>): the best <c>dpal</c> local alignment of the oligo with its
    /// reverse complement, scored +1 per complementary base pair, −1 per mismatch, −0.25 against N (any non-ACGT
    /// character is scored as N), −2 per single-base gap (max gap 1), never below 0 — a line-by-line port of
    /// dpal.c <c>_dpal_long_nopath_maxgap1_local</c> (the routine Primer3 runs in DPM_FAST mode) with Primer3's
    /// <c>align</c> rule that a second sequence shorter than 3 scores its length. No length limit (dpal's
    /// score-only routine is linear in memory). Bit-exact to dpal.c / primer3-py 2.3.1.
    /// </summary>
    /// <param name="oligo">Oligo (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>self_any</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerSelfAnyComplementarity(string oligo)
    {
        if (string.IsNullOrEmpty(oligo))
            return 0;
        string p = oligo.ToUpperInvariant();
        return DpalLocalScore(p, ReverseComplementPrimer3(p));
    }

    /// <summary>
    /// Primer3 alignment-mode pair complementarity <c>compl_any</c> (<c>libprimer3.cc</c> <c>characterize_pair</c>,
    /// <c>PRIMER_PAIR_COMPL_ANY</c> with <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; default limit
    /// <c>PRIMER_PAIR_MAX_COMPL_ANY = 8.00</c>): <c>align(s1, s2, DPAL_LOCAL)</c> where s1 is the left primer and s2 the
    /// right primer's top-strand copy, i.e. the dpal local alignment of <paramref name="leftPrimer"/> with
    /// revcomp(<paramref name="rightPrimer"/>) (scoring as <see cref="CalculatePrimerSelfAnyComplementarity"/>).
    /// Primer3 evaluates this one orientation only.
    /// </summary>
    /// <param name="leftPrimer">Left (forward) primer, 5′→3′, case-insensitive.</param>
    /// <param name="rightPrimer">Right (reverse) primer, 5′→3′, case-insensitive.</param>
    /// <returns>The Primer3 <c>compl_any</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerDimerAnyComplementarity(string leftPrimer, string rightPrimer)
    {
        if (string.IsNullOrEmpty(leftPrimer) || string.IsNullOrEmpty(rightPrimer))
            return 0;
        return DpalLocalScore(leftPrimer.ToUpperInvariant(), ReverseComplementPrimer3(rightPrimer.ToUpperInvariant()));
    }
}
