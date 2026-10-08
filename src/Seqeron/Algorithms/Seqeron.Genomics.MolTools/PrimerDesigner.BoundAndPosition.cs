namespace Seqeron.Genomics.MolTools;

// Primer3 fraction bound at the annealing temperature (PRIMER_ANNEALING_TEMP, PRIMER_MIN/MAX/OPT_BOUND,
// PRIMER_WT_BOUND_GT/LT and the PRIMER_INTERNAL_* counterparts) and the position penalty relative to the target
// (PRIMER_INSIDE_PENALTY / PRIMER_OUTSIDE_PENALTY, PRIMER_WT_POS_PENALTY) — primer3-py 2.3.1 libprimer3.c
// calc_and_check_oligo_features / p_obj_fn / compute_position_penalty, oligotm.c oligotm.
public static partial class PrimerDesigner
{
    /// <summary>Primer3 OLIGOTM_ERROR (<c>oligotm.h</c>, −999999.9999): the <c>bound</c> value of an oligo for which
    /// oligotm computed none (PRIMER_ANNEALING_TEMP ≤ 0, more than 36 bases, a non-ACGT base).</summary>
    internal const double Primer3OligoTmError = -999999.9999;

    /// <summary>Primer3's default PRIMER_ANNEALING_TEMP (−10 °C, <c>pr_set_default_global_args_1</c>): a value ≤ 0
    /// disables the fraction-bound computation, checks and penalty terms of the primers.</summary>
    public const double Primer3DefaultAnnealingTemperature = -10.0;

    /// <summary>Primer3's default PRIMER_OPT_BOUND / PRIMER_INTERNAL_OPT_BOUND (97 %).</summary>
    public const double Primer3OptBound = 97.0;

    /// <summary>Primer3's default PRIMER_MIN_BOUND / PRIMER_INTERNAL_MIN_BOUND (−10 %).</summary>
    public const double Primer3MinBound = -10.0;

    /// <summary>Primer3's default PRIMER_MAX_BOUND / PRIMER_INTERNAL_MAX_BOUND (110 %).</summary>
    public const double Primer3MaxBound = 110.0;

    /// <summary>Primer3 PR_INFINITE_POSITION_PENALTY (−1) = the default PRIMER_INSIDE_PENALTY (<c>libprimer3.h</c>
    /// PR_DEFAULT_INSIDE_PENALTY): together with PRIMER_OUTSIDE_PENALTY = 0 primers may not overlap the target.</summary>
    public const double Primer3DefaultInsidePenalty = -1.0;

    /// <summary>Primer3's default PRIMER_OUTSIDE_PENALTY (0, PR_DEFAULT_OUTSIDE_PENALTY).</summary>
    public const double Primer3DefaultOutsidePenalty = 0.0;

    /// <summary>Primer3's default PRIMER_WT_POS_PENALTY (1, <c>pr_set_default_global_args_1</c>).</summary>
    public const double Primer3WeightPositionPenalty = 1.0;

    /// <summary>
    /// Primer3's fraction (%) of an oligo bound to its target at the annealing temperature (<c>oligotm.c</c>
    /// <c>oligotm</c>, PRIMER_ANNEALING_TEMP; primer3-py 2.3.1 PRIMER_LEFT/RIGHT/INTERNAL_n_BOUND), on the scale of
    /// <see cref="CalculateMeltingTemperaturePrimer3"/> (SantaLucia 1998 nearest-neighbour ΔH/ΔS with the SantaLucia
    /// salt correction ΔS += 0.368·(N − 1)·ln([Mon]_eq/1000)):
    /// ΔG = ΔH − (Ta + 273.15)·ΔS, K = exp(−ΔG / (1.987·(Ta + 273.15))),
    /// bound = 100 / (1 + √(1 / ((C/x)·K))), with C the oligo concentration in nM and x = 4·10⁹ (10⁹ for a
    /// self-complementary oligo).
    /// </summary>
    /// <param name="primer">Oligo sequence (case-insensitive).</param>
    /// <param name="annealingTemperature">Annealing temperature Ta, °C (0 &lt; Ta ≤ 100, Primer3 <c>_pr_data_control</c>).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>The bound fraction in percent (0–100), or <c>double.NaN</c> when Primer3 computes none: fewer than
    /// 2 bases, a non-ACGT base, or more than <see cref="Primer3MaxNnTmLength"/> bases (Primer3 then uses
    /// <c>long_seq_tm</c> for the Tm and leaves <c>bound</c> = OLIGOTM_ERROR, which fails any PRIMER_MIN_BOUND ≥ −999999).</returns>
    /// <exception cref="ArgumentOutOfRangeException">Ta ≤ 0, Ta &gt; 100 or not finite, or illegal concentrations.</exception>
    public static double CalculateFractionBoundPrimer3(
        string primer,
        double annealingTemperature,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar)
    {
        if (!(annealingTemperature > 0 && annealingTemperature <= 100))
            throw new ArgumentOutOfRangeException(nameof(annealingTemperature), annealingTemperature,
                "Annealing temperature must be in (0, 100] °C (Primer3 computes no fraction bound at ≤ 0; \"Annealing temperature higher than 100 C\").");
        ValidatePrimer3Conditions(monovalentMillimolar, divalentMillimolar, dntpMillimolar, dnaConcentrationNanomolar,
            nameof(dnaConcentrationNanomolar));
        double monovalentEq = Primer3MonovalentEquivalent(monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        double bound = Primer3SeqTm(primer, dnaConcentrationNanomolar, monovalentEq, annealingTemperature).Bound;
        return bound == Primer3OligoTmError ? double.NaN : bound;
    }

    /// <summary>
    /// Primer3 <c>compute_position_penalty</c> (non-default PRIMER_INSIDE_PENALTY / PRIMER_OUTSIDE_PENALTY, one
    /// SEQUENCE_TARGET): the position penalty of a left (<paramref name="isForward"/>) or right primer relative to the
    /// target [<paramref name="targetStart"/>, <paramref name="targetEnd"/>) (0-based, end exclusive). The primer's
    /// 3′ base is <c>position + length − 1</c> (left) or <c>position</c> (right; the leftmost top-strand base). A left
    /// primer whose 3′ base lies before the target gets <c>outside·(targetStart − 3′ − 1)</c>, one whose 3′ base lies
    /// inside it <c>inside·(3′ − targetStart + 1)</c>, and one whose 3′ base lies past the target end is not allowed
    /// (returns <c>null</c>, Primer3's infinite position penalty / "overlap target"); a right primer mirrors this
    /// (outside: <c>3′ − (targetEnd − 1) − 1</c>, inside: <c>targetEnd − 1 − 3′ + 1</c>, not allowed when 3′ &lt;
    /// targetStart). As in Primer3 the multiplier is applied as given, so the default inside penalty −1 makes an
    /// inside position negative when only the outside penalty is changed.
    /// </summary>
    /// <param name="position">0-based leftmost top-strand base of the primer.</param>
    /// <param name="length">Primer length.</param>
    /// <param name="isForward">True for a left (forward) primer.</param>
    /// <param name="targetStart">0-based inclusive target start.</param>
    /// <param name="targetEnd">0-based exclusive target end.</param>
    /// <param name="insidePenalty">PRIMER_INSIDE_PENALTY.</param>
    /// <param name="outsidePenalty">PRIMER_OUTSIDE_PENALTY.</param>
    /// <returns>The position penalty, or <c>null</c> when the 3′ end lies beyond the target.</returns>
    /// <exception cref="ArgumentOutOfRangeException">length &lt; 1 or targetStart ≥ targetEnd.</exception>
    public static double? CalculatePositionPenaltyPrimer3(
        int position, int length, bool isForward, int targetStart, int targetEnd,
        double insidePenalty, double outsidePenalty)
    {
        if (length < 1)
            throw new ArgumentOutOfRangeException(nameof(length), length, "Primer length must be ≥ 1.");
        if (targetStart >= targetEnd)
            throw new ArgumentOutOfRangeException(nameof(targetEnd), targetEnd, "Target must satisfy targetStart < targetEnd.");
        int targetBegin = targetStart, targetLast = targetEnd - 1;
        int threePrime = isForward ? position + length - 1 : position;
        double penalty;
        bool inside;
        if (isForward)
        {
            if (threePrime > targetLast)
                return null;
            inside = threePrime >= targetBegin;
            penalty = inside ? threePrime - targetBegin + 1 : targetBegin - threePrime - 1;
        }
        else
        {
            if (threePrime < targetBegin)
                return null;
            inside = threePrime <= targetLast;
            penalty = inside ? targetLast - threePrime + 1 : threePrime - targetLast - 1;
        }
        return penalty * (inside ? insidePenalty : outsidePenalty);
    }

    // Primer3 _PR_DEFAULT_POSITION_PENALTIES: exactly the default inside (−1) and outside (0) penalties.
    internal static bool IsDefaultPositionPenalties(double insidePenalty, double outsidePenalty) =>
        insidePenalty == Primer3DefaultInsidePenalty && outsidePenalty == Primer3DefaultOutsidePenalty;

    // _pr_data_control: PRIMER_OPT_BOUND within [MIN, MAX]; PRIMER_ANNEALING_TEMP ≤ 100 when > 0.
    internal static void ValidatePrimer3Bound(double annealingTemperature, double minBound, double maxBound, double optBound,
        bool internalOligo, string paramName)
    {
        if (double.IsNaN(annealingTemperature) || double.IsNaN(minBound) || double.IsNaN(maxBound) || double.IsNaN(optBound))
            throw new ArgumentOutOfRangeException(paramName, "Fraction-bound settings must not be NaN.");
        if (optBound < minBound || optBound > maxBound)
            throw new ArgumentOutOfRangeException(paramName, internalOligo
                ? "Optimum internal oligo fraction binding lower than minimum or higher than maximum (Primer3 _pr_data_control)."
                : "Optimum primer fraction binding lower than minimum or higher than maximum (Primer3 _pr_data_control).");
        if (annealingTemperature > 100.0)
            throw new ArgumentOutOfRangeException(paramName, "Annealing temperature higher than 100 C (Primer3 _pr_data_control).");
    }
}
