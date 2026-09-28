namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Constants for thermodynamic calculations (melting temperature, etc.).
/// Centralizes magic numbers used in Tm formulas across the codebase.
/// </summary>
public static class ThermoConstants
{
    #region Wallace Rule (Short Oligos < 14 bp)

    /// <summary>
    /// Contribution of A-T base pairs to Tm in Wallace rule.
    /// Tm = 2*(A+T) + 4*(G+C)
    /// </summary>
    public const int WallaceAtContribution = 2;

    /// <summary>
    /// Contribution of G-C base pairs to Tm in Wallace rule.
    /// </summary>
    public const int WallaceGcContribution = 4;

    /// <summary>
    /// Maximum primer length for Wallace rule application.
    /// </summary>
    public const int WallaceMaxLength = 14;

    #endregion

    #region Basic GC Formula (Marmur-Doty)

    /// <summary>
    /// Base temperature constant in the "basic" GC formula for oligos of 14 nt or more,
    /// Tm = 64.9 + 41*(GC - 16.4)/N, as published in OligoCalc (Kibbe 2007, NAR 35:W43)
    /// and customarily attributed to Marmur &amp; Doty (1962). The formula has a fixed
    /// ionic condition built in: OligoCalc states it assumes 50 nM primer, 50 mM Na+ and pH 7.0
    /// (algebraically 64.9 + 0.41·%GC − 672.4/N, i.e. the 81.5 + 0.41·%GC − 675/N family with the
    /// 16.6·log10[Na+] term already folded into the constant). It must therefore NOT be combined
    /// with an absolute 16.6·log10[Na+] correction (that would count the salt term twice).
    /// </summary>
    public const double MarmurDotyBase = 64.9;

    /// <summary>
    /// GC coefficient in Marmur-Doty formula.
    /// </summary>
    public const double MarmurDotyGcCoefficient = 41.0;

    /// <summary>
    /// GC offset constant in Marmur-Doty formula.
    /// </summary>
    public const double MarmurDotyGcOffset = 16.4;

    #endregion

    #region Salt-Adjusted Formula

    /// <summary>
    /// Base temperature for the salt-adjusted GC formula
    /// Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/N  ([Na+] in M).
    /// This is Biopython <c>MeltingTemp.Tm_GC</c> valueset 7 ("used by Primer3Plus to calculate
    /// the product Tm"; salt term = Schildkraut &amp; Lifson 1965, Biopython salt_correction method 1).
    /// </summary>
    public const double SaltAdjustedBase = 81.5;

    /// <summary>
    /// Salt concentration coefficient.
    /// </summary>
    public const double SaltCoefficient = 16.6;

    /// <summary>
    /// GC percentage coefficient in salt-adjusted formula.
    /// </summary>
    public const double SaltAdjustedGcCoefficient = 41.0;

    /// <summary>
    /// Length correction factor in salt-adjusted formula.
    /// </summary>
    public const double SaltAdjustedLengthFactor = 600.0;

    /// <summary>
    /// Default Na+ concentration in M (0.05 = 50mM).
    /// </summary>
    public const double DefaultNaConcentration = 0.05;

    /// <summary>
    /// OligoCalc salt-adjusted Tm for oligos of 14 nt or more (Kibbe 2007, NAR 35:W43):
    /// Tm = 100.5 + 41*(G+C)/N − 820/N + 16.6*log10([Na+]) ([Na+] in M). Base constant.
    /// </summary>
    public const double OligoCalcSaltAdjustedBase = 100.5;

    /// <summary>
    /// Length term (820/N) of the OligoCalc salt-adjusted Tm for oligos of 14 nt or more.
    /// </summary>
    public const double OligoCalcSaltAdjustedLengthFactor = 820.0;

    /// <summary>
    /// Reference [Na+] (M) at which the basic Wallace rule is defined in OligoCalc's salt-adjusted
    /// short-oligo formula: Tm = 2(A+T) + 4(G+C) − 16.6*log10(0.050) + 16.6*log10([Na+]).
    /// </summary>
    public const double WallaceReferenceNaMolar = 0.050;

    #endregion

    #region Calculation Methods

    /// <summary>
    /// Calculates Tm using Wallace rule for short oligonucleotides.
    /// </summary>
    /// <param name="countAT">Number of A and T nucleotides.</param>
    /// <param name="countGC">Number of G and C nucleotides.</param>
    /// <returns>Melting temperature in °C.</returns>
    public static double CalculateWallaceTm(int countAT, int countGC) =>
        WallaceAtContribution * countAT + WallaceGcContribution * countGC;

    /// <summary>
    /// Calculates Tm using Marmur-Doty formula for longer primers.
    /// </summary>
    /// <param name="gcCount">Number of G and C nucleotides.</param>
    /// <param name="length">Total sequence length.</param>
    /// <returns>Melting temperature in °C.</returns>
    public static double CalculateMarmurDotyTm(int gcCount, int length)
    {
        if (length == 0) return 0;
        return MarmurDotyBase + MarmurDotyGcCoefficient * (gcCount - MarmurDotyGcOffset) / length;
    }

    /// <summary>
    /// Calculates salt-adjusted Tm.
    /// </summary>
    /// <param name="gcFraction">GC content as fraction (0-1).</param>
    /// <param name="length">Sequence length.</param>
    /// <param name="naConcentration">Na+ concentration in M (default 0.05 = 50mM).</param>
    /// <returns>Melting temperature in °C.</returns>
    public static double CalculateSaltAdjustedTm(double gcFraction, int length, double naConcentration = DefaultNaConcentration)
    {
        if (length == 0) return 0;
        return SaltAdjustedBase + SaltCoefficient * Math.Log10(naConcentration) +
               SaltAdjustedGcCoefficient * gcFraction - SaltAdjustedLengthFactor / length;
    }

    /// <summary>
    /// OligoCalc salt-adjusted Tm (Kibbe 2007, NAR 35:W43), the salt-aware counterpart of the
    /// basic Wallace / 64.9-formula pair:
    /// <list type="bullet">
    /// <item>fewer than <see cref="WallaceMaxLength"/> bases: Tm = 2(A+T) + 4(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+]);</item>
    /// <item>otherwise: Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+]).</item>
    /// </list>
    /// [Na+] in mol/L. Cross-check: OligoCalc reports 78 °C (salt adjusted) and 67.6 °C (basic)
    /// for the 39-mer GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCC at 50 mM Na+.
    /// </summary>
    /// <param name="countAT">Number of A and T nucleotides.</param>
    /// <param name="countGC">Number of G and C nucleotides.</param>
    /// <param name="naConcentrationMolar">Na+ concentration in M (must be &gt; 0).</param>
    /// <returns>Melting temperature in °C (0 when there are no counted bases).</returns>
    public static double CalculateOligoCalcSaltAdjustedTm(int countAT, int countGC, double naConcentrationMolar)
    {
        int length = countAT + countGC;
        if (length == 0) return 0;
        double saltTerm = SaltCoefficient * Math.Log10(naConcentrationMolar);
        if (length < WallaceMaxLength)
            return CalculateWallaceTm(countAT, countGC) - SaltCoefficient * Math.Log10(WallaceReferenceNaMolar) + saltTerm;
        return OligoCalcSaltAdjustedBase + MarmurDotyGcCoefficient * countGC / length
               - OligoCalcSaltAdjustedLengthFactor / length + saltTerm;
    }

    /// <summary>
    /// Absolute Schildkraut &amp; Lifson (1965) salt term 16.6·log10([Na+]) with [Na+] given in mM
    /// (Biopython <c>salt_correction(method=1)</c>; −21.60 °C at 50 mM). Only valid for Tm formulas whose
    /// constant is referenced to 1 M Na+ (e.g. the 81.5 family); do not add it to the Wallace rule or to
    /// the 64.9 basic formula, which already assume 50 mM Na+.
    /// </summary>
    /// <param name="naConcentrationMM">Na+ concentration in mM.</param>
    /// <returns>Salt correction in °C.</returns>
    public static double CalculateSaltCorrection(double naConcentrationMM) =>
        SaltCoefficient * Math.Log10(naConcentrationMM / 1000.0);

    #endregion
}
