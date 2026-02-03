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
    /// Base temperature constant in Marmur-Doty formula.
    /// Tm = 64.9 + 41*(GC - 16.4)/N
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
    /// Base temperature for salt-adjusted Tm formula.
    /// Tm = 81.5 + 16.6*log10([Na+]) + 41*(%GC) - 600/length
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
    /// Calculates salt correction factor for Tm.
    /// </summary>
    /// <param name="naConcentrationMM">Na+ concentration in mM.</param>
    /// <returns>Salt correction in °C.</returns>
    public static double CalculateSaltCorrection(double naConcentrationMM) =>
        SaltCoefficient * Math.Log10(naConcentrationMM / 1000.0);

    #endregion
}
