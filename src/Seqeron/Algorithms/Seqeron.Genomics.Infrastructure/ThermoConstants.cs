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

    /// <summary>
    /// The canonical "basic" oligonucleotide melting temperature (OligoCalc basic Tm, Kibbe 2007,
    /// NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for fewer than
    /// <see cref="WallaceMaxLength"/> counted bases, otherwise Tm = 64.9 + 41·(G+C − 16.4)/N
    /// (<see cref="CalculateMarmurDotyTm"/>), at OligoCalc's fixed conditions (50 nM primer, 50 mM Na+).
    /// <para>
    /// Counting: A, C, G, T and U, case-insensitive; U is read as T, as Biopython
    /// <c>MeltingTemp._check</c> back-transcribes RNA before <c>Tm_Wallace</c>/<c>Tm_GC</c>. Every other
    /// character (N, IUPAC codes, gaps, whitespace) is ignored, and N — the length that selects the
    /// formula and divides the GC term — is the number of counted bases (OligoCalc divides by
    /// N = wA+xT+yG+zC). For ≥ 14 counted bases the GC formula is ≥ 64.9 − 672.4/14 &gt; 16.8 °C, so the
    /// result is never negative.
    /// </para>
    /// </summary>
    /// <param name="sequence">Oligonucleotide (DNA or RNA).</param>
    /// <returns>Tm in °C; 0 for null/empty input or when no A/C/G/T/U base is present.</returns>
    public static double CalculateBasicTm(string? sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        var (at, gc) = CountBasicTmBases(sequence);
        int length = at + gc;
        if (length == 0)
            return 0;

        return length < WallaceMaxLength
            ? CalculateWallaceTm(at, gc)
            : CalculateMarmurDotyTm(gc, length);
    }

    /// <summary>
    /// Counts the bases used by the basic / OligoCalc Tm formulas: (A+T+U, G+C), case-insensitive,
    /// U read as T (Biopython <c>_check</c> back-transcription); every other character is ignored.
    /// </summary>
    /// <param name="sequence">Oligonucleotide (DNA or RNA).</param>
    /// <returns>(A+T+U count, G+C count); (0, 0) for null/empty input.</returns>
    public static (int CountAT, int CountGC) CountBasicTmBases(string? sequence)
    {
        int at = 0, gc = 0;
        if (string.IsNullOrEmpty(sequence))
            return (0, 0);
        foreach (char ch in sequence)
        {
            switch (ch)
            {
                case 'A': case 'a': case 'T': case 't': case 'U': case 'u': at++; break;
                case 'G': case 'g': case 'C': case 'c': gc++; break;
            }
        }
        return (at, gc);
    }

    #endregion

    #region Nearest-neighbour Tm core (Biopython Bio.SeqUtils.MeltingTemp.Tm_NN)

    // One canonical DNA/DNA nearest-neighbour engine, a line-by-line port of Biopython 1.88
    // Bio.SeqUtils.MeltingTemp (Tm_NN, salt_correction, _check, the DNA_NN1..4, DNA_IMM1, DNA_TMM1 and
    // DNA_DE1 tables). It lives in Infrastructure so that both MolTools (PrimerDesigner) and Analysis
    // (SequenceStatistics) can call it.
    //
    // Tables (kcal/mol, cal/(K·mol), 1 M NaCl) — verbatim from Biopython 1.88 and cross-checked
    // entry by entry against the MELTING 5 data files (Dumousseau et al. 2012; mirrored in
    // aravind-j/rmelting inst/extdata/Data): Breslauer1986nn, Sugimoto1996nn, AllawiSantalucia1997nn,
    // Santalucia2004nn, AllawiSantaluciaPeyret1997_1998_1999mm (+ …tanmm tandem G·T),
    // Santalucia2005inomn (inosine), Bommarito2000de — all identical except MELTING's
    // Santalucia2004nn TA/AT ΔS° −20.4, a MELTING transcription error (SantaLucia & Hicks 2004
    // Table 1 / SantaLucia 1998 Table 2: TA/AT −7.2 / −21.3, ΔG°37 −0.58). The terminal-mismatch
    // table (SantaLucia & Peyret 2001, patent WO 01/94611) agrees with Primer3's tstack2 parameters
    // in 47/48 entries (AC/TC ΔS° +0.5 here, −0.5 in Primer3); Biopython's value is kept.

    private sealed class NnTable
    {
        public required (double H, double S) Init { get; init; }
        public required (double H, double S) InitAT { get; init; }
        public required (double H, double S) InitGC { get; init; }
        public required (double H, double S) InitOneGC { get; init; }
        public required (double H, double S) InitAllAT { get; init; }
        public required (double H, double S) Init5TA { get; init; }
        public required (double H, double S) Sym { get; init; }
        public required Dictionary<string, (double H, double S)> Stacks { get; init; }
    }

    // Biopython DNA_NN1 — Breslauer et al. (1986) PNAS 83:3746.
    private static readonly NnTable DnaNn1 = new()
    {
        Init = (0, 0), InitAT = (0, 0), InitGC = (0, 0),
        InitOneGC = (0, -16.8), InitAllAT = (0, -20.1), Init5TA = (0, 0), Sym = (0, -1.3),
        Stacks = new()
        {
            ["AA/TT"] = (-9.1, -24.0), ["AT/TA"] = (-8.6, -23.9), ["TA/AT"] = (-6.0, -16.9),
            ["CA/GT"] = (-5.8, -12.9), ["GT/CA"] = (-6.5, -17.3), ["CT/GA"] = (-7.8, -20.8),
            ["GA/CT"] = (-5.6, -13.5), ["CG/GC"] = (-11.9, -27.8), ["GC/CG"] = (-11.1, -26.7),
            ["GG/CC"] = (-11.0, -26.6)
        }
    };

    // Biopython DNA_NN2 — Sugimoto et al. (1996) NAR 24:4501.
    private static readonly NnTable DnaNn2 = new()
    {
        Init = (0.6, -9.0), InitAT = (0, 0), InitGC = (0, 0),
        InitOneGC = (0, 0), InitAllAT = (0, 0), Init5TA = (0, 0), Sym = (0, -1.4),
        Stacks = new()
        {
            ["AA/TT"] = (-8.0, -21.9), ["AT/TA"] = (-5.6, -15.2), ["TA/AT"] = (-6.6, -18.4),
            ["CA/GT"] = (-8.2, -21.0), ["GT/CA"] = (-9.4, -25.5), ["CT/GA"] = (-6.6, -16.4),
            ["GA/CT"] = (-8.8, -23.5), ["CG/GC"] = (-11.8, -29.0), ["GC/CG"] = (-10.5, -26.4),
            ["GG/CC"] = (-10.9, -28.4)
        }
    };

    // Biopython DNA_NN3 — Allawi & SantaLucia (1997) Biochemistry 36:10581 (= SantaLucia 1998 PNAS
    // 95:1460 unified set with per-terminal-pair initiation).
    private static readonly NnTable DnaNn3 = new()
    {
        Init = (0, 0), InitAT = (2.3, 4.1), InitGC = (0.1, -2.8),
        InitOneGC = (0, 0), InitAllAT = (0, 0), Init5TA = (0, 0), Sym = (0, -1.4),
        Stacks = new()
        {
            ["AA/TT"] = (-7.9, -22.2), ["AT/TA"] = (-7.2, -20.4), ["TA/AT"] = (-7.2, -21.3),
            ["CA/GT"] = (-8.5, -22.7), ["GT/CA"] = (-8.4, -22.4), ["CT/GA"] = (-7.8, -21.0),
            ["GA/CT"] = (-8.2, -22.2), ["CG/GC"] = (-10.6, -27.2), ["GC/CG"] = (-9.8, -24.4),
            ["GG/CC"] = (-8.0, -19.9)
        }
    };

    // Biopython DNA_NN4 — SantaLucia & Hicks (2004) Annu Rev Biophys Biomol Struct 33:415, Table 1
    // (duplex initiation + terminal A·T penalty; AA/TT −7.6/−21.3 — differs from 1998's −7.9/−22.2).
    private static readonly NnTable DnaNn4 = new()
    {
        Init = (0.2, -5.7), InitAT = (2.2, 6.9), InitGC = (0, 0),
        InitOneGC = (0, 0), InitAllAT = (0, 0), Init5TA = (0, 0), Sym = (0, -1.4),
        Stacks = new()
        {
            ["AA/TT"] = (-7.6, -21.3), ["AT/TA"] = (-7.2, -20.4), ["TA/AT"] = (-7.2, -21.3),
            ["CA/GT"] = (-8.5, -22.7), ["GT/CA"] = (-8.4, -22.4), ["CT/GA"] = (-7.8, -21.0),
            ["GA/CT"] = (-8.2, -22.2), ["CG/GC"] = (-10.6, -27.2), ["GC/CG"] = (-9.8, -24.4),
            ["GG/CC"] = (-8.0, -19.9)
        }
    };

    // Biopython DNA_IMM1 — internal single mismatches and inosine: Allawi & SantaLucia (1997)
    // Biochemistry 36:10581 (G·T, incl. tandem G·T); (1998) Biochemistry 37:9435 (G·A), 37:2170 (C·T),
    // NAR 26:2694 (A·C); Peyret et al. (1999) Biochemistry 38:3468 (A·A, C·C, G·G, T·T); Watkins &
    // SantaLucia (2005) NAR 33:6258 (inosine). Key "top/bottom", bottom written 3'→5'.
    private static readonly Dictionary<string, (double H, double S)> DnaImm1 = new()
    {
        ["AG/TT"] = (1.0, 0.9), ["AT/TG"] = (-2.5, -8.3), ["CG/GT"] = (-4.1, -11.7),
        ["CT/GG"] = (-2.8, -8.0), ["GG/CT"] = (3.3, 10.4), ["GG/TT"] = (5.8, 16.3),
        ["GT/CG"] = (-4.4, -12.3), ["GT/TG"] = (4.1, 9.5), ["TG/AT"] = (-0.1, -1.7),
        ["TG/GT"] = (-1.4, -6.2), ["TT/AG"] = (-1.3, -5.3), ["AA/TG"] = (-0.6, -2.3),
        ["AG/TA"] = (-0.7, -2.3), ["CA/GG"] = (-0.7, -2.3), ["CG/GA"] = (-4.0, -13.2),
        ["GA/CG"] = (-0.6, -1.0), ["GG/CA"] = (0.5, 3.2), ["TA/AG"] = (0.7, 0.7),
        ["TG/AA"] = (3.0, 7.4),
        ["AC/TT"] = (0.7, 0.2), ["AT/TC"] = (-1.2, -6.2), ["CC/GT"] = (-0.8, -4.5),
        ["CT/GC"] = (-1.5, -6.1), ["GC/CT"] = (2.3, 5.4), ["GT/CC"] = (5.2, 13.5),
        ["TC/AT"] = (1.2, 0.7), ["TT/AC"] = (1.0, 0.7),
        ["AA/TC"] = (2.3, 4.6), ["AC/TA"] = (5.3, 14.6), ["CA/GC"] = (1.9, 3.7),
        ["CC/GA"] = (0.6, -0.6), ["GA/CC"] = (5.2, 14.2), ["GC/CA"] = (-0.7, -3.8),
        ["TA/AC"] = (3.4, 8.0), ["TC/AA"] = (7.6, 20.2),
        ["AA/TA"] = (1.2, 1.7), ["CA/GA"] = (-0.9, -4.2), ["GA/CA"] = (-2.9, -9.8),
        ["TA/AA"] = (4.7, 12.9), ["AC/TC"] = (0.0, -4.4), ["CC/GC"] = (-1.5, -7.2),
        ["GC/CC"] = (3.6, 8.9), ["TC/AC"] = (6.1, 16.4), ["AG/TG"] = (-3.1, -9.5),
        ["CG/GG"] = (-4.9, -15.3), ["GG/CG"] = (-6.0, -15.8), ["TG/AG"] = (1.6, 3.6),
        ["AT/TT"] = (-2.7, -10.8), ["CT/GT"] = (-5.0, -15.8), ["GT/CT"] = (-2.2, -8.4),
        ["TT/AT"] = (0.2, -1.5),
        ["AI/TC"] = (-8.9, -25.5), ["TI/AC"] = (-5.9, -17.4), ["AC/TI"] = (-8.8, -25.4),
        ["TC/AI"] = (-4.9, -13.9), ["CI/GC"] = (-5.4, -13.7), ["GI/CC"] = (-6.8, -19.1),
        ["CC/GI"] = (-8.3, -23.8), ["GC/CI"] = (-5.0, -12.6),
        ["AI/TA"] = (-8.3, -25.0), ["TI/AA"] = (-3.4, -11.2), ["AA/TI"] = (-0.7, -2.6),
        ["TA/AI"] = (-1.3, -4.6), ["CI/GA"] = (2.6, 8.9), ["GI/CA"] = (-7.8, -21.1),
        ["CA/GI"] = (-7.0, -20.0), ["GA/CI"] = (-7.6, -20.2),
        ["AI/TT"] = (0.49, -0.7), ["TI/AT"] = (-6.5, -22.0), ["AT/TI"] = (-5.6, -18.7),
        ["TT/AI"] = (-0.8, -4.3), ["CI/GT"] = (-1.0, -2.4), ["GI/CT"] = (-3.5, -10.6),
        ["CT/GI"] = (0.1, -1.0), ["GT/CI"] = (-4.3, -12.1),
        ["AI/TG"] = (-4.9, -15.8), ["TI/AG"] = (-1.9, -8.5), ["AG/TI"] = (0.1, -1.8),
        ["TG/AI"] = (1.0, 1.0), ["CI/GG"] = (7.1, 21.3), ["GI/CG"] = (-1.1, -3.2),
        ["CG/GI"] = (5.8, 16.9), ["GG/CI"] = (-7.6, -22.0),
        ["AI/TI"] = (-3.3, -11.9), ["TI/AI"] = (0.1, -2.3), ["CI/GI"] = (1.3, 3.0),
        ["GI/CI"] = (-0.5, -1.3)
    };

    // Biopython DNA_TMM1 — terminal mismatches, SantaLucia & Peyret (2001) WO 01/94611.
    private static readonly Dictionary<string, (double H, double S)> DnaTmm1 = new()
    {
        ["AA/TA"] = (-3.1, -7.8), ["TA/AA"] = (-2.5, -6.3), ["CA/GA"] = (-4.3, -10.7),
        ["GA/CA"] = (-8.0, -22.5),
        ["AC/TC"] = (-0.1, 0.5), ["TC/AC"] = (-0.7, -1.3), ["CC/GC"] = (-2.1, -5.1),
        ["GC/CC"] = (-3.9, -10.6),
        ["AG/TG"] = (-1.1, -2.1), ["TG/AG"] = (-1.1, -2.7), ["CG/GG"] = (-3.8, -9.5),
        ["GG/CG"] = (-0.7, -19.2),
        ["AT/TT"] = (-2.4, -6.5), ["TT/AT"] = (-3.2, -8.9), ["CT/GT"] = (-6.1, -16.9),
        ["GT/CT"] = (-7.4, -21.2),
        ["AA/TC"] = (-1.6, -4.0), ["AC/TA"] = (-1.8, -3.8), ["CA/GC"] = (-2.6, -5.9),
        ["CC/GA"] = (-2.7, -6.0), ["GA/CC"] = (-5.0, -13.8), ["GC/CA"] = (-3.2, -7.1),
        ["TA/AC"] = (-2.3, -5.9), ["TC/AA"] = (-2.7, -7.0),
        ["AC/TT"] = (-0.9, -1.7), ["AT/TC"] = (-2.3, -6.3), ["CC/GT"] = (-3.2, -8.0),
        ["CT/GC"] = (-3.9, -10.6), ["GC/CT"] = (-4.9, -13.5), ["GT/CC"] = (-3.0, -7.8),
        ["TC/AT"] = (-2.5, -6.3), ["TT/AC"] = (-0.7, -1.2),
        ["AA/TG"] = (-1.9, -4.4), ["AG/TA"] = (-2.5, -5.9), ["CA/GG"] = (-3.9, -9.6),
        ["CG/GA"] = (-6.0, -15.5), ["GA/CG"] = (-4.3, -11.1), ["GG/CA"] = (-4.6, -11.4),
        ["TA/AG"] = (-2.0, -4.7), ["TG/AA"] = (-2.4, -5.8),
        ["AG/TT"] = (-3.2, -8.7), ["AT/TG"] = (-3.5, -9.4), ["CG/GT"] = (-3.8, -9.0),
        ["CT/GG"] = (-6.6, -18.7), ["GG/CT"] = (-5.7, -15.9), ["GT/CG"] = (-5.9, -16.1),
        ["TG/AT"] = (-3.9, -10.5), ["TT/AG"] = (-3.6, -9.8)
    };

    // Biopython DNA_DE1 — single dangling ends, Bommarito, Peyret & SantaLucia (2000) NAR 28:1929.
    private static readonly Dictionary<string, (double H, double S)> DnaDe1 = new()
    {
        ["AA/.T"] = (0.2, 2.3), ["AC/.G"] = (-6.3, -17.1), ["AG/.C"] = (-3.7, -10.0),
        ["AT/.A"] = (-2.9, -7.6), ["CA/.T"] = (0.6, 3.3), ["CC/.G"] = (-4.4, -12.6),
        ["CG/.C"] = (-4.0, -11.9), ["CT/.A"] = (-4.1, -13.0), ["GA/.T"] = (-1.1, -1.6),
        ["GC/.G"] = (-5.1, -14.0), ["GG/.C"] = (-3.9, -10.9), ["GT/.A"] = (-4.2, -15.0),
        ["TA/.T"] = (-6.9, -20.0), ["TC/.G"] = (-4.0, -10.9), ["TG/.C"] = (-4.9, -13.8),
        ["TT/.A"] = (-0.2, -0.5),
        [".A/AT"] = (-0.7, -0.8), [".C/AG"] = (-2.1, -3.9), [".G/AC"] = (-5.9, -16.5),
        [".T/AA"] = (-0.5, -1.1), [".A/CT"] = (4.4, 14.9), [".C/CG"] = (-0.2, -0.1),
        [".G/CC"] = (-2.6, -7.4), [".T/CA"] = (4.7, 14.2), [".A/GT"] = (-1.6, -3.6),
        [".C/GG"] = (-3.9, -11.2), [".G/GC"] = (-3.2, -10.4), [".T/GA"] = (-4.1, -13.1),
        [".A/TT"] = (2.9, 10.4), [".C/TG"] = (-4.4, -13.1), [".G/TC"] = (-5.2, -15.0),
        [".T/TA"] = (-3.8, -12.6)
    };

    /// <summary>Gas constant used by Biopython <c>Tm_NN</c> (R = 1.987 cal/(K·mol); SantaLucia 1998).</summary>
    public const double NnGasConstantBiopython = 1.987;

    /// <summary>Gas constant of SantaLucia &amp; Hicks (2004) Eq. 3 (R = 1.9872 cal/(K·mol)).</summary>
    public const double NnGasConstantSantaLuciaHicks2004 = 1.9872;

    /// <summary>Biopython <c>salt_correction</c> method 7: dNTP·Mg²⁺ association constant Ka (M⁻¹).</summary>
    public const double NnDntpMagnesiumAssociationConstant = 3e4;

    private static NnTable GetNnTable(NnParameterSet set) => set switch
    {
        NnParameterSet.Breslauer1986 => DnaNn1,
        NnParameterSet.Sugimoto1996 => DnaNn2,
        NnParameterSet.AllawiSantaLucia1997 => DnaNn3,
        NnParameterSet.SantaLuciaHicks2004 => DnaNn4,
        _ => throw new ArgumentOutOfRangeException(nameof(set), set, "Unknown nearest-neighbour parameter set.")
    };

    /// <summary>
    /// Biopython <c>MeltingTemp._check(seq, "Tm_NN")</c>: removes whitespace, upper-cases, back-transcribes
    /// (U → T) and keeps only the bases the NN tables can score (A, C, G, T and inosine I).
    /// </summary>
    /// <param name="sequence">Raw sequence.</param>
    /// <returns>The normalised sequence (empty for null input).</returns>
    public static string NormalizeForNearestNeighbor(string? sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return string.Empty;
        var sb = new System.Text.StringBuilder(sequence.Length);
        foreach (char ch in sequence)
        {
            char b = char.ToUpperInvariant(ch);
            if (b == 'U') b = 'T';
            if (b is 'A' or 'C' or 'G' or 'T' or 'I')
                sb.Append(b);
        }
        return sb.ToString();
    }

    // Biopython Seq.complement (Bio.Data.IUPACData.ambiguous_dna_complement, U → A, case kept), same
    // left-to-right order; every other character (I, N, S, W, X, '.', …) maps to itself.
    private static string NnComplement(string seq) =>
        string.Create(seq.Length, seq, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
            {
                char ch = src[i];
                char up = char.ToUpperInvariant(ch);
                char co = up switch
                {
                    'A' => 'T', 'T' => 'A', 'U' => 'A', 'C' => 'G', 'G' => 'C',
                    'M' => 'K', 'K' => 'M', 'R' => 'Y', 'Y' => 'R',
                    'V' => 'B', 'B' => 'V', 'H' => 'D', 'D' => 'H',
                    _ => up
                };
                dest[i] = char.IsLower(ch) ? char.ToLowerInvariant(co) : co;
            }
        });

    // Biopython SeqUtils.gc_fraction(seq, "ignore"): (G + C + S, either case) / len(seq).
    private static double NnGcFraction(string seq)
    {
        if (seq.Length == 0) return 0;
        int gc = 0;
        foreach (char c in seq)
            if (c is 'G' or 'C' or 'S' or 'g' or 'c' or 's') gc++;
        return (double)gc / seq.Length;
    }

    private static string ReverseString(string s)
    {
        var a = s.ToCharArray();
        Array.Reverse(a);
        return new string(a);
    }

    /// <summary>
    /// Biopython <c>MeltingTemp.salt_correction</c> (1.88), all seven methods:
    /// <list type="bullet">
    /// <item>[Mon] = [Na⁺] + [K⁺] + [Tris]/2 (mM); for methods 1–6, when any of K⁺/Mg²⁺/Tris/dNTPs is &gt; 0
    /// and [dNTPs] &lt; [Mg²⁺], the von Ahsen et al. (2001) equivalent [Mon] += 120·√([Mg²⁺] − [dNTPs]);</item>
    /// <item>1: 16.6·log10[Mon] (Schildkraut &amp; Lifson 1965); 2: 16.6·log10([Mon]/(1 + 0.7[Mon])) (Wetmur 1991);
    /// 3: 12.5·log10[Mon] (SantaLucia 1996); 4: 11.7·log10[Mon] (SantaLucia 1998) — added to Tm (°C);</item>
    /// <item>5: 0.368·(N − 1)·ln[Mon] (SantaLucia 1998) — added to ΔS° (cal/(K·mol));</item>
    /// <item>6: (4.29·f(GC) − 3.95)·10⁻⁵·ln[Mon] + 9.40·10⁻⁶·ln²[Mon] (Owczarzy et al. 2004) — added to 1/Tm (K⁻¹);</item>
    /// <item>7: Owczarzy et al. (2008) divalent model (free Mg²⁺ after dNTP binding, Ka = 3·10⁴; R = √[Mg²⁺]/[Mon]
    /// &lt; 0.22 → method 6 form; 0.22 ≤ R &lt; 6 → a, d, g re-parameterised by [Mon]) — added to 1/Tm.</item>
    /// </list>
    /// Concentrations in mM; ions enter the formulas in mol/L. f(GC) and N are taken from
    /// <paramref name="sequence"/> as given (Biopython passes the checked sequence).
    /// </summary>
    /// <param name="method">Correction method; <see cref="NnSaltCorrection.None"/> returns 0.</param>
    /// <param name="sodium">[Na⁺], mM.</param>
    /// <param name="potassium">[K⁺], mM.</param>
    /// <param name="tris">[Tris], mM.</param>
    /// <param name="magnesium">[Mg²⁺], mM.</param>
    /// <param name="dntps">[dNTPs], mM.</param>
    /// <param name="sequence">Sequence (required for methods 5–7).</param>
    /// <returns>The correction term (unit depends on the method, see above).</returns>
    /// <exception cref="ArgumentException">Methods 5–7 without a sequence; a zero total ion concentration
    /// (methods 1–6); a ln of a non-positive concentration; method 7 on a 1-nt sequence.</exception>
    public static double CalculateNnSaltCorrection(
        NnSaltCorrection method,
        double sodium = 0, double potassium = 0, double tris = 0, double magnesium = 0, double dntps = 0,
        string? sequence = null)
    {
        int m = (int)method;
        if (m is < 0 or > 7)
            throw new ArgumentOutOfRangeException(nameof(method), method, "Allowed salt-correction methods are 0–7.");
        if (m is 5 or 6 or 7 && string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence is missing (needed for the GC content or sequence length).", nameof(sequence));
        if (m == 0)
            return 0;

        double monMm = sodium + potassium + tris / 2.0;
        double mg = magnesium * 1e-3;
        if ((potassium + magnesium + tris + dntps) > 0 && m != 7 && dntps < magnesium)
            monMm += 120 * Math.Sqrt(magnesium - dntps);
        double mon = monMm * 1e-3;
        if (m <= 6 && mon == 0)
            throw new ArgumentException("Total ion concentration of zero is not allowed in this method.", nameof(sodium));

        double corr = m switch
        {
            1 => 16.6 * Math.Log10(Positive(mon)),
            2 => 16.6 * Math.Log10(Positive(mon / (1.0 + 0.7 * mon))),
            3 => 12.5 * Math.Log10(Positive(mon)),
            4 => 11.7 * Math.Log10(Positive(mon)),
            5 => 0.368 * (sequence!.Length - 1) * Math.Log(Positive(mon)),
            6 => (4.29 * NnGcFraction(sequence!) - 3.95) * 1e-5 * Math.Log(Positive(mon))
                 + 9.40e-6 * Math.Pow(Math.Log(Positive(mon)), 2),
            _ => 0
        };
        if (m != 7)
            return corr;

        double a = 3.92, b = -0.911, c = 6.26, d = 1.42, e = -48.2, f = 52.5, g = 8.31;
        if (dntps > 0)
        {
            double dn = dntps * 1e-3;
            double ka = NnDntpMagnesiumAssociationConstant;
            mg = (-(ka * dn - ka * mg + 1.0)
                  + Math.Sqrt(Math.Pow(ka * dn - ka * mg + 1.0, 2) + 4.0 * ka * mg)) / (2.0 * ka);
        }
        if (monMm > 0)
        {
            double r = Math.Sqrt(mg) / mon;
            double lnMon = Math.Log(Positive(mon));
            if (r < 0.22)
                return (4.29 * NnGcFraction(sequence!) - 3.95) * 1e-5 * lnMon + 9.40e-6 * Math.Pow(lnMon, 2);
            if (r < 6.0)
            {
                a = 3.92 * (0.843 - 0.352 * Math.Sqrt(mon) * lnMon);
                d = 1.42 * (1.279 - 4.03e-3 * lnMon - 8.03e-3 * Math.Pow(lnMon, 2));
                g = 8.31 * (0.486 - 0.258 * lnMon + 5.25e-3 * Math.Pow(lnMon, 3));
            }
        }
        if (sequence!.Length < 2)
            throw new ArgumentException("The Owczarzy (2008) correction needs at least 2 bases.", nameof(sequence));
        double lnMg = Math.Log(Positive(mg));
        return (a + b * lnMg + NnGcFraction(sequence) * (c + d * lnMg)
                + (1 / (2.0 * (sequence.Length - 1))) * (e + f * lnMg + g * Math.Pow(lnMg, 2))) * 1e-5;

        static double Positive(double x) => x > 0
            ? x
            : throw new ArgumentException("A salt concentration entering a logarithm must be > 0.");
    }

    /// <summary>
    /// Duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl, before any salt correction) by the
    /// nearest-neighbour model, exactly as Biopython <c>Tm_NN</c> accumulates them: dangling ends
    /// (<c>de_table</c> DNA_DE1), terminal mismatches (<c>tmm_table</c> DNA_TMM1), the initiation terms of
    /// the chosen table (general, allA/T or oneG/C, 5'-T/3'-A, per terminal A·T / G·C pair of
    /// <paramref name="sequence"/>'s ends), the stacks (internal mismatch / inosine table DNA_IMM1 first,
    /// then the Watson–Crick table, each key also tried reversed) and, for a self-complementary duplex,
    /// the symmetry term. See <see cref="CalculateNearestNeighborDuplex"/> for the parameters.
    /// </summary>
    /// <returns>(ΔH°, ΔS°).</returns>
    public static (double DeltaH, double DeltaS) CalculateNearestNeighborThermodynamics(
        string sequence,
        string? complement = null,
        int shift = 0,
        NnParameterSet parameterSet = NnParameterSet.AllawiSantaLucia1997,
        bool selfComplementary = false,
        bool check = true,
        bool strict = true)
    {
        var (dh, ds, _) = NnSum(sequence, complement, shift, parameterSet, selfComplementary, check, strict);
        return (dh, ds);
    }

    /// <summary>
    /// The canonical nearest-neighbour duplex model — a line-by-line port of Biopython 1.88
    /// <c>Bio.SeqUtils.MeltingTemp.Tm_NN</c> (parameter names and defaults kept):
    /// Tm = 1000·ΔH° / (ΔS° + R·ln k) − 273.15 with k = (dnac1 − dnac2/2)·10⁻⁹ M, or k = dnac1·10⁻⁹ M for a
    /// self-complementary duplex; salt per <see cref="CalculateNnSaltCorrection"/> (method 5 corrects ΔS°,
    /// methods 1–4 add to Tm, methods 6–7 add to 1/Tm). Mismatched duplexes (<paramref name="complement"/>)
    /// and dangling ends (<paramref name="shift"/> / unequal lengths) are corrected automatically.
    /// </summary>
    /// <param name="sequence">Primer/probe 5'→3'.</param>
    /// <param name="complement">Template/target written 3'→5' (complement direction, base i under base i of
    /// <paramref name="sequence"/> at shift 0); <c>null</c> = the perfect complement.</param>
    /// <param name="shift">Offset of <paramref name="sequence"/> on <paramref name="complement"/>
    /// (&gt; 0: complement overhangs on the left; &lt; 0: sequence overhangs on the left); a shorter strand is
    /// padded on the right. Overhangs longer than one base are trimmed (only one dangling base is scored).</param>
    /// <param name="parameterSet">Watson–Crick NN table (default DNA_NN3, as Biopython).</param>
    /// <param name="dnac1">Concentration of the higher-concentrated strand, nM (default 25).</param>
    /// <param name="dnac2">Concentration of the lower-concentrated strand, nM (default 25). Equal dnac1 = dnac2
    /// = C_T/2 gives the SantaLucia k = C_T/4.</param>
    /// <param name="selfComplementary">Self-complementary duplex (adds the symmetry term, k = dnac1); the caller's
    /// statement, not inferred (as in Biopython).</param>
    /// <param name="sodium">[Na⁺], mM (default 50).</param>
    /// <param name="potassium">[K⁺], mM.</param>
    /// <param name="tris">[Tris], mM.</param>
    /// <param name="magnesium">[Mg²⁺], mM.</param>
    /// <param name="dntps">[dNTPs], mM.</param>
    /// <param name="saltCorrection">Salt-correction method (default 5, as Biopython).</param>
    /// <param name="check">Normalise both strands with <see cref="NormalizeForNearestNeighbor"/> (default true).</param>
    /// <param name="strict">Throw on a neighbour pair with no parameter (default true); otherwise skip it
    /// (Biopython only warns).</param>
    /// <param name="gasConstant">R in cal/(K·mol) (default Biopython's 1.987).</param>
    /// <returns>ΔH° and ΔS° (1 M NaCl), ΔS° after the method-5 correction (= ΔS° for other methods) and Tm (°C).</returns>
    /// <exception cref="ArgumentException">Empty sequence after normalisation, a missing parameter with
    /// <paramref name="strict"/>, a non-positive k, or an invalid salt input.</exception>
    public static NnDuplexResult CalculateNearestNeighborDuplex(
        string sequence,
        string? complement = null,
        int shift = 0,
        NnParameterSet parameterSet = NnParameterSet.AllawiSantaLucia1997,
        double dnac1 = 25,
        double dnac2 = 25,
        bool selfComplementary = false,
        double sodium = 50,
        double potassium = 0,
        double tris = 0,
        double magnesium = 0,
        double dntps = 0,
        NnSaltCorrection saltCorrection = NnSaltCorrection.SantaLucia1998Entropy,
        bool check = true,
        bool strict = true,
        double gasConstant = NnGasConstantBiopython)
    {
        var (deltaH, deltaS, seq) = NnSum(sequence, complement, shift, parameterSet, selfComplementary, check, strict);
        return CalculateNearestNeighborTmFromThermodynamics(deltaH, deltaS, seq, dnac1, dnac2, selfComplementary,
            sodium, potassium, tris, magnesium, dntps, saltCorrection, gasConstant);
    }

    /// <summary>
    /// The final step of Biopython <c>Tm_NN</c> for a duplex whose 1 M NaCl ΔH° (kcal/mol) and ΔS°
    /// (cal/(K·mol)) are already known (e.g. an NN sum extended with LNA increments):
    /// k = (dnac1 − dnac2/2)·10⁻⁹ M (k = dnac1·10⁻⁹ M when self-complementary), salt term per
    /// <see cref="CalculateNnSaltCorrection"/> evaluated on <paramref name="saltSequence"/> (method 5 → ΔS°,
    /// 1–4 → Tm, 6–7 → 1/Tm), Tm = 1000·ΔH° / (ΔS° + R·ln k) − 273.15.
    /// </summary>
    /// <param name="deltaH">ΔH° (kcal/mol), 1 M NaCl.</param>
    /// <param name="deltaS">ΔS° (cal/(K·mol)), 1 M NaCl, including the symmetry term if any.</param>
    /// <param name="saltSequence">Sequence whose length / GC fraction the salt correction uses (Biopython: the checked primer).</param>
    /// <param name="dnac1">Higher strand concentration, nM.</param>
    /// <param name="dnac2">Lower strand concentration, nM.</param>
    /// <param name="selfComplementary">Self-complementary duplex (k = dnac1).</param>
    /// <param name="sodium">[Na⁺], mM.</param>
    /// <param name="potassium">[K⁺], mM.</param>
    /// <param name="tris">[Tris], mM.</param>
    /// <param name="magnesium">[Mg²⁺], mM.</param>
    /// <param name="dntps">[dNTPs], mM.</param>
    /// <param name="saltCorrection">Salt-correction method.</param>
    /// <param name="gasConstant">R, cal/(K·mol).</param>
    /// <returns>ΔH°, ΔS°, salt-corrected ΔS° (method 5) and Tm (°C).</returns>
    /// <exception cref="ArgumentException">A non-positive k, an invalid salt input or a degenerate duplex.</exception>
    public static NnDuplexResult CalculateNearestNeighborTmFromThermodynamics(
        double deltaH,
        double deltaS,
        string saltSequence,
        double dnac1 = 25,
        double dnac2 = 25,
        bool selfComplementary = false,
        double sodium = 50,
        double potassium = 0,
        double tris = 0,
        double magnesium = 0,
        double dntps = 0,
        NnSaltCorrection saltCorrection = NnSaltCorrection.SantaLucia1998Entropy,
        double gasConstant = NnGasConstantBiopython)
    {
        double k = (dnac1 - (dnac2 / 2.0)) * 1e-9;
        if (selfComplementary)
            k = dnac1 * 1e-9;
        if (!(k > 0))
            throw new ArgumentException("The effective strand concentration k must be > 0 (check dnac1/dnac2).", nameof(dnac1));

        double corr = 0;
        if (saltCorrection != NnSaltCorrection.None)
            corr = CalculateNnSaltCorrection(saltCorrection, sodium, potassium, tris, magnesium, dntps, saltSequence);
        double deltaSCorrected = deltaS;
        if (saltCorrection == NnSaltCorrection.SantaLucia1998Entropy)
            deltaSCorrected += corr;

        double denominator = deltaSCorrected + (gasConstant * Math.Log(k));
        if (denominator == 0)
            throw new ArgumentException("Degenerate duplex: ΔS° + R·ln k is zero (Biopython raises ZeroDivisionError).", nameof(deltaS));
        double tm = (1000 * deltaH) / denominator - 273.15;
        if (saltCorrection is NnSaltCorrection.SchildkrautLifson1965 or NnSaltCorrection.Wetmur1991
            or NnSaltCorrection.SantaLucia1996 or NnSaltCorrection.SantaLucia1998Tm)
            tm += corr;
        if (saltCorrection is NnSaltCorrection.Owczarzy2004 or NnSaltCorrection.Owczarzy2008)
        {
            double tmKelvin = tm + 273.15;
            if (tmKelvin == 0 || 1 / tmKelvin + corr == 0)
                throw new ArgumentException("Degenerate duplex: Tm = 0 K cannot be salt-corrected in 1/Tm form (Biopython raises ZeroDivisionError).", nameof(deltaH));
            tm = 1 / (1 / tmKelvin + corr) - 273.15;
        }

        return new NnDuplexResult(deltaH, deltaS, deltaSCorrected, tm);
    }

    /// <summary>
    /// Watson–Crick stack parameters (ΔH° kcal/mol, ΔS° cal/(K·mol)) of the 5'→3' top-strand dinucleotide
    /// <paramref name="dinucleotide"/> (A/C/G/T, upper case) paired with its complement, from the given table
    /// (key "XY/X'Y'" or its reverse, as <c>Tm_NN</c> looks it up).
    /// </summary>
    /// <param name="parameterSet">NN table.</param>
    /// <param name="dinucleotide">Two upper-case bases.</param>
    /// <param name="parameters">The stack parameters when found.</param>
    /// <returns><c>true</c> when the table has the stack.</returns>
    public static bool TryGetNearestNeighborStack(
        NnParameterSet parameterSet, string dinucleotide, out (double DeltaH, double DeltaS) parameters)
    {
        parameters = default;
        if (dinucleotide is null || dinucleotide.Length != 2)
            return false;
        var stacks = GetNnTable(parameterSet).Stacks;
        string key = dinucleotide + "/" + NnComplement(dinucleotide);
        if (stacks.TryGetValue(key, out var p) || stacks.TryGetValue(ReverseString(key), out p))
        {
            parameters = (p.H, p.S);
            return true;
        }
        return false;
    }

    /// <summary>
    /// Parameters (ΔH° kcal/mol, ΔS° cal/(K·mol)) of one internal nearest-neighbour step of a possibly
    /// mismatched duplex, looked up exactly as Biopython <c>Tm_NN</c>'s zipping loop does: key
    /// "<paramref name="top"/>/<paramref name="bottom"/>" first in the internal-mismatch / inosine table
    /// (DNA_IMM1; Allawi &amp; SantaLucia 1997–1998, Peyret et al. 1999), then in the Watson–Crick table of
    /// <paramref name="parameterSet"/>, each key also tried reversed.
    /// </summary>
    /// <param name="parameterSet">Watson–Crick NN table.</param>
    /// <param name="top">Top-strand dinucleotide 5'→3' (upper case).</param>
    /// <param name="bottom">The two opposite bottom-strand bases, written 3'→5' (upper case).</param>
    /// <param name="parameters">The step parameters when found.</param>
    /// <returns><c>true</c> when either table has the step.</returns>
    public static bool TryGetNearestNeighborDuplexStep(
        NnParameterSet parameterSet, string top, string bottom, out (double DeltaH, double DeltaS) parameters)
    {
        parameters = default;
        if (top is null || bottom is null || top.Length != 2 || bottom.Length != 2)
            return false;
        string key = top + "/" + bottom;
        string reversed = ReverseString(key);
        var stacks = GetNnTable(parameterSet).Stacks;
        if (DnaImm1.TryGetValue(key, out var p) || DnaImm1.TryGetValue(reversed, out p)
            || stacks.TryGetValue(key, out p) || stacks.TryGetValue(reversed, out p))
        {
            parameters = (p.H, p.S);
            return true;
        }
        return false;
    }

    /// <summary>
    /// Initiation-type terms of a nearest-neighbour table (ΔH° kcal/mol, ΔS° cal/(K·mol)) — Biopython keys
    /// <c>init</c>, <c>init_A/T</c>, <c>init_G/C</c>, <c>init_oneG/C</c>, <c>init_allA/T</c>, <c>init_5T/A</c>, <c>sym</c>.
    /// </summary>
    /// <param name="parameterSet">NN table.</param>
    /// <returns>The initiation terms.</returns>
    public static NnInitiationTerms GetNearestNeighborInitiation(NnParameterSet parameterSet)
    {
        var t = GetNnTable(parameterSet);
        return new NnInitiationTerms(t.Init, t.InitAT, t.InitGC, t.InitOneGC, t.InitAllAT, t.Init5TA, t.Sym);
    }

    /// <summary>
    /// Nearest-neighbour Tm (°C) — <see cref="CalculateNearestNeighborDuplex"/>'s
    /// <see cref="NnDuplexResult.MeltingTemperature"/> (Biopython <c>Tm_NN</c>).
    /// </summary>
    public static double CalculateNearestNeighborTm(
        string sequence,
        string? complement = null,
        int shift = 0,
        NnParameterSet parameterSet = NnParameterSet.AllawiSantaLucia1997,
        double dnac1 = 25,
        double dnac2 = 25,
        bool selfComplementary = false,
        double sodium = 50,
        double potassium = 0,
        double tris = 0,
        double magnesium = 0,
        double dntps = 0,
        NnSaltCorrection saltCorrection = NnSaltCorrection.SantaLucia1998Entropy,
        bool check = true,
        bool strict = true,
        double gasConstant = NnGasConstantBiopython) =>
        CalculateNearestNeighborDuplex(sequence, complement, shift, parameterSet, dnac1, dnac2, selfComplementary,
            sodium, potassium, tris, magnesium, dntps, saltCorrection, check, strict, gasConstant).MeltingTemperature;

    // Tm_NN up to (and including) the symmetry term; returns the checked sequence for the salt correction.
    private static (double DeltaH, double DeltaS, string Seq) NnSum(
        string sequence, string? complement, int shift, NnParameterSet parameterSet,
        bool selfComplementary, bool check, bool strict)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var table = GetNnTable(parameterSet);

        // Biopython: c_seq defaults to the complement of the raw input; then both strands are checked.
        string seq = sequence;
        string cSeq = string.IsNullOrEmpty(complement) ? NnComplement(sequence) : complement;
        if (check)
        {
            seq = NormalizeForNearestNeighbor(seq);
            cSeq = NormalizeForNearestNeighbor(cSeq);
        }
        if (seq.Length == 0)
            throw new ArgumentException("The sequence contains no A/C/G/T/I base.", nameof(sequence));

        string tmpSeq = seq, tmpCseq = cSeq;
        double dH = 0, dS = 0;

        // Dangling ends?
        if (shift != 0 || seq.Length != cSeq.Length)
        {
            if (shift > 0) tmpSeq = new string('.', shift) + seq;
            if (shift < 0) tmpCseq = new string('.', -shift) + cSeq;
            if (tmpCseq.Length > tmpSeq.Length) tmpSeq += new string('.', tmpCseq.Length - tmpSeq.Length);
            if (tmpCseq.Length < tmpSeq.Length) tmpCseq += new string('.', tmpSeq.Length - tmpCseq.Length);
            // Remove 'over-dangling' ends.
            while (tmpSeq.StartsWith("..", StringComparison.Ordinal) || tmpCseq.StartsWith("..", StringComparison.Ordinal))
            {
                tmpSeq = Drop1(tmpSeq);
                tmpCseq = Drop1(tmpCseq);
            }
            while (tmpSeq.EndsWith("..", StringComparison.Ordinal) || tmpCseq.EndsWith("..", StringComparison.Ordinal))
            {
                tmpSeq = DropLast(tmpSeq);
                tmpCseq = DropLast(tmpCseq);
            }
            if (tmpSeq.StartsWith('.') || tmpCseq.StartsWith('.'))
            {
                string leftDe = Head(tmpSeq) + "/" + Head(tmpCseq);
                if (DnaDe1.TryGetValue(leftDe, out var de)) { dH += de.H; dS += de.S; }
                else KeyError(leftDe, strict);
                tmpSeq = Drop1(tmpSeq);
                tmpCseq = Drop1(tmpCseq);
            }
            if (tmpSeq.EndsWith('.') || tmpCseq.EndsWith('.'))
            {
                string rightDe = ReverseString(Tail(tmpCseq)) + "/" + ReverseString(Tail(tmpSeq));
                if (DnaDe1.TryGetValue(rightDe, out var de)) { dH += de.H; dS += de.S; }
                else KeyError(rightDe, strict);
                tmpSeq = DropLast(tmpSeq);
                tmpCseq = DropLast(tmpCseq);
            }
        }

        // Terminal mismatches.
        string leftTmm = ReverseString(Head(tmpCseq)) + "/" + ReverseString(Head(tmpSeq));
        if (DnaTmm1.TryGetValue(leftTmm, out var lt))
        {
            dH += lt.H; dS += lt.S;
            tmpSeq = Drop1(tmpSeq);
            tmpCseq = Drop1(tmpCseq);
        }
        string rightTmm = Tail(tmpSeq) + "/" + Tail(tmpCseq);
        if (DnaTmm1.TryGetValue(rightTmm, out var rt))
        {
            dH += rt.H; dS += rt.S;
            tmpSeq = DropLast(tmpSeq);
            tmpCseq = DropLast(tmpCseq);
        }

        // Initiation.
        dH += table.Init.H; dS += table.Init.S;
        if (NnGcFraction(seq) == 0) { dH += table.InitAllAT.H; dS += table.InitAllAT.S; }
        else { dH += table.InitOneGC.H; dS += table.InitOneGC.S; }
        if (seq.StartsWith('T')) { dH += table.Init5TA.H; dS += table.Init5TA.S; }
        if (seq.EndsWith('A')) { dH += table.Init5TA.H; dS += table.Init5TA.S; }
        char e0 = seq[0], e1 = seq[^1];
        int at = (e0 is 'A' or 'T' ? 1 : 0) + (e1 is 'A' or 'T' ? 1 : 0);
        int gc = (e0 is 'G' or 'C' ? 1 : 0) + (e1 is 'G' or 'C' ? 1 : 0);
        dH += table.InitAT.H * at; dS += table.InitAT.S * at;
        dH += table.InitGC.H * gc; dS += table.InitGC.S * gc;

        // Zipping.
        for (int i = 0; i < tmpSeq.Length - 1; i++)
        {
            string neighbors = Slice2(tmpSeq, i) + "/" + Slice2(tmpCseq, i);
            string reversed = ReverseString(neighbors);
            if (DnaImm1.TryGetValue(neighbors, out var p) || DnaImm1.TryGetValue(reversed, out p)
                || table.Stacks.TryGetValue(neighbors, out p) || table.Stacks.TryGetValue(reversed, out p))
            {
                dH += p.H; dS += p.S;
            }
            else
            {
                KeyError(neighbors, strict);
            }
        }

        if (selfComplementary)
        {
            dH += table.Sym.H; dS += table.Sym.S;
        }
        return (dH, dS, seq);

        // Python slices never throw: s[:2], s[-2:], s[i:i+2] of a short string are shorter.
        static string Head(string s) => s.Length <= 2 ? s : s[..2];
        static string Tail(string s) => s.Length <= 2 ? s : s[^2..];
        static string Drop1(string s) => s.Length == 0 ? s : s[1..];
        static string DropLast(string s) => s.Length == 0 ? s : s[..^1];
        static string Slice2(string s, int i) => s.Substring(i, Math.Min(2, s.Length - i));
        static void KeyError(string neighbors, bool strict)
        {
            if (strict)
                throw new ArgumentException($"No thermodynamic data for neighbors '{neighbors}' available.");
        }
    }

    #endregion
}

/// <summary>DNA/DNA Watson–Crick nearest-neighbour parameter sets (Biopython <c>MeltingTemp</c> tables).</summary>
public enum NnParameterSet
{
    /// <summary>Breslauer et al. (1986) PNAS 83:3746 — Biopython <c>DNA_NN1</c>.</summary>
    Breslauer1986,

    /// <summary>Sugimoto et al. (1996) NAR 24:4501 — Biopython <c>DNA_NN2</c>.</summary>
    Sugimoto1996,

    /// <summary>Allawi &amp; SantaLucia (1997) Biochemistry 36:10581 — Biopython <c>DNA_NN3</c> (Tm_NN default;
    /// SantaLucia 1998 unified stacks, initiation per terminal A·T (+2.3/+4.1) / G·C (+0.1/−2.8) pair).</summary>
    AllawiSantaLucia1997,

    /// <summary>SantaLucia &amp; Hicks (2004) Annu Rev Biophys Biomol Struct 33:415, Table 1 — Biopython
    /// <c>DNA_NN4</c> (duplex initiation +0.2/−5.7, terminal A·T penalty +2.2/+6.9, AA/TT −7.6/−21.3).</summary>
    SantaLuciaHicks2004
}

/// <summary>Salt-correction methods of Biopython <c>MeltingTemp.salt_correction</c> (numeric values = Biopython method numbers).</summary>
public enum NnSaltCorrection
{
    /// <summary>No correction (method 0).</summary>
    None = 0,

    /// <summary>16.6·log10[Mon], Schildkraut &amp; Lifson (1965) — added to Tm.</summary>
    SchildkrautLifson1965 = 1,

    /// <summary>16.6·log10([Mon]/(1 + 0.7[Mon])), Wetmur (1991) — added to Tm.</summary>
    Wetmur1991 = 2,

    /// <summary>12.5·log10[Mon], SantaLucia et al. (1996) — added to Tm.</summary>
    SantaLucia1996 = 3,

    /// <summary>11.7·log10[Mon], SantaLucia (1998) — added to Tm.</summary>
    SantaLucia1998Tm = 4,

    /// <summary>0.368·(N − 1)·ln[Mon], SantaLucia (1998) — added to ΔS° (Tm_NN default).</summary>
    SantaLucia1998Entropy = 5,

    /// <summary>Owczarzy et al. (2004) Biochemistry 43:3537 — added to 1/Tm.</summary>
    Owczarzy2004 = 6,

    /// <summary>Owczarzy et al. (2008) Biochemistry 47:5336 (Mg²⁺, dNTPs) — added to 1/Tm.</summary>
    Owczarzy2008 = 7
}

/// <summary>
/// Initiation-type terms of a nearest-neighbour table, (ΔH° kcal/mol, ΔS° cal/(K·mol)) each: general duplex
/// initiation, per terminal A·T pair, per terminal G·C pair, duplex with ≥ 1 G·C pair, all-A·T duplex, per
/// 5'-T / 3'-A terminus, and the symmetry correction of a self-complementary duplex.
/// </summary>
public readonly record struct NnInitiationTerms(
    (double DeltaH, double DeltaS) Initiation,
    (double DeltaH, double DeltaS) TerminalAT,
    (double DeltaH, double DeltaS) TerminalGC,
    (double DeltaH, double DeltaS) OneGC,
    (double DeltaH, double DeltaS) AllAT,
    (double DeltaH, double DeltaS) FiveTerminalTA,
    (double DeltaH, double DeltaS) Symmetry);

/// <summary>
/// Result of <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/>: ΔH° (kcal/mol) and ΔS°
/// (cal/(K·mol)) at 1 M NaCl, ΔS° after the SantaLucia (1998) method-5 entropy correction (equal to
/// <see cref="DeltaS"/> for every other method) and the melting temperature (°C).
/// </summary>
public readonly record struct NnDuplexResult(
    double DeltaH,
    double DeltaS,
    double SaltCorrectedDeltaS,
    double MeltingTemperature);

