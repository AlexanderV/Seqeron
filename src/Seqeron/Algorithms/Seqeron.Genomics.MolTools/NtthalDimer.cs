using System;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Faithful port of the Primer3 <c>ntthal</c> dimer thermodynamic-alignment engine
/// (oligo–oligo hybridisation, <c>type == 1</c> / mode ANY) for PRIMER-TM-001.
/// <para>
/// Computes the most stable intermolecular DNA duplex between two oligonucleotides via the
/// full ntthal dynamic program — matched nearest-neighbour stacks, single internal mismatches
/// (1×1 via the SantaLucia internal-mismatch table), internal loops (≥ 2 unpaired bases per
/// side via the interior-loop length parameters + the <c>tstack</c> terminal table + the
/// −0.3 kcal/mol asymmetry correction), single- and multi-base bulges (bulge-loop length
/// parameters; 1-base bulges add the intervening NN stack), and terminal overhangs / dangling
/// ends (the <c>tstack2</c> terminal-stacking table and the 5′/3′ dangling-end tables). It then
/// returns the ntthal ΔH° (cal/mol), ΔS° (cal/(K·mol), salt-corrected), ΔG°37 (cal/mol) and the
/// bimolecular melting temperature (°C).
/// </para>
/// <para>
/// This is a line-for-line translation of <c>thal.c</c> (<c>fillMatrix</c>, <c>LSH</c>,
/// <c>RSH</c>, <c>maxTM</c>, <c>calc_bulge_internal</c>, <c>traceback</c>, <c>calcDimer</c>) from
/// the primer3-py vendored libprimer3, retrieved this session from
/// https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c .
/// All thermodynamic tables (<c>stack</c>, <c>stackmm</c>, <c>tstack2</c>, <c>tstack</c>,
/// <c>dangle</c>, interior/bulge loop lengths) are taken verbatim from the primer3 config files
/// (<c>primer3_config/*.dh</c>, <c>*.ds</c>), the same authoritative parameter set ntthal loads.
/// </para>
/// <para>
/// Cross-checked against primer3-py 2.3.0 <c>calc_homodimer</c> / <c>calc_heterodimer</c>
/// (mv = 50 mM, dv = 0, dntp = 0, dna_conc = 50 nM): this engine reproduces ntthal's ΔH°, ΔS°,
/// ΔG° and Tm to machine precision for contiguous Watson–Crick duplexes <b>and</b> for dimers
/// whose optimal structure contains an internal mismatch, an internal loop, a bulge, or a
/// terminal overhang.
/// </para>
/// <para>
/// Sources (retrieved &amp; extracted this session):
/// SantaLucia &amp; Hicks (2004) Annu Rev Biophys 33:415-440 (unified NN parameters; loop model);
/// Untergasser et al. (2012) Nucleic Acids Res 40:e115 (Primer3 2.0 / ntthal engine);
/// primer3 <c>thal.c</c> + <c>primer3_config</c> tables (URLs above).
/// </para>
/// </summary>
internal static class NtthalDimer
{
    // ---- physical constants (thal.c lines 125-135, 588-589) -----------------------------
    // Several of these constants/tables are shared verbatim with the ntthal HAIRPIN engine
    // (NtthalHairpin.cs); those are marked internal so both engines use a single source of truth.
    internal const double Inf = double.PositiveInfinity;
    internal const double R = 1.9872;                 // gas constant, cal/(K·mol)
    internal const double IlAs = -300.0 / 310.15;     // internal-loop entropy asymmetry, cal/(K·mol) per |Δn|
    internal const double IlAh = 0.0;                 // internal-loop enthalpy asymmetry
    internal const double AtH = 2200.0;               // terminal-A·T penalty ΔH, cal/mol
    internal const double AtS = 6.9;                  // terminal-A·T penalty ΔS, cal/(K·mol)
    internal const double MinEntropyCutoff = -2500.0; // filters out non-existing entropies
    internal const double MinEntropy = -3224.0;       // initiation sentinel
    internal const double TempKelvin = 310.15;        // 37 °C reference for the internal ΔG ranking
    internal const double AbsoluteZero = 273.15;
    internal const int MaxLoop = 30;                  // maximum loop length the loop tables cover
    private const double DplxInitH = 200.0;          // duplex initiation ΔH, cal/mol
    private const double DplxInitS = -5.7;           // duplex initiation ΔS, cal/(K·mol)
    internal const double AtPenaltySEntry = 1e-11;    // tableStartATS default (non-A·T) entropy
    private const double Equal = 1e-6;               // traceback equality tolerance (thal.c equal())

    // bp index matrix BPI[5][5] (A,C,G,T,N): 1 = Watson-Crick pair, 0 = none (thal.c lines 140-145).
    internal static readonly int[,] Bpi =
    {
        { 0, 0, 0, 1, 0 },
        { 0, 0, 1, 0, 0 },
        { 0, 1, 0, 0, 0 },
        { 1, 0, 0, 0, 0 },
        { 0, 0, 0, 0, 0 },
    };

    internal static int Str2Int(char c) => c switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' => 3,
        _ => 4,
    };

    internal static bool IsFinite(double x) => !double.IsInfinity(x);

    /// <summary>The most stable dimer's ntthal thermodynamics (native ntthal units).</summary>
    /// <param name="DeltaH">Dimer ΔH° in cal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the N·saltCorrection term.</param>
    /// <param name="DeltaG37">Dimer ΔG°37 = ΔH° − 310.15·ΔS° in cal/mol.</param>
    /// <param name="TmCelsius">Bimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">Number of paired bases N+1 in the optimal structure.</param>
    /// <param name="Strand1End">1-based 3′-most paired index on strand 1 (ntthal align_end_1).</param>
    /// <param name="Strand2End">1-based paired index on the reversed strand 2 (ntthal align_end_2).</param>
    internal readonly record struct Result(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius,
        int BasePairs, int Strand1End, int Strand2End);

    /// <summary>
    /// Runs the full ntthal dimer DP for two oligos (5′→3′). Returns <c>null</c> when no duplex
    /// can be formed (no finite terminal base pair, i.e. ntthal <c>no_structure</c>).
    /// </summary>
    /// <param name="oligo1">Strand 1 (5′→3'), ACGT only.</param>
    /// <param name="oligo2">Strand 2 (5′→3'), ACGT only.</param>
    /// <param name="mvMolar">Monovalent cation concentration in mol/L (ntthal mv is in mM).</param>
    /// <param name="dnaConcMolar">Total strand concentration in mol/L (ntthal dna_conc in nM).</param>
    internal static Result? Run(string oligo1, string oligo2, double mvMolar, double dnaConcMolar)
    {
        // ntthal mv is in mM, dna_conc in nM; convert from the SI (mol/L) the caller passes.
        double mv = mvMolar * 1000.0;
        double dnaConc = dnaConcMolar * 1e9;

        // numSeq1 = oligo1 (5'->3'); numSeq2 = REVERSE of oligo2 so the DP runs 3'->5' (thal.c 602).
        int len1 = oligo1.Length;
        int len2 = oligo2.Length;
        // 1-indexed sequences with N (=4) sentinels at indices 0 and len+1 (thal.c line 634).
        var a = new int[len1 + 2];
        var b = new int[len2 + 2];
        a[0] = a[len1 + 1] = 4;
        b[0] = b[len2 + 1] = 4;
        for (int i = 0; i < len1; i++) a[i + 1] = Str2Int(oligo1[i]);
        for (int j = 0; j < len2; j++) b[j + 1] = Str2Int(oligo2[len2 - 1 - j]); // reversed

        // RC = R·ln(dna_conc / x), x=1e9 if both palindromic (symmetry_thermo) else 4e9 (thal.c 590-593).
        bool symmetric = IsSymmetric(oligo1) && IsSymmetric(oligo2);
        double rc = symmetric ? R * Math.Log(dnaConc / 1e9) : R * Math.Log(dnaConc / 4e9);
        // saltCorrectS (thal.c 1042); dv=dntp=0 here so the divalent term vanishes.
        double saltCorrection = 0.368 * Math.Log(mv / 1000.0);

        // A·T penalty tables (thal.c tableStartATH/ATS): only A·T (0,3)/(3,0) carry the penalty.
        double AtPenaltyH(int x, int y) => (x == 0 && y == 3) || (x == 3 && y == 0) ? AtH : 0.0;
        double AtPenaltyS(int x, int y) => (x == 0 && y == 3) || (x == 3 && y == 0) ? AtS : AtPenaltySEntry;

        int Bp(int x, int y) => Bpi[x, y];

        // DP tables (1-indexed). EnthalpyDPT/EntropyDPT.
        var enH = new double[len1 + 2, len2 + 2];
        var enS = new double[len1 + 2, len2 + 2];

        // ---- table accessors (4-D flat arrays indexed [i][ii][j][jj] -> i*125+ii*25+j*5+jj) ----
        static double T4(double[] t, int i, int ii, int j, int jj) => t[((i * 5 + ii) * 5 + j) * 5 + jj];
        static double T3(double[] t, int i, int j, int k) => t[(i * 5 + j) * 5 + k];

        // initMatrix (thal.c 1547-1562).
        for (int i = 1; i <= len1; i++)
        {
            for (int j = 1; j <= len2; j++)
            {
                if (Bp(a[i], b[j]) == 0) { enH[i, j] = Inf; enS[i, j] = -1.0; }
                else { enH[i, j] = 0.0; enS[i, j] = MinEntropy; }
            }
        }

        // Ss / Hs — stack ΔS/ΔH for the extension step (k=1 path used by maxTM/traceback) (thal.c 1985-2034).
        double Ss(int i, int j) => T4(StackS, a[i], a[i + 1], b[j], b[j + 1]);
        double Hs(int i, int j) => T4(StackH, a[i], a[i + 1], b[j], b[j + 1]);

        // RSH — right terminal stack (3'-side): tstack2 terminal-mismatch / dangling-end (thal.c 1857-1983).
        (double S, double H) Rsh(int i, int j)
        {
            if (Bp(a[i], b[j]) == 0) return (-1.0, Inf);
            double s1 = AtPenaltyS(a[i], b[j]) + T4(Tstack2S, a[i], a[i + 1], b[j], b[j + 1]);
            double h1 = AtPenaltyH(a[i], b[j]) + T4(Tstack2H, a[i], a[i + 1], b[j], b[j + 1]);
            double g1 = h1 - TempKelvin * s1;
            if (!IsFinite(h1) || g1 > 0) { h1 = Inf; s1 = -1.0; g1 = 1.0; }
            double s2, h2, t2;
            double t1;

            bool unpaired = Bp(a[i + 1], b[j + 1]) == 0;
            bool d3 = IsFinite(T3(Dangle3H, a[i], a[i + 1], b[j]));
            bool d5 = IsFinite(T3(Dangle5H, a[i], b[j], b[j + 1]));
            if (unpaired && d3 && d5)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle3S, a[i], a[i + 1], b[j]) + T3(Dangle5S, a[i], b[j], b[j + 1]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle3H, a[i], a[i + 1], b[j]) + T3(Dangle5H, a[i], b[j], b[j + 1]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (unpaired && d3)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle3S, a[i], a[i + 1], b[j]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle3H, a[i], a[i + 1], b[j]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (unpaired && d5)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle5S, a[i], b[j], b[j + 1]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle5H, a[i], b[j], b[j + 1]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            // bare A·T-penalty alternative (no terminal stack).
            s2 = AtPenaltyS(a[i], b[j]);
            h2 = AtPenaltyH(a[i], b[j]);
            t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
            if (IsFinite(h1))
            {
                t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                return t1 < t2 ? (s2, h2) : (s1, h1);
            }
            return (s2, h2);
        }

        // LSH — left terminal stack (5'-side) (thal.c 1741-1855).
        (double S, double H)? Lsh(int i, int j)
        {
            if (Bp(a[i], b[j]) == 0) { enS[i, j] = -1.0; enH[i, j] = Inf; return null; }
            double s1 = AtPenaltyS(a[i], b[j]) + T4(Tstack2S, b[j], b[j - 1], a[i], a[i - 1]);
            double h1 = AtPenaltyH(a[i], b[j]) + T4(Tstack2H, b[j], b[j - 1], a[i], a[i - 1]);
            double g1 = h1 - TempKelvin * s1;
            if (!IsFinite(h1) || g1 > 0) { h1 = Inf; s1 = -1.0; g1 = 1.0; }
            double s2, h2, t2, t1;

            bool notPaired = Bp(a[i - 1], b[j - 1]) != 1;
            bool d3 = IsFinite(T3(Dangle3H, b[j], b[j - 1], a[i]));
            bool d5 = IsFinite(T3(Dangle5H, b[j], a[i], a[i - 1]));
            if (notPaired && d3 && d5)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle3S, b[j], b[j - 1], a[i]) + T3(Dangle5S, b[j], a[i], a[i - 1]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle3H, b[j], b[j - 1], a[i]) + T3(Dangle5H, b[j], a[i], a[i - 1]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (notPaired && d3)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle3S, b[j], b[j - 1], a[i]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle3H, b[j], b[j - 1], a[i]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (notPaired && d5)
            {
                s2 = AtPenaltyS(a[i], b[j]) + T3(Dangle5S, b[j], a[i], a[i - 1]);
                h2 = AtPenaltyH(a[i], b[j]) + T3(Dangle5H, b[j], a[i], a[i - 1]);
                double g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
                if (IsFinite(h1) && g1 < 0)
                {
                    t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                    if (t1 < t2 && g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            s2 = AtPenaltyS(a[i], b[j]);
            h2 = AtPenaltyH(a[i], b[j]);
            t2 = (h2 + DplxInitH) / (s2 + DplxInitS + rc);
            if (IsFinite(h1))
            {
                t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc);
                return t1 < t2 ? (s2, h2) : (s1, h1);
            }
            return (s2, h2);
        }

        // maxTM — stack extension from (i-1,j-1) vs current (thal.c 1662-1702).
        void MaxTm(int i, int j)
        {
            double s0 = enS[i, j], h0 = enH[i, j];
            var (s, h) = Rsh(i, j);
            double t0 = (h0 + DplxInitH + h) / (s0 + DplxInitS + s + rc);
            double s1, h1, t1;
            if (IsFinite(enH[i - 1, j - 1]) && IsFinite(Hs(i - 1, j - 1)))
            {
                s1 = enS[i - 1, j - 1] + Ss(i - 1, j - 1);
                h1 = enH[i - 1, j - 1] + Hs(i - 1, j - 1);
                t1 = (h1 + DplxInitH + h) / (s1 + DplxInitS + s + rc);
            }
            else { s1 = -1.0; h1 = Inf; t1 = (h1 + DplxInitH) / (s1 + DplxInitS + rc); }
            if (s1 < MinEntropyCutoff) { s1 = MinEntropy; h1 = 0.0; }
            if (s0 < MinEntropyCutoff) { s0 = MinEntropy; h0 = 0.0; }
            if (t1 > t0) { enS[i, j] = s1; enH[i, j] = h1; }
            else { enS[i, j] = s0; enH[i, j] = h0; }
        }

        // calc_bulge_internal — bulges + internal loops closing at (i,j) from inner pair (ii,jj)
        // (thal.c 2149-2304). traceback=false returns the candidate only when it lowers ΔG.
        (double S, double H)? CalcBulgeInternal(int i, int j, int ii, int jj, bool traceback)
        {
            int loopSize1 = ii - i - 1;
            int loopSize2 = jj - j - 1;
            int loopSize = loopSize1 + loopSize2 - 1;
            var (rs, rh) = Rsh(ii, jj);
            double s, h, g1, g2;

            if ((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0))
            {
                if (loopSize2 == 1 || loopSize1 == 1) // bulge of size one: add intervening NN stack
                {
                    h = Inf; s = -1.0;
                    if ((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1))
                    {
                        h = BulgeH[loopSize] + T4(StackH, a[i], a[ii], b[j], b[jj]);
                        s = BulgeS[loopSize] + T4(StackS, a[i], a[ii], b[j], b[jj]);
                    }
                    if (h > 0 || s > 0) { h = Inf; s = -1.0; }
                    h += enH[i, j]; s += enS[i, j];
                    if (!IsFinite(h)) { h = Inf; s = -1.0; }
                    g1 = h + rh - TempKelvin * (s + rs);
                    g2 = enH[ii, jj] + rh - TempKelvin * (enS[ii, jj] + rs);
                    if (g1 < g2 || traceback) return (s, h);
                }
                else // larger bulge: two terminal-A·T penalties
                {
                    h = BulgeH[loopSize] + AtPenaltyH(a[i], b[j]) + AtPenaltyH(a[ii], b[jj]) + enH[i, j];
                    s = BulgeS[loopSize] + AtPenaltyS(a[i], b[j]) + AtPenaltyS(a[ii], b[jj]) + enS[i, j];
                    if (!IsFinite(h)) { h = Inf; s = -1.0; }
                    if (h > 0 && s > 0) { h = Inf; s = -1.0; }
                    g1 = h + rh - TempKelvin * (s + rs);
                    g2 = enH[ii, jj] + rh - TempKelvin * (enS[ii, jj] + rs);
                    if (g1 < g2 || traceback) return (s, h);
                }
            }
            else if (loopSize1 == 1 && loopSize2 == 1) // single internal mismatch (1x1): stackmm both sides
            {
                s = T4(Int2S, a[i], a[i + 1], b[j], b[j + 1]) + T4(Int2S, b[jj], b[jj - 1], a[ii], a[ii - 1]) + enS[i, j];
                h = T4(Int2H, a[i], a[i + 1], b[j], b[j + 1]) + T4(Int2H, b[jj], b[jj - 1], a[ii], a[ii - 1]) + enH[i, j];
                if (!IsFinite(h)) { h = Inf; s = -1.0; }
                if (h > 0 && s > 0) { h = Inf; s = -1.0; }
                g1 = h + rh - TempKelvin * (s + rs);
                g2 = enH[ii, jj] + rh - TempKelvin * (enS[ii, jj] + rs);
                if (g1 < g2 || traceback) return (s, h);
            }
            else // general internal loop: interior-loop length + tstack terminal both sides + asymmetry
            {
                h = InteriorH[loopSize] + T4(TstackH, a[i], a[i + 1], b[j], b[j + 1]) +
                    T4(TstackH, b[jj], b[jj - 1], a[ii], a[ii - 1]) + IlAh * Math.Abs(loopSize1 - loopSize2) + enH[i, j];
                s = InteriorS[loopSize] + T4(TstackS, a[i], a[i + 1], b[j], b[j + 1]) +
                    T4(TstackS, b[jj], b[jj - 1], a[ii], a[ii - 1]) + IlAs * Math.Abs(loopSize1 - loopSize2) + enS[i, j];
                if (!IsFinite(h)) { h = Inf; s = -1.0; }
                if (h > 0 && s > 0) { h = Inf; s = -1.0; }
                g1 = h + rh - TempKelvin * (s + rs);
                g2 = enH[ii, jj] + rh - TempKelvin * (enS[ii, jj] + rs);
                if (g1 < g2 || traceback) return (s, h);
            }
            return null;
        }

        // fillMatrix (thal.c 1581-1629).
        for (int i = 1; i <= len1; i++)
        {
            for (int j = 1; j <= len2; j++)
            {
                if (!IsFinite(enH[i, j])) continue;
                var lsh = Lsh(i, j);
                if (lsh is { } l && IsFinite(l.H)) { enS[i, j] = l.S; enH[i, j] = l.H; }
                if (i > 1 && j > 1)
                {
                    MaxTm(i, j);
                    for (int d = 3; d <= MaxLoop + 2; d++)
                    {
                        int ii = i - 1;
                        int jj = -ii - d + (j + i);
                        if (jj < 1) { ii -= Math.Abs(jj - 1); jj = 1; }
                        for (; ii > 0 && jj < j; --ii, ++jj)
                        {
                            if (!IsFinite(enH[ii, jj])) continue;
                            var r = CalcBulgeInternal(ii, jj, i, j, traceback: false);
                            if (r is { } rr)
                            {
                                double cs = rr.S, ch = rr.H;
                                if (cs < MinEntropyCutoff) { cs = MinEntropy; ch = 0.0; }
                                if (IsFinite(ch)) { enH[i, j] = ch; enS[i, j] = cs; }
                            }
                        }
                    }
                }
            }
        }

        // Best terminal base pair over all (i,j) (type==1 / ANY) (thal.c 710-723).
        double bestG = Inf;
        int bestI = 0, bestJ = 0;
        for (int i = 1; i <= len1; i++)
        {
            for (int j = 1; j <= len2; j++)
            {
                var (s, h) = Rsh(i, j);
                double g = (enH[i, j] + h + DplxInitH) - TempKelvin * (enS[i, j] + s + DplxInitS);
                if (g < bestG) { bestG = g; bestI = i; bestJ = j; }
            }
        }
        if (!IsFinite(bestG)) return null; // ntthal no_structure

        var (bs, bh) = Rsh(bestI, bestJ);
        double dH = enH[bestI, bestJ] + bh + DplxInitH;
        double dS = enS[bestI, bestJ] + bs + DplxInitS;

        // traceback to count paired bases N (thal.c 2957-3003).
        bool Eq(double x, double y)
        {
            if (double.IsInfinity(x) && double.IsInfinity(y)) return (x > 0) == (y > 0);
            if (double.IsInfinity(x) || double.IsInfinity(y)) return false;
            return Math.Abs(x - y) < Equal;
        }

        var ps1 = new int[len1];
        var ps2 = new int[len2];
        {
            int i = bestI, j = bestJ;
            ps1[i - 1] = j; ps2[j - 1] = i;
            while (true)
            {
                var lsh = Lsh(i, j);
                if (lsh is { } l && Eq(enS[i, j], l.S) && Eq(enH[i, j], l.H)) break;
                bool done = false;
                if (i > 1 && j > 1 &&
                    Eq(enS[i, j], Ss(i - 1, j - 1) + enS[i - 1, j - 1]) &&
                    Eq(enH[i, j], Hs(i - 1, j - 1) + enH[i - 1, j - 1]))
                {
                    i -= 1; j -= 1; ps1[i - 1] = j; ps2[j - 1] = i; done = true;
                }
                for (int d = 3; !done && d <= MaxLoop + 2; d++)
                {
                    int ii = i - 1;
                    int jj = -ii - d + (j + i);
                    if (jj < 1) { ii -= Math.Abs(jj - 1); jj = 1; }
                    for (; !done && ii > 0 && jj < j; --ii, ++jj)
                    {
                        var r = CalcBulgeInternal(ii, jj, i, j, traceback: true);
                        if (r is { } rr && Eq(enS[i, j], rr.S) && Eq(enH[i, j], rr.H))
                        {
                            i = ii; j = jj; ps1[i - 1] = j; ps2[j - 1] = i; done = true; break;
                        }
                    }
                }
            }
        }

        int n = 0;
        for (int i = 0; i < len1; i++) if (ps1[i] > 0) n++;
        for (int j = 0; j < len2; j++) if (ps2[j] > 0) n++;
        n = n / 2 - 1; // number of NN stacks (thal.c calcDimer line 3027)

        double tm = dH / (dS + n * saltCorrection + rc) - AbsoluteZero;
        double dsOut = dS + n * saltCorrection;
        double dg = dH - TempKelvin * dsOut;
        return new Result(dH, dsOut, dg, tm, n + 1, bestI, bestJ);
    }

    /// <summary>Reverse-complement palindrome test (thal.c symmetry_thermo).</summary>
    private static bool IsSymmetric(string oligo)
    {
        int len = oligo.Length;
        if (len % 2 != 0) return false;
        for (int i = 0; i < len; i++)
        {
            char x = char.ToUpperInvariant(oligo[i]);
            char y = char.ToUpperInvariant(oligo[len - 1 - i]);
            bool wc = (x == 'A' && y == 'T') || (x == 'T' && y == 'A') ||
                      (x == 'C' && y == 'G') || (x == 'G' && y == 'C');
            if (!wc) return false;
        }
        return true;
    }

    // ===================== thermodynamic parameter tables =====================
    // All values copied verbatim from the primer3 config files (primer3_config/*.dh, *.ds),
    // the authoritative parameter set loaded by ntthal, ordered [i][ii][j][jj] (A,C,G,T,N=0..4)
    // exactly as getStack/getStackint2/getTstack2/getTstack read them. 4-D tables are flattened
    // row-major (index = ((i*5+ii)*5+j)*5+jj); 3-D dangle tables as i*25 + col*5 + col2 per getDangle.
    // INF entries are non-existent pairings/stacks (thal.c marks them _INFINITY/-1).
    internal static readonly double[] StackH = {
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-7900.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-8400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-7800.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-7200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-8500.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-8000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-10600.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-7800.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-8200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,-9800.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,-8000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,-8400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,-7200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-8200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-8500.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -7900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] StackS = {
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-22.2d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-22.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-21.0d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-20.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-22.7d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-19.9d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-27.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-21.0d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-22.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-24.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-19.9d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-22.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-21.3d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-22.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-22.7d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -22.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] Int2H = {
        Inf,Inf,Inf,4700.0d,Inf,Inf,Inf,Inf,7600.0d,Inf,Inf,Inf,Inf,3000.0d,Inf,1200.0d,2300.0d,-600.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-2900.0d,Inf,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,500.0d,Inf,Inf,5300.0d,-10.0d,Inf,700.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-900.0d,Inf,Inf,Inf,Inf,600.0d,Inf,Inf,Inf,Inf,-4000.0d,Inf,Inf,Inf,-700.0d,Inf,-3100.0d,1000.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        1200.0d,Inf,Inf,Inf,Inf,5300.0d,Inf,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,Inf,-1200.0d,-2500.0d,-2700.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,3400.0d,Inf,Inf,Inf,Inf,6100.0d,Inf,-900.0d,1900.0d,-700.0d,Inf,Inf,Inf,Inf,Inf,1000.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,5200.0d,Inf,Inf,Inf,Inf,3600.0d,Inf,Inf,600.0d,-1500.0d,Inf,-800.0d,Inf,Inf,Inf,5200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,1900.0d,Inf,Inf,Inf,Inf,-1500.0d,Inf,Inf,Inf,-4000.0d,Inf,-4900.0d,-4100.0d,Inf,Inf,-1500.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        2300.0d,Inf,Inf,Inf,Inf,-10.0d,Inf,Inf,Inf,Inf,Inf,-1500.0d,-2800.0d,-5000.0d,Inf,-1200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,700.0d,Inf,-2900.0d,5200.0d,-600.0d,Inf,Inf,Inf,Inf,Inf,1600.0d,Inf,Inf,Inf,Inf,-1300.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-600.0d,Inf,Inf,-700.0d,3600.0d,Inf,2300.0d,Inf,Inf,Inf,-6000.0d,Inf,Inf,Inf,Inf,-4400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-700.0d,Inf,Inf,Inf,500.0d,Inf,-6000.0d,3300.0d,Inf,Inf,-4900.0d,Inf,Inf,Inf,Inf,-2800.0d,Inf,5800.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        -600.0d,Inf,Inf,Inf,Inf,Inf,5200.0d,-4400.0d,-2200.0d,Inf,-3100.0d,Inf,Inf,Inf,Inf,-2500.0d,Inf,4100.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,3400.0d,700.0d,Inf,Inf,Inf,Inf,Inf,1200.0d,Inf,Inf,Inf,Inf,-100.0d,Inf,Inf,Inf,Inf,200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        7600.0d,6100.0d,Inf,1200.0d,Inf,Inf,Inf,2300.0d,Inf,Inf,Inf,Inf,3300.0d,Inf,Inf,Inf,Inf,-2200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        3000.0d,Inf,1600.0d,-100.0d,Inf,Inf,-800.0d,Inf,Inf,Inf,Inf,-4100.0d,Inf,-1400.0d,Inf,Inf,-5000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,1000.0d,-1300.0d,200.0d,Inf,700.0d,Inf,Inf,Inf,Inf,1000.0d,Inf,5800.0d,Inf,Inf,-2700.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] Int2S = {
        -1d,-1d,-1d,12.9d,-1d,-1d,-1d,-1d,20.2d,-1d,-1d,-1d,-1d,7.4d,-1d,1.7d,4.6d,-2.3d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-9.8d,-1d,-1d,-1d,-1d,-3.8d,-1d,-1d,-1d,-1d,3.2d,-1d,-1d,14.6d,-4.4d,-1d,0.2d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-4.2d,-1d,-1d,-1d,-1d,-0.6d,-1d,-1d,-1d,-1d,-13.2d,-1d,-1d,-1d,-2.3d,-1d,-9.5d,0.9d,-1d,-1d,-1d,-1d,-1d,-1d,
        1.7d,-1d,-1d,-1d,-1d,14.6d,-1d,-1d,-1d,-1d,-2.3d,-1d,-1d,-1d,-1d,-1d,-6.2d,-8.3d,-10.8d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,8.0d,-1d,-1d,-1d,-1d,16.4d,-1d,-4.2d,3.7d,-2.3d,-1d,-1d,-1d,-1d,-1d,0.7d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,14.2d,-1d,-1d,-1d,-1d,8.9d,-1d,-1d,-0.6d,-7.2d,-1d,-4.5d,-1d,-1d,-1d,13.5d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,3.7d,-1d,-1d,-1d,-1d,-7.2d,-1d,-1d,-1d,-13.2d,-1d,-15.3d,-11.7d,-1d,-1d,-6.1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        4.6d,-1d,-1d,-1d,-1d,-4.4d,-1d,-1d,-1d,-1d,-1d,-6.1d,-8.0d,-15.8d,-1d,-6.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,0.7d,-1d,-9.8d,14.2d,-1d,-1d,-1d,-1d,-1d,-1d,3.6d,-1d,-1d,-1d,-1d,-5.3d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-3.8d,8.9d,-1d,5.4d,-1d,-1d,-1d,-15.8d,-1d,-1d,-1d,-1d,-12.3d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-2.3d,-1d,-1d,-1d,3.2d,-1d,-15.8d,10.4d,-1d,-1d,-15.3d,-1d,-1d,-1d,-1d,-8.0d,-1d,16.3d,-1d,-1d,-1d,-1d,-1d,-1d,
        -2.3d,-1d,-1d,-1d,-1d,-1d,13.5d,-12.3d,-8.4d,-1d,-9.5d,-1d,-1d,-1d,-1d,-8.3d,-1d,9.5d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,8.0d,0.7d,-1d,-1d,-1d,-1d,-1d,0.7d,-1d,-1d,-1d,-1d,-1.7d,-1d,-1d,-1d,-1d,-1.5d,-1d,-1d,-1d,-1d,-1d,-1d,
        20.2d,16.4d,-1d,0.7d,-1d,-1d,-1d,5.4d,-1d,-1d,-1d,-1d,10.4d,-1d,-1d,-1d,-1d,-8.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        7.4d,-1d,3.6d,-1.7d,-1d,-1d,-4.5d,-1d,-1d,-1d,-1d,-11.7d,-1d,-6.2d,-1d,-1d,-15.8d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,0.7d,-5.3d,-1.5d,-1d,0.2d,-1d,-1d,-1d,-1d,0.9d,-1d,16.3d,-1d,-1d,-10.8d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] Tstack2H = {
        Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,-2700.0d,Inf,Inf,Inf,Inf,-2400.0d,Inf,-3100.0d,-1600.0d,-1900.0d,-5000.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-8000.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-4600.0d,Inf,Inf,-1800.0d,-100.0d,-6000.0d,-900.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-4300.0d,Inf,Inf,Inf,Inf,-2700.0d,Inf,Inf,Inf,Inf,-6000.0d,Inf,Inf,Inf,-2500.0d,-6000.0d,-1100.0d,-3200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        -3100.0d,Inf,Inf,Inf,Inf,-1800.0d,Inf,Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,-5000.0d,-2300.0d,-3500.0d,-2400.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,-2300.0d,Inf,Inf,Inf,Inf,-700.0d,Inf,-4300.0d,-2600.0d,-3900.0d,-6000.0d,Inf,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-5000.0d,Inf,Inf,Inf,Inf,-3900.0d,Inf,Inf,-2700.0d,-2100.0d,-7000.0d,-3200.0d,Inf,Inf,Inf,-3000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-2600.0d,Inf,Inf,Inf,Inf,-2100.0d,Inf,Inf,Inf,-6000.0d,-7000.0d,-3800.0d,-3800.0d,Inf,Inf,-3900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -1600.0d,Inf,Inf,Inf,Inf,-100.0d,Inf,Inf,Inf,Inf,-6000.0d,-3900.0d,-6600.0d,-6100.0d,Inf,-2300.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,-2000.0d,Inf,-8000.0d,-5000.0d,-4300.0d,-6000.0d,Inf,Inf,Inf,Inf,-1100.0d,Inf,Inf,Inf,Inf,-3600.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-4300.0d,Inf,Inf,-3200.0d,-3900.0d,-7000.0d,-4900.0d,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,-5900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-3900.0d,Inf,Inf,Inf,-4600.0d,-7000.0d,-700.0d,-5700.0d,Inf,Inf,-3800.0d,Inf,Inf,Inf,Inf,-6600.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -1900.0d,Inf,Inf,Inf,Inf,-6000.0d,-3000.0d,-5900.0d,-7400.0d,Inf,-1100.0d,Inf,Inf,Inf,Inf,-3500.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -2500.0d,-2300.0d,-2000.0d,-5000.0d,Inf,Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,-3900.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        -2700.0d,-700.0d,-6000.0d,-2500.0d,Inf,Inf,Inf,-4900.0d,Inf,Inf,Inf,Inf,-5700.0d,Inf,Inf,Inf,Inf,-7400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -2400.0d,-6000.0d,-1100.0d,-3900.0d,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-3800.0d,Inf,Inf,Inf,Inf,-6100.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -5000.0d,-700.0d,-3600.0d,-3200.0d,Inf,-900.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-2400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] Tstack2S = {
        -1d,-1d,-1d,-6.3d,-1d,-1d,-1d,-1d,-7.0d,-1d,-1d,-1d,-1d,-5.8d,-1d,-7.8d,-4.0d,-4.4d,-13.5d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-22.5d,-1d,-1d,-1d,-1d,-7.1d,-1d,-1d,-1d,-1d,-11.4d,-1d,-1d,-3.8d,-0.5d,-16.1d,-1.7d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-10.7d,-1d,-1d,-1d,-1d,-6.0d,-1d,-1d,-1d,-1d,-15.5d,-1d,-1d,-1d,-5.9d,-16.1d,-2.1d,-8.7d,-1d,-1d,-1d,-1d,-1d,-1d,
        -7.8d,-1d,-1d,-1d,-1d,-3.8d,-1d,-1d,-1d,-1d,-5.9d,-1d,-1d,-1d,-1d,-13.6d,-6.3d,-9.4d,-6.5d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-5.9d,-1d,-1d,-1d,-1d,-1.3d,-1d,-10.7d,-5.9d,-9.6d,-16.1d,-1d,-1d,-1d,-1d,-1.2d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-13.8d,-1d,-1d,-1d,-1d,-10.6d,-1d,-1d,-6.0d,-5.1d,-19.3d,-8.0d,-1d,-1d,-1d,-7.8d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-5.9d,-1d,-1d,-1d,-1d,-5.1d,-1d,-1d,-1d,-15.5d,-19.3d,-9.5d,-9.0d,-1d,-1d,-10.6d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -4.0d,-1d,-1d,-1d,-1d,-0.5d,-1d,-1d,-1d,-1d,-16.1d,-10.6d,-18.7d,-16.9d,-1d,-6.3d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-4.7d,-1d,-22.5d,-13.8d,-11.1d,-16.1d,-1d,-1d,-1d,-1d,-2.7d,-1d,-1d,-1d,-1d,-9.8d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-11.1d,-1d,-1d,-7.1d,-10.6d,-19.3d,-13.5d,-1d,-1d,-1d,-19.2d,-1d,-1d,-1d,-1d,-16.1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-9.6d,-1d,-1d,-1d,-11.4d,-19.3d,-19.2d,-15.9d,-1d,-1d,-9.5d,-1d,-1d,-1d,-1d,-18.7d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -4.4d,-1d,-1d,-1d,-1d,-16.1d,-7.8d,-16.1d,-21.2d,-1d,-2.1d,-1d,-1d,-1d,-1d,-9.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -6.3d,-5.9d,-4.7d,-14.2d,-1d,-1d,-1d,-1d,-6.3d,-1d,-1d,-1d,-1d,-10.5d,-1d,-1d,-1d,-1d,-8.9d,-1d,-1d,-1d,-1d,-1d,-1d,
        -7.0d,-1.3d,-16.1d,-6.3d,-1d,-1d,-1d,-13.5d,-1d,-1d,-1d,-1d,-15.9d,-1d,-1d,-1d,-1d,-21.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -5.8d,-16.1d,-2.7d,-10.5d,-1d,-1d,-8.0d,-1d,-1d,-1d,-1d,-9.0d,-1d,-1d,-1d,-1d,-16.9d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -13.5d,-1.2d,-9.8d,-8.9d,-1d,-1.7d,-1d,-1d,-1d,-1d,-8.7d,-1d,-1d,-1d,-1d,-6.5d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] TstackH = {
        Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,-2700.0d,Inf,Inf,Inf,Inf,-2400.0d,Inf,-3100.0d,-1600.0d,-1900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-8000.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-4600.0d,Inf,Inf,-1800.0d,-100.0d,Inf,-900.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-4300.0d,Inf,Inf,Inf,Inf,-2700.0d,Inf,Inf,Inf,Inf,-6000.0d,Inf,Inf,Inf,-2500.0d,Inf,-1100.0d,-3200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        -3100.0d,Inf,Inf,Inf,Inf,-1800.0d,Inf,Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,Inf,-2300.0d,-3500.0d,-2400.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,-2300.0d,Inf,Inf,Inf,Inf,-700.0d,Inf,-4300.0d,-2600.0d,-3900.0d,Inf,Inf,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-5000.0d,Inf,Inf,Inf,Inf,-3900.0d,Inf,Inf,-2700.0d,-2100.0d,Inf,-3200.0d,Inf,Inf,Inf,-3000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-2600.0d,Inf,Inf,Inf,Inf,-2100.0d,Inf,Inf,Inf,-6000.0d,Inf,-3800.0d,-3800.0d,Inf,Inf,-3900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -1600.0d,Inf,Inf,Inf,Inf,-100.0d,Inf,Inf,Inf,Inf,Inf,-3900.0d,-6600.0d,-6100.0d,Inf,-2300.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,-2000.0d,Inf,-8000.0d,-5000.0d,-4300.0d,Inf,Inf,Inf,Inf,Inf,-1100.0d,Inf,Inf,Inf,Inf,-3600.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-4300.0d,Inf,Inf,-3200.0d,-3900.0d,Inf,-4900.0d,Inf,Inf,Inf,-700.0d,Inf,Inf,Inf,Inf,-5900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-3900.0d,Inf,Inf,Inf,-4600.0d,Inf,-700.0d,-5700.0d,Inf,Inf,-3800.0d,Inf,Inf,Inf,Inf,-6600.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -1900.0d,Inf,Inf,Inf,Inf,Inf,-3000.0d,-5900.0d,-7400.0d,Inf,-1100.0d,Inf,Inf,Inf,Inf,-3500.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -2500.0d,-2300.0d,-2000.0d,Inf,Inf,Inf,Inf,Inf,-2500.0d,Inf,Inf,Inf,Inf,-3900.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        -2700.0d,-700.0d,Inf,-2500.0d,Inf,Inf,Inf,-4900.0d,Inf,Inf,Inf,Inf,-5700.0d,Inf,Inf,Inf,Inf,-7400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -2400.0d,Inf,-1100.0d,-3900.0d,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-3800.0d,Inf,Inf,Inf,Inf,-6100.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-700.0d,-3600.0d,-3200.0d,Inf,-900.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-2400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] TstackS = {
        -1d,-1d,-1d,-6.3d,-1d,-1d,-1d,-1d,-7.0d,-1d,-1d,-1d,-1d,-5.8d,-1d,-7.8d,-4.0d,-4.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-22.5d,-1d,-1d,-1d,-1d,-7.1d,-1d,-1d,-1d,-1d,-11.4d,-1d,-1d,-3.8d,-0.5d,-1d,-1.7d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-10.7d,-1d,-1d,-1d,-1d,-6.0d,-1d,-1d,-1d,-1d,-15.5d,-1d,-1d,-1d,-5.9d,-1d,-2.1d,-8.7d,-1d,-1d,-1d,-1d,-1d,-1d,
        -7.8d,-1d,-1d,-1d,-1d,-3.8d,-1d,-1d,-1d,-1d,-5.9d,-1d,-1d,-1d,-1d,-1d,-6.3d,-9.4d,-6.5d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-5.9d,-1d,-1d,-1d,-1d,-1.3d,-1d,-10.7d,-5.9d,-9.6d,-1d,-1d,-1d,-1d,-1d,-1.2d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-13.8d,-1d,-1d,-1d,-1d,-10.6d,-1d,-1d,-6.0d,-5.1d,-1d,-8.0d,-1d,-1d,-1d,-7.8d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-5.9d,-1d,-1d,-1d,-1d,-5.1d,-1d,-1d,-1d,-15.5d,-1d,-9.5d,-9.0d,-1d,-1d,-10.6d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -4.0d,-1d,-1d,-1d,-1d,-0.5d,-1d,-1d,-1d,-1d,-1d,-10.6d,-18.7d,-16.9d,-1d,-6.3d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-4.7d,-1d,-22.5d,-13.8d,-11.1d,-1d,-1d,-1d,-1d,-1d,-2.7d,-1d,-1d,-1d,-1d,-9.8d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-11.1d,-1d,-1d,-7.1d,-10.6d,-1d,-13.5d,-1d,-1d,-1d,-19.2d,-1d,-1d,-1d,-1d,-16.1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-9.6d,-1d,-1d,-1d,-11.4d,-1d,-19.2d,-15.9d,-1d,-1d,-9.5d,-1d,-1d,-1d,-1d,-18.7d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -4.4d,-1d,-1d,-1d,-1d,-1d,-7.8d,-16.1d,-21.2d,-1d,-2.1d,-1d,-1d,-1d,-1d,-9.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -6.3d,-5.9d,-4.7d,-1d,-1d,-1d,-1d,-1d,-6.3d,-1d,-1d,-1d,-1d,-10.5d,-1d,-1d,-1d,-1d,-8.9d,-1d,-1d,-1d,-1d,-1d,-1d,
        -7.0d,-1.3d,-1d,-6.3d,-1d,-1d,-1d,-13.5d,-1d,-1d,-1d,-1d,-15.9d,-1d,-1d,-1d,-1d,-21.2d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -5.8d,-1d,-2.7d,-10.5d,-1d,-1d,-8.0d,-1d,-1d,-1d,-1d,-9.0d,-1d,-1d,-1d,-1d,-16.9d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1.2d,-9.8d,-8.9d,-1d,-1.7d,-1d,-1d,-1d,-1d,-8.7d,-1d,-1d,-1d,-1d,-6.5d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] Dangle3H = {
        Inf,Inf,Inf,-500.0d,Inf,Inf,Inf,Inf,4700.0d,Inf,Inf,Inf,Inf,-4100.0d,Inf,Inf,Inf,Inf,-3800.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,-5900.0d,Inf,Inf,Inf,Inf,-2600.0d,Inf,Inf,Inf,Inf,-3200.0d,Inf,Inf,Inf,Inf,-5200.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,-2100.0d,Inf,Inf,Inf,Inf,-200.0d,Inf,Inf,Inf,Inf,-3900.0d,Inf,Inf,Inf,Inf,-4400.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        -700.0d,Inf,Inf,Inf,Inf,4400.0d,Inf,Inf,Inf,Inf,-1600.0d,Inf,Inf,Inf,Inf,2900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] Dangle3S = {
        -1d,-1d,-1d,-1.1d,-1d,-1d,-1d,-1d,14.2d,-1d,-1d,-1d,-1d,-13.1d,-1d,-1d,-1d,-1d,-12.6d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-16.5d,-1d,-1d,-1d,-1d,-7.4d,-1d,-1d,-1d,-1d,-10.4d,-1d,-1d,-1d,-1d,-15.0d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-3.9d,-1d,-1d,-1d,-1d,-0.1d,-1d,-1d,-1d,-1d,-11.2d,-1d,-1d,-1d,-1d,-13.1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -0.8d,-1d,-1d,-1d,-1d,14.9d,-1d,-1d,-1d,-1d,-3.6d,-1d,-1d,-1d,-1d,10.4d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] Dangle5H = {
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-2900.0d,-4100.0d,-4200.0d,-200.0d,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,-3700.0d,-4000.0d,-3900.0d,-4900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,-6300.0d,-4400.0d,-5100.0d,-4000.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        200.0d,600.0d,-1100.0d,-6900.0d,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,
        Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf
    };
    internal static readonly double[] Dangle5S = {
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-7.6d,-13.0d,-15.0d,-0.5d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-10.0d,-11.9d,-10.9d,-13.8d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-17.1d,-12.6d,-14.0d,-10.9d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        2.3d,3.3d,-1.6d,-20.0d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,
        -1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d,-1d
    };
    internal static readonly double[] InteriorH = { Inf,Inf,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d };
    internal static readonly double[] InteriorS = { -1d,-1d,-10.31d,-11.6d,-12.89d,-14.18d,-14.83d,-15.47d,-15.79d,-15.79d,-16.26d,-16.76d,-17.15d,-17.41d,-17.74d,-18.05d,-18.34d,-18.7d,-18.96d,-19.02d,-19.25d,-19.48d,-19.7d,-19.9d,-20.31d,-20.5d,-20.68d,-20.86d,-21.03d,-21.28d };
    internal static readonly double[] BulgeH = { 0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d,0.0d };
    internal static readonly double[] BulgeS = { -12.89d,-9.35d,-9.99d,-10.31d,-10.64d,-11.28d,-11.92d,-12.57d,-13.21d,-13.86d,-14.32d,-14.5d,-14.89d,-15.47d,-15.81d,-16.12d,-16.41d,-16.76d,-17.02d,-17.08d,-17.32d,-17.55d,-17.76d,-17.97d,-18.05d,-18.24d,-18.42d,-18.6d,-18.77d,-19.02d };
}
