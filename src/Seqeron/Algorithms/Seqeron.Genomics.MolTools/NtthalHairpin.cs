namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Primer3 <c>ntthal</c> intramolecular-hairpin (monomer, <c>type==4</c>) thermodynamic dynamic
/// program — a bit-faithful port of primer3-py 2.3.1 <c>thal.c</c> (<c>thal</c> type-4 path,
/// <c>initMatrix2</c>, <c>fillMatrix2</c>, <c>maxTM2</c>, <c>CBI</c>, <c>calc_bulge_internal2</c>,
/// <c>calc_hairpin</c>, <c>RSH</c>, <c>Ss</c>/<c>Hs</c>, <c>calc_terminal_bp</c>,
/// <c>END5_1..4</c>, <c>tracebacku</c>, <c>equal</c>, <c>calcHairpin</c>, <c>drawHairpin</c>) as
/// used by <c>primer3.calc_hairpin</c>.
/// <para>
/// The model folds one oligo onto itself, scoring Watson-Crick stacks, internal mismatches/loops
/// and bulges, the hairpin-loop initiation (<c>loops</c> hairpin column), the terminal mismatch
/// (<c>tstack2</c>) of loops ≥ 4 nt, the closing-A·T penalty of 3-nt loops, the special
/// <b>triloop</b>/<b>tetraloop</b> bonuses (<c>triloop.dh/.ds</c>, <c>tetraloop.dh/.ds</c>, keyed on
/// the loop string including the closing pair) and, on the exterior loop, the A·T penalty plus
/// dangling-end / terminal-mismatch contributions (<c>END5_1..4</c>).
/// </para>
/// <para>
/// thal.c quirks reproduced deliberately (they change results for ≈5 % of random oligos):
/// <c>calc_hairpin</c> calls the dimer <c>RSH(i, j)</c> on the non-reversed oligo, i.e. with the
/// neighbours <c>s[i+1]</c>/<c>s[j+1]</c>, and that RSH keeps its running Tm at −∞ unless a
/// dangling-end branch is entered (<see cref="NtthalDimer.TerminalPair"/>);
/// <c>calc_terminal_bp</c> accepts an exterior candidate only when ΔG at the analysis temperature
/// <c>temp</c> is below the global <c>G2 = 0</c>; <c>equal()</c> uses a 1e−5 tolerance and never
/// equates non-finite values; the 1×1 internal-mismatch candidate is skipped unless its Tm exceeds
/// the cell's by at least <c>SMALL_NON_ZERO</c> (<c>DBL_EQ</c>); the traceback consults
/// <c>CBI(…, traceback = 2)</c>'s last qualifying candidate before searching a loop.
/// </para>
/// <para>
/// All stacking / terminal-mismatch / dangle / interior / bulge tables and the physical
/// constants are reused verbatim from <see cref="NtthalDimer"/> (the same primer3
/// <c>primer3_config/*.dh,*.ds</c> values). The hairpin-loop length table and the special
/// tri/tetraloop bonus tables are embedded here.
/// </para>
/// </summary>
internal static class NtthalHairpin
{
    private const double Inf = NtthalDimer.Inf;
    private const double IlAs = NtthalDimer.IlAs;
    private const double IlAh = NtthalDimer.IlAh;
    private const double MinEntropyCutoff = NtthalDimer.MinEntropyCutoff;
    private const double MinEntropy = NtthalDimer.MinEntropy;
    private const double TempKelvin = NtthalDimer.TempKelvin;
    private const double AbsoluteZero = NtthalDimer.AbsoluteZero;
    private const int MaxLoop = NtthalDimer.MaxLoop;
    private const int MinHrpnLoop = 3;                  // thal.c MIN_HRPN_LOOP
    private const double EqualTolerance = 1e-5;         // thal.c equal()
    private const double SmallNonZero = 0.000001;       // thal.c SMALL_NON_ZERO (DBL_EQ)
    private const double G2 = 0.0;                      // thal.c global G2 (calc_terminal_bp)

    // For a unimolecular structure ntthal sets dplx_init_H = 0, dplx_init_S = -1e-11, RC = 0
    // (thal.c thal(), type 4). There is NO strand-concentration term in a hairpin.
    private const double DplxInitH = 0.0;
    private const double DplxInitS = -1e-11;
    private const double Rc = 0.0;

    // thal.c isFinite == C isfinite (false for ±∞ and NaN).
    private static bool IsFinite(double x) => double.IsFinite(x);

    /// <summary>The most stable hairpin's ntthal thermodynamics (native ntthal units).</summary>
    /// <param name="DeltaH">Hairpin ΔH° in cal/mol (salt-independent).</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol), including the (N/2−1)·saltCorrection term.</param>
    /// <param name="DeltaG37">Hairpin ΔG = ΔH° − temp·ΔS° in cal/mol at the analysis temperature
    /// (310.15 K unless another <c>temp</c> is passed; primer3-py <c>temp_c</c>).</param>
    /// <param name="TmCelsius">Unimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">ntthal N/2 (half the number of paired positions counted over
    /// bp[0..len−2]); the ntthal salt-correction base-pair count, not necessarily the literal
    /// stem length.</param>
    /// <param name="AsciiStructure">thal.c <c>drawHairpin</c> lines ("SEQ\t…", "STR\t…"; primer3-py
    /// <c>ascii_structure_lines</c>) when requested, else null.</param>
    /// <param name="PairPartners">thal.c traceback <c>bp</c>: for each 0-based position the 1-based index of its
    /// partner, 0 when unpaired (the table <c>drawHairpin</c> draws).</param>
    internal readonly record struct Result(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs,
        string[]? AsciiStructure = null, int[]? PairPartners = null);

    // Hairpin-loop ΔS by size (loops.ds hairpin column, sizes 1..30). ΔH = 0 for sizes 3..30
    // (loops.dh hairpin column; sizes 1-2 are never used, MIN_HRPN_LOOP = 3). thal.c indexes
    // hairpinLoop*[loopSize - 1] and uses index 29 for loops longer than 30.
    private static readonly double[] HairpinLoopS =
    {
        -1.0, -1.0, -11.28, -11.28, -10.64, -12.89, -13.54, -13.86, -14.5, -14.83,
        -15.29, -16.12, -16.5, -16.44, -16.77, -17.08, -17.38, -17.73, -17.99, -18.37,
        -18.61, -18.84, -19.05, -19.26, -19.66, -19.85, -20.04, -20.21, -20.38, -20.31,
    };

    // table accessors (4-D flat arrays indexed [i][ii][j][jj] -> ((i*5+ii)*5+j)*5+jj) — same
    // layout NtthalDimer uses. 3-D dangle: i*25 + col*5 + col2.
    private static double T4(double[] t, int i, int ii, int j, int jj) => t[((i * 5 + ii) * 5 + j) * 5 + jj];
    private static double T3(double[] t, int i, int j, int k) => t[(i * 5 + j) * 5 + k];

    /// <summary>
    /// Runs the full ntthal hairpin DP on one oligo (5′→3′, ACGT only). Returns <c>null</c> when
    /// no hairpin can form (ntthal <c>no_structure</c> — e.g. a homopolymer).
    /// </summary>
    /// <param name="oligo">The DNA oligo (5′→3').</param>
    /// <param name="mvMolar">Monovalent cation concentration in mol/L (ntthal mv is in mM).</param>
    internal static Result? Run(string oligo, double mvMolar) => Run(oligo, mvMolar, 0.0, 0.0);

    /// <summary>
    /// Runs the ntthal hairpin DP with divalent-cation and dNTP concentrations (mol/L), which enter
    /// only through ntthal's <c>saltCorrectS</c> (see <see cref="NtthalDimer.SaltCorrectS"/>).
    /// </summary>
    /// <param name="tempKelvin">ntthal <c>temp</c> (K): enters the exterior-loop acceptance test
    /// (<c>calc_terminal_bp</c>: ΔH − temp·ΔS &lt; 0) and the reported ΔG (<c>calcHairpin</c>).</param>
    /// <param name="maxLoop">ntthal <c>maxLoop</c>: maximum internal-loop / bulge size (0–30).</param>
    /// <param name="withStructure">Also return the thal.c <c>drawHairpin</c> ASCII lines.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <param name="maxAlign">thal.h <c>THAL_MAX_ALIGN</c> (default <see cref="NtthalDimer.ThalMaxAlign"/> = 60):
    /// a compile-time guard only, so a larger value folds longer oligos with the unchanged recursions (ntthal built
    /// with <c>-DTHAL_MAX_ALIGN=…</c>). 1 ≤ maxAlign ≤ 10 000.</param>
    /// <exception cref="ArgumentException">The oligo is longer than <paramref name="maxAlign"/> nt (thal.c
    /// <c>THAL_MAX_ALIGN</c>: both "sequences" of a hairpin are the oligo).</exception>
    internal static Result? Run(
        string oligo, double mvMolar, double dvMolar, double dntpMolar,
        double tempKelvin = TempKelvin, int maxLoop = MaxLoop, bool withStructure = false,
        int maxAlign = NtthalDimer.ThalMaxAlign)
    {
        if (maxLoop < 0 || maxLoop > MaxLoop)
            throw new ArgumentOutOfRangeException(nameof(maxLoop), maxLoop, "ntthal max_loop must be in 0..30.");
        NtthalDimer.CheckMaxAlign(maxAlign);
        if (oligo.Length > maxAlign)
            throw new ArgumentException(NtthalDimer.MaxAlignMessage(maxAlign));

        int len1 = oligo.Length;
        int len2 = len1; // monomer: numSeq2 == numSeq1 (NOT reversed for type 4)

        // 1-indexed numeric sequence with N (=4) sentinels at 0 and len+1 (thal.c thal()).
        var s = new int[len1 + 2];
        s[0] = s[len1 + 1] = 4;
        for (int k = 0; k < len1; k++) s[k + 1] = NtthalDimer.Str2Int(oligo[k]);

        double saltCorrection = NtthalDimer.SaltCorrectS(mvMolar * 1000.0, dvMolar * 1000.0, dntpMolar * 1000.0);

        static int Bp(int x, int y) => NtthalDimer.Bpi[x, y];
        static double AtPenaltyH(int x, int y) => NtthalDimer.AtPenaltyHOf(x, y);
        static double AtPenaltyS(int x, int y) => NtthalDimer.AtPenaltySOf(x, y);

        var enH = new double[len1 + 2, len2 + 2];
        var enS = new double[len1 + 2, len2 + 2];

        // thal.c Ss/Hs (k == 2): inward stack of pair (i,j) on (i+1,j-1), with the range guards.
        double Ss2(int i, int j)
        {
            if (i >= j || i == len1 || j == len2 + 1) return -1.0;
            return T4(NtthalDimer.StackS, s[i], s[i + 1], s[j], s[j - 1]);
        }
        double Hs2(int i, int j)
        {
            if (i >= j || i == len1 || j == len2 + 1) return Inf;
            double v = T4(NtthalDimer.StackH, s[i], s[i + 1], s[j], s[j - 1]);
            return IsFinite(v) ? v : Inf;
        }

        // dangle / tstack2 helpers (thal.c Sd5/Hd5/Sd3/Hd3/Ststack/Htstack, all on numSeq1).
        double Sd5(int i, int j) => T3(NtthalDimer.Dangle5S, s[i], s[j], s[j - 1]);
        double Hd5(int i, int j) => T3(NtthalDimer.Dangle5H, s[i], s[j], s[j - 1]);
        double Sd3(int i, int j) => T3(NtthalDimer.Dangle3S, s[i], s[i + 1], s[j]);
        double Hd3(int i, int j) => T3(NtthalDimer.Dangle3H, s[i], s[i + 1], s[j]);
        double Ststack(int i, int j) => T4(NtthalDimer.Tstack2S, s[i], s[i + 1], s[j], s[j - 1]);
        double Htstack(int i, int j) => T4(NtthalDimer.Tstack2H, s[i], s[i + 1], s[j], s[j - 1]);

        static bool Equal(double a, double b) =>
            IsFinite(a) && IsFinite(b) && Math.Abs(a - b) < EqualTolerance;

        static double Tm(double h, double sv) => (h + DplxInitH) / (sv + DplxInitS + Rc);

        // ----- calc_hairpin (thal.c) -----
        // traceback == 0 (fill): the existing cell value is kept when it is more stable with the
        // same RSH terminal; traceback == 1: the bare hairpin value is returned.
        void CalcHairpin(int i, int j, ref double eS, ref double eH, bool traceback)
        {
            int loopSize = j - i - 1;
            if (loopSize < MinHrpnLoop) { eS = -1.0; eH = Inf; return; }
            if (loopSize <= 30) { eH = 0.0; eS = HairpinLoopS[loopSize - 1]; }
            else { eH = 0.0; eS = HairpinLoopS[29]; }

            if (loopSize > 3)
            {
                eH += T4(NtthalDimer.Tstack2H, s[i], s[i + 1], s[j], s[j - 1]);
                eS += T4(NtthalDimer.Tstack2S, s[i], s[i + 1], s[j], s[j - 1]);
            }
            else if (loopSize == 3)
            {
                eH += AtPenaltyH(s[i], s[j]);
                eS += AtPenaltyS(s[i], s[j]);
            }

            if (loopSize == 3)
            {
                if (TriloopBonus(s, i, out double th, out double ts)) { eH += th; eS += ts; }
            }
            else if (loopSize == 4 && TetraloopBonus(s, i, out double th, out double ts)) { eH += th; eS += ts; }

            if (!IsFinite(eH)) { eH = Inf; eS = -1.0; }
            if (eH > 0 && eS > 0 && (!(enH[i, j] > 0) || !(enS[i, j] > 0))) { eH = Inf; eS = -1.0; }

            var (rs, rh) = NtthalDimer.RightTerminalPair(s, s, i, j, DplxInitH, DplxInitS, Rc);
            double g1 = eH + rh - TempKelvin * (eS + rs);
            double g2 = enH[i, j] + rh - TempKelvin * (enS[i, j] + rs);
            if (g2 < g1 && !traceback) { eS = enS[i, j]; eH = enH[i, j]; }
        }

        // ----- maxTM2 (thal.c) -----
        void MaxTm2(int i, int j)
        {
            double s0 = enS[i, j], h0 = enH[i, j];
            double t0 = Tm(h0, s0);
            double s1, h1;
            if (IsFinite(enH[i, j])) { s1 = enS[i + 1, j - 1] + Ss2(i, j); h1 = enH[i + 1, j - 1] + Hs2(i, j); }
            else { s1 = -1.0; h1 = Inf; }
            double t1 = Tm(h1, s1);
            if (s1 < MinEntropyCutoff) { s1 = MinEntropy; h1 = 0.0; }
            if (s0 < MinEntropyCutoff) { s0 = MinEntropy; h0 = 0.0; }
            if (t1 > t0) { enS[i, j] = s1; enH[i, j] = h1; }
            else { enS[i, j] = s0; enH[i, j] = h0; }
        }

        // ----- calc_bulge_internal2 (thal.c) -----
        // traceback: 0 = fill, 1 = traceback loop search (bare loop term, always written),
        // 2 = CBI traceback probe (cell-inclusive, written when T1 >= T2).
        void CalcBulgeInternal2(int i, int j, int ii, int jj, ref double eS, ref double eH, int traceback)
        {
            int loopSize1 = ii - i - 1;
            int loopSize2 = j - jj - 1;
            if (loopSize1 + loopSize2 > maxLoop) { eS = -1.0; eH = Inf; return; }
            int loopSize = loopSize1 + loopSize2 - 1;
            double sv = MinEntropy, h = 0.0, t1, t2;
            bool tb = traceback != 0;
            if ((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) // bulge
            {
                if (loopSize2 == 1 || loopSize1 == 1)
                {
                    if ((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1))
                    {
                        h = NtthalDimer.BulgeH[loopSize] + T4(NtthalDimer.StackH, s[i], s[ii], s[j], s[jj]);
                        sv = NtthalDimer.BulgeS[loopSize] + T4(NtthalDimer.StackS, s[i], s[ii], s[j], s[jj]);
                    }
                }
                else
                {
                    h = NtthalDimer.BulgeH[loopSize] + AtPenaltyH(s[i], s[j]) + AtPenaltyH(s[ii], s[jj]);
                    sv = NtthalDimer.BulgeS[loopSize] + AtPenaltyS(s[i], s[j]) + AtPenaltyS(s[ii], s[jj]);
                }
                if (traceback != 1) { h += enH[ii, jj]; sv += enS[ii, jj]; }
                if (!IsFinite(h)) { h = Inf; sv = -1.0; }
                t1 = Tm(h, sv);
                t2 = Tm(enH[i, j], enS[i, j]);
                if (t1 > t2 || (tb && t1 >= t2) || traceback == 1) { eS = sv; eH = h; }
            }
            else if (loopSize1 == 1 && loopSize2 == 1) // single internal mismatch (1×1)
            {
                sv = T4(NtthalDimer.Int2S, s[i], s[i + 1], s[j], s[j - 1]) +
                     T4(NtthalDimer.Int2S, s[jj], s[jj + 1], s[ii], s[ii - 1]);
                if (traceback != 1) sv += enS[ii, jj];
                h = T4(NtthalDimer.Int2H, s[i], s[i + 1], s[j], s[j - 1]) +
                    T4(NtthalDimer.Int2H, s[jj], s[jj + 1], s[ii], s[ii - 1]);
                if (traceback != 1) h += enH[ii, jj];
                if (!IsFinite(h)) { h = Inf; sv = -1.0; }
                t1 = Tm(h, sv);
                t2 = Tm(enH[i, j], enS[i, j]);
                // DBL_EQ(T1,T2) == 2  <=>  !((T1 - T2) < SMALL_NON_ZERO)
                if ((!(t1 - t2 < SmallNonZero) || tb) && (t1 > t2 || (tb && t1 >= t2) || traceback == 1))
                {
                    eS = sv; eH = h;
                }
            }
            else // general internal loop
            {
                h = NtthalDimer.InteriorH[loopSize] + T4(NtthalDimer.TstackH, s[i], s[i + 1], s[j], s[j - 1]) +
                    T4(NtthalDimer.TstackH, s[jj], s[jj + 1], s[ii], s[ii - 1]) + IlAh * Math.Abs(loopSize1 - loopSize2);
                if (traceback != 1) h += enH[ii, jj];
                sv = NtthalDimer.InteriorS[loopSize] + T4(NtthalDimer.TstackS, s[i], s[i + 1], s[j], s[j - 1]) +
                     T4(NtthalDimer.TstackS, s[jj], s[jj + 1], s[ii], s[ii - 1]) + IlAs * Math.Abs(loopSize1 - loopSize2);
                if (traceback != 1) sv += enS[ii, jj];
                if (!IsFinite(h)) { h = Inf; sv = -1.0; }
                t1 = Tm(h, sv);
                t2 = Tm(enH[i, j], enS[i, j]);
                if (t1 > t2 || (tb && t1 >= t2) || traceback == 1) { eS = sv; eH = h; }
            }
        }

        // ----- CBI (thal.c) -----
        void Cbi(int i, int j, ref double eS, ref double eH, int traceback)
        {
            for (int d = j - i - 3; d >= MinHrpnLoop + 1 && d >= j - i - 2 - maxLoop; --d)
            {
                for (int ii = i + 1; ii < j - d && ii <= len1; ++ii)
                {
                    int jj = d + ii;
                    if (traceback == 0) { eS = -1.0; eH = Inf; }
                    if (IsFinite(enH[ii, jj]) && IsFinite(enH[i, j]))
                    {
                        CalcBulgeInternal2(i, j, ii, jj, ref eS, ref eH, traceback);
                        if (IsFinite(eH))
                        {
                            if (eS < MinEntropyCutoff) { eS = MinEntropy; eH = 0.0; }
                            if (traceback == 0) { enH[i, j] = eH; enS[i, j] = eS; }
                        }
                    }
                }
            }
        }

        // ----- initMatrix2 (thal.c) -----
        for (int i = 1; i <= len1; ++i)
            for (int j = i; j <= len2; ++j)
            {
                if (j - i < MinHrpnLoop + 1 || Bp(s[i], s[j]) == 0) { enH[i, j] = Inf; enS[i, j] = -1.0; }
                else { enH[i, j] = 0.0; enS[i, j] = MinEntropy; }
            }

        // ----- fillMatrix2 (thal.c) -----
        for (int j = 2; j <= len2; ++j)
        {
            for (int i = j - MinHrpnLoop - 1; i >= 1; --i)
            {
                if (!IsFinite(enH[i, j])) continue;
                double sh0 = -1.0, sh1 = Inf;
                MaxTm2(i, j);
                Cbi(i, j, ref sh0, ref sh1, 0);
                sh0 = -1.0; sh1 = Inf;
                CalcHairpin(i, j, ref sh0, ref sh1, traceback: false);
                if (IsFinite(sh1))
                {
                    if (sh0 < MinEntropyCutoff) { sh0 = MinEntropy; sh1 = 0.0; }
                    enS[i, j] = sh0; enH[i, j] = sh1;
                }
            }
        }

        // ----- calc_terminal_bp (thal.c) — exterior 5' loop DP -----
        var send5 = new double[len1 + 1];
        var hend5 = new double[len1 + 1];
        send5[0] = -1.0; hend5[0] = Inf;
        if (len1 >= 1) { send5[1] = -1.0; hend5[1] = Inf; }
        for (int i = 2; i <= len1; i++) { send5[i] = MinEntropy; hend5[i] = 0.0; }

        // END5_1..4 (thal.c). Returns (H_max, S_max).
        static void End5Candidate(double h, double sv, ref double hMax, ref double sMax, ref double maxTm)
        {
            if (!IsFinite(h) || h > 0 || sv > 0) { h = Inf; sv = -1.0; }
            double t = Tm(h, sv);
            if (maxTm < t && sv > MinEntropyCutoff) { hMax = h; sMax = sv; maxTm = t; }
        }
        // thal.c: the previous exterior value (k) is carried only when its Tm >= that of an empty
        // exterior (0/(0 + dplx_init_S) = -0.0); otherwise the candidate starts from 0.
        (double H, double S) Base(int k) =>
            Tm(hend5[k], send5[k]) >= Tm(0, 0) ? (hend5[k], send5[k]) : (0.0, 0.0);
        (double H, double S) End5_1(int i)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 2; ++k)
            {
                var (bh, bs) = Base(k);
                End5Candidate(bh + AtPenaltyH(s[k + 1], s[i]) + enH[k + 1, i],
                              bs + AtPenaltyS(s[k + 1], s[i]) + enS[k + 1, i], ref hMax, ref sMax, ref maxTm);
            }
            return (hMax, sMax);
        }
        (double H, double S) End5_2(int i)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
            {
                var (bh, bs) = Base(k);
                End5Candidate(bh + AtPenaltyH(s[k + 2], s[i]) + Hd5(i, k + 2) + enH[k + 2, i],
                              bs + AtPenaltyS(s[k + 2], s[i]) + Sd5(i, k + 2) + enS[k + 2, i], ref hMax, ref sMax, ref maxTm);
            }
            return (hMax, sMax);
        }
        (double H, double S) End5_3(int i)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
            {
                var (bh, bs) = Base(k);
                End5Candidate(bh + AtPenaltyH(s[k + 1], s[i - 1]) + Hd3(i - 1, k + 1) + enH[k + 1, i - 1],
                              bs + AtPenaltyS(s[k + 1], s[i - 1]) + Sd3(i - 1, k + 1) + enS[k + 1, i - 1], ref hMax, ref sMax, ref maxTm);
            }
            return (hMax, sMax);
        }
        (double H, double S) End5_4(int i)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 4; ++k)
            {
                var (bh, bs) = Base(k);
                End5Candidate(bh + AtPenaltyH(s[k + 2], s[i - 1]) + Htstack(i - 1, k + 2) + enH[k + 2, i - 1],
                              bs + AtPenaltyS(s[k + 2], s[i - 1]) + Ststack(i - 1, k + 2) + enS[k + 2, i - 1], ref hMax, ref sMax, ref maxTm);
            }
            return (hMax, sMax);
        }

        for (int i = 2; i <= len1; ++i)
        {
            var e1 = End5_1(i); var e2 = End5_2(i); var e3 = End5_3(i); var e4 = End5_4(i);
            double t1 = Tm(hend5[i - 1], send5[i - 1]);
            int max = Max5(t1, Tm(e1.H, e1.S), Tm(e2.H, e2.S), Tm(e3.H, e3.S), Tm(e4.H, e4.S));
            (double H, double S) pick = max switch { 2 => e1, 3 => e2, 4 => e3, 5 => e4, _ => (hend5[i - 1], send5[i - 1]) };
            if (max != 1 && !(pick.H - tempKelvin * pick.S < G2))
                pick = (hend5[i - 1], send5[i - 1]);
            send5[i] = pick.S; hend5[i] = pick.H;
        }

        double mh = hend5[len1];
        double ms = send5[len1];
        if (!IsFinite(mh)) return null; // ntthal no_structure

        // ----- tracebacku (thal.c) -----
        var bp = new int[len1];
        var stack = new System.Collections.Generic.Stack<(int i, int j, int mtrx)>();
        stack.Push((len1, 0, 1));
        while (stack.Count > 0)
        {
            var (i, j, mtrx) = stack.Pop();
            if (mtrx == 1)
            {
                while (i >= 1 && Equal(send5[i], send5[i - 1]) && Equal(hend5[i], hend5[i - 1])) --i;
                if (i == 0) continue;
                var e1 = End5_1(i);
                if (Equal(send5[i], e1.S) && Equal(hend5[i], e1.H))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 2; ++k)
                    {
                        double ds = AtPenaltyS(s[k + 1], s[i]) + enS[k + 1, i];
                        double dh = AtPenaltyH(s[k + 1], s[i]) + enH[k + 1, i];
                        if (Equal(send5[i], ds) && Equal(hend5[i], dh)) { stack.Push((k + 1, i, 0)); break; }
                        if (Equal(send5[i], send5[k] + ds) && Equal(hend5[i], hend5[k] + dh))
                        { stack.Push((k + 1, i, 0)); stack.Push((k, 0, 1)); break; }
                    }
                    continue;
                }
                var e2 = End5_2(i);
                if (Equal(send5[i], e2.S) && Equal(hend5[i], e2.H))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
                    {
                        double ds = AtPenaltyS(s[k + 2], s[i]) + Sd5(i, k + 2) + enS[k + 2, i];
                        double dh = AtPenaltyH(s[k + 2], s[i]) + Hd5(i, k + 2) + enH[k + 2, i];
                        if (Equal(send5[i], ds) && Equal(hend5[i], dh)) { stack.Push((k + 2, i, 0)); break; }
                        if (Equal(send5[i], send5[k] + ds) && Equal(hend5[i], hend5[k] + dh))
                        { stack.Push((k + 2, i, 0)); stack.Push((k, 0, 1)); break; }
                    }
                    continue;
                }
                var e3 = End5_3(i);
                if (Equal(send5[i], e3.S) && Equal(hend5[i], e3.H))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
                    {
                        double ds = AtPenaltyS(s[k + 1], s[i - 1]) + Sd3(i - 1, k + 1) + enS[k + 1, i - 1];
                        double dh = AtPenaltyH(s[k + 1], s[i - 1]) + Hd3(i - 1, k + 1) + enH[k + 1, i - 1];
                        if (Equal(send5[i], ds) && Equal(hend5[i], dh)) { stack.Push((k + 1, i - 1, 0)); break; }
                        if (Equal(send5[i], send5[k] + ds) && Equal(hend5[i], hend5[k] + dh))
                        { stack.Push((k + 1, i - 1, 0)); stack.Push((k, 0, 1)); break; }
                    }
                    continue;
                }
                var e4 = End5_4(i);
                if (Equal(send5[i], e4.S) && Equal(hend5[i], e4.H))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 4; ++k)
                    {
                        double ds = AtPenaltyS(s[k + 2], s[i - 1]) + Ststack(i - 1, k + 2) + enS[k + 2, i - 1];
                        double dh = AtPenaltyH(s[k + 2], s[i - 1]) + Htstack(i - 1, k + 2) + enH[k + 2, i - 1];
                        if (Equal(send5[i], ds) && Equal(hend5[i], dh)) { stack.Push((k + 2, i - 1, 0)); break; }
                        if (Equal(send5[i], send5[k] + ds) && Equal(hend5[i], hend5[k] + dh))
                        { stack.Push((k + 2, i - 1, 0)); stack.Push((k, 0, 1)); break; }
                    }
                }
            }
            else // mtrx == 0: pair (i, j)
            {
                bp[i - 1] = j;
                bp[j - 1] = i;
                double sh1S = -1.0, sh1H = Inf;
                CalcHairpin(i, j, ref sh1S, ref sh1H, traceback: true);
                double sh2S = -1.0, sh2H = Inf;
                Cbi(i, j, ref sh2S, ref sh2H, 2);
                if (Equal(enS[i, j], Ss2(i, j) + enS[i + 1, j - 1]) &&
                    Equal(enH[i, j], Hs2(i, j) + enH[i + 1, j - 1]))
                {
                    stack.Push((i + 1, j - 1, 0));
                }
                else if (Equal(enS[i, j], sh1S) && Equal(enH[i, j], sh1H))
                {
                    // hairpin loop closed by (i, j): end of this branch
                }
                else if (Equal(enS[i, j], sh2S) && Equal(enH[i, j], sh2H))
                {
                    bool done = false;
                    for (int d = j - i - 3; d >= MinHrpnLoop + 1 && d >= j - i - 2 - maxLoop && !done; --d)
                    {
                        for (int ii = i + 1; ii < j - d; ++ii)
                        {
                            int jj = d + ii;
                            double es = -1.0, eh = Inf;
                            CalcBulgeInternal2(i, j, ii, jj, ref es, ref eh, 1);
                            if (Equal(enS[i, j], es + enS[ii, jj]) && Equal(enH[i, j], eh + enH[ii, jj]))
                            {
                                stack.Push((ii, jj, 0));
                                done = true;
                                break;
                            }
                        }
                    }
                }
            }
        }

        // ----- calcHairpin / drawHairpin (thal.c): N over bp[0 .. len1-2] -----
        int n = 0;
        for (int i = 1; i < len1; i++) if (bp[i - 1] > 0) n++;
        int half = n / 2; // integer division, as in thal.c
        double dsOut = ms + (half - 1) * saltCorrection;
        double tm = mh / dsOut - AbsoluteZero;
        double dg = mh - tempKelvin * dsOut;
        string[]? structure = withStructure ? DrawHairpin(oligo.ToUpperInvariant(), bp) : null;
        return new Result(mh, dsOut, dg, tm, half, structure, bp);
    }

    /// <summary>
    /// thal.c <c>drawHairpin</c> ASCII lines (primer3-py <c>ascii_structure_lines</c>): "SEQ\t" + one
    /// character per base ('-' unpaired; for a pair (i, j), i &lt; j, '/' is written at i and '\' at
    /// j) and "STR\t" + the oligo.
    /// </summary>
    private static string[] DrawHairpin(string oligo, int[] bp)
    {
        int len1 = oligo.Length;
        var row = new char[len1];
        for (int i = 1; i < len1 + 1; ++i)
        {
            if (bp[i - 1] == 0) row[i - 1] = '-';
            else if (bp[i - 1] > i - 1) row[bp[i - 1] - 1] = '\\';
            else row[bp[i - 1] - 1] = '/';
        }
        return new[] { "SEQ\t" + new string(row), "STR\t" + oligo };
    }

    // thal.c max5, verbatim: 1 only when T1 is strictly greater than all others; otherwise the first
    // of T2..T4 strictly greater than every LATER value, else 5 (NaN comparisons are false).
    private static int Max5(double a, double b, double c, double d, double e)
    {
        if (a > b && a > c && a > d && a > e) return 1;
        if (b > c && b > d && b > e) return 2;
        if (c > d && c > e) return 3;
        if (d > e) return 4;
        return 5;
    }

    // ----- special-loop bonus lookups (verbatim primer3 triloop/tetraloop tables, below) -----
    // The key is the full loop string INCLUDING the closing base pair: for a triloop the 5 bases
    // s[i],s[i+1],s[i+2],s[i+3],s[i+4] (closing-5' + 3 loop + closing-3'); tetraloop = 6 bases.
    private static bool TriloopBonus(int[] s, int i, out double dh, out double ds)
    {
        int key = Enc(s, i, 5); // closing-5' + 3 loop + closing-3'
        for (int k = 0; k < TriloopKeys.Length; k++)
        {
            if (TriloopKeys[k] == key) { dh = TriloopDh[k]; ds = TriloopDs[k]; return true; }
        }
        dh = 0.0; ds = 0.0; return false;
    }

    private static bool TetraloopBonus(int[] s, int i, out double dh, out double ds)
    {
        int key = Enc(s, i, 6); // closing-5' + 4 loop + closing-3'
        for (int k = 0; k < TetraloopKeys.Length; k++)
        {
            if (TetraloopKeys[k] == key) { dh = TetraloopDh[k]; ds = TetraloopDs[k]; return true; }
        }
        dh = 0.0; ds = 0.0; return false;
    }

    // Encode a loop string of <paramref name="n"/> bases (A=0,C=1,G=2,T=3) as a 4-bit-nibble int.
    private static int Enc(string loop)
    {
        int v = 0;
        foreach (char c in loop) v = (v << 4) | NtthalDimer.Str2Int(c);
        return v;
    }

    private static int Enc(int[] s, int i, int n)
    {
        int v = 0;
        for (int k = 0; k < n; k++) v = (v << 4) | s[i + k];
        return v;
    }

    // ===================== special-loop bonus tables (verbatim primer3) =====================
    // Triloop bonuses: primer3_config/triloop.dh (ΔH, cal/mol) + triloop.ds (ΔS, cal/(K·mol), all 0).
    // Tetraloop bonuses: primer3_config/tetraloop.dh (ΔH) + tetraloop.ds (ΔS). Keys are the full
    // loop string incl. the closing base pair (5-char triloop, 6-char tetraloop). Provenance: the
    // libprimer3 thermodynamic parameter files vendored in primer3-py (GPL-2.0); values transcribed
    // verbatim. Source: SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 special hairpin loops.
    private static readonly (string Loop, double Dh, double Ds)[] TriloopTable =
    {
        ("AGAAT", -1500, 0), ("AGCAT", -1500, 0), ("AGGAT", -1500, 0), ("AGTAT", -1500, 0),
        ("CGAAG", -2000, 0), ("CGCAG", -2000, 0), ("CGGAG", -2000, 0), ("CGTAG", -2000, 0),
        ("GGAAC", -2000, 0), ("GGCAC", -2000, 0), ("GGGAC", -2000, 0), ("GGTAC", -2000, 0),
        ("TGAAA", -1500, 0), ("TGCAA", -1500, 0), ("TGGAA", -1500, 0), ("TGTAA", -1500, 0),
    };

    private static readonly (string Loop, double Dh, double Ds)[] TetraloopTable =
    {
        ("AAAAAT", 500, -650), ("AAAACT", 700, 1610), ("AAACAT", 1000, 1610), ("ACTTGT", 0, 4190),
        ("AGAAAT", -1100, 1610), ("AGAGAT", -1100, 1610), ("AGATAT", -1500, 1610), ("AGCAAT", -1600, 1610),
        ("AGCGAT", -1100, 1610), ("AGCTTT", 200, 1610), ("AGGAAT", -1100, 1610), ("AGGGAT", -1100, 1610),
        ("AGGGGT", 500, 640), ("AGTAAT", -1600, 1610), ("AGTGAT", -1100, 1610), ("AGTTCT", 800, 1610),
        ("ATTCGT", -200, 1610), ("ATTTGT", 0, 1610), ("ATTTTT", -500, 1610), ("CAAAAG", 500, -1290),
        ("CAAACG", 700, 0), ("CAACAG", 1000, 0), ("CAACCG", 0, 0), ("CCTTGG", 0, 2570),
        ("CGAAAG", -1100, 0), ("CGAGAG", -1100, 0), ("CGATAG", -1500, 0), ("CGCAAG", -1600, 0),
        ("CGCGAG", -1100, 0), ("CGCTTG", 200, 0), ("CGGAAG", -1100, 0), ("CGGGAG", -1000, 0),
        ("CGGGGG", 500, -970), ("CGTAAG", -1600, 0), ("CGTGAG", -1100, 0), ("CGTTCG", 800, 0),
        ("CTTCGG", -200, 0), ("CTTTGG", 0, 0), ("CTTTTG", -500, 0), ("GAAAAC", 500, -3230),
        ("GAAACC", 700, 0), ("GAACAC", 1000, 0), ("GCTTGC", 0, 2570), ("GGAAAC", -1100, 0),
        ("GGAGAC", -1100, 0), ("GGATAC", -1600, 0), ("GGCAAC", -1600, 0), ("GGCGAC", -1100, 0),
        ("GGCTTC", 200, 0), ("GGGAAC", -1100, 0), ("GGGGAC", -1100, 0), ("GGGGGC", 500, -970),
        ("GGTAAC", -1600, 0), ("GGTGAC", -1100, 0), ("GGTTCC", 800, 0), ("GTTCGC", -200, 0),
        ("GTTTGC", 0, 0), ("GTTTTC", -500, 0), ("TAAAAA", 500, 320), ("TAAACA", 700, 1610),
        ("TAACAA", 1000, 1610), ("TCTTGA", 0, 4190), ("TGAAAA", -1100, 1610), ("TGAGAA", -1100, 1610),
        ("TGATAA", -1600, 1610), ("TGCAAA", -1600, 1610), ("TGCGAA", -1100, 1610), ("TGCTTA", 200, 1610),
        ("TGGAAA", -1100, 1610), ("TGGGAA", -1100, 1610), ("TGGGGA", 500, 640), ("TGTAAA", -1600, 1610),
        ("TGTGAA", -1100, 1610), ("TGTTCA", 800, 1610), ("TTTCGA", -200, 1610), ("TTTTGA", 0, 1610),
        ("TTTTTA", -500, 1610),
    };

    private static readonly int[] TriloopKeys = TriloopTable.Select(t => Enc(t.Loop)).ToArray();
    private static readonly double[] TriloopDh = TriloopTable.Select(t => t.Dh).ToArray();
    private static readonly double[] TriloopDs = TriloopTable.Select(t => t.Ds).ToArray();
    private static readonly int[] TetraloopKeys = TetraloopTable.Select(t => Enc(t.Loop)).ToArray();
    private static readonly double[] TetraloopDh = TetraloopTable.Select(t => t.Dh).ToArray();
    private static readonly double[] TetraloopDs = TetraloopTable.Select(t => t.Ds).ToArray();
}
