using System;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Full Primer3 <c>ntthal</c> intramolecular-hairpin (monomer, <c>type==4</c>) thermodynamic
/// dynamic program — a line-for-line translation of <c>thal.c</c> (<c>initMatrix2</c>,
/// <c>fillMatrix2</c>, <c>maxTM2</c>, <c>CBI</c>/<c>calc_bulge_internal2</c>, <c>calc_hairpin</c>,
/// <c>calc_terminal_bp</c>, <c>END5_1..4</c>, <c>tracebacku</c>, <c>calcHairpin</c>) for the
/// hairpin (single-sequence) folding path used by <c>primer3.calc_hairpin</c>.
/// <para>
/// The model folds one oligo onto itself, scoring a single Watson-Crick stem, internal
/// mismatches/loops and bulges, the terminal-mismatch (<c>tstack2</c>) / dangling-end
/// (<c>dangle</c>) terminal contributions, the size-keyed hairpin-loop initiation
/// (<c>loops</c>), the closing-A·T penalty (3-nt loops) and — the addition this file bundles —
/// the sequence-specific special <b>triloop</b> (3-nt) and <b>tetraloop</b> (4-nt) stability
/// bonuses (the primer3 <c>triloop.dh/.ds</c> + <c>tetraloop.dh/.ds</c> tables), keyed on the
/// full loop string including the closing base pair.
/// </para>
/// <para>
/// All stacking / terminal-mismatch / dangle / interior / bulge tables and the physical
/// constants are reused verbatim from <see cref="NtthalDimer"/> (the same primer3
/// <c>primer3_config/*.dh,*.ds</c> values <c>ntthal</c> loads). The hairpin-loop length table
/// (<c>loops</c> hairpin column) and the special tri/tetraloop bonus tables are embedded here.
/// </para>
/// </summary>
internal static class NtthalHairpin
{
    private const double Inf = NtthalDimer.Inf;
    private const double R = NtthalDimer.R;
    private const double IlAs = NtthalDimer.IlAs;
    private const double IlAh = NtthalDimer.IlAh;
    private const double AtH = NtthalDimer.AtH;
    private const double AtS = NtthalDimer.AtS;
    private const double AtPenaltySEntry = NtthalDimer.AtPenaltySEntry;
    private const double MinEntropyCutoff = NtthalDimer.MinEntropyCutoff;
    private const double MinEntropy = NtthalDimer.MinEntropy;
    private const double TempKelvin = NtthalDimer.TempKelvin;
    private const double AbsoluteZero = NtthalDimer.AbsoluteZero;
    private const int MaxLoop = NtthalDimer.MaxLoop;
    private const int MinHrpnLoop = 3;                  // thal.c MIN_HRPN_LOOP
    private const double Equal = 1e-6;                  // traceback equality tolerance

    // For a unimolecular structure ntthal sets dplx_init_H = 0, dplx_init_S = -1e-11, RC = 0
    // (thal.c lines 583-585). There is NO strand-concentration term in a hairpin.
    private const double DplxInitH = 0.0;
    private const double DplxInitS = -1e-11;
    private const double Rc = 0.0;

    private static bool IsFinite(double x) => NtthalDimer.IsFinite(x);

    /// <summary>The most stable hairpin's ntthal thermodynamics (native ntthal units).</summary>
    /// <param name="DeltaH">Hairpin ΔH° in cal/mol (salt-independent).</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol), including the (N/2−1)·saltCorrection term.</param>
    /// <param name="DeltaG37">Hairpin ΔG°37 = ΔH° − 310.15·ΔS° in cal/mol.</param>
    /// <param name="TmCelsius">Unimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">ntthal N/2 (half the number of paired positions counted over
    /// bp[0..len−2]); the ntthal salt-correction base-pair count, not necessarily the literal
    /// stem length.</param>
    internal readonly record struct Result(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

    // Hairpin-loop ΔS by size (loops.ds hairpin column, sizes 1..30). ΔH = 0 for all sizes
    // (loops.dh hairpin column). thal.c indexes hairpinLoop*[loopSize - 1].
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
    internal static Result? Run(string oligo, double mvMolar)
    {
        double mv = mvMolar * 1000.0;
        int len1 = oligo.Length;
        int len2 = len1; // monomer: oligo1 == oligo2 (NOT reversed for type==4)

        // 1-indexed sequence with N (=4) sentinels at 0 and len+1 (thal.c line 634).
        var s = new int[len1 + 2];
        s[0] = s[len1 + 1] = 4;
        for (int k = 0; k < len1; k++) s[k + 1] = NtthalDimer.Str2Int(oligo[k]);

        // saltCorrectS (thal.c 1042); dv = dntp = 0 so the divalent term vanishes.
        double saltCorrection = 0.368 * Math.Log(mv / 1000.0);

        int Bp(int x, int y) => NtthalDimer.Bpi[x, y];
        double AtPenaltyH(int x, int y) => (x == 0 && y == 3) || (x == 3 && y == 0) ? AtH : 0.0;
        double AtPenaltyS(int x, int y) => (x == 0 && y == 3) || (x == 3 && y == 0) ? AtS : AtPenaltySEntry;

        var enH = new double[len1 + 2, len2 + 2];
        var enS = new double[len1 + 2, len2 + 2];

        // Stack ΔS/ΔH for the inward-extension step (Ss/Hs mode k==2): numSeq2 == s, so
        // stack[ s[i] ][ s[i+1] ][ s[j] ][ s[j-1] ] (thal.c 2002/2026).
        double Ss2(int i, int j) => T4(NtthalDimer.StackS, s[i], s[i + 1], s[j], s[j - 1]);
        double Hs2(int i, int j) => T4(NtthalDimer.StackH, s[i], s[i + 1], s[j], s[j - 1]);

        // dangle / tstack2 helpers (thal.c Sd5/Hd5/Sd3/Hd3/Ststack/Htstack). numSeq2 == s.
        double Sd5(int i, int j) => T3(NtthalDimer.Dangle5S, s[i], s[j], s[j - 1]);
        double Hd5(int i, int j) => T3(NtthalDimer.Dangle5H, s[i], s[j], s[j - 1]);
        double Sd3(int i, int j) => T3(NtthalDimer.Dangle3S, s[i], s[i + 1], s[j]);
        double Hd3(int i, int j) => T3(NtthalDimer.Dangle3H, s[i], s[i + 1], s[j]);
        double Ststack(int i, int j) => T4(NtthalDimer.Tstack2S, s[i], s[i + 1], s[j], s[j - 1]);
        double Htstack(int i, int j) => T4(NtthalDimer.Tstack2H, s[i], s[i + 1], s[j], s[j - 1]);

        // ----- calc_hairpin (thal.c 2067-2146) -----
        // Returns (S,H) for the hairpin closed by pair (i,j); loop = j-i-1 unpaired bases.
        (double S, double H) CalcHairpin(int i, int j)
        {
            int loopSize = j - i - 1;
            if (loopSize < MinHrpnLoop) return (-1.0, Inf);

            double h, eS;
            if (loopSize <= 30) { h = 0.0; eS = HairpinLoopS[loopSize - 1]; }
            else { h = 0.0; eS = HairpinLoopS[29]; }

            if (loopSize > 3)
            {
                // terminal mismatch (tstack2) for loops of 4+ (thal.c 2099-2100).
                h += T4(NtthalDimer.Tstack2H, s[i], s[i + 1], s[j], s[j - 1]);
                eS += T4(NtthalDimer.Tstack2S, s[i], s[i + 1], s[j], s[j - 1]);
            }
            else // loopSize == 3: closing-A·T penalty (thal.c 2102-2103).
            {
                h += AtPenaltyH(s[i], s[j]);
                eS += AtPenaltyS(s[i], s[j]);
            }

            if (loopSize == 3) // triloop bonus, keyed on s[i..i+4] (closing pair + 3 loop) (thal.c 2106-2117).
            {
                if (TriloopBonus(s, i, out double th, out double ts)) { h += th; eS += ts; }
            }
            else if (loopSize == 4) // tetraloop bonus, keyed on s[i..i+5] (thal.c 2118-2127).
            {
                if (TetraloopBonus(s, i, out double th, out double ts)) { h += th; eS += ts; }
            }

            if (!IsFinite(h)) { h = Inf; eS = -1.0; }
            // both S and H positive and the cell isn't already positive → reject (thal.c 2133-2136).
            if (h > 0 && eS > 0 && (!(enH[i, j] > 0) || !(enS[i, j] > 0))) { h = Inf; eS = -1.0; }

            // RSH for monomer: 3'-side terminal stack of the closing pair.
            var (rs, rh) = Rsh(i, j);
            double g1 = h + rh - TempKelvin * (eS + rs);
            double g2 = enH[i, j] + rh - TempKelvin * (enS[i, j] + rs);
            if (g2 < g1) return (enS[i, j], enH[i, j]); // keep the existing (stack-extension) value
            return (eS, h);
        }

        // ----- RSH (thal.c 1857-1983) right terminal stack (3'-side) for the closing pair (i,j) -----
        (double S, double H) Rsh(int i, int j)
        {
            if (Bp(s[i], s[j]) == 0) return (-1.0, Inf);
            double s1 = AtPenaltyS(s[i], s[j]) + T4(NtthalDimer.Tstack2S, s[i], s[i + 1], s[j], s[j - 1]);
            double h1 = AtPenaltyH(s[i], s[j]) + T4(NtthalDimer.Tstack2H, s[i], s[i + 1], s[j], s[j - 1]);
            double g1 = h1 - TempKelvin * s1;
            if (!IsFinite(h1) || g1 > 0) { h1 = Inf; s1 = -1.0; g1 = 1.0; }

            bool unpaired = Bp(s[i + 1], s[j - 1]) == 0;
            bool d3 = IsFinite(Hd3(i, j));
            bool d5 = IsFinite(Hd5(i, j));
            double s2, h2, g2;
            if (unpaired && d3 && d5)
            {
                s2 = AtPenaltyS(s[i], s[j]) + Sd3(i, j) + Sd5(i, j);
                h2 = AtPenaltyH(s[i], s[j]) + Hd3(i, j) + Hd5(i, j);
                g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                if (IsFinite(h1) && g1 < 0)
                {
                    if (h1 - TempKelvin * s1 < (h2 - TempKelvin * s2)) { /* keep 1 */ }
                    else if (g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (unpaired && d3)
            {
                s2 = AtPenaltyS(s[i], s[j]) + Sd3(i, j);
                h2 = AtPenaltyH(s[i], s[j]) + Hd3(i, j);
                g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                if (IsFinite(h1) && g1 < 0)
                {
                    if (h1 - TempKelvin * s1 < (h2 - TempKelvin * s2)) { /* keep 1 */ }
                    else if (g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else if (unpaired && d5)
            {
                s2 = AtPenaltyS(s[i], s[j]) + Sd5(i, j);
                h2 = AtPenaltyH(s[i], s[j]) + Hd5(i, j);
                g2 = h2 - TempKelvin * s2;
                if (!IsFinite(h2) || g2 > 0) { h2 = Inf; s2 = -1.0; g2 = 1.0; }
                if (IsFinite(h1) && g1 < 0)
                {
                    if (h1 - TempKelvin * s1 < (h2 - TempKelvin * s2)) { /* keep 1 */ }
                    else if (g2 < 0) { s1 = s2; h1 = h2; }
                }
                else if (g2 < 0) { s1 = s2; h1 = h2; }
            }
            else
            {
                s2 = AtPenaltyS(s[i], s[j]);
                h2 = AtPenaltyH(s[i], s[j]);
                if (IsFinite(h1) && g1 < 0) { /* keep 1 */ }
                else { s1 = s2; h1 = h2; }
            }
            return (s1, h1);
        }

        // ----- maxTM2 (thal.c 1704-1738) -----
        void MaxTm2(int i, int j)
        {
            double s0 = enS[i, j], h0 = enH[i, j];
            double t0 = (h0 + DplxInitH) / (s0 + DplxInitS + Rc);
            double s1, h1;
            if (IsFinite(enH[i, j])) { s1 = enS[i + 1, j - 1] + Ss2(i, j); h1 = enH[i + 1, j - 1] + Hs2(i, j); }
            else { s1 = -1.0; h1 = Inf; }
            double t1 = (h1 + DplxInitH) / (s1 + DplxInitS + Rc);
            if (s1 < MinEntropyCutoff) { s1 = MinEntropy; h1 = 0.0; }
            if (s0 < MinEntropyCutoff) { s0 = MinEntropy; h0 = 0.0; }
            if (t1 > t0) { enS[i, j] = s1; enH[i, j] = h1; }
            else { enS[i, j] = s0; enH[i, j] = h0; }
        }

        // ----- calc_bulge_internal2 (thal.c 2307-2472) -----
        // For traceback==false the cell (i,j) is updated only when the candidate has higher Tm.
        (double S, double H)? CalcBulgeInternal2(int i, int j, int ii, int jj, bool traceback)
        {
            int loopSize1 = ii - i - 1;
            int loopSize2 = j - jj - 1;
            if (loopSize1 + loopSize2 > MaxLoop) return null;
            int loopSize = loopSize1 + loopSize2 - 1;
            double h, eS, t1, t2;

            if ((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) // bulge
            {
                if (loopSize2 == 1 || loopSize1 == 1)
                {
                    h = Inf; eS = MinEntropy;
                    if ((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1))
                    {
                        h = NtthalDimer.BulgeH[loopSize] + T4(NtthalDimer.StackH, s[i], s[ii], s[j], s[jj]);
                        eS = NtthalDimer.BulgeS[loopSize] + T4(NtthalDimer.StackS, s[i], s[ii], s[j], s[jj]);
                    }
                    if (!traceback) { h += enH[ii, jj]; eS += enS[ii, jj]; }
                    if (!IsFinite(h)) { h = Inf; eS = -1.0; }
                    t1 = (h + DplxInitH) / (eS + DplxInitS + Rc);
                    t2 = (enH[i, j] + DplxInitH) / (enS[i, j] + DplxInitS + Rc);
                    if (t1 > t2 || traceback) return (eS, h);
                }
                else
                {
                    h = NtthalDimer.BulgeH[loopSize] + AtPenaltyH(s[i], s[j]) + AtPenaltyH(s[ii], s[jj]);
                    eS = NtthalDimer.BulgeS[loopSize] + AtPenaltyS(s[i], s[j]) + AtPenaltyS(s[ii], s[jj]);
                    if (!traceback) { h += enH[ii, jj]; eS += enS[ii, jj]; }
                    if (!IsFinite(h)) { h = Inf; eS = -1.0; }
                    t1 = (h + DplxInitH) / (eS + DplxInitS + Rc);
                    t2 = (enH[i, j] + DplxInitH) / (enS[i, j] + DplxInitS + Rc);
                    if (t1 > t2 || traceback) return (eS, h);
                }
            }
            else if (loopSize1 == 1 && loopSize2 == 1) // single internal mismatch (1×1)
            {
                eS = T4(NtthalDimer.Int2S, s[i], s[i + 1], s[j], s[j - 1]) +
                     T4(NtthalDimer.Int2S, s[jj], s[jj + 1], s[ii], s[ii - 1]);
                h = T4(NtthalDimer.Int2H, s[i], s[i + 1], s[j], s[j - 1]) +
                    T4(NtthalDimer.Int2H, s[jj], s[jj + 1], s[ii], s[ii - 1]);
                if (!traceback) { h += enH[ii, jj]; eS += enS[ii, jj]; }
                if (!IsFinite(h)) { h = Inf; eS = -1.0; }
                t1 = (h + DplxInitH) / (eS + DplxInitS + Rc);
                t2 = (enH[i, j] + DplxInitH) / (enS[i, j] + DplxInitS + Rc);
                if (t1 > t2 || traceback) return (eS, h);
                return null;
            }
            else // general internal loop
            {
                h = NtthalDimer.InteriorH[loopSize] + T4(NtthalDimer.TstackH, s[i], s[i + 1], s[j], s[j - 1]) +
                    T4(NtthalDimer.TstackH, s[jj], s[jj + 1], s[ii], s[ii - 1]) + IlAh * Math.Abs(loopSize1 - loopSize2);
                eS = NtthalDimer.InteriorS[loopSize] + T4(NtthalDimer.TstackS, s[i], s[i + 1], s[j], s[j - 1]) +
                     T4(NtthalDimer.TstackS, s[jj], s[jj + 1], s[ii], s[ii - 1]) + IlAs * Math.Abs(loopSize1 - loopSize2);
                if (!traceback) { h += enH[ii, jj]; eS += enS[ii, jj]; }
                if (!IsFinite(h)) { h = Inf; eS = -1.0; }
                t1 = (h + DplxInitH) / (eS + DplxInitS + Rc);
                t2 = (enH[i, j] + DplxInitH) / (enS[i, j] + DplxInitS + Rc);
                if (t1 > t2 || traceback) return (eS, h);
            }
            return null;
        }

        // ----- CBI (thal.c 2036-2064) -----
        void Cbi(int i, int j)
        {
            for (int d = j - i - 3; d >= MinHrpnLoop + 1 && d >= j - i - 2 - MaxLoop; --d)
            {
                for (int ii = i + 1; ii < j - d && ii <= len1; ++ii)
                {
                    int jj = d + ii;
                    if (IsFinite(enH[ii, jj]) && IsFinite(enH[i, j]))
                    {
                        var r = CalcBulgeInternal2(i, j, ii, jj, traceback: false);
                        if (r is { } rr && IsFinite(rr.H))
                        {
                            double cs = rr.S, ch = rr.H;
                            if (cs < MinEntropyCutoff) { cs = MinEntropy; ch = 0.0; }
                            enH[i, j] = ch; enS[i, j] = cs;
                        }
                    }
                }
            }
        }

        // ----- initMatrix2 (thal.c 1565-1580) -----
        for (int i = 1; i <= len1; ++i)
            for (int j = i; j <= len2; ++j)
            {
                if (j - i < MinHrpnLoop + 1 || Bp(s[i], s[j]) == 0) { enH[i, j] = Inf; enS[i, j] = -1.0; }
                else { enH[i, j] = 0.0; enS[i, j] = MinEntropy; }
            }

        // ----- fillMatrix2 (thal.c 1631-1659) -----
        for (int j = 2; j <= len2; ++j)
        {
            for (int i = j - MinHrpnLoop - 1; i >= 1; --i)
            {
                if (!IsFinite(enH[i, j])) continue;
                MaxTm2(i, j);
                Cbi(i, j);
                var (sh0, sh1) = CalcHairpin(i, j);
                if (IsFinite(sh1))
                {
                    if (sh0 < MinEntropyCutoff) { sh0 = MinEntropy; sh1 = 0.0; }
                    enS[i, j] = sh0; enH[i, j] = sh1;
                }
            }
        }

        // ----- calc_terminal_bp (thal.c 2475-2551) — exterior 5' loop DP -----
        var send5 = new double[len1 + 1];
        var hend5 = new double[len1 + 1];
        double Send5(int i) => send5[i];
        double Hend5(int i) => hend5[i];
        send5[0] = -1.0; hend5[0] = Inf;
        if (len1 >= 1) { send5[1] = -1.0; hend5[1] = Inf; }
        for (int i = 2; i <= len1; i++) { send5[i] = MinEntropy; hend5[i] = 0.0; }

        // END5_1..4 (thal.c 2553-2730). hs: 1 = enthalpy, 2 = entropy.
        double End5_1(int i, int hs)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 2; ++k)
            {
                double t1 = (Hend5(k) + DplxInitH) / (Send5(k) + DplxInitS + Rc);
                double t2 = (0 + DplxInitH) / (0 + DplxInitS + Rc);
                double h, es;
                double bh = (t1 >= t2) ? Hend5(k) : 0.0;
                double bs = (t1 >= t2) ? Send5(k) : 0.0;
                h = bh + AtPenaltyH(s[k + 1], s[i]) + enH[k + 1, i];
                es = bs + AtPenaltyS(s[k + 1], s[i]) + enS[k + 1, i];
                if (!IsFinite(h) || h > 0 || es > 0) { h = Inf; es = -1.0; }
                double t = (h + DplxInitH) / (es + DplxInitS + Rc);
                if (maxTm < t && es > MinEntropyCutoff) { hMax = h; sMax = es; maxTm = t; }
            }
            return hs == 1 ? hMax : sMax;
        }
        double End5_2(int i, int hs)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
            {
                double t1 = (Hend5(k) + DplxInitH) / (Send5(k) + DplxInitS + Rc);
                double t2 = (0 + DplxInitH) / (0 + DplxInitS + Rc);
                double bh = (t1 >= t2) ? Hend5(k) : 0.0;
                double bs = (t1 >= t2) ? Send5(k) : 0.0;
                double h = bh + AtPenaltyH(s[k + 2], s[i]) + Hd5(i, k + 2) + enH[k + 2, i];
                double es = bs + AtPenaltyS(s[k + 2], s[i]) + Sd5(i, k + 2) + enS[k + 2, i];
                if (!IsFinite(h) || h > 0 || es > 0) { h = Inf; es = -1.0; }
                double t = (h + DplxInitH) / (es + DplxInitS + Rc);
                if (maxTm < t && es > MinEntropyCutoff) { hMax = h; sMax = es; maxTm = t; }
            }
            return hs == 1 ? hMax : sMax;
        }
        double End5_3(int i, int hs)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
            {
                double t1 = (Hend5(k) + DplxInitH) / (Send5(k) + DplxInitS + Rc);
                double t2 = (0 + DplxInitH) / (0 + DplxInitS + Rc);
                double bh = (t1 >= t2) ? Hend5(k) : 0.0;
                double bs = (t1 >= t2) ? Send5(k) : 0.0;
                double h = bh + AtPenaltyH(s[k + 1], s[i - 1]) + Hd3(i - 1, k + 1) + enH[k + 1, i - 1];
                double es = bs + AtPenaltyS(s[k + 1], s[i - 1]) + Sd3(i - 1, k + 1) + enS[k + 1, i - 1];
                if (!IsFinite(h) || h > 0 || es > 0) { h = Inf; es = -1.0; }
                double t = (h + DplxInitH) / (es + DplxInitS + Rc);
                if (maxTm < t && es > MinEntropyCutoff) { hMax = h; sMax = es; maxTm = t; }
            }
            return hs == 1 ? hMax : sMax;
        }
        double End5_4(int i, int hs)
        {
            double hMax = Inf, sMax = -1.0, maxTm = double.NegativeInfinity;
            for (int k = 0; k <= i - MinHrpnLoop - 4; ++k)
            {
                double t1 = (Hend5(k) + DplxInitH) / (Send5(k) + DplxInitS + Rc);
                double t2 = (0 + DplxInitH) / (0 + DplxInitS + Rc);
                double bh = (t1 >= t2) ? Hend5(k) : 0.0;
                double bs = (t1 >= t2) ? Send5(k) : 0.0;
                double h = bh + AtPenaltyH(s[k + 2], s[i - 1]) + Htstack(i - 1, k + 2) + enH[k + 2, i - 1];
                double es = bs + AtPenaltyS(s[k + 2], s[i - 1]) + Ststack(i - 1, k + 2) + enS[k + 2, i - 1];
                if (!IsFinite(h) || h > 0 || es > 0) { h = Inf; es = -1.0; }
                double t = (h + DplxInitH) / (es + DplxInitS + Rc);
                if (maxTm < t && es > MinEntropyCutoff) { hMax = h; sMax = es; maxTm = t; }
            }
            return hs == 1 ? hMax : sMax;
        }

        for (int i = 2; i <= len1; ++i)
        {
            double t1 = (Hend5(i - 1) + DplxInitH) / (Send5(i - 1) + DplxInitS + Rc);
            double t2 = (End5_1(i, 1) + DplxInitH) / (End5_1(i, 2) + DplxInitS + Rc);
            double t3 = (End5_2(i, 1) + DplxInitH) / (End5_2(i, 2) + DplxInitS + Rc);
            double t4 = (End5_3(i, 1) + DplxInitH) / (End5_3(i, 2) + DplxInitS + Rc);
            double t5 = (End5_4(i, 1) + DplxInitH) / (End5_4(i, 2) + DplxInitS + Rc);
            int max = Max5(t1, t2, t3, t4, t5);
            double g, g2 = -1.0; // thal.c uses a global G2 == -1.0 here (no extra structure better than nothing)
            switch (max)
            {
                case 1: send5[i] = Send5(i - 1); hend5[i] = Hend5(i - 1); break;
                case 2:
                    g = End5_1(i, 1) - TempKelvin * End5_1(i, 2);
                    if (g < g2) { send5[i] = End5_1(i, 2); hend5[i] = End5_1(i, 1); }
                    else { send5[i] = Send5(i - 1); hend5[i] = Hend5(i - 1); }
                    break;
                case 3:
                    g = End5_2(i, 1) - TempKelvin * End5_2(i, 2);
                    if (g < g2) { send5[i] = End5_2(i, 2); hend5[i] = End5_2(i, 1); }
                    else { send5[i] = Send5(i - 1); hend5[i] = Hend5(i - 1); }
                    break;
                case 4:
                    g = End5_3(i, 1) - TempKelvin * End5_3(i, 2);
                    if (g < g2) { send5[i] = End5_3(i, 2); hend5[i] = End5_3(i, 1); }
                    else { send5[i] = Send5(i - 1); hend5[i] = Hend5(i - 1); }
                    break;
                case 5:
                    g = End5_4(i, 1) - TempKelvin * End5_4(i, 2);
                    if (g < g2) { send5[i] = End5_4(i, 2); hend5[i] = End5_4(i, 1); }
                    else { send5[i] = Send5(i - 1); hend5[i] = Hend5(i - 1); }
                    break;
            }
        }

        double mh = Hend5(len1);
        double ms = Send5(len1);
        if (!IsFinite(mh) || !IsFinite(ms)) return null; // ntthal no_structure

        // ----- tracebacku (thal.c 2820-2954) — count paired bases N -----
        var bp = new int[len1];

        bool Eq(double x, double y)
        {
            if (double.IsInfinity(x) && double.IsInfinity(y)) return (x > 0) == (y > 0);
            if (double.IsInfinity(x) || double.IsInfinity(y)) return false;
            return Math.Abs(x - y) < Equal;
        }

        // explicit stack of (i, j, mtrx)
        var stack = new System.Collections.Generic.Stack<(int i, int j, int mtrx)>();
        stack.Push((len1, 0, 1));
        while (stack.Count > 0)
        {
            var (i, j, mtrx) = stack.Pop();
            if (mtrx == 1)
            {
                while (i >= 1 && Eq(Send5(i), Send5(i - 1)) && Eq(Hend5(i), Hend5(i - 1))) --i;
                if (i == 0) continue;
                if (Eq(Send5(i), End5_1(i, 2)) && Eq(Hend5(i), End5_1(i, 1)))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 2; ++k)
                    {
                        if (Eq(Send5(i), AtPenaltyS(s[k + 1], s[i]) + enS[k + 1, i]) &&
                            Eq(Hend5(i), AtPenaltyH(s[k + 1], s[i]) + enH[k + 1, i]))
                        { stack.Push((k + 1, i, 0)); break; }
                        if (Eq(Send5(i), Send5(k) + AtPenaltyS(s[k + 1], s[i]) + enS[k + 1, i]) &&
                            Eq(Hend5(i), Hend5(k) + AtPenaltyH(s[k + 1], s[i]) + enH[k + 1, i]))
                        { stack.Push((k + 1, i, 0)); stack.Push((k, 0, 1)); break; }
                    }
                }
                else if (Eq(Send5(i), End5_2(i, 2)) && Eq(Hend5(i), End5_2(i, 1)))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
                    {
                        if (Eq(Send5(i), AtPenaltyS(s[k + 2], s[i]) + Sd5(i, k + 2) + enS[k + 2, i]) &&
                            Eq(Hend5(i), AtPenaltyH(s[k + 2], s[i]) + Hd5(i, k + 2) + enH[k + 2, i]))
                        { stack.Push((k + 2, i, 0)); break; }
                        if (Eq(Send5(i), Send5(k) + AtPenaltyS(s[k + 2], s[i]) + Sd5(i, k + 2) + enS[k + 2, i]) &&
                            Eq(Hend5(i), Hend5(k) + AtPenaltyH(s[k + 2], s[i]) + Hd5(i, k + 2) + enH[k + 2, i]))
                        { stack.Push((k + 2, i, 0)); stack.Push((k, 0, 1)); break; }
                    }
                }
                else if (Eq(Send5(i), End5_3(i, 2)) && Eq(Hend5(i), End5_3(i, 1)))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 3; ++k)
                    {
                        if (Eq(Send5(i), AtPenaltyS(s[k + 1], s[i - 1]) + Sd3(i - 1, k + 1) + enS[k + 1, i - 1]) &&
                            Eq(Hend5(i), AtPenaltyH(s[k + 1], s[i - 1]) + Hd3(i - 1, k + 1) + enH[k + 1, i - 1]))
                        { stack.Push((k + 1, i - 1, 0)); break; }
                        if (Eq(Send5(i), Send5(k) + AtPenaltyS(s[k + 1], s[i - 1]) + Sd3(i - 1, k + 1) + enS[k + 1, i - 1]) &&
                            Eq(Hend5(i), Hend5(k) + AtPenaltyH(s[k + 1], s[i - 1]) + Hd3(i - 1, k + 1) + enH[k + 1, i - 1]))
                        { stack.Push((k + 1, i - 1, 0)); stack.Push((k, 0, 1)); break; }
                    }
                }
                else if (Eq(Send5(i), End5_4(i, 2)) && Eq(Hend5(i), End5_4(i, 1)))
                {
                    for (int k = 0; k <= i - MinHrpnLoop - 4; ++k)
                    {
                        if (Eq(Send5(i), AtPenaltyS(s[k + 2], s[i - 1]) + Ststack(i - 1, k + 2) + enS[k + 2, i - 1]) &&
                            Eq(Hend5(i), AtPenaltyH(s[k + 2], s[i - 1]) + Htstack(i - 1, k + 2) + enH[k + 2, i - 1]))
                        { stack.Push((k + 2, i - 1, 0)); break; }
                        if (Eq(Send5(i), Send5(k) + AtPenaltyS(s[k + 2], s[i - 1]) + Ststack(i - 1, k + 2) + enS[k + 2, i - 1]) &&
                            Eq(Hend5(i), Hend5(k) + AtPenaltyH(s[k + 2], s[i - 1]) + Htstack(i - 1, k + 2) + enH[k + 2, i - 1]))
                        { stack.Push((k + 2, i - 1, 0)); stack.Push((k, 0, 1)); break; }
                    }
                }
            }
            else // mtrx == 0
            {
                bp[i - 1] = j;
                bp[j - 1] = i;
                var (sh10, sh11) = CalcHairpin(i, j);
                bool stackStep = Eq(enS[i, j], Ss2(i, j) + enS[i + 1, j - 1]) && Eq(enH[i, j], Hs2(i, j) + enH[i + 1, j - 1]);
                if (stackStep) { stack.Push((i + 1, j - 1, 0)); }
                else if (Eq(enS[i, j], sh10) && Eq(enH[i, j], sh11)) { /* hairpin closed: terminal */ }
                else
                {
                    bool done = false;
                    for (int d = j - i - 3; d >= MinHrpnLoop + 1 && d >= j - i - 2 - MaxLoop && !done; --d)
                    {
                        for (int ii = i + 1; ii < j - d; ++ii)
                        {
                            int jj = d + ii;
                            var r = CalcBulgeInternal2(i, j, ii, jj, traceback: true);
                            if (r is { } rr && Eq(enS[i, j], rr.S + enS[ii, jj]) && Eq(enH[i, j], rr.H + enH[ii, jj]))
                            { stack.Push((ii, jj, 0)); done = true; break; }
                        }
                    }
                }
            }
        }

        // thal.c calcHairpin counts paired positions over bp[i-1] for i = 1 .. len1-1
        // (i.e. bp[0 .. len1-2]; the last index bp[len1-1] is intentionally excluded).
        int n = 0;
        for (int i = 1; i < len1; i++) if (bp[i - 1] > 0) n++;

        // ----- calcHairpin Tm (thal.c 3229-3266): t = mh/(ms + ((N/2 - 1)*saltCorrection)) - 273.15 -----
        int half = n / 2; // integer division, as in thal.c
        double dsOut = ms + (half - 1) * saltCorrection;
        double tm = mh / dsOut - AbsoluteZero;
        double dg = mh - TempKelvin * dsOut;
        return new Result(mh, dsOut, dg, tm, half);
    }

    // thal.c max5 (returns 1..5 for the largest of T1..T5).
    private static int Max5(double t1, double t2, double t3, double t4, double t5)
    {
        int max = 1; double m = t1;
        if (t2 > m) { m = t2; max = 2; }
        if (t3 > m) { m = t3; max = 3; }
        if (t4 > m) { m = t4; max = 4; }
        if (t5 > m) { max = 5; }
        return max;
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

    private static readonly int[] TriloopKeys;
    private static readonly double[] TriloopDh;
    private static readonly double[] TriloopDs;
    private static readonly int[] TetraloopKeys;
    private static readonly double[] TetraloopDh;
    private static readonly double[] TetraloopDs;

    static NtthalHairpin()
    {
        TriloopKeys = new int[TriloopTable.Length];
        TriloopDh = new double[TriloopTable.Length];
        TriloopDs = new double[TriloopTable.Length];
        for (int k = 0; k < TriloopTable.Length; k++)
        {
            TriloopKeys[k] = Enc(TriloopTable[k].Loop);
            TriloopDh[k] = TriloopTable[k].Dh;
            TriloopDs[k] = TriloopTable[k].Ds;
        }
        TetraloopKeys = new int[TetraloopTable.Length];
        TetraloopDh = new double[TetraloopTable.Length];
        TetraloopDs = new double[TetraloopTable.Length];
        for (int k = 0; k < TetraloopTable.Length; k++)
        {
            TetraloopKeys[k] = Enc(TetraloopTable[k].Loop);
            TetraloopDh[k] = TetraloopTable[k].Dh;
            TetraloopDs[k] = TetraloopTable[k].Ds;
        }
    }
}
