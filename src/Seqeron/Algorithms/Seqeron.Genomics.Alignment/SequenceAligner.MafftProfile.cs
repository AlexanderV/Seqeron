namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: A__align (Salignmm.c), the group-group dynamic program of the progressive stage

    /// <summary>
    /// The state of one run of MAFFT's alignment binaries (<c>disttbfast</c>, later <c>tbfast</c>, <c>dvtditr</c>, <c>pairlocalalign</c>): the global variables (<see cref="MafftGlobals"/>) and the
    /// function-local <c>static</c> work arrays of the dynamic programs, which the program re-uses from one call to the next (and which it reads at one place before they are written: see
    /// <see cref="AAlign"/>, <c>ogcp1[1]</c>), so the port keeps them across calls exactly as the program does (single-thread run).
    /// </summary>
    internal sealed partial class MafftRun
    {
        internal readonly MafftGlobals G;
        internal System.Text.StringBuilder Log => G.Log;

        internal MafftRun(MafftGlobals globals)
        {
            G = globals;
        }

        // ---- the static work arrays of A__align (Salignmm.c), re-used and zeroed only when a sequence is longer than any before
        private int _aOrlgth1, _aOrlgth2;
        private double[] _aOgcp1 = Array.Empty<double>(), _aOgcp1o = Array.Empty<double>(), _aOgcp2 = Array.Empty<double>(), _aOgcp2o = Array.Empty<double>();
        private double[] _aFgcp1 = Array.Empty<double>(), _aFgcp1o = Array.Empty<double>(), _aFgcp2 = Array.Empty<double>(), _aFgcp2o = Array.Empty<double>();
        private double[] _aGapfreq1 = Array.Empty<double>(), _aGapfreq2 = Array.Empty<double>();
        private double[][] _aCpmx1 = Array.Empty<double[]>(), _aCpmx2 = Array.Empty<double[]>();

        // -------------------------------------------------------------------- tddis.c / mltaln9.c helpers

        /// <summary><c>cpmx_calc_new</c>: the frequency of every letter in every column, weighted by <paramref name="eff"/>.</summary>
        private void CpmxCalcNew(string[] seq, double[][] cpmx, double[] eff, int lgth, int clus)
        {
            for (int j = 0; j < G.Nalphabets; j++)
                Array.Clear(cpmx[j], 0, lgth);
            for (int k = 0; k < clus; k++)
            {
                double feff = eff[k];
                string s = seq[k];
                for (int j = 0; j < lgth; j++)
                    cpmx[G.AminoN[s[j]]][j] += feff;
            }
        }

        /// <summary><c>gapcountf</c>: the weighted gap frequency of every column.</summary>
        private static void GapCountF(double[] freq, string[] seq, int nseq, double[] eff, int lgth)
        {
            for (int i = 0; i < lgth; i++)
            {
                double fr = 0.0;
                for (int j = 0; j < nseq; j++)
                    if (seq[j][i] == '-')
                        fr += eff[j];
                freq[i] = fr;
            }
        }

        /// <summary><c>st_OpeningGapCount</c>: the weighted number of gap openings at every column (a gap that starts here).</summary>
        private static void StOpeningGapCount(double[] ogcp, int clus, string[] seq, double[] eff, int len)
        {
            for (int i = 0; i < len; i++)
                ogcp[i] = 0.0;
            for (int j = 0; j < clus; j++)
            {
                double feff = eff[j];
                string spt = seq[j];
                int gc = 0;
                for (int i = 0; i < len; i++)
                {
                    int gb = gc;
                    gc = spt[i] == '-' ? 1 : 0;
                    if (gb == 0 && gc != 0)
                        ogcp[i] += feff;
                }
            }
            ogcp[len] = 0.0;
        }

        /// <summary><c>st_FinalGapCount</c>: the weighted number of gap closings at every column (a gap that ends here).</summary>
        private static void StFinalGapCount(double[] fgcp, int clus, string[] seq, double[] eff, int len)
        {
            for (int i = 0; i < len; i++)
                fgcp[i] = 0.0;
            for (int j = 0; j < clus; j++)
            {
                double feff = eff[j];
                string spt = seq[j];
                int gc = spt[0] == '-' ? 1 : 0;
                for (int i = 0; i < len; i++)
                {
                    int gb = gc;
                    char following = i + 1 < spt.Length ? spt[i + 1] : '\0';
                    gc = following == '-' ? 1 : 0;
                    if (gb != 0 && gc == 0)
                        fgcp[i] += feff;
                }
            }
        }

        /// <summary><c>new_OpeningGapCount</c>: as <see cref="StOpeningGapCount"/> but the first column is compared with the last column of the segment before (<paramref name="sgappat"/>).</summary>
        private static void NewOpeningGapCount(double[] ogcp, int clus, string[] seq, double[] eff, int len, string sgappat)
        {
            for (int i = 0; i < len; i++)
                ogcp[i] = 0.0;
            for (int j = 0; j < clus; j++)
            {
                double feff = eff[j];
                string spt = seq[j];
                int gc = sgappat[j] == '-' ? 1 : 0;
                for (int i = 0; i < len; i++)
                {
                    int gb = gc;
                    gc = spt[i] == '-' ? 1 : 0;
                    if (gb == 0 && gc != 0)
                        ogcp[i] += feff;
                }
            }
        }

        /// <summary><c>new_FinalGapCount</c>: as <see cref="StFinalGapCount"/> (the program compares the last column with the terminating NUL, so the pattern of the following segment has no effect).</summary>
        private static void NewFinalGapCount(double[] fgcp, int clus, string[] seq, double[] eff, int len)
        {
            StFinalGapCount(fgcp, clus, seq, eff, len);
        }

        /// <summary><c>outgapcount</c>: the weight of the sequences that have a gap at the boundary.</summary>
        private static double OutGapCount(int nseq, string gappat, double[] eff)
        {
            double fr = 0.0;
            for (int j = 0; j < nseq; j++)
                if (gappat[j] == '-')
                    fr += eff[j];
            return fr;
        }

        /// <summary><c>match_calc</c> (FASTMATCHCALC): the score of the letters of column <paramref name="i1"/> of the first profile with every column of the second.</summary>
        private void MatchCalc(double[][] dyn, double[] match, double[][] cpmx1, double[][] cpmx2, int i1, int lgth2, double[][] doublework, int[][] intwork, bool initialize)
        {
            int na = G.Nalphabets;
            if (initialize)
            {
                for (int j = 0; j < lgth2; j++)
                {
                    int count = 0;
                    for (int l = 0; l < na; l++)
                    {
                        if (cpmx2[l][j] != 0.0)
                        {
                            doublework[j][count] = cpmx2[l][j];
                            intwork[j][count] = l;
                            count++;
                        }
                    }
                    intwork[j][count] = -1;
                }
            }
            var scarr = new double[na];
            for (int l = 0; l < na; l++)
            {
                scarr[l] = 0.0;
                for (int j = 0; j < na; j++)
                    scarr[l] += dyn[j][l] * cpmx1[j][i1];
            }
            for (int j = 0; j < lgth2; j++)
            {
                double m = 0.0;
                int[] idx = intwork[j];
                double[] val = doublework[j];
                for (int k = 0; idx[k] > -1; k++)
                    m += scarr[idx[k]] * val[k];
                match[j] = m;
            }
        }

        // -------------------------------------------------------------------- A__align

        /// <summary>
        /// <c>A__align</c> of <c>Salignmm.c</c>: the dynamic program that aligns two groups of aligned sequences with the gap penalties scaled by the gap frequencies of the columns
        /// (<c>legacygapcost == 0</c>) of the groups. <paramref name="seq1"/> / <paramref name="seq2"/> hold the rows of the groups and receive the aligned rows (in place, as the program does).
        /// <paramref name="hist"/> (the <c>cpmxhist</c> of <c>treebase</c>) with <paramref name="child0"/>, <paramref name="child1"/> and <paramref name="resultIndex"/> is the memory of the
        /// profiles of the children that <c>treebase</c> hands over (only for groups of more than 20 sequences); <paramref name="sgap1"/> ... <paramref name="egap2"/> are the boundary gaps of the
        /// segment that <c>Falign</c> hands over (<c>null</c> for <c>treebase</c>).
        /// </summary>
        internal double AAlign(double[][] dyn, int penaltyL, int penaltyExL, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, int alloclen,
            string? sgap1, string? sgap2, string? egap1, string? egap2, int headgp, int tailgp, double[][][]? hist, int child0, int child1, int resultIndex, double orieff1, double orieff2, double[][]? imp = null)
        {
            int na = G.Nalphabets;
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            if (lgth1 == 0 && lgth2 == 0)
                return 0.0;
            if (lgth1 == 0)
            {
                for (int i = 0; i < icyc; i++)
                    seq1[i] = new string(G.NewGapChar, lgth2);
                return 0.0;
            }
            if (lgth2 == 0)
            {
                for (int i = 0; i < jcyc; i++)
                    seq2[i] = new string(G.NewGapChar, lgth1);
                return 0.0;
            }
            if (G.TryWarp)
                throw new NotSupportedException("--allowshift / --shiftpenalty (the warp of A__align) is not available in this port.");

            if (lgth1 > _aOrlgth1 || lgth2 > _aOrlgth2)
            {
                int ll1 = Math.Max((int)(1.3 * lgth1), _aOrlgth1) + 100;
                int ll2 = Math.Max((int)(1.3 * lgth2), _aOrlgth2) + 100;
                _aOgcp1 = new double[ll1 + 2]; _aOgcp1o = new double[ll1 + 2];
                _aOgcp2 = new double[ll2 + 2]; _aOgcp2o = new double[ll2 + 2];
                _aFgcp1 = new double[ll1 + 2]; _aFgcp1o = new double[ll1 + 2];
                _aFgcp2 = new double[ll2 + 2]; _aFgcp2o = new double[ll2 + 2];
                _aCpmx1 = new double[na][];
                _aCpmx2 = new double[na][];
                for (int i = 0; i < na; i++)
                {
                    _aCpmx1[i] = new double[ll1 + 2];
                    _aCpmx2[i] = new double[ll2 + 2];
                }
                _aGapfreq1 = new double[ll1 + 2];
                _aGapfreq2 = new double[ll2 + 2];
                _aOrlgth1 = ll1 - 100;
                _aOrlgth2 = ll2 - 100;
            }

            double fpenalty = penaltyL;
            double fpenaltyEx = penaltyExL;
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            var lastverticalw = new double[lgth1 + 2];
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var ijp = new int[lgth1 + 1][];
            for (int i = 0; i <= lgth1; i++)
                ijp[i] = new int[lgth2 + 1];
            var doublework = new double[Math.Max(lgth1, lgth2) + 2][];
            var intwork = new int[Math.Max(lgth1, lgth2) + 2][];
            for (int i = 0; i < doublework.Length; i++)
            {
                doublework[i] = new double[na];
                intwork[i] = new int[na + 1];
            }

            double[][] cpmx1pt, cpmx2pt;
            double[] gapfreq1pt, gapfreq2pt, ogcp1opt, ogcp2opt, fgcp1opt, fgcp2opt;
            double headgapfreq1, headgapfreq2;
            bool history1 = false, history2 = false;
            double[][]? cpmx1hist = null, cpmx2hist = null;

            if (hist != null && G.SpecificityConsideration == 0.0)
            {
                if (sgap1 != null)
                    throw new InvalidOperationException("The combination of sgap1 and cpmxhit is not supported. See Salignmm.c");
                if (child0 >= 0 && hist[child0] != null)
                {
                    cpmx1hist = hist[child0];
                    history1 = true;
                    cpmx1pt = cpmx1hist!;
                    gapfreq1pt = cpmx1hist![na];
                    ogcp1opt = cpmx1hist[na + 1];
                    fgcp1opt = cpmx1hist[na + 2];
                }
                else
                {
                    cpmx1pt = _aCpmx1;
                    CpmxCalcNew(seq1, cpmx1pt, eff1, lgth1, icyc);
                    gapfreq1pt = _aGapfreq1;
                    GapCountF(gapfreq1pt, seq1, icyc, eff1, lgth1);
                    for (int i = 0; i < lgth1; i++)
                        gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
                    ogcp1opt = _aOgcp1o;
                    fgcp1opt = _aFgcp1o;
                    StOpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1);
                    StFinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1);
                }
                if (child1 >= 0 && hist[child1] != null)
                {
                    cpmx2hist = hist[child1];
                    history2 = true;
                    cpmx2pt = cpmx2hist!;
                    gapfreq2pt = cpmx2hist![na];
                    ogcp2opt = cpmx2hist[na + 1];
                    fgcp2opt = cpmx2hist[na + 2];
                }
                else
                {
                    cpmx2pt = _aCpmx2;
                    CpmxCalcNew(seq2, cpmx2pt, eff2, lgth2, jcyc);
                    gapfreq2pt = _aGapfreq2;
                    GapCountF(gapfreq2pt, seq2, jcyc, eff2, lgth2);
                    for (int i = 0; i < lgth2; i++)
                        gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
                    ogcp2opt = _aOgcp2o;
                    fgcp2opt = _aFgcp2o;
                    StOpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2);
                    StFinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2);
                }
                headgapfreq1 = 1.0;
                headgapfreq2 = 1.0;
                gapfreq1pt[lgth1] = 1.0;
                gapfreq2pt[lgth2] = 1.0;
            }
            else
            {
                cpmx1pt = _aCpmx1;
                cpmx2pt = _aCpmx2;
                gapfreq1pt = _aGapfreq1;
                gapfreq2pt = _aGapfreq2;
                ogcp1opt = _aOgcp1o;
                ogcp2opt = _aOgcp2o;
                fgcp1opt = _aFgcp1o;
                fgcp2opt = _aFgcp2o;
                CpmxCalcNew(seq1, cpmx1pt, eff1, lgth1, icyc);
                CpmxCalcNew(seq2, cpmx2pt, eff2, lgth2, jcyc);
                if (sgap1 != null)
                {
                    NewOpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1, sgap1);
                    NewFinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1);
                    NewOpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2, sgap2!);
                    NewFinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2);
                    headgapfreq1 = OutGapCount(icyc, sgap1, eff1);
                    headgapfreq2 = OutGapCount(jcyc, sgap2!, eff2);
                    gapfreq1pt[lgth1] = OutGapCount(icyc, egap1!, eff1);
                    gapfreq2pt[lgth2] = OutGapCount(jcyc, egap2!, eff2);
                }
                else
                {
                    StOpeningGapCount(ogcp1opt, icyc, seq1, eff1, lgth1);
                    StFinalGapCount(fgcp1opt, icyc, seq1, eff1, lgth1);
                    StOpeningGapCount(ogcp2opt, jcyc, seq2, eff2, lgth2);
                    StFinalGapCount(fgcp2opt, jcyc, seq2, eff2, lgth2);
                    headgapfreq1 = 0.0;
                    headgapfreq2 = 0.0;
                    gapfreq1pt[lgth1] = 0.0;
                    gapfreq2pt[lgth2] = 0.0;
                }
                if (G.LegacyGapCost == 0)
                {
                    GapCountF(gapfreq1pt, seq1, icyc, eff1, lgth1);
                    GapCountF(gapfreq2pt, seq2, jcyc, eff2, lgth2);
                    for (int i = 0; i < lgth1 + 1; i++)
                        gapfreq1pt[i] = 1.0 - gapfreq1pt[i];
                    for (int i = 0; i < lgth2 + 1; i++)
                        gapfreq2pt[i] = 1.0 - gapfreq2pt[i];
                    headgapfreq1 = 1.0 - headgapfreq1;
                    headgapfreq2 = 1.0 - headgapfreq2;
                }
                else
                {
                    for (int i = 0; i < lgth1 + 1; i++)
                        gapfreq1pt[i] = 1.0;
                    for (int i = 0; i < lgth2 + 1; i++)
                        gapfreq2pt[i] = 1.0;
                    headgapfreq1 = 1.0;
                    headgapfreq2 = 1.0;
                }
            }

            double[] ogcp1 = _aOgcp1, ogcp2 = _aOgcp2, fgcp1 = _aFgcp1, fgcp2 = _aFgcp2;
            for (int i = 0; i < lgth1; i++)
            {
                ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * gapfreq1pt[i];
                fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * gapfreq1pt[i];
            }
            for (int i = 0; i < lgth2; i++)
            {
                ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * gapfreq2pt[i];
                fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * gapfreq2pt[i];
            }

            double[] currentw = w1, previousw = w2;
            MatchCalc(dyn, initverticalw, cpmx2pt, cpmx1pt, 0, lgth1, doublework, intwork, true);
            if (imp != null)
                for (int i = 0; i < lgth1; i++)
                    initverticalw[i] += imp[i][0];
            MatchCalc(dyn, currentw, cpmx1pt, cpmx2pt, 0, lgth2, doublework, intwork, true);
            if (imp != null)
                for (int j = 0; j < lgth2; j++)
                    currentw[j] += imp[0][j];
            if (headgp == 1)
            {
                for (int i = 1; i < lgth1 + 1; i++)
                {
                    initverticalw[i] += ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2pt[0];
                    initverticalw[i] += fpenaltyEx * i;
                }
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    currentw[j] += ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1pt[0];
                    currentw[j] += fpenaltyEx * j;
                }
            }
            else
            {
                for (int i = 1; i < lgth1 + 1; i++)
                {
                    initverticalw[i] += fpenalty * 0.0;
                    initverticalw[i] += fpenaltyEx * i * 0.0;
                }
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    currentw[j] += fpenalty * 0.0;
                    currentw[j] += fpenaltyEx * j * 0.0;
                }
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2pt[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int lasti = tailgp != 0 ? lgth1 + 1 : lgth1;
            int lastj = lgth2 + 1;
            double wm = 0.0;
            double g;
            double mi;
            int mpi;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                MatchCalc(dyn, currentw, cpmx1pt, cpmx2pt, i, lgth2, doublework, intwork, false);
                if (imp != null)
                {
                    var impRow = imp[i];
                    for (int j = 0; j < lgth2; j++)
                        currentw[j] += impRow[j];
                }
                currentw[0] = initverticalw[i];
                mi = previousw[0] + ogcp2[1] * gapfreq1pt[i - 1];
                mpi = 0;
                var ijpRow = ijp[i];
                double fgcp1va = fgcp1[i - 1];
                double ogcp1va = ogcp1[i];
                double gf1va = gapfreq1pt[i];
                double gf1vapre = gapfreq1pt[i - 1];
                for (int j = 1; j < lastj; j++)
                {
                    wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    if ((g = mi + fgcp2[j - 1] * gf1va) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = previousw[j - 1] + ogcp2[j] * gf1vapre) >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    mi += fpenaltyEx;
                    if ((g = m[j] + fgcp1va * gapfreq2pt[j]) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    if ((g = previousw[j - 1] + ogcp1va * gapfreq2pt[j - 1]) >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    m[j] += fpenaltyEx;
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
            }

            var (gt1, gt2, mseq1, mseq2, wmo) = Atracking(G.NewGapChar, currentw, lastverticalw, fpenalty, fpenaltyEx, seq1, seq2, ijp, icyc, jcyc, tailgp, imp != null);
            if (tailgp == 0 && imp == null)
                wm = wmo;

            if (hist != null && G.SpecificityConsideration == 0.0)
            {
                if (icyc + jcyc > 20)
                {
                    int alen = gt1.Length;
                    int limk = (alen + 1000) / 1000 * 1000;
                    _ = limk;
                    double totaleff1 = 0.0, totaleff2 = 0.0;
                    for (int i = 0; i < icyc; i++)
                        totaleff1 += eff1[i];
                    for (int j = 0; j < jcyc; j++)
                        totaleff2 += eff2[j];
                    if (Math.Abs(totaleff1 - 1.0) > 0.001 || Math.Abs(totaleff2 - 1.0) > 0.001)
                        throw new InvalidOperationException("MAFFT stops: the weights of a group do not add up to 1.");
                    totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
                    totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);
                    var res = new double[na + 3][];
                    CreateCpmxResult(res, alen, totaleff1, totaleff2, cpmx1pt, cpmx2pt, gt1, gt2);
                    res[na] = CreateGapFreqResult(alen, totaleff1, totaleff2, gapfreq1pt, gapfreq2pt, gt1, gt2);
                    res[na + 1] = CreateOgResult(alen, totaleff1, totaleff2, ogcp1opt, ogcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2);
                    res[na + 2] = CreateFgResult(alen, totaleff1, totaleff2, fgcp1opt, fgcp2opt, gapfreq1pt, gapfreq2pt, gt1, gt2);
                    hist[resultIndex] = res;
                    if (history1)
                        hist[child0] = null!;
                    if (history2)
                        hist[child1] = null!;
                }
                else
                {
                    hist[resultIndex] = null!;
                }
            }

            int resultlen = mseq1[0].Length;
            if (alloclen < resultlen)
                throw new InvalidOperationException("LENGTH OVER!");
            for (int i = 0; i < icyc; i++)
                seq1[i] = mseq1[i];
            for (int j = 0; j < jcyc; j++)
                seq2[j] = mseq2[j];
            return wm;
        }

        /// <summary><c>Atracking</c> of <c>Salignmm.c</c> (no shift): the trace back, the gap tables of the two groups ('o' = a column of the group, '-' = a new gap) and the aligned rows.</summary>
        private static (string Gt1, string Gt2, string[] Mseq1, string[] Mseq2, double Wm) Atracking(char gap, double[] lasthorizontalw, double[] lastverticalw, double fpenalty, double fpenaltyEx,
            string[] seq1, string[] seq2, int[][] ijp, int icyc, int jcyc, int tailgp, bool localhom = false)
        {
            int lgth1 = seq1[0].Length, lgth2 = seq2[0].Length;
            double wm = 0.0;
            if (tailgp != 1 && localhom)
            {
                // Atracking_localhom: the end of the alignment is the cell of the last row / column with the largest value (no end-gap term)
                wm = lastverticalw[0];
                for (int i = 0; i < lgth1; i++)
                {
                    if (lastverticalw[i] >= wm)
                    {
                        wm = lastverticalw[i];
                        ijp[lgth1][lgth2] = +(lgth1 - i);
                    }
                }
                for (int j = 0; j < lgth2; j++)
                {
                    if (lasthorizontalw[j] >= wm)
                    {
                        wm = lasthorizontalw[j];
                        ijp[lgth1][lgth2] = -(lgth2 - j);
                    }
                }
            }
            else if (tailgp != 1)
            {
                double g;
                wm = lasthorizontalw[lgth2 - 1] - 1.0;
                for (int j = lgth2 - 2; j >= 0; j--)
                {
                    if ((g = lasthorizontalw[j] + (fpenalty * 0.0 + fpenaltyEx * (lgth2 - 1 - j) * 0.0)) > wm)
                    {
                        wm = g;
                        ijp[lgth1][lgth2] = -(lgth2 - j);
                    }
                }
                for (int i = lgth1 - 2; i >= 0; i--)
                {
                    if ((g = lastverticalw[i] + (fpenalty * 0.0 + fpenaltyEx * (lgth1 - 1 - i) * 0.0)) > wm)
                    {
                        wm = g;
                        ijp[lgth1][lgth2] = +(lgth1 - i);
                    }
                }
                if (lasthorizontalw[lgth2 - 1] > wm)
                {
                    wm = lasthorizontalw[lgth2 - 1];
                    ijp[lgth1][lgth2] = 0;
                }
            }
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = i + 1;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = -(j + 1);

            var gt1 = new char[lgth1 + lgth2 + 1];
            var gt2 = new char[lgth1 + lgth2 + 1];
            int p1 = lgth1 + lgth2, p2 = lgth1 + lgth2;
            int iin = lgth1, jin = lgth2;
            int limk = lgth1 + lgth2 + 1;
            int k = 0;
            while (k < limk)
            {
                int ifi, jfi;
                if (ijp[iin][jin] < 0)
                {
                    ifi = iin - 1;
                    jfi = jin + ijp[iin][jin];
                }
                else if (ijp[iin][jin] > 0)
                {
                    ifi = iin - ijp[iin][jin];
                    jfi = jin - 1;
                }
                else
                {
                    ifi = iin - 1;
                    jfi = jin - 1;
                }
                int l = iin - ifi;
                while (--l > 0)
                {
                    gt1[--p1] = 'o';
                    gt2[--p2] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    gt1[--p1] = '-';
                    gt2[--p2] = 'o';
                    k++;
                }
                if (iin <= 0 || jin <= 0)
                    break;
                gt1[--p1] = 'o';
                gt2[--p2] = 'o';
                k++;
                iin = ifi;
                jin = jfi;
                k++;
            }
            string gaptable1 = new string(gt1, p1, lgth1 + lgth2 - p1);
            string gaptable2 = new string(gt2, p2, lgth1 + lgth2 - p2);
            var mseq1 = new string[icyc];
            var mseq2 = new string[jcyc];
            bool ngap1 = gaptable1.Contains('-');
            bool ngap2 = gaptable2.Contains('-');
            for (int i = 0; i < icyc; i++)
                mseq1[i] = ngap1 ? GapIreru(seq1[i], gaptable1, gap) : seq1[i][..Math.Min(gaptable1.Length, seq1[i].Length)];
            for (int j = 0; j < jcyc; j++)
                mseq2[j] = ngap2 ? GapIreru(seq2[j], gaptable2, gap) : seq2[j][..Math.Min(gaptable2.Length, seq2[j].Length)];
            return (gaptable1, gaptable2, mseq1, mseq2, wm);
        }

        /// <summary><c>gapireru</c>: the residues of <paramref name="ori"/> in the columns marked 'o' of the gap table, a gap in the others.</summary>
        internal static string GapIreru(string ori, string gt, char gap = '-')
        {
            var sb = new System.Text.StringBuilder(gt.Length);
            int k = 0;
            foreach (char g in gt)
                sb.Append(g == '-' ? gap : ori[k++]);
            return sb.ToString();
        }

        private void CreateCpmxResult(double[][] res, int alen, double eff1, double eff2, double[][] cpmx1, double[][] cpmx2, string gaptable1, string gaptable2)
        {
            for (int i = 0; i < G.Nalphabets; i++)
            {
                var r = new double[alen + 1];
                int p = 0;
                for (int j = 0; j < alen; j++)
                    if (gaptable1[j] != '-')
                        r[j] += cpmx1[i][p++] * eff1;
                p = 0;
                for (int j = 0; j < alen; j++)
                    if (gaptable2[j] != '-')
                        r[j] += cpmx2[i][p++] * eff2;
                res[i] = r;
            }
        }

        private static double[] CreateGapFreqResult(int alen, double eff1, double eff2, double[] gapf1, double[] gapf2, string gaptable1, string gaptable2)
        {
            var r = new double[alen + 1];
            int p = 0;
            for (int j = 0; j < alen + 1; j++)
                if (j >= alen || gaptable1[j] != '-')
                    r[j] += gapf1[p++] * eff1;
            p = 0;
            for (int j = 0; j < alen; j++)
                if (gaptable2[j] != '-')
                    r[j] += gapf2[p++] * eff2;
            r[alen] = 1.0;
            return r;
        }

        private static double[] CreateOgResult(int alen, double eff1, double eff2, double[] ori1, double[] ori2, double[] gf1, double[] gf2, string gaptable1, string gaptable2)
        {
            var r = new double[alen + 1];
            int p = 0;
            for (int j = 0; j < alen; j++)
            {
                if (gaptable1[j] == '-')
                {
                    if (j == 0)
                        r[j] += 1.0 * eff1;
                    else if (gaptable1[j - 1] != '-')
                        r[j] += gf1[p - 1] * eff1;
                }
                else
                {
                    if (j == 0 || gaptable1[j - 1] != '-')
                        r[j] += ori1[p] * eff1;
                    p++;
                }
            }
            p = 0;
            for (int j = 0; j < alen; j++)
            {
                if (gaptable2[j] == '-')
                {
                    if (j == 0)
                        r[j] += 1.0 * eff2;
                    else if (gaptable2[j - 1] != '-')
                        r[j] += gf2[p - 1] * eff2;
                }
                else
                {
                    if (j == 0 || gaptable2[j - 1] != '-')
                        r[j] += ori2[p] * eff2;
                    p++;
                }
            }
            return r;
        }

        private static double[] CreateFgResult(int alen, double eff1, double eff2, double[] ori1, double[] ori2, double[] gf1, double[] gf2, string gaptable1, string gaptable2)
        {
            var r = new double[alen + 1];
            int p = 0;
            for (int j = 0; j < alen; j++)
            {
                if (gaptable1[j] == '-')
                {
                    if (j == alen - 1)
                        r[j] += eff1;
                    else if (gaptable1[j + 1] != '-')
                        r[j] += gf1[p] * eff1;
                }
                else
                {
                    if (j + 1 >= alen || gaptable1[j + 1] != '-')
                        r[j] += ori1[p] * eff1;
                    p++;
                }
            }
            p = 0;
            for (int j = 0; j < alen; j++)
            {
                if (gaptable2[j] == '-')
                {
                    if (j == alen - 1)
                        r[j] += eff2;
                    else if (gaptable2[j + 1] != '-')
                        r[j] += gf2[p] * eff2;
                }
                else
                {
                    if (j + 1 >= alen || gaptable2[j + 1] != '-')
                        r[j] += ori2[p] * eff2;
                    p++;
                }
            }
            return r;
        }
    }

    #endregion
}
