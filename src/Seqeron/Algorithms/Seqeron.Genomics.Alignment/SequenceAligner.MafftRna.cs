#pragma warning disable S3776, S107, S138, S1541, S1854, S3358, S1066, S125, CA1814, S2368, S1450, IDE0059, S2234
using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    /// <summary><c>RNApair</c> of <c>mltaln.h</c>: a partner position and its base-pairing probability.</summary>
    internal sealed class RnaPair
    {
        internal int Pos;
        internal double Score;

        internal RnaPair(int pos, double score)
        {
            Pos = pos;
            Score = score;
        }
    }

    internal sealed partial class MafftRun
    {
        // -------------------------------------------------------------------- the RNA branch of tbfast (rnakozo, rnaprediction 'm'): hat4, foldrna, Lalignmm_hmout

        /// <summary>The text of the file <c>hat4</c> that <c>tbfast</c> reads when <see cref="MafftGlobals.RnaKozo"/> is set (written by <c>mccaskillwrap</c> / <c>contrafoldwrap</c>); <c>null</c>: "Make hat4 using mccaskill.".</summary>
        internal string? Hat4Text { get; set; }

        /// <summary><c>singlerna</c> of <c>tbfast.c</c>: for every sequence and every residue the list of its partners.</summary>
        private List<RnaPair>[][]? _singleRna;

        /// <summary><c>tbfast.c</c>: "Loading 'hat4' ... " (<c>readmccaskill</c> for every sequence).</summary>
        private void LoadHat4(string[] bseq)
        {
            if (Hat4Text == null)
                throw new MafftProgramFailureException(Log + "Make hat4 using mccaskill.\n");
            Report("Loading 'hat4' ... ");
            var lines = Hat4Text.Split('\n');
            int pos = 0;
            int njob = bseq.Length;
            _singleRna = new List<RnaPair>[njob][];
            for (int i = 0; i < njob; i++)
            {
                int nogaplen = bseq[i].Length;
                var lists = new List<RnaPair>[nogaplen];
                for (int j = 0; j < nogaplen; j++)
                    lists[j] = new List<RnaPair>();
                _singleRna[i] = lists;
                // readmccaskill
                if (pos >= lines.Length || lines[pos].Length == 0 || lines[pos][0] != '>')
                    throw new MafftProgramFailureException(Log + "format error in hat4 - 1\n");
                pos++;
                while (pos < lines.Length && lines[pos].Length > 0 && lines[pos][0] != '>')
                {
                    var f = lines[pos++].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (f.Length < 3)
                        continue;
                    int left = MafftGlobals.CAtoi(f[0]);
                    int right = MafftGlobals.CAtoi(f[1]);
                    double prob = MafftGlobals.CAtof(f[2]);
                    if (left >= nogaplen || right >= nogaplen)
                        throw new MafftProgramFailureException(Log + "format error in hat4 - 2\n");
                    if (prob < 0.01)
                        continue;
                    if (left != right && prob > 0.0)
                    {
                        lists[left].Add(new RnaPair(right, prob));
                        lists[right].Add(new RnaPair(left, prob));
                    }
                }
            }
            Report("\ndone.\n");
        }

        /// <summary><c>imp_rna</c> / <c>imp_rnaD</c>: the structure term of one step of the progressive alignment is added to the importance matrix of the pair library (<c>foldrna</c>, <c>rna.c</c>).</summary>
        private void ImpRna(int clus1, int clus2, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int[] members1, int[] members2, double[][] impmtx)
        {
            var group1 = new List<RnaPair>[clus1][];
            var group2 = new List<RnaPair>[clus2][];
            for (int k = 0; k < clus1; k++)
                group1[k] = _singleRna![members1[k]];
            for (int k = 0; k < clus2; k++)
                group2[k] = _singleRna![members2[k]];
            FoldRna(clus1, clus2, seq1, seq2, eff1, eff2, group1, group2, impmtx);
        }

        private static string Utot(string s)
        {
            var c = new char[s.Length];
            for (int i = 0; i < c.Length; i++)
                c[i] = s[i] switch { 'a' => 'a', 't' => 't', 'u' => 't', 'g' => 'g', 'c' => 'c', '-' => '-', _ => 'n' };
            return new string(c);
        }

        /// <summary><c>commongappick_record</c> for one sequence: the sequence without its gaps and the column of every residue.</summary>
        private static (string Residues, int[] Map) PickGaps(string s)
        {
            var sb = new System.Text.StringBuilder(s.Length);
            var map = new List<int>(s.Length);
            for (int i = 0; i < s.Length; i++)
            {
                if (s[i] == '-')
                    continue;
                sb.Append(s[i]);
                map.Add(i);
            }
            return (sb.ToString(), map.ToArray());
        }

        /// <summary><c>mccaskillextract</c>: the weighted base-pairing probabilities of the sequences of a group, by column.</summary>
        private static List<RnaPair>[] McCaskillExtract(int lgth, string[] nogap, int nseq, List<RnaPair>[][] single, int[][] sgapmap, double[] eff)
        {
            var pairprob = new List<RnaPair>[lgth];
            for (int i = 0; i < lgth; i++)
                pairprob[i] = new List<RnaPair>();
            for (int i = 0; i < nseq; i++)
            {
                int nogaplgth = nogap[i].Length;
                for (int j = 0; j < nogaplgth; j++)
                {
                    foreach (var pt in single[i][j])
                    {
                        int left = sgapmap[i][j];
                        int right = sgapmap[i][pt.Pos];
                        double prob = pt.Score;
                        RnaPair? pt2 = null;
                        foreach (var cand in pairprob[left])
                        {
                            if (cand.Pos == right)
                            {
                                pt2 = cand;
                                break;
                            }
                        }
                        if (pt2 == null)
                        {
                            pt2 = new RnaPair(right, 0.0);
                            pairprob[left].Add(pt2);
                        }
                        pt2.Score += prob * eff[i];
                    }
                }
            }
            return pairprob;
        }

        /// <summary><c>foldrna</c> of <c>rna.c</c> (<c>rnaprediction</c> 'm', <c>RNAscoremtx</c> 'n'): the importance of the pair of columns (i, j) grows with the probability that the partners of i and j are aligned too (four-way consistency).</summary>
        private void FoldRna(int nseq1, int nseq2, string[] seq1, string[] seq2, double[] eff1, double[] eff2, List<RnaPair>[][] grouprna1, List<RnaPair>[][] grouprna2, double[][] impmtx)
        {
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            var useq1 = new string[nseq1];
            var useq2 = new string[nseq2];
            var sgapmap1 = new int[nseq1][];
            var sgapmap2 = new int[nseq2][];
            var oseq1 = new string[nseq1];
            var oseq2 = new string[nseq2];
            for (int i = 0; i < nseq1; i++)
            {
                (useq1[i], sgapmap1[i]) = PickGaps(seq1[i]);
                oseq1[i] = Utot(seq1[i]);
            }
            for (int i = 0; i < nseq2; i++)
            {
                (useq2[i], sgapmap2[i]) = PickGaps(seq2[i]);
                oseq2[i] = Utot(seq2[i]);
            }
            var pairprob1 = McCaskillExtract(lgth1, useq1, nseq1, grouprna1, sgapmap1, eff1);
            var pairprob2 = McCaskillExtract(lgth2, useq2, nseq2, grouprna2, sgapmap2, eff2);
            var map = new double[lgth1][];
            for (int i = 0; i < lgth1; i++)
                map[i] = new double[lgth2];
            LalignmmHmout(oseq1, oseq2, eff1, eff2, nseq1, nseq2, map);
            var impmtx2 = new double[lgth1][];
            for (int i = 0; i < lgth1; i++)
                impmtx2[i] = new double[lgth2];
            for (int i = 0; i < lgth1; i++)
            {
                foreach (var pairpt1 in pairprob1[i])
                {
                    for (int j = 0; j < lgth2; j++)
                    {
                        foreach (var pairpt2 in pairprob2[j])
                        {
                            int uido = pairpt1.Pos;
                            int ujdo = pairpt2.Pos;
                            double prob = pairpt1.Score * pairpt2.Score;
                            if (uido > -1 && ujdo > -1 && (i > uido && j > ujdo || i < uido && j < ujdo))
                                impmtx2[i][j] += Math.Max(0.0, map[uido][ujdo]) * G.ConsweightRna * 600 * prob;
                        }
                    }
                }
            }
            for (int i = 0; i < lgth1; i++)
                for (int j = 0; j < lgth2; j++)
                    impmtx[i][j] += impmtx2[i][j];
        }

        /// <summary><c>match_calc</c> of <c>Lalignmm.c</c>: the score of column <paramref name="i1"/> of profile 1 against every column of profile 2 (matrix scores minus <c>RNAthr</c>); <paramref name="initialize"/> makes the sparse lists of profile 2.</summary>
        private void RnaMatchCalc(double[] match, double[][] cpmx1, double[][] cpmx2, int i1, int lgth2, double[][] doublework, int[][] intwork, bool initialize)
        {
            int na = G.Nalphabets;
            if (initialize)
            {
                for (int j = 0; j < lgth2; j++)
                {
                    int count = 0;
                    for (int l = 0; l < na; l++)
                    {
                        if (cpmx2[j][l] != 0.0)
                        {
                            doublework[j][count] = cpmx2[j][l];
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
                for (int k = 0; k < na; k++)
                    scarr[l] += (G.NDis[k][l] - G.RnaThr) * cpmx1[i1][k];
            }
            for (int j = 0; j < lgth2; j++)
            {
                match[j] = 0.0;
                for (int k = 0, cpkd; (cpkd = intwork[j][k]) > -1; k++)
                    match[j] += scarr[cpkd] * doublework[j][k];
            }
        }

        private void MsCpmxCalcNew(string[] seq, double[][] cpmx, double[] eff, int lgth, int clus)
        {
            for (int j = 0; j < lgth; j++)
                Array.Clear(cpmx[j], 0, cpmx[j].Length);
            for (int k = 0; k < clus; k++)
            {
                double feff = eff[k];
                string s = seq[k];
                for (int j = 0; j < lgth; j++)
                    cpmx[j][G.AminoN[s[j]]] += feff;
            }
        }

        /// <summary>
        /// <c>Lalignmm_hmout</c> (<c>Lalignmm.c</c>): the forward and backward pass of the profile alignment (the Myers-Miller scan of <c>MSalignmm_rec</c> with the penalties <c>RNApenalty</c> and the scores of <c>match_calc</c>) and, for every pair of columns,
        /// the score of the best alignment through it divided by the score of the best alignment (<paramref name="map"/>). Profiles shorter than ten columns (<c>DPTANNI</c>) give no map (the matrix stays zero).
        /// </summary>
        private void LalignmmHmout(string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, double[][] map)
        {
            double fpenalty = G.RnaPenalty;
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            int ll1 = lgth1 + 100;
            int ll2 = lgth2 + 100;
            int na = G.Nalphabets;
            var cpmx1 = new double[ll1 + 2][];
            var cpmx2 = new double[ll2 + 2][];
            for (int i = 0; i < cpmx1.Length; i++)
                cpmx1[i] = new double[na + 1];
            for (int i = 0; i < cpmx2.Length; i++)
                cpmx2[i] = new double[na + 1];
            var ogcp1 = new double[ll1 + 2];
            var ogcp2 = new double[ll2 + 2];
            var fgcp1 = new double[ll1 + 2];
            var fgcp2 = new double[ll2 + 2];
            MsCpmxCalcNew(seq1, cpmx1, eff1, lgth1, icyc);
            MsCpmxCalcNew(seq2, cpmx2, eff2, lgth2, jcyc);
            StOpeningGapCount(ogcp1, icyc, seq1, eff1, lgth1);
            StOpeningGapCount(ogcp2, jcyc, seq2, eff2, lgth2);
            StFinalGapCount(fgcp1, icyc, seq1, eff1, lgth1);
            StFinalGapCount(fgcp2, jcyc, seq2, eff2, lgth2);
            for (int i = 0; i < lgth1; i++)
            {
                ogcp1[i] = 0.5 * (1.0 - ogcp1[i]) * fpenalty;
                fgcp1[i] = 0.5 * (1.0 - fgcp1[i]) * fpenalty;
            }
            for (int i = 0; i < lgth2; i++)
            {
                ogcp2[i] = 0.5 * (1.0 - ogcp2[i]) * fpenalty;
                fgcp2[i] = 0.5 * (1.0 - fgcp2[i]) * fpenalty;
            }
            // MSalignmm_rec( ..., 0, lgth1-1, 0, lgth2-1, ... )
            if (lgth2 <= 0 || lgth1 < 10 || lgth2 < 10)
                return;
            var w1 = new double[ll2 + 2];
            var w2 = new double[ll2 + 2];
            var midw = new double[ll2 + 2];
            var midn = new double[ll2 + 2];
            var midm = new double[ll2 + 2];
            var jumpbacki = new int[ll2 + 2];
            var jumpbackj = new int[ll2 + 2];
            var jumpforwi = new int[ll2 + 2];
            var jumpforwj = new int[ll2 + 2];
            var initverticalw = new double[ll1 + 2];
            var lastverticalw = new double[ll1 + 2];
            var m = new double[ll2 + 2];
            var mp = new int[ll2 + 2];
            int wl = Math.Max(ll1, ll2) + 2;
            var doublework = new double[wl][];
            var intwork = new int[wl][];
            for (int i = 0; i < wl; i++)
            {
                doublework[i] = new double[na + 1];
                intwork[i] = new int[na + 1];
            }
            var wmmtx = new double[ll1][];
            var wmmtx2 = new double[ll1][];
            for (int i = 0; i < ll1; i++)
            {
                wmmtx[i] = new double[ll2];
                wmmtx2[i] = new double[ll2];
            }
            double[] currentw = w1, previousw = w2, wtmp;
            RnaMatchCalc(initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, true);
            RnaMatchCalc(currentw, cpmx1, cpmx2, 0, lgth2, doublework, intwork, true);
            for (int i = 1; i < lgth1 + 1; i++)
                initverticalw[i] += ogcp1[0] + fgcp1[i - 1];
            for (int j = 1; j < lgth2 + 1; j++)
                currentw[j] += ogcp2[0] + fgcp2[j - 1];
            wmmtx[0][0] = initverticalw[0];
            for (int i = 1; i < lgth1 + 1; i++)
                wmmtx[i][0] = initverticalw[i];
            for (int j = 1; j < lgth2 + 1; j++)
                wmmtx[0][j] = currentw[j];
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1] + ogcp1[1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int imid = (int)(lgth1 * 0.5);
            int jumpi = 0, jumpj = 0, jmid = 0;
            int lasti = lgth1 + 1;
            double maxwm = 0.0;
            for (int i = 1; i < lasti; i++)
            {
                wtmp = previousw;
                previousw = currentw;
                currentw = wtmp;
                previousw[0] = initverticalw[i - 1];
                RnaMatchCalc(currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, false);
                currentw[0] = initverticalw[i];
                m[0] = ogcp1[i];
                if (i == imid)
                    midm[0] = m[0];
                double mi = previousw[0] + ogcp2[1];
                int mpi = 0;
                int lastj = lgth2 + 1;
                for (int j = 1; j < lastj; j++)
                {
                    double wm = previousw[j - 1];
                    double g = mi + fgcp2[j - 1];
                    if (g > wm)
                        wm = g;
                    g = previousw[j - 1] + ogcp2[j];
                    if (g >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    g = m[j] + fgcp1[i - 1];
                    if (g > wm)
                        wm = g;
                    g = previousw[j - 1] + ogcp1[i];
                    if (g >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    currentw[j] += wm;
                    wmmtx[i][j] = currentw[j];
                    wmmtx2[i][j] = m[j];
                    if (i == imid)
                    {
                        jumpbackj[j] = mp[j];
                        jumpbacki[j] = mpi;
                        midw[j] = currentw[j];
                        midm[j] = m[j];
                        midn[j] = mi;
                    }
                }
                lastverticalw[i] = currentw[lgth2 - 1];
                wmmtx2[i][lgth2] = m[lgth2 - 1];
            }
            RnaMatchCalc(initverticalw, cpmx2, cpmx1, lgth2 - 1, lgth1, doublework, intwork, true);
            RnaMatchCalc(currentw, cpmx1, cpmx2, lgth1 - 1, lgth2, doublework, intwork, true);
            for (int i = 0; i < lgth1 - 1; i++)
                initverticalw[i] += fgcp1[lgth1 - 1] + ogcp1[i + 1];
            for (int j = 0; j < lgth2 - 1; j++)
                currentw[j] += fgcp2[lgth2 - 1] + ogcp2[j + 1];
            for (int i = 0; i < lgth1 - 1; i++)
                wmmtx[i][lgth2 - 1] += fgcp1[lgth1 - 1] + ogcp1[i + 1];
            for (int j = 0; j < lgth2 - 1; j++)
                wmmtx[lgth1 - 1][j] += fgcp2[lgth2 - 1] + ogcp2[j + 1];
            for (int j = lgth2 - 1; j > 0; --j)
            {
                m[j - 1] = currentw[j] + fgcp2[lgth2 - 2];
                mp[j] = lgth1 - 1;
            }
            double firstm = -9999999.9;
            int firstmp = lgth1 - 1;
            for (int i = lgth1 - 2; i > -1; i--)
            {
                wtmp = previousw;
                previousw = currentw;
                currentw = wtmp;
                previousw[lgth2 - 1] = initverticalw[i + 1];
                RnaMatchCalc(currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, false);
                currentw[lgth2 - 1] = initverticalw[i];
                double mi = previousw[lgth2 - 1] + fgcp2[lgth2 - 2];
                int mpi = lgth2 - 1;
                int j;
                for (j = lgth2 - 2; j > -1; j--)
                {
                    double wm = previousw[j + 1];
                    int ijpi = i + 1;
                    int ijpj = j + 1;
                    double g = mi + ogcp2[j + 1];
                    if (g > wm)
                    {
                        wm = g;
                        ijpj = mpi;
                        ijpi = i + 1;
                    }
                    g = previousw[j + 1] + fgcp2[j];
                    if (g >= mi)
                    {
                        mi = g;
                        mpi = j + 1;
                    }
                    g = m[j] + ogcp1[i + 1];
                    if (g > wm)
                    {
                        wm = g;
                        ijpi = mp[j];
                        ijpj = j + 1;
                    }
                    g = previousw[j + 1] + fgcp1[i];
                    if (g >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i + 1;
                    }
                    if (i == jumpi || i == imid - 1)
                    {
                        jumpforwi[j] = ijpi;
                        jumpforwj[j] = ijpj;
                    }
                    if (i == imid)
                    {
                        midw[j] += wm;
                        midm[j + 1] += m[j];
                    }
                    if (i == imid - 1)
                        midn[j] += mi;
                    wmmtx[i][j] += wm;
                    wmmtx2[i][j] += currentw[j];
                    currentw[j] += wm;
                }
                // here j == -1: the pointers of the program stand one cell in front of the row (prept = previousw + j + 1 - 1 + 1)
                double gg = previousw[0] + fgcp1[i];
                if (firstm < gg)
                {
                    firstm = gg;
                    firstmp = i + 1;
                }
                if (i == imid)
                    midm[j + 1] += firstm;
                if (i == imid - 1)
                {
                    maxwm = midw[1];
                    jmid = 0;
                    for (j = 2; j < lgth2 - 1; j++)
                    {
                        double wm = midw[j];
                        if (wm > maxwm)
                        {
                            jmid = j;
                            maxwm = wm;
                        }
                    }
                    for (j = 0; j < lgth2 + 1; j++)
                    {
                        double wm = midm[j];
                        if (wm > maxwm)
                        {
                            jmid = j;
                            maxwm = wm;
                        }
                    }
                    double wmid = midw[jmid];
                    jumpi = imid - 1;
                    jumpj = jmid - 1;
                    if (jmid > 0 && midn[jmid - 1] > wmid)
                    {
                        jumpi = imid - 1;
                        jumpj = jumpbacki[jmid];
                        wmid = midn[jmid - 1];
                    }
                    if (midm[jmid] > wmid)
                    {
                        jumpi = jumpbackj[jmid];
                        jumpj = jmid - 1;
                        wmid = midm[jmid];
                    }
                }
                if (i == jumpi)
                {
                    if (jmid == 0)
                    {
                        jumpj = 0;
                        jmid = 1;
                        jumpi = firstmp - 1;
                        imid = firstmp;
                    }
                    else if (jmid >= lgth2)
                    {
                        jumpi = imid - 1;
                        jmid = lgth2;
                        jumpj = lgth2 - 1;
                    }
                    else
                    {
                        imid = jumpforwi[jumpj];
                        jmid = jumpforwj[jumpj];
                    }
                }
            }
            for (int i = 0; i < lgth1; i++)
                for (int j = 0; j < lgth2; j++)
                    map[i][j] = wmmtx[i][j] / maxwm;
        }
    }
}
