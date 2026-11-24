#pragma warning disable IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854, S2325, S1481, S3260
namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: MSalignmm (MSalignmm.c), the memory-saving group-group dynamic program (--memsave, groups of more than 30 000 columns)

    internal sealed partial class MafftRun
    {
        /// <summary>A C <c>char *</c> into a buffer: <c>B[O]</c> is <c>*p</c>; the strings of the program end with a NUL (0).</summary>
        private readonly struct CPtr
        {
            internal readonly char[] B;
            internal readonly int O;

            internal CPtr(char[] b, int o)
            {
                B = b;
                O = o;
            }

            internal CPtr Plus(int n) => new(B, O + n);

            internal int Len()
            {
                int k = O;
                while (B[k] != 0)
                    k++;
                return k - O;
            }

            internal void Append(string src, int start, int count)
            {
                int at = O + Len();
                for (int k = 0; k < count; k++)
                    B[at + k] = src[start + k];
                B[at + count] = (char)0;
            }

            internal void AppendChar(char c, int count)
            {
                int at = O + Len();
                for (int k = 0; k < count; k++)
                    B[at + k] = c;
                B[at + count] = (char)0;
            }

            internal string Text() => new(B, O, Len());
        }

        private static CPtr[] Plus(CPtr[] rows, int n)
        {
            var r = new CPtr[rows.Length];
            for (int i = 0; i < rows.Length; i++)
                r[i] = rows[i].Plus(n);
            return r;
        }

        /// <summary><c>match_calc_alphabet_seq</c> (MSalignmm.c): the scores of the letters of column <paramref name="i1"/> of <paramref name="cpmx1"/> with the columns <c>start2 .. start2 + lgth2 - 1</c> of <paramref name="cpmx2"/>.</summary>
        private void MatchCalcAlphabetSeq(double[][] dyn, double[] match, double[][] cpmx1, double[][] cpmx2, int i1, int start2, int lgth2, double[][] doublework, int[][] intwork, bool initialize)
        {
            int na = G.Nalphabets;
            if (initialize)
            {
                for (int j = 0, p = start2; j < lgth2; j++, p++)
                {
                    int count = 0;
                    for (int l = 0; l < na; l++)
                    {
                        if (cpmx2[l][p] != 0.0)
                        {
                            doublework[j][count] = cpmx2[l][p];
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

        /// <summary><c>Atracking</c> of MSalignmm.c: the trace back of one block; the gap tables go to <paramref name="mgt1"/> / <paramref name="mgt2"/>, the aligned rows to <paramref name="mseq1"/> / <paramref name="mseq2"/>.</summary>
        private void MsAtracking(double[] lasthorizontalw, double[] lastverticalw, string[] seq1, string[] seq2, CPtr[] mseq1, CPtr[] mseq2, int[][] ijp, int icyc, int jcyc, int ist, int ien, int jst, int jen, int fulllen1, int fulllen2, int tailgp, CPtr mgt1, CPtr mgt2)
        {
            int lgth1 = ien - ist + 1;
            int lgth2 = jen - jst + 1;
            var gt1 = new char[lgth1 + lgth2 + 1];
            var gt2 = new char[lgth1 + lgth2 + 1];
            int p1 = lgth1 + lgth2, p2 = lgth1 + lgth2;
            if (tailgp != 1 && (ien == fulllen1 - 1 || jen == fulllen2 - 1))
            {
                double wm = lastverticalw[0];
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
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = i + 1;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = -(j + 1);
            int iin = lgth1, jin = lgth2;
            int klim = lgth1 + lgth2;
            for (int k = 0; k <= klim; k++)
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
                while (--l != 0)
                {
                    gt1[--p1] = 'o';
                    gt2[--p2] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l != 0)
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
            }
            int klen = lgth1 + lgth2 - p1;
            string table1 = new(gt1, p1, klen);
            string table2 = new(gt2, p2, lgth1 + lgth2 - p2);
            bool ngap1 = table1.Contains('-');
            bool ngap2 = table2.Contains('-');
            for (int i = 0; i < icyc; i++)
                WriteRow(mseq1[i], ngap1 ? GapIreru(seq1[i].Substring(ist), table1, G.NewGapChar) : seq1[i].Substring(ist, klen));
            for (int j = 0; j < jcyc; j++)
                WriteRow(mseq2[j], ngap2 ? GapIreru(seq2[j].Substring(jst), table2, G.NewGapChar) : seq2[j].Substring(jst, klen));
            WriteRow(mgt1, table1);
            WriteRow(mgt2, table2);
        }

        private static void WriteRow(CPtr row, string s)
        {
            for (int k = 0; k < s.Length; k++)
                row.B[row.O + k] = s[k];
            row.B[row.O + s.Length] = (char)0;
        }

        /// <summary><c>MSalignmm_tanni</c>: the plain dynamic program of a block that is short in one direction.</summary>
        private double MsTanni(double[][] dyn, int icyc, int jcyc, double[] eff1, double[] eff2, string[] seq1, string[] seq2, double[][] cpmx1pt, double[][] cpmx2pt, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2,
            CPtr[] mseq1, CPtr[] mseq2, CPtr mgt1, CPtr mgt2, double[][] gapinfo, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g)
        {
            int na = G.Nalphabets;
            var ogcp1 = gapinfo[0].AsSpan(ist);
            var fgcp1 = gapinfo[1].AsSpan(ist);
            var ogcp2 = gapinfo[2].AsSpan(jst);
            var fgcp2 = gapinfo[3].AsSpan(jst);
            var gapfreq1f = gapinfo[4].AsSpan(ist);
            var gapfreq2f = gapinfo[5].AsSpan(jst);
            double headgapfreq1 = ist > 0 ? gapinfo[4][ist - 1] : headgapfreq1_g;
            double headgapfreq2 = jst > 0 ? gapinfo[5][jst - 1] : headgapfreq2_g;
            int lgth1 = ien - ist + 1;
            int lgth2 = jen - jst + 1;
            int ll1 = lgth1 + 100;
            int ll2 = lgth2 + 100;
            var w1 = new double[ll2 + 2];
            var w2 = new double[ll2 + 2];
            var initverticalw = new double[ll1 + 2];
            var lastverticalw = new double[ll1 + 2];
            var m = new double[ll2 + 2];
            var mp = new int[ll2 + 2];
            int dim = Math.Max(ll1, ll2) + 2;
            var doublework = new double[dim][];
            var intwork = new int[dim][];
            for (int i = 0; i < dim; i++)
            {
                doublework[i] = new double[na + 1];
                intwork[i] = new int[na + 1];
            }
            var ijp = new int[ll1 + 1][];
            for (int i = 0; i < ijp.Length; i++)
                ijp[i] = new int[ll2 + 1];
            double[] currentw = w1, previousw = w2;
            MatchCalcAlphabetSeq(dyn, initverticalw, cpmx2pt, cpmx1pt, jst, ist, lgth1, doublework, intwork, true);
            MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist, jst, lgth2, doublework, intwork, true);
            if (headgp != 0 || ist != 0)
            {
                for (int i = 1; i < lgth1 + 1; i++)
                    initverticalw[i] += ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2f[0];
            }
            if (headgp != 0 || jst != 0)
            {
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] += ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1f[0];
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2f[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int lasti = tailgp != 0 || jen != fulllen2 - 1 ? lgth1 + 1 : lgth1;
            double wm = 0.0;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist + i, jst, lgth2, doublework, intwork, false);
                currentw[0] = initverticalw[i];
                double mi = previousw[0] + ogcp2[1] * gapfreq1f[i - 1];
                int mpi = 0;
                var ijpRow = ijp[i];
                int lastj = lgth2 + 1;
                for (int j = 1; j < lastj; j++)
                {
                    wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    double g = mi + fgcp2[j - 1] * gapfreq1f[i];
                    if (g > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    g = previousw[j - 1] + ogcp2[j] * gapfreq1f[i - 1];
                    if (g >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    g = m[j] + fgcp1[i - 1] * gapfreq2f[j];
                    if (g > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    g = previousw[j - 1] + ogcp1[i] * gapfreq2f[j - 1];
                    if (g >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
            }
            MsAtracking(currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, ist, ien, jst, jen, fulllen1, fulllen2, tailgp, mgt1, mgt2);
            return wm;
        }

        /// <summary><c>MSalignmm_rec</c>: the recursive (Hirschberg-type) dynamic program; the block is cut at the middle row, the alignments of the two halves are joined.</summary>
        private double MsRec(double[][] dyn, int icyc, int jcyc, double[] eff1, double[] eff2, string[] seq1, string[] seq2, double[][] cpmx1pt, double[][] cpmx2pt, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2,
            CPtr[] mseq1, CPtr[] mseq2, CPtr mgt1, CPtr mgt2, int depth, double[][] gapinfo, int headgp, int tailgp, double headgapfreq1_g, double headgapfreq2_g)
        {
            int na = G.Nalphabets;
            double value;
            var ogcp1 = gapinfo[0].AsSpan(ist);
            var fgcp1 = gapinfo[1].AsSpan(ist);
            var ogcp2 = gapinfo[2].AsSpan(jst);
            var fgcp2 = gapinfo[3].AsSpan(jst);
            var gapfreq1f = gapinfo[4].AsSpan(ist);
            var gapfreq2f = gapinfo[5].AsSpan(jst);
            double headgapfreq1 = ist > 0 ? gapinfo[4][ist - 1] : headgapfreq1_g;
            double headgapfreq2 = jst > 0 ? gapinfo[5][jst - 1] : headgapfreq2_g;
            depth++;
            int lgth1 = ien - ist + 1;
            int lgth2 = jen - jst + 1;
            if (lgth2 <= 0)
            {
                for (int i = 0; i < icyc; i++)
                    WriteRow(mseq1[i], seq1[i].Substring(ist, lgth1));
                var gapRow = new string(G.NewGapChar, lgth1);
                for (int i = 0; i < jcyc; i++)
                    WriteRow(mseq2[i], gapRow);
                WriteRow(mgt1, new string('o', lgth1));
                WriteRow(mgt2, new string('-', lgth1));
                return 0.0;
            }
            var aseq1 = (CPtr[])mseq1.Clone();
            var aseq2 = (CPtr[])mseq2.Clone();
            CPtr agt1 = mgt1, agt2 = mgt2;
            if (lgth1 < DpTanni || lgth2 < DpTanni)
                return MsTanni(dyn, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist, ien, jst, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g);
            int ll1 = lgth1 + 100;
            int ll2 = lgth2 + 100;
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
            int dim = Math.Max(ll1, ll2) + 2;
            var doublework = new double[dim][];
            var intwork = new int[dim][];
            for (int i = 0; i < dim; i++)
            {
                doublework[i] = new double[na];
                intwork[i] = new int[na + 1];
            }
            double[] currentw = w1, previousw = w2;
            MatchCalcAlphabetSeq(dyn, initverticalw, cpmx2pt, cpmx1pt, jst, ist, lgth1, doublework, intwork, true);
            MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist, jst, lgth2, doublework, intwork, true);
            for (int i = 1; i < lgth1 + 1; i++)
                initverticalw[i] += ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2f[0];
            for (int j = 1; j < lgth2 + 1; j++)
                currentw[j] += ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1f[0];
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2f[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int imid = (int)(lgth1 * 0.5);
            int jumpi = 0;
            int jumpj = 0;
            int jmid = 0;
            double wm = 0.0, g;
            for (int i = 1; i <= imid; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist + i, jst, lgth2, doublework, intwork, false);
                currentw[0] = initverticalw[i];
                m[0] = ogcp1[i];
                if (i == imid)
                    midm[0] = m[0];
                double mi = previousw[0] + ogcp2[1] * gapfreq1f[i - 1];
                int mpi = 0;
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    wm = previousw[j - 1];
                    g = mi + fgcp2[j - 1] * gapfreq1f[i];
                    if (g > wm)
                        wm = g;
                    g = previousw[j - 1] + ogcp2[j] * gapfreq1f[i - 1];
                    if (g >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    g = m[j] + fgcp1[i - 1] * gapfreq2f[j];
                    if (g > wm)
                        wm = g;
                    g = previousw[j - 1] + ogcp1[i] * gapfreq2f[j - 1];
                    if (g >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    currentw[j] += wm;
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
            }
            MatchCalcAlphabetSeq(dyn, initverticalw, cpmx2pt, cpmx1pt, jst + lgth2 - 1, ist, lgth1, doublework, intwork, true);
            MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist + lgth1 - 1, jst, lgth2, doublework, intwork, true);
            for (int i = 0; i < lgth1 - 1; i++)
                initverticalw[i] += fgcp1[lgth1 - 1] * gapfreq2f[lgth2] + ogcp1[i + 1] * gapfreq2f[lgth2 - 1];
            for (int j = 0; j < lgth2 - 1; j++)
                currentw[j] += fgcp2[lgth2 - 1] * gapfreq1f[lgth1] + ogcp2[j + 1] * gapfreq1f[lgth1 - 1];
            for (int j = lgth2 - 1; j > -1; --j)
            {
                m[j] = currentw[j + 1] + fgcp1[lgth1 - 2] * gapfreq2f[j + 1];
                mp[j] = lgth1 - 1;
            }
            double firstm = -9999999.9;
            int firstmp = lgth1;
            for (int i = lgth1 - 2; i > -1; i--)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[lgth2 - 1] = initverticalw[i + 1];
                MatchCalcAlphabetSeq(dyn, currentw, cpmx1pt, cpmx2pt, ist + i, jst, lgth2, doublework, intwork, false);
                currentw[lgth2 - 1] = initverticalw[i];
                double mi = previousw[lgth2 - 1] + fgcp2[lgth2 - 2] * gapfreq1f[i + 1];
                int mpi = lgth2 - 1;
                int j;
                for (j = lgth2 - 2; j > -1; j--)
                {
                    wm = previousw[j + 1];
                    int ijpi = i + 1;
                    int ijpj = j + 1;
                    g = mi + ogcp2[j + 1] * gapfreq1f[i];
                    if (g > wm)
                    {
                        wm = g;
                        ijpj = mpi;
                        ijpi = i + 1;
                    }
                    g = previousw[j + 1] + fgcp2[j] * gapfreq1f[i + 1];
                    if (g >= mi)
                    {
                        mi = g;
                        mpi = j + 1;
                    }
                    g = m[j] + ogcp1[i + 1] * gapfreq2f[j];
                    if (g > wm)
                    {
                        wm = g;
                        ijpi = mp[j];
                        ijpj = j + 1;
                    }
                    g = previousw[j + 1] + fgcp1[i] * gapfreq2f[j + 1];
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
                    currentw[j] += wm;
                }
                // the pointers prept / curpt stand on column 0 after the loop
                g = previousw[0] + fgcp1[i];
                if (firstm < g)
                {
                    firstm = g;
                    firstmp = i + 1;
                }
                if (i == imid)
                    midm[j + 1] += firstm;
                if (i == imid - 1)
                {
                    double maxwm = midw[1];
                    jmid = 0;
                    for (j = 2; j < lgth2 - 1; j++)
                    {
                        wm = midw[j];
                        if (wm > maxwm)
                        {
                            jmid = j;
                            maxwm = wm;
                        }
                    }
                    for (j = 0; j < lgth2 + 1; j++)
                    {
                        wm = midm[j];
                        if (wm > maxwm)
                        {
                            jmid = j;
                            maxwm = wm;
                        }
                    }
                    wm = midw[jmid];
                    jumpi = imid - 1;
                    jumpj = jmid - 1;
                    if (jmid > 0 && midn[jmid - 1] > wm)
                    {
                        jumpi = imid - 1;
                        jumpj = jumpbacki[jmid];
                        wm = midn[jmid - 1];
                    }
                    if (midm[jmid] > wm)
                    {
                        jumpi = jumpbackj[jmid];
                        jumpj = jmid - 1;
                        wm = midm[jmid];
                    }
                }
                if (i == jumpi)
                {
                    if (jmid == 0)
                    {
                        jumpj = 0;
                        jmid = 1;
                        if (imid < firstmp - 1)
                        {
                            jumpi = firstmp;
                            imid = firstmp + 1;
                        }
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
                        if (imid == jumpi)
                            jumpi = imid - 1;
                    }
                    break;
                }
            }
            value = MsRec(dyn, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist, ist + jumpi, jst, jst + jumpj, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, depth, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g);
            int len = mseq1[0].Len();
            int l = jmid - jumpj - 1;
            if (l > 0)
            {
                for (int i = 0; i < icyc; i++)
                    mseq1[i].AppendChar(G.NewGapChar, l);
                for (int j = 0; j < jcyc; j++)
                    mseq2[j].Append(seq2[j], jst + jumpj + 1, l);
                for (int k = 0; k < l; k++)
                    mgt1.B[mgt1.O + len + k] = '-';
                mgt1.B[mgt1.O + len + l] = (char)0;
                for (int k = 0; k < l; k++)
                    mgt2.B[mgt2.O + len + k] = 'o';
                // the program ends mgt1 again here (mgt1[j] = 0) and leaves mgt2 as it is
                value += ogcp2[jumpj + 1] + fgcp2[jmid - 1];
            }
            len = mseq1[0].Len();
            l = imid - jumpi - 1;
            if (l > 0)
            {
                for (int i = 0; i < icyc; i++)
                    mseq1[i].Append(seq1[i], ist + jumpi + 1, l);
                for (int j = 0; j < jcyc; j++)
                    mseq2[j].AppendChar(G.NewGapChar, l);
                for (int k = 0; k < l; k++)
                    mgt1.B[mgt1.O + len + k] = 'o';
                mgt1.B[mgt1.O + len + l] = (char)0;
                for (int k = 0; k < l; k++)
                    mgt2.B[mgt2.O + len + k] = '-';
                value += ogcp1[jumpi + 1] + fgcp1[imid - 1];
            }
            int alnlen = aseq1[0].Len();
            aseq1 = Plus(aseq1, alnlen);
            aseq2 = Plus(aseq2, alnlen);
            agt1 = agt1.Plus(alnlen);
            agt2 = agt2.Plus(alnlen);
            value += MsRec(dyn, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1pt, cpmx2pt, ist + imid, ien, jst + jmid, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, agt1, agt2, depth, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g);
            return value;
        }

        private const int DpTanni = 100;

        /// <summary>
        /// <c>MSalignmm</c> (<c>MSalignmm.c</c>): the memory-saving version of <see cref="AAlign"/> (the same gap penalties from the profiles of the groups, the same history of profiles); <paramref name="seq1"/> / <paramref name="seq2"/> receive the aligned rows.
        /// </summary>
        internal double MsAlign(double[][] dyn, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, int alloclen, string? sgap1, string? sgap2, string? egap1, string? egap2, int headgp, int tailgp,
            double[][][]? hist, int child0, int child1, int resultIndex, double orieff1, double orieff2)
        {
            int na = G.Nalphabets;
            double fpenalty = G.Penalty;
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            int ll1 = lgth1 + 100;
            int ll2 = lgth2 + 100;
            var mseqBuf1 = new CPtr[icyc];
            var mseqBuf2 = new CPtr[jcyc];
            for (int i = 0; i < icyc; i++)
                mseqBuf1[i] = new CPtr(new char[ll1 + ll2 + 1], 0);
            for (int j = 0; j < jcyc; j++)
                mseqBuf2[j] = new CPtr(new char[ll1 + ll2 + 1], 0);
            var mgt1 = new CPtr(new char[ll1 + ll2 + 1], 0);
            var mgt2 = new CPtr(new char[ll1 + ll2 + 1], 0);
            var ogcp1 = new double[ll1 + 2];
            var ogcp2 = new double[ll2 + 2];
            var fgcp1 = new double[ll1 + 2];
            var fgcp2 = new double[ll2 + 2];
            var work = new ProfileWork
            {
                Cpmx1 = NewMatrix(na, ll1 + 2), Cpmx2 = NewMatrix(na, ll2 + 2),
                Gapfreq1 = new double[ll1 + 2], Gapfreq2 = new double[ll2 + 2],
                Ogcp1o = new double[ll1 + 2], Ogcp2o = new double[ll2 + 2], Fgcp1o = new double[ll1 + 2], Fgcp2o = new double[ll2 + 2],
            };
            var prof = PrepareProfiles(work, seq1, seq2, eff1, eff2, icyc, jcyc, lgth1, lgth2, sgap1, sgap2, egap1, egap2, hist, child0, child1);
            for (int i = 0; i < lgth1; i++)
            {
                ogcp1[i] = 0.5 * (1.0 - prof.Ogcp1o[i]) * fpenalty * prof.Gapfreq1[i];
                fgcp1[i] = 0.5 * (1.0 - prof.Fgcp1o[i]) * fpenalty * prof.Gapfreq1[i];
            }
            for (int i = 0; i < lgth2; i++)
            {
                ogcp2[i] = 0.5 * (1.0 - prof.Ogcp2o[i]) * fpenalty * prof.Gapfreq2[i];
                fgcp2[i] = 0.5 * (1.0 - prof.Fgcp2o[i]) * fpenalty * prof.Gapfreq2[i];
            }
            var gapinfo = new[] { ogcp1, fgcp1, ogcp2, fgcp2, prof.Gapfreq1, prof.Gapfreq2 };
            double wm = MsRec(dyn, icyc, jcyc, eff1, eff2, seq1, seq2, prof.Cpmx1, prof.Cpmx2, 0, lgth1 - 1, 0, lgth2 - 1, alloclen, lgth1, lgth2, mseqBuf1, mseqBuf2, mgt1, mgt2, 0, gapinfo, headgp, tailgp, prof.Head1, prof.Head2);
            if (hist != null && G.SpecificityConsideration == 0.0)
            {
                if (icyc + jcyc > 20)
                {
                    double totaleff1 = 0.0, totaleff2 = 0.0;
                    for (int i = 0; i < icyc; i++)
                        totaleff1 += eff1[i];
                    for (int j = 0; j < jcyc; j++)
                        totaleff2 += eff2[j];
                    if (Math.Abs(totaleff1 - 1.0) > 0.001 || Math.Abs(totaleff2 - 1.0) > 0.001)
                        throw new MafftProgramFailureException("Warning: rounding error may be large.  totaleff1 = " + totaleff1.ToString("F40", System.Globalization.CultureInfo.InvariantCulture).PadLeft(50) + "\nWarning: rounding error may be large.  totaleff2 = " + totaleff2.ToString("F40", System.Globalization.CultureInfo.InvariantCulture).PadLeft(50) + "\n");
                    totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
                    totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);
                    // the program takes the arrays ogcp1o / fgcp1o of its own (zero when a child's profile was used), not the counts of the children (ogcp1opt)
                    string t1 = mgt1.Text(), t2 = mgt2.Text();
                    int alen = t1.Length;
                    t2 = t2.Length >= alen ? t2 : t2.PadRight(alen, '\0');
                    var res = new double[na + 3][];
                    CreateCpmxResult(res, alen, totaleff1, totaleff2, prof.Cpmx1, prof.Cpmx2, t1, t2);
                    res[na] = CreateGapFreqResult(alen, totaleff1, totaleff2, prof.Gapfreq1, prof.Gapfreq2, t1, t2);
                    res[na + 1] = CreateOgResult(alen, totaleff1, totaleff2, work.Ogcp1o, work.Ogcp2o, prof.Gapfreq1, prof.Gapfreq2, t1, t2);
                    res[na + 2] = CreateFgResult(alen, totaleff1, totaleff2, work.Fgcp1o, work.Fgcp2o, prof.Gapfreq1, prof.Gapfreq2, t1, t2);
                    hist[resultIndex] = res;
                }
                else
                {
                    hist[resultIndex] = null!;
                }
                if (prof.History1)
                    hist[child0] = null!;
                if (prof.History2)
                    hist[child1] = null!;
            }
            for (int i = 0; i < icyc; i++)
                seq1[i] = mseqBuf1[i].Text();
            for (int j = 0; j < jcyc; j++)
                seq2[j] = mseqBuf2[j].Text();
            return wm;
        }

        private static double[][] NewMatrix(int rows, int cols)
        {
            var r = new double[rows][];
            for (int i = 0; i < rows; i++)
                r[i] = new double[cols];
            return r;
        }
    }

    #endregion
}
