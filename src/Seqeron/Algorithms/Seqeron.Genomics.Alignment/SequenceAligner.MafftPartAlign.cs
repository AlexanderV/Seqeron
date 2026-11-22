using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: partA__align and Falign_localhom (the dynamic program of the iterative refinement with a pair library)

    internal sealed partial class MafftRun
    {
        // ---- the static work arrays of partA__align (partSalignmm.c): separate from those of A__align (Salignmm.c), re-used and zeroed only when a sequence is longer than any before
        private int _pOrlgth1, _pOrlgth2;
        private double[] _pOgcp1 = Array.Empty<double>(), _pOgcp2 = Array.Empty<double>(), _pFgcp1 = Array.Empty<double>(), _pFgcp2 = Array.Empty<double>();
        private double[] _pGapfreq1 = Array.Empty<double>(), _pGapfreq2 = Array.Empty<double>();
        private double[][] _pCpmx1 = Array.Empty<double[]>(), _pCpmx2 = Array.Empty<double[]>();

        /// <summary>The value of <c>impmatch</c> that the last call of <c>partA__align</c> that set it left (<c>Falign_localhom</c> adds its local variable even after a call that returned early).</summary>
        private double _lastImpMatch;

        /// <summary>
        /// <c>partA__align</c> (<c>partSalignmm.c</c>) with <c>constraint</c>: the dynamic program of two groups with the importance matrix of the pair library added to the match scores. The group rows were freed of their common gaps
        /// (<paramref name="gapmap1"/> / <paramref name="gapmap2"/> give the column of the uncompacted segment of every column) and <paramref name="start1"/> / <paramref name="start2"/> are the first columns of the segment in the
        /// alignment the importance matrix was built from. The rows are replaced by the aligned rows; <paramref name="impmatch"/> is the sum of the importance of the cells of the path.
        /// </summary>
        private void PartAAlign(string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, int alloclen, out double impmatch, int start1, int start2, int[] gapmap1, int[] gapmap2,
            string sgap1, string sgap2, string egap1, string egap2)
        {
            impmatch = _lastImpMatch;
            int na = G.Nalphabets;
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            if (lgth1 == 0 && lgth2 == 0)
                return;
            if (lgth1 == 0)
            {
                for (int i = 0; i < icyc; i++)
                    seq1[i] = new string('-', lgth2);
                return;
            }
            if (lgth2 == 0)
            {
                for (int i = 0; i < jcyc; i++)
                    seq2[i] = new string('-', lgth1);
                return;
            }
            if (G.TryWarp)
                throw new NotSupportedException("--allowshift / --shiftpenalty (the warp of partA__align) is not available in this port.");
            if (lgth1 > _pOrlgth1 || lgth2 > _pOrlgth2)
            {
                int ll1 = Math.Max((int)(1.3 * lgth1), _pOrlgth1) + 100;
                int ll2 = Math.Max((int)(1.3 * lgth2), _pOrlgth2) + 100;
                _pOgcp1 = new double[ll1 + 2];
                _pOgcp2 = new double[ll2 + 2];
                _pFgcp1 = new double[ll1 + 2];
                _pFgcp2 = new double[ll2 + 2];
                _pCpmx1 = new double[na][];
                _pCpmx2 = new double[na][];
                for (int i = 0; i < na; i++)
                {
                    _pCpmx1[i] = new double[ll1 + 2];
                    _pCpmx2[i] = new double[ll2 + 2];
                }
                _pGapfreq1 = new double[ll1 + 2];
                _pGapfreq2 = new double[ll2 + 2];
                _pOrlgth1 = ll1 - 100;
                _pOrlgth2 = ll2 - 100;
            }
            double fpenalty = G.Penalty;
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
            double[][] cpmx1pt = _pCpmx1, cpmx2pt = _pCpmx2;
            double[] ogcp1 = _pOgcp1, ogcp2 = _pOgcp2, fgcp1 = _pFgcp1, fgcp2 = _pFgcp2, gapfreq1 = _pGapfreq1, gapfreq2 = _pGapfreq2;
            CpmxCalcNew(seq1, cpmx1pt, eff1, lgth1, icyc);
            CpmxCalcNew(seq2, cpmx2pt, eff2, lgth2, jcyc);
            double headgapfreq1, headgapfreq2;
            NewOpeningGapCount(ogcp1, icyc, seq1, eff1, lgth1, sgap1);
            NewOpeningGapCount(ogcp2, jcyc, seq2, eff2, lgth2, sgap2);
            NewFinalGapCount(fgcp1, icyc, seq1, eff1, lgth1);
            NewFinalGapCount(fgcp2, jcyc, seq2, eff2, lgth2);
            headgapfreq1 = OutGapCount(icyc, sgap1, eff1);
            headgapfreq2 = OutGapCount(jcyc, sgap2, eff2);
            gapfreq1[lgth1] = OutGapCount(icyc, egap1, eff1);
            gapfreq2[lgth2] = OutGapCount(jcyc, egap2, eff2);
            if (G.LegacyGapCost == 0)
            {
                GapCountF(gapfreq1, seq1, icyc, eff1, lgth1);
                GapCountF(gapfreq2, seq2, jcyc, eff2, lgth2);
                for (int i = 0; i < lgth1 + 1; i++)
                    gapfreq1[i] = 1.0 - gapfreq1[i];
                for (int i = 0; i < lgth2 + 1; i++)
                    gapfreq2[i] = 1.0 - gapfreq2[i];
                headgapfreq1 = 1.0 - headgapfreq1;
                headgapfreq2 = 1.0 - headgapfreq2;
            }
            else
            {
                for (int i = 0; i < lgth1 + 1; i++)
                    gapfreq1[i] = 1.0;
                for (int i = 0; i < lgth2 + 1; i++)
                    gapfreq2[i] = 1.0;
                headgapfreq1 = 1.0;
                headgapfreq2 = 1.0;
            }
            for (int i = 0; i < lgth1; i++)
            {
                ogcp1[i] = 0.5 * (1.0 - ogcp1[i]) * fpenalty * gapfreq1[i];
                fgcp1[i] = 0.5 * (1.0 - fgcp1[i]) * fpenalty * gapfreq1[i];
            }
            for (int i = 0; i < lgth2; i++)
            {
                ogcp2[i] = 0.5 * (1.0 - ogcp2[i]) * fpenalty * gapfreq2[i];
                fgcp2[i] = 0.5 * (1.0 - fgcp2[i]) * fpenalty * gapfreq2[i];
            }
            var dyn = G.NDisConsweightMulti;
            var impmtx = _partImpMtx!;
            double[] currentw = w1, previousw = w2;
            MatchCalc(dyn, initverticalw, cpmx2pt, cpmx1pt, 0, lgth1, doublework, intwork, true);
            {
                int j1 = gapmap2[0] + start2;
                for (int i = 0; i < lgth1; i++)
                    initverticalw[i] += impmtx[start1 + gapmap1[i]][j1];
            }
            MatchCalc(dyn, currentw, cpmx1pt, cpmx2pt, 0, lgth2, doublework, intwork, true);
            {
                var row = impmtx[gapmap1[0] + start1];
                for (int j = 0; j < lgth2; j++)
                    currentw[j] += row[start2 + gapmap2[j]];
            }
            int outgap = G.Outgap;
            if (outgap == 1)
            {
                for (int i = 1; i < lgth1 + 1; i++)
                    initverticalw[i] += ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2[0];
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] += ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1[0];
            }
            else
            {
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] -= G.Offset * j / 2.0;
                for (int i = 1; i < lgth1 + 1; i++)
                    initverticalw[i] -= G.Offset * i / 2.0;
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1] + ogcp1[1] * gapfreq2[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int lasti = outgap != 0 ? lgth1 + 1 : lgth1;
            int lastj = lgth2 + 1;
            double wm = 0.0;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                MatchCalc(dyn, currentw, cpmx1pt, cpmx2pt, i, lgth2, doublework, intwork, false);
                {
                    var row = impmtx[gapmap1[i] + start1];
                    for (int j = 0; j < lgth2; j++)
                        currentw[j] += row[start2 + gapmap2[j]];
                }
                currentw[0] = initverticalw[i];
                double mi = previousw[0] + ogcp2[1] * gapfreq1[i - 1];
                int mpi = 0;
                var ijpRow = ijp[i];
                double fgcp1va = fgcp1[i - 1];
                double ogcp1va = ogcp1[i];
                double gf1va = gapfreq1[i];
                double gf1vapre = gapfreq1[i - 1];
                for (int j = 1; j < lastj; j++)
                {
                    double g;
                    wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    if ((g = mi + fgcp2[j - 1] * gf1va) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = previousw[j - 1] + ogcp2[j] * gf1vapre) > mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    if ((g = m[j] + fgcp1va * gapfreq2[j]) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    if ((g = previousw[j - 1] + ogcp1va * gapfreq2[j - 1]) > m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
            }
            if (outgap == 0)
            {
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] -= G.Offset * (lgth2 - j) / 2.0;
                for (int i = 1; i < lgth1 + 1; i++)
                    lastverticalw[i] -= G.Offset * (lgth1 - i / 2.0);
            }

            // Atracking_localhom of partSalignmm.c
            if (outgap != 1)
            {
                double wmt = lastverticalw[0];
                for (int i = 0; i < lgth1; i++)
                {
                    if (lastverticalw[i] >= wmt)
                    {
                        wmt = lastverticalw[i];
                        ijp[lgth1][lgth2] = +(lgth1 - i);
                    }
                }
                for (int j = 0; j < lgth2; j++)
                {
                    if (currentw[j] >= wmt)
                    {
                        wmt = currentw[j];
                        ijp[lgth1][lgth2] = -(lgth2 - j);
                    }
                }
            }
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = i + 1;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = -(j + 1);
            var a1 = new char[icyc][];
            var a2 = new char[jcyc][];
            for (int i = 0; i < icyc; i++)
                a1[i] = new char[lgth1 + lgth2 + 1];
            for (int j = 0; j < jcyc; j++)
                a2[j] = new char[lgth1 + lgth2 + 1];
            int p = lgth1 + lgth2;
            int iin = lgth1, jin = lgth2;
            double impwm = 0.0;
            int limk = lgth1 + lgth2 + 1;
            int kk = 0;
            while (kk < limk)
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
                    p--;
                    for (int i = 0; i < icyc; i++)
                        a1[i][p] = seq1[i][ifi + l];
                    for (int j = 0; j < jcyc; j++)
                        a2[j][p] = '-';
                    kk++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    p--;
                    for (int i = 0; i < icyc; i++)
                        a1[i][p] = '-';
                    for (int j = 0; j < jcyc; j++)
                        a2[j][p] = seq2[j][jfi + l];
                    kk++;
                }
                if (iin != lgth1 && jin != lgth2)
                    impwm += impmtx[gapmap1[iin] + start1][gapmap2[jin] + start2];
                if (iin <= 0 || jin <= 0)
                    break;
                p--;
                for (int i = 0; i < icyc; i++)
                    a1[i][p] = seq1[i][ifi];
                for (int j = 0; j < jcyc; j++)
                    a2[j][p] = seq2[j][jfi];
                kk++;
                iin = ifi;
                jin = jfi;
                kk++;
            }
            impmatch = impwm;
            _lastImpMatch = impwm;
            int resultlen = lgth1 + lgth2 - p;
            if (alloclen < resultlen)
                throw new InvalidOperationException("LENGTH OVER!");
            for (int i = 0; i < icyc; i++)
                seq1[i] = new string(a1[i], p, resultlen);
            for (int j = 0; j < jcyc; j++)
                seq2[j] = new string(a2[j], p, resultlen);
        }

        /// <summary><c>commongappick_record</c>: the columns that are not a gap in every row of the group, and the column of the original of each (the terminating NUL is a column too).</summary>
        private static int[] CommonGapPickRecord(string[] seq)
        {
            int nseq = seq.Length;
            int len = seq[0].Length;
            var map = new List<int>(len + 1);
            var sbs = new System.Text.StringBuilder[nseq];
            for (int j = 0; j < nseq; j++)
                sbs[j] = new System.Text.StringBuilder(len);
            for (int i = 0; i < len; i++)
            {
                bool keep = false;
                for (int j = 0; j < nseq; j++)
                {
                    if (seq[j][i] != '-')
                    {
                        keep = true;
                        break;
                    }
                }
                if (keep)
                {
                    for (int j = 0; j < nseq; j++)
                        sbs[j].Append(seq[j][i]);
                    map.Add(i);
                }
            }
            map.Add(len);
            for (int j = 0; j < nseq; j++)
                seq[j] = sbs[j].ToString();
            return map.ToArray();
        }

        /// <summary>
        /// <c>Falign_localhom</c> (<c>Falign_localhom.c</c>) for <c>kobetsubunkatsu</c> (the refinement: the groups' own columns as anchors, no FFT): the alignment is cut at the centres of the alignable regions, every piece is freed of its
        /// common gaps and aligned by <see cref="PartAAlign"/>; <paramref name="totalimpmatch"/> is the sum of the importance of the paths.
        /// </summary>
        private void FalignLocalHom(string[] seq1, string[] seq2, double[] eff1, double[] eff2, int clus1, int clus2, int alloclen, out double totalimpmatch)
        {
            if (G.Kobetsubunkatsu == 0)
                throw new NotSupportedException("The FFT variant of Falign_localhom is not available in this port.");
            int len1 = seq1[0].Length, len2 = seq2[0].Length;
            int count = 0;
            int tmpint = AlignableRegion(clus1, clus2, seq1, 0, seq2, 0, eff1, eff2, _fSegment, count);
            if (count + tmpint > MaxSeg - 3)
                throw new InvalidOperationException("TOO MANY SEGMENTS.");
            while (tmpint-- > 0)
            {
                var sg = _fSegment[count];
                _fSegment1[count].Start = sg.Start;
                _fSegment1[count].End = sg.End;
                _fSegment1[count].Center = sg.Center;
                _fSegment1[count].Score = sg.Score;
                _fSegment2[count].Start = sg.Start;
                _fSegment2[count].End = sg.End;
                _fSegment2[count].Center = sg.Center;
                _fSegment2[count].Score = sg.Score;
                _fSegment1[count].Pair = _fSegment2[count];
                _fSegment2[count].Pair = _fSegment1[count];
                count++;
            }
            var sorted1 = new FftSegment[count];
            var sorted2 = new FftSegment[count];
            for (int i = 0; i < count; i++)
            {
                sorted1[i] = _fSegment1[i];
                sorted2[i] = _fSegment2[i];
            }
            MyMergeSort(0, count - 1, sorted1);
            MyMergeSort(0, count - 1, sorted2);
            for (int i = 0; i < count; i++)
                sorted1[i].Number = i;
            for (int i = 0; i < count; i++)
                sorted2[i].Number = i;
            int[] cut1 = _fCut1, cut2 = _fCut2;
            for (int i = 0; i < count; i++)
            {
                cut1[i + 1] = sorted1[i].Center;
                cut2[i + 1] = sorted2[i].Center;
            }
            cut1[0] = 0;
            cut2[0] = 0;
            cut1[count + 1] = len1;
            cut2[count + 1] = len2;
            count += 2;
            var result1 = new System.Text.StringBuilder[clus1];
            var result2 = new System.Text.StringBuilder[clus2];
            for (int j = 0; j < clus1; j++)
                result1[j] = new System.Text.StringBuilder();
            for (int j = 0; j < clus2; j++)
                result2[j] = new System.Text.StringBuilder();
            int totallen = 0;
            totalimpmatch = 0.0;
            for (int i = 0; i < count - 1; i++)
            {
                var sgap1 = new char[clus1];
                var egap1 = new char[clus1];
                var sgap2 = new char[clus2];
                var egap2 = new char[clus2];
                if (cut1[i] != 0)
                {
                    for (int j = 0; j < clus1; j++)
                        sgap1[j] = seq1[j][cut1[i] - 1];
                    for (int j = 0; j < clus2; j++)
                        sgap2[j] = seq2[j][cut2[i] - 1];
                }
                else
                {
                    for (int j = 0; j < clus1; j++)
                        sgap1[j] = 'o';
                    for (int j = 0; j < clus2; j++)
                        sgap2[j] = 'o';
                }
                if (cut1[i + 1] != len1)
                {
                    for (int j = 0; j < clus1; j++)
                        egap1[j] = seq1[j][cut1[i + 1]];
                    for (int j = 0; j < clus2; j++)
                        egap2[j] = seq2[j][cut2[i + 1]];
                }
                else
                {
                    for (int j = 0; j < clus1; j++)
                        egap1[j] = 'o';
                    for (int j = 0; j < clus2; j++)
                        egap2[j] = 'o';
                }
                var tmpres1 = new string[clus1];
                var tmpres2 = new string[clus2];
                for (int j = 0; j < clus1; j++)
                    tmpres1[j] = seq1[j].Substring(cut1[i], cut1[i + 1] - cut1[i]);
                int[] gapmap1 = CommonGapPickRecord(tmpres1);
                for (int j = 0; j < clus2; j++)
                    tmpres2[j] = seq2[j].Substring(cut2[i], cut2[i + 1] - cut2[i]);
                int[] gapmap2 = CommonGapPickRecord(tmpres2);
                if (G.Alg != 'A')
                    throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                PartAAlign(tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, out double impmatch, cut1[i], cut2[i], gapmap1, gapmap2,
                    new string(sgap1), new string(sgap2), new string(egap1), new string(egap2));
                totalimpmatch += impmatch;
                int nl = tmpres1[0].Length;
                if (totallen + nl > alloclen)
                    throw new InvalidOperationException("LENGTH OVER in Falign");
                for (int j = 0; j < clus1; j++)
                    result1[j].Append(tmpres1[j]);
                for (int j = 0; j < clus2; j++)
                    result2[j].Append(tmpres2[j]);
                totallen += nl;
            }
            for (int j = 0; j < clus1; j++)
                seq1[j] = result1[j].ToString();
            for (int j = 0; j < clus2; j++)
                seq2[j] = result2[j].ToString();
        }
    }

    #endregion
}
