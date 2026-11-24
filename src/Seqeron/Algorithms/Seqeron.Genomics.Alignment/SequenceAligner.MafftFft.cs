using System.Globalization;
namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: Falign.c, fftFunctions.c and fft.c (the FFT anchors of the progressive stage)

    internal sealed partial class MafftRun
    {
        private const int NKouho = 20;
        private const int MaxSeg = 100000;
        private const int SegmentSize = 150;

        /// <summary><c>Segment</c> of <c>fft.h</c>: a window of the alignment of two groups with a high correlation score.</summary>
        internal sealed class FftSegment
        {
            internal int Start, End, Center, Number;
            internal double Score;
            internal FftSegment? Pair;
        }

        // ---- fft.c: the sine table and the bit reversal are re-computed when the size changes (static in the program)
        private int _fftLastN;
        private double[] _fftSinTbl = Array.Empty<double>();
        private int[] _fftBitrev = Array.Empty<int>();

        private static void MakeSinTbl(int n, double[] sintbl)
        {
            int n2 = n / 2, n4 = n / 4, n8 = n / 8;
            double t = Math.Sin(Math.PI / n);
            double dc = 2 * t * t;
            double ds = Math.Sqrt(dc * (2 - dc));
            t = 2 * dc;
            double c = sintbl[n4] = 1;
            double s = sintbl[0] = 0;
            for (int i = 1; i < n8; i++)
            {
                c -= dc;
                dc += t * c;
                s += ds;
                ds -= t * s;
                sintbl[i] = s;
                sintbl[n4 - i] = c;
            }
            if (n8 != 0)
                sintbl[n8] = Math.Sqrt(0.5);
            for (int i = 0; i < n4; i++)
                sintbl[n2 - i] = sintbl[i];
            for (int i = 0; i < n2 + n4; i++)
                sintbl[i + n2] = -sintbl[i];
        }

        private static void MakeBitrev(int n, int[] bitrev)
        {
            int n2 = n / 2, i = 0, j = 0;
            while (true)
            {
                bitrev[i] = j;
                if (++i >= n)
                    break;
                int k = n2;
                while (k <= j)
                {
                    j -= k;
                    k /= 2;
                }
                j += k;
            }
        }

        /// <summary><c>fft</c> of <c>fft.c</c> (Okumura): the in-place radix-2 transform; a negative <paramref name="n"/> is the inverse transform; the forward transform divides by n.</summary>
        private void Fft(int n, double[] xr, double[] xi)
        {
            bool inverse;
            if (n < 0)
            {
                n = -n;
                inverse = true;
            }
            else
            {
                inverse = false;
            }
            int n4 = n / 4;
            if (n != _fftLastN || n == 0)
            {
                _fftLastN = n;
                var st = new double[n + n4];
                var br = new int[n];
                MakeSinTbl(n, st);
                MakeBitrev(n, br);
                _fftSinTbl = st;
                _fftBitrev = br;
            }
            double[] sintbl = _fftSinTbl;
            int[] bitrev = _fftBitrev;
            for (int i = 0; i < n; i++)
            {
                int j = bitrev[i];
                if (i < j)
                {
                    double t = xr[i]; xr[i] = xr[j]; xr[j] = t;
                    t = xi[i]; xi[i] = xi[j]; xi[j] = t;
                }
            }
            int k2;
            for (int k = 1; k < n; k = k2)
            {
                int h = 0;
                k2 = k + k;
                int d = n / k2;
                for (int j = 0; j < k; j++)
                {
                    double c = sintbl[h + n4];
                    double s = inverse ? -sintbl[h] : sintbl[h];
                    for (int i = j; i < n; i += k2)
                    {
                        int ik = i + k;
                        double dR = s * xi[ik] + c * xr[ik];
                        double dI = c * xi[ik] - s * xr[ik];
                        xr[ik] = xr[i] - dR;
                        xr[i] += dR;
                        xi[ik] = xi[i] - dI;
                        xi[i] += dI;
                    }
                    h += d;
                }
            }
            if (!inverse)
            {
                for (int i = 0; i < n; i++)
                {
                    xr[i] /= n;
                    xi[i] /= n;
                }
            }
        }

        // ---- fftFunctions.c

        /// <summary><c>getKouho</c>: the <paramref name="nkouho"/> highest peaks of the correlation (lags).</summary>
        private static void GetKouho(int[] kouho, int nkouho, double[] soukan, int nlen2)
        {
            int nlen4 = nlen2 / 2;
            int ikouho = 0;
            for (int j = 0; j < nkouho; j++)
            {
                double max = -9999.9;
                for (int i = 0; i < nlen2; i++)
                {
                    double tmp;
                    if ((tmp = soukan[i]) > max)
                    {
                        ikouho = i;
                        max = tmp;
                    }
                }
                soukan[ikouho] = -9999.9;
                kouho[j] = ikouho - nlen4;
            }
        }

        private double[] _arStra = Array.Empty<double>();
        private int _arAlloclen;
        private double _arThreshold;

        /// <summary><c>alignableReagion</c>: the windows of the shifted groups (<paramref name="off1"/>, <paramref name="off2"/> characters skipped) whose average column score exceeds the threshold.</summary>
        private int AlignableRegion(int clus1, int clus2, string[] seq1, int off1, string[] seq2, int off2, double[] eff1, double[] eff2, SegList seg, int segBase)
        {
            int na = G.Nalphabets;
            var prf1 = new double[na];
            var prf2 = new double[na];
            var hat1 = new int[na + 1];
            var hat2 = new int[na + 1];
            int len1 = seq1[0].Length - off1, len2 = seq2[0].Length - off2;
            int len = Math.Min(len1, len2);
            int maxlen = Math.Max(len1, len2) + G.FftWinSize;
            if (_arAlloclen < maxlen)
            {
                if (_arAlloclen == 0)
                    _arThreshold = (double)G.FftThreshold / 100.0 * 600.0 * G.FftWinSize;
                _arStra = new double[maxlen];
                _arAlloclen = maxlen;
            }
            double[] stra = _arStra;
            double totaleff = 0.0;
            for (int i = 0; i < clus1; i++)
                for (int j = 0; j < clus2; j++)
                    totaleff += eff1[i] * eff2[j];
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < na; j++)
                {
                    prf1[j] = 0.0;
                    prf2[j] = 0.0;
                }
                for (int j = 0; j < clus1; j++)
                    prf1[G.AminoN[seq1[j][off1 + i]]] += eff1[j];
                for (int j = 0; j < clus2; j++)
                    prf2[G.AminoN[seq2[j][off2 + i]]] += eff2[j];
                int pre1 = na, pre2 = na;
                for (int j = 25; j >= 0; j--)
                {
                    if (prf1[j] != 0.0)
                    {
                        hat1[pre1] = j;
                        pre1 = j;
                    }
                    if (prf2[j] != 0.0)
                    {
                        hat2[pre2] = j;
                        pre2 = j;
                    }
                }
                hat1[pre1] = -1;
                hat2[pre2] = -1;
                stra[i] = 0.0;
                for (int k = hat1[na]; k != -1; k = hat1[k])
                    for (int j = hat2[na]; j != -1; j = hat2[j])
                        stra[i] += G.NDisFft[k][j] * prf1[k] * prf2[j];
                stra[i] /= totaleff;
            }

            int winSize = G.FftWinSize;
            int status = 0, starttmp = 0, length = 0, value = 0;
            double cumscore = 0.0, score = 0.0;
            int s = segBase;
            for (int j = 0; j < winSize; j++)
                score += stra[j];
            int ii;
            for (ii = 1; ii < len - winSize; ii++)
            {
                score = score - stra[ii - 1] + stra[ii + winSize - 1];
                if (score > _arThreshold)
                {
                    if (status == 0)
                    {
                        status = 1;
                        starttmp = ii;
                        length = 0;
                        cumscore = 0.0;
                    }
                    length++;
                    cumscore += score;
                }
                if ((score <= _arThreshold || length > SegmentSize) && status != 0)
                {
                    if (length > winSize)
                    {
                        var sg = seg[s];
                        sg.Start = starttmp;
                        sg.End = ii;
                        sg.Center = (sg.Start + sg.End + winSize) / 2;
                        sg.Score = cumscore;
                        value++;
                        s++;
                    }
                    length = 0;
                    cumscore = 0.0;
                    status = 0;
                    starttmp = ii;
                    if (value > MaxSeg - 3)
                        throw new InvalidOperationException("TOO MANY SEGMENTS!");
                }
            }
            if (status != 0 && length > winSize)
            {
                var sg = seg[s];
                sg.End = ii;
                sg.Start = starttmp;
                sg.Center = (starttmp + ii + winSize) / 2;
                sg.Score = cumscore;
                value++;
            }
            return value;
        }

        // ---- blockAlign2 (fftFunctions.c): statics of the program that survive from call to call
        private double _baMaxj, _baMaxi;
        private int _baPointj, _baPointi;
        private readonly int[] _baResult1 = new int[MaxSeg], _baResult2 = new int[MaxSeg], _baOcut1 = new int[MaxSeg], _baOcut2 = new int[MaxSeg];
        private int _baCrossscoresize;
        private int[][] _baTrack = Array.Empty<int[]>();
        private double[][] _baCrossscore = Array.Empty<double[]>();

        /// <summary><c>blockAlign2</c>: the highest-scoring chain of the anchors (<paramref name="ocrossscore"/>). <paramref name="ncut"/> receives the number of cuts that remain.</summary>
        private void BlockAlign2(int[] cut1, int[] cut2, double[][] ocrossscore, ref int ncut)
        {
            if (_baCrossscoresize < ncut + 2)
            {
                _baCrossscoresize = ncut + 2;
                _baTrack = new int[_baCrossscoresize][];
                _baCrossscore = new double[_baCrossscoresize][];
                for (int i = 0; i < _baCrossscoresize; i++)
                {
                    _baTrack[i] = new int[_baCrossscoresize];
                    _baCrossscore[i] = new double[_baCrossscoresize];
                }
            }
            int[][] track = _baTrack;
            double[][] crossscore = _baCrossscore;
            int[] ocut1 = _baOcut1, ocut2 = _baOcut2, result1 = _baResult1, result2 = _baResult2;
            for (int i = 0; i < ncut; i++)
                for (int j = 0; j < ncut; j++)
                    crossscore[i][j] = ocrossscore[i][j];
            for (int i = 0; i < ncut; i++)
            {
                ocut1[i] = cut1[i];
                ocut2[i] = cut2[i];
            }
            int penalty = G.Penalty;
            for (int i = 1; i < ncut; i++)
            {
                for (int j = 1; j < ncut; j++)
                {
                    _baPointi = 0;
                    _baMaxi = 0.0;
                    int klim = j - 2;
                    for (int k = 0; k < klim; k++)
                    {
                        // permit() always returns 0 in the program, so for k > 0 the candidate is skipped while j < ncut-1
                        if (k != 0 && k < ncut - 1 && j < ncut - 1)
                            continue;
                        if (crossscore[i - 1][k] > _baMaxj)
                        {
                            _baPointi = k;
                            _baMaxi = crossscore[i - 1][k];
                        }
                    }
                    _baPointj = 0;
                    _baMaxj = 0.0;
                    klim = i - 2;
                    for (int k = 0; k < klim; k++)
                    {
                        if (k != 0 && k < ncut - 1 && i < ncut - 1)
                            continue;
                        if (crossscore[k][j - 1] > _baMaxj)
                        {
                            _baPointj = k;
                            _baMaxj = crossscore[k][j - 1];
                        }
                    }
                    _baMaxi += penalty;
                    _baMaxj += penalty;
                    double maximum = crossscore[i - 1][j - 1];
                    track[i][j] = 0;
                    if (maximum < _baMaxi)
                    {
                        maximum = _baMaxi;
                        track[i][j] = j - _baPointi;
                    }
                    if (maximum < _baMaxj)
                    {
                        maximum = _baMaxj;
                        track[i][j] = _baPointj - i;
                    }
                    crossscore[i][j] += maximum;
                }
            }
            result1[MaxSeg - 1] = ncut - 1;
            result2[MaxSeg - 1] = ncut - 1;
            int ii;
            for (ii = MaxSeg - 1; ii >= 1; ii--)
            {
                int cur1 = result1[ii];
                int cur2 = result2[ii];
                if (cur1 == 0 || cur2 == 0)
                    break;
                int shift = track[cur1][cur2];
                if (shift == 0)
                {
                    result1[ii - 1] = cur1 - 1;
                    result2[ii - 1] = cur2 - 1;
                }
                else if (shift > 0)
                {
                    result1[ii - 1] = cur1 - 1;
                    result2[ii - 1] = cur2 - shift;
                }
                else
                {
                    result1[ii - 1] = cur1 + shift;
                    result2[ii - 1] = cur2 - 1;
                }
            }
            int count = 0;
            for (int j = ii; j < MaxSeg; j++)
            {
                if (ocrossscore[result1[j]][result2[j]] == 0.0)
                    continue;
                if ((result1[j] == result1[j - 1] || result2[j] == result2[j - 1]) && ocrossscore[result1[j]][result2[j]] > ocrossscore[result1[j - 1]][result2[j - 1]])
                    count--;
                cut1[count] = ocut1[result1[j]];
                cut2[count] = ocut2[result2[j]];
                count++;
            }
            ncut = count;
        }

        /// <summary><c>mymergesort</c> of <c>Falign.c</c>: the stable merge sort of the segments by their centre.</summary>
        private static void MyMergeSort(int first, int last, FftSegment[] seg)
        {
            if (first < last)
            {
                int middle = (first + last) / 2;
                MyMergeSort(first, middle, seg);
                MyMergeSort(middle + 1, last, seg);
                var work = new FftSegment[middle - first + 1];
                int p = 0;
                for (int i = first; i <= middle; i++)
                    work[p++] = seg[i];
                int ii = middle + 1, j = 0, k = first;
                while (ii <= last && j < p)
                {
                    if (work[j].Center <= seg[ii].Center)
                        seg[k++] = work[j++];
                    else
                        seg[k++] = seg[ii++];
                }
                while (j < p)
                    seg[k++] = work[j++];
            }
        }

        // ---- Falign (Falign.c)

        private readonly SegList _fSegment = new(), _fSegment1 = new(), _fSegment2 = new();
        private int _fCrossscoresize;
        private double[][] _fCrossscore = Array.Empty<double[]>();
        private readonly int[] _fCut1 = new int[MaxSeg], _fCut2 = new int[MaxSeg];
        private readonly int[] _fKouho = new int[NKouhoLong];
        private const int NKouhoLong = 500;

        /// <summary><c>*fftlog</c> after the last <see cref="Falign"/>: the number of segments of the alignment minus 2 (-1 for a single stretch).</summary>
        internal int FftlogResult;

        /// <summary>
        /// <c>Falign</c> of <c>Falign.c</c>: the alignment of two groups of aligned sequences by the fast Fourier transform (Katoh et al. 2002): the correlation of the polarity / volume (amino acids) or
        /// of the letter 'a' (nucleotides) of the two groups gives the 20 best shifts, the windows with a high score are the anchors, the chain of anchors cuts the problem into segments, which are
        /// aligned by <c>G__align11</c> (two single sequences) or <see cref="AAlign"/>. The rows of the groups are replaced by the aligned rows.
        /// With <paramref name="memsave"/> it is <c>Falign_udpari_long</c> (alg M): 500 shifts, the segments are aligned by <see cref="MsAlign"/>, the progress text "DP n / m" is written.
        /// </summary>
        internal double Falign(double[][] dyn, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int clus1, int clus2, int alloclen, bool memsave = false)
        {
            int len1 = seq1[0].Length, len2 = seq2[0].Length;
            int nlentmp = Math.Max(len1, len2);
            int nlen = 1;
            while (nlentmp >= nlen)
                nlen <<= 1;
            int nlen2 = nlen / 2;
            var tmpseq1 = new string[clus1];
            var tmpseq2 = new string[clus2];
            for (int j = 0; j < clus1; j++)
                tmpseq1[j] = seq1[j];
            for (int j = 0; j < clus2; j++)
                tmpseq2[j] = seq2[j];

            if ((G.CodonPos || G.CodonScore) && !_codonReady)
                InitCodon(seq1[0]);
            bool kobetsu = G.Kobetsubunkatsu != 0;
            int[] kouho = _fKouho;
            if (!kobetsu)
            {
                var vec1R = new double[nlen + 1];
                var vec1I = new double[nlen + 1];
                var vec2R = new double[nlen + 1];
                var vec2I = new double[nlen + 1];
                bool scoreVec = G.FftScore != 0 && G.Scoremtx != -1;
                // Falign_udpari_long (alg M) codes a nucleotide as a complex number (seq_vec_4: a = +1, t = -1, g = +i, c = -i); Falign uses the indicator of a (seq_vec_3, one vector)
                bool complexNucleotides = memsave && G.Scoremtx == -1;
                for (int i = 0; i < clus1; i++)
                    SeqVec(vec1R, vec1I, eff1[i], tmpseq1[i], scoreVec, complexNucleotides);
                for (int i = 0; i < clus2; i++)
                    SeqVec(vec2R, vec2I, eff2[i], tmpseq2[i], scoreVec, complexNucleotides);
                Fft(nlen, vec2R, vec2I);
                Fft(nlen, vec1R, vec1I);
                var nR = new double[nlen];
                var nI = new double[nlen];
                for (int l = 0; l < nlen; l++)
                {
                    double x = vec1R[l], xi = vec1I[l], y = vec2R[l], yi = vec2I[l];
                    double vr = x * y + xi * yi;
                    double vi = -x * yi + xi * y;
                    nR[l] = 0.0 + vr;
                    nI[l] = 0.0 + vi;
                }
                Fft(-nlen, nR, nI);
                var soukan = new double[nlen + 1];
                for (int m = 0; m <= nlen2; m++)
                    soukan[m] = nR[nlen2 - m];
                for (int m = nlen2 + 1; m < nlen; m++)
                    soukan[m] = nR[nlen + nlen2 - m];
                GetKouho(kouho, memsave ? NKouhoLong : NKouho, soukan, nlen);

            }

            int count = 0;
            var tmpptr1Off = 0;
            var tmpptr2Off = 0;
            int maxk = memsave ? NKouhoLong : NKouho;
            if (kobetsu)
            {
                maxk = 1;
                kouho[0] = 0;
            }
            for (int k = 0; k < maxk; k++)
            {
                int lag = kouho[k];
                if (lag <= -len1 || len2 <= lag)
                    continue;
                if (lag > 0)
                {
                    tmpptr1Off = 0;
                    tmpptr2Off = lag;
                }
                else
                {
                    tmpptr1Off = -lag;
                    tmpptr2Off = 0;
                }
                int tmpint = AlignableRegion(clus1, clus2, seq1, tmpptr1Off, seq2, tmpptr2Off, eff1, eff2, _fSegment, count);
                if (count + tmpint > MaxSeg - 3)
                    throw new InvalidOperationException("TOO MANY SEGMENTS.");
                if (tmpint == 0 && !memsave)
                    break;
                while (tmpint-- > 0)
                {
                    var sg = _fSegment[count];
                    if (lag > 0)
                    {
                        _fSegment1[count].Start = sg.Start;
                        _fSegment1[count].End = sg.End;
                        _fSegment1[count].Center = sg.Center;
                        _fSegment1[count].Score = sg.Score;
                        _fSegment2[count].Start = sg.Start + lag;
                        _fSegment2[count].End = sg.End + lag;
                        _fSegment2[count].Center = sg.Center + lag;
                        _fSegment2[count].Score = sg.Score;
                    }
                    else
                    {
                        _fSegment1[count].Start = sg.Start - lag;
                        _fSegment1[count].End = sg.End - lag;
                        _fSegment1[count].Center = sg.Center - lag;
                        _fSegment1[count].Score = sg.Score;
                        _fSegment2[count].Start = sg.Start;
                        _fSegment2[count].End = sg.End;
                        _fSegment2[count].Center = sg.Center;
                        _fSegment2[count].Score = sg.Score;
                    }
                    _fSegment1[count].Pair = _fSegment2[count];
                    _fSegment2[count].Pair = _fSegment1[count];
                    count++;
                }
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
            if (kobetsu)
            {
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
            }
            else if (count >= 5000)
            {
                // Falign_udpari_long: many anchors, the upper-level DP is skipped
                Report("\nMany anchors were found. The upper-level DP is skipped.\n\n");
                cut1[0] = 0;
                cut2[0] = 0;
                int count0 = 0;
                for (int i = 0; i < count; i++)
                {
                    if (sorted1[i].Center > cut1[count0] && sorted1[i].Pair!.Center > cut2[count0])
                    {
                        count0++;
                        cut1[count0] = sorted1[i].Center;
                        cut2[count0] = sorted1[i].Pair!.Center;
                    }
                    else if (i != 0 && sorted1[i].Score > sorted1[i - 1].Score
                        && sorted1[i].Center > cut1[count0 - 1] && sorted1[i].Pair!.Center > cut2[count0 - 1])
                    {
                        cut1[count0] = sorted1[i].Center;
                        cut2[count0] = sorted1[i].Pair!.Center;
                    }
                }
                cut1[count0 + 1] = len1;
                cut2[count0 + 1] = len2;
                count = count0 + 2;
            }
            else
            {
                if (_fCrossscoresize < count + 2)
                {
                    _fCrossscoresize = count + 2;
                    _fCrossscore = new double[_fCrossscoresize][];
                    for (int i = 0; i < _fCrossscoresize; i++)
                        _fCrossscore[i] = new double[_fCrossscoresize];
                }
                double[][] crossscore = _fCrossscore;
                for (int i = 0; i < count + 2; i++)
                    for (int j = 0; j < count + 2; j++)
                        crossscore[i][j] = 0.0;
                for (int i = 0; i < count; i++)
                {
                    crossscore[_fSegment1[i].Number + 1][_fSegment1[i].Pair!.Number + 1] = _fSegment1[i].Score;
                    cut1[i + 1] = sorted1[i].Center;
                    cut2[i + 1] = sorted2[i].Center;
                }
                crossscore[0][0] = 10000000.0;
                cut1[0] = 0;
                cut2[0] = 0;
                crossscore[count + 1][count + 1] = 10000000.0;
                cut1[count + 1] = len1;
                cut2[count + 1] = len2;
                count += 2;
                BlockAlign2(cut1, cut2, crossscore, ref count);

            }

            var result1 = new System.Text.StringBuilder[clus1];
            var result2 = new System.Text.StringBuilder[clus2];
            for (int j = 0; j < clus1; j++)
                result1[j] = new System.Text.StringBuilder();
            for (int j = 0; j < clus2; j++)
                result2[j] = new System.Text.StringBuilder();
            int totallen = 0;
            double totalscore = 0.0;
            FftlogResult = -1; // *fftlog = -1 and one more per segment: tbfast keeps it per sequence and tries the FFT again only for groups whose last FFT alignment had segments
            for (int i = 0; i < count - 1; i++)
            {
                FftlogResult++;
                int headgp = i == 0 ? G.Outgap : 1;
                int tailgp = i == count - 2 ? G.Outgap : 1;
                var sgap1 = new char[clus1];
                var egap1 = new char[clus1];
                var sgap2 = new char[clus2];
                var egap2 = new char[clus2];
                if (cut1[i] != 0)
                    for (int j = 0; j < clus1; j++)
                        sgap1[j] = seq1[j][cut1[i] - 1];
                else
                    for (int j = 0; j < clus1; j++)
                        sgap1[j] = 'o';
                if (cut2[i] != 0)
                    for (int j = 0; j < clus2; j++)
                        sgap2[j] = seq2[j][cut2[i] - 1];
                else
                    for (int j = 0; j < clus2; j++)
                        sgap2[j] = 'o';
                if (cut1[i + 1] != len1)
                    for (int j = 0; j < clus1; j++)
                        egap1[j] = seq1[j][cut1[i + 1]];
                else
                    for (int j = 0; j < clus1; j++)
                        egap1[j] = 'o';
                if (cut2[i + 1] != len2)
                    for (int j = 0; j < clus2; j++)
                        egap2[j] = seq2[j][cut2[i + 1]];
                else
                    for (int j = 0; j < clus2; j++)
                        egap2[j] = 'o';
                var tmpres1 = new string[clus1];
                var tmpres2 = new string[clus2];
                for (int j = 0; j < clus1; j++)
                    tmpres1[j] = seq1[j].Substring(cut1[i], cut1[i + 1] - cut1[i]);
                for (int j = 0; j < clus2; j++)
                    tmpres2[j] = seq2[j].Substring(cut2[i], cut2[i + 1] - cut2[i]);
                if (memsave)
                    Report("DP " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(5, '0') + " / " + (count - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5, '0') + " \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                if (kobetsu && G.FftKeika != 0)
                {
                    CommonGapPick(tmpres1);
                    CommonGapPick(tmpres2);
                }
                if (memsave || G.Alg == 'M')
                {
                    totalscore += MsAlign(dyn, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, new string(sgap1), new string(sgap2), new string(egap1), new string(egap2), headgp, tailgp, null, -1, -1, -1, 0.0, 0.0);
                }
                else if (G.Alg == 'A')
                {
                    if (clus1 == 1 && clus2 == 1)
                    {
                        if (G.CodonPos || G.CodonScore)
                            totalscore += G__align11(dyn, tmpres1, tmpres2, alloclen, headgp, tailgp, false, new MafftCodonDp(_codonMatrix, _codonGStart!, _codonGEnd!, cut1[i]));
                        else
                            totalscore += G__align11(dyn, tmpres1, tmpres2, alloclen, headgp, tailgp);
                    }
                    else
                    {
                        if (G.CodonPos)
                            throw new MafftProgramFailureException("\n\ncodonpos will be soon supported for a reference MSA. For now, use a single sequence as reference.\n\n\n");
                        totalscore += AAlign(dyn, G.Penalty, G.PenaltyEx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen,
                            new string(sgap1), new string(sgap2), new string(egap1), new string(egap2), headgp, tailgp, null, -1, -1, -1, 0.0, 0.0);
                    }
                }
                else if (G.Alg == 'd')
                {
                    if (clus1 == 1 && clus2 == 1)
                    {
                        totalscore += G__align11(dyn, tmpres1, tmpres2, alloclen, headgp, tailgp);
                    }
                    else
                    {
                        double dumdb = 0.0;
                        totalscore += DAlign(dyn, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, null, ref dumdb, headgp, tailgp);
                    }
                }
                else
                {
                    throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                }
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
            return totalscore;
        }

        /// <summary>The segment arrays of the program (<c>calloc( MAXSEG, sizeof( Segment ) )</c>), grown on demand.</summary>
        internal sealed class SegList
        {
            private FftSegment[] _a = new FftSegment[64];

            internal FftSegment this[int i]
            {
                get
                {
                    if (i >= _a.Length)
                        Array.Resize(ref _a, Math.Max(i + 1, _a.Length * 2));
                    return _a[i] ??= new FftSegment();
                }
            }
        }

        /// <summary><c>seq_vec_5</c> / <c>seq_vec_3</c>: adds the weighted physico-chemical value (polarity, volume) of every residue (amino acids) or the indicator of the letter 'a' (nucleotides) to the vectors.</summary>
        private void SeqVec(double[] resR, double[] resI, double incr, string seq, bool scoreVec, bool complexNucleotides)
        {
            if (complexNucleotides)
            {
                for (int i = 0; i < seq.Length; i++)
                {
                    char c = seq[i];
                    if (c == 'a')
                        resR[i] += incr;
                    else if (c == 't')
                        resR[i] -= incr;
                    else if (c == 'g')
                        resI[i] += incr;
                    else if (c == 'c')
                        resI[i] -= incr;
                }
            }
            else if (scoreVec)
            {
                for (int i = 0; i < seq.Length; i++)
                {
                    int n = G.AminoN[seq[i]];
                    if (n > 20)
                        continue;
                    resR[i] += incr * G.Polarity[n];
                    resI[i] += incr * G.Volume[n];
                }
            }
            else
            {
                for (int i = 0; i < seq.Length; i++)
                {
                    int n = G.AminoN[seq[i]];
                    if (n < 1 && n >= 0)
                        resR[i] += incr;
                }
            }
        }
    }

    #endregion
}
