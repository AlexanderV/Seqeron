using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region FASTA 3.4 (fasta34): the optimised score that splittbfast -S reads

    /// <summary>
    /// The part of FASTA 3.4 (W. R. Pearson, release 34.26.5: <c>comp_lib.c</c>, <c>dropnfa.c</c> (<c>init_work</c>, <c>do_fasta</c>, <c>savemax</c>, <c>spam</c>, <c>sconn</c>, <c>dmatch</c> / <c>FLOCAL_ALIGN</c>), <c>mshowbest.c</c>,
    /// <c>scaleswn.c</c> (<c>-z 3</c>: <c>proc_hist_a</c>, <c>find_za</c>, <c>zs_to_E</c>, <c>scale_scores</c>), <c>compacc.c</c> (<c>sortbeste</c>), <c>apam.c</c> / <c>upam.h</c> (BLOSUM50, +5/-4), <c>initfa.c</c>) that <c>splittbfast -S</c> of MAFFT 7.505
    /// (<c>mafft --fastaparttree</c>) uses: the command <c>fasta34 -z3 -m10 -p|-n -Q -H -b n -E n -d0 query library 2|6</c> and the <c>opt</c> that <c>ReadFasta34m10_scoreonly</c> reads from the list of best scores.
    /// </summary>
    internal static class MafftFasta34
    {
        private const int Na = 123, Nmap = 33, MaxSav = 10, BigNum = 1000000000;
        private const double EulerG = 0.57721566490153286060, PiSqrt6 = 1.28254983016186409554, ZsMax = 400.0;

        // ------------------------------------------------------------------ tables (upam.h, uascii.h)
        private static readonly int[] Bl50 =
        {
            5,-2,7,-1,-1,7,-2,-2,2,8,-1,-4,-2,-4,13,-1,1,0,0,-3,7,-1,0,0,2,-3,2,6,0,-3,0,-1,-3,-2,-3,8,-2,0,1,-1,
            -3,1,0,-2,10,-1,-4,-3,-4,-2,-3,-4,-4,-4,5,-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-1,-2,
            -2,-4,-2,0,-2,-3,-1,2,3,-2,7,-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,
            1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-3,-4,-5,-5,-1,-3,
            -3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,0,-3,-3,-4,-1,-3,-3,-4,-4,4,
            1,-3,1,-1,-3,-2,0,-3,-1,5,-2,-1,4,5,-3,0,1,-1,0,-4,-4,0,-3,-4,-2,0,0,-5,-3,-4,5,-1,0,0,1,-3,4,5,-2,0,
            -3,-3,1,-1,-4,-1,0,-1,-2,-2,-3,2,5,-1,-1,-1,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,0,-3,-1,-1,-1,-1,-1,-1,-1,-1,-1,
            -2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,0,-3,-1,-1,-1,-1,-1,7,
        };

        private static readonly int[] Nt54 =
        {
            5,-4,5,-4,-4,5,-4,-4,-4,5,-4,-4,-4,5,5,2,-1,2,-1,-1,2,-1,2,-1,2,2,-2,2,2,2,-1,-1,-1,-1,-1,2,2,-1,-1,2,
            2,1,1,1,2,-1,2,2,-1,-1,1,1,1,-1,2,-1,-1,2,2,2,1,1,-1,1,1,2,1,-2,1,1,1,1,-1,-1,1,-1,1,1,1,1,
            -2,1,1,-1,1,1,1,-1,-1,-1,1,1,1,1,-2,-2,1,-1,1,-1,1,-1,-1,-1,1,-2,1,1,1,1,-1,1,-1,-1,1,1,-1,-1,-1,1,
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
        };

        private const string AaOrder = "ARNDCQEGHILKMFPSTWYVBZX*";

        private sealed class Molecule
        {
            internal int Nsq;
            internal int[][] Pam2 = null!;
            internal int[] Hsq = null!;
            internal int[] Ascii = new int[128];
            internal int PamL, PamH, Ktup, GDel, GGap, BestScale, BestOff, BkFact, BkTup, BestMax, ScFact, PamFact, OptWid;
            internal double Lambda, K, H;
            internal bool Nucleotide;
            internal int[] CNt = Array.Empty<int>();
        }

        private static readonly Lazy<Molecule> Protein = new(() => BuildProtein());
        private static readonly Lazy<Molecule> Dna = new(() => BuildDna());

        private static int[][] LowerTriangle(int[] tab, int nsq)
        {
            var pam2 = new int[nsq + 1][];
            for (int i = 0; i <= nsq; i++)
            {
                pam2[i] = new int[nsq + 1];
                pam2[i][0] = -BigNum;
            }
            for (int j = 0; j <= nsq; j++)
                pam2[0][j] = -BigNum;
            int k = 0;
            for (int i = 1; i <= nsq; i++)
                for (int j = 1; j <= i; j++)
                    pam2[j][i] = pam2[i][j] = tab[k++];
            return pam2;
        }

        private static Molecule BuildProtein()
        {
            var m = new Molecule { Nsq = 24 };
            m.Pam2 = LowerTriangle(Bl50, 24);
            m.Hsq = new[] { Nmap, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 3, 7, Nmap, Nmap };
            for (int i = 0; i < 128; i++)
                m.Ascii[i] = Na;
            for (int i = 0; i < AaOrder.Length; i++)
            {
                m.Ascii[AaOrder[i]] = i + 1;
                if (AaOrder[i] is >= 'A' and <= 'Z')
                    m.Ascii[char.ToLowerInvariant(AaOrder[i])] = i + 1;
            }
            m.Ascii['U'] = m.Ascii['u'] = 23;
            m.PamL = 1;
            m.PamH = -1;
            for (int i = 1; i <= 24; i++)
                for (int j = 1; j <= i; j++)
                {
                    m.PamL = Math.Min(m.PamL, m.Pam2[i][j]);
                    m.PamH = Math.Max(m.PamH, m.Pam2[i][j]);
                }
            m.Ktup = 2;
            m.GDel = -10;
            m.GGap = -2;
            m.BestScale = 300; m.BestOff = 36; m.BkFact = 6; m.BkTup = 2; m.BestMax = 50; m.ScFact = 3; m.PamFact = 1; m.OptWid = 16;
            m.Lambda = 0.158; m.K = 0.019; m.H = 0.10; // bl50_p {12, 2}
            return m;
        }

        private static Molecule BuildDna()
        {
            var m = new Molecule { Nsq = 17, Nucleotide = true };
            m.Pam2 = LowerTriangle(Nt54, 17);
            m.Hsq = new[] { Nmap, 0, 1, 2, 3, 3, 0, 1, 0, 0, 1, 2, 0, 0, 0, 1, Nmap, Nmap };
            const string nt = "ACGTURYMWSKDHVBNX";
            for (int i = 0; i < 128; i++)
                m.Ascii[i] = Na;
            // nascii of uascii.h: the letters of the nucleotide alphabet; P and Q are mapped like R and Y, as in the table of the program
            for (int i = 0; i < nt.Length; i++)
            {
                m.Ascii[nt[i]] = i + 1;
                m.Ascii[char.ToLowerInvariant(nt[i])] = i + 1;
            }
            m.Ascii['P'] = m.Ascii['p'] = 6;
            m.Ascii['Q'] = m.Ascii['q'] = 7;
            m.PamL = 1;
            m.PamH = -1;
            for (int i = 1; i <= 17; i++)
                for (int j = 1; j <= i; j++)
                {
                    m.PamL = Math.Min(m.PamL, m.Pam2[i][j]);
                    m.PamH = Math.Max(m.PamH, m.Pam2[i][j]);
                }
            m.Ktup = 6;
            m.GDel = -12;
            m.GGap = -4;
            m.BestScale = 80; m.BestOff = 45; m.BkFact = 5; m.BkTup = 6; m.BestMax = 80; m.ScFact = 1; m.PamFact = 0; m.OptWid = 16;
            m.Lambda = 0.192; m.K = 0.177; m.H = 0.36; // nt54_p {16, 4}
            m.CNt = new[] { 0, 4, 3, 2, 1, 1, 7, 6, 11, 9, 10, 8, 13, 12, 15, 14, 16, 16 };
            return m;
        }

        // ------------------------------------------------------------------ the scores of one query against a library

        /// <summary>One query frame (<c>init_work</c>): the hash table of the query, the diagonal arrays.</summary>
        private sealed class FStr
        {
            internal int[] Aa0 = null!; // 1 + n0 + 1: [0] = 0, [1..n0] the residues, [n0 + 1] = 0
            internal int N0, Kshft, Hmask, Ktup, Kt1, Noff;
            internal int[] Harr = null!, Link = null!, Pamh1 = null!, Pamh2 = null!;
            internal int OptCut, Cgap, Pgap;
            // diagonal arrays
            internal int[] DScore = Array.Empty<int>(), DStart = Array.Empty<int>(), DStop = Array.Empty<int>(), DDmax = Array.Empty<int>();
            internal readonly int[] VScore = new int[MaxSav], VDp = new int[MaxSav], VStart = new int[MaxSav], VStop = new int[MaxSav];
            internal int LowMax, LowScor;
            internal int[] Vptr = new int[MaxSav];
            internal int[] BssCc = Array.Empty<int>(), BssDd = Array.Empty<int>();
        }

        private static int Shscore(int[] aa0, int n0, int[][] pam2)
        {
            int sum = 0;
            for (int i = 0; i < n0; i++)
                sum += pam2[aa0[i + 1]][aa0[i + 1]];
            return sum;
        }

        private static FStr InitWork(Molecule m, int[] aa0p, int n0)
        {
            var f = new FStr { Aa0 = aa0p, N0 = n0, Ktup = m.Ktup };
            int ktup = m.Ktup;
            int btemp = 2 * m.BestOff / 3 + n0 / m.BestScale + m.BkFact * (m.BkTup - ktup);
            if (m.Nucleotide)
                btemp = btemp * m.PamH / 5;
            btemp = Math.Min(btemp, m.BestMax);
            if (btemp > 3 * n0)
                btemp = 3 * Shscore(aa0p, n0, m.Pam2) / 5;
            f.Cgap = btemp + m.BestOff / 3;
            f.OptCut = btemp;
            f.Pgap = m.GDel + 2 * m.GGap;
            int pamfact = m.PamFact;
            int fact = m.ScFact * ktup;
            int mhv = -1;
            for (int i0 = 1; i0 <= m.Nsq; i0++)
                if (m.Hsq[i0] < Nmap && m.Hsq[i0] > mhv)
                    mhv = m.Hsq[i0];
            f.Kshft = 0;
            for (; mhv > 0; mhv /= 2)
                f.Kshft++;
            int kt1 = ktup - 1;
            f.Kt1 = kt1;
            int hv = 1;
            for (int i0 = 0; i0 < ktup; i0++)
                hv <<= f.Kshft;
            int hmax = hv;
            f.Hmask = (hmax >> f.Kshft) - 1;
            f.Harr = new int[hmax];
            f.Pamh1 = new int[m.Nsq + 1];
            f.Pamh2 = new int[hmax];
            f.Link = new int[n0];
            Array.Fill(f.Harr, -1);
            Array.Fill(f.Link, -1);
            int[] hsq = m.Hsq;
            int phv = hv = 0;
            int lkt = kt1;
            int i;
            for (i = 0; i < Math.Min(lkt, n0); i++)
            {
                if (hsq[aa0p[i + 1]] >= Nmap)
                {
                    hv = phv = 0;
                    lkt = i + ktup;
                    continue;
                }
                hv = (hv << f.Kshft) + hsq[aa0p[i + 1]];
                phv += m.Pam2[aa0p[i + 1]][aa0p[i + 1]] * ktup;
            }
            while (i < n0)
            {
                if (hsq[aa0p[i + 1]] >= Nmap)
                {
                    hv = phv = 0;
                    lkt = i + kt1;
                    while ((i < lkt || hsq[aa0p[i + 1]] >= Nmap) && i < n0)
                    {
                        if (hsq[aa0p[i + 1]] >= Nmap)
                        {
                            hv = phv = 0;
                            lkt = i + ktup;
                        }
                        else
                        {
                            hv = (hv << f.Kshft) + hsq[aa0p[i + 1]];
                            phv += m.Pam2[aa0p[i + 1]][aa0p[i + 1]] * ktup;
                        }
                        i++;
                    }
                }
                if (i >= n0)
                    break;
                hv = ((hv & f.Hmask) << f.Kshft) + hsq[aa0p[i + 1]];
                f.Link[i] = f.Harr[hv];
                f.Harr[hv] = i;
                if (pamfact != 0)
                {
                    phv += m.Pam2[aa0p[i + 1]][aa0p[i + 1]] * ktup;
                    f.Pamh2[hv] = phv;
                    if (hsq[aa0p[i - kt1 + 1]] < Nmap)
                        phv -= m.Pam2[aa0p[i - kt1 + 1]][aa0p[i - kt1 + 1]] * ktup;
                }
                else
                {
                    f.Pamh2[hv] = fact * ktup;
                }
                i++;
            }
            if (pamfact != 0)
                for (int i0 = 1; i0 <= m.Nsq; i0++)
                    f.Pamh1[i0] = m.Pam2[i0][i0] * ktup;
            else
                for (int i0 = 1; i0 <= m.Nsq; i0++)
                    f.Pamh1[i0] = fact;
            f.Noff = n0 - 1;
            f.BssCc = new int[m.OptWid * 2 + 8];
            f.BssDd = new int[m.OptWid * 2 + 8];
            return f;
        }

        private static void SaveMax(FStr f, int dpos)
        {
            int vm = f.DDmax[dpos];
            int i;
            if (vm != -1 && f.VDp[vm] == dpos && f.VStart[vm] == f.DStart[dpos])
            {
                f.VStop[vm] = f.DStop[dpos];
                if ((i = f.DScore[dpos]) <= f.VScore[vm])
                    return;
                f.VScore[vm] = i;
                if (vm != f.LowMax)
                    return;
            }
            else
            {
                i = f.VScore[f.LowMax] = f.DScore[dpos];
                f.VDp[f.LowMax] = dpos;
                f.VStart[f.LowMax] = f.DStart[dpos];
                f.VStop[f.LowMax] = f.DStop[dpos];
                f.DDmax[dpos] = f.LowMax;
            }
            for (int v = 0; v < MaxSav; v++)
            {
                if (f.VScore[v] < i)
                {
                    i = f.VScore[v];
                    f.LowMax = v;
                }
            }
            f.LowScor = i;
        }

        private static int Spam(int[] aa0, int[] aa1, FStr f, int vm, int[][] pam2, int xdrop, int noff, bool doExtend)
        {
            int start = f.VStart[vm], stop = f.VStop[vm], dp = f.VDp[vm];
            int lpos = start;
            int p1 = lpos + 1;               // aa1[lpos] in the padded array
            int p0 = lpos - dp + noff + 1;   // aa0[lpos - dp + noff]
            int curStart = lpos, curStop = 0, curScore = 0, maxStart = 0, maxStop = 0, maxScore = 0, tot = 0;
            for (; lpos <= stop; lpos++)
            {
                tot += pam2[aa0[p0++]][aa1[p1++]];
                if (tot > curScore)
                {
                    curStop = lpos;
                    curScore = tot;
                }
                else if (tot < 0)
                {
                    if (curScore > maxScore)
                    {
                        maxStart = curStart;
                        maxStop = curStop;
                        maxScore = curScore;
                    }
                    tot = curScore = 0;
                    curStart = lpos + 1;
                    if (lpos >= stop)
                        break;
                }
            }
            if (curScore > maxScore)
            {
                maxStart = curStart;
                maxStop = curStop;
                maxScore = curScore;
            }
            if (doExtend && maxScore > xdrop)
            {
                if (maxStop == stop)
                {
                    tot = maxScore;
                    int dropThresh = maxScore - xdrop;
                    p1 = stop + 1;
                    lpos = stop;
                    p0 = lpos - dp + noff + 1;
                    while (tot > dropThresh)
                    {
                        ++lpos;
                        tot += pam2[aa0[++p0]][aa1[++p1]];
                        if (tot > maxScore)
                        {
                            maxStart = lpos;
                            maxScore = tot;
                            dropThresh = tot - xdrop;
                        }
                    }
                }
                if (maxStart == start)
                {
                    tot = maxScore;
                    int dropThresh = maxScore - xdrop;
                    lpos = start;
                    p1 = lpos + 1;
                    p0 = lpos - dp + noff + 1;
                    while (tot > dropThresh)
                    {
                        --lpos;
                        tot += pam2[aa0[--p0]][aa1[--p1]];
                        if (tot > maxScore)
                        {
                            maxStart = lpos;
                            maxScore = tot;
                            dropThresh = tot - xdrop;
                        }
                    }
                }
            }
            f.VStart[vm] = maxStart;
            f.VStop[vm] = maxStop;
            return maxScore;
        }

        private static void KpSort(int[] v, int n, int[] start)
        {
            for (int gap = n / 2; gap > 0; gap /= 2)
                for (int i = gap; i < n; i++)
                    for (int j = i - gap; j >= 0; j -= gap)
                    {
                        if (start[v[j]] <= start[v[j + gap]])
                            break;
                        (v[j], v[j + gap]) = (v[j + gap], v[j]);
                    }
        }

        private static int Sconn(FStr f, int[] v, int n, int cgap, int pgap, int noff)
        {
            KpSort(v, n, f.VStart);
            var sScore = new int[MaxSav];
            var sVp = new int[MaxSav];
            var sNext = new int[MaxSav];
            int start = -1;
            int si = 0;
            for (int i = 0; i < n; i++)
            {
                if (f.VScore[v[i]] < cgap)
                    continue;
                int lstart = f.VStart[v[i]];
                int tstart = lstart - f.VDp[v[i]] + noff;
                sVp[si] = v[i];
                sScore[si] = f.VScore[v[i]];
                sNext[si] = -1;
                for (int sl = start; sl != -1; sl = sNext[sl])
                {
                    int plstop = f.VStop[sVp[sl]];
                    int ptstop = plstop - f.VDp[sVp[sl]] + noff;
                    if (plstop < lstart && ptstop < tstart)
                    {
                        sScore[si] = sScore[sl] + f.VScore[v[i]] + pgap;
                        break;
                    }
                }
                if (start == -1)
                {
                    start = si;
                }
                else
                {
                    int so = -1;
                    for (int sj = start; sj != -1; sj = sNext[sj])
                    {
                        if (sScore[si] > sScore[sj])
                        {
                            sNext[si] = sj;
                            if (so != -1)
                                sNext[so] = si;
                            else
                                start = si;
                            break;
                        }
                        so = sj;
                    }
                }
                si++;
            }
            return start != -1 ? sScore[start] : 0;
        }

        /// <summary><c>FLOCAL_ALIGN</c> (the banded local alignment of <c>dmatch</c>): <paramref name="a"/> and <paramref name="b"/> are the padded arrays (a[i] = residue i of the 1-based sequence).</summary>
        private static int FlocalAlign(FStr f, int[] a, int[] b, int mm, int nn, int low, int up, int[][] w, int g, int h)
        {
            var cc = f.BssCc;
            var dd = f.BssDd;
            const int minInt = -9999999;
            int m = g + h;
            low = Math.Max(-mm, low);
            up = Math.Min(nn, up);
            if (nn <= 0)
                return 0;
            if (mm <= 0)
                return 0;
            int band = up - low + 1;
            if (band < 1)
                return 0;
            int leftd;
            if (low > 0)
                leftd = 1;
            else if (up < 0)
                leftd = band;
            else
                leftd = 1 - low;
            int rightd = band;
            int si = Math.Max(0, -up);
            int ei = Math.Min(mm, nn - low);
            // the program indexes bss[-1 .. band + 1]; here bss[x] is stored at x + 1
            cc[leftd + 1] = 0;
            for (int j = leftd + 1; j <= rightd; j++)
            {
                cc[j + 1] = 0;
                dd[j + 1] = -g;
            }
            cc[rightd + 2] = minInt;
            dd[rightd + 2] = minInt;
            int bestScore = 0;
            cc[leftd] = minInt;
            dd[leftd + 1] = -g;
            int c = 0, d, e;
            for (int i = si + 1; i <= ei; i++)
            {
                if (i > nn - up)
                    rightd--;
                if (leftd > 1)
                    leftd--;
                int[] wa = w[a[i]];
                if ((c = cc[leftd + 2] - m) > (d = dd[leftd + 2] - h))
                    d = c;
                int ib;
                if ((ib = leftd + low - 1 + i) > 0)
                    c = cc[leftd + 1] + wa[b[ib]];
                if (d > c)
                    c = d;
                if (c < 0)
                    c = 0;
                e = c - g;
                dd[leftd + 1] = d;
                cc[leftd + 1] = c;
                if (c > bestScore)
                    bestScore = c;
                for (int curd = leftd + 1; curd <= rightd; curd++)
                {
                    if ((c = c - m) > (e = e - h))
                        e = c;
                    if ((c = cc[curd + 2] - m) > (d = dd[curd + 2] - h))
                        d = c;
                    c = cc[curd + 1] + wa[b[curd + low - 1 + i]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    if (c < 0)
                        c = 0;
                    cc[curd + 1] = c;
                    dd[curd + 1] = d;
                    if (c > bestScore)
                        bestScore = c;
                }
            }
            return bestScore;
        }

        /// <summary><c>do_fasta</c>: the three scores (initn, init1, opt) of the query against one library sequence.</summary>
        private static int DoFasta(Molecule m, FStr f, int[] aa1, int n1)
        {
            int ktup = m.Ktup, kt1 = ktup - 1;
            int[] hsq = m.Hsq;
            int n0 = f.N0;
            int[] aa0 = f.Aa0;
            int xdrop = -m.PamL;
            bool doExtend = !m.Nucleotide || n0 < 50 || n1 < 50;
            if (n1 < ktup)
                return 0;
            int nd = n0 + n1;
            if (f.DScore.Length < nd)
            {
                f.DScore = new int[nd];
                f.DStart = new int[nd];
                f.DStop = new int[nd];
                f.DDmax = new int[nd];
            }
            for (int d = 0; d < nd; d++)
            {
                f.DStop[d] = -1;
                f.DDmax[d] = -1;
                f.DScore[d] = 0;
            }
            for (int v = 0; v < MaxSav; v++)
                f.VScore[v] = 0;
            f.LowMax = 0;
            f.LowScor = 0;
            int lhval = 0;
            int lkt = kt1;
            int lpos;
            for (lpos = 0; (lpos < lkt || hsq[aa1[lpos + 1]] >= Nmap) && lpos < n1; lpos++)
            {
                if (hsq[aa1[lpos + 1]] >= Nmap)
                {
                    lhval = 0;
                    lkt = lpos + ktup;
                    continue;
                }
                lhval = ((lhval & f.Hmask) << f.Kshft) + hsq[aa1[lpos + 1]];
            }
            int noff = f.Noff;
            int diagp = noff + lkt;
            while (lpos < n1)
            {
                if (hsq[aa1[lpos + 1]] >= Nmap)
                {
                    lpos++;
                    diagp++;
                    while (lpos < n1 && hsq[aa1[lpos + 1]] >= Nmap)
                    {
                        lpos++;
                        diagp++;
                    }
                    if (lpos >= n1)
                        break;
                    lhval = 0;
                }
                lhval = ((lhval & f.Hmask) << f.Kshft) + hsq[aa1[lpos + 1]];
                for (int tpos = f.Harr[lhval]; tpos >= 0; tpos = f.Link[tpos])
                {
                    int dptr = diagp - tpos;
                    int tscor = f.DStop[dptr];
                    if (tscor >= 0)
                    {
                        tscor += ktup;
                        tscor -= lpos;
                        if (tscor <= 0)
                        {
                            int scor = f.DScore[dptr];
                            int kfact = f.Pamh2[lhval];
                            tscor += kfact;
                            if (tscor < 0 && f.LowScor < scor)
                                SaveMax(f, dptr);
                            tscor += scor;
                            if (tscor >= kfact)
                            {
                                f.DScore[dptr] = tscor;
                                f.DStop[dptr] = lpos;
                            }
                            else
                            {
                                f.DScore[dptr] = kfact;
                                f.DStop[dptr] = lpos;
                                f.DStart[dptr] = lpos - kt1;
                            }
                        }
                        else
                        {
                            f.DScore[dptr] += f.Pamh1[aa0[tpos + 1]];
                            f.DStop[dptr] = lpos;
                        }
                    }
                    else
                    {
                        f.DScore[dptr] = f.Pamh2[lhval];
                        f.DStop[dptr] = lpos;
                        f.DStart[dptr] = lpos - kt1;
                    }
                }
                lpos++;
                diagp++;
            }
            for (int d = 0; d < nd; d++)
            {
                if (f.DScore[d] > f.LowScor)
                    SaveMax(f, d);
                f.DStop[d] = -1;
                f.DDmax[d] = -1;
                f.DScore[d] = 0;
            }
            int nsave = 0;
            for (int v = 0; v < MaxSav; v++)
            {
                if (f.VScore[v] > 0)
                {
                    f.VScore[v] = Spam(aa0, aa1, f, v, m.Pam2, xdrop, noff, doExtend);
                    f.Vptr[nsave++] = v;
                }
            }
            if (nsave <= 0)
                return 0;
            int scorN = Sconn(f, f.Vptr, nsave, f.Cgap, f.Pgap, noff);
            int best = f.Vptr[0];
            for (int ib = 1; ib < nsave; ib++)
                if (f.VScore[f.Vptr[ib]] > f.VScore[best])
                    best = f.Vptr[ib];
            int init1 = f.VScore[best];
            int initn = Math.Max(scorN, init1);
            int opt = initn;
            if (initn > f.OptCut)
            {
                int hoff = noff - f.VDp[best];
                int window = Math.Min(n1, m.OptWid);
                int low = -window / 2 - hoff;
                int up = low + window;
                opt = FlocalAlign(f, aa0, aa1, n0, n1, low, up, m.Pam2, -m.GDel, -m.GGap);
            }
            return opt;
        }

        private static int[] Encode(string s, int[] ascii)
        {
            var list = new List<int> { 0 };
            foreach (char ch in s)
            {
                int c = ch < 128 ? ascii[ch] : Na;
                if (c < Na)
                    list.Add(c);
            }
            list.Add(0);
            return list.ToArray();
        }

        private sealed class Entry
        {
            internal int Seq, Frame, N1, Opt;
            internal double Z, E;
        }

        private static double FindZa(Molecule m, double a_n0, double a_n0f, int opt, int n1)
        {
            double a_n1 = n1;
            double a_n1f = Math.Log(a_n1) / m.H;
            double mp = a_n0 - a_n0f - a_n1f;
            double np = a_n1 - a_n0f - a_n1f;
            if (np < 1.0)
                np = 1.0;
            if (mp < 1.0)
                mp = 1.0;
            double z = m.Lambda * opt - Math.Log(m.K * np * mp);
            z = -z + EulerG;
            z /= -PiSqrt6;
            return 50.0 + z * 10.0;
        }

        private static double ZsToE(double zs, int n1, bool dna, long entries, long dbLength, long dbEntries)
        {
            double z = (zs - 50.0) / 10.0;
            if (z > ZsMax)
                return 0.0;
            if (entries < 1)
                entries = dbEntries;
            double k;
            if (dna)
                k = (double)dbLength / n1;
            else
                k = entries;
            if (k < 1.0)
                k = 1.0;
            z *= PiSqrt6;
            z += EulerG;
            double e = Math.Exp(-z);
            return k * (e > .01 ? 1.0 - Math.Exp(-e) : e);
        }

        private static void SortBestZ(List<Entry> b, int n)
        {
            int[] incs = { 1391376, 463792, 198768, 86961, 33936, 13776, 4592, 1968, 861, 336, 112, 48, 21, 7, 3, 1 };
            foreach (int gap in incs)
            {
                for (int i = gap; i < n; i++)
                {
                    var tmp = b[i];
                    int j = i;
                    double v = b[i].Z;
                    while (j >= gap && b[j - gap].Z < v)
                    {
                        b[j] = b[j - gap];
                        j -= gap;
                    }
                    b[j] = tmp;
                }
            }
        }

        private static void SortBestE(List<Entry> b, int n)
        {
            int[] incs = { 1391376, 463792, 198768, 86961, 33936, 13776, 4592, 1968, 861, 336, 112, 48, 21, 7, 3, 1 };
            foreach (int gap in incs)
            {
                for (int i = gap; i < n; i++)
                {
                    int j = i;
                    var tmp = b[i];
                    double v = tmp.E;
                    while (j >= gap && b[j - gap].E > v)
                    {
                        b[j] = b[j - gap];
                        j -= gap;
                    }
                    b[j] = tmp;
                }
            }
            int jj = 0;
            while (jj < n && b[jj].E <= 2.0 * 2.2250738585072014e-308)
                jj++;
            if (jj > 1)
                SortBestZ(b, jj);
        }

        /// <summary>
        /// <c>fasta34 -z3 -m10 -p|-n -Q -H -b n -E n -d0 query library 2|6</c> read by <c>ReadFasta34m10_scoreonly[_nuc]</c>: the <c>opt</c> score of every library sequence that is on the list of the best scores (0 for those that are not).
        /// The library is the sequences of <paramref name="library"/> in this order.
        /// </summary>
        internal static double[] Search(string query, IReadOnlyList<string> library, bool nucleotide, int mshow = -1)
        {
            var m = nucleotide ? Dna.Value : Protein.Value;
            int nin = library.Count;
            if (mshow < 0)
                mshow = nin;
            double eCut = nin;
            var res = new double[nin];
            int[] q0 = Encode(query, m.Ascii);
            int n0 = q0.Length - 2;
            if (n0 <= 0)
                throw new InvalidOperationException("Query sequence length <= 0");
            var frames = new List<FStr> { InitWork(m, q0, n0) };
            if (nucleotide)
            {
                var rc = new int[n0 + 2];
                for (int i = 0; i < n0; i++)
                    rc[i + 1] = m.CNt[q0[n0 - i]];
                frames.Add(InitWork(m, rc, n0));
            }
            var entries = new List<Entry>();
            long dbLength = 0;
            int dbEntries = 0;
            for (int s = 0; s < nin; s++)
            {
                int[] aa1 = Encode(library[s], m.Ascii);
                int n1 = aa1.Length - 2;
                if (n1 <= 1)
                    continue;
                dbEntries++;
                dbLength += n1;
                for (int itt = 0; itt < frames.Count; itt++)
                {
                    int opt = DoFasta(m, frames[itt], aa1, n1);
                    entries.Add(new Entry { Seq = s, Frame = itt, N1 = n1, Opt = opt });
                }
            }
            if (entries.Count == 0)
                return res;
            // process_hist (-z 3: Altschul / Gish parameters) and scale_scores
            double a_n0 = n0;
            double a_n0f = Math.Log(m.K * a_n0) / m.H;
            foreach (var e in entries)
            {
                e.Z = FindZa(m, a_n0, a_n0f, e.Opt, e.N1);
                e.E = ZsToE(e.Z, e.N1, nucleotide, dbEntries, dbLength, dbEntries);
            }
            SortBestE(entries, entries.Count);
            // showbest (quiet, -b n -E n): the list is printed in batches of 20, 30, 40 ... entries until -b entries are shown
            int nbest = entries.Count;
            if (!(entries[0].E < eCut))
                return res;
            int nshow = Math.Min(mshow, 20);
            int istart = 0;
            var seen = new bool[nin];
            while (true)
            {
                int istop = Math.Min(nbest, nshow);
                int ib;
                Entry? last = null;
                for (ib = istart; ib < istop; ib++)
                {
                    var bbp = entries[ib];
                    last = bbp;
                    if (bbp.E > eCut)
                        return res;
                    if (bbp.E < 0.0)
                        continue;
                    if ((!nucleotide || bbp.Frame == 0) && !seen[bbp.Seq])
                    {
                        res[bbp.Seq] += bbp.Opt;
                        seen[bbp.Seq] = true;
                    }
                }
                if (ib < nbest && last != null && last.E < eCut)
                {
                    if (istop >= mshow)
                        break;
                    istart = istop;
                    nshow += 10;
                    continue;
                }
                break;
            }
            return res;
        }

        /// <summary><c>getfastascoremtx</c>: the diagonal of the scoring matrix of FASTA (the score of six equal residues, divided by six) for protein; the fixed scores of the nucleotide alphabet (no call of <c>fasta34</c>).</summary>
        internal static int[,] ScoreMatrix(bool nucleotide, Action<string> report)
        {
            var t = new int[128, 128];
            if (nucleotide)
            {
                t['a', 'a'] = 5;
                t['g', 'g'] = 5;
                t['c', 'c'] = 5;
                t['t', 't'] = 5;
                t['n', 'n'] = -1;
                return t;
            }
            report("done.\n");
            string tmpseq = string.Concat("ACDEFGHIKLMNPQRSTVWY".Select(c => new string(c, 6) + "XXXXXX"));
            const string amino = "ARNDCQEGHILKMFPSTWYV";
            for (int i = 0; i < 20; i++)
            {
                char aa = amino[i];
                double resvec = Search(new string(aa, 6), new[] { tmpseq }, false, 500000)[0];
                report(((char)('A' + i)).ToString() + ": " + (resvec / 6).ToString("F6", CultureInfo.InvariantCulture) + "\n");
                if (((int)resvec % 6) > 0.0)
                    throw new InvalidOperationException("Error in blast, *resvec=" + resvec.ToString("F6", CultureInfo.InvariantCulture));
                t[aa, aa] = (int)(resvec / 6);
            }
            t['X', 'X'] = -1;
            return t;
        }
    }

    #endregion
}
