#pragma warning disable S3776, S107, S1541, S138, S1854, S1121, S127, S3358, S1066, S3241
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region FASTA 3.4 (fasta34): the alignments that dndfast7 reads (fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup)

    internal static partial class MafftFasta34
    {
        // ------------------------------------------------------------------ C formatting helpers

        /// <summary>C <c>%.{prec}g</c> of a finite double (no <c>#</c> flag).</summary>
        internal static string FormatG(double v, int prec)
        {
            if (v == 0.0)
                return "0";
            string e = v.ToString("E" + (prec - 1).ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
            int x = int.Parse(e[(e.IndexOf('E') + 1)..], CultureInfo.InvariantCulture);
            if (x < -4 || x >= prec)
            {
                string mant = e[..e.IndexOf('E')];
                if (mant.Contains('.'))
                    mant = mant.TrimEnd('0').TrimEnd('.');
                return mant + "e" + (x < 0 ? "-" : "+") + Math.Abs(x).ToString("D2", CultureInfo.InvariantCulture);
            }

            string fixedText = v.ToString("F" + (prec - 1 - x).ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
            if (fixedText.Contains('.'))
                fixedText = fixedText.TrimEnd('0').TrimEnd('.');
            return fixedText;
        }

        private static string F(double v, int decimals) => v.ToString("F" + decimals.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);

        private static string I(long v) => v.ToString(CultureInfo.InvariantCulture);

        private static string Iw(long v, int width) => v.ToString(CultureInfo.InvariantCulture).PadLeft(width);

        // ------------------------------------------------------------------ the alignment arrays (dropnfa.c)

        /// <summary><c>a_res_str</c>: the span of the alignment and its edit script (<c>nres</c> columns).</summary>
        private sealed class AlnRes
        {
            internal int Min0, Max0, Min1, Max1, Nres;
            internal int[] Res = Array.Empty<int>();
        }

        /// <summary>The script of <c>nw_align</c> / <c>bg_align</c>: <c>*sapp</c> (index) and <c>*last</c>; <c>DEL</c>, <c>INS</c> and <c>REP</c> of the program.</summary>
        private sealed class Script
        {
            internal int[] S;
            internal int Pos, Last;

            internal Script(int[] s) => S = s;

            private void Need(int n)
            {
                if (n >= S.Length)
                    Array.Resize(ref S, Math.Max(S.Length * 2, n + 16));
            }

            internal void Del(int k)
            {
                if (Last < 0)
                {
                    S[Pos - 1] -= k;
                    Last = S[Pos - 1];
                }
                else
                {
                    Need(Pos + 1);
                    S[Pos] = -k;
                    Last = S[Pos];
                    Pos++;
                }
            }

            internal void Ins(int k)
            {
                if (Last > 0)
                {
                    S[Pos - 1] += k;
                    Last = S[Pos - 1];
                }
                else
                {
                    Need(Pos + 1);
                    S[Pos] = k;
                    Last = S[Pos];
                    Pos++;
                }
            }

            internal void Rep()
            {
                Need(Pos + 1);
                S[Pos] = 0;
                Last = 0;
                Pos++;
            }
        }

        private static int Gap(int k, int g, int h) => k <= 0 ? 0 : g + h * k;

        /// <summary><c>nw_align</c>: the optimum alignment of <c>A[1..M]</c> and <c>B[1..N]</c> in linear space (Myers and Miller); <c>A[i]</c> is <c>a[aoff + i]</c>.</summary>
        private static int NwAlign(int[] a, int aoff, int[] b, int boff, int mm, int nn, int tb, int te, int[][] w, int g, int h, int[] fH, int[] fE, int[] rH, int[] rE, Script sc)
        {
            int m = g + h;
            if (nn <= 0)
            {
                if (mm > 0)
                    sc.Del(mm);
                return -Gap(mm, g, h);
            }

            if (mm <= 1)
            {
                if (mm <= 0)
                {
                    sc.Ins(nn);
                    return -Gap(nn, g, h);
                }

                if (tb < te)
                    tb = te;
                int midc0 = (tb - h) - Gap(nn, g, h);
                int midj0 = 0;
                int[] wa0 = w[a[aoff + 1]];
                for (int j = 1; j <= nn; j++)
                {
                    int c0 = -Gap(j - 1, g, h) + wa0[b[boff + j]] - Gap(nn - j, g, h);
                    if (c0 > midc0)
                    {
                        midc0 = c0;
                        midj0 = j;
                    }
                }

                if (midj0 == 0)
                {
                    sc.Del(1);
                    sc.Ins(nn);
                }
                else
                {
                    if (midj0 > 1)
                        sc.Ins(midj0 - 1);
                    sc.Rep();
                    if (midj0 < nn)
                        sc.Ins(nn - midj0);
                }

                return midc0;
            }

            int midi = mm / 2;
            fH[0] = 0;
            int t = -g;
            for (int j = 1; j <= nn; j++)
            {
                fH[j] = t = t - h;
                fE[j] = t - g;
            }

            t = tb;
            int c, e, d, s;
            for (int i = 1; i <= midi; i++)
            {
                s = fH[0];
                fH[0] = c = t = t - h;
                e = t - g;
                int[] wa = w[a[aoff + i]];
                for (int j = 1; j <= nn; j++)
                {
                    if ((c = c - m) > (e = e - h))
                        e = c;
                    if ((c = fH[j] - m) > (d = fE[j] - h))
                        d = c;
                    c = s + wa[b[boff + j]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    s = fH[j];
                    fH[j] = c;
                    fE[j] = d;
                }
            }

            fE[0] = fH[0];
            rH[nn] = 0;
            t = -g;
            for (int j = nn - 1; j >= 0; j--)
            {
                rH[j] = t = t - h;
                rE[j] = t - g;
            }

            t = te;
            for (int i = mm - 1; i >= midi; i--)
            {
                s = rH[nn];
                rH[nn] = c = t = t - h;
                e = t - g;
                int[] wa = w[a[aoff + i + 1]];
                for (int j = nn - 1; j >= 0; j--)
                {
                    if ((c = c - m) > (e = e - h))
                        e = c;
                    if ((c = rH[j] - m) > (d = rE[j] - h))
                        d = c;
                    c = s + wa[b[boff + j + 1]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    s = rH[j];
                    rH[j] = c;
                    rE[j] = d;
                }
            }

            rE[nn] = rH[nn];
            int midc = fH[0] + rH[0];
            int midj = 0;
            int type = 1;
            for (int j = 0; j <= nn; j++)
            {
                if ((c = fH[j] + rH[j]) >= midc)
                {
                    if (c > midc || (fH[j] != fE[j] && rH[j] == rE[j]))
                    {
                        midc = c;
                        midj = j;
                    }
                }
            }

            for (int j = nn; j >= 0; j--)
            {
                if ((c = fE[j] + rE[j] + g) > midc)
                {
                    midc = c;
                    midj = j;
                    type = 2;
                }
            }

            if (type == 1)
            {
                NwAlign(a, aoff, b, boff, midi, midj, tb, -g, w, g, h, fH, fE, rH, rE, sc);
                NwAlign(a, aoff + midi, b, boff + midj, mm - midi, nn - midj, -g, te, w, g, h, fH, fE, rH, rE, sc);
            }
            else
            {
                NwAlign(a, aoff, b, boff, midi - 1, midj, tb, 0, w, g, h, fH, fE, rH, rE, sc);
                sc.Del(2);
                NwAlign(a, aoff + midi + 1, b, boff + midj, mm - midi - 1, nn - midj, 0, te, w, g, h, fH, fE, rH, rE, sc);
            }

            return midc;
        }

        /// <summary><c>CHECK_SCORE</c> / <c>BCHECK_SCORE</c>: the score of the script and its column count (<c>nres</c>).</summary>
        private static int CheckScore(int[] a, int aoff, int[] b, int boff, int mm, int nn, int[] s, int[][] w, int g, int h, out int nres)
        {
            int score = 0, i = 0, j = 0, nc = 0, p = 0;
            while (i < mm || j < nn)
            {
                int op = s[p++];
                if (op == 0)
                {
                    score = w[a[aoff + ++i]][b[boff + ++j]] + score;
                    nc++;
                }
                else if (op > 0)
                {
                    score -= g + op * h;
                    j += op;
                    nc += op;
                }
                else
                {
                    score -= g - op * h;
                    i -= op;
                    nc -= op;
                }
            }

            nres = nc;
            return score;
        }

        /// <summary><c>ALIGN</c>: <c>nw_align</c> followed by <c>CHECK_SCORE</c> (which also counts the columns).</summary>
        private static void Align(int[] a, int aoff, int[] b, int boff, int mm, int nn, int[][] w, int g, int h, AlnRes res, StringBuilder err)
        {
            var sc = new Script(res.Res);
            var fH = new int[nn + 2];
            var fE = new int[nn + 2];
            var rH = new int[nn + 2];
            var rE = new int[nn + 2];
            int c = NwAlign(a, aoff, b, boff, mm, nn, -g, -g, w, g, h, fH, fE, rH, rE, sc);
            res.Res = sc.S;
            int ck = CheckScore(a, aoff, b, boff, mm, nn, res.Res, w, g, h, out int nres);
            res.Nres = nres;
            if (c != ck)
                err.Append("Check_score error ").Append(I(c)).Append(" != ").Append(I(ck)).Append('\n');
        }

        /// <summary><c>sw_walign</c>: the Smith-Waterman alignment of the whole sequences (forward pass for the score and the end, backward pass for the start, <c>ALIGN</c> for the path). <paramref name="aa0"/> and <paramref name="aa1"/> are the padded arrays.</summary>
        private static int SwWalign(Molecule m, int[] aa0, int n0, int[] aa1, int n1, AlnRes ar, StringBuilder err)
        {
            int[][] pam = m.Pam2;
            int q = -m.GDel, r = -m.GGap, mm = q + r;
            var ssH = new int[n0 + 4];
            var ssE = new int[n0 + 4];
            for (int j = 0; j < n0; j++)
            {
                ssH[j] = 0;
                ssE[j] = -q;
            }

            int score = 0, bi = 0, bj = 0;
            for (int i = 0; i < n1; i++)
            {
                int h = 0, p = 0, f = -q, e;
                int code = aa1[i + 1];
                for (int j = 0; j < n0; j++)
                {
                    if ((h = h - mm) > (f = f - r))
                        f = h;
                    if ((h = ssH[j] - mm) > (e = ssE[j] - r))
                        e = h;
                    h = p + pam[aa0[j + 1]][code];
                    if (h < 0)
                        h = 0;
                    if (h < f)
                        h = f;
                    if (h < e)
                        h = e;
                    p = ssH[j];
                    ssH[j] = h;
                    ssE[j] = e;
                    if (h > score)
                    {
                        score = h;
                        bi = i;
                        bj = j;
                    }
                }
            }

            if (score <= 0)
                return 0;
            int cost = 0, kk = 0, ll = 0;
            for (int j = bj; j >= 0; j--)
            {
                ssH[j] = -1;
                ssE[j] = -1;
            }

            bool found = false;
            for (int i = bi; i >= 0 && !found; i--)
            {
                int h = -1, f = -1;
                int p = i == bi ? 0 : -1;
                for (int j = bj; j >= 0; j--)
                {
                    f = Math.Max(f, h - q) - r;
                    ssE[j] = Math.Max(ssE[j], ssH[j] - q) - r;
                    h = Math.Max(Math.Max(ssE[j], f), p + pam[aa0[j + 1]][aa1[i + 1]]);
                    p = ssH[j];
                    ssH[j] = h;
                    if (h > cost)
                    {
                        cost = h;
                        kk = i;
                        ll = j;
                        if (cost >= score)
                        {
                            found = true;
                            break;
                        }
                    }
                }
            }

            ar.Max0 = bj + 1;
            ar.Min0 = ll;
            ar.Max1 = bi + 1;
            ar.Min1 = kk;
            // ALIGN(&aa0[L-1], &aa1[K-1], J-L+1, I-K+1, ...): A[i] = aa0[L - 1 + i] of the program = aa0[L + i] of the padded array
            Align(aa0, ll, aa1, kk, bj - ll + 1, bi - kk + 1, pam, q, r, ar, err);
            return score;
        }

        // ------------------------------------------------------------------ banded alignment (LOCAL_ALIGN, bg_align, B_ALIGN)

        private const int MinInt = -9999999;

        /// <summary><c>LOCAL_ALIGN</c>: the best local alignment inside the band; the arrays are padded (A[i] = <c>a[i]</c>). Returns the score, or -1 for an empty band.</summary>
        private static int LocalAlign(int[] a, int[] b, int mm, int nn, int low, int up, int[][] w, int gg, int hh, out int psi, out int psj, out int pei, out int pej, int[] cc, int[] dd, StringBuilder err)
        {
            psi = psj = pei = pej = 0;
            int m = gg + hh;
            low = Math.Max(-mm, low);
            up = Math.Min(nn, up);
            if (nn <= 0 || mm <= 0)
                return 0;
            int band = up - low + 1;
            if (band < 1)
            {
                err.Append("low > up is unacceptable!: M: ").Append(I(mm)).Append(" N: ").Append(I(nn)).Append(" l/u: ").Append(I(low)).Append('/').Append(I(up)).Append('\n');
                return -1;
            }

            // cc[x] / dd[x] hold bss[x - 1]
            int leftd = low > 0 ? 1 : up < 0 ? band : 1 - low;
            int rightd = band;
            int si = Math.Max(0, -up);
            int ei = Math.Min(mm, nn - low);
            cc[leftd + 1] = 0;
            for (int j = leftd + 1; j <= rightd; j++)
            {
                cc[j + 1] = 0;
                dd[j + 1] = -gg;
            }

            cc[rightd + 2] = MinInt;
            dd[rightd + 2] = MinInt;
            int bestScore = 0;
            int endi = si;
            int endj = si + low;
            cc[leftd] = MinInt;
            dd[leftd + 1] = -gg;
            int c, d, e;
            for (int i = si + 1; i <= ei; i++)
            {
                if (i > nn - up)
                    rightd--;
                if (leftd > 1)
                    leftd--;
                int[] wa = w[a[i]];
                if ((c = cc[leftd + 2] - m) > (d = dd[leftd + 2] - hh))
                    d = c;
                int ib;
                if ((ib = leftd + low - 1 + i) > 0)
                    c = cc[leftd + 1] + wa[b[ib]];
                if (d > c)
                    c = d;
                if (c < 0)
                    c = 0;
                e = c - gg;
                dd[leftd + 1] = d;
                cc[leftd + 1] = c;
                if (c > bestScore)
                {
                    bestScore = c;
                    endi = i;
                    endj = ib;
                }

                for (int curd = leftd + 1; curd <= rightd; curd++)
                {
                    if ((c = c - m) > (e = e - hh))
                        e = c;
                    if ((c = cc[curd + 2] - m) > (d = dd[curd + 2] - hh))
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
                    {
                        bestScore = c;
                        endi = i;
                        endj = curd + low - 1 + i;
                    }
                }
            }

            leftd = Math.Max(1, -endi - low + 1);
            rightd = band - (up - (endj - endi));
            cc[rightd + 1] = 0;
            int t = -gg;
            for (int j = rightd - 1; j >= leftd; j--)
            {
                cc[j + 1] = t = t - hh;
                dd[j + 1] = t - gg;
            }

            for (int j = rightd + 1; j <= band; ++j)
                cc[j + 1] = MinInt;
            cc[leftd] = dd[leftd] = MinInt;
            dd[rightd + 1] = -gg;
            bool flag = false;
            int starti = 0, startj = 0;
            for (int i = endi; i >= 1; i--)
            {
                if (i + low <= 0)
                    leftd++;
                if (rightd < band)
                    rightd++;
                int[] wa = w[a[i]];
                if ((c = cc[rightd] - m) > (d = dd[rightd] - hh))
                    d = c;
                int ib;
                if ((ib = rightd + low - 1 + i) <= nn)
                    c = cc[rightd + 1] + wa[b[ib]];
                if (d > c)
                    c = d;
                e = c - gg;
                dd[rightd + 1] = d;
                cc[rightd + 1] = c;
                if (c == bestScore)
                {
                    starti = i;
                    startj = ib;
                    flag = true;
                    break;
                }

                for (int curd = rightd - 1; curd >= leftd; curd--)
                {
                    if ((c = c - m) > (e = e - hh))
                        e = c;
                    if ((c = cc[curd] - m) > (d = dd[curd] - hh))
                        d = c;
                    c = cc[curd + 1] + wa[b[curd + low - 1 + i]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    cc[curd + 1] = c;
                    dd[curd + 1] = d;
                    if (c == bestScore)
                    {
                        starti = i;
                        startj = curd + low - 1 + i;
                        flag = true;
                        break;
                    }
                }

                if (flag)
                    break;
            }

            psi = starti;
            psj = startj;
            pei = endi;
            pej = endj;
            return bestScore;
        }

        /// <summary>The state of <c>B_ALIGN</c> (the file-level arrays <c>MP</c>, <c>MT</c>, <c>FP</c>, <c>FT</c>, <c>IP</c> of <c>dropnfa.c</c>) and the band arrays of <c>f_str-&gt;bss</c>.</summary>
        private sealed class BandState
        {
            internal int Ip;
            internal int[][] Mp = new int[3][];
            internal int[] Fp = Array.Empty<int>();
            internal char[][] Mt = new char[3][];
            internal char[] Ft = Array.Empty<char>();
            internal int[] CC = Array.Empty<int>(), DD = Array.Empty<int>(), CP = Array.Empty<int>(), DP = Array.Empty<int>();
        }

        /// <summary><c>bg_align</c> (bss[x] of the program is stored at x).</summary>
        private static int BgAlign(int[] a, int aoff, int[] b, int boff, int mm, int nn, int low, int up, int tb, int te, int[][] w, int g, int h, BandState bs, Script sc)
        {
            int rmid, k, l, r, v, kt, t1, t2, t3;
            {
                int band, midd, leftd, rightd, i, j, c = 0, d = 0, e = 0, t, fr, ib, m;
                int[] wa;
                var CC = bs.CC;
                var DD = bs.DD;
                var CP = bs.CP;
                var DP = bs.DP;
                if (nn <= 0)
                {
                    if (mm > 0)
                        sc.Del(mm);
                    return 0;
                }

                if (mm <= 0)
                {
                    sc.Ins(nn);
                    return 0;
                }

                if ((band = up - low + 1) <= 1)
                {
                    for (i = 1; i <= mm; i++)
                        sc.Rep();
                    return 0;
                }

                m = g + h;
                midd = band / 2 + 1;
                rmid = low + midd - 1;
                leftd = 1 - low;
                rightd = up - low + 1;
                if (leftd < midd)
                {
                    fr = -1;
                    for (j = 0; j < midd; j++)
                        CP[j] = DP[j] = -1;
                    for (j = midd; j <= rightd; j++)
                        CP[j] = DP[j] = 0;
                    bs.Mp[0][0] = -1;
                    bs.Mp[1][0] = -1;
                    bs.Mp[2][0] = -1;
                    bs.Mt[0][0] = bs.Mt[1][0] = bs.Mt[2][0] = (char)0;
                }
                else if (leftd > midd)
                {
                    fr = leftd - midd;
                    for (j = 0; j <= midd; j++)
                        CP[j] = DP[j] = fr;
                    for (j = midd + 1; j <= rightd; j++)
                        CP[j] = DP[j] = -1;
                    bs.Mp[0][fr] = -1;
                    bs.Mp[1][fr] = -1;
                    bs.Mp[2][fr] = -1;
                    bs.Mt[0][fr] = bs.Mt[1][fr] = bs.Mt[2][fr] = (char)0;
                }
                else
                {
                    fr = 0;
                    for (j = 0; j < midd; j++)
                        CP[j] = DP[j] = 0;
                    for (j = midd; j <= rightd; j++)
                        CP[j] = DP[j] = 0;
                    bs.Mp[0][0] = -1;
                    bs.Mp[1][0] = -1;
                    bs.Mp[2][0] = -1;
                    bs.Mt[0][0] = bs.Mt[1][0] = bs.Mt[2][0] = (char)0;
                }

                CC[leftd] = 0;
                t = tb == 2 ? 0 : -g;
                for (j = leftd + 1; j <= rightd; j++)
                {
                    CC[j] = t = t - h;
                    DD[j] = t - g;
                }

                CC[rightd + 1] = MinInt;
                DD[rightd + 1] = MinInt;
                DD[leftd] = tb == 1 ? 0 : -g;
                CC[leftd - 1] = MinInt;
                for (i = 1; i <= mm; i++)
                {
                    if (i > nn - up)
                        rightd--;
                    if (leftd > 1)
                        leftd--;
                    wa = w[a[aoff + i]];
                    if ((c = CC[leftd + 1] - m) > (d = DD[leftd + 1] - h))
                    {
                        d = c;
                        DP[leftd] = CP[leftd + 1];
                    }
                    else
                    {
                        DP[leftd] = DP[leftd + 1];
                    }

                    if ((ib = leftd + low - 1 + i) > 0)
                        c = CC[leftd] + wa[b[boff + ib]];
                    if (d > c || ib <= 0)
                    {
                        c = d;
                        CP[leftd] = DP[leftd];
                    }

                    e = c - g;
                    DD[leftd] = d;
                    CC[leftd] = c;
                    bs.Ip = CP[leftd];
                    if (leftd == midd)
                        CP[leftd] = DP[leftd] = bs.Ip = i;
                    for (int curd = leftd + 1; curd <= rightd; curd++)
                    {
                        if (curd != midd)
                        {
                            if ((c = c - m) > (e = e - h))
                            {
                                e = c;
                                bs.Ip = CP[curd - 1];
                            }

                            if ((c = CC[curd + 1] - m) > (d = DD[curd + 1] - h))
                            {
                                d = c;
                                DP[curd] = CP[curd + 1];
                            }
                            else
                            {
                                DP[curd] = DP[curd + 1];
                            }

                            c = CC[curd] + wa[b[boff + curd + low - 1 + i]];
                            if (c < d || c < e)
                            {
                                if (e > d)
                                {
                                    c = e;
                                    CP[curd] = bs.Ip;
                                }
                                else
                                {
                                    c = d;
                                    CP[curd] = DP[curd];
                                }
                            }

                            CC[curd] = c;
                            DD[curd] = d;
                        }
                        else
                        {
                            if ((c = c - m) > (e = e - h))
                            {
                                e = c;
                                bs.Mp[1][i] = CP[curd - 1];
                                bs.Mt[1][i] = (char)2;
                            }
                            else
                            {
                                bs.Mp[1][i] = bs.Ip;
                                bs.Mt[1][i] = (char)2;
                            }

                            if ((c = CC[curd + 1] - m) > (d = DD[curd + 1] - h))
                            {
                                d = c;
                                bs.Mp[2][i] = CP[curd + 1];
                                bs.Mt[2][i] = (char)1;
                            }
                            else
                            {
                                bs.Mp[2][i] = DP[curd + 1];
                                bs.Mt[2][i] = (char)1;
                            }

                            c = CC[curd] + wa[b[boff + curd + low - 1 + i]];
                            if (c < d || c < e)
                            {
                                if (e > d)
                                {
                                    c = e;
                                    bs.Mp[0][i] = bs.Mp[1][i];
                                    bs.Mt[0][i] = (char)2;
                                }
                                else
                                {
                                    c = d;
                                    bs.Mp[0][i] = bs.Mp[2][i];
                                    bs.Mt[0][i] = (char)1;
                                }
                            }
                            else
                            {
                                bs.Mp[0][i] = i - 1;
                                bs.Mt[0][i] = (char)0;
                            }

                            if (c - g > e)
                            {
                                bs.Mp[1][i] = bs.Mp[0][i];
                                bs.Mt[1][i] = bs.Mt[0][i];
                            }

                            if (c - g > d)
                            {
                                bs.Mp[2][i] = bs.Mp[0][i];
                                bs.Mt[2][i] = bs.Mt[0][i];
                            }

                            CP[curd] = DP[curd] = bs.Ip = i;
                            CC[curd] = c;
                            DD[curd] = d;
                        }
                    }
                }

                if (te == 1 && d + g > c)
                {
                    k = DP[rightd];
                    l = 2;
                }
                else if (te == 2 && e + g > c)
                {
                    k = bs.Ip;
                    l = 1;
                }
                else
                {
                    k = CP[rightd];
                    l = 0;
                }

                if (rmid > nn - mm)
                    l = 2;
                else if (rmid < nn - mm)
                    l = 1;
                v = c;
            }

            r = -1;
            for (; k > -1; r = k, k = bs.Mp[l][r], l = bs.Mt[l][r])
            {
                bs.Fp[k] = r;
                bs.Ft[k] = (char)l;
            }

            if (r == -1)
            {
                if (rmid < 0)
                    BgAlign(a, aoff, b, boff, mm, nn, rmid + 1, up, tb, te, w, g, h, bs, sc);
                else
                    BgAlign(a, aoff, b, boff, mm, nn, low, rmid - 1, tb, te, w, g, h, bs, sc);
            }
            else
            {
                k = r;
                l = bs.Fp[k];
                kt = bs.Ft[k];
                if (rmid < 0)
                {
                    BgAlign(a, aoff, b, boff, r - 1, r + rmid, rmid + 1, Math.Min(up, r + rmid), tb, 1, w, g, h, bs, sc);
                    sc.Del(1);
                }
                else if (rmid > 0)
                {
                    BgAlign(a, aoff, b, boff, r, r + rmid - 1, Math.Max(-r, low), rmid - 1, tb, 2, w, g, h, bs, sc);
                    sc.Ins(1);
                }

                t2 = up - rmid - 1;
                t3 = low - rmid + 1;
                for (; l > -1; k = l, l = bs.Fp[k], kt = bs.Ft[k])
                {
                    if (kt == 0)
                    {
                        sc.Rep();
                    }
                    else if (kt == 1)
                    {
                        sc.Ins(1);
                        t1 = l - k - 1;
                        BgAlign(a, aoff + k, b, boff + k + rmid + 1, t1, t1, 0, Math.Min(t1, t2), 2, 1, w, g, h, bs, sc);
                        sc.Del(1);
                    }
                    else
                    {
                        sc.Del(1);
                        t1 = l - k - 1;
                        BgAlign(a, aoff + k + 1, b, boff + k + rmid, t1, t1, Math.Max(-t1, t3), 0, 1, 2, w, g, h, bs, sc);
                        sc.Ins(1);
                    }
                }

                if (nn - mm > rmid)
                {
                    sc.Ins(1);
                    t1 = k + rmid + 1;
                    BgAlign(a, aoff + k, b, boff + t1, mm - k, nn - t1, 0, Math.Min(nn - t1, t2), 2, te, w, g, h, bs, sc);
                }
                else if (nn - mm < rmid)
                {
                    sc.Del(1);
                    t1 = mm - (k + 1);
                    BgAlign(a, aoff + k + 1, b, boff + k + rmid, t1, nn - (k + rmid), Math.Max(-t1, t3), 0, 1, te, w, g, h, bs, sc);
                }
            }

            return v;
        }

        /// <summary><c>B_ALIGN</c>: the global alignment inside the band.</summary>
        private static void BAlign(int[] a, int aoff, int[] b, int boff, int mm, int nn, int low, int up, int[][] w, int g, int h, AlnRes ar, int mx, int bssSize)
        {
            var sc = new Script(ar.Res);
            low = Math.Min(Math.Max(-mm, low), Math.Min(nn - mm, 0));
            up = Math.Max(Math.Min(nn, up), Math.Max(nn - mm, 0));
            if (nn <= 0)
            {
                if (mm > 0)
                    sc.Del(mm);
                _ = -Gap(mm, g, h);
            }
            else if (mm <= 0)
            {
                sc.Ins(nn);
                _ = -Gap(nn, g, h);
            }
            else if (up - low + 1 <= 1)
            {
                for (int i = 1; i <= mm; i++)
                {
                    sc.Rep();
                }

            }
            else
            {
                var bs = new BandState();
                int len = mx + 1 + 8;
                for (int x = 0; x < 3; x++)
                {
                    bs.Mp[x] = new int[len];
                    bs.Mt[x] = new char[len];
                }

                bs.Fp = new int[len];
                bs.Ft = new char[len];
                int size = bssSize + 8;
                bs.CC = new int[size];
                bs.DD = new int[size];
                bs.CP = new int[size];
                bs.DP = new int[size];
                BgAlign(a, aoff, b, boff, mm, nn, low, up, 0, 0, w, g, h, bs, sc);
                ar.Res = sc.S;
                CheckScore(a, aoff, b, boff, mm, nn, ar.Res, w, g, h, out int nres);
                ar.Nres = nres;
                return;
            }

            ar.Res = sc.S;
        }

        /// <summary><c>bd_walign</c>: the banded alignment around the diagonal that <c>do_fasta</c> found (the <c>optwid</c> window).</summary>
        private static int BdWalign(Molecule m, int[] aa0, int n0, int[] aa1, int n1, int optWid, int hoff, AlnRes ar, FStr f, StringBuilder err)
        {
            int window = Math.Min(n1, optWid);
            int low = -window / 2 - hoff;
            int up = low + window;
            int score = LocalAlign(aa0, aa1, n0, n1, low, up, m.Pam2, -m.GDel, -m.GGap, out int min0, out int min1, out int max0, out int max1, f.BssCc, f.BssDd, err);
            if (score <= 0)
            {
                err.Append("n0/n1: ").Append(I(n0)).Append('/').Append(I(n1)).Append(" hoff: ").Append(I(hoff)).Append(" window: ").Append(I(window)).Append('\n');
                return 0;
            }

            ar.Min0 = min0 - 1;
            ar.Min1 = min1 - 1;
            ar.Max0 = max0;
            ar.Max1 = max1;
            // B_ALIGN(aa0 - 1 + min0 - 1, ...): A[i] = padded aa0[min0 - 1 + i]
            BAlign(aa0, min0 - 1, aa1, min1 - 1, max0 - min0 + 1, max1 - min1 + 1, low - (min1 - min0), up - (min1 - min0), m.Pam2, -m.GDel, -m.GGap, ar, n0, optWid * 4 + 16);
            return score;
        }

        /// <summary><c>do_walign</c>: the alignment of the query with one library sequence; sw_flag decides between the Smith-Waterman alignment and the banded one.</summary>
        private static int DoWalign(Molecule m, FStr f, int[] aa1, int n1, bool swFlag, AlnRes ar, StringBuilder err)
        {
            ar.Res = new int[Math.Max(3 * f.N0 / 2, 2000)];
            ar.Nres = 0;
            if (swFlag)
                return SwWalign(m, f.Aa0, f.N0, aa1, n1, ar, err);
            int optWid = m.OptWid * 2;
            var sc = DoFastaFull(m, f, aa1, n1, 0, optWid);
            if (sc.Initn > 0)
                return BdWalign(m, f.Aa0, f.N0, aa1, n1, optWid, sc.Hoff, ar, f, err);
            ar.Nres = 0;
            return 0;
        }

        // ------------------------------------------------------------------ the consensus of the alignment and its display (calcons, discons)

        private sealed class Aln
        {
            internal int Llen = 60, Llcntx = 30;
            internal int Smin0, Smin1, Amin0, Amax0, Amin1, Amax1;
            internal int Nident, Nsim, NgapQ, NgapL;
            internal bool Qlrev;
            internal long DStart0, DStop0, DStart1, DStop1;
        }

        /// <summary><c>aancpy</c>: the letters of <paramref name="count"/> residues (padded array, <paramref name="from"/> is the index of the first residue of the C array plus one) into a C string at <paramref name="to"/>.</summary>
        private static void Aancpy(char[] to, int toPos, int[] aa, int from, int count, string sq, int nsq)
        {
            int tp = toPos;
            int p = from;
            while (count-- > 0 && p < aa.Length && aa[p] != 0)
            {
                to[tp++] = aa[p] <= nsq ? sq[aa[p]] : (char)aa[p];
                p++;
            }

            to[tp] = '\0';
        }

        /// <summary><c>calcons</c>: the two aligned rows with their context (30 residues before, the rest after up to the line length), the alignment codes and the counts of identities. Returns the length of the rows.</summary>
        private static int Calcons(Molecule m, int[] aa0, int n0, int[] aa1, int n1, out int lc, Aln aln, AlnRes res, char[] seqc0, char[] seqc1, char[] seqca)
        {
            string sq = m.Nucleotide ? "\0ACGTURYMWSKDHVBNX" : "\0ARNDCQEGHILKMFPSTWYVBZX*";
            int nsq = m.Nsq;
            int mins, smins;
            aln.Amin0 = res.Min0;
            aln.Amax0 = res.Max0;
            aln.Amin1 = res.Min1;
            aln.Amax1 = res.Max1;
            if (Math.Min(res.Min0, res.Min1) < aln.Llen)
            {
                if (res.Min0 >= res.Min1)
                {
                    smins = 0;
                    mins = Math.Min(res.Min0, aln.Llcntx);
                    Aancpy(seqc0, 0, aa0, res.Min0 - mins + 1, mins, sq, nsq);
                    aln.Smin0 = res.Min0 - mins;
                    if (mins - res.Min1 > 0)
                    {
                        for (int x = 0; x < mins - res.Min1; x++)
                            seqc1[x] = ' ';
                        Aancpy(seqc1, mins - res.Min1, aa1, 1, res.Min1, sq, nsq);
                        aln.Smin1 = 0;
                    }
                    else
                    {
                        Aancpy(seqc1, 0, aa1, res.Min1 - mins + 1, mins, sq, nsq);
                        aln.Smin1 = res.Min1 - mins;
                    }
                }
                else
                {
                    smins = 0;
                    mins = Math.Min(res.Min1, aln.Llcntx);
                    Aancpy(seqc1, 0, aa1, res.Min1 - mins + 1, mins, sq, nsq);
                    aln.Smin1 = res.Min1 - mins;
                    if (mins - res.Min0 > 0)
                    {
                        for (int x = 0; x < mins - res.Min0; x++)
                            seqc0[x] = ' ';
                        Aancpy(seqc0, mins - res.Min0, aa0, 1, res.Min0, sq, nsq);
                        aln.Smin0 = 0;
                    }
                    else
                    {
                        Aancpy(seqc0, 0, aa0, res.Min0 - mins + 1, mins, sq, nsq);
                        aln.Smin0 = res.Min0 - mins;
                    }
                }
            }
            else
            {
                mins = Math.Min(aln.Llcntx, Math.Min(res.Min0, res.Min1));
                smins = mins;
                aln.Smin0 = res.Min0 - mins;
                aln.Smin1 = res.Min1 - mins;
                Aancpy(seqc0, 0, aa0, res.Min0 - mins + 1, mins, sq, nsq);
                Aancpy(seqc1, 0, aa1, res.Min1 - mins + 1, mins, sq, nsq);
            }

            _ = smins;
            for (int x = 0; x < mins; x++)
                seqca[x] = (char)0;
            int spa = mins, sp0 = mins, sp1 = mins;
            int rp = 0;
            int lenc = 0, op = 0;
            aln.Nident = aln.Nsim = aln.NgapQ = aln.NgapL = 0;
            int i0 = res.Min0, i1 = res.Min1;
            while (i0 < res.Max0 || i1 < res.Max1)
            {
                if (op == 0 && res.Res[rp] == 0)
                {
                    op = res.Res[rp++];
                    lenc++;
                    int itmp = m.Pam2[aa0[i0 + 1]][aa1[i1 + 1]];
                    if (itmp < 0)
                        seqca[spa] = (char)1;
                    else if (itmp == 0)
                        seqca[spa] = (char)2;
                    else
                        seqca[spa] = (char)3;
                    if (seqca[spa] == 3 || seqca[spa] == 2)
                        aln.Nsim++;
                    seqc0[sp0] = aa0[i0 + 1] <= nsq ? sq[aa0[i0 + 1]] : (char)aa0[i0 + 1];
                    i0++;
                    seqc1[sp1] = aa1[i1 + 1] <= nsq ? sq[aa1[i1 + 1]] : (char)aa1[i1 + 1];
                    i1++;
                    char u0 = char.ToUpperInvariant(seqc0[sp0]), u1 = char.ToUpperInvariant(seqc1[sp1]);
                    if (u0 == u1)
                    {
                        aln.Nident++;
                        seqca[spa] = (char)4;
                    }
                    else if (m.Nucleotide)
                    {
                        if ((u0 == 'T' && u1 == 'U') || (u0 == 'U' && u1 == 'T'))
                        {
                            aln.Nident++;
                            seqca[spa] = (char)4;
                        }
                        else if (u0 == 'N')
                        {
                            aln.NgapQ++;
                        }
                        else if (u1 == 'N')
                        {
                            aln.NgapL++;
                        }
                    }

                    sp0++;
                    sp1++;
                    spa++;
                }
                else
                {
                    if (op == 0)
                        op = res.Res[rp++];
                    if (op > 0)
                    {
                        seqc0[sp0++] = '-';
                        seqc1[sp1++] = aa1[i1 + 1] <= nsq ? sq[aa1[i1 + 1]] : (char)aa1[i1 + 1];
                        i1++;
                        seqca[spa++] = (char)5;
                        op--;
                        lenc++;
                        aln.NgapQ++;
                    }
                    else
                    {
                        seqc0[sp0++] = aa0[i0 + 1] <= nsq ? sq[aa0[i0 + 1]] : (char)aa0[i0 + 1];
                        i0++;
                        seqc1[sp1++] = '-';
                        seqca[spa++] = (char)5;
                        op++;
                        lenc++;
                        aln.NgapL++;
                    }
                }
            }

            lc = lenc;
            seqca[spa] = '\0';
            int ns = mins + lenc + aln.Llen;
            int itmp2;
            ns -= itmp2 = ns % aln.Llen;
            if (itmp2 > aln.Llen / 2)
                ns += aln.Llen;
            int nd = ns - (mins + lenc);
            if (nd > Math.Max(n0 - res.Max0, n1 - res.Max1))
                nd = Math.Max(n0 - res.Max0, n1 - res.Max1);
            if (nd - (n0 - res.Max0) > 0)
            {
                Aancpy(seqc0, mins + lenc, aa0, res.Max0 + 1, n0 - res.Max0, sq, nsq);
                for (int x = 0; x < nd - (n0 - res.Max0); x++)
                    seqc0[mins + lenc + n0 - res.Max0 + x] = ' ';
            }
            else
            {
                Aancpy(seqc0, mins + lenc, aa0, res.Max0 + 1, nd, sq, nsq);
            }

            if (nd - (n1 - res.Max1) > 0)
            {
                Aancpy(seqc1, mins + lenc, aa1, res.Max1 + 1, n1 - res.Max1, sq, nsq);
                for (int x = 0; x < nd - (n1 - res.Max1); x++)
                    seqc1[mins + lenc + n1 - res.Max1 + x] = ' ';
            }
            else
            {
                Aancpy(seqc1, mins + lenc, aa1, res.Max1 + 1, nd, sq, nsq);
            }

            return mins + lenc + nd;
        }

        /// <summary><c>cal_coord</c> (<c>sq0off</c> 1, the offset of the library sequence 0).</summary>
        private static void CalCoord(int n0, int n1, Aln aln)
        {
            long qoffset;
            int qlsgn, qfx0, qfxn;
            if (aln.Qlrev)
            {
                qoffset = n0;
                qlsgn = -1;
                qfx0 = 0;
                qfxn = 1;
            }
            else
            {
                qoffset = 0;
                qlsgn = 1;
                qfx0 = 1;
                qfxn = 0;
            }

            aln.DStart0 = qoffset + qlsgn * aln.Amin0 + qfx0;
            aln.DStop0 = qoffset + qlsgn * aln.Amax0 + qfxn;
            aln.DStart1 = 0 + aln.Amin1 + 1;
            aln.DStop1 = 0 + aln.Amax1;
        }

        private static string CString(char[] a)
        {
            int n = Array.IndexOf(a, '\0');
            return new string(a, 0, n < 0 ? a.Length : n);
        }

        /// <summary><c>discons</c> of <c>c_dispn.c</c> for <c>-m 10</c>: the two records with their coordinates and the aligned rows (50 residues per line).</summary>
        private static void Discons(StringBuilder o, char[] seqc0, char[] seqc1, int nc, int n0, int n1, string name0, string name1, string sqtype, Aln aln)
        {
            long qoffset = aln.Qlrev ? n0 : 0;
            int qlsgn = aln.Qlrev ? -1 : 1;
            int qfx0 = aln.Qlrev ? 0 : 1;
            int ioff0 = aln.Smin0, ioff1 = aln.Smin1;
            o.Append('>').Append(name0).Append(" ..\n");
            o.Append("; sq_len: ").Append(I(n0)).Append('\n');
            o.Append("; sq_offset: 1\n");
            o.Append("; sq_type: ").Append(sqtype[0]).Append('\n');
            o.Append("; al_start: ").Append(I(aln.DStart0)).Append('\n');
            o.Append("; al_stop: ").Append(I(aln.DStop0)).Append('\n');
            o.Append("; al_display_start: ").Append(I(qoffset + qlsgn * (long)ioff0 * 1 + qfx0)).Append('\n');
            Rows(o, seqc0, nc);
            o.Append('>').Append(name1).Append(" ..\n");
            o.Append("; sq_len: ").Append(I(n1)).Append('\n');
            o.Append("; sq_type: ").Append(sqtype[0]).Append('\n');
            o.Append("; al_start: ").Append(I(aln.DStart1)).Append('\n');
            o.Append("; al_stop: ").Append(I(aln.DStop1)).Append('\n');
            o.Append("; al_display_start: ").Append(I(0 + ioff1 + 1)).Append('\n');
            Rows(o, seqc1, nc);
        }

        private static void Rows(StringBuilder o, char[] seqc, int nc)
        {
            bool haveRes = false;
            int i;
            for (i = 0; i < nc && seqc[i] != '\0'; i++)
            {
                if (!haveRes && seqc[i] == ' ')
                {
                    o.Append('-');
                }
                else if (seqc[i] == ' ')
                {
                    break;
                }
                else
                {
                    haveRes = true;
                    o.Append(seqc[i]);
                }

                if (i % 50 == 49)
                    o.Append('\n');
            }

            if ((i - 1) % 50 != 49 || seqc[i - 1] == ' ')
                o.Append('\n');
        }
    }

    #endregion
}
