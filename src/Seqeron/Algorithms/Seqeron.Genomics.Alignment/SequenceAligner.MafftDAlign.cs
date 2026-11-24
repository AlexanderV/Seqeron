#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854, S125, S1871, S907, S2234
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --sp, D__align and D__align_variousdist (Dalignmm.c): the group-group dynamic program with the gap penalty of every cell from the distribution of the gap lengths of the two groups

    /// <summary>
    /// One entry of the gap-length table of a column (<c>Gaplen</c> of <c>Dalignmm.c</c>): a run of <see cref="Len"/> gaps of weight <see cref="Freq"/> that ends in this column (<see cref="Relend"/> 0) or covers it
    /// (<see cref="Relend"/> columns before its end); <see cref="Idatnext"/> / <see cref="Idatprev"/> are the indices of the entry of the same run in the next / the previous column (the terminator of a column has <c>Idatnext</c> -1).
    /// </summary>
    internal struct Gaplen
    {
        internal int Relend, Idatnext, Idatprev, Npat, Len;
        internal double Freq;
    }

    internal sealed partial class MafftRun
    {
        private readonly DpCapacity _dCap = new();

        /// <summary>The marker <c>(Gaplen *)1</c> that ends a row of column tables (the loops of the program stop there).</summary>
        private static readonly Gaplen[] GlOne = Array.Empty<Gaplen>();

        private static Gaplen[]?[] AllocGlRow(int count)
        {
            var row = new Gaplen[]?[count + 1];
            row[count] = GlOne;
            return row;
        }

        private static Gaplen Terminator() => new Gaplen { Len = -1, Idatnext = -1, Idatprev = -1, Relend = -1, Freq = 0.0, Npat = -1 };

        /// <summary><c>strralpha</c>: the number of gaps that directly precede the position <paramref name="j"/> of the row.</summary>
        private static int Strralpha(string s, int j)
        {
            int v = 0;
            for (int p = j - 1; p >= 0; p--)
            {
                if (s[p] != '-')
                    return v;
                v++;
            }
            return v;
        }

        private static char CharAtNul(string s, int j) => j < s.Length ? s[j] : '\0';

        /// <summary><c>gapvariety</c>: the number of (sequence, column) pairs that have a residue (or the end of the row) behind a run of gaps.</summary>
        private static int GapVariety(int n, int l, string[] seq)
        {
            int val = 0;
            for (int j = 0; j <= l; j++)
            {
                int nknown = 0;
                for (int i = 0; i < n; i++)
                {
                    if (CharAtNul(seq[i], j) == '-')
                        continue;
                    if (Strralpha(seq[i], j) > 0)
                        nknown++;
                }
                val += nknown;
            }
            return val;
        }

        /// <summary><c>fillgaplen</c>: the entries of the columns that a run of gaps covers (before the column where it ends), linked by <c>idatnext</c> / <c>idatprev</c>.</summary>
        private static void FillGaplen(Gaplen[]?[] mtx, int l)
        {
            for (int i = 0; i <= l; i++)
            {
                if (mtx[i] == null)
                    continue;
                for (int n = 0; mtx[i]![n].Len != -1; n++)
                {
                    int len = mtx[i]![n].Len;
                    double freq = mtx[i]![n].Freq;
                    int idatnext = n;
                    for (int j = 0; j < len; j++)
                    {
                        int pos = i - 1 - j;
                        int k;
                        if (mtx[pos] == null)
                        {
                            mtx[pos] = new Gaplen[2];
                            mtx[pos]![0] = Terminator();
                            k = 0;
                        }
                        else
                        {
                            k = mtx[pos]![0].Npat;
                        }
                        var row = mtx[pos]!;
                        Array.Resize(ref row, k + 2);
                        mtx[pos] = row;
                        row[k].Len = len;
                        row[k].Relend = j + 1;
                        row[k].Freq = freq;
                        row[k].Idatnext = idatnext;
                        row[k].Idatprev = -1;
                        row[k].Npat = -1;
                        row[k + 1].Len = -1;
                        row[k + 1].Idatnext = -1;
                        row[k + 1].Relend = -1;
                        row[k + 1].Freq = 0.0;
                        row[k + 1].Npat = -1;
                        row[0].Npat = k + 1;
                        mtx[pos + 1]![idatnext].Idatprev = k;
                        idatnext = k;
                    }
                }
            }
        }

        /// <summary><c>gaplencount</c>: for every column of a group the lengths of the runs of gaps that end there (<c>strralpha</c>), their weights, then <see cref="FillGaplen"/>.</summary>
        private void GaplenCount(int n, int l, Gaplen[]?[] mtx, string[] seq, double[] eff)
        {
            for (int j = 0; j <= l; j++)
            {
                if (ReferenceEquals(mtx[j], GlOne))
                    break;
                mtx[j] = null;
            }
            var known = new int[l + 1];
            for (int j = 0; j <= l; j++)
            {
                for (int i = 0; i < j; i++)
                    known[i] = 0;
                int nknown = 0;
                for (int i = 0; i < n; i++)
                {
                    if (CharAtNul(seq[i], j) == '-')
                        continue;
                    int gl = Strralpha(seq[i], j);
                    if (gl > 0)
                    {
                        if (known[gl] != 0)
                        {
                            int k;
                            for (k = 0; mtx[j]![k].Len != -1; k++)
                                if (mtx[j]![k].Len == gl)
                                    break;
                            if (mtx[j]![k].Len == -1)
                                throw new MafftProgramFailureException("Unexpected error!\n");
                            mtx[j]![k].Freq += eff[i];
                        }
                        else
                        {
                            var row = mtx[j] ?? Array.Empty<Gaplen>();
                            Array.Resize(ref row, nknown + 2);
                            mtx[j] = row;
                            row[nknown].Len = gl;
                            row[nknown].Relend = 0;
                            row[nknown].Freq = eff[i];
                            row[nknown].Idatnext = -2;
                            row[nknown + 1].Len = -1;
                            row[nknown + 1].Idatnext = -1;
                            row[nknown + 1].Relend = -1;
                            row[nknown + 1].Freq = 0.0;
                            row[nknown + 1].Npat = -1;
                            known[gl]++;
                            nknown++;
                            row[0].Npat = nknown;
                        }
                    }
                }
            }
            FillGaplen(mtx, l);
        }

        /// <summary><c>pairgapcount</c>: the number of runs of gaps in the two rows after the columns that both have a gap in are removed.</summary>
        private static int PairGapCount(string s1, string s2)
        {
            int v = 0;
            bool st1 = false, st2 = false;
            for (int i = 0; i < s1.Length; i++)
            {
                if (s1[i] == '-' && s2[i] == '-')
                    continue;
                if (s1[i] == '-')
                {
                    if (!st1)
                        v++;
                    st1 = true;
                }
                else
                {
                    st1 = false;
                }
                if (s2[i] == '-')
                {
                    if (!st2)
                        v++;
                    st2 = true;
                }
                else
                {
                    st2 = false;
                }
            }
            return v;
        }

        // ---- the probabilities of a new gap, an extension and a gap that ends / starts in a cell (the weight of the penalty)

        private static double CalcPfacGapNoIdatEnd(Gaplen[]?[] gaplen1, Gaplen[]?[] gaplen2, int newgaplen, int i, int j)
        {
            Gaplen[]? gaplen2j = gaplen2[j];
            Gaplen[]? gaplen1i = gaplen1[i];
            double pfac = 0.0, pfac1 = 0.0, pfac10 = 0.0, pfac2;
            if (gaplen1i != null)
            {
                for (int k = 0; gaplen1i[k].Idatnext != -1; k++)
                {
                    var g1 = gaplen1i[k];
                    int pos1 = g1.Relend;
                    if (pos1 != 0)
                    {
                        pfac2 = 0.0;
                        if (gaplen2j != null)
                        {
                            for (int l = 0; gaplen2j[l].Idatnext != -1; l++)
                            {
                                var g2 = gaplen2j[l];
                                int pos2 = g2.Relend;
                                if (pos2 == 0 && g2.Len > g1.Len - pos1 + newgaplen)
                                    pfac2 += g2.Freq;
                            }
                        }
                        pfac += pfac2 * g1.Freq;
                        pfac1 += g1.Freq;
                    }
                    else
                    {
                        pfac2 = 1.0;
                        if (gaplen2j != null)
                        {
                            for (int l = 0; gaplen2j[l].Idatnext != -1; l++)
                            {
                                var g2 = gaplen2j[l];
                                int pos2 = g2.Relend;
                                if (pos2 == 0 && g2.Len == g1.Len + newgaplen)
                                    pfac2 -= g2.Freq;
                                if (pos2 != 0 && g2.Len - (pos2 - 1) > g1.Len + newgaplen)
                                    pfac2 -= g2.Freq;
                            }
                        }
                        pfac += pfac2 * g1.Freq;
                        pfac10 += g1.Freq;
                    }
                }
            }
            pfac1 = 1.0 - pfac1 - pfac10;
            pfac2 = 1.0;
            if (gaplen2j != null)
            {
                for (int l = 0; gaplen2j[l].Idatnext != -1; l++)
                {
                    var g2 = gaplen2j[l];
                    int pos2 = g2.Relend;
                    if (pos2 == 0 && g2.Len == newgaplen)
                        pfac2 -= g2.Freq;
                    if (pos2 != 0 && g2.Len - (pos2 - 1) > newgaplen)
                        pfac2 -= g2.Freq;
                }
            }
            pfac += pfac1 * pfac2;
            return pfac;
        }

        private static double CalcPfacGapexNoIdatEnd(Gaplen[]?[] gaplen1, Gaplen[]?[] gaplen2, int i, int j, int newgaplen)
        {
            Gaplen[]? gaplen1i = gaplen1[i];
            Gaplen[]? gaplen2j = gaplen2[j];
            double pfac = 0.0;
            if (gaplen2j != null)
            {
                for (int k = 0; gaplen2j[k].Idatnext != -1; k++)
                {
                    var g2 = gaplen2j[k];
                    int pos2 = g2.Relend;
                    if (pos2 != 0)
                        continue;
                    double pfac1 = 1.0;
                    double pfac10 = 1.0;
                    if (gaplen1i != null)
                    {
                        for (int l = 0; gaplen1i[l].Idatnext != -1; l++)
                        {
                            var g1 = gaplen1i[l];
                            int pos1 = g1.Relend;
                            pfac10 -= g1.Freq;
                            if (newgaplen + g1.Len - pos1 > g2.Len - pos2)
                                pfac1 -= g1.Freq;
                        }
                    }
                    pfac += pfac1 * g2.Freq;
                    if (newgaplen >= g2.Len - (pos2 - 1))
                        pfac -= pfac10 * g2.Freq;
                }
            }
            return pfac;
        }

        private static double CalcPfacNoIdatEnd(Gaplen[]?[] gaplen1, Gaplen[]?[] gaplen2, int i, int j)
        {
            Gaplen[]? gaplen1i = gaplen1[i];
            Gaplen[]? gaplen2j = gaplen2[j];
            double pfac1 = 0.0, pfac2 = 0.0;
            if (gaplen1i != null)
                for (int k = 0; gaplen1i[k].Idatnext != -1; k++)
                    if (gaplen1i[k].Relend == 0)
                        pfac1 += gaplen1i[k].Freq;
            if (gaplen2j != null)
                for (int l = 0; gaplen2j[l].Idatnext != -1; l++)
                    if (gaplen2j[l].Relend == 0)
                        pfac2 += gaplen2j[l].Freq;
            double pfac = pfac1 * pfac2 + pfac1 * (1 - pfac2) + pfac2 * (1 - pfac1);
            if (gaplen1i != null)
            {
                for (int k = 0; gaplen1i[k].Idatnext != -1; k++)
                {
                    var g1 = gaplen1i[k];
                    int pos1 = g1.Relend;
                    if (gaplen2j != null)
                    {
                        for (int l = 0; gaplen2j[l].Idatnext != -1; l++)
                        {
                            var g2 = gaplen2j[l];
                            int pos2 = g2.Relend;
                            if (pos1 == 0 && pos2 == 0 && g1.Len == g2.Len)
                                pfac -= g1.Freq * g2.Freq;
                            else if (pos1 == 0 && pos2 != 0 && g2.Len - (pos2 - 1) > g1.Len)
                                pfac -= g1.Freq * g2.Freq;
                            else if (pos1 != 0 && pos2 == 0 && g1.Len - (pos1 - 1) > g2.Len)
                                pfac -= g1.Freq * g2.Freq;
                        }
                    }
                }
            }
            return pfac;
        }

        // ---- copies of the column tables (the program keeps one table per cell it may come back to)

        /// <summary><c>extendgaplencompactx</c>: the table of the column <paramref name="start"/> of <paramref name="orig"/> (from the offset <paramref name="origOff"/>) is copied when <paramref name="cpy"/> has none, and the lengths of the runs that
        /// started in the column before are taken from the table of that column.</summary>
        private static void ExtendGaplenCompactX(Gaplen[]?[] cpy, Gaplen[]?[] orig, int origOff, int start)
        {
            var os = orig[origOff + start];
            if (os == null)
            {
                if (cpy[start] != null)
                    cpy[start] = null;
                return;
            }
            if (cpy[start] == null)
            {
                int l = os[0].Npat;
                var c = new Gaplen[l + 2];
                int q = 0;
                for (; os[q].Idatnext != -1; q++)
                    c[q] = os[q];
                c[q] = os[q];
                cpy[start] = c;
            }
            if (start == 0)
                return;
            var prev = cpy[start - 1];
            if (prev == null)
                return;
            var cs = cpy[start]!;
            for (int l = 0; os[l].Idatnext != -1; l++)
            {
                int id = os[l].Idatprev;
                if (id == -1)
                    continue;
                cs[l].Len = prev[id].Len;
            }
        }

        /// <summary><c>duplicategaplencompactx</c>: the tables of the columns <paramref name="start"/> ... <paramref name="end"/> of <paramref name="orig"/> (from the offset <paramref name="origOff"/>) are copied into <paramref name="cpy"/> (the end of the row is the marker).</summary>
        private static void DuplicateGaplenCompactX(Gaplen[]?[] cpy, Gaplen[]?[] orig, int origOff, int start, int end)
        {
            if (start < 0)
                start = 0;
            for (int i = start; i <= end; i++)
            {
                if (ReferenceEquals(cpy[i], GlOne))
                {
                    end = i - 1;
                    break;
                }
                cpy[i] = null;
            }
            for (int i = start; i <= end; i++)
            {
                var o = orig[origOff + i];
                if (o == null)
                {
                    cpy[i] = null;
                    continue;
                }
                int l = o[0].Npat;
                var c = new Gaplen[l + 2];
                for (l = 0; o[l].Idatnext != -1; l++)
                    c[l] = o[l];
                c[l] = o[l];
                cpy[i] = c;
            }
        }

        /// <summary><c>gaplenextendnoidatend</c>: the runs of the table of <paramref name="gapstartpos"/> become <paramref name="insertionlen"/> longer (the entry of the run in the column before is changed too).</summary>
        private static void GaplenExtendNoIdatEnd(Gaplen[]?[] cpy, int gapstartpos, int insertionlen)
        {
            for (int l = 0; cpy[gapstartpos] != null && cpy[gapstartpos]![l].Idatnext != -1; l++)
            {
                int len = cpy[gapstartpos]![l].Len + insertionlen;
                int pos = gapstartpos, idn = l;
                if (cpy[pos] != null && cpy[pos]![idn].Relend != 0)
                {
                    cpy[pos]![idn].Len = len;
                    idn = cpy[pos]![idn].Idatnext;
                }
                idn = cpy[gapstartpos]![l].Idatprev;
                if (gapstartpos != 0 && idn != -1 && cpy[gapstartpos - 1] != null)
                    cpy[gapstartpos - 1]![idn].Len = len;
            }
        }

        /// <summary><c>copygaplencompactx</c>: the lengths of the table <paramref name="posinori"/> of <paramref name="orig"/> are written into the table <paramref name="posincopy"/> of <paramref name="cpy"/>, then <see cref="GaplenExtendNoIdatEnd"/>.</summary>
        private static void CopyGaplenCompactX(Gaplen[]?[] cpy, Gaplen[]?[] orig, int gapstartpos, int insertionlen, int posincopy, int posinori)
        {
            var o = orig[posinori];
            if (o == null)
                return;
            var c = cpy[posincopy];
            for (int q = 0; o[q].Idatnext != -1; q++)
                c![q].Len = o[q].Len;
            if (gapstartpos == -1)
                gapstartpos = posincopy;
            GaplenExtendNoIdatEnd(cpy, gapstartpos, insertionlen);
        }

        /// <summary><c>freegaplenpartly</c>.</summary>
        private static void FreeGaplenPartly(Gaplen[]?[] mtx, int startpos, int endpos)
        {
            if (startpos < 0)
                startpos = 0;
            for (int i = startpos; i <= endpos; i++)
            {
                if (ReferenceEquals(mtx[i], GlOne))
                    break;
                mtx[i] = null;
            }
        }

        /// <summary><c>naivepairscore11_dynmtx</c>: the score of two aligned rows with the matrix <paramref name="mtx"/>, a run of gaps in one row costs <paramref name="penal"/> once.</summary>
        private double NaivePairScore11Dyn(double[][] mtx, string seq1, string seq2, int penal)
        {
            var s1 = new StringBuilder(seq1.Length);
            var s2 = new StringBuilder(seq1.Length);
            for (int i = 0; i < seq1.Length; i++)
            {
                if (seq1[i] == '-' && seq2[i] == '-')
                    continue;
                s1.Append(seq1[i]);
                s2.Append(seq2[i]);
            }
            double vali = 0.0;
            int a = 0, n = s1.Length;
            while (a < n)
            {
                if (s1[a] == '-')
                {
                    vali += penal;
                    while (a < n && s1[a] == '-')
                        a++;
                    continue;
                }
                if (s2[a] == '-')
                {
                    vali += penal;
                    while (a < n && s2[a] == '-')
                        a++;
                    continue;
                }
                int c1 = G.AminoN[s1[a]], c2 = G.AminoN[s2[a]];
                if (c1 < 0)
                    throw UndefinedBehaviour("naivepairscore11_dynmtx reads mtx[-1] for the character '" + s1[a] + "' (the new gap character of --add that is not a letter of the alphabet): the program reads a pointer from the header of the heap block and stops with a segmentation fault");
                if (c2 >= 0) // mtx[c1][-1] is the size field of the heap block read as a double (a denormal number, 0 in the sum)
                    vali += mtx[c1][c2];
                a++;
            }
            return vali;
        }

        /// <summary><c>D__align_ls</c>: the group with more (sequence, column) pairs behind a run of gaps is the first one of <see cref="DAlign"/> (the rows of both groups are replaced in place).</summary>
        private void DAlignLs(double[][] dyn, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, int alloclen, ref double impmatch, int headgp, int tailgp)
        {
            int v1 = GapVariety(icyc, seq1[0].Length, seq1);
            int v2 = GapVariety(jcyc, seq2[0].Length, seq2);
            Report(" v1,v2 = " + v1.ToString(CultureInfo.InvariantCulture) + "," + v2.ToString(CultureInfo.InvariantCulture) + "\n");
            if (v1 >= v2)
                DAlign(dyn, seq1, seq2, eff1, eff2, icyc, jcyc, alloclen, null, ref impmatch, headgp, tailgp);
            else
                DAlign(dyn, seq2, seq1, eff2, eff1, jcyc, icyc, alloclen, null, ref impmatch, headgp, tailgp);
        }

        /// <summary>
        /// <c>D__align</c> (<c>--sp</c>) and, with <paramref name="various"/>, <c>D__align_variousdist</c> (<c>Dalignmm.c</c>): the dynamic program of two groups of aligned rows in which the penalty of a gap is not the constant of the scoring model but
        /// the penalty weighted by the probability that the gap is new (it does not continue a gap of the group that the other group's gaps cover); the probabilities come from the table of the gap lengths of every column of both groups
        /// (<see cref="GaplenCount"/>), which the program copies for every cell it may come back to. <paramref name="imp"/> is the importance matrix of a pair library (<c>constraint</c>); <paramref name="impmatch"/> receives the
        /// sum of the importance on the path then. The program prints the progress, the score and a check of the score by the pairwise scores of the result to the standard error, and stops with the groups on the standard output when the check fails.
        /// </summary>
        internal double DAlign(double[][] dyn, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int icyc, int jcyc, int alloclen, double[][]? imp, ref double impmatch, int headgp, int tailgp, VariousDist? various = null)
        {
            bool constraint = imp != null;
            double fpenalty = G.Penalty;
            double fpenaltyShift = G.PenaltyShift;
            int na = G.Nalphabets;
            ClassProfiles? classes = null;
            if (various != null)
            {
                bool complex = false;
                for (int c = 0; c < various.MaxDistClass && !complex; c++)
                    for (int i = 0; i < icyc && !complex; i++)
                        for (int j = 0; j < jcyc; j++)
                            if (various.Eff1s[c][i] * various.Eff2s[c][j] != 0.0 && c != various.Which[i][j])
                            {
                                complex = true;
                                break;
                            }
                if (complex)
                    Report("Found a complex grouping. This step may be a bit slow.\n");
            }
            int lgth1 = seq1[0].Length;
            int lgth2 = seq2[0].Length;
            Report(icyc.ToString(CultureInfo.InvariantCulture) + " x " + jcyc.ToString(CultureInfo.InvariantCulture) + " sequences, len=" + lgth1.ToString(CultureInfo.InvariantCulture) + ", " + lgth2.ToString(CultureInfo.InvariantCulture) + "\n");
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
            int warpbase = lgth1 + lgth2;
            bool trywarp = G.TryWarp;
            WarpState? warp = null;
            if (trywarp)
            {
                if (headgp == 0 || tailgp == 0)
                    throw new MafftProgramFailureException("At present, headgp and tailgp must be 1 to allow shift.\n");
                warp = new WarpState(lgth1, lgth2, G.PenaltyShift);
            }
            _dCap.Fit(lgth1, lgth2);
            var ijp = G.Common.IpFor(_dCap.Orlgth1, _dCap.Orlgth2);

            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            var lastverticalw = new double[lgth1 + 2];
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var cpmx1 = new double[na][];
            var cpmx2 = new double[na][];
            for (int i = 0; i < na; i++)
            {
                cpmx1[i] = new double[lgth1 + 2];
                cpmx2[i] = new double[lgth2 + 2];
            }
            int ll = Math.Max(lgth1, lgth2) + 2;
            double[][]? doublework = null;
            int[][]? intwork = null;
            if (various == null)
            {
                doublework = new double[ll][];
                intwork = new int[ll][];
                for (int i = 0; i < ll; i++)
                {
                    doublework[i] = new double[na];
                    intwork[i] = new int[na + 1];
                }
            }
            else
            {
                classes = new ClassProfiles(this, various, seq1, seq2, icyc, jcyc, lgth1, lgth2, false);
            }

            // the tables of the gap lengths of the columns of both groups, and the copies that the dynamic program keeps
            var gaplen1 = AllocGlRow(lgth1 + 1);
            var gaplen2 = AllocGlRow(lgth2 + 1);
            var gaplen1ibestkamo = new Gaplen[]?[lgth1 + 2][];
            var gaplen2ibestkamo = new Gaplen[]?[lgth1 + 2][];
            var gaplen1icurr = new Gaplen[]?[lgth1 + 2][];
            var gaplen2icurr = new Gaplen[]?[lgth1 + 2][];
            for (int i = 0; i < lgth1 + 1; i++)
            {
                gaplen1ibestkamo[i] = AllocGlRow(lgth1 + 1 - i);
                gaplen2ibestkamo[i] = AllocGlRow(lgth2 + 1);
                gaplen1icurr[i] = AllocGlRow(lgth1 + 1 - i);
                gaplen2icurr[i] = AllocGlRow(lgth2 + 1);
            }
            var gaplen1jbestkamo = new Gaplen[]?[lgth2 + 2][];
            var gaplen2jbestkamo = new Gaplen[]?[lgth2 + 2][];
            var gaplen1jbest = new Gaplen[]?[lgth2 + 2][];
            var gaplen2jbest = new Gaplen[]?[lgth2 + 2][];
            var gaplen1jcurr = new Gaplen[]?[lgth2 + 2][];
            var gaplen2jcurr = new Gaplen[]?[lgth2 + 2][];
            var gaplen1jprev = new Gaplen[]?[lgth2 + 2][];
            var gaplen2jprev = new Gaplen[]?[lgth2 + 2][];
            for (int i = 0; i < lgth2 + 1; i++)
            {
                gaplen1jbestkamo[i] = AllocGlRow(lgth1 + 1);
                gaplen2jbestkamo[i] = AllocGlRow(lgth2 + 1 - i);
                gaplen1jbest[i] = AllocGlRow(lgth1 + 1);
                gaplen2jbest[i] = AllocGlRow(lgth2 + 1 - i);
                gaplen1jcurr[i] = AllocGlRow(lgth1 + 1);
                gaplen2jcurr[i] = AllocGlRow(lgth2 + 1 - i);
                gaplen1jprev[i] = AllocGlRow(lgth1 + 1);
                gaplen2jprev[i] = AllocGlRow(lgth2 + 1 - i);
            }

            if (various == null)
            {
                CpmxCalcNew(seq1, cpmx1, eff1, lgth1, icyc);
                CpmxCalcNew(seq2, cpmx2, eff2, lgth2, jcyc);
            }
            GaplenCount(icyc, lgth1, gaplen1, seq1, eff1);
            GaplenCount(jcyc, lgth2, gaplen2, seq2, eff2);

            // i = 0
            DuplicateGaplenCompactX(gaplen1icurr[0], gaplen1, 0, 0, 1);
            DuplicateGaplenCompactX(gaplen2icurr[0], gaplen2, 0, 0, 0);
            CopyGaplenCompactX(gaplen2icurr[0], gaplen2, 0, 0, 0, 0);
            DuplicateGaplenCompactX(gaplen1ibestkamo[0], gaplen1, 0, 0, 0);
            DuplicateGaplenCompactX(gaplen2ibestkamo[0], gaplen2, 0, 0, 0);
            CopyGaplenCompactX(gaplen2ibestkamo[0], gaplen2, 0, 0, 0, 0);
            for (int j = 0; j < lgth2 + 1; j++)
            {
                DuplicateGaplenCompactX(gaplen1jcurr[j], gaplen1, 0, 0, 0);
                DuplicateGaplenCompactX(gaplen2jcurr[j], gaplen2, j, 0, 0);
                CopyGaplenCompactX(gaplen1jcurr[j], gaplen1, 0, j, 0, 0);
            }
            for (int j = 0; j < lgth2 + 1; j++)
            {
                DuplicateGaplenCompactX(gaplen1jprev[j], gaplen1, 0, 0, 0);
                DuplicateGaplenCompactX(gaplen2jprev[j], gaplen2, j, 0, 0);
                CopyGaplenCompactX(gaplen1jprev[j], gaplen1, 0, j, 0, 0);
            }
            for (int j = 0; j < lgth2 + 1; j++)
            {
                DuplicateGaplenCompactX(gaplen1jbestkamo[j], gaplen1, 0, 0, 0);
                DuplicateGaplenCompactX(gaplen2jbestkamo[j], gaplen2, j, 0, 1);
                CopyGaplenCompactX(gaplen1jbestkamo[j], gaplen1, 0, j, 0, 0);
                DuplicateGaplenCompactX(gaplen1jbest[j], gaplen1, 0, 0, 0);
                DuplicateGaplenCompactX(gaplen2jbest[j], gaplen2, j, 0, 1);
                CopyGaplenCompactX(gaplen1jbest[j], gaplen1, 0, j, 0, 0);
            }

            double[] currentw = w1, previousw = w2;
            if (various == null)
            {
                MatchCalc(dyn, initverticalw, cpmx2, cpmx1, 0, lgth1, doublework!, intwork!, true);
            }
            else
            {
                classes!.Match(initverticalw, true, 0, lgth1, true, seq1, seq2, eff1, eff2);
            }
            if (constraint)
                for (int i = 0; i < lgth1; i++)
                    initverticalw[i] += imp![i][0];
            if (various == null)
                MatchCalc(dyn, currentw, cpmx1, cpmx2, 0, lgth2, doublework!, intwork!, true);
            else
                classes!.Match(currentw, false, 0, lgth2, true, seq1, seq2, eff1, eff2);
            if (constraint)
                for (int j = 0; j < lgth2; j++)
                    currentw[j] += imp![0][j];
            for (int j = 1; j < lgth2 + 1; j++)
            {
                double pfac = CalcPfacGapNoIdatEnd(gaplen1, gaplen2, j, 0, j);
                currentw[j] += fpenalty * pfac;
            }
            for (int i = 1; i < lgth1 + 1; i++)
            {
                double pfac = CalcPfacGapNoIdatEnd(gaplen2, gaplen1, i, 0, i);
                initverticalw[i] += fpenalty * pfac;
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                double pfac = CalcPfacGapexNoIdatEnd(gaplen2, gaplen1, j, 1, j);
                m[j] = currentw[j - 1] + fpenalty * pfac;
                mp[j] = 0;
            }
            lastverticalw[0] = lgth2 == 0 ? 0.0 : currentw[lgth2 - 1];
            int lasti = tailgp != 0 ? lgth1 + 1 : lgth1;
            int lastj = lgth2 + 1;
            double wm = 0.0, g, mi;
            int mpi;
            double pfacL;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                (gaplen1jprev, gaplen1jcurr) = (gaplen1jcurr, gaplen1jprev);
                (gaplen2jprev, gaplen2jcurr) = (gaplen2jcurr, gaplen2jprev);
                DuplicateGaplenCompactX(gaplen1icurr[i], gaplen1, i, 0, 1);
                DuplicateGaplenCompactX(gaplen2icurr[i], gaplen2, 0, 0, 0);
                CopyGaplenCompactX(gaplen2icurr[i], gaplen2, 0, i, 0, 0);
                DuplicateGaplenCompactX(gaplen1ibestkamo[i], gaplen1, i, 0, 1);
                DuplicateGaplenCompactX(gaplen2ibestkamo[i], gaplen2, 0, 0, 0);
                CopyGaplenCompactX(gaplen2ibestkamo[i], gaplen2, 0, i, 0, 0);
                ExtendGaplenCompactX(gaplen1jprev[0], gaplen1, 0, i);
                CopyGaplenCompactX(gaplen1jprev[0], gaplen1icurr[i - 1], -1, 0, i, 1);
                CopyGaplenCompactX(gaplen2jprev[0], gaplen2icurr[i - 1], -1, 0, 0, 0);
                if (various == null)
                    MatchCalc(dyn, currentw, cpmx1, cpmx2, i, lgth2, doublework!, intwork!, false);
                else
                    classes!.Match(currentw, false, i, lgth2, false, seq1, seq2, eff1, eff2);
                if (constraint)
                {
                    var impRow = imp![i];
                    for (int j = 0; j < lgth2; j++)
                        currentw[j] += impRow[j];
                }
                currentw[0] = initverticalw[i];
                pfacL = CalcPfacGapexNoIdatEnd(gaplen1, gaplen2, i, 1, i);
                mi = previousw[0] + fpenalty * pfacL;
                mpi = 0;
                var ijpRow = ijp[i];
                for (int j = 1; j < lastj; j++)
                {
                    ExtendGaplenCompactX(gaplen1icurr[i - 1], gaplen1, i - 1, 1);
                    ExtendGaplenCompactX(gaplen2icurr[i - 1], gaplen2, 0, j);
                    ExtendGaplenCompactX(gaplen1jprev[j - 1], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen1jcurr[j - 1], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen2jprev[j - 1], gaplen2, j - 1, 1);
                    ExtendGaplenCompactX(gaplen2jcurr[j - 1], gaplen2, j - 1, 1);
                    ExtendGaplenCompactX(gaplen1jbest[j - 1], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen1jbestkamo[j - 1], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen1jprev[mpi], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen2jprev[mpi], gaplen2, mpi, j - mpi);
                    ExtendGaplenCompactX(gaplen1ibestkamo[i - 1], gaplen1, i - 1, 1);
                    ExtendGaplenCompactX(gaplen2ibestkamo[i - 1], gaplen2, 0, j);

                    double pfac = CalcPfacNoIdatEnd(gaplen1jprev[j - 1], gaplen2jprev[j - 1], i, 1);
                    wm = previousw[j - 1] + fpenalty * pfac;
                    ijpRow[j] = 0;
                    int newgaplen = j - mpi - 1;
                    pfac = CalcPfacGapNoIdatEnd(gaplen1ibestkamo[i - 1], gaplen2ibestkamo[i - 1], newgaplen, 1, j);
                    if ((g = mi + fpenalty * pfac) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = previousw[j - 1]) >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                        FreeGaplenPartly(gaplen2ibestkamo[i - 1], j - 3, j - 2);
                        CopyGaplenCompactX(gaplen1ibestkamo[i - 1], gaplen1jprev[j - 1], -1, 0, 1, i);
                        CopyGaplenCompactX(gaplen2ibestkamo[i - 1], gaplen2jprev[j - 1], -1, 0, j, 1);
                    }
                    pfac = CalcPfacGapexNoIdatEnd(gaplen1ibestkamo[i - 1], gaplen2ibestkamo[i - 1], 1, j, j - mpi);
                    mi += fpenalty * pfac;
                    newgaplen = i - mp[j] - 1;
                    pfac = CalcPfacGapNoIdatEnd(gaplen2jbestkamo[j - 1], gaplen1jbestkamo[j - 1], newgaplen, 1, i);
                    if ((g = m[j] + fpenalty * pfac) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                        FreeGaplenPartly(gaplen1jbest[j - 1], i - 3, i - 2);
                        CopyGaplenCompactX(gaplen1jbest[j - 1], gaplen1jbestkamo[j - 1], -1, 0, i, i);
                        CopyGaplenCompactX(gaplen2jbest[j - 1], gaplen2jbestkamo[j - 1], -1, 0, 1, 1);
                    }
                    if ((g = previousw[j - 1]) >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                        FreeGaplenPartly(gaplen1jbestkamo[j - 1], i - 3, i - 2);
                        CopyGaplenCompactX(gaplen1jbestkamo[j - 1], gaplen1jprev[j - 1], -1, 0, i, i);
                        CopyGaplenCompactX(gaplen2jbestkamo[j - 1], gaplen2jprev[j - 1], -1, 0, 1, 1);
                    }
                    pfac = CalcPfacGapexNoIdatEnd(gaplen2jbestkamo[j - 1], gaplen1jbestkamo[j - 1], 1, i, i - mp[j]);
                    m[j] += fpenalty * pfac;
                    if (warp != null)
                        wm = warp.StepD(i, j, wm, currentw[j], fpenaltyShift, ijpRow);
                    ExtendGaplenCompactX(gaplen1jcurr[j], gaplen1, 0, i);
                    ExtendGaplenCompactX(gaplen2jcurr[j], gaplen2, j, 0);
                    if (ijpRow[j] < 0)
                    {
                        CopyGaplenCompactX(gaplen1jcurr[j], gaplen1jprev[j + ijpRow[j]], i, -ijpRow[j] - 1, i, i);
                        CopyGaplenCompactX(gaplen2jcurr[j], gaplen2jprev[j + ijpRow[j]], -1, 0, 0, -ijpRow[j]);
                    }
                    else if (ijpRow[j] > 0)
                    {
                        CopyGaplenCompactX(gaplen1jcurr[j], gaplen1jbest[j - 1], -1, 0, i, i);
                        CopyGaplenCompactX(gaplen2jcurr[j], gaplen2jbest[j - 1], 0, ijpRow[j] - 1, 0, 1);
                    }
                    else
                    {
                        CopyGaplenCompactX(gaplen1jcurr[j], gaplen1jprev[j - 1], -1, 0, i, i);
                        CopyGaplenCompactX(gaplen2jcurr[j], gaplen2jprev[j - 1], -1, 0, 0, 1);
                    }
                    FreeGaplenPartly(gaplen1jcurr[j - 1], i - 3, i - 2);
                    FreeGaplenPartly(gaplen1jbestkamo[j - 1], i - 3, i - 2);
                    FreeGaplenPartly(gaplen1jbest[j - 1], i - 3, i - 2);
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
                FreeGaplenPartly(gaplen1icurr[i - 1], 0, lgth1 - i);
                FreeGaplenPartly(gaplen2icurr[i - 1], 0, lgth2);
                FreeGaplenPartly(gaplen1ibestkamo[i - 1], 0, lgth1 - i);
                FreeGaplenPartly(gaplen2ibestkamo[i - 1], 0, lgth2);
                warp?.EndRow(lastj);
            }

            string gt1, gt2;
            double impwm = 0.0;
            try
            {
                (gt1, gt2, impwm) = DTracking(constraint, imp, currentw, lastverticalw, seq1, seq2, ijp, tailgp, warp, warpbase);
            }
            catch (IndexOutOfRangeException)
            {
                throw StaleTraceBack();
            }
            if (constraint)
                impmatch = impwm;
            var mseq1 = new string[icyc];
            var mseq2 = new string[jcyc];
            for (int i = 0; i < icyc; i++)
                mseq1[i] = GapIreru(seq1[i], gt1, G.NewGapChar);
            for (int j = 0; j < jcyc; j++)
                mseq2[j] = GapIreru(seq2[j], gt2, G.NewGapChar);
            int resultlen = mseq1[0].Length;
            if (alloclen < resultlen || resultlen > 5000000)
                throw new MafftProgramFailureException("alloclen=" + alloclen.ToString(CultureInfo.InvariantCulture) + ", resultlen=" + resultlen.ToString(CultureInfo.InvariantCulture) + ", N=5000000\nLENGTH OVER!\n");
            for (int i = 0; i < icyc; i++)
                seq1[i] = mseq1[i];
            for (int j = 0; j < jcyc; j++)
                seq2[j] = mseq2[j];

            Report("wm = " + MafftGlobals.F(wm - impmatch, 0, 6) + "\n");
            Report("*impmatch = " + MafftGlobals.F(impmatch, 0, 6) + "\n");
            int kenzan = 0;
            for (int i = 0; i < icyc; i++)
                for (int j = 0; j < jcyc; j++)
                    kenzan += PairGapCount(mseq1[i], mseq2[j]);
            Report("kenzan = " + kenzan.ToString(CultureInfo.InvariantCulture) + " -> " + MafftGlobals.F((double)kenzan / (icyc * jcyc), 0, 6) + "\n");
            double pairscore = 0.0, nogappairscore = 0.0;
            for (int i = 0; i < icyc; i++)
            {
                for (int j = 0; j < jcyc; j++)
                {
                    var pseq = new[] { seq1[i], seq2[j] };
                    CommonGapPick(pseq);
                    var mtx = various == null ? dyn : various.Matrices[various.Which[i][j]];
                    pairscore += eff1[i] * eff2[j] * NaivePairScore11Dyn(mtx, pseq[0], pseq[1], G.Penalty);
                    nogappairscore += eff1[i] * eff2[j] * NaivePairScore11Dyn(mtx, pseq[0], pseq[1], 0);
                }
            }
            double diff;
            if (various == null)
            {
                Report("pairscore = " + MafftGlobals.F(pairscore, 0, 6) + "\n");
                Report("pairscore-nogappairscore = " + MafftGlobals.F(pairscore - nogappairscore, 0, 6) + "\n");
                Report("pairscore-nogappairscore / penalty = " + MafftGlobals.F((pairscore - nogappairscore) / fpenalty, 0, 6) + "\n");
                diff = (pairscore - wm + impmatch) / fpenalty;
                Report("diff = " + MafftGlobals.F(diff, 0, 6) + "\n\n");
                if (!trywarp && Math.Abs(pairscore - wm + impmatch) > 0.01)
                    throw DAlignCheckFailure();
            }
            else
            {
                diff = (pairscore - wm + impmatch) / seq1[0].Length;
                Report("pairscore = " + MafftGlobals.F(pairscore, 0, 6) + "\n");
                Report("pairscore-nogappairscore = " + MafftGlobals.F(pairscore - nogappairscore, 0, 6) + "\n");
                Report("pairscore-nogappairscore / penalty = " + MafftGlobals.F((pairscore - nogappairscore) / fpenalty, 0, 6) + "\n");
                Report("diff = " + MafftGlobals.F(diff, 0, 6) + "\n\n");
                if (!trywarp && Math.Abs(diff) > 0.01)
                    throw DAlignCheckFailure();
            }
            return wm;
        }

        /// <summary>The check of <c>D__align</c> failed (the score of the dynamic program is not the sum of the pairwise scores of the result): the groups are written to the standard output and the program stops with exit status 1.</summary>
        private static MafftProgramFailureException DAlignCheckFailure() => new MafftProgramFailureException("", 1);

        /// <summary><c>Atracking</c> / <c>Atracking_localhom</c> of <c>Dalignmm.c</c>: the trace back; returns the gap tables of the two groups and, with the importance matrix, the sum of the importance of the cells of the path.</summary>
        private (string Gt1, string Gt2, double Imp) DTracking(bool localhom, double[][]? imp, double[] lasthorizontalw, double[] lastverticalw, string[] seq1, string[] seq2, int[][] ijp, int tailgp, WarpState? warp, int warpbase)
        {
            int lgth1 = seq1[0].Length, lgth2 = seq2[0].Length;
            if (tailgp != 1)
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
            var gt1 = new char[lgth1 + lgth2 + 1];
            var gt2 = new char[lgth1 + lgth2 + 1];
            int p1 = lgth1 + lgth2, p2 = lgth1 + lgth2;
            int iin = lgth1, jin = lgth2;
            int limk = lgth1 + lgth2 + 1;
            double impwm = 0.0;
            for (int k = 0; k < limk; k++)
            {
                int ifi, jfi;
                if (warp != null && ijp[iin][jin] >= warpbase)
                {
                    ifi = warp.Is[ijp[iin][jin] - warpbase];
                    jfi = warp.Js[ijp[iin][jin] - warpbase];
                }
                else if (ijp[iin][jin] < 0)
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
                if (warp != null && ifi == -warpbase && jfi == -warpbase)
                {
                    int rest = iin;
                    while (--rest >= 0)
                    {
                        gt1[--p1] = 'o';
                        gt2[--p2] = '-';
                        k++;
                    }
                    rest = jin;
                    while (--rest >= 0)
                    {
                        gt1[--p1] = '-';
                        gt2[--p2] = 'o';
                    }
                    break;
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
                if (localhom && iin != lgth1 && jin != lgth2)
                    impwm += imp![iin][jin];
                if (iin <= 0 || jin <= 0)
                    break;
                gt1[--p1] = 'o';
                gt2[--p2] = 'o';
                k++;
                iin = ifi;
                jin = jfi;
            }
            return (new string(gt1, p1, lgth1 + lgth2 - p1), new string(gt2, p2, lgth1 + lgth2 - p2), impwm);
        }
    }

    #endregion
}
