using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: disttbfast (the progressive stage for three or more sequences)

    internal sealed partial class MafftRun
    {
        /// <summary>The result of one run of <c>disttbfast</c>: the aligned rows in the order of the input, the order in which the guide tree visits the sequences (the file <c>order</c> of <c>--reorder</c>) and the messages of the program (stderr).</summary>
        internal sealed record Result(string[] Rows, int[] Order);

        /// <summary>The stack size that <c>setstacksize</c> reports (the limit of the process: 8192 kb by default; the message is the only place where it shows).</summary>
        internal int StackSizeKb = 8192;

        private int _njob;
        private int _nlenmax;
        private int[] _nlen = Array.Empty<int>();
        private int[] _nogaplen = Array.Empty<int>();
        private int _nunknown;

        private void Report(string s) => Log.Append(s);

        // -------------------------------------------------------------------- the front end of disttbfast

        /// <summary>
        /// <c>disttbfast</c> (<c>main</c> of <c>disttbfast.c</c> for the strategies of the plain progressive method FFT-NS-1 / -2 and NW-NS-1 / -2): the 6-mer (10-mer for nucleotides) distance, the UPGMA guide tree, the
        /// progressive alignment (<c>treebase</c>), and for <c>-E 2</c> the distances of the first alignment, a second tree and a second progressive alignment.
        /// <paramref name="seq"/> are the sequences as the program holds them after <c>replaceu</c>.
        /// </summary>
        internal Result Disttbfast(IReadOnlyList<string> seq0)
        {
            int njob = seq0.Count;
            _njob = njob;
            var seq = new string[njob];
            for (int i = 0; i < njob; i++)
                seq[i] = seq0[i];
            _nlenmax = 0;
            for (int i = 0; i < njob; i++)
                _nlenmax = Math.Max(_nlenmax, seq[i].Length);
            _nlen = new int[njob];
            for (int i = 0; i < njob; i++)
                _nlen[i] = seq[i].Length;
            _nogaplen = new int[njob];
            if (G.Dorp == 'd')
                for (int i = 0; i < njob; i++)
                    seq[i] = seq[i].ToLowerInvariant();

            Report("stacksize: " + StackSizeKb.ToString(CultureInfo.InvariantCulture) + " kb\n");
            G.Constants(seq);
            Report("Gap Penalty = " + MafftGlobals.F((double)G.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.Poffset / 1000, 5, 2, true) + "\n");
            for (int i = 0; i < njob; i++)
            {
                for (int k = 0; k < seq[i].Length; k++)
                {
                    if (G.AminoN[seq[i][k]] == -1)
                        throw new ArgumentException($"Alphabet '{seq[i][k]}' is unknown (MAFFT stops: try --anysymbol).");
                }
            }
            Report("\n");
            if (G.Dorp == 'p' && G.TupleSize != 6)
                throw new ArgumentException("tuplesize must be 6 for aa sequence");
            if (G.Dorp == 'd' && G.TupleSize != 6 && G.TupleSize != 10)
                throw new ArgumentException("tuplesize must be 6 or 10 for dna sequence");
            if (G.TreeIn != 0 || G.NoAlign != 0 || G.DistOut != 0 || G.TreeOut != 0)
                throw new NotSupportedException("--treein / --treeout / --distout of MAFFT are not available in this port.");

            double[][] mtx = SixmerDistanceMatrix(seq);
            int nguidetree = G.NGuideTree;
            int[] order = new int[njob];
            string[]? bseq = null;
            for (int iguidetree = 0; iguidetree < nguidetree; iguidetree++)
            {
                bool calcpairdists = iguidetree != nguidetree - 1;
                var topol0 = new int[njob][];
                var topol1 = new int[njob][];
                var len = new double[njob][];
                for (int i = 0; i < njob; i++)
                    len[i] = new double[2];
                var dep = new Treedep[njob];
                Report("Constructing a UPGMA tree (efffree=" + (calcpairdists ? 0 : 1) + ") ... ");
                MuscleSupg(njob, mtx, topol0, topol1, len, dep);
                Report("\ndone.\n\n");
                var tree = new GuideTree(topol0, topol1, dep, len);
                int[] lm = tree.Order(njob - 2, 2);
                Array.Copy(lm, order, njob);
                double[] eff = new double[njob];
                if (G.Tbrweight != 0)
                {
                    G.Weight = 3;
                    CountEffSimple(njob, tree, eff);
                }
                else
                {
                    for (int i = 0; i < njob; i++)
                        eff[i] = 1.0;
                }
                bseq = new string[njob];
                for (int i = 0; i < njob; i++)
                    bseq[i] = seq[i].Replace("-", "");
                int alloclen = _nlenmax * 2 + 1;
                int[] selfscore = Array.Empty<int>();
                if (calcpairdists)
                {
                    selfscore = new int[njob];
                    for (int i = 0; i < njob; i++)
                        selfscore[i] = (int)NaivePairScore11(seq[i], seq[i], G.PenaltyDist);
                }
                Report("Progressive alignment " + (iguidetree + 1) + "/" + nguidetree + "... \n");
                TreeBase(bseq, tree, eff, ref alloclen);
                Report("\ndone.\n\n");
                if (calcpairdists)
                {
                    Report("Making a distance matrix from msa.. \n");
                    mtx = MsaDistanceMatrix(njob, selfscore, bseq);
                    Report("\ndone.\n\n");
                }
            }
            Report("disttbfast (" + (G.Dorp == 'd' ? "nuc" : "aa") + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new Result(bseq!, order);
        }

        // -------------------------------------------------------------------- distances

        private static int[] MakePointTable(int[] n, int count, int tuple, int[] mult, int baseN)
        {
            // pointt = the code of every window of 'tuple' consecutive groups (base 'baseN'); empty when the sequence has fewer than 'tuple' groups
            if (count < tuple)
                return Array.Empty<int>();
            var res = new int[count - tuple + 1];
            int point = 0;
            for (int i = 0; i < tuple; i++)
                point += n[i] * mult[i];
            res[0] = point;
            int p = 0;
            for (int k = tuple; k < count; k++)
            {
                point -= n[p++] * mult[0];
                point *= baseN;
                point += n[k];
                res[k - tuple + 1] = point;
            }
            return res;
        }

        private double[][] SixmerDistanceMatrix(string[] seq)
        {
            int njob = _njob;
            bool dna = G.Dorp == 'd';
            int tuple = G.TupleSize;
            int baseN = dna ? 4 : 6;
            int tsize = dna ? (int)Math.Pow(4, tuple) : (int)Math.Pow(6, 6);
            int[] mult;
            if (dna)
            {
                mult = new int[tuple];
                int v = 1;
                for (int i = tuple - 1; i >= 0; i--)
                {
                    mult[i] = v;
                    v *= 4;
                }
            }
            else
            {
                mult = new[] { 7776, 1296, 216, 36, 6, 1 };
            }
            double lenfaca, lenfacb, lenfacc, lenfacd;
            if (dna && tuple == 6)
            {
                lenfaca = 0.01; lenfacb = 2500; lenfacc = 2500; lenfacd = 0.1;
            }
            else if (dna && tuple == 10)
            {
                lenfaca = 0.01; lenfacb = 1000000; lenfacc = 1000000; lenfacd = 0.0;
            }
            else
            {
                lenfaca = 0.01; lenfacb = 10000; lenfacc = 10000; lenfacd = 0.1;
            }
            Report("\n\nMaking a distance matrix ..\n");
            var pointt = new int[njob][];
            int maxl = 0;
            for (int i = 0; i < njob; i++)
            {
                string tmpseq = seq[i].Replace("-", "");
                _nogaplen[i] = tmpseq.Length;
                if (_nogaplen[i] > maxl)
                    maxl = _nogaplen[i];
                var grp = new List<int>(tmpseq.Length);
                int limit = dna ? 4 : 6;
                foreach (char c in tmpseq)
                {
                    int tmp = G.AminoGrp[c];
                    if (tmp < limit)
                        grp.Add(tmp);
                    else
                        _nunknown++;
                }
                pointt[i] = grp.Count < tuple ? Array.Empty<int>() : MakePointTable(grp.ToArray(), grp.Count, tuple, mult, baseN);
            }
            if (_nunknown != 0)
                Report("\nThere are " + _nunknown + " ambiguous characters.\n");
            var mtx = new double[njob][];
            for (int i = 0; i < njob; i++)
                mtx[i] = new double[njob - i];
            var memo = new int[tsize];
            for (int i = 0; i < njob; i++)
            {
                var table1 = new int[tsize];
                if (i % 100 == 0)
                    Report("\r" + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + " / " + njob);
                foreach (int point in pointt[i])
                    table1[point]++;
                for (int j = i; j < njob; j++)
                    mtx[i][j - i] = CommonSextet(table1, pointt[j], memo);
            }
            Report("\ndone.\n\n");
            int ien = njob - 1;
            for (int i = 0; i < ien; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    double longer, shorter;
                    if (_nogaplen[i] > _nogaplen[j])
                    {
                        longer = _nogaplen[i];
                        shorter = _nogaplen[j];
                    }
                    else
                    {
                        longer = _nogaplen[j];
                        shorter = _nogaplen[i];
                    }
                    double lenfac = 1.0 / (shorter / longer * lenfacd + lenfacb / (longer + lenfacc) + lenfaca);
                    double bunbo = Math.Min(mtx[i][0], mtx[j][0]);
                    if (bunbo == 0.0)
                        mtx[i][j - i] = 2.0;
                    else
                        mtx[i][j - i] = (1.0 - mtx[i][j - i] / bunbo) * lenfac * 2.0;
                }
            }
            return mtx;
        }

        /// <summary><c>commonsextet_p</c>: the number of 6-mers (10-mers) that two sequences share, a k-mer counted as often as it occurs in both.</summary>
        private static int CommonSextet(int[] table, int[] pointt, int[] memo)
        {
            if (pointt.Length == 0)
                return 0;
            int value = 0;
            foreach (int point in pointt)
            {
                int tmp = memo[point]++;
                if (tmp < table[point])
                    value++;
            }
            foreach (int point in pointt)
                memo[point] = 0;
            return value;
        }

        /// <summary><c>naivepairscore11</c>: the score of a pair of aligned sequences, a run of gaps costs <paramref name="penal"/> once.</summary>
        private double NaivePairScore11(string seq1, string seq2, int penal)
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
            return PairScore(s1.ToString(), s2.ToString(), penal);
        }

        private double PairScore(string p1, string p2, int penal)
        {
            double vali = 0.0;
            int a = 0;
            int n = p1.Length;
            while (a < n)
            {
                if (p1[a] == '-')
                {
                    vali += penal;
                    while (a < n && p1[a] == '-')
                        a++;
                    continue;
                }
                if (p2[a] == '-')
                {
                    vali += penal;
                    while (a < n && p2[a] == '-')
                        a++;
                    continue;
                }
                vali += G.AminoDis[p1[a]][p2[a]];
                a++;
            }
            return vali;
        }

        /// <summary><c>naivepairscorefast</c> with the skip tables of <c>makeskiptable</c>: the same score as <see cref="NaivePairScore11"/> (the columns with a gap in both sequences are dropped).</summary>
        private double NaivePairScoreFast(string seq1, string seq2, int[] skip1, int[] skip2, int penal)
        {
            var r1 = new StringBuilder(seq1.Length);
            var r2 = new StringBuilder(seq1.Length);
            int skipped1 = 0, skipped2 = 0;
            int a = 0, k1 = 0, k2 = 0;
            while (true)
            {
                int skip = Math.Min(skip1[k1] - skipped1, skip2[k2] - skipped2);
                a += skip;
                skipped1 += skip;
                skipped2 += skip;
                if (a >= seq1.Length)
                    break;
                if (seq1[a] != '-')
                {
                    skipped1 = 0;
                    k1++;
                }
                else
                {
                    skipped1++;
                }
                if (seq2[a] != '-')
                {
                    skipped2 = 0;
                    k2++;
                }
                else
                {
                    skipped2++;
                }
                r1.Append(seq1[a]);
                r2.Append(seq2[a]);
                a++;
            }
            return PairScore(r1.ToString(), r2.ToString(), penal);
        }

        private static int[] MakeSkipTable(string seq)
        {
            int nogaplen = 0;
            foreach (char c in seq)
                if (c != '-')
                    nogaplen++;
            var skip = new int[nogaplen + 1];
            int posinseq = 0;
            foreach (char c in seq)
            {
                if (c == '-')
                    skip[posinseq]++;
                else
                    posinseq++;
            }
            return skip;
        }

        /// <summary><c>msadistmtxthread</c>: the distance of every pair of rows of the alignment, 2 (1 - score / min(self scores)), at most 10.</summary>
        private double[][] MsaDistanceMatrix(int njob, int[] selfscore, string[] bseq)
        {
            var skip = new int[njob][];
            for (int i = 0; i < njob; i++)
                skip[i] = MakeSkipTable(bseq[i]);
            var iscore = new double[njob][];
            for (int i = 0; i < njob; i++)
                iscore[i] = new double[njob - i];
            int nlim = njob - 1;
            for (int i = 0; i < nlim; i++)
            {
                if (i % 100 == 0)
                    Report("\r" + i.ToString(CultureInfo.InvariantCulture).PadLeft(5) + " / " + njob);
                double ssi = selfscore[i];
                for (int j = i + 1; j < njob; j++)
                {
                    double ssj = selfscore[j];
                    double bunbo = Math.Min(ssi, ssj);
                    double iscoretmp;
                    if (bunbo == 0.0)
                    {
                        iscoretmp = 2.0;
                    }
                    else
                    {
                        iscoretmp = (1.0 - NaivePairScoreFast(bseq[i], bseq[j], skip[i], skip[j], G.PenaltyDist) / bunbo) * 2.0;
                        if (iscoretmp > 10)
                            iscoretmp = 10.0;
                    }
                    if (iscoretmp < 0.0)
                    {
                        Report("WARNING: negative distance, iscoretmp = " + iscoretmp.ToString("F6", CultureInfo.InvariantCulture) + "\n");
                        iscoretmp = 0.0;
                    }
                    iscore[i][j - i] = iscoretmp;
                }
            }
            return iscore;
        }

        // -------------------------------------------------------------------- the guide tree

        internal struct Treedep
        {
            internal int Child0, Child1;
            internal double DistFromTip;
        }

        /// <summary>The guide tree of the program: <c>topol[k][0/1]</c> (the representative of the two groups joined at step k), <c>dep</c> (the steps that produced the groups) and the branch lengths.</summary>
        internal sealed class GuideTree
        {
            internal readonly int[] Topol0, Topol1;
            internal readonly Treedep[] Dep;
            internal readonly double[][] Len;

            internal GuideTree(int[][] t0, int[][] t1, Treedep[] dep, double[][] len)
            {
                int n = t0.Length;
                Topol0 = new int[n];
                Topol1 = new int[n];
                for (int i = 0; i < n; i++)
                {
                    Topol0[i] = t0[i] is { Length: > 0 } ? t0[i][0] : -1;
                    Topol1[i] = t1[i] is { Length: > 0 } ? t1[i][0] : -1;
                }
                Dep = dep;
                Len = len;
            }

            /// <summary><c>topolorderz</c>: the leaves below the group <paramref name="nchild"/> (0, 1 or 2 = both) of step <paramref name="pos"/> in the order of the tree.</summary>
            internal int[] Order(int pos, int nchild)
            {
                var res = new List<int>();
                if (nchild == 0 || nchild == 2)
                {
                    int child = Dep[pos].Child0;
                    if (child == -1)
                        res.Add(Topol0[pos]);
                    else
                        res.AddRange(Order(child, 2));
                }
                if (nchild == 1 || nchild == 2)
                {
                    int child = Dep[pos].Child1;
                    if (child == -1)
                        res.Add(Topol1[pos]);
                    else
                        res.AddRange(Order(child, 2));
                }
                return res.ToArray();
            }
        }

        /// <summary><c>fixed_musclesupg_double_realloc_nobk_halfmtx_memsave</c>: the UPGMA (mixed linkage, <c>-X</c>) of the half matrix <paramref name="eff"/> (MUSCLE's nearest-neighbour bookkeeping); step k joins the two closest groups.</summary>
        private void MuscleSupg(int nseq, double[][] eff, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool fullMembership = false)
        {
            double sueff1 = 1 - G.SueffGlobal;
            double sueff05 = G.SueffGlobal * 0.5;
            Func<double, double, double> cluster = G.TreeMethod switch
            {
                'X' => (d1, d2) => Math.Min(d1, d2) * sueff1 + (d1 + d2) * sueff05,
                'E' => (d1, d2) => (d1 + d2) * 0.5,
                'q' => (d1, d2) => Math.Min(d1, d2),
                _ => throw new InvalidOperationException("Unknown treemethod"),
            };
            var next = new int[nseq];
            var prev = new int[nseq];
            for (int i = 0; i < nseq; i++)
            {
                next[i] = i + 1;
                prev[i] = i - 1;
            }
            next[nseq - 1] = -1;
            var hist = new int[nseq];
            var tmptmplen = new double[nseq];
            var nmemar = new int[nseq];
            var mindisfrom = new double[nseq];
            var nearest = new int[nseq];
            for (int i = 0; i < nseq; i++)
                SetNearest(eff, next, mindisfrom, nearest, i);
            for (int i = 0; i < nseq; i++)
            {
                hist[i] = -1;
                nmemar[i] = 1;
            }
            Report("\n");
            int im = -1, jm = -1;
            for (int k = 0; k < nseq - 1; k++)
            {
                if (k % 10 == 0)
                    Report("\r" + k.ToString(CultureInfo.InvariantCulture).PadLeft(5) + " / " + nseq);
                double minscore = 999.9;
                for (int a = 0; next[a] != -1; a = next[a])
                {
                    int i = a;
                    if (mindisfrom[i] < minscore)
                    {
                        im = i;
                        minscore = mindisfrom[i];
                    }
                }
                jm = nearest[im];
                if (jm < im)
                    (im, jm) = (jm, im);
                int prevnode = hist[im];
                dep[k].Child0 = prevnode;
                int nmemim = nmemar[im];
                topol0[k] = GroupMembers(prevnode, im, topol0, topol1, fullMembership);
                prevnode = hist[jm];
                dep[k].Child1 = prevnode;
                int nmemjm = nmemar[jm];
                topol1[k] = GroupMembers(prevnode, jm, topol0, topol1, fullMembership);
                minscore *= 0.5;
                len[k][0] = minscore - tmptmplen[im];
                len[k][1] = minscore - tmptmplen[jm];
                dep[k].DistFromTip = minscore;
                tmptmplen[im] = minscore;
                hist[im] = k;
                nmemar[im] = nmemim + nmemjm;
                mindisfrom[im] = 999.9;
                for (int a = 0; a != -1; a = next[a])
                {
                    int i = a;
                    if (i != im && i != jm)
                    {
                        int miniim, maxiim, minijm, maxijm;
                        if (i < im)
                        {
                            miniim = i; maxiim = im; minijm = i; maxijm = jm;
                        }
                        else if (i < jm)
                        {
                            miniim = im; maxiim = i; minijm = i; maxijm = jm;
                        }
                        else
                        {
                            miniim = im; maxiim = i; minijm = jm; maxijm = i;
                        }
                        double eff0 = eff[miniim][maxiim - miniim];
                        double eff1 = eff[minijm][maxijm - minijm];
                        double tmpdouble = eff[miniim][maxiim - miniim] = cluster(eff0, eff1);
                        if (tmpdouble < mindisfrom[i])
                        {
                            mindisfrom[i] = tmpdouble;
                            nearest[i] = im;
                        }
                        if (tmpdouble < mindisfrom[im])
                        {
                            mindisfrom[im] = tmpdouble;
                            nearest[im] = i;
                        }
                        if (nearest[i] == jm)
                            nearest[i] = im;
                    }
                }
                int acjmprev = prev[jm];
                int acjmnext = next[jm];
                next[acjmprev] = acjmnext;
                if (acjmnext != -1)
                    prev[acjmnext] = acjmprev;
                for (int a = 0; a != -1; a = next[a])
                {
                    int i = a;
                    if (nearest[i] == im)
                    {
                        int miniim, maxiim;
                        if (i < im)
                        {
                            miniim = i; maxiim = im;
                        }
                        else
                        {
                            miniim = im; maxiim = i;
                        }
                        if (eff[miniim][maxiim - miniim] > mindisfrom[i])
                            SetNearest(eff, next, mindisfrom, nearest, i);
                    }
                }
            }
        }

        /// <summary>The members of the group formed at step <paramref name="prevnode"/> (the two groups, the one whose first member is smaller first); only the first member with the memory-saving tree of <c>disttbfast</c>.</summary>
        private static int[] GroupMembers(int prevnode, int leaf, int[][] topol0, int[][] topol1, bool full)
        {
            if (prevnode == -1)
                return new[] { leaf };
            int[] pt1 = topol0[prevnode], pt2 = topol1[prevnode];
            if (!full)
                return new[] { Math.Min(pt1[0], pt2[0]) };
            int[] a = pt1, b = pt2;
            if (pt1[0] > pt2[0])
            {
                a = pt2;
                b = pt1;
            }
            var res = new int[a.Length + b.Length];
            a.CopyTo(res, 0);
            b.CopyTo(res, a.Length);
            return res;
        }

        /// <summary><c>setnearest</c>: the group nearest to <paramref name="pos"/> among the groups that are still in the list.</summary>
        private static void SetNearest(double[][] eff, int[] next, double[] mindisfrom, int[] nearest, int pos)
        {
            double mindis = 999.9;
            int near = -1;
            for (int j = next[pos]; j != -1; j = next[j])
            {
                double tmpdouble;
                if ((tmpdouble = eff[pos][j - pos]) < mindis)
                {
                    mindis = tmpdouble;
                    near = j;
                }
            }
            for (int j = 0; j != -1 && j != pos; j = next[j])
            {
                double tmpdouble;
                if ((tmpdouble = eff[j][pos - j]) < mindis)
                {
                    mindis = tmpdouble;
                    near = j;
                }
            }
            mindisfrom[pos] = mindis;
            nearest[pos] = near;
        }

        /// <summary><c>counteff_simple_double_nostatic_memsave</c>: the weight of every sequence from the branch lengths of the tree (every branch is shared equally among the sequences below it), +0.001, normalised.</summary>
        private static void CountEffSimple(int nseq, GuideTree tree, double[] node)
        {
            var len = tree.Len;
            for (int i = 0; i < nseq; i++)
            {
                if (len[i][0] < 0.0)
                    len[i][0] = 0.0;
                if (len[i][1] < 0.0)
                    len[i][1] = 0.0;
            }
            var rootnode = new double[nseq];
            var eff = new double[nseq];
            var memhist = new int[nseq - 1][];
            for (int i = 0; i < nseq; i++)
            {
                rootnode[i] = 0.0;
                eff[i] = 1.0;
            }
            for (int i = 0; i < nseq - 1; i++)
            {
                int[] m0 = tree.Dep[i].Child0 == -1 ? new[] { tree.Topol0[i] } : memhist[tree.Dep[i].Child0];
                int[] m1 = tree.Dep[i].Child1 == -1 ? new[] { tree.Topol1[i] } : memhist[tree.Dep[i].Child1];
                var all = new int[m0.Length + m1.Length];
                m0.CopyTo(all, 0);
                m1.CopyTo(all, m0.Length);
                memhist[i] = all;
                foreach (int s1 in m0)
                {
                    rootnode[s1] += len[i][0] * eff[s1];
                    eff[s1] *= 0.5;
                }
                foreach (int s2 in m1)
                {
                    rootnode[s2] += len[i][1] * eff[s2];
                    eff[s2] *= 0.5;
                }
            }
            for (int i = 0; i < nseq; i++)
                rootnode[i] += 0.001;
            double total = 0.0;
            for (int i = 0; i < nseq; i++)
                total += rootnode[i];
            for (int i = 0; i < nseq; i++)
                node[i] = rootnode[i] / total;
        }

        // -------------------------------------------------------------------- treebase

        /// <summary><c>treebase</c> (<c>disttbfast.c</c>, all groups aligned: <c>mergeoralign[l] == 'a'</c>): joins the groups in the order of the guide tree.</summary>
        private void TreeBase(string[] aseq, GuideTree tree, double[] effarr, ref int alloclen, MafftLocalHom[][]? localhom = null)
        {
            int njob = _njob;
            var memhist = new int[njob - 1][];
            var cpmxhist = new double[njob - 1][][];
            if (localhom != null)
                CalcImportanceHalf(njob, effarr, aseq, localhom, alloclen);
            for (int l = 0; l < njob - 1; l++)
            {
                int m1 = tree.Topol0[l];
                int m2 = tree.Topol1[l];
                int child0 = tree.Dep[l].Child0, child1 = tree.Dep[l].Child1;
                int[] localmem0 = child0 == -1 ? new[] { m1 } : memhist[child0];
                int[] localmem1 = child1 == -1 ? new[] { m2 } : memhist[child1];
                int clus1 = localmem0.Length, clus2 = localmem1.Length;
                if (l != njob - 2)
                {
                    var h = new int[clus1 + clus2];
                    localmem0.CopyTo(h, 0);
                    localmem1.CopyTo(h, clus1);
                    memhist[l] = h;
                }
                // tbfast makes the matrix of every step with makedynamicmtx (disttbfast only when --amax > 0): the matrix is shifted by the (negative) offset of the distance of the step from the tips, which is
                // zero for ordinary distances and not for the negative distances that a pair library can give
                double[][] dynamicmtx = localhom != null ? MakeDynamicMtx(G.NDisConsweightMulti, tree.Dep[l].DistFromTip) : G.NDisConsweightMulti;
                int len1 = aseq[m1].Length;
                int len2 = aseq[m2].Length;
                if (G.NeverMemSave == 0 && G.Alg != 'M' && (len1 > 30000 || len2 > 30000))
                    throw new NotSupportedException("MAFFT switches to its memory-saving mode (MSalignmm) for groups longer than 30 000 columns; that mode is not available in this port.");
                if (alloclen < len1 + len2)
                {
                    Report("\nReallocating..");
                    alloclen = len1 + len2 + 1000;
                    Report("done. *alloclen = " + alloclen + "\n");
                }
                var mseq1 = new string[clus1];
                var mseq2 = new string[clus2];
                var effarr1 = new double[clus1];
                var effarr2 = new double[clus2];
                double minimumWeight = localhom != null ? G.MinimumWeight : 0.0;
                double orieff1 = FastConjunction(localmem0, aseq, mseq1, effarr1, effarr, minimumWeight);
                double orieff2 = FastConjunction(localmem1, aseq, mseq2, effarr2, effarr, minimumWeight);
                if (localhom != null)
                {
                    // tbfast, constraint == 2: the importance of the blocks of the pairwise alignments enters the match scores of A__align
                    if (l < 1000 || l % 100 == 0)
                        Report("\rSTEP " + (l + 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + " /" + (njob - 1) + " ");
                    var shrink = ShrinkLocalHom(localmem0, localmem1, localhom);
                    var imp = FillImp(clus1, clus2, mseq1[0].Length, mseq2[0].Length, mseq1, mseq2, effarr1, effarr2, shrink, localmem0, localmem1);
                    AAlign(dynamicmtx, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap,
                        cpmxhist, child0, child1, l, orieff1, orieff2, imp);
                    _nlen[m1] = (int)(0.5 * (_nlen[m1] + _nlen[m2]));
                    for (int k = 0; k < clus1; k++)
                        aseq[localmem0[k]] = mseq1[k];
                    for (int k = 0; k < clus2; k++)
                        aseq[localmem1[k]] = mseq2[k];
                    continue;
                }
                if (l < 500 || l % 100 == 0)
                    Report("\rSTEP " + (l + 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + " / " + (njob - 1) + " ");
                bool ffttry = _nlen[m1] > clus1 && _nlen[m2] > clus2;
                if (G.ForceFft || (G.UseFft && ffttry))
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" f\b\b");
                    Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                }
                else
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" d\b\b");
                    switch (G.Alg)
                    {
                        case 'A':
                            if (clus1 == 1 && clus2 == 1)
                            {
                                G__align11(dynamicmtx, mseq1, mseq2, alloclen, G.Outgap, G.Outgap);
                            }
                            else
                            {
                                if ((l < 500 || l % 100 == 0) && ((child1 != -1 && cpmxhist[child1] != null) || (child0 != -1 && cpmxhist[child0] != null)))
                                    Report(" h");
                                AAlign(dynamicmtx, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap,
                                    cpmxhist, child0, child1, l, orieff1, orieff2);
                            }
                            break;
                        default:
                            throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                    }
                }
                _nlen[m1] = (int)(0.5 * (_nlen[m1] + _nlen[m2]));
                for (int k = 0; k < clus1; k++)
                    aseq[localmem0[k]] = mseq1[k];
                for (int k = 0; k < clus2; k++)
                    aseq[localmem1[k]] = mseq2[k];
            }
            FreeStatics();
        }

        /// <summary><c>makedynamicmtx</c> with <c>dist2offset</c> (<c>mltaln9.c</c>): a copy of the matrix, the letters' scores raised by <c>offset x 600</c> where the offset is half the distance (at most 0) minus <c>--amax</c>.</summary>
        private double[][] MakeDynamicMtx(double[][] inp, double distfromtip)
        {
            double offset = (distfromtip * 2.0) * 0.5 - G.SpecificityConsideration;
            if (offset > 0.0)
                offset = 0.0;
            int n = G.Nalphabets;
            var res = new double[n][];
            for (int i = 0; i < n; i++)
            {
                res[i] = new double[n];
                Array.Copy(inp[i], res[i], n);
            }
            if (offset == 0.0)
                return res;
            for (int i = 0; i < n; i++)
            {
                if (G.Amino[i] == '-')
                    continue;
                for (int j = 0; j < n; j++)
                {
                    if (G.Amino[j] == '-')
                        continue;
                    res[i][j] = inp[i][j] + offset * 600;
                }
            }
            return res;
        }

        /// <summary>The end of <c>treebase</c> / <c>TreeDependentIteration</c>: <c>A__align( NULL ... )</c>, <c>Falign( NULL ... )</c> release the work arrays of the dynamic programs; the next call allocates (and clears) them again.</summary>
        private void FreeStatics()
        {
            _aOrlgth1 = 0;
            _aOrlgth2 = 0;
            _arAlloclen = 0;
            _impMtx = null;
        }

        /// <summary><c>fastconjuction_noname</c>: the rows of a group and their weights, normalised to 1; returns the total weight before normalising.</summary>
        private static double FastConjunction(int[] memlist, string[] seq, string[] aseq, double[] peff, double[] eff, double mineff = 0.0)
        {
            double total = 0.0;
            for (int k = 0; k < memlist.Length; k++)
            {
                int m = memlist[k];
                aseq[k] = seq[m];
                peff[k] = eff[m] < mineff ? mineff : eff[m];
                total += peff[k];
            }
            for (int m = 0; m < memlist.Length; m++)
                peff[m] /= total;
            return total;
        }

        // -------------------------------------------------------------------- Galign11.c

        /// <summary><c>G__align11</c>: the alignment of two single sequences (Gotoh with the penalties of the program, free end gaps when <paramref name="headgp"/> / <paramref name="tailgp"/> are 0); in place.</summary>
        internal double G__align11(double[][] dynamicmtx, string[] seq1, string[] seq2, int alloclen, int headgp, int tailgp)
        {
            string s1 = seq1[0], s2 = seq2[0];
            int lgth1 = s1.Length, lgth2 = s2.Length;
            if (lgth1 == 0 && lgth2 == 0)
                return 0.0;
            if (lgth1 == 0)
            {
                seq1[0] = new string('-', lgth2);
                return 0.0;
            }
            if (lgth2 == 0)
            {
                seq2[0] = new string('-', lgth1);
                return 0.0;
            }
            if (G.TryWarp)
                throw new NotSupportedException("--allowshift / --shiftpenalty (the warp of G__align11) is not available in this port.");
            var amino = new double[256][];
            for (int i = 0; i < 256; i++)
                amino[i] = new double[256];
            for (int i = 0; i < G.Nalphabets; i++)
                for (int j = 0; j < G.Nalphabets; j++)
                    amino[G.Amino[i]][G.Amino[j]] = dynamicmtx[i][j];
            double fpenalty = G.Penalty;
            double fpenaltyEx = G.PenaltyEx;
            const double termGapFac = 0.0, termGapFacEx = 0.0;
            RequireCells(((long)lgth1 + 1) * ((long)lgth2 + 1), MafftMaxPairCells, nameof(MafftMaxPairCells), "MAFFT", "keeps an int trace-back matrix ijp[len1 + 1][len2 + 1] (G__align11)",
                string.Create(CultureInfo.InvariantCulture, $"{lgth1:N0} x {lgth2:N0} residues = {((long)lgth1 + 1) * ((long)lgth2 + 1):N0} cells"), "The program has no limit and is killed when the allocation (4 bytes per cell plus the row arrays) fails", "use Kalign for sequences of this length");
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            var lastverticalw = new double[lgth1 + 2];
            var ijp = new int[lgth1 + 1][];
            for (int i = 0; i <= lgth1; i++)
                ijp[i] = new int[lgth2 + 1];
            double[] currentw = w1, previousw = w2;
            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = amino[s2[0]][s1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = amino[s1[0]][s2[j]];
            if (headgp == 1)
            {
                for (int i = 1; i < lgth1 + 1; i++)
                    initverticalw[i] += fpenalty;
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] += fpenalty;
            }
            else
            {
                for (int i = 1; i < lgth1 + 1; i++)
                {
                    initverticalw[i] += fpenalty * termGapFac;
                    initverticalw[i] += fpenaltyEx * i * termGapFacEx;
                }
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    currentw[j] += fpenalty * termGapFac;
                    currentw[j] += fpenaltyEx * j * termGapFacEx;
                }
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = currentw[lgth2 - 1];
            int lasti = tailgp != 0 ? lgth1 + 1 : lgth1;
            int lastj = lgth2 + 1;
            double mi;
            int mpi;
            double wm = 0.0;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                double[] aminoRow = amino[i < lgth1 ? s1[i] : '\0']; // the terminating NUL of the sequence scores 0
                for (int j = 0; j < lgth2; j++)
                    currentw[j] = aminoRow[s2[j]];
                currentw[0] = initverticalw[i];
                mi = previousw[0];
                mpi = 0;
                double fpenaltyExI = i < lgth1 ? fpenaltyEx : 0.0;
                var ijpRow = ijp[i];
                for (int j = 1; j < lastj; j++)
                {
                    wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    double g;
                    if ((g = mi + fpenalty) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = previousw[j - 1]) >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    mi += fpenaltyExI;
                    if ((g = m[j] + fpenalty) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    if ((g = previousw[j - 1]) >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    if (j < lgth2)
                        m[j] += fpenaltyEx;
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
            }
            var (a, b, wmo) = Atracking11(currentw, lastverticalw, s1, s2, ijp, tailgp);
            if (tailgp == 0)
                wm = wmo;
            if (alloclen < a.Length)
                throw new InvalidOperationException("LENGTH OVER!");
            seq1[0] = a;
            seq2[0] = b;
            return wm;
        }

        private static (string, string, double) Atracking11(double[] lasthorizontalw, double[] lastverticalw, string seq1, string seq2, int[][] ijp, int tailgp)
        {
            int lgth1 = seq1.Length, lgth2 = seq2.Length;
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = i + 1;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = -(j + 1);
            double wm = 0.0;
            if (tailgp != 1)
            {
                wm = lasthorizontalw[lgth2 - 1] - 1.0;
                double g;
                for (int j = lgth2 - 2; j >= 0; j--)
                {
                    if ((g = lasthorizontalw[j] + 0.0) > wm)
                    {
                        wm = g;
                        ijp[lgth1][lgth2] = -(lgth2 - j);
                    }
                }
                for (int i = lgth1 - 2; i >= 0; i--)
                {
                    if ((g = lastverticalw[i] + 0.0) > wm)
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
            var a = new char[lgth1 + lgth2 + 1];
            var b = new char[lgth1 + lgth2 + 1];
            int pa = lgth1 + lgth2, pb = lgth1 + lgth2;
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
                    a[--pa] = seq1[ifi + l];
                    b[--pb] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    a[--pa] = '-';
                    b[--pb] = seq2[jfi + l];
                    k++;
                }
                if (iin <= 0 || jin <= 0)
                    break;
                a[--pa] = seq1[ifi];
                b[--pb] = seq2[jfi];
                k++;
                iin = ifi;
                jin = jfi;
                k++;
            }
            return (new string(a, pa, lgth1 + lgth2 - pa), new string(b, pb, lgth1 + lgth2 - pb), wm);
        }
    }

    #endregion
}
