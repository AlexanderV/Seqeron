#pragma warning disable IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854, S2325, S1481, S3260
using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: nodepair (--large / --memsavetree / --stepadd with --globalpair, --localpair, --genafpair)

    /// <summary>What <c>nodepair</c> leaves in the directory <c>hat3dir</c> for <c>tbfast</c>: the guide tree (<c>hat3dir/tree</c>), the sequences whose pairs were aligned (<c>hat3dir/uselh</c>) and the ungapped blocks of the pairwise alignments of every node (<c>hat3node-n-subid</c>).</summary>
    internal sealed class NodePairData
    {
        internal int[][] Topol0 = Array.Empty<int[]>(), Topol1 = Array.Empty<int[]>();
        internal double[][] Len = Array.Empty<double[]>();
        internal SequenceAligner.MafftRun.Treedep[] Dep = Array.Empty<SequenceAligner.MafftRun.Treedep>();
        internal int[] Uselh = Array.Empty<int>();
        internal int[] FilesForNode = Array.Empty<int>();

        /// <summary>For every node the records in the order the program writes them: the positions i, j of the two sequences in the member lists of the node and the blocks (head of the list).</summary>
        internal List<(int I, int J, MafftLocalHom Head)>[] Records = Array.Empty<List<(int, int, MafftLocalHom)>>();

        internal string? TreeFile;
        internal string Log = "";
    }

    internal sealed partial class MafftGlobals
    {
        internal int LhLimit = int.MaxValue;
        internal string NodeInput = "";

        /// <summary><c>arguments()</c> of <c>nodepair.c</c>: the options of <c>pairlocalalign</c> and <c>-i file</c>, <c>-q n</c>, <c>-U</c>, <c>-t</c>, <c>-y</c>, <c>-T</c>, <c>-B</c>, <c>-F</c>; the others are "illegal option" and stop the program.</summary>
        internal static MafftGlobals ParseNodepair(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Alg = 'A', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            g.FastaThreshold = 2.7;
            g.Constraint = 2;
            int k = 0;
            int argc = argv.Count + 1;
            string Next()
            {
                argc--;
                return k < argv.Count ? argv[k++] : "";
            }
            while (--argc > 0 && k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'i': g.NodeInput = Next(); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'O': g.PpenaltyOp = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'E': g.PpenaltyExGenaf = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; rest = true; break;
                        case 'l': g.FastaThreshold = CAtof(Next()); g.Constraint = 2; rest = true; break;
                        case 'C':
                            g.Nthread = CAtoi(Next());
                            if (g.Nthread == 0)
                                g.Nthread = 1;
                            rest = true;
                            break;
                        case 'u': g.SpecificityConsideration = CAtof(Next()); rest = true; break;
                        case 'K': break;
                        case 'a': g.Fmodel = 1; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'L': g.Alg = 'L'; break;
                        case 'N': g.Alg = 'N'; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'Z': g.UseNaiveScore = true; break;
                        case 'B': break;
                        case 'T': break;
                        case 'F': g.UseFft = true; break;
                        case 'U': g.TreeIn = 1; break;
                        case 't': g.TreeOut = 1; break;
                        case 'y': g.DivPairScore = 1; break;
                        case ':': g.NWildcard = 1; break;
                        case 'q': g.LhLimit = CAtoi(Next()); rest = true; break;
                        default:
                            g.Report("illegal option " + c + "\n");
                            argc = 0;
                            rest = true;
                            k = argv.Count;
                            break;
                    }
                }
            }
            if (argc == 1)
                argc--;
            if (argc != 0)
                throw new MafftProgramFailureException(g.Log + "pairlocalalign options: Check source file !\n");
            return g;
        }

        internal int DivPairScore;
    }

    internal sealed partial class MafftRun
    {
        private const int LhBlockFactor = 2;
        private const double MinBlockLen2 = 1000000000.0;

        /// <summary>The conversion of a <c>double</c> to <c>int</c> as the program (x86-64) does it: NaN and values out of range give the smallest <c>int</c>.</summary>
        private static int CInt(double x) => double.IsNaN(x) || x >= 2147483648.0 || x <= -2147483649.0 ? int.MinValue : (int)x;

        /// <summary><c>limitlh</c> (<c>--lhlimit n</c>): the <paramref name="limit"/> longest sequences (shuffled within each half of the order) are used for the local homologies.</summary>
        private static int[] LimitLh(int[] lengths, int limit)
        {
            int size = lengths.Length;
            var uselh = new int[size];
            var items = Enumerable.Range(0, size).Select(i => (Len: lengths[i], Num: i)).OrderByDescending(x => x.Len).ToArray();
            var rand = new GlibcRandom();
            void Shuffle(int start, int count)
            {
                for (int i = 0; i < count; i++)
                {
                    int j = rand.Next() % count;
                    (items[start + i], items[start + j]) = (items[start + j], items[start + i]);
                }
            }
            Shuffle(0, size / 2);
            Shuffle(size / 2, size - size / 2);
            if (limit > size)
                limit = size;
            for (int i = 0; i < limit; i++)
                uselh[items[i].Num] = 1;
            for (int i = limit; i < size; i++)
                uselh[items[i].Num] = 0;
            return uselh;
        }

        private double DistFuncNoAlign(string s1, string s2, double selfscore1, double selfscore2)
        {
            double v;
            if (G.Alg == 'A')
            {
                v = G__align11Noalign(G.NDisConsweightMulti, s1, s2);
            }
            else if (G.Alg == 'L')
            {
                v = L__align11(G.NDisConsweightMulti, 0.0, new[] { s1 }, new[] { s2 }, _nlenmax * 2, false, out _, out _);
            }
            else
            {
                var t1 = new[] { s1 };
                var t2 = new[] { s2 };
                v = GenL__align11(G.NDisConsweightMulti, t1, t2, _nlenmax * 2, out _, out _);
                if (G.UseNaiveScore)
                    v = NaivePairScore11(t1[0], t2[0], 0);
            }
            return Score2Dist(v, selfscore1, selfscore2);
        }

        private void DistFunc(double[][] mtx, string s1, string s2, MafftLocalHom head)
        {
            var a = new[] { s1 };
            var b = new[] { s2 };
            int off1 = 0, off2 = 0;
            if (G.Alg == 'A')
                G__align11(mtx, a, b, _nlenmax * 2, 1, 1);
            else if (G.Alg == 'L')
                L__align11(mtx, 0.0, a, b, _nlenmax * 2, true, out off1, out off2);
            else
                GenL__align11(mtx, a, b, _nlenmax * 2, out off1, out off2);
            PutLocalHom2(a[0], b[0], head, off1, off2);
        }

        /// <summary>
        /// <c>nodepair</c> (<c>nodepair.c</c> with <c>pairalign_node</c> of <c>pairlocalalign.c</c>): the guide tree from the dynamic-programming distance of every sequence to the earlier ones (<c>compacttreedpdist</c>, the tree of
        /// <see cref="CompactTreeGivenDist"/>) or the tree of the user, then the pairwise alignment of every pair of sequences at the node that joins them (only pairs with a sequence of <c>uselh</c>): their blocks are what <c>tbfast</c> reads.
        /// </summary>
        internal NodePairData NodePair(IReadOnlyList<string> seq0, IReadOnlyList<string>? names)
        {
            int njob = seq0.Count;
            _njob = njob;
            if (G.Alg is not ('A' or 'L' or 'N'))
                throw new MafftProgramFailureException("alg " + G.Alg + " is not yet supported\n");
            if (G.Alg != 'N' && G.UseNaiveScore)
                throw new MafftProgramFailureException("The combination of usenaivescoreinsteadofalignmentscore and alg " + G.Alg + " is not yet supported\n");
            int lhlimit = G.LhLimit;
            int constraint = G.Constraint;
            if (G.FastaThreshold < 0.0001)
            {
                constraint = 0;
                lhlimit = 0;
            }
            _nlenmax = 0;
            var seq = new string[njob];
            for (int i = 0; i < njob; i++)
            {
                seq[i] = seq0[i];
                _nlenmax = Math.Max(_nlenmax, seq[i].Length);
            }
            if (njob < 2)
                throw new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
            Report(StackSizeText(njob));
            G.Constants(seq);
            if (CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            var fileNames = names ?? DefaultNames(njob);
            var topol0 = new int[njob][];
            var topol1 = new int[njob][];
            var len = new double[njob][];
            for (int i = 0; i < njob; i++)
                len[i] = new double[2];
            var dep = new Treedep[njob];
            bool treeout = G.TreeOut != 0;
            bool treegiven = G.TreeIn != 0;
            if (treegiven)
                LoadTree(njob, topol0, topol1, len, dep, treeout, fileNames);

            // ---- pairalign_node
            var bseq = new string[njob];
            var dseq = new string[njob];
            for (int i = 0; i < njob; i++)
            {
                bseq[i] = seq[i].Replace("-", "");
                dseq[i] = bseq[i].Replace("X", "").Replace("x", "");
            }
            var lengths = bseq.Select(x => x.Length).ToArray();
            int[] uselh = LimitLh(lengths, lhlimit);
            var selfscore = new double[njob];
            for (int i = 0; i < njob; i++)
                selfscore[i] = NaivePairScore11(seq[i], seq[i], 0);
            if (!treegiven)
            {
                // calcnearest_para
                var neighbors = new int[njob];
                var mindists = new double[njob];
                int nlim = njob - 2;
                for (int pos = njob - 1; pos >= 1; pos--)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    if ((nlim - pos) % 100 == 0)
                    {
                        int progress = CInt(((double)(nlim - pos) * nlim - (nlim - pos) * ((nlim - pos) - 1.0) * 0.5) / (nlim * (nlim - 1.0) * 0.5) * 100);
                        Report("Step " + (nlim - pos).ToString(CultureInfo.InvariantCulture) + " (" + progress.ToString(CultureInfo.InvariantCulture) + "%), thread 0       \r");
                    }
                    double mindist = 999.9;
                    int neighbor = 0;
                    for (int i = pos - 1; i > -1; i--)
                    {
                        if (dseq[pos].Length == 0 || dseq[i].Length == 0)
                            throw UndefinedBehaviour("an empty sequence: the program reads the byte before the sequence");
                        double tmpdist = DistFuncNoAlign(dseq[pos], dseq[i], selfscore[pos], selfscore[i]);
                        if (mindist > tmpdist)
                        {
                            mindist = tmpdist;
                            neighbor = i;
                        }
                    }
                    mindists[pos] = mindist;
                    neighbors[pos] = neighbor;
                }
                CompactTreeGivenDist(njob, mindists, neighbors, topol0, topol1, len, dep, treeout, fileNames);
            }
            else
            {
                Report("treegiven\n");
            }

            // ---- recalcpairs_para4
            var data = new NodePairData { Topol0 = topol0, Topol1 = topol1, Len = len, Dep = dep, Uselh = uselh, FilesForNode = new int[njob - 1], Records = new List<(int, int, MafftLocalHom)>[njob - 1], TreeFile = TreeFile };
            double lenav = 0.0;
            for (int n = 0; n < njob; n++)
                lenav += selfscore[n];
            lenav /= njob * 600;
            var memhist = new int[njob - 1][];
            var n0s = new int[njob - 1];
            var n1s = new int[njob - 1];
            var mem0s = new int[njob - 1][];
            var mem1s = new int[njob - 1][];
            var mergeoralign = new char[njob - 1];
            double sizeav = 0.0;
            for (int n = 0; n < njob - 1; n++)
            {
                int[] mem0 = dep[n].Child0 == -1 ? new[] { topol0[n][0] } : memhist[dep[n].Child0];
                int[] mem1 = dep[n].Child1 == -1 ? new[] { topol1[n][0] } : memhist[dep[n].Child1];
                memhist[n] = mem0.Concat(mem1).ToArray();
                mem0s[n] = mem0;
                mem1s[n] = mem1;
                n0s[n] = mem0.Length;
                n1s[n] = mem1.Length;
                sizeav += (double)mem0.Length * mem1.Length;
                mergeoralign[n] = 'a';
                bool i0 = mem0.Any(x => uselh[x] != 0);
                bool i1 = mem1.Any(x => uselh[x] != 0);
                if (!i0 && !i1)
                    mergeoralign[n] = 'n';
            }
            sizeav /= njob - 1;
            Report("sizeav = " + MafftGlobals.F(sizeav, 0, 6) + ", BLOCKSIZE = " + ((int)(sizeav * LhBlockFactor)).ToString(CultureInfo.InvariantCulture) + ", npairs=" + ((long)n0s[0] * n1s[0]).ToString(CultureInfo.InvariantCulture) + "-"
                + ((long)n0s[njob - 2] * n1s[njob - 2]).ToString(CultureInfo.InvariantCulture) + "\n");
            var jobs = new List<(int Node, int Start, int End, int SubId, long NPairs)>();
            int nnodesdivided = 0;
            for (int k = njob - 2; k >= 0; k--)
            {
                int n = k;
                if (mergeoralign[n] == 'n')
                    continue;
                int n0 = n0s[n], n1 = n1s[n];
                int blocksize0;
                if ((double)n0 * n1 * lenav * lenav > MinBlockLen2)
                {
                    blocksize0 = (int)(MinBlockLen2 / n1 / lenav / lenav);
                    if (blocksize0 == 0)
                        blocksize0 = 1;
                }
                else
                {
                    blocksize0 = n0;
                }
                if (n0 > blocksize0)
                    nnodesdivided++;
                int subid = 0;
                for (int b = 0; b < n0; b += blocksize0)
                {
                    int blocksize = Math.Min(blocksize0, n0 - b);
                    jobs.Add((n, b, b + blocksize, subid, (long)n1 * blocksize));
                    subid++;
                }
                data.FilesForNode[n] = subid;
            }
            Report("numjob=" + jobs.Count.ToString(CultureInfo.InvariantCulture) + ", nnodesdivided=" + nnodesdivided.ToString(CultureInfo.InvariantCulture) + ", njob=" + njob.ToString(CultureInfo.InvariantCulture) + "\n");
            Report("Divided " + nnodesdivided.ToString(CultureInfo.InvariantCulture) + " nodes to " + (jobs.Count - (njob - 1)).ToString(CultureInfo.InvariantCulture) + " jobs\n");
            for (int n = 0; n < njob - 1; n++)
                data.Records[n] = new List<(int, int, MafftLocalHom)>();
            double[][]? dynamicmtx = G.SpecificityConsideration > 0.0 ? null : G.NDisConsweightMulti;
            ulong done = 0;
            ulong totalpairs = (ulong)njob * (ulong)(njob - 1) / 2;
            for (int step = 0; step < jobs.Count; step++)
            {
                AlignmentCancellation.ThrowIfRequested();
                var job = jobs[step];
                int n = job.Node;
                done += (ulong)job.NPairs;
                double[][] mtxptr = dynamicmtx ?? MakeDynamicMtx(G.NDisConsweightMulti, dep[n].DistFromTip);
                if (step % 100 == 0)
                {
                    int progress = CInt((double)done / totalpairs * 100);
                    Report("Node " + n.ToString("D6", CultureInfo.InvariantCulture) + "-" + job.SubId.ToString("D3", CultureInfo.InvariantCulture) + " (" + progress.ToString(CultureInfo.InvariantCulture) + "%)          \r");
                }
                int[] mem0 = mem0s[n], mem1 = mem1s[n];
                for (int i = job.Start; i < job.End; i++)
                {
                    int m00 = mem0[i];
                    for (int j = 0; j < mem1.Length; j++)
                    {
                        int m11 = mem1[j];
                        if (uselh[m11] == 0 && uselh[m00] == 0)
                            continue;
                        int m0 = Math.Min(m00, m11), m1 = Math.Max(m00, m11);
                        if (bseq[m0].Length == 0 || bseq[m1].Length == 0)
                            throw UndefinedBehaviour("an empty sequence: the program reads the byte before the aligned sequence");
                        var head = new MafftLocalHom();
                        DistFunc(mtxptr, bseq[m0], bseq[m1], head);
                        // readlocalhomtable2_single_bin_noseek: the blocks as tbfast reads them back (opt and importance are the mean score per column, rescaled)
                        for (var t = head; t != null; t = t.Next)
                        {
                            t.Importance = (t.Opt + 0.00) / 5.8 * 600;
                            t.Opt = (t.Opt + 0.00) / 5.8 * 600;
                        }
                        data.Records[n].Add((i, j, head));
                    }
                }
            }
            Report("nodepair (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return data;
        }

        /// <summary>The records of node <paramref name="node"/> as the table that <see cref="FillImp"/> takes (<c>fillimp_file</c> of the program): the entry (i, j) of the two groups is the blocks of that pair.</summary>
        private static MafftLocalHom?[][] NodeShrink(NodePairData data, int node, int clus1, int clus2)
        {
            var shrink = new MafftLocalHom?[clus1][];
            for (int i = 0; i < clus1; i++)
                shrink[i] = new MafftLocalHom?[clus2];
            foreach (var (i, j, head) in data.Records[node])
                shrink[i][j] = head;
            return shrink;
        }
    }

    #endregion
}
