using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: splittbfast (--parttree, --dpparttree, --fastaparttree)

    /// <summary>The command line of <c>splittbfast</c> (<c>arguments()</c> of <c>splittbfast.c</c>, the options that the script can produce).</summary>
    internal sealed partial class MafftGlobals
    {
        /// <summary><c>picksize</c>, <c>classsize</c>: the number of sequences that represent a group (<c>-p</c>), the size of a group that is aligned as a whole (<c>-s</c>).</summary>
        internal int PickSize = NotSpecified, ClassSize = NotSpecified;
        /// <summary><c>doalign</c>: 0 = 6-mer distance, 1 = dynamic programming (<c>-U</c>), 'f' = fasta34 (<c>-S</c>).</summary>
        internal int DoAlign;
        /// <summary><c>fromaln</c> (<c>-Z</c>): the distances come from the alignment that the previous pass made.</summary>
        internal bool FromAln;
        /// <summary><c>reorder</c> (<c>-x</c> switches it off).</summary>
        internal bool Reorder = true;
        /// <summary><c>treeout</c> (<c>-t</c>).</summary>
        internal bool SplitTreeOut;
        internal double Tokyoripara = NotSpecified;
        internal string? InputFile;

        internal static MafftGlobals ParseSplittbfast(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = -1530, PpenaltyEx = NotSpecified, Poffset = -123, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Alg = 'A', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : throw new ArgumentException("A MAFFT option is missing its argument.");
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'p': g.PickSize = CAtoi(Next()); g.Report("picksize = " + g.PickSize.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 's': g.ClassSize = CAtoi(Next()); g.Report("groupsize = " + g.ClassSize.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'i': g.InputFile = Next(); g.Report("inputfile = " + g.InputFile + "\n"); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); g.Report("ppenalty_ex = " + g.PpenaltyEx.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); g.Report("kimuraR = " + g.KimuraR.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; g.Report("jtt " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; g.Report("tm " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'T': g.Tokyoripara = CAtof(Next()); rest = true; break;
                        case 'n': g.OutNumber = 1; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'S': g.DoAlign = 'f'; break;
                        case 'Z': g.FromAln = true; break;
                        case 'U': g.DoAlign = 1; break;
                        case 'x': g.Reorder = false; break;
                        case 't': g.SplitTreeOut = true; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'L': g.LegacyGapCost = 1; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'N': g.NeverMemSave = 1; break;
                        case ':': g.NWildcard = 1; break;
                        case 'X': g.TreeMethod = 'X'; g.SueffGlobal = CAtof(Next()); g.Report("sueff_global = " + g.SueffGlobal.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'E': g.TreeMethod = 'E'; break;
                        case 'q': g.TreeMethod = 'q'; break;
                        default: throw new NotSupportedException($"The MAFFT option -{c} of splittbfast is not available in this port.");
                    }
                }
            }
            return g;
        }
    }

    internal sealed partial class MafftRun
    {
        // -------------------------------------------------------------------- splittbfast.c

        /// <summary>The result of one run of <c>splittbfast</c>: the rows in the order of the output, the input index of every row, the guide tree file (<c>-t</c>), the label of the group of every row.</summary>
        internal sealed record SplitResult(string[] Rows, int[] Order, string? Tree, string[] GroupLabels);

        private struct SplitScore
        {
            internal double Score;
            internal int SelfScore;
            internal int OriLen;
            internal int[] Pointt;
            internal int NumInSeq;
            internal int Shimon;
        }

        private const double SplitTokyoriPara = 0.70;
        private const double SplitMax6Dist = 10.0;

        private int _spNjob, _spNlenmax, _spAlloclen, _spClassSize, _spPickSize, _spBranchId, _spGroupId, _spTsize, _spOrderPos;
        private string[] _spSeq = Array.Empty<string>();
        private string[]? _spOriAln;
        private int[] _spNlen = Array.Empty<int>(), _spOrder = Array.Empty<int>(), _spWhichGroup = Array.Empty<int>();
        private short[]? _spMemo;
        private double _spLenfaca, _spLenfacb, _spLenfacc, _spLenfacd;
        private int _spDepth;
        private GlibcRandom _spRandom = new();
        private bool _spTreeOut;

        /// <summary>
        /// <c>main</c> of <c>splittbfast.c</c>: the sequences are classified by the distance to a few representatives (the 6-mer distance, the dynamic programming score of <c>-U</c>, the score of <c>fasta34</c> of <c>-S</c>, or the distance of the
        /// alignment of the previous pass of <c>-Z</c>), the guide tree of the representatives (UPGMA) is followed recursively and the groups are aligned in the order of the tree.
        /// <paramref name="input"/> are the sequences as the file <c>infile</c> / <c>pre</c> holds them (aligned for <c>-Z</c>).
        /// </summary>
        internal SplitResult Splittbfast(IReadOnlyList<string> input, IReadOnlyList<string>? names = null)
        {
            int njob = input.Count;
            _spNjob = njob;
            if (njob < 2)
                throw new ArgumentException("At least 2 sequences should be input!");
            _spNlenmax = 0;
            foreach (var s in input)
                _spNlenmax = Math.Max(_spNlenmax, s.Length);
            if (G.Dorp == '\0')
            {
                int atgc = 0, total = 0;
                foreach (var s in input)
                {
                    if (total >= 1000000)
                        break;
                    int n = 0, a = 0;
                    foreach (char ch in s)
                    {
                        if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')))
                            continue;
                        n++;
                        char c = char.ToLowerInvariant(ch);
                        if (c is 'a' or 't' or 'g' or 'c' or 'u' or 'n')
                            a++;
                    }
                    atgc += a;
                    total += n;
                }
                G.Dorp = (double)atgc / total > 0.75 ? 'd' : 'p';
            }
            _spPickSize = G.PickSize;
            _spClassSize = G.ClassSize;
            if (_spPickSize == MafftGlobals.NotSpecified || _spPickSize < 2)
                _spPickSize = 50;
            if (_spClassSize == MafftGlobals.NotSpecified || _spClassSize < 0)
                _spClassSize = njob + 1;
            int doalign = G.DoAlign;
            if (G.FromAln)
                doalign = 1;
            double tokyoripara = G.Tokyoripara;
            if (tokyoripara == MafftGlobals.NotSpecified)
                tokyoripara = SplitTokyoriPara;
            if (_spPickSize > njob)
                tokyoripara = 0.0;
            _spAlloclen = _spNlenmax * 2;
            Report("alloclen = " + _spAlloclen.ToString(CultureInfo.InvariantCulture) + " in main\n");
            _spWhichGroup = new int[njob];
            _spOrder = new int[njob + 1];
            for (int i = 0; i < njob; i++)
                _spOrder[i] = -1;
            _spSeq = new string[njob];
            for (int i = 0; i < njob; i++)
                _spSeq[i] = G.Dorp == 'd' ? input[i].ToLowerInvariant() : input[i];
            if (G.FromAln)
            {
                _spOriAln = new string[njob];
                for (int i = 0; i < njob; i++)
                {
                    if (_spSeq[i].Length != _spNlenmax)
                        throw new ArgumentException("Input sequences must be aligned");
                    _spOriAln[i] = _spSeq[i];
                }
            }
            G.Constants(_spSeq);
            _spTsize = G.Dorp == 'd' ? 4096 : 46656;
            if (G.Dorp == 'd')
            {
                _spLenfaca = 0.01; _spLenfacb = 2500; _spLenfacc = 2500; _spLenfacd = 0.1;
            }
            else
            {
                _spLenfaca = 0.01; _spLenfacb = 10000; _spLenfacc = 10000; _spLenfacd = 0.1;
            }
            _spNlen = new int[njob];
            var scores = new SplitScore[njob];
            int[] mult;
            if (G.Dorp == 'd')
                mult = new[] { 1024, 256, 64, 16, 4, 1 };
            else
                mult = new[] { 7776, 1296, 216, 36, 6, 1 };
            int limit = G.Dorp == 'd' ? 4 : 6;
            int nunknown = 0;
            for (int i = 0; i < njob; i++)
            {
                string tmpseq = _spSeq[i].Replace("-", "");
                _spNlen[i] = tmpseq.Length;
                _spSeq[i] = tmpseq;
                if (_spNlen[i] < 6)
                    throw new MafftProgramFailureException("Seq " + (i + 1).ToString(CultureInfo.InvariantCulture) + ", too short, " + _spNlen[i].ToString(CultureInfo.InvariantCulture) + " characters\nname = " + (names?[i] ?? "=seq" + (i + 1).ToString(CultureInfo.InvariantCulture)) + "\nseq = " + _spSeq[i] + "\n");
                var grp = new List<int>(tmpseq.Length);
                foreach (char c in tmpseq)
                {
                    int tmp = G.AminoGrp[c];
                    if (tmp < limit)
                        grp.Add(tmp);
                    else
                        nunknown++;
                }
                if (grp.Count < 6)
                    throw new MafftProgramFailureException("Seq " + (i + 1).ToString(CultureInfo.InvariantCulture) + ", too short.\nname = " + (names?[i] ?? "=seq" + (i + 1).ToString(CultureInfo.InvariantCulture)) + "\nseq = " + _spSeq[i] + "\n");
                var core = MakePointTable(grp.ToArray(), grp.Count, 6, mult, limit);
                var pointt = new int[_spNlen[i] + 1];
                Array.Copy(core, pointt, core.Length);
                pointt[core.Length] = -1;
                scores[i].Pointt = pointt;
            }
            if (nunknown != 0)
                Report("\nThere are " + nunknown.ToString(CultureInfo.InvariantCulture) + " ambiguous characters\n");
            // WriteOptions( trap_g ): only the gap penalty line goes to the standard error
            Report("Gap Penalty = " + MafftGlobals.F((double)G.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.Poffset / 1000, 5, 2, true) + "\n");
            if (CheckAlphabet(_spSeq, "Illeagal") is { } unknownLetter)
                throw unknownLetter;
            for (int i = 0; i < njob; i++)
            {
                scores[i].NumInSeq = i;
                scores[i].OriLen = _spSeq[i].Length;
            }
            int[,]? tmpaminodis = null;
            if (doalign == 'f')
            {
                Report("fastapath=fasta34\n");
                tmpaminodis = MafftFasta34.ScoreMatrix(G.Dorp == 'd', Report);
            }
            _spMemo = new short[_spTsize];
            for (int i = 0; i < njob; i++)
            {
                int[] pt = scores[i].Pointt;
                scores[i].Shimon = pt[0] + pt[1] + pt[scores[i].OriLen - 6];
                if (doalign != 0)
                {
                    Report("\r " + i.ToString("D5", CultureInfo.InvariantCulture) + "/" + njob.ToString("D5", CultureInfo.InvariantCulture) + "   ");
                    scores[i].Pointt = Array.Empty<int>();
                    if (doalign == 'f')
                    {
                        int pscore = 0;
                        var sb = new StringBuilder(_spSeq[i]);
                        if (G.Dorp == 'd')
                        {
                            int st = 1, en = 0;
                            for (int q = 0; q < sb.Length; q++)
                            {
                                if (sb[q] == 'u')
                                    sb[q] = 't';
                                int aan = G.AminoN[sb[q]];
                                if (aan < 0 || aan >= 4)
                                    sb[q] = 'n';
                                if (sb[q] == 'n')
                                {
                                    en++;
                                    if (st != 0)
                                        continue;
                                    pscore += tmpaminodis![sb[q], sb[q]];
                                }
                                else
                                {
                                    st = 0;
                                    en = 0;
                                    pscore += tmpaminodis![sb[q], sb[q]];
                                }
                            }
                            scores[i].SelfScore = pscore - en * tmpaminodis!['n', 'n'];
                        }
                        else
                        {
                            int st = 1, en = 0;
                            for (int q = 0; q < sb.Length; q++)
                            {
                                int aan = G.AminoN[sb[q]];
                                if (aan < 0 || aan >= 20)
                                    sb[q] = 'X';
                                if (sb[q] == 'X')
                                {
                                    en++;
                                    if (st != 0)
                                        continue;
                                    pscore += tmpaminodis![sb[q], sb[q]];
                                }
                                else
                                {
                                    st = 0;
                                    en = 0;
                                    pscore += tmpaminodis![sb[q], sb[q]];
                                }
                            }
                            scores[i].SelfScore = pscore - en * tmpaminodis!['X', 'X'];
                        }
                        _spSeq[i] = sb.ToString();
                    }
                    else
                    {
                        int pscore = 0;
                        foreach (char c in _spSeq[i])
                            pscore += G.AminoDis[c][c];
                        scores[i].SelfScore = pscore;
                    }
                }
                else
                {
                    var table1 = new short[_spTsize];
                    foreach (int p in pt)
                    {
                        if (p == -1)
                            break;
                        table1[p]++;
                    }
                    scores[i].SelfScore = LocalCommonSextet(table1, pt);
                }
            }
            _spOrderPos = 0;
            _spTreeOut = G.SplitTreeOut;
            _spDepth = 0;
            _spBranchId = 0;
            _spGroupId = 0;
            _spRandom = new GlibcRandom();
            string? tree = null;
            SplitNode(scores, njob, 0, out tree, doalign, tokyoripara, -1, tmpaminodis);
            Report("\nDone.\n\n");
            string? treeFile = null;
            if (_spTreeOut)
                treeFile = tree + "\n";
            // the groups: the sequences that belong to a group of the tree that was aligned get the label Group-n, the others Group-para
            int groupnum = 0, groupid = -1;
            var labels = new string[njob];
            var paramem = new List<int>();
            for (int i = 0; i < njob; i++)
            {
                int pos = _spOrder[i];
                if (_spWhichGroup[pos] != groupid)
                {
                    groupnum++;
                    groupid = _spWhichGroup[pos];
                }
                if (_spWhichGroup[pos] != 0)
                {
                    FlushParalogs(paramem);
                    labels[i] = "Group-" + groupnum.ToString(CultureInfo.InvariantCulture);
                }
                else
                {
                    paramem.Add(pos);
                    labels[i] = "Group-para";
                }
            }
            FlushParalogs(paramem);
            var rows = new string[njob];
            var order = new int[njob];
            for (int i = 0; i < njob; i++)
            {
                int k = G.Reorder ? _spOrder[i] : i;
                rows[i] = _spSeq[k];
                order[i] = k;
            }
            string[] outLabels = labels;
            if (!G.Reorder)
            {
                // writeData_pointer: the names are those of the input order; the label of a sequence is the one of its place in the order of the tree
                outLabels = new string[njob];
                for (int i = 0; i < njob; i++)
                    outLabels[_spOrder[i]] = labels[i];
            }
            if (!(_spClassSize > 1 && _spClassSize <= njob))
                outLabels = new string[njob]; // the names are only changed with --groupsize 2 ... n
            Report("\n\n");
            Report("----------------------------------------------------------------------------\n");
            Report("\n");
            if (_spClassSize == 1)
            {
                Report("nseq = " + njob.ToString(CultureInfo.InvariantCulture) + "\n");
                Report("groupsize = " + _spClassSize.ToString(CultureInfo.InvariantCulture) + ", picksize=" + _spPickSize.ToString(CultureInfo.InvariantCulture) + "\n");
                Report("The input sequences have been sorted so that similar sequences are close.\n");
                if (G.Reorder)
                    Report("The order of sequences has been changed according to estimated similarity.\n");
                if (_spTreeOut)
                    Report("\nA guide tree is in the 'infile.tree' file.\n");
                if (doalign == 0)
                    Report("\nmafft --dpparttree might give a better result, although slow.\nmafft --fastaparttree is also available if you have fasta34.\n");
            }
            else if (groupnum > 1)
            {
                Report("groupsize = " + _spClassSize.ToString(CultureInfo.InvariantCulture) + ", picksize=" + _spPickSize.ToString(CultureInfo.InvariantCulture) + "\n");
                Report("The input sequences have been classified into " + groupnum.ToString(CultureInfo.InvariantCulture) + " groups + some paraphyletic groups\n");
                Report("Note that the alignment is not completed.\n");
                if (G.Reorder)
                    Report("The order of sequences has been changed according to estimated similarity.\n");
                if (_spTreeOut)
                    Report("\nA guide tree is in the 'infile.tree' file.\n");
                if (doalign == 0)
                    Report("\nmafft --dpparttree might give a better result, although slow.\nmafft --fastaparttree is also available if you have fasta34.\n");
            }
            else
            {
                Report("nseq = " + njob.ToString(CultureInfo.InvariantCulture) + "\n");
                Report("groupsize = " + _spClassSize.ToString(CultureInfo.InvariantCulture) + ", partsize=" + _spPickSize.ToString(CultureInfo.InvariantCulture) + "\n");
                if (G.Reorder)
                    Report("The order of sequences has been changed according to estimated similarity.\n");
                if (_spTreeOut)
                    Report("\nA guide tree is in the 'infile.tree' file.\n");
                if (doalign == 0 || G.FromAln)
                    Report("\nmafft --dpparttree might give a better result, although slow.\nmafft --fastaparttree is also available if you have fasta34.\n");
            }
            Report("\n");
            Report("----------------------------------------------------------------------------\n");
            if (doalign == 'f')
                Report("rm -f /tmp/data-<pid>*  /tmp/query-<pid> /tmp/fasta-<pid>\n");
            Report("splittbfast (" + (G.Dorp == 'd' ? "nuc" : "aa") + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new SplitResult(rows, order, treeFile, outLabels);

            void FlushParalogs(List<int> mem)
            {
                if (mem.Count > 0)
                {
                    if (mem.Count > 1 && _spClassSize > 2)
                    {
                        mem.Sort();
                        AlignParaphiles(mem.ToArray());
                    }
                    mem.Clear();
                }
            }
        }

        /// <summary><c>alignparaphiles</c>: the sequences that belong to no group are added one by one to the growing alignment.</summary>
        private void AlignParaphiles(int[] members)
        {
            int nmem = members.Length;
            for (int i = 0; i < nmem - 1; i++)
            {
                var mem1 = new int[i + 1];
                Array.Copy(members, mem1, i + 1);
                SplitPairAlign(mem1, new[] { members[i + 1] });
            }
        }

        /// <summary><c>localcommonsextet_p</c>.</summary>
        private int LocalCommonSextet(short[] table, int[] pointt)
        {
            var memo = _spMemo!;
            int value = 0;
            var touched = new List<int>();
            foreach (int point in pointt)
            {
                if (point == -1)
                    break;
                short tmp = memo[point]++;
                if (tmp < table[point])
                    value++;
                if (tmp == 0)
                    touched.Add(point);
            }
            foreach (int p in touched)
                memo[p] = 0;
            return value;
        }

        private double LenFac(double a, double b)
        {
            double longer = a > b ? a : b;
            double shorter = a > b ? b : a;
            return 1.0 / (shorter / longer * _spLenfacd + _spLenfacb / (longer + _spLenfacc) + _spLenfaca);
        }

        /// <summary><c>pairalign</c>: the two groups are aligned as profiles (<c>A__align</c>, <c>G__align11</c> for two single sequences); the rows of both are replaced in <c>seq</c>.</summary>
        private void SplitPairAlign(int[] mem1, int[] mem2)
        {
            int m1 = mem1[0], m2 = mem2[0];
            int len1 = _spSeq[m1].Length, len2 = _spSeq[m2].Length;
            if (_spAlloclen < len1 + len2)
            {
                Report("\nReallocating..");
                _spAlloclen = len1 + len2 + 1000;
                Report("done. *alloclen = " + _spAlloclen.ToString(CultureInfo.InvariantCulture) + "\n");
            }
            int clus1 = mem1.Length, clus2 = mem2.Length;
            var mseq1 = new string[clus1];
            var mseq2 = new string[clus2];
            var eff1 = new double[clus1];
            var eff2 = new double[clus2];
            for (int k = 0; k < clus1; k++)
            {
                mseq1[k] = _spSeq[mem1[k]];
                eff1[k] = 1.0;
            }
            for (int k = 0; k < clus1; k++)
                eff1[k] /= clus1;
            for (int k = 0; k < clus2; k++)
            {
                mseq2[k] = _spSeq[mem2[k]];
                eff2[k] = 1.0;
            }
            for (int k = 0; k < clus2; k++)
                eff2[k] /= clus2;
            if (G.NeverMemSave == 0 && G.Alg != 'M' && (len1 > 10000 || len2 > 10000))
                throw new NotSupportedException("MAFFT switches to its memory-saving mode (MSalignmm) for groups longer than 10 000 columns; that mode is not available in this port.");
            Report("\bd");
            switch (G.Alg)
            {
                case 'A':
                    if (clus1 == 1 && clus2 == 1)
                        G__align11(G.NDisConsweightMulti, mseq1, mseq2, _spAlloclen, G.Outgap, G.Outgap);
                    else
                        AAlign(G.NDisConsweightMulti, G.Penalty, G.PenaltyEx, mseq1, mseq2, eff1, eff2, clus1, clus2, _spAlloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, 0.0, 0.0);
                    break;
                default:
                    throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
            }
            _spNlen[m1] = (int)(0.5 * (_spNlen[m1] + _spNlen[m2]));
            for (int k = 0; k < clus1; k++)
                _spSeq[mem1[k]] = mseq1[k];
            for (int k = 0; k < clus2; k++)
                _spSeq[mem2[k]] = mseq2[k];
        }

        /// <summary><c>G__align11_noalign( n_disLN, -1200, -60, ... )</c>: the score of the global alignment of two sequences with the penalties of <c>splittbfast</c>, without the alignment.</summary>
        private double G__align11NoalignLN(string s1, string s2)
        {
            const double fpenalty = -1200.0, fpenaltyEx = -60.0;
            int lgth1 = s1.Length, lgth2 = s2.Length;
            var amino = new double[256][];
            for (int i = 0; i < 256; i++)
                amino[i] = new double[256];
            for (int i = 0; i < G.Nalphabets; i++)
                for (int j = 0; j < G.Nalphabets; j++)
                    amino[G.Amino[i]][G.Amino[j]] = G.NDisLN[i][j];
            var m = new double[lgth2 + 2];
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            double[] currentw = w1, previousw = w2;
            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = amino[s2[0]][s1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = amino[s1[0]][s2[j]];
            for (int i = 1; i < lgth1 + 1; i++)
                initverticalw[i] += fpenalty;
            for (int j = 1; j < lgth2 + 1; j++)
                currentw[j] += fpenalty;
            for (int j = 1; j < lgth2 + 1; ++j)
                m[j] = currentw[j - 1];
            int lasti = lgth1 + 1;
            double wm = 0.0;
            for (int i = 1; i < lasti; i++)
            {
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                double[] row = amino[i < lgth1 ? s1[i] : '\0'];
                for (int j = 0; j < lgth2; j++)
                    currentw[j] = row[s2[j]];
                currentw[0] = initverticalw[i];
                double mi = previousw[0];
                double fpenaltyExI = i < lgth1 ? fpenaltyEx : 0.0;
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    wm = previousw[j - 1];
                    double g;
                    if ((g = mi + fpenalty) > wm)
                        wm = g;
                    if ((g = previousw[j - 1]) > mi)
                        mi = g;
                    mi += fpenaltyExI;
                    if ((g = m[j] + fpenalty) > wm)
                        wm = g;
                    if ((g = previousw[j - 1]) > m[j])
                        m[j] = g;
                    if (j < lgth2)
                        m[j] += fpenaltyEx;
                    currentw[j] += wm;
                }
            }
            return wm;
        }

        private static List<T> StableSorted<T>(T[] arr, int n, Comparison<T> cmp)
        {
            var idx = Enumerable.Range(0, n).ToArray();
            var sorted = idx.OrderBy(i => i, Comparer<int>.Create((a, b) =>
            {
                int c = cmp(arr[a], arr[b]);
                return c;
            })).ToList();
            var res = new List<T>(n);
            foreach (int i in sorted)
                res.Add(arr[i]);
            return res;
        }

        private static int DCompare(SplitScore a, SplitScore b)
        {
            if (a.Score > b.Score)
                return 1;
            if (a.Score < b.Score)
                return -1;
            if (a.SelfScore < b.SelfScore)
                return 1;
            if (a.SelfScore > b.SelfScore)
                return -1;
            if (a.OriLen < b.OriLen)
                return 1;
            if (a.OriLen > b.OriLen)
                return -1;
            return 0;
        }

        private bool SameSequence(SplitScore a, SplitScore b) => a.Shimon == b.Shimon && string.Equals(_spSeq[a.NumInSeq], _spSeq[b.NumInSeq], StringComparison.Ordinal);

        /// <summary><c>splitseq_mq</c>: the recursion. <paramref name="scores"/> is changed in place (the sequences of the group, sorted by the distance to its representative).</summary>
        private int SplitNode(SplitScore[] scores, int nin, int uniform, out string? tree, int doalign, double tokyoripara, int qinoya, int[,]? tmpaminodis)
        {
            tree = null;
            if (nin == 0)
            {
                if (_spTreeOut)
                    tree = "";
                return 1;
            }
            if (nin < 2 || uniform == -1)
            {
                Report("\nLeaf  " + (++_spBranchId).ToString(CultureInfo.InvariantCulture) + " / " + _spNjob.ToString(CultureInfo.InvariantCulture) + "                ");
                if (_spTreeOut)
                {
                    var sb = new StringBuilder("\n");
                    if (nin > 1)
                        sb.Append('(');
                    for (int j = 0; j < nin - 1; j++)
                        sb.Append(scores[j].NumInSeq + 1).Append(',');
                    sb.Append(scores[nin - 1].NumInSeq + 1);
                    if (nin > 1)
                        sb.Append(")\n");
                    else
                        sb.Append('\n');
                    tree = sb.ToString();
                }
                for (int j = 0; j < nin; j++)
                    _spOrder[_spOrderPos++] = scores[j].NumInSeq;
                return 1;
            }
            // the representative: the sequence with the highest self score
            int selfscore0 = scores[0].SelfScore;
            int belongto = 0;
            for (int i = 0; i < nin; i++)
            {
                if (scores[i].SelfScore > selfscore0)
                {
                    selfscore0 = scores[i].SelfScore;
                    belongto = i;
                }
            }
            if (belongto != 0)
                (scores[0], scores[belongto]) = (scores[belongto], scores[0]);

            double[]? fastaRes = null;
            short[]? table1 = null;
            if (qinoya != scores[0].NumInSeq)
            {
                var lib = new List<string>(nin);
                if (doalign != 0)
                {
                    if (doalign == 'f')
                    {
                        for (int i = 0; i < nin; i++)
                            lib.Add(_spSeq[scores[i].NumInSeq]);
                        fastaRes = MafftFasta34.Search(_spSeq[scores[0].NumInSeq], lib, G.Dorp == 'd');
                        if (scores[0].SelfScore != (int)fastaRes[0])
                            Report("\n\nWARNING1: selfscore\nscores->numinseq = " + (scores[0].NumInSeq + 1) + "\nscores->orilen = " + scores[0].OriLen + "\nscores->selfscore = " + scores[0].SelfScore + ", but blastresults[0] = " + fastaRes[0].ToString("F6", CultureInfo.InvariantCulture) + "\n");
                    }
                }
                else
                {
                    table1 = new short[_spTsize];
                    foreach (int p in scores[0].Pointt)
                    {
                        if (p == -1)
                            break;
                        table1[p]++;
                    }
                }
                selfscore0 = scores[0].SelfScore;
                for (int i = 0; i < nin; i++)
                {
                    double lenfac = LenFac(scores[0].OriLen, scores[i].OriLen);
                    if (doalign != 0)
                    {
                        if (doalign == 'f')
                        {
                            scores[i].Score = (1.0 - fastaRes![i] / Math.Min(scores[0].SelfScore, scores[i].SelfScore)) * 1;
                            if (scores[i].Score < 0.0)
                                scores[i].Score = 0.0;
                        }
                        else if (G.FromAln)
                        {
                            scores[i].Score = (1.0 - NaivePairScore11(_spOriAln![scores[i].NumInSeq], _spOriAln[scores[0].NumInSeq], G.Penalty) / Math.Min(selfscore0, scores[i].SelfScore)) * 1;
                        }
                        else
                        {
                            if (_spDepth == 0)
                                Report("\r" + i + " / " + nin + "   ");
                            scores[i].Score = (1.0 - G__align11NoalignLN(_spSeq[scores[0].NumInSeq], _spSeq[scores[i].NumInSeq]) / Math.Min(selfscore0, scores[i].SelfScore)) * 1;
                        }
                    }
                    else
                    {
                        scores[i].Score = (1.0 - (double)LocalCommonSextet(table1!, scores[i].Pointt) / Math.Min(selfscore0, scores[i].SelfScore)) * lenfac;
                        if (scores[i].Score > SplitMax6Dist)
                            scores[i].Score = SplitMax6Dist;
                    }
                }
            }
            var sortedAll = StableSorted(scores, nin, DCompare);
            for (int i = 0; i < nin; i++)
                scores[i] = sortedAll[i];
            double maxdist = scores[nin - 1].Score;
            if (G.FromAln && SameSequence(scores[0], scores[nin - 1]))
                maxdist = 0.0;

            // the candidates for the representatives of the classes
            var pickkouho = new int[nin + 1];
            for (int i = 0; i < nin; i++)
                pickkouho[i] = i + 1;
            int nkouho = nin - 1;
            var picks = new List<int> { 0 };
            if (nkouho > 0)
            {
                int picktmp = pickkouho[nkouho - 1];
                nkouho--;
                if (!SameSequence(scores[picktmp], scores[0]))
                    picks.Add(picktmp);
            }
            bool first = true;
            while (picks.Count < _spPickSize && nkouho > 0)
            {
                int rn;
                if (first)
                {
                    first = false;
                    rn = (int)(nkouho * 0.5);
                }
                else
                {
                    rn = (int)((1.0 / 2147483648.0) * _spRandom.Next() * nkouho);
                }
                int picktmp = pickkouho[rn];
                nkouho--;
                pickkouho[rn] = pickkouho[nkouho];
                int j;
                for (j = 0; j < picks.Count; j++)
                    if (SameSequence(scores[picktmp], scores[picks[j]]))
                        break;
                if (j == picks.Count)
                    picks.Add(picktmp);
            }
            picks.Sort();
            int npick = picks.Count;
            var picksArr = picks.ToArray();
            var sPMap = new int[nin + 1];
            for (int i = 0; i < nin; i++)
                sPMap[i] = -1;
            for (int i = 0; i < npick; i++)
                sPMap[picksArr[i]] = i;
            var pickmtx = new double[npick][];
            for (int i = 0; i < npick; i++)
                pickmtx[i] = new double[npick - i];
            for (int j = 0; j < nin; j++)
                if (sPMap[j] != -1)
                    pickmtx[0][sPMap[j]] = scores[j].Score;
            for (int j = 1; j < npick; j++)
            {
                double[]? res = null;
                short[]? t1 = null;
                string? q1 = null;
                if (doalign != 0)
                {
                    if (doalign == 'f')
                    {
                        var lib = new List<string>(npick);
                        for (int i = 0; i < npick; i++)
                            lib.Add(_spSeq[scores[picksArr[i]].NumInSeq]);
                        res = MafftFasta34.Search(_spSeq[scores[picksArr[j]].NumInSeq], lib, G.Dorp == 'd');
                        if (scores[picksArr[j]].SelfScore != (int)res[j])
                            Report("\n\nWARNING2: selfscore j=" + j + "/" + npick + "\nscores[picks[j]].numinseq = " + (scores[picksArr[j]].NumInSeq + 1) + "\nscores[picks[j]].orilen = " + scores[picksArr[j]].OriLen + "\nscores[picks[j]].selfscore = " + scores[picksArr[j]].SelfScore + ", but blastresults[j] = " + res[j].ToString("F6", CultureInfo.InvariantCulture) + "\n");
                    }
                    else
                    {
                        q1 = _spSeq[scores[picksArr[j]].NumInSeq];
                    }
                }
                else
                {
                    t1 = new short[_spTsize];
                    foreach (int p in scores[picksArr[j]].Pointt)
                    {
                        if (p == -1)
                            break;
                        t1[p]++;
                    }
                }
                int sj = scores[picksArr[j]].SelfScore;
                pickmtx[j][0] = 0.0;
                for (int i = j + 1; i < npick; i++)
                {
                    double lenfac = LenFac(scores[picksArr[j]].OriLen, scores[picksArr[i]].OriLen);
                    int si = scores[picksArr[i]].SelfScore;
                    if (doalign != 0)
                    {
                        if (doalign == 'f')
                        {
                            pickmtx[j][i - j] = (1.0 - res![i] / Math.Min(sj, si)) * 1;
                            if (pickmtx[j][i - j] < 0.0)
                                pickmtx[j][i - j] = 0.0;
                        }
                        else if (G.FromAln)
                        {
                            Report(j + "-" + i + "/" + npick + "\r");
                            pickmtx[j][i - j] = (1.0 - NaivePairScore11(_spOriAln![scores[picksArr[i]].NumInSeq], _spOriAln[scores[picksArr[j]].NumInSeq], G.Penalty) / Math.Min(sj, si)) * 1;
                        }
                        else
                        {
                            pickmtx[j][i - j] = (1.0 - G__align11NoalignLN(q1!, _spSeq[scores[picksArr[i]].NumInSeq]) / Math.Min(sj, si)) * 1;
                        }
                    }
                    else
                    {
                        pickmtx[j][i - j] = (1.0 - (double)LocalCommonSextet(t1!, scores[picksArr[i]].Pointt) / Math.Min(sj, si)) * lenfac;
                        if (pickmtx[j][i - j] > SplitMax6Dist)
                            pickmtx[j][i - j] = SplitMax6Dist;
                    }
                }
            }
            var dfromcp = new double[npick][];
            for (int i = 0; i < npick; i++)
            {
                dfromcp[i] = new double[nin];
                for (int j = 0; j < nin; j++)
                    dfromcp[i][j] = -0.5;
            }
            for (int j = 0; j < nin; j++)
                dfromcp[0][j] = scores[j].Score;
            for (int i = 0; i < npick; i++)
                for (int j = i; j < npick; j++)
                    dfromcp[i][picksArr[j]] = dfromcp[j][picksArr[i]] = pickmtx[i][j - i];

            var tsukau = new bool[npick];
            for (int i = 0; i < npick; i++)
                tsukau[i] = true;
            if (npick > 2)
            {
                for (int i = 0; i < npick - 1; i++)
                {
                    if (!tsukau[i])
                        continue;
                    for (int j = i + 1; j < npick; j++)
                    {
                        double kijun = maxdist * tokyoripara;
                        if (!tsukau[j])
                            continue;
                        if (pickmtx[i][j - i] < kijun)
                            tsukau[j] = false;
                    }
                }
            }
            var dfromc = new List<double[]>();
            var yukos = new List<int>();
            var keep = new List<int>();
            for (int i = 0; i < npick; i++)
            {
                if (tsukau[i])
                {
                    dfromc.Add(dfromcp[i]);
                    yukos.Add(picksArr[i]);
                    keep.Add(i);
                }
            }
            int nyuko = yukos.Count;
            var yukomtx = new double[nyuko][];
            for (int ii = 0; ii < nyuko; ii++)
            {
                yukomtx[ii] = new double[nyuko - ii];
                for (int jj = ii; jj < nyuko; jj++)
                    yukomtx[ii][jj - ii] = pickmtx[keep[ii]][keep[jj] - keep[ii]];
            }
            var sYMap = new int[nin + 1];
            for (int i = 0; i < nin; i++)
                sYMap[i] = -1;
            var yOMap = new int[nyuko + 1];
            for (int i = 0; i < nyuko; i++)
            {
                sYMap[yukos[i]] = i;
                yOMap[i] = scores[yukos[i]].NumInSeq;
            }
            string[]? children = _spTreeOut ? new string[nyuko + 1] : null;
            var outs = new List<SplitScore>[nyuko];
            for (int i = 0; i < nyuko; i++)
                outs[i] = new List<SplitScore>();
            var treeorder = new int[nyuko + 1];
            var topol0 = new int[Math.Max(nyuko, 1)][];
            var topol1 = new int[Math.Max(nyuko, 1)][];
            if (nyuko == 1)
            {
                if (npick != 1)
                    throw new InvalidOperationException($"okashii, nyuko = 1, shikashi npick = {npick}");
                uniform = -1;
                for (int j = 0; j < nin; j++)
                    outs[0].Add(scores[j]);
            }
            else
            {
                Report("\n\n" + nyuko + "x" + nin + " distance matrix\n");
                for (int i = 1; i < nyuko; i++)
                {
                    Report(i + " / " + nyuko + " \r");
                    double[]? res = null;
                    short[]? t1 = null;
                    string? q1 = null;
                    if (doalign != 0)
                    {
                        if (doalign == 'f')
                        {
                            var lib = new List<string>(nin);
                            for (int a = 0; a < nin; a++)
                                lib.Add(_spSeq[scores[a].NumInSeq]);
                            res = MafftFasta34.Search(_spSeq[scores[yukos[i]].NumInSeq], lib, G.Dorp == 'd');
                        }
                        else
                        {
                            q1 = _spSeq[scores[yukos[i]].NumInSeq];
                        }
                    }
                    else
                    {
                        t1 = new short[_spTsize];
                        foreach (int p in scores[yukos[i]].Pointt)
                        {
                            if (p == -1)
                                break;
                            t1[p]++;
                        }
                    }
                    int sy = scores[yukos[i]].SelfScore;
                    for (int j = 0; j < nin; j++)
                    {
                        double lenfac = LenFac(scores[yukos[i]].OriLen, scores[j].OriLen);
                        if (dfromc[i][j] == -0.5)
                        {
                            if (doalign != 0)
                            {
                                if (doalign == 'f')
                                {
                                    dfromc[i][j] = (1.0 - res![j] / Math.Min(sy, scores[j].SelfScore)) * 1;
                                    if (dfromc[i][j] < 0.0)
                                        dfromc[i][j] = 0.0;
                                }
                                else if (G.FromAln)
                                {
                                    dfromc[i][j] = (1.0 - NaivePairScore11(_spOriAln![scores[j].NumInSeq], _spOriAln[scores[yukos[i]].NumInSeq], G.Penalty) / Math.Min(sy, scores[j].SelfScore)) * 1;
                                }
                                else
                                {
                                    dfromc[i][j] = (1.0 - G__align11NoalignLN(q1!, _spSeq[scores[j].NumInSeq]) / Math.Min(sy, scores[j].SelfScore)) * 1;
                                }
                            }
                            else
                            {
                                dfromc[i][j] = (1.0 - (double)LocalCommonSextet(t1!, scores[j].Pointt) / Math.Min(sy, scores[j].SelfScore)) * lenfac;
                                if (dfromc[i][j] > SplitMax6Dist)
                                    dfromc[i][j] = SplitMax6Dist;
                            }
                        }
                    }
                }
                Report("                \r");
                for (int j = 0; j < nin; j++)
                {
                    int bt = sYMap[j];
                    if (bt == -1)
                    {
                        bt = 0;
                        double minscore = dfromc[0][j];
                        for (int i = 0; i < nyuko; i++)
                        {
                            if (SameSequence(scores[j], scores[yukos[i]]))
                            {
                                bt = i;
                                break;
                            }
                            if (dfromc[i][j] < minscore)
                            {
                                minscore = dfromc[i][j];
                                bt = i;
                            }
                        }
                    }
                    outs[bt].Add(scores[j]);
                }
                var len = new double[nyuko][];
                for (int i = 0; i < nyuko; i++)
                    len[i] = new double[2];
                if (nyuko > 2)
                {
                    Report("upgma       ");
                    var dep = new Treedep[nyuko];
                    MuscleSupg(nyuko, yukomtx, topol0, topol1, len, dep, true);
                    Report("\r                      \r");
                }
                else
                {
                    topol0[0] = new[] { 0 };
                    topol1[0] = new[] { 1 };
                }
                int pos = 0;
                foreach (int v in topol0[nyuko - 2])
                    treeorder[pos++] = v;
                foreach (int v in topol1[nyuko - 2])
                    treeorder[pos++] = v;
            }
            int aligned = 1;
            for (int i = 0; i < nyuko; i++)
            {
                int ii = treeorder[i];
                var childScores = outs[ii].ToArray();
                string? ctree;
                aligned *= SplitNode(childScores, childScores.Length, uniform, out ctree, doalign, tokyoripara, scores[0].NumInSeq, tmpaminodis);
                if (children != null)
                    children[ii] = ctree ?? "";
            }
            for (int i = 0; i < nyuko; i++)
                if (outs[i].Count == 0)
                    throw new InvalidOperationException($"i={i}/{nyuko}, ERROR!");
            int val = (nin >= _spClassSize || aligned == 0) ? 0 : 1;
            if (nyuko > 1)
            {
                string[]? parttree = null;
                if (_spTreeOut)
                {
                    parttree = new string[nyuko];
                    for (int i = 0; i < nyuko; i++)
                        parttree[i] = children![i];
                }
                int v1 = 0;
                int nlim = nyuko - 1;
                int[] mem1 = Array.Empty<int>(), mem2 = Array.Empty<int>();
                for (int l = 0; l < nlim; l++)
                {
                    var m1 = new List<int>();
                    foreach (int idx in topol0[l])
                        foreach (var sc in outs[idx])
                            m1.Add(sc.NumInSeq);
                    var m2 = new List<int>();
                    foreach (int idx in topol1[l])
                        foreach (var sc in outs[idx])
                            m2.Add(sc.NumInSeq);
                    m1.Sort();
                    m2.Sort();
                    mem1 = m1.ToArray();
                    mem2 = m2.ToArray();
                    if (val != 0)
                    {
                        Report("\r  Alignment " + mem1.Length + "-" + mem2.Length + "                                 \r");
                        int[] firstGroup = mem1[0] < mem2[0] ? mem1 : mem2;
                        int[] secondGroup = mem1[0] < mem2[0] ? mem2 : mem1;
                        SplitPairAlign(firstGroup, secondGroup);
                    }
                    if (_spTreeOut)
                    {
                        v1 = topol0[l][0];
                        int v2 = topol1[l][0];
                        if (v1 > v2)
                            (v1, v2) = (v2, v1);
                        tree = "(" + parttree![v1] + "," + parttree[v2] + ")";
                        parttree[v1] = tree;
                        parttree[v2] = null!;
                    }
                }
                if (val != 0)
                {
                    _spGroupId++;
                    foreach (int x in mem1)
                        _spWhichGroup[x] = _spGroupId;
                    foreach (int x in mem2)
                        _spWhichGroup[x] = _spGroupId;
                }
                _spDepth++;
            }
            else if (_spTreeOut)
            {
                tree = children![0];
            }
            return val;
        }

        #endregion
    }
}
