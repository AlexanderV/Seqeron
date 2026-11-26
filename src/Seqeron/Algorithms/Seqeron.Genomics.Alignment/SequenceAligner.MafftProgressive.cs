using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: disttbfast (the progressive stage for three or more sequences)

    internal sealed partial class MafftRun
    {
        /// <summary>The result of one run of <c>disttbfast</c>: the aligned rows in the order of the input, the order in which the guide tree visits the sequences (the file <c>order</c> of <c>--reorder</c>) and the messages of the program (stderr).</summary>
        internal sealed record Result(string[] Rows, int[] Order)
        {
            /// <summary>--sparsepickup: the input index of every row of <see cref="Rows"/> (the other sequences are in <see cref="NotUsed"/>), otherwise null.</summary>
            internal int[]? Picked { get; init; }

            /// <summary>The file <c>notused</c> of --sparsepickup.</summary>
            internal string? NotUsed { get; init; }

            /// <summary>--keeplength: the residues of every added sequence that were deleted (the file <c>_deletelist</c>).</summary>
            internal List<(int Row, int Pos, int Len)>? Deleted { get; init; }

            /// <summary>--mapout: the file <c>_deletemap</c>.</summary>
            internal string? DeleteMap { get; init; }

            /// <summary>--compactmapout: the program names of the added sequences after "Nins:...|" was put in front (null entries: unchanged).</summary>
            internal string?[]? NewNames { get; init; }

            /// <summary>--merge: the file <c>GuideTree</c> that <c>system( "cp infile.tree GuideTree" )</c> left in the directory of the run when a subalignment is no group of the guide tree (the script copies it to the current directory).</summary>
            internal string? GuideTreeCopy { get; init; }
        }

        /// <summary>The stack size that <c>setstacksize</c> reports (the limit of the process: 8192 kb by default; the message is the only place where it shows).</summary>
        internal int StackSizeKb = 8192;

        /// <summary>
        /// <c>setstacksize( 200 * njob )</c> (io.c): the limit of the process is raised when it is smaller than 200 bytes per sequence, and the message then shows both sizes ("8192 kb->19531 kb", truncated to kb).
        /// </summary>
        private string StackSizeText(int njob) => StackSizeMessage(StackSizeKb, njob);

        /// <summary>The text of <c>setstacksize( 200 * njob )</c> for a process whose stack limit is <paramref name="stackKb"/> kb.</summary>
        internal static string StackSizeMessage(int stackKb, int njob)
        {
            long wanted = 200L * njob;
            string current = stackKb.ToString(CultureInfo.InvariantCulture);
            return (long)stackKb * 1024 < wanted
                ? "stacksize: " + current + " kb->" + (wanted / 1024).ToString(CultureInfo.InvariantCulture) + " kb\n"
                : "stacksize: " + current + " kb\n";
        }

        /// <summary>The file <c>hat2</c> is written again by a later program of the script (<c>dndpre</c> of the iterative methods): a guide tree without distances then still leaves the file that <c>cp</c> copies.</summary>
        internal bool DistFileRewrittenLater { get; init; }

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
        internal Result Disttbfast(IReadOnlyList<string> seq0, IReadOnlyList<string>? names = null)
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

            Report(StackSizeText(njob));
            MafftSubalignments? merge = G.Subalignment ? ReadSubalignmentsTable(njob) : null;
            if (G.UseExternalAnchors)
                PrepareExternalAnchors(njob, seq);
            G.Constants(seq);
            Report("Gap Penalty = " + MafftGlobals.F((double)G.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.Poffset / 1000, 5, 2, true) + "\n");
            if (CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            Report("\n");
            if (G.Dorp == 'p' && G.TupleSize != 6)
                throw new MafftProgramFailureException("tuplesize must be 6 for aa sequence\n");
            if (G.Dorp == 'd' && G.TupleSize != 6 && G.TupleSize != 10)
                throw new MafftProgramFailureException("tuplesize must be 6 or 10 for dna sequence\n");
            int nguidetree = G.NGuideTree;
            int randomseed = 0;
            char treein = '\0';
            if (G.TreeIn != 0)
            {
                var (kind, seed, npickx, limitram) = CheckGuideTreeFile();
                randomseed = seed;
                _maxDistMtxSize = limitram;
                treein = kind;
                if (treein == 't')
                {
                    VarPairScore(njob, npickx, seq, randomseed, out _);
                    throw new MafftProgramFailureException("", 1);
                }
                if (treein == 'c')
                {
                    _compactTree = 1;
                    treein = '\0';
                }
                else if (treein == 'Y')
                {
                    _compactTree = 4;
                    treein = '\0';
                }
                else if (treein is 'S' or 'C')
                {
                    _compactTree = 2;
                    treein = '\0';
                }
                else if (treein == 'a')
                {
                    if (njob > 200000)
                    {
                        Report("Chain?\n");
                        treein = 's';
                        nguidetree = 1;
                    }
                    else if (njob < 100 || VarPairScoreChoice(njob, npickx, seq, randomseed) == 't')
                    {
                        Report("Tree!\n");
                        treein = '\0';
                        nguidetree = 2;
                    }
                    else
                    {
                        Report("Chain!\n");
                        treein = 's';
                        nguidetree = 1;
                    }
                }
                else if (treein != '\0')
                {
                    nguidetree = 1;
                }
            }
            if (_compactTree == 1 && _maxDistMtxSize > (double)njob * (njob - 1) * sizeof(double) / 2)
            {
                Report("Use conventional tree.\n");
                _compactTree = 0;
            }
            if (G.DistOut != 0 && treein != '\0' && !DistFileRewrittenLater)
                throw new NotSupportedException("--distout with --pileup / --randomchain / --treein: the program writes no distance matrix and the message of the script's cp names a random temporary directory (cp: cannot stat '/tmp/mafft.XXXXXXXXXX/hat2').");
            if (G.DistOut != 0 && _compactTree != 0)
                throw new NotSupportedException("--distout with --large / --memsavetree / --stepadd / --youngestlinkage: the program writes the distance matrix that it does not hold and stops with a segmentation fault whose message names a process id.");
            double[][]? mtx = null;
            int[][]? pointt = null;
            int[]? mindistfrom = null;
            double[]? mindist = null;
            int[]? compactSelfscore = null;
            if (treein == '\0')
            {
                if (_compactTree != 0)
                {
                    pointt = PrepareKmerTables(seq);
                    Report("Compact tree, step 1\n");
                    mindistfrom = new int[njob];
                    mindist = new double[njob];
                    compactSelfscore = new int[njob];
                    PreparePartMtx(njob);
                    var memo = new int[_tsize];
                    for (int i = 0; i < njob; i++)
                    {
                        var table1 = new int[_tsize];
                        foreach (int point in pointt[i])
                            table1[point]++;
                        compactSelfscore[i] = CommonSextet(table1, pointt[i], memo);
                    }
                    for (int i = 0; i < njob; i++)
                    {
                        mindist[i] = 999.9;
                        mindistfrom[i] = -1;
                    }
                    CompactDistancesKmer(njob, pointt, compactSelfscore, mindist, mindistfrom);
                    for (int i = 0; i < njob; i++)
                        mindist[i] -= PreferenceVal(i, mindistfrom[i], njob);
                    Report("\ndone.\n");
                }
                else
                {
                    mtx = SixmerDistanceMatrix(seq);
                }
            }
            int[] order = new int[njob];
            string[]? bseq = null;
            string? originalGaps = null;
            string[]? addbk = null;
            if (G.Nadd > 0 && G.KeepLength != 0)
            {
                originalGaps = RecordOriginalGaps(njob - G.Nadd, seq);
                if (G.MapOut != 0)
                    addbk = seq.Skip(njob - G.Nadd).Select(x => x.Replace("-", "", StringComparison.Ordinal)).ToArray();
            }
            GuideTree? lastTree = null;
            double[]? lastEff = null;
            int lastAlloclen = 0;
            for (int iguidetree = 0; iguidetree < nguidetree; iguidetree++)
            {
                bool calcpairdists = iguidetree != nguidetree - 1;
                var topol0 = new int[njob][];
                var topol1 = new int[njob][];
                var len = new double[njob][];
                for (int i = 0; i < njob; i++)
                    len[i] = new double[2];
                var dep = new Treedep[njob];
                bool treeOut = G.TreeOut != 0;
                if (treein != '\0')
                {
                    calcpairdists = false;
                    if (treein == 'l')
                    {
                        LoadTree(njob, topol0, topol1, len, dep, treeOut, names ?? DefaultNames(njob));
                    }
                    else if (treein is 's' or 'p')
                    {
                        CreateChain(njob, topol0, topol1, len, dep, treeOut, treein == 's', randomseed, names ?? DefaultNames(njob));
                        G.NthreadPair = G.Nthread = 0;
                        G.Weight = 0;
                        G.Tbrweight = 0;
                        G.NthreadTb = 0;
                    }
                    else
                    {
                        throw new MafftProgramFailureException("Error. treein = " + (int)treein + " or " + treein + "\n");
                    }
                }
                else
                {
                    if (G.DistOut != 0)
                    {
                        DistFile = WriteHat2(mtx!, names ?? DefaultNames(njob));
                        if (G.TreeOut == 0 && G.NoAlign != 0)
                            return NoAlignResult(seq, njob);
                    }
                    if (merge != null)
                    {
                        Report("Constructing a UPGMA tree ... ");
                        ConstrainedSupg(njob, mtx!, topol0, topol1, len, dep, merge, names ?? DefaultNames(njob));
                    }
                    else if (_compactTree != 0)
                    {
                        Report("Constructing a tree ... ");
                        if (_compactTree == 4)
                            CompactTreeSelectable(njob, mindistfrom!, mindist!, iguidetree == 0 ? pointt : null, compactSelfscore!, iguidetree == 0 ? null : bseq, topol0, topol1, len, dep, treeOut, names ?? DefaultNames(njob), 0);
                        else
                            CompactTreeGivenDist(njob, mindist!, mindistfrom!, topol0, topol1, len, dep, treeOut, names ?? DefaultNames(njob));
                        mindistfrom = null;
                        mindist = null;
                        pointt = null;
                    }
                    else if (treeOut)
                    {
                        Report("Constructing a UPGMA tree (treeout, efffree=" + (calcpairdists ? 0 : 1) + ") ... ");
                        MuscleSupg(njob, mtx!, topol0, topol1, len, dep, false, names ?? DefaultNames(njob));
                    }
                    else
                    {
                        Report("Constructing a UPGMA tree (efffree=" + (calcpairdists ? 0 : 1) + ") ... ");
                        MuscleSupg(njob, mtx!, topol0, topol1, len, dep);
                    }
                }
                Report("\ndone.\n\n");
                if (G.SparsePickup != 0 && iguidetree == nguidetree - 1)
                {
                    Report("Sparsepickup! \n");
                    var picked = Pickup(njob, seq, topol0, topol1, names ?? DefaultNames(njob), treein != '\0');
                    Report("done. \n");
                    Report("disttbfast (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                        + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
                    return picked;
                }
                var tree = new GuideTree(topol0, topol1, dep, len);
                int[] lm;
                if (treein is 's' or 'p')
                {
                    lm = new int[njob];
                    lm[0] = tree.Topol0[0];
                    for (int i = 1; i < njob; i++)
                        lm[i] = tree.Topol1[i - 1];
                }
                else
                {
                    lm = tree.Order(njob - 2, 2);
                }
                Array.Copy(lm, order, njob);
                if ((G.TreeOut != 0 || G.DistOut != 0) && G.NoAlign != 0)
                    return NoAlignResult(seq, njob);
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
                MafftAddState? addState = null;
                if (G.Nadd > 0)
                {
                    addState = PrepareAdd(njob, G.Nadd, seq, bseq, tree, names ?? DefaultNames(njob));
                    if (addState == null)
                        return new Result(Array.Empty<string>(), Array.Empty<int>()) { Picked = Array.Empty<int>() };
                }
                else if (merge != null)
                {
                    addState = PrepareMerge(njob, seq, bseq, topol0, topol1, merge, names ?? DefaultNames(njob), tbfastMessages: false);
                    if (addState == null)
                        return new Result(Array.Empty<string>(), Array.Empty<int>()) { Picked = Array.Empty<int>(), GuideTreeCopy = GuideTreeCopy };
                }
                else
                {
                    for (int i = 0; i < njob; i++)
                        bseq[i] = seq[i].Replace("-", "");
                }
                int alloclen = _nlenmax * 2 + 1;
                int[] selfscore = Array.Empty<int>();
                if (calcpairdists)
                {
                    selfscore = new int[njob];
                    for (int i = 0; i < njob; i++)
                        selfscore[i] = (int)NaivePairScore11(seq[i], seq[i], G.PenaltyDist);
                }
                Report("Progressive alignment " + (iguidetree + 1) + "/" + nguidetree + "... \n");
                TreeBase(bseq, tree, eff, ref alloclen, add: addState);
                Report("\ndone.\n\n");
                lastTree = tree;
                lastEff = eff;
                lastAlloclen = alloclen;
                if (calcpairdists && _compactTree == 0)
                {
                    Report("Making a distance matrix from msa.. \n");
                    mtx = MsaDistanceMatrix(njob, selfscore, bseq);
                    Report("\ndone.\n\n");
                }
                else if (calcpairdists)
                {
                    Report("Making a compact tree from msa, step 1.. \n");
                    mindistfrom = new int[njob];
                    mindist = new double[njob];
                    PreparePartMtx(njob);
                    for (int i = 0; i < njob; i++)
                    {
                        mindist[i] = 999.9;
                        mindistfrom[i] = -1;
                    }
                    CompactDistancesMsa(njob, bseq, selfscore, mindist, mindistfrom);
                    for (int i = 0; i < njob; i++)
                        mindist[i] -= PreferenceVal(i, mindistfrom[i], njob);
                    Report("\n");
                    compactSelfscore = selfscore;
                }
            }
            if (G.OneIteration)
            {
                if (G.TextMode)
                    throw new NotSupportedException("--text with --oneiteration: dooneiteration of the program stops with a segmentation fault (a stack overflow: valgrind 'Invalid write of size 8', the message of the script names a process id) on every text input.");
                Report("Iterative refinement (one vs others)\n");
                DoOneIteration(bseq!, lastTree!, lastEff!, ref lastAlloclen);
            }
            var (deletedRuns, deleteMap, newNames, ndeleted) = KeepLength(njob, seq, bseq!, originalGaps, addbk, names);
            if (G.SpScoreOut != 0)
                ReportSumOfPairs(bseq!);
            Report("disttbfast (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            if (ndeleted > 0)
            {
                Report("\nTo keep the alignment length, " + ndeleted.ToString(CultureInfo.InvariantCulture) + " letters were DELETED.\n");
                Report(G.MapOut != 0 ? "The deleted letters are shown in the (filename).map file.\n" : "To know the positions of deleted letters, rerun the same command with the --mapout option.\n");
            }
            return new Result(bseq!, order) { Deleted = deletedRuns, DeleteMap = deleteMap, NewNames = newNames };
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
            var pointt = PrepareKmerTables(seq);
            int tsize = _tsize;
            var mtx = new double[njob][];
            for (int i = 0; i < njob; i++)
                mtx[i] = new double[njob - i];
            var memo = new int[tsize];
            for (int i = 0; i < njob; i++)
            {
                var table1 = new int[tsize];
                if (i % 100 == 0)
                    Report("\r" + MafftGlobals.D5(i + 1) + " / " + njob);
                foreach (int point in pointt[i])
                    table1[point]++;
                for (int j = i; j < njob; j++)
                    mtx[i][j - i] = CommonSextet(table1, pointt[j], memo);
            }
            Report("\ndone.\n\n");
            int ien = njob - 1;
            double lenfaca = _lenfaca, lenfacb = _lenfacb, lenfacc = _lenfacc, lenfacd = _lenfacd;
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

        private int _tsize;
        private double _lenfaca, _lenfacb, _lenfacc, _lenfacd;

        /// <summary>The first half of the distance stage of <c>disttbfast</c>: "Making a distance matrix ..", the k-mer code tables <c>pointt</c> of every sequence (the lengths without gaps go to <c>nogaplen</c>), the count of the ambiguous letters, and the length factors of the distance.</summary>
        private int[][] PrepareKmerTables(string[] seq)
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
            _tsize = tsize;
            _lenfaca = lenfaca; _lenfacb = lenfacb; _lenfacc = lenfacc; _lenfacd = lenfacd;
            return pointt;
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

        /// <summary><c>sumofpairsscore</c> (<c>--scoreout</c>, <c>-S</c>): the message "Unweighted sum-of-pairs score".</summary>
        private void ReportSumOfPairs(IReadOnlyList<string> seq)
        {
            double v = 0;
            for (int i = 1; i < seq.Count; i++)
                for (int j = 0; j < i; j++)
                    v += NaivePairScore11(seq[i], seq[j], G.Penalty) / 600;
            Report("Unweighted sum-of-pairs score = " + MafftGlobals.F(v, 10, 5) + "\n");
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
                    Report("\r" + MafftGlobals.D5(i) + " / " + njob);
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
        private void MuscleSupg(int nseq, double[][] eff, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool fullMembership = false, IReadOnlyList<string>? treeNames = null, bool progress = true)
        {
            string[]? tree = null;
            double[]? density = null;
            if (treeNames != null)
            {
                tree = new string[nseq];
                for (int i = 0; i < nseq; i++)
                    tree[i] = "\n" + (i + 1).ToString(CultureInfo.InvariantCulture) + "_" + TreeLabel(treeNames[i]) + "\n";
            }
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
            if (tree != null && G.TreeOut == 2)
            {
                density = new double[nseq];
                for (int i = 0; i < nseq; i++)
                    density[i] = SetDensity(eff, next, nseq, i);
            }
            for (int i = 0; i < nseq; i++)
            {
                hist[i] = -1;
                nmemar[i] = 1;
            }
            if (progress)
                Report("\n");
            int im = -1, jm = -1;
            for (int k = 0; k < nseq - 1; k++)
            {
                if (progress && k % 10 == 0)
                    Report("\r" + MafftGlobals.D5(k) + " / " + nseq);
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
                if (tree != null)
                {
                    tree[im] = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                    tree[jm] = null!;
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
            if (tree != null)
            {
                // a single sequence: im is -1 and the program prints the pointer in front of its array of trees, which is NULL ("(null)" of glibc)
                var file = new StringBuilder(im < 0 ? "(null)" : tree[im]);
                file.Append(";\n");
                if (G.TreeOut == 2)
                {
                    var gt = new GuideTree(topol0, topol1, dep, len);
                    file.Append("\nDensity:");
                    for (int k = 0; k < nseq; k++)
                        file.Append("\nSequence ").Append(k + 1).Append(", ").Append(MafftGlobals.F(density![k], 7, 4));
                    file.Append("\n\nNode info:");
                    for (int k = 0; k < nseq - 1; k++)
                    {
                        file.Append("\nNode ").Append(k + 1).Append(", Height=").Append(dep[k].DistFromTip.ToString("F6", CultureInfo.InvariantCulture)).Append('\n');
                        foreach (int child in new[] { 0, 1 })
                        {
                            int[] mem = gt.Order(k, child);
                            file.Append(Densest(mem, density!) + 1).Append(':');
                            foreach (int mm in mem)
                                file.Append(' ').Append(mm + 1);
                            file.Append('\n');
                        }
                    }
                }
                TreeFile = file.ToString();
            }
        }

        /// <summary>The text of the file <c>infile.tree</c> (<c>--treeout</c> / <c>--nodeout</c>) after the last guide tree, and the distance file <c>hat2</c> (<c>--distout</c>) written for the last guide tree.</summary>
        internal string? TreeFile, DistFile;

        private static IReadOnlyList<string> DefaultNames(int n) => Enumerable.Range(0, n).Select(i => "=seq" + (i + 1).ToString(CultureInfo.InvariantCulture)).ToArray();

        /// <summary>The name of a sequence in the tree file: letters, digits, '/', '=', '-', '{' and '}' stay, everything else is replaced by '_'; the first character (the '=' of the program's name) and a prefix up to <c>_oe_</c> are cut.</summary>
        private string TreeLabel(string name)
        {
            var sb = new StringBuilder();
            foreach (char c in name.Length > 999 ? name[..999] : name)
                sb.Append(char.IsAsciiLetterOrDigit(c) || c is '/' or '=' or '-' or '{' or '}' ? c : '_');
            string all = sb.ToString();
            string nameptr = sb.Length > 0 ? all[1..] : "";
            if (G.OutNumber != 0)
            {
                // outnumber: nameptr = strstr( nametmp, "_numo_e" ) + 8
                int numo = all.IndexOf("_numo_e", StringComparison.Ordinal);
                nameptr = numo >= 0 && numo + 8 <= all.Length ? all[(numo + 8)..] : "";
            }
            int oe = nameptr.IndexOf("_oe_", StringComparison.Ordinal);
            if (oe >= 0)
                nameptr = nameptr[(oe + 4)..];
            return nameptr.Length > 900 ? nameptr[..900] : nameptr;
        }

        /// <summary><c>setdensity</c>: the sum of 2 - distance over the sequences that are closer than 1.</summary>
        private static double SetDensity(double[][] eff, int[] next, int nseq, int pos)
        {
            double d = 0.0;
            for (int j = next[pos]; j != -1 && j < nseq; j = next[j])
                if (eff[pos][j - pos] < 1.0)
                    d += 2.0 - eff[pos][j - pos];
            for (int j = 0; j != -1 && j != pos; j = next[j])
                if (eff[j][pos - j] < 1.0)
                    d += 2.0 - eff[j][pos - j];
            return d;
        }

        private static int Densest(int[] members, double[] density)
        {
            double dmax = -100.0;
            int pmax = -1;
            foreach (int m in members)
            {
                if (density[m] > dmax)
                {
                    dmax = density[m];
                    pmax = m;
                }
            }
            return pmax;
        }

        /// <summary><c>WriteFloatHat2_pointer_halfmtx</c>: the distance file of <c>--distout</c>.</summary>
        private static string WriteHat2(double[][] mtx, IReadOnlyList<string> names)
        {
            int n = names.Count;
            double max = 0.0;
            for (int i = 0; i < n - 1; i++)
                for (int j = 1; j < n - i; j++)
                    if (mtx[i][j] > max)
                        max = mtx[i][j];
            var sb = new StringBuilder();
            sb.Append(MafftGlobals.F(1, 5, 0)).Append('\n').Append(n.ToString(CultureInfo.InvariantCulture).PadLeft(5)).Append('\n');
            sb.Append(' ').Append(Hat2Number(max * 2.5)).Append('\n');
            for (int i = 0; i < n; i++)
                sb.Append((i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4)).Append(". ").Append(names[i]).Append('\n');
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    sb.Append(Hat2Number(mtx[i][j - i]));
                    int ijsa = j - i;
                    if (ijsa % 12 == 0 || ijsa == n - i - 1)
                        sb.Append('\n');
                }
            }
            return sb.ToString();
        }

        /// <summary><c>%#6.3f</c>.</summary>
        private static string Hat2Number(double v) => v.ToString("F3", CultureInfo.InvariantCulture).PadLeft(6);

        /// <summary>The output of <c>disttbfast -T</c> (<c>--retree 0</c>, or the guide tree / distances only): the sequences as they were read, not aligned.</summary>
        private Result NoAlignResult(string[] seq, int njob, string program = "disttbfast")
        {
            Report("\n");
            Report(program + " (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new Result((string[])seq.Clone(), Enumerable.Range(0, njob).ToArray());
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
        private void CountEffSimple(int nseq, GuideTree tree, double[] node)
        {
            var memhist = new int[nseq - 1][];
            CountEffCore(nseq, tree.Len, node, i =>
            {
                int[] m0 = tree.Dep[i].Child0 == -1 ? new[] { tree.Topol0[i] } : memhist[tree.Dep[i].Child0];
                int[] m1 = tree.Dep[i].Child1 == -1 ? new[] { tree.Topol1[i] } : memhist[tree.Dep[i].Child1];
                var all = new int[m0.Length + m1.Length];
                m0.CopyTo(all, 0);
                m1.CopyTo(all, m0.Length);
                memhist[i] = all;
                return (m0, m1);
            });
        }

        /// <summary><c>counteff_simple_double_nostatic</c> (the groups of every step are given in full, <c>addsingle</c>): the same weights as <see cref="CountEffSimple(int, GuideTree, double[])"/>.</summary>
        private void CountEffSimple(int nseq, int[][] groups0, int[][] groups1, double[][] len, double[] node) =>
            CountEffCore(nseq, len, node, i => (groups0[i], groups1[i]));

        private void CountEffCore(int nseq, double[][] len, double[] node, Func<int, (int[] Group0, int[] Group1)> groups)
        {
            for (int i = 0; i < nseq; i++)
            {
                if (len[i][0] < 0.0)
                {
                    Report("WARNING: negative branch length " + len[i][0].ToString("F6", CultureInfo.InvariantCulture) + ", step " + i + "-0\n");
                    len[i][0] = 0.0;
                }
                if (len[i][1] < 0.0)
                {
                    Report("WARNING: negative branch length " + len[i][1].ToString("F6", CultureInfo.InvariantCulture) + ", step " + i + "-1\n");
                    len[i][1] = 0.0;
                }
            }
            var rootnode = new double[nseq];
            var eff = new double[nseq];
            for (int i = 0; i < nseq; i++)
            {
                rootnode[i] = 0.0;
                eff[i] = 1.0;
            }
            for (int i = 0; i < nseq - 1; i++)
            {
                var (m0, m1) = groups(i);
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

        /// <summary><c>treebase</c> (<c>disttbfast.c</c>, all groups aligned: <c>mergeoralign[l] == 'a'</c>): joins the groups in the order of the guide tree. <paramref name="threadsTb"/> &gt; 0: the program runs <c>treebasethread</c> of <c>tbfast.c</c>
        /// (<c>nthreadtb &gt; 0 &amp;&amp; nadd == 0</c>), whose standard error differs from <c>treebase</c>: "STEP k /m (thread t) " and "Reallocating (by thread t) ..", the buffer grows when <c>alloclen &lt;= len1 + len2</c> (<c>treebase</c>: <c>&lt;</c>). The threads take the steps from a queue in the order of
        /// the steps; which thread takes which step is up to the scheduler (the port names thread <c>l % threadsTb</c>, one of the possible schedules; with one thread the program prints thread 0 every time); the alignment is the same.</summary>
        private void TreeBase(string[] aseq, GuideTree tree, double[] effarr, ref int alloclen, MafftLocalHom[][]? localhom = null, bool tbfast = false, MafftAddState? add = null, NodePairData? nodePair = null, double[]? effarrKozo = null, MafftLocalHom[][]? seedHom = null, int nseed = 0, int threadsTb = 0)
        {
            bool constraint = localhom != null || nodePair != null && G.Constraint == 2;
            double dumdb = 0.0;
            int njob = _njob;
            // tbfast without constraints (--weighti 0): the FFT is tried again for a group only when the last FFT alignment of its first member had segments (fftlog)
            int[]? fftlog = tbfast ? Enumerable.Repeat(1, njob).ToArray() : null;
            var memhist = new int[njob - 1][];
            var cpmxhist = new double[njob - 1][][];
            if (localhom != null)
                CalcImportanceHalf(njob, effarr, aseq, localhom, alloclen);
            int[]? seedIn1 = null, seedIn2 = null;
            if (constraint && nodePair != null && seedHom != null && nseed > 0)
            {
                // compacttree == 3 with seeds: dontcalcimportance_half( nseed, ... ) (the importance of a block is its score), the positions of the seeds of both groups in the table of the seeds
                for (int i = 0; i < nseed; i++)
                    for (int j = i + 1; j < nseed; j++)
                        for (var t = seedHom[i][j - i]; t != null; t = t.Next)
                            t.Importance = t.Opt;
                seedIn1 = new int[njob];
                seedIn2 = new int[njob];
            }
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
                if (add != null && add.MergeOrAlign[l] == 'n')
                    continue; // both groups are old sequences that are aligned already
                // tbfast makes the matrix of every step with makedynamicmtx (disttbfast only when --amax > 0): the matrix is shifted by the (negative) offset of the distance of the step from the tips, which is
                // zero for ordinary distances and not for the negative distances that a pair library can give
                double[][] dynamicmtx = constraint || localhom != null || tbfast ? MakeDynamicMtx(G.NDisConsweightMulti, tree.Dep[l].DistFromTip) : G.NDisConsweightMulti;
                int len1 = aseq[m1].Length;
                int len2 = aseq[m2].Length;
                if (threadsTb > 0 ? alloclen <= len1 + len2 : alloclen < len1 + len2)
                {
                    Report(threadsTb > 0 ? "\nReallocating (by thread " + (l % threadsTb).ToString(CultureInfo.InvariantCulture) + ") .." : "\nReallocating..");
                    alloclen = len1 + len2 + 1000;
                    Report("done. *alloclen = " + alloclen + "\n");
                }
                var mseq1 = new string[clus1];
                var mseq2 = new string[clus2];
                var effarr1 = new double[clus1];
                var effarr2 = new double[clus2];
                double minimumWeight = constraint || tbfast ? G.MinimumWeight : 0.0;
                double[]? effarr1Kozo = null, effarr2Kozo = null;
                double orieff1, orieff2;
                if (effarrKozo != null)
                {
                    // structures (blocks of kind 'k'): fastconjuction_noname_kozo, no minimum weight, no reuse of the profiles of the children (cpmxchild0 = cpmxchild1 = NULL; orieff1 / orieff2 are not set and not used)
                    effarr1Kozo = new double[clus1];
                    effarr2Kozo = new double[clus2];
                    FastConjunctionKozo(localmem0, aseq, mseq1, effarr1, effarr, effarr1Kozo, effarrKozo);
                    FastConjunctionKozo(localmem1, aseq, mseq2, effarr2, effarr, effarr2Kozo, effarrKozo);
                    orieff1 = orieff2 = 1.0;
                    child0 = child1 = -1;
                }
                else
                {
                    orieff1 = FastConjunction(localmem0, aseq, mseq1, effarr1, effarr, minimumWeight);
                    orieff2 = FastConjunction(localmem1, aseq, mseq2, effarr2, effarr, minimumWeight);
                }
                int len1nocommongap = len1;
                int[]? gapmap = null;
                char mergeOrAlign = add?.MergeOrAlign[l] ?? 'a';
                G.NewGapChar = mergeOrAlign is '1' or '2' ? '=' : '-';
                if (mergeOrAlign == '1')
                {
                    gapmap = FindCommonGaps(clus2, mseq2);
                    CommonGapPickN(clus2, mseq2);
                }
                else if (mergeOrAlign == '2')
                {
                    gapmap = FindCommonGaps(clus1, mseq1);
                    CommonGapPickN(clus1, mseq1);
                    len1nocommongap = mseq1[0].Length;
                }
                if (constraint)
                {
                    // tbfast, constraint == 2: the importance of the blocks of the pairwise alignments enters the match scores of A__align
                    if (l < 1000 || l % 100 == 0)
                        Report("\rSTEP " + MafftGlobals.D5(l + 1) + " /" + (njob - 1) + (threadsTb > 0 ? " (thread " + D4(l % threadsTb) + ") " : " "));
                    MafftLocalHom?[][]? seedShrink = null;
                    if (seedIn1 != null)
                    {
                        seedShrink = ShrinkSeedHom(localmem0, localmem1, nseed, seedIn1, seedIn2!, seedHom!);
                        for (int i = 0; i < njob; i++)
                            Report("seedinlh1[" + i.ToString(CultureInfo.InvariantCulture) + "]=" + seedIn1[i].ToString(CultureInfo.InvariantCulture) + "\n");
                        for (int i = 0; i < njob; i++)
                            Report("seedinlh2[" + i.ToString(CultureInfo.InvariantCulture) + "]=" + seedIn2![i].ToString(CultureInfo.InvariantCulture) + "\n");
                    }
                    if (G.Alg == 'M')
                        throw new MafftProgramFailureException("\n\nMemory saving mode is not supported.\n\n");
                    var shrink = nodePair != null ? NodeShrink(nodePair, l, clus1, clus2) : ShrinkLocalHom(localmem0, localmem1, localhom!);
                    var imp = FillImp(clus1, clus2, mseq1[0].Length, mseq2[0].Length, mseq1, mseq2, effarr1, effarr2, shrink, localmem0, localmem1, false, effarr1Kozo, effarr2Kozo, seedShrink, seedIn1, seedIn2, G.Alg == 'd');
                    if (G.Alg == 'd')
                        DAlign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, imp, ref dumdb, G.Outgap, G.Outgap);
                    else
                        AAlign(dynamicmtx, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap,
                            cpmxhist, child0, child1, l, orieff1, orieff2, imp);
                    _nlen[m1] = (int)(0.5 * (_nlen[m1] + _nlen[m2]));
                    MergeAfterAlignment(aseq, njob, add, mergeOrAlign, gapmap, len1, len1nocommongap, clus1, clus2, mseq1, mseq2, localmem0, localmem1, alloclen);
                    continue;
                }
                if (tbfast)
                {
                    if (l < 1000 || l % 100 == 0)
                        Report("\rSTEP " + MafftGlobals.D5(l + 1) + " /" + (njob - 1) + (threadsTb > 0 ? " (thread " + D4(l % threadsTb) + ") " : " "));
                    bool tryFft = fftlog![m1] != 0 && fftlog[m2] != 0 && _nlen[m1] > clus1 && _nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000;
                    if (G.ForceFft || (G.UseFft && tryFft))
                    {
                        Report(" f\b\b");
                        Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                        fftlog[m1] = FftlogResult;
                    }
                    else
                    {
                        Report(" d\b\b");
                        fftlog[m1] = 0;
                        if (G.Alg == 'd')
                        {
                            double noImp = 0.0;
                            DAlign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, ref noImp, G.Outgap, G.Outgap);
                        }
                        else if (G.Alg != 'A')
                            throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                        else
                            AAlign(dynamicmtx, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap,
                                cpmxhist, child0, child1, l, orieff1, orieff2);
                    }
                    _nlen[m1] = (int)(0.5 * (_nlen[m1] + _nlen[m2]));
                    MergeAfterAlignment(aseq, njob, add, mergeOrAlign, gapmap, len1, len1nocommongap, clus1, clus2, mseq1, mseq2, localmem0, localmem1, alloclen);
                    continue;
                }
                if (l < 500 || l % 100 == 0)
                    Report("\rSTEP " + MafftGlobals.D5(l + 1) + " / " + (njob - 1) + " ");
                if (G.NeverMemSave == 0 && G.Alg != 'M' && (len1 > 30000 || len2 > 30000))
                {
                    Report("\nlen1=" + len1.ToString(CultureInfo.InvariantCulture) + ", len2=" + len2.ToString(CultureInfo.InvariantCulture) + ", Switching to the memsave mode\n");
                    G.Alg = 'M';
                    G.Common.FreeIp();
                }
                bool ffttry = _nlen[m1] > clus1 && _nlen[m2] > clus2;
                if (G.UseExternalAnchors)
                {
                    var pairanch = PickPairAnch(clus1, clus2, localmem0, localmem1, mseq1, mseq2);
                    FalignGivenAnchors(pairanch, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                }
                else if (G.ForceFft || (G.UseFft && ffttry))
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" f\b\b");
                    if (G.Alg == 'M')
                    {
                        if (l < 500 || l % 100 == 0)
                            Report("m");
                        Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, memsave: true);
                    }
                    else
                    {
                        Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                    }
                }
                else
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" d\b\b");
                    switch (G.Alg)
                    {
                        case 'M':
                            if (l < 500 || l % 100 == 0)
                                Report("m");
                            if ((l < 500 || l % 100 == 0) && ((child1 != -1 && cpmxhist[child1] != null) || (child0 != -1 && cpmxhist[child0] != null)))
                                Report(" h");
                            MsAlign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, cpmxhist, child0, child1, l, orieff1, orieff2);
                            break;
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
                        case 'd':
                            if (clus1 == 1 && clus2 == 1)
                            {
                                G__align11(dynamicmtx, mseq1, mseq2, alloclen, G.Outgap, G.Outgap);
                            }
                            else
                            {
                                double noImp = 0.0;
                                DAlignLs(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, ref noImp, G.Outgap, G.Outgap);
                            }
                            break;
                        default:
                            throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                    }
                }
                _nlen[m1] = (int)(0.5 * (_nlen[m1] + _nlen[m2]));
                MergeAfterAlignment(aseq, njob, add, mergeOrAlign, gapmap, len1, len1nocommongap, clus1, clus2, mseq1, mseq2, localmem0, localmem1, alloclen);
            }
            G.NewGapChar = '-';
            FreeStatics();
        }

        /// <summary><c>--keeplength</c> / <c>--mapout</c> of <c>disttbfast</c> and <c>tbfast</c> (<c>deletenewinsertions_whole</c>, <c>restoreoriginalgaps</c>, <c>reconstructdeletemap</c>): the letters of the added sequences that would insert columns into the
        /// alignment are deleted, the columns that the old alignment had only gaps in are put back, and the map of the deleted letters is written.</summary>
        private (List<(int Row, int Pos, int Len)>? DeletedRuns, string? DeleteMap, string?[]? NewNames, int NDeleted) KeepLength(int njob, string[] seq, string[] bseq, string? originalGaps, string[]? addbk, IReadOnlyList<string>? names)
        {
            List<(int Row, int Pos, int Len)>? deletedRuns = null;
            string? deleteMap = null;
            string?[]? newNames = null;
            int ndeleted = 0;
            if (G.Nadd > 0 && G.KeepLength != 0)
            {
                int nadd = G.Nadd;
                var oldRows = bseq.Take(njob - nadd).ToArray();
                var addRows = bseq.Skip(njob - nadd).ToArray();
                var deletelist = DeleteNewInsertions(njob - nadd, nadd, oldRows, addRows);
                deletedRuns = new List<(int Row, int Pos, int Len)>();
                for (int i = 0; i < nadd; i++)
                    foreach (var d in deletelist[i])
                        deletedRuns.Add((njob - nadd + i, d.Pos, d.Len));
                ndeleted = bseq[0].Length - oldRows[0].Length;
                for (int i = 0; i < njob - nadd; i++)
                    bseq[i] = oldRows[i];
                for (int i = 0; i < nadd; i++)
                    bseq[njob - nadd + i] = addRows[i];
                RestoreOriginalGaps(njob, bseq, originalGaps!);
                if (G.MapOut != 0)
                {
                    var allNames = names ?? DefaultNames(njob);
                    var addNames = allNames.Skip(njob - nadd).ToArray();
                    var realn = G.MapOut == 1 ? bseq.Skip(njob - nadd).ToArray() : seq.Skip(njob - nadd).ToArray();
                    if (G.MapOut == 1)
                    {
                        deleteMap = ReconstructDeleteMap(nadd, addbk!, deletelist, realn, addNames);
                    }
                    else
                    {
                        var changed = new string?[nadd];
                        deleteMap = ReconstructDeleteMapCompact(nadd, addbk!, deletelist, realn, addNames, changed!);
                        newNames = changed;
                    }
                }
            }
            return (deletedRuns, deleteMap, newNames, ndeleted);
        }

        /// <summary>The end of a step of <c>treebase</c> (<c>disttbfast.c</c>, <c>tbfast.c</c>): with <c>--add</c> the columns that only gaps of the old group had are restored and the new columns are put into the other old sequences
        /// (<c>mergeoralign</c> '2'); then the rows of both groups go back into the alignment.</summary>
        private void MergeAfterAlignment(string[] aseq, int njob, MafftAddState? add, char mergeOrAlign, int[]? gapmap, int len1, int len1nocommongap, int clus1, int clus2, string[] mseq1, string[] mseq2, int[] localmem0, int[] localmem1, int alloclen,
            bool addSingle = false, int len2 = 0, int len2nocommongap = 0)
        {
            if (mergeOrAlign == '1' && !addSingle)
                throw new MafftProgramFailureException(Log + "Check source!!!\n");
            if (mergeOrAlign == '1')
            {
                // addsingle.c: the mirror image of the branch '2' (the first group is the one that was added, the second one the old alignment whose common gaps were taken out)
                (mseq1, mseq2) = (mseq2, mseq1);
                (localmem0, localmem1) = (localmem1, localmem0);
                (clus1, clus2) = (clus2, clus1);
                (len1, len1nocommongap) = (len2, len2nocommongap);
                mergeOrAlign = '2';
            }
            if (mergeOrAlign == '2')
            {
                // restore the columns that only gaps of the old group had, put the new columns into the other old sequences
                int gapmaplen = mseq1[0].Length - len1nocommongap + len1;
                gapmap = AdjustGapMap(gapmaplen, gapmap!, mseq1[0]);
                for (int k = 0; k < clus1; k++)
                    aseq[localmem0[k]] = mseq1[k];
                for (int k = 0; k < clus2; k++)
                    aseq[localmem1[k]] = mseq2[k];
                bool smoothing = !addSingle && G.Smoothing != 0;
                gapmap = smoothing ? RestoreCommonGapsSmoothly(njob, njob - (clus1 + clus2), aseq, localmem0, localmem1, gapmap, '-') : RestoreCommonGaps(njob, njob - (clus1 + clus2), aseq, localmem0, localmem1, gapmap, '-');
                for (int k = 0; k < clus1; k++)
                    mseq1[k] = aseq[localmem0[k]];
                for (int k = 0; k < clus2; k++)
                    mseq2[k] = aseq[localmem1[k]];
                var gaplen = FindNewGaps(mseq1, 0);
                if (smoothing)
                    InsertNewGapsBothOrders(njob, add!.AlreadyAligned, aseq, localmem0, localmem1, gaplen, gapmap, gapmaplen, alloclen, '-');
                else
                    InsertNewGaps(njob, add!.AlreadyAligned, aseq, localmem0, localmem1, gaplen, gapmap, alloclen, '-');
                for (int k = 0; k < clus1; k++)
                    mseq1[k] = aseq[localmem0[k]];
                for (int k = 0; k < clus2; k++)
                    mseq2[k] = aseq[localmem1[k]];
                if (!addSingle)
                {
                    // disttbfast / tbfast resolve the new columns ('=') at the end of the step; addsingle keeps them: the old sequences show where the added sequence made new columns
                    Eq2DashMatome(mseq1, clus1);
                    Eq2DashMatome(mseq2, clus2);
                }
                foreach (int m in localmem1)
                    add.AlreadyAligned[m] = true;
            }
            for (int k = 0; k < clus1; k++)
                aseq[localmem0[k]] = mseq1[k];
            for (int k = 0; k < clus2; k++)
                aseq[localmem1[k]] = mseq2[k];
        }

        /// <summary>
        /// The part of <c>main</c> of <c>disttbfast.c</c> for <c>-I n</c> (<c>--add</c> with <c>-K</c>, <c>--addprofile</c> without): the checks, which step of the tree does what (<c>mergeoralign</c>) and the rows of <paramref name="bseq"/>.
        /// Returns null when the program stops after its message (the original sequences are not aligned).
        /// </summary>
        private MafftAddState? PrepareAdd(int njob, int nadd, string[] seq, string[] bseq, GuideTree tree, IReadOnlyList<string> names, bool tbfastMessages = false)
        {
            _ = names;
            int alignmentlength = seq[0].Length;
            for (int i = 0; i < njob - nadd; i++)
            {
                if (seq[i].Length != alignmentlength)
                {
                    if (tbfastMessages)
                        throw new MafftProgramFailureException("#################################################################################\n# ERROR!                                                                        #\n# The original" + (njob - nadd).ToString(CultureInfo.InvariantCulture).PadLeft(4)
                            + " sequences must be aligned                                    #\n#################################################################################\n");
                    Report("#################################################################################\n# ERROR!\n# The original " + (njob - nadd).ToString(CultureInfo.InvariantCulture) + " sequences must be aligned\n# alignmentlength = "
                        + alignmentlength.ToString(CultureInfo.InvariantCulture) + ", but strlen(seq[" + i.ToString(CultureInfo.InvariantCulture) + "])=" + seq[i].Length.ToString(CultureInfo.InvariantCulture) + "\n#################################################################################\n");
                    return null;
                }
            }
            var mergeoralign = new char[njob - 1];
            var addSet = new HashSet<int>(Enumerable.Range(njob - nadd, nadd));
            if (G.AddProfile)
            {
                alignmentlength = seq[njob - nadd].Length;
                for (int i = njob - nadd; i < njob; i++)
                {
                    if (seq[i].Length != alignmentlength)
                    {
                        if (tbfastMessages)
                            throw new MafftProgramFailureException("###############################################################################\n# ERROR!                                                                      #\n# The" + nadd.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                + " additional sequences must be aligned                                #\n# Otherwise, try the '--add' option, instead of '--addprofile' option.        #\n###############################################################################\n");
                        throw new MafftProgramFailureException("###############################################################################\n# ERROR!\n# The additional sequences must be aligned\n".Replace("additional sequences", nadd.ToString(CultureInfo.InvariantCulture) + " additional sequences", StringComparison.Ordinal)
                            + "# Otherwise, try the '--add' option, instead of '--addprofile' option.\n###############################################################################\n");
                    }
                }
                bool found = false;
                for (int i = 0; i < njob - 1; i++)
                {
                    var m0 = tree.Order(i, 0);
                    var m1 = tree.Order(i, 1);
                    if (m0.Length == nadd && m0.All(addSet.Contains))
                    {
                        mergeoralign[i] = '1';
                        found = true;
                    }
                    else if (m1.Length == nadd && m1.All(addSet.Contains))
                    {
                        mergeoralign[i] = '2';
                        found = true;
                    }
                    else
                    {
                        mergeoralign[i] = 'n';
                    }
                }
                if (!found && tbfastMessages)
                    throw new MafftProgramFailureException("###############################################################################\n# ERROR!                                                                      #\n# There is no appropriate position to add the" + nadd.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                        + " sequences in the guide tree.#\n# Check whether the" + nadd.ToString(CultureInfo.InvariantCulture).PadLeft(4) + " sequences form a monophyletic cluster.                #\n# If not, try the '--add' option, instead of the '--addprofile' option.       #\n############################################################################### \n");
                if (!found)
                    throw new MafftProgramFailureException("###############################################################################\n# ERROR!\n# There is no appropriate position to add the " + nadd.ToString(CultureInfo.InvariantCulture)
                        + " sequences in the guide tree.\n# Check whether the " + nadd.ToString(CultureInfo.InvariantCulture) + " sequences form a monophyletic cluster.\n# If not, try the '--add' option, instead of the '--addprofile' option.\n############################################################################### \n");
                var added = new string[nadd];
                for (int i = 0; i < nadd; i++)
                    added[i] = seq[njob - nadd + i];
                CommonGapPickN(nadd, added);
                for (int i = 0; i < nadd; i++)
                    bseq[njob - nadd + i] = added[i];
            }
            else
            {
                for (int i = 0; i < njob - 1; i++)
                {
                    mergeoralign[i] = 'n';
                    bool in0 = tree.Order(i, 0).All(addSet.Contains);
                    bool in1 = tree.Order(i, 1).All(addSet.Contains);
                    if (in0 && in1)
                        mergeoralign[i] = 'w';
                    else if (in0)
                        mergeoralign[i] = '1';
                    else if (in1)
                        mergeoralign[i] = '2';
                }
                for (int i = njob - nadd; i < njob; i++)
                    bseq[i] = seq[i].Replace("-", "", StringComparison.Ordinal);
            }
            var old = new string[njob - nadd];
            for (int i = 0; i < old.Length; i++)
                old[i] = seq[i];
            CommonGapPickN(njob - nadd, old);
            for (int i = 0; i < old.Length; i++)
            {
                seq[i] = old[i];
                bseq[i] = old[i];
            }
            var alreadyAligned = new bool[njob];
            for (int i = 0; i < njob - nadd; i++)
                alreadyAligned[i] = true;
            return new MafftAddState(mergeoralign, alreadyAligned);
        }

        /// <summary>
        /// <c>pickup</c> (<c>disttbfast -q n</c>, <c>--sparsepickup n</c>): the longest sequence of each side of the last steps of the guide tree (the root and its descendants until <paramref name="seq"/> n steps are used) are written, the others go to the file
        /// <c>notused</c>; nothing is aligned.
        /// </summary>
        private Result Pickup(int n, string[] seq, int[][] topol0, int[][] topol1, IReadOnlyList<string> names, bool zeroLengths = false)
        {
            // with a guide tree that is not computed from the distances the program has not measured the sequences: nogaplen is all zero
            var nogaplen = seq.Select(x => zeroLengths ? 0 : x.Replace("-", "", StringComparison.Ordinal).Length).ToArray();
            var longest = new int[n - 1, 2];
            for (int i = 0; i < n - 1; i++)
            {
                int best0 = -1, bestLen0 = -1, best1 = -1, bestLen1 = -1;
                foreach (int m in topol0[i])
                {
                    if (nogaplen[m] > bestLen0)
                    {
                        bestLen0 = nogaplen[m];
                        best0 = m;
                    }
                }
                foreach (int m in topol1[i])
                {
                    if (nogaplen[m] > bestLen1)
                    {
                        bestLen1 = nogaplen[m];
                        best1 = m;
                    }
                }
                longest[i, 0] = best0;
                longest[i, 1] = best1;
            }
            var select = new bool[n];
            int count = 1;
            for (int i = n - 2; i > -1; i--)
            {
                select[longest[i, 0]] = true;
                select[longest[i, 1]] = true;
                count++;
                if (count >= G.SparsePickup)
                    break;
            }
            var outIdx = new List<int>();
            var notIdx = new List<int>();
            for (int i = 0; i < n; i++)
                { if (select[i]) outIdx.Add(i); else notIdx.Add(i); }
            if (outIdx.Count != count)
                throw new NotSupportedException("The program writes " + count + " sequences but picked " + outIdx.Count + ": it reads names that do not exist and stops with a segmentation fault (the message names a process id).");
            string Strip(int i) => names[i][1..];
            string notUsed = MafftProgramIo.Fasta(notIdx.Select(Strip).ToList(), notIdx.Select(i => seq[i]).ToList());
            return new Result(outIdx.Select(i => seq[i]).ToArray(), Enumerable.Range(0, outIdx.Count).ToArray()) { Picked = outIdx.ToArray(), NotUsed = notUsed };
        }

        /// <summary>
        /// <c>dooneiteration</c> (<c>disttbfast -r</c>, <c>--oneiteration</c>): twice around the sequences, every sequence in turn is cut out of the alignment and aligned to the others again (<c>A__align</c> / <c>Falign</c>); the new alignment is kept
        /// when its sum-of-pairs score is not lower.
        /// </summary>
        private void DoOneIteration(string[] aseq, GuideTree tree, double[] effarr, ref int alloclen)
        {
            int njob = _njob;
            const int IterativeCycle = 2;
            for (int ll = 0; ll < njob * IterativeCycle; ll++)
            {
                AlignmentCancellation.ThrowIfRequested();
                int l = ll % njob;
                int[] localmem0 = { l };
                var localmem1 = new int[njob - 1];
                for (int i = 0, j = 0; i < njob; i++)
                    if (i != l)
                        localmem1[j++] = i;
                int clus1 = 1, clus2 = njob - 1;
                int m1 = l, m2 = localmem1[0];
                double distFromTip = l < tree.Dep.Length ? tree.Dep[l].DistFromTip : 0.0;
                double[][] dynamicmtx = G.SpecificityConsideration != 0.0 ? MakeDynamicMtx(G.NDisConsweightMulti, distFromTip) : G.NDisConsweightMulti;
                int len1 = aseq[m1].Length, len2 = aseq[m2].Length;
                if (alloclen < len1 + len2)
                {
                    Report("\nReallocating..");
                    alloclen = len1 + len2 + 1000;
                }
                var mseq1 = new string[clus1];
                var mseq2 = new string[clus2];
                var effarr1 = new double[clus1];
                var effarr2 = new double[clus2];
                double orieff1 = FastConjunction(localmem0, aseq, mseq1, effarr1, effarr, 0.0);
                double orieff2 = FastConjunction(localmem1, aseq, mseq2, effarr2, effarr, 0.0);
                double oscore = IntergroupScore(mseq1, mseq2, effarr1, effarr2, clus1, clus2, mseq1[0].Length);
                var oseq1 = (string[])mseq1.Clone();
                var oseq2 = (string[])mseq2.Clone();
                CommonGapPick(mseq2);
                CommonGapPick(mseq1);
                if (l < 500 || l % 100 == 0)
                    Report("\rIteration " + MafftGlobals.D5(ll + 1) + " / " + (njob * IterativeCycle).ToString(CultureInfo.InvariantCulture) + " ");
                if (G.NeverMemSave == 0 && G.Alg != 'M' && (len1 > 30000 || len2 > 30000))
                {
                    Report("\nlen1=" + len1.ToString(CultureInfo.InvariantCulture) + ", len2=" + len2.ToString(CultureInfo.InvariantCulture) + ", Switching to the memsave mode\n");
                    G.Alg = 'M';
                    G.Common.FreeIp();
                }
                bool ffttry = _nlen[m1] > clus1 && _nlen[m2] > clus2;
                if (G.UseExternalAnchors)
                {
                    var pairanch = PickPairAnch(clus1, clus2, localmem0, localmem1, mseq1, mseq2);
                    FalignGivenAnchors(pairanch, dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                }
                else if (G.ForceFft || (G.UseFft && ffttry))
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" f\b\b");
                    if (G.Alg == 'M')
                    {
                        if (l < 500 || l % 100 == 0)
                            Report("m");
                        Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, memsave: true);
                    }
                    else
                    {
                        Falign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                    }
                }
                else
                {
                    if (l < 500 || l % 100 == 0)
                        Report(" d\b\b");
                    if (G.Alg == 'M')
                    {
                        if (l < 500 || l % 100 == 0)
                            Report("m");
                        MsAlign(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, orieff1, orieff2);
                    }
                    else if (G.Alg == 'd')
                    {
                        if (clus1 == 1 && clus2 == 1)
                        {
                            G__align11(dynamicmtx, mseq1, mseq2, alloclen, G.Outgap, G.Outgap);
                        }
                        else
                        {
                            double noImp = 0.0;
                            DAlignLs(dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, ref noImp, G.Outgap, G.Outgap);
                        }
                    }
                    else if (G.Alg != 'A')
                        throw new NotSupportedException("Only the default algorithm (A) of MAFFT is available in this port.");
                    else if (clus1 == 1 && clus2 == 1)
                        G__align11(dynamicmtx, mseq1, mseq2, alloclen, G.Outgap, G.Outgap);
                    else
                        AAlign(dynamicmtx, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, orieff1, orieff2);
                }
                double nscore = IntergroupScore(mseq1, mseq2, effarr1, effarr2, clus1, clus2, mseq1[0].Length);
                if (nscore < oscore)
                {
                    for (int i = 0; i < clus1; i++)
                        mseq1[i] = oseq1[i];
                    for (int i = 0; i < clus2; i++)
                        mseq2[i] = oseq2[i];
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
            if (!_isTbfast)
                _gCap.Reset(); // disttbfast calls G__align11( NULL ... ) with A__align( NULL ... ); tbfast does not
            _dCap.Reset();
            _arAlloclen = 0;
            _impMtx = null;
            _dImpMtx = null;
        }

        /// <summary><c>fastconjuction_noname_kozo</c>: the rows of a group, their weights (normalised to 1) and the structure weights (normalised to 1, plus the weight of the sequence when positive; all 0 when no member has a structure weight).</summary>
        private static void FastConjunctionKozo(int[] memlist, string[] seq, string[] aseq, double[] peff, double[] eff, double[] peffKozo, double[] effKozo)
        {
            double total = 0.0, totalKozo = 0.0;
            for (int k = 0; k < memlist.Length; k++)
            {
                int m = memlist[k];
                aseq[k] = seq[m];
                peff[k] = eff[m];
                peffKozo[k] = effKozo[m];
                total += peff[k];
                totalKozo += peffKozo[k];
            }
            for (int m = 0; m < memlist.Length; m++)
                peff[m] /= total;
            if (totalKozo != 0.0)
            {
                for (int m = 0; m < memlist.Length; m++)
                {
                    peffKozo[m] /= totalKozo;
                    if (peffKozo[m] > 0.0)
                        peffKozo[m] += peff[m];
                }
            }
            else
            {
                for (int m = 0; m < memlist.Length; m++)
                    peffKozo[m] = 0.0;
            }
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
        internal double G__align11(double[][] dynamicmtx, string[] seq1, string[] seq2, int alloclen, int headgp, int tailgp, bool noalign = false, MafftCodonDp? codon = null)
        {
            string s1 = seq1[0], s2 = seq2[0];
            int lgth1 = s1.Length, lgth2 = s2.Length;
            if (lgth1 == 0 && lgth2 == 0)
                return 0.0;
            if (lgth1 == 0)
            {
                seq1[0] = new string(G.NewGapChar, lgth2);
                return 0.0;
            }
            if (lgth2 == 0)
            {
                seq2[0] = new string(G.NewGapChar, lgth1);
                return 0.0;
            }
            // G__align11_noalign (the distance only) has no warp
            bool trywarp = G.TryWarp && !noalign;
            if (trywarp && (headgp == 0 || tailgp == 0))
                throw new MafftProgramFailureException("At present, headgp and tailgp must be 1.\n");
            var warp = trywarp ? new WarpState(lgth1, lgth2, G.PenaltyShift) : null;
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
            int[][] ijp;
            if (noalign)
            {
                // G__align11_noalign keeps no trace-back matrix (and does not touch commonIP)
                ijp = new int[lgth1 + 1][];
                for (int i = 0; i <= lgth1; i++)
                    ijp[i] = new int[lgth2 + 1];
            }
            else
            {
                _gCap.Fit(lgth1, lgth2);
                ijp = G.Common.IpFor(_gCap.Orlgth1, _gCap.Orlgth2);
            }
            double[] currentw = w1, previousw = w2;
            int[]? codonId1 = null, codonId2 = null;
            if (codon?.Matrix != null)
            {
                // G__align11psg: the id of the codon that ends at every position (the first two positions have none)
                codonId1 = new int[lgth1 + 1];
                codonId2 = new int[lgth2 + 1];
                for (int i = 0; i < 3 && i < lgth1; i++)
                    codonId1[i] = -1;
                for (int i = 0; i < 3 && i < lgth2; i++)
                    codonId2[i] = -1;
                for (int i = 3; i < lgth1; i++)
                    codonId1[i] = MafftCodon.Id(s1, i - 2);
                for (int i = 3; i < lgth2; i++)
                    codonId2[i] = MafftCodon.Id(s2, i - 2);
            }

            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = amino[s2[0]][s1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = amino[s1[0]][s2[j]];
            if (codon?.Matrix != null && codon.GStart[codon.Offset] == 0.5 && codon.GEnd[codon.Offset] > 0.5)
                AddCodonMatch(codon.Matrix, currentw, codonId1![0], codonId2!, lgth2);
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
                if (codon?.Matrix != null && codon.GStart[codon.Offset + i] == 0.5 && codon.GEnd[codon.Offset + i] > 0.5)
                    AddCodonMatch(codon.Matrix, currentw, codonId1![i], codonId2!, lgth2);
                currentw[0] = initverticalw[i];
                mi = previousw[0];
                mpi = 0;
                double fpenaltyExI = i < lgth1 ? fpenaltyEx : 0.0;
                var ijpRow = ijp[i];
                double penaltyEnd = codon != null ? fpenalty * codon.GEnd[codon.Offset + i] : fpenalty;
                double penaltyStart = codon != null ? fpenalty * codon.GStart[codon.Offset + i - 1] : fpenalty;
                for (int j = 1; j < lastj; j++)
                {
                    wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    double g;
                    if ((g = mi + penaltyEnd) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = codon != null ? previousw[j - 1] + penaltyStart : previousw[j - 1]) >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    mi += fpenaltyExI;
                    if ((g = m[j] + penaltyEnd) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    if ((g = codon != null ? previousw[j - 1] + penaltyStart : previousw[j - 1]) >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    if (j < lgth2)
                        m[j] += fpenaltyEx;
                    if (warp != null)
                        wm = warp.Step(i, j, wm, currentw[j], fpenaltyEx, ijpRow);
                    currentw[j] += wm;
                }
                warp?.EndRow(lastj);
                lastverticalw[i] = currentw[lgth2 - 1];
            }
            (string, string, double) tracked;
            try
            {
                tracked = Atracking11(G.NewGapChar, currentw, lastverticalw, s1, s2, ijp, tailgp, warp);
            }
            catch (IndexOutOfRangeException)
            {
                throw StaleTraceBack();
            }
            var (a, b, wmo) = tracked;
            if (tailgp == 0)
                wm = wmo;
            if (alloclen < a.Length)
                throw new InvalidOperationException("LENGTH OVER!");
            seq1[0] = a;
            seq2[0] = b;
            return wm;
        }

        /// <summary><c>match_calc_mtx_codon</c>: the score of the codon pair (scaled by 600) is added to the match score of the position of the first sequence against every position of the second.</summary>
        private static void AddCodonMatch(double[][] codonmtx, double[] match, int codonid1, int[] codonSeq2, int lgth2)
        {
            for (int j = 0; j < lgth2; j++)
            {
                int codonid2 = codonSeq2[j];
                double codonmatch = 0.0;
                if (codonid1 > -1 && codonid2 > -1)
                    codonmatch = codonmtx[codonid1][codonid2] * 600.0;
                match[j] += codonmatch;
            }
        }

        private static (string, string, double) Atracking11(char gap, double[] lasthorizontalw, double[] lastverticalw, string seq1, string seq2, int[][] ijp, int tailgp, WarpState? warp = null)
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
                if (warp != null && ijp[iin][jin] >= warp.Base)
                {
                    ifi = warp.Is[ijp[iin][jin] - warp.Base];
                    jfi = warp.Js[ijp[iin][jin] - warp.Base];
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
                if (warp != null && ifi == -warp.Base && jfi == -warp.Base)
                {
                    int rest = iin;
                    while (--rest >= 0)
                    {
                        a[--pa] = seq1[rest];
                        b[--pb] = gap;
                        k++;
                    }
                    rest = jin;
                    while (--rest >= 0)
                    {
                        a[--pa] = gap;
                        b[--pb] = seq2[rest];
                        k++;
                    }
                    break;
                }
                int l = iin - ifi;
                while (--l > 0)
                {
                    a[--pa] = seq1[ifi + l];
                    b[--pb] = gap;
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    a[--pa] = gap;
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
