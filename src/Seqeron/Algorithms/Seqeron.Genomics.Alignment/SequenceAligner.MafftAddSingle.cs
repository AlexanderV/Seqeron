#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325, S138, S1541
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: addsingle (--addfragments, --addfull, --addlong, --addtotop)

    internal sealed partial class MafftGlobals
    {
        /// <summary><c>allowlongadds</c> (<c>-V</c>), <c>addtotop</c> (<c>-x</c>), <c>multidist</c> (<c>-d</c>), <c>smoothing</c> (<c>-p</c>) of <c>addsingle.c</c>; they are not the <c>smoothing</c> of <c>disttbfast</c>.</summary>
        internal bool AllowLongAdds, AddToTop, MultiDist, SingleSmoothing, DiffOut, CodonPos, CodonScore;

        /// <summary><c>arguments()</c> of <c>addsingle.c</c> for the options the script <c>mafft</c> passes to it (the script's <c>$codonposopt $codonscoreopt -Q 100 $legacygapopt -W n -O $outnum $addsinglearg $addarg $add2ndhalfarg -C n ...</c>).</summary>
        internal static MafftGlobals ParseAddsingle(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Alg = 'A', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            g.TupleSize = -1;
            g.AddProfile = true;
            g.Nadd = 0;
            g.Weight = 3;
            g.Tbrweight = 3;
            g.TreeMethod = 'X';
            g.SueffGlobal = 0.1;
            int k = 0;
            bool illegal = false;
            string Next() => k < argv.Count ? argv[k++] : "";
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-' && !illegal)
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'I': g.Nadd = CAtoi(Next()); g.Report("nadd = " + g.Nadd.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'e':
                        case 'o':
                        case 'r':
                            CAtof(Next());
                            rest = true;
                            break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); g.Report("ppenalty_ex = " + g.PpenaltyEx.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); g.Report("kappa = " + g.KimuraR.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; g.Report("blosum " + g.Nblosum.ToString(CultureInfo.InvariantCulture) + " / kimura 200\n"); rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; g.Report("jtt/kimura " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; g.Report("tm " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'l': g.FastaThreshold = CAtof(Next()); g.Constraint = 2; rest = true; break;
                        case 'c': g.ConsweightMulti = CAtof(Next()); rest = true; break;
                        case 'C': g.Nthread = CAtoi(Next()); g.Report("nthread = " + g.Nthread.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'R': g.CodonPos = true; break;
                        case 'S': g.CodonScore = true; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'K': g.AddProfile = false; break;
                        case 'y': g.DistOut = 1; break;
                        case 't': g.TreeOut = 1; break;
                        case 'T': g.NoAlign = 1; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'O': g.Outgap = 0; break;
                        case 'X': g.TreeMethod = 'X'; g.SueffGlobal = CAtof(Next()); g.Report("sueff_global = " + g.SueffGlobal.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'E': g.TreeMethod = 'E'; break;
                        case 'q': g.TreeMethod = 'q'; break;
                        case 'n': g.OutNumber = 1; break;
                        case 'H': g.Alg = 'H'; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'M': g.Alg = 'M'; break;
                        case 'N': g.NeverMemSave = 1; break;
                        case 'B': break;
                        case 'F': g.UseFft = true; break;
                        case 'G': g.ForceFft = true; g.UseFft = true; break;
                        case 'U': g.TreeIn = 1; break;
                        case 'x': g.AddToTop = true; break;
                        case 'V': g.AllowLongAdds = true; break;
                        case 'p': g.SingleSmoothing = true; break;
                        case 'u': g.Tbrweight = 0; g.Weight = 0; break;
                        case 'v': g.Tbrweight = 3; break;
                        case 'd': g.MultiDist = true; break;
                        case 'W': g.TupleSize = CAtoi(Next()); rest = true; break;
                        case 'J': g.TbUtree = false; break;
                        case 'w': g.FftWinSize = CAtoi(Next()); rest = true; break;
                        case 'L': g.LegacyGapCost = 1; break;
                        case 'Y': g.KeepLength = 1; break;
                        case 'z': g.MapOut = 2; break;
                        case 'Z': g.MapOut = 1; break;
                        case '%': g.DiffOut = true; break;
                        case ':': g.NWildcard = 1; break;
                        case '@': g.Report("illegal option @\n"); illegal = true; break; // --sp: addsingle does not know -@
                        default:
                            throw new NotSupportedException($"The MAFFT option -{c} of addsingle is not available in this port.");
                    }
                }
            }
            if (illegal)
                throw new MafftProgramFailureException(g.Log + "options: Check source file !\n");
            if (g.FastaThreshold < 0.0001)
                g.Constraint = 0;
            return g;
        }
    }

    /// <summary>The block of <c>infile.tree</c> that <c>addonetip</c> records for one added sequence (<c>Addtree</c>).</summary>
    internal sealed class MafftAddTree
    {
        internal int Nearest;
        internal double Dist1, Dist2;
        internal string Neighbors = "";
    }

    internal sealed partial class MafftRun
    {
        private const int AddSingleMaxDist = 2;

        /// <summary>What <c>addsingle</c> reads besides the sequences: the distance files of the pair stage of the script (<c>hat2</c>, or <c>hat2i</c> / <c>hat2n</c> with <c>-d</c>) and the file <c>hat3</c> (<c>-l</c>).</summary>
        internal sealed record AddSingleFiles(string? Hat2 = null, string? Hat2I = null, string? Hat2N = null, string? Hat3 = null, string? CodonPos = null, string? CodonScore = null);

        /// <summary><c>strlen</c>-free helpers of <c>addsingle.c</c>: the sequence without its gaps (<c>gappick0</c>).</summary>
        private static string GapPick0(string s) => s.Replace("-", "", StringComparison.Ordinal);

        private static string Eq2Dash(string s) => s.Replace('=', '-');

        /// <summary><c>addsingle</c>: adds the last <c>nadd</c> sequences, one by one, to the alignment that the others form (<c>--addfragments</c>, <c>--addfull</c>, <c>--addlong</c>, <c>--addtotop</c>).</summary>
        /// <param name="seq0">The sequences as the program holds them (the alignment first, the new sequences after it).</param>
        /// <param name="names">The titles as the program holds them (with the leading character of the file).</param>
        /// <param name="files">The files of the pair stage.</param>
        internal Result AddSingle(IReadOnlyList<string> seq0, IReadOnlyList<string> names, AddSingleFiles? files = null)
        {
            files ??= new AddSingleFiles();
            G.CodonPosText = files.CodonPos;
            G.CodonScoreText = files.CodonScore;
            int njob = seq0.Count;
            _njob = njob;
            int nadd = G.Nadd;
            if (G.Constraint != 0 && G.FastaThreshold < 0.0001)
                G.Constraint = 0;
            if (njob < 2)
                throw new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
            int norg = njob - nadd;
            int njobc = norg + 1;
            Report("norg = " + norg.ToString(CultureInfo.InvariantCulture) + "\n");
            Report("njobc = " + njobc.ToString(CultureInfo.InvariantCulture) + "\n");
            _nlenmax = 0;
            foreach (var s in seq0)
                _nlenmax = Math.Max(_nlenmax, s.Length);
            int alloclen = _nlenmax * 4 + 1;
            int fullseqlen = alloclen;
            var seq = new string[njob];
            for (int i = 0; i < njob; i++)
                seq[i] = seq0[i];
            var nlen = new int[njob];
            for (int i = 0; i < njob; i++)
                nlen[i] = seq[i].Length;
            _nlen = nlen;
            bool multidistOrTuple = G.MultiDist || G.TupleSize > 0;
            double[][] iscore;
            double[][]? nscore;
            if (multidistOrTuple)
            {
                iscore = HalfMtx(norg);
                nscore = new double[norg][];
                for (int i = 0; i < norg; i++)
                    nscore[i] = new double[nadd];
            }
            else
            {
                iscore = HalfMtx(njob);
                nscore = null;
            }

            MafftLocalHom[][]? localhomtable = null;
            if (G.Constraint != 0)
            {
                Report("Loading 'hat3' ... ");
                if (files.Hat3 == null)
                    throw new MafftProgramFailureException("Make hat3.\n");
                var kozoarivec = new bool[njob];
                localhomtable = ReadHat3ForAddsingle(files.Hat3, norg, nadd, G.MultiDist, kozoarivec);
                Report("\ndone.\n");
                if (kozoarivec.Any(x => x))
                    throw new NotSupportedException("addsingle with blocks of kind 'k' (structures) in hat3 is not ported yet.");
            }

            G.Constants(seq);
            if (CheckAlphabet(seq, "Illegal") is { } unknownLetter)
                throw unknownLetter;
            int alignmentlength = seq[0].Length;
            for (int i = 0; i < norg; i++)
            {
                if (alignmentlength != seq[i].Length)
                    throw new MafftProgramFailureException("#################################################################################\n" +
                        "# ERROR!                                                                        #\n" +
                        "# The original" + (njob - nadd).ToString(CultureInfo.InvariantCulture).PadLeft(4) + " sequences must be aligned                                    #\n" +
                        "#################################################################################\n");
            }
            if (G.AddProfile)
                throw new MafftProgramFailureException("Not supported!\n");

            if (G.TupleSize > 0)
            {
                KtupleDistanceMatrix(njob, norg, seq, iscore, nscore!);
                if (G.MultiDist)
                    AddSingleDndpre(norg, seq, iscore);
            }
            else if (G.MultiDist)
            {
                Report("Loading 'hat2n' (aligned sequences - new sequences) ... ");
                if (files.Hat2N == null)
                    throw new MafftProgramFailureException("Make hat2n.\n");
                ReadHat2PartText(files.Hat2N, njob, nadd, nscore!);
                Report("done.\n");
                Report("Loading 'hat2i' (aligned sequences) ... ");
                if (files.Hat2I == null)
                    throw new MafftProgramFailureException("Make hat2i.\n");
                ReadHat2DoubleHalf(files.Hat2I, norg, iscore);
                Report("done.\n");
            }
            else
            {
                Report("Loading 'hat2' ... ");
                if (files.Hat2 == null)
                    throw new MafftProgramFailureException("Make hat2.\n");
                ReadHat2DoubleHalf(files.Hat2, njob, iscore);
                Report("done.\n");
            }

            if (G.DistOut != 0)
            {
                Report("Writing distances between new sequences and existing msa.\n");
                var sb = new StringBuilder();
                for (int iadd = 0; iadd < nadd; iadd++)
                {
                    sb.Append("Distance between new sequence ").Append(iadd + 1).Append(" and ").Append(norg).Append(" sequences in existing msa\n");
                    for (int i = 0; i < norg; i++)
                    {
                        double v = multidistOrTuple ? nscore![i][iadd] : iscore[i][norg + iadd - i];
                        sb.Append(MafftGlobals.F(v, 5, 3)).Append(' ');
                        if ((i + 1) % 12 == 0)
                            sb.Append('\n');
                    }

                    sb.Append("\n\n");
                }

                DistFile = sb.ToString();
            }

            if (G.DiffOut)
                throw new NotSupportedException("addsingle -% (diffout) is marked 'Incomplete. Not used' in the program and stops with 'Not yet written'.");
            string? originalGaps = null;
            List<GapPos>[]? deletelist = null;
            if (G.KeepLength != 0)
            {
                originalGaps = RecordOriginalGaps(norg, seq);
                deletelist = new List<GapPos>[nadd];
                for (int i = 0; i < nadd; i++)
                    deletelist[i] = new List<GapPos>();
            }

            var old = seq.Take(norg).ToArray();
            CommonGapPickN(norg, old);
            for (int i = 0; i < norg; i++)
                seq[i] = old[i];
            int lenfull = seq[0].Length;
            string[]? addbk = null;
            if (G.KeepLength != 0 && G.MapOut != 0)
            {
                addbk = new string[nadd];
                for (int i = 0; i < nadd; i++)
                    addbk[i] = GapPick0(seq[norg + i]);
            }

            var newgaplistO = new int[nadd][];
            for (int i = 0; i < nadd; i++)
                newgaplistO[i] = new int[lenfull * 2];
            var newgaplistCompact = new int[lenfull * 2];
            var istherenewgap = new bool[nadd];
            var follower = new List<int>[norg];
            for (int i = 0; i < norg; i++)
                follower[i] = new List<int>();
            var follows = new int[nadd];
            var dep = new Treedep[norg];
            var topol0 = new int[norg][];
            var topol1 = new int[norg][];
            var len = new double[norg][];
            for (int i = 0; i < norg; i++)
                len[i] = new double[2];
            var iscoreo = HalfMtx(norg);
            for (int i = 0; i < norg - 1; i++)
                for (int j = i + 1; j < norg; j++)
                    iscoreo[i][j - i] = iscore[i][j - i];
            MafftAddTree[]? addtree = G.TreeOut != 0 ? Enumerable.Range(0, nadd).Select(_ => new MafftAddTree()).ToArray() : null;
            if (G.TreeIn != 0)
            {
                Report("Loading a tree ... ");
                LoadTop(norg, iscoreo, topol0, topol1, len, names, dep);
                Report("\ndone.\n\n");
            }
            else if (G.TreeOut != 0)
            {
                MuscleSupg(norg, iscoreo, topol0, topol1, len, dep, true, names.Take(norg).ToArray());
            }
            else
            {
                MuscleSupg(norg, iscoreo, topol0, topol1, len, dep, true, null, progress: false);
            }

            int[] ordertable;
            if (norg > 1)
                ordertable = topol0[norg - 2].Concat(topol1[norg - 2]).ToArray();
            else
                ordertable = new[] { 0 };

            // one thread per added sequence in turn: the program takes the next sequence of the shared counter, the first thread does all of them in order when the sequences are short (the model that the installed program shows)
            AddSingleThread(norg, njob, nadd, nlen, names, seq, iscore, nscore, localhomtable, istherenewgap, newgaplistO, follows, dep, topol0, topol1, len, addtree, deletelist, alloclen, ref _addSingleNdeleted);

            string? treeText = null;
            if (G.TreeOut != 0)
            {
                var sb = new StringBuilder(TreeFile ?? "");
                for (int i = 0; i < nadd; i++)
                {
                    sb.Append('\n');
                    sb.Append((norg + i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append(": ").Append(names[norg + i][1..]).Append('\n');
                    sb.Append("          nearest sequence: ").Append(addtree![i].Nearest + 1).Append('\n');
                    sb.Append("          approximate distance: ").Append(addtree[i].Dist1.ToString("F6", CultureInfo.InvariantCulture)).Append('\n');
                    sb.Append("          sister group: ").Append(addtree[i].Neighbors).Append('\n');
                    sb.Append("          approximate distance: ").Append(addtree[i].Dist2.ToString("F6", CultureInfo.InvariantCulture)).Append('\n');
                }

                treeText = sb.ToString();
                TreeFile = treeText;
            }

            for (int iadd = 0; iadd < nadd; iadd++)
                follower[follows[iadd]].Add(iadd);
            var order = new List<int>();
            foreach (int o in ordertable)
            {
                order.Add(o);
                foreach (int f in follower[o])
                    order.Add(f + norg);
            }

            var posmap = new int[lenfull + 2];
            var realign = new Blocktorealign[lenfull + 2];
            for (int i = 0; i < lenfull + 1; i++)
                posmap[i] = i;
            Report("\n\nCombining ..\n");
            int tmpseqlen = alloclen * 100;
            for (int iadd = 0; iadd < nadd; iadd++)
            {
                MakeGapListCompact(lenfull, posmap, newgaplistCompact, newgaplistO[iadd]);
                if (iadd == 0 || istherenewgap[iadd])
                {
                    string tmpseq1 = GapList2Alnx(lenfull, seq[0], newgaplistCompact, posmap, tmpseqlen);
                    int tmplen;
                    if ((tmplen = tmpseq1.Length) >= fullseqlen)
                    {
                        fullseqlen = tmplen * 2 + 1;
                        Report("reallocating...");
                        Report("done.\n");
                    }

                    seq[0] = tmpseq1;
                    int ien = norg + iadd;
                    if (G.Nthread > 0 && ien > 500)
                    {
                        Report(iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + " (threads 0-" + G.Nthread.ToString(CultureInfo.InvariantCulture) + ")\r");
                        for (int i = 1; i < ien; i++)
                            seq[i] = GapList2Alnx(lenfull, seq[i], newgaplistCompact, posmap, tmpseqlen);
                    }
                    else
                    {
                        Report(iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + "\r");
                        for (int i = 1; i < ien; i++)
                        {
                            if (i == 1)
                                Report(" " + iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + "\r");
                            seq[i] = GapList2Alnx(lenfull, seq[i], newgaplistCompact, posmap, tmpseqlen);
                        }
                    }
                }

                string merged = InsertGapsByOtherFragmentsCompact(lenfull, seq[norg + iadd], newgaplistO[iadd], posmap);
                seq[norg + iadd] = Eq2Dash(merged);
                AdjustPosMap(lenfull, posmap, newgaplistCompact);
                CountNewRes(lenfull, realign, posmap, newgaplistO[iadd]);
            }

            Report("\r   done.                      \n\n");
            if (G.SingleSmoothing)
            {
                SmoothingOne(njob, nadd, lenfull, seq, realign);
                SmoothingTwo(njob, nadd, lenfull, seq, realign);
            }

            for (int i = 0; i < lenfull + 1; i++)
            {
                if (realign[i].NNewRes > 1)
                {
                    Report("\rRealigning " + i.ToString(CultureInfo.InvariantCulture) + "/" + lenfull.ToString(CultureInfo.InvariantCulture) + "           \r");
                    int zure = DoRealignmentTree(realign[i], seq, ref fullseqlen, norg, topol0, topol1, follows, njob);
                    for (int j = i + 1; j < lenfull + 1; j++)
                    {
                        if (realign[j].NNewRes != 0)
                        {
                            realign[j].Start -= zure;
                            realign[j].End -= zure;
                        }
                    }
                }
            }

            Report("\r   done.                      \n\n");
            List<(int Row, int Pos, int Len)>? deleted = null;
            string? deleteMap = null;
            string?[]? newNames = null;
            if (G.KeepLength != 0)
            {
                RestoreOriginalGaps(njob, seq, originalGaps!);
                deleted = new List<(int Row, int Pos, int Len)>();
                for (int i = 0; i < nadd; i++)
                    foreach (var d in deletelist![i])
                        deleted.Add((njob - nadd + i, d.Pos, d.Len));
                if (G.MapOut != 0)
                {
                    var addNames = names.Skip(njob - nadd).ToArray();
                    var realn = seq.Skip(njob - nadd).ToArray();
                    if (G.MapOut == 1)
                    {
                        deleteMap = ReconstructDeleteMap(nadd, addbk!, deletelist!, realn, addNames);
                    }
                    else
                    {
                        var changed = new string?[nadd];
                        deleteMap = ReconstructDeleteMapCompact(nadd, addbk!, deletelist!, realn, addNames, changed!);
                        newNames = changed;
                    }
                }
            }

            FreeStatics();
            Report("addsingle (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            if (_addSingleNdeleted > 0)
            {
                Report("\nTo keep the alignment length, " + _addSingleNdeleted.ToString(CultureInfo.InvariantCulture) + " letters were DELETED.\n");
                Report(G.MapOut != 0 ? "The deleted letters are shown in the (filename).map file.\n" : "To know the positions of deleted letters, rerun the same command with the --mapout option.\n");
            }

            return new Result(seq, order.ToArray()) { Deleted = deleted, DeleteMap = deleteMap, NewNames = newNames };
        }

        private int _addSingleNdeleted;

        private static double[][] HalfMtx(int n)
        {
            var m = new double[n][];
            for (int i = 0; i < n; i++)
                m[i] = new double[n - i];
            return m;
        }

        /// <summary><c>loadtop</c> (<c>addsingle -U</c>): the guide tree of the old sequences from the file <c>_guidetree</c> (one line per node: two sequence numbers and two branch lengths, the groups keep their full membership); a branch length of -1 is
        /// replaced by half the mean distance of the two groups minus the height of the node. The tree goes to <c>infile.tree</c> ("#by loadtop").</summary>
        private void LoadTop(int nseq, double[][] mtx, int[][] topol0, int[][] topol1, double[][] len, IReadOnlyList<string> names, Treedep[] dep)
        {
            if (GuideTreeText == null)
                throw new MafftProgramFailureException("cannot open _guidetree\n");
            var fp = new FgetsReader(GuideTreeText);
            var hist = new int[nseq];
            var nmemar = new int[nseq];
            var height = new double[nseq];
            var tree = TreeLeafStrings(nseq, names);
            string? treetmp = null;
            for (int i = 0; i < nseq; i++)
            {
                hist[i] = -1;
                nmemar[i] = 1;
            }

            Report("\n");
            for (int k = 0; k < nseq - 1; k++)
            {
                if (k % 10 == 0)
                    Report("\r" + MafftGlobals.D5(k) + " / " + nseq.ToString(CultureInfo.InvariantCulture));
                len[k][0] = len[k][1] = -1.0;
                var (im, jm) = LoadTreeOneLine(fp, len[k]);
                if (im > nseq - 1 || jm > nseq - 1 || im >= 0 && tree[im] == null || jm >= 0 && tree[jm] == null)
                    throw new MafftProgramFailureException("\n\nCheck the guide tree.\nim=" + (im + 1).ToString(CultureInfo.InvariantCulture) + ", jm=" + (jm + 1).ToString(CultureInfo.InvariantCulture)
                        + "\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n\n");
                if (im < 0)
                    throw UndefinedBehaviour("a node names the sequence 0: the program reads hist[-1]");
                int prevnode = hist[im];
                dep[k].Child0 = prevnode;
                int nmemim = nmemar[im];
                topol0[k] = GroupMembers(prevnode, im, topol0, topol1, true);
                int nmemjm = nmemar[jm];
                prevnode = hist[jm];
                dep[k].Child1 = prevnode;
                topol1[k] = GroupMembers(prevnode, jm, topol0, topol1, true);
                hist[im] = k;
                nmemar[im] = nmemim + nmemjm;
                if (len[k][0] == -1 || len[k][1] == -1)
                {
                    Report("Re-computing the length of branch " + k.ToString(CultureInfo.InvariantCulture) + "..\n");
                    double clusterdist = 0.0;
                    int mpair = 0;
                    foreach (int mi in topol0[k])
                    {
                        foreach (int mj in topol1[k])
                        {
                            int minijm = Math.Min(mi, mj), maxijm = Math.Max(mi, mj);
                            clusterdist += mtx[minijm][maxijm - minijm];
                            mpair++;
                        }
                    }

                    clusterdist /= mpair;
                    Report("clusterdist = " + clusterdist.ToString("F6", CultureInfo.InvariantCulture) + "\n");
                    if (len[k][0] == -1)
                        len[k][0] = clusterdist / 2.0 - height[im];
                    if (len[k][1] == -1)
                        len[k][1] = clusterdist / 2.0 - height[im];
                    Report("len0 = " + len[k][0].ToString("F6", CultureInfo.InvariantCulture) + "\n");
                    Report("len1 = " + len[k][1].ToString("F6", CultureInfo.InvariantCulture) + "\n\n");
                }

                height[im] += len[k][0];
                dep[k].DistFromTip = height[im];
                treetmp = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                tree[im] = treetmp;
                tree[jm] = null!;
            }

            TreeFile = (treetmp ?? "(null)") + ";\n#by loadtop\n";
        }

        // -------------------------------------------------------------------- the distances

        /// <summary><c>ktupledistancematrix</c> of <c>addsingle.c</c>: the 6-mer (10-mer) distances of the old sequences among themselves (<paramref name="imtx"/>) and to the new ones (<paramref name="nmtx"/>); the new sequences are not compared with each other.</summary>
        private void KtupleDistanceMatrix(int nseq, int norg, string[] seq, double[][] imtx, double[][] nmtx)
        {
            bool dna = G.Dorp == 'd';
            int tuple = G.TupleSize;
            Report("\n\nMaking a distance matrix ..\n");
            int tsize = dna ? (int)Math.Pow(4, tuple) : 46656;
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

            var nogaplen = new int[nseq];
            var pointt = new int[nseq][];
            int nunknown = 0;
            int[] mult;
            int baseN;
            if (dna)
            {
                if (tuple is not (6 or 10))
                    throw new MafftProgramFailureException("tuplesize=" + tuple.ToString(CultureInfo.InvariantCulture) + ": not supported\n");
                mult = new int[tuple];
                int v = 1;
                for (int i = tuple - 1; i >= 0; i--)
                {
                    mult[i] = v;
                    v *= 4;
                }

                baseN = 4;
            }
            else
            {
                mult = new[] { 7776, 1296, 216, 36, 6, 1 };
                baseN = 6;
                tuple = 6;
            }

            for (int i = 0; i < nseq; i++)
            {
                string tmpseq = GapPick0(seq[i]);
                nogaplen[i] = tmpseq.Length;
                var grp = new List<int>(tmpseq.Length);
                int limit = dna ? 4 : 6;
                foreach (char c in tmpseq)
                {
                    int tmp = G.AminoGrp[c];
                    if (tmp < limit)
                        grp.Add(tmp);
                    else
                        nunknown++;
                }

                pointt[i] = grp.Count < tuple ? Array.Empty<int>() : MakePointTable(grp.ToArray(), grp.Count, tuple, mult, baseN);
            }

            if (nunknown != 0)
                Report("\nThere are " + nunknown.ToString(CultureInfo.InvariantCulture) + " ambiguous characters\n");
            var memo = new int[tsize];
            var selfscore = new double[nseq];
            for (int i = 0; i < nseq; i++)
            {
                var table1 = new int[tsize];
                foreach (int point in pointt[i])
                    table1[point]++;
                selfscore[i] = CommonSextet(table1, pointt[i], memo);
            }

            for (int i = 0; i < norg; i++)
            {
                var table1 = new int[tsize];
                if (i % 100 == 0)
                    Report("\r" + MafftGlobals.D5(i + 1) + " / " + norg.ToString(CultureInfo.InvariantCulture) + (G.Nthread > 0 ? " (thread " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(4) + ")" : ""));
                foreach (int point in pointt[i])
                    table1[point]++;
                for (int j = i + 1; j < nseq; j++)
                {
                    double mtxv = CommonSextet(table1, pointt[j], memo);
                    double longer, shorter;
                    if (nogaplen[i] > nogaplen[j])
                    {
                        longer = nogaplen[i];
                        shorter = nogaplen[j];
                    }
                    else
                    {
                        longer = nogaplen[j];
                        shorter = nogaplen[i];
                    }

                    double lenfac = 1.0 / (shorter / longer * lenfacd + lenfacb / (longer + lenfacc) + lenfaca);
                    double bunbo = Math.Min(selfscore[i], selfscore[j]);
                    double d = bunbo == 0.0 ? AddSingleMaxDist : AddSingleMaxDist * (1.0 - mtxv / bunbo) * lenfac;
                    if (j < norg)
                        imtx[i][j - i] = d;
                    else
                        nmtx[i][j - norg] = d;
                }
            }

            Report("\ndone.\n\n");
        }

        /// <summary><c>dndpre</c> of <c>addsingle.c</c> (used with <c>-d -W n</c> only): the distances of the old sequences from the alignment (the naive score of the pairs).</summary>
        private void AddSingleDndpre(int nseq, string[] seq, double[][] mtx)
        {
            var selfscore = new double[nseq];
            for (int i = 0; i < nseq; i++)
                selfscore[i] = NaivePairScore11(seq[i], seq[i], 0);
            if (G.Nthread > 0)
            {
                // the threads take the pairs in turn; the first one does them all, and each of the others that arrives after the last pair prints the progress of the last row once more
                if (nseq == 1)
                    return;
                int i = 0, j = 0;
                int thread = 0;
                while (true)
                {
                    j++;
                    if (j == nseq)
                    {
                        Report(D4(i + 1) + "/" + D4(nseq) + " (thread " + D4(thread) + "), dndpre\r");
                        i++;
                        j = i + 1;
                        if (i == nseq - 1)
                            break;
                    }

                    AddSingleDndpreDistance(i, j, seq, selfscore, mtx);
                }

                for (int t = 1; t < G.Nthread; t++)
                    Report(D4(nseq - 1) + "/" + D4(nseq) + " (thread " + D4(t) + "), dndpre\r");
                return;
            }

            for (int i = 0; i < nseq - 1; i++)
            {
                Report(D4(i + 1) + "/" + D4(nseq) + "\r");
                for (int j = i + 1; j < nseq; j++)
                    AddSingleDndpreDistance(i, j, seq, selfscore, mtx);
            }
        }

        private static string D4(int v) => v.ToString(CultureInfo.InvariantCulture).PadLeft(4);

        private void AddSingleDndpreDistance(int i, int j, string[] seq, double[] selfscore, double[][] mtx)
        {
            double bunbo = Math.Min(selfscore[i], selfscore[j]);
            double mtxv = bunbo == 0.0 ? AddSingleMaxDist : AddSingleMaxDist * (1.0 - NaivePairScore11(seq[i], seq[j], G.Penalty * 10) / bunbo);
            if (mtxv < 0.0)
            {
                Report("WARNING: distance " + i.ToString(CultureInfo.InvariantCulture) + "-" + j.ToString(CultureInfo.InvariantCulture) + " is strange, " + mtxv.ToString("F6", CultureInfo.InvariantCulture) + ".\n");
                mtxv = 0.0;
            }

            if (mtxv > 9.9)
            {
                Report("WARNING: distance " + i.ToString(CultureInfo.InvariantCulture) + "-" + j.ToString(CultureInfo.InvariantCulture) + " is strange, " + mtxv.ToString("F6", CultureInfo.InvariantCulture) + ".\n");
                mtxv = 9.9;
            }

            mtx[i][j - i] = mtxv;
        }

        private static void ReadHat2DoubleHalf(string text, int n, double[][] mtx)
        {
            var full = ReadHat2Text(text, n);
            for (int i = 0; i < mtx.Length; i++)
                for (int j = i + 1; j < Math.Min(n, mtx[i].Length + i); j++)
                    mtx[i][j - i] = full[i][j];
        }

        // -------------------------------------------------------------------- adding one sequence

        /// <summary><c>addsinglethread</c>: for every new sequence the place in the guide tree of the old sequences (<c>addonetip</c>), the alignment with its sister group (<c>treebase</c>) and the list of the new columns that it made.</summary>
        private void AddSingleThread(int norg, int njob, int nadd, int[] nlen, IReadOnlyList<string> name, string[] seq, double[][] iscore, double[][]? nscore, MafftLocalHom[][]? localhomtable, bool[] istherenewgap, int[][] newgaplist, int[] follows, Treedep[] dep,
            int[][] topol0, int[][] topol1, double[][] len, MafftAddTree[]? addtree, List<GapPos>[]? deletelist, int alloclen, ref int ndeleted)
        {
            int njobc = norg + 1;
            bool multidistOrTuple = G.MultiDist || G.TupleSize > 0;
            var nlenc = new int[njobc];
            for (int i = 0; i < norg; i++)
                nlenc[i] = nlen[i];
            var nogaplenjusttodecideaddhereornot = new int[njobc];
            if (!G.AllowLongAdds)
                for (int i = 0; i < norg; i++)
                    nogaplenjusttodecideaddhereornot[i] = GapPick0(seq[i]).Length;
            MafftLocalHom[][]? localhomtablec = null;
            if (G.Constraint != 0)
            {
                localhomtablec = new MafftLocalHom[njobc][];
                for (int i = 0; i < njobc; i++)
                {
                    localhomtablec[i] = new MafftLocalHom[G.MultiDist ? 1 : njobc];
                    for (int j = 0; j < localhomtablec[i].Length; j++)
                        localhomtablec[i][j] = new MafftLocalHom();
                }

                if (!G.MultiDist)
                    for (int i = 0; i < norg; i++)
                        for (int j = 0; j < norg; j++)
                            localhomtablec[i][j] = CloneHead(localhomtable![i][j]);
            }

            var bseq = new string[njobc];
            var alnleninnode = new int[norg];
            for (int i = 0; i < norg; i++)
                bseq[i] = seq[i];
            if (norg == 1)
            {
                alnleninnode[0] = bseq[0].Length;
            }
            else
            {
                for (int i = norg - 2; i >= 0; i--)
                {
                    var rows = new List<string>();
                    foreach (int m in topol0[i])
                        rows.Add(bseq[m]);
                    foreach (int m in topol1[i])
                        rows.Add(bseq[m]);
                    var arr = rows.ToArray();
                    CommonGapPickN(arr.Length, arr);
                    alnleninnode[i] = arr[0].Length;
                }
            }

            for (int iadd = 0; iadd < nadd; iadd++)
            {
                if (G.Nthread > 0)
                {
                    if (iadd < 500)
                        Report("\rSTEP " + iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + " (thread 0)                    \r");
                    else if (iadd % 100 == 0)
                        Report("\nSTEP " + iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + " (thread 0)                    \n");
                }
                else
                {
                    if (iadd < 500)
                        Report("\rSTEP " + iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + "                    \r");
                    else if (iadd % 100 == 0)
                        Report("\nSTEP " + iadd.ToString(CultureInfo.InvariantCulture) + " / " + nadd.ToString(CultureInfo.InvariantCulture) + "                    \n");
                }

                string targetseq = seq[norg + iadd];
                for (int i = 0; i < norg; i++)
                    bseq[i] = seq[i];
                bseq[norg] = GapPick0(targetseq);
                nogaplenjusttodecideaddhereornot[norg] = G.AllowLongAdds ? 0 : bseq[norg].Length;

                var iscorec = HalfMtx(njobc);
                for (int i = 0; i < norg - 1; i++)
                    for (int j = i + 1; j < norg; j++)
                        iscorec[i][j - i] = iscore[i][j - i];
                for (int i = 0; i < norg; i++)
                    iscorec[i][norg - i] = multidistOrTuple ? nscore![i][iadd] : iscore[i][norg + iadd - i];
                nlenc[norg] = nlen[norg + iadd];
                if (G.Constraint != 0)
                {
                    for (int i = 0; i < norg; i++)
                    {
                        if (G.MultiDist)
                        {
                            localhomtablec![i][0] = CloneHead(localhomtable![i][iadd]);
                        }
                        else
                        {
                            localhomtablec![i][norg] = CloneHead(localhomtable![i][norg + iadd]);
                            localhomtablec[norg][i] = CloneHead(localhomtable[norg + iadd][i]);
                        }
                    }
                }


                var topolc0 = new int[njobc][];
                var topolc1 = new int[njobc][];
                var lenc = new double[njobc][];
                for (int i = 0; i < njobc; i++)
                    lenc[i] = new double[2];
                int neighbor = G.AddToTop
                    ? AddOneTip(njobc, topolc0, topolc1, lenc, iscorec, topol0, topol1, len, dep, addtree, iadd, alnleninnode, nogaplenjusttodecideaddhereornot, true)
                    : AddOneTip(njobc, topolc0, topolc1, lenc, iscorec, topol0, topol1, len, dep, addtree, iadd, alnleninnode, nogaplenjusttodecideaddhereornot, false);
                var effc = new double[njobc];
                if (G.Tbrweight != 0)
                {
                    G.Weight = 3;
                    CountEffSimple(njobc, topolc0, topolc1, lenc, effc);
                }
                else
                {
                    for (int i = 0; i < njobc; i++)
                        effc[i] = 1.0;
                }

                if (G.NoAlign != 0)
                    continue;
                follows[iadd] = neighbor;
                var mergeoralign = new char[njobc];
                for (int i = 0; i < njobc - 1; i++)
                    mergeoralign[i] = 'n';
                int newid = njobc - 1;
                for (int i = 0; i < njobc - 1; i++)
                {
                    if (topolc0[i].Length == 1 && topolc0[i][0] == newid)
                        mergeoralign[i] = mergeoralign[i] != 'n' ? 'w' : '1';
                    else if (topolc1[i].Length == 1 && topolc1[i][0] == newid)
                        mergeoralign[i] = mergeoralign[i] != 'n' ? 'w' : '2';
                }

                for (int i = 0; i < njobc - 1; i++)
                {
                    bool in0 = topolc0[i].All(m => m == newid);
                    bool in1 = topolc1[i].All(m => m == newid);
                    if (in0 && in1)
                        mergeoralign[i] = 'w';
                    else if (in0)
                        mergeoralign[i] = '1';
                    else if (in1)
                        mergeoralign[i] = '2';
                }

                int allocBefore = alloclen;
                AddSingleTreeBase(njobc, nlenc, bseq, mergeoralign, topolc0, topolc1, effc, ref allocBefore, localhomtablec);
                if (G.KeepLength != 0)
                {
                    var oldRows = bseq.Take(norg).ToArray();
                    var addRow = new[] { bseq[norg] };
                    int lenBefore = oldRows[0].Length;
                    var deleted = DeleteNewInsertions(norg, 1, oldRows, addRow, true);
                    deletelist![iadd] = deleted[0];
                    for (int i = 0; i < norg; i++)
                        bseq[i] = oldRows[i];
                    bseq[norg] = addRow[0];
                    ndeleted += lenBefore - oldRows[0].Length;
                    targetseq = bseq[norg];
                    istherenewgap[iadd] = false;
                }
                else
                {
                    targetseq = bseq[norg];
                    int rep = -1;
                    int i;
                    for (i = 0; i < norg; i++)
                        if (bseq[i].Contains('='))
                            break;
                    if (i == norg)
                    {
                        istherenewgap[iadd] = false;
                    }
                    else
                    {
                        rep = i;
                        istherenewgap[iadd] = true;
                        MakeNewGapList(newgaplist[iadd], bseq[rep]);
                    }
                }

                seq[norg + iadd] = Eq2Dash(targetseq);
            }
        }

        /// <summary><c>makenewgaplist</c>: the number of inserted ('=') columns in front of every column of the old alignment, ended by -1.</summary>
        private static void MakeNewGapList(int[] l, string a)
        {
            int li = 0, ai = 0;
            while (true)
            {
                while (ai < a.Length && a[ai] == '=')
                {
                    ai++;
                    l[li]++;
                }

                l[++li] = 0;
                if (ai >= a.Length)
                    break;
                ai++;
            }

            l[li] = -1;
        }

        /// <summary><c>addonetip</c> / <c>addonetip2top</c> (<c>mltaln9.c</c>): the guide tree of the old sequences plus the new one, which joins the tree at the height of its distance to the nearest old sequence (<c>addonetip2top</c>: at the first node whose height is not below half the distance 0, i.e. above the first step).</summary>
        private int AddOneTip(int njobc, int[][] topolc0, int[][] topolc1, double[][] lenc, double[][] iscorec, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, MafftAddTree[]? addtree, int iadd, int[] alnleninnode, int[] nogaplen, bool toTop)
        {
            int norg = njobc - 1;
            int nstep = njobc - 2;
            double sueff1 = 1.0 - G.SueffGlobal;
            double sueff05 = G.SueffGlobal * 0.5;
            var leaf2node = new int[norg];
            var alive = new bool[norg];
            for (int i = 0; i < norg; i++)
            {
                alive[i] = true;
                leaf2node[i] = -1;
            }

            double minscore;
            int nearest;
            if (toTop)
            {
                nearest = 0;
                minscore = 0.0;
            }
            else
            {
                minscore = 9999.9;
                nearest = -1;
                for (int i = 0; i < norg; i++)
                {
                    double tmpmin = iscorec[i][norg - i];
                    if (minscore > tmpmin)
                    {
                        minscore = tmpmin;
                        nearest = i;
                    }
                }
            }

            int nearesto = nearest;
            double minscoreo = minscore;
            int posinnew = 0;
            int repnorg = -1;
            int nogaplentoadd = nogaplen[norg];
            double addedlen = 0.0;
            int neighbor = 0;
            int i0;
            for (i0 = 0; i0 < nstep; i0++)
            {
                int mem0 = topol0[i0][0];
                int mem1 = topol1[i0][0];
                int nearestnode = leaf2node[nearest];
                int reflen = nearestnode == -1 ? nogaplen[nearest] : alnleninnode[nearestnode];
                bool seqlengthcondition = G.NoAlign != 0 || nogaplentoadd <= reflen;
                bool take = toTop ? repnorg == -1 && dep[i0].DistFromTip * 2 >= minscore : repnorg == -1 && dep[i0].DistFromTip * 2 > minscore && seqlengthcondition;
                if (take)
                {
                    if (nearestnode == -1)
                    {
                        topolc0[posinnew] = new[] { nearest };
                        addedlen = lenc[posinnew][0] = minscore / 2;
                    }
                    else
                    {
                        topolc0[posinnew] = topol0[nearestnode].Concat(topol1[nearestnode]).ToArray();
                        addedlen = lenc[posinnew][0] = dep[i0].DistFromTip - minscore / 2;
                    }

                    neighbor = topolc0[posinnew][^1];
                    if (G.TreeOut != 0)
                        RecordAddTree(addtree![iadd], nearesto, minscoreo, minscore, topolc0[posinnew]);
                    topolc1[posinnew] = new[] { norg };
                    lenc[posinnew][1] = minscore / 2;
                    repnorg = nearest;
                    posinnew++;
                }

                double eff0 = iscorec[mem0][norg - mem0];
                double eff1 = iscorec[mem1][norg - mem1];
                iscorec[mem0][norg - mem0] = Math.Min(eff0, eff1) * sueff1 + (eff0 + eff1) * sueff05;
                iscorec[mem1][norg - mem1] = 9999.9;
                alive[mem1] = false;
                if (nearest == mem1 || nearest == mem0)
                {
                    minscore = 9999.9;
                    for (int j = 0; j < norg; j++)
                    {
                        if (!alive[j])
                            continue;
                        double tmpmin = iscorec[j][norg - j];
                        if (minscore > tmpmin)
                        {
                            minscore = tmpmin;
                            nearest = j;
                        }
                    }
                }

                if (topol0[i0][0] == repnorg)
                {
                    topolc0[posinnew] = topol0[i0].Concat(new[] { norg }).ToArray();
                    lenc[posinnew][0] = len[i0][0] - addedlen;
                    addedlen = 0.0;
                }
                else
                {
                    topolc0[posinnew] = topol0[i0].ToArray();
                    lenc[posinnew][0] = len[i0][0];
                }

                if (topol1[i0][0] == repnorg)
                {
                    topolc1[posinnew] = topol1[i0].Concat(new[] { norg }).ToArray();
                    lenc[posinnew][1] = len[i0][1] - addedlen;
                    addedlen = 0.0;
                    repnorg = topolc0[posinnew][0];
                }
                else
                {
                    topolc1[posinnew] = topol1[i0].ToArray();
                    lenc[posinnew][1] = len[i0][1];
                }

                foreach (int m in topol0[i0])
                    leaf2node[m] = i0;
                foreach (int m in topol1[i0])
                    leaf2node[m] = i0;
                posinnew++;
            }

            int[] topolo0, topolo1;
            if (nstep > 0)
            {
                topolo0 = topol0[nstep - 1];
                topolo1 = topol1[nstep - 1];
            }
            else
            {
                topolo0 = new[] { 0 };
                topolo1 = Array.Empty<int>();
            }

            if (repnorg == -1)
            {
                topolc0[posinnew] = topolo0.Concat(topolo1).ToArray();
                lenc[posinnew][0] = nstep > 0 ? minscore / 2 - dep[nstep - 1].DistFromTip : minscore / 2;
                topolc1[posinnew] = new[] { norg };
                lenc[posinnew][1] = minscore / 2;
                neighbor = norg - 1;
                if (G.TreeOut != 0)
                    RecordAddTree(addtree![iadd], nearesto, minscoreo, minscore, topolc0[posinnew]);
            }

            return neighbor;
        }

        private static void RecordAddTree(MafftAddTree t, int nearesto, double minscoreo, double minscore, int[] group)
        {
            t.Nearest = nearesto;
            t.Dist1 = minscoreo;
            t.Dist2 = minscore;
            var sb = new StringBuilder();
            foreach (int m in group)
                sb.Append((m + 1).ToString(CultureInfo.InvariantCulture)).Append(' ');
            t.Neighbors = sb.ToString();
        }

        /// <summary><c>treebase</c> of <c>addsingle.c</c>: only the step at which the new sequence joins is aligned (<c>mergeoralign</c> '1' / '2'); the other steps join old sequences that are aligned already.</summary>
        private void AddSingleTreeBase(int nseq, int[] nlen, string[] aseq, char[] mergeoralign, int[][] topol0, int[][] topol1, double[] effarr, ref int alloclen, MafftLocalHom[][]? localhomtable)
        {
            var fftlog = Enumerable.Repeat(1, nseq).ToArray();
            var alreadyaligned = new bool[nseq];
            for (int i = 0; i < nseq - 1; i++)
                alreadyaligned[i] = true;
            alreadyaligned[nseq - 1] = false;
            var add = new MafftAddState(mergeoralign, alreadyaligned);
            if (G.Constraint != 0)
            {
                if (G.MultiDist)
                    DontCalcImportanceFirstOne(nseq, localhomtable!);
                else
                    CalcImportanceFull(nseq, effarr, aseq, localhomtable!);
            }

            for (int l = 0; l < nseq - 1; l++)
            {
                if (mergeoralign[l] == 'n')
                    continue;
                int m1 = topol0[l][0];
                int m2 = topol1[l][0];
                int len1 = aseq[m1].Length;
                int len2 = aseq[m2].Length;
                if (alloclen < len1 + len2)
                    throw new MafftProgramFailureException("Length over!\n");
                int clus1 = topol0[l].Length, clus2 = topol1[l].Length;
                var mseq1 = new string[clus1];
                var mseq2 = new string[clus2];
                var effarr1 = new double[clus1];
                var effarr2 = new double[clus2];
                FastConjunction(topol0[l], aseq, mseq1, effarr1, effarr, 0.0);
                FastConjunction(topol1[l], aseq, mseq2, effarr2, effarr, 0.0);
                char mode = mergeoralign[l];
                G.NewGapChar = mode is '1' or '2' ? '=' : '-';
                int len1nocommongap = len1;
                int len2nocommongap = len2;
                int[]? gapmap = null;
                if (mode == '1')
                {
                    gapmap = FindCommonGaps(clus2, mseq2);
                    CommonGapPickN(clus2, mseq2);
                    len2nocommongap = mseq2[0].Length;
                }
                else if (mode == '2')
                {
                    gapmap = FindCommonGaps(clus1, mseq1);
                    CommonGapPickN(clus1, mseq1);
                    len1nocommongap = mseq1[0].Length;
                }

                if (G.NeverMemSave == 0 && G.Constraint != 2 && G.Alg != 'M' && (len1 > 50000 || len2 > 50000))
                {
                    Report("\nlen1=" + len1.ToString(CultureInfo.InvariantCulture) + ", len2=" + len2.ToString(CultureInfo.InvariantCulture) + ", Switching to the memsave mode.\n");
                    G.Alg = 'M';
                    G.Common.FreeIp();
                }

                if (G.Alg == 'M')
                    G.Common.FreeIp();
                bool ffttry = fftlog[m1] != 0 && fftlog[m2] != 0 && nlen[m1] > clus1 && nlen[m2] > clus2;
                var dyn = G.NDisConsweightMulti;
                if (G.Constraint == 2)
                {
                    if (G.CodonPos)
                        throw new MafftProgramFailureException("\n\nThe --codonpos option is supported only for --6merpair --addfragments at this point. Add the --6merpair flag, for now.\n\n\n");
                    if (G.Alg == 'M')
                        throw new MafftProgramFailureException("\n\nMemory saving mode is not supported.\n\n");
                    Report("c");
                    if (G.Alg == 'A')
                    {
                        var shrink = G.MultiDist ? ShrinkLocalHomOne(topol0[l], topol1[l], nseq - 1, localhomtable!) : ShrinkLocalHomFull(topol0[l], topol1[l], localhomtable!);
                        var imp = FillImp(clus1, clus2, mseq1[0].Length, mseq2[0].Length, mseq1, mseq2, effarr1, effarr2, shrink, topol0[l], topol1[l]);
                        _ = AAlign(dyn, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, 1.0, 1.0, imp);
                    }
                    else if (G.Alg == 'Q')
                    {
                        throw new MafftProgramFailureException("Q has been disabled.\n");
                    }
                }
                else if (G.ForceFft || (G.UseFft && ffttry))
                {
                    Report("f");
                    if (G.Alg == 'M')
                    {
                        Report("m");
                        _ = Falign(dyn, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, memsave: true);
                    }
                    else
                    {
                        _ = Falign(dyn, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                    }

                    fftlog[m1] = FftlogResult;
                }
                else
                {
                    Report("d");
                    fftlog[m1] = 0;
                    switch (G.Alg)
                    {
                        case 'M':
                            Report("m");
                            _ = MsAlign(dyn, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, 1.0, 1.0);
                            break;
                        case 'A':
                            _ = AAlign(dyn, G.Penalty, G.PenaltyEx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, null, null, null, null, G.Outgap, G.Outgap, null, -1, -1, 0, 1.0, 1.0);
                            break;
                        default:
                            throw new MafftProgramFailureException("ERROR IN SOURCE FILE\n");
                    }
                }

                nlen[m1] = (int)(0.5 * (nlen[m1] + nlen[m2]));
                if (mode is '1' or '2')
                    MergeAfterAlignment(aseq, nseq, add, mode, gapmap, len1, len1nocommongap, clus1, clus2, mseq1, mseq2, topol0[l], topol1[l], alloclen, true, len2, len2nocommongap);
                else
                {
                    for (int k = 0; k < clus1; k++)
                        aseq[topol0[l][k]] = mseq1[k];
                    for (int k = 0; k < clus2; k++)
                        aseq[topol1[l][k]] = mseq2[k];
                }
            }

            G.NewGapChar = '-';
        }

        private static MafftLocalHom CloneHead(MafftLocalHom h) => new()
        {
            Next = h.Next, Start1 = h.Start1, End1 = h.End1, Start2 = h.Start2, End2 = h.End2, Opt = h.Opt, Overlapaa = h.Overlapaa, Importance = h.Importance, Rimportance = h.Rimportance, Korh = h.Korh,
        };

        /// <summary><c>fastshrinklocalhom</c>: the blocks of every pair (member of group 1, member of group 2) of the full table; null when the pair has none.</summary>
        private static MafftLocalHom?[][] ShrinkLocalHomFull(int[] mem1, int[] mem2, MafftLocalHom[][] localhom)
        {
            var shrink = new MafftLocalHom?[mem1.Length][];
            for (int k1 = 0; k1 < mem1.Length; k1++)
            {
                shrink[k1] = new MafftLocalHom?[mem2.Length];
                for (int k2 = 0; k2 < mem2.Length; k2++)
                    shrink[k1][k2] = localhom[mem1[k1]][mem2[k2]].Opt == -1 ? null : localhom[mem1[k1]][mem2[k2]];
            }

            return shrink;
        }

        /// <summary><c>fastshrinklocalhom_one</c>: the second group is the new sequence (number <paramref name="norg"/>); its column of the table is column 0.</summary>
        private static MafftLocalHom?[][] ShrinkLocalHomOne(int[] mem1, int[] mem2, int norg, MafftLocalHom[][] localhom)
        {
            var shrink = new MafftLocalHom?[mem1.Length][];
            for (int k1 = 0; k1 < mem1.Length; k1++)
            {
                shrink[k1] = new MafftLocalHom?[mem2.Length];
                for (int k2 = 0; k2 < mem2.Length; k2++)
                {
                    if (mem2[k2] != norg)
                        throw new MafftProgramFailureException("ERROR! *intpt2 = " + mem2[k2].ToString(CultureInfo.InvariantCulture) + "\n");
                    shrink[k1][k2] = localhom[mem1[k1]][0].Opt == -1 ? null : localhom[mem1[k1]][0];
                }
            }

            return shrink;
        }

        /// <summary><c>dontcalcimportance_firstone</c>: the importance of a block is half its score.</summary>
        private static void DontCalcImportanceFirstOne(int nseq, MafftLocalHom[][] localhom)
        {
            for (int i = 0; i < nseq - 1; i++)
                for (var t = localhom[i][0]; t != null; t = t.Next)
                    t.Importance = t.Opt * 0.5;
        }

        /// <summary><c>calcimportance</c> (the full table: both directions of every pair): the importance of a block = the mean over its residues of the weighted number of sequences that have a block there, times its score; the two directions are averaged.</summary>
        private void CalcImportanceFull(int nseq, double[] eff, string[] seq, MafftLocalHom[][] localhom)
        {
            var importance = new double[Math.Max(_nlenmax, 1) + 2];
            var ieff = new double[nseq];
            double totaleff = 0.0;
            for (int i = 0; i < nseq; i++)
            {
                int nogaplen = seq[i].Count(c => c != '-');
                ieff[i] = nogaplen == 0 ? 0.0 : eff[i];
                totaleff += ieff[i];
            }

            for (int i = 0; i < nseq; i++)
                ieff[i] /= totaleff;
            for (int i = 0; i < nseq; i++)
            {
                Array.Clear(importance);
                for (int j = 0; j < nseq; j++)
                {
                    if (i == j)
                        continue;
                    for (var t = localhom[i][j]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1)
                            continue;
                        for (int pos = t.Start1; pos <= t.End1; pos++)
                            importance[pos] += ieff[j];
                    }
                }

                for (int j = 0; j < nseq; j++)
                {
                    if (i == j)
                        continue;
                    if (localhom[i][j].Opt == -1.0)
                        continue;
                    for (var t = localhom[i][j]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0)
                            continue;
                        double tmpdouble = 0.0;
                        int len = 0;
                        for (int pos = t.Start1; pos <= t.End1; pos++)
                        {
                            tmpdouble += importance[pos];
                            len++;
                        }

                        tmpdouble /= len;
                        t.Importance = tmpdouble * t.Opt;
                    }
                }
            }

            for (int i = 0; i < nseq - 1; i++)
            {
                for (int j = i + 1; j < nseq; j++)
                {
                    for (MafftLocalHom? t1 = localhom[i][j], t2 = localhom[j][i]; t1 != null && t2 != null; t1 = t1.Next, t2 = t2.Next)
                    {
                        if (t1.Opt == -1.0 || t2.Opt == -1.0)
                            continue;
                        double imp = 0.5 * (t1.Importance + t2.Importance);
                        t1.Importance = t2.Importance = imp;
                    }
                }
            }
        }

        // -------------------------------------------------------------------- the new columns of all new sequences

        /// <summary><c>makegaplistcompact</c>.</summary>
        private static void MakeGapListCompact(int len, int[] p, int[] c, int[] l)
        {
            int prep = -1;
            for (int i = 0; i < len + 2; i++)
            {
                int pg;
                if ((pg = p[i] - prep - 1) > 0 && l[i] > 0)
                    c[i] = pg < l[i] ? l[i] - pg : 0;
                else
                    c[i] = l[i];
                prep = p[i];
            }
        }

        /// <summary><c>gaplist2alnx</c>: the sequence <paramref name="s"/> with <c>l[k]</c> gap columns in front of its column <c>p[k]</c>.</summary>
        private static string GapList2Alnx(int len, string s, int[] l, int[] p, int lenlimit)
        {
            var a = new StringBuilder(s.Length + 16);
            int prevp = -1;
            int reslen = 0;
            for (int k = 0; k < len; k++)
            {
                int gaplen = l[k];
                int pi = p[k];
                if ((reslen += gaplen) > lenlimit)
                    throw new MafftProgramFailureException("Length over. Please recompile!\n");
                a.Append('-', gaplen);
                int pos = prevp + 1;
                int posl = pi - pos;
                if (posl != 0)
                {
                    if ((reslen += posl) > lenlimit)
                        throw new MafftProgramFailureException("Length over. Please recompile\n");
                    a.Append(s, pos, posl);
                }

                if (reslen++ > lenlimit)
                    throw new MafftProgramFailureException("Length over. Please recompile\n");
                a.Append(s[pos + posl]);
                prevp = pi;
            }

            int gaplenLast = l[len];
            int piLast = p[len];
            if ((reslen += gaplenLast) > lenlimit)
                throw new MafftProgramFailureException("Length over. Please recompile!\n");
            a.Append('-', gaplenLast);
            int posLast = prevp + 1;
            int poslLast = piLast - posLast;
            if (poslLast != 0)
            {
                if ((reslen += poslLast) > lenlimit)
                    throw new MafftProgramFailureException("Length over. Please recompile\n");
                a.Append(s, posLast, poslLast);
            }

            return a.ToString();
        }

        /// <summary><c>insertgapsbyotherfragments_compact</c>: the new sequence, in the coordinates of the old alignment, with '=' where the other new sequences made columns.</summary>
        private static string InsertGapsByOtherFragmentsCompact(int len, string s, int[] l, int[] p)
        {
            // the program reads the buffer of the sequence: behind its end are NUL characters (the buffer is zeroed), and the C string ends at the first of them
            char At(int k) => k < s.Length ? s[k] : '\0';
            var a = new StringBuilder(s.Length + 16);
            int si = 0;
            int prevp = -1;
            int i;
            for (i = 0; i < len; i++)
            {
                int gaplen = l[i];
                int pi = p[i];
                int pos = prevp + 1;
                while (gaplen-- > 0)
                {
                    pos++;
                    a.Append(At(si++));
                }

                while (pos++ < pi)
                    a.Append('=');
                a.Append(At(si++));
                prevp = pi;
            }

            int gaplenLast = l[i];
            int piLast = p[i];
            int posLast = prevp + 1;
            while (gaplenLast-- > 0)
            {
                posLast++;
                a.Append(At(si++));
            }

            while (posLast++ < piLast)
                a.Append('=');
            string res = a.ToString();
            int nul = res.IndexOf('\0');
            return nul >= 0 ? res[..nul] : res;
        }

        /// <summary><c>adjustposmap</c>: the positions of the columns of the old alignment after the new columns of <paramref name="gaplist"/>.</summary>
        private static void AdjustPosMap(int len, int[] posmap, int[] gaplist)
        {
            var newposmap = new int[len + 2];
            int zure = 0;
            for (int i = 0; i < len; i++)
            {
                zure += gaplist[i];
                newposmap[i] = posmap[i] + zure;
            }

            zure += gaplist[len];
            newposmap[len] = posmap[len] + zure;
            for (int i = 0; i < len + 1; i++)
                posmap[i] = newposmap[i];
        }

        internal struct Blocktorealign
        {
            internal int Start, End, NNewRes;
        }

        /// <summary><c>countnewres</c>: the blocks of the alignment in which new sequences have columns of their own.</summary>
        private static void CountNewRes(int len, Blocktorealign[] realign, int[] posmap, int[] gaplist)
        {
            int regstart = 0;
            int len1 = len + 1;
            for (int i = 0; i < len1; i++)
            {
                if (realign[i].NNewRes != 0 || gaplist[i] != 0)
                {
                    int regend = posmap[i] - 1;
                    realign[i].Start = regstart;
                    realign[i].End = regend;
                }

                if (gaplist[i] != 0)
                    realign[i].NNewRes++;
                regstart = posmap[i] + 1;
            }
        }

        private static int SearchLet(string s, int from, int to)
        {
            for (int p = from; p <= to; p++)
                if (s[p] != '-')
                    return p - from;
            return -1;
        }

        /// <summary><c>smoothing1</c> (<c>-p</c>): neighbouring blocks are joined.</summary>
        private void SmoothingOne(int njob, int nadd, int lenfull, string[] seq, Blocktorealign[] realign)
        {
            int norg = njob - nadd;
            Report("Smoothing1\n");
            var rows = seq.Select(x => x.ToCharArray()).ToArray();
            for (int i = 1; i < lenfull + 1; i++)
            {
                if (realign[i - 1].NNewRes != 0 && realign[i].NNewRes != 0)
                {
                    int postoshiftto = realign[i - 1].Start;
                    int postoshiftfrom = realign[i - 1].End + 1;
                    if (postoshiftfrom != realign[i].Start - 1)
                    {
                        Report("Unexpected pattern??? i=" + i + ", realign[i-1].end=" + realign[i - 1].End + ", realign[i].start=" + realign[i].Start + "\n");
                        throw new MafftProgramFailureException("");
                    }

                    realign[i].NNewRes += realign[i - 1].NNewRes;
                    realign[i].Start = realign[i - 1].Start + 1;
                    realign[i - 1].NNewRes = 0;
                    realign[i - 1].Start = realign[i - 1].End = 0;
                    for (int j = 0; j < norg; j++)
                    {
                        if (rows[j][postoshiftto] != '-')
                        {
                            Report("Unexpected pattern 2???\n");
                            throw new MafftProgramFailureException("");
                        }

                        rows[j][postoshiftto] = rows[j][postoshiftfrom];
                        rows[j][postoshiftfrom] = '-';
                    }

                    for (int j = norg; j < njob; j++)
                    {
                        if (rows[j][postoshiftto] == '-')
                        {
                            int shiftonadd = SearchLet(new string(rows[j]), postoshiftto, postoshiftfrom);
                            if (shiftonadd != -1)
                            {
                                rows[j][postoshiftto] = rows[j][postoshiftto + shiftonadd];
                                rows[j][postoshiftto + shiftonadd] = '-';
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < njob; i++)
                seq[i] = new string(rows[i]);
        }

        /// <summary><c>smoothing2</c> (<c>-p</c>): blocks that are two columns apart are joined.</summary>
        private void SmoothingTwo(int njob, int nadd, int lenfull, string[] seq, Blocktorealign[] realign)
        {
            int norg = njob - nadd;
            Report("Smoothing2\n");
            var rows = seq.Select(x => x.ToCharArray()).ToArray();
            for (int i = 2; i < lenfull + 1; i++)
            {
                if (realign[i - 2].NNewRes != 0 && realign[i - 1].NNewRes == 0 && realign[i].NNewRes != 0)
                {
                    int postoshiftto = realign[i - 2].Start;
                    int postoshiftfrom = realign[i - 2].End + 1;
                    if (postoshiftfrom != realign[i].Start - 2)
                    {
                        Report("Unexpected pattern??? i=" + i + ", realign[i-1].end=" + realign[i - 1].End + ", realign[i].start=" + realign[i].Start + "\n");
                        throw new MafftProgramFailureException("");
                    }

                    realign[i].NNewRes += realign[i - 2].NNewRes;
                    realign[i].Start = realign[i - 2].Start + 2;
                    realign[i - 2].NNewRes = 0;
                    realign[i - 2].Start = realign[i - 2].End = 0;
                    for (int j = 0; j < norg; j++)
                    {
                        if (rows[j][postoshiftto] != '-')
                        {
                            Report("Unexpected pattern 2. i=" + i + ", postoshiftto=" + postoshiftto + ", postoshiftfrom=" + postoshiftfrom + ", seq[" + j + "][" + postoshiftto + "]=" + rows[j][postoshiftto] + "\n");
                            throw new MafftProgramFailureException("");
                        }

                        rows[j][postoshiftto] = rows[j][postoshiftfrom];
                        rows[j][postoshiftfrom] = '-';
                        if (rows[j][postoshiftto + 1] != '-')
                        {
                            Report("Unexpected pattern 3???\n");
                            throw new MafftProgramFailureException("");
                        }

                        rows[j][postoshiftto + 1] = rows[j][postoshiftfrom + 1];
                        rows[j][postoshiftfrom + 1] = '-';
                    }

                    for (int j = norg; j < njob; j++)
                    {
                        if (rows[j][postoshiftto] == '-')
                        {
                            int shiftonadd = SearchLet(new string(rows[j]), postoshiftto, postoshiftfrom);
                            if (shiftonadd != -1)
                            {
                                rows[j][postoshiftto] = rows[j][postoshiftto + shiftonadd];
                                rows[j][postoshiftto + shiftonadd] = '-';
                            }
                        }

                        if (rows[j][postoshiftto + 1] == '-')
                        {
                            int shiftonadd = SearchLet(new string(rows[j]), postoshiftto + 1, postoshiftfrom + 1);
                            if (shiftonadd != -1)
                            {
                                rows[j][postoshiftto + 1] = rows[j][postoshiftto + 1 + shiftonadd];
                                rows[j][postoshiftto + 1 + shiftonadd] = '-';
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < njob; i++)
                seq[i] = new string(rows[i]);
        }

        // -------------------------------------------------------------------- the realignment of the new residues that share a place

        /// <summary><c>profilealignment2</c> (<c>addfunctions.c</c>): the profile of <paramref name="aln0"/> against the profile of <paramref name="aln2"/>, both rows replaced by the aligned rows.</summary>
        private void ProfileAlignment2(int n0, int n2, string[] aln0, string[] aln2, int alloclen)
        {
            if (aln0[0].Length <= 1 && aln2[0].Length <= 1)
                return;
            CommonGapPickN(n0, aln0);
            CommonGapPickN(n2, aln2);
            var effarr0 = new double[n0];
            var effarr2 = new double[n2];
            int alcount0 = 0;
            var allgap0 = new bool[n0];
            for (int i = 0; i < n0; i++)
            {
                if (aln0[i].All(c => c == '-'))
                    allgap0[i] = true;
                else
                    alcount0++;
            }

            int alcount2 = 0;
            var allgap2 = new bool[n2];
            for (int i = 0; i < n2; i++)
            {
                if (aln2[i].All(c => c == '-'))
                    allgap2[i] = true;
                else
                    alcount2++;
            }

            for (int i = 0; i < n0; i++)
                if (!allgap0[i])
                    effarr0[i] = 1.0 / alcount0;
            for (int i = 0; i < n2; i++)
                if (!allgap2[i])
                    effarr2[i] = 1.0 / alcount2;
            G.NewGapChar = '-';
            if (G.Alg == 'M')
                MsAlign(G.NDisConsweightMulti, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, null, null, null, null, 1, 1, null, -1, -1, 0, 1.0, 1.0);
            else
                AAlign(G.NDisConsweightMulti, G.Penalty, G.PenaltyEx, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, null, null, null, null, 1, 1, null, -1, -1, 0, 1.0, 1.0);
        }

        /// <summary><c>dorealignment_tree</c>: the residues that several new sequences put into the same place of the alignment are aligned again, the sequences that follow the same old sequence first, the groups then in the order of the guide tree of the old sequences; returns the number of columns by which the block became shorter.</summary>
        private int DoRealignmentTree(Blocktorealign block, string[] fullseq, ref int fullseqlen, int norg, int[][] topol0, int[][] topol1, int[] follows, int njob)
        {
            int localloclen = 4 * (block.End - block.Start + 1);
            var tmpaln0 = new List<string>();
            var iinf0 = new List<int>();
            var iinf1 = new List<int>();
            var nearest = new List<int>();
            for (int i = 0; i < njob; i++)
            {
                string tmp = GapPick0(fullseq[i].Substring(block.Start, block.End - block.Start + 1));
                if (tmp.Length != 0)
                {
                    if (i < norg)
                        throw new MafftProgramFailureException("BUG!!!!\n");
                    tmpaln0.Add(tmp);
                    iinf0.Add(i);
                    nearest.Add(follows[i - norg]);
                }
                else
                {
                    iinf1.Add(i);
                }
            }

            int n0 = tmpaln0.Count, n1 = iinf1.Count;
            var aln0 = tmpaln0.ToArray();
            var group = new int[n0];
            var g2n = new int[n0];
            var members = new List<List<int>>();
            int ngroup = 0;
            for (int i = 0; i < n0; i++)
            {
                int j;
                for (j = 0; j < i; j++)
                    if (nearest[j] == nearest[i])
                        break;
                if (j == i)
                {
                    group[i] = ngroup++;
                    members.Add(new List<int>());
                }
                else
                {
                    group[i] = group[j];
                }

                members[group[i]].Add(i);
                g2n[group[i]] = nearest[i];
            }

            for (int g = 0; g < ngroup; g++)
                members[g] = members[g].OrderByDescending(m => aln0[m].Length).ToList(); // qsort( lencomp ): the longer first (glibc's qsort is a stable merge sort)

            void Align(IReadOnlyList<int> rows0, IReadOnlyList<int> rows2)
            {
                var a = rows0.Select(m => aln0[m]).ToArray();
                var b = rows2.Select(m => aln0[m]).ToArray();
                ProfileAlignment2(a.Length, b.Length, a, b, localloclen);
                for (int k = 0; k < rows0.Count; k++)
                    aln0[rows0[k]] = a[k];
                for (int k = 0; k < rows2.Count; k++)
                    aln0[rows2[k]] = b[k];
            }

            for (int g = 0; g < ngroup; g++)
            {
                for (int j = 1; j < members[g].Count; j++)
                    Align(new[] { members[g][j] }, members[g].Take(j).ToList());
            }

            var topolpick0 = new List<int>[ngroup];
            var topolpick1 = new List<int>[ngroup];
            var pickhistory = new int[ngroup];
            for (int i = 0; i < ngroup; i++)
            {
                topolpick0[i] = new List<int>();
                topolpick1[i] = new List<int>();
                pickhistory[i] = -1;
            }

            int nhit = 0;
            for (int i = 0; i < norg - 1; i++)
            {
                int hit1 = FindGroup(topol0[i], g2n, ngroup);
                if (hit1 < 0)
                    continue;
                int hit2 = FindGroup(topol1[i], g2n, ngroup);
                if (hit2 < 0)
                    continue;
                topolpick0[nhit] = pickhistory[hit1] == -1 ? new List<int> { hit1 } : topolpick0[pickhistory[hit1]].Concat(topolpick1[pickhistory[hit1]]).ToList();
                topolpick1[nhit] = pickhistory[hit2] == -1 ? new List<int> { hit2 } : topolpick0[pickhistory[hit2]].Concat(topolpick1[pickhistory[hit2]]).ToList();
                pickhistory[hit1] = nhit;
                pickhistory[hit2] = nhit;
                nhit++;
            }

            for (int i = 0; i < ngroup - 1; i++)
            {
                var prof1 = topolpick0[i].SelectMany(g => members[g]).ToList();
                var prof2 = topolpick1[i].SelectMany(g => members[g]).ToList();
                Align(prof1, prof2);
            }

            int newlen = aln0[0].Length;
            int zure = block.End - block.Start + 1 - newlen;
            if (fullseqlen < fullseq[0].Length - (block.End - block.Start + 1) + newlen + 1)
            {
                fullseqlen = fullseq[0].Length * 2;
                Report("reallocating...");
                Report("done.\n");
            }

            int posinold = block.End + 1;
            for (int i = 0; i < n0; i++)
                fullseq[iinf0[i]] = string.Concat(fullseq[iinf0[i]].AsSpan(0, block.Start), aln0[i].AsSpan(0, newlen), fullseq[iinf0[i]].AsSpan(posinold));
            for (int i = 0; i < n1; i++)
                fullseq[iinf1[i]] = string.Concat(fullseq[iinf1[i]].AsSpan(0, block.Start), new string('-', newlen), fullseq[iinf1[i]].AsSpan(posinold));
            return zure;
        }

        private static int FindGroup(int[] members, int[] g2n, int ngroup)
        {
            foreach (int e in members)
                for (int k = 0; k < ngroup; k++)
                    if (e == g2n[k])
                        return k;
            return -1;
        }
    }

    #endregion
}
