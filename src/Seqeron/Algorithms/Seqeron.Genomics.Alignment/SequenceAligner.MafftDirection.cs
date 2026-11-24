#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325, S138, S1541
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --adjustdirection / --adjustdirectionaccurately (makedirectionlist.c, setdirection.c)

    internal sealed partial class MafftGlobals
    {
        /// <summary><c>thresholdtorev</c>, <c>dodp</c> (<c>-d</c>), <c>addfragment</c> (<c>-F</c>), <c>mode</c> (<c>-o</c>), <c>reflim</c> (<c>-r</c>), <c>contrastsort</c> (<c>-c</c> turns it off) of <c>makedirectionlist</c>.</summary>
        internal double ThresholdToRev;
        internal bool DoDp, AddFragment, ContrastSort = true;
        internal char DirectionMode = '2';
        internal int RefLim = 1000;

        /// <summary><c>arguments()</c> of <c>makedirectionlist.c</c> for the options the script <c>mafft</c> passes (<c>-F -C n -m -I nadd -i infile -t 0.00 -r 5000 -o a [-d]</c>).</summary>
        internal static MafftGlobals ParseMakedirectionlist(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = 2, PamN = 200, Alg = 'm' };
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : "";
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'i': g.Report("inputfile = " + Next() + "\n"); rest = true; break;
                        case 'I': g.Nadd = CAtoi(Next()); g.Report("nadd = " + g.Nadd.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'C': g.Nthread = CAtoi(Next()); g.Report("nthread = " + g.Nthread.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); g.Report("ppenalty_ex = " + g.PpenaltyEx.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); g.Report("kappa = " + g.KimuraR.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; g.Report("jtt/kimura " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 't': g.ThresholdToRev = CAtof(Next()); g.Report("thresholdtorev = " + g.ThresholdToRev.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'o': { string v = Next(); g.DirectionMode = v.Length > 0 ? v[0] : '\0'; g.Report("mode = " + g.DirectionMode + "\n"); rest = true; break; }
                        case 'r': g.RefLim = CAtoi(Next()); g.Report("reflim = " + g.RefLim.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'c': g.ContrastSort = false; break;
                        case 'd': g.DoDp = true; break;
                        case 'F': g.AddFragment = true; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'S': g.Alg = 'S'; break;
                        case 'M': g.Alg = 'M'; break;
                        case 'm': g.Alg = 'm'; break;
                        case 'G': g.Alg = 'G'; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        default:
                            throw new NotSupportedException($"The MAFFT option -{c} of makedirectionlist is not available in this port.");
                    }
                }
            }
            return g;
        }
    }

    internal sealed partial class MafftRun
    {
        /// <summary><c>creverse</c> of <c>io.c</c>: the complement of a letter (IUPAC codes, both cases); other characters stay.</summary>
        private static char CReverse(char f) => f switch
        {
            'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'U' => 'A', 'M' => 'K', 'R' => 'Y', 'W' => 'W', 'S' => 'S', 'Y' => 'R', 'K' => 'M', 'V' => 'B', 'H' => 'D', 'D' => 'H', 'B' => 'V', 'N' => 'N',
            'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a', 'u' => 'a', 'm' => 'k', 'r' => 'y', 'w' => 'w', 's' => 's', 'y' => 'r', 'k' => 'm', 'v' => 'b', 'h' => 'd', 'd' => 'h', 'b' => 'v', 'n' => 'n',
            _ => f,
        };

        /// <summary><c>sreverse</c> of <c>io.c</c>: the reverse complement; a sequence that has more U than T gets U for the T it produces.</summary>
        internal static string SReverse(string s)
        {
            int numt = s.Count(c => c is 't' or 'T');
            int numu = s.Count(c => c is 'u' or 'U');
            var r = new char[s.Length];
            for (int i = 0; i < s.Length; i++)
                r[s.Length - 1 - i] = CReverse(s[i]);
            if (numu > numt)
            {
                for (int i = 0; i < r.Length; i++)
                {
                    if (r[i] == 't')
                        r[i] = 'u';
                    else if (r[i] == 'T')
                        r[i] = 'U';
                }
            }

            return new string(r);
        }

        /// <summary>What <c>makedirectionlist</c> wrote: the lines of the file <c>_direction</c> (standard output) and its exit status.</summary>
        internal sealed record DirectionList(string[] Lines, int ExitCode);

        /// <summary><c>_F_</c> / <c>_R_</c> and the first ten characters of the title.</summary>
        private static string DirectionLabel(char dir, string name)
        {
            string tail = name.Length > 1 ? name[1..] : "";
            return "_" + dir + "_" + (tail.Length > 10 ? tail[..10] : tail);
        }

        private static int ShortCommon(short[] table, int[] pointt, int maxl, int tsize)
        {
            // localcommonsextet_p2: the number of 6-mers of the second sequence that the first one has (a 6-mer is counted as often as the first one holds it)
            if (pointt.Length == 0)
                return 0;
            var memo = new short[tsize];
            int value = 0;
            foreach (int point in pointt)
            {
                short tmp = memo[point]++;
                if (tmp < table[point])
                    value++;
            }

            return value;
        }

        private static int[] KmerPoints(string s, int[] grpBuf, int[] mult)
        {
            // seq_grp_nuc + makepointtable_nuc (the 6-mers of the groups a, c, g, t; other letters are skipped); fewer than six groups: no table
            var grp = new List<int>(s.Length);
            foreach (char c in s)
            {
                int tmp = grpBuf[c];
                if (tmp < 4)
                    grp.Add(tmp);
            }

            return grp.Count < 6 ? Array.Empty<int>() : MakePointTable(grp.ToArray(), grp.Count, 6, mult, 4);
        }

        /// <summary>glibc's <c>qsort</c> (the merge sort of <c>msort.c</c>): stable, with the comparison of the program (which may not be a total order).</summary>
        private static void MergeSortLikeGlibc<T>(T[] a, Comparison<T> cmp)
        {
            var tmp = new T[a.Length];
            Sort(0, a.Length);

            void Sort(int lo, int n)
            {
                if (n <= 1)
                    return;
                int n1 = n / 2, n2 = n - n1;
                Sort(lo, n1);
                Sort(lo + n1, n2);
                int i1 = lo, i2 = lo + n1, k = 0;
                int e1 = lo + n1, e2 = lo + n;
                while (i1 < e1 && i2 < e2)
                    tmp[k++] = cmp(a[i1], a[i2]) <= 0 ? a[i1++] : a[i2++];
                while (i1 < e1)
                    tmp[k++] = a[i1++];
                Array.Copy(tmp, 0, a, lo, k);
            }
        }

        /// <summary><c>makedirectionlist</c> (<c>-m</c>, the 6-mer method; <c>-d</c>: the local alignment score instead): for every sequence of the input (every new one with <c>-I n</c>) whether it is closer to the others as it is (_F_) or reverse complemented (_R_).</summary>
        internal DirectionList MakeDirectionList(IReadOnlyList<string> rows, IReadOnlyList<string> names)
        {
            int njob = rows.Count;
            int nadd = G.Nadd;
            int nlenmax = rows.Max(r => r.Length);
            _njob = njob;
            _nlenmax = nlenmax;
            char dorp = MafftEngine.GuessNucleotide(rows, new MafftOptions()) ? 'd' : 'p';
            var name = names.ToArray();
            if (dorp != 'd')
            {
                Report("Not necessary!\n");
                return new DirectionList(Enumerable.Range(0, njob).Select(i => "_F_" + (name[i].Length > 1 ? name[i][1..].PadRight(10)[..10] : "          ")).ToArray(), 1);
            }

            G.Dorp = 'd';
            var seq = rows.Select(r => r.ToLowerInvariant()).ToArray();
            G.Constants(seq);
            if (CheckAlphabet(seq, "Illegal") is { } unknown)
                throw unknown;
            Report("\n");
            if (G.Alg != 'm')
                throw new NotSupportedException("makedirectionlist -" + G.Alg + ": the script passes -m only.");
            bool dodp = G.DoDp;
            int tsize = 4096;
            var mult = new[] { 1024, 256, 64, 16, 4, 1 };
            int[]? grpBuf = null;
            if (!dodp)
            {
                grpBuf = new int[256];
                Array.Fill(grpBuf, 99);
                for (int c = 0; c < 256; c++)
                    grpBuf[c] = G.AminoGrp[c];
            }

            int maxl = 0;
            var tmpseq = "";
            for (int i = 0; i < njob; i++)
                maxl = Math.Max(maxl, GapPick0(seq[i]).Length);
            Report("Step 1/2\n");
            var pointt = new int[njob][];
            var pointtRev = new int[njob][];
            var spointt = new int[njob][];
            if (!dodp)
            {
                int iend0 = nadd > 0 ? njob - nadd : 0;
                for (int i = 0; i < iend0; i++)
                {
                    tmpseq = GapPick0(seq[i]);
                    seq[i] = tmpseq;
                    pointt[i] = KmerPoints(tmpseq, grpBuf!, mult);
                    spointt[i] = pointt[i];
                }

                int istart0 = nadd > 0 ? njob - nadd : 0;
                for (int i = istart0; i < njob; i++)
                {
                    tmpseq = GapPick0(seq[i]);
                    seq[i] = tmpseq;
                    string rev = SReverse(tmpseq);
                    pointt[i] = KmerPoints(tmpseq, grpBuf!, mult);
                    pointtRev[i] = KmerPoints(rev, grpBuf!, mult);
                    spointt[i] = pointt[i];
                }
            }

            var contrastorder = new int[njob];
            if (G.ContrastSort)
            {
                int istart, iend;
                if (nadd > 0)
                {
                    for (int i = 0; i < njob - nadd; i++)
                        contrastorder[i] = i;
                    istart = njob - nadd;
                    iend = nadd;
                }
                else
                {
                    istart = 0;
                    iend = njob;
                }

                var res = new double[iend];
                if (dodp)
                {
                    // selfdpthread: the score of the sequence with itself minus the score with its reverse complement
                    int j = -1;
                    while (true)
                    {
                        j++;
                        if (j % 100 == 0)
                            Report(G.Nthread > 0 ? j.ToString(CultureInfo.InvariantCulture) + " / " + iend.ToString(CultureInfo.InvariantCulture) + " (thread 0)   \r" : j.ToString(CultureInfo.InvariantCulture) + " / " + iend.ToString(CultureInfo.InvariantCulture) + "      \r");
                        if (j == iend)
                            break;
                        string s = seq[istart + j];
                        string rev = SReverse(s);
                        res[j] = L__align11(G.NDisConsweightMulti, 0.0, new[] { s }, new[] { s }, nlenmax * 2, false, out _, out _);
                        res[j] -= L__align11(G.NDisConsweightMulti, 0.0, new[] { s }, new[] { rev }, nlenmax * 2, false, out _, out _);
                    }
                }
                else
                {
                    for (int i = 0; i < iend; i++)
                    {
                        if (i % 100 == 1)
                            Report(i.ToString(CultureInfo.InvariantCulture) + "   \r");
                        int ii = istart + i;
                        var table1 = new short[tsize];
                        foreach (int point in pointt[ii])
                            table1[point]++;
                        res[i] = ShortCommon(table1, pointt[ii], maxl, tsize);
                        var table1Rev = new short[tsize];
                        foreach (int point in pointtRev[ii])
                            table1Rev[point]++;
                        res[i] -= ShortCommon(table1Rev, pointt[ii], maxl, tsize);
                    }
                }

                var arr = Enumerable.Range(0, iend).Select(i => (Pos: i, Dif: res[i])).ToArray();
                MergeSortLikeGlibc(arr, (x, y) => (int)(y.Dif - x.Dif));
                for (int i = 0; i < iend; i++)
                    contrastorder[istart + i] = arr[i].Pos + istart;
            }
            else
            {
                for (int i = 0; i < njob; i++)
                    contrastorder[i] = i;
            }

            int iendRef = nadd > 0 ? njob - nadd : 1;
            for (int i = 0; i < iendRef; i++)
            {
                int ic = contrastorder[i];
                seq[ic] = GapPick0(seq[ic]);
                name[ic] = DirectionLabel('F', name[ic]);
            }

            Report("\n\nStep 2/2\n");
            int istart2 = nadd > 0 ? njob - nadd : 1;
            var resf = new int[njob];
            var resr = new int[njob];
            var map = new int[njob];
            for (int i = istart2; i < njob; i++)
            {
                int ic = contrastorder[i];
                string fwd = GapPick0(seq[ic]);
                seq[ic] = fwd;
                string revseq = SReverse(fwd);
                short[]? table1 = null, table1Rev = null;
                if (!dodp)
                {
                    table1 = new short[tsize];
                    table1Rev = new short[tsize];
                    foreach (int point in pointt[ic])
                        table1[point]++;
                    foreach (int point in pointtRev[ic])
                        table1Rev[point]++;
                }

                int iend;
                if (nadd > 0 && G.AddFragment)
                    iend = njob - nadd;
                else
                    iend = i;
                for (int j = 0; j < iend; j++)
                    map[j] = contrastorder[j];
                if (iend > G.RefLim)
                    iend = G.RefLim;
                if (i % 100 == 1)
                    Report(G.Nthread > 0 ? " " + i.ToString(CultureInfo.InvariantCulture) + " / " + njob.ToString(CultureInfo.InvariantCulture) + " (" + G.Nthread.ToString(CultureInfo.InvariantCulture) + " threads)   \r" : " " + i.ToString(CultureInfo.InvariantCulture) + " / " + njob.ToString(CultureInfo.InvariantCulture) + "   \r");
                for (int j = 0; j < iend; j++)
                {
                    resf[j] = dodp ? (int)L__align11(G.NDisConsweightMulti, 0.0, new[] { fwd }, new[] { seq[map[j]] }, nlenmax * 2, false, out _, out _) : ShortCommon(table1!, spointt[map[j]], maxl, tsize);
                }

                for (int j = 0; j < iend; j++)
                {
                    resr[j] = dodp ? (int)L__align11(G.NDisConsweightMulti, 0.0, new[] { revseq }, new[] { seq[map[j]] }, nlenmax * 2, false, out _, out _) : ShortCommon(table1Rev!, spointt[map[j]], maxl, tsize);
                }

                double resForward, resReverse;
                char mode = G.DirectionMode;
                if (mode == '2')
                {
                    resForward = BestTwoMean(resf, iend);
                    resReverse = BestTwoMean(resr, iend);
                }
                else if (mode == '1')
                {
                    resReverse = 0.0;
                    for (int j = 0; j < iend; j++)
                        if (resReverse < resr[j])
                            resReverse = resr[j];
                    resForward = 0.0;
                    for (int j = 0; j < iend; j++)
                        if (resForward < resf[j])
                            resForward = resf[j];
                }
                else if (mode == 'd')
                {
                    resReverse = 0.0;
                    for (int j = 0; j < iend; j++)
                        if (resReverse < resr[j] - resf[j])
                            resReverse = resr[j] - resf[j];
                    resForward = 0.0;
                    for (int j = 0; j < iend; j++)
                        if (resForward < resf[j] - resr[j])
                            resForward = resf[j] - resr[j];
                }
                else if (mode == 'a')
                {
                    resReverse = 0.0;
                    for (int j = 0; j < iend; j++)
                        resReverse += resr[j];
                    resReverse /= iend;
                    resForward = 0.0;
                    for (int j = 0; j < iend; j++)
                        resForward += resf[j];
                    resForward /= iend;
                }
                else
                {
                    Report("Unknown mode!\n");
                    throw new MafftProgramFailureException("");
                }

                if (resReverse > resForward)
                {
                    seq[ic] = revseq;
                    name[ic] = DirectionLabel('R', name[ic]);
                    if (!dodp)
                        spointt[ic] = pointtRev[ic];
                }
                else
                {
                    name[ic] = DirectionLabel('F', name[ic]);
                    if (!dodp)
                        spointt[ic] = pointt[ic];
                }
            }

            if (name[0].Length > 1 && name[0][1] == 'R')
            {
                for (int j = 0; j < njob; j++)
                {
                    var chars = name[j].ToCharArray();
                    chars[1] = chars[1] == 'R' ? 'F' : 'R';
                    name[j] = new string(chars);
                }
            }

            Report("\n");
            Report("makedirectionlist (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n" + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new DirectionList(name, 0);
        }

        private static double BestTwoMean(int[] r, int iend)
        {
            int mres = 0, mres2 = 0;
            for (int j = 0; j < iend; j++)
            {
                int ires = r[j];
                if (ires > mres2)
                {
                    if (ires > mres)
                    {
                        mres2 = mres;
                        mres = ires;
                    }
                    else
                    {
                        mres2 = ires;
                    }
                }
            }

            return (double)(mres + mres2) / 2;
        }

        /// <summary>
        /// <c>setdirection -d _direction -i file [-H n] [-r]</c>: the sequences that the file <c>_direction</c> marks _R_ are reverse complemented and, unless <paramref name="originalFile"/> (<c>-r</c>: the file <c>orig</c> of <c>--anysymbol</c>), get the prefix _R_
        /// in their title; with <c>-H</c> (<c>--merge</c>) the sequences of one subalignment must have the same direction. The messages go to <see cref="MafftGlobals.Log"/>.
        /// </summary>
        internal (string[] Rows, string[] Names) SetDirectionProgram(IReadOnlyList<string> rows, IReadOnlyList<string> names, IReadOnlyList<string> directionLines, bool originalFile, bool subalignment)
        {
            int njob = rows.Count;
            Report("directionfile = _direction\n");
            Report("inputfile = " + (originalFile ? "orig" : "infile") + "\n");
            Report("subalignment = " + (subalignment ? 1 : 0).ToString(CultureInfo.InvariantCulture) + "\n");
            Report("subalignmentoffset = " + (subalignment ? G.SubalignmentOffset : 0).ToString(CultureInfo.InvariantCulture) + "\n");
            int nlenmax = rows.Max(r => r.Length);
            int nlenmin = rows.Min(r => r.Length);
            char dorp = MafftEngine.GuessNucleotide(rows, new MafftOptions()) ? 'd' : 'p';
            Report(njob.ToString(CultureInfo.InvariantCulture) + " x " + nlenmax.ToString(CultureInfo.InvariantCulture) + " - " + nlenmin.ToString(CultureInfo.InvariantCulture) + " " + dorp + "\n");
            var seq = rows.ToArray();
            var name = names.ToArray();
            var directions = new char[njob];
            for (int i = 0; i < njob; i++)
            {
                string line = i < directionLines.Count ? directionLines[i] : "";
                if (line.Length < 1 || line[0] != '_')
                    throw new MafftProgramFailureException("Format error!\n");
                if (line.Length > 1 && line[1] == 'R')
                {
                    directions[i] = 'R';
                    seq[i] = SReverse(seq[i]);
                    string t = name[i].Length > 1 ? name[i][1..] : "";
                    if (t.Length > MafftProgramIo.NameLimit - 1)
                        t = t[..(MafftProgramIo.NameLimit - 1)];
                    name[i] = originalFile ? "=" + t : "=_R_" + t;
                }
                else if (line.Length > 1 && line[1] == 'F')
                {
                    directions[i] = 'F';
                }
                else
                {
                    throw new MafftProgramFailureException("Format error!\n");
                }
            }

            if (subalignment)
            {
                var sub = ReadSubalignmentsTable(njob);
                for (int j = 0; j < sub.Groups.Count; j++)
                {
                    Report("Checking directions of sequences in subalignment" + j.ToString(CultureInfo.InvariantCulture) + "\n");
                    char firstdir = directions[sub.Groups[j][0]];
                    Report("firstdir = " + firstdir + "\n");
                    foreach (int member in sub.Groups[j])
                    {
                        if (directions[member] != firstdir)
                            throw new MafftProgramFailureException("\n\n#############################################################################\n\nDirection of nucleotide sequences seems to be inconsistent.\nPlease check the following two sequences:\n\tSequece no."
                                + (sub.Groups[j][0] + 1).ToString(CultureInfo.InvariantCulture) + " (" + name[sub.Groups[j][0]] + ")\n\tSequece no." + (member + 1).ToString(CultureInfo.InvariantCulture) + " (" + name[member] + ")\n\nThese sequences are in sub alignment no."
                                + (j + 1).ToString(CultureInfo.InvariantCulture) + " in your setting of --merge,\nbut their directions seem to be different.\n\n#############################################################################\n\n\n\n");
                    }

                    Report("OK!\n");
                }
            }

            return (seq, name);
        }
    }

    #endregion
}
