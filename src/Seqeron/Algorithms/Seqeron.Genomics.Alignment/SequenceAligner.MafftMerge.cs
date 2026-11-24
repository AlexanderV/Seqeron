#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --merge (readsubalignmentstable, the constrained UPGMA tree and the subalignment branch of disttbfast / tbfast / dvtditr)

    internal sealed partial class MafftGlobals
    {
        [ThreadStatic]
        private static string? _ambientSubalignmentsTable;

        /// <summary>The text of the file <c>_subalignmentstable</c> of the run in progress on this thread (<c>mafft --merge table</c>; the script copies the file there, the programs read it).</summary>
        internal string? SubalignmentsTableText = _ambientSubalignmentsTable;

        /// <summary><c>-H n</c> (<c>subalignment = 1</c>, <c>subalignmentoffset = n</c>).</summary>
        internal bool Subalignment;
        internal int SubalignmentOffset;

        /// <summary>Makes <paramref name="text"/> the file <c>_subalignmentstable</c> for every program that is started on this thread until the scope is disposed.</summary>
        internal static IDisposable UseSubalignmentsTable(string? text)
        {
            var previous = _ambientSubalignmentsTable;
            _ambientSubalignmentsTable = text;
            return new Scope(() => _ambientSubalignmentsTable = previous);
        }
    }

    /// <summary>The subalignments of <c>--merge</c>: the sequences (0-based) of every group, whether the group keeps the columns that only gaps have (<c>preservegaps</c>, a negative number in the table), and the size of the largest group.</summary>
    internal sealed record MafftSubalignments(List<int[]> Groups, List<bool> PreserveGaps, int MaxMem);

    internal sealed partial class MafftRun
    {
        /// <summary><c>atoi</c> of the C library applied to the text that starts at <paramref name="start"/> (leading blanks, a sign, digits; 0 when there are none).</summary>
        private static int AtoiAt(string s, int start) => MafftGlobals.CAtoi(s[start..]);

        /// <summary><c>countspace</c> (<c>io.c</c>): the number of numbers on a line of the table; a token that does not start with a number stops the program.</summary>
        private static int CountSpace(string line)
        {
            int v = 0;
            char status = 's';
            for (int p = 0; p < line.Length; p++)
            {
                char c = line[p];
                if (c == ' ' || c == '\t')
                {
                    status = 's';
                }
                else if (status == 's')
                {
                    if (c == '\n' || c == '#')
                        break;
                    v++;
                    status = 'n';
                    if (AtoiAt(line, p) == 0)
                        throw new MafftProgramFailureException("Format error? Sequences should be specified as 1, 2, 3...\n");
                }
            }
            return v;
        }

        /// <summary><c>readasubalignment</c> (<c>io.c</c>): the numbers of one line, 0-based; a negative number turns <c>preservegaps</c> on.</summary>
        private static List<int> ReadASubalignment(string line, out bool preserveGaps)
        {
            var t = new List<int>();
            preserveGaps = false;
            char status = 's';
            for (int p = 0; p < line.Length; p++)
            {
                char c = line[p];
                if (c == ' ' || c == '\t')
                {
                    status = 's';
                }
                else if (status == 's')
                {
                    if (c == '\n' || c == '#')
                        break;
                    status = 'n';
                    int v = AtoiAt(line, p);
                    if (v == 0)
                        throw new MafftProgramFailureException("Format error? Sequences must be specified as 1, 2, 3...\n");
                    if (v < 0)
                        preserveGaps = true;
                    t.Add(Math.Abs(v) - 1);
                }
            }
            return t;
        }

        /// <summary>
        /// <c>readsubalignmentstable</c> (<c>io.c</c>), both passes: the lines of the file <c>_subalignmentstable</c> that start with a number (comment lines and lines that <c>atoi</c> reads as 0 are skipped); stops with the program's messages
        /// for a number that is not a sequence number, a sequence in two groups and a number above the number of sequences. The program reads <c>tab01[p]</c> before it tests <c>p</c>: a number that is only a little above
        /// the number of sequences lands in the unused end of the allocated block (zero), a larger one writes into the next heap chunk (undefined: refused).
        /// </summary>
        private MafftSubalignments ReadSubalignmentsTable(int nseq)
        {
            string? text = G.SubalignmentsTableText;
            if (text == null)
                throw new MafftProgramFailureException("Cannot open _subalignmentstable\n");
            var lines = new List<string>();
            foreach (string raw in SplitFgetsLines(text, 999999))
            {
                if (raw.Length == 0 || raw[^1] != '\n')
                    throw new MafftProgramFailureException("too long line? \n");
                lines.Add(raw);
            }
            int nsub = 0, maxmem = 0;
            foreach (string line in lines)
            {
                if (line[0] == '#' || AtoiAt(line, 0) == 0)
                    continue;
                int nmem = CountSpace(line);
                maxmem = Math.Max(maxmem, nmem);
                nsub++;
            }
            Report("nsubalignments = " + nsub.ToString(CultureInfo.InvariantCulture) + "\nmaxmem = " + maxmem.ToString(CultureInfo.InvariantCulture) + "\n");
            // calloc( nseq, sizeof( int ) ): the block has room for the slots that glibc's chunk rounding leaves
            long usable = ((long)nseq * 4 + 8 + 15) / 16 * 16;
            if (usable < 32)
                usable = 32;
            int slots = (int)((usable - 8) / 4);
            var tab01 = new bool[slots];
            var groups = new List<int[]>();
            var preserve = new List<bool>();
            foreach (string line in lines)
            {
                if (line[0] == '#' || AtoiAt(line, 0) == 0)
                    continue;
                var members = ReadASubalignment(line, out bool pg);
                foreach (int p in members)
                {
                    if (p >= slots)
                        throw new NotSupportedException("--merge: the table names a sequence far above the number of sequences; the program writes behind its block 'tab01' into the next heap chunk (undefined behaviour).");
                    if (tab01[p])
                        throw new MafftProgramFailureException("\nSequence " + (p + 1).ToString(CultureInfo.InvariantCulture) + " appears in different groups.\nHierarchical grouping is not supported.\n\n");
                    tab01[p] = true;
                    if (p > nseq - 1)
                        throw new MafftProgramFailureException("Sequence " + (p + 1).ToString(CultureInfo.InvariantCulture) + " does not exist in the input sequence file.\n");
                }
                groups.Add(members.ToArray());
                preserve.Add(pg);
            }
            return new MafftSubalignments(groups, preserve, maxmem);
        }

        /// <summary><c>fgets( line, linelen - 1, fp )</c> over the whole text: each piece holds at most <paramref name="linelen"/> - 2 characters, a piece that ends without a newline is a piece that was cut.</summary>
        private static IEnumerable<string> SplitFgetsLines(string text, int linelen)
        {
            int pos = 0;
            while (pos < text.Length)
            {
                int nl = text.IndexOf('\n', pos);
                int end = nl < 0 ? text.Length : nl + 1;
                int take = Math.Min(end - pos, linelen - 2);
                yield return text.Substring(pos, take);
                pos += take;
            }
        }

        /// <summary><c>increaseintergroupdistanceshalfmtx</c> (<c>mltaln9.c</c>): the distances between a member of a group and a sequence outside it are raised to the largest distance inside the group, so that the group is joined first.</summary>
        private void IncreaseInterGroupDistances(double[][] eff, MafftSubalignments subs, int nseq)
        {
            int nwarned = 0;
            var tft = new bool[nseq];
            var others = new List<int>();
            Report("\n");
            for (int i = 0; i < subs.Groups.Count; i++)
            {
                var g = subs.Groups[i];
                if (g.Length < 2)
                    continue;
                Array.Clear(tft);
                foreach (int s in g)
                    tft[s] = true;
                others.Clear();
                for (int m = 0; m < nseq; m++)
                    if (!tft[m])
                        others.Add(m);
                double maxdist = 0.0;
                for (int m = 1; m < g.Length; m++)
                {
                    for (int k = 0; k < m; k++)
                    {
                        int s2 = g[m], s1 = g[k];
                        int sl = s2 > s1 ? s2 : s1, ss = s2 > s1 ? s1 : s2;
                        double dtmp = eff[ss][sl - ss];
                        if (dtmp > maxdist)
                            maxdist = dtmp;
                    }
                }
                foreach (int s2 in g)
                {
                    foreach (int s1 in others)
                    {
                        int sl = s2 > s1 ? s2 : s1, ss = s2 > s1 ? s1 : s2;
                        if (eff[ss][sl - ss] < maxdist)
                        {
                            if (eff[ss][sl - ss] < 0.5 && nwarned++ < 100)
                                Report("# Sequences " + (s2 + 1).ToString(CultureInfo.InvariantCulture) + " and " + (s1 + 1).ToString(CultureInfo.InvariantCulture) + " seem to be closely related, but are not in the same sub MSA (" + (i + 1).ToString(CultureInfo.InvariantCulture) + ") in your setting.\n");
                            eff[ss][sl - ss] = maxdist;
                        }
                    }
                }
            }
            if (nwarned > 100)
                Report("# Sequenc.... (more pairs)\n");
        }

        /// <summary>
        /// <c>fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained</c> (disttbfast, tbfast) and <c>fixed_supg_double_treeout_constrained</c> (dvtditr): the UPGMA tree of the half matrix <paramref name="eff"/> in which every subalignment is
        /// a monophyletic group: a join that would break a group is marked inconsistent and the closest consistent pair is taken (<c>inconsistent</c>, <c>inconsistentpairlist</c>). The groups keep their full membership in
        /// <paramref name="topol0"/> / <paramref name="topol1"/>; the file <c>infile.tree</c> is always written.
        /// </summary>
        private void ConstrainedSupg(int nseq, double[][] eff, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, MafftSubalignments subs, IReadOnlyList<string> treeNames)
        {
            IncreaseInterGroupDistances(eff, subs, nseq);
            double sueff1 = 1 - G.SueffGlobal;
            double sueff05 = G.SueffGlobal * 0.5;
            Func<double, double, double> cluster = G.TreeMethod switch
            {
                'X' => (d1, d2) => Math.Min(d1, d2) * sueff1 + (d1 + d2) * sueff05,
                'E' => (d1, d2) => (d1 + d2) * 0.5,
                'q' => (d1, d2) => Math.Min(d1, d2),
                _ => throw new MafftProgramFailureException("Unknown treemethod, " + G.TreeMethod + "\n"),
            };
            var tree = new string[nseq];
            for (int i = 0; i < nseq; i++)
                tree[i] = "\n" + (i + 1).ToString(CultureInfo.InvariantCulture) + "_" + TreeLabel(treeNames[i]) + "\n";
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
            var inconsistent = new bool[nseq, nseq];
            var inconsistentPairs = new List<(int, int)>();
            int maxInconsistentPairs = 1;
            var warned = new bool[subs.Groups.Count];
            for (int i = 0; i < nseq; i++)
                SetNearest(eff, next, mindisfrom, nearest, i);
            for (int i = 0; i < nseq; i++)
            {
                tmptmplen[i] = 0.0;
                hist[i] = -1;
                nmemar[i] = 1;
            }
            Report("\n");
            int im = -1, jm = -1;
            double minscore = 999.9;
            string treetmp = "";
            for (int k = 0; k < nseq - 1; k++)
            {
                if (k % 10 == 0)
                    Report("\r" + MafftGlobals.D5(k) + " / " + nseq);
                foreach (var (a, b) in inconsistentPairs)
                    inconsistent[a, b] = false;
                inconsistentPairs.Clear();
                bool firsttime = true;
                while (true)
                {
                    if (firsttime)
                    {
                        firsttime = false;
                        minscore = 999.9;
                        for (int a = 0; next[a] != -1; a = next[a])
                        {
                            if (mindisfrom[a] < minscore)
                            {
                                im = a;
                                minscore = mindisfrom[a];
                            }
                        }
                        jm = nearest[im];
                        if (jm < im)
                            (im, jm) = (jm, im);
                    }
                    else
                    {
                        minscore = 999.9;
                        for (int a = 0; next[a] != -1; a = next[a])
                        {
                            int i = a;
                            for (int b = next[a]; b != -1; b = next[b])
                            {
                                int j = b;
                                double d;
                                if (!inconsistent[i, j] && (d = eff[i][j - i]) < minscore)
                                {
                                    minscore = d;
                                    im = i;
                                    jm = j;
                                }
                            }
                            for (int b = 0; b != -1 && b != i; b = next[b])
                            {
                                int j = b;
                                double d;
                                if (!inconsistent[j, i] && (d = eff[j][i - j]) < minscore)
                                {
                                    minscore = d;
                                    im = j;
                                    jm = i;
                                }
                            }
                        }
                    }
                    bool allinconsistent = true;
                    for (int a = 0; next[a] != -1 && allinconsistent; a = next[a])
                        for (int b = next[a]; b != -1; b = next[b])
                            if (!inconsistent[a, b])
                            {
                                allinconsistent = false;
                                break;
                            }
                    if (allinconsistent)
                        throw new MafftProgramFailureException("\n\n\nPlease check whether the grouping is possible.\n\n\n");
                    var testtopol = new List<int>();
                    int prevnode = hist[im];
                    if (prevnode == -1)
                    {
                        testtopol.Add(im);
                    }
                    else
                    {
                        testtopol.AddRange(topol0[prevnode]);
                        testtopol.AddRange(topol1[prevnode]);
                    }
                    prevnode = hist[jm];
                    if (prevnode == -1)
                    {
                        testtopol.Add(jm);
                    }
                    else
                    {
                        testtopol.AddRange(topol0[prevnode]);
                        testtopol.AddRange(topol1[prevnode]);
                    }
                    int gi;
                    for (gi = 0; gi < subs.Groups.Count; gi++)
                    {
                        var grp = subs.Groups[gi];
                        if (grp.Any(testtopol.Contains))
                        {
                            if (!testtopol.All(grp.Contains) && !grp.All(testtopol.Contains))
                            {
                                if (!warned[gi])
                                {
                                    warned[gi] = true;
                                    Report("\n###################################################################\n# WARNING: Group " + (gi + 1).ToString(CultureInfo.InvariantCulture) + " is forced to be a monophyletic cluster.\n###################################################################\n");
                                }
                                inconsistent[im, jm] = true;
                                if (maxInconsistentPairs < inconsistentPairs.Count + 1)
                                {
                                    maxInconsistentPairs = inconsistentPairs.Count + 1;
                                    Report("Reallocated inconsistentpairlist, size=" + maxInconsistentPairs.ToString(CultureInfo.InvariantCulture) + "\n");
                                }
                                inconsistentPairs.Add((im, jm));
                                break;
                            }
                        }
                    }
                    if (gi == subs.Groups.Count)
                        break;
                }
                int pn = hist[im];
                dep[k].Child0 = pn;
                int nmemim = nmemar[im];
                topol0[k] = GroupMembers(pn, im, topol0, topol1, true);
                pn = hist[jm];
                dep[k].Child1 = pn;
                int nmemjm = nmemar[jm];
                topol1[k] = GroupMembers(pn, jm, topol0, topol1, true);
                minscore *= 0.5;
                len[k][0] = minscore - tmptmplen[im];
                len[k][1] = minscore - tmptmplen[jm];
                if (len[k][0] < 0.0)
                    len[k][0] = 0.0;
                if (len[k][1] < 0.0)
                    len[k][1] = 0.0;
                dep[k].DistFromTip = minscore;
                tmptmplen[im] = minscore;
                hist[im] = k;
                nmemar[im] = nmemim + nmemjm;
                mindisfrom[im] = 999.9;
                eff[im][jm - im] = 999.9;
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
                treetmp = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                tree[im] = treetmp;
                tree[jm] = null!;
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
                        int miniim = i < im ? i : im, maxiim = i < im ? im : i;
                        if (eff[miniim][maxiim - miniim] > mindisfrom[i])
                            SetNearest(eff, next, mindisfrom, nearest, i);
                    }
                }
            }
            TreeFile = treetmp + ";\n";
        }

        /// <summary>The file <c>GuideTree</c> that <c>cp infile.tree GuideTree</c> made (null when it failed or was not run).</summary>
        internal string? GuideTreeCopy;

        private static bool IncludeMemberOf(int[] mem, int[] cand) => mem.All(cand.Contains);

        private static bool SameMemberAs(int[] mem, int[] cand) => mem.Length == cand.Length && mem.All(cand.Contains);

        /// <summary>
        /// The <c>else if( subalignment )</c> branch of <c>main</c> of <c>disttbfast.c</c> / <c>tbfast.c</c> (<c>--merge</c>): every step of the tree is aligned ('a') except the steps inside a subalignment ('n'); a subalignment must consist of
        /// equally long rows and be a group of the tree. The rows of the subalignments keep their columns (minus the columns of gaps only unless <c>preservegaps</c>), the other rows lose their gaps.
        /// Returns null when disttbfast stops after its message (<c>goto chudan</c>, exit status 0); <paramref name="tbfastMessages"/> makes the messages those of tbfast, which stops with <c>exit( 1 )</c>.
        /// </summary>
        private MafftAddState? PrepareMerge(int njob, string[] seq, string[] bseq, int[][] topol0, int[][] topol1, MafftSubalignments subs, IReadOnlyList<string> names, bool tbfastMessages)
        {
            var mergeoralign = new char[njob - 1];
            Array.Fill(mergeoralign, 'a');
            var insubtable = new bool[njob];
            string Nm(int i) => FormatNamePrefix(names[i].Length > 0 ? names[i][1..] : "", 10);
            string OffsetNote()
            {
                int off = G.SubalignmentOffset;
                if (off == 0)
                    return "";
                string o = off.ToString(CultureInfo.InvariantCulture);
                return "#\n# You specified seed alignment(s) consisting of " + o + " sequences.\n# In this case, the rule of numbering is:\n#   The aligned seed sequences are numbered as 1 .. " + o + "\n#   The input sequences to be aligned are numbered as "
                    + (off + 1).ToString(CultureInfo.InvariantCulture) + " .. " + (off + njob).ToString(CultureInfo.InvariantCulture) + "\n";
            }
            for (int i = 0; i < subs.Groups.Count; i++)
            {
                var sub = subs.Groups[i];
                Report("Checking subalignment " + (i + 1).ToString(CultureInfo.InvariantCulture) + ":\n");
                int alignmentlength = seq[sub[0]].Length;
                foreach (int m in sub)
                {
                    if (m >= njob)
                        throw new MafftProgramFailureException("No such sequence, " + (m + 1).ToString(CultureInfo.InvariantCulture) + ".\n");
                    if (seq[m].Length != alignmentlength)
                    {
                        string text = "\n###############################################################################\n# ERROR!\n# Subalignment " + (i + 1).ToString(CultureInfo.InvariantCulture) + " must be aligned.\n# Please check the alignment lengths of following sequences.\n#\n"
                            + "# " + (sub[0] + 1).ToString(CultureInfo.InvariantCulture) + ". " + Nm(sub[0]) + " -> " + alignmentlength.ToString(CultureInfo.InvariantCulture) + " letters (including gaps)\n"
                            + "# " + (m + 1).ToString(CultureInfo.InvariantCulture) + ". " + Nm(m) + " -> " + seq[m].Length.ToString(CultureInfo.InvariantCulture) + " letters (including gaps)\n#\n"
                            + "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" + OffsetNote()
                            + "###############################################################################\n\n";
                        if (tbfastMessages)
                            throw new MafftProgramFailureException(text);
                        Report(text);
                        return null;
                    }
                    insubtable[m] = true;
                }
                for (int j = 0; j < njob - 1; j++)
                    if (IncludeMemberOf(topol0[j], sub) && IncludeMemberOf(topol1[j], sub))
                        mergeoralign[j] = 'n';
                bool found = false;
                for (int j = 0; j < njob - 1; j++)
                {
                    if (SameMemberAs(topol0[j], sub) || SameMemberAs(topol1[j], sub))
                    {
                        found = true;
                        Report(" -> OK\n");
                        break;
                    }
                }
                if (!found)
                {
                    // system( "cp infile.tree GuideTree" ): the copy fails (and says so on the standard error of the program) when no program has written the tree file
                    string cpMessage = TreeFile != null ? "" : "cp: cannot stat 'infile.tree': No such file or directory\n";
                    Report(cpMessage);
                    GuideTreeCopy = TreeFile;
                    string text = "\n###############################################################################\n# ERROR!\n# Subalignment " + (i + 1).ToString(CultureInfo.InvariantCulture)
                        + (tbfastMessages ? " does not form a monophyletic cluster\n" : " does not seem to form a monophyletic cluster\n")
                        + "# in the guide tree ('GuideTree' in this directory) internally computed.\n# If you really want to use this subalignment, pelase give a tree with --treein \n# http://mafft.cbrc.jp/alignment/software/treein.html\n# http://mafft.cbrc.jp/alignment/software/merge.html\n"
                        + OffsetNote() + "############################################################################### \n\n";
                    if (tbfastMessages)
                        throw new MafftProgramFailureException(text);
                    Report(text);
                    return null;
                }
            }
            for (int i = 0; i < njob; i++)
                bseq[i] = insubtable[i] ? seq[i] : seq[i].Replace("-", "", StringComparison.Ordinal);
            for (int i = 0; i < subs.Groups.Count; i++)
            {
                var rows = subs.Groups[i].Select(m => bseq[m]).ToArray();
                if (!subs.PreserveGaps[i])
                    CommonGapPickN(rows.Length, rows);
                for (int j = 0; j < rows.Length; j++)
                    bseq[subs.Groups[i][j]] = rows[j];
            }
            var alreadyAligned = new bool[njob];
            return new MafftAddState(mergeoralign, alreadyAligned);
        }

        /// <summary><c>%-10.10s</c>: at most ten characters, padded with blanks to ten.</summary>
        private static string FormatNamePrefix(string s, int width) => (s.Length > width ? s[..width] : s).PadRight(width);
    }

    #endregion
}
