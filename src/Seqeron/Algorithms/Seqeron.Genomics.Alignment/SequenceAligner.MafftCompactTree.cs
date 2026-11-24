#pragma warning disable IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854, S3260
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: the guide trees of disttbfast that are not the UPGMA of the full matrix (--treein, --pileup, --randomchain, --large / --memsavetree, --memsavetreex, --stepadd, --youngestlinkage)

    internal sealed partial class MafftRun
    {
        /// <summary>The text of the file <c>_guidetree</c> that the script writes for <c>-U</c> (<c>--treein</c>, <c>--pileup</c>, <c>--randomchain</c>, <c>--large</c> ...): the first line names the method, the other lines are the nodes of a tree given by the user.</summary>
        internal string GuideTreeText = "";

        /// <summary>The text that the program wrote to standard output before the alignment.</summary>
        internal string StdoutPrefix = "";

        /// <summary><c>compacttree</c> of <c>disttbfast.c</c>: 0 the UPGMA of the full matrix, 1 the compact tree with a limit of RAM, 2 the compact tree with the nearest earlier neighbour, 4 the youngest linkage.</summary>
        private int _compactTree;

        /// <summary><c>maxdistmtxsize</c>: the limit of RAM of <c>compact</c> (bytes).</summary>
        private double _maxDistMtxSize = 10.0 * 1000 * 1000 * 1000;

        /// <summary>The UB of the program that the port cannot copy: it reads memory that does not belong to the data (the reason is the message).</summary>
        private static NotSupportedException UndefinedBehaviour(string what) => new("The program has undefined behaviour here and its result is not reproducible: " + what);

        // -------------------------------------------------------------------- scanf and fgets as the programs use them

        /// <summary><c>fgets( buf, 999, fp )</c> over a text: the next line with its newline, at most 998 characters; null at the end of the file.</summary>
        private sealed class FgetsReader
        {
            private readonly string _text;
            private int _pos;

            internal FgetsReader(string text) => _text = text;

            internal string? Next()
            {
                if (_pos >= _text.Length)
                    return null;
                int end = _pos;
                int count = 0;
                while (end < _text.Length && count < 998)
                {
                    char c = _text[end++];
                    count++;
                    if (c == '\n')
                        break;
                }
                string line = _text.Substring(_pos, end - _pos);
                _pos = end;
                return line;
            }
        }

        /// <summary>The conversions <c>%d</c> and <c>%lf</c> of glibc <c>sscanf</c> over a string: leading white space is skipped, the characters that can belong to the number are collected greedily (<c>vfscanf</c>), then <c>strtol</c> / <c>strtod</c> convert the longest valid prefix of them (the rest is consumed and ignored: "1e" is read as 1).</summary>
        internal sealed class CScanner
        {
            private readonly string _s;
            private int _pos;

            internal CScanner(string s, int start = 0)
            {
                _s = s;
                _pos = start;
            }

            private void SkipSpace()
            {
                while (_pos < _s.Length && _s[_pos] is ' ' or '\t' or '\n' or '\r' or '\v' or '\f')
                    _pos++;
            }

            private char Peek => _pos < _s.Length ? _s[_pos] : '\0';

            private bool More => _pos < _s.Length;

            /// <summary><c>%d</c>: false when no integer follows.</summary>
            internal bool Int(out int value)
            {
                value = 0;
                SkipSpace();
                bool neg = false;
                if (More && Peek is '+' or '-')
                {
                    neg = Peek == '-';
                    _pos++;
                }
                int digits = 0;
                System.Numerics.BigInteger v = 0;
                while (More && Peek is >= '0' and <= '9')
                {
                    v = v * 10 + (Peek - '0');
                    digits++;
                    _pos++;
                }
                if (digits == 0)
                    return false;
                if (neg)
                    v = -v;
                long l = v > long.MaxValue ? long.MaxValue : v < long.MinValue ? long.MinValue : (long)v;
                value = unchecked((int)l);
                return true;
            }

            /// <summary><c>%lf</c>: false when no number follows.</summary>
            internal bool Double(out double value)
            {
                value = 0;
                SkipSpace();
                var buf = new StringBuilder();
                if (More && Peek is '+' or '-')
                {
                    buf.Append(Peek);
                    _pos++;
                }
                bool hex = false;
                // nan / inf
                if (More && (Peek is 'n' or 'N'))
                {
                    if (!Word("nan"))
                        return false;
                    if (More && Peek == '(')
                    {
                        int k = _pos + 1;
                        while (k < _s.Length && (char.IsAsciiLetterOrDigit(_s[k]) || _s[k] == '_'))
                            k++;
                        if (k < _s.Length && _s[k] == ')')
                            _pos = k + 1;
                    }
                    value = BitConverter.Int64BitsToDouble(0x7FF8000000000000L); // strtod("nan") is a positive NaN
                    if (buf.Length > 0 && buf[0] == '-')
                        value = -value;
                    return true;
                }
                if (More && (Peek is 'i' or 'I'))
                {
                    if (!Word("inf"))
                        return false;
                    int save = _pos;
                    if (!Word("inity"))
                        _pos = save;
                    value = buf.Length > 0 && buf[0] == '-' ? double.NegativeInfinity : double.PositiveInfinity;
                    return true;
                }
                if (More && Peek == '0' && _pos + 1 < _s.Length && _s[_pos + 1] is 'x' or 'X')
                {
                    hex = true;
                    buf.Append("0x");
                    _pos += 2;
                }
                bool gotDigit = false, gotDot = false, gotE = false;
                char expChar = hex ? 'p' : 'e';
                while (More)
                {
                    char c = Peek;
                    if (c is >= '0' and <= '9' || (hex && Uri.IsHexDigit(c) && !gotE))
                    {
                        buf.Append(c);
                        gotDigit = true;
                    }
                    else if (gotE && buf.Length > 0 && char.ToLowerInvariant(buf[^1]) == expChar && c is '-' or '+')
                    {
                        buf.Append(c);
                    }
                    else if (gotDigit && !gotE && char.ToLowerInvariant(c) == expChar)
                    {
                        buf.Append(c);
                        gotE = gotDot = true;
                    }
                    else if (!gotDot && c == '.')
                    {
                        buf.Append(c);
                        gotDot = true;
                    }
                    else
                    {
                        break;
                    }
                    _pos++;
                }
                return StrToD(buf.ToString(), out value);
            }

            private bool Word(string w)
            {
                for (int i = 0; i < w.Length; i++)
                {
                    if (_pos >= _s.Length || char.ToLowerInvariant(_s[_pos]) != w[i])
                        return false;
                    _pos++;
                }
                return true;
            }

            /// <summary><c>strtod</c> of the collected characters: the longest prefix that is a number; false when there is none.</summary>
            private static bool StrToD(string t, out double value)
            {
                value = 0;
                int p = 0;
                bool neg = false;
                if (p < t.Length && t[p] is '+' or '-')
                {
                    neg = t[p] == '-';
                    p++;
                }
                if (p + 1 < t.Length && t[p] == '0' && t[p + 1] is 'x' or 'X')
                {
                    int h = p + 2;
                    double mant = 0;
                    int hexDigits = 0;
                    int exp2 = 0;
                    while (h < t.Length && Uri.IsHexDigit(t[h]))
                    {
                        mant = mant * 16 + Convert.ToInt32(t[h].ToString(), 16);
                        hexDigits++;
                        h++;
                    }
                    if (h < t.Length && t[h] == '.')
                    {
                        int h2 = h + 1;
                        int fd = 0;
                        while (h2 < t.Length && Uri.IsHexDigit(t[h2]))
                        {
                            mant = mant * 16 + Convert.ToInt32(t[h2].ToString(), 16);
                            exp2 -= 4;
                            fd++;
                            h2++;
                        }
                        if (hexDigits + fd > 0)
                        {
                            hexDigits += fd;
                            h = h2;
                        }
                    }
                    if (hexDigits == 0)
                    {
                        // "0x" without digits: strtod reads the "0"
                        value = neg ? -0.0 : 0.0;
                        return true;
                    }
                    if (h < t.Length && t[h] is 'p' or 'P')
                    {
                        int k = h + 1;
                        bool eneg = false;
                        if (k < t.Length && t[k] is '+' or '-')
                        {
                            eneg = t[k] == '-';
                            k++;
                        }
                        int ed = 0;
                        int ev = 0;
                        while (k < t.Length && t[k] is >= '0' and <= '9')
                        {
                            ev = Math.Min(ev * 10 + (t[k] - '0'), 100000);
                            ed++;
                            k++;
                        }
                        if (ed > 0)
                            exp2 += eneg ? -ev : ev;
                    }
                    double r = mant * Math.Pow(2, exp2);
                    value = neg ? -r : r;
                    return true;
                }
                int d = p;
                int nd = 0;
                while (d < t.Length && t[d] is >= '0' and <= '9')
                {
                    d++;
                    nd++;
                }
                if (d < t.Length && t[d] == '.')
                {
                    int d2 = d + 1;
                    int fd = 0;
                    while (d2 < t.Length && t[d2] is >= '0' and <= '9')
                    {
                        d2++;
                        fd++;
                    }
                    if (nd + fd > 0)
                    {
                        d = d2;
                        nd += fd;
                    }
                }
                if (nd == 0)
                    return false;
                if (d < t.Length && t[d] is 'e' or 'E')
                {
                    int k = d + 1;
                    if (k < t.Length && t[k] is '+' or '-')
                        k++;
                    int ed = 0;
                    while (k < t.Length && t[k] is >= '0' and <= '9')
                    {
                        k++;
                        ed++;
                    }
                    if (ed > 0)
                        d = k;
                }
                value = double.Parse(t.AsSpan(p, d - p), NumberStyles.Float, CultureInfo.InvariantCulture);
                if (neg)
                    value = -value;
                return true;
            }
        }

        // -------------------------------------------------------------------- check_guidetreefile

        /// <summary><c>check_guidetreefile</c> of <c>mltaln9.c</c>: the method named by the first line of <c>_guidetree</c>.</summary>
        private (char Kind, int Seed, int NPick, double LimitRam) CheckGuideTreeFile()
        {
            int seed = 0, npick = 200;
            double limitram = 10.0 * 1000 * 1000 * 1000;
            var fp = new FgetsReader(GuideTreeText);
            string? line = fp.Next();
            if (line == null)
            {
                // fgets returns NULL and leaves string[100] as it was: stack content that matches none of the keywords (repeated runs of the installed program: "loadtree." and the format error of the empty file, 5 / 5)
                Report("loadtree.\n");
                return ('l', seed, npick, limitram);
            }
            if (line.Length > 99)
                throw UndefinedBehaviour("the first line of _guidetree is longer than 98 characters: fgets( string, 999, fp ) overflows the buffer string[100] and the program aborts with a message that names a process id");
            string str = line;
            if (str.StartsWith("shuffle", StringComparison.Ordinal))
            {
                new CScanner(str, 7).Int(out seed);
                Report("shuffle, seed=" + seed.ToString(CultureInfo.InvariantCulture) + "\n");
                return ('s', seed, npick, limitram);
            }
            if (str.StartsWith("pileup", StringComparison.Ordinal))
            {
                Report("pileup.\n");
                return ('p', seed, npick, limitram);
            }
            if (str.StartsWith("auto", StringComparison.Ordinal))
            {
                var sc = new CScanner(str, 4);
                if (sc.Int(out int a))
                {
                    seed = a;
                    if (sc.Int(out int b))
                        npick = b;
                }
                Report("auto, seed=" + seed.ToString(CultureInfo.InvariantCulture) + ", npick=" + npick.ToString(CultureInfo.InvariantCulture) + "\n");
                if (npick < 2)
                    throw new MafftProgramFailureException("Check npick\n");
                return ('a', seed, npick, limitram);
            }
            if (str.StartsWith("test", StringComparison.Ordinal))
            {
                var sc = new CScanner(str, 4);
                if (sc.Int(out int a))
                {
                    seed = a;
                    if (sc.Int(out int b))
                        npick = b;
                }
                Report("calc, seed=" + seed.ToString(CultureInfo.InvariantCulture) + ", npick=" + npick.ToString(CultureInfo.InvariantCulture) + "\n");
                if (npick < 2)
                    throw new MafftProgramFailureException("Check npick\n");
                return ('t', seed, npick, limitram);
            }
            if (str.StartsWith("compact", StringComparison.Ordinal))
            {
                string sizestring = str[7..];
                Report("sizestring = " + sizestring + "\n");
                double tanni;
                if (sizestring.Contains('k'))
                    tanni = 1.0 * 1000;
                else if (sizestring.Contains('M') || sizestring.Contains('m'))
                    tanni = 1.0 * 1000 * 1000;
                else if (sizestring.Contains('G') || sizestring.Contains('g'))
                    tanni = 1.0 * 1000 * 1000 * 1000;
                else if (sizestring.Contains('T') || sizestring.Contains('t'))
                    tanni = 1.0 * 1000 * 1000 * 1000 * 1000;
                else
                    throw new MafftProgramFailureException("\nSpecify initial ram usage by '--initialramusage xGB'\n\n\n");
                double tmpd = 0.0;
                new CScanner(sizestring).Double(out tmpd);
                limitram = tmpd * tanni;
                Report("Initial RAM usage = " + MafftGlobals.F(limitram / 1000 / 1000 / 1000, 10, 3) + "GB\n");
                return ('c', seed, npick, limitram);
            }
            if (str.StartsWith("very compact", StringComparison.Ordinal))
            {
                Report("very compact.\n");
                return ('C', seed, npick, limitram);
            }
            if (str.StartsWith("stepadd", StringComparison.Ordinal))
            {
                Report("stepwise addition (disttbfast).\n");
                return ('S', seed, npick, limitram);
            }
            if (str.StartsWith("youngestlinkage", StringComparison.Ordinal))
            {
                Report("youngest linkage (disttbfast).\n");
                return ('Y', seed, npick, limitram);
            }
            if (str.StartsWith("nodepair", StringComparison.Ordinal))
            {
                Report("Use nodepair.\n");
                return ('n', seed, npick, limitram);
            }
            Report("loadtree.\n");
            return ('l', seed, npick, limitram);
        }

        // -------------------------------------------------------------------- the leaf names of the tree files

        private string[] TreeLeafStrings(int nseq, IReadOnlyList<string> names)
        {
            var tree = new string[nseq];
            for (int i = 0; i < nseq; i++)
                tree[i] = "\n" + (i + 1).ToString(CultureInfo.InvariantCulture) + "_" + TreeLabel(names[i]) + "\n";
            return tree;
        }

        private char VarPairScoreChoice(int nseq, int npick, string[] seq, int seed)
        {
            char c = VarPairScore(nseq, npick, seq, seed, out string text);
            StdoutPrefix = text;
            return c;
        }

        /// <summary><c>varpairscore</c> of <c>disttbfast.c</c> (the guide tree file <c>auto</c> / <c>test</c>): the scores of the pairs of <paramref name="npick"/> random sequences decide between a tree ('t') and a chain ('s', when the coefficient of variation is below 0.2).</summary>
        private char VarPairScore(int nseq, int npick, string[] seq, int seed, out string stdout)
        {
            var sb = new StringBuilder();
            Report("nseq = " + nseq.ToString(CultureInfo.InvariantCulture) + ", nlenmax=" + _nlenmax.ToString(CultureInfo.InvariantCulture) + ", seed=" + seed.ToString(CultureInfo.InvariantCulture) + "\n");
            if (npick > nseq)
                throw UndefinedBehaviour("more sequences are picked than there are: the program reads the list of sequences beyond its end");
            var rand = new GlibcRandom(unchecked((uint)seed));
            var slist = new int[nseq];
            for (int i = 0; i < nseq; i++)
                slist[i] = i;
            StringShuffle(slist, nseq, rand);
            var pickseq = new string[npick];
            for (int i = 0; i < npick; i++)
                pickseq[i] = seq[slist[i]].Replace("-", "");
            if (_nlenmax % 16 == 8 && pickseq.Any(x => x.Length == _nlenmax))
                throw UndefinedBehaviour("the program writes the terminating NUL of a picked sequence of the longest length one byte behind its buffer (calloc( nlenmax )): the next block of the heap is damaged and the program may abort with a message that names a process id");
            double scoreav = 0.0, scoreav2 = 0.0;
            int npair = npick * (npick - 1) / 2;
            for (int i = 1; i < npick; i++)
            {
                Report(i.ToString(CultureInfo.InvariantCulture) + " / " + npick.ToString(CultureInfo.InvariantCulture) + "\r");
                for (int j = 0; j < i; j++)
                {
                    double score = G__align11NoalignWith(G.NDisConsweightMulti, -1200.0, -60.0, pickseq[i], pickseq[j]);
                    scoreav += score;
                    scoreav2 += score * score;
                    sb.Append("score = ").Append(CInt(score).ToString(CultureInfo.InvariantCulture)).Append('\n');
                }
            }
            scoreav /= npair;
            scoreav2 /= npair;
            double scorevar = (scoreav2 - scoreav * scoreav) * npair / (npair - 1);
            double scorestd = Math.Sqrt(scorevar);
            sb.Append("av = ").Append(MafftGlobals.F(scoreav, 0, 6)).Append('\n');
            sb.Append("stddev = ").Append(MafftGlobals.F(scorestd, 0, 6)).Append('\n');
            sb.Append("cv = ").Append(MafftGlobals.F(scorestd / scoreav, 0, 6)).Append('\n');
            stdout = sb.ToString();
            return scorestd / scoreav < 0.2 ? 's' : 't';
        }

        // -------------------------------------------------------------------- loadtree

        private static (int A, int B) LoadTreeOneLine(FgetsReader fp, double[] len)
        {
            const string Hint = "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n";
            string? gett = fp.Next();
            if (gett == null)
                throw new MafftProgramFailureException("\n\nFormat error (1) in the tree?  It has to be a bifurcated and rooted tree.\n" + Hint);
            var sc = new CScanner(gett);
            int res = 0;
            int a = 0, b = 0;
            double l0 = len[0], l1 = len[1];
            if (sc.Int(out a))
            {
                res = 1;
                if (sc.Int(out b))
                {
                    res = 2;
                    if (sc.Double(out double d0))
                    {
                        l0 = d0;
                        res = 3;
                        if (sc.Double(out double d1))
                        {
                            l1 = d1;
                            res = 4;
                        }
                    }
                }
            }
            len[0] = l0;
            len[1] = l1;
            if (res != 4)
                throw new MafftProgramFailureException("\n\nFormat error (2) in the tree?  It has to be a bifurcated and rooted tree.\n" + Hint);
            a--;
            b--;
            if (a >= b)
                throw new MafftProgramFailureException("\n\nIncorrect guide tree\n" + Hint);
            return (a, b);
        }

        /// <summary><c>loadtree</c> (<c>--treein</c>): the tree of <c>_guidetree</c> (one line per node: two sequence numbers and two branch lengths); the groups keep their full membership.</summary>
        private void LoadTree(int nseq, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool treeout, IReadOnlyList<string> names)
        {
            Report("Loading a tree\n");
            var fp = new FgetsReader(GuideTreeText);
            var hist = new int[nseq];
            var nmemar = new int[nseq];
            var height = new double[nseq];
            string[]? tree = treeout ? TreeLeafStrings(nseq, names) : null;
            string treetmp = "";
            var next = new int[nseq];
            var prev = new int[nseq];
            for (int i = 0; i < nseq; i++)
            {
                next[i] = i + 1;
                prev[i] = i - 1;
                hist[i] = -1;
                nmemar[i] = 1;
            }
            next[nseq - 1] = -1;
            Report("\n");
            for (int k = 0; k < nseq - 1; k++)
            {
                if (k % 10 == 0)
                    Report("\r" + MafftGlobals.D5(k) + " / " + nseq);
                len[k][0] = len[k][1] = -1.0;
                var (im, jm) = LoadTreeOneLine(fp, len[k]);
                if (im > nseq - 1 || jm > nseq - 1)
                    throw new MafftProgramFailureException("\n\nCheck the guide tree.\nim=" + (im + 1).ToString(CultureInfo.InvariantCulture) + ", jm=" + (jm + 1).ToString(CultureInfo.InvariantCulture)
                        + "\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n\n");
                if (len[k][0] == -1.0 || len[k][1] == -1.0)
                    throw new MafftProgramFailureException("\n\nERROR: Branch length is not given.\n");
                if (len[k][0] < 0.0)
                    len[k][0] = 0.0;
                if (len[k][1] < 0.0)
                    len[k][1] = 0.0;
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
                if (tree != null)
                {
                    treetmp = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                    tree[im] = treetmp;
                    tree[jm] = null!;
                }
                int acjmprev = prev[jm];
                int acjmnext = next[jm];
                if (acjmprev >= 0)
                    next[acjmprev] = acjmnext;
                if (acjmnext != -1)
                    prev[acjmnext] = acjmprev;
                height[im] += len[k][0];
                dep[k].DistFromTip = height[im];
            }
            if (tree != null)
                TreeFile = treetmp + ";\n#by loadtree\n";
        }

        /// <summary><c>veryfastsupg_double_loadtree</c> (<c>dvtditr -U</c>): the tree of <c>_guidetree</c> for the iterative refinement; the half matrix <paramref name="eff"/> is updated as the program does (the values are not used afterwards), the groups keep their full membership, the file <c>infile.tree</c> is always written.</summary>
        private void VeryFastSupgLoadTree(int nseq, double[][] eff, int[][] topol0, int[][] topol1, double[][] len, IReadOnlyList<string> names)
        {
            var fp = new FgetsReader(GuideTreeText);
            string[] tree = TreeLeafStrings(nseq, names);
            var hist = new int[nseq];
            var next = new int[nseq];
            var prev = new int[nseq];
            for (int i = 0; i < nseq; i++)
            {
                next[i] = i + 1;
                prev[i] = i - 1;
                hist[i] = -1;
            }
            next[nseq - 1] = -1;
            string treetmp = "";
            Report("\n");
            var lenfl = new double[2];
            for (int k = 0; k < nseq - 1; k++)
            {
                if (k % 10 == 0)
                    Report(k.ToString(CultureInfo.InvariantCulture) + " / " + nseq.ToString(CultureInfo.InvariantCulture) + "\r");
                lenfl[0] = lenfl[1] = -1.0;
                var (im, jm) = LoadTreeOneLine(fp, lenfl);
                if (im < 0 || im > nseq - 1)
                    throw UndefinedBehaviour("a node names a sequence outside the matrix of distances: the program reads eff[" + (im + 1).ToString(CultureInfo.InvariantCulture) + "]");
                if (jm > nseq - 1 || tree[im] == null || jm < 0 || tree[jm] == null)
                    throw new MafftProgramFailureException("\n\nCheck the guide tree.\nim=" + (im + 1).ToString(CultureInfo.InvariantCulture) + ", jm=" + (jm + 1).ToString(CultureInfo.InvariantCulture)
                        + "\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n\n");
                if (lenfl[0] == -1.0 || lenfl[1] == -1.0)
                    throw new MafftProgramFailureException("\n\nWARNING: Branch length is not given.\n");
                if (lenfl[0] < 0.0)
                    lenfl[0] = 0.0;
                if (lenfl[1] < 0.0)
                    lenfl[1] = 0.0;
                topol0[k] = GroupMembers(hist[im], im, topol0, topol1, true);
                topol1[k] = GroupMembers(hist[jm], jm, topol0, topol1, true);
                len[k][0] = lenfl[0];
                len[k][1] = lenfl[1];
                hist[im] = k;
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
                        eff[miniim][maxiim - miniim] = Math.Min(eff0, eff1) * (1.0 - G.SueffGlobal) + (eff0 + eff1) * 0.5 * G.SueffGlobal;
                    }
                }
                int acjmprev = prev[jm];
                int acjmnext = next[jm];
                next[acjmprev] = acjmnext;
                if (acjmnext != -1)
                    prev[acjmnext] = acjmprev;
                treetmp = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                tree[im] = treetmp;
                tree[jm] = null!;
            }
            TreeFile = treetmp + ";\n";
            Report("\n");
        }

        // -------------------------------------------------------------------- createchain

        /// <summary><c>stringshuffle</c>: <c>rand() % size</c> swaps.</summary>
        private static void StringShuffle(int[] ary, int size, GlibcRandom rand)
        {
            for (int i = 0; i < size; i++)
            {
                int j = rand.Next() % size;
                (ary[i], ary[j]) = (ary[j], ary[i]);
            }
        }

        /// <summary><c>createchain</c> (<c>--pileup</c> in the order of the input, <c>--randomchain</c> after <c>stringshuffle</c>): the sequences are added one after the other to one group; the tree is a chain. The groups keep one member (MEMSAVE).</summary>
        private void CreateChain(int nseq, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool treeout, bool shuffle, int seed, IReadOnlyList<string> names)
        {
            string[]? tree = treeout ? TreeLeafStrings(nseq, names) : null;
            var instant = new StringBuilder();
            if (treeout)
                instant.Append('(', nseq - 1);
            var order = new int[nseq];
            for (int i = 0; i < nseq; i++)
                order[i] = i;
            var rand = new GlibcRandom(unchecked((uint)seed));
            if (shuffle)
                StringShuffle(order, nseq, rand);
            double l, ll;
            ll = l = 2.0 / nseq;
            int im = order[0];
            int jm = order[1];
            int mm;
            if (im < jm)
            {
                topol0[0] = new[] { im };
                topol1[0] = new[] { jm };
                mm = im;
            }
            else
            {
                topol0[0] = new[] { jm };
                topol1[0] = new[] { im };
                mm = jm;
            }
            len[0][0] = len[0][1] = l;
            dep[0].Child1 = -1;
            dep[0].Child0 = -1;
            dep[0].DistFromTip = l;
            ll += l;
            if (treeout)
                instant.Append(tree![im]).Append(':').Append(MafftGlobals.F(len[0][0], 7, 5)).Append(',');
            for (int i = 1; i < nseq - 1; i++)
            {
                im = order[i];
                jm = order[i + 1];
                if (mm < jm)
                {
                    topol0[i] = new[] { mm };
                    topol1[i] = new[] { jm };
                    len[i][0] = l;
                    len[i][1] = ll;
                    dep[i].Child0 = i - 1;
                    dep[i].Child1 = -1;
                    dep[i].DistFromTip = ll;
                }
                else
                {
                    topol1[i] = new[] { mm };
                    topol0[i] = new[] { jm };
                    mm = jm;
                    len[i][1] = l;
                    len[i][0] = ll;
                    dep[i].Child1 = i - 1;
                    dep[i].Child0 = -1;
                    dep[i].DistFromTip = ll;
                }
                if (treeout)
                    instant.Append(tree![im]).Append(':').Append(MafftGlobals.F(ll - l, 7, 5)).Append("):").Append(MafftGlobals.F(l, 7, 5)).Append(',');
                ll += l;
            }
            if (treeout)
            {
                instant.Append(tree![jm]).Append(':').Append(MafftGlobals.F(ll - l, 7, 5)).Append(')');
                TreeFile = instant + ";\n";
            }
            // the program writes the chain as a tree to _guidetree (CANONICALTREEFORMAT), which the next program of the script (dvtditr -U) reads
            var canonical = new StringBuilder();
            for (int i = 0; i < nseq - 1; i++)
                canonical.Append(topol0[i][0] + 1).Append(' ').Append(topol1[i][0] + 1).Append(' ').Append(MafftGlobals.F(len[i][0], 0, 6)).Append(' ').Append(MafftGlobals.F(len[i][1], 0, 6)).Append('\n');
            GuideTreeText = canonical.ToString();
        }

        // -------------------------------------------------------------------- the compact trees

        /// <summary><c>preferenceval</c>: a tie-break that grows with the distance of the index, in the cyclic order.</summary>
        private static double PreferenceVal(int ori, int pos, int max)
        {
            pos -= ori;
            if (pos < 0)
                pos += max;
            return 0.00000000000001 * pos;
        }

        /// <summary><c>distcompact</c>: the distance of two sequences from their common k-mers (table1 of the first, pointt of the second) and the self scores.</summary>
        private double DistCompact(int len1, int len2, int[] table1, int[] point2, int ss1, int ss2, int[] memo)
        {
            double longer, shorter;
            if (len1 > len2)
            {
                longer = len1;
                shorter = len2;
            }
            else
            {
                longer = len2;
                shorter = len1;
            }
            double lenfac = 1.0 / (shorter / longer * _lenfacd + _lenfacb / (longer + _lenfacc) + _lenfaca);
            if (ss1 == 0 || ss2 == 0)
                return 2.0;
            return (1.0 - (double)CommonSextet(table1, point2, memo) / Math.Min(ss1, ss2)) * lenfac * 2.0;
        }

        /// <summary><c>distcompact_msa</c>: the distance of two rows of an alignment from their score and the self scores.</summary>
        private double DistCompactMsa(string seq1, string seq2, int ss1, int ss2)
        {
            int bunbo = Math.Min(ss1, ss2);
            if (bunbo == 0)
                return 2.0;
            double value = (1.0 - NaivePairScoreStream(seq1, seq2, G.PenaltyDist) / bunbo) * 2.0;
            if (value > 10)
                value = 10.0;
            return value;
        }

        /// <summary>The score of <see cref="NaivePairScore11"/> (the columns where both rows have a gap are left out, a run of gaps of one row costs <paramref name="penal"/> once) without the copies of the rows.</summary>
        private double NaivePairScoreStream(string seq1, string seq2, int penal)
        {
            double vali = 0.0;
            int mode = 0;
            var dis = G.AminoDis;
            for (int a = 0; a < seq1.Length; a++)
            {
                char x = seq1[a], y = seq2[a];
                if (x == '-')
                {
                    if (y == '-')
                        continue;
                    if (mode != 1)
                    {
                        vali += penal;
                        mode = 1;
                    }
                }
                else if (y == '-')
                {
                    if (mode != 2)
                    {
                        vali += penal;
                        mode = 2;
                    }
                }
                else
                {
                    vali += dis[x][y];
                    mode = 0;
                }
            }
            return vali;
        }

        /// <summary><c>preparepartmtx</c>: the rows of the matrix of distances that fit into <c>maxdistmtxsize</c> (only the messages matter: the matrix is a cache); </summary>
        private void PreparePartMtx(int nseq)
        {
            if (_compactTree == 1)
            {
                double size = 0;
                int i;
                for (i = 0; i < nseq; i++)
                {
                    size += (double)sizeof(double) * nseq;
                    if (size > _maxDistMtxSize)
                    {
                        Report("\n\nThe size of full distance matrix is estimated to exceed " + MafftGlobals.F(_maxDistMtxSize / 1000 / 1000 / 1000, 0, 2) + "GB.\n");
                        Report("Will try the calculation using a " + nseq.ToString(CultureInfo.InvariantCulture) + " x " + i.ToString(CultureInfo.InvariantCulture) + " matrix.\n");
                        Report("This calculation will be slow due to the limited RAM space.\n");
                        Report("To avoid the slowdown, please try '--initialramusage xGB' (x>>" + MafftGlobals.F(_maxDistMtxSize / 1000 / 1000 / 1000, 0, 2) + "),\n");
                        Report("if larger RAM space is available.\n");
                        Report("Note that xGB is NOT the upper limit of RAM usage.\n");
                        Report("Two to three times larger space may be used for building a guide tree.\n");
                        Report("Memory usage of the MSA stage depends on similarity of input sequences.\n\n");
                        Report("The '--memsavetree' option uses smaller RAM space.\n");
                        Report("If tree-like relationship can be ignored, try '--pileup' or '--randomchain'.\n\n");
                        Report("The result of --initialramusage xGB is almost identical to the default, except for rounding differences.\n");
                        Report("In the cases of --memsavetree, --pileup and --randomchain, the result differs from the default.\n\n");
                        break;
                    }
                }
                if (i == nseq)
                    Report("The full matrix will be used.\n");
            }
        }

        /// <summary>The nearest earlier (or, for the youngest linkage, any) sequence of every sequence: <c>compactdisthalfmtxthread</c> / <c>ylcompactdisthalfmtxthread</c> of <c>disttbfast.c</c>.</summary>
        private void CompactDistancesKmer(int njob, int[][] pointt, int[] selfscore, double[] mindist, int[] mindistfrom)
        {
            var memo = new int[_tsize];
            bool yl = _compactTree == 4;
            int jobpos = yl ? 0 : njob - 1;
            while (true)
            {
                AlignmentCancellation.ThrowIfRequested();
                int i = jobpos;
                if (yl)
                {
                    if (i == njob - 1)
                        return;
                    jobpos = i + 1;
                }
                else
                {
                    if (i == -1)
                        return;
                    jobpos = i - 1;
                }
                var table1 = new int[_tsize];
                if (i % 100 == 0)
                    Report("\r" + MafftGlobals.D5(yl ? i + 1 : njob - i) + " / " + njob);
                foreach (int point in pointt[i])
                    table1[point]++;
                if (yl)
                {
                    for (int j = i + 1; j < njob; j++)
                    {
                        double tmpdist = DistCompact(_nogaplen[i], _nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j], memo);
                        double tmpdistx = tmpdist + PreferenceVal(i, j, njob);
                        if (tmpdistx < mindist[i])
                        {
                            mindist[i] = tmpdistx;
                            mindistfrom[i] = j;
                        }
                        double tmpdisty = tmpdist + PreferenceVal(j, i, njob);
                        if (tmpdisty < mindist[j])
                        {
                            mindist[j] = tmpdisty;
                            mindistfrom[j] = i;
                        }
                    }
                }
                else
                {
                    for (int j = i - 1; j > -1; j--)
                    {
                        double tmpdist = DistCompact(_nogaplen[i], _nogaplen[j], table1, pointt[j], selfscore[i], selfscore[j], memo);
                        double tmpdistx = tmpdist + PreferenceVal(i, j, njob);
                        if (tmpdistx < mindist[i])
                        {
                            mindist[i] = tmpdistx;
                            mindistfrom[i] = j;
                        }
                    }
                }
            }
        }

        /// <summary><c>msacompactdisthalfmtxthread</c> / <c>ylmsacompactdisthalfmtxthread</c>: the same from the rows of an alignment.</summary>
        private void CompactDistancesMsa(int njob, string[] seq, int[] selfscore, double[] mindist, int[] mindistfrom)
        {
            bool yl = _compactTree == 4;
            int jobpos = yl ? 0 : njob - 1;
            while (true)
            {
                AlignmentCancellation.ThrowIfRequested();
                int i = jobpos;
                if (yl)
                {
                    if (i == njob - 1)
                        return;
                    jobpos = i + 1;
                }
                else
                {
                    if (i == -1)
                        return;
                    jobpos = i - 1;
                }
                if (i % 100 == 0)
                    Report("\r" + MafftGlobals.D5(i) + " / " + njob);
                if (yl)
                {
                    for (int j = i + 1; j < njob; j++)
                    {
                        double tmpdist = DistCompactMsa(seq[i], seq[j], selfscore[i], selfscore[j]);
                        double tmpdistx = tmpdist + PreferenceVal(i, j, njob);
                        if (tmpdistx < mindist[i])
                        {
                            mindist[i] = tmpdistx;
                            mindistfrom[i] = j;
                        }
                        double tmpdisty = tmpdist + PreferenceVal(j, i, njob);
                        if (tmpdisty < mindist[j])
                        {
                            mindist[j] = tmpdisty;
                            mindistfrom[j] = i;
                        }
                    }
                }
                else
                {
                    for (int j = i - 1; j > -1; j--)
                    {
                        double tmpdist = DistCompactMsa(seq[i], seq[j], selfscore[i], selfscore[j]);
                        double tmpdistx = tmpdist + PreferenceVal(i, j, njob);
                        if (tmpdistx < mindist[i])
                        {
                            mindist[i] = tmpdistx;
                            mindistfrom[i] = j;
                        }
                    }
                }
            }
        }

        private sealed class TreePoint
        {
            internal int Child0 = -1, Child1 = -1, Parent = -1;
            internal int Rep0, Rep1 = -1;
            internal double Len0, Len1, Height;
        }

        /// <summary><c>compacttreegivendist</c>: the sequences are added one after the other (from the second) above the nearest earlier sequence at the height of their distance; <c>reformattree</c> turns the tree into the steps of the progressive alignment.</summary>
        private void CompactTreeGivenDist(int njob, double[] mindists, int[] neighbors, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool treeout, IReadOnlyList<string> names)
        {
            var treept = new TreePoint[njob * 2];
            for (int i = 0; i < treept.Length; i++)
                treept[i] = new TreePoint();
            for (int i = 0; i < njob; i++)
                treept[i].Rep0 = i;
            int n = njob;
            double tmpdist = mindists[1];
            treept[0].Parent = n;
            treept[1].Parent = n;
            treept[n].Child0 = 0;
            treept[n].Child1 = 1;
            treept[n].Height = tmpdist * 0.5;
            treept[n].Len0 = tmpdist * 0.5;
            treept[n].Len1 = tmpdist * 0.5;
            treept[n].Parent = -1;
            treept[n].Rep0 = 0;
            treept[n].Rep1 = 1;
            int root = n;
            for (int i = 2; i < njob; i++)
            {
                n++;
                int neighbor = neighbors[i];
                double mindist = mindists[i];
                int b = neighbor;
                int p = treept[b].Parent;
                for (; p != -1; b = p, p = treept[p].Parent)
                {
                    if (treept[p].Height > mindist * 0.5)
                        break;
                }
                if (p == -1)
                {
                    treept[n].Parent = -1;
                    root = n;
                }
                else if (treept[p].Child0 == b)
                {
                    treept[p].Child0 = n;
                    treept[p].Len0 = treept[p].Height - mindist * 0.5;
                    treept[n].Parent = p;
                }
                else if (treept[p].Child1 == b)
                {
                    treept[p].Child1 = n;
                    treept[p].Len1 = treept[p].Height - mindist * 0.5;
                    treept[n].Parent = p;
                }
                else
                {
                    throw new MafftProgramFailureException("okashii\n");
                }
                treept[i].Parent = n;
                treept[b].Parent = n;
                treept[n].Child0 = b;
                treept[n].Child1 = i;
                treept[n].Height = mindist * 0.5;
                treept[n].Rep0 = treept[b].Rep0;
                treept[n].Rep1 = treept[i].Rep0;
                treept[n].Len0 = mindist * 0.5 - treept[b].Height;
                treept[n].Len1 = mindist * 0.5;
            }
            ReformatTree(root, treept, njob, topol0, topol1, len, dep, names, treeout);
        }

        /// <summary><c>reformattree</c> with <c>reformat_younger0_rec</c>: the nodes in the order children first (child0 before child1); a group is named by its smaller representative first.</summary>
        private void ReformatTree(int root, TreePoint[] ori, int n, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, IReadOnlyList<string> names, bool treeout)
        {
            if (treeout)
            {
                string[] leaf = TreeLeafStrings(n, names);
                var sb = new StringBuilder();
                // reformat_rec_newick: (child of the smaller representative first)
                var stack = new Stack<(int Node, int State)>();
                stack.Push((root, 0));
                var lenText = new Stack<string>();
                while (stack.Count > 0)
                {
                    var (node, _) = stack.Pop();
                    if (node < 0)
                    {
                        sb.Append(lenText.Pop());
                        continue;
                    }
                    var pt = ori[node];
                    if (pt.Rep1 == -1)
                    {
                        sb.Append(leaf[pt.Rep0]);
                        continue;
                    }
                    int first, second;
                    double l1, l2;
                    if (pt.Rep0 < pt.Rep1)
                    {
                        first = pt.Child0; second = pt.Child1; l1 = pt.Len0; l2 = pt.Len1;
                    }
                    else
                    {
                        first = pt.Child1; second = pt.Child0; l1 = pt.Len1; l2 = pt.Len0;
                    }
                    sb.Append('(');
                    // emit: first ":" l1 "," second ":" l2 ")"  (pushed in reverse order)
                    stack.Push((-1, 0));
                    lenText.Push(":" + MafftGlobals.F(l2, 7, 5) + ")");
                    stack.Push((second, 0));
                    stack.Push((-1, 0));
                    lenText.Push(":" + MafftGlobals.F(l1, 7, 5) + ",");
                    stack.Push((first, 0));
                }
                TreeFile = sb + ";\n";
            }
            var lastappear = new int[n];
            for (int i = 0; i < n; i++)
                lastappear[i] = -1;
            int pos = 0;
            // post-order walk (child0, child1, node)
            var work = new Stack<(int Node, bool Visited)>();
            work.Push((root, false));
            while (work.Count > 0)
            {
                var (node, visited) = work.Pop();
                var pt = ori[node];
                if (pt.Rep1 == -1)
                    continue;
                if (!visited)
                {
                    work.Push((node, true));
                    if (pt.Child1 != -1)
                        work.Push((pt.Child1, false));
                    if (pt.Child0 != -1)
                        work.Push((pt.Child0, false));
                    continue;
                }
                if (pt.Rep0 < pt.Rep1)
                {
                    topol0[pos] = new[] { pt.Rep0 };
                    topol1[pos] = new[] { pt.Rep1 };
                    len[pos][0] = pt.Len0;
                    len[pos][1] = pt.Len1;
                    dep[pos].Child0 = lastappear[pt.Rep0];
                    dep[pos].Child1 = lastappear[pt.Rep1];
                }
                else
                {
                    topol1[pos] = new[] { pt.Rep0 };
                    topol0[pos] = new[] { pt.Rep1 };
                    len[pos][1] = pt.Len0;
                    len[pos][0] = pt.Len1;
                    dep[pos].Child1 = lastappear[pt.Rep0];
                    dep[pos].Child0 = lastappear[pt.Rep1];
                }
                lastappear[pt.Rep0] = pos;
                lastappear[pt.Rep1] = pos;
                dep[pos].DistFromTip = pt.Height;
                pos++;
            }
        }

        /// <summary><c>compacttree_memsaveselectable</c> with <c>howcompact == 2</c>, <c>memsave == 1</c> (<c>--youngestlinkage</c>): the UPGMA (mixed linkage) over the nearest neighbours only; the distances to a new group are computed from the sequences again (the k-mers when <paramref name="seq"/> is null, else the rows of the alignment).</summary>
        private void CompactTreeSelectable(int nseq, int[] nearest, double[] mindist, int[][]? pointt, int[] tselfscore, string[]? seq, int[][] topol0, int[][] topol1, double[][] len, Treedep[] dep, bool treeout, IReadOnlyList<string> names, int numfilled)
        {
            double sueff1 = 1 - G.SueffGlobal;
            double sueff05 = G.SueffGlobal * 0.5;
            if (G.TreeMethod != 'X')
                throw new MafftProgramFailureException("Unknown treemethod, " + G.TreeMethod + "\n");
            double Cluster(double d1, double d2) => Math.Min(d1, d2) * sueff1 + (d1 + d2) * sueff05;
            var tmptmplen = new double[nseq];
            var hist = new int[nseq];
            var nmemar = new int[nseq];
            var next = new int[nseq];
            var prev = new int[nseq];
            string[]? tree = treeout ? TreeLeafStrings(nseq, names) : null;
            string treetmp = "";
            for (int i = 0; i < nseq; i++)
            {
                next[i] = i + 1;
                prev[i] = i - 1;
                tmptmplen[i] = 0.0;
                hist[i] = -1;
                nmemar[i] = 1;
            }
            next[nseq - 1] = -1;
            var memo = new int[_tsize];
            var joblist = new int[nseq];
            Report("\n");
            int im = -1, jm = -1;
            for (int k = 0; k < nseq - 1; k++)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (k % 100 == 0)
                    Report("\r" + MafftGlobals.D5(k) + " / " + nseq);
                double minscore = 999.9;
                for (int a = 0; next[a] != -1; a = next[a])
                {
                    int i = a;
                    if (mindist[i] < minscore)
                    {
                        im = i;
                        minscore = mindist[i];
                    }
                }
                jm = nearest[im];
                if (jm < im)
                    (im, jm) = (jm, im);
                int prevnode = hist[im];
                dep[k].Child0 = prevnode;
                int nmemim = nmemar[im];
                topol0[k] = GroupMembers(prevnode, im, topol0, topol1, false);
                prevnode = hist[jm];
                dep[k].Child1 = prevnode;
                int nmemjm = nmemar[jm];
                topol1[k] = GroupMembers(prevnode, jm, topol0, topol1, false);
                minscore *= 0.5;
                len[k][0] = minscore - tmptmplen[im];
                len[k][1] = minscore - tmptmplen[jm];
                dep[k].DistFromTip = minscore;
                tmptmplen[im] = minscore;
                hist[im] = k;
                nmemar[im] = nmemim + nmemjm;
                mindist[im] = 999.9;
                int[]? table1 = null, table2 = null;
                if (pointt != null)
                {
                    table1 = new int[_tsize];
                    table2 = new int[_tsize];
                    foreach (int point in pointt[im])
                        table1[point]++;
                    foreach (int point in pointt[jm])
                        table2[point]++;
                }
                var newarr = new double[nseq];
                int nactive = 0;
                for (int a = 0; a != -1; a = next[a])
                    joblist[nactive++] = a;
                if (k % 100 == 0)
                    Report(" (serial, nactive=" + nactive.ToString(CultureInfo.InvariantCulture) + ", nfilled=" + numfilled.ToString(CultureInfo.InvariantCulture) + ")             \r");
                for (int p = 0; p < nactive; p++)
                {
                    int i = joblist[p];
                    if (i == im || i == jm)
                        continue;
                    double tmpdist1, tmpdist2;
                    if (pointt != null)
                    {
                        tmpdist1 = DistCompact(_nogaplen[im], _nogaplen[i], table1!, pointt[i], tselfscore[im], tselfscore[i], memo);
                        tmpdist2 = DistCompact(_nogaplen[jm], _nogaplen[i], table2!, pointt[i], tselfscore[jm], tselfscore[i], memo);
                    }
                    else
                    {
                        tmpdist1 = DistCompactMsa(seq![im], seq[i], tselfscore[im], tselfscore[i]);
                        tmpdist2 = DistCompactMsa(seq[jm], seq[i], tselfscore[jm], tselfscore[i]);
                    }
                    double tmpdouble = Cluster(tmpdist1, tmpdist2);
                    newarr[i] = tmpdouble;
                    if (tmpdouble < mindist[i])
                    {
                        mindist[i] = tmpdouble;
                        nearest[i] = im;
                    }
                    if (nearest[i] == jm)
                        nearest[i] = im;
                }
                for (int a = 0; a != -1; a = next[a])
                {
                    int i = a;
                    if (i != im && i != jm && newarr[i] < mindist[im])
                    {
                        mindist[im] = newarr[i];
                        nearest[im] = i;
                    }
                }
                if (tree != null)
                {
                    treetmp = "(" + tree[im] + ":" + MafftGlobals.F(len[k][0], 7, 5) + "," + tree[jm] + ":" + MafftGlobals.F(len[k][1], 7, 5) + ")";
                    tree[im] = treetmp;
                    tree[jm] = null!;
                }
                int acjmprev = prev[jm];
                int acjmnext = next[jm];
                next[acjmprev] = acjmnext;
                if (acjmnext != -1)
                    prev[acjmnext] = acjmprev;
            }
            if (tree != null)
                TreeFile = treetmp + ";\n";
        }
    }

    #endregion
}
