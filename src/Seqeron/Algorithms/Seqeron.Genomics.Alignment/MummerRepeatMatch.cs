using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The switches of <c>repeat-match</c> (MUMmer 3.23).</summary>
public sealed record MummerRepeatMatchOptions
{
    /// <summary><c>-n #</c>: minimum length of an exact match (default 20; at least 1).</summary>
    public int MinMatchLength { get; init; } = 20;

    /// <summary><c>-f</c>: the forward strand only, the reverse complement is not considered.</summary>
    public bool ForwardOnly { get; init; }

    /// <summary><c>-t</c>: only tandem repeats (matches whose two copies touch: <c>start1 + length &gt;= start2</c>); implies <see cref="ForwardOnly"/>.</summary>
    public bool TandemOnly { get; init; }

    /// <summary><c>-E</c>: find the matches by the exhaustive (quadratic) search of the program instead of the suffix tree.</summary>
    public bool Exhaustive { get; init; }

    /// <summary><c>-V #</c>: verbose level; above 0 the sizes of the program's structures are printed, above 1 the construction and the traversal of the suffix tree.</summary>
    public int Verbose { get; init; }
}

/// <summary>The MUMmer program stops with an exit status (<c>exit (EXIT_FAILURE)</c> after a message): raised inside the port and turned into <see cref="MummerRepeatMatchResult.ExitCode"/>.</summary>
/// <param name="code">The exit status of the program.</param>
public sealed class MummerExitException(int code) : Exception("The program exits with status " + code.ToString(CultureInfo.InvariantCulture))
{
    /// <summary>The exit status.</summary>
    public int Code { get; } = code;
}

/// <summary>Output of <see cref="MummerRepeatMatch.Run"/>.</summary>
/// <param name="Text">What the program writes to standard output.</param>
/// <param name="Messages">What the program writes to standard error.</param>
/// <param name="ExitCode">The exit status: 0, or 255 when the program stops (the genome is an exact palindrome, an inconsistent suffix tree).</param>
public sealed record MummerRepeatMatchResult(string Text, string Messages, int ExitCode);

/// <summary>
/// <c>repeat-match</c> of MUMmer 3.23 (<c>src/tigr/repeat-match.cc</c>, Delcher 2002): the maximal exact repeats (longer than <c>-n</c>) of the first sequence of a FASTA file, on the forward strand and, unless <c>-f</c>, between
/// the forward and the reverse-complement strand (<c>r</c> after the second position). The program builds a suffix tree of the sequence and of its reverse complement (McCreight style with suffix links; a new child is put in
/// front of its siblings) and lists the matches while it walks the tree; the order of the listing is therefore the order of that tree and is reproduced by the same construction.
/// </summary>
/// <remarks>
/// An empty file or a file without a <c>&gt;</c> line leaves the program's sequence buffer uninitialised (valgrind: use of uninitialised values; the output changes under <c>MALLOC_PERTURB_</c>): refused with
/// <see cref="FormatException"/>, as is a header whose first word has 500 characters or more (the program's name buffer overflows: <c>*** buffer overflow detected ***</c>, abort).
/// </remarks>
public static class MummerRepeatMatch
{
    /// <summary>Runs <c>repeat-match [options] genome-file</c>.</summary>
    /// <param name="genomeFasta">The FASTA text; only its first record is used.</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <exception cref="ArgumentException"><c>-n</c> is below 1 (the program prints <c>ERROR:  Illegal min match length</c> and exits).</exception>
    /// <exception cref="FormatException">The program fails or reads uninitialised memory (see the remarks).</exception>
    public static MummerRepeatMatchResult Run(string genomeFasta, MummerRepeatMatchOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(genomeFasta);
        var o = options ?? new MummerRepeatMatchOptions();
        if (o.MinMatchLength < 1)
            throw new ArgumentException("ERROR:  Illegal min match length \"" + o.MinMatchLength.ToString(CultureInfo.InvariantCulture) + "\"");
        var reader = new MummerFastaReader(genomeFasta);
        var warn = new StringBuilder();
        bool ok = reader.Read(out string name, out byte[] seq, out long len, warn);
        if (!ok)
            throw new FormatException("repeat-match: no sequence record: the program works on an uninitialised buffer (undefined; valgrind: use of uninitialised values).");
        if (name.Length >= 500)
            throw new FormatException("repeat-match: *** buffer overflow detected ***: the first word of the header has 500 characters or more (the 500 byte name buffer overflows; the program aborts).");
        var engine = new Tree(o, seq, (int)len);
        if (o.Verbose > 0)
            engine.Messages.Append("Verbose = ").Append(o.Verbose.ToString(CultureInfo.InvariantCulture)).Append('\n');
        engine.Messages.Append(warn);
        engine.Execute();
        return new MummerRepeatMatchResult(engine.Out.ToString(), engine.Messages.ToString(), engine.ExitCode);
    }

    private const int Nil = 0;
    private const byte Dollar = (byte)'$';
    private const byte StartChar = (byte)'%';

    private sealed class Tree
    {
        private readonly MummerRepeatMatchOptions _o;
        private readonly byte[] _data;
        private readonly int _genomeLen, _stringSeparator;
        private readonly bool _forwardOnly, _tandemOnly;
        private int _nextAvailNode = 1;

        // nodes
        private readonly int[] _nLo, _nChild, _nSibling, _nLink, _nLen, _nSubtreeSize;
        private readonly bool[] _nChildIsLeaf, _nSiblingIsLeaf, _nShouldSkip;

        // leaves
        private readonly int[] _lLo, _lSibling, _lLen;
        private readonly bool[] _lSiblingIsLeaf, _lIsDuplicate;
        private readonly int[] _nextLeaf;

        public StringBuilder Out { get; } = new();
        public StringBuilder Messages { get; } = new();
        public int ExitCode { get; private set; }

        public Tree(MummerRepeatMatchOptions o, byte[] seq, int genomeLen)
        {
            _o = o;
            _forwardOnly = o.ForwardOnly || o.TandemOnly;
            _tandemOnly = o.TandemOnly;
            _genomeLen = genomeLen;
            int dataLen = 3 + 2 * genomeLen;
            _stringSeparator = 1 + genomeLen;
            _data = new byte[dataLen + 2];
            _data[0] = StartChar;
            Array.Copy(seq, 1, _data, 1, genomeLen);
            Array.Copy(seq, 1, _data, genomeLen + 2, genomeLen);
            long lo = genomeLen + 2, hi = 2L * genomeLen + 1;
            while (lo <= hi)
            {
                byte ch = _data[hi];
                _data[hi] = MummerSequence.Complement(_data[lo]);
                _data[lo] = MummerSequence.Complement(ch);
                lo++;
                hi--;
            }

            _data[1 + genomeLen] = Dollar;
            _data[2 + 2 * genomeLen] = Dollar;
            _data[dataLen] = 0;
            int n = dataLen;
            _nLo = new int[n]; _nChild = new int[n]; _nSibling = new int[n]; _nLink = new int[n]; _nLen = new int[n]; _nSubtreeSize = new int[n];
            _nChildIsLeaf = new bool[n]; _nSiblingIsLeaf = new bool[n]; _nShouldSkip = new bool[n];
            _lLo = new int[n]; _lSibling = new int[n]; _lLen = new int[n];
            _lSiblingIsLeaf = new bool[n]; _lIsDuplicate = new bool[n];
            _nextLeaf = new int[n];
        }

        private void P(string s) => Out.Append(s);

        private static string D(long v) => v.ToString(CultureInfo.InvariantCulture);

        public void Execute()
        {
            try
            {
                Run();
            }
            catch (MummerExitException e)
            {
                ExitCode = e.Code;
            }
        }

        private void Run()
        {
            if (_o.Verbose > 0)
            {
                P("Node size = 24\nLeaf size = 12\nsizeof (int) = 4\nsizeof (size_t) = 8\n");
            }

            int root = BuildSuffixTree(1);
            if (!_forwardOnly && AddString(2 + _genomeLen, root) != 0)
                throw new MummerExitException(255);

            Messages.Append("Genome Length = ").Append(D(_genomeLen)).Append("   Used ").Append(D(_nextAvailNode)).Append(" internal nodes\n");
            SetSubtreeSize(root);
            MarkSkipableNodes(root);

            if (_o.Exhaustive)
            {
                _data[_genomeLen + 1] = (byte)'#';
                P("Exhaustive Exact Matches:\n");
                P(Header());
                for (int i = 1; i <= _genomeLen - _o.MinMatchLength; i++)
                {
                    for (int j = i + 1; j <= 1 + _genomeLen - _o.MinMatchLength; j++)
                    {
                        if (_data[i - 1] == _data[j - 1])
                            continue;
                        int match = LongestPrefixMatch(i, j);
                        if (match >= _o.MinMatchLength)
                            P(Line(i, j, ' ', match));
                    }
                }

                for (int i = 1; i <= _genomeLen + 1 - _o.MinMatchLength; i++)
                {
                    for (int j = _genomeLen + 2; j <= 2 * _genomeLen - i; j++)
                    {
                        if (_data[i - 1] == _data[j - 1])
                            continue;
                        int match = LongestPrefixMatch(i, j);
                        if (match >= _o.MinMatchLength)
                            P(Line(i, 2 * _genomeLen + 2 - j, 'r', match));
                    }
                }

                _data[_genomeLen + 1] = Dollar;
            }
            else
            {
                P("Long Exact Matches:\n");
                P(Header());
                ListMaximalMatches(root);
            }

            if (_o.Verbose > 1)
                Messages.Append("Global_Skip_Ct = ").Append(D(_skipCt)).Append("\nGlobal_Non_Skip_Ct = ").Append(D(_nonSkipCt)).Append('\n');
        }

        private int _skipCt, _nonSkipCt;

        private static string Header() => "Start1".PadLeft(9) + " " + "Start2".PadLeft(10) + "  " + "Length".PadLeft(8) + "\n";

        private static string Line(long a, long b, char c, long n) => D(a).PadLeft(9) + " " + D(b).PadLeft(10) + c + " " + D(n).PadLeft(8) + "\n";

        private int LongestPrefixMatch(int p, int q)
        {
            int i = 0;
            while (_data[p + i] == _data[q + i])
                i++;
            return i;
        }

        private int NewNode() => _nextAvailNode++;

        // ------------------------------------------------------------------ leaves and nodes: one place to link a new leaf in front of the children of a node

        private void PutLeaf(int place, int leaf, int lo, int len)
        {
            _lLo[leaf] = lo;
            _lSibling[leaf] = _nChild[place];
            _lSiblingIsLeaf[leaf] = _nChildIsLeaf[place];
            _nChild[place] = leaf;
            _nChildIsLeaf[place] = true;
            _lLen[leaf] = len;
        }

        // ------------------------------------------------------------------ construction

        private int BuildSuffixTree(int start)
        {
            int end = start;
            while (_data[end] != Dollar)
                end++;
            int root = NewNode();
            int leaf = start;
            _nLo[root] = 0;
            _nChild[root] = leaf;
            _nSibling[root] = Nil;
            _nLink[root] = root;
            _nLen[root] = 0;
            _nSiblingIsLeaf[root] = false;
            _nChildIsLeaf[root] = true;
            _lLo[leaf] = start;
            _lSibling[leaf] = Nil;
            _lLen[leaf] = 1 + end - start;
            _lSiblingIsLeaf[leaf] = false;

            int lastParent = root, lastParentDepth = 0, grandparent = 0, newPlace;
            while (++start <= end)
            {
                if (_o.Verbose > 1)
                    P("Build_Suffix_Tree:  Start = " + D(start) + "\n");
                int segmentStart, segmentLen, linkDepth, matched;
                if (_nLink[lastParent] != Nil)
                {
                    if (lastParent == root)
                    {
                        segmentStart = start;
                        segmentLen = 1 + end - start;
                        linkDepth = lastParentDepth;
                    }
                    else
                    {
                        segmentStart = start + lastParentDepth - 1;
                        segmentLen = 2 + end - start - lastParentDepth;
                        linkDepth = lastParentDepth - 1;
                    }

                    newPlace = StepDown(_nLink[lastParent], linkDepth, segmentStart, segmentLen, out matched, false, out _, ref grandparent);
                    if (matched >= segmentLen)
                    {
                        P("Ooops:  Suffix can't appear twice.\n");
                        throw new MummerExitException(255);
                    }

                    leaf = start;
                    PutLeaf(newPlace, leaf, segmentStart + matched, segmentLen - matched);
                    lastParent = newPlace;
                    lastParentDepth = linkDepth + matched;
                }
                else
                {
                    int leafLen = _lLen[leaf];
                    if (grandparent == root)
                    {
                        segmentStart = start;
                        segmentLen = _nLen[lastParent] - 1;
                        linkDepth = 0;
                    }
                    else
                    {
                        segmentStart = start + lastParentDepth - 1 - _nLen[lastParent];
                        segmentLen = _nLen[lastParent];
                        linkDepth = lastParentDepth - 1 - _nLen[lastParent];
                    }

                    newPlace = JumpDown(_nLink[grandparent], linkDepth, segmentStart, segmentLen, out bool madeNewNode, ref grandparent);
                    if (madeNewNode)
                    {
                        leaf = start;
                        PutLeaf(newPlace, leaf, segmentStart + segmentLen, leafLen);
                        _nLink[lastParent] = newPlace;
                        lastParent = newPlace;
                        lastParentDepth = linkDepth + segmentLen;
                    }
                    else
                    {
                        _nLink[lastParent] = newPlace;
                        segmentStart += segmentLen;
                        linkDepth += segmentLen;
                        segmentLen = leafLen;
                        newPlace = StepDown(newPlace, linkDepth, segmentStart, segmentLen, out matched, false, out _, ref grandparent);
                        if (matched >= segmentLen)
                        {
                            P("Ooops:  Suffix can't appear twice.\n");
                            throw new MummerExitException(255);
                        }

                        leaf = start;
                        PutLeaf(newPlace, leaf, segmentStart + matched, segmentLen - matched);
                        lastParent = newPlace;
                        lastParentDepth = linkDepth + matched;
                    }
                }
            }

            return root;
        }

        private int AddString(int start, int root)
        {
            int end = start;
            while (_data[end] != Dollar)
                end++;
            int grandparent = 0, matched;
            int segmentStart = start, segmentLen = 1 + end - start;
            int newPlace = StepDown(root, 0, segmentStart, segmentLen, out matched, true, out bool newPlaceIsLeaf, ref grandparent);
            if (segmentLen == matched)
            {
                P("*** Genome is exact palindrome ***\n");
                return -1;
            }

            if (1 + matched == segmentLen)
            {
                if (newPlaceIsLeaf)
                    Messages.Append("ERROR:  Unexpected leaf\n");
                P("String 1 is a substring of previous string.\n");
                return 0;
            }

            int leaf = start;
            PutLeaf(newPlace, leaf, segmentStart + matched, segmentLen - matched);
            int lastParent = newPlace, lastParentDepth = matched;

            while (++start <= end)
            {
                int linkDepth;
                if (_nLink[lastParent] != Nil)
                {
                    if (lastParent == root)
                    {
                        segmentStart = start;
                        segmentLen = 1 + end - start;
                        linkDepth = lastParentDepth;
                    }
                    else
                    {
                        segmentStart = start + lastParentDepth - 1;
                        segmentLen = 2 + end - start - lastParentDepth;
                        linkDepth = lastParentDepth - 1;
                    }

                    newPlace = StepDown(_nLink[lastParent], linkDepth, segmentStart, segmentLen, out matched, false, out _, ref grandparent);
                    if (matched == segmentLen)
                        return AddDuplicates(start, end, newPlace);
                    leaf = start;
                    PutLeaf(newPlace, leaf, segmentStart + matched, segmentLen - matched);
                    lastParent = newPlace;
                    lastParentDepth = linkDepth + matched;
                }
                else
                {
                    int leafLen = _lLen[leaf];
                    if (grandparent == root)
                    {
                        segmentStart = start;
                        segmentLen = _nLen[lastParent] - 1;
                        linkDepth = 0;
                    }
                    else
                    {
                        segmentStart = start + lastParentDepth - 1 - _nLen[lastParent];
                        segmentLen = _nLen[lastParent];
                        linkDepth = lastParentDepth - 1 - _nLen[lastParent];
                    }

                    newPlace = JumpDown(_nLink[grandparent], linkDepth, segmentStart, segmentLen, out bool madeNewNode, ref grandparent);
                    if (madeNewNode)
                    {
                        leaf = start;
                        PutLeaf(newPlace, leaf, segmentStart + segmentLen, leafLen);
                        _nLink[lastParent] = newPlace;
                        lastParent = newPlace;
                        lastParentDepth = linkDepth + segmentLen;
                    }
                    else
                    {
                        _nLink[lastParent] = newPlace;
                        segmentStart += segmentLen;
                        linkDepth += segmentLen;
                        segmentLen = leafLen;
                        newPlace = StepDown(newPlace, linkDepth, segmentStart, segmentLen, out matched, false, out _, ref grandparent);
                        if (matched >= segmentLen)
                            return AddDuplicates(start, end, newPlace);
                        leaf = start;
                        PutLeaf(newPlace, leaf, segmentStart + matched, segmentLen - matched);
                        lastParent = newPlace;
                        lastParentDepth = linkDepth + matched;
                    }
                }
            }

            return 0;
        }

        private int AddDuplicates(int start, int end, int leaf)
        {
            int j = leaf;
            for (int i = start; i <= end; i++, j++)
                _lIsDuplicate[j] = true;
            return 0;
        }

        /// <summary><c>New_Find_Child</c>: the child of <paramref name="node"/> whose string starts with <paramref name="ch"/> (0 = none), with its kind and its predecessor among the siblings.</summary>
        private int FindChild(int node, byte ch, out bool isLeaf, out int pred, out bool predIsLeaf)
        {
            pred = Nil;
            predIsLeaf = false;
            int i = _nChild[node];
            isLeaf = _nChildIsLeaf[node];
            byte startCh = StartCh(i, isLeaf);
            while (i != Nil && startCh != ch)
            {
                pred = i;
                predIsLeaf = isLeaf;
                if (isLeaf)
                {
                    isLeaf = _lSiblingIsLeaf[i];
                    i = _lSibling[i];
                }
                else
                {
                    isLeaf = _nSiblingIsLeaf[i];
                    i = _nSibling[i];
                }

                startCh = StartCh(i, isLeaf);
            }

            return i;
        }

        private byte StartCh(int i, bool isLeaf)
        {
            int lo = isLeaf ? _lLo[i] : _nLo[i];
            int len = isLeaf ? _lLen[i] : _nLen[i];
            return len > 0 ? _data[lo] : MummerSequence.Complement(_data[lo]);
        }

        /// <summary>Splits the edge to <paramref name="p"/> after <paramref name="len"/> characters (the common code of <c>New_Jump_Down</c> and <c>New_Step_Down</c>): the new node becomes the child of <paramref name="node"/> where <paramref name="p"/> was.</summary>
        private int Split(int node, int p, bool pIsLeaf, int pSib, bool pSibIsLeaf, int pred, bool predIsLeaf, int lo, int len)
        {
            int q = NewNode();
            _nLo[q] = lo;
            _nLen[q] = len;
            _nChild[q] = p;
            _nSibling[q] = pSib;
            _nSiblingIsLeaf[q] = pSibIsLeaf;
            _nLink[q] = Nil;
            _nChildIsLeaf[q] = pIsLeaf;
            if (pred == Nil)
            {
                _nChild[node] = q;
                _nChildIsLeaf[node] = false;
            }
            else if (predIsLeaf)
            {
                _lSibling[pred] = q;
                _lSiblingIsLeaf[pred] = false;
            }
            else
            {
                _nSibling[pred] = q;
                _nSiblingIsLeaf[pred] = false;
            }

            if (!pIsLeaf)
            {
                _nLo[p] += len;
                _nLen[p] -= len;
                _nSibling[p] = Nil;
                _nSiblingIsLeaf[p] = false;
            }
            else
            {
                _lLo[p] += len;
                _lLen[p] -= len;
                _lSibling[p] = Nil;
                _lSiblingIsLeaf[p] = false;
            }

            return q;
        }

        private int JumpDown(int node, int nodeDepth, int lo, int len, out bool madeNewNode, ref int parNewNode)
        {
            madeNewNode = false;
            if (len == 0)
                return node;
            if (_o.Verbose > 1)
                P("Jump_Down:  Node = " + D(node) + "  Depth = " + D(nodeDepth) + "  Lo = " + D(lo) + "  Len = " + D(len) + "\n");
            int p = FindChild(node, _data[lo], out bool pIsLeaf, out int pred, out bool predIsLeaf);
            while (p != Nil)
            {
                int d, pSib;
                bool pSibIsLeaf;
                if (pIsLeaf)
                {
                    d = _lLen[p];
                    pSib = _lSibling[p];
                    pSibIsLeaf = _lSiblingIsLeaf[p];
                }
                else
                {
                    d = _nLen[p];
                    pSib = _nSibling[p];
                    pSibIsLeaf = _nSiblingIsLeaf[p];
                }

                if (len == d)
                    return p;
                if (d < len)
                {
                    lo += d;
                    len -= d;
                    node = p;
                    nodeDepth += d;
                    p = FindChild(node, _data[lo], out pIsLeaf, out pred, out predIsLeaf);
                }
                else
                {
                    int q = Split(node, p, pIsLeaf, pSib, pSibIsLeaf, pred, predIsLeaf, lo, len);
                    parNewNode = node;
                    madeNewNode = true;
                    return q;
                }
            }

            P("Ooops:  Couldn't find appropriate child node\n");
            throw new MummerExitException(255);
        }

        private int StepDown(int node, int nodeDepth, int lo, int len, out int depth, bool doingPrefix, out bool returnIsLeaf, ref int parNewNode)
        {
            depth = 0;
            returnIsLeaf = false;
            if (len == 0)
                return Nil;
            int p = FindChild(node, _data[lo], out bool pIsLeaf, out int pred, out bool predIsLeaf);
            while (p != Nil)
            {
                int i, d, pSib;
                bool pSibIsLeaf;
                if (pIsLeaf)
                {
                    i = _lLo[p];
                    d = _lLen[p];
                    pSib = _lSibling[p];
                    pSibIsLeaf = _lSiblingIsLeaf[p];
                }
                else
                {
                    i = _nLo[p];
                    d = _nLen[p];
                    pSib = _nSibling[p];
                    pSibIsLeaf = _nSiblingIsLeaf[p];
                }

                int j = 1;
                while (j < len && j < d && _data[i + j] == _data[lo + j])
                    j++;
                depth += j;

                if (j < d)
                {
                    if (doingPrefix && 1 + j == len)
                    {
                        returnIsLeaf = pIsLeaf;
                        return p;
                    }

                    int q = Split(node, p, pIsLeaf, pSib, pSibIsLeaf, pred, predIsLeaf, lo, j);
                    parNewNode = node;
                    returnIsLeaf = false;
                    return q;
                }

                if (pIsLeaf || j == len)
                {
                    returnIsLeaf = pIsLeaf;
                    return p;
                }

                lo += j;
                len -= j;
                node = p;
                nodeDepth += j;
                p = FindChild(node, _data[lo], out pIsLeaf, out pred, out predIsLeaf);
            }

            returnIsLeaf = false;
            return node;
        }

        // ------------------------------------------------------------------ the walks over the tree (recursion of the program, with an explicit stack: the tree of a repeat of 10^5 bases is that deep)

        private struct Frame(int root, bool isLeaf, int parent, int parentDepth)
        {
            public int Root = root, Parent = parent, ParentDepth = parentDepth, State;
            public bool IsLeaf = isLeaf;
        }

        /// <summary>
        /// <c>Set_Subtree_Size</c>, <c>Mark_Skipable_Nodes</c> and <c>List_Maximal_Matches</c> all recurse the same way: a leaf (action <paramref name="leafAction"/>) continues with its sibling; an internal node first visits its child
        /// subtree, then its sibling, and only then does its own action <paramref name="nodeAction"/> (<c>root, parent, parentDepth, depth</c>).
        /// </summary>
        private void Walk(int root, Action<int, int, int>? leafAction, Action<int, int, int, int> nodeAction)
        {
            var stack = new Stack<Frame>();
            stack.Push(new Frame(root, false, Nil, 0));
            while (stack.Count > 0)
            {
                var f = stack.Pop();
                if (f.Root == Nil)
                    continue;
                if (f.IsLeaf)
                {
                    leafAction?.Invoke(f.Root, f.Parent, f.ParentDepth);
                    stack.Push(new Frame(_lSibling[f.Root], _lSiblingIsLeaf[f.Root], f.Parent, f.ParentDepth));
                    continue;
                }

                int depth = f.ParentDepth + _nLen[f.Root];
                switch (f.State)
                {
                    case 0:
                        f.State = 1;
                        stack.Push(f);
                        stack.Push(new Frame(_nChild[f.Root], _nChildIsLeaf[f.Root], f.Root, depth));
                        break;
                    case 1:
                        f.State = 2;
                        stack.Push(f);
                        stack.Push(new Frame(_nSibling[f.Root], _nSiblingIsLeaf[f.Root], f.Parent, f.ParentDepth));
                        break;
                    default:
                        nodeAction(f.Root, f.Parent, f.ParentDepth, depth);
                        break;
                }
            }
        }

        private void SetSubtreeSize(int root) => Walk(root,
            (leaf, parent, _) => _nSubtreeSize[parent] += _lIsDuplicate[leaf] ? 1 : 2,
            (node, parent, parentDepth, _) =>
            {
                if (_o.Verbose > 1)
                    P("Set_Subtree_Size:  Root = " + D(node) + "  Parent = " + D(parent) + "  Parent_Depth = " + D(parentDepth) + "\n");
                if (parent != Nil)
                    _nSubtreeSize[parent] += _nSubtreeSize[node];
            });

        private void MarkSkipableNodes(int root) => Walk(root, null, (node, parent, parentDepth, _) =>
        {
            if (_o.Verbose > 1)
                P("Mark_Skipable_Nodes:  Root = " + D(node) + "  Parent = " + D(parent) + "  Parent_Depth = " + D(parentDepth) + "\n");
            int link = _nLink[node];
            if (link != Nil && _nSubtreeSize[link] == _nSubtreeSize[node])
                _nShouldSkip[link] = true;
        });

        private void ListMaximalMatches(int root) => Walk(root, null, (node, parent, parentDepth, depth) =>
        {
            if (_o.Verbose > 1)
            {
                P("List_Maximal_Matches:  Root = " + D(node) + "  Parent = " + D(parent) + "  Parent_Depth = " + D(parentDepth) + "\n   Depth = " + D(depth) + "  Subtree_Size = " + D(_nSubtreeSize[node]) + "\n");
            }

            if (depth >= _o.MinMatchLength)
            {
                if (_nShouldSkip[node])
                    _skipCt++;
                else
                    _nonSkipCt++;
            }

            if (depth >= _o.MinMatchLength && !_nShouldSkip[node])
            {
                bool isLeaf = _nChildIsLeaf[node];
                for (int j = _nChild[node]; j != Nil;)
                {
                    if (_o.Verbose > 1)
                        P("  Child = " + D(j) + " " + (isLeaf ? "leaf" : "node") + "\n");
                    int list1, k;
                    bool kIsLeaf;
                    if (isLeaf)
                    {
                        list1 = j;
                        k = _lSibling[j];
                        kIsLeaf = _lSiblingIsLeaf[j];
                        isLeaf = _lSiblingIsLeaf[j];
                        j = _lSibling[j];
                    }
                    else
                    {
                        list1 = _nLink[j];
                        k = _nSibling[j];
                        kIsLeaf = _nSiblingIsLeaf[j];
                        isLeaf = _nSiblingIsLeaf[j];
                        j = _nSibling[j];
                    }

                    while (k != Nil)
                    {
                        if (kIsLeaf)
                        {
                            ListMatches(list1, k, depth);
                            kIsLeaf = _lSiblingIsLeaf[k];
                            k = _lSibling[k];
                        }
                        else
                        {
                            ListMatches(list1, _nLink[k], depth);
                            kIsLeaf = _nSiblingIsLeaf[k];
                            k = _nSibling[k];
                        }
                    }
                }
            }

            int start = 0, end = 0;
            bool leafKind = _nChildIsLeaf[node];
            for (int j = _nChild[node]; j != Nil;)
            {
                if (leafKind)
                {
                    if (start == 0)
                    {
                        start = end = j;
                    }
                    else
                    {
                        _nextLeaf[end] = j;
                        end = j;
                    }

                    leafKind = _lSiblingIsLeaf[j];
                    j = _lSibling[j];
                }
                else
                {
                    if (start == 0)
                        start = end = _nLink[j];
                    else
                        _nextLeaf[end] = _nLink[j];
                    while (_nextLeaf[end] != Nil)
                        end = _nextLeaf[end];
                    leafKind = _nSiblingIsLeaf[j];
                    j = _nSibling[j];
                }
            }

            _nLink[node] = start;
        });

        private void ListMatches(int a, int b, int n)
        {
            for (int i = a; i != Nil; i = _nextLeaf[i])
            {
                for (int j = b; j != Nil; j = _nextLeaf[j])
                {
                    if (_data[i - 1] == _data[j - 1] || _data[i + n] == _data[j + n] || (i > _stringSeparator && j > _stringSeparator))
                        continue;
                    bool reversed = false;
                    int l, r;
                    if (j > _stringSeparator)
                    {
                        int k = _genomeLen - (j - _stringSeparator) - n + 2;
                        if (k < i)
                            continue;
                        l = i;
                        r = k + n - 1;
                        reversed = true;
                    }
                    else if (i > _stringSeparator)
                    {
                        int k = _genomeLen - (i - _stringSeparator) - n + 2;
                        if (k < j)
                            continue;
                        l = j;
                        r = k + n - 1;
                        reversed = true;
                    }
                    else if (i < j)
                    {
                        l = i;
                        r = j;
                    }
                    else
                    {
                        l = j;
                        r = i;
                    }

                    if (_o.Verbose > 1)
                        VerifyMatch(l, r, n, reversed);
                    if (_tandemOnly && l + n < r)
                        continue;
                    P(Line(l, r, reversed ? 'r' : ' ', n));
                }
            }
        }

        private void VerifyMatch(int a, int b, int n, bool reverse)
        {
            string rev = reverse ? "t" : "f";
            string Pre(string what) => "Verify_Match:  a=" + D(a) + "  b=" + D(b) + "  rev=" + rev + "  " + what + " chars match\n";
            if (reverse)
            {
                if (_data[a - 1] == MummerSequence.Complement(_data[b + 1]))
                    P(Pre("preceding"));
                if (_data[a + n] == MummerSequence.Complement(_data[b - n]))
                    P(Pre("following"));
                for (int i = 0; i < n; i++)
                {
                    if (_data[a + i] != MummerSequence.Complement(_data[b - i]))
                    {
                        P("Verify_Match:  a=" + D(a) + "  b=" + D(b) + "  rev=t  mismatch at " + D(a + i) + "[" + (char)_data[a + i] + "]:" + D(b - i) + "[" + (char)MummerSequence.Complement(_data[b - i]) + "]\n");
                        throw new MummerExitException(255);
                    }
                }
            }
            else
            {
                if (_data[a - 1] == _data[b - 1])
                    P(Pre("preceding"));
                if (_data[a + n] == _data[b + n])
                    P(Pre("following"));
                for (int i = 0; i < n; i++)
                {
                    if (_data[a + i] != _data[b + i])
                    {
                        P("Verify_Match:  a=" + D(a) + "  b=" + D(b) + "  rev=f  mismatch at " + D(a + i) + "[" + (char)_data[a + i] + "]:" + D(b + i) + "[" + (char)_data[b + i] + "]\n");
                        throw new MummerExitException(255);
                    }
                }
            }
        }
    }
}
