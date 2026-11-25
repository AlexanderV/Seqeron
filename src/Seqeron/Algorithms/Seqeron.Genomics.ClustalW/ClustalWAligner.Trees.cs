// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

public static partial class ClustalWAligner
{
    // ClustalW 2.1 guide-tree files (audit A2-14b, B09 F104): the writer of the NJ guide tree (`-NEWTREE`, ClusterTreeOutput::printPhylipTree / twoWaySplit; the two-sequence form of
    // UnRootedClusterTree::treeFromDistMatrix) and the reader of a user guide tree (`-USETREE`, Tree::readTree / createTree / insertNode / reRoot / createSets / groupSeqs /
    // calcSeqWeights / calcSimilarities of the clustalw 2.1+lgpl-7 Debian source: src/tree/{Tree,ClusterTreeOutput,UnRootedClusterTree,TreeInterface}.cpp).

    private const int CwTreeMaxNames = 150;   // MAXNAMES of clustalw.h

    /// <summary>The result of reading a user guide tree: the rooted tree (<c>Tree::root</c>), the leaf node of every sequence of the group (<c>olptr</c>, group-relative), all parsed leaves (<c>lptr</c>) and <c>distanceTree</c>.</summary>
    private sealed record CwUserTree(CwNode Root, CwNode[] Leaves, IReadOnlyList<CwNode> ParsedLeaves, bool DistanceTree);

    /// <summary>
    /// <c>Tree::readTree</c>: parses the Newick text (<c>createTree</c>: a hand-written reader with the quirks of the program), checks it against the group of sequences first .. first + count - 1 of
    /// <paramref name="names"/>, re-roots an unrooted tree and assigns the sequences to the leaves by name.
    /// </summary>
    private static CwUserTree ClustalWReadTree(string text, IReadOnlyList<string> names, int first, int count)
    {
        var reader = new CwTreeCharStream(text);
        reader.SkipSpace();
        int c = reader.Get();
        if (c != '(')
            throw new ClustalWTreeFileException(ClustalWTreeFileProblem.WrongFormat, "Wrong format in tree file: the program requires the first character of the tree to be '('.");

        reader = new CwTreeCharStream(text);   // file.seekg(0, ios::beg)
        var ctx = new CwTreeParse(reader, count);
        var seqTree = new CwNode();
        ctx.DistanceTree = true;
        ctx.CreateTree(seqTree, null);

        if (ctx.NumSeq != count)
            throw new ClustalWTreeFileException(ClustalWTreeFileProblem.NotCompatible, $"tree not compatible with alignment ({count} sequences in alignment and {ctx.NumSeq} in tree).");
        if (ctx.Ptrs.Count > 3 * (count + 1))
            throw new ClustalWTreeFileException(ClustalWTreeFileProblem.Undefined, "tree has more nodes than the program's arrays can hold (the program overruns its heap buffer and crashes).");

        if (!ctx.DistanceTree && !ctx.RootedTree)
            throw new ClustalWTreeFileException(ClustalWTreeFileProblem.UnrootedWithoutDistances, "input tree is unrooted and has no distances. Cannot align sequences.");

        CwNode root = ctx.RootedTree ? seqTree : ReRoot(seqTree, ctx.Ptrs, ctx.Lptr, count + 1);

        // Tree::clearTree deletes the nodes below the root: a node reachable twice (the trees of a nested trichotomy are no trees) is freed twice and the program aborts
        var visited = new HashSet<CwNode>();
        bool shared = false;
        void Visit(CwNode? n)
        {
            if (n == null || shared)
                return;
            if (!visited.Add(n))
            {
                shared = true;
                return;
            }

            Visit(n.Left);
            Visit(n.Right);
        }

        Visit(root);
        if (shared)
            throw new ClustalWTreeFileException(ClustalWTreeFileProblem.Undefined, "the guide tree shares a node between two branches (nested trichotomy): the program frees it twice and aborts or loops.");

        var olptr = new CwNode[count];
        for (int i = 0; i < count; i++)
        {
            string name = names[first + i];
            if (name.Length > CwTreeMaxNames)
                name = name[..CwTreeMaxNames];
            var name2 = new StringBuilder();
            foreach (char ch in name)
            {
                char lowered = ch is > '@' and < '[' ? (char)(ch | 0x20) : ch;
                name2.Append(ch == ' ' ? '_' : lowered);
            }

            CwNode? found = null;
            foreach (var leaf in ctx.Lptr)
            {
                string raw = leaf.Name.Length > CwTreeMaxNames ? leaf.Name[..CwTreeMaxNames] : leaf.Name;
                var name1 = new StringBuilder();
                foreach (char ch in raw)
                    name1.Append(ch is > '@' and < '[' ? (char)(ch | 0x20) : ch);
                if (name1.ToString() == name2.ToString())
                    found = leaf;
            }

            if (found == null)
                throw new ClustalWTreeFileException(ClustalWTreeFileProblem.NameNotFound, $"tree not compatible with alignment: {name2} not found.");
            olptr[i] = found;
        }

        return new CwUserTree(root, olptr, ctx.Lptr, ctx.DistanceTree);
    }

    /// <summary>The state of <c>Tree::createTree</c> (the nodes in the order the program stores them, the counters and flags).</summary>
    private sealed class CwTreeParse
    {
        private readonly CwTreeCharStream _file;
        private readonly int _capacity;
        private int _c;   // Tree::charFromFile

        public CwTreeParse(CwTreeCharStream file, int count)
        {
            _file = file;
            _capacity = count;
        }

        public List<CwNode> Ptrs { get; } = [];
        public List<CwNode> Lptr { get; } = [];
        public int NumSeq => Lptr.Count;
        public bool RootedTree { get; private set; } = true;
        public bool DistanceTree { get; set; } = true;

        private static void CreateNode(CwNode pptr, CwNode? parent)
        {
            pptr.Parent = parent;
            pptr.Left = new CwNode();
            pptr.Right = new CwNode();
        }

        private static CwNode InsertNode(CwNode pptr)
        {
            var newnode = new CwNode();
            CreateNode(newnode, pptr.Parent);
            newnode.Left = pptr;
            pptr.Parent = newnode;
            // setInfo(newnode, pptr->parent, NODE, "", 0.0): the parent was just set to the new node's own parent - which is overwritten by pptr->parent (now the new node itself)
            newnode.Parent = pptr.Parent;
            newnode.Dist = 0f;
            newnode.Order = 0;
            newnode.Leaf = false;
            return newnode;
        }

        private void Register(CwNode node)
        {
            // the program's arrays hold 3 * (n + 1) nodes and n + 1 leaves; a larger tree overruns them (heap corruption: the outcome depends on the heap, the program mostly prints its
            // "tree not compatible" error and sometimes crashes). The reader goes on far enough to count the leaves, but not into absurd sizes (recursion depth).
            if (Ptrs.Count >= 3 * (_capacity + 1) + 10000)
                throw new ClustalWTreeFileException(ClustalWTreeFileProblem.Undefined, "tree has far more nodes than the program's arrays can hold (the program overruns its heap buffer and crashes).");
            Ptrs.Add(node);
        }

        public void CreateTree(CwNode ptree, CwNode? parent)
        {
            bool isNode;
            string name = "";
            _file.SkipSpace();
            _c = _file.Get();
            if (_c == '(')
            {
                isNode = true;
                Register(ptree);
                CreateNode(ptree, parent);

                var p = ptree.Left!;
                CreateTree(p, ptree);

                if (_c == ',')
                {
                    p = ptree.Right!;
                    CreateTree(p, ptree);
                    if (_c == ',')
                    {
                        ptree = InsertNode(ptree);
                        Register(ptree);
                        p = ptree.Right!;
                        CreateTree(p, ptree);
                        RootedTree = false;
                    }
                }

                _file.SkipSpace();
                _c = _file.Get();
            }
            else
            {
                isNode = false;
                Register(ptree);
                Lptr.Add(ptree);
                var sb = new StringBuilder();
                sb.Append((char)_c);
                _c = _file.Get();
                int i = 1;
                while (_c != ':' && _c != ',' && _c != ')')
                {
                    if (_c == CwTreeCharStream.Eof && _file.Failed)
                        throw new ClustalWTreeFileException(ClustalWTreeFileProblem.Undefined, "the tree ends inside a sequence name: the program reads past the end of the file and never terminates.");
                    if (i < CwTreeMaxNames)
                    {
                        sb.Append((char)_c);
                        i++;
                    }

                    _c = _file.Get();
                }

                name = sb.ToString();
                if (_c != ':')
                    DistanceTree = false;
            }

            float dist = 0f;
            if (_c == ':')
            {
                _file.SkipSpace();
                dist = _file.ReadFloat();
                _file.SkipSpace();
                _c = _file.Get();
            }

            // setInfo(ptree, parent, type, name, dist)
            ptree.Parent = parent;
            ptree.Leaf = !isNode;
            ptree.Dist = dist;
            ptree.Order = 0;
            ptree.Name = name;
            if (ptree.Leaf)
            {
                ptree.Left = null;
                ptree.Right = null;
            }
        }
    }

    /// <summary>
    /// The <c>std::ifstream</c> the tree reader of the program uses, as far as it matters: <c>get</c> (EOF after the end sets failbit and every later read fails), <c>putback</c> (a no-op in the fail state),
    /// and <c>operator&gt;&gt;(float)</c> of libstdc++'s <c>num_get</c> (the characters it accepts, <c>strtof</c> of them, failbit on a failed conversion).
    /// </summary>
    private sealed class CwTreeCharStream
    {
        public const int Eof = -1;
        private readonly string _text;
        private int _pos;

        public CwTreeCharStream(string text) => _text = text;

        /// <summary>failbit (and eofbit) set: every later read fails.</summary>
        public bool Failed { get; private set; }

        public int Get()
        {
            if (Failed)
                return Eof;
            if (_pos >= _text.Length)
            {
                Failed = true;
                return Eof;
            }

            return _text[_pos++];
        }

        // Tree::skipSpace: read characters until one is not a space, then put it back
        public void SkipSpace()
        {
            int c;
            do
            {
                c = Get();
            }
            while (c != Eof && IsSpace(c));

            if (!Failed && c != Eof)
                _pos--;
        }

        private static bool IsSpace(int c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';

        /// <summary><c>file &gt;&gt; dist</c> for a <c>float</c> (libstdc++ <c>num_get::_M_extract_float</c> + <c>__convert_to_v</c>); returns 0 and sets failbit when the conversion fails.</summary>
        public float ReadFloat()
        {
            if (Failed)
                return 0f;
            // the sentry skips white space (the reader has already done so)
            while (_pos < _text.Length && IsSpace(_text[_pos]))
                _pos++;

            var x = new StringBuilder();
            int p = _pos;
            bool EofAt() => p >= _text.Length;
            if (!EofAt() && (_text[p] == '+' || _text[p] == '-'))
            {
                if (_text[p] == '-')
                    x.Append('-');
                p++;
            }

            bool foundMantissa = false;
            while (!EofAt())
            {
                char ch = _text[p];
                if (ch == '.')
                    break;
                if (ch == '0')
                {
                    if (!foundMantissa)
                    {
                        x.Append('0');
                        foundMantissa = true;
                    }

                    p++;
                }
                else
                {
                    break;
                }
            }

            bool foundDec = false, foundSci = false;
            while (!EofAt())
            {
                char ch = _text[p];
                if (ch is >= '0' and <= '9')
                {
                    x.Append(ch);
                    foundMantissa = true;
                }
                else if (ch == '.' && !foundDec && !foundSci)
                {
                    x.Append('.');
                    foundDec = true;
                }
                else if ((ch == 'e' || ch == 'E') && !foundSci && foundMantissa)
                {
                    x.Append('e');
                    foundSci = true;
                    p++;
                    if (!EofAt())
                    {
                        char s = _text[p];
                        if (s == '+' || s == '-')
                            x.Append(s);
                        else
                            continue;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }

                p++;
            }

            bool atEnd = EofAt();
            _pos = p;
            string s2 = x.ToString();
            float value = 0f;
            bool ok = s2.Length > 0 && float.TryParse(s2, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
            if (ok)
            {
                if (float.IsInfinity(value))
                {
                    value = value > 0 ? float.MaxValue : float.MinValue;
                    ok = false;
                }
            }
            else
            {
                value = 0f;
            }

            if (atEnd)
                Failed = true;   // eofbit while reading: the next read fails; the extraction itself still succeeds when the conversion did
            if (!ok)
                Failed = true;
            return value;
        }
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Tree::createSets / groupSeqs
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// <summary><c>Tree::createSets(0, n)</c>: one alignment step per internal node (post-order, left subtree = 1, right subtree = 2), 1-based vectors of n + 1 entries (<c>AlignmentSteps::saveSet</c>).</summary>
    private static List<int[]> ClustalWCreateSets(CwNode root, CwNode[] olptr)
    {
        int n = olptr.Length;
        var steps = new List<int[]>();
        void MarkGroup1(CwNode p, int[] groups)
        {
            for (int i = 0; i < n; i++)
                groups[i] = olptr[i] == p ? 1 : 0;
        }

        void MarkGroup2(CwNode p, int[] groups)
        {
            for (int i = 0; i < n; i++)
            {
                if (olptr[i] == p)
                    groups[i] = 2;
                else if (groups[i] != 0)
                    groups[i] = 1;
            }
        }

        int depth = 0;
        void GroupSeqs(CwNode p, int[] nextGroups)
        {
            // a node below itself (nested trichotomies leave such tree): the program recurses until its stack is exhausted
            if (++depth > 10000)
                throw new ClustalWTreeFileException(ClustalWTreeFileProblem.Undefined, "the guide tree contains a cycle: the program recurses until it crashes.");
            var tmp = new int[n + 1];
            if (p.Left != null)
            {
                if (!p.Left.Leaf)
                {
                    GroupSeqs(p.Left, nextGroups);
                    for (int i = 0; i < n; i++)
                    {
                        if (nextGroups[i] != 0)
                            tmp[i] = 1;
                    }
                }
                else
                {
                    MarkGroup1(p.Left, tmp);
                }
            }

            if (p.Right != null)
            {
                if (!p.Right.Leaf)
                {
                    GroupSeqs(p.Right, nextGroups);
                    for (int i = 0; i < n; i++)
                    {
                        if (nextGroups[i] != 0)
                            tmp[i] = 2;
                    }
                }
                else
                {
                    MarkGroup2(p.Right, tmp);
                }

                var saved = new int[n + 1];
                for (int i = 1; i <= n; i++)
                    saved[i] = tmp[i - 1];
                steps.Add(saved);
            }

            for (int i = 0; i < n; i++)
                nextGroups[i] = tmp[i];
            depth--;
        }

        GroupSeqs(root, new int[n + 1]);
        return steps;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // The guide-tree file of the NJ clustering (ClusterTreeOutput::printPhylipTree)
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// <summary><c>ostream &lt;&lt; fixed &lt;&lt; setprecision(digits) &lt;&lt; value</c> (printf <c>%.*f</c>: the exact binary value rounded half to even; see <see cref="CFormat.Fixed"/>).</summary>
    internal static string ClustalWFixed(double value, int digits) => CFormat.Fixed(value, digits);

    /// <summary>
    /// The one-pair guide tree of <c>UnRootedClusterTree::treeFromDistMatrix</c>: <c>float dist = distMat(first, first + 1) / 2.0</c> written with <c>setprecision(5)</c> (<c>%.5g</c> of the float) and the
    /// value the program reads back from the file as a <c>float</c>.
    /// </summary>
    private static (string Text, float ReadBack) ClustalWPairTreeDistance(double distance)
    {
        string text = CppGeneral((double)(float)(distance / 2.0), 5);
        return (text, float.Parse(text, CultureInfo.InvariantCulture));
    }

    /// <summary>
    /// The guide-tree file <c>clustalw</c> writes for the neighbour-joining clustering (<c>UnRootedClusterTree::treeFromDistMatrix</c> + <c>ClusterTreeOutput::printPhylipTree</c> / <c>twoWaySplit</c>) of the sequences
    /// <paramref name="names"/> from their distances: <c>(a:d,b:d);</c> for two sequences (d = half the distance, five significant digits), else the unrooted tree as three subtrees under the first parenthesis,
    /// every node on its own lines, branch lengths with five decimals.
    /// </summary>
    /// <param name="names">Sequence names (file order).</param>
    /// <param name="distances">ClustalW's pairwise distances (0-based square matrix).</param>
    public static string ClustalWGuideTreeText(IReadOnlyList<string> names, double[,] distances)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(distances);
        int k = names.Count;
        if (k < 2)
            throw new ArgumentException("A guide tree needs at least two sequences.", nameof(names));
        if (distances.GetLength(0) != k || distances.GetLength(1) != k)
            throw new ArgumentException("The distance matrix must be square with one row per name.", nameof(distances));
        return ClustalWNjTreeText(names, k == 2 ? null : ClustalWNeighborJoining(k, distances), k == 2 ? distances[0, 1] : 0.0);
    }

    private static string ClustalWNjTreeText(IReadOnlyList<string> names, PhylogeneticAnalyzer.PhyloNode? unrooted, double pairDistance)
    {
        if (unrooted == null)
        {
            string d = ClustalWPairTreeDistance(pairDistance).Text;
            return $"({names[0]}:{d},{names[1]}:{d});\n";
        }

        var sb = new StringBuilder();
        sb.Append("(\n");
        var top = unrooted.Children.OrderBy(MinLeafIndex).ToList();
        for (int i = 0; i < 3; i++)
        {
            var child = top[i];
            WriteSubtree(sb, child, names, top: true);
            sb.Append(':').Append(ClustalWFixed(child.BranchLength, 5));
            sb.Append(i < 2 ? ",\n" : ")");
        }

        sb.Append(";\n");
        return sb.ToString();
    }

    // twoWaySplit: a leaf prints its name; a node "(\n" first-child [":" len] ",\n" second-child [":" len] ")\n" (children by lowest leaf index); the caller appends the branch length
    private static void WriteSubtree(StringBuilder sb, PhylogeneticAnalyzer.PhyloNode node, IReadOnlyList<string> names, bool top)
    {
        if (node.IsLeaf)
        {
            sb.Append(names[int.Parse(node.Name, CultureInfo.InvariantCulture)]);
            return;
        }

        var kids = node.Children.OrderBy(MinLeafIndex).ToList();
        sb.Append("(\n");
        WriteSubtree(sb, kids[0], names, false);
        sb.Append(':').Append(ClustalWFixed(kids[0].BranchLength, 5)).Append(",\n");
        WriteSubtree(sb, kids[1], names, false);
        sb.Append(':').Append(ClustalWFixed(kids[1].BranchLength, 5)).Append(")\n");
    }
}
