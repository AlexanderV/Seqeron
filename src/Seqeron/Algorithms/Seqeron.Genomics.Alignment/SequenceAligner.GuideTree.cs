using System.Globalization;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region Guide trees of the progressive family (UPGMA / neighbour-joining)

    // The guide tree of the progressive, iterative and consistency aligners is built by the canonical tree builders of
    // Seqeron.Genomics.Phylogenetics (PhylogeneticAnalyzer.BuildTreeFromMatrix: UPGMA, Sokal & Michener 1958;
    // neighbour-joining, Saitou & Nei 1987) — this file only adapts their result to the binary, rooted tree that a
    // progressive alignment needs. (Until A1-12 this class carried a private copy of UPGMA, justified by a "circular"
    // project reference Phylogenetics -> Alignment that no code used; the reference is gone, the copy is gone.)
    //
    //   * MsaGuideTree.Upgma  — rooted and binary by construction; ties (distances equal within ProfileScoreEpsilon) are
    //     resolved to the lowest (i, j) index pair, exactly the rule of the former private copy.
    //   * MsaGuideTree.NeighborJoining — ClustalW's progressive guide tree (Thompson, Higgins & Gibson 1994; ClustalW 2.1
    //     `-CLUSTERING=NJ`, the default). NJ is an UNROOTED tree whose centre is a trifurcation, so ClustalW (a) prints the
    //     unrooted tree (`.dnd`: branch lengths with 5 decimals, |b| < 1e-4 set to 0, children ordered by their lowest
    //     leaf index), (b) re-reads it with single-precision branch lengths and (c) roots it with Tree::reRoot (the node
    //     whose root position gives the smallest mean root-to-leaf difference, "midpoint-like" rooting), after which
    //     Tree::createSets lists the merges in post-order (left subtree first). ClustalWRoot below is a line-by-line port
    //     of (b)-(c) (clustalw 2.1+lgpl-7 src/tree/Tree.cpp: createTree, insertNode, reRoot, calcMean, insertRoot,
    //     groupSeqs; src/tree/ClusterTreeOutput.cpp: printPhylipTree/twoWaySplit for the child order).

    /// <summary>
    /// A node of the progressive-alignment guide tree. Leaves carry a single input-sequence index; internal nodes carry
    /// two children (<see cref="Left"/> is merged first / listed first, as in ClustalW's alignment steps).
    /// </summary>
    internal sealed class ProgressiveGuideNode
    {
        public int LeafIndex = -1; // >= 0 for leaves, -1 for internal nodes
        public ProgressiveGuideNode? Left;
        public ProgressiveGuideNode? Right;
        public bool IsLeaf => Left == null && Right == null;
    }

    /// <summary>
    /// Builds the binary rooted guide tree of <paramref name="k"/> sequences from the symmetric distance matrix.
    /// </summary>
    internal static ProgressiveGuideNode BuildProgressiveGuideTree(int k, double[,] distance, MsaGuideTree method = MsaGuideTree.Upgma)
    {
        if (!Enum.IsDefined(method))
            throw new ArgumentOutOfRangeException(nameof(method), method, "Unknown guide-tree method.");

        return method == MsaGuideTree.NeighborJoining
            ? BuildClustalWNeighborJoiningGuideTree(k, distance).Rooted
            : BuildUpgmaGuideTree(k, distance);
    }

    /// <summary>Leaf names of the canonical builders: the decimal input index.</summary>
    private static List<string> GuideTaxa(int k)
    {
        var taxa = new List<string>(k);
        for (int i = 0; i < k; i++)
            taxa.Add(i.ToString(CultureInfo.InvariantCulture));
        return taxa;
    }

    /// <summary>
    /// UPGMA guide tree from <see cref="PhylogeneticAnalyzer.BuildTreeFromMatrix(IReadOnlyList{string}, double[,], PhylogeneticAnalyzer.TreeMethod, PhylogeneticAnalyzer.TreeBuildOptions?)"/>
    /// (proportional averaging, Sokal &amp; Michener 1958). Tie-breaking is deterministic: the lowest (i, j) index pair wins
    /// (distances equal within <see cref="ProfileScoreEpsilon"/> tie), so hand-derived test cases have a unique tree.
    /// </summary>
    private static ProgressiveGuideNode BuildUpgmaGuideTree(int k, double[,] distance)
    {
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            GuideTaxa(k), distance, PhylogeneticAnalyzer.TreeMethod.UPGMA,
            new PhylogeneticAnalyzer.TreeBuildOptions { TieTolerance = ProfileScoreEpsilon, TieOrder = PhylogeneticAnalyzer.PairScanOrder.FirstIndexMajor });
        return ToGuideNode(tree.Root);
    }

    private static ProgressiveGuideNode ToGuideNode(PhylogeneticAnalyzer.PhyloNode node)
    {
        if (node.IsLeaf)
            return new ProgressiveGuideNode { LeafIndex = int.Parse(node.Name, CultureInfo.InvariantCulture) };

        // UPGMA nodes are binary.
        return new ProgressiveGuideNode { Left = ToGuideNode(node.Children[0]), Right = ToGuideNode(node.Children[1]) };
    }

    /// <summary>The unrooted neighbour-joining tree and the rooted binary guide tree ClustalW derives from it.</summary>
    internal sealed record ClustalWGuideTreeResult(PhylogeneticAnalyzer.PhyloNode Unrooted, ProgressiveGuideNode Rooted);

    /// <summary>
    /// ClustalW's guide tree of <paramref name="k"/> sequences: canonical neighbour-joining
    /// (<see cref="PhylogeneticAnalyzer.TreeMethod.NeighborJoining"/>) with ClustalW's pair-scan order, then the ClustalW printing / rooting steps.
    /// </summary>
    internal static ClustalWGuideTreeResult BuildClustalWNeighborJoiningGuideTree(int k, double[,] distance)
    {
        if (k < 2)
            throw new ArgumentOutOfRangeException(nameof(k), k, "A guide tree needs at least two sequences.");

        var taxa = GuideTaxa(k);
        var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(
            taxa, distance, PhylogeneticAnalyzer.TreeMethod.NeighborJoining,
            new PhylogeneticAnalyzer.TreeBuildOptions { TieTolerance = ProfileScoreEpsilon, TieOrder = PhylogeneticAnalyzer.PairScanOrder.SecondIndexMajor });

        return new ClustalWGuideTreeResult(tree.Root, ClustalWRoot(tree.Root, k));
    }

    // ---- ClustalW Tree.cpp port --------------------------------------------------------------------------------

    /// <summary>TreeNode of ClustalW's Tree.cpp (single-precision distances).</summary>
    private sealed class CwNode
    {
        public CwNode? Left, Right, Parent;
        public float Dist;
        public bool Leaf;
        public int LeafIndex = -1;
    }

    /// <summary>
    /// Prints the branch length like ClustalW (<c>fixed &lt;&lt; setprecision(5)</c>, |b| &lt; 1e-4 reset to 0 first) and
    /// reads it back as a single-precision number (<c>file &gt;&gt; float</c>).
    /// </summary>
    private static float ClustalWBranch(double b)
    {
        if (Math.Abs(b) < 0.0001)
            b = 0.0;
        return float.Parse(b.ToString("F5", CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
    }

    private static int MinLeafIndex(PhylogeneticAnalyzer.PhyloNode node)
    {
        if (node.IsLeaf)
            return int.Parse(node.Name, CultureInfo.InvariantCulture);
        int min = int.MaxValue;
        foreach (var c in node.Children)
            min = Math.Min(min, MinLeafIndex(c));
        return min;
    }

    /// <summary>
    /// ClustalW <c>Tree::readTree</c> + <c>Tree::createSets</c> on the unrooted NJ tree it just wrote: rebuilds the node
    /// structure (<c>createTree</c>, root trifurcation turned into a root with a left (A,B) node and a right C by
    /// <c>insertNode</c>), re-roots it (<c>reRoot</c>) and returns the binary tree whose post-order (left first) is ClustalW's
    /// list of alignment steps. Two sequences: the single step (0, 1).
    /// </summary>
    private static ProgressiveGuideNode ClustalWRoot(PhylogeneticAnalyzer.PhyloNode unrooted, int k)
    {
        if (k == 2)
        {
            return new ProgressiveGuideNode
            {
                Left = new ProgressiveGuideNode { LeafIndex = 0 },
                Right = new ProgressiveGuideNode { LeafIndex = 1 },
            };
        }

        var ptrs = new List<CwNode>();
        var lptr = new List<CwNode>();

        CwNode Build(PhylogeneticAnalyzer.PhyloNode src, CwNode? parent)
        {
            var node = new CwNode { Parent = parent, Dist = ClustalWBranch(src.BranchLength) };
            ptrs.Add(node);
            if (src.IsLeaf)
            {
                node.Leaf = true;
                node.LeafIndex = int.Parse(src.Name, CultureInfo.InvariantCulture);
                lptr.Add(node);
                return node;
            }

            // printPhylipTree / twoWaySplit: children are listed by their lowest leaf index.
            var kids = src.Children.OrderBy(MinLeafIndex).ToList();
            node.Left = Build(kids[0], node);
            node.Right = Build(kids[1], node);
            return node;
        }

        // NJ centre: (c1, c2, c3) listed by lowest leaf index.
        var centre = unrooted.Children.OrderBy(MinLeafIndex).ToList();

        // createTree on the file's outermost '(': the node X = (c1, c2) ...
        var x = new CwNode(); // seqTree: dist 0, parent set below
        ptrs.Add(x);
        x.Left = Build(centre[0], x);
        x.Right = Build(centre[1], x);

        // ... a third ',' triggers insertNode: a new root N above X; X keeps dist 0, N gets c3 as right child.
        var n = new CwNode { Left = x, Parent = null, Dist = 0f };
        x.Parent = n;
        ptrs.Add(n);
        n.Right = Build(centre[2], n);

        CwNode root = ReRoot(x, ptrs, lptr, k);
        return ToGuideNode(root);
    }

    private static ProgressiveGuideNode ToGuideNode(CwNode p)
    {
        if (p.Leaf)
            return new ProgressiveGuideNode { LeafIndex = p.LeafIndex };
        return new ProgressiveGuideNode { Left = ToGuideNode(p.Left!), Right = ToGuideNode(p.Right!) };
    }

    /// <summary>Port of <c>Tree::reRoot</c>.</summary>
    private static CwNode ReRoot(CwNode ptree, List<CwNode> ptrs, List<CwNode> lptr, int nseqs)
    {
        float minDiff = 0.0f, minDepth = 1.0f;
        bool first = true;

        // find the difference between the means of leaf->node distances on the left and on the right of each node
        CwNode rootPtr = ptree;
        foreach (var p in ptrs)
        {
            if (p.Parent == null)
                continue; // AW bug 94 guard of the C++ source

            float diff = CalcMean(p, lptr, nseqs, out float maxDist);
            if ((diff == 0 || (diff > 0 && diff < 2 * p.Dist)) && (maxDist < minDepth || first))
            {
                first = false;
                rootPtr = p;
                minDepth = maxDist;
                minDiff = diff;
            }
        }

        // insert a new node as the ancestor of the node which produces the shallowest tree
        if (rootPtr == ptree)
        {
            minDiff = rootPtr.Left!.Dist + rootPtr.Right!.Dist;
            rootPtr = rootPtr.Right!;
        }

        return InsertRoot(rootPtr, minDiff);
    }

    /// <summary>Port of <c>Tree::calcMean</c> (RIGHT = descendants of <paramref name="nptr"/>, LEFT = the rest).</summary>
    private static float CalcMean(CwNode nptr, List<CwNode> lptr, int nSeqs, out float maxDist)
    {
        float leftSum = 0.0f, rightSum = 0.0f;
        int depth = 0;
        int n = 0;
        int numLeft = 0, numRight = 0;
        var pathToRoot = new CwNode[nSeqs + 2];
        var distToNode = new float[nSeqs + 2];

        // all nodes between the selected node and the root
        maxDist = 0.0f;
        float dist = 0.0f;
        CwNode? p = nptr;
        while (p != null)
        {
            if (depth == pathToRoot.Length)
            {
                Array.Resize(ref pathToRoot, depth * 2);
                Array.Resize(ref distToNode, depth * 2);
            }

            pathToRoot[depth] = p;
            dist += p.Dist;
            distToNode[depth] = dist;
            p = p.Parent;
            depth++;
        }

        for (int i = 0; i < lptr.Count; i++)
        {
            p = lptr[i];
            bool right;
            if (p == nptr)
            {
                right = true;
                dist = 0.0f;
            }
            else
            {
                right = false;
                dist = 0.0f;

                // find the common ancestor
                bool found = false;
                n = 0;
                while (!found && p!.Parent != null)
                {
                    for (int j = 0; j < depth; j++)
                    {
                        if (p.Parent == pathToRoot[j])
                        {
                            found = true;
                            n = j;
                        }
                    }

                    dist += p.Dist;
                    p = p.Parent;
                }

                if (p == nptr)
                    right = true;
            }

            if (!right)
            {
                leftSum += dist;
                leftSum += distToNode[n - 1];
                numLeft++;
            }
            else
            {
                rightSum += dist;
                numRight++;
            }

            if (dist > maxDist)
                maxDist = dist;
        }

        float leftMean = leftSum / numLeft;
        float rightMean = rightSum / numRight;
        return leftMean - rightMean;
    }

    /// <summary>Port of <c>Tree::insertRoot</c>: a new root above <paramref name="p"/>, the old path re-oriented, the old root removed.</summary>
    private static CwNode InsertRoot(CwNode p, float diff)
    {
        var newp = new CwNode();

        CwNode t = p.Parent!;
        float prevDist = t.Dist;

        p.Parent = newp;
        float dist = p.Dist;
        p.Dist = diff / 2;
        if (p.Dist < 0.0f)
            p.Dist = 0.0f;
        if (p.Dist > dist)
            p.Dist = dist;

        t.Dist = dist - p.Dist;

        newp.Left = t;
        newp.Right = p;
        newp.Parent = null;
        newp.Dist = 0.0f;
        newp.Leaf = false;

        if (t.Left == p)
            t.Left = t.Parent;
        else
            t.Right = t.Parent;

        CwNode prev = t;
        CwNode? q = t.Parent;
        t.Parent = newp;

        while (q != null)
        {
            float td;
            if (q.Left == prev)
            {
                q.Left = q.Parent;
                q.Parent = prev;
                td = q.Dist;
                q.Dist = prevDist;
                prevDist = td;
                prev = q;
                q = q.Left;
            }
            else
            {
                q.Right = q.Parent;
                q.Parent = prev;
                td = q.Dist;
                q.Dist = prevDist;
                prevDist = td;
                prev = q;
                q = q.Right;
            }
        }

        // remove the old root node
        q = prev;
        if (q.Left == null)
        {
            float d = q.Dist;
            q = q.Right!;
            q.Dist += d;
            q.Parent = prev.Parent;
            if (prev.Parent!.Left == prev)
                prev.Parent.Left = q;
            else
                prev.Parent.Right = q;
            prev.Right = null;
        }
        else
        {
            float d = q.Dist;
            q = q.Left;
            q.Dist += d;
            q.Parent = prev.Parent;
            if (prev.Parent!.Left == prev)
                prev.Parent.Left = q;
            else
                prev.Parent.Right = q;
            prev.Left = null;
        }

        return newp;
    }

    /// <summary>
    /// ClustalW's progressive-alignment steps: the internal nodes of the rooted guide tree in post-order, left subtree first
    /// (<c>Tree::groupSeqs</c>), each as the number of sequences it merges — the "Group n: Sequences: s" lines of
    /// <c>clustalw -ALIGN</c>.
    /// </summary>
    internal static List<int> GuideTreeStepSizes(ProgressiveGuideNode root)
    {
        var sizes = new List<int>();

        int Visit(ProgressiveGuideNode node)
        {
            if (node.IsLeaf)
                return 1;
            int s = Visit(node.Left!) + Visit(node.Right!);
            sizes.Add(s);
            return s;
        }

        Visit(root);
        return sizes;
    }

    /// <summary>
    /// ClustalW's guide-tree distance of two sequences in single precision exactly as <c>FullPairwiseAlign</c> stores it:
    /// <c>mmScore = 100f·count; mmScore /= (float)MIN(len1, len2); d = (100f − mmScore) / 100f</c> (an empty sequence: 1).
    /// </summary>
    internal static double ClustalWSinglePrecisionDistance(string a, string b)
    {
        if (a.Length == 0 || b.Length == 0)
            return 1.0;

        var r = ClustalWPairwiseAlign(a, b);
        float mm = 100.0f * r.Identities;
        mm /= Math.Min(a.Length, b.Length);
        return (100.0f - mm) / 100.0f;
    }

    #endregion
}
