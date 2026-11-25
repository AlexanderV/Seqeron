// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;

namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

public static partial class ClustalWAligner
{
    // ClustalW 2.1 `-CLUSTERING=UPGMA` (audit A2-14, B09 F92): src/tree/UPGMA/{UPGMAAlgorithm,Node,RootedGuideTree}.cpp and TreeInterface::getWeightsAndStepsFromDistMatUPGMA, ported with the pointer
    // arithmetic of the C++ turned into offsets of the packed lower-triangular distance array (row r starts at r(r+1)/2, entry j < r), the linked list of live clusters into `Next` references and the
    // `Node**` of the min-distance search into the PREDECESSOR of the node (the location the C++ pointer points at).

    internal sealed class UpgmaNode
    {
        public UpgmaNode? Next, Left, Right;
        public int Size = 1;
        public int SeqNum;          // 1-based sequence number, -1 for an internal node
        public int Order;           // number of leaves below (RootedGuideTree::orderNodes)
        public int Row = -1;        // offset of the node's row in the distance array, -1 = none (the first node)
        public int NumDists;
        public int IndexToMinDist = -1;
        public double MinDist = double.MaxValue;
        public double Height;
        public List<int> AllElements = [];

        public UpgmaNode Copy() => new()
        {
            Next = Next, Left = Left, Right = Right, Size = Size, SeqNum = SeqNum, Order = Order, Row = Row, NumDists = NumDists,
            IndexToMinDist = IndexToMinDist, MinDist = MinDist, Height = Height, AllElements = [.. AllElements],
        };
    }

    /// <summary>The result of <c>UPGMAAlgorithm::generateTree</c>: the rooted tree and the progressive alignment steps (1-based group vectors: 1 / 2 = the two groups of the join, 0 = not involved).</summary>
    internal sealed record UpgmaResult(UpgmaNode Root, List<int[]> Steps);

    private const double UpgmaBlankDist = -1;

    // UPGMAAlgorithm::generateTree + doUPGMA (distances: the n x n matrix of the pairwise stage, 0-based)
    internal static UpgmaResult ClustalWUpgmaTree(double[,] dist, int numSeqs)
    {
        // DistMatrix::elements: entry (i, j), i > j (0-based) at i(i+1)/2 + j; the diagonal slot i(i+3)/2 and slot 0 stay 0
        var elements = new double[((numSeqs + 1) * (numSeqs + 2)) / 2];
        for (int i = 1; i < numSeqs; i++)
            for (int j = 0; j < i; j++)
                elements[i * (i + 1) / 2 + j] = dist[i, j];

        // initialiseNodes
        var nodes = new UpgmaNode[numSeqs];
        nodes[0] = new UpgmaNode { SeqNum = 1, AllElements = [1] };
        for (int e = 1; e < numSeqs; e++)
        {
            var node = new UpgmaNode { SeqNum = e + 1, Row = e * (e + 1) / 2, NumDists = e, AllElements = [e + 1] };
            FindMinDist(node, elements);
            nodes[e - 1].Next = node;
            nodes[e] = node;
        }

        var steps = new List<int[]>();
        int orderNode1 = 0, orderNode2 = 0, orderNew = 0;
        double CalcNewDist(double d1, double d2) => ((d1 * orderNode1) + (d2 * orderNode2)) / orderNew;
        void MovePastUnused(ref int p)
        {
            while (elements[p] < 0)
                p++;
        }

        while (nodes[0].Next != null)
        {
            // getNodeWithMinDist: the predecessor of the node with the smallest minDist (the first one wins a tie)
            UpgmaNode? predOfMin = null;
            double minDistance = double.MaxValue;
            for (var pred = nodes[0]; pred.Next != null; pred = pred.Next)
            {
                if (pred.Next.MinDist < minDistance)
                {
                    minDistance = pred.Next.MinDist;
                    predOfMin = pred;
                }
            }

            var node2 = predOfMin!.Next!;
            int indexToMinDist = node2.IndexToMinDist;
            var node1 = nodes[indexToMinDist];
            orderNode1 = node1.Size;
            orderNode2 = node2.Size;
            orderNew = orderNode1 + orderNode2;

            // STEP 2: recompute the distance row of node1
            int it2 = node2.Row;
            if (indexToMinDist != 0)
            {
                int it1 = node1.Row;
                elements[it1] = CalcNewDist(elements[it1], elements[it2]);
                int minIndex1 = it1;
                it1++;
                it2++;
                for (int num = node1.NumDists - 1; num > 0; num--)
                {
                    if (elements[it1] >= 0)
                    {
                        elements[it1] = CalcNewDist(elements[it1], elements[it2]);
                        if (elements[it1] < elements[minIndex1])
                            minIndex1 = it1;
                    }

                    it1++;
                    it2++;
                }

                node1.MinDist = elements[minIndex1];
                node1.IndexToMinDist = minIndex1 - node1.Row;
            }

            // STEP 3A: computeDistsUpToNodeToJoin2
            int idx2 = node2.IndexToMinDist;
            MovePastUnused(ref it2);
            for (var iter = node1.Next; iter != node2; iter = iter!.Next)
            {
                it2++;
                MovePastUnused(ref it2);
                double distToNode = elements[iter!.Row + idx2];
                double newDist = CalcNewDist(distToNode, elements[it2]);
                elements[iter.Row + idx2] = newDist;
                if (newDist < iter.MinDist)
                {
                    iter.MinDist = newDist;
                    iter.IndexToMinDist = idx2;
                }
                else if (newDist > iter.MinDist && iter.IndexToMinDist == idx2)
                    FindMinDist(iter, elements);
            }

            // STEP 3B: computeDistsForNodesAfterNode2
            int indexToNode2 = node2.NumDists;
            const double compareRelEpsilon = 1E-14;
            for (var iter = node2.Next; iter != null; iter = iter.Next)
            {
                int slot = iter.Row + idx2;
                double distUpdate = ((elements[slot] * orderNode1) + (elements[iter.Row + indexToNode2] * orderNode2)) / orderNew;
                elements[slot] = distUpdate;
                if (distUpdate < iter.MinDist && (iter.MinDist - distUpdate) / iter.MinDist > compareRelEpsilon)
                {
                    iter.MinDist = distUpdate;
                    iter.IndexToMinDist = idx2;
                }
                else if ((distUpdate > iter.MinDist && (distUpdate - iter.MinDist) / distUpdate > compareRelEpsilon && iter.IndexToMinDist == idx2) || iter.IndexToMinDist == indexToNode2)
                {
                    elements[iter.Row + indexToNode2] = UpgmaBlankDist;
                    FindMinDist(iter, elements);
                }
                else
                    elements[iter.Row + indexToNode2] = UpgmaBlankDist;
            }

            // STEP 4: addAlignmentStep (the groups are indexed by the 1-based sequence numbers)
            var groups = new int[numSeqs + 1];
            foreach (int s in node1.AllElements.Where(s => s < groups.Length))
                groups[s] = 1;
            foreach (int s in node2.AllElements.Where(s => s < groups.Length))
                groups[s] = 2;
            steps.Add(groups);

            // STEP 5: Node::merge (height = minDist / 2; both children, whatever their own height, get the height of the join)
            double height = node2.MinDist / 2.0;
            var left = node1.Copy();
            left.Row = -1;
            node1.Left = left;
            node1.Right = node2;
            node1.Size = left.Size + node2.Size;
            node1.SeqNum = -1;
            node1.Height = height;
            left.Height = height;
            node2.Height = height;
            node1.AllElements.AddRange(node2.AllElements);
            node2.AllElements.Clear();
            predOfMin.Next = node2.Next;
        }

        return new UpgmaResult(nodes[0], steps);
    }

    // RootedClusterTree::treeFromDistMatrix + RootedTreeOutput::printPhylipTree: the guide-tree file `-NEWTREE` of `-CLUSTERING=UPGMA`. Two sequences: "(a:d,b:d);" with d = distance / 2 as a float printed with
    // setprecision(5) (general format); otherwise the rooted tree, "(\n" left ",\n" right "):height" for an internal node, "name:height" for a leaf (heights are doubles in the default stream
    // format, six significant digits), the outermost ":height" included, then ";\n".
    internal static string ClustalWUpgmaGuideTreeText(double[,] dist, IReadOnlyList<string> names)
    {
        int n = names.Count;
        var sb = new System.Text.StringBuilder();
        if (n == 2)
        {
            float d = (float)(dist[0, 1] / 2.0);
            string t = CppGeneral(d, 5);
            return $"({names[0]}:{t},{names[1]}:{t});\n";
        }

        var tree = ClustalWUpgmaTree(dist, n);
        void Write(UpgmaNode t)
        {
            if (t.Left == null && t.Right == null)
            {
                sb.Append(names[t.SeqNum - 1]).Append(':').Append(CppGeneral(t.Height, 6));
                return;
            }

            sb.Append("(\n");
            Write(t.Left!);
            sb.Append(",\n");
            Write(t.Right!);
            sb.Append("):").Append(CppGeneral(t.Height, 6));
        }

        Write(tree.Root);
        return sb.Append(";\n").ToString();
    }

    // the default (%g-like) formatting of a C++ ostream with the given precision
    private static string CppGeneral(double v, int precision)
    {
        if (v == 0)
            return "0";
        string s = v.ToString("G" + precision, CultureInfo.InvariantCulture);
        int e = s.IndexOf('E');
        if (e < 0)
            return s;
        char sign = s[e + 1];
        string digits = s[(e + 2)..].TrimStart('0').PadLeft(2, '0');
        return s[..e] + "e" + sign + digits;
    }

    /// <summary>The <c>-NEWTREE</c> text of ClustalW's DNA <c>-CLUSTERING=UPGMA</c> run (guide tree of the pairwise-stage distances; sequences are named by <paramref name="names"/>).</summary>
    internal static string ClustalWUpgmaGuideTreeText(string[] raw, ClustalWOptions options, IReadOnlyList<string> names) =>
        ClustalWUpgmaGuideTreeText(ClustalWPairwiseDistanceMatrix(raw, options), names);

    /// <summary>The <c>-NEWTREE</c> text of ClustalW's PROTEIN <c>-CLUSTERING=UPGMA</c> run.</summary>
    internal static string ClustalWUpgmaGuideTreeText(string[] raw, ClustalWProteinOptions options, IReadOnlyList<string> names) =>
        ClustalWUpgmaGuideTreeText(ClustalWProteinDistanceMatrix(raw.Select(s => s.ToUpperInvariant()).ToArray(), options), names);

    // Node::findMinDist: the first entry is the initial minimum even when it is negative; the others replace it when they are non-negative and strictly smaller
    private static void FindMinDist(UpgmaNode node, double[] elements)
    {
        int min = node.Row;
        for (int i = 1; i < node.NumDists; i++)
        {
            int p = node.Row + i;
            if (elements[p] >= 0 && elements[p] < elements[min])
                min = p;
        }

        node.MinDist = elements[min];
        node.IndexToMinDist = min - node.Row;
    }

    // RootedGuideTree::calcSeqWeights(0, n, ...) for n >= 2 with distance tree and weights on: orderNodes, doWeightCalc (float arithmetic), normalisation to INT_SCALE_FACTOR
    internal static int[] ClustalWUpgmaSeqWeights(UpgmaNode? root, int n)
    {
        var weights = new float[n];
        if (root != null)
        {
            int Order(UpgmaNode t)
            {
                t.Order = t.Left == null && t.Right == null ? 1 : Order(t.Left!) + Order(t.Right!);
                return t.Order;
            }

            void Walk(float soFar, UpgmaNode t)
            {
                if (t.Left == null && t.Right == null)
                {
                    weights[t.SeqNum - 1] = soFar;
                    return;
                }

                float w = (float)((double)soFar + (t.Height / t.Order));
                Walk(w, t.Left!);
                Walk(w, t.Right!);
            }

            Order(root);
            Walk(0.0f, root);
        }

        var seqWeight = new int[n];
        for (int i = 0; i < n; i++)
            seqWeight[i] = (int)(weights[i] * 100);
        int sum = 0;
        for (int i = 0; i < n; i++)
            sum += seqWeight[i];
        if (sum == 0)
        {
            for (int i = 0; i < n; i++)
                seqWeight[i] = 1;
            sum = n;
        }

        for (int i = 0; i < n; i++)
        {
            seqWeight[i] = seqWeight[i] * CwIntScaleFactor / sum;
            if (seqWeight[i] < 1)
                seqWeight[i] = 1;
        }

        return seqWeight;
    }
}
