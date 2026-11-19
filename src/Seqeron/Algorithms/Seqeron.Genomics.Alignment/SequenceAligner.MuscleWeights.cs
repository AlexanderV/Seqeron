namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MuscleEngine
    {

        // ---------------------------------------------------------------------------------------------------------------------------------
        // setblosumweights.cpp / setgscweights.cpp / cluster.cpp ClusterTree::Create (single linkage, the cluster tree of the HMMer weights)

        private sealed class ClusterNode
        {
            public double Weight, Weight2;
            public ClusterNode? Left, Right, Parent, NextDisjoint, PrevDisjoint;
            public int Index;

            public int ClusterSize => Left is null && Right is null ? 1 : (Left?.ClusterSize ?? 0) + (Right?.ClusterSize ?? 0);

            public double ClusterWeight => (Left?.ClusterWeight ?? 0.0) + (Right?.ClusterWeight ?? 0.0) + Weight;
        }

        private static ClusterNode[] CreateClusterTree(float[,] dist, int leafCount)
        {
            int nodeCount = 2 * leafCount - 1;
            var nodes = new ClusterNode[nodeCount];
            for (int i = 0; i < nodeCount; i++)
                nodes[i] = new ClusterNode { Index = i };
            for (int i = 0; i < leafCount - 1; i++)
                nodes[i].NextDisjoint = nodes[i + 1];
            for (int i = 1; i < leafCount; i++)
                nodes[i].PrevDisjoint = nodes[i - 1];
            ClusterNode? disjoints = nodes[0];
            var cd = new float[nodeCount, nodeCount];
            for (int i = 0; i < leafCount; i++)
                for (int j = 0; j < leafCount; j++)
                    cd[i, j] = dist[i, j];

            void Delete(ClusterNode node)
            {
                var prev = node.PrevDisjoint;
                var next = node.NextDisjoint;
                if (prev is not null)
                    prev.NextDisjoint = next;
                else
                    disjoints = next;
                if (next is not null)
                    next.PrevDisjoint = prev;
            }

            for (int join = leafCount; join < nodeCount; join++)
            {
                int closest1 = 0, closest2 = 0;
                double distClosest = 9e99;
                for (var n1 = disjoints; n1 is not null; n1 = n1.NextDisjoint)
                {
                    for (var n2 = n1.NextDisjoint; n2 is not null; n2 = n2.NextDisjoint)
                    {
                        double d = cd[n1.Index, n2.Index];
                        if (d < distClosest)
                        {
                            distClosest = d;
                            closest1 = n1.Index;
                            closest2 = n2.Index;
                        }
                    }
                }

                var jn = nodes[join];
                var c1 = nodes[closest1];
                var c2 = nodes[closest2];
                jn.Left = c1;
                jn.Right = c2;
                jn.Weight = distClosest;
                c1.Parent = jn;
                c2.Parent = jn;
                Delete(c1);
                Delete(c2);
                jn.NextDisjoint = disjoints;
                jn.PrevDisjoint = null;
                if (disjoints is not null)
                    disjoints.PrevDisjoint = jn;
                disjoints = jn;

                for (var p = disjoints; p is not null; p = p.NextDisjoint)
                {
                    int idx = p.Index;
                    float d1 = cd[idx, closest1];
                    float d2 = cd[idx, closest2];
                    float d = d1 < d2 ? d1 : d2;
                    cd[join, idx] = d;
                    cd[idx, join] = d;
                }
            }

            return nodes;
        }

        private static void SetSubtreeGscWeight(ClusterNode? node)
        {
            if (node is null)
                return;
            var parent = node.Parent!;
            double parentWeight2 = parent.Weight2;
            double parentClusterWeight = parent.ClusterWeight;
            if (0.0 == parentClusterWeight)
            {
                double thisSize = node.ClusterSize;
                double parentSize = parent.ClusterSize;
                node.Weight2 = parentWeight2 * thisSize / parentSize;
            }
            else
            {
                double thisClusterWeight = node.ClusterWeight;
                double parentWeight = parent.Weight;
                double num = thisClusterWeight + parentWeight;
                double denom = parentClusterWeight + parentWeight;
                node.Weight2 = parentWeight2 * (num / denom);
            }

            SetSubtreeGscWeight(node.Left);
            SetSubtreeGscWeight(node.Right);
        }

        private static void CopyGscWeights(ClusterNode node, float[] weights)
        {
            if (node.Right is null && node.Left is null)
            {
                weights[node.Index] = (float)node.Weight2;
                return;
            }

            if (node.Left is not null)
                CopyGscWeights(node.Left, weights);
            if (node.Right is not null)
                CopyGscWeights(node.Right, weights);
        }

        private static void SetGscWeights(MuscleMsa m)
        {
            int n = m.SeqCount;
            var dist = new float[n, n];
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                {
                    double d = PctIdentityPair(m.Rows[i], m.Rows[j]);
                    dist[i, j] = (float)(1.0 - d);
                    dist[j, i] = dist[i, j];
                }

            var nodes = CreateClusterTree(dist, n);
            var root = nodes[2 * n - 2];
            root.Weight2 = 1.0;
            SetSubtreeGscWeight(root.Left);
            SetSubtreeGscWeight(root.Right);
            CopyGscWeights(root, m.Weights);
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // threewaywt.cpp (CalcThreeWayWeights / SetThreeWayWeightsMuscle). The rooted tree is read as unrooted (the root deleted); the edge weights the program computes are not used
        // for the sequence weights, only the edge lengths.

        private static int[] Neighbors(MuscleTree t, int node) => [t.Parent[node], t.IsLeaf(node) ? -1 : t.Left[node], t.IsLeaf(node) ? -1 : t.Right[node]];

        private static bool IsEdge(MuscleTree t, int a, int b) => t.Parent[a] == b || t.Parent[b] == a;

        private static int FirstNeighbor(MuscleTree t, int node, int exclude)
        {
            foreach (int nb in Neighbors(t, node))
                if (nb != -1 && nb != exclude)
                    return nb;
            return -1;
        }

        private static int SecondNeighbor(MuscleTree t, int node, int exclude)
        {
            bool foundOne = false;
            foreach (int nb in Neighbors(t, node))
            {
                if (nb != -1 && nb != exclude)
                {
                    if (foundOne)
                        return nb;
                    foundOne = true;
                }
            }

            return -1;
        }

        private static int FirstNeighborUnrooted(MuscleTree t, int node, int other)
        {
            if (!IsEdge(t, node, other))
                return FirstNeighbor(t, node, t.Root);
            int nb = FirstNeighbor(t, node, other);
            return nb == t.Root ? FirstNeighbor(t, nb, node) : nb;
        }

        private static int SecondNeighborUnrooted(MuscleTree t, int node, int other)
        {
            if (!IsEdge(t, node, other))
                return SecondNeighbor(t, node, t.Root);
            int nb = SecondNeighbor(t, node, other);
            return nb == t.Root ? FirstNeighbor(t, nb, node) : nb;
        }

        private static double EdgeLengthUnrooted(MuscleTree t, int n1, int n2)
        {
            if (!IsEdge(t, n1, n2))
                return t.EdgeLength[n1] + t.EdgeLength[n2];
            return t.Parent[n1] == n2 ? t.EdgeLength[n1] : t.EdgeLength[n2];
        }

        private static void SetThreeWaySeqWeights(MuscleTree t, int from, int to, double pathWeight, float[] weights)
        {
            double thisLength = EdgeLengthUnrooted(t, from, to);
            if (t.IsLeaf(to))
            {
                weights[t.LeafId[to]] = (float)(pathWeight + thisLength);
                return;
            }

            int first = FirstNeighborUnrooted(t, to, from);
            int second = SecondNeighborUnrooted(t, to, from);
            pathWeight *= thisLength;
            SetThreeWaySeqWeights(t, to, first, pathWeight, weights);
            SetThreeWaySeqWeights(t, to, second, pathWeight, weights);
        }

        private void SetThreeWayWeights(MuscleMsa m)
        {
            if (_splitNode1 < 0 || _splitNode2 < 0)
            {
                SetHenikoffWeightsPb(m);
                return;
            }

            var t = _muscleTree ?? throw new InvalidOperationException("g_ptrMuscleTree = 0");
            int n1 = _splitNode1, n2 = _splitNode2;
            if (n1 == t.Root)
                n1 = FirstNeighbor(t, n1, n2);
            else if (n2 == t.Root)
                n2 = FirstNeighbor(t, n2, n1);
            var weights = new float[t.LeafCount];
            SetThreeWaySeqWeights(t, n1, n2, 0.0, weights);
            SetThreeWaySeqWeights(t, n2, n1, 0.0, weights);
            for (int i = 0; i < m.SeqCount; i++)
            {
                int id = m.Ids[i];
                if (id >= weights.Length)
                    throw new InvalidOperationException("SetThreeWayWeightsMuscle: id out of range");
                m.Weights[i] = weights[id];
            }

            NormalizeWeights(m, 1.0f);
        }
    }
}
