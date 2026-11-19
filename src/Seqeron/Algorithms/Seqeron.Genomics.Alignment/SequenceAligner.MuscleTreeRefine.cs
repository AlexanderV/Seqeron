namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MuscleEngine
    {
        // ---------------------------------------------------------------------------------------------------------------------------------
        // treefrommsa.cpp TreeFromMSA / distcalc.cpp DistCalcMSA: the guide tree of a (draft) alignment

        private float[,] MsaDistances(IReadOnlyList<char[]> rows, MuscleDistance distance, bool neighborJoining)
        {
            int n = rows.Count;
            var d = new float[n, n];
            for (int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    float v;
                    switch (distance)
                    {
                        // TreeFromMSA_NJ computes the distances with MSADist in double precision (msadist.h), TreeFromMSA_UPGMA with DistCalcMSA from the single-precision identity
                        case MuscleDistance.PctIdKimura:
                            v = (float)KimuraDist(neighborJoining ? PctIdentityPair(rows[i], rows[j]) : (float)PctIdentityPair(rows[i], rows[j]));
                            break;
                        case MuscleDistance.PctIdLog:
                            v = (float)PctIdToMafftDist(neighborJoining ? PctIdentityPair(rows[i], rows[j]) : (float)PctIdentityPair(rows[i], rows[j]));
                            break;
                        case MuscleDistance.Edit when !neighborJoining:
                            {
                                float pctId = (float)PctIdentityPair(rows[i], rows[j]);
                                if (pctId > 1.0)
                                    throw new InvalidOperationException("Internal error, DISTANCE_Edit");
                                v = (float)1.0 - pctId;
                                break;
                            }
                        case MuscleDistance.ScoreDist:
                            v = (float)ScoreDist(rows[i], rows[j]);
                            break;
                        default:
                            throw new ArgumentException($"{distance} is not a distance of an aligned pair (use pctid_kimura, pctid_log, edit or scoredist for -distance2).");
                    }

                    d[i, j] = v;
                    d[j, i] = v;
                }

            return d;
        }

        private static double PctIdToMafftDist(double pctId)
        {
            if (pctId < 0.05)
                pctId = 0.05;
            return -Math.Log(pctId);
        }

        private MuscleTree TreeFromMsaRows(int[] ids, IReadOnlyList<char[]> rows, MuscleRoot? root = null)
        {
            var distance = _o.Distance2 == MuscleDistance.Default ? MuscleDistance.PctIdKimura : _o.Distance2;
            var tree = TreeFromMsaRows(ids, rows, _o.Cluster2, distance, root ?? _o.Root2, out string text);
            Dist2Text = text;
            return tree;
        }

        // TreeFromMSA(msa, tree, Cluster, Distance, Root)
        private MuscleTree TreeFromMsaRows(int[] ids, IReadOnlyList<char[]> rows, MuscleClusterMethod cluster, MuscleDistance distance, MuscleRoot root, out string distanceText)
        {
            bool nj = cluster == MuscleClusterMethod.NeighborJoining;
            var dist = MsaDistances(rows, distance, nj);
            distanceText = DistanceMatrixText(dist, rows.Count, nj, ids);
            return TreeFromDistances(dist, ids, cluster, root);
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // refinetree.cpp RefineTree (the pipeline of a weighting other than ClustalW; g_uMaxTreeRefineIters rounds)

        private void RefineTreeMsa(MuscleMsa msa, ref MuscleTree tree)
        {
            int seqCount = msa.SeqCount;
            if (tree.LeafCount != seqCount)
                throw new InvalidOperationException("Refine tree, tree has different number of nodes");
            if (seqCount < 3)
                return;

            int diffsCount = seqCount;
            for (int iter = 0; iter < _o.MaxTrees; iter++)
            {
                var tree2 = TreeFromMsaRows(msa.Ids, msa.Rows);
                var (diffs, idToDiffsLeaf) = DiffTrees(tree2, tree);
                tree = tree2;

                int newDiffsNodeCount = diffs.Count;
                int newDiffsCount = (newDiffsNodeCount - 1) / 2;
                if (0 == newDiffsCount || newDiffsCount >= diffsCount)
                    break;
                diffsCount = newDiffsCount;

                var msa2 = RealignDiffs(msa, diffs, idToDiffsLeaf);
                msa.CopyFrom(msa2);
            }
        }

        // difftrees.cpp DiffTrees: the diffs tree (node 0 = root, children appended in pairs) and, for every sequence id, the diffs leaf it belongs to
        internal static (List<(int Left, int Right)> Diffs, int[] IdToDiffsLeaf) DiffTrees(MuscleTree tree1, MuscleTree tree2)
        {
            int nodeCount = tree1.NodeCount;
            int leafCount = tree1.LeafCount;
            if (nodeCount != tree2.NodeCount)
                throw new InvalidOperationException("DiffTrees: different node counts");
            var nodeIndexToId1 = new int[nodeCount];
            var idToNodeIndex2 = new int[nodeCount];
            var bachelor1 = new bool[nodeCount];
            var diff1 = new bool[nodeCount];
            for (int i = 0; i < nodeCount; i++)
            {
                nodeIndexToId1[i] = nodeCount;
                idToNodeIndex2[i] = nodeCount;
            }

            for (int i = 0; i < nodeCount; i++)
            {
                if (tree1.IsLeaf(i))
                {
                    int id = tree1.LeafId[i];
                    if (id >= nodeCount)
                        throw new InvalidOperationException("Diff trees requires existing leaf ids in range 0 .. (N-1)");
                    nodeIndexToId1[i] = id;
                }

                if (tree2.IsLeaf(i))
                {
                    int id = tree2.LeafId[i];
                    if (id >= nodeCount)
                        throw new InvalidOperationException("Diff trees requires existing leaf ids in range 0 .. (N-1)");
                    idToNodeIndex2[id] = i;
                }
            }

            for (int id = 0; id < leafCount; id++)
                if (idToNodeIndex2[id] == nodeCount)
                    throw new InvalidOperationException("DiffTrees, check 2");

            int internalId = leafCount;
            foreach (int n1 in tree1.PostOrder())
            {
                if (tree1.IsLeaf(n1) || bachelor1[n1])
                    continue;
                int left1 = tree1.Left[n1];
                if (bachelor1[left1])
                {
                    bachelor1[n1] = true;
                    continue;
                }

                int right1 = tree1.Right[n1];
                if (bachelor1[right1])
                {
                    bachelor1[n1] = true;
                    continue;
                }

                int idLeft = nodeIndexToId1[left1], idRight = nodeIndexToId1[right1];
                if (idLeft == nodeCount || idRight == nodeCount)
                    throw new InvalidOperationException("DiffTrees, check 5");
                int left2 = idToNodeIndex2[idLeft], right2 = idToNodeIndex2[idRight];
                if (left2 == nodeCount || right2 == nodeCount)
                    throw new InvalidOperationException("DiffTrees, check 6");
                int parentLeft2 = tree2.Parent[left2], parentRight2 = tree2.Parent[right2];
                if (parentLeft2 == parentRight2)
                {
                    nodeIndexToId1[n1] = internalId;
                    idToNodeIndex2[internalId] = parentLeft2;
                    internalId++;
                }
                else
                {
                    bachelor1[n1] = true;
                }
            }

            for (int n = 0; n < nodeCount; n++)
            {
                if (bachelor1[n])
                    continue;
                if (n == tree1.Root)
                {
                    diff1[n] = true;
                    continue;
                }

                if (bachelor1[tree1.Parent[n]])
                    diff1[n] = true;
            }

            var diffs = new List<(int Left, int Right)> { (-1, -1) };
            var idToDiffsLeaf = new int[leafCount];
            Array.Fill(idToDiffsLeaf, nodeCount);

            void Build(int treeNode, int diffsNode)
            {
                if (diff1[treeNode])
                {
                    var leaves = new List<int>();
                    GetLeaves(tree1, treeNode, leaves);
                    foreach (int leaf in leaves)
                    {
                        int id = tree1.LeafId[leaf];
                        if (id >= tree1.LeafCount)
                            throw new InvalidOperationException("BuildDiffs, id out of range");
                        idToDiffsLeaf[id] = diffsNode;
                    }

                    return;
                }

                if (tree1.IsLeaf(treeNode))
                    throw new InvalidOperationException("BuildDiffs: should never reach leaf");
                int dl = diffs.Count, dr = dl + 1;
                diffs.Add((-1, -1));
                diffs.Add((-1, -1));
                diffs[diffsNode] = (dl, dr);
                Build(tree1.Left[treeNode], dl);
                Build(tree1.Right[treeNode], dr);
            }

            Build(tree1.Root, 0);
            for (int id = 0; id < leafCount; id++)
                if (idToDiffsLeaf[id] == nodeCount)
                    throw new InvalidOperationException("TreeDiffs check 7");
            return (diffs, idToDiffsLeaf);
        }

        // realigndiffs.cpp RealignDiffs: progressive alignment of the diffs tree (the leaves are the unchanged subtrees, kept as they are with their gap-only columns removed)
        private MuscleMsa RealignDiffs(MuscleMsa msaIn, List<(int Left, int Right)> diffs, int[] idToDiffsLeaf)
        {
            var nodes = new MuscleMsa?[diffs.Count];
            void Visit(int node)
            {
                var (l, r) = diffs[node];
                if (l < 0)
                {
                    var ids = new List<int>();
                    for (int id = 0; id < msaIn.SeqCount; id++)
                        if (idToDiffsLeaf[id] == node)
                            ids.Add(id);
                    if (ids.Count == 0)
                        throw new InvalidOperationException("MakeNode: no seqs in diff");
                    var m = msaIn.SubsetByIds(ids);
                    m.DeleteGappedColumns();
                    nodes[node] = m;
                    return;
                }

                Visit(l);
                Visit(r);
                nodes[node] = AlignTwoMsas(nodes[l]!, nodes[r]!, false, false).Merged;
            }

            Visit(0);
            return nodes[0]!;
        }
    }
}
