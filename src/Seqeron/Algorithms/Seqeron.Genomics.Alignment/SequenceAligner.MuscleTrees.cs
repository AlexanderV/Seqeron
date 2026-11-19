namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    internal sealed partial class MuscleEngine
    {
        // ---------------------------------------------------------------------------------------------------------------------------------
        // distfunc.cpp DistUnaligned: the distance of the first guide tree (TreeFromSeqVect, fastclust.cpp)

        private static readonly int[,] MuscleBlosum62 =
        {
            //  A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
            { 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2 },
            { 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 },
            { -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3 },
            { -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2 },
            { -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3 },
            { 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3 },
            { -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2 },
            { -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1 },
            { -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2 },
            { -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1 },
            { -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1 },
            { -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 },
            { -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3 },
            { -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1 },
            { -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2 },
            { 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2 },
            { 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2 },
            { 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1 },
            { -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2 },
            { -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7 },
        };

        private const double Blosum62Expected = -0.5209;

        // scoredist.cpp (Sonnhammer and Hollich 2005, BMC Bioinformatics 6:108), JTT calibration factor
        private double Sigma(char[] a, char[] b, out uint length)
        {
            uint len = 0;
            double score = 0;
            for (int c = 0; c < a.Length; c++)
            {
                int l1 = LetterEx(a[c]), l2 = LetterEx(b[c]);
                if (l1 >= 20 || l2 >= 20)
                    continue;
                len++;
                score += MuscleBlosum62[l1, l2];
            }

            length = len;
            return score;
        }

        internal double ScoreDist(char[] a, char[] b)
        {
            if (_nuc)
                throw new InvalidOperationException("Scoredist is only for amino acid sequences");
            double score = Sigma(a, b, out uint length);
            double sN = score - length * Blosum62Expected;
            double s11 = Sigma(a, a, out _);
            double s22 = Sigma(b, b, out length);
            double sUn = (s11 + s22) / 2 - length * Blosum62Expected;
            double d = 0.0;
            if (sUn != 0)
            {
                double ratio = sN / sUn;
                d = ratio < 0.001 ? 3.0 : -Math.Log(ratio);
            }

            return d * 1.2873;
        }

        internal float[,] UnalignedDistances()
        {
            switch (_distance1)
            {
                case MuscleDistance.Kmer6_6:
                case MuscleDistance.Kmer4_6:
                    return ComputeKmerDistances();
                case MuscleDistance.Kmer20_3:
                    return DistKmer20_3();
                case MuscleDistance.Kmer20_4:
                    return DistKmer20_4();
                case MuscleDistance.Kbit20_3:
                    return DistKbit20_3();
                case MuscleDistance.PwKimura:
                case MuscleDistance.PwScoreDist:
                    return DistPairwiseProfile(_distance1 == MuscleDistance.PwScoreDist);
                default:
                    throw new ArgumentException($"{_distance1} is not a distance of unaligned sequences (-distance1 kmer6_6, kmer20_3, kmer20_4, kbit20_3, kmer4_6, pwkimura, pwscoredist).");
            }
        }

        // fastdistjones.cpp DistKmer20_3: the value stored is the number of shared triples per minimum length (a similarity, as the program does)
        internal float[,] DistKmer20_3()
        {
            int n = _n;
            const int tripleCount = 20 * 20 * 20;
            var counts = new uint[tripleCount][];
            var seqCountOf = new int[tripleCount];
            var d = new float[n, n];
            for (int s = 0; s < n; s++)
            {
                var seq = _seq[s];
                for (int pos = 0; pos + 2 < seq.Length; pos++)
                {
                    int l1 = LetterEx(seq[pos]);
                    if (l1 >= 20)
                        continue;
                    int l2 = LetterEx(seq[pos + 1]);
                    if (l2 >= 20)
                        continue;
                    int l3 = LetterEx(seq[pos + 2]);
                    if (l3 >= 20)
                        continue;
                    int word = l1 + l2 * 20 + l3 * 20 * 20;
                    counts[word] ??= new uint[n];
                    if (counts[word][s] == 0)
                        seqCountOf[word]++;
                    counts[word][s]++;
                }
            }

            for (int word = 0; word < tripleCount; word++)
            {
                if (seqCountOf[word] == 0)
                    continue;
                var list = new List<int>();
                for (int s = 0; s < n; s++)
                {
                    if (counts[word]![s] > 0)
                    {
                        list.Add(s);
                        if (list.Count == seqCountOf[word])
                            break;
                    }
                }

                for (int a = 0; a < list.Count; a++)
                {
                    uint c1 = counts[word]![list[a]];
                    for (int b = 0; b < a; b++)
                    {
                        uint c2 = counts[word]![list[b]];
                        uint min = c1 < c2 ? c1 : c2;
                        double cur = d[list[a], list[b]];
                        d[list[a], list[b]] = (float)(cur + min);
                        d[list[b], list[a]] = d[list[a], list[b]];
                    }
                }
            }

            for (int s1 = 0; s1 < n; s1++)
            {
                for (int s2 = 0; s2 < s1; s2++)
                {
                    uint minLength = (uint)Math.Min(_seq[s1].Length, _seq[s2].Length);
                    if (minLength < 3)
                    {
                        d[s1, s2] = d[s2, s1] = 1.0f;
                        continue;
                    }

                    double tripleSum = d[s1, s2];
                    if (tripleSum == 0)
                    {
                        d[s1, s2] = d[s2, s1] = 1.0f;
                        continue;
                    }

                    d[s1, s2] = d[s2, s1] = (float)(tripleSum / (minLength - 2));
                }
            }

            return d;
        }

        // fastdistkmer.cpp FastDistKmer ("kmer20_4"): 4-mers of the 20-letter alphabet (byte counters, as the program), distance 1 - F; sequences shorter than 4 have no k-mers
        internal float[,] DistKmer20_4()
        {
            const int K = 4;
            const int table = 20 * 20 * 20 * 20;
            int n = _n;
            var d = new float[n, n];
            var letters = new int[n][];
            for (int s = 0; s < n; s++)
                letters[s] = _seq[s].Select(c => _isWildcard[c] ? _charToLetterEx['A'] : _charToLetterEx[c]).ToArray();

            static IEnumerable<int> Kmers(int[] l)
            {
                for (int p = 0; p + K <= l.Length; p++)
                    yield return ((l[p] * 20 + l[p + 1]) * 20 + l[p + 2]) * 20 + l[p + 3];
            }

            var counts1 = new byte[table];
            var counts2 = new byte[table];
            for (int s1 = 0; s1 < n - 1; s1++)
            {
                Array.Clear(counts1);
                foreach (int km in Kmers(letters[s1]))
                    unchecked { counts1[km]++; }
                for (int s2 = s1 + 1; s2 < n; s2++)
                {
                    Array.Clear(counts2);
                    foreach (int km in Kmers(letters[s2]))
                        unchecked { counts2[km]++; }
                    uint common = 0;
                    foreach (int km in Kmers(letters[s1]))
                    {
                        common += Math.Min(counts1[km], counts2[km]);
                        counts2[km] = 0;
                    }

                    uint minLength = (uint)Math.Min(letters[s1].Length, letters[s2].Length);
                    double f = (double)common / unchecked(minLength - K + 1);
                    if (0.0 == f)
                        f = 0.01;
                    d[s1, s2] = d[s2, s1] = (float)(1 - f);
                }
            }

            return d;
        }

        // fastdistkbit.cpp DistKbit20_3 ("kbit20_3"): bit vectors of the 3-mers, every second residue is skipped (the program's loop increments the position twice), the "common" count of the
        // program counts the bits set in either vector plus those of the second vector again (b & 0x101, see CommonBitCount below)
        internal float[,] DistKbit20_3()
        {
            int n = _n;
            var bits = new byte[n][];
            for (int s = 0; s < n; s++)
            {
                bits[s] = new byte[1000];
                SetKmerBitVector(_seq[s], bits[s]);
            }

            var d = new float[n, n];
            for (int s1 = 0; s1 < n; s1++)
            {
                for (int s2 = 0; s2 < s1; s2++)
                {
                    int count = 0;
                    for (int b = 0; b < 1000; b++)
                    {
                        // CommonBitCount: b = *p1 | (*p2 << 8); while (b) { if (b & 0x101) ++count; b >>= 1; }: the shifts 0..7 count the bit positions set in either byte,
                        // the shifts 8..15 find the bits of *p2 again at position 0 (position 8 is empty by then) - a set bit of the second sequence is counted twice
                        count += System.Numerics.BitOperations.PopCount((uint)(bits[s1][b] | bits[s2][b])) + System.Numerics.BitOperations.PopCount((uint)bits[s2][b]);
                    }
                    float dist = (float)count / (uint)Math.Min(_seq[s1].Length, _seq[s2].Length);
                    d[s1, s2] = d[s2, s1] = dist;
                }
            }

            return d;
        }

        private void SetKmerBitVector(char[] s, byte[] bits)
        {
            const uint k = 3;
            uint len = (uint)s.Length;
            uint i = 0, c = 0, h = 0;
            for (uint j = 0; j < k - 1; ++j)
            {
                if (i >= len)
                    break;   // (the program reads past the end of a sequence shorter than 2 here)
                uint x = (uint)_charToLetterEx[s[(int)i]];
                i++;
                if (x <= 19)
                {
                    c = c * 20 + x;
                }
                else
                {
                    c = 0;
                    h = j + 1;
                }
            }

            while (i < len)
            {
                uint x = (uint)_charToLetterEx[s[(int)i]];
                i++;
                if (x <= 19)
                {
                    c = (c * 20 + x) % 8000;
                }
                else
                {
                    c = 0;
                    h = i + k;
                }

                if (i >= h)
                    bits[c / 8] |= (byte)(1 << (int)(c % 8));
                i++;
            }
        }

        // distpwkimura.cpp / scoredist.cpp DistPWScoreDist: every pair is aligned (profile-profile with the Henikoff weights), the distance is the Kimura / score distance of the alignment
        private float[,] DistPairwiseProfile(bool scoreDist)
        {
            var save = _curWeight;
            _curWeight = MuscleSeqWeight.Henikoff;
            int n = _n;
            var d = new float[n, n];
            for (int s1 = 0; s1 < n; s1++)
            {
                var m1 = new MuscleMsa([s1], [(char[])_seq[s1].Clone()]);
                for (int s2 = 0; s2 < s1; s2++)
                {
                    var m2 = new MuscleMsa([s2], [(char[])_seq[s2].Clone()]);
                    var merged = AlignTwoMsas(m1, m2, false, false).Merged;
                    float f = scoreDist
                        ? (float)ScoreDist(merged.Rows[0], merged.Rows[1])
                        : (float)KimuraDist(PctIdentityPair(merged.Rows[0], merged.Rows[1]));
                    d[s1, s2] = d[s2, s1] = f;
                }
            }

            _curWeight = save;
            return d;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // clust.cpp Clust (join style neighbour joining) + phyfromclust.cpp Tree::FromClust

        private sealed class NjCluster
        {
            public readonly int LeafCount, NodeCount;
            public readonly int[] Left, Right, Size;
            public readonly float[] Length;
            public readonly float[] Dist;
            public int ClusterCount;
            private readonly int[] _next, _prev;
            private int _head = -1;

            public NjCluster(int leafCount)
            {
                LeafCount = leafCount;
                NodeCount = 2 * leafCount - 1;
                Left = new int[NodeCount];
                Right = new int[NodeCount];
                Size = new int[NodeCount];
                Length = new float[NodeCount];
                Dist = new float[NodeCount * (NodeCount - 1) / 2 + 1];
                _next = new int[NodeCount];
                _prev = new int[NodeCount];
                Array.Fill(Left, -1);
                Array.Fill(Right, -1);
                Array.Fill(_next, -1);
                Array.Fill(_prev, -1);
            }

            public static int Vi(int a, int b) => a >= b ? b + a * (a - 1) / 2 : a + b * (b - 1) / 2;

            public float GetDist(int a, int b) => Dist[Vi(a, b)];

            public void SetDist(int a, int b, float d) => Dist[Vi(a, b)] = d;

            public int First => _head;

            public int Next(int i) => _next[i];

            public void Add(int i)
            {
                if (_head != -1)
                    _prev[_head] = i;
                _next[i] = _head;
                _prev[i] = -1;
                _head = i;
            }

            public void Delete(int i)
            {
                int prev = _prev[i], next = _next[i];
                if (next != -1)
                    _prev[next] = prev;
                if (prev == -1)
                    _head = next;
                else
                    _next[prev] = next;
                _next[i] = -1;
                _prev[i] = -1;
            }

            public float CalcR(int node)
            {
                if (ClusterCount == 2)
                    return 0;
                float sum = 0;
                for (int i = First; i != -1; i = Next(i))
                {
                    if (i == node)
                        continue;
                    sum += GetDist(node, i);
                }

                return sum / (uint)(ClusterCount - 2);
            }

            public float Metric(int i, int j)
            {
                float ri = CalcR(i), rj = CalcR(j);
                float dij = GetDist(i, j);
                return dij - (ri + rj);
            }
        }

        internal static MuscleTree NeighborJoiningTree(float[,] dist, int[] ids)
        {
            int n = ids.Length;
            if (n <= 1)
                throw new InvalidOperationException("Clust::Create: no leaves");
            var c = new NjCluster(n);
            for (int i = 0; i < c.NodeCount; i++)
            {
                if (i < n)
                {
                    c.Size[i] = 1;
                    c.Add(i);
                }
            }

            for (int i = 0; i < n; i++)
                for (int j = 0; j < i; j++)
                    c.SetDist(i, j, dist[i, j]);

            c.ClusterCount = n;
            for (int node = n; node < c.NodeCount; node++)
            {
                // ChooseJoinNeighborJoining
                int minLeft = -1, minRight = -1;
                float minMetric = float.PositiveInfinity;
                for (int l = c.First; l != -1; l = c.Next(l))
                {
                    for (int r = c.Next(l); r != -1; r = c.Next(r))
                    {
                        float metric = c.Metric(l, r);
                        if (metric < minMetric)
                        {
                            minMetric = metric;
                            minLeft = l;
                            minRight = r;
                        }
                    }
                }

                float dLr = c.GetDist(minLeft, minRight);
                float rL = c.CalcR(minLeft), rR = c.CalcR(minRight);
                float leftLen = (dLr + rL - rR) / 2;
                float rightLen = (dLr - rL + rR) / 2;
                int newNode = c.NodeCount - c.ClusterCount + 1;
                c.Length[minLeft] = leftLen;
                c.Length[minRight] = rightLen;
                c.Left[newNode] = minLeft;
                c.Right[newNode] = minRight;
                c.Size[newNode] = c.Size[minLeft] + c.Size[minRight];
                c.Delete(minLeft);
                c.Delete(minRight);
                c.Add(newNode);
                c.ClusterCount--;
                for (int other = c.First; other != -1; other = c.Next(other))
                {
                    if (other == minLeft || other == minRight || other == newNode)
                        continue;
                    float dL = c.GetDist(minLeft, other);
                    float dR = c.GetDist(minRight, other);
                    float dLR = c.GetDist(minLeft, minRight);
                    c.SetDist(newNode, other, (dL + dR - dLR) / 2);
                }
            }

            var t = new MuscleTree { LeafCount = n };
            int nodes = c.NodeCount;
            t.Left = new int[nodes];
            t.Right = new int[nodes];
            t.Parent = new int[nodes];
            t.EdgeLength = new double[nodes];
            t.LeafId = new int[nodes];
            Array.Fill(t.Left, -1);
            Array.Fill(t.Right, -1);
            Array.Fill(t.Parent, -1);
            for (int i = 0; i < n; i++)
                t.LeafId[i] = ids[i];
            for (int node = n; node < nodes; node++)
            {
                int l = c.Left[node], r = c.Right[node];
                t.Left[node] = l;
                t.Right[node] = r;
                t.Parent[l] = node;
                t.Parent[r] = node;
                t.EdgeLength[l] = c.Length[l];
                t.EdgeLength[r] = c.Length[r];
            }

            t.Root = nodes - 1;
            return t;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // phy3.cpp FixRoot / FindRoot, phy.cpp UnrootByDeletingRoot / RootUnrootedTree / OrientParent

        private sealed class UnrootedTree
        {
            public int Count;
            public int[][] Nb = Array.Empty<int[]>();       // neighbour slots 1, 2, 3 of the program (parent, left, right of a rooted tree), -1 = none
            public double[][] Len = Array.Empty<double[]>();
            public int Root = -1;

            public int NeighborCount(int n) => Nb[n].Count(x => x != -1);

            public bool IsLeaf(int n) => NeighborCount(n) == 1;

            public int Slot(int n, int nb) => Array.IndexOf(Nb[n], nb);

            public double EdgeLength(int a, int b) => Len[a][Slot(a, b)];

            public void SetEdgeLength(int a, int b, double d)
            {
                Len[a][Slot(a, b)] = d;
                Len[b][Slot(b, a)] = d;
            }

            public void OrientParent(int node, int parent)
            {
                if (node == -1)
                    return;
                if (Nb[node][0] != parent && Nb[node][1] == parent)
                {
                    (Nb[node][1], Nb[node][0]) = (Nb[node][0], parent);
                    (Len[node][1], Len[node][0]) = (Len[node][0], Len[node][1]);
                }
                else if (Nb[node][0] != parent)
                {
                    (Nb[node][2], Nb[node][0]) = (Nb[node][0], parent);
                    (Len[node][2], Len[node][0]) = (Len[node][0], Len[node][2]);
                }

                OrientParent(Nb[node][1], node);
                OrientParent(Nb[node][2], node);
            }
        }

        private struct EdgeInfo
        {
            public int Node1, Node2;
            public double MaxDistToLeaf, TotalDistToLeaves;
            public int MaxStep, MostDistantLeaf, LeafCount;
        }

        internal static MuscleTree FixRoot(MuscleTree t, MuscleRoot method)
        {
            if (method == MuscleRoot.Pseudo)
                return t;

            // UnrootByDeletingRoot (the root is the last node)
            int count = t.NodeCount;
            var u = new UnrootedTree { Count = count - 1, Nb = new int[count - 1][], Len = new double[count - 1][] };
            for (int i = 0; i < count - 1; i++)
            {
                bool leaf = t.IsLeaf(i);
                u.Nb[i] = [t.Parent[i], leaf ? -1 : t.Left[i], leaf ? -1 : t.Right[i]];
                u.Len[i] = [t.EdgeLength[i], leaf ? 0 : t.EdgeLength[t.Left[i]], leaf ? 0 : t.EdgeLength[t.Right[i]]];
            }

            int rl = t.Left[t.Root], rr = t.Right[t.Root];
            u.Nb[rl][0] = rr;
            u.Nb[rr][0] = rl;
            double joined = t.EdgeLength[rl] + t.EdgeLength[rr];
            u.Len[rl][0] = joined;
            u.Len[rr][0] = joined;

            var (node1, node2, length1, length2) = FindRoot(u, method);

            // RootUnrootedTree
            int root = u.Count;
            u.Count++;
            u.Nb = u.Nb.Append([-1, node1, node2]).ToArray();
            u.Len = u.Len.Append([0, 0, 0]).ToArray();
            int s1 = u.Slot(node1, node2);
            u.Nb[node1][s1] = root;
            int s2 = u.Slot(node2, node1);
            u.Nb[node2][s2] = root;
            u.OrientParent(node1, root);
            u.OrientParent(node2, root);
            u.SetEdgeLength(root, node1, length1);
            u.SetEdgeLength(root, node2, length2);

            var r = new MuscleTree { LeafCount = t.LeafCount, Root = root, LeafId = t.LeafId };
            r.Parent = new int[u.Count];
            r.Left = new int[u.Count];
            r.Right = new int[u.Count];
            r.EdgeLength = new double[u.Count];
            for (int i = 0; i < u.Count; i++)
            {
                r.Parent[i] = u.Nb[i][0];
                r.Left[i] = u.Nb[i][1];
                r.Right[i] = u.Nb[i][2];
                r.EdgeLength[i] = u.Len[i][0];
            }

            return r;
        }

        private static (int Node1, int Node2, double Length1, double Length2) FindRoot(UnrootedTree u, MuscleRoot method)
        {
            int nodeCount = u.Count;
            if (nodeCount < 2)
                throw new InvalidOperationException("Root: don't support trees with < 2 edges");
            var memo = new Dictionary<(int, int), EdgeInfo>();

            EdgeInfo Info(int from, int to)
            {
                if (memo.TryGetValue((from, to), out var cached))
                    return cached;
                var ei = new EdgeInfo { Node1 = from, Node2 = to };
                if (u.IsLeaf(to))
                {
                    ei.MaxDistToLeaf = 0;
                    ei.TotalDistToLeaves = 0;
                    ei.MaxStep = -1;
                    ei.MostDistantLeaf = to;
                    ei.LeafCount = 1;
                                        memo[(from, to)] = ei;
                    return ei;
                }

                double maxDist = -1e29, total = 0.0;
                int leafCount = 0, mostDistantLeaf = -1, maxStep = -1;
                for (int sub = 0; sub < 3; sub++)
                {
                    int node3 = u.Nb[to][sub];
                    if (node3 == -1 || node3 == from)
                        continue;
                    var next = Info(to, node3);
                    leafCount += next.LeafCount;
                    double edge = u.EdgeLength(to, node3);
                    double totalDist = next.TotalDistToLeaves + next.LeafCount * edge;
                    total += totalDist;
                    double dd = next.MaxDistToLeaf + edge;
                    if (dd > maxDist)
                    {
                        maxDist = dd;
                        mostDistantLeaf = next.MostDistantLeaf;
                        maxStep = node3;
                    }
                }

                if (maxStep == -1 || mostDistantLeaf == -1 || leafCount == 0)
                    throw new InvalidOperationException("CalcInfo: internal error 2");
                ei.MaxDistToLeaf = maxDist;
                ei.TotalDistToLeaves = total;
                ei.MaxStep = maxStep;
                ei.MostDistantLeaf = mostDistantLeaf;
                ei.LeafCount = leafCount;
                                memo[(from, to)] = ei;
                return ei;
            }

            const double veryLarge = 1e20;
            if (method == MuscleRoot.MidLongestSpan)
            {
                int leaf1 = -1;
                double maxDist = -veryLarge;
                for (int node = 0; node < nodeCount; node++)
                {
                    if (!u.IsLeaf(node))
                        continue;
                    int node2 = u.Nb[node][0];
                    if (node2 == -1)
                        throw new InvalidOperationException("RootByMidLongestSpan: internal error 0");
                    double edge = u.EdgeLength(node, node2);
                    var ei = Info(node, node2);
                    double span = edge + ei.MaxDistToLeaf;
                    if (span > maxDist)
                    {
                        maxDist = span;
                        leaf1 = node;
                    }
                }

                if (leaf1 == -1)
                    throw new InvalidOperationException("RootByMidLongestSpan: internal error 3");
                double treeHeight = maxDist / 2.0;
                int n1 = leaf1, n2 = u.Nb[leaf1][0];
                double accum = 0;
                while (true)
                {
                    double edge = u.EdgeLength(n1, n2);
                    if (accum + edge >= treeHeight)
                    {
                        double l1 = treeHeight - accum;
                        return (n1, n2, l1, edge - l1);
                    }

                    if (u.IsLeaf(n2))
                        throw new InvalidOperationException("RootByMidLongestSpan: internal error 4");
                    accum += edge;
                    var ei = Info(n1, n2);
                    n1 = n2;
                    n2 = ei.MaxStep;
                }
            }

            // RootByMinAvgLeafDist
            int leafCountTotal = Enumerable.Range(0, nodeCount).Count(u.IsLeaf);
            int best1 = -1, best2 = -1;
            double minHeight = veryLarge, bestLen1 = veryLarge, bestLen2 = veryLarge;
            for (int node = 0; node < nodeCount; node++)
            {
                int nbCount = u.NeighborCount(node);
                for (int sub = 0; sub < nbCount; sub++)
                {
                    int nb = u.Nb[node][sub];
                    if (nb < node)
                        continue;
                    var ei = Info(node, nb);
                    var rev = Info(nb, node);
                    if (leafCountTotal != ei.LeafCount + rev.LeafCount)
                        throw new InvalidOperationException("RootByMinAvgLeafDist, internal error 4");
                    double edge = u.EdgeLength(node, nb);
                    double adl1 = ei.TotalDistToLeaves / ei.LeafCount;
                    double adr2 = rev.TotalDistToLeaves / rev.LeafCount;
                    double x = (adr2 - adl1 + edge) / 2.0;
                    if (x >= 0 && x <= edge)
                    {
                        double length1 = x, length2 = edge - x;
                        double height1 = ei.MaxDistToLeaf + length1;
                        double height2 = rev.MaxDistToLeaf + length2;
                        double height = height1 >= height2 ? height1 : height2;
                        if (height < minHeight)
                        {
                            best1 = node;
                            best2 = nb;
                            bestLen1 = length1;
                            bestLen2 = length2;
                            minHeight = height;
                        }
                    }
                }
            }

            if (best1 == -1 || best2 == -1)
                throw new InvalidOperationException("RootByMinAvgLeafDist, internal error 6");
            return (best1, best2, bestLen1, bestLen2);
        }
    }
}
