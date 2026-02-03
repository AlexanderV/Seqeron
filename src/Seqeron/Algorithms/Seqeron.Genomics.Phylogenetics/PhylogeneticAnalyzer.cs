using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Phylogenetics;

/// <summary>
/// Provides phylogenetic tree construction and analysis algorithms.
/// Supports UPGMA, Neighbor-Joining, and distance matrix calculations.
/// </summary>
public static class PhylogeneticAnalyzer
{
    /// <summary>
    /// Represents a node in a phylogenetic tree.
    /// </summary>
    public class PhyloNode
    {
        public string Name { get; set; } = "";
        public double BranchLength { get; set; }
        public PhyloNode? Left { get; set; }
        public PhyloNode? Right { get; set; }
        public bool IsLeaf => Left == null && Right == null;
        public List<string> Taxa { get; set; } = new();

        public PhyloNode() { }

        public PhyloNode(string name)
        {
            Name = name;
            Taxa = new List<string> { name };
        }
    }

    /// <summary>
    /// Result of phylogenetic tree construction.
    /// </summary>
    public readonly record struct PhylogeneticTree(
        PhyloNode Root,
        IReadOnlyList<string> Taxa,
        double[,] DistanceMatrix,
        string Method);

    /// <summary>
    /// Distance calculation methods.
    /// </summary>
    public enum DistanceMethod
    {
        /// <summary>Proportion of differing sites (p-distance).</summary>
        PDistance,
        /// <summary>Jukes-Cantor corrected distance.</summary>
        JukesCantor,
        /// <summary>Kimura 2-parameter distance.</summary>
        Kimura2Parameter,
        /// <summary>Number of differing positions (raw count).</summary>
        Hamming
    }

    /// <summary>
    /// Tree construction methods.
    /// </summary>
    public enum TreeMethod
    {
        /// <summary>Unweighted Pair Group Method with Arithmetic Mean.</summary>
        UPGMA,
        /// <summary>Neighbor-Joining algorithm.</summary>
        NeighborJoining
    }

    /// <summary>
    /// Builds a phylogenetic tree from aligned sequences.
    /// </summary>
    /// <param name="sequences">Named aligned sequences (must be same length).</param>
    /// <param name="distanceMethod">Method for calculating distances.</param>
    /// <param name="treeMethod">Method for tree construction.</param>
    /// <returns>The constructed phylogenetic tree.</returns>
    public static PhylogeneticTree BuildTree(
        IReadOnlyDictionary<string, string> sequences,
        DistanceMethod distanceMethod = DistanceMethod.JukesCantor,
        TreeMethod treeMethod = TreeMethod.UPGMA)
    {
        if (sequences == null || sequences.Count < 2)
            throw new ArgumentException("At least 2 sequences required.", nameof(sequences));

        var taxa = sequences.Keys.ToList();
        var seqs = sequences.Values.ToList();

        // Validate alignment
        int length = seqs[0].Length;
        if (seqs.Any(s => s.Length != length))
            throw new ArgumentException("All sequences must have the same length (aligned).");

        // Calculate distance matrix
        var distMatrix = CalculateDistanceMatrix(seqs, distanceMethod);

        // Build tree
        PhyloNode root = treeMethod switch
        {
            TreeMethod.UPGMA => BuildUPGMA(taxa, distMatrix),
            TreeMethod.NeighborJoining => BuildNeighborJoining(taxa, distMatrix),
            _ => BuildUPGMA(taxa, distMatrix)
        };

        return new PhylogeneticTree(root, taxa, distMatrix, treeMethod.ToString());
    }

    /// <summary>
    /// Calculates a distance matrix for aligned sequences.
    /// </summary>
    public static double[,] CalculateDistanceMatrix(
        IReadOnlyList<string> alignedSequences,
        DistanceMethod method = DistanceMethod.JukesCantor)
    {
        int n = alignedSequences.Count;
        var matrix = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double dist = CalculatePairwiseDistance(
                    alignedSequences[i], alignedSequences[j], method);
                matrix[i, j] = dist;
                matrix[j, i] = dist;
            }
        }

        return matrix;
    }

    /// <summary>
    /// Calculates pairwise distance between two aligned sequences.
    /// </summary>
    public static double CalculatePairwiseDistance(
        string seq1, string seq2, DistanceMethod method = DistanceMethod.JukesCantor)
    {
        if (seq1.Length != seq2.Length)
            throw new ArgumentException("Sequences must have the same length.");

        int differences = 0;
        int transitions = 0;
        int transversions = 0;
        int comparableSites = 0;

        for (int i = 0; i < seq1.Length; i++)
        {
            char c1 = char.ToUpperInvariant(seq1[i]);
            char c2 = char.ToUpperInvariant(seq2[i]);

            // Skip gaps
            if (c1 == '-' || c2 == '-') continue;
            comparableSites++;

            if (c1 != c2)
            {
                differences++;
                if (IsTransition(c1, c2))
                    transitions++;
                else
                    transversions++;
            }
        }

        if (comparableSites == 0) return 0;

        double p = (double)differences / comparableSites;

        return method switch
        {
            DistanceMethod.Hamming => differences,
            DistanceMethod.PDistance => p,
            DistanceMethod.JukesCantor => JukesCantorDistance(p),
            DistanceMethod.Kimura2Parameter => Kimura2ParameterDistance(
                (double)transitions / comparableSites,
                (double)transversions / comparableSites),
            _ => p
        };
    }

    private static bool IsTransition(char c1, char c2)
    {
        // Purines: A, G; Pyrimidines: C, T
        bool bothPurines = (c1 == 'A' || c1 == 'G') && (c2 == 'A' || c2 == 'G');
        bool bothPyrimidines = (c1 == 'C' || c1 == 'T') && (c2 == 'C' || c2 == 'T');
        return bothPurines || bothPyrimidines;
    }

    private static double JukesCantorDistance(double p)
    {
        // JC69 correction: d = -3/4 * ln(1 - 4p/3)
        double arg = 1 - (4 * p / 3);
        if (arg <= 0) return double.PositiveInfinity;
        return -0.75 * Math.Log(arg);
    }

    private static double Kimura2ParameterDistance(double s, double v)
    {
        // K80: d = -0.5 * ln((1 - 2S - V) * sqrt(1 - 2V))
        double arg1 = 1 - 2 * s - v;
        double arg2 = 1 - 2 * v;
        if (arg1 <= 0 || arg2 <= 0) return double.PositiveInfinity;
        return -0.5 * Math.Log(arg1 * Math.Sqrt(arg2));
    }

    /// <summary>
    /// Builds a tree using UPGMA (Unweighted Pair Group Method with Arithmetic Mean).
    /// </summary>
    private static PhyloNode BuildUPGMA(List<string> taxa, double[,] distMatrix)
    {
        int n = taxa.Count;

        // Use a dictionary to track nodes and their sizes by original index
        var nodes = new Dictionary<int, PhyloNode>();
        var clusterSizes = new Dictionary<int, int>();

        // Initialize with leaf nodes
        for (int i = 0; i < n; i++)
        {
            nodes[i] = new PhyloNode(taxa[i]);
            clusterSizes[i] = 1;
        }

        // Working distance matrix as dictionary for sparse access
        var dist = new Dictionary<(int, int), double>();
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double d = distMatrix[i, j];
                dist[(i, j)] = d;
                dist[(j, i)] = d;
            }
        }

        var active = new HashSet<int>(Enumerable.Range(0, n));

        // Merge clusters until only one remains
        while (active.Count > 1)
        {
            // Find minimum distance pair
            double minDist = double.MaxValue;
            int minI = -1, minJ = -1;

            var activeList = active.ToList();
            for (int ii = 0; ii < activeList.Count; ii++)
            {
                for (int jj = ii + 1; jj < activeList.Count; jj++)
                {
                    int i = activeList[ii];
                    int j = activeList[jj];
                    var key = i < j ? (i, j) : (j, i);
                    if (dist.TryGetValue(key, out double d) && d < minDist)
                    {
                        minDist = d;
                        minI = i;
                        minJ = j;
                    }
                    else if (!dist.ContainsKey(key) && 0 < minDist)
                    {
                        // Zero distance (identical sequences)
                        minDist = 0;
                        minI = i;
                        minJ = j;
                    }
                }
            }

            if (minI == -1 || minJ == -1)
            {
                // Fallback: just pick first two
                minI = activeList[0];
                minJ = activeList[1];
                minDist = 0;
            }

            // Create new node
            var newNode = new PhyloNode
            {
                Name = $"({nodes[minI].Name},{nodes[minJ].Name})",
                Left = nodes[minI],
                Right = nodes[minJ],
                Taxa = nodes[minI].Taxa.Concat(nodes[minJ].Taxa).ToList()
            };

            // Set branch lengths (distance to midpoint)
            double height = minDist / 2;
            nodes[minI].BranchLength = Math.Max(0, height);
            nodes[minJ].BranchLength = Math.Max(0, height);

            // Update distances using UPGMA formula
            int newSize = clusterSizes[minI] + clusterSizes[minJ];
            foreach (int k in active)
            {
                if (k != minI && k != minJ)
                {
                    var keyIK = minI < k ? (minI, k) : (k, minI);
                    var keyJK = minJ < k ? (minJ, k) : (k, minJ);
                    double dIK = dist.GetValueOrDefault(keyIK, 0);
                    double dJK = dist.GetValueOrDefault(keyJK, 0);
                    double newDist = (dIK * clusterSizes[minI] + dJK * clusterSizes[minJ]) / newSize;
                    var newKey = minI < k ? (minI, k) : (k, minI);
                    dist[newKey] = newDist;
                    dist[(k, minI)] = newDist;
                    dist[(minI, k)] = newDist;
                }
            }

            // Replace minI with new cluster, remove minJ
            nodes[minI] = newNode;
            clusterSizes[minI] = newSize;
            active.Remove(minJ);
        }

        return nodes[active.First()];
    }

    /// <summary>
    /// Builds a tree using Neighbor-Joining algorithm.
    /// </summary>
    private static PhyloNode BuildNeighborJoining(List<string> taxa, double[,] distMatrix)
    {
        int n = taxa.Count;
        var nodes = new List<PhyloNode>();

        // Initialize with leaf nodes
        for (int i = 0; i < n; i++)
        {
            nodes.Add(new PhyloNode(taxa[i]));
        }

        // Working distance matrix
        var dist = new double[n, n];
        Array.Copy(distMatrix, dist, distMatrix.Length);
        var active = Enumerable.Range(0, n).ToList();

        while (active.Count > 2)
        {
            int m = active.Count;

            // Calculate r values (sum of distances)
            var r = new double[n];
            foreach (int i in active)
            {
                r[i] = 0;
                foreach (int j in active)
                {
                    if (i != j) r[i] += dist[i, j];
                }
            }

            // Find minimum Q value
            double minQ = double.MaxValue;
            int minI = -1, minJ = -1;

            for (int ii = 0; ii < m; ii++)
            {
                for (int jj = ii + 1; jj < m; jj++)
                {
                    int i = active[ii];
                    int j = active[jj];
                    double q = (m - 2) * dist[i, j] - r[i] - r[j];
                    if (q < minQ)
                    {
                        minQ = q;
                        minI = i;
                        minJ = j;
                    }
                }
            }

            // Calculate branch lengths
            double distIJ = dist[minI, minJ];
            double branchI = (distIJ / 2) + (r[minI] - r[minJ]) / (2 * (m - 2));
            double branchJ = distIJ - branchI;

            // Create new node
            var newNode = new PhyloNode
            {
                Name = $"({nodes[minI].Name},{nodes[minJ].Name})",
                Left = nodes[minI],
                Right = nodes[minJ],
                Taxa = nodes[minI].Taxa.Concat(nodes[minJ].Taxa).ToList()
            };

            nodes[minI].BranchLength = Math.Max(0, branchI);
            nodes[minJ].BranchLength = Math.Max(0, branchJ);

            // Update distances
            foreach (int k in active)
            {
                if (k != minI && k != minJ)
                {
                    double newDist = (dist[minI, k] + dist[minJ, k] - dist[minI, minJ]) / 2;
                    dist[minI, k] = newDist;
                    dist[k, minI] = newDist;
                }
            }

            nodes[minI] = newNode;
            active.Remove(minJ);
        }

        // Join last two nodes
        if (active.Count == 2)
        {
            int i = active[0];
            int j = active[1];
            double branchLen = dist[i, j] / 2;

            var root = new PhyloNode
            {
                Name = $"({nodes[i].Name},{nodes[j].Name})",
                Left = nodes[i],
                Right = nodes[j],
                Taxa = nodes[i].Taxa.Concat(nodes[j].Taxa).ToList()
            };

            nodes[i].BranchLength = branchLen;
            nodes[j].BranchLength = branchLen;

            return root;
        }

        return nodes[active[0]];
    }

    /// <summary>
    /// Converts a tree to Newick format.
    /// </summary>
    public static string ToNewick(PhyloNode node, bool includeBranchLengths = true)
    {
        if (node == null) return "";

        var sb = new StringBuilder();
        ToNewickRecursive(node, sb, includeBranchLengths, isRoot: true);
        sb.Append(';');
        return sb.ToString();
    }

    private static void ToNewickRecursive(PhyloNode node, StringBuilder sb, bool includeBranchLengths, bool isRoot)
    {
        if (node.IsLeaf)
        {
            sb.Append(node.Name);
        }
        else
        {
            sb.Append('(');
            if (node.Left != null)
            {
                ToNewickRecursive(node.Left, sb, includeBranchLengths, isRoot: false);
                if (includeBranchLengths)
                    sb.Append($":{node.Left.BranchLength:F4}");
            }
            sb.Append(',');
            if (node.Right != null)
            {
                ToNewickRecursive(node.Right, sb, includeBranchLengths, isRoot: false);
                if (includeBranchLengths)
                    sb.Append($":{node.Right.BranchLength:F4}");
            }
            sb.Append(')');
        }
    }

    /// <summary>
    /// Parses a Newick format tree string.
    /// </summary>
    public static PhyloNode ParseNewick(string newick)
    {
        if (string.IsNullOrWhiteSpace(newick))
            throw new ArgumentException("Newick string is empty.");

        newick = newick.Trim();
        if (newick.EndsWith(";"))
            newick = newick[..^1];

        int pos = 0;
        return ParseNewickRecursive(newick, ref pos);
    }

    private static PhyloNode ParseNewickRecursive(string newick, ref int pos)
    {
        var node = new PhyloNode();

        if (pos < newick.Length && newick[pos] == '(')
        {
            pos++; // skip '('

            // Parse left child
            node.Left = ParseNewickRecursive(newick, ref pos);
            node.Taxa.AddRange(node.Left.Taxa);

            // Parse branch length for left child
            if (pos < newick.Length && newick[pos] == ':')
            {
                pos++;
                node.Left.BranchLength = ParseNumber(newick, ref pos);
            }

            // Skip comma
            if (pos < newick.Length && newick[pos] == ',')
                pos++;

            // Parse right child
            node.Right = ParseNewickRecursive(newick, ref pos);
            node.Taxa.AddRange(node.Right.Taxa);

            // Parse branch length for right child
            if (pos < newick.Length && newick[pos] == ':')
            {
                pos++;
                node.Right.BranchLength = ParseNumber(newick, ref pos);
            }

            // Skip ')'
            if (pos < newick.Length && newick[pos] == ')')
                pos++;

            // Parse internal node name (if any)
            if (pos < newick.Length && newick[pos] != ':' && newick[pos] != ',' &&
                newick[pos] != ')' && newick[pos] != ';')
            {
                node.Name = ParseLabel(newick, ref pos);
            }
        }
        else
        {
            // Leaf node
            node.Name = ParseLabel(newick, ref pos);
            node.Taxa.Add(node.Name);
        }

        return node;
    }

    private static string ParseLabel(string newick, ref int pos)
    {
        var sb = new StringBuilder();
        while (pos < newick.Length &&
               newick[pos] != ':' && newick[pos] != ',' &&
               newick[pos] != ')' && newick[pos] != '(' &&
               newick[pos] != ';')
        {
            sb.Append(newick[pos]);
            pos++;
        }
        return sb.ToString();
    }

    private static double ParseNumber(string newick, ref int pos)
    {
        var sb = new StringBuilder();
        while (pos < newick.Length &&
               (char.IsDigit(newick[pos]) || newick[pos] == '.' ||
                newick[pos] == '-' || newick[pos] == 'e' || newick[pos] == 'E'))
        {
            sb.Append(newick[pos]);
            pos++;
        }
        return double.TryParse(sb.ToString(), out double val) ? val : 0;
    }

    /// <summary>
    /// Gets all leaf nodes (taxa) from the tree.
    /// </summary>
    public static IEnumerable<PhyloNode> GetLeaves(PhyloNode root)
    {
        if (root == null) yield break;

        if (root.IsLeaf)
        {
            yield return root;
        }
        else
        {
            if (root.Left != null)
                foreach (var leaf in GetLeaves(root.Left))
                    yield return leaf;
            if (root.Right != null)
                foreach (var leaf in GetLeaves(root.Right))
                    yield return leaf;
        }
    }

    /// <summary>
    /// Calculates the total tree length (sum of all branch lengths).
    /// </summary>
    public static double CalculateTreeLength(PhyloNode root)
    {
        if (root == null) return 0;

        double length = root.BranchLength;
        if (root.Left != null) length += CalculateTreeLength(root.Left);
        if (root.Right != null) length += CalculateTreeLength(root.Right);

        return length;
    }

    /// <summary>
    /// Gets the depth (height) of the tree.
    /// </summary>
    public static int GetTreeDepth(PhyloNode root)
    {
        if (root == null || root.IsLeaf) return 0;

        int leftDepth = root.Left != null ? GetTreeDepth(root.Left) : 0;
        int rightDepth = root.Right != null ? GetTreeDepth(root.Right) : 0;

        return 1 + Math.Max(leftDepth, rightDepth);
    }

    /// <summary>
    /// Calculates Robinson-Foulds distance between two trees.
    /// </summary>
    public static int RobinsonFouldsDistance(PhyloNode tree1, PhyloNode tree2)
    {
        var splits1 = GetSplits(tree1);
        var splits2 = GetSplits(tree2);

        int symmetricDiff = splits1.Except(splits2).Count() + splits2.Except(splits1).Count();
        return symmetricDiff;
    }

    private static HashSet<string> GetSplits(PhyloNode root)
    {
        var splits = new HashSet<string>();
        var allTaxa = GetLeaves(root).Select(l => l.Name).OrderBy(n => n).ToList();

        CollectSplits(root, splits, allTaxa);
        return splits;
    }

    private static List<string> CollectSplits(PhyloNode node, HashSet<string> splits, List<string> allTaxa)
    {
        if (node == null) return new List<string>();

        if (node.IsLeaf)
        {
            return new List<string> { node.Name };
        }

        var leftTaxa = node.Left != null ? CollectSplits(node.Left, splits, allTaxa) : new List<string>();
        var rightTaxa = node.Right != null ? CollectSplits(node.Right, splits, allTaxa) : new List<string>();
        var subtreeTaxa = leftTaxa.Concat(rightTaxa).OrderBy(n => n).ToList();

        // Create split representation (smaller side)
        if (subtreeTaxa.Count > 0 && subtreeTaxa.Count < allTaxa.Count)
        {
            var complement = allTaxa.Except(subtreeTaxa).OrderBy(n => n).ToList();
            var smaller = subtreeTaxa.Count <= complement.Count ? subtreeTaxa : complement;
            splits.Add(string.Join("|", smaller));
        }

        return subtreeTaxa;
    }

    /// <summary>
    /// Finds the most recent common ancestor of two taxa.
    /// </summary>
    public static PhyloNode? FindMRCA(PhyloNode root, string taxon1, string taxon2)
    {
        if (root == null) return null;

        if (root.IsLeaf)
        {
            return root.Name == taxon1 || root.Name == taxon2 ? root : null;
        }

        var leftResult = FindMRCA(root.Left!, taxon1, taxon2);
        var rightResult = FindMRCA(root.Right!, taxon1, taxon2);

        if (leftResult != null && rightResult != null)
            return root;

        return leftResult ?? rightResult;
    }

    /// <summary>
    /// Calculates the patristic distance (tree path length) between two taxa.
    /// </summary>
    public static double PatristicDistance(PhyloNode root, string taxon1, string taxon2)
    {
        var mrca = FindMRCA(root, taxon1, taxon2);
        if (mrca == null) return double.NaN;

        double dist1 = DistanceToTaxon(mrca, taxon1);
        double dist2 = DistanceToTaxon(mrca, taxon2);

        return dist1 + dist2;
    }

    private static double DistanceToTaxon(PhyloNode node, string taxon)
    {
        if (node == null) return double.NaN;

        if (node.IsLeaf)
            return node.Name == taxon ? 0 : double.NaN;

        double leftDist = DistanceToTaxon(node.Left!, taxon);
        if (!double.IsNaN(leftDist))
            return leftDist + (node.Left?.BranchLength ?? 0);

        double rightDist = DistanceToTaxon(node.Right!, taxon);
        if (!double.IsNaN(rightDist))
            return rightDist + (node.Right?.BranchLength ?? 0);

        return double.NaN;
    }

    /// <summary>
    /// Bootstrap analysis - builds multiple trees from resampled alignments.
    /// </summary>
    public static IReadOnlyDictionary<string, double> Bootstrap(
        IReadOnlyDictionary<string, string> sequences,
        int replicates = 100,
        DistanceMethod distanceMethod = DistanceMethod.JukesCantor,
        TreeMethod treeMethod = TreeMethod.UPGMA)
    {
        var taxa = sequences.Keys.ToList();
        var seqs = sequences.Values.ToList();
        int alignmentLength = seqs[0].Length;

        // Build reference tree
        var refTree = BuildTree(sequences, distanceMethod, treeMethod);
        var refSplits = GetSplits(refTree.Root);

        // Count support for each split
        var supportCounts = new Dictionary<string, int>();
        foreach (var split in refSplits)
            supportCounts[split] = 0;

        var random = new Random(42);

        for (int rep = 0; rep < replicates; rep++)
        {
            // Resample columns with replacement
            var resampledSeqs = new Dictionary<string, string>();
            var columns = new int[alignmentLength];
            for (int i = 0; i < alignmentLength; i++)
                columns[i] = random.Next(alignmentLength);

            for (int t = 0; t < taxa.Count; t++)
            {
                var sb = new StringBuilder();
                foreach (int col in columns)
                    sb.Append(seqs[t][col]);
                resampledSeqs[taxa[t]] = sb.ToString();
            }

            // Build tree from resampled data
            var bootTree = BuildTree(resampledSeqs, distanceMethod, treeMethod);
            var bootSplits = GetSplits(bootTree.Root);

            // Count matching splits
            foreach (var split in refSplits)
            {
                if (bootSplits.Contains(split))
                    supportCounts[split]++;
            }
        }

        // Convert to proportions
        return supportCounts.ToDictionary(
            kvp => kvp.Key,
            kvp => (double)kvp.Value / replicates);
    }
}
