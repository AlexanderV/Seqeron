using System;
using System.Collections.Generic;
using System.Globalization;
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
    /// <remarks>
    /// <para>
    /// The tree model is <b>N-ary (multifurcating)</b>: an internal node holds an ordered
    /// collection of <see cref="Children"/> rather than a fixed left/right pair. A node with
    /// 0 children is a leaf (taxon); 2 children is the common bifurcation; ≥3 children is a
    /// <i>multifurcation</i> / <i>polytomy</i> — an unresolved or hard polytomy
    /// (Newick format; a node "is multifurcating if it has three or more immediate descendant lineages").
    /// This enables genuine Newick multifurcations <c>(A,B,C);</c> and the natural unrooted
    /// trifurcation of Neighbor-Joining (Saitou &amp; Nei 1987).
    /// </para>
    /// <para>
    /// <see cref="Left"/> and <see cref="Right"/> are convenience accessors over the first two
    /// children for the common binary case; reading them on a node with fewer children yields
    /// <c>null</c>, and assigning them rewrites <see cref="Children"/>[0]/[1] (growing the list
    /// as needed). Code that must handle multifurcations should traverse <see cref="Children"/>.
    /// </para>
    /// </remarks>
    public class PhyloNode
    {
        public string Name { get; set; } = "";
        public double BranchLength { get; set; }

        /// <summary>
        /// The ordered child subtrees. Empty for a leaf, 2 for a bifurcation, ≥3 for a multifurcation.
        /// </summary>
        public List<PhyloNode> Children { get; set; } = new();

        /// <summary>
        /// Convenience accessor for the first child (binary-tree compatibility).
        /// Returns <c>null</c> when there are no children. Setting to a non-null value places
        /// the node at <see cref="Children"/>[0]; setting to <c>null</c> removes the first child.
        /// </summary>
        public PhyloNode? Left
        {
            get => Children.Count > 0 ? Children[0] : null;
            set => SetChildAt(0, value);
        }

        /// <summary>
        /// Convenience accessor for the second child (binary-tree compatibility).
        /// Returns <c>null</c> when there are fewer than two children. Setting to a non-null value
        /// places the node at <see cref="Children"/>[1]; setting to <c>null</c> removes the second child.
        /// </summary>
        public PhyloNode? Right
        {
            get => Children.Count > 1 ? Children[1] : null;
            set => SetChildAt(1, value);
        }

        private void SetChildAt(int index, PhyloNode? value)
        {
            if (value == null)
            {
                if (Children.Count > index)
                    Children.RemoveAt(index);
                return;
            }
            while (Children.Count <= index)
                Children.Add(value);
            Children[index] = value;
        }

        /// <summary>A leaf is a node with no children (a terminal taxon/OTU).</summary>
        public bool IsLeaf => Children.Count == 0;

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
    /// Builds a phylogenetic tree from a pre-computed distance matrix.
    /// Allows direct testing against known reference matrices (e.g., Wikipedia examples).
    /// </summary>
    /// <param name="taxa">Taxon names in the same order as the distance matrix.</param>
    /// <param name="distanceMatrix">Pre-computed symmetric distance matrix.</param>
    /// <param name="treeMethod">Method for tree construction.</param>
    /// <returns>The constructed phylogenetic tree.</returns>
    public static PhylogeneticTree BuildTreeFromMatrix(
        IReadOnlyList<string> taxa,
        double[,] distanceMatrix,
        TreeMethod treeMethod = TreeMethod.UPGMA)
    {
        if (taxa == null || taxa.Count < 2)
            throw new ArgumentException("At least 2 taxa required.", nameof(taxa));
        if (distanceMatrix == null)
            throw new ArgumentException("Distance matrix is required.", nameof(distanceMatrix));
        if (distanceMatrix.GetLength(0) != taxa.Count || distanceMatrix.GetLength(1) != taxa.Count)
            throw new ArgumentException("Distance matrix dimensions must match the number of taxa.");

        PhyloNode root = treeMethod switch
        {
            TreeMethod.UPGMA => BuildUPGMA(taxa.ToList(), distanceMatrix),
            TreeMethod.NeighborJoining => BuildNeighborJoining(taxa.ToList(), distanceMatrix),
            _ => BuildUPGMA(taxa.ToList(), distanceMatrix)
        };

        return new PhylogeneticTree(root, taxa, distanceMatrix, treeMethod.ToString());
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
        if (seq1 == null) throw new ArgumentNullException(nameof(seq1));
        if (seq2 == null) throw new ArgumentNullException(nameof(seq2));
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

            // Skip gaps and ambiguous bases (only compare A, C, G, T)
            if (c1 == '-' || c2 == '-') continue;
            if (!IsStandardBase(c1) || !IsStandardBase(c2)) continue;
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

    private static bool IsStandardBase(char c) =>
        c == 'A' || c == 'C' || c == 'G' || c == 'T';

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
    /// Branch lengths are computed as incremental heights per the UPGMA definition:
    /// height(new_cluster) = d(i,j)/2; branch_length(child) = height(new) - height(child).
    /// Source: Wikipedia UPGMA, Sokal &amp; Michener (1958).
    /// </summary>
    private static PhyloNode BuildUPGMA(List<string> taxa, double[,] distMatrix)
    {
        int n = taxa.Count;

        var nodes = new Dictionary<int, PhyloNode>();
        var clusterSizes = new Dictionary<int, int>();
        var clusterHeights = new Dictionary<int, double>();

        // Initialize with leaf nodes (height 0)
        for (int i = 0; i < n; i++)
        {
            nodes[i] = new PhyloNode(taxa[i]);
            clusterSizes[i] = 1;
            clusterHeights[i] = 0.0;
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
                    int a = activeList[ii];
                    int b = activeList[jj];
                    var key = a < b ? (a, b) : (b, a);
                    if (dist.TryGetValue(key, out double d) && d < minDist)
                    {
                        minDist = d;
                        minI = a;
                        minJ = b;
                    }
                    else if (!dist.ContainsKey(key) && 0 < minDist)
                    {
                        minDist = 0;
                        minI = a;
                        minJ = b;
                    }
                }
            }

            if (minI == -1 || minJ == -1)
            {
                minI = activeList[0];
                minJ = activeList[1];
                minDist = 0;
            }

            // New cluster height = d(i,j) / 2 (ultrametric property)
            double newHeight = minDist / 2;

            // Incremental branch lengths: height(new) - height(child)
            nodes[minI].BranchLength = Math.Max(0, newHeight - clusterHeights[minI]);
            nodes[minJ].BranchLength = Math.Max(0, newHeight - clusterHeights[minJ]);

            var newNode = new PhyloNode
            {
                Name = $"({nodes[minI].Name},{nodes[minJ].Name})",
                Left = nodes[minI],
                Right = nodes[minJ],
                Taxa = nodes[minI].Taxa.Concat(nodes[minJ].Taxa).ToList()
            };

            // Update distances using UPGMA formula (proportional averaging)
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
            clusterHeights[minI] = newHeight;
            active.Remove(minJ);
        }

        return nodes[active.First()];
    }

    /// <summary>
    /// Builds a tree using Neighbor-Joining algorithm.
    /// </summary>
    /// <remarks>
    /// Neighbor-Joining infers an <b>unrooted</b> tree, so the algorithm stops when three OTUs
    /// remain and connects them to a single central node — the characteristic
    /// <i>trifurcation</i> at the centre of the unrooted tree (Saitou &amp; Nei 1987; the
    /// three final branch lengths are obtained from the three pairwise distances via the
    /// standard additive system <c>L_ix = (d_ij + d_ik − d_jk)/2</c>). Under the N-ary
    /// <see cref="PhyloNode"/> model this final node carries three children rather than two,
    /// so an NJ tree round-trips through Newick as a genuine trifurcation.
    /// </remarks>
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

        // Stop at three remaining OTUs: NJ yields an unrooted tree whose centre is a
        // trifurcation of the last three nodes (Saitou & Nei 1987).
        while (active.Count > 3)
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

            // NJ may produce negative branch lengths (Wikipedia: INV-N02)
            nodes[minI].BranchLength = branchI;
            nodes[minJ].BranchLength = branchJ;

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

        // Final step: connect the last three OTUs to one central node (NJ trifurcation).
        // The three branch lengths solve the additive 3-taxon system (Saitou & Nei 1987):
        //   L_i = (d_ij + d_ik − d_jk)/2, and symmetrically for j and k.
        if (active.Count == 3)
        {
            int i = active[0];
            int j = active[1];
            int k = active[2];

            double dij = dist[i, j];
            double dik = dist[i, k];
            double djk = dist[j, k];

            nodes[i].BranchLength = (dij + dik - djk) / 2;
            nodes[j].BranchLength = (dij + djk - dik) / 2;
            nodes[k].BranchLength = (dik + djk - dij) / 2;

            var root = new PhyloNode
            {
                Name = $"({nodes[i].Name},{nodes[j].Name},{nodes[k].Name})",
                Children = new List<PhyloNode> { nodes[i], nodes[j], nodes[k] },
                Taxa = nodes[i].Taxa.Concat(nodes[j].Taxa).Concat(nodes[k].Taxa).ToList()
            };

            return root;
        }

        // Two-taxon degenerate input: a single edge between the two leaves.
        if (active.Count == 2)
        {
            int i = active[0];
            int j = active[1];
            double branchLen = dist[i, j] / 2;

            var root = new PhyloNode
            {
                Name = $"({nodes[i].Name},{nodes[j].Name})",
                Children = new List<PhyloNode> { nodes[i], nodes[j] },
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
            // N-ary: emit every child (≥2), comma-separated, per the Newick grammar
            // Internal → "(" BranchSet ")" Name where BranchSet is a comma-separated
            // list of Branch entries. A node with ≥3 children is a multifurcation.
            sb.Append('(');
            for (int c = 0; c < node.Children.Count; c++)
            {
                if (c > 0) sb.Append(',');
                var child = node.Children[c];
                ToNewickRecursive(child, sb, includeBranchLengths, isRoot: false);
                if (includeBranchLengths)
                {
                    sb.Append(':');
                    sb.Append(child.BranchLength.ToString("F4", CultureInfo.InvariantCulture));
                }
            }
            sb.Append(')');

            // Newick grammar: Internal → "(" BranchSet ")" Name
            // Emit internal node name if it is a valid unquoted Newick label.
            // Names containing Newick metacharacters (e.g. auto-generated UPGMA/NJ names)
            // are omitted per Olsen spec: unquoted labels may not contain
            // blanks, parentheses, square brackets, single quotes, colons, semicolons, or commas.
            if (IsValidUnquotedNewickLabel(node.Name))
                sb.Append(node.Name);
        }
    }

    /// <summary>
    /// Checks whether a label is a valid unquoted Newick label per the Olsen specification.
    /// Unquoted labels may not contain blanks, parentheses, square brackets,
    /// single quotes, colons, semicolons, or commas.
    /// Source: Olsen (1990), Wikipedia Newick format § Notes.
    /// </summary>
    private static bool IsValidUnquotedNewickLabel(string label)
    {
        if (string.IsNullOrEmpty(label)) return false;
        foreach (char c in label)
        {
            if (c == ' ' || c == '(' || c == ')' || c == '[' || c == ']' ||
                c == '\'' || c == ':' || c == ';' || c == ',')
                return false;
        }
        return true;
    }

    /// <summary>
    /// Maximum parenthesis-nesting depth accepted by <see cref="ParseNewick"/>.
    /// The parser is recursive descent, so each level of nested parentheses consumes one
    /// stack frame; an unbounded depth (e.g. a malicious 50&#160;000&#215;<c>(</c> input) would
    /// otherwise overflow the call stack and crash the process — an uncatchable
    /// <see cref="StackOverflowException"/> and a denial-of-service hazard. The recursion is
    /// therefore capped at this depth and a pathologically nested input is rejected with a
    /// <see cref="FormatException"/> rather than being allowed to overflow. The limit is far
    /// larger than any biologically meaningful tree depth (a binary tree of depth 1000 has up
    /// to 2^1000 leaves) yet stays comfortably within a default 1&#160;MB thread stack, on which
    /// the unguarded recursion was measured to overflow at roughly 3000 frames.
    /// </summary>
    public const int MaxParseDepth = 1000;

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
        var root = ParseNewickRecursive(newick, ref pos, depth: 0);

        // Grammar alternative: Tree → Branch ";" (Olsen: tree ==> descendant_list [root_label] [:branch_length] ;)
        // Handle optional root branch length after the main subtree.
        if (pos < newick.Length && newick[pos] == ':')
        {
            pos++;
            root.BranchLength = ParseNumber(newick, ref pos);
        }

        // Any remaining input (the terminal ';' and surrounding whitespace were already
        // stripped above) is unconsumed trailing garbage and indicates a malformed tree.
        // Throw rather than silently ignoring it.
        if (pos < newick.Length)
            throw new FormatException(
                $"Malformed Newick string: unexpected trailing input '{newick[pos..]}' at position {pos}.");

        return root;
    }

    private static PhyloNode ParseNewickRecursive(string newick, ref int pos, int depth)
    {
        // Depth guard: a recursive-descent parser consumes one stack frame per level of
        // parenthesis nesting, so an unbounded nesting depth would overflow the call stack
        // (an uncatchable StackOverflowException that terminates the process — a DoS hazard).
        // Reject pathologically nested input with a catchable FormatException instead.
        if (depth > MaxParseDepth)
            throw new FormatException(
                $"Malformed Newick string: parenthesis nesting exceeds the maximum supported " +
                $"depth of {MaxParseDepth} (possible malformed or malicious input).");

        var node = new PhyloNode();

        if (pos < newick.Length && newick[pos] == '(')
        {
            pos++; // skip '('

            // N-ary descendant list: parse one or more children separated by commas.
            // Newick grammar: BranchSet → Branch | Branch "," BranchSet. A node with three or
            // more children is a genuine multifurcation (polytomy) and is parsed faithfully.
            while (true)
            {
                var child = ParseNewickRecursive(newick, ref pos, depth + 1);

                // Parse branch length for this child
                if (pos < newick.Length && newick[pos] == ':')
                {
                    pos++;
                    child.BranchLength = ParseNumber(newick, ref pos);
                }

                node.Children.Add(child);
                node.Taxa.AddRange(child.Taxa);

                if (pos < newick.Length && newick[pos] == ',')
                {
                    pos++; // consume separator, parse the next child
                    continue;
                }
                break;
            }

            // Require the matching ')'. The Newick grammar (Olsen; Wikipedia) defines an
            // internal node as Internal → "(" BranchSet ")" Name: the closing parenthesis is
            // mandatory. An opened descendant list that is never closed (e.g. "(A,B" or
            // "((A,B);") is a malformed tree with unbalanced parentheses and must be rejected
            // rather than silently accepted as a (truncated/degenerate) tree.
            if (pos < newick.Length && newick[pos] == ')')
            {
                pos++;
            }
            else
            {
                throw new FormatException(
                    "Malformed Newick string: unbalanced parentheses — an opening '(' has no " +
                    $"matching ')' (unexpected end of descendant list at position {pos}).");
            }

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
                newick[pos] == '-' || newick[pos] == '+' ||
                newick[pos] == 'e' || newick[pos] == 'E'))
        {
            sb.Append(newick[pos]);
            pos++;
        }
        // Use InvariantCulture: Newick format always uses '.' as decimal separator
        // regardless of the system locale (Wikipedia grammar: Length → ":" number).
        return double.TryParse(sb.ToString(), NumberStyles.Float, CultureInfo.InvariantCulture, out double val) ? val : 0;
    }

    /// <summary>
    /// Gets all leaf (terminal) nodes of the tree in left-to-right (pre-order) traversal order.
    /// A leaf is a node with no children (<see cref="PhyloNode.IsLeaf"/>); in a phylogenetic tree
    /// the leaves are the taxa/operational taxonomic units.
    /// Source: a leaf is "a vertex with no children" (Tree (graph theory)); Biopython
    /// <c>Tree.get_terminals</c> returns "all of this tree's terminal (leaf) nodes".
    /// </summary>
    /// <param name="root">The root of the (sub)tree. A <c>null</c> root yields no leaves.</param>
    /// <returns>The leaf nodes; empty when <paramref name="root"/> is null.</returns>
    public static IEnumerable<PhyloNode> GetLeaves(PhyloNode root)
    {
        if (root == null) yield break;

        if (root.IsLeaf)
        {
            yield return root;
        }
        else
        {
            foreach (var child in root.Children)
                foreach (var leaf in GetLeaves(child))
                    yield return leaf;
        }
    }

    /// <summary>
    /// Calculates the total tree length: the sum of the branch lengths of every node in the
    /// (sub)tree rooted at <paramref name="root"/>. This is the quantity minimized by the
    /// minimum-evolution criterion (Rzhetsky &amp; Nei 1992).
    /// Source: DendroPy <c>Tree.length()</c> — "the sum of edge lengths of self"; Biopython
    /// <c>Tree.total_branch_length</c> — "the sum of all the branch lengths in this tree".
    /// </summary>
    /// <param name="root">The root of the (sub)tree.</param>
    /// <returns>
    /// The sum of all branch lengths; <c>0</c> for a null tree (no edges) and for a tree whose
    /// branch lengths are all zero. (DendroPy treats undefined edge lengths as 0; here the default
    /// <see cref="PhyloNode.BranchLength"/> is 0.)
    /// </returns>
    public static double CalculateTreeLength(PhyloNode root)
    {
        if (root == null) return 0;

        double length = root.BranchLength;
        foreach (var child in root.Children)
            length += CalculateTreeLength(child);

        return length;
    }

    /// <summary>Height of an empty tree (no vertices), by convention.</summary>
    /// <remarks>
    /// "Conventionally, an empty tree (a tree with no vertices, if such are allowed) has depth and
    /// height −1." (Tree (graph theory) / Tree (abstract data type)).
    /// </remarks>
    private const int EmptyTreeHeight = -1;

    /// <summary>
    /// Gets the height (topological depth) of the tree: the number of edges on the longest
    /// downward path from <paramref name="root"/> to a leaf. The height of the tree is the height
    /// of its root.
    /// Source: "The height of a node is the length of the longest downward path to a leaf from that
    /// node. The height of the root is the height of the tree." (Tree (abstract data type)).
    /// </summary>
    /// <param name="root">The root of the tree.</param>
    /// <returns>
    /// The height in edges: <c>0</c> for a single-node (leaf-only) tree — "a tree with only a single
    /// node ... has depth and height zero" — and <c>-1</c> for a null/empty tree by the cited
    /// convention.
    /// </returns>
    public static int GetTreeDepth(PhyloNode root)
    {
        // Empty tree (null): height -1 by convention; a single leaf node has height 0.
        if (root == null) return EmptyTreeHeight;
        if (root.IsLeaf) return 0;

        int maxChildDepth = 0;
        foreach (var child in root.Children)
            maxChildDepth = Math.Max(maxChildDepth, GetTreeDepth(child));

        return 1 + maxChildDepth;
    }

    /// <summary>
    /// Calculates Robinson-Foulds distance between two rooted trees.
    /// Uses clade (cluster) comparison: each internal node defines a clade
    /// (the set of taxa in its subtree). RF = |clades(T1) △ clades(T2)|.
    /// </summary>
    public static int RobinsonFouldsDistance(PhyloNode tree1, PhyloNode tree2)
    {
        var clades1 = GetClades(tree1);
        var clades2 = GetClades(tree2);

        int symmetricDiff = clades1.Except(clades2).Count() + clades2.Except(clades1).Count();
        return symmetricDiff;
    }

    private static HashSet<string> GetClades(PhyloNode root)
    {
        var clades = new HashSet<string>();
        int totalLeaves = GetLeaves(root).Count();

        CollectClades(root, clades, totalLeaves);
        return clades;
    }

    /// <summary>
    /// Calculates the <b>unrooted</b> Robinson–Foulds distance between two trees,
    /// as the symmetric difference of their sets of non-trivial <i>bipartitions</i>
    /// (splits) — the original Robinson &amp; Foulds (1981) metric.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is an additive alternative to <see cref="RobinsonFouldsDistance"/>, which
    /// compares rooted <i>clades</i>. Removing any internal edge of a tree partitions
    /// the leaf set into two complementary subsets; that unordered pair is a
    /// <i>bipartition</i> (split). The unrooted RF distance is
    /// <c>RF = |A\B| + |B\A|</c> over the bipartition sets A = splits(T1), B = splits(T2).
    /// </para>
    /// <para>
    /// <b>Root invariance.</b> A bipartition is an <i>unordered</i> {S, complement} pair,
    /// so it does not depend on where the tree is rooted. Two binary rooted trees that
    /// describe the same unrooted topology but are rooted on different edges therefore
    /// have unrooted RF = 0, even though their rooted-clade RF can be non-zero. This is the
    /// distinguishing behaviour of this metric.
    /// </para>
    /// <para>
    /// <b>Trivial splits are excluded.</b> Splits induced by terminal (leaf) edges —
    /// a single leaf versus all-but-one — are shared by every tree on the same taxa and
    /// carry no topological information, so only splits with at least 2 leaves on each side
    /// are counted. The binary <see cref="PhyloNode"/> model is read as unrooted by ignoring
    /// the root placement when forming bipartitions.
    /// </para>
    /// </remarks>
    /// <param name="tree1">First tree (root of a binary <see cref="PhyloNode"/> tree).</param>
    /// <param name="tree2">Second tree.</param>
    /// <returns>The unrooted RF distance (a non-negative, even integer).</returns>
    /// <exception cref="ArgumentNullException">Either tree is null.</exception>
    /// <exception cref="ArgumentException">
    /// Either tree has fewer than 3 leaves (no non-trivial split is possible), or the two
    /// trees do not have identical leaf sets (the RF metric is only defined on a common taxon set).
    /// </exception>
    public static int CalculateUnrootedRobinsonFoulds(PhyloNode tree1, PhyloNode tree2)
    {
        if (tree1 == null) throw new ArgumentNullException(nameof(tree1));
        if (tree2 == null) throw new ArgumentNullException(nameof(tree2));

        var leaves1 = GetLeaves(tree1).Select(l => l.Name).ToList();
        var leaves2 = GetLeaves(tree2).Select(l => l.Name).ToList();

        if (leaves1.Count < 3 || leaves2.Count < 3)
            throw new ArgumentException(
                "Unrooted Robinson–Foulds distance requires at least 3 leaves per tree " +
                "(no non-trivial bipartition exists below that).");

        var taxa1 = new HashSet<string>(leaves1);
        var taxa2 = new HashSet<string>(leaves2);
        if (!taxa1.SetEquals(taxa2))
            throw new ArgumentException(
                "Unrooted Robinson–Foulds distance is only defined for trees over the same leaf set.");

        var splits1 = GetBipartitions(tree1, taxa1);
        var splits2 = GetBipartitions(tree2, taxa2);

        return splits1.Except(splits2).Count() + splits2.Except(splits1).Count();
    }

    /// <summary>
    /// Calculates the <b>normalized</b> unrooted Robinson–Foulds distance:
    /// <c>RF / (2n − 6)</c>, where <c>n</c> is the number of leaves. The denominator is the
    /// maximum possible unrooted RF for two fully-resolved (binary) trees on <c>n</c> leaves
    /// (each has <c>n − 3</c> internal edges, so the symmetric difference is at most
    /// <c>2(n − 3) = 2n − 6</c>). The result is in [0, 1] for binary trees.
    /// </summary>
    /// <param name="tree1">First tree.</param>
    /// <param name="tree2">Second tree.</param>
    /// <returns>Normalized RF in [0, 1] for binary trees; 0 when n &lt; 4 (no internal split possible).</returns>
    /// <exception cref="ArgumentNullException">Either tree is null.</exception>
    /// <exception cref="ArgumentException">Fewer than 3 leaves, or mismatched leaf sets.</exception>
    public static double CalculateNormalizedUnrootedRobinsonFoulds(PhyloNode tree1, PhyloNode tree2)
    {
        int rf = CalculateUnrootedRobinsonFoulds(tree1, tree2);

        int n = GetLeaves(tree1).Count();
        int denominator = 2 * n - 6;
        // n == 3 ⇒ denominator 0: only one unrooted topology exists, RF is always 0.
        if (denominator <= 0) return 0.0;

        return (double)rf / denominator;
    }

    /// <summary>
    /// Collects the set of non-trivial bipartitions (splits) of a tree, read as unrooted.
    /// Each internal edge induces a split of the leaf set; the split is canonicalised to its
    /// smaller side (ties broken lexicographically) so that {S, complement} is represented once,
    /// independent of root placement.
    /// </summary>
    private static HashSet<string> GetBipartitions(PhyloNode root, HashSet<string> allTaxa)
    {
        var splits = new HashSet<string>();
        CollectBipartitions(root, splits, allTaxa);
        return splits;
    }

    private static List<string> CollectBipartitions(
        PhyloNode? node, HashSet<string> splits, HashSet<string> allTaxa)
    {
        if (node == null) return new List<string>();

        if (node.IsLeaf)
            return new List<string> { node.Name };

        // Each internal node corresponds to one edge (incident to its parent); the leaves
        // beneath it are one side of the induced split. A multifurcation simply produces one
        // such side (one candidate split) instead of the two an internal binary node would —
        // i.e. an unresolved node contributes fewer non-trivial bipartitions (Robinson–Foulds:
        // collapsing/contracting an internal edge removes its split from Σ(T)).
        var subtreeTaxa = new List<string>();
        foreach (var child in node.Children)
            subtreeTaxa.AddRange(CollectBipartitions(child, splits, allTaxa));

        int total = allTaxa.Count;
        int side = subtreeTaxa.Count;
        // Non-trivial split: at least 2 leaves on each side. A side of size 0/1 (leaf edge)
        // or total-1 (the complementary leaf edge) is trivial and shared by all trees.
        if (side >= 2 && total - side >= 2)
        {
            var key = CanonicalSplitKey(subtreeTaxa, allTaxa);
            splits.Add(key);
        }

        return subtreeTaxa;
    }

    /// <summary>
    /// Builds a root-invariant key for a bipartition: the lexicographically sorted leaf names
    /// of whichever side is smaller (lexicographically smaller side on a size tie), so the same
    /// {S, complement} pair maps to one key regardless of which side the edge removal exposed.
    /// </summary>
    private static string CanonicalSplitKey(List<string> oneSide, HashSet<string> allTaxa)
    {
        var sideA = oneSide.OrderBy(n => n, StringComparer.Ordinal).ToList();
        var setA = new HashSet<string>(oneSide);
        var sideB = allTaxa.Where(t => !setA.Contains(t))
                           .OrderBy(n => n, StringComparer.Ordinal).ToList();

        // Choose the smaller side (by count, then lexicographically) for a canonical, unordered key.
        bool aFirst = sideA.Count < sideB.Count ||
                      (sideA.Count == sideB.Count &&
                       string.CompareOrdinal(string.Join("|", sideA), string.Join("|", sideB)) <= 0);

        return aFirst ? string.Join("|", sideA) : string.Join("|", sideB);
    }

    private static List<string> CollectClades(PhyloNode node, HashSet<string> clades, int totalLeaves)
    {
        if (node == null) return new List<string>();

        if (node.IsLeaf)
        {
            return new List<string> { node.Name };
        }

        // N-ary: a clade is the set of taxa beneath an internal node, gathered over all children.
        // A multifurcation contributes a single clade (its own subtree) rather than the nested
        // clades a fully-resolved binary subtree would, so an unresolved node yields fewer clades.
        var collected = new List<string>();
        foreach (var child in node.Children)
            collected.AddRange(CollectClades(child, clades, totalLeaves));
        var subtreeTaxa = collected.OrderBy(n => n).ToList();

        // Non-trivial clade: more than one taxon (not a leaf) and fewer than all (not root).
        // Each clade is represented as the sorted, joined taxon names of the subtree.
        if (subtreeTaxa.Count > 1 && subtreeTaxa.Count < totalLeaves)
        {
            clades.Add(string.Join("|", subtreeTaxa));
        }

        return subtreeTaxa;
    }

    /// <summary>
    /// Finds the most recent common ancestor of two taxa.
    /// Returns null when the root is null or either taxon does not exist in the tree.
    /// </summary>
    public static PhyloNode? FindMRCA(PhyloNode root, string taxon1, string taxon2)
    {
        if (root == null) return null;

        var result = FindMRCAInternal(root, taxon1, taxon2);

        // When taxon1 == taxon2, a leaf result is correct (self-MRCA).
        // When taxon1 != taxon2, both taxa found ⇒ MRCA is always an internal node.
        // A leaf result with different taxa means only one was found → other missing.
        if (result != null && result.IsLeaf && taxon1 != taxon2)
            return null;

        return result;
    }

    private static PhyloNode? FindMRCAInternal(PhyloNode node, string taxon1, string taxon2)
    {
        if (node == null) return null;

        if (node.IsLeaf)
        {
            return node.Name == taxon1 || node.Name == taxon2 ? node : null;
        }

        // N-ary: recurse into every child. The MRCA is the deepest node whose subtree contains
        // both taxa. If two (or more) distinct children each return a hit, the split between the
        // taxa happens here, so this node is the MRCA. If exactly one child returns a hit, the
        // MRCA lies within (or is) that child's result and is propagated upward.
        PhyloNode? singleHit = null;
        int childrenWithHit = 0;

        foreach (var child in node.Children)
        {
            var childResult = FindMRCAInternal(child, taxon1, taxon2);
            if (childResult != null)
            {
                childrenWithHit++;
                if (childrenWithHit >= 2)
                    return node;
                singleHit = childResult;
            }
        }

        return singleHit;
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

        // N-ary: descend into whichever child subtree contains the taxon, accumulating that
        // child's branch length onto the distance returned from below it.
        foreach (var child in node.Children)
        {
            double childDist = DistanceToTaxon(child, taxon);
            if (!double.IsNaN(childDist))
                return childDist + child.BranchLength;
        }

        return double.NaN;
    }

    /// <summary>Default number of bootstrap replicates (Felsenstein 1985 uses ≈100+ replicates).</summary>
    private const int DefaultBootstrapReplicates = 100;

    /// <summary>
    /// Default RNG seed for the column-resampling step. Fixed so that, for a given
    /// alignment and parameters, the returned support values are reproducible.
    /// </summary>
    private const int DefaultBootstrapSeed = 42;

    /// <summary>
    /// Felsenstein's phylogenetic bootstrap: estimates clade support by resampling
    /// alignment columns (sites) with replacement and rebuilding the tree on each replicate.
    /// </summary>
    /// <remarks>
    /// Procedure (Felsenstein 1985; Lemoine et al. 2018; Biopython <c>Bio.Phylo.Consensus</c>):
    /// keep all taxa, resample the alignment columns with replacement to a pseudo-alignment of the
    /// <em>same length</em> as the original, rebuild a tree, and for every non-trivial clade of the
    /// reference (original-data) tree count the proportion of replicate trees that contain a clade
    /// with the identical set of leaf names. The returned support is that proportion in [0,1]
    /// (multiply by 100 for the published percentage).
    /// </remarks>
    /// <param name="sequences">Named aligned sequences (≥2, equal length).</param>
    /// <param name="replicates">Number of bootstrap replicates (≥1).</param>
    /// <param name="distanceMethod">Distance method used to build each tree.</param>
    /// <param name="treeMethod">Tree-construction method used for the reference and replicate trees.</param>
    /// <param name="seed">RNG seed for column resampling; fixed value makes results reproducible.</param>
    /// <returns>Map from clade (sorted, '|'-joined leaf names) to bootstrap support in [0,1].</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">When fewer than 2 sequences are supplied, or <paramref name="replicates"/> &lt; 1.</exception>
    public static IReadOnlyDictionary<string, double> Bootstrap(
        IReadOnlyDictionary<string, string> sequences,
        int replicates = DefaultBootstrapReplicates,
        DistanceMethod distanceMethod = DistanceMethod.JukesCantor,
        TreeMethod treeMethod = TreeMethod.UPGMA,
        int seed = DefaultBootstrapSeed)
    {
        if (sequences == null)
            throw new ArgumentNullException(nameof(sequences));
        if (sequences.Count < 2)
            throw new ArgumentException("At least 2 sequences required.", nameof(sequences));
        if (replicates < 1)
            throw new ArgumentException("At least 1 replicate required.", nameof(replicates));

        var taxa = sequences.Keys.ToList();
        var seqs = sequences.Values.ToList();
        int alignmentLength = seqs[0].Length;

        // Build reference tree
        var refTree = BuildTree(sequences, distanceMethod, treeMethod);
        var refClades = GetClades(refTree.Root);

        // Count support for each clade
        var supportCounts = new Dictionary<string, int>();
        foreach (var clade in refClades)
            supportCounts[clade] = 0;

        var random = new Random(seed);

        for (int rep = 0; rep < replicates; rep++)
        {
            // Resample alignmentLength columns with replacement (Felsenstein 1985:
            // pseudo-alignments are the same size as the original; Biopython bootstrap_trees).
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
            var bootClades = GetClades(bootTree.Root);

            // Count matching clades
            foreach (var clade in refClades)
            {
                if (bootClades.Contains(clade))
                    supportCounts[clade]++;
            }
        }

        // Convert to proportions
        return supportCounts.ToDictionary(
            kvp => kvp.Key,
            kvp => (double)kvp.Value / replicates);
    }
}
