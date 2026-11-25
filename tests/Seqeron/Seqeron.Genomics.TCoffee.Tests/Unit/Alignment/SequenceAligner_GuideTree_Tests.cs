// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
using System.Globalization;
using static Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-12, B09 F53): the guide tree of the progressive / iterative / consistency aligners is built by the
/// canonical <see cref="PhylogeneticAnalyzer"/> builders (the private UPGMA copy of <c>SequenceAligner</c> is gone) and
/// <see cref="MsaGuideTree.NeighborJoining"/> selects ClustalW's neighbour-joining guide tree.
///
/// Oracles (TestData/Alignment/*.tsv):
///  * clustalw_guidetree_oracle.tsv (scripts/alignment/generate_clustalw_guidetree_oracle.py): the installed <c>clustalw</c> 2.1
///    (<c>-TYPE=DNA -ALIGN -MAXDIV=0 -NEWTREE</c>) on 420 DNA sets — the unrooted NJ tree it writes to the <c>.dnd</c> file (5-decimal
///    branch lengths) and the sizes of its progressive merges ("Group n: Sequences: s") = the post-order of its rooted tree.
///  * biopython_guidetree_oracle.tsv (scripts/alignment/generate_guidetree_biopython_oracle.py): Biopython 1.88
///    <c>DistanceTreeConstructor().nj</c> and scipy <c>linkage(method='average')</c> (= UPGMA, size-weighted) on 400 random real
///    distance matrices (no ties). Biopython's own <c>upgma()</c> is not used: it averages without cluster sizes (WPGMA).
///  * the former private UPGMA of SequenceAligner, kept below as <see cref="ReferenceUpgma"/> (verbatim copy of the removed code),
///    on 900 random / tie-heavy matrices.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_GuideTree_Tests
{
    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    #region Tree helpers

    /// <summary>Neutral tree: leaves carry the taxon index, edges their branch length.</summary>
    private sealed class TNode
    {
        public int Leaf = -1;
        public double Len;
        public List<TNode> Kids = new();
    }

    private static TNode FromPhylo(PhyloNode node, Func<string, int> leafId)
    {
        var t = new TNode { Len = node.BranchLength };
        if (node.IsLeaf)
            t.Leaf = leafId(node.Name);
        else
            t.Kids = node.Children.Select(c => FromPhylo(c, leafId)).ToList();
        return t;
    }

    /// <summary>Minimal Newick reader (labels, optional internal labels, ':length'); leaf label = prefix + index.</summary>
    private static TNode ParseNewick(string s, Func<string, int> leafId)
    {
        int pos = 0;

        TNode Parse()
        {
            var node = new TNode();
            if (s[pos] == '(')
            {
                pos++;
                node.Kids.Add(Parse());
                while (s[pos] == ',')
                {
                    pos++;
                    node.Kids.Add(Parse());
                }

                Assert.That(s[pos], Is.EqualTo(')'));
                pos++;
            }

            int start = pos;
            while (pos < s.Length && s[pos] != ':' && s[pos] != ',' && s[pos] != ')' && s[pos] != ';')
                pos++;
            string label = s[start..pos];
            if (node.Kids.Count == 0)
                node.Leaf = leafId(label);
            if (pos < s.Length && s[pos] == ':')
            {
                pos++;
                int ls = pos;
                while (pos < s.Length && s[pos] != ',' && s[pos] != ')' && s[pos] != ';')
                    pos++;
                node.Len = double.Parse(s[ls..pos], CultureInfo.InvariantCulture);
            }

            return node;
        }

        return Parse();
    }

    private static List<int> Leaves(TNode n) =>
        n.Kids.Count == 0 ? new List<int> { n.Leaf } : n.Kids.SelectMany(Leaves).ToList();

    /// <summary>Edge table: leaf set below every non-root node (unrooted: the side without taxon 0) → total branch length.</summary>
    private static Dictionary<string, double> Edges(TNode root, int n, bool unrooted)
    {
        var edges = new Dictionary<string, double>();

        void Visit(TNode node, bool isRoot)
        {
            if (!isRoot)
            {
                var set = Leaves(node);
                if (unrooted && set.Contains(0))
                    set = Enumerable.Range(0, n).Where(i => !set.Contains(i)).ToList();
                string key = string.Join(",", set.OrderBy(x => x));
                edges[key] = edges.GetValueOrDefault(key) + node.Len;
            }

            foreach (var k in node.Kids)
                Visit(k, false);
        }

        Visit(root, true);
        return edges;
    }

    private static string Describe(Dictionary<string, double> e) =>
        string.Join(" ", e.OrderBy(kv => kv.Key.Length).ThenBy(kv => kv.Key).Select(kv => $"[{kv.Key}]={kv.Value:F5}"));

    private static bool SameEdges(Dictionary<string, double> a, Dictionary<string, double> b, double tol)
    {
        if (a.Count != b.Count)
            return false;
        foreach (var kv in a)
            if (!b.TryGetValue(kv.Key, out double v) || Math.Abs(v - kv.Value) > tol)
                return false;
        return true;
    }

    private static string Nested(Seqeron.Genomics.ClustalW.ClustalWAligner.ProgressiveGuideNode node) =>
        node.IsLeaf ? node.LeafIndex.ToString(CultureInfo.InvariantCulture) : $"({Nested(node.Left!)},{Nested(node.Right!)})";

    #endregion

    #region Former private UPGMA (test oracle) vs the canonical builder

    /// <summary>The removed <c>Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree</c>, verbatim (dictionary UPGMA, lowest index wins within 1e-9).</summary>
    private static string ReferenceUpgma(int k, double[,] distance)
    {
        const double eps = 1e-9;
        var nodes = new Dictionary<int, string>();
        var sizes = new Dictionary<int, int>();
        var dist = new Dictionary<(int, int), double>();
        for (int i = 0; i < k; i++)
        {
            nodes[i] = i.ToString(CultureInfo.InvariantCulture);
            sizes[i] = 1;
        }

        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                dist[(i, j)] = distance[i, j];
                dist[(j, i)] = distance[i, j];
            }

        var active = Enumerable.Range(0, k).ToList();
        while (active.Count > 1)
        {
            double minDist = double.MaxValue;
            int minI = -1, minJ = -1;
            for (int ii = 0; ii < active.Count; ii++)
                for (int jj = ii + 1; jj < active.Count; jj++)
                {
                    int a = active[ii];
                    int b = active[jj];
                    var key = a < b ? (a, b) : (b, a);
                    double d = dist.GetValueOrDefault(key, 0);
                    if (d < minDist - eps)
                    {
                        minDist = d;
                        minI = a;
                        minJ = b;
                    }
                }

            int newSize = sizes[minI] + sizes[minJ];
            string newNode = $"({nodes[minI]},{nodes[minJ]})";
            foreach (int c in active)
            {
                if (c == minI || c == minJ)
                    continue;
                var keyIK = minI < c ? (minI, c) : (c, minI);
                var keyJK = minJ < c ? (minJ, c) : (c, minJ);
                double newDist = (dist.GetValueOrDefault(keyIK, 0) * sizes[minI] + dist.GetValueOrDefault(keyJK, 0) * sizes[minJ]) / newSize;
                dist[minI < c ? (minI, c) : (c, minI)] = newDist;
            }

            nodes[minI] = newNode;
            sizes[minI] = newSize;
            nodes.Remove(minJ);
            sizes.Remove(minJ);
            active.Remove(minJ);
        }

        return nodes[active[0]];
    }

    private static double[,] RandomMatrix(Random rng, int k, int kind)
    {
        var m = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
            {
                double d = kind switch
                {
                    0 => rng.NextDouble(),                       // generic reals
                    1 => rng.Next(0, 5),                         // integers 0..4: massive exact ties
                    2 => rng.Next(0, 11) / 10.0,                 // tenths: ties that only differ by rounding noise (0.3+0.4 vs 0.7)
                    _ => 1.0 - rng.Next(0, 21) / 20.0,           // ClustalW-like 1 - identities/len
                };
                m[i, j] = m[j, i] = d;
            }

        return m;
    }

    [Test]
    public void UpgmaGuideTree_EqualsFormerPrivateImplementation_On900RandomAndTieHeavyMatrices()
    {
        var rng = new Random(12345);
        int checkedCount = 0;
        for (int kind = 0; kind < 4; kind++)
            for (int t = 0; t < 225; t++)
            {
                int k = rng.Next(2, 15);
                var d = RandomMatrix(rng, k, kind);
                string expected = ReferenceUpgma(k, d);
                string actual = Nested(Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree(k, d));
                Assert.That(actual, Is.EqualTo(expected), $"kind {kind}, k {k}");
                checkedCount++;
            }

        Assert.That(checkedCount, Is.EqualTo(900));
    }

    [Test]
    public void UpgmaGuideTree_HandRun_TieGoesToLowestIndexPair()
    {
        // d(0,3) = d(1,2) = 1 are the minimum: the lowest (i, j) pair (0,3) merges first, then (1,2), then the two clusters.
        var d = new double[4, 4];
        double[,] raw =
        {
            { 0, 5, 5, 1 },
            { 5, 0, 1, 5 },
            { 5, 1, 0, 5 },
            { 1, 5, 5, 0 },
        };
        Array.Copy(raw, d, raw.Length);
        Assert.That(Nested(Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree(4, d)), Is.EqualTo("((0,3),(1,2))"));
    }

    [Test]
    public void UpgmaGuideTree_ValueFromRoundingNoiseTies()
    {
        // 0.3 + 0.4 vs 0.7 cannot be told apart beyond rounding noise: the tie must go to the lowest index pair.
        double[,] d =
        {
            { 0.0, 0.7, 0.3, 0.9 },
            { 0.7, 0.0, 0.9, 0.3 },
            { 0.3, 0.9, 0.0, 0.8 },
            { 0.9, 0.3, 0.8, 0.0 },
        };
        Assert.That(Nested(Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree(4, d)), Is.EqualTo(ReferenceUpgma(4, d)));
    }

    #endregion

    #region Canonical builders vs Biopython 1.88 (NJ) / scipy average linkage (UPGMA)

    private static readonly Lazy<List<string[]>> BiopythonRows = new(() => ReadTsv("biopython_guidetree_oracle.tsv").ToList());

    private static double[,] ParseMatrix(string text, int n)
    {
        var m = new double[n, n];
        var rows = text.Split(';');
        for (int i = 0; i < n; i++)
        {
            var cells = rows[i].Split(',');
            for (int j = 0; j < n; j++)
                m[i, j] = double.Parse(cells[j], CultureInfo.InvariantCulture);
        }

        return m;
    }

    [TestCase(false)]
    [TestCase(true)]
    public void CanonicalBuilders_MatchBiopythonNj_AndScipyAverageLinkageUpgma_OnRandomMatrices(bool neighborJoining)
    {
        var rows = BiopythonRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(400));
        int bad = 0;
        string? first = null;
        foreach (var f in rows)
        {
            int n = int.Parse(f[0]);
            var m = ParseMatrix(f[1], n);
            var taxa = Enumerable.Range(0, n).Select(i => "t" + i).ToList();
            var tree = BuildTreeFromMatrix(taxa, m, neighborJoining ? TreeMethod.NeighborJoining : TreeMethod.UPGMA);
            int Id(string name) => int.Parse(name.AsSpan(1), CultureInfo.InvariantCulture);
            var actual = Edges(FromPhylo(tree.Root, Id), n, neighborJoining);
            var expected = Edges(ParseNewick(neighborJoining ? f[3] : f[2], Id), n, neighborJoining);
            if (!SameEdges(actual, expected, 1e-6) && ++bad == 1)
                first = $"n={n}\nactual   {Describe(actual)}\nexpected {Describe(expected)}";
        }

        Assert.That(bad, Is.Zero, first);
    }

    #endregion

    #region ClustalW 2.1 neighbour-joining guide tree

    private static readonly Lazy<List<string[]>> ClustalWRows = new(() => ReadTsv("clustalw_guidetree_oracle.tsv").ToList());

    private static double[,] ClustalWDistances(string[] seqs)
    {
        int k = seqs.Length;
        var d = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                d[i, j] = d[j, i] = Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWSinglePrecisionDistance(seqs[i], seqs[j]);
        return d;
    }

    [Test]
    public void ClustalWNeighborJoining_UnrootedTreeAndStepOrder_MatchClustalW21Oracle()
    {
        var rows = ClustalWRows.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(400));
        int badTree = 0, badSteps = 0;
        string? firstTree = null, firstSteps = null;
        foreach (var f in rows)
        {
            var seqs = f[0].Split(';');
            int k = seqs.Length;
            var result = Seqeron.Genomics.ClustalW.ClustalWAligner.BuildClustalWNeighborJoiningGuideTree(k, ClustalWDistances(seqs));

            // The .dnd file: ClustalW prints |b| < 1e-4 as 0 and the others with 5 decimals.
            int Id(string name) => int.Parse(name.AsSpan(1), CultureInfo.InvariantCulture);
            var expected = Edges(ParseNewick(f[1], Id), k, unrooted: true);
            var actualTree = FromPhylo(result.Unrooted, name => int.Parse(name, CultureInfo.InvariantCulture));
            RoundLikeClustalW(actualTree);
            var actual = Edges(actualTree, k, unrooted: true);
            if (!SameEdges(actual, expected, 1.5e-5) && ++badTree == 1)
                firstTree = $"{f[0]}\nactual   {Describe(actual)}\nexpected {Describe(expected)}";

            var steps = Seqeron.Genomics.ClustalW.ClustalWAligner.GuideTreeStepSizes(result.Rooted);
            var expectedSteps = f[2].Split(',').Select(int.Parse).ToList();
            if (!steps.SequenceEqual(expectedSteps) && ++badSteps == 1)
                firstSteps = $"{f[0]}\nactual {string.Join(",", steps)} expected {f[2]}\n{Nested(result.Rooted)}";
        }

        Assert.Multiple(() =>
        {
            Assert.That(badTree, Is.Zero, firstTree);
            Assert.That(badSteps, Is.Zero, firstSteps);
        });
    }

    private static void RoundLikeClustalW(TNode node)
    {
        double b = Math.Abs(node.Len) < 0.0001 ? 0.0 : node.Len;
        node.Len = double.Parse(b.ToString("F5", CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
        foreach (var k in node.Kids)
            RoundLikeClustalW(k);
    }

    [Test]
    public void ClustalWNeighborJoining_OracleCoversTiesNegativeBranchesAndAllSizes()
    {
        var rows = ClustalWRows.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(f => f[0].Split(';').Length).Distinct().Min(), Is.EqualTo(2));
            Assert.That(rows.Select(f => f[0].Split(';').Length).Max(), Is.GreaterThanOrEqualTo(10));
            Assert.That(rows.Count(f => f[1].Contains(":-")), Is.GreaterThan(20), "negative NJ branch lengths");
            Assert.That(rows.Count(f => f[0].Split(';').Distinct().Count() < f[0].Split(';').Length), Is.GreaterThan(30), "duplicate sequences (tied distances)");
        });
    }

    [Test]
    public void ClustalWNeighborJoining_TwoSequences_SingleStep()
    {
        var d = new double[2, 2];
        d[0, 1] = d[1, 0] = 0.4;
        var r = Seqeron.Genomics.ClustalW.ClustalWAligner.BuildClustalWNeighborJoiningGuideTree(2, d);
        Assert.That(Nested(r.Rooted), Is.EqualTo("(0,1)"));
    }

    [Test]
    public void ClustalWNeighborJoining_FewerThanTwoSequences_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.BuildClustalWNeighborJoiningGuideTree(1, new double[1, 1]));
    }

    #endregion

    #region Canonical tie options

    private static string EdgeKeys(PhylogeneticTree tree, int n) =>
        string.Join(" ", Edges(FromPhylo(tree.Root, name => int.Parse(name, CultureInfo.InvariantCulture)), n, unrooted: true)
            .Keys.Where(k => k.Contains(',') && k.Split(',').Length < n - 1).OrderBy(k => k));

    [Test]
    public void TreeBuildOptions_TieOrder_FirstIndexMajorIsDefault_SecondIndexMajorIsClustalWOrder()
    {
        // Neighbour-joining on an integer matrix where the disjoint pairs (0,4) and (1,3) tie for the smallest Q. Choosing one first
        // changes the Q of the other afterwards, so the two scan orders give different trees:
        //   lexicographic (a, b) scan, default   : (0,4) joins first -> splits {0,4} and {2,3}   i.e. ((0,4),1)|2|3
        //   (b, a) scan = ClustalW NJTree loops : (1,3) joins first -> splits {1,3} and {2,4}   i.e. (0,(1,3))|2|4
        // (hand-run Saitou-Nei steps below; the matrix was found by exhaustive search over 400 000 random integer matrices).
        double[,] d =
        {
            { 0, 2, 3, 3, 1 },
            { 2, 0, 3, 2, 2 },
            { 3, 3, 0, 2, 1 },
            { 3, 2, 2, 0, 3 },
            { 1, 2, 1, 3, 0 },
        };
        var taxa = new[] { "0", "1", "2", "3", "4" };
        var def = BuildTreeFromMatrix(taxa, d, TreeMethod.NeighborJoining);
        var first = BuildTreeFromMatrix(taxa, d, TreeMethod.NeighborJoining, new TreeBuildOptions { TieTolerance = 1e-9, TieOrder = PairScanOrder.FirstIndexMajor });
        var second = BuildTreeFromMatrix(taxa, d, TreeMethod.NeighborJoining, new TreeBuildOptions { TieTolerance = 1e-9, TieOrder = PairScanOrder.SecondIndexMajor });
        Assert.Multiple(() =>
        {
            Assert.That(EdgeKeys(def, 5), Is.EqualTo("1,2,3 2,3"));   // {0,4} (normalised to the side without taxon 0) and {2,3}
            Assert.That(EdgeKeys(first, 5), Is.EqualTo("1,2,3 2,3"));
            Assert.That(EdgeKeys(second, 5), Is.EqualTo("1,3 2,4"));
        });
    }

    [Test]
    public void TreeBuildOptions_TieTolerance_TreatsRoundingNoiseAsTie()
    {
        // UPGMA: d(0,1) = 0.3 + 1e-12 is a hair above d(1,2) = 0.3. Exact comparison merges (1,2) first -> (a,(b,c));
        // tolerance 1e-9 makes it a tie, and the lowest index pair (0,1) wins -> ((a,b),c).
        double[,] d =
        {
            { 0, 0.3 + 1e-12, 0.9 },
            { 0.3 + 1e-12, 0, 0.3 },
            { 0.9, 0.3, 0 },
        };
        var taxa = new[] { "a", "b", "c" };
        var exact = BuildTreeFromMatrix(taxa, d, TreeMethod.UPGMA);
        var tolerant = BuildTreeFromMatrix(taxa, d, TreeMethod.UPGMA, new TreeBuildOptions { TieTolerance = 1e-9 });
        Assert.Multiple(() =>
        {
            Assert.That(exact.Root.Children[0].Taxa, Is.EqualTo(new[] { "a" }));
            Assert.That(exact.Root.Children[1].Taxa, Is.EqualTo(new[] { "b", "c" }));
            Assert.That(tolerant.Root.Children[0].Taxa, Is.EqualTo(new[] { "a", "b" }));
            Assert.That(tolerant.Root.Children[1].Taxa, Is.EqualTo(new[] { "c" }));
        });
    }

    [Test]
    public void TreeBuildOptions_Default_ReproducesTheParameterlessOverload()
    {
        var rng = new Random(7);
        for (int t = 0; t < 100; t++)
        {
            int k = rng.Next(2, 12);
            var m = RandomMatrix(rng, k, t % 4);
            var taxa = Enumerable.Range(0, k).Select(i => "x" + i).ToList();
            foreach (var method in new[] { TreeMethod.UPGMA, TreeMethod.NeighborJoining })
            {
                var a = BuildTreeFromMatrix(taxa, m, method);
                var b = BuildTreeFromMatrix(taxa, m, method, new TreeBuildOptions());
                var c = BuildTreeFromMatrix(taxa, m, method, null);
                Assert.That(ToNewick(b.Root), Is.EqualTo(ToNewick(a.Root)));
                Assert.That(ToNewick(c.Root), Is.EqualTo(ToNewick(a.Root)));
            }
        }
    }

    #endregion

    #region ConsistencyMsaOptions.GuideTree through the aligners

    private static List<DnaSequence> RandomFamily(Random rng, int n)
    {
        string baseSeq = new(Enumerable.Range(0, rng.Next(12, 30)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var list = new List<DnaSequence>();
        for (int i = 0; i < n; i++)
        {
            var sb = new System.Text.StringBuilder();
            foreach (char c in baseSeq)
            {
                double r = rng.NextDouble();
                if (r < 0.15)
                    sb.Append("ACGT"[rng.Next(4)]);
                else if (r < 0.2)
                    continue;
                else
                    sb.Append(c);
            }

            if (sb.Length == 0)
                sb.Append('A');
            list.Add(new DnaSequence(sb.ToString()));
        }

        return list;
    }

    private static void AssertValidMsa(IReadOnlyList<DnaSequence> input, MultipleAlignmentResult r)
    {
        Assert.That(r.AlignedSequences.Length, Is.EqualTo(input.Count));
        int w = r.AlignedSequences[0].Length;
        for (int i = 0; i < input.Count; i++)
        {
            Assert.That(r.AlignedSequences[i].Length, Is.EqualTo(w));
            Assert.That(r.AlignedSequences[i].Replace("-", ""), Is.EqualTo(input[i].Sequence));
        }
    }

    [Test]
    public void GuideTreeOption_DefaultIsUpgma_ExplicitUpgmaIsIdentical()
    {
        Assert.That(new ConsistencyMsaOptions().GuideTree, Is.EqualTo(MsaGuideTree.Upgma));
        var rng = new Random(99);
        for (int t = 0; t < 40; t++)
        {
            var seqs = RandomFamily(rng, rng.Next(3, 8));
            var a = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna);
            var b = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna, new ConsistencyMsaOptions { GuideTree = MsaGuideTree.Upgma });
            Assert.That(b.AlignedSequences, Is.EqualTo(a.AlignedSequences));
            Assert.That(b.TotalScore, Is.EqualTo(a.TotalScore));
        }
    }

    [TestCase(MsaDistanceMode.Linear)]
    [TestCase(MsaDistanceMode.ClustalWAffine)]
    public void NeighborJoiningGuideTree_ProducesValidAlignments_OnAllThreeGuideTreeAligners(MsaDistanceMode distance)
    {
        var rng = new Random(2024);
        var options = new ConsistencyMsaOptions { GuideTree = MsaGuideTree.NeighborJoining, Distance = distance };
        for (int t = 0; t < 40; t++)
        {
            var seqs = RandomFamily(rng, rng.Next(2, 9));
            var prog = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna, options);
            var iter = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(seqs, SequenceAligner.SimpleDna, options, 3);
            var cons = Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(seqs, SequenceAligner.SimpleDna, options);
            AssertValidMsa(seqs, prog);
            AssertValidMsa(seqs, iter);
            AssertValidMsa(seqs, cons);
            Assert.That(iter.TotalScore, Is.GreaterThanOrEqualTo(prog.TotalScore), "refinement never lowers the sum-of-pairs score");
        }
    }

    [Test]
    public void NeighborJoiningGuideTree_ChangesTheMergeOrder_OnANonUltrametricSet()
    {
        // Saitou & Nei's point: UPGMA assumes a molecular clock. Taxa 0,1 are close; 2 sits on a long branch next to 1 in
        // distance terms; the NJ and UPGMA trees differ in the rooted structure on this matrix.
        double[,] d =
        {
            { 0.0, 0.3, 0.9, 1.0, 1.0 },
            { 0.3, 0.0, 0.8, 0.9, 0.9 },
            { 0.9, 0.8, 0.0, 0.4, 0.5 },
            { 1.0, 0.9, 0.4, 0.0, 0.1 },
            { 1.0, 0.9, 0.5, 0.1, 0.0 },
        };
        string upgma = Nested(Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree(5, d, MsaGuideTree.Upgma));
        var nj = Seqeron.Genomics.ClustalW.ClustalWAligner.BuildProgressiveGuideTree(5, d, MsaGuideTree.NeighborJoining);
        Assert.Multiple(() =>
        {
            Assert.That(upgma, Is.EqualTo("((0,1),(2,(3,4)))"));
            // NJ: {0,1} and {3,4} are the cherries; ClustalW's reRoot rule puts the root between {0,1} and {2,3,4}.
            Assert.That(Seqeron.Genomics.ClustalW.ClustalWAligner.GuideTreeStepSizes(nj), Is.EqualTo(new[] { 2, 3, 2, 5 }));
            Assert.That(Nested(nj), Is.EqualTo("(((4,3),2),(0,1))"));
        });
    }

    [Test]
    public void InvalidGuideTree_Throws()
    {
        var seqs = new List<DnaSequence> { new("ACGT"), new("ACGA"), new("ACCT") };
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna, new ConsistencyMsaOptions { GuideTree = (MsaGuideTree)42 }));
    }

    #endregion

    #region Layering (A1-12: the Phylogenetics -> Alignment reference was unused and is gone)

    [Test]
    public void Phylogenetics_DoesNotReferenceAlignment_AlignmentReferencesPhylogenetics()
    {
        var phylo = typeof(PhylogeneticAnalyzer).Assembly.GetReferencedAssemblies().Select(a => a.Name).ToList();
        var alignment = typeof(Seqeron.Genomics.ClustalW.ClustalWAligner).Assembly.GetReferencedAssemblies().Select(a => a.Name).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(phylo, Does.Not.Contain("Seqeron.Genomics.Alignment"));
            Assert.That(alignment, Does.Contain("Seqeron.Genomics.Phylogenetics"));
        });
    }

    #endregion
}
