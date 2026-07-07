// 08_DIFFERENTIAL_TESTING rows 39-42 (Phylogenetics). Independent oracles: closed-form JC69/K80 with
// independent site counting, both-trees-valid property checks, a hand-built Newick string, and a brute
// bipartition symmetric-difference for the Robinson-Foulds distance.

using Node = Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.PhyloNode;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class PhylogeneticsDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 39: PHYLO-DIST-001 — JC / K2P vs closed-form with independent site counting ----

    private static (double jc, double k2p) DistanceOracle(string a, string b)
    {
        int comp = 0, diff = 0, ts = 0, tv = 0;
        const string std = "ACGT";
        for (int i = 0; i < a.Length; i++)
        {
            char c1 = char.ToUpperInvariant(a[i]), c2 = char.ToUpperInvariant(b[i]);
            if (!std.Contains(c1) || !std.Contains(c2)) continue;
            comp++;
            if (c1 == c2) continue;
            diff++;
            bool purines = (c1 is 'A' or 'G') && (c2 is 'A' or 'G');
            bool pyrimidines = (c1 is 'C' or 'T') && (c2 is 'C' or 'T');
            if (purines || pyrimidines) ts++; else tv++;
        }
        if (comp == 0) return (0, 0);
        double p = (double)diff / comp, s = (double)ts / comp, v = (double)tv / comp;
        double jc = -0.75 * Math.Log(1 - 4 * p / 3);
        double k2p = -0.5 * Math.Log((1 - 2 * s - v) * Math.Sqrt(1 - 2 * v));
        return (jc, k2p);
    }

    [Test]
    [Category("PHYLO-DIST-001")]
    [TestCase("ACGTACGTAC", "ACGAACGTTC")]
    [TestCase("AAAACCCCGGGGTTTT", "AAGACCTCGCGGTATT")]
    [TestCase("ACGTACGT", "ACGTACGT")]   // identical -> both 0
    public void PairwiseDistance_JcAndK2p_MatchClosedForm(string a, string b)
    {
        var (jc, k2p) = DistanceOracle(a, b);
        Assert.That(PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.JukesCantor),
            Is.EqualTo(jc).Within(Tol), "JC69");
        Assert.That(PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter),
            Is.EqualTo(k2p).Within(Tol), "K80");
        // Per the checklist relation: with transversions present, K2P >= JC >= 0.
        Assert.That(jc, Is.GreaterThanOrEqualTo(0));
        Assert.That(k2p, Is.GreaterThanOrEqualTo(jc - Tol));
    }

    // ---- Row 40: PHYLO-TREE-001 — UPGMA and NJ both produce valid trees over the same taxa ----

    private static List<string> LeafNames(Node n)
    {
        if (n.IsLeaf) return new List<string> { n.Name };
        return n.Children.SelectMany(LeafNames).ToList();
    }

    [Test]
    [Category("PHYLO-TREE-001")]
    public void UpgmaAndNj_BothProduceValidTreesOverInputTaxa()
    {
        var seqs = new Dictionary<string, string>
        {
            ["A"] = "AAAACCCCGG",
            ["B"] = "AAAACCCCGT",
            ["C"] = "AAGACCTCGG",
            ["D"] = "TAGACCTCGA",
        };
        var taxa = seqs.Keys.OrderBy(x => x).ToList();

        foreach (var method in new[] { PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining })
        {
            var tree = PhylogeneticAnalyzer.BuildTree(seqs, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, method);
            var leaves = LeafNames(tree.Root).OrderBy(x => x).ToList();
            Assert.That(leaves, Is.EqualTo(taxa), $"{method}: leaf set equals input taxa (no dups/drops)");
            Assert.That(tree.Taxa.OrderBy(x => x), Is.EqualTo(taxa), $"{method}: Taxa property");
        }
    }

    // ---- Row 41: PHYLO-NEWICK-001 — ToNewick vs hand-built string ----

    private static Node Leaf(string name, double bl) => new(name) { BranchLength = bl };
    private static Node Internal(double bl, params Node[] children)
    {
        var n = new Node { BranchLength = bl };
        n.Children.AddRange(children);
        return n;
    }

    [Test]
    [Category("PHYLO-NEWICK-001")]
    public void ToNewick_MatchesHandBuiltString()
    {
        // root( A:0.1 , inner( C:0.2 , D:0.25 ):0.3 )
        var root = Internal(0, Leaf("A", 0.1), Internal(0.3, Leaf("C", 0.2), Leaf("D", 0.25)));
        Assert.That(PhylogeneticAnalyzer.ToNewick(root),
            Is.EqualTo("(A:0.1000,(C:0.2000,D:0.2500):0.3000);"));
    }

    // ---- Row 42: PHYLO-COMP-001 — unrooted RF vs brute bipartition symmetric difference ----

    private static Node Pair(string a, string b) => Internal(1, Leaf(a, 1), Leaf(b, 1));

    private static HashSet<string> BipartitionOracle(Node root)
    {
        var allLeaves = LeafNames(root).ToHashSet();
        int n = allLeaves.Count;
        var splits = new HashSet<string>();

        List<string> Collect(Node node)
        {
            if (node.IsLeaf) return new List<string> { node.Name };
            var sub = node.Children.SelectMany(Collect).ToList();
            int side = sub.Count;
            if (side >= 2 && side <= n - 2) // non-trivial split
            {
                var smaller = side <= n - side ? sub : allLeaves.Except(sub).ToList();
                var other = allLeaves.Except(smaller).ToList();
                string ca = string.Join(",", smaller.OrderBy(x => x, StringComparer.Ordinal));
                string cb = string.Join(",", other.OrderBy(x => x, StringComparer.Ordinal));
                // Canonicalise to the lexicographically smaller side (handles the even split tie).
                splits.Add(string.CompareOrdinal(ca, cb) <= 0 ? ca : cb);
            }
            return sub;
        }
        Collect(root);
        return splits;
    }

    private static int RfOracle(Node t1, Node t2)
    {
        var s1 = BipartitionOracle(t1);
        var s2 = BipartitionOracle(t2);
        return s1.Except(s2).Count() + s2.Except(s1).Count();
    }

    [Test]
    [Category("PHYLO-COMP-001")]
    public void RobinsonFoulds_MatchesBruteBipartitionDifference()
    {
        var ab_cd = Internal(0, Pair("A", "B"), Pair("C", "D"));
        var ab_cd2 = Internal(0, Pair("A", "B"), Pair("C", "D"));
        var ac_bd = Internal(0, Pair("A", "C"), Pair("B", "D"));

        // Same topology -> RF 0; the two distinct 4-taxon topologies -> RF 2.
        Assert.That(PhylogeneticAnalyzer.CalculateUnrootedRobinsonFoulds(ab_cd, ab_cd2), Is.EqualTo(0));
        Assert.That(PhylogeneticAnalyzer.CalculateUnrootedRobinsonFoulds(ab_cd, ac_bd), Is.EqualTo(2));

        // ... and the production value equals the independent bipartition symmetric difference.
        Assert.That(PhylogeneticAnalyzer.CalculateUnrootedRobinsonFoulds(ab_cd, ac_bd), Is.EqualTo(RfOracle(ab_cd, ac_bd)));
        Assert.That(PhylogeneticAnalyzer.CalculateUnrootedRobinsonFoulds(ab_cd, ab_cd2), Is.EqualTo(RfOracle(ab_cd, ab_cd2)));
    }
}
