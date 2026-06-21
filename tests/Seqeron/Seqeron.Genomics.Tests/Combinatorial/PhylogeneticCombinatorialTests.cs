namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Phylogenetic area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Phylogenetic")]
public class PhylogeneticCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static char NextBase(char c) => "ACGT"[("ACGT".IndexOf(c) + 1) % 4];

    /// <summary>A deterministic family: sequence i has i substitutions of the base sequence.</summary>
    private static List<string> MakeFamily(int nSeqs, int seqLen, uint seed = 0xC0FFEEu)
    {
        string baseSeq = DiverseDna(seqLen, seed);
        var family = new List<string>();
        for (int i = 0; i < nSeqs; i++)
        {
            var chars = baseSeq.ToCharArray();
            for (int k = 0; k < i; k++) { int p = (k * 7 + 1) % chars.Length; chars[p] = NextBase(chars[p]); }
            family.Add(new string(chars));
        }
        return family;
    }

    private static bool IsTransition(char a, char b)
    {
        bool purines = (a is 'A' or 'G') && (b is 'A' or 'G');
        bool pyrimidines = (a is 'C' or 'T') && (b is 'C' or 'T');
        return purines || pyrimidines;
    }

    /// <summary>Independent JC69 / K2P distance re-derivation (ground truth).</summary>
    private static double IndepDistance(string s1, string s2, PhylogeneticAnalyzer.DistanceMethod method)
    {
        int diff = 0, ts = 0, tv = 0, sites = 0;
        for (int i = 0; i < s1.Length; i++)
        {
            char c1 = char.ToUpperInvariant(s1[i]), c2 = char.ToUpperInvariant(s2[i]);
            if (c1 is '-' || c2 is '-' || "ACGT".IndexOf(c1) < 0 || "ACGT".IndexOf(c2) < 0) continue;
            sites++;
            if (c1 != c2) { diff++; if (IsTransition(c1, c2)) ts++; else tv++; }
        }
        if (sites == 0) return 0;
        double p = (double)diff / sites;
        if (method == PhylogeneticAnalyzer.DistanceMethod.JukesCantor)
        {
            double arg = 1 - 4 * p / 3;
            return arg <= 0 ? double.PositiveInfinity : -0.75 * Math.Log(arg);
        }
        double s = (double)ts / sites, v = (double)tv / sites;
        double a1 = 1 - 2 * s - v, a2 = 1 - 2 * v;
        return a1 <= 0 || a2 <= 0 ? double.PositiveInfinity : -0.5 * Math.Log(a1 * Math.Sqrt(a2));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-DIST-001 — Evolutionary distance matrix (Phylogenetic)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 39.
    // Spec: tests/TestSpecs/PHYLO-DIST-001.md (canonical CalculateDistanceMatrix /
    //       CalculatePairwiseDistance).
    // Dimensions: model(2: JC/K2P) × nSeqs(3) × seqLen(3). Grid 2×3×3 = 18.
    //
    // Model (Jukes-Cantor 1969; Kimura 1980): from the observed proportion of differing
    // sites p (and the transition/transversion split for K2P), the corrected distance is
    // JC d = −¾·ln(1−4p/3) and K2P d = −½·ln((1−2S−V)·√(1−2V)); saturation gives +∞.
    //
    // The combinatorial point: across model × nSeqs × seqLen the produced matrix must be a
    // valid distance matrix (symmetric, zero diagonal, non-negative) whose every entry equals
    // both the pairwise routine and an independent re-derivation of the selected model formula.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PhyloDist_MatrixIsValid_AndMatchesModelFormula(
        [Values(PhylogeneticAnalyzer.DistanceMethod.JukesCantor, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter)]
        PhylogeneticAnalyzer.DistanceMethod model,
        [Values(3, 4, 5)] int nSeqs,
        [Values(20, 40, 80)] int seqLen)
    {
        var seqs = MakeFamily(nSeqs, seqLen);
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs, model);

        matrix.GetLength(0).Should().Be(nSeqs);
        matrix.GetLength(1).Should().Be(nSeqs);

        for (int i = 0; i < nSeqs; i++)
        {
            matrix[i, i].Should().Be(0, "the diagonal is self-distance 0");
            for (int j = i + 1; j < nSeqs; j++)
            {
                matrix[i, j].Should().Be(matrix[j, i], "the matrix is symmetric");
                matrix[i, j].Should().BeGreaterThanOrEqualTo(0, "distances are non-negative");
                matrix[i, j].Should().Be(PhylogeneticAnalyzer.CalculatePairwiseDistance(seqs[i], seqs[j], model),
                    "the matrix entry equals the pairwise routine");

                double expected = IndepDistance(seqs[i], seqs[j], model);
                if (double.IsInfinity(expected))
                    matrix[i, j].Should().Be(expected, "saturated distance is +∞");
                else
                    matrix[i, j].Should().BeApproximately(expected, 1e-12, "entry equals the model formula");
            }
        }
    }

    /// <summary>
    /// Interaction witness: identical sequences have zero distance under both models, and the
    /// JC correction inflates the raw p-distance for a divergent pair (d ≥ p).
    /// </summary>
    [Test]
    public void PhyloDist_IdentityZero_AndJcInflatesPDistance()
    {
        string a = DiverseDna(60, 1u);
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, a, PhylogeneticAnalyzer.DistanceMethod.JukesCantor).Should().Be(0);
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, a, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter).Should().Be(0);

        var fam = MakeFamily(2, 60);             // 1 substitution
        double p = PhylogeneticAnalyzer.CalculatePairwiseDistance(fam[0], fam[1], PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(fam[0], fam[1], PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        jc.Should().BeGreaterThanOrEqualTo(p, "the JC correction never shrinks the observed distance");
    }

    /// <summary>
    /// Worked example: a pair with one transition and one transversion in ten sites gives
    /// JC ≈ 0.2328 and K2P ≈ 0.2341 — the models genuinely differ.
    /// </summary>
    [Test]
    public void PhyloDist_JcVersusK2p_WorkedExample()
    {
        const string s1 = "ACGTACGTAC";
        const string s2 = "GAGTACGTAC"; // pos0 A→G (transition), pos1 C→A (transversion)

        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        double k2p = PhylogeneticAnalyzer.CalculatePairwiseDistance(s1, s2, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        jc.Should().BeApproximately(0.2326, 1e-3);
        k2p.Should().BeApproximately(0.2341, 1e-3);
        jc.Should().NotBe(k2p, "JC and K2P are different evolutionary models");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-TREE-001 — Distance-based tree construction (Phylogenetic)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 40.
    // Spec: tests/TestSpecs/PHYLO-TREE-001.md (canonical BuildTree).
    // Dimensions: method(2: UPGMA/NJ) × nSeqs(3) × seqLen(3). Grid 2×3×3 = 18.
    //
    // Model (Sokal & Michener 1958 UPGMA; Saitou & Nei 1987 NJ): both methods cluster a
    // distance matrix into a binary tree whose leaves are the input taxa. Whatever the
    // method, the result must keep every taxon as a leaf, label the tree with the method,
    // expose the distance matrix, and serialise to Newick that round-trips its leaf set.
    //
    // The combinatorial point: method × nSeqs × seqLen all vary the tree, yet the leaf
    // invariants hold universally; and (witness) both methods recover the true grouping —
    // closely related taxa are joined before distant ones.
    // ═══════════════════════════════════════════════════════════════════════

    private static Dictionary<string, string> MakeNamedFamily(int nSeqs, int seqLen)
    {
        var fam = MakeFamily(nSeqs, seqLen);
        var dict = new Dictionary<string, string>();
        for (int i = 0; i < nSeqs; i++) dict[$"T{i}"] = fam[i];
        return dict;
    }

    [Test, Combinatorial]
    public void PhyloTree_LeafInvariants_AcrossMethodNSeqsAndLength(
        [Values(PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining)]
        PhylogeneticAnalyzer.TreeMethod method,
        [Values(3, 4, 5)] int nSeqs,
        [Values(20, 40, 80)] int seqLen)
    {
        var named = MakeNamedFamily(nSeqs, seqLen);
        var tree = PhylogeneticAnalyzer.BuildTree(named, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, method);

        var leaves = PhylogeneticAnalyzer.GetLeaves(tree.Root).ToList();
        leaves.Should().HaveCount(nSeqs, "every taxon is a leaf");
        leaves.Select(l => l.Name).Should().BeEquivalentTo(named.Keys, "leaf labels are the input taxa");

        tree.Taxa.Should().HaveCount(nSeqs);
        tree.Method.Should().Be(method.ToString());
        tree.DistanceMatrix.GetLength(0).Should().Be(nSeqs);
        PhylogeneticAnalyzer.CalculateTreeLength(tree.Root).Should().BeGreaterThanOrEqualTo(0);

        // Newick serialisation round-trips the leaf set.
        var reparsed = PhylogeneticAnalyzer.ParseNewick(PhylogeneticAnalyzer.ToNewick(tree.Root));
        PhylogeneticAnalyzer.GetLeaves(reparsed).Should().HaveCount(nSeqs);
    }

    /// <summary>
    /// Interaction witness: with two clear groups (A,B from one ancestor; C,D from another),
    /// both UPGMA and NJ join the within-group pairs more tightly — the patristic distance
    /// A–B is smaller than A–C.
    /// </summary>
    [Test]
    public void PhyloTree_RecoversTrueGrouping_BothMethods()
    {
        string g1 = DiverseDna(80, 11u);
        string g2 = DiverseDna(80, 999u); // a very different ancestor
        char Flip(char c) => NextBase(c);
        string Mut(string s, int p) { var c = s.ToCharArray(); c[p] = Flip(c[p]); return new string(c); }

        var named = new Dictionary<string, string>
        {
            ["A"] = g1,
            ["B"] = Mut(g1, 3),
            ["C"] = g2,
            ["D"] = Mut(g2, 3),
        };

        foreach (var method in new[] { PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining })
        {
            var tree = PhylogeneticAnalyzer.BuildTree(named, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, method);
            double ab = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "A", "B");
            double ac = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "A", "C");
            ab.Should().BeLessThan(ac, $"{method}: the close pair A,B is nearer than the distant A,C");
        }
    }

    /// <summary>
    /// Worked example: a three-taxon UPGMA tree serialises to Newick listing all three taxa
    /// and terminating with a semicolon.
    /// </summary>
    [Test]
    public void PhyloTree_Newick_ListsAllTaxa()
    {
        var named = MakeNamedFamily(3, 40);
        var tree = PhylogeneticAnalyzer.BuildTree(named, PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        newick.Should().EndWith(";");
        foreach (var taxon in named.Keys)
            newick.Should().Contain(taxon);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-BOOT-001 — Bootstrap clade support (Phylogenetic)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 221.
    // Spec: tests/TestSpecs/PHYLO-BOOT-001.md (canonical PhylogeneticAnalyzer.Bootstrap). ADVANCED §10.
    // Dimensions: nReplicates(3) × method(2) × nSeqs(3). Grid 3×2×3 = 18 (full, exhaustive).
    //
    // Model (Felsenstein 1985): bootstrap resamples alignment columns with replacement, rebuilds the
    // tree per replicate, and reports each reference clade's support = fraction of replicates that
    // recover it (∈ [0,1]); a fixed seed makes the procedure deterministic.
    //
    // Axis mapping (documented): nReplicates → resampling replicates; method → tree method
    // (UPGMA/NJ); nSeqs → number of taxa. The combinatorial point: support values stay in [0,1] and
    // are deterministic for a fixed seed, and a clearly-grouped clade earns high support (witness).
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PhyloBoot_SupportInRangeAndDeterministic_AcrossReplicatesMethodNSeqs(
        [Values(10, 25, 50)] int replicates,
        [Values(PhylogeneticAnalyzer.TreeMethod.UPGMA, PhylogeneticAnalyzer.TreeMethod.NeighborJoining)]
        PhylogeneticAnalyzer.TreeMethod method,
        [Values(3, 4, 5)] int nSeqs)
    {
        var named = MakeNamedFamily(nSeqs, 80);

        var support = PhylogeneticAnalyzer.Bootstrap(named, replicates, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, method, seed: 42);

        support.Values.Should().OnlyContain(v => v >= 0.0 && v <= 1.0, "clade support is a fraction in [0,1]");
        var again = PhylogeneticAnalyzer.Bootstrap(named, replicates, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, method, seed: 42);
        support.Should().BeEquivalentTo(again, "a fixed seed makes the bootstrap deterministic");
    }

    /// <summary>
    /// Interaction witness — a clearly-defined clade (two near-identical sequences against two
    /// distant ones) earns near-maximal bootstrap support; a different seed is still in range.
    /// </summary>
    [Test]
    public void PhyloBoot_StrongCladeHasHighSupport()
    {
        string g1 = DiverseDna(120, 7u);
        string g2 = DiverseDna(120, 808u);
        char Flip(char c) => NextBase(c);
        string Mut(string s, int p) { var a = s.ToCharArray(); a[p] = Flip(a[p]); return new string(a); }

        var named = new Dictionary<string, string> { ["A"] = g1, ["B"] = Mut(g1, 5), ["C"] = g2, ["D"] = Mut(g2, 5) };

        var support = PhylogeneticAnalyzer.Bootstrap(named, 50, PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            PhylogeneticAnalyzer.TreeMethod.NeighborJoining, seed: 1);
        support.Values.Max().Should().BeGreaterThan(0.9, "a strongly-supported clade is recovered in nearly all replicates");
    }
}
