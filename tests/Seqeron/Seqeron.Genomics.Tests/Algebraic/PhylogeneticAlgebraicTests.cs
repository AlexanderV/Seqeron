using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

using PhyloNode = PhylogeneticAnalyzer.PhyloNode;
using DistanceMethod = PhylogeneticAnalyzer.DistanceMethod;

/// <summary>
/// Algebraic-law tests for the Phylogenetic area (genetic distance, Newick I/O,
/// Robinson–Foulds tree comparison).
///
/// Algebraic testing pins the metric-space axioms of the distance functions and
/// the parse∘serialize round-trip of the Newick representation.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 39, 41, 42.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Phylogenetic")]
public class PhylogeneticAlgebraicTests
{
    /// <summary>Three gap-free A/C/G/T sequences of equal length, the domain on
    /// which the per-site distances are genuine metrics.</summary>
    private static Arbitrary<(string A, string B, string C)> EqualLengthTriple() =>
        (from n in Gen.Choose(1, 40)
         from a in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         from b in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         from c in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         select (new string(a), new string(b), new string(c)))
        .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-DIST-001 — Pairwise genetic distance (Phylogenetic)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 39.
    //
    // Model: a per-site distance over aligned sequences. The uncorrected
    //        distances (raw Hamming difference count, and p-distance = differing
    //        sites / comparable sites) are genuine metrics on equal-length gap-free
    //        sequences. The JC69 / K80 evolutionary corrections are monotone
    //        nonlinear transforms of p and are NOT required to satisfy the
    //        triangle inequality, so TRI is asserted only on the uncorrected
    //        distances; ID and COMM hold for every method.
    //   — docs/algorithms/Phylogenetics; PhylogeneticAnalyzer.CalculatePairwiseDistance.
    //
    // Laws under test (checklist row 39):
    //   • ID   — d(x, x) = 0 (identity of indiscernibles).
    //   • COMM — d(a, b) = d(b, a) (symmetry).
    //   • TRI  — d(a, c) ≤ d(a, b) + d(b, c) (triangle inequality).
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly DistanceMethod[] AllMethods =
    {
        DistanceMethod.Hamming, DistanceMethod.PDistance, DistanceMethod.JukesCantor
    };

    /// <summary>
    /// ID: d(x, x) = 0 for every distance method.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Distance_Identity_SelfDistanceIsZero()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            bool ok = AllMethods.All(m =>
                PhylogeneticAnalyzer.CalculatePairwiseDistance(t.A, t.A, m) == 0.0);
            return ok.Label($"d(x,x) != 0 for \"{t.A}\"");
        });
    }

    /// <summary>
    /// COMM: d(a, b) = d(b, a) for every distance method.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Distance_Commutative_Symmetric()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            bool ok = AllMethods.All(m =>
                PhylogeneticAnalyzer.CalculatePairwiseDistance(t.A, t.B, m)
                == PhylogeneticAnalyzer.CalculatePairwiseDistance(t.B, t.A, m));
            return ok.Label($"d(a,b) != d(b,a) for \"{t.A}\"/\"{t.B}\"");
        });
    }

    /// <summary>
    /// TRI: d(a, c) ≤ d(a, b) + d(b, c) for the genuine metrics (Hamming, p-distance).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Distance_TriangleInequality()
    {
        return Prop.ForAll(EqualLengthTriple(), t =>
        {
            foreach (var m in new[] { DistanceMethod.Hamming, DistanceMethod.PDistance })
            {
                double ac = PhylogeneticAnalyzer.CalculatePairwiseDistance(t.A, t.C, m);
                double ab = PhylogeneticAnalyzer.CalculatePairwiseDistance(t.A, t.B, m);
                double bc = PhylogeneticAnalyzer.CalculatePairwiseDistance(t.B, t.C, m);
                if (ac > ab + bc + 1e-9)
                    return false.Label($"{m}: d(a,c)={ac} > {ab}+{bc}");
            }
            return true.ToProperty();
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-NEWICK-001 — Newick tree I/O (Phylogenetic)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 41.
    //
    // Model: ToNewick serializes a PhyloNode tree to the Newick grammar and
    //        ParseNewick parses it back. On a branch-length-free canonical form the
    //        serialization is the tree's normal form, so parse∘serialize is the
    //        identity isomorphism: ToNewick(ParseNewick(s)) = s.
    //   — docs/algorithms/Phylogenetics; PhylogeneticAnalyzer.ToNewick / ParseNewick.
    //
    // Laws under test (checklist row 41):
    //   • RT — ToNewick(ParseNewick(s)) = s for canonical (no-branch-length) Newick,
    //          and the tree is structurally invariant across a parse-reserialize cycle.
    //   • ID — a single leaf serializes to the simple Newick "name;".
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] CanonicalNewicks =
    {
        "(A,B);",
        "((A,B),C);",
        "((A,B),(C,D));",
        "(A,(B,(C,D)));",
        "((A,C),(B,D));",
        "((A,B),(C,(D,E)));"
    };

    private static bool SameTopology(PhyloNode x, PhyloNode y)
    {
        if (x.Name != y.Name) return false;
        if (x.Children.Count != y.Children.Count) return false;
        for (int i = 0; i < x.Children.Count; i++)
            if (!SameTopology(x.Children[i], y.Children[i])) return false;
        return true;
    }

    /// <summary>
    /// RT: ToNewick(ParseNewick(s)) = s on the canonical no-branch-length form, and
    /// the parsed tree is structurally identical before and after a reserialize cycle.
    /// </summary>
    [Test]
    public void Newick_RoundTrip_ParseSerializeIsIdentity()
    {
        foreach (var s in CanonicalNewicks)
        {
            var tree = PhylogeneticAnalyzer.ParseNewick(s);
            string serialized = PhylogeneticAnalyzer.ToNewick(tree, includeBranchLengths: false);
            serialized.Should().Be(s, $"round-trip of \"{s}\"");

            var reparsed = PhylogeneticAnalyzer.ParseNewick(serialized);
            SameTopology(tree, reparsed).Should().BeTrue($"topology drift on \"{s}\"");
        }
    }

    /// <summary>
    /// ID: a single leaf serializes to the simple Newick "name;" and round-trips.
    /// </summary>
    [Test]
    public void Newick_Identity_LeafIsSimpleNewick()
    {
        var leaf = PhylogeneticAnalyzer.ParseNewick("A;");
        leaf.IsLeaf.Should().BeTrue();
        leaf.Name.Should().Be("A");
        PhylogeneticAnalyzer.ToNewick(leaf, includeBranchLengths: false).Should().Be("A;");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PHYLO-COMP-001 — Robinson–Foulds tree comparison (Phylogenetic)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 42.
    //
    // Model: the rooted Robinson–Foulds distance is the cardinality of the
    //        symmetric difference of the two trees' clade sets,
    //        RF(a,b) = |clades(a) △ clades(b)|. Set symmetric difference is a
    //        metric, so RF inherits the metric axioms.
    //   — docs/algorithms/Phylogenetics; PhylogeneticAnalyzer.RobinsonFouldsDistance.
    //
    // Laws under test (checklist row 42):
    //   • ID   — RF(t, t) = 0.
    //   • COMM — RF(a, b) = RF(b, a).
    //   • TRI  — RF(a, c) ≤ RF(a, b) + RF(b, c) (symmetric-difference triangle ineq.).
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] DistinctTopologies =
    {
        "((A,B),(C,(D,E)));",
        "((A,C),(B,(D,E)));",
        "((A,B),((C,D),E));",
        "(((A,B),C),(D,E));",
        "((A,D),(B,(C,E)));"
    };

    private static PhyloNode[] Trees() =>
        DistinctTopologies.Select(PhylogeneticAnalyzer.ParseNewick).ToArray();

    /// <summary>
    /// ID: RF(t, t) = 0 for every tree.
    /// </summary>
    [Test]
    public void RobinsonFoulds_Identity_SelfDistanceIsZero()
    {
        foreach (var t in Trees())
            PhylogeneticAnalyzer.RobinsonFouldsDistance(t, t).Should().Be(0);
    }

    /// <summary>
    /// COMM: RF(a, b) = RF(b, a) for every pair.
    /// </summary>
    [Test]
    public void RobinsonFoulds_Commutative_Symmetric()
    {
        var trees = Trees();
        foreach (var a in trees)
            foreach (var b in trees)
                PhylogeneticAnalyzer.RobinsonFouldsDistance(a, b)
                    .Should().Be(PhylogeneticAnalyzer.RobinsonFouldsDistance(b, a));
    }

    /// <summary>
    /// TRI: RF(a, c) ≤ RF(a, b) + RF(b, c) over all triples.
    /// </summary>
    [Test]
    public void RobinsonFoulds_TriangleInequality()
    {
        var trees = Trees();
        foreach (var a in trees)
            foreach (var b in trees)
                foreach (var c in trees)
                {
                    int ac = PhylogeneticAnalyzer.RobinsonFouldsDistance(a, c);
                    int ab = PhylogeneticAnalyzer.RobinsonFouldsDistance(a, b);
                    int bc = PhylogeneticAnalyzer.RobinsonFouldsDistance(b, c);
                    ac.Should().BeLessThanOrEqualTo(ab + bc);
                }
    }
}
