// PHYLO-BOOT-001 — Phylogenetic Bootstrap Analysis
// Evidence: docs/Evidence/PHYLO-BOOT-001-Evidence.md
// TestSpec: tests/TestSpecs/PHYLO-BOOT-001.md
// Source: Felsenstein J (1985). Evolution 39(4):783-791. doi:10.1111/j.1558-5646.1985.tb00420.x
//         Lemoine et al. (2018). Nature 556:452-456. PMC6030568.
//         Biopython Bio.Phylo.Consensus (bootstrap, get_support).

using NUnit.Framework;
using Seqeron.Genomics.Phylogenetics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for PHYLO-BOOT-001: Phylogenetic Bootstrap Analysis (Felsenstein's bootstrap proportions).
/// Covers <see cref="PhylogeneticAnalyzer.Bootstrap"/>. Randomized resampling is made deterministic
/// via a fixed, documented seed (default 42; see Evidence Dataset 1).
/// </summary>
[TestFixture]
[Category("PHYLO-BOOT-001")]
public class PhylogeneticAnalyzer_Bootstrap_Tests
{
    // Fixed documented seed: makes the column-resampling reproducible (Evidence §Test Datasets).
    private const int Seed = 42;

    // Two well-separated groups: A=B (all A), C=D (all G). For any resampled multiset of these
    // columns, d(A,B)=d(C,D)=0 and the A/B-vs-C/D distances stay saturated, so every replicate
    // recovers the same {A,B},{C,D} topology (Evidence Dataset 1).
    private static Dictionary<string, string> TwoGroupAlignment() => new()
    {
        ["A"] = "AAAAAAAAAA",
        ["B"] = "AAAAAAAAAA",
        ["C"] = "GGGGGGGGGG",
        ["D"] = "GGGGGGGGGG",
    };

    #region Bootstrap — support values (M1, M2, M3, M6)

    [Test]
    [Description("M1: a clade recovered in every replicate has support exactly 1.0 (Felsenstein 1985: 100% -> P=1).")]
    public void Bootstrap_TwoSeparatedGroups_SupportIsExactlyOne()
    {
        // Arrange
        var sequences = TwoGroupAlignment();

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(
            sequences, replicates: 100, distanceMethod: PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.UPGMA, seed: Seed);

        // Assert — clade keys are sorted, '|'-joined leaf names (GetClades convention).
        Assert.Multiple(() =>
        {
            Assert.That(support.ContainsKey("A|B"), Is.True, "Reference clade {A,B} must be scored.");
            Assert.That(support.ContainsKey("C|D"), Is.True, "Reference clade {C,D} must be scored.");
            Assert.That(support["A|B"], Is.EqualTo(1.0).Within(1e-10),
                "{A,B} is recovered in every replicate (distances invariant under column resampling).");
            Assert.That(support["C|D"], Is.EqualTo(1.0).Within(1e-10),
                "{C,D} is recovered in every replicate (distances invariant under column resampling).");
        });
    }

    [Test]
    [Description("M2 (INV-1): every reported support value lies in [0,1] (support = count/replicates).")]
    public void Bootstrap_MixedAlignment_AllSupportsInUnitInterval()
    {
        // Arrange — partly informative alignment so not every clade is 0 or 1.
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAA",
            ["C"] = "TGCATGCATG",
            ["D"] = "TGCATGCATT",
        };

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 50, seed: Seed);

        // Assert
        Assert.That(support, Is.Not.Empty, "Reference tree must define at least one non-trivial clade.");
        Assert.That(support.Values, Has.All.InRange(0.0, 1.0),
            "Support is a proportion count/replicates, hence within [0,1] (INV-1).");
    }

    [Test]
    [Description("M3 (INV-2): every support value equals k/replicates for some integer k in [0, replicates].")]
    public void Bootstrap_KnownReplicateCount_SupportsAreQuantizedToCountOverReplicates()
    {
        // Arrange
        const int replicates = 20;
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAA",
            ["C"] = "TGCATGCATG",
            ["D"] = "TGCATGCATT",
        };

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: replicates, seed: Seed);

        // Assert — value * replicates must be a non-negative integer not exceeding replicates.
        Assert.Multiple(() =>
        {
            foreach (var kvp in support)
            {
                double scaled = kvp.Value * replicates;
                double rounded = Math.Round(scaled);
                Assert.That(scaled, Is.EqualTo(rounded).Within(1e-9),
                    $"Support for {kvp.Key} must be an integer count over {replicates} replicates (INV-2).");
                Assert.That(rounded, Is.InRange(0.0, replicates),
                    $"Count for {kvp.Key} must be within [0,{replicates}].");
            }
        });
    }

    [Test]
    [Description("M6 (INV-5): all-identical sequences -> every reported clade has support 1.0 (one topology each replicate).")]
    public void Bootstrap_AllIdenticalSequences_EveryCladeHasFullSupport()
    {
        // Arrange — zero-distance matrix reproduces an identical topology in every replicate.
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "ACGTACGT",
            ["C"] = "ACGTACGT",
        };

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 30, seed: Seed);

        // Assert
        Assert.That(support, Is.Not.Empty, "A 3-taxon tree has one non-trivial clade.");
        Assert.That(support.Values, Has.All.EqualTo(1.0).Within(1e-10),
            "Identical sequences give identical reference and replicate topologies (INV-5).");
    }

    [Test]
    [Description("M7 (N-ary NJ): NeighborJoining infers an UNROOTED tree whose centre is a trifurcation (Saitou & Nei 1987). For the 4-taxon two-group dataset the NJ root is the trifurcation ((A,B),C,D): {A,B} is the single non-trivial rooted CLADE (recovered in every replicate, support 1.0), while {C,D} sits across the unrooted trifurcation and is therefore NOT a rooted clade. The bootstrap proportion procedure is tree-method-agnostic (Felsenstein 1985).")]
    public void Bootstrap_NeighborJoining_TwoSeparatedGroups_SupportIsExactlyOne()
    {
        // Arrange — same two-group dataset, exercising the NeighborJoining branch.
        // d(A,B)=d(C,D)=0; the cross distances are saturated, so every NJ replicate recovers the
        // same unrooted topology ((A,B),C,D). Hand-derived: with the N-ary model NJ stops at three
        // OTUs {AB, C, D} and connects them to one central node — a trifurcation, NOT a binary root.
        // Under a rooted-clade reading the only non-trivial clade is {A,B}; C and D hang directly off
        // the trifurcation so no internal node groups exactly {C,D}.
        var sequences = TwoGroupAlignment();

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(
            sequences, replicates: 50, distanceMethod: PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            treeMethod: PhylogeneticAnalyzer.TreeMethod.NeighborJoining, seed: Seed);

        // Assert
        Assert.Multiple(() =>
        {
            // Exactly one non-trivial clade for the trifurcating root.
            Assert.That(support.Keys, Is.EquivalentTo(new[] { "A|B" }),
                "NJ trifurcation ((A,B),C,D) has exactly one non-trivial rooted clade: {A,B}.");
            Assert.That(support["A|B"], Is.EqualTo(1.0).Within(1e-10),
                "{A,B} is recovered in every NJ replicate (distances invariant under column resampling).");
            Assert.That(support.ContainsKey("C|D"), Is.False,
                "{C,D} crosses the unrooted NJ trifurcation, so it is not a rooted clade.");
        });
    }

    #endregion

    #region Bootstrap — structure and determinism (M4, M5, S4)

    [Test]
    [Description("M4 (INV-3): result keys equal exactly the non-trivial clades of the reference tree.")]
    public void Bootstrap_ResultKeys_MatchReferenceTreeNonTrivialClades()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAA",
            ["C"] = "GGGGGGGGGG",
            ["D"] = "GGGGGGGGGG",
        };

        // Act
        var support = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 10, seed: Seed);

        // Reference (original-data) tree's non-trivial clades, computed the same way Bootstrap does.
        var refTree = PhylogeneticAnalyzer.BuildTree(sequences);
        var expectedClades = NonTrivialCladeKeys(refTree.Root);

        // Assert
        Assert.That(support.Keys.OrderBy(k => k), Is.EqualTo(expectedClades.OrderBy(k => k)),
            "Bootstrap scores exactly the non-trivial clades of the reference tree (INV-3).");
    }

    [Test]
    [Description("M5 (INV-4): identical inputs and seed produce identical results across runs.")]
    public void Bootstrap_SameSeed_IsDeterministic()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTACGTAA",
            ["C"] = "TGCATGCATG",
            ["D"] = "TGCATGCATT",
        };

        // Act
        var first = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 25, seed: Seed);
        var second = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 25, seed: Seed);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(second.Keys.OrderBy(k => k), Is.EqualTo(first.Keys.OrderBy(k => k)),
                "Same seed must yield the same clade keys (INV-4).");
            foreach (var key in first.Keys)
                Assert.That(second[key], Is.EqualTo(first[key]).Within(1e-12),
                    $"Same seed must yield identical support for {key} (INV-4).");
        });
    }

    [Test]
    [Description("S4: different seeds may differ but every result still satisfies the [0,1] / quantization invariants.")]
    public void Bootstrap_DifferentSeeds_RemainValidProportions()
    {
        // Arrange
        const int replicates = 40;
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTAC",
            ["B"] = "ACGTTCGTAA",
            ["C"] = "TGCATGCATG",
            ["D"] = "TGCATCCATT",
        };

        // Act
        var s1 = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: replicates, seed: 1);
        var s7 = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: replicates, seed: 7);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(s1.Values, Has.All.InRange(0.0, 1.0), "Seed 1 results must be valid proportions.");
            Assert.That(s7.Values, Has.All.InRange(0.0, 1.0), "Seed 7 results must be valid proportions.");
        });
    }

    #endregion

    #region Bootstrap — input validation (S1, S2, S3)

    [Test]
    [Description("S1: null sequences throw ArgumentNullException.")]
    public void Bootstrap_NullSequences_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => PhylogeneticAnalyzer.Bootstrap(null!),
            "Null input is invalid.");
    }

    [Test]
    [Description("S2: fewer than 2 sequences throw ArgumentException (a tree needs >=2 taxa).")]
    public void Bootstrap_SingleSequence_Throws()
    {
        var sequences = new Dictionary<string, string> { ["A"] = "ACGT" };
        Assert.Throws<ArgumentException>(() => PhylogeneticAnalyzer.Bootstrap(sequences),
            "At least 2 sequences are required to build a tree.");
    }

    [Test]
    [Description("S3: replicates < 1 throws ArgumentException (the denominator must be >=1).")]
    public void Bootstrap_ZeroReplicates_Throws()
    {
        var sequences = TwoGroupAlignment();
        Assert.Throws<ArgumentException>(() => PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 0),
            "At least 1 replicate is required.");
    }

    [Test]
    [Description("S3b: negative replicates throw ArgumentException (boundary below the >=1 contract).")]
    public void Bootstrap_NegativeReplicates_Throws()
    {
        var sequences = TwoGroupAlignment();
        Assert.Throws<ArgumentException>(() => PhylogeneticAnalyzer.Bootstrap(sequences, replicates: -5),
            "A negative replicate count is invalid (denominator must be >=1).");
    }

    [Test]
    [Description("S5: unequal-length sequences throw ArgumentException (bootstrap resamples columns of an alignment; Felsenstein 1985 requires aligned characters). Contract: surfaces from BuildTree.")]
    public void Bootstrap_UnequalLengthSequences_Throws()
    {
        // Bootstrap resamples alignment columns; sequences of different lengths are not an alignment.
        var sequences = new Dictionary<string, string> { ["A"] = "ACGT", ["B"] = "ACG" };
        Assert.Throws<ArgumentException>(() => PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 10),
            "Sequences must be aligned (equal length) to resample columns.");
    }

    #endregion

    #region Helpers

    // Mirrors PhylogeneticAnalyzer.GetClades/CollectClades: non-trivial clade = sorted '|'-joined
    // leaf names of a subtree with >1 and < all leaves.
    private static HashSet<string> NonTrivialCladeKeys(PhylogeneticAnalyzer.PhyloNode root)
    {
        var clades = new HashSet<string>();
        int total = PhylogeneticAnalyzer.GetLeaves(root).Count();
        Collect(root, clades, total);
        return clades;

        static List<string> Collect(PhylogeneticAnalyzer.PhyloNode? node, HashSet<string> acc, int total)
        {
            if (node == null) return new List<string>();
            if (node.IsLeaf) return new List<string> { node.Name };
            var taxa = Collect(node.Left, acc, total)
                .Concat(Collect(node.Right, acc, total))
                .OrderBy(n => n).ToList();
            if (taxa.Count > 1 && taxa.Count < total)
                acc.Add(string.Join("|", taxa));
            return taxa;
        }
    }

    #endregion
}
