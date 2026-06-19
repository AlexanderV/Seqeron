using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Phylogenetics;
using static Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Phylogenetic area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-DIST-001 — pairwise evolutionary distance (Phylogenetic).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 39.
///
/// API under test (PhylogeneticAnalyzer.CalculatePairwiseDistance):
///   Over the columns where both aligned sequences carry a standard base (A/C/G/T; gaps and
///   ambiguous bases skipped), let p = differences / comparableSites. Then
///     Hamming         = raw difference count
///     PDistance       = p
///     JukesCantor     = −¾·ln(1 − 4p/3)                          (JC69)
///     Kimura2Parameter= −½·ln((1 − 2S − V)·√(1 − 2V))            (K80; S,V = transition,
///                                                                  transversion proportions)
///   comparableSites = 0 ⇒ distance 0.
///
/// Relations (derived from these formulas, NOT from output):
///   • INV (zero self-distance): every column of d(x,x) is identical, so differences = 0 and
///          p = 0; each formula maps p = 0 (and S = V = 0) to exactly 0.
///   • SYM (symmetry): difference, transition and transversion counts are symmetric in the two
///          sequences, so d(a,b) = d(b,a) for every method.
///   • MON (more mutations ⇒ larger distance): adding one substitution raises the difference
///          count (and S or V), and each method is strictly increasing in those — so the
///          distance strictly increases.
///   • COMP (triangle inequality): Hamming and PDistance are genuine metrics on equal-length
///          gap-free sequences (a≠c at a site implies a≠b or b≠c there), so
///          d(a,c) ≤ d(a,b) + d(b,c). The corrected distances (JukesCantor, Kimura2Parameter)
///          are CONVEX in p and designed for additivity along a tree, NOT metricity, and can
///          violate the triangle inequality — so this relation is asserted only for the two
///          raw metrics, rather than rubber-stamping a property the corrections lack.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class PhylogeneticMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>A base guaranteed to differ from <paramref name="c"/>, so flipping yields a difference.</summary>
    private static char Flip(char c) => c == 'A' ? 'C' : 'A';

    /// <summary>Returns a copy of <paramref name="seq"/> with exactly <paramref name="count"/> spread-out positions mutated.</summary>
    private static string MutatePositions(string seq, int count)
    {
        var chars = seq.ToCharArray();
        for (int n = 0; n < count; n++)
        {
            int pos = (int)((long)n * seq.Length / count);
            chars[pos] = Flip(chars[pos]);
        }
        return new string(chars);
    }

    private static readonly DistanceMethod[] AllMethods =
    {
        DistanceMethod.PDistance,
        DistanceMethod.JukesCantor,
        DistanceMethod.Kimura2Parameter,
        DistanceMethod.Hamming,
    };

    /// <summary>Methods that are genuine metrics on equal-length gap-free sequences.</summary>
    private static readonly DistanceMethod[] MetricMethods =
    {
        DistanceMethod.PDistance,
        DistanceMethod.Hamming,
    };

    /// <summary>Reorders taxa and the symmetric distance matrix by the given permutation, consistently.</summary>
    private static (List<string> Taxa, double[,] Matrix) Permute(
        IReadOnlyList<string> taxa, double[,] matrix, int[] perm)
    {
        int n = perm.Length;
        var newTaxa = perm.Select(p => taxa[p]).ToList();
        var newMatrix = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                newMatrix[i, j] = matrix[perm[i], perm[j]];
        return (newTaxa, newMatrix);
    }

    #endregion

    #region INV — d(x,x) = 0

    [Test]
    [Description("INV: a sequence has zero distance to itself under every method, because no column differs.")]
    public void PairwiseDistance_SelfDistance_IsZero()
    {
        foreach (var method in AllMethods)
        {
            foreach (int len in new[] { 10, 40, 100 })
            {
                string x = RandomDna(len);
                PhylogeneticAnalyzer.CalculatePairwiseDistance(x, x, method)
                    .Should().BeApproximately(0.0, 1e-12,
                        because: $"{method}: identical sequences differ at no comparable site, so p = 0 maps to distance 0");
            }
        }
    }

    #endregion

    #region SYM — d(a,b) = d(b,a)

    [Test]
    [Description("SYM: distance is symmetric, since difference/transition/transversion counts do not depend on argument order.")]
    public void PairwiseDistance_IsSymmetric()
    {
        foreach (var method in AllMethods)
        {
            for (int trial = 0; trial < 20; trial++)
            {
                int len = 60;
                string a = RandomDna(len);
                string b = MutatePositions(a, Rng.Next(1, 10));   // a few substitutions

                double ab = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
                double ba = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, a, method);

                ba.Should().BeApproximately(ab, 1e-12,
                    because: $"{method}: swapping the arguments leaves the per-column comparison counts unchanged");
            }
        }
    }

    #endregion

    #region MON — more substitutions strictly increase the distance

    [Test]
    [Description("MON: each additional substitution raises the difference count, and every method is strictly increasing in it, so the distance strictly increases.")]
    public void PairwiseDistance_MoreMutations_StrictlyIncreases()
    {
        const int len = 40;   // keep mutation fraction well below the JC/K80 divergence limits

        foreach (var method in AllMethods)
        {
            string baseSeq = RandomDna(len);
            double previous = double.NegativeInfinity;

            for (int d = 0; d <= 6; d++)
            {
                string variant = MutatePositions(baseSeq, d);
                double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(baseSeq, variant, method);

                dist.Should().BeGreaterThan(previous,
                    because: $"{method}: {d} substitutions differ at one more site than {d - 1}, and the distance is strictly increasing in the difference count");
                previous = dist;
            }
        }
    }

    #endregion

    #region COMP — triangle inequality for the raw metrics

    [Test]
    [Description("COMP: on equal-length gap-free sequences Hamming and PDistance satisfy the triangle inequality d(a,c) ≤ d(a,b) + d(b,c).")]
    public void PairwiseDistance_RawMetrics_SatisfyTriangleInequality()
    {
        foreach (var method in MetricMethods)
        {
            for (int trial = 0; trial < 50; trial++)
            {
                int len = 50;
                string a = RandomDna(len);
                string b = RandomDna(len);
                string c = RandomDna(len);

                double ab = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
                double bc = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, c, method);
                double ac = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, c, method);

                ac.Should().BeLessThanOrEqualTo(ab + bc + 1e-9,
                    because: $"{method}: a≠c at a site forces a≠b or b≠c there, so per-site disagreements are subadditive");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PHYLO-TREE-001 — UPGMA tree construction (Phylogenetic).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 40.
    //
    // API under test (PhylogeneticAnalyzer.BuildTreeFromMatrix, TreeMethod.UPGMA):
    //   UPGMA repeatedly merges the closest pair of clusters; the new cluster's height is
    //   d(i,j)/2 and each child's branch length is height(new) − height(child). The result is
    //   an ultrametric rooted tree. PatristicDistance sums the branch lengths between two
    //   leaves; RobinsonFouldsDistance is the rooted-clade symmetric difference (0 ⇔ identical
    //   topology on the same taxa).
    //
    // Relations (derived from the UPGMA recurrence, NOT from output):
    //   • INV (input order ⇒ same topology): UPGMA depends only on the SET of pairwise
    //          distances, not the taxon ordering. With all distances distinct (no merge ties),
    //          the rooted topology is uniquely determined, so permuting the taxa+matrix yields
    //          an RF distance of 0 and an unchanged total tree length.
    //   • MON (closer ⇒ shorter branches): the first-merged cherry {X,Y} reproduces its input
    //          distance exactly as a patristic distance (d/2 + d/2 = d), so making X,Y closer
    //          shortens their branches; and the cherry's patristic distance stays below the
    //          patristic distance to any outgroup taxon.
    // ───────────────────────────────────────────────────────────────────────────

    #region INV — permuting the input order preserves the UPGMA topology

    [Test]
    [Description("INV: UPGMA depends only on the set of distances, so permuting the taxa+matrix gives an identical rooted topology (RF = 0) and the same total tree length, when distances are distinct.")]
    public void BuildUpgma_PermuteInputOrder_PreservesTopology()
    {
        var taxa = new List<string> { "A", "B", "C", "D" };
        // All pairwise distances distinct ⇒ no merge ties ⇒ unique UPGMA topology ((A,B),(C,D)).
        var matrix = new double[,]
        {
            { 0,  2, 10, 11 },
            { 2,  0, 12, 13 },
            { 10, 12, 0,  4 },
            { 11, 13, 4,  0 },
        };

        var reference = PhylogeneticAnalyzer.BuildTreeFromMatrix(taxa, matrix, TreeMethod.UPGMA);

        foreach (var perm in new[]
                 {
                     new[] { 3, 2, 1, 0 },
                     new[] { 1, 3, 0, 2 },
                     new[] { 2, 0, 3, 1 },
                 })
        {
            var (pTaxa, pMatrix) = Permute(taxa, matrix, perm);
            var permuted = PhylogeneticAnalyzer.BuildTreeFromMatrix(pTaxa, pMatrix, TreeMethod.UPGMA);

            PhylogeneticAnalyzer.RobinsonFouldsDistance(reference.Root, permuted.Root).Should().Be(0,
                because: "UPGMA uses only the (distinct) pairwise distances, so reordering the taxa cannot change which clusters merge — the rooted topology is identical");
            PhylogeneticAnalyzer.CalculateTreeLength(permuted.Root).Should().BeApproximately(
                PhylogeneticAnalyzer.CalculateTreeLength(reference.Root), 1e-9,
                because: "the same merges at the same heights produce the same branch lengths regardless of input order");
        }
    }

    #endregion

    #region MON — a closer cherry has shorter branches than a more distant one

    [Test]
    [Description("MON: in UPGMA the first-merged cherry's patristic distance equals its input distance, so a smaller input distance gives strictly shorter branches — and stays below the distance to the outgroup.")]
    public void BuildUpgma_CloserPair_HasShorterBranches()
    {
        const double outgroup = 10.0;   // fixed distance from each of X,Y to Z
        var taxa = new List<string> { "X", "Y", "Z" };

        double previous = double.NegativeInfinity;

        foreach (double dxy in new[] { 1.0, 2.0, 4.0, 6.0, 8.0 })
        {
            var matrix = new double[,]
            {
                { 0,        dxy,      outgroup },
                { dxy,      0,        outgroup },
                { outgroup, outgroup, 0 },
            };

            var tree = PhylogeneticAnalyzer.BuildTreeFromMatrix(taxa, matrix, TreeMethod.UPGMA);

            double cherry = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "X", "Y");
            double toOut = PhylogeneticAnalyzer.PatristicDistance(tree.Root, "X", "Z");

            cherry.Should().BeApproximately(dxy, 1e-9,
                because: "UPGMA merges X,Y first at height dxy/2, so their patristic distance reconstructs the input distance dxy exactly");
            cherry.Should().BeGreaterThan(previous,
                because: "a larger input distance between the cherry taxa gives proportionally longer branches");
            cherry.Should().BeLessThan(toOut,
                because: "the closest pair sits in the lowest cherry, so its patristic distance is below the distance to the outgroup");
            previous = cherry;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PHYLO-NEWICK-001 — Newick serialization round-trip (Phylogenetic).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 41.
    //
    // API under test (PhylogeneticAnalyzer.ToNewick / ParseNewick):
    //   ToNewick emits the tree in Newick format with child branch lengths formatted "F4";
    //   internal labels containing Newick metacharacters (e.g. auto-generated UPGMA names) are
    //   omitted per the Olsen spec. ParseNewick reads it back, trimming leading/trailing
    //   whitespace and the terminal ';'.
    //
    // Relations (derived from the grammar, NOT from output):
    //   • COMP (round-trip identity): parsing a serialized tree and re-serializing reproduces
    //          the string exactly (re-serialization is idempotent): leaf names, child order and
    //          the F4 branch lengths all survive a parse∘serialize cycle, and the reparsed tree
    //          has the identical topology (RF = 0) and leaf set.
    //   • INV (surrounding whitespace irrelevant): ParseNewick trims leading/trailing
    //          whitespace, so padding the string with spaces/tabs/newlines yields the same tree.
    //          NOTE: this parser does NOT ignore INTERNAL whitespace — it would become part of a
    //          label — and the Olsen spec forbids blanks in unquoted labels, so only SURROUNDING
    //          whitespace invariance is asserted, not a false "all whitespace" claim.
    // ───────────────────────────────────────────────────────────────────────────

    #region Helpers (Newick)

    private static PhylogeneticTree SampleTree(int which) => which switch
    {
        0 => PhylogeneticAnalyzer.BuildTreeFromMatrix(
                new List<string> { "A", "B", "C", "D" },
                new double[,] { { 0, 2, 10, 11 }, { 2, 0, 12, 13 }, { 10, 12, 0, 4 }, { 11, 13, 4, 0 } },
                TreeMethod.UPGMA),
        _ => PhylogeneticAnalyzer.BuildTreeFromMatrix(
                new List<string> { "Homo", "Pan", "Gorilla", "Pongo", "Macaca" },
                new double[,]
                {
                    { 0, 3, 7, 13, 19 },
                    { 3, 0, 8, 14, 20 },
                    { 7, 8, 0, 15, 21 },
                    { 13, 14, 15, 0, 22 },
                    { 19, 20, 21, 22, 0 },
                },
                TreeMethod.UPGMA),
    };

    private static List<string> SortedLeafNames(PhyloNode root) =>
        PhylogeneticAnalyzer.GetLeaves(root).Select(l => l.Name).OrderBy(s => s).ToList();

    /// <summary>Builds a right-nested "ladder" rooted tree from a leaf order, e.g. (1,(2,(3,(4,(5,6))))).</summary>
    private static PhyloNode Ladder(params string[] order)
    {
        string s = order[^1];
        for (int i = order.Length - 2; i >= 0; i--)
            s = $"({order[i]},{s})";
        return PhylogeneticAnalyzer.ParseNewick(s + ";");
    }

    #endregion

    #region COMP — parse∘serialize round-trips the tree

    [Test]
    [Description("COMP: ParseNewick(ToNewick(t)) reproduces the tree — re-serialization is idempotent, the topology is preserved (RF = 0) and the leaf set is unchanged.")]
    public void Newick_RoundTrip_PreservesTree()
    {
        foreach (int which in new[] { 0, 1 })
        {
            var root = SampleTree(which).Root;

            string newick = PhylogeneticAnalyzer.ToNewick(root);
            var reparsed = PhylogeneticAnalyzer.ParseNewick(newick);

            PhylogeneticAnalyzer.ToNewick(reparsed).Should().Be(newick,
                because: "leaf names, child order and the F4-formatted branch lengths all survive a parse∘serialize cycle, so re-serialization is idempotent");
            PhylogeneticAnalyzer.RobinsonFouldsDistance(root, reparsed).Should().Be(0,
                because: "the round-trip must preserve the tree topology exactly");
            SortedLeafNames(reparsed).Should().Equal(SortedLeafNames(root),
                because: "every taxon must reappear after the round-trip");
        }
    }

    #endregion

    #region INV — surrounding whitespace does not affect the parse

    [Test]
    [Description("INV: padding a Newick string with leading/trailing whitespace yields the same tree, because ParseNewick trims surrounding whitespace.")]
    public void Newick_SurroundingWhitespace_DoesNotAffectParse()
    {
        foreach (int which in new[] { 0, 1 })
        {
            string newick = PhylogeneticAnalyzer.ToNewick(SampleTree(which).Root);
            string canonical = PhylogeneticAnalyzer.ToNewick(PhylogeneticAnalyzer.ParseNewick(newick));

            foreach (string padded in new[] { "  " + newick, newick + "\n", "\t " + newick + "  \r\n" })
            {
                string reSerialized = PhylogeneticAnalyzer.ToNewick(PhylogeneticAnalyzer.ParseNewick(padded));
                reSerialized.Should().Be(canonical,
                    because: "ParseNewick trims leading/trailing whitespace, so surrounding padding cannot change the parsed tree");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PHYLO-COMP-001 — tree comparison / Robinson–Foulds distance (Phylogenetic).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 42.
    //
    // API under test (PhylogeneticAnalyzer.RobinsonFouldsDistance):
    //   RF(t1,t2) = the symmetric difference of the two trees' sets of NON-TRIVIAL rooted
    //   clades (subtree leaf-sets of size 2..n−1). For two fully-resolved trees with the same
    //   clade count c, RF = 2·(c − shared).
    //
    // Relations (derived from the symmetric-difference definition, NOT from output):
    //   • SYM (symmetry): the symmetric difference |A\B| + |B\A| is symmetric, so RF(a,b) = RF(b,a).
    //   • COMP (self-distance zero): a tree has the same clade set as itself, so RF(t,t) = 0.
    //   • MON (more rearrangements ⇒ higher RF): on a ladder, the non-trivial clades are the
    //          nested suffix-sets of the leaf order. The σ_k family below corrupts the leaf
    //          order so that σ_k keeps exactly the (4−k) largest suffix-clades of the reference
    //          and breaks the rest, giving RF(R, σ_k) = 2k — strictly increasing with the number
    //          of rearrangements. (Hand-verified suffix sets; see σ_k construction.)
    // ───────────────────────────────────────────────────────────────────────────

    #region SYM — RF(a,b) = RF(b,a)

    [Test]
    [Description("SYM: Robinson–Foulds distance is symmetric, being a symmetric set difference of clades.")]
    public void RobinsonFoulds_IsSymmetric()
    {
        var trees = new[]
        {
            Ladder("1", "2", "3", "4", "5", "6"),
            Ladder("1", "2", "3", "5", "6", "4"),
            Ladder("6", "5", "4", "3", "2", "1"),
            Ladder("2", "4", "6", "1", "3", "5"),
        };

        for (int i = 0; i < trees.Length; i++)
            for (int j = 0; j < trees.Length; j++)
                PhylogeneticAnalyzer.RobinsonFouldsDistance(trees[i], trees[j])
                    .Should().Be(PhylogeneticAnalyzer.RobinsonFouldsDistance(trees[j], trees[i]),
                        because: "the symmetric difference of the two clade sets does not depend on argument order");
    }

    #endregion

    #region COMP — RF(t,t) = 0

    [Test]
    [Description("COMP: a tree has Robinson–Foulds distance 0 to itself, having an identical clade set.")]
    public void RobinsonFoulds_SelfDistance_IsZero()
    {
        var trees = new[]
        {
            Ladder("1", "2", "3", "4", "5", "6"),
            Ladder("6", "5", "4", "3", "2", "1"),
            SampleTree(0).Root,
            SampleTree(1).Root,
        };

        foreach (var t in trees)
            PhylogeneticAnalyzer.RobinsonFouldsDistance(t, t).Should().Be(0,
                because: "a tree shares every clade with itself, so the symmetric difference is empty");
    }

    #endregion

    #region MON — more rearrangements give a strictly larger RF distance

    [Test]
    [Description("MON: progressively corrupting more of a ladder's suffix-clades strictly increases the RF distance from the reference (RF(R, σ_k) = 2k).")]
    public void RobinsonFoulds_MoreRearrangements_IncreasesDistance()
    {
        var reference = Ladder("1", "2", "3", "4", "5", "6");

        // σ_k keeps the (4−k) largest suffix-clades of R and breaks the rest ⇒ RF = 2k.
        var family = new[]
        {
            Ladder("1", "2", "3", "4", "5", "6"),   // k=0: shares all 4 clades  → RF 0
            Ladder("1", "2", "3", "5", "6", "4"),   // k=1: shares 3 → RF 2
            Ladder("1", "2", "4", "5", "6", "3"),   // k=2: shares 2 → RF 4
            Ladder("1", "3", "4", "5", "6", "2"),   // k=3: shares 1 → RF 6
            Ladder("2", "3", "4", "5", "6", "1"),   // k=4: shares 0 → RF 8
        };

        int previous = -1;
        for (int k = 0; k < family.Length; k++)
        {
            int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(reference, family[k]);

            rf.Should().Be(2 * k,
                because: $"σ_{k} retains the {4 - k} largest suffix-clades of the reference and breaks the rest, so RF = 2·{k}");
            rf.Should().BeGreaterThan(previous,
                because: "each further rearrangement breaks one more shared clade, strictly increasing the Robinson–Foulds distance");
            previous = rf;
        }
    }

    #endregion
}
