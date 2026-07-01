using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Phylogenetics;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for MIRNA-PCT-001 — the Friedman et al. (2009) / TargetScan PCT
/// (<i>probability of conserved targeting</i>) branch-length conservation surface.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// the code NEVER fails in an undisciplined way: no hang, no unhandled runtime
/// exception, no out-of-contract output. PCT is a logistic (sigmoid) of a branch-length
/// score, so the headline hazards are: (a) <c>Math.Exp</c> OVERFLOW producing ±Inf/NaN
/// from extreme params or extreme branch lengths, and (b) a PCT value escaping its
/// theoretical band. Every input must resolve to EITHER a well-defined, theory-correct
/// finite result OR a documented validation exception. — §8.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-PCT-001 — PCT from branch-length conservation
/// Checklist: docs/checklists/03_FUZZING.md, row 253.
/// Strategies (per §Description): MC = Malformed Content, BE = Boundary Exploitation.
/// Fuzz targets for row 253: "empty tree, single-species alignment, negative branch
/// length, missing sigmoid params".
/// Source doc: docs/algorithms/MiRNA/Target_Site_Prediction.md
///   (Bls definition + the published logistic + the truncate-at-0 rule).
/// Source: src/.../Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs
///   • PctFromBranchLength(double bls, PctSigmoidParameters p)  (~line 1645)
///   • ComputeBranchLengthScore(PhyloNode tree, IReadOnlyCollection&lt;string&gt;)  (~1574)
///   • record struct PctSigmoidParameters(B0,B1,B2,B3)  (~589)
///   • record struct PctConservation(Tree, SpeciesWithSite, SigmoidParameters)  (~601)
///
/// ───────────────────────────────────────────────────────────────────────────
/// The PCT contract under test (independently derived from the doc, NOT read off code)
/// ───────────────────────────────────────────────────────────────────────────
/// The library maps a branch-length score b (Bls) to PCT via the published TargetScan
/// logistic (targetscan_70_BL_PCT.pl, calculatePCTthisBL):
///
///     PCT(b) = B0 + B1 / (1 + e^(−B2·b + B3)),   then truncated at 0 (negative ⇒ 0).
///
/// Pinned invariants (theory, re-derived here):
///   INV-PCT-1  (FLOOR / saturation): the logistic term B1/(1+e^…) ∈ (0, B1) for any
///              finite b, so PCT ∈ (B0, B0+B1) before truncation. As b→+∞ the term →B1
///              (PCT→B0+B1, the ceiling); as b→−∞ it →0 (PCT→B0, the floor). For the
///              PUBLISHED-STYLE band B0≥0, B1>0, B0+B1≤1 this forces PCT ∈ [0,1] for
///              ALL b — including the degenerate b=0 (empty / single-species tree), the
///              boundary b=0, and NEGATIVE b.
///   INV-PCT-2  (MONOTONIC): for B1>0, B2>0, PCT is non-decreasing in b (dPCT/db =
///              B1·B2·e^(…)/(1+e^(…))² ≥ 0). More conservation ⇒ never-lower PCT.
///   INV-PCT-3  (HAND PIN): for B0=0,B1=1,B2=1,B3=0 → PCT(0)=1/(1+e^0)=0.5 exactly;
///              PCT(2)=1/(1+e^−2)=0.880797…  (computed independently below).
///   INV-PCT-4  (FINITE / no overflow escape): with FINITE params and ANY finite b
///              (incl. ±1e300 and huge B2) the result is finite and ≥0 — Math.Exp may
///              return +Inf, but B1/(1+Inf)=0 saturates to the floor B0, NOT NaN/Inf.
///   INV-PCT-5  (Bls semantics, Friedman 2009): Bls = total branch length of the minimal
///              subtree connecting the species-with-site. EMPTY species set or a SINGLE
///              species ⇒ Bls = 0 (no connecting path). NEGATIVE edge lengths flow
///              through additively (the score is a plain sum), and the resulting b feeds
///              the same logistic — so even a negative Bls yields a finite PCT≥0.
///
/// MC/BE boundaries fuzzed: empty tree (single leaf), single-species alignment (Bls=0),
/// negative branch length (negative b into the logistic), all-zero / degenerate /
/// extreme sigmoid params, b = 0 / ±huge, and a hand-traced end-to-end Bls→PCT pin.
/// If a source-derived assertion and the code disagree, the CODE is wrong (fix per doc);
/// a PCT &lt; 0, or any NaN/Inf from exp overflow, would be a real bug.
///
/// LimitationPolicy: these are pure static functions; the assembly module-initializer
/// (_LimitationPolicyTestBootstrap) already runs under LimitationMode.Permissive.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaPctFuzzTests
{
    private const double Tol = 1e-12;

    #region Helpers

    /// <summary>Independent reference logistic — the doc's formula re-implemented here,
    /// deliberately NOT calling the library, so the pins test theory not tautology.</summary>
    private static double ReferencePct(double b, double b0, double b1, double b2, double b3)
    {
        double pct = b0 + b1 / (1.0 + Math.Exp(-b2 * b + b3));
        return pct < 0.0 ? 0.0 : pct;
    }

    /// <summary>A leaf node with a name and a branch length.</summary>
    private static PhylogeneticAnalyzer.PhyloNode Leaf(string name, double bl) =>
        new(name) { BranchLength = bl };

    /// <summary>An internal node with the given branch length and children.</summary>
    private static PhylogeneticAnalyzer.PhyloNode Inner(double bl, params PhylogeneticAnalyzer.PhyloNode[] kids)
    {
        var n = new PhylogeneticAnalyzer.PhyloNode { BranchLength = bl };
        n.Children.AddRange(kids);
        return n;
    }

    #endregion

    #region MIRNA-PCT-001 — hand-derived sigmoid pin (INV-PCT-3)

    // PCT(0) with B0=0,B1=1,B2=1,B3=0 = 0 + 1/(1+e^0) = 1/(1+1) = 0.5 EXACTLY.
    [Test]
    public void Pct_AtZeroBranchLength_IsTheLogisticFloor_HalfForCanonicalParams()
    {
        var p = new PctSigmoidParameters(0.0, 1.0, 1.0, 0.0);

        double pct = PctFromBranchLength(0.0, p);

        pct.Should().BeApproximately(0.5, Tol,
            "PCT(0)=B0+B1/(1+e^0)=0+1/2 for B0=0,B1=1,B3=0 — the documented b=0 floor");
    }

    // PCT(2) with B0=0,B1=1,B2=1,B3=0 = 1/(1+e^-2) = 0.8807970779778823… (hand value).
    [Test]
    public void Pct_HandDerivedValue_MatchesPublishedLogistic()
    {
        var p = new PctSigmoidParameters(0.0, 1.0, 1.0, 0.0);
        const double expected = 0.8807970779778823; // 1/(1+e^-2), computed independently

        double pct = PctFromBranchLength(2.0, p);

        pct.Should().BeApproximately(expected, Tol,
            "PCT(2)=1/(1+e^-2) for the canonical logistic — pinned from the doc formula, not the code");
        pct.Should().BeApproximately(ReferencePct(2.0, 0, 1, 1, 0), Tol,
            "and must equal the independent reference logistic");
    }

    #endregion

    #region MIRNA-PCT-001 — PCT ∈ [0,1] for published-style params, ALL branch lengths (INV-PCT-1)

    // The headline invariant: across degenerate, boundary, negative and huge b — and a
    // randomised sweep of PUBLISHED-STYLE params (B0≥0, B1>0, B0+B1≤1) — PCT never leaves [0,1].
    [Test]
    public void Pct_PublishedStyleParams_StayWithinUnitInterval_ForEveryBranchLength([Values(20260626, 42, 7)] int seed)
    {
        var rng = new Random(seed);
        // Branch lengths spanning empty/zero, tiny, typical, huge, and NEGATIVE corners.
        double[] branchLengths =
        {
            0.0, 1e-12, 0.5, 1.0, 2.0, 3.5, 12.7, 1e6, 1e300,
            -1e-12, -0.5, -3.0, -1e6, -1e300,
        };

        for (int i = 0; i < 400; i++)
        {
            double b0 = rng.NextDouble() * 0.3;                 // floor ∈ [0,0.3]
            double b1 = rng.NextDouble() * (1.0 - b0);          // ceiling B0+B1 ≤ 1
            double b2 = rng.NextDouble() * 5.0;                 // slope ≥ 0
            double b3 = (rng.NextDouble() - 0.5) * 8.0;         // offset, either sign
            var p = new PctSigmoidParameters(b0, b1, b2, b3);

            foreach (double b in branchLengths)
            {
                double pct = PctFromBranchLength(b, p);

                double.IsNaN(pct).Should().BeFalse($"PCT must not be NaN (b={b}, params={b0},{b1},{b2},{b3})");
                double.IsInfinity(pct).Should().BeFalse($"PCT must not be ±Inf (b={b}, params={b0},{b1},{b2},{b3})");
                pct.Should().BeGreaterThanOrEqualTo(0.0, "the logistic is truncated at 0");
                pct.Should().BeLessThanOrEqualTo(1.0 + Tol,
                    $"PCT ≤ B0+B1 ≤ 1 for published-style params (b={b})");
                pct.Should().BeApproximately(ReferencePct(b, b0, b1, b2, b3), 1e-9,
                    "library PCT must equal the independent reference logistic");
            }
        }
    }

    #endregion

    #region MIRNA-PCT-001 — monotonic non-decreasing in branch length (INV-PCT-2)

    [Test]
    public void Pct_IsMonotonicNonDecreasing_InBranchLength_ForPositiveSlope([Values(1, 99, 20260626)] int seed)
    {
        var rng = new Random(seed);

        for (int i = 0; i < 200; i++)
        {
            double b0 = rng.NextDouble() * 0.3;
            double b1 = 0.05 + rng.NextDouble() * (1.0 - b0 - 0.05); // B1 > 0
            double b2 = 0.05 + rng.NextDouble() * 5.0;               // B2 > 0
            double b3 = (rng.NextDouble() - 0.5) * 6.0;
            var p = new PctSigmoidParameters(b0, b1, b2, b3);

            double prev = double.NegativeInfinity;
            // Increasing branch lengths, including the negative→zero→positive range.
            foreach (double b in new[] { -10.0, -3.0, -1.0, 0.0, 0.5, 1.0, 2.0, 5.0, 20.0, 1e6 })
            {
                double pct = PctFromBranchLength(b, p);
                pct.Should().BeGreaterThanOrEqualTo(prev - 1e-9,
                    $"PCT is non-decreasing in b for B1,B2>0 (b={b}, params={b0},{b1},{b2},{b3})");
                prev = pct;
            }
        }
    }

    #endregion

    #region MIRNA-PCT-001 — MC: degenerate / extreme sigmoid params never overflow (INV-PCT-4)

    // "missing sigmoid params" = all-zero params, and extreme/adversarial params. With FINITE
    // params and ANY finite b the result must be finite and ≥0 (exp overflow saturates, not NaN/Inf).
    [Test]
    public void Pct_DegenerateAndExtremeParams_AlwaysProduceFiniteNonNegativePct()
    {
        var cases = new (double b, PctSigmoidParameters p, string why)[]
        {
            // All-zero params ("missing"): PCT = 0 + 0/(1+e^0) = 0.
            (0.0,   new PctSigmoidParameters(0, 0, 0, 0), "all-zero params ⇒ PCT=0"),
            (5.0,   new PctSigmoidParameters(0, 0, 0, 0), "all-zero params, b≠0 ⇒ still 0"),
            (-5.0,  new PctSigmoidParameters(0, 0, 0, 0), "all-zero params, negative b ⇒ still 0"),
            // Huge slope × huge +b: -B2·b → −∞, exp→0, PCT→B0+B1 (ceiling), finite.
            (1e6,   new PctSigmoidParameters(0.0, 0.9, 1e6, 0), "huge +arg ⇒ ceiling, finite"),
            // Huge slope × huge -b: -B2·b → +∞, exp→+Inf, B1/(1+Inf)=0 ⇒ PCT→B0 (floor), finite.
            (-1e6,  new PctSigmoidParameters(0.2, 0.7, 1e6, 0), "huge -arg ⇒ floor, exp saturates not NaN"),
            (1e300, new PctSigmoidParameters(0.1, 0.8, 1e300, 0), "double-overflow magnitude ⇒ finite floor/ceiling"),
            // Large B3 offset both directions.
            (1.0,   new PctSigmoidParameters(0.0, 1.0, 1.0, 1e3), "huge +B3 ⇒ exp +Inf ⇒ PCT floor B0"),
            (1.0,   new PctSigmoidParameters(0.0, 1.0, 1.0, -1e3), "huge -B3 ⇒ exp 0 ⇒ PCT ceiling B0+B1"),
            // Negative B1 with the logistic ≤1 ⇒ PCT could dip below 0 ⇒ truncated to 0.
            (-50.0, new PctSigmoidParameters(0.0, 1.0, 1.0, 0.0), "very negative b ⇒ near-0, truncated, ≥0"),
            // Negative B0 (out-of-band param) ⇒ truncation-at-0 must still hold the floor.
            (-50.0, new PctSigmoidParameters(-0.5, 0.4, 1.0, 0.0), "negative B0 below 0 ⇒ truncated to 0"),
        };

        foreach (var (b, p, why) in cases)
        {
            double pct = PctFromBranchLength(b, p);

            double.IsNaN(pct).Should().BeFalse($"PCT NaN-free: {why}");
            double.IsInfinity(pct).Should().BeFalse($"PCT finite: {why}");
            pct.Should().BeGreaterThanOrEqualTo(0.0, $"PCT truncated at 0: {why}");
        }
    }

    // Explicit saturation pins: the exp-overflow corners resolve to the documented floor/ceiling.
    [Test]
    public void Pct_ExpOverflowCorners_SaturateToFloorAndCeiling()
    {
        // -B2·b = +∞ ⇒ exp = +Inf ⇒ B1/(1+Inf)=0 ⇒ PCT = B0 (floor).
        PctFromBranchLength(-1e6, new PctSigmoidParameters(0.25, 0.6, 1e6, 0.0))
            .Should().BeApproximately(0.25, 1e-9, "huge +exp-arg saturates to the floor B0");

        // -B2·b = −∞ ⇒ exp = 0 ⇒ B1/1 = B1 ⇒ PCT = B0+B1 (ceiling).
        PctFromBranchLength(1e6, new PctSigmoidParameters(0.25, 0.6, 1e6, 0.0))
            .Should().BeApproximately(0.85, 1e-9, "huge -exp-arg saturates to the ceiling B0+B1");
    }

    #endregion

    #region MIRNA-PCT-001 — BE: Bls boundaries — empty tree, single species, negative edges (INV-PCT-5)

    // Empty tree (a lone leaf) and a single-species alignment both connect 0 or 1 species ⇒ Bls = 0.
    [Test]
    public void Bls_EmptyTreeAndSingleSpecies_AreZero()
    {
        var singleLeaf = Leaf("A", 3.0);

        // No species-with-site at all ⇒ no connecting subtree ⇒ Bls 0.
        ComputeBranchLengthScore(singleLeaf, Array.Empty<string>())
            .Should().Be(0.0, "empty species-with-site set ⇒ Bls 0 (Friedman 2009)");

        // Exactly one species-with-site ⇒ no path between two sites ⇒ Bls 0.
        ComputeBranchLengthScore(singleLeaf, new[] { "A" })
            .Should().Be(0.0, "a single conserved species has no connecting branch ⇒ Bls 0");

        // A richer tree, but only one species marked ⇒ still 0.
        var tree = Inner(0.0, Inner(0.5, Leaf("A", 1.0), Leaf("B", 2.0)),
                              Inner(1.5, Leaf("C", 3.0), Leaf("D", 4.0)));
        ComputeBranchLengthScore(tree, new[] { "C" })
            .Should().Be(0.0, "single-species alignment ⇒ Bls 0");

        // A species not present in the tree at all ⇒ ignored ⇒ Bls 0.
        ComputeBranchLengthScore(tree, new[] { "ZZZ", "QQQ" })
            .Should().Be(0.0, "unknown leaf names contribute nothing ⇒ Bls 0");
    }

    // End-to-end hand-traced Bls→PCT pin:
    // Tree ((A:1.0,B:2.0):0.5,(C:3.0,D:4.0):1.5); species-with-site {A,C}.
    // Connecting subtree edges: A(1.0)+innerAB(0.5)+C(3.0)+innerCD(1.5) = 6.0.
    [Test]
    public void Bls_ThenPct_EndToEnd_MatchesHandTrace()
    {
        var tree = Inner(0.0,
            Inner(0.5, Leaf("A", 1.0), Leaf("B", 2.0)),
            Inner(1.5, Leaf("C", 3.0), Leaf("D", 4.0)));

        double bls = ComputeBranchLengthScore(tree, new[] { "A", "C" });
        bls.Should().BeApproximately(6.0, Tol, "Bls = 1.0+0.5+3.0+1.5 along the A↔C connecting subtree");

        var p = new PctSigmoidParameters(0.0, 0.9, 2.0, 1.0);
        double pct = PctFromBranchLength(bls, p);

        pct.Should().BeApproximately(ReferencePct(6.0, 0.0, 0.9, 2.0, 1.0), Tol,
            "PCT(6.0) must equal the independent reference logistic for these params");
        pct.Should().BeInRange(0.0, 1.0, "published-style params keep PCT in [0,1]");
    }

    // NEGATIVE branch length: the score is a plain additive sum, so a negative edge yields a
    // negative Bls — which is a legal real-number input to the logistic and must NOT crash and
    // must still produce a finite PCT ≥ 0.
    [Test]
    public void Bls_NegativeEdgeLengths_FlowThroughToFinitePct()
    {
        // A and C connected; the inner CD edge is NEGATIVE.
        var tree = Inner(0.0,
            Inner(0.5, Leaf("A", 1.0), Leaf("B", 2.0)),
            Inner(-2.5, Leaf("C", -1.0), Leaf("D", 4.0)));

        double bls = ComputeBranchLengthScore(tree, new[] { "A", "C" });
        // 1.0 + 0.5 + (-1.0) + (-2.5) = -2.0
        bls.Should().BeApproximately(-2.0, Tol, "negative edges sum additively into Bls");

        var p = new PctSigmoidParameters(0.0, 0.9, 2.0, 1.0);
        double pct = PctFromBranchLength(bls, p);

        double.IsNaN(pct).Should().BeFalse("negative Bls must not yield NaN");
        double.IsInfinity(pct).Should().BeFalse("negative Bls must not yield Inf");
        pct.Should().BeGreaterThanOrEqualTo(0.0, "PCT truncated at 0 even for negative Bls");
        pct.Should().BeApproximately(ReferencePct(-2.0, 0.0, 0.9, 2.0, 1.0), Tol,
            "negative Bls feeds the same logistic — pinned to the reference");
    }

    // Randomised malformed trees (random topology, random ± branch lengths, random marked
    // species, random params) must NEVER throw, hang, or emit a non-finite / out-of-floor PCT.
    [Test]
    public void BlsToPct_RandomMalformedTrees_NeverThrow_AndProduceFinitePct([Values(20260626, 13, 555)] int seed)
    {
        var rng = new Random(seed);
        string[] taxa = { "A", "B", "C", "D", "E", "F", "G", "H" };

        for (int iter = 0; iter < 150; iter++)
        {
            // Build a random small bifurcating tree over a random subset of taxa.
            var leaves = taxa.OrderBy(_ => rng.Next()).Take(2 + rng.Next(taxa.Length - 1)).ToList();
            var nodes = leaves.Select(t => Leaf(t, (rng.NextDouble() - 0.3) * 6.0)).ToList<PhylogeneticAnalyzer.PhyloNode>();
            while (nodes.Count > 1)
            {
                var a = nodes[rng.Next(nodes.Count)]; nodes.Remove(a);
                var b = nodes[rng.Next(nodes.Count)]; nodes.Remove(b);
                nodes.Add(Inner((rng.NextDouble() - 0.3) * 6.0, a, b));
            }
            var root = nodes[0];

            // Random marked species (possibly empty, possibly unknown names).
            var marked = taxa.Concat(new[] { "UNKNOWN" })
                             .Where(_ => rng.NextDouble() < 0.5)
                             .ToArray();

            double bls = 0;
            Action computeBls = () => bls = ComputeBranchLengthScore(root, marked);
            computeBls.Should().NotThrow("malformed/degenerate trees must not crash Bls");
            double.IsNaN(bls).Should().BeFalse("Bls must be finite");
            double.IsInfinity(bls).Should().BeFalse("Bls must be finite");

            // Feed Bls through a random published-style sigmoid.
            double b0 = rng.NextDouble() * 0.2;
            double b1 = rng.NextDouble() * (1.0 - b0);
            double b2 = rng.NextDouble() * 4.0;
            double b3 = (rng.NextDouble() - 0.5) * 6.0;
            var p = new PctSigmoidParameters(b0, b1, b2, b3);

            double pct = PctFromBranchLength(bls, p);
            double.IsNaN(pct).Should().BeFalse("PCT must be finite for any tree-derived Bls");
            double.IsInfinity(pct).Should().BeFalse("PCT must be finite for any tree-derived Bls");
            pct.Should().BeGreaterThanOrEqualTo(0.0, "PCT truncated at 0");
            pct.Should().BeLessThanOrEqualTo(1.0 + Tol, "published-style params keep PCT ≤ 1");
        }
    }

    // Null guards are part of the documented contract (ArgumentNullException, not a crash).
    [Test]
    public void Bls_NullArguments_ThrowArgumentNull()
    {
        var leaf = Leaf("A", 1.0);
        Action nullTree = () => ComputeBranchLengthScore(null!, new[] { "A" });
        Action nullSpecies = () => ComputeBranchLengthScore(leaf, null!);

        nullTree.Should().Throw<ArgumentNullException>("a null tree is a documented validation error");
        nullSpecies.Should().Throw<ArgumentNullException>("a null species set is a documented validation error");
    }

    #endregion
}
