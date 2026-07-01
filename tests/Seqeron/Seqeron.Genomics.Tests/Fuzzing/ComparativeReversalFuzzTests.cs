using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Reversal Distance
/// (COMPGEN-REVERSAL-001), the unsigned breakpoint lower bound
/// <see cref="ComparativeGenomics.CalculateReversalDistance"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the O(n) breakpoint scan must always terminate), no state corruption, no
/// nonsense output (a negative distance, a value above the documented maximum, a
/// non-zero distance for the identity), and no *unhandled* runtime exception
/// (IndexOutOfRange on a singleton/empty permutation, off-by-one in the adjacency
/// walk, DivideByZero). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* validation exception
/// (ArgumentException when the two orders differ in length — contract §3.3, §6.1).
/// A raw runtime exception, a hang, a wrong distance, or distance ≠ 0 for the
/// identity is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-REVERSAL-001 — Reversal Distance (Breakpoint Lower Bound)
/// Checklist: docs/checklists/03_FUZZING.md, row 138.
/// Algorithm doc: docs/algorithms/Comparative_Genomics/Reversal_Distance.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row:
///          – IDENTITY PERMUTATION (perm1 == perm2, or any order vs itself → the
///            identity is the unique permutation with 0 breakpoints, so the
///            distance is 0 — the canonical no-op — §2.4 INV-01, §6.1).
///          – FULL REVERSAL (the completely reversed order [n..1] → [1..n]): the
///            extended unsigned permutation has exactly b = 2 breakpoints (the two
///            sentinel flanks; every internal pair is consecutive), so the
///            documented distance is ⌈2/2⌉ = 1, INDEPENDENT of n — §6.1 row
///            "fully reversed [4,3,2,1] → [1,2,3,4] ⇒ 1".
///          – SINGLETON (a length-1 permutation, and the empty permutation):
///            n ≤ 1 has no internal adjacency ⇒ 0, with NO IndexOutOfRange /
///            DivideByZero — §3.3, §6.1 "empty / single element ⇒ 0".
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope relative to COMPGEN-REARR-001 (row 137)
/// ───────────────────────────────────────────────────────────────────────────
/// Row 137 covers the breakpoint-model rearrangement *detector*
/// (DetectRearrangements / ClassifyRearrangement) — it COUNTS and classifies
/// disrupted signed adjacencies of two gene lists. THIS row is the reversal
/// *distance* METRIC over two integer permutations: the unsigned breakpoint lower
/// bound ⌈b/2⌉ on the minimum number of reversals (§2.5). The breakpoint count is
/// only a lower bound on the reversal distance (d ≥ b/2 — §2.2), so the two units
/// measure different quantities; this file is scoped to CalculateReversalDistance.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Reversal_Distance.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Inputs are UNSIGNED permutations of the same marker set. permutation2 is
/// relabelled to the identity 0..n−1 and permutation1 to its relative permutation;
/// the extended permutation is (−1, relative…, n). A consecutive pair is a
/// BREAKPOINT iff its two values are not consecutive integers (|Δ| ≠ 1), counting
/// the left flank (relative[0] ≠ 0) and right flank (relative[n−1] ≠ n−1) (§2.2,
/// §4.1). The result is ⌈b/2⌉ = (b + 1) / 2 (§2.2, §4.1 step 4). Invariants:
///   INV-01 d(π,π) = 0 (identity is the unique 0-breakpoint permutation);
///   INV-02 result ≥ 0;
///   INV-03 d(α,β) = d(β,α) (symmetric);
///   INV-04 result = ⌈b/2⌉ over the extended relative permutation;
///   INV-05 result is a LOWER BOUND on the true reversal distance (§2.4).
/// b lies in [0, n+1] (n+1 extended pairs), so the result lies in [0, ⌈(n+1)/2⌉].
/// Unequal-length inputs ⇒ ArgumentException; empty/single ⇒ 0 (§3.3, §6.1).
///   ComparativeGenomics.CalculateReversalDistance(
///       IReadOnlyList&lt;int&gt; permutation1, IReadOnlyList&lt;int&gt; permutation2) → int
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeReversalFuzzTests
{
    #region Helpers

    /// <summary>The completely reversed order of 1..n, i.e. [n, n-1, …, 1].</summary>
    private static int[] DescendingOneToN(int n)
        => Enumerable.Range(1, n).Reverse().ToArray();

    /// <summary>The sorted order 1..n.</summary>
    private static int[] AscendingOneToN(int n)
        => Enumerable.Range(1, n).ToArray();

    /// <summary>
    /// Asserts a distance is well-formed per the documented contract:
    /// non-negative (INV-02), within the documented ceiling ⌈(n+1)/2⌉ (b ∈ [0, n+1],
    /// result = ⌈b/2⌉ — INV-04), and equal to 0 iff the two orders are identical
    /// (the identity is the unique 0-breakpoint permutation — INV-01).
    /// </summary>
    private static void AssertWellFormed(int distance, int[] perm1, int[] perm2)
    {
        int n = perm1.Length;
        distance.Should().BeGreaterThanOrEqualTo(0, "INV-02: ⌈b/2⌉ with b ≥ 0 is non-negative");

        // The extended permutation of n markers has n+1 internal pairs ⇒ b ∈ [0, n+1],
        // and result = ⌈b/2⌉, so the maximum possible distance is ⌈(n+1)/2⌉.
        int max = n <= 1 ? 0 : (n + 1 + 1) / 2;
        distance.Should().BeLessThanOrEqualTo(max,
            $"INV-04: ⌈b/2⌉ cannot exceed ⌈(n+1)/2⌉ for n={n} markers");

        bool identical = perm1.SequenceEqual(perm2);
        if (identical)
            distance.Should().Be(0, "INV-01: identical orders have 0 breakpoints ⇒ distance 0");
        else if (n >= 2)
            distance.Should().BeGreaterThan(0,
                "a non-identity permutation of ≥2 markers has ≥1 breakpoint ⇒ distance ≥ 1");
    }

    #endregion

    #region COMPGEN-REVERSAL-001 — Reversal Distance (BE: identity permutation, full reversal, singleton)

    #region BE — Boundary: identity permutation (the canonical no-op, distance 0 — INV-01)

    // Identity at every size: an order compared to ITSELF has 0 breakpoints ⇒ distance 0.
    [Test]
    public void CalculateReversalDistance_IdentityOrderVsItself_AlwaysZero()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 13, 21 })
        {
            var order = AscendingOneToN(n);

            int d = ComparativeGenomics.CalculateReversalDistance(order, order);

            d.Should().Be(0, $"INV-01: an order vs itself has 0 breakpoints (n={n})");
            AssertWellFormed(d, order, order);
        }
    }

    // Identity holds for ARBITRARY (non-contiguous, negative, large) distinct marker ids: the
    // markers are relabelled to the relative permutation, so any order vs an equal copy ⇒ 0.
    [Test]
    public void CalculateReversalDistance_IdentityArbitraryMarkerIds_StillZero()
    {
        var orders = new[]
        {
            new[] { 42, 7, 1000, -5, 0 },
            new[] { int.MinValue, 0, int.MaxValue },
            new[] { -3, -2, -1 },
            new[] { 100 },
        };

        foreach (var order in orders)
        {
            // A fresh copy: equal VALUES but a distinct reference, so the no-op contract is exercised
            // through the relabelling path, not reference identity.
            var copy = order.ToArray();

            int d = ComparativeGenomics.CalculateReversalDistance(order, copy);

            d.Should().Be(0, "any order equals itself ⇒ 0 breakpoints ⇒ distance 0 (INV-01)");
            AssertWellFormed(d, order, copy);
        }
    }

    #endregion

    #region BE — Boundary: full reversal (b = 2 ⇒ distance 1, independent of n — §6.1)

    // The documented edge-case row: fully reversed [n..1] → [1..n] has exactly b = 2 breakpoints
    // (the two sentinel flanks; every internal pair (k+1,k) is consecutive |Δ|=1), so the distance
    // is ⌈2/2⌉ = 1 for EVERY n ≥ 2 — a single reversal sorts the whole block (§6.1, §2.2).
    [Test]
    public void CalculateReversalDistance_FullReversal_AlwaysOne()
    {
        foreach (int n in new[] { 2, 3, 4, 5, 8, 12, 50, 200 })
        {
            var reversed = DescendingOneToN(n);
            var sorted = AscendingOneToN(n);

            int forward = ComparativeGenomics.CalculateReversalDistance(reversed, sorted);
            int backward = ComparativeGenomics.CalculateReversalDistance(sorted, reversed);

            forward.Should().Be(1,
                $"a full reversal has b=2 ⇒ ⌈2/2⌉=1, independent of n (n={n}) — §6.1");
            backward.Should().Be(1, $"INV-03: reversal distance is symmetric (n={n})");
            AssertWellFormed(forward, reversed, sorted);
            AssertWellFormed(backward, sorted, reversed);
        }
    }

    // The exact documented edge-case literal: [4,3,2,1] → [1,2,3,4] ⇒ 1 (§6.1).
    [Test]
    public void CalculateReversalDistance_DocEdgeCaseLiteral_FullReversalIsOne()
    {
        int d = ComparativeGenomics.CalculateReversalDistance(
            new[] { 4, 3, 2, 1 }, new[] { 1, 2, 3, 4 });

        d.Should().Be(1, "doc §6.1: fully reversed [4,3,2,1] → [1,2,3,4] ⇒ ⌈2/2⌉ = 1");
    }

    #endregion

    #region BE — Boundary: singleton & empty (no internal adjacency ⇒ 0, no crash — §6.1)

    // A length-1 permutation: no internal adjacency to evaluate ⇒ 0, and critically NO
    // IndexOutOfRange / DivideByZero in the breakpoint walk (relative[n-1] with n=1).
    [Test]
    public void CalculateReversalDistance_Singleton_ReturnsZeroNoCrash()
    {
        foreach (var single in new[] { new[] { 0 }, new[] { 1 }, new[] { -99 }, new[] { int.MaxValue } })
        {
            Action act = () =>
            {
                int d = ComparativeGenomics.CalculateReversalDistance(single, single.ToArray());
                d.Should().Be(0, "a single marker has no internal adjacency ⇒ 0 (§6.1)");
            };
            act.Should().NotThrow("a length-1 permutation must not crash the breakpoint walk");
        }
    }

    // The empty permutation: n = 0 ⇒ 0, no exception (the n ≤ 1 short-circuit covers it).
    [Test]
    public void CalculateReversalDistance_Empty_ReturnsZeroNoCrash()
    {
        Action act = () =>
        {
            int d = ComparativeGenomics.CalculateReversalDistance(
                Array.Empty<int>(), Array.Empty<int>());
            d.Should().Be(0, "an empty permutation has no breakpoints ⇒ 0 (§6.1)");
        };

        act.Should().NotThrow("the empty permutation is a documented degenerate input, not an error");
    }

    #endregion

    #region BE — Boundary: unequal lengths (documented validation exception — §3.3, §6.1)

    // Different-length inputs ⇒ ArgumentException: the distance is undefined across different
    // marker sets (§3.3, §6.1, §5.4 deviation #2). Includes the empty-vs-nonempty boundary.
    [Test]
    public void CalculateReversalDistance_UnequalLengths_ThrowArgumentException()
    {
        var pairs = new (int[] a, int[] b)[]
        {
            (new[] { 1, 2, 3 }, new[] { 1, 2 }),
            (new[] { 1 }, Array.Empty<int>()),
            (Array.Empty<int>(), new[] { 1, 2 }),
            (new[] { 1, 2, 3, 4, 5 }, new[] { 1, 2, 3 }),
        };

        foreach (var (a, b) in pairs)
        {
            Action act = () => ComparativeGenomics.CalculateReversalDistance(a, b);
            act.Should().Throw<ArgumentException>(
                "unequal-length orders have an undefined distance (§3.3)");
        }
    }

    #endregion

    #region Positive sanity & metric correctness (the documented breakpoint contract)

    // The doc's worked example (§7.1): [2,3,1,6,5,4] → [1..6] has b=4 ⇒ ⌈4/2⌉ = 2.
    // A POSITIVE sanity check that the metric computes a real breakpoint count, not merely 0/1.
    [Test]
    public void CalculateReversalDistance_DocWorkedExample_ReturnsTwo()
    {
        int d = ComparativeGenomics.CalculateReversalDistance(
            new[] { 2, 3, 1, 6, 5, 4 }, new[] { 1, 2, 3, 4, 5, 6 });

        d.Should().Be(2,
            "doc §7.1: relative [1,2,0,5,4,3] has b=4 ⇒ ⌈4/2⌉ = 2");
    }

    // Hand-derived small case: [1,3,2] → [1,2,3]. Relative = [0,2,1]; extended (−1,0,2,1,3):
    // (−1,0)|Δ|=1, (0,2)|Δ|=2 ✔, (2,1)|Δ|=1, (1,3)|Δ|=2 ✔ ⇒ b=2 ⇒ ⌈2/2⌉ = 1 (one adjacent swap).
    [Test]
    public void CalculateReversalDistance_SingleAdjacentTransposition_ReturnsOne()
    {
        int d = ComparativeGenomics.CalculateReversalDistance(
            new[] { 1, 3, 2 }, new[] { 1, 2, 3 });

        d.Should().Be(1, "[1,3,2] has b=2 ⇒ ⌈2/2⌉ = 1");
    }

    // Hand-derived: [3,1,2] → [1,2,3]. Relative = [2,0,1]; extended (−1,2,0,1,3):
    // (−1,2)✔, (2,0)✔, (0,1)|Δ|=1, (1,3)✔ ⇒ b=3 ⇒ ⌈3/2⌉ = 2. The odd-breakpoint rounding boundary.
    [Test]
    public void CalculateReversalDistance_OddBreakpointCount_RoundsUp()
    {
        int d = ComparativeGenomics.CalculateReversalDistance(
            new[] { 3, 1, 2 }, new[] { 1, 2, 3 });

        d.Should().Be(2, "[3,1,2] has b=3 ⇒ ⌈3/2⌉ = 2 (rounds up)");
    }

    // Symmetry (INV-03): d(α,β) = d(β,α) for an arbitrary hand-picked pair.
    [Test]
    public void CalculateReversalDistance_IsSymmetric()
    {
        var a = new[] { 2, 3, 1, 6, 5, 4 };
        var b = new[] { 1, 2, 3, 4, 5, 6 };

        int forward = ComparativeGenomics.CalculateReversalDistance(a, b);
        int backward = ComparativeGenomics.CalculateReversalDistance(b, a);

        backward.Should().Be(forward, "INV-03: reversal distance is symmetric");
    }

    // Metric discrimination: on the SAME n, identity ⇒ 0, full reversal ⇒ 1, and a maximally
    // scrambled order ⇒ a strictly larger distance — so the three outcomes are visibly distinct.
    [Test]
    public void CalculateReversalDistance_DiscriminatesIdentityReversalAndScramble()
    {
        const int n = 6;
        var sorted = AscendingOneToN(n);

        int identity = ComparativeGenomics.CalculateReversalDistance(sorted, sorted);
        int fullReversal = ComparativeGenomics.CalculateReversalDistance(DescendingOneToN(n), sorted);
        // An alternating "every adjacency broken" order: [2,4,6,...,1,3,5,...] maximises breakpoints.
        var scrambled = new[] { 2, 4, 6, 1, 3, 5 };
        int scramble = ComparativeGenomics.CalculateReversalDistance(scrambled, sorted);

        identity.Should().Be(0, "identity ⇒ 0 (INV-01)");
        fullReversal.Should().Be(1, "full reversal ⇒ ⌈2/2⌉ = 1 (§6.1)");
        scramble.Should().BeGreaterThan(fullReversal,
            "a maximally broken order has more breakpoints than a single full reversal");
    }

    #endregion

    #region Robustness — randomized permutations always well-formed (no hang, no crash)

    // Random permutations of varying size: the result must ALWAYS be well-formed (non-negative,
    // ≤ ⌈(n+1)/2⌉, 0 iff identity), regardless of the random order, with no hang (O(n) scan).
    [Test]
    [CancelAfter(30000)]
    public void CalculateReversalDistance_RandomPermutations_AlwaysWellFormed()
    {
        var rng = new Random(138); // locally-seeded, deterministic.

        for (int trial = 0; trial < 500; trial++)
        {
            int n = rng.Next(0, 24);
            var target = AscendingOneToN(n);

            // Random permutation of the marker set via a Fisher–Yates shuffle.
            var perm = target.ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (perm[i], perm[j]) = (perm[j], perm[i]);
            }

            int d = ComparativeGenomics.CalculateReversalDistance(perm, target);

            AssertWellFormed(d, perm, target);
            // INV-03 cross-check on every trial: symmetric distance.
            ComparativeGenomics.CalculateReversalDistance(target, perm)
                .Should().Be(d, "INV-03: distance is symmetric for any random pair");
        }
    }

    #endregion

    #endregion
}
