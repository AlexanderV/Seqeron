// COMPGEN-REVERSAL-001 — Reversal Distance (unsigned breakpoint lower bound)
// Evidence: docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-REVERSAL-001.md
// Source: Bafna V, Pevzner PA (1998). Sorting by Transpositions. SIAM J. Discrete Math. 11(2):224-240, §2.
//         Hunter College CompBio Lecture 16 — sorting by reversals (extended permutation, d >= b/2).
//         Hübotter J (2020). On Sorting by Reversals (unsigned breakpoint |Δ| != 1).

using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using System;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ComparativeGenomics_CalculateReversalDistance_Tests
{
    #region CalculateReversalDistance

    // M1 — Identity: identity is the only permutation with 0 breakpoints (Bafna & Pevzner 1998, §2),
    //      so the reversal distance to itself is 0.
    [Test]
    public void CalculateReversalDistance_IdenticalOrders_ReturnsZero()
    {
        var order = new[] { 1, 2, 3, 4, 5 };

        int d = ComparativeGenomics.CalculateReversalDistance(order, order);

        Assert.That(d, Is.EqualTo(0),
            "Identity permutation has 0 breakpoints (Bafna & Pevzner 1998 §2) ⇒ distance 0.");
    }

    // M2 — Hunter Lecture 16 worked example (unsigned specialization). Signed α=(−2,−3,+1,+6,−5,−4)
    //      has b=6; the unsigned magnitudes [2,3,1,6,5,4] vs identity give the extended permutation
    //      (−1,1,2,0,5,4,3,6) with breakpoints at (−1,1),(2,0),(0,5),(3,6) ⇒ b=4 ⇒ ⌈4/2⌉ = 2.
    [Test]
    public void CalculateReversalDistance_HunterWorkedExampleUnsigned_ReturnsTwo()
    {
        var perm1 = new[] { 2, 3, 1, 6, 5, 4 };
        var perm2 = new[] { 1, 2, 3, 4, 5, 6 };

        int d = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(d, Is.EqualTo(2),
            "Unsigned breakpoint count b=4 (Hunter Lecture 16 example, unsigned) ⇒ ⌈b/2⌉ = 2.");
    }

    // M3 — Fully reversed sequence. Extended [0,4,3,2,1,5]: breakpoints 0→4 and 1→5 ⇒ b=2 ⇒ ⌈2/2⌉ = 1.
    //      A single reversal of the whole block sorts it, so the bound is exact here.
    [Test]
    public void CalculateReversalDistance_FullyReversed_ReturnsOne()
    {
        var perm1 = new[] { 4, 3, 2, 1 };
        var perm2 = new[] { 1, 2, 3, 4 };

        int d = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(d, Is.EqualTo(1),
            "Fully reversed order has b=2 breakpoints (only the two end boundaries) ⇒ ⌈2/2⌉ = 1.");
    }

    // M4 — Single adjacent swap. relative=[0,1,3,2]; extended (−1,0,1,3,2,4): breakpoints (1,3),(2,4)
    //      ⇒ b=2 ⇒ ⌈2/2⌉ = 1 (TestSpec §7.2 derivation).
    [Test]
    public void CalculateReversalDistance_SingleAdjacentSwap_ReturnsOne()
    {
        var perm1 = new[] { 1, 2, 4, 3 };
        var perm2 = new[] { 1, 2, 3, 4 };

        int d = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(d, Is.EqualTo(1),
            "Swapping the last two markers yields b=2 breakpoints ⇒ ⌈2/2⌉ = 1.");
    }

    // M5 — Lower-bound property: the result must never exceed the number of reversals actually used
    //      to construct perm1 (Hunter Lecture 16: b(α) ≤ 2t ⇒ d(α) ≥ b/2, i.e. bound ≤ true distance).
    //      Build [1..8] by applying 3 reversals, then assert the bound is between 0 and 3.
    [Test]
    public void CalculateReversalDistance_IsLowerBoundOnAppliedReversals()
    {
        var target = new[] { 1, 2, 3, 4, 5, 6, 7, 8 };
        var perm = (int[])target.Clone();
        // Reversal 1: indices [1..4]; Reversal 2: [0..2]; Reversal 3: [4..7].
        Reverse(perm, 1, 4);
        Reverse(perm, 0, 2);
        Reverse(perm, 4, 7);
        const int appliedReversals = 3;

        int d = ComparativeGenomics.CalculateReversalDistance(perm, target);

        Assert.Multiple(() =>
        {
            Assert.That(d, Is.GreaterThanOrEqualTo(0),
                "A reversal-count lower bound is non-negative.");
            Assert.That(d, Is.LessThanOrEqualTo(appliedReversals),
                "The breakpoint bound ⌈b/2⌉ never exceeds the number of reversals actually applied.");
        });
    }

    // S1 — Empty inputs: n=0 ⇒ no internal adjacency ⇒ 0.
    [Test]
    public void CalculateReversalDistance_EmptyInputs_ReturnsZero()
    {
        int d = ComparativeGenomics.CalculateReversalDistance(Array.Empty<int>(), Array.Empty<int>());

        Assert.That(d, Is.EqualTo(0), "An empty permutation has no breakpoints ⇒ distance 0.");
    }

    // S2 — Single element: n=1 ⇒ no internal adjacency ⇒ 0.
    [Test]
    public void CalculateReversalDistance_SingleElement_ReturnsZero()
    {
        var single = new[] { 7 };

        int d = ComparativeGenomics.CalculateReversalDistance(single, single);

        Assert.That(d, Is.EqualTo(0), "A single-marker permutation has no breakpoints ⇒ distance 0.");
    }

    // S3 — Unequal lengths: distance is undefined across different marker sets ⇒ ArgumentException.
    [Test]
    public void CalculateReversalDistance_UnequalLengths_Throws()
    {
        var perm1 = new[] { 1, 2 };
        var perm2 = new[] { 1 };

        Assert.Throws<ArgumentException>(
            () => ComparativeGenomics.CalculateReversalDistance(perm1, perm2),
            "Reversal distance is undefined when the two orders cover different marker sets.");
    }

    // S4 — Symmetry: d(α, β) = d(β, α) (Hunter Lecture 16: d_β(α) = d_α(β)).
    [Test]
    public void CalculateReversalDistance_IsSymmetric()
    {
        var a = new[] { 2, 3, 1, 6, 5, 4 };
        var b = new[] { 1, 2, 3, 4, 5, 6 };

        int dab = ComparativeGenomics.CalculateReversalDistance(a, b);
        int dba = ComparativeGenomics.CalculateReversalDistance(b, a);

        Assert.That(dab, Is.EqualTo(dba),
            "Reversal distance is symmetric: d(α,β) = d(β,α) (Hunter Lecture 16).");
    }

    // C1 — Arbitrary (non 1..n) labels: relabelling to the relative permutation must give the same
    //      breakpoint structure. [30,10,20] vs [10,20,30] ⇒ relative [2,0,1]; extended (−1,2,0,1,3):
    //      breakpoints (−1,2),(2,0),(1,3) ⇒ b=3 ⇒ ⌈3/2⌉ = 2.
    [Test]
    public void CalculateReversalDistance_ArbitraryLabels_UsesRelativeOrder()
    {
        var perm1 = new[] { 30, 10, 20 };
        var perm2 = new[] { 10, 20, 30 };

        int d = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(d, Is.EqualTo(2),
            "Cyclic shift [30,10,20] over target [10,20,30] has b=3 breakpoints ⇒ ⌈3/2⌉ = 2.");
    }

    private static void Reverse(int[] a, int i, int j)
    {
        while (i < j)
        {
            (a[i], a[j]) = (a[j], a[i]);
            i++;
            j--;
        }
    }

    #endregion
}
