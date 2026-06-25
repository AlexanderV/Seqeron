using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology clonal/subclonal classification area — ONCO-CLONAL-001.
/// The units under test are the deterministic clonality entry points implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
/// <see cref="OncologyAnalyzer.IdentifyClonalMutations(IEnumerable{double})"/>
/// — the point-estimate classifier that selects the indices of CCF values exceeding the
/// clonal threshold (CCF &gt; 0.95, strict) — and
/// <see cref="OncologyAnalyzer.ClassifyClonality(IEnumerable{OncologyAnalyzer.ClonalityVariant}, double)"/>
/// — the posterior-grid classifier producing per-variant calls, counts and the clonal fraction.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault. Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome —
/// here, an ArgumentNullException for null inputs, or an ArgumentException for a
/// CCF value that is NaN or outside [0, 1]. There is NO silent clamping in
/// IdentifyClonalMutations: a CCF &gt; 1 (overshoot from sampling noise / copy-number
/// error) is REJECTED, not coerced, because the point-estimate contract requires
/// each CCF ∈ [0, 1] (Clonal_Subclonal_Classification.md §3.1, §3.3). A CCF = 0
/// (no cancer cell carries the mutation) is a valid in-range value that is simply
/// not clonal — no divide-by-zero. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CLONAL-001 — clonal vs subclonal mutation classification (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 108.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 108): "CCF at threshold, CCF&gt;1, CCF=0".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The targets map onto the documented contract as follows:
///   • CCF at threshold  ⇔ CCF EXACTLY 0.95 ⇒ NOT clonal — the rule is strict
///                          (CCF &gt; 0.95), so the cutoff itself is subclonal
///                          (Clonal_Subclonal_Classification.md §4.2, §6.1, INV-06);
///   • CCF &gt; 1          ⇔ a CCF exceeding 1.0 (noise/CN overshoot) ⇒ REJECTED
///                          with ArgumentException; the point-estimate contract is
///                          CCF ∈ [0, 1], not silently clamped (§3.1, §3.3);
///   • CCF = 0           ⇔ no cancer cells carry the mutation ⇒ valid, NOT clonal,
///                          no divide-by-zero (§3.1 ccfValues each ∈ [0, 1]).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Clonal_Subclonal_Classification.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • IdentifyClonalMutations selects index i ⇔ ccf[i] &gt; 0.95 (STRICT)   (INV-06, §4.2, §6.1)
///   • ClonalCcfThreshold == 0.95                                          (§4.2)
///   • ClassifyClonality: clonal ⇔ P(CCF &gt; 0.95) &gt; 0.5                  (INV-05, §2.2)
///   • ClonalCount + SubclonalCount == total                              (INV-01)
///   • ClonalFraction == ClonalCount / total (0 if empty)                 (INV-02)
///   • CCF point estimate ∈ [0.01, 1]                                     (INV-03)
///   • ProbabilityClonal ∈ [0, 1]                                         (INV-04)
///   • Empty input ⇒ empty calls, counts 0, ClonalFraction 0             (§3.3, §6.1)
///   • null variants / null ccfValues ⇒ ArgumentNullException            (§3.3)
///   • purity NaN or ∉ (0, 1] ⇒ ArgumentOutOfRangeException              (§3.3)
///   • a CCF NaN or ∉ [0, 1] ⇒ ArgumentException                         (§3.3)
///   • an invalid variant (N&lt;1, a∉[0,N], q&lt;1, M∉[1,q]) ⇒ ArgumentException (§3.3)
///
/// No source bug was found; no test was weakened.
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyClonalityFuzzTests
{
    private const double Threshold = ClonalCcfThreshold; // 0.95 (Landau et al. 2013)

    // ── Well-formed-result assertion helpers ─────────────────────────────────
    // Pin the documented numeric contract on EVERY returned result so a fuzz
    // test cannot rubber-stamp a malformed result (mis-summed counts, out-of-range
    // CCF, NaN probability) green.
    private static void AssertWellFormedResult(ClonalityResult result, int expectedTotal)
    {
        result.Calls.Should().NotBeNull("a result always carries its per-variant call list");
        result.Calls.Count.Should().Be(expectedTotal, "calls are returned one per input variant, in order");
        result.ClonalCount.Should().BeGreaterThanOrEqualTo(0, "a count is non-negative");
        result.SubclonalCount.Should().BeGreaterThanOrEqualTo(0, "a count is non-negative");
        (result.ClonalCount + result.SubclonalCount).Should().Be(expectedTotal,
            "every variant is classified into exactly one of two states (INV-01)");

        double expectedFraction = expectedTotal == 0 ? 0.0 : (double)result.ClonalCount / expectedTotal;
        result.ClonalFraction.Should().BeApproximately(expectedFraction, 1e-12,
            "ClonalFraction == ClonalCount / total, 0 if empty (INV-02)");

        int clonalSeen = 0;
        foreach (ClonalityCall call in result.Calls)
        {
            call.Ccf.Should().BeInRange(CcfGridLowerBoundForTest, 1.0,
                "the CCF point estimate is the posterior mean over the grid c ∈ [0.01, 1] (INV-03)");
            double.IsNaN(call.Ccf).Should().BeFalse("a posterior-mean CCF is finite (INV-03)");
            call.ProbabilityClonal.Should().BeInRange(0.0, 1.0,
                "ProbabilityClonal is normalised posterior mass (INV-04)");
            double.IsNaN(call.ProbabilityClonal).Should().BeFalse("a normalised probability is finite (INV-04)");

            bool isClonal = call.Status == ClonalityStatus.Clonal;
            isClonal.Should().Be(call.ProbabilityClonal > 0.5,
                "clonal ⇔ P(CCF > 0.95) > 0.5 (INV-05)");
            if (isClonal)
            {
                clonalSeen++;
            }
        }

        clonalSeen.Should().Be(result.ClonalCount, "ClonalCount equals the number of clonal calls (INV-01)");
    }

    // The grid lower bound (0.01) is private in the source; the documented INV-03 range
    // is [0.01, 1]. We replicate it here for the well-formed assertion only.
    private const double CcfGridLowerBoundForTest = 0.01;

    #region ONCO-CLONAL-001 — positive sanity (high CCF clonal, low CCF subclonal, threshold strict)

    // IdentifyClonalMutations: high CCF (~1.0) → clonal, low CCF (~0.3) → subclonal,
    // and the threshold value 0.95 itself → NOT clonal (strict >).
    [Test]
    public void IdentifyClonalMutations_PositiveSanity_SelectsOnlyAboveThreshold()
    {
        double[] ccf = { 0.30, 0.95, 0.951, 0.99, 1.0, 0.0, 0.5 };

        IReadOnlyList<int> clonal = IdentifyClonalMutations(ccf);

        // Indices 3 (0.99), 4 (1.0) and 2 (0.951) exceed 0.95; 0.95 itself (index 1) does not.
        clonal.Should().Equal(2, 3, 4);
    }

    // ClassifyClonality: a deep, high-VAF clonal variant vs a low-VAF subclonal variant.
    [Test]
    public void ClassifyClonality_PositiveSanity_HighVafClonal_LowVafSubclonal()
    {
        // From Clonal_Subclonal_Classification.md §7.1 worked example.
        var variants = new[]
        {
            new ClonalityVariant(altReads: 400, totalReads: 1000, localCopyNumber: 2), // clonal
            new ClonalityVariant(altReads: 240, totalReads: 1000, localCopyNumber: 2), // subclonal
        };

        ClonalityResult result = ClassifyClonality(variants, purity: 0.8);

        AssertWellFormedResult(result, expectedTotal: 2);
        result.Calls[0].Status.Should().Be(ClonalityStatus.Clonal, "VAF 0.40 at purity 0.8, q=2 ⇒ CCF ≈ 1.0");
        result.Calls[1].Status.Should().Be(ClonalityStatus.Subclonal, "VAF 0.24 ⇒ CCF ≈ 0.6");
        result.ClonalCount.Should().Be(1);
        result.SubclonalCount.Should().Be(1);
        result.ClonalFraction.Should().BeApproximately(0.5, 1e-12);
    }

    #endregion

    #region ONCO-CLONAL-001 — BE: CCF at threshold (CCF == 0.95 is NOT clonal, strict >)

    [Test]
    public void IdentifyClonalMutations_CcfExactlyThreshold_IsNotClonal()
    {
        IReadOnlyList<int> clonal = IdentifyClonalMutations(new[] { Threshold });

        clonal.Should().BeEmpty("the clonal rule is strict (CCF > 0.95); the cutoff itself is subclonal (INV-06, §6.1)");
    }

    [Test]
    public void IdentifyClonalMutations_JustAboveThreshold_IsClonal()
    {
        // The smallest representable value above 0.95 must flip to clonal — no off-by-one slack.
        double justAbove = Math.BitIncrement(Threshold);

        IReadOnlyList<int> clonal = IdentifyClonalMutations(new[] { justAbove });

        clonal.Should().Equal(new[] { 0 }, "any value strictly greater than 0.95 is clonal (INV-06)");
    }

    [Test]
    public void IdentifyClonalMutations_JustBelowThreshold_IsNotClonal()
    {
        double justBelow = Math.BitDecrement(Threshold);

        IReadOnlyList<int> clonal = IdentifyClonalMutations(new[] { justBelow });

        clonal.Should().BeEmpty("a value below the cutoff is subclonal (INV-06)");
    }

    // Fuzz a band around the cutoff: each value is clonal IFF it is strictly > 0.95, never otherwise.
    [Test]
    public void IdentifyClonalMutations_FuzzAroundThreshold_StrictRuleHolds()
    {
        var rng = new Random(108_001);
        for (int trial = 0; trial < 2000; trial++)
        {
            // Values within ±0.05 of the threshold, clamped to [0, 1].
            double ccf = Math.Clamp(Threshold + (rng.NextDouble() - 0.5) * 0.10, 0.0, 1.0);

            IReadOnlyList<int> clonal = IdentifyClonalMutations(new[] { ccf });

            bool selected = clonal.Count == 1;
            selected.Should().Be(ccf > Threshold,
                $"clonal ⇔ CCF > 0.95 (strict); ccf={ccf:R} (trial {trial})");
        }
    }

    #endregion

    #region ONCO-CLONAL-001 — BE: CCF > 1 (overshoot is REJECTED, not silently clamped)

    [Test]
    public void IdentifyClonalMutations_CcfJustAboveOne_Throws()
    {
        double overshoot = Math.BitIncrement(1.0);

        Action act = () => IdentifyClonalMutations(new[] { overshoot });

        act.Should().Throw<ArgumentException>("a CCF must be in [0, 1]; > 1 is not silently clamped (§3.1, §3.3)")
            .WithMessage("*[0, 1]*");
    }

    [Test]
    public void IdentifyClonalMutations_CcfFarAboveOne_Throws()
    {
        // CNAqc reports e.g. 1.06 under sampling noise; the point-estimate contract still rejects it.
        Action act = () => IdentifyClonalMutations(new[] { 0.5, 1.06 });

        act.Should().Throw<ArgumentException>("an out-of-range CCF is a contract violation (§3.3)");
    }

    [Test]
    public void IdentifyClonalMutations_FuzzAboveOne_AlwaysThrows()
    {
        var rng = new Random(108_002);
        for (int trial = 0; trial < 1000; trial++)
        {
            double overshoot = 1.0 + rng.NextDouble() * 1e6 + double.Epsilon;

            Action act = () => IdentifyClonalMutations(new[] { overshoot });

            act.Should().Throw<ArgumentException>($"CCF > 1 is rejected; ccf={overshoot:R} (trial {trial})");
        }
    }

    [Test]
    public void IdentifyClonalMutations_PositiveInfinityAndNaN_Throw()
    {
        Action infinity = () => IdentifyClonalMutations(new[] { double.PositiveInfinity });
        Action nan = () => IdentifyClonalMutations(new[] { double.NaN });

        infinity.Should().Throw<ArgumentException>("+Inf is not in [0, 1]");
        nan.Should().Throw<ArgumentException>("NaN is explicitly rejected (§3.3)");
    }

    #endregion

    #region ONCO-CLONAL-001 — BE: CCF = 0 (valid, not clonal, no divide-by-zero)

    [Test]
    public void IdentifyClonalMutations_CcfZero_IsValidAndNotClonal()
    {
        IReadOnlyList<int> clonal = IdentifyClonalMutations(new[] { 0.0 });

        clonal.Should().BeEmpty("CCF = 0 (no cancer cell carries the mutation) is in-range and not clonal (§3.1)");
    }

    [Test]
    public void IdentifyClonalMutations_NegativeCcf_Throws()
    {
        Action act = () => IdentifyClonalMutations(new[] { -Math.BitIncrement(0.0) });

        act.Should().Throw<ArgumentException>("a CCF below 0 is out of range (§3.3)");
    }

    [Test]
    public void IdentifyClonalMutations_AllZeros_ReturnsNoneAndDoesNotThrow()
    {
        double[] ccf = Enumerable.Repeat(0.0, 50).ToArray();

        IReadOnlyList<int> clonal = IdentifyClonalMutations(ccf);

        clonal.Should().BeEmpty("no zero-CCF mutation is clonal; iterating zeros must not divide by zero");
    }

    // ClassifyClonality with zero alternate reads ⇒ posterior mass near c=0.01 ⇒ subclonal, finite.
    [Test]
    public void ClassifyClonality_ZeroAltReads_IsSubclonalAndFinite()
    {
        var variants = new[] { new ClonalityVariant(altReads: 0, totalReads: 500, localCopyNumber: 2) };

        ClonalityResult result = ClassifyClonality(variants, purity: 0.8);

        AssertWellFormedResult(result, expectedTotal: 1);
        result.Calls[0].Status.Should().Be(ClonalityStatus.Subclonal, "no alt reads ⇒ CCF ≈ 0 ⇒ subclonal");
        result.ClonalCount.Should().Be(0);
    }

    #endregion

    #region ONCO-CLONAL-001 — BE: empty input & documented exceptions

    [Test]
    public void IdentifyClonalMutations_EmptyInput_ReturnsEmpty()
    {
        IdentifyClonalMutations(Array.Empty<double>()).Should().BeEmpty("nothing to classify (§6.1)");
    }

    [Test]
    public void ClassifyClonality_EmptyInput_ReturnsZeroCountsAndZeroFraction()
    {
        ClonalityResult result = ClassifyClonality(Array.Empty<ClonalityVariant>(), purity: 0.5);

        AssertWellFormedResult(result, expectedTotal: 0);
        result.ClonalCount.Should().Be(0);
        result.SubclonalCount.Should().Be(0);
        result.ClonalFraction.Should().Be(0.0, "clonal fraction is 0 for an empty set (INV-02, §6.1)");
    }

    [Test]
    public void IdentifyClonalMutations_NullInput_Throws()
    {
        Action act = () => IdentifyClonalMutations(null!);
        act.Should().Throw<ArgumentNullException>("null ccfValues is rejected (§3.3)");
    }

    [Test]
    public void ClassifyClonality_NullVariants_Throws()
    {
        Action act = () => ClassifyClonality(null!, purity: 0.5);
        act.Should().Throw<ArgumentNullException>("null variants is rejected (§3.3)");
    }

    [TestCase(0.0)]
    [TestCase(-0.1)]
    [TestCase(1.0001)]
    [TestCase(double.NaN)]
    public void ClassifyClonality_PurityOutOfRange_Throws(double purity)
    {
        Action act = () => ClassifyClonality(Array.Empty<ClonalityVariant>(), purity);
        act.Should().Throw<ArgumentOutOfRangeException>("purity must be in (0, 1] (§3.3)");
    }

    [Test]
    public void ClassifyClonality_PurityExactlyOne_IsAccepted()
    {
        // Boundary: purity = 1 is the inclusive upper end of (0, 1].
        var variants = new[] { new ClonalityVariant(altReads: 480, totalReads: 1000, localCopyNumber: 2) };

        Action act = () => ClassifyClonality(variants, purity: 1.0);

        act.Should().NotThrow("purity = 1 is the inclusive upper bound of (0, 1] (§3.1)");
    }

    [TestCase(0, 0, 2, 1)]   // TotalReads < 1
    [TestCase(5, 4, 2, 1)]   // AltReads > TotalReads
    [TestCase(-1, 10, 2, 1)] // AltReads < 0
    [TestCase(5, 10, 0, 1)]  // LocalCopyNumber < 1
    [TestCase(5, 10, 2, 3)]  // Multiplicity > LocalCopyNumber
    [TestCase(5, 10, 2, 0)]  // Multiplicity < 1
    public void ClassifyClonality_InvalidVariant_Throws(int alt, int total, int cn, int mult)
    {
        var variants = new[] { new ClonalityVariant(alt, total, cn, mult) };

        Action act = () => ClassifyClonality(variants, purity: 0.7);

        act.Should().Throw<ArgumentException>("an invalid variant is rejected (§3.3)");
    }

    #endregion

    #region ONCO-CLONAL-001 — robustness: large / random sweeps (no hang, contract holds)

    [Test]
    [CancelAfter(20_000)]
    public void IdentifyClonalMutations_RandomInRangeSweep_MatchesStrictRule()
    {
        var rng = new Random(108_010);
        for (int trial = 0; trial < 500; trial++)
        {
            int n = rng.Next(0, 64);
            double[] ccf = new double[n];
            for (int i = 0; i < n; i++)
            {
                ccf[i] = rng.NextDouble(); // in [0, 1)
            }

            IReadOnlyList<int> clonal = IdentifyClonalMutations(ccf);

            int[] expected = Enumerable.Range(0, n).Where(i => ccf[i] > Threshold).ToArray();
            clonal.Should().Equal(expected, "selection is exactly the strict CCF > 0.95 rule, in order (INV-06)");
        }
    }

    [Test]
    [CancelAfter(20_000)]
    public void ClassifyClonality_RandomVariantsSweep_ProducesWellFormedResults()
    {
        var rng = new Random(108_011);
        for (int trial = 0; trial < 300; trial++)
        {
            int n = rng.Next(0, 16);
            var variants = new ClonalityVariant[n];
            for (int i = 0; i < n; i++)
            {
                int total = rng.Next(1, 5000);
                int alt = rng.Next(0, total + 1);
                int cn = rng.Next(1, 8);
                int mult = rng.Next(1, cn + 1);
                variants[i] = new ClonalityVariant(alt, total, cn, mult);
            }

            double purity = Math.Clamp(rng.NextDouble(), 1e-6, 1.0);

            ClonalityResult result = ClassifyClonality(variants, purity);

            AssertWellFormedResult(result, expectedTotal: n);
        }
    }

    [Test]
    [CancelAfter(20_000)]
    public void ClassifyClonality_ExtremeReadDepth_StaysFinite()
    {
        // Very deep coverage stresses the log-space binomial kernel for underflow.
        var variants = new[]
        {
            new ClonalityVariant(altReads: int.MaxValue / 2, totalReads: int.MaxValue, localCopyNumber: 2),
            new ClonalityVariant(altReads: 0, totalReads: int.MaxValue, localCopyNumber: 2),
            new ClonalityVariant(altReads: int.MaxValue, totalReads: int.MaxValue, localCopyNumber: 2),
        };

        ClonalityResult result = ClassifyClonality(variants, purity: 0.9);

        AssertWellFormedResult(result, expectedTotal: 3);
    }

    #endregion
}
