using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Epigenetics area — epigenetic-age estimation (EPIGEN-AGE-001).
/// The unit under test is the Horvath-style DNA-methylation clock entry point in
/// <see cref="EpigeneticsAnalyzer"/>:
///   • <see cref="EpigeneticsAnalyzer.CalculateEpigeneticAge"/> — the linear predictor
///     <c>Y = intercept + Σ coef_i · β_i</c> over the clock CpGs, mapped to years by
///     the Horvath inverse calibration F⁻¹.
///   • <see cref="EpigeneticsAnalyzer.HorvathAntiTransform"/> — the published two-branch
///     <c>anti.trafo</c> (helper exercised indirectly and directly for the transform).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the code
/// NEVER fails in an undisciplined way: no hang, no state corruption, no nonsense
/// output, and no *unhandled* runtime exception (NullReference / DivideByZero /
/// Overflow / NaN). Every input must resolve to EITHER a well-defined, theory-correct
/// value OR a *documented, intentional* outcome. For the epigenetic clock the headline
/// hazards (this row's BE targets) are:
///   • NO CLOCK SITES — none of the input CpGs match the coefficient table: the
///     weighted sum collapses to the intercept; the result must be exactly F⁻¹(intercept),
///     never NaN and never a DivideByZero from a phantom "site-count" normalisation.
///   • ALL-METHYLATED — every β = 1.0: the predictor is intercept + Σ coef_i; the result
///     must be finite (no overflow / Infinity / NaN in the exp/linear transform).
///   • EMPTY — an empty methylation map: documented as valid, result = F⁻¹(intercept);
///     a null map / null or empty coefficient table is a documented throw.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-AGE-001 — Epigenetic age estimation (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 181.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///     Targets (checklist row 181): "no clock sites, all-methylated, empty".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The epigenetic-clock contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Linear predictor (Horvath 2013 reference R code, StepwiseAnalysis.R):
///     Y = intercept + Σ_i (coef_i · β_i)   over CpGs present in BOTH maps.
///   CpGs absent from the coefficient table do NOT enter the sum (INV-05).
///   — docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md §2.2, §4.1.
///
/// Inverse calibration F⁻¹ = anti.trafo, adult.age = 20 (horvath2013.R):
///     F⁻¹(Y) = 21·exp(Y) − 1   if Y < 0   (INV-02)
///     F⁻¹(Y) = 21·Y + 20       if Y ≥ 0   (INV-01)
///     F⁻¹(0) = 20.0 exactly                (INV-03)
///   F⁻¹ is strictly increasing and continuous; the branches meet at (0, 20) (INV-04).
///   — docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md §2.2, §2.4, §6.1.
///
/// Preconditions / documented guards (§3.3, §6.1):
///   • methylationAtClockCpGs == null      → ArgumentNullException.
///   • coefficients == null                → ArgumentNullException.
///   • coefficients.Count == 0             → ArgumentException (empty clock undefined).
///   • EMPTY methylation map               → VALID; result = F⁻¹(intercept).
///   • CpG not in the coefficient table    → ignored (only clock CpGs contribute).
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticAgeFuzzTests
{
    // ── Reference reimplementation of the documented contract ────────────────
    // Independent recomputation of F⁻¹ (anti.trafo) used to cross-check the unit
    // without trusting its own arithmetic. This is the published two-branch form.
    private static double ExpectedAntiTransform(double y)
        => y < 0 ? 21.0 * Math.Exp(y) - 1.0 : 21.0 * y + 20.0;

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the discipline contract: a predicted age must be a real, finite number
    // (never NaN / Infinity) and, for any finite linear predictor, must lie in the
    // mathematically reachable range of anti.trafo: [ -1, +∞ ) — the exponential
    // branch 21·exp(Y)−1 asymptotes to −1 from above as Y→−∞ (and, once exp(Y)
    // underflows to 0.0 in IEEE-754 double, reaches exactly −1.0), while the linear
    // branch is unbounded above. So a valid age is ALWAYS ≥ −1 and never below it.
    // This is what stops a test from rubber-stamping a NaN/Infinity or a sub-(−1)
    // artefact green.
    private static void AssertWellFormedAge(double age)
    {
        double.IsNaN(age).Should().BeFalse("predicted age must never be NaN");
        double.IsInfinity(age).Should().BeFalse("predicted age must never be Infinity");
        age.Should().BeGreaterThanOrEqualTo(-1.0,
            because: "anti.trafo asymptotes to −1 from above; age can never drop below it");
    }

    #region EPIGEN-AGE-001 — Epigenetic age (Horvath clock)

    // ── BE: EMPTY methylation map (documented valid → F⁻¹(intercept)) ─────────

    [Test]
    public void Empty_MethylationMap_PositiveIntercept_ReturnsAntiTransformOfIntercept()
    {
        // Docs §3.3 / §6.1: an empty methylation map is VALID — no CpG contributions,
        // so Y = intercept and age = F⁻¹(intercept). With intercept = 1.0 ≥ 0 the
        // linear branch gives 21·1 + 20 = 41 exactly.
        var coefficients = new Dictionary<string, double> { ["cg00000029"] = 0.5 };
        var methylation = new Dictionary<string, double>(); // empty

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 1.0);

        age.Should().BeApproximately(41.0, 1e-9);
        age.Should().BeApproximately(ExpectedAntiTransform(1.0), 1e-12);
        AssertWellFormedAge(age);
    }

    [Test]
    public void Empty_MethylationMap_ZeroIntercept_ReturnsExactly20()
    {
        // INV-03: Y = intercept = 0 → linear branch → 21·0 + 20 = 20.0 exactly.
        var coefficients = new Dictionary<string, double> { ["cg1"] = -3.7, ["cg2"] = 9.1 };
        var methylation = new Dictionary<string, double>();

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        age.Should().Be(20.0);
        AssertWellFormedAge(age);
    }

    [Test]
    public void Empty_MethylationMap_NegativeIntercept_TakesExponentialBranch()
    {
        // INV-02: Y = intercept = −0.5 < 0 → 21·exp(−0.5) − 1 ≈ 11.737, finite, > −1.
        var coefficients = new Dictionary<string, double> { ["cg1"] = 1.0 };
        var methylation = new Dictionary<string, double>();

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: -0.5);

        age.Should().BeApproximately(ExpectedAntiTransform(-0.5), 1e-12);
        age.Should().BeLessThan(20.0, "an exponential-branch (Y<0) result is below adult age");
        AssertWellFormedAge(age);
    }

    // ── BE: NO CLOCK SITES (none of the input CpGs match → intercept only) ────

    [Test]
    public void NoClockSites_NoneOfInputMatchesTable_CollapsesToInterceptOnly()
    {
        // INV-05: the methylation map is non-empty but NONE of its CpGs are in the
        // coefficient table, so the weighted sum adds nothing → Y = intercept.
        // No DivideByZero / NaN from a phantom normalisation by "matched site count".
        var coefficients = new Dictionary<string, double>
        {
            ["cgClockA"] = 5.0,
            ["cgClockB"] = -2.0,
        };
        var methylation = new Dictionary<string, double>
        {
            ["cgOtherX"] = 0.9,
            ["cgOtherY"] = 0.1,
            ["cgOtherZ"] = 0.5,
        };

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.3);

        // None matched → identical to the empty-map result at the same intercept.
        double emptyAge = CalculateEpigeneticAge(
            new Dictionary<string, double>(), coefficients, intercept: 0.3);
        age.Should().BeApproximately(emptyAge, 1e-12);
        age.Should().BeApproximately(ExpectedAntiTransform(0.3), 1e-12);
        AssertWellFormedAge(age);
    }

    [Test]
    public void NoClockSites_ExtremeBetaValues_StillIgnored_NoOverflow()
    {
        // Even with absurd β-values on the non-clock CpGs, none contribute, so the
        // transform sees only the intercept and cannot overflow from the inputs.
        var coefficients = new Dictionary<string, double> { ["cgClock"] = 1.0 };
        var methylation = new Dictionary<string, double>
        {
            ["cgNoiseA"] = 1e9,
            ["cgNoiseB"] = -1e9,
        };

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        age.Should().Be(20.0, "no clock site matched → Y = intercept = 0 → 20");
        AssertWellFormedAge(age);
    }

    // ── BE: ALL-METHYLATED (every β = 1.0 → finite, no overflow) ──────────────

    [Test]
    public void AllMethylated_EveryBetaIsOne_PredictorIsInterceptPlusSumOfCoefficients()
    {
        // β = 1 for every clock CpG ⇒ Σ coef_i · 1 = Σ coef_i, so Y = intercept + Σ coef_i.
        // coefficients sum = 0.0127 − 0.0312 + 0.0245 = 0.006; intercept 0.5 → Y = 0.506.
        var coefficients = new Dictionary<string, double>
        {
            ["cg00000029"] = 0.0127,
            ["cg00000165"] = -0.0312,
            ["cg00000363"] = 0.0245,
        };
        var methylation = coefficients.Keys.ToDictionary(k => k, _ => 1.0);

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.5);

        double expectedY = 0.5 + coefficients.Values.Sum();
        age.Should().BeApproximately(ExpectedAntiTransform(expectedY), 1e-12);
        age.Should().BeApproximately(21.0 * expectedY + 20.0, 1e-9);
        AssertWellFormedAge(age);
    }

    [Test]
    public void AllMethylated_LargePositiveCoefficients_FiniteNoOverflow()
    {
        // A large but realistic positive predictor stays on the LINEAR branch, so even
        // big Y produces a finite, large age — never Infinity/NaN. Y = 0 + 100·1 = 100.
        var coefficients = new Dictionary<string, double> { ["cg1"] = 100.0 };
        var methylation = new Dictionary<string, double> { ["cg1"] = 1.0 };

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        age.Should().BeApproximately(21.0 * 100.0 + 20.0, 1e-6); // 2120
        AssertWellFormedAge(age);
    }

    [Test]
    public void AllMethylated_LargeNegativePredictor_ExponentialBranchSaturatesAboveMinusOne()
    {
        // A strongly NEGATIVE Y (all β=1 with a big negative coefficient) takes the
        // exponential branch, where 21·exp(Y) − 1 → −1⁺ as Y → −∞. The result must be
        // finite, just above −1, and NEVER NaN/−Infinity. Y = 0 + (−1000)·1 = −1000.
        var coefficients = new Dictionary<string, double> { ["cg1"] = -1000.0 };
        var methylation = new Dictionary<string, double> { ["cg1"] = 1.0 };

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        age.Should().BeApproximately(-1.0, 1e-9, "exp(−1000) underflows to 0 ⇒ age = −1, the asymptote");
        AssertWellFormedAge(age); // ≥ −1, finite, never NaN/−Infinity
    }

    // ── POSITIVE sanity: the documented worked example (§7.1) ─────────────────

    [Test]
    public void DocumentedWorkedExample_ReproducesPublishedAge()
    {
        // docs §7.1: coefficients {0.0127, −0.0312, 0.0245}, β {0.5, 0.8, 0.3},
        // intercept 0.695507258 → Y = 0.684247258 (≥0) → age = 21·Y + 20 = 34.369192418.
        var coefficients = new Dictionary<string, double>
        {
            ["cg00000029"] = 0.0127,
            ["cg00000165"] = -0.0312,
            ["cg00000363"] = 0.0245,
        };
        var methylation = new Dictionary<string, double>
        {
            ["cg00000029"] = 0.5,
            ["cg00000165"] = 0.8,
            ["cg00000363"] = 0.3,
        };

        double age = CalculateEpigeneticAge(methylation, coefficients, intercept: 0.695507258);

        age.Should().BeApproximately(34.369192418, 1e-6);
        AssertWellFormedAge(age);
    }

    // ── Documented throws: null map, null/empty coefficient table ─────────────

    [Test]
    public void NullMethylationMap_ThrowsArgumentNull()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = 1.0 };
        Action act = () => CalculateEpigeneticAge(null!, coefficients);
        act.Should().Throw<ArgumentNullException>()
            .WithParameterName("methylationAtClockCpGs");
    }

    [Test]
    public void NullCoefficients_ThrowsArgumentNull()
    {
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };
        Action act = () => CalculateEpigeneticAge(methylation, null!);
        act.Should().Throw<ArgumentNullException>()
            .WithParameterName("coefficients");
    }

    [Test]
    public void EmptyCoefficientTable_ThrowsArgument()
    {
        // Docs §3.3 / §6.1: an empty clock has no defined output → ArgumentException.
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };
        var emptyCoefficients = new Dictionary<string, double>();
        Action act = () => CalculateEpigeneticAge(methylation, emptyCoefficients);
        act.Should().Throw<ArgumentException>()
            .WithParameterName("coefficients");
    }

    [Test]
    public void EmptyCoefficientTable_TakesPrecedenceOverEmptyMethylation()
    {
        // Both maps empty: the empty-coefficient guard must still fire (no silent 20.0).
        Action act = () => CalculateEpigeneticAge(
            new Dictionary<string, double>(), new Dictionary<string, double>());
        act.Should().Throw<ArgumentException>().WithParameterName("coefficients");
    }

    // ── INV-05 reinforcement: extra non-clock CpG never changes the age ───────

    [Test]
    public void NonClockCpG_AddedToInput_DoesNotChangeAge()
    {
        var coefficients = new Dictionary<string, double>
        {
            ["cgClock"] = 0.3, // Y contribution 0.3·0.7 = 0.21
        };
        var baseMeth = new Dictionary<string, double> { ["cgClock"] = 0.7 };
        var withExtra = new Dictionary<string, double>
        {
            ["cgClock"] = 0.7,
            ["cgNotInTable"] = 0.9, // must be ignored
        };

        double ageBase = CalculateEpigeneticAge(baseMeth, coefficients, intercept: 0.1);
        double ageExtra = CalculateEpigeneticAge(withExtra, coefficients, intercept: 0.1);

        ageBase.Should().BeApproximately(ageExtra, 1e-12);
        // Exact: Y = 0.1 + 0.3·0.7 = 0.31 → 21·0.31 + 20 = 26.51.
        ageBase.Should().BeApproximately(26.51, 1e-9);
    }

    // ── INV-04: anti.trafo is strictly increasing and continuous at Y=0 ───────

    [Test]
    public void HorvathAntiTransform_IsStrictlyIncreasing_AndContinuousAtZero()
    {
        // Strict monotonicity across the branch boundary, plus continuity at Y=0.
        double prev = double.NegativeInfinity;
        for (double y = -5.0; y <= 5.0; y += 0.01)
        {
            double age = HorvathAntiTransform(y);
            AssertWellFormedAge(age);
            age.Should().BeGreaterThan(prev, $"anti.trafo must be strictly increasing at Y={y}");
            prev = age;
        }

        // Branches meet at (0, 20): left limit and the exact value agree.
        HorvathAntiTransform(0.0).Should().Be(20.0);
        HorvathAntiTransform(-1e-12).Should().BeApproximately(20.0, 1e-9);
        HorvathAntiTransform(1e-12).Should().BeApproximately(20.0, 1e-9);
    }

    // ── BE / robustness: random fuzz — never crash, always well-formed ────────

    [Test]
    [CancelAfter(30000)]
    public void RandomClocks_NeverCrash_AgeAlwaysFiniteAndAboveMinusOne()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);

            // Random clock: 0..10 coefficients (count==0 must throw, handled below).
            int nCoef = rng.Next(0, 11);
            var coefficients = new Dictionary<string, double>();
            for (int i = 0; i < nCoef; i++)
                coefficients[$"cg{i}"] = (rng.NextDouble() - 0.5) * 20.0; // [-10, 10]

            // Random methylation: a mix of clock and non-clock CpGs, β in [0,1] plus
            // occasional out-of-range extremes to stress the predictor.
            int nMeth = rng.Next(0, 15);
            var methylation = new Dictionary<string, double>();
            for (int i = 0; i < nMeth; i++)
            {
                bool clockSite = rng.NextDouble() < 0.5 && nCoef > 0;
                string key = clockSite ? $"cg{rng.Next(nCoef)}" : $"other{i}";
                double beta = rng.Next(4) switch
                {
                    0 => 0.0,
                    1 => 1.0,                       // all-methylated boundary
                    2 => rng.NextDouble(),          // in-range
                    _ => (rng.NextDouble() - 0.5) * 4.0, // out-of-range extreme
                };
                methylation[key] = beta;
            }

            double intercept = (rng.NextDouble() - 0.5) * 10.0;

            if (nCoef == 0)
            {
                Action bad = () => CalculateEpigeneticAge(methylation, coefficients, intercept);
                bad.Should().Throw<ArgumentException>($"seed={seed}: empty clock");
                continue;
            }

            double age = 0;
            Action act = () => age = CalculateEpigeneticAge(methylation, coefficients, intercept);
            act.Should().NotThrow($"seed={seed}");
            AssertWellFormedAge(age);

            // Cross-check against an independent recomputation of Y and F⁻¹.
            double y = intercept;
            foreach (var (cpg, beta) in methylation)
                if (coefficients.TryGetValue(cpg, out double c))
                    y += c * beta;
            age.Should().BeApproximately(ExpectedAntiTransform(y), 1e-6, $"seed={seed}");
        }
    }

    #endregion
}
