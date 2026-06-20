using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology gene-EXPRESSION area — ONCO-EXPR-001.
/// The units under test are the cBioPortal-style expression z-score / outlier /
/// signature-activity entry points implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.CalculateExpressionZScore"/> — z = (r − μ)/σ.
///   • <see cref="OncologyAnalyzer.IdentifyOutlierGenes"/>      — strict ±t outlier rule.
///   • <see cref="OncologyAnalyzer.CalculateSignatureScore"/>   — combined z = (Σz)/√k.
/// (ONCO-IMMUNE-001, the ESTIMATE/ssGSEA immune-signature scorer, is covered
/// separately in OncologyImmuneFuzzTests.cs; this file is scoped to general
/// expression normalization / z-score / differential analysis.)
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts the code
/// NEVER fails in an undisciplined way: no hang, no state corruption, no nonsense
/// output, and no *unhandled* runtime fault (DivideByZero / Inf / silently
/// NaN-corrupted result). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a *documented, intentional* outcome (here: a typed
/// ArgumentException). For z-score normalization the headline hazards are:
///   • a zero-variance / all-equal reference cohort → σ = 0 → z = (r − μ)/0,
///     i.e. DivideByZero / ±Infinity / NaN. The DOCUMENTED guard converts this
///     into an ArgumentException (mirroring the reference NormalizeExpressionLevels
///     .java fatal error) — NO Inf/NaN may leak out. THIS is the headline contract.
///   • a single-sample cohort (n = 1) → sample SD (divisor n − 1) divides by zero
///     → variance undefined. The DOCUMENTED guard throws ArgumentException; there
///     must be no DivideByZero or variance-of-one crash.
///   • a NaN expression value / NaN cohort entry → must not crash, hang, or be
///     mistaken for a valid finite z-score (per IEEE-754 it propagates as NaN;
///     the σ = 0 guard must NOT misfire on a NaN σ and silently emit a finite z).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-EXPR-001 — Tumor gene expression outlier / signature score (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 120.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення (0, single sample, all-equal,
///     NaN). Targets (checklist row 120): "zero variance, single sample,
///     all-equal expression, NaN".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// docs/algorithms/Oncology/Tumor_Gene_Expression_Outlier.md:
///   • z = (r − μ)/σ, μ = (Σrᵢ)/n, σ = √(Σ(rᵢ−μ)²/(n−1)) — sample SD,
///       divisor (n − 1)                                            (§2.2, INV-01..03)
///   • z = 0 when value = cohort mean                               (INV-01, §6.1)
///   • Outlier iff z > +t (Over) or z < −t (Under), strict;
///       |z| = t is NOT an outlier; default t = 2                   (INV-04, §6.1)
///   • Combined signature score a = (Σz)/√k; k member z-scores      (§2.2, INV-05)
///   • k = 1 ⇒ a = z₁                                               (INV-06, §6.1)
///   • DEGENERATE handling (§3.3 / §6.1):
///       – zero-spread (σ = 0) cohort      → ArgumentException
///       – cohort size n < 2 (single)      → ArgumentException (sample SD undefined)
///       – null cohort / dictionaries      → ArgumentNullException
///       – non-positive threshold          → ArgumentOutOfRangeException
///       – sampled gene with no cohort     → ArgumentException
///       – empty signature (k = 0)         → ArgumentException
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyExpressionFuzzTests
{
    private const double Tol = 1e-9;

    // The canonical worked-example cohort from the algorithm doc §7.1:
    //   {2,2,4,6,6}: μ = 4, Σ(rᵢ−μ)² = 16, sample var = 16/(5−1) = 4, σ = 2.
    private static IReadOnlyList<double> WorkedCohort() => new[] { 2.0, 2.0, 4.0, 6.0, 6.0 };

    // ── Well-formed-result helper ─────────────────────────────────────────────
    // A z-score that the contract DEFINES must be finite (never NaN/±Inf): that is
    // exactly what proves the σ = 0 guard fired and no DivideByZero leaked. Used on
    // every path where the inputs are finite and non-degenerate.
    private static void AssertFiniteZScore(double z)
    {
        double.IsNaN(z).Should().BeFalse("a z-score over a finite, non-degenerate cohort must never be NaN");
        double.IsInfinity(z).Should().BeFalse("a z-score over a finite, non-degenerate cohort must be finite");
    }

    #region ONCO-EXPR-001 — Expression z-score / outlier / signature score

    // ── POSITIVE sanity: documented z-scores, fold-change, ordering ───────────

    [Test]
    public void ZScore_WorkedExample_MatchesHandComputedValues()
    {
        // Doc §7.1: cohort {2,2,4,6,6} → μ=4, σ=2 (sample SD).
        //   value 10 → z = (10−4)/2 = 3.0 (Over)
        //   value  4 → z = 0          (= mean, INV-01)
        //   value −1 → z = (−1−4)/2 = −2.5 (Under)
        //   value  8 → z = (8−4)/2 = 2.0   (boundary → not an outlier)
        var cohort = WorkedCohort();

        double zOver = CalculateExpressionZScore(10.0, cohort);
        double zMean = CalculateExpressionZScore(4.0, cohort);
        double zUnder = CalculateExpressionZScore(-1.0, cohort);
        double zBoundary = CalculateExpressionZScore(8.0, cohort);

        AssertFiniteZScore(zOver);
        AssertFiniteZScore(zMean);
        AssertFiniteZScore(zUnder);
        AssertFiniteZScore(zBoundary);

        zOver.Should().BeApproximately(3.0, Tol);
        zMean.Should().BeApproximately(0.0, Tol, "value = cohort mean ⇒ z = 0 (INV-01)");
        zUnder.Should().BeApproximately(-2.5, Tol);
        zBoundary.Should().BeApproximately(2.0, Tol);

        // Sign of the z-score encodes over/under-expression direction (INV-02 monotone).
        zOver.Should().BePositive("a value above the mean has a positive z-score (over-expression)");
        zUnder.Should().BeNegative("a value below the mean has a negative z-score (under-expression)");
    }

    [Test]
    public void ZScore_UsesSampleStandardDeviation_NotPopulation()
    {
        // Sample SD divisor (n−1) is the contract (cBioPortal NormalizeExpressionLevels.java).
        // For {2,2,4,6,6}, value 6: sample SD = 2 → z = (6−4)/2 = 1.0.
        // A population-SD impl (divisor n) would give σ = √(16/5) ≈ 1.78885 and z ≈ 1.118.
        double z = CalculateExpressionZScore(6.0, WorkedCohort());

        z.Should().BeApproximately(1.0, Tol,
            "sample SD (divisor n−1) is the documented contract; population SD would give ≈1.118");
    }

    [Test]
    public void SignatureScore_WorkedExample_MatchesDocumentedCombinedZ()
    {
        // Doc §7.1: z = {3,1,−1,1} → a = (3+1−1+1)/√4 = 4/2 = 2.0.
        double a = CalculateSignatureScore(new[] { 3.0, 1.0, -1.0, 1.0 });

        a.Should().BeApproximately(2.0, Tol);
    }

    [Test]
    public void OutlierRule_StrictThreshold_DocumentedDirections()
    {
        // Doc §6.1 / INV-04: |z| = t is NOT an outlier; z > +t Over, z < −t Under.
        // Cohort {2,2,4,6,6} (μ=4, σ=2). values: BOUNDARY z=2.0 (not outlier),
        // OVER z=3.0, UNDER z=−2.5, NEUTRAL z=0.
        var sample = new Dictionary<string, double>
        {
            ["BOUNDARY"] = 8.0,   // z = 2.0  → exactly threshold → excluded
            ["OVER"] = 10.0,      // z = 3.0  → Over
            ["UNDER"] = -1.0,     // z = −2.5 → Under
            ["NEUTRAL"] = 4.0,    // z = 0    → excluded
        };
        var cohorts = new Dictionary<string, IReadOnlyList<double>>
        {
            ["BOUNDARY"] = WorkedCohort(),
            ["OVER"] = WorkedCohort(),
            ["UNDER"] = WorkedCohort(),
            ["NEUTRAL"] = WorkedCohort(),
        };

        IReadOnlyList<ExpressionOutlier> outliers = IdentifyOutlierGenes(sample, cohorts);

        outliers.Select(o => o.Gene).Should().BeEquivalentTo(new[] { "OVER", "UNDER" },
            "|z| = threshold (BOUNDARY) and z = 0 (NEUTRAL) are not outliers (strict rule, INV-04)");

        outliers.Single(o => o.Gene == "OVER").Direction.Should().Be(ExpressionDirection.Over);
        outliers.Single(o => o.Gene == "UNDER").Direction.Should().Be(ExpressionDirection.Under);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE — Zero variance / all-equal expression (HEADLINE: σ = 0 ⇒ throw, no Inf/NaN)
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void BE_AllEqualCohort_ThrowsArgumentException_NoInfNaNLeak()
    {
        // Every cohort value identical ⇒ σ = 0 ⇒ z = (r − μ)/0 hazard. The
        // documented guard (§6.1, mirrors NormalizeExpressionLevels.java) converts
        // this to ArgumentException — it must NOT return ±Infinity or NaN.
        var cohort = new[] { 5.0, 5.0, 5.0, 5.0 };

        Action act = () => CalculateExpressionZScore(9.0, cohort);

        act.Should().Throw<ArgumentException>(
            "a zero-variance (all-equal) cohort has σ = 0 and no defined z-score");
    }

    [Test]
    [CancelAfter(30_000)]
    public void BE_RandomAllEqualCohorts_NeverLeakInfOrNaN()
    {
        // Fuzz a range of constant cohorts (including 0, ±large, negative) and
        // assert EVERY one is rejected by the σ = 0 guard rather than producing a
        // non-finite z. The σ = 0 case is the single most dangerous BE input.
        var rng = new Random(120_001);

        for (int i = 0; i < 400; i++)
        {
            double constant = (rng.NextDouble() - 0.5) * Math.Pow(10, rng.Next(0, 9));
            int n = rng.Next(2, 12);
            var cohort = Enumerable.Repeat(constant, n).ToArray();
            double value = (rng.NextDouble() - 0.5) * 1e6;

            Action act = () => CalculateExpressionZScore(value, cohort);

            act.Should().Throw<ArgumentException>(
                "constant cohort {0}×{1} ⇒ σ=0 ⇒ must throw, never emit Inf/NaN", constant, n);
        }
    }

    [Test]
    public void BE_AllEqualValueEqualsConstant_StillThrows_NoZeroOverZeroNaN()
    {
        // Adversarial: value EQUALS the constant cohort, so the numerator (r − μ)
        // is 0 too — a naive impl would compute 0/0 = NaN and quietly return it.
        // The guard must still throw (σ = 0 is checked BEFORE the division).
        var cohort = new[] { 7.0, 7.0, 7.0 };

        Action act = () => CalculateExpressionZScore(7.0, cohort);

        act.Should().Throw<ArgumentException>("σ = 0 is rejected before any (0/0) division can occur");
    }

    [Test]
    public void BE_AllEqualCohort_InOutlierDetection_PropagatesThrow()
    {
        // IdentifyOutlierGenes delegates to CalculateExpressionZScore, so a
        // degenerate cohort must surface as the same documented ArgumentException
        // (§3.3 "or a cohort is degenerate"), not a swallowed Inf/NaN outlier.
        var sample = new Dictionary<string, double> { ["G1"] = 3.0 };
        var cohorts = new Dictionary<string, IReadOnlyList<double>>
        {
            ["G1"] = new[] { 2.0, 2.0, 2.0 },
        };

        Action act = () => IdentifyOutlierGenes(sample, cohorts);

        act.Should().Throw<ArgumentException>();
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE — Single sample (n = 1 / n = 0): sample SD undefined ⇒ throw
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void BE_SingleSampleCohort_ThrowsArgumentException_NoVarianceOfOneCrash()
    {
        // n = 1 ⇒ sample-SD divisor (n − 1) = 0 ⇒ variance undefined. Documented
        // guard throws ArgumentException (§6.1 "cohort size 1"); must not
        // DivideByZero or compute a spurious variance-of-one.
        Action act = () => CalculateExpressionZScore(42.0, new[] { 42.0 });

        act.Should().Throw<ArgumentException>("a single-sample cohort has an undefined sample SD");
    }

    [Test]
    public void BE_EmptyCohort_ThrowsArgumentException()
    {
        // n = 0 ⇒ no mean, no SD. Same n < 2 guard.
        Action act = () => CalculateExpressionZScore(1.0, Array.Empty<double>());

        act.Should().Throw<ArgumentException>("an empty cohort cannot define a z-score");
    }

    [Test]
    [CancelAfter(30_000)]
    public void BE_RandomSingleSampleCohorts_AlwaysThrow()
    {
        var rng = new Random(120_002);

        for (int i = 0; i < 200; i++)
        {
            double only = (rng.NextDouble() - 0.5) * 1e7;
            double value = (rng.NextDouble() - 0.5) * 1e7;

            Action act = () => CalculateExpressionZScore(value, new[] { only });

            act.Should().Throw<ArgumentException>("single-sample cohort {0} ⇒ sample SD undefined", only);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE — NaN inputs: never crash/hang; never a finite z masquerading for NaN
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void BE_NaNValue_OverValidCohort_PropagatesNaN_NotFinite()
    {
        // A NaN sample value over a well-formed cohort: (NaN − μ)/σ = NaN per
        // IEEE-754. The contract does not silently coerce it to a finite z; it
        // must NOT throw a non-Argument runtime fault and must NOT return a
        // plausible-looking finite number.
        double z = CalculateExpressionZScore(double.NaN, WorkedCohort());

        double.IsNaN(z).Should().BeTrue("a NaN expression value yields a NaN z-score (IEEE-754 propagation)");
    }

    [Test]
    public void BE_NaNInCohort_DoesNotMisfireZeroVarianceGuard_PropagatesNaN()
    {
        // A NaN cohort entry makes μ and σ NaN. The σ = 0 guard compares (NaN == 0)
        // which is FALSE, so the method does NOT throw — it returns (value − NaN)/NaN
        // = NaN. The danger this pins: the guard must not be tricked into emitting a
        // FINITE z from a NaN-poisoned cohort. The result is NaN, never a finite number.
        var cohort = new[] { 2.0, double.NaN, 6.0, 4.0 };

        double z = CalculateExpressionZScore(5.0, cohort);

        double.IsNaN(z).Should().BeTrue("a NaN-poisoned cohort produces a NaN z-score, never a finite one");
        double.IsInfinity(z).Should().BeFalse();
    }

    [Test]
    public void BE_NaNZScore_InSignatureScore_PropagatesNaN()
    {
        // Combined z a = (Σz)/√k. A NaN member z-score propagates to a NaN
        // activity; it must not be silently dropped to yield a finite (wrong) score.
        double a = CalculateSignatureScore(new[] { 1.0, double.NaN, -1.0 });

        double.IsNaN(a).Should().BeTrue("a NaN member z-score yields a NaN combined signature score");
    }

    [Test]
    public void BE_NaNValue_InOutlierDetection_NotClassifiedAsOutlier()
    {
        // A NaN z-score satisfies neither z > +t nor z < −t (all comparisons with
        // NaN are false), so a NaN-valued gene is silently excluded — it is NOT a
        // spurious Over/Under outlier. (No crash; finite genes still classified.)
        var sample = new Dictionary<string, double>
        {
            ["NAN"] = double.NaN, // z = NaN  → excluded (NaN > t and NaN < −t both false)
            ["OVER"] = 10.0,      // z = 3.0  → Over
        };
        var cohorts = new Dictionary<string, IReadOnlyList<double>>
        {
            ["NAN"] = WorkedCohort(),
            ["OVER"] = WorkedCohort(),
        };

        IReadOnlyList<ExpressionOutlier> outliers = IdentifyOutlierGenes(sample, cohorts);

        outliers.Select(o => o.Gene).Should().BeEquivalentTo(new[] { "OVER" },
            "a NaN z-score is not > +t nor < −t, so it is not reported as an outlier");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE — Other documented degenerate guards (null / threshold / missing / empty)
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void BE_NullCohort_ThrowsArgumentNullException()
    {
        Action act = () => CalculateExpressionZScore(1.0, null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void BE_NullDictionaries_ThrowArgumentNullException()
    {
        Action sampleNull = () => IdentifyOutlierGenes(
            null!, new Dictionary<string, IReadOnlyList<double>>());
        Action cohortsNull = () => IdentifyOutlierGenes(
            new Dictionary<string, double>(), null!);

        sampleNull.Should().Throw<ArgumentNullException>();
        cohortsNull.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void BE_NullSignature_ThrowsArgumentNullException()
    {
        Action act = () => CalculateSignatureScore(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void BE_EmptySignature_ThrowsArgumentException()
    {
        // k = 0 ⇒ a = (Σz)/√0 undefined (and 0/0). Documented guard throws (§6.1).
        Action act = () => CalculateSignatureScore(Array.Empty<double>());
        act.Should().Throw<ArgumentException>("an empty signature has no defined combined z-score");
    }

    [Test]
    public void BE_SingletonSignature_EqualsThatZScore()
    {
        // INV-06: k = 1 ⇒ a = z₁/√1 = z₁. The √0 hazard only bites at k = 0.
        double a = CalculateSignatureScore(new[] { -1.7 });
        a.Should().BeApproximately(-1.7, Tol);
        AssertFiniteZScore(a);
    }

    [Test]
    [TestCase(0.0)]
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    public void BE_NonPositiveThreshold_ThrowsArgumentOutOfRange(double threshold)
    {
        // threshold ≤ 0 (or NaN, since NaN ≤ 0 is false but the contract requires
        // a positive cutoff) — boundary value 0 and -1 must be rejected (§3.3).
        var sample = new Dictionary<string, double> { ["G"] = 10.0 };
        var cohorts = new Dictionary<string, IReadOnlyList<double>> { ["G"] = WorkedCohort() };

        Action act = () => IdentifyOutlierGenes(sample, cohorts, threshold);

        if (double.IsNaN(threshold))
        {
            // NaN ≤ 0 is false, so the guard does not fire; the run must still not
            // crash and must yield no outliers (NaN comparisons are all false).
            IReadOnlyList<ExpressionOutlier> outliers = IdentifyOutlierGenes(sample, cohorts, threshold);
            outliers.Should().BeEmpty("a NaN threshold makes every z > t / z < −t comparison false");
        }
        else
        {
            act.Should().Throw<ArgumentOutOfRangeException>("the outlier threshold must be positive");
        }
    }

    [Test]
    public void BE_SampledGeneWithoutCohort_ThrowsArgumentException()
    {
        var sample = new Dictionary<string, double> { ["KNOWN"] = 10.0, ["MISSING"] = 5.0 };
        var cohorts = new Dictionary<string, IReadOnlyList<double>> { ["KNOWN"] = WorkedCohort() };

        Action act = () => IdentifyOutlierGenes(sample, cohorts);

        act.Should().Throw<ArgumentException>("every sampled gene must have a reference cohort (§3.3)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE — Boundary-magnitude fuzz on well-formed cohorts: z stays finite & exact
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    [CancelAfter(30_000)]
    public void BE_RandomNonDegenerateCohorts_ZScoreFiniteAndMatchesFormula()
    {
        // Random non-degenerate cohorts (guaranteed spread) across a wide magnitude
        // range. Every z must be finite and equal the independently re-derived
        // z = (r − μ)/σ with σ the sample SD. This is the "z is well-formed when
        // defined" invariant from the doc, fuzzed.
        var rng = new Random(120_003);

        for (int i = 0; i < 500; i++)
        {
            int n = rng.Next(2, 16);
            double scale = Math.Pow(10, rng.Next(-3, 6));
            var cohort = new double[n];
            for (int j = 0; j < n; j++)
            {
                cohort[j] = (rng.NextDouble() - 0.5) * scale;
            }

            // Force non-zero spread (avoid the σ = 0 path, which is tested separately).
            cohort[0] += scale + 1.0;

            double mean = cohort.Average();
            double ss = cohort.Sum(x => (x - mean) * (x - mean));
            double sd = Math.Sqrt(ss / (n - 1));
            if (sd == 0.0)
            {
                continue; // extremely unlikely; skip degenerate
            }

            double value = (rng.NextDouble() - 0.5) * scale;
            double expected = (value - mean) / sd;

            double z = CalculateExpressionZScore(value, cohort);

            AssertFiniteZScore(z);
            z.Should().BeApproximately(expected, Math.Max(1e-6, Math.Abs(expected) * 1e-9),
                "z must equal (value − μ)/σ with σ the sample SD (divisor n−1)");
        }
    }

    [Test]
    public void BE_ZeroValuedSample_OverNonDegenerateCohort_FiniteZScore()
    {
        // value = 0 boundary (a real expression floor) over the worked cohort:
        // z = (0 − 4)/2 = −2.0 (exactly the threshold → not an outlier).
        double z = CalculateExpressionZScore(0.0, WorkedCohort());

        AssertFiniteZScore(z);
        z.Should().BeApproximately(-2.0, Tol);
    }

    #endregion
}
