// EPIGEN-AGE-001 — Epigenetic Age Estimation (Horvath DNA methylation clock)
// Evidence: docs/Evidence/EPIGEN-AGE-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-AGE-001.md
// Source: Horvath S (2013). DNA methylation age of human tissues and cell types. Genome Biology 14:R115.
//         Reference R: aldringsvitenskap/epigeneticclock horvath2013.R (anti.trafo, adult.age=20)
//         and StepwiseAnalysis.R (predictedAge = anti.trafo(intercept + meth·coef)).

namespace Seqeron.Genomics.Tests;

using System;
using System.Collections.Generic;
using Seqeron.Genomics.Annotation;

[TestFixture]
public class EpigeneticsAnalyzer_CalculateEpigeneticAge_Tests
{
    #region CalculateEpigeneticAge

    // M1 — linear branch: intercept + 3 CpGs. Y = 0.695507258 + 0.0127*0.5 - 0.0312*0.8 + 0.0245*0.3
    //      = 0.684247258 (>= 0) → anti.trafo = 21*Y + 20 = 34.369192418 (StepwiseAnalysis.R + anti.trafo).
    [Test]
    public void CalculateEpigeneticAge_LinearBranchWithIntercept_Returns21YPlus20()
    {
        const double intercept = 0.695507258;
        var coefficients = new Dictionary<string, double>
        {
            ["cg1"] = 0.0127,
            ["cg2"] = -0.0312,
            ["cg3"] = 0.0245,
        };
        var methylation = new Dictionary<string, double>
        {
            ["cg1"] = 0.5,
            ["cg2"] = 0.8,
            ["cg3"] = 0.3,
        };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, intercept);

        Assert.That(age, Is.EqualTo(34.369192418).Within(1e-9),
            "Y=0.684247258 (>=0) → linear branch 21*Y+20 per Horvath anti.trafo (sources #2,#3)");
    }

    // M2 — boundary Y = 0: single coefficient times beta 0 with zero intercept → Y=0 → 21*0+20 = 20.
    [Test]
    public void CalculateEpigeneticAge_LinearPredictorZero_ReturnsAdultAge20()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = 0.5 };
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.0 };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        Assert.That(age, Is.EqualTo(20.0).Within(1e-10),
            "At Y=0 the x<0 test is false → linear branch yields exactly adult.age = 20 (anti.trafo, source #2)");
    }

    // M3 — negative branch: coef -2.0 * beta 0.5 = -1.0, intercept 0 → Y=-1.0 → 21*e^-1 - 1.
    [Test]
    public void CalculateEpigeneticAge_NegativeLinearPredictor_ReturnsExponentialBranch()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = -2.0 };
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);

        Assert.That(age, Is.EqualTo(6.7254682646002895).Within(1e-12),
            "Y=-1.0 (<0) → exponential branch 21*exp(-1)-1 per Horvath anti.trafo (source #2)");
    }

    // M4 — a CpG present in the methylation map but absent from the coefficient table is ignored:
    //      only clock CpGs enter the weighted sum (StepwiseAnalysis.R matrix product, source #3).
    [Test]
    public void CalculateEpigeneticAge_NonClockCpG_IsIgnored()
    {
        const double intercept = 0.1;
        var coefficients = new Dictionary<string, double> { ["cg_clock"] = 0.5 };
        var withExtra = new Dictionary<string, double>
        {
            ["cg_clock"] = 0.4,
            ["cg_not_in_clock"] = 0.9, // huge value that must not affect the result
        };
        var withoutExtra = new Dictionary<string, double> { ["cg_clock"] = 0.4 };

        double ageWithExtra = EpigeneticsAnalyzer.CalculateEpigeneticAge(withExtra, coefficients, intercept);
        double ageWithoutExtra = EpigeneticsAnalyzer.CalculateEpigeneticAge(withoutExtra, coefficients, intercept);

        // Exact sourced value: only cg_clock enters the sum → Y = 0.1 + 0.5*0.4 = 0.3 (>=0)
        // → linear branch 21*0.3 + 20 = 26.3 (anti.trafo, sources #2,#3). Asserting the exact value
        // (not just equality of the two calls) prevents a consistent-but-wrong predictor from passing.
        Assert.Multiple(() =>
        {
            Assert.That(ageWithExtra, Is.EqualTo(26.3).Within(1e-12),
                "Y=0.3 (>=0) → 21*0.3+20=26.3; non-clock CpG must not change the sourced age (sources #2,#3)");
            Assert.That(ageWithExtra, Is.EqualTo(ageWithoutExtra).Within(1e-12),
                "CpGs without a clock coefficient contribute nothing (source #3); both calls must match exactly");
        });
    }

    // S1 — empty methylation map: no CpG contributions, so Y = intercept only. intercept 1.0 → 21*1+20 = 41.
    [Test]
    public void CalculateEpigeneticAge_EmptyMethylationMap_ReturnsAntiTransformOfIntercept()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = 0.5 };
        var methylation = new Dictionary<string, double>();

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, intercept: 1.0);

        Assert.That(age, Is.EqualTo(41.0).Within(1e-10),
            "No CpGs → Y=intercept=1.0 → 21*1+20=41 (intercept still applies, source #3)");
    }

    // S3 — monotonicity (INV-04): increasing the linear predictor never decreases age, across the Y=0 boundary.
    [Test]
    public void CalculateEpigeneticAge_IncreasingPredictor_IsMonotonic()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = 1.0 };
        // beta values chosen so Y spans negative, zero, and positive (intercept 0).
        double[] betas = { -1.0, -0.25, 0.0, 0.25, 1.0 };
        double previous = double.NegativeInfinity;

        Assert.Multiple(() =>
        {
            foreach (double beta in betas)
            {
                var methylation = new Dictionary<string, double> { ["cg1"] = beta };
                double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients, intercept: 0.0);
                Assert.That(age, Is.GreaterThan(previous),
                    $"anti.trafo is strictly increasing in Y; age at Y={beta} must exceed the previous (INV-04)");
                previous = age;
            }
        });
    }

    // M6 — null methylation map → ArgumentNullException.
    [Test]
    public void CalculateEpigeneticAge_NullMethylation_Throws()
    {
        var coefficients = new Dictionary<string, double> { ["cg1"] = 0.5 };

        Assert.That(() => EpigeneticsAnalyzer.CalculateEpigeneticAge(null!, coefficients),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null methylation map is invalid input (contract)");
    }

    // M7 — null coefficient table → ArgumentNullException.
    [Test]
    public void CalculateEpigeneticAge_NullCoefficients_Throws()
    {
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };

        Assert.That(() => EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null coefficient table is invalid input; defaults are never fabricated (contract)");
    }

    // M8 — empty coefficient table → ArgumentException (an empty clock has no defined output).
    [Test]
    public void CalculateEpigeneticAge_EmptyCoefficients_Throws()
    {
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };
        var empty = new Dictionary<string, double>();

        Assert.That(() => EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, empty),
            NUnit.Framework.Throws.ArgumentException,
            "An empty clock coefficient table has no defined age output (contract)");
    }

    #endregion

    #region HorvathAntiTransform

    // M5 — published anti.trafo in isolation, negative branch: x=-2.5 → 21*exp(-2.5)-1.
    [Test]
    public void HorvathAntiTransform_NegativeValue_ReturnsExponentialBranch()
    {
        double age = EpigeneticsAnalyzer.HorvathAntiTransform(-2.5);

        Assert.That(age, Is.EqualTo(0.7237849711018749).Within(1e-12),
            "x=-2.5 (<0) → 21*exp(-2.5)-1 per Horvath anti.trafo (source #2)");
    }

    // S2 — anti.trafo boundary at x=0 → exactly 20 (INV-03).
    [Test]
    public void HorvathAntiTransform_Zero_ReturnsAdultAge20()
    {
        double age = EpigeneticsAnalyzer.HorvathAntiTransform(0.0);

        Assert.That(age, Is.EqualTo(20.0).Within(1e-10),
            "x=0 → linear branch 21*0+20 = 20 (adult.age boundary, source #2)");
    }

    // C1 — anti.trafo positive branch: x=1.0 → 21*1+20 = 41.
    [Test]
    public void HorvathAntiTransform_PositiveValue_ReturnsLinearBranch()
    {
        double age = EpigeneticsAnalyzer.HorvathAntiTransform(1.0);

        Assert.That(age, Is.EqualTo(41.0).Within(1e-10),
            "x=1.0 (>=0) → linear branch 21*1+20 = 41 (source #2)");
    }

    #endregion
}
