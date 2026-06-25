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

    #region Embedded Horvath 2013 multi-tissue clock

    // E1 — table integrity: the embedded clock has exactly 353 CpG coefficients and the
    //      published intercept. Source: Horvath (2013) Additional file 3 (CoefficientTraining).
    [Test]
    public void HorvathMultiTissueClock_TableIntegrity_Has353CpGsAndPublishedIntercept()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.HorvathMultiTissueCoefficients.Count, Is.EqualTo(353),
                "Horvath 2013 multi-tissue clock selects exactly 353 CpGs (Additional file 3)");
            Assert.That(EpigeneticsAnalyzer.HorvathMultiTissueIntercept, Is.EqualTo(0.695507258).Within(1e-12),
                "Embedded intercept must equal CoefficientTraining[1] = 0.695507258 (Additional file 3)");
        });
    }

    // E2 — spot-check four named probe coefficients against the published table
    //      (verified byte-identical across Springer supplement + GitHub mirror).
    [Test]
    public void HorvathMultiTissueClock_NamedProbes_MatchPublishedCoefficients()
    {
        var c = EpigeneticsAnalyzer.HorvathMultiTissueCoefficients;
        Assert.Multiple(() =>
        {
            Assert.That(c["cg00075967"], Is.EqualTo(0.12933661).Within(1e-12),
                "cg00075967 CoefficientTraining = 0.12933661 (Additional file 3)");
            Assert.That(c["cg00374717"], Is.EqualTo(0.005017857).Within(1e-12),
                "cg00374717 CoefficientTraining = 0.005017857 (Additional file 3)");
            Assert.That(c["cg00864867"], Is.EqualTo(1.59976405).Within(1e-12),
                "cg00864867 CoefficientTraining = 1.59976405 (Additional file 3)");
            Assert.That(c["cg09809672"], Is.EqualTo(-0.391318905).Within(1e-12),
                "cg09809672 (EDARADD) CoefficientTraining = -0.391318905 (Additional file 3)");
            Assert.That(c["cg27544190"], Is.EqualTo(-0.869124446).Within(1e-12),
                "cg27544190 CoefficientTraining = -0.869124446 (last supplement row, Additional file 3)");
        });
    }

    // E3 — default overload, empty methylation map → DNAmAge = anti.trafo(intercept).
    //      LP = 0.695507258 (>=0) → 21*LP + 20 = 34.605652418 (Horvath anti.trafo, adult.age=20).
    [Test]
    public void CalculateEpigeneticAge_BuiltInClock_EmptyMethylation_ReturnsAntiTransformOfIntercept()
    {
        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(new Dictionary<string, double>());

        Assert.That(age, Is.EqualTo(34.605652418).Within(1e-9),
            "No CpGs → LP = intercept 0.695507258 → 21*LP+20 = 34.605652418 (built-in clock + anti.trafo)");
    }

    // E4 — default overload, linear branch: one clock CpG with β=1 contributes its exact
    //      published coefficient. LP = 0.695507258 + 1.59976405*1 = 2.295271308 (>=0)
    //      → 21*LP + 20 = 68.200697468.
    [Test]
    public void CalculateEpigeneticAge_BuiltInClock_SingleProbe_ReturnsExactLinearBranch()
    {
        var methylation = new Dictionary<string, double> { ["cg00864867"] = 1.0 };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);

        Assert.That(age, Is.EqualTo(68.200697468).Within(1e-9),
            "LP = 0.695507258 + 1.59976405 = 2.295271308 → 21*LP+20 = 68.200697468 (built-in coef)");
    }

    // E5 — default overload, exponential branch: two negative-weight clock CpGs at β=1 push
    //      LP below 0. LP = 0.695507258 + (-0.391318905) + (-0.869124446) = -0.564936093 (<0)
    //      → 21*exp(LP) - 1 = 10.936325872311789.
    [Test]
    public void CalculateEpigeneticAge_BuiltInClock_NegativePredictor_ReturnsExponentialBranch()
    {
        var methylation = new Dictionary<string, double>
        {
            ["cg09809672"] = 1.0,
            ["cg27544190"] = 1.0,
        };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);

        Assert.That(age, Is.EqualTo(10.936325872311789).Within(1e-9),
            "LP = -0.564936093 (<0) → 21*exp(LP)-1 = 10.9363... (built-in coef, exponential branch)");
    }

    // E6 — default overload equals the explicit overload called with the built-in table+intercept
    //      (the parameterless form is a thin wrapper; it must not alter the result).
    [Test]
    public void CalculateEpigeneticAge_BuiltInClock_EqualsExplicitOverloadWithSameTable()
    {
        var methylation = new Dictionary<string, double>
        {
            ["cg00075967"] = 0.42,
            ["cg00864867"] = 0.61,
            ["cg09809672"] = 0.33,
        };

        double viaDefault = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);
        double viaExplicit = EpigeneticsAnalyzer.CalculateEpigeneticAge(
            methylation,
            EpigeneticsAnalyzer.HorvathMultiTissueCoefficients,
            EpigeneticsAnalyzer.HorvathMultiTissueIntercept);

        Assert.That(viaDefault, Is.EqualTo(viaExplicit).Within(1e-12),
            "Parameterless overload must delegate to the explicit one with the built-in table/intercept");
    }

    // E7 — default overload rejects a null methylation map (same contract as the explicit overload).
    [Test]
    public void CalculateEpigeneticAge_BuiltInClock_NullMethylation_Throws()
    {
        Assert.That(() => EpigeneticsAnalyzer.CalculateEpigeneticAge(null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null methylation map is invalid input for the built-in clock (contract)");
    }

    #endregion

    #region Embedded Horvath 2018 skin & blood clock

    // SB1 — table integrity: 391 CpGs and the published intercept -0.447119319.
    //       Source: Horvath et al. (2018) Aging 10(7):1758-1775, Supplementary Dataset 2 (Coef).
    [Test]
    public void HorvathSkinBloodClock_TableIntegrity_Has391CpGsAndPublishedIntercept()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.HorvathSkinBloodCoefficients.Count, Is.EqualTo(391),
                "Horvath 2018 skin & blood clock selects exactly 391 CpGs (Supplementary Dataset 2)");
            Assert.That(EpigeneticsAnalyzer.HorvathSkinBloodIntercept, Is.EqualTo(-0.447119319).Within(1e-12),
                "Embedded intercept must equal the '(Intercept)' Coef = -0.447119319 (Supplementary Dataset 2)");
        });
    }

    // SB2 — spot-check named probe coefficients (verified numerically identical: aging-us SD5 + biolearn Horvath2.csv).
    [Test]
    public void HorvathSkinBloodClock_NamedProbes_MatchPublishedCoefficients()
    {
        var c = EpigeneticsAnalyzer.HorvathSkinBloodCoefficients;
        Assert.Multiple(() =>
        {
            Assert.That(c["cg12140144"], Is.EqualTo(0.36318117).Within(1e-12),
                "cg12140144 Coef = 0.36318117 (first CpG row, Supplementary Dataset 2)");
            Assert.That(c["cg26933021"], Is.EqualTo(-0.09050009).Within(1e-12),
                "cg26933021 Coef = -0.09050009 (Supplementary Dataset 2)");
            Assert.That(c["cg00431549"], Is.EqualTo(8.83E-06).Within(1e-12),
                "cg00431549 Coef = 8.83E-06 (scientific-notation entry, Supplementary Dataset 2)");
            Assert.That(c["cg22982767"], Is.EqualTo(-0.41219696).Within(1e-12),
                "cg22982767 Coef = -0.41219696 (Supplementary Dataset 2)");
            Assert.That(c["cg01892695"], Is.EqualTo(-0.250618874).Within(1e-12),
                "cg01892695 Coef = -0.250618874 (last supplement row, Supplementary Dataset 2)");
        });
    }

    // SB3 — built-in clock, empty methylation → DNAmAge = anti.trafo(intercept).
    //       intercept = -0.447119319 (<0) → 21*exp(-0.447119319)-1 = 12.428819664840216.
    [Test]
    public void CalculateSkinBloodAge_BuiltInClock_EmptyMethylation_ReturnsAntiTransformOfIntercept()
    {
        double age = EpigeneticsAnalyzer.CalculateSkinBloodAge(new Dictionary<string, double>());

        Assert.That(age, Is.EqualTo(12.428819664840216).Within(1e-9),
            "No CpGs → LP = intercept -0.447119319 (<0) → 21*exp(LP)-1 = 12.4288... (anti.trafo, adult.age=20)");
    }

    // SB4 — built-in clock, single probe cg12140144 (0.36318117) at β=1.
    //       LP = -0.447119319 + 0.36318117 = -0.08393814900000002 (<0) → 21*exp(LP)-1 = 18.309250637525345.
    [Test]
    public void CalculateSkinBloodAge_BuiltInClock_SingleProbe_ReturnsExactExponentialBranch()
    {
        var methylation = new Dictionary<string, double> { ["cg12140144"] = 1.0 };

        double age = EpigeneticsAnalyzer.CalculateSkinBloodAge(methylation);

        Assert.That(age, Is.EqualTo(18.309250637525345).Within(1e-9),
            "LP = -0.447119319 + 0.36318117 = -0.083938149 (<0) → 21*exp(LP)-1 = 18.30925... (built-in coef)");
    }

    // SB5 — built-in clock delegates to the anti.trafo overload with the built-in table/intercept.
    [Test]
    public void CalculateSkinBloodAge_BuiltInClock_EqualsExplicitOverloadWithSameTable()
    {
        var methylation = new Dictionary<string, double>
        {
            ["cg12140144"] = 0.42,
            ["cg26933021"] = 0.61,
            ["cg01892695"] = 0.33,
        };

        double viaDefault = EpigeneticsAnalyzer.CalculateSkinBloodAge(methylation);
        double viaExplicit = EpigeneticsAnalyzer.CalculateEpigeneticAge(
            methylation,
            EpigeneticsAnalyzer.HorvathSkinBloodCoefficients,
            EpigeneticsAnalyzer.HorvathSkinBloodIntercept);

        Assert.That(viaDefault, Is.EqualTo(viaExplicit).Within(1e-12),
            "Skin & blood overload must delegate to CalculateEpigeneticAge with the built-in table/intercept (anti.trafo path)");
    }

    // SB6 — built-in clock rejects a null methylation map (same contract as the explicit overload).
    [Test]
    public void CalculateSkinBloodAge_BuiltInClock_NullMethylation_Throws()
    {
        Assert.That(() => EpigeneticsAnalyzer.CalculateSkinBloodAge(null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "A null methylation map is invalid input for the built-in skin & blood clock (contract)");
    }

    #endregion

    #region Embedded Levine 2018 PhenoAge clock

    // PA1 — table integrity: 513 CpGs and the published intercept 60.664.
    //       Source: Levine et al. (2018) Aging 10(4):573-591, Supplementary Dataset 2 (Weight).
    [Test]
    public void PhenoAgeClock_TableIntegrity_Has513CpGsAndPublishedIntercept()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.PhenoAgeCoefficients.Count, Is.EqualTo(513),
                "Levine 2018 PhenoAge clock uses exactly 513 CpGs (Supplementary Dataset 2)");
            Assert.That(EpigeneticsAnalyzer.PhenoAgeIntercept, Is.EqualTo(60.664).Within(1e-12),
                "Embedded intercept must equal the 'Intercept' Weight = 60.664 (Supplementary Dataset 2)");
        });
    }

    // PA2 — spot-check named probe weights (verified numerically identical: aging-us SD2 + biolearn PhenoAge.csv).
    [Test]
    public void PhenoAgeClock_NamedProbes_MatchPublishedWeights()
    {
        var c = EpigeneticsAnalyzer.PhenoAgeCoefficients;
        Assert.Multiple(() =>
        {
            Assert.That(c["cg15611364"], Is.EqualTo(63.12415047).Within(1e-12),
                "cg15611364 Weight = 63.12415047 (largest positive weight, Supplementary Dataset 2)");
            Assert.That(c["cg17605084"], Is.EqualTo(-44.00939313).Within(1e-12),
                "cg17605084 Weight = -44.00939313 (Supplementary Dataset 2)");
            Assert.That(c["cg26382071"], Is.EqualTo(40.42085373).Within(1e-12),
                "cg26382071 Weight = 40.42085373 (Supplementary Dataset 2)");
            Assert.That(c["cg00503840"], Is.EqualTo(0.002679625).Within(1e-12),
                "cg00503840 Weight = 0.002679625 (Supplementary Dataset 2)");
            Assert.That(c["cg15381313"], Is.EqualTo(0.002078841).Within(1e-12),
                "cg15381313 Weight = 0.002078841 (last supplement row, Supplementary Dataset 2)");
        });
    }

    // PA3 — built-in clock, empty methylation → age = intercept (NO transform). Result = 60.664.
    [Test]
    public void CalculatePhenoAge_BuiltInClock_EmptyMethylation_ReturnsIntercept()
    {
        double age = EpigeneticsAnalyzer.CalculatePhenoAge(new Dictionary<string, double>());

        Assert.That(age, Is.EqualTo(60.664).Within(1e-12),
            "No CpGs and no transform → DNAm PhenoAge = intercept = 60.664 (Levine 2018 linear predictor)");
    }

    // PA4 — built-in clock, single probe cg15611364 (63.12415047) at β=1 → 60.664 + 63.12415047 = 123.78815047.
    [Test]
    public void CalculatePhenoAge_BuiltInClock_SingleProbe_ReturnsExactLinearPredictor()
    {
        var methylation = new Dictionary<string, double> { ["cg15611364"] = 1.0 };

        double age = EpigeneticsAnalyzer.CalculatePhenoAge(methylation);

        Assert.That(age, Is.EqualTo(123.78815047).Within(1e-9),
            "age = 60.664 + 63.12415047*1 = 123.78815047, returned untransformed (Levine 2018)");
    }

    // PA5 — built-in clock, two probes at β=0.5 → 60.664 + 0.5*63.12415047 + 0.5*(-44.00939313) = 70.22137867.
    [Test]
    public void CalculatePhenoAge_BuiltInClock_TwoProbes_ReturnsExactLinearPredictor()
    {
        var methylation = new Dictionary<string, double>
        {
            ["cg15611364"] = 0.5,
            ["cg17605084"] = 0.5,
        };

        double age = EpigeneticsAnalyzer.CalculatePhenoAge(methylation);

        Assert.That(age, Is.EqualTo(70.22137867).Within(1e-9),
            "age = 60.664 + 0.5*63.12415047 + 0.5*(-44.00939313) = 70.22137867 (no transform, Levine 2018)");
    }

    // PA6 — PhenoAge applies NO anti.trafo: the same linear predictor must differ from the Horvath path.
    //       For a non-clock-only input with LP=intercept, PhenoAge returns 60.664 whereas anti.trafo(60.664)
    //       would be 21*60.664+20 = 1293.944 — confirming no transform is applied.
    [Test]
    public void CalculatePhenoAge_DoesNotApplyAntiTransform()
    {
        var empty = new Dictionary<string, double>();
        double phenoAge = EpigeneticsAnalyzer.CalculatePhenoAge(empty);
        double withTransform = EpigeneticsAnalyzer.HorvathAntiTransform(EpigeneticsAnalyzer.PhenoAgeIntercept);

        Assert.Multiple(() =>
        {
            Assert.That(phenoAge, Is.EqualTo(60.664).Within(1e-12),
                "PhenoAge returns the linear predictor directly (= intercept here), no transform (Levine 2018)");
            Assert.That(phenoAge, Is.Not.EqualTo(withTransform).Within(1e-6),
                "Applying anti.trafo would give 21*60.664+20 = 1293.944; PhenoAge must NOT do this");
        });
    }

    // PA7 — caller-supplied overload: non-clock CpG ignored; exact linear value with no transform.
    [Test]
    public void CalculatePhenoAge_CallerSupplied_NonClockCpGIgnored_NoTransform()
    {
        const double intercept = 10.0;
        var coefficients = new Dictionary<string, double> { ["cg_clock"] = 2.0 };
        var methylation = new Dictionary<string, double>
        {
            ["cg_clock"] = 0.5,
            ["cg_not_in_clock"] = 99.0, // must be ignored
        };

        double age = EpigeneticsAnalyzer.CalculatePhenoAge(methylation, coefficients, intercept);

        Assert.That(age, Is.EqualTo(11.0).Within(1e-12),
            "age = 10 + 2.0*0.5 = 11.0 returned untransformed; non-clock CpG ignored (Levine 2018 formula)");
    }

    // PA8 — built-in PhenoAge delegates to the caller-supplied overload with the built-in table/intercept.
    [Test]
    public void CalculatePhenoAge_BuiltInClock_EqualsExplicitOverloadWithSameTable()
    {
        var methylation = new Dictionary<string, double>
        {
            ["cg15611364"] = 0.42,
            ["cg17605084"] = 0.61,
            ["cg00503840"] = 0.33,
        };

        double viaDefault = EpigeneticsAnalyzer.CalculatePhenoAge(methylation);
        double viaExplicit = EpigeneticsAnalyzer.CalculatePhenoAge(
            methylation,
            EpigeneticsAnalyzer.PhenoAgeCoefficients,
            EpigeneticsAnalyzer.PhenoAgeIntercept);

        Assert.That(viaDefault, Is.EqualTo(viaExplicit).Within(1e-12),
            "Parameterless PhenoAge overload must delegate to the explicit one with the built-in table/intercept");
    }

    // PA9 — null / empty contract for the caller-supplied PhenoAge overload.
    [Test]
    public void CalculatePhenoAge_NullOrEmpty_Throws()
    {
        var methylation = new Dictionary<string, double> { ["cg1"] = 0.5 };
        var coefficients = new Dictionary<string, double> { ["cg1"] = 1.0 };
        Assert.Multiple(() =>
        {
            Assert.That(() => EpigeneticsAnalyzer.CalculatePhenoAge(null!),
                NUnit.Framework.Throws.ArgumentNullException,
                "Null methylation map is invalid (built-in PhenoAge contract)");
            Assert.That(() => EpigeneticsAnalyzer.CalculatePhenoAge(methylation, null!),
                NUnit.Framework.Throws.ArgumentNullException,
                "Null coefficient table is invalid (PhenoAge contract)");
            Assert.That(() => EpigeneticsAnalyzer.CalculatePhenoAge(methylation, new Dictionary<string, double>()),
                NUnit.Framework.Throws.ArgumentException,
                "Empty coefficient table has no defined PhenoAge output (contract)");
        });
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
