// ONCO-SIG-002 — Mutational Signature Fitting / Refitting (NNLS + cosine similarity)
// Evidence: docs/Evidence/ONCO-SIG-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-002.md
// Source: Blokzijl F. et al. (2018). MutationalPatterns. Genome Medicine 10:33. https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/
//         Rosenthal R. et al. (2016). deconstructSigs. Genome Biology 17:31. https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
//         Lawson C.L. & Hanson R.J. (1974). Solving Least Squares Problems, Ch.23 (NNLS). https://en.wikipedia.org/wiki/Non-negative_least_squares
//         Pan W., Wang X. (2020). iMutSig. https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_FitSignatures_Tests
{
    // 1/sqrt(2): cosine of [1,1] vs [1,0] = 1 / (sqrt(2) * 1). Hand-derived from the formula.
    private const double OneOverSqrt2 = 0.70710678118654752440;

    #region CosineSimilarity

    // M1 — Identical vectors: sim = 14/(sqrt(14)*sqrt(14)) = 1 (Blokzijl 2018; iMutSig).
    [Test]
    public void CosineSimilarity_IdenticalVectors_ReturnsOne()
    {
        var a = new double[] { 1, 2, 3 };
        var b = new double[] { 1, 2, 3 };

        double sim = OncologyAnalyzer.CosineSimilarity(a, b);

        Assert.That(sim, Is.EqualTo(1.0).Within(1e-10),
            "Identical non-zero vectors are parallel, so cosine similarity must be exactly 1 (sim=14/(√14·√14)).");
    }

    // M2 — Orthogonal vectors: dot product 0 -> sim 0 (Blokzijl 2018).
    [Test]
    public void CosineSimilarity_OrthogonalVectors_ReturnsZero()
    {
        var a = new double[] { 1, 0 };
        var b = new double[] { 0, 1 };

        double sim = OncologyAnalyzer.CosineSimilarity(a, b);

        Assert.That(sim, Is.EqualTo(0.0).Within(1e-10),
            "Disjoint-support vectors have dot product 0, so cosine similarity must be exactly 0.");
    }

    // M3 — General case: [1,1] vs [1,0] = 1/(√2·1) = 1/√2 (formula, hand-derived).
    [Test]
    public void CosineSimilarity_GeneralVectors_ReturnsExactRatio()
    {
        var a = new double[] { 1, 1 };
        var b = new double[] { 1, 0 };

        double sim = OncologyAnalyzer.CosineSimilarity(a, b);

        Assert.That(sim, Is.EqualTo(OneOverSqrt2).Within(1e-10),
            "cos([1,1],[1,0]) = 1/(√2·1) = 1/√2 ≈ 0.70710678 by the dot-product/norm formula.");
    }

    // M4 — Scale invariance: [3,4] vs [6,8] = 50/(5·10) = 1 (cosine of angle, iMutSig).
    [Test]
    public void CosineSimilarity_PositivelyScaledVector_IsScaleInvariant()
    {
        var a = new double[] { 3, 4 };
        var b = new double[] { 6, 8 };

        double sim = OncologyAnalyzer.CosineSimilarity(a, b);

        Assert.That(sim, Is.EqualTo(1.0).Within(1e-10),
            "b is 2·a (same direction); cosine is scale-invariant so sim = 50/(5·10) = 1.");
    }

    // S2 — Zero-norm vector: cosine undefined (÷0) -> documented 0.0 (Assumption 2).
    [Test]
    public void CosineSimilarity_ZeroVector_ReturnsZero()
    {
        var a = new double[] { 0, 0 };
        var b = new double[] { 1, 1 };

        double sim = OncologyAnalyzer.CosineSimilarity(a, b);

        Assert.That(sim, Is.EqualTo(0.0).Within(1e-12),
            "A zero-norm vector has no direction; the documented degenerate result is 0.0.");
    }

    // C1 — Property: cos(a, k·a) = 1 for any k>0 (scale invariance, INV-03).
    [Test]
    public void CosineSimilarity_VectorAgainstPositiveMultipleOfItself_IsAlwaysOne()
    {
        var a = new double[] { 2, 5, 1, 7 };

        Assert.Multiple(() =>
        {
            foreach (double k in new[] { 0.5, 1.0, 3.0, 100.0 })
            {
                var scaled = a.Select(v => v * k).ToArray();
                Assert.That(OncologyAnalyzer.CosineSimilarity(a, scaled), Is.EqualTo(1.0).Within(1e-10),
                    $"cos(a, {k}·a) must be 1 because positive scaling preserves direction.");
            }
        });
    }

    [Test]
    public void CosineSimilarity_NullInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.CosineSimilarity(null!, new double[] { 1 }),
                "Null first vector must throw ArgumentNullException.");
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.CosineSimilarity(new double[] { 1 }, null!),
                "Null second vector must throw ArgumentNullException.");
        });
    }

    [Test]
    public void CosineSimilarity_LengthMismatch_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CosineSimilarity(new double[] { 1, 2 }, new double[] { 1 }),
            "Vectors of different lengths cannot be compared and must throw ArgumentException.");
    }

    // Cosine similarity is "undefined for empty vectors" (no components to sum) -> ArgumentException.
    [Test]
    public void CosineSimilarity_EmptyVectors_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CosineSimilarity(Array.Empty<double>(), Array.Empty<double>()),
            "Cosine similarity has no components to sum for empty vectors and must throw ArgumentException.");
    }

    #endregion

    #region FitSignatures

    // M5 — Identity signature matrix recovers exposures exactly: x = d (unconstrained LS already ≥ 0).
    [Test]
    public void FitSignatures_IdentitySignatures_RecoversExposures()
    {
        var catalog = new double[] { 3, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 },
            new double[] { 0, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(3.0).Within(1e-9),
                "With orthonormal signatures S=I, the NNLS solution equals d, so exposure[0] = 3.");
            Assert.That(fit.Exposures[1], Is.EqualTo(5.0).Within(1e-9),
                "With S=I the NNLS solution equals d, so exposure[1] = 5.");
        });
    }

    // M6 / M10 — Constraint binds: S=[[1,1],[0,1]], d=[0,1]; unconstrained x1=-1 -> clamp to 0, refit -> [0,0.5].
    [Test]
    public void FitSignatures_NegativeUnconstrainedCoefficient_ClampsToZeroAndRefits()
    {
        var catalog = new double[] { 0, 1 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 }, // signature 1
            new double[] { 1, 1 }  // signature 2
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(0.0).Within(1e-9),
                "Signature 1's unconstrained weight is -1 (<0); NNLS clamps it to exactly 0.");
            Assert.That(fit.Exposures[1], Is.EqualTo(0.5).Within(1e-9),
                "Refitting signature 2 alone gives ([1,1]·[0,1])/([1,1]·[1,1]) = 1/2.");
            Assert.That(fit.Exposures.Min(), Is.GreaterThanOrEqualTo(0.0),
                "INV-04: all NNLS exposures must be non-negative.");
        });
    }

    // M7 / M8 — Reconstruction S·x and reconstruction cosine for an exactly representable catalog.
    [Test]
    public void FitSignatures_ExactlyRepresentableCatalog_ReconstructsWithCosineOne()
    {
        var catalog = new double[] { 3, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 },
            new double[] { 0, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Reconstruction[0], Is.EqualTo(3.0).Within(1e-9),
                "Reconstruction S·x must reproduce channel 0 = 3.");
            Assert.That(fit.Reconstruction[1], Is.EqualTo(5.0).Within(1e-9),
                "Reconstruction S·x must reproduce channel 1 = 5.");
            Assert.That(fit.ReconstructionCosineSimilarity, Is.EqualTo(1.0).Within(1e-10),
                "An exactly representable catalog reconstructs to itself, so cosine = 1 (Blokzijl 2018).");
        });
    }

    // M9 — Normalised exposures are proportions summing to 1: [3,5] -> [0.375, 0.625] (deconstructSigs).
    [Test]
    public void FitSignatures_NormalizedExposures_AreProportionsSummingToOne()
    {
        var catalog = new double[] { 3, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 },
            new double[] { 0, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.NormalizedExposures[0], Is.EqualTo(0.375).Within(1e-9),
                "3 / (3+5) = 0.375 (deconstructSigs normalises weights to proportions).");
            Assert.That(fit.NormalizedExposures[1], Is.EqualTo(0.625).Within(1e-9),
                "5 / (3+5) = 0.625.");
            Assert.That(fit.NormalizedExposures.Sum(), Is.EqualTo(1.0).Within(1e-10),
                "INV-06: normalised exposures must sum to 1 when the total is positive.");
        });
    }

    // S1 — Zero catalog: all exposures 0, reconstruction 0, proportions 0 (degenerate NNLS minimiser).
    [Test]
    public void FitSignatures_ZeroCatalog_YieldsAllZeroFit()
    {
        var catalog = new double[] { 0, 0 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 },
            new double[] { 0, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures, Is.All.EqualTo(0.0),
                "The only feasible minimiser of ‖S·x‖² with x≥0 for d=0 is x=0.");
            Assert.That(fit.Reconstruction, Is.All.EqualTo(0.0),
                "Reconstruction of an all-zero exposure vector is the zero catalog.");
            Assert.That(fit.NormalizedExposures, Is.All.EqualTo(0.0),
                "INV-06: with Σexposures = 0 the normalised exposures are all 0 (no division).");
        });
    }

    // S3 / INV-05 — Fit residual SSE never exceeds the SSE of the all-zero fit (‖d‖²).
    [Test]
    public void FitSignatures_ResidualSse_DoesNotExceedZeroFitSse()
    {
        var catalog = new double[] { 4, 1, 7 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 1, 0 },
            new double[] { 0, 1, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        double fitSse = 0.0;
        for (int k = 0; k < catalog.Length; k++)
        {
            double diff = catalog[k] - fit.Reconstruction[k];
            fitSse += diff * diff;
        }

        double zeroFitSse = catalog.Sum(v => v * v); // ‖d‖² is the SSE of the feasible x = 0.

        Assert.That(fitSse, Is.LessThanOrEqualTo(zeroFitSse + 1e-9),
            "INV-05: x=0 is feasible, so the NNLS minimiser's residual SSE must be ≤ ‖d‖².");
    }

    // Imperfect (under-determined) fit: a single flat signature [1,1,1] cannot represent d=[3,0,0].
    // NNLS exposure = (s·d)/(s·s) = 3/3 = 1; reconstruction = [1,1,1]; reconstruction cosine =
    // cos([3,0,0],[1,1,1]) = 3/(3·√3) = 1/√3 ≈ 0.5773502691896258 (Blokzijl 2018 reconstruction-quality
    // measure; below the 0.95 "successful reconstruction" threshold). Cross-checked vs scipy.optimize.nnls.
    [Test]
    public void FitSignatures_ImperfectFit_ReportsExactSubUnityReconstructionCosine()
    {
        var catalog = new double[] { 3, 0, 0 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 1, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(1.0).Within(1e-9),
                "NNLS exposure of the flat signature = (s·d)/(s·s) = 3/3 = 1.");
            Assert.That(fit.Reconstruction[0], Is.EqualTo(1.0).Within(1e-9), "S·x channel 0 = 1·1 = 1.");
            Assert.That(fit.Reconstruction[1], Is.EqualTo(1.0).Within(1e-9), "S·x channel 1 = 1·1 = 1.");
            Assert.That(fit.Reconstruction[2], Is.EqualTo(1.0).Within(1e-9), "S·x channel 2 = 1·1 = 1.");
            Assert.That(fit.ReconstructionCosineSimilarity, Is.EqualTo(0.57735026918962584).Within(1e-12),
                "cos([3,0,0],[1,1,1]) = 3/(3·√3) = 1/√3 (reconstruction quality below the 0.95 threshold).");
        });
    }

    [Test]
    public void FitSignatures_NullCatalog_Throws()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1 } };
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.FitSignatures(null!, signatures),
            "Null catalog must throw ArgumentNullException.");
    }

    [Test]
    public void FitSignatures_DimensionMismatch_Throws()
    {
        var catalog = new double[] { 1, 2, 3 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 } };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.FitSignatures(catalog, signatures),
            "A catalog length differing from the signature channel count must throw ArgumentException.");
    }

    [Test]
    public void FitSignatures_NoSignatures_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.FitSignatures(new double[] { 1 }, Array.Empty<IReadOnlyList<double>>()),
            "Fitting requires at least one reference signature.");
    }

    #endregion

    #region ReconstructCatalog

    // M7 (direct) — S·x with explicit exposures.
    [Test]
    public void ReconstructCatalog_TwoSignatures_ComputesMatrixVectorProduct()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0 },
            new double[] { 1, 1 }
        };
        var exposures = new double[] { 2, 3 };

        var reconstruction = OncologyAnalyzer.ReconstructCatalog(signatures, exposures);

        Assert.Multiple(() =>
        {
            Assert.That(reconstruction[0], Is.EqualTo(5.0).Within(1e-10),
                "Channel 0 = 1·2 + 1·3 = 5 (S·x).");
            Assert.That(reconstruction[1], Is.EqualTo(3.0).Within(1e-10),
                "Channel 1 = 0·2 + 1·3 = 3 (S·x).");
        });
    }

    [Test]
    public void ReconstructCatalog_ExposureCountMismatch_Throws()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 } };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ReconstructCatalog(signatures, new double[] { 1, 2 }),
            "Exposure count must equal the signature count; otherwise ArgumentException.");
    }

    [Test]
    public void ReconstructCatalog_NullArguments_Throw()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.ReconstructCatalog(null!, new double[] { 1 }),
                "Null signatures must throw ArgumentNullException.");
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.ReconstructCatalog(signatures, null!),
                "Null exposures must throw ArgumentNullException.");
        });
    }

    #endregion
}
