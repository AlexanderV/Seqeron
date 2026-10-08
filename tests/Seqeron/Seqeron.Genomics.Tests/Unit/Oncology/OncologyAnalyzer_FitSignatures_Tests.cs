// ONCO-SIG-002 — Mutational Signature Fitting / Refitting (NNLS + cosine similarity)
// Evidence: docs/Evidence/ONCO-SIG-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-002.md
// Source: Blokzijl F. et al. (2018). MutationalPatterns. Genome Medicine 10:33. https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/
//         Rosenthal R. et al. (2016). deconstructSigs. Genome Biology 17:31. https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
//         Lawson C.L. & Hanson R.J. (1974). Solving Least Squares Problems, Ch.23 (NNLS). https://en.wikipedia.org/wiki/Non-negative_least_squares
//         Pan W., Wang X. (2020). iMutSig. https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/

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

    #region Reference cross-checks (scipy.optimize.nnls / Lawson-Hanson NNLS) — campaign 2026-09

    // scipy.optimize.nnls docstring worked example (SciPy 1.17, _nnls.py "Examples"):
    // A=[[1,0],[1,0],[0,1]], b=[2,1,1] -> x=[1.5, 1.0], ‖Ax−b‖ = 0.7071067811865476.
    // Columns of A are the signatures: s0=[1,1,0], s1=[0,0,1].
    [Test]
    public void FitSignatures_ScipyDocstringExample_MatchesReferenceSolution()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 1, 0 }, new double[] { 0, 0, 1 } };

        var fit = OncologyAnalyzer.FitSignatures(new double[] { 2, 1, 1 }, signatures);

        double rnorm = Math.Sqrt(
            Math.Pow(fit.Reconstruction[0] - 2, 2) + Math.Pow(fit.Reconstruction[1] - 1, 2) +
            Math.Pow(fit.Reconstruction[2] - 1, 2));
        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(1.5).Within(1e-12), "scipy.optimize.nnls example: x0 = 1.5.");
            Assert.That(fit.Exposures[1], Is.EqualTo(1.0).Within(1e-12), "scipy.optimize.nnls example: x1 = 1.0.");
            Assert.That(rnorm, Is.EqualTo(0.7071067811865476).Within(1e-12), "scipy rnorm = 1/√2.");
        });
    }

    // scipy.optimize.nnls docstring worked example: b=[-1,-1,-1] -> x=[0,0], rnorm=√3 (every dual coefficient
    // w = Aᵀb is negative, so the Kuhn-Tucker conditions hold at x = 0).
    [Test]
    public void FitSignatures_ScipyDocstringNegativeTarget_ReturnsZero()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 1, 0 }, new double[] { 0, 0, 1 } };

        var fit = OncologyAnalyzer.FitSignatures(new double[] { -1, -1, -1 }, signatures);

        Assert.That(fit.Exposures, Is.EqualTo(new[] { 0.0, 0.0 }), "scipy.optimize.nnls example: x = [0, 0].");
    }

    // Regression (campaign 2026-09, F1): the former solver stopped when max(w_R) ≤ 1e-12 (absolute), so a catalog
    // whose dual coefficients are below 1e-12 returned the all-zero fit. Lawson-Hanson NNLS uses the sign test
    // w_j > 0 (no absolute threshold); scipy.optimize.nnls(I, [3e-13, 5e-13]) = [3e-13, 5e-13], rnorm 0.
    [Test]
    public void FitSignatures_TinyMagnitudeCatalog_IsNotTruncatedToZero()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };

        var fit = OncologyAnalyzer.FitSignatures(new double[] { 3e-13, 5e-13 }, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(3e-13).Within(1e-27), "scipy: x0 = 3e-13.");
            Assert.That(fit.Exposures[1], Is.EqualTo(5e-13).Within(1e-27), "scipy: x1 = 5e-13.");
            Assert.That(fit.NormalizedExposures[0], Is.EqualTo(0.375).Within(1e-12), "Proportions are scale-free.");
        });
    }

    // NNLS is positively homogeneous: x*(c·d) = c·x*(d) for c > 0 (the feasible cone and the objective scale
    // together). Constraint-binding case M6: S=[[1,1],[0,1]] (columns), d=[0,1] -> x=[0,0.5] (scipy).
    [TestCase(1e-12)]
    [TestCase(1e-6)]
    [TestCase(1.0)]
    [TestCase(1e6)]
    [TestCase(1e12)]
    public void FitSignatures_ScaledCatalog_ExposuresScaleLinearly(double scale)
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 1, 1 } };

        var fit = OncologyAnalyzer.FitSignatures(new[] { 0.0, 1.0 * scale }, signatures);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Exposures[0], Is.EqualTo(0.0), "The binding coefficient stays exactly 0 at every scale.");
            Assert.That(fit.Exposures[1], Is.EqualTo(0.5 * scale).Within(1e-12 * scale), "x1 = 0.5·scale.");
        });
    }

    // Non-finite input: scipy.optimize.nnls rejects NaN/Inf (np.asarray_chkfinite -> ValueError). Previously a NaN
    // catalog silently produced the all-zero fit (every NaN comparison is false).
    [Test]
    public void FitSignatures_NonFiniteInput_Throws()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        var badSignatures = new IReadOnlyList<double>[] { new double[] { 1, double.PositiveInfinity }, new double[] { 0, 1 } };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.FitSignatures(new[] { double.NaN, 1.0 }, signatures), "NaN catalog value.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.FitSignatures(new[] { 1.0, double.PositiveInfinity }, signatures), "Infinite catalog value.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.FitSignatures(new[] { 1.0, 1.0 }, badSignatures), "Infinite signature value.");
        });
    }

    // Duplicate (collinear) signature columns: Lawson-Hanson rejects a candidate column that is numerically dependent
    // on the passive set, so the duplicate keeps exposure 0 and the fit is unique in its first occurrence.
    // scipy.optimize.nnls([[1,1,0],[0,0,1]]ᵀ-with-duplicate, [2,1,1]) = [1.5, 0, 1].
    [Test]
    public void FitSignatures_DuplicateSignature_AssignsExposureToFirstCopyOnly()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 1, 0 }, new double[] { 1, 1, 0 }, new double[] { 0, 0, 1 }
        };

        var fit = OncologyAnalyzer.FitSignatures(new double[] { 2, 1, 1 }, signatures);

        Assert.That(fit.Exposures.ToArray(), Is.EqualTo(new[] { 1.5, 0.0, 1.0 }).Within(1e-12),
            "scipy.optimize.nnls returns [1.5, 0, 1] (duplicate column rejected as linearly dependent).");
    }

    private const string CosmicResource = "Seqeron.Genomics.Tests.TestData.Cosmic.COSMIC_v3.4_SBS_GRCh37.txt";

    // Sample catalog (SBS96, COSMIC file row order): Poisson draw of 300·SBS1 + 800·SBS5 + 200·SBS32
    // (numpy default_rng(7)); total 1345 mutations.
    private static readonly double[] CosmicSampleCatalog =
    {
        18, 10, 0, 13, 4, 6, 0, 5, 46, 41, 126, 39, 6, 8, 7, 7, 38, 12, 29, 32, 2, 1, 9, 5, 5, 7, 4, 12, 10, 5, 4, 9,
        20, 10, 68, 28, 3, 6, 5, 2, 12, 10, 19, 10, 3, 5, 4, 6, 8, 4, 3, 9, 3, 4, 0, 12, 21, 24, 88, 24, 6, 5, 3, 5,
        17, 6, 21, 12, 1, 1, 5, 2, 7, 11, 2, 15, 8, 6, 2, 13, 32, 25, 55, 33, 4, 6, 3, 6, 20, 16, 11, 21, 4, 8, 6, 16
    };

    // scipy.optimize.nnls(COSMIC_v3.4 (96×86), CosmicSampleCatalog): the 38 non-zero exposures (all others 0),
    // ‖Sx − d‖ = 21.901597834456467, cos(d, Sx) = 0.9953491259911486.
    private static readonly (int Index, double Exposure)[] CosmicScipyExposures =
    {
        (0, 304.0688771932464), (1, 23.73211791015279), (2, 53.86002360763532), (4, 311.23237916325286),
        (6, 1.2438724289510348), (12, 0.9025016760233332), (13, 23.513232263685637), (17, 88.67694064659317),
        (18, 7.968045767608636), (19, 17.716670821732798), (20, 40.218825827644956), (21, 37.976441064275285),
        (23, 2.1268139033159748), (25, 2.1731210882910674), (27, 6.138432736752313), (35, 4.146875469314331),
        (37, 14.375339604630556), (39, 191.62931731138946), (40, 1.2447328629706593), (46, 22.163770905276937),
        (48, 24.11656401977636), (53, 20.35475176897135), (57, 0.2353436147162052), (58, 0.22641770305395006),
        (59, 17.71897831104093), (60, 9.29127287421852), (63, 30.429149364272057), (64, 3.0423928376482268),
        (66, 11.854094560342078), (68, 2.484737117036816), (70, 2.05630013505543), (74, 20.782638412503584),
        (75, 9.756729641081455), (79, 1.0553900904886457), (82, 4.343798767503239), (83, 13.630725698197786),
        (84, 3.8763612238757204), (85, 14.554738516180588)
    };

    private static IReadOnlyList<IReadOnlyList<double>> LoadCosmicSignatures()
    {
        var asm = typeof(OncologyAnalyzer_FitSignatures_Tests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(CosmicResource)
            ?? throw new InvalidOperationException($"Embedded resource '{CosmicResource}' not found.");
        using var reader = new StreamReader(stream);
        string[] header = reader.ReadLine()!.Split('\t');
        int signatureCount = header.Length - 1;
        var columns = new List<double>[signatureCount];
        for (int j = 0; j < signatureCount; j++)
        {
            columns[j] = new List<double>(96);
        }

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (line.Length == 0)
            {
                continue;
            }

            string[] fields = line.Split('\t');
            for (int j = 0; j < signatureCount; j++)
            {
                columns[j].Add(double.Parse(fields[j + 1], System.Globalization.CultureInfo.InvariantCulture));
            }
        }

        return columns.Select(c => (IReadOnlyList<double>)c.ToArray()).ToArray();
    }

    // Refit to the full COSMIC v3.4 SBS set (MutationalPatterns fit_to_signatures usage: pracma::lsqnonneg over all
    // reference signatures). Every exposure matches scipy.optimize.nnls (Lawson-Hanson) to 1e-9 relative.
    // scale = 1e6 is the regression for F1: the former solver did not terminate within 90 s on this input
    // (absolute 1e-12 dual tolerance + no ztest/independence safeguards -> add/drop cycling).
    [TestCase(1.0)]
    [TestCase(1e6)]
    public void FitSignatures_FullCosmicV34Refit_MatchesScipyNnls(double scale)
    {
        IReadOnlyList<IReadOnlyList<double>> cosmic = LoadCosmicSignatures();
        double[] catalog = CosmicSampleCatalog.Select(v => v * scale).ToArray();

        var task = Task.Run(() => OncologyAnalyzer.FitSignatures(catalog, cosmic));
        Assert.That(task.Wait(TimeSpan.FromSeconds(60)), Is.True, "NNLS must terminate (Lawson-Hanson is finite).");
        var fit = task.Result;

        var expected = new double[cosmic.Count];
        foreach ((int index, double exposure) in CosmicScipyExposures)
        {
            expected[index] = exposure * scale;
        }

        double rnorm = Math.Sqrt(catalog.Select((d, k) => Math.Pow(d - fit.Reconstruction[k], 2)).Sum());
        Assert.Multiple(() =>
        {
            Assert.That(cosmic.Count, Is.EqualTo(86), "COSMIC v3.4 SBS has 86 signatures.");
            for (int j = 0; j < cosmic.Count; j++)
            {
                Assert.That(fit.Exposures[j], Is.EqualTo(expected[j]).Within(1e-9 * scale * 311.23),
                    $"Exposure of COSMIC column {j} must equal scipy.optimize.nnls.");
            }

            Assert.That(rnorm, Is.EqualTo(21.901597834456467 * scale).Within(1e-9 * scale), "scipy rnorm.");
            Assert.That(fit.ReconstructionCosineSimilarity, Is.EqualTo(0.9953491259911486).Within(1e-12),
                "Reconstruction cosine (MutationalPatterns cos_sim) of the scipy solution.");
        });
    }

    #endregion
}
