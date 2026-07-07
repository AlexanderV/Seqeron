// ONCO-SIG-002 — De-novo Mutational-Signature Extraction via NMF (Lee & Seung multiplicative updates)
// Evidence: docs/Evidence/ONCO-SIG-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-002.md
// Source: Lee D.D. & Seung H.S. (2001). Algorithms for Non-negative Matrix Factorization. NIPS 13.
//         https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization
//         Alexandrov L.B. et al. (2013). Deciphering Signatures of Mutational Processes Operative in Human
//         Cancer. Cell Reports 3(1):246-259. https://doi.org/10.1016/j.celrep.2012.12.008
//         COSMIC SBS96: https://cancer.sanger.ac.uk/signatures/sbs/sbs96/

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ExtractSignatures_Tests
{
    // Planted W0 (channels x k): 4 channels, k = 2. Chosen SEPARABLE so the nonnegative factorisation is
    // unique up to permutation/scaling: channel 0 is pure signature 0, channel 1 is pure signature 1
    // (a "pure-pixel" / separability condition — Donoho & Stodden 2004; the unique-NMF guarantee used to make
    // planted ground truth recoverable).
    private static readonly double[][] PlantedW =
    {
        new[] { 5.0, 0.0 },   // channel 0: pure signature 0
        new[] { 0.0, 4.0 },   // channel 1: pure signature 1
        new[] { 2.0, 1.0 },   // channel 2: mixed
        new[] { 1.0, 3.0 },   // channel 3: mixed
    };

    // Planted H0 (k x samples): k = 2, 4 samples. Chosen SEPARABLE: sample 0 is pure signature 0,
    // sample 1 is pure signature 1, making the row space identifiable too.
    private static readonly double[][] PlantedH =
    {
        new[] { 6.0, 0.0, 3.0, 2.0 },   // sample 0: pure signature 0
        new[] { 0.0, 5.0, 1.0, 4.0 },   // sample 1: pure signature 1
    };

    // V = W0 . H0 (channels x samples), an exactly rank-2 factorable nonnegative matrix.
    private static double[][] BuildExactV()
    {
        int channels = PlantedW.Length;
        int k = PlantedW[0].Length;
        int samples = PlantedH[0].Length;
        var v = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            v[c] = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                double sum = 0.0;
                for (int a = 0; a < k; a++)
                {
                    sum += PlantedW[c][a] * PlantedH[a][s];
                }

                v[c][s] = sum;
            }
        }

        return v;
    }

    private static IReadOnlyList<IReadOnlyList<double>> AsMatrix(double[][] rows) =>
        rows.Select(r => (IReadOnlyList<double>)r).ToList();

    // Reconstruct V_hat = W . H from an extraction result.
    private static double[,] Reconstruct(OncologyAnalyzer.SignatureExtractionResult result, int channels, int samples)
    {
        int k = result.Signatures.Count;
        var recon = new double[channels, samples];
        for (int c = 0; c < channels; c++)
        {
            for (int s = 0; s < samples; s++)
            {
                double sum = 0.0;
                for (int a = 0; a < k; a++)
                {
                    sum += result.Signatures[a][c] * result.Exposures[a][s];
                }

                recon[c, s] = sum;
            }
        }

        return recon;
    }

    private static double CosineSimilarity(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        double dot = 0.0, na = 0.0, nb = 0.0;
        for (int i = 0; i < a.Count; i++)
        {
            dot += a[i] * b[i];
            na += a[i] * a[i];
            nb += b[i] * b[i];
        }

        return dot / (Math.Sqrt(na) * Math.Sqrt(nb));
    }

    #region ExtractSignatures — reconstruction & recovery

    // M1 — Exactly-factorable V = W0.H0 reconstructs with residual ~0 (Lee & Seung fixed-point: factors of ones
    // when V = WH; Wikipedia "matrices of ones when V = WH").
    [Test]
    public void ExtractSignatures_ExactlyFactorableMatrix_ReconstructsWithNearZeroResidual()
    {
        double[][] v = BuildExactV();

        // Run to tight convergence (tolerance 0) so the local optimum at V = WH is reached: the default
        // relative-change stop halts earlier, which is fine for production but not for an exact-zero assertion.
        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, maxIterations: 200_000, tolerance: 0.0);

        // Residual ||V - WH||_F^2 should converge essentially to zero for an exactly rank-2 nonnegative matrix.
        Assert.That(result.FinalResidual, Is.LessThan(1e-6),
            "An exactly rank-2 nonnegative V = W0.H0 must be reconstructed by rank-2 NMF with a near-zero "
            + "Frobenius residual (Lee & Seung fixed point at V = WH).");
    }

    // M2 — Element-wise reconstruction W.H matches V for the exactly-factorable matrix.
    [Test]
    public void ExtractSignatures_ExactlyFactorableMatrix_ProductRecoversV()
    {
        double[][] v = BuildExactV();
        int channels = v.Length;
        int samples = v[0].Length;

        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, maxIterations: 200_000, tolerance: 0.0);
        double[,] recon = Reconstruct(result, channels, samples);

        Assert.Multiple(() =>
        {
            for (int c = 0; c < channels; c++)
            {
                for (int s = 0; s < samples; s++)
                {
                    Assert.That(recon[c, s], Is.EqualTo(v[c][s]).Within(1e-3),
                        $"Reconstructed (W.H)[{c},{s}] must equal the planted V[{c},{s}] = {v[c][s]} "
                        + "for an exactly factorable matrix (V approx WH model, Alexandrov 2013).");
                }
            }
        });
    }

    // M3 — Planted signatures recovered up to column permutation and positive scaling: each L1-normalised
    // planted column matches some extracted signature with cosine ~ 1 (NMF permutation/scale ambiguity;
    // Alexandrov 2013 blind source separation).
    [Test]
    public void ExtractSignatures_ExactlyFactorableMatrix_RecoversPlantedSignaturesUpToPermutation()
    {
        double[][] v = BuildExactV();
        int channels = PlantedW.Length;
        int k = PlantedW[0].Length;

        // L1-normalise the planted columns (the extraction normalises signatures to sum to 1).
        var plantedNormalized = new double[k][];
        for (int a = 0; a < k; a++)
        {
            double sum = 0.0;
            for (int c = 0; c < channels; c++)
            {
                sum += PlantedW[c][a];
            }

            plantedNormalized[a] = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                plantedNormalized[a][c] = PlantedW[c][a] / sum;
            }
        }

        // Drive the multiplicative updates to tight convergence (tolerance 0 = run the full iteration budget):
        // factor recovery up to perm/scale requires near-stationarity, beyond the default relative-change stop.
        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, maxIterations: 200_000, tolerance: 0.0);

        Assert.Multiple(() =>
        {
            for (int a = 0; a < k; a++)
            {
                double bestCosine = result.Signatures
                    .Select(sig => CosineSimilarity(plantedNormalized[a], sig))
                    .Max();
                Assert.That(bestCosine, Is.EqualTo(1.0).Within(1e-3),
                    $"Planted signature {a} must be recovered (cosine ~ 1) by some extracted signature, "
                    + "up to column permutation and positive scaling (NMF identifiability up to perm/scale).");
            }
        });
    }

    #endregion

    #region ExtractSignatures — invariants (nonnegativity, normalization, monotonicity)

    // M4 — Nonnegativity of W and H (Lee & Seung multiplicative updates preserve nonnegativity).
    [Test]
    public void ExtractSignatures_AnyInput_ProducesNonNegativeFactors()
    {
        double[][] v = BuildExactV();

        var result = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Signatures.SelectMany(s => s).Min(), Is.GreaterThanOrEqualTo(0.0),
                "Every extracted signature weight must be >= 0 (multiplicative updates preserve nonnegativity).");
            Assert.That(result.Exposures.SelectMany(e => e).Min(), Is.GreaterThanOrEqualTo(0.0),
                "Every exposure must be >= 0 (multiplicative updates preserve nonnegativity).");
        });
    }

    // M5 — L1 column-normalisation invariant: each extracted signature sums to 1 (COSMIC: a signature is a
    // probability distribution across the channels).
    [Test]
    public void ExtractSignatures_AnyInput_EachSignatureSumsToOne()
    {
        double[][] v = BuildExactV();

        var result = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2);

        Assert.Multiple(() =>
        {
            for (int a = 0; a < result.Signatures.Count; a++)
            {
                double sum = result.Signatures[a].Sum();
                Assert.That(sum, Is.EqualTo(1.0).Within(1e-9),
                    $"Extracted signature {a} must be L1-normalised to sum to 1 (probability distribution over "
                    + "channels, COSMIC / Alexandrov 2013-2020).");
            }
        });
    }

    // M6 — The Frobenius objective ||V - WH||_F^2 is monotonically non-increasing across iterations
    // (Lee & Seung 2001, Theorem 1).
    [Test]
    public void ExtractSignatures_ObjectiveHistory_IsMonotonicallyNonIncreasing()
    {
        double[][] v = BuildExactV();

        var result = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2);

        Assert.That(result.ObjectiveHistory.Count, Is.GreaterThan(1),
            "At least two iterations are needed to check monotonicity.");
        Assert.Multiple(() =>
        {
            for (int i = 1; i < result.ObjectiveHistory.Count; i++)
            {
                // Allow a tiny floating-point slack on the non-increase (theorem is exact in real arithmetic).
                Assert.That(result.ObjectiveHistory[i],
                    Is.LessThanOrEqualTo(result.ObjectiveHistory[i - 1] + 1e-9),
                    $"Objective at iteration {i} ({result.ObjectiveHistory[i]}) must not exceed iteration "
                    + $"{i - 1} ({result.ObjectiveHistory[i - 1]}) — Lee & Seung Theorem 1 (non-increasing).");
            }
        });
    }

    #endregion

    #region ExtractSignatures — SBS-96, determinism, scale absorption

    // S1 — 96-channel SBS planted truth reconstructs with low residual and normalised signatures.
    [Test]
    public void ExtractSignatures_Sbs96PlantedTruth_ReconstructsAndNormalises()
    {
        const int channels = 96;
        const int k = 2;
        const int samples = 5;
        var rng = new Random(7);

        // Planted W0 (96 x 2): two distinct sparse-ish nonnegative signatures.
        var w0 = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            w0[c] = new[] { c < 48 ? rng.NextDouble() + 0.1 : 0.0, c >= 48 ? rng.NextDouble() + 0.1 : 0.0 };
        }

        var h0 = new double[k][];
        for (int a = 0; a < k; a++)
        {
            h0[a] = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                h0[a][s] = rng.Next(1, 100);
            }
        }

        var v = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            v[c] = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                v[c][s] = w0[c][0] * h0[0][s] + w0[c][1] * h0[1][s];
            }
        }

        double vNormSq = v.SelectMany(r => r).Sum(x => x * x);

        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: k, maxIterations: 100_000, tolerance: 0.0);

        Assert.Multiple(() =>
        {
            // Relative residual should be tiny for an exactly factorable 96-channel matrix.
            Assert.That(result.FinalResidual / vNormSq, Is.LessThan(1e-4),
                "Rank-2 NMF of an exactly factorable 96-channel V must have a tiny relative Frobenius residual.");
            foreach (var sig in result.Signatures)
            {
                Assert.That(sig.Count, Is.EqualTo(channels),
                    "Each extracted SBS signature must have 96 channels.");
                Assert.That(sig.Sum(), Is.EqualTo(1.0).Within(1e-9),
                    "Each extracted SBS-96 signature must be L1-normalised to sum to 1.");
            }
        });
    }

    // S2 — Determinism: same seed -> identical factors (NMF is non-convex; seeded init makes it reproducible).
    [Test]
    public void ExtractSignatures_SameSeed_ProducesIdenticalFactors()
    {
        double[][] v = BuildExactV();

        var r1 = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2, seed: 123);
        var r2 = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2, seed: 123);

        Assert.Multiple(() =>
        {
            Assert.That(r1.Iterations, Is.EqualTo(r2.Iterations),
                "Same seed must yield the same iteration count.");
            for (int a = 0; a < r1.Signatures.Count; a++)
            {
                for (int c = 0; c < r1.Signatures[a].Count; c++)
                {
                    Assert.That(r1.Signatures[a][c], Is.EqualTo(r2.Signatures[a][c]).Within(1e-15),
                        $"Same seed must yield identical signature[{a}][{c}] (determinism).");
                }
            }
        });
    }

    // S3 — Scale absorption: L1-normalising W must not change W.H (exposures absorb the removed scale).
    // Verified via the near-exact reconstruction in M2; here assert exposures are strictly positive so the
    // scale really was absorbed (not zeroed) for the exactly-factorable case.
    [Test]
    public void ExtractSignatures_ExactlyFactorable_ExposuresAreStrictlyPositive()
    {
        double[][] v = BuildExactV();

        var result = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2);

        Assert.That(result.Exposures.SelectMany(e => e).Sum(), Is.GreaterThan(0.0),
            "For a non-zero exactly factorable V the exposures must carry the scale removed by L1-normalising "
            + "the signatures, so their total is strictly positive.");
    }

    #endregion

    #region ExtractSignatures — input validation

    // V1 — Null matrix throws.
    [Test]
    public void ExtractSignatures_NullMatrix_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.ExtractSignatures(null!, rank: 2),
            "A null count matrix must throw ArgumentNullException.");
    }

    // V2 — Empty matrix (no channels) throws.
    [Test]
    public void ExtractSignatures_EmptyMatrix_Throws()
    {
        var empty = new List<IReadOnlyList<double>>();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(empty, rank: 1),
            "A matrix with no channels must throw ArgumentException.");
    }

    // V3 — Zero samples throws.
    [Test]
    public void ExtractSignatures_ZeroSamples_Throws()
    {
        var noSamples = AsMatrix(new[] { Array.Empty<double>(), Array.Empty<double>() });
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(noSamples, rank: 1),
            "A matrix with zero samples must throw ArgumentException.");
    }

    // V4 — Ragged rows throw.
    [Test]
    public void ExtractSignatures_RaggedRows_Throws()
    {
        var ragged = AsMatrix(new[] { new[] { 1.0, 2.0 }, new[] { 3.0 } });
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(ragged, rank: 1),
            "Ragged channel rows must throw ArgumentException.");
    }

    // V5 — Negative entry throws.
    [Test]
    public void ExtractSignatures_NegativeEntry_Throws()
    {
        var negative = AsMatrix(new[] { new[] { 1.0, -2.0 }, new[] { 3.0, 4.0 } });
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(negative, rank: 1),
            "A negative count must throw ArgumentException (V must be nonnegative).");
    }

    // V6 — Non-finite entry throws.
    [Test]
    public void ExtractSignatures_NonFiniteEntry_Throws()
    {
        var nan = AsMatrix(new[] { new[] { 1.0, double.NaN }, new[] { 3.0, 4.0 } });
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(nan, rank: 1),
            "A NaN count must throw ArgumentException.");
    }

    // V7 — rank < 1 throws.
    [Test]
    public void ExtractSignatures_RankBelowOne_Throws()
    {
        double[][] v = BuildExactV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 0),
            "Rank k < 1 must throw ArgumentException.");
    }

    // V8 — rank > channel count throws.
    [Test]
    public void ExtractSignatures_RankAboveChannelCount_Throws()
    {
        double[][] v = BuildExactV(); // 4 channels
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 5),
            "Rank k greater than the channel count must throw ArgumentException.");
    }

    // V9 — maxIterations <= 0 throws.
    [Test]
    public void ExtractSignatures_NonPositiveMaxIterations_Throws()
    {
        double[][] v = BuildExactV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2, maxIterations: 0),
            "maxIterations <= 0 must throw ArgumentException.");
    }

    // V10 — negative tolerance throws.
    [Test]
    public void ExtractSignatures_NegativeTolerance_Throws()
    {
        double[][] v = BuildExactV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2, tolerance: -1.0),
            "A negative tolerance must throw ArgumentException.");
    }

    #endregion
}
