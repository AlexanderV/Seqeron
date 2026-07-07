// ONCO-SIG-002 — Automatic NMF rank selection + Poisson/KL objective + COSMIC cosine matching
// Evidence: docs/Evidence/ONCO-SIG-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-002.md
// Source: Lee D.D. & Seung H.S. (2001). Algorithms for Non-negative Matrix Factorization. NIPS 13.
//         https://arxiv.org/html/2501.11341v1 (verbatim Theorem 1/2 update rules)
//         Brunet J-P. et al. (2004). Metagenes and molecular pattern discovery using matrix factorization.
//         PNAS 101(12):4164-4169. https://doi.org/10.1073/pnas.0308531101 (consensus / cophenetic)
//         Alexandrov L.B. et al. (2013). Cell Reports 3(1):246-259. https://doi.org/10.1016/j.celrep.2012.12.008
//         Islam S.M.A. et al. (2022). SigProfilerExtractor. Cell Genomics 2(11):100179
//         https://doi.org/10.1016/j.xgen.2022.100179 (silhouette stability >= 0.80; cosine matching)

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_SelectRank_Tests
{
    // Planted truth: k0 = 2 known signatures over 4 channels (L1-normalised probability distributions),
    // sampled into 6 samples so each true signature dominates >= 1 sample (Alexandrov 2013 blind source
    // separation framing). V = W0 . H0 is exactly rank-2 factorable. The signatures are SEPARABLE (channel 0 is
    // pure signature 0, channel 2 is pure signature 1 — a "pure-channel" condition, Donoho & Stodden 2004) so
    // the nonnegative factorisation is identifiable up to permutation/scaling and the planted truth is exactly
    // recoverable.
    private static readonly double[] TrueSignature0 = { 0.70, 0.20, 0.00, 0.10 };
    private static readonly double[] TrueSignature1 = { 0.00, 0.10, 0.70, 0.20 };

    // H0 (k x samples): samples 0-2 dominated by signature 0, samples 3-5 by signature 1, with mild mixing.
    private static readonly double[][] PlantedH =
    {
        new[] { 100.0, 80.0, 60.0, 5.0, 10.0, 0.0 },   // signature 0 exposures
        new[] { 0.0, 10.0, 5.0, 70.0, 90.0, 120.0 },   // signature 1 exposures
    };

    private static double[][] BuildPlantedV()
    {
        const int channels = 4;
        int samples = PlantedH[0].Length;
        var v = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            v[c] = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                v[c][s] = TrueSignature0[c] * PlantedH[0][s] + TrueSignature1[c] * PlantedH[1][s];
            }
        }

        return v;
    }

    private static IReadOnlyList<IReadOnlyList<double>> AsMatrix(double[][] rows) =>
        rows.Select(r => (IReadOnlyList<double>)r).ToList();

    private static IReadOnlyList<IReadOnlyList<double>> AsSignatures(params double[][] sigs) =>
        sigs.Select(s => (IReadOnlyList<double>)s).ToList();

    private static double Cosine(IReadOnlyList<double> a, IReadOnlyList<double> b)
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

    #region KL / Poisson objective extraction (Lee & Seung Theorem 2)

    // M-KL1 — The KL divergence D(V||WH) is monotonically non-increasing across iterations (Lee & Seung 2001,
    // Theorem 2: "the divergence D(V||WH) is nonincreasing under the update rules").
    [Test]
    public void ExtractSignatures_KullbackLeibler_DivergenceHistoryIsMonotonicallyNonIncreasing()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.KullbackLeibler,
            maxIterations: 5_000, tolerance: 0.0, seed: 1);

        Assert.That(result.ObjectiveHistory.Count, Is.GreaterThan(1),
            "At least two iterations are needed to check monotonicity.");
        Assert.Multiple(() =>
        {
            for (int i = 1; i < result.ObjectiveHistory.Count; i++)
            {
                // Allow a tiny floating-point slack on the non-increase (theorem is exact in real arithmetic).
                Assert.That(result.ObjectiveHistory[i],
                    Is.LessThanOrEqualTo(result.ObjectiveHistory[i - 1] + 1e-7),
                    $"KL divergence at iteration {i} ({result.ObjectiveHistory[i]}) must not exceed iteration "
                    + $"{i - 1} ({result.ObjectiveHistory[i - 1]}) — Lee & Seung Theorem 2 (non-increasing).");
            }
        });
    }

    // M-KL2 — KL extraction of the planted k0=2 data recovers BOTH true signatures with cosine ~ 1
    // (Alexandrov 2013: NMF de novo signatures recover the latent mutational processes).
    [Test]
    public void ExtractSignatures_KullbackLeibler_RecoversPlantedSignatures()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.KullbackLeibler,
            maxIterations: 50_000, tolerance: 0.0, seed: 1);

        double bestForS0 = result.Signatures.Select(sig => Cosine(TrueSignature0, sig)).Max();
        double bestForS1 = result.Signatures.Select(sig => Cosine(TrueSignature1, sig)).Max();

        Assert.Multiple(() =>
        {
            Assert.That(bestForS0, Is.GreaterThan(0.999),
                "KL-NMF must recover planted signature 0 (cosine ~ 1) up to permutation (Alexandrov 2013).");
            Assert.That(bestForS1, Is.GreaterThan(0.999),
                "KL-NMF must recover planted signature 1 (cosine ~ 1) up to permutation (Alexandrov 2013).");
        });
    }

    // M-KL3 — KL extraction produces non-negative, L1-normalised signatures (COSMIC: probability distributions).
    [Test]
    public void ExtractSignatures_KullbackLeibler_ProducesNormalisedNonNegativeSignatures()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.KullbackLeibler, seed: 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.Signatures.SelectMany(s => s).Min(), Is.GreaterThanOrEqualTo(0.0),
                "KL multiplicative updates must preserve signature non-negativity.");
            Assert.That(result.Exposures.SelectMany(e => e).Min(), Is.GreaterThanOrEqualTo(0.0),
                "KL multiplicative updates must preserve exposure non-negativity.");
            foreach (var sig in result.Signatures)
            {
                Assert.That(sig.Sum(), Is.EqualTo(1.0).Within(1e-9),
                    "Each KL-extracted signature must be L1-normalised to sum to 1 (COSMIC convention).");
            }
        });
    }

    // M-KL4 — Determinism: same seed -> identical KL factors (NMF is non-convex; seeded init is reproducible).
    [Test]
    public void ExtractSignatures_KullbackLeibler_SameSeedProducesIdenticalFactors()
    {
        double[][] v = BuildPlantedV();

        var r1 = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.KullbackLeibler, seed: 99);
        var r2 = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.KullbackLeibler, seed: 99);

        Assert.Multiple(() =>
        {
            Assert.That(r1.Iterations, Is.EqualTo(r2.Iterations), "Same seed must yield the same iteration count.");
            for (int a = 0; a < r1.Signatures.Count; a++)
            {
                for (int c = 0; c < r1.Signatures[a].Count; c++)
                {
                    Assert.That(r1.Signatures[a][c], Is.EqualTo(r2.Signatures[a][c]).Within(1e-15),
                        $"Same seed must yield identical KL signature[{a}][{c}] (determinism).");
                }
            }
        });
    }

    // M-KL5 — The default ExtractSignatures overload still uses the Frobenius objective (backward-compat):
    // its objective value equals the Frobenius residual of the explicit-Frobenius overload.
    [Test]
    public void ExtractSignatures_DefaultOverload_MatchesExplicitFrobenius()
    {
        double[][] v = BuildPlantedV();

        var def = OncologyAnalyzer.ExtractSignatures(AsMatrix(v), rank: 2, seed: 5);
        var frob = OncologyAnalyzer.ExtractSignatures(
            AsMatrix(v), rank: 2, OncologyAnalyzer.NmfObjective.Frobenius, seed: 5);

        Assert.That(def.FinalResidual, Is.EqualTo(frob.FinalResidual).Within(1e-15),
            "The default ExtractSignatures overload must be identical to the explicit Frobenius objective "
            + "(backward compatibility — the original Frobenius path is unchanged).");
    }

    #endregion

    #region Automatic rank selection (Brunet 2004 cophenetic; Alexandrov 2013 / SigProfiler stability)

    // M-RS1 — On planted k0=2 data the automatic selection picks rank 2 (high stability + acceptable error),
    // and reports per-rank diagnostics (Brunet 2004; Alexandrov 2013 stability/error trade-off).
    [Test]
    public void SelectRank_PlantedTwoSignatureData_SelectsTrueRankTwo()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 1, maxRank: 3,
            objective: OncologyAnalyzer.NmfObjective.KullbackLeibler,
            runs: 10, maxIterations: 3_000, tolerance: 1e-8, seed: 7);

        Assert.Multiple(() =>
        {
            Assert.That(result.SelectedRank, Is.EqualTo(2),
                "Data synthesised from k0 = 2 known signatures must yield selected rank 2 "
                + "(SigProfiler stability/error trade-off; Brunet 2004).");
            Assert.That(result.PerRank.Count, Is.EqualTo(3),
                "Per-rank diagnostics must be reported for every candidate rank 1..3 (auditability).");
            Assert.That(result.PerRank.Select(p => p.Rank), Is.EqualTo(new[] { 1, 2, 3 }),
                "Per-rank diagnostics must be in ascending rank order.");
        });
    }

    // M-RS2 — Rank 1 has cophenetic correlation = 1.0 exactly: with k=1 every sample is in the single cluster,
    // so every connectivity matrix is all-ones and the consensus is all-ones (Brunet 2004 perfect consensus).
    [Test]
    public void SelectRank_RankOne_HasCopheneticCorrelationOne()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 1, maxRank: 2,
            objective: OncologyAnalyzer.NmfObjective.KullbackLeibler,
            runs: 5, maxIterations: 1_000, tolerance: 1e-6, seed: 7);

        var rank1 = result.PerRank.Single(p => p.Rank == 1);
        Assert.That(rank1.CopheneticCorrelation, Is.EqualTo(1.0).Within(1e-12),
            "At rank 1 the consensus matrix is all-ones (one cluster), so the cophenetic correlation is 1.0 "
            + "by construction (Brunet 2004 perfect consensus).");
    }

    // M-RS3 — Reconstruction error decreases (or stays flat) as rank increases on factorable data: a richer
    // model fits at least as well (Alexandrov 2013 reconstruction-error component).
    [Test]
    public void SelectRank_IncreasingRank_DoesNotIncreaseMeanReconstructionError()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 1, maxRank: 3,
            objective: OncologyAnalyzer.NmfObjective.Frobenius,
            runs: 5, maxIterations: 3_000, tolerance: 1e-9, seed: 7);

        double err1 = result.PerRank.Single(p => p.Rank == 1).MeanReconstructionError;
        double err2 = result.PerRank.Single(p => p.Rank == 2).MeanReconstructionError;

        Assert.That(err2, Is.LessThanOrEqualTo(err1 + 1e-6),
            "Rank 2 must fit the planted rank-2 data at least as well as rank 1 (lower/equal Frobenius error).");
    }

    // M-RS4 — Determinism: same seed -> identical selection report (fixed derived seed sequence).
    [Test]
    public void SelectRank_SameSeed_ProducesIdenticalReport()
    {
        double[][] v = BuildPlantedV();

        var r1 = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 1, maxRank: 3, runs: 5, maxIterations: 1_000, tolerance: 1e-6, seed: 13);
        var r2 = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 1, maxRank: 3, runs: 5, maxIterations: 1_000, tolerance: 1e-6, seed: 13);

        Assert.Multiple(() =>
        {
            Assert.That(r1.SelectedRank, Is.EqualTo(r2.SelectedRank), "Same seed must select the same rank.");
            for (int i = 0; i < r1.PerRank.Count; i++)
            {
                Assert.That(r1.PerRank[i].AverageStability,
                    Is.EqualTo(r2.PerRank[i].AverageStability).Within(1e-15),
                    $"Same seed must yield identical average stability at rank {r1.PerRank[i].Rank}.");
                Assert.That(r1.PerRank[i].CopheneticCorrelation,
                    Is.EqualTo(r2.PerRank[i].CopheneticCorrelation).Within(1e-15),
                    $"Same seed must yield identical cophenetic correlation at rank {r1.PerRank[i].Rank}.");
            }
        });
    }

    // S-RS5 — Single candidate rank: minRank == maxRank returns that rank with one diagnostic row.
    [Test]
    public void SelectRank_SingleCandidate_ReturnsThatRank()
    {
        double[][] v = BuildPlantedV();

        var result = OncologyAnalyzer.SelectRank(
            AsMatrix(v), minRank: 2, maxRank: 2, runs: 3, maxIterations: 1_000, tolerance: 1e-6, seed: 7);

        Assert.Multiple(() =>
        {
            Assert.That(result.SelectedRank, Is.EqualTo(2), "A single candidate rank must be the selected rank.");
            Assert.That(result.PerRank.Count, Is.EqualTo(1), "Exactly one diagnostic row for one candidate rank.");
        });
    }

    #endregion

    #region SelectRank — input validation

    [Test]
    public void SelectRank_NullMatrix_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.SelectRank(null!, minRank: 1, maxRank: 2),
            "A null count matrix must throw ArgumentNullException.");
    }

    [Test]
    public void SelectRank_MinRankBelowOne_Throws()
    {
        double[][] v = BuildPlantedV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.SelectRank(AsMatrix(v), minRank: 0, maxRank: 2),
            "minRank < 1 must throw ArgumentException.");
    }

    [Test]
    public void SelectRank_MaxRankBelowMinRank_Throws()
    {
        double[][] v = BuildPlantedV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.SelectRank(AsMatrix(v), minRank: 3, maxRank: 2),
            "maxRank < minRank must throw ArgumentException.");
    }

    [Test]
    public void SelectRank_MaxRankAboveChannelCount_Throws()
    {
        double[][] v = BuildPlantedV(); // 4 channels
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.SelectRank(AsMatrix(v), minRank: 1, maxRank: 5),
            "maxRank greater than the channel count must throw ArgumentException.");
    }

    [Test]
    public void SelectRank_RunsBelowOne_Throws()
    {
        double[][] v = BuildPlantedV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.SelectRank(AsMatrix(v), minRank: 1, maxRank: 2, runs: 0),
            "runs < 1 must throw ArgumentException.");
    }

    [Test]
    public void SelectRank_StabilityThresholdOutOfRange_Throws()
    {
        double[][] v = BuildPlantedV();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.SelectRank(AsMatrix(v), minRank: 1, maxRank: 2, stabilityThreshold: 1.5),
            "stabilityThreshold outside [0,1] must throw ArgumentException.");
    }

    #endregion

    #region MatchToReferenceSignatures — cosine matching (Alexandrov 2013/2020; SigProfiler)

    // M-MT1 — A positively-scaled copy of a reference matches that reference with cosine 1.0 (cosine is
    // scale-invariant), and an UNRELATED signature matches it with low cosine (SigProfiler cosine matching).
    [Test]
    public void MatchToReferenceSignatures_ScaledCopyAndUnrelated_MatchCorrectly()
    {
        var references = AsSignatures(TrueSignature0, TrueSignature1);

        // Extracted[0] = 5 * TrueSignature0 (a positively scaled copy of reference 0).
        var scaledCopyOfRef0 = TrueSignature0.Select(x => 5.0 * x).ToArray();
        // Extracted[1] = TrueSignature1 itself (exact copy of reference 1).
        var extracted = AsSignatures(scaledCopyOfRef0, TrueSignature1);

        var matches = OncologyAnalyzer.MatchToReferenceSignatures(extracted, references);

        Assert.Multiple(() =>
        {
            Assert.That(matches[0].ReferenceIndex, Is.EqualTo(0),
                "A scaled copy of reference 0 must match reference 0 (cosine is scale-invariant).");
            Assert.That(matches[0].CosineSimilarity, Is.EqualTo(1.0).Within(1e-12),
                "Cosine of a positively-scaled copy to its reference must be exactly 1.0.");
            Assert.That(matches[1].ReferenceIndex, Is.EqualTo(1),
                "Reference 1's exact copy must match reference 1.");
            Assert.That(matches[1].CosineSimilarity, Is.EqualTo(1.0).Within(1e-12),
                "Cosine of an exact copy to its reference must be exactly 1.0.");
            // The scaled copy of ref0 must NOT match ref1 closely: cos(TrueSignature0, TrueSignature1) is low.
            double crossCosine = Cosine(TrueSignature0, TrueSignature1);
            Assert.That(crossCosine, Is.LessThan(0.5),
                "The two unrelated planted signatures must have low cosine (distinct mutational processes).");
        });
    }

    // M-MT2 — A channel-permuted copy of a reference matches that permuted form, demonstrating that matching is
    // by cosine direction; here we verify the best match is the reference equal to the (un-permuted) query.
    [Test]
    public void MatchToReferenceSignatures_OneExtractedPerMatch_LabelsEachWithClosest()
    {
        var references = AsSignatures(
            new[] { 0.6, 0.2, 0.1, 0.1 },
            new[] { 0.1, 0.1, 0.2, 0.6 });

        // Extracted is exactly reference 1; must be labelled with reference index 1, cosine 1.0.
        var extracted = AsSignatures(new[] { 0.1, 0.1, 0.2, 0.6 });

        var matches = OncologyAnalyzer.MatchToReferenceSignatures(extracted, references);

        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1), "One match per extracted signature.");
            Assert.That(matches[0].ExtractedIndex, Is.EqualTo(0), "Match must reference its extracted index.");
            Assert.That(matches[0].ReferenceIndex, Is.EqualTo(1), "Extracted == reference 1 must match reference 1.");
            Assert.That(matches[0].CosineSimilarity, Is.EqualTo(1.0).Within(1e-12),
                "An exact copy must match with cosine 1.0.");
        });
    }

    [Test]
    public void MatchToReferenceSignatures_NullExtracted_Throws()
    {
        var references = AsSignatures(TrueSignature0);
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.MatchToReferenceSignatures(null!, references),
            "Null extracted signatures must throw ArgumentNullException.");
    }

    [Test]
    public void MatchToReferenceSignatures_ChannelCountMismatch_Throws()
    {
        var extracted = AsSignatures(new[] { 0.5, 0.5 });           // 2 channels
        var references = AsSignatures(new[] { 0.25, 0.25, 0.5 });   // 3 channels
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.MatchToReferenceSignatures(extracted, references),
            "Differing channel counts must throw ArgumentException.");
    }

    [Test]
    public void MatchToReferenceSignatures_EmptyReferences_Throws()
    {
        var extracted = AsSignatures(TrueSignature0);
        var emptyRefs = new List<IReadOnlyList<double>>();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.MatchToReferenceSignatures(extracted, emptyRefs),
            "An empty reference set must throw ArgumentException.");
    }

    #endregion
}
