using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for Oncology mutational-signature FITTING / REFITTING — ONCO-SIG-002.
/// The unit under test is the deterministic NNLS signature refit
/// <see cref="OncologyAnalyzer.FitSignatures"/> (deconvolve an observed catalog d against a caller-supplied
/// reference signature matrix S, returning per-signature exposures x ≥ 0, the proportion-normalised
/// exposures, the reconstruction S·x, and the reconstruction cosine), together with its helpers
/// <see cref="OncologyAnalyzer.ReconstructCatalog"/> (S·x) and
/// <see cref="OncologyAnalyzer.CosineSimilarity"/>, implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This is SIG-002 (the fit / exposure estimation), distinct from SIG-001 (the SBS-96 context CATALOGUE
/// builder, OncologySignatureContextFuzzTests). SIG-001 produces the catalogue d; SIG-002 deconvolves it.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts that the code NEVER fails
/// in an undisciplined way: no hang (the NNLS active-set loop must terminate), no NaN / Infinity leaking
/// out of a normalisation-by-total or a matrix solve, no exception from "inverting" a matrix that has no
/// inverse, and no result that violates the model's non-negativity constraint. Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a documented, intentional outcome (an
/// <see cref="ArgumentNullException"/> / <see cref="ArgumentException"/> for null / empty / ragged /
/// dimension-mismatched inputs). — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SIG-002 — mutational-signature fitting (NNLS refit), Oncology
/// Checklist: docs/checklists/03_FUZZING.md, row 97.
/// Fuzz strategy exercised: BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty).
///   Targets (checklist row 97): "zero catalogue, singular signature matrix, negative counts".
///     • ZERO CATALOGUE — an all-zero observed catalog d = 0. The only feasible minimiser of ‖S·x‖² with
///       x ≥ 0 is x = 0, so exposures, reconstruction and normalised exposures are ALL zero (§6.1). The
///       proportion normalisation must NOT divide by a zero total (no NaN / DivideByZero, INV-06), and the
///       reconstruction cosine of a zero-norm vector must be the documented 0.0 (§3.3), not NaN.
///     • SINGULAR SIGNATURE MATRIX — linearly dependent reference signatures (identical columns, an
///       all-equal matrix, a zero column, or scalar multiples). The passive-set normal-equations matrix
///       SᵀS is then rank-deficient / non-invertible. The solver must NOT throw from inverting a matrix
///       that has no inverse (§5.2: "a singular passive-set matrix … leaves the affected component at 0
///       rather than throwing"); it must still converge to a non-negative solution with NO Inf / NaN, and
///       still reconstruct a catalog that lies in the signatures' span (INV-04, INV-05).
///     • NEGATIVE COUNTS — a malformed catalog with negative entries. The fit is unconstrained in d but the
///       exposures are constrained x ≥ 0; the result must therefore NEVER contain a negative exposure
///       (INV-04, "negative contributions make no biological sense" [2]) and must stay finite — the
///       non-negativity constraint clamps the contribution rather than letting it go negative.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Mutational_Signature_Fitting.md (docs/algorithms/Oncology/Mutational_Signature_Fitting.md):
///   • Objective: minₓ ‖S·x − d‖₂², subject to x ≥ 0 (NNLS); S columns = reference signatures.   (§2.2)
///   • Solved by the Lawson-Hanson active-set algorithm; ε = 1e-12; finite termination.            (§4.2 [3])
///   • INV-04: every fitted exposure ≥ 0 (the NNLS constraint).                                     (§2.4)
///   • INV-05: ‖S·x − d‖² ≤ ‖d‖² (x = 0 is feasible, so the minimiser is no worse).                (§2.4)
///   • INV-06: normalised exposures sum to 1 when Σ>0, else ALL 0 (no ÷0).                          (§2.4)
///   • INV-07: NNLS = unconstrained LS when the latter is already non-negative.                     (§2.4)
///   • Reconstruction R = S·x; reconstruction cosine = sim(d, S·x) (≥0.95 ⇒ good fit).              (§2.2)
///   • Cosine of a zero-norm vector is 0.0 (÷0 undefined; treated as no shared direction).          (§3.3)
///   • Zero catalog d = 0 ⇒ exposures all 0, reconstruction 0, proportions all 0.                   (§6.1)
///   • Singular (collinear) passive-set matrix ⇒ affected component left at 0, no throw.            (§5.2, §6.2)
///   • null / empty / ragged / dimension-mismatched ⇒ ArgumentNullException / ArgumentException.    (§3.3)
///   • Worked example (§7.1): S=[[1,1],[0,1]] (sig1=[1,0], sig2=[1,1]), d=[0,1] ⇒ x=[0, 0.5].
///   • API example (§7.1): d=[3,5], S=I₂ ⇒ x=[3,5], reconstruction cosine = 1.0.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySignatureFittingFuzzTests
{
    // Tolerance for exact-arithmetic NNLS results on small, well-conditioned matrices.
    private const double Eps = 1e-9;

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY successful fit: the exposure / normalised / recon
    // vectors have the documented lengths; every value is FINITE (no NaN / ±Infinity from a divide-by-zero
    // normalisation or a singular matrix solve); every exposure is ≥ 0 (INV-04); the reconstruction cosine
    // is finite and in the documented [0,1] band (or exactly 0 for a degenerate zero-norm vector). This is
    // what stops a fuzz test from rubber-stamping a NaN / negative-exposure result green.
    private static void AssertWellFormedFit(SignatureFitResult fit, int signatureCount, int channelCount)
    {
        fit.Exposures.Should().HaveCount(signatureCount, "one exposure per reference signature (§3.2)");
        fit.NormalizedExposures.Should().HaveCount(signatureCount, "one normalised exposure per signature");
        fit.Reconstruction.Should().HaveCount(channelCount, "the reconstruction S·x spans all channels (§3.2)");

        foreach (double e in fit.Exposures)
        {
            double.IsFinite(e).Should().BeTrue("an exposure must never be NaN/Infinity (no ÷0, no bad solve)");
            e.Should().BeGreaterThanOrEqualTo(-Eps, "every fitted exposure is ≥ 0 — NNLS constraint (INV-04)");
        }

        foreach (double e in fit.NormalizedExposures)
        {
            double.IsFinite(e).Should().BeTrue("a normalised exposure must never be NaN/Infinity (no ÷0, INV-06)");
            e.Should().BeGreaterThanOrEqualTo(-Eps, "normalised exposures are non-negative proportions");
        }

        foreach (double r in fit.Reconstruction)
        {
            double.IsFinite(r).Should().BeTrue("a reconstruction entry must never be NaN/Infinity");
        }

        double cosine = fit.ReconstructionCosineSimilarity;
        double.IsFinite(cosine).Should().BeTrue("the reconstruction cosine must never be NaN/Infinity (§3.3)");
        cosine.Should().BeInRange(-Eps, 1.0 + Eps, "cosine of non-negative vectors lies in [0,1] (INV-01)");

        // INV-06: normalised exposures sum to 1 when the exposure total is positive, else are all zero.
        double exposureTotal = fit.Exposures.Sum();
        double normalisedSum = fit.NormalizedExposures.Sum();
        if (exposureTotal > Eps)
        {
            normalisedSum.Should().BeApproximately(1.0, 1e-6,
                "normalised exposures sum to 1 when the total is positive (INV-06)");
        }
        else
        {
            normalisedSum.Should().BeApproximately(0.0, 1e-9,
                "a zero exposure total ⇒ all-zero proportions, never a ÷0 NaN (INV-06)");
        }
    }

    // Builds a random non-negative reference signature matrix: `count` signatures, each of length `channels`.
    private static IReadOnlyList<double>[] RandomSignatures(Random rng, int count, int channels, double scale = 1.0)
    {
        var signatures = new IReadOnlyList<double>[count];
        for (int j = 0; j < count; j++)
        {
            var sig = new double[channels];
            for (int k = 0; k < channels; k++)
            {
                sig[k] = rng.NextDouble() * scale;
            }

            signatures[j] = sig;
        }

        return signatures;
    }

    // d = Σ_j weights[j] · signatures[j]  (a known non-negative combination — the synthetic ground truth).
    private static double[] Combine(IReadOnlyList<IReadOnlyList<double>> signatures, IReadOnlyList<double> weights)
    {
        int channels = signatures[0].Count;
        var d = new double[channels];
        for (int j = 0; j < signatures.Count; j++)
        {
            for (int k = 0; k < channels; k++)
            {
                d[k] += signatures[j][k] * weights[j];
            }
        }

        return d;
    }

    #region ONCO-SIG-002 — positive sanity (exact documented recovery of known exposures)

    // The documented API worked example (§7.1): catalog [3,5] against the identity signature matrix recovers
    // exposures [3,5] exactly with a perfect (1.0) reconstruction cosine.
    [Test]
    public void FitSignatures_IdentityBasis_RecoversExposuresExactly_WorkedExample()
    {
        var catalog = new double[] { 3, 5 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 2);
        fit.Exposures[0].Should().BeApproximately(3.0, Eps, "identity basis ⇒ exposure = catalog count (§7.1)");
        fit.Exposures[1].Should().BeApproximately(5.0, Eps);
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-9,
            "an exactly representable catalog reconstructs perfectly (§7.1)");
    }

    // The documented NNLS clamp worked example (§7.1): S = [[1,1],[0,1]], d = [0,1]. The unconstrained
    // optimum x₁ = −1 is infeasible, so sig1 is clamped to 0 and the fit is x = [0, 0.5] — proving the
    // non-negativity constraint actively clamps a would-be-negative exposure.
    [Test]
    public void FitSignatures_NnlsClamp_NegativeUnconstrainedCoefficient_ClampedToZero_WorkedExample()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 1, 1 } };
        var catalog = new double[] { 0, 1 };

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 2);
        fit.Exposures[0].Should().BeApproximately(0.0, Eps,
            "the unconstrained coefficient is negative ⇒ clamped to 0, never returned negative (§7.1, INV-04)");
        fit.Exposures[1].Should().BeApproximately(0.5, Eps, "sig2 is refit alone to 1/2 (§7.1)");
        fit.Reconstruction[0].Should().BeApproximately(0.5, Eps);
        fit.Reconstruction[1].Should().BeApproximately(0.5, Eps);
    }

    // A catalog synthesised as a known non-negative mixture 0.7·SigA + 0.3·SigB of two LINEARLY INDEPENDENT
    // signatures must be recovered with exposures close to [0.7, 0.3] and a (near-)perfect reconstruction.
    [Test]
    public void FitSignatures_KnownMixture_RecoversComponentExposures()
    {
        // Two clearly independent profiles over 4 channels.
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 }, // SigA
            new double[] { 0.1, 0.1, 0.3, 0.5 }, // SigB
        };
        var catalog = Combine(signatures, new[] { 0.7, 0.3 });

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 4);
        fit.Exposures[0].Should().BeApproximately(0.7, 1e-6, "0.7·SigA recovered (NNLS exact convex optimum)");
        fit.Exposures[1].Should().BeApproximately(0.3, 1e-6, "0.3·SigB recovered");
        fit.NormalizedExposures[0].Should().BeApproximately(0.7, 1e-6, "proportions = 0.7 / (0.7+0.3) (INV-06)");
        fit.NormalizedExposures[1].Should().BeApproximately(0.3, 1e-6);
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-9,
            "a catalog in the span reconstructs perfectly (INV-05 with residual ≈ 0)");
    }

    // A single-signature catalog (d = c·SigA, c > 0) recovers exactly that one signature's exposure and a
    // zero exposure for the other, present-but-unused signature.
    [Test]
    public void FitSignatures_SingleSignatureCatalog_RecoversThatSignatureOnly()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.5, 0.3, 0.2, 0.0 }, // SigA (the source)
            new double[] { 0.0, 0.1, 0.4, 0.5 }, // SigB (independent, should be unused)
        };
        var catalog = Combine(signatures, new[] { 2.0, 0.0 }); // pure SigA, weight 2

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 4);
        fit.Exposures[0].Should().BeApproximately(2.0, 1e-6, "the source signature is recovered with its weight");
        fit.Exposures[1].Should().BeApproximately(0.0, 1e-6, "the unused signature gets a zero exposure (INV-04)");
        fit.NormalizedExposures[0].Should().BeApproximately(1.0, 1e-6, "100% of the catalog is SigA (INV-06)");
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-9);
    }

    // INV-07: when the unconstrained least-squares solution is already non-negative, NNLS returns exactly it.
    // Fuzz over random catalogs against the identity basis — exposures must equal the catalog counts.
    [Test]
    [CancelAfter(20_000)]
    public void FitSignatures_NonNegativeUnconstrainedCase_EqualsUnconstrainedLs_IdentityBasis()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(2, 6);
            var signatures = new IReadOnlyList<double>[n];
            var catalog = new double[n];
            for (int i = 0; i < n; i++)
            {
                var sig = new double[n];
                sig[i] = 1.0; // identity column ⇒ exposures = catalog
                signatures[i] = sig;
                catalog[i] = rng.NextDouble() * 100.0; // non-negative counts
            }

            var fit = FitSignatures(catalog, signatures);
            AssertWellFormedFit(fit, n, n);
            for (int i = 0; i < n; i++)
            {
                fit.Exposures[i].Should().BeApproximately(catalog[i], 1e-6,
                    $"seed {seed}: identity basis ⇒ NNLS = unconstrained LS = catalog (INV-07)");
            }
        }
    }

    #endregion

    #region ONCO-SIG-002 — BE: zero catalogue (d = 0 ⇒ all-zero fit, no ÷0 NaN)

    // The all-zero catalog: the only feasible minimiser of ‖S·x‖² with x ≥ 0 is x = 0, so exposures,
    // reconstruction and normalised exposures are ALL exactly zero (§6.1). The proportion normalisation must
    // NOT produce NaN from 0/0 (INV-06), and the reconstruction cosine of a zero-norm vector is exactly 0.0
    // (§3.3) — never NaN.
    [Test]
    public void FitSignatures_ZeroCatalog_AllZeroExposures_NoNaN_NoDivideByZero()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new double[4]; // all zero

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 4);
        fit.Exposures.Should().OnlyContain(e => e == 0.0, "x = 0 is the only feasible minimiser of ‖Sx‖² (§6.1)");
        fit.Reconstruction.Should().OnlyContain(r => r == 0.0, "S·0 = 0 (§6.1)");
        fit.NormalizedExposures.Should().OnlyContain(e => e == 0.0,
            "a zero exposure total ⇒ all-zero proportions, NOT a 0/0 NaN (INV-06)");
        fit.ReconstructionCosineSimilarity.Should().Be(0.0,
            "the cosine of a zero-norm reconstruction is the documented 0.0, never NaN (§3.3)");
    }

    // Fuzz: an all-zero catalog against random (varying size / scale) NON-degenerate signature matrices is
    // always the all-zero fit with no NaN / DivideByZero, regardless of the signature shape.
    [Test]
    [CancelAfter(20_000)]
    public void FitSignatures_ZeroCatalog_RandomSignatureMatrices_AlwaysAllZero_NoNaN()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(1, 6);          // 1..5 signatures
            int channels = rng.Next(1, 12);  // 1..11 channels
            var signatures = RandomSignatures(rng, k, channels, scale: rng.NextDouble() * 100.0 + 0.01);
            var catalog = new double[channels]; // all zero

            SignatureFitResult fit = default;
            FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
                .Should().NotThrow($"seed {seed}: a zero catalog is a valid degenerate input");

            AssertWellFormedFit(fit, k, channels);
            fit.Exposures.Should().OnlyContain(e => e == 0.0, $"seed {seed}: x = 0 is the only minimiser");
            fit.NormalizedExposures.Should().OnlyContain(e => e == 0.0, $"seed {seed}: no 0/0 NaN (INV-06)");
            fit.ReconstructionCosineSimilarity.Should().Be(0.0, $"seed {seed}: zero-norm recon cosine = 0.0 (§3.3)");
        }
    }

    #endregion

    #region ONCO-SIG-002 — BE: singular signature matrix (collinear/rank-deficient ⇒ no throw, no Inf/NaN)

    // IDENTICAL signatures: S has two identical columns, so SᵀS is singular (non-invertible). The solver must
    // NOT throw from "inverting" it (§5.2) and must still converge to a non-negative solution that
    // reconstructs the catalog (which lies exactly in the degenerate span).
    [Test]
    public void FitSignatures_IdenticalSignatures_SingularMatrix_NoThrow_NonNegative_Reconstructs()
    {
        var sig = new double[] { 0.4, 0.3, 0.2, 0.1 };
        var signatures = new IReadOnlyList<double>[] { (double[])sig.Clone(), (double[])sig.Clone() };
        var catalog = new double[] { 0.8, 0.6, 0.4, 0.2 }; // = 2·sig — in the (degenerate) span

        SignatureFitResult fit = default;
        FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
            .Should().NotThrow("a singular (collinear) signature matrix must not throw from a non-existent inverse (§5.2)");

        AssertWellFormedFit(fit, 2, 4);
        // The decomposition is non-unique (§6.2) but the RECONSTRUCTION is determined: S·x must equal the
        // catalog (cosine 1) since it lies in the span — that is the testable, well-defined contract.
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-6,
            "the catalog lies in the (collinear) span ⇒ reconstructed exactly despite the singular matrix (INV-05)");
        (fit.Exposures[0] + fit.Exposures[1]).Should().BeApproximately(2.0, 1e-6,
            "the total contribution along the shared direction is determined (= 2), even if the split is not");
    }

    // A ZERO signature column (a degenerate / empty mutational process) makes SᵀS singular. The fit must not
    // throw, the zero signature must receive a (non-negative, finite) exposure that contributes nothing, and
    // the real signature must still be recovered.
    [Test]
    public void FitSignatures_ZeroSignatureColumn_SingularMatrix_NoThrow_RecoversRealSignature()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.5, 0.3, 0.2, 0.0 }, // real
            new double[4],                       // all-zero column ⇒ singular SᵀS
        };
        var catalog = Combine(signatures, new[] { 3.0, 0.0 });

        SignatureFitResult fit = default;
        FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
            .Should().NotThrow("a zero signature column makes the matrix singular but must not throw (§5.2)");

        AssertWellFormedFit(fit, 2, 4);
        fit.Exposures[0].Should().BeApproximately(3.0, 1e-6, "the real signature is recovered with its weight");
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-6,
            "the catalog reconstructs exactly; the zero column contributes nothing (INV-05)");
    }

    // SCALAR-MULTIPLE (collinear) signatures: column 2 = 2·column 1 — rank-1 matrix, singular SᵀS. No throw,
    // non-negative finite solution, catalog (in the shared direction) reconstructed.
    [Test]
    public void FitSignatures_ScalarMultipleSignatures_RankDeficient_NoThrow_Reconstructs()
    {
        var baseSig = new double[] { 0.3, 0.2, 0.4, 0.1 };
        var signatures = new IReadOnlyList<double>[]
        {
            (double[])baseSig.Clone(),
            baseSig.Select(v => 2.0 * v).ToArray(), // collinear (2×)
        };
        var catalog = baseSig.Select(v => 5.0 * v).ToArray(); // 5·baseSig — in the span

        SignatureFitResult fit = default;
        FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
            .Should().NotThrow("a rank-deficient (scalar-multiple) signature matrix must not throw (§5.2)");

        AssertWellFormedFit(fit, 2, 4);
        fit.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-6,
            "the catalog is in the (rank-1) span ⇒ reconstructed exactly (INV-05)");
    }

    // An ALL-EQUAL signature matrix (every column the same constant vector) — maximally collinear. No throw,
    // no Inf/NaN, non-negative finite solution across a fuzzed range of sizes.
    [Test]
    [CancelAfter(30_000)]
    public void FitSignatures_AllEqualOrCollinearMatrices_Fuzz_NoThrow_NoInfNaN_NonNegative()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(2, 6);          // ≥2 signatures so the matrix can be rank-deficient
            int channels = rng.Next(1, 12);

            // Build a random base column, then make EVERY signature a (random non-negative) scalar multiple
            // of it ⇒ a rank-1, maximally collinear matrix ⇒ singular SᵀS.
            var baseCol = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                baseCol[c] = rng.NextDouble() * 10.0;
            }

            var signatures = new IReadOnlyList<double>[k];
            for (int j = 0; j < k; j++)
            {
                double mult = rng.NextDouble() * 5.0 + 0.001;
                signatures[j] = baseCol.Select(v => v * mult).ToArray();
            }

            // A catalog along the shared direction (in the span) plus, sometimes, an out-of-span component.
            var catalog = new double[channels];
            double catScale = rng.NextDouble() * 7.0;
            for (int c = 0; c < channels; c++)
            {
                catalog[c] = baseCol[c] * catScale + (rng.Next(3) == 0 ? rng.NextDouble() : 0.0);
            }

            SignatureFitResult fit = default;
            FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
                .Should().NotThrow($"seed {seed}: a singular/collinear matrix must never throw (§5.2)");

            AssertWellFormedFit(fit, k, channels);
        }
    }

    #endregion

    #region ONCO-SIG-002 — BE: negative counts (malformed catalog ⇒ exposures still ≥ 0, finite)

    // A catalog with NEGATIVE entries is malformed ("negative contributions make no biological sense" [2]).
    // The fit is unconstrained in d, but the non-negativity constraint is on the EXPOSURES x ≥ 0 — so the
    // result must NEVER contain a negative exposure (INV-04) and must stay finite (no NaN/Inf). A fully
    // negative catalog pulls the unconstrained optimum into the negative orthant; NNLS must clamp it to 0.
    [Test]
    public void FitSignatures_FullyNegativeCatalog_ExposuresClampedToZero_NeverNegative()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new double[] { -1.0, -2.0, -3.0, -4.0 }; // wholly negative ⇒ best non-neg fit is 0

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 4); // already asserts every exposure ≥ 0 and finite (INV-04)
        fit.Exposures.Should().OnlyContain(e => e <= Eps,
            "a wholly negative catalog ⇒ the non-negativity constraint clamps every exposure to 0 (INV-04)");
        fit.NormalizedExposures.Should().OnlyContain(e => e == 0.0,
            "a zero exposure total ⇒ all-zero proportions, no ÷0 NaN (INV-06)");
    }

    // A MIXED-SIGN catalog (some entries negative, some positive). Exposures must still be non-negative and
    // finite — the malformed negative entries must not drive any exposure negative (INV-04).
    [Test]
    public void FitSignatures_MixedSignCatalog_ExposuresNonNegative_Finite()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.5, 0.3, 0.2, 0.0 },
            new double[] { 0.0, 0.1, 0.4, 0.5 },
        };
        var catalog = new double[] { 4.0, -1.0, 2.0, -3.0 };

        var fit = FitSignatures(catalog, signatures);

        AssertWellFormedFit(fit, 2, 4);
        // No additional shape assertion beyond the well-formed helper: the KEY guarantee under fuzzing is
        // that malformed negative counts never break the non-negativity constraint or produce NaN/Inf.
    }

    // Fuzz: random catalogs with arbitrary (possibly negative) entries against random signature matrices —
    // the result is ALWAYS well-formed (finite, non-negative exposures, INV-06 proportions) and the NNLS
    // loop ALWAYS terminates. This is the core BE "negative counts" robustness sweep.
    [Test]
    [CancelAfter(30_000)]
    public void FitSignatures_RandomNegativeAndPositiveCatalogs_Fuzz_AlwaysNonNegativeFinite()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(1, 6);
            int channels = rng.Next(1, 12);
            var signatures = RandomSignatures(rng, k, channels, scale: rng.NextDouble() * 50.0 + 0.01);

            var catalog = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                // Centred around zero so ~half the entries are negative (the malformed-count regime).
                catalog[c] = (rng.NextDouble() - 0.5) * 200.0;
            }

            SignatureFitResult fit = default;
            FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
                .Should().NotThrow($"seed {seed}: a negative-count catalog is malformed but must not crash");

            AssertWellFormedFit(fit, k, channels); // finite + every exposure ≥ 0 (INV-04) + INV-06
        }
    }

    // Combined adversary: NEGATIVE counts AND a SINGULAR (collinear) signature matrix together — the two
    // BE hazards at once. Still no throw, no Inf/NaN, non-negative exposures.
    [Test]
    [CancelAfter(30_000)]
    public void FitSignatures_NegativeCatalog_WithSingularMatrix_Fuzz_StillRobust()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(2, 5);
            int channels = rng.Next(1, 10);

            var baseCol = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                baseCol[c] = rng.NextDouble() * 10.0;
            }

            var signatures = new IReadOnlyList<double>[k];
            for (int j = 0; j < k; j++)
            {
                double mult = rng.NextDouble() * 4.0 + 0.001;
                signatures[j] = baseCol.Select(v => v * mult).ToArray(); // collinear ⇒ singular
            }

            var catalog = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                catalog[c] = (rng.NextDouble() - 0.5) * 100.0; // mixed sign
            }

            SignatureFitResult fit = default;
            FluentActions.Invoking(() => fit = FitSignatures(catalog, signatures))
                .Should().NotThrow($"seed {seed}: negative counts + singular matrix must not crash (§5.2, INV-04)");

            AssertWellFormedFit(fit, k, channels);
        }
    }

    #endregion

    #region ONCO-SIG-002 — BE: structural boundaries (empty/ragged/dimension-mismatch ⇒ documented throws)

    // null catalog / null signatures ⇒ ArgumentNullException (§3.3).
    [Test]
    public void FitSignatures_NullArguments_ThrowArgumentNull()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };

        FluentActions.Invoking(() => FitSignatures(null!, signatures))
            .Should().Throw<ArgumentNullException>("a null catalog is a documented guard (§3.3)");
        FluentActions.Invoking(() => FitSignatures(new double[] { 1, 2 }, null!))
            .Should().Throw<ArgumentNullException>("a null signature matrix is a documented guard (§3.3)");
    }

    // Empty signature matrix (no signatures at all) ⇒ ArgumentException (§3.3) — the boundary "no columns".
    [Test]
    public void FitSignatures_NoSignatures_ThrowsArgumentException()
    {
        FluentActions.Invoking(() => FitSignatures(new double[] { 1, 2 }, Array.Empty<IReadOnlyList<double>>()))
            .Should().Throw<ArgumentException>("at least one reference signature is required (§3.3)");
    }

    // Empty (zero-channel) signatures ⇒ ArgumentException — the boundary "no rows".
    [Test]
    public void FitSignatures_EmptySignatureVectors_ThrowsArgumentException()
    {
        var signatures = new IReadOnlyList<double>[] { Array.Empty<double>() };
        FluentActions.Invoking(() => FitSignatures(Array.Empty<double>(), signatures))
            .Should().Throw<ArgumentException>("signature vectors cannot be empty (§3.3)");
    }

    // Ragged signatures (unequal lengths) ⇒ ArgumentException (§3.3).
    [Test]
    public void FitSignatures_RaggedSignatures_ThrowsArgumentException()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0 },
            new double[] { 0, 1 }, // shorter ⇒ ragged
        };
        FluentActions.Invoking(() => FitSignatures(new double[] { 1, 2, 3 }, signatures))
            .Should().Throw<ArgumentException>("all signatures must have the same length (§3.3)");
    }

    // Catalog length ≠ signature channel count ⇒ ArgumentException (§3.3) — the dimension-mismatch boundary.
    [Test]
    public void FitSignatures_CatalogLengthMismatch_ThrowsArgumentException()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        FluentActions.Invoking(() => FitSignatures(new double[] { 1, 2, 3 }, signatures))
            .Should().Throw<ArgumentException>("catalog length must equal the signature channel count (§3.3)");
    }

    // Fuzz: random null / ragged / dimension-mismatch corruptions always raise the documented Argument*
    // exception family — never an undisciplined IndexOutOfRange / NullReference crash.
    [Test]
    [CancelAfter(20_000)]
    public void FitSignatures_RandomMalformedStructure_AlwaysDocumentedArgumentException()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int channels = rng.Next(1, 8);
            int k = rng.Next(1, 5);
            var signatures = RandomSignatures(rng, k, channels);

            int defect = rng.Next(3);
            Action act;
            switch (defect)
            {
                case 0: // catalog length mismatch
                    int wrongLen = channels + (rng.Next(2) == 0 ? 1 : -1);
                    wrongLen = Math.Max(0, wrongLen);
                    if (wrongLen == channels)
                    {
                        wrongLen = channels + 1;
                    }

                    var badCatalog = new double[wrongLen];
                    act = () => FitSignatures(badCatalog, signatures);
                    break;

                case 1: // ragged signatures (corrupt one signature's length)
                    var ragged = signatures.ToArray();
                    int corrupt = rng.Next(k);
                    ragged[corrupt] = new double[channels + 1];
                    var cat1 = new double[channels];
                    act = () => FitSignatures(cat1, ragged);
                    break;

                default: // null signature vector
                    var withNull = signatures.ToArray();
                    withNull[rng.Next(k)] = null!;
                    var cat2 = new double[channels];
                    act = () => FitSignatures(cat2, withNull);
                    break;
            }

            FluentActions.Invoking(act).Should().Throw<ArgumentException>(
                $"seed {seed}: malformed structure (defect {defect}) ⇒ documented ArgumentException, not a crash");
        }
    }

    #endregion
}
