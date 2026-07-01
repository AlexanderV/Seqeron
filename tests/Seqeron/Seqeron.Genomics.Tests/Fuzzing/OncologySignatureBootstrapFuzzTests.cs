using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for Oncology mutational-signature exposure BOOTSTRAP confidence intervals — ONCO-SIG-003.
/// The unit under test is the parametric (multinomial) bootstrap
/// <see cref="OncologyAnalyzer.BootstrapExposures"/>: an observed integer mutational catalog is repeatedly
/// resampled as a draw of N = Σ catalog mutations from the multinomial distribution with per-channel
/// probabilities pₖ = catalogₖ / N, each resampled catalog is refit to the caller-supplied reference
/// signatures by NNLS (<see cref="OncologyAnalyzer.FitSignatures"/>, ONCO-SIG-002), and a two-sided
/// percentile confidence interval is taken per signature from the resulting bootstrap exposure
/// distribution. Implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This is SIG-003 (resampling-based exposure CONFIDENCE / stability), distinct from SIG-002 (the
/// deterministic NNLS point fit, OncologySignatureFittingFuzzTests) and SIG-001 (the SBS-96 context
/// CATALOGUE builder, OncologySignatureContextFuzzTests). SIG-003 calls SIG-002's fit once per replicate;
/// these tests target the bootstrap machinery (resampling, percentile CI, determinism), not the NNLS fit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts that the code NEVER fails in
/// an undisciplined way: no hang (the per-replicate resampling + NNLS loop must terminate), no DivideByZero
/// when averaging / taking percentiles (e.g. over zero replicates, or over a total-mutation count of 0 or
/// 1), no NaN / Infinity leaking out of the mean / quantile, no overflow on an extreme seed, no inverted
/// (lower > upper) or negative interval. Every input must resolve to EITHER a well-defined, theory-correct
/// value OR a documented, intentional outcome (an Argument*Exception per the contract). Because the method
/// is Monte-Carlo, the headline contract is DETERMINISM UNDER A FIXED SEED: identical (catalog, signatures,
/// replicates, confidence, seed) ⇒ identical intervals. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SIG-003 — mutational-signature exposure bootstrap CIs, Oncology
/// Checklist: docs/checklists/03_FUZZING.md, row 98.
/// Fuzz strategy exercised: BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty).
///   Targets (checklist row 98): "0 bootstrap reps, 1 mutation, fixed-seed extremes".
///     • 0 BOOTSTRAP REPS — replicates &lt; 1 is the documented guard boundary. The mean / percentile would
///       average / index over an empty replicate distribution (a DivideByZero / empty-array hazard); the
///       contract (§3.3) is a clean ArgumentOutOfRangeException, NOT a crash or a NaN result.
///     • 1 MUTATION — a catalog with total N = 1 (a single mutation in one channel). The multinomial draw
///       over N = 1 is near-degenerate (it must place the single mutation in a channel; with one non-zero
///       channel every replicate is identical), and pₖ = dₖ/N has a denominator of 1 — there must be no
///       DivideByZero / NaN on the tiny total, the CI must stay well-defined, and the single-non-zero-channel
///       case is exactly deterministic (every replicate exposure = the point estimate, §6.1).
///     • FIXED-SEED EXTREMES — with a fixed seed the bootstrap is FULLY REPRODUCIBLE: the same seed ⇒
///       byte-identical intervals (INV-03); different seeds may differ but stay within the documented
///       non-negative, lower ≤ upper bounds. Extreme rep counts (replicates = 1, the §6.1 single-replicate
///       case) and extreme / edge seeds (0, -1, int.MinValue, int.MaxValue) must not overflow or crash —
///       Random(seed) must accept them and produce a well-formed, reproducible result.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Mutational_Signature_Exposure_Bootstrap.md (docs/algorithms/Oncology/Mutational_Signature_Exposure_Bootstrap.md):
///   • Procedure: resample d⁽ᵇ⁾ ~ Multinomial(N, p), refit by NNLS, summarise per signature.          (§2.2, §4.1)
///   • Output: one ExposureConfidenceInterval per signature in signature order; PointEstimate =
///     NNLS exposure of the un-resampled observed catalog; Mean over replicates; Lower/Upper =
///     [(1−c)/2, 1−(1−c)/2] type-7 percentiles of the replicate exposures.                            (§3.2, INV-04)
///   • INV-01: all exposures and bounds ≥ 0 (NNLS x ≥ 0; multinomial counts ≥ 0).                     (§2.4)
///   • INV-02: lower_j ≤ upper_j (the lower percentile ≤ the upper percentile).                       (§2.4)
///   • INV-03: determinism for fixed (d, S, R, c, seed) — the RNG is seeded from a fixed value.       (§2.4)
///   • INV-05: N = 0 (all-zero catalog) ⇒ every interval is [0,0] with point 0 and mean 0.            (§2.4, §6.1)
///   • §6.1 single non-zero channel ⇒ every replicate exposure = point ⇒ lower = upper = mean = point.
///   • §6.1 replicates = 1 ⇒ lower = upper = mean = the single replicate exposure.
///   • §3.3 validation: null catalog / signatures (or a null signature vector) ⇒ ArgumentNullException;
///     empty / ragged signatures, catalog length ≠ channel count, or a negative count ⇒ ArgumentException;
///     replicates &lt; 1 or confidence ∉ (0,1) ⇒ ArgumentOutOfRangeException.
///   • §7.1 walk-through: catalog [10] vs signature [[1.0]] ⇒ every draw = [10] ⇒ lower=upper=mean=point=10.
///   • Defaults: DefaultBootstrapReplicates = 1000, DefaultBootstrapConfidence = 0.95, DefaultBootstrapSeed = 42.
///
/// All randomness for BUILDING test catalogs is LOCALLY seeded (new Random(seed)); the bootstrap's own RNG
/// seed is an explicit argument passed to BootstrapExposures (never a shared static Rng).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySignatureBootstrapFuzzTests
{
    // Tolerance for exact-arithmetic NNLS / percentile results on small catalogs.
    private const double Eps = 1e-9;

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY successful bootstrap result: one interval per
    // signature (§3.2); every PointEstimate / Mean / Lower / Upper is FINITE (no NaN / ±Infinity from a
    // mean-over-zero, a 0/0 probability, or a bad quantile) and ≥ 0 (INV-01); Lower ≤ Upper (INV-02, the CI
    // is never inverted); and the recorded Confidence is the one requested. (The Mean is NOT bounded by the
    // percentile band — at a narrow confidence the [lower, upper] band can sit on one side of a skewed mean —
    // so the mean is only asserted finite and non-negative, never "within [lower, upper]".) This is what stops
    // a fuzz test from rubber-stamping a NaN / negative / inverted interval green.
    private static void AssertWellFormedIntervals(
        IReadOnlyList<ExposureConfidenceInterval> intervals, int signatureCount, double confidence)
    {
        intervals.Should().HaveCount(signatureCount, "one interval per reference signature (§3.2)");

        foreach (ExposureConfidenceInterval ci in intervals)
        {
            double.IsFinite(ci.PointEstimate).Should().BeTrue("the point estimate must never be NaN/Infinity");
            double.IsFinite(ci.Mean).Should().BeTrue("the replicate mean must never be NaN/Infinity (no ÷0 over reps)");
            double.IsFinite(ci.Lower).Should().BeTrue("the lower bound must never be NaN/Infinity (no bad quantile)");
            double.IsFinite(ci.Upper).Should().BeTrue("the upper bound must never be NaN/Infinity");

            ci.PointEstimate.Should().BeGreaterThanOrEqualTo(-Eps, "exposures are NNLS-constrained x ≥ 0 (INV-01)");
            ci.Mean.Should().BeGreaterThanOrEqualTo(-Eps, "the replicate mean of ≥ 0 exposures is ≥ 0 (INV-01)");
            ci.Lower.Should().BeGreaterThanOrEqualTo(-Eps, "the lower bound is a percentile of ≥ 0 exposures (INV-01)");
            ci.Upper.Should().BeGreaterThanOrEqualTo(-Eps, "the upper bound is a percentile of ≥ 0 exposures (INV-01)");

            ci.Lower.Should().BeLessThanOrEqualTo(ci.Upper + Eps,
                "the lower percentile must not exceed the upper percentile (INV-02, never inverted)");

            ci.Confidence.Should().BeApproximately(confidence, Eps, "the recorded confidence is the requested level (§3.2)");
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

    // Builds a random non-negative INTEGER catalog of length `channels`, each count in [0, maxCount].
    private static int[] RandomCatalog(Random rng, int channels, int maxCount)
    {
        var catalog = new int[channels];
        for (int k = 0; k < channels; k++)
        {
            catalog[k] = rng.Next(0, maxCount + 1);
        }

        return catalog;
    }

    #region ONCO-SIG-003 — positive sanity (a clear dominant signature ⇒ reproducible CI brackets the point)

    // The documented numerical walk-through (§7.1): catalog [10] against the single signature [[1.0]]. N = 10,
    // p = [1] ⇒ every multinomial draw assigns all 10 mutations to the single channel ⇒ d⁽ᵇ⁾ = [10] for all
    // replicates ⇒ NNLS exposure 10 every time. The bootstrap distribution is the constant {10}, so
    // lower = upper = mean = point = 10 at any confidence level.
    [Test]
    public void BootstrapExposures_SingleChannelSingleSignature_ConstantDistribution_WorkedExample()
    {
        var catalog = new[] { 10 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1.0 } };

        var intervals = BootstrapExposures(catalog, signatures, replicates: 200, confidence: 0.95, seed: 42);

        AssertWellFormedIntervals(intervals, 1, 0.95);
        ExposureConfidenceInterval ci = intervals[0];
        ci.PointEstimate.Should().BeApproximately(10.0, Eps, "NNLS of [10] onto [[1.0]] is 10 (§7.1)");
        ci.Mean.Should().BeApproximately(10.0, Eps, "every replicate draw = [10] ⇒ exposure 10 (§7.1)");
        ci.Lower.Should().BeApproximately(10.0, Eps, "the constant {10} distribution has lower = 10 (§7.1)");
        ci.Upper.Should().BeApproximately(10.0, Eps, "the constant {10} distribution has upper = 10 (§7.1)");
    }

    // A catalog with a clear DOMINANT signature: d is built almost entirely from SigA. The bootstrap CI for
    // SigA must BRACKET its point estimate (lower ≤ point ≤ upper) and the result must be reproducible across
    // two runs with the same seed; the dominant signature's interval must sit well above the minor one's.
    [Test]
    [CancelAfter(30_000)]
    public void BootstrapExposures_DominantSignature_CiBracketsPoint_AndReproducible()
    {
        // Two clearly independent profiles over 4 channels.
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 6, 2, 1, 1 }, // SigA (dominant source)
            new double[] { 1, 1, 3, 5 }, // SigB (minor)
        };
        // Catalog ≈ 50·SigA + 3·SigB, as integer counts.
        var catalog = new[] { 6 * 50 + 3, 2 * 50 + 3, 1 * 50 + 9, 1 * 50 + 15 };

        var first = BootstrapExposures(catalog, signatures, replicates: 200, confidence: 0.95, seed: 7);
        var second = BootstrapExposures(catalog, signatures, replicates: 200, confidence: 0.95, seed: 7);

        AssertWellFormedIntervals(first, 2, 0.95);

        // The CI brackets the point estimate for the dominant signature.
        first[0].Lower.Should().BeLessThanOrEqualTo(first[0].PointEstimate + Eps,
            "the bootstrap lower bound brackets the point estimate from below (INV-02 + percentile method)");
        first[0].Upper.Should().BeGreaterThanOrEqualTo(first[0].PointEstimate - Eps,
            "the bootstrap upper bound brackets the point estimate from above");
        first[0].PointEstimate.Should().BeGreaterThan(first[1].PointEstimate,
            "the dominant signature carries far more exposure than the minor one");

        // Determinism (INV-03): same seed ⇒ byte-identical intervals.
        for (int j = 0; j < first.Count; j++)
        {
            first[j].Should().Be(second[j],
                $"signature {j}: a fixed seed makes the bootstrap fully reproducible (INV-03)");
        }
    }

    // Different seeds MAY produce different intervals, but each remains well-formed (non-negative,
    // lower ≤ upper) and brackets the same deterministic point estimate (the point does not depend on seed).
    [Test]
    [CancelAfter(30_000)]
    public void BootstrapExposures_DifferentSeeds_DifferentCis_ButPointEstimateSeedIndependent()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new[] { 40, 18, 14, 18 };

        var a = BootstrapExposures(catalog, signatures, replicates: 200, confidence: 0.95, seed: 1);
        var b = BootstrapExposures(catalog, signatures, replicates: 200, confidence: 0.95, seed: 2);

        AssertWellFormedIntervals(a, 2, 0.95);
        AssertWellFormedIntervals(b, 2, 0.95);

        for (int j = 0; j < a.Count; j++)
        {
            a[j].PointEstimate.Should().BeApproximately(b[j].PointEstimate, Eps,
                $"signature {j}: the point estimate is the NNLS fit of the observed catalog — seed-independent (INV-04)");
        }
    }

    #endregion

    #region ONCO-SIG-003 — BE: 0 bootstrap reps (replicates < 1 ⇒ documented guard, no ÷0 / empty crash)

    // replicates = 0: averaging the mean and taking percentiles would operate over an EMPTY replicate
    // distribution — a DivideByZero / empty-array hazard. The documented guard (§3.3) is a clean
    // ArgumentOutOfRangeException, never a crash, NaN, or a degenerate result.
    [Test]
    public void BootstrapExposures_ZeroReplicates_ThrowsArgumentOutOfRange()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        var catalog = new[] { 5, 3 };

        FluentActions.Invoking(() => BootstrapExposures(catalog, signatures, replicates: 0))
            .Should().Throw<ArgumentOutOfRangeException>()
            .WithMessage("*replicate*",
                "0 replicates would average over an empty distribution ⇒ documented guard, never a ÷0 (§3.3)");
    }

    // Negative replicate counts are equally invalid and must raise the SAME documented guard — including the
    // extreme int.MinValue, which must not be (mis)read as a huge positive loop bound.
    [Test]
    public void BootstrapExposures_NegativeOrExtremeNegativeReplicates_ThrowsArgumentOutOfRange()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        var catalog = new[] { 5, 3 };

        foreach (int reps in new[] { -1, -1000, int.MinValue })
        {
            FluentActions.Invoking(() => BootstrapExposures(catalog, signatures, replicates: reps))
                .Should().Throw<ArgumentOutOfRangeException>(
                    $"replicates = {reps} < 1 is the documented out-of-range boundary (§3.3)");
        }
    }

    #endregion

    #region ONCO-SIG-003 — BE: 1 mutation (total N = 1 ⇒ no ÷0 on tiny total, CI well-defined)

    // A catalog whose total is a SINGLE mutation in one channel: N = 1, p = [..,1,..]. The multinomial draw
    // over N = 1 with one non-zero channel deterministically places that mutation in the channel every time,
    // so every replicate is identical and the per-signature distribution is constant ⇒ lower = upper = mean.
    // The denominator pₖ = dₖ/N = 1/1 must not trigger any DivideByZero / NaN on the tiny total.
    [Test]
    public void BootstrapExposures_SingleMutationSingleChannel_Deterministic_NoDivideByZero()
    {
        // One mutation in channel 0; identity signatures so exposure = the resampled count.
        var catalog = new[] { 1, 0, 0 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0 },
            new double[] { 0, 1, 0 },
            new double[] { 0, 0, 1 },
        };

        var intervals = BootstrapExposures(catalog, signatures, replicates: 100, confidence: 0.95, seed: 13);

        AssertWellFormedIntervals(intervals, 3, 0.95);
        // The single mutation always lands in channel 0 ⇒ exposure of sig0 is constantly 1, the rest 0.
        intervals[0].PointEstimate.Should().BeApproximately(1.0, Eps, "the one mutation maps to sig0 exposure 1");
        intervals[0].Lower.Should().BeApproximately(1.0, Eps, "constant {1} distribution ⇒ lower = 1 (§6.1)");
        intervals[0].Upper.Should().BeApproximately(1.0, Eps, "constant {1} distribution ⇒ upper = 1 (§6.1)");
        intervals[0].Mean.Should().BeApproximately(1.0, Eps);
        for (int j = 1; j < 3; j++)
        {
            intervals[j].PointEstimate.Should().BeApproximately(0.0, Eps, $"sig{j} is unused ⇒ exposure 0");
            intervals[j].Upper.Should().BeApproximately(0.0, Eps, $"sig{j} exposure is constantly 0");
        }
    }

    // A single mutation spread across a multinomial with SEVERAL non-zero channels (each pₖ = dₖ/N where the
    // counts are all 1, so the "total" of the resampled draw is still 1 distributed across channels). The
    // single drawn mutation lands in exactly one channel per replicate, so the distribution is non-trivial
    // but every interval stays well-formed (finite, ≥ 0, lower ≤ upper) with no ÷0 on the tiny total.
    [Test]
    [CancelAfter(20_000)]
    public void BootstrapExposures_SingleMutationMultiChannel_WellFormed_NoDivideByZero()
    {
        // Total N = 1 but the distributing probabilities are spread over 4 equal channels (each count 1 over
        // a catalog whose Σ = 4 is NOT N=1) — so build a true N=1 catalog: one channel = 1, rest 0, but make
        // the signatures NON-identity so the resampled single mutation maps through NNLS to several exposures.
        var catalog = new[] { 1, 0, 0, 0 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };

        var intervals = BootstrapExposures(catalog, signatures, replicates: 150, confidence: 0.95, seed: 99);

        AssertWellFormedIntervals(intervals, 2, 0.95);
        // N = 1 with a single non-zero channel ⇒ every draw is identical ⇒ each interval is a constant
        // (lower = upper = mean = point), and finite — the key BE guarantee is "no ÷0 on total = 1".
        foreach (ExposureConfidenceInterval ci in intervals)
        {
            ci.Lower.Should().BeApproximately(ci.Upper, Eps,
                "a single-channel N=1 catalog gives a constant resample ⇒ degenerate (lower = upper) CI (§6.1)");
            ci.Mean.Should().BeApproximately(ci.PointEstimate, Eps, "constant distribution ⇒ mean = point");
        }
    }

    // Fuzz: random TINY catalogs (total N ∈ {0, 1, 2, 3}) against random signature matrices — every result is
    // well-formed with no DivideByZero / NaN on the tiny total, across many shapes and seeds.
    [Test]
    [CancelAfter(30_000)]
    public void BootstrapExposures_TinyTotalCatalogs_Fuzz_AlwaysWellFormed_NoDivideByZero()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(1, 5);          // 1..4 signatures
            int channels = rng.Next(1, 8);   // 1..7 channels
            var signatures = RandomSignatures(rng, k, channels, scale: rng.NextDouble() * 10.0 + 0.01);

            // Build a catalog with a tiny total N ∈ {0,1,2,3}.
            int targetTotal = rng.Next(0, 4);
            var catalog = new int[channels];
            for (int i = 0; i < targetTotal; i++)
            {
                catalog[rng.Next(channels)]++;
            }

            IReadOnlyList<ExposureConfidenceInterval> intervals = null!;
            FluentActions.Invoking(() =>
                    intervals = BootstrapExposures(catalog, signatures, replicates: 30, confidence: 0.9, seed: 5))
                .Should().NotThrow($"seed {seed}: a tiny-total catalog is a valid degenerate input (no ÷0)");

            AssertWellFormedIntervals(intervals, k, 0.9);
        }
    }

    #endregion

    #region ONCO-SIG-003 — BE: N = 0 (all-zero catalog ⇒ every interval [0,0], INV-05)

    // The all-zero catalog has total N = 0. The multinomial over zero mutations is empty ⇒ every resample is
    // the zero vector ⇒ NNLS(0) = 0 ⇒ every interval is [0,0] with point 0 and mean 0 (INV-05, §6.1). The
    // p = d/N = 0/0 must NOT produce a NaN probability or a DivideByZero.
    [Test]
    public void BootstrapExposures_ZeroCatalog_AllIntervalsZero_NoNaN_NoDivideByZero()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new int[4]; // all zero ⇒ N = 0

        var intervals = BootstrapExposures(catalog, signatures, replicates: 50, confidence: 0.95, seed: 42);

        AssertWellFormedIntervals(intervals, 2, 0.95);
        foreach (ExposureConfidenceInterval ci in intervals)
        {
            ci.PointEstimate.Should().Be(0.0, "NNLS(0) = 0 (INV-05, §6.1)");
            ci.Mean.Should().Be(0.0, "every resample of an empty multinomial is 0 (INV-05)");
            ci.Lower.Should().Be(0.0, "the [0,0] interval lower bound (INV-05)");
            ci.Upper.Should().Be(0.0, "the [0,0] interval upper bound (INV-05)");
        }
    }

    #endregion

    #region ONCO-SIG-003 — BE: fixed-seed extremes (determinism INV-03, extreme reps & seeds)

    // INV-03 headline contract: a fixed seed makes the WHOLE bootstrap reproducible across independent calls
    // — byte-identical intervals — over a fuzzed range of catalog shapes, rep counts, confidences, and seeds.
    [Test]
    [CancelAfter(30_000)]
    public void BootstrapExposures_FixedSeed_FullyReproducible_Fuzz()
    {
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int k = rng.Next(1, 4);
            int channels = rng.Next(1, 8);
            var signatures = RandomSignatures(rng, k, channels, scale: rng.NextDouble() * 20.0 + 0.01);
            var catalog = RandomCatalog(rng, channels, maxCount: rng.Next(1, 40));
            int reps = rng.Next(1, 60);
            double conf = 0.80 + rng.NextDouble() * 0.18; // (0.80, 0.98)
            int bootSeed = rng.Next(int.MinValue, int.MaxValue);

            var first = BootstrapExposures(catalog, signatures, reps, conf, bootSeed);
            var second = BootstrapExposures(catalog, signatures, reps, conf, bootSeed);

            AssertWellFormedIntervals(first, k, conf);
            for (int j = 0; j < k; j++)
            {
                first[j].Should().Be(second[j],
                    $"seed {seed}, sig {j}: identical (d,S,R,c,seed) ⇒ identical intervals (INV-03)");
            }
        }
    }

    // §6.1 single-replicate boundary: replicates = 1 ⇒ the percentile of a 1-element sample is that element,
    // so lower = upper = mean = the single replicate exposure (a finite, non-negative, non-inverted interval).
    [Test]
    public void BootstrapExposures_SingleReplicate_LowerEqualsUpperEqualsMean()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new[] { 30, 12, 10, 8 };

        var intervals = BootstrapExposures(catalog, signatures, replicates: 1, confidence: 0.95, seed: 3);

        AssertWellFormedIntervals(intervals, 2, 0.95);
        foreach (ExposureConfidenceInterval ci in intervals)
        {
            ci.Lower.Should().BeApproximately(ci.Upper, Eps,
                "a 1-element percentile distribution ⇒ lower = upper (§6.1)");
            ci.Mean.Should().BeApproximately(ci.Lower, Eps,
                "the mean of a single replicate equals that replicate ⇒ mean = lower = upper (§6.1)");
        }
    }

    // Extreme / edge seeds (0, -1, int.MinValue, int.MaxValue) must be accepted by Random(seed) with NO
    // overflow / crash, and each must produce a well-formed AND reproducible result.
    [Test]
    [CancelAfter(30_000)]
    public void BootstrapExposures_ExtremeSeeds_NoOverflow_WellFormed_Reproducible()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.5, 0.3, 0.2 },
            new double[] { 0.2, 0.3, 0.5 },
        };
        var catalog = new[] { 25, 15, 20 };

        foreach (int extremeSeed in new[] { 0, -1, int.MinValue, int.MaxValue })
        {
            IReadOnlyList<ExposureConfidenceInterval> first = null!;
            IReadOnlyList<ExposureConfidenceInterval> second = null!;

            FluentActions.Invoking(() =>
            {
                first = BootstrapExposures(catalog, signatures, replicates: 80, confidence: 0.95, seed: extremeSeed);
                second = BootstrapExposures(catalog, signatures, replicates: 80, confidence: 0.95, seed: extremeSeed);
            }).Should().NotThrow($"seed {extremeSeed}: an extreme seed must not overflow or crash Random(seed)");

            AssertWellFormedIntervals(first, 2, 0.95);
            for (int j = 0; j < first.Count; j++)
            {
                first[j].Should().Be(second[j],
                    $"seed {extremeSeed}, sig {j}: extreme seeds are still fully reproducible (INV-03)");
            }
        }
    }

    // Extreme replicate counts: replicates = 1 (minimal) up through a modest-but-large count must terminate
    // (no hang — the per-replicate NNLS loop is bounded) and stay well-formed. Kept modest under CancelAfter.
    [Test]
    [CancelAfter(60_000)]
    public void BootstrapExposures_ExtremeReplicateCounts_Terminate_WellFormed()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new[] { 35, 14, 11, 10 };

        foreach (int reps in new[] { 1, 2, 5, 50, 500, 2000 })
        {
            IReadOnlyList<ExposureConfidenceInterval> intervals = null!;
            FluentActions.Invoking(() =>
                    intervals = BootstrapExposures(catalog, signatures, replicates: reps, confidence: 0.95, seed: 42))
                .Should().NotThrow($"replicates = {reps}: the bootstrap loop must terminate, not hang");

            AssertWellFormedIntervals(intervals, 2, 0.95);
        }
    }

    // Extreme confidence levels near the open-interval boundary (very wide / very narrow) must still produce
    // well-formed, non-inverted intervals, while exactly 0 / 1 / NaN confidence is the documented out-of-range
    // guard (§3.3).
    [Test]
    public void BootstrapExposures_ExtremeConfidenceLevels_BoundaryAndGuard()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 0.6, 0.2, 0.1, 0.1 },
            new double[] { 0.1, 0.1, 0.3, 0.5 },
        };
        var catalog = new[] { 30, 12, 10, 8 };

        // Near-boundary VALID confidences must work and stay well-formed.
        foreach (double conf in new[] { 0.001, 0.5, 0.999 })
        {
            var intervals = BootstrapExposures(catalog, signatures, replicates: 60, confidence: conf, seed: 42);
            AssertWellFormedIntervals(intervals, 2, conf);
        }

        // Out-of-range / NaN confidences are the documented guard (§3.3).
        foreach (double bad in new[] { 0.0, 1.0, -0.1, 1.1, double.NaN })
        {
            FluentActions.Invoking(() =>
                    BootstrapExposures(catalog, signatures, replicates: 10, confidence: bad, seed: 42))
                .Should().Throw<ArgumentOutOfRangeException>(
                    $"confidence = {bad} ∉ (0,1) is the documented out-of-range guard (§3.3)");
        }
    }

    #endregion

    #region ONCO-SIG-003 — BE: structural boundaries (null / empty / ragged / mismatch / negative ⇒ throws)

    // null catalog / null signatures ⇒ ArgumentNullException; a null signature VECTOR is caught by the shared
    // ValidateSignatures helper as an ArgumentException ("Signature vectors cannot be null.") — consistent
    // with ONCO-SIG-002's FitSignatures. (ArgumentNullException ⊂ ArgumentException, so both are within the
    // documented §3.3 Argument* guard family; the null-vector case lands on the base ArgumentException.)
    [Test]
    public void BootstrapExposures_NullArguments_ThrowArgumentExceptionFamily()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };

        FluentActions.Invoking(() => BootstrapExposures(null!, signatures))
            .Should().Throw<ArgumentNullException>("a null catalog is a documented guard (§3.3)");
        FluentActions.Invoking(() => BootstrapExposures(new[] { 1, 2 }, null!))
            .Should().Throw<ArgumentNullException>("a null signature matrix is a documented guard (§3.3)");

        var withNull = new IReadOnlyList<double>[] { new double[] { 1, 0 }, null! };
        FluentActions.Invoking(() => BootstrapExposures(new[] { 1, 2 }, withNull))
            .Should().Throw<ArgumentException>("a null signature vector is a documented Argument* guard (§3.3)");
    }

    // Empty signature matrix / empty (zero-channel) signatures / ragged signatures ⇒ ArgumentException (§3.3).
    [Test]
    public void BootstrapExposures_EmptyOrRaggedSignatures_ThrowArgumentException()
    {
        FluentActions.Invoking(() => BootstrapExposures(new[] { 1, 2 }, Array.Empty<IReadOnlyList<double>>()))
            .Should().Throw<ArgumentException>("at least one reference signature is required (§3.3)");

        var emptyVectors = new IReadOnlyList<double>[] { Array.Empty<double>() };
        FluentActions.Invoking(() => BootstrapExposures(Array.Empty<int>(), emptyVectors))
            .Should().Throw<ArgumentException>("signature vectors cannot be empty (§3.3)");

        var ragged = new IReadOnlyList<double>[] { new double[] { 1, 0, 0 }, new double[] { 0, 1 } };
        FluentActions.Invoking(() => BootstrapExposures(new[] { 1, 2, 3 }, ragged))
            .Should().Throw<ArgumentException>("all signatures must have the same length (§3.3)");
    }

    // Catalog length ≠ signature channel count ⇒ ArgumentException (§3.3) — the dimension-mismatch boundary.
    [Test]
    public void BootstrapExposures_CatalogLengthMismatch_ThrowsArgumentException()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        FluentActions.Invoking(() => BootstrapExposures(new[] { 1, 2, 3 }, signatures))
            .Should().Throw<ArgumentException>("catalog length must equal the signature channel count (§3.3)");
    }

    // A NEGATIVE catalog count is malformed ⇒ ArgumentException (§3.3) — the bootstrap requires non-negative
    // integer counts (the multinomial sample size N = Σ catalog is otherwise ill-defined).
    [Test]
    public void BootstrapExposures_NegativeCount_ThrowsArgumentException()
    {
        var signatures = new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };
        FluentActions.Invoking(() => BootstrapExposures(new[] { 5, -3 }, signatures))
            .Should().Throw<ArgumentException>("a negative mutation count is malformed (§3.3)");
    }

    // Fuzz: random structural corruptions (length mismatch / ragged / null vector / negative count) always
    // raise the documented Argument* family — never an undisciplined IndexOutOfRange / NullReference crash.
    [Test]
    [CancelAfter(20_000)]
    public void BootstrapExposures_RandomMalformedStructure_AlwaysDocumentedArgumentException()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int channels = rng.Next(1, 8);
            int k = rng.Next(1, 5);
            var signatures = RandomSignatures(rng, k, channels);

            int defect = rng.Next(4);
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

                    act = () => BootstrapExposures(new int[wrongLen], signatures);
                    break;

                case 1: // ragged signatures
                    var ragged = signatures.ToArray();
                    ragged[rng.Next(k)] = new double[channels + 1];
                    act = () => BootstrapExposures(new int[channels], ragged);
                    break;

                case 2: // null signature vector
                    var withNull = signatures.ToArray();
                    withNull[rng.Next(k)] = null!;
                    act = () => BootstrapExposures(new int[channels], withNull);
                    break;

                default: // negative count
                    var negCatalog = new int[channels];
                    negCatalog[rng.Next(channels)] = -(rng.Next(1, 100));
                    act = () => BootstrapExposures(negCatalog, signatures);
                    break;
            }

            FluentActions.Invoking(act).Should().Throw<ArgumentException>(
                $"seed {seed}: malformed structure (defect {defect}) ⇒ documented ArgumentException, not a crash");
        }
    }

    #endregion
}
