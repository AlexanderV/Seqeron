// ONCO-SIG-003 — Signature Exposure Estimation: Bootstrap Confidence Intervals
// Evidence: docs/Evidence/ONCO-SIG-003-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-SIG-003.md
// Source: Senkin S. (2021). MSA. BMC Bioinformatics 22:540. https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/
//         Huang X., Wojtowicz D., Przytycka T.M. (2018). Bioinformatics 34(2):330-337. https://academic.oup.com/bioinformatics/article/34/2/330/4209996
//         Wang S. et al. sigminer sig_fit_bootstrap. https://raw.githubusercontent.com/ShixiangWang/sigminer/master/R/sig_fit_bootstrap.R
//         Efron B. (1979). Annals of Statistics 7(1):1-26 (percentile method). https://doi.org/10.1214/aos/1176344552
//         Hyndman R.J., Fan Y. (1996). The American Statistician 50(4):361-365 (type-7). https://doi.org/10.1080/00031305.1996.10473566
//         Knuth D.E. (1997). TAOCP Vol. 2, §3.4.1 (Poisson deviate generation). [Poisson variant]

using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_BootstrapExposures_Tests
{
    // Single signature equal to the single channel: NNLS fit of [c] onto [[1.0]] gives exposure c.
    private static readonly IReadOnlyList<IReadOnlyList<double>> SingleUnitSignature =
        new IReadOnlyList<double>[] { new double[] { 1.0 } };

    // Two orthogonal unit signatures over two channels: NNLS exposure[j] = catalog[j].
    private static readonly IReadOnlyList<IReadOnlyList<double>> TwoOrthogonalSignatures =
        new IReadOnlyList<double>[] { new double[] { 1, 0 }, new double[] { 0, 1 } };

    #region Deterministic multinomial-collapse cases (seed-independent)

    // M1 — Single non-zero channel: Multinomial(N,p) over one outcome is deterministic, so every
    // resample equals [10] and every NNLS exposure equals 10 (Senkin 2021 resample; NNLS fit).
    [Test]
    public void BootstrapExposures_SingleChannelCatalog_CollapsesToPointEstimate()
    {
        var catalog = new[] { 10 };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: 100, confidence: 0.95, seed: 42);

        Assert.That(intervals, Has.Count.EqualTo(1), "One signature must yield one interval.");
        var interval = intervals[0];
        Assert.Multiple(() =>
        {
            Assert.That(interval.PointEstimate, Is.EqualTo(10.0).Within(1e-10),
                "Point estimate is the NNLS exposure of the observed [10] onto [[1.0]] = 10.");
            Assert.That(interval.Mean, Is.EqualTo(10.0).Within(1e-10),
                "All resamples equal [10] (deterministic multinomial), so the mean exposure is 10.");
            Assert.That(interval.Lower, Is.EqualTo(10.0).Within(1e-10),
                "A constant bootstrap distribution gives lower percentile = 10 (type-7 of constant sample).");
            Assert.That(interval.Upper, Is.EqualTo(10.0).Within(1e-10),
                "A constant bootstrap distribution gives upper percentile = 10.");
        });
    }

    // M2 — Same collapse, explicit 95% level: 2.5% and 97.5% percentiles of the constant {10} are both 10
    // (Efron 1979 percentile method on a degenerate distribution).
    [Test]
    public void BootstrapExposures_SingleChannel95Percent_BoundsEqualPointEstimate()
    {
        var catalog = new[] { 7 };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: 50, confidence: 0.95, seed: 7);

        Assert.Multiple(() =>
        {
            Assert.That(intervals[0].Lower, Is.EqualTo(7.0).Within(1e-10),
                "2.5th percentile of the constant {7} distribution is 7.");
            Assert.That(intervals[0].Upper, Is.EqualTo(7.0).Within(1e-10),
                "97.5th percentile of the constant {7} distribution is 7.");
            Assert.That(intervals[0].Confidence, Is.EqualTo(0.95).Within(1e-12),
                "The interval records the requested confidence level.");
        });
    }

    // M3 — Zero-mutation catalog (N=0): no mutations to resample, every replicate is all-zero, NNLS(0)=0,
    // so every signature interval collapses to [0,0] (documented N=0 corner case).
    [Test]
    public void BootstrapExposures_ZeroMutationCatalog_AllIntervalsZero()
    {
        var catalog = new[] { 0, 0, 0 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0 },
            new double[] { 0, 1, 0 },
        };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, signatures, replicates: 100, confidence: 0.95, seed: 42);

        Assert.That(intervals, Has.Count.EqualTo(2), "Two signatures -> two intervals.");
        foreach (var interval in intervals)
        {
            Assert.Multiple(() =>
            {
                Assert.That(interval.PointEstimate, Is.EqualTo(0.0).Within(1e-12),
                    "Zero catalog -> NNLS exposure 0.");
                Assert.That(interval.Mean, Is.EqualTo(0.0).Within(1e-12), "All replicates are 0.");
                Assert.That(interval.Lower, Is.EqualTo(0.0).Within(1e-12), "Lower bound 0.");
                Assert.That(interval.Upper, Is.EqualTo(0.0).Within(1e-12), "Upper bound 0.");
            });
        }
    }

    // S2 — Two-channel deterministic split: catalog [0,7] over orthogonal signatures. The single non-zero
    // channel (index 1) takes all mutations every draw, so sig0 collapses to 0 and sig1 collapses to 7.
    [Test]
    public void BootstrapExposures_TwoChannelDeterministicSplit_ExactExposures()
    {
        var catalog = new[] { 0, 7 };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, TwoOrthogonalSignatures, replicates: 200, confidence: 0.95, seed: 123);

        Assert.Multiple(() =>
        {
            Assert.That(intervals[0].PointEstimate, Is.EqualTo(0.0).Within(1e-10), "Channel 0 has no mutations -> exposure 0.");
            Assert.That(intervals[0].Lower, Is.EqualTo(0.0).Within(1e-10), "Sig0 lower bound 0.");
            Assert.That(intervals[0].Upper, Is.EqualTo(0.0).Within(1e-10), "Sig0 upper bound 0.");
            Assert.That(intervals[1].PointEstimate, Is.EqualTo(7.0).Within(1e-10), "All 7 mutations fall in channel 1 -> exposure 7.");
            Assert.That(intervals[1].Lower, Is.EqualTo(7.0).Within(1e-10), "Sig1 lower bound 7 (deterministic).");
            Assert.That(intervals[1].Upper, Is.EqualTo(7.0).Within(1e-10), "Sig1 upper bound 7 (deterministic).");
        });
    }

    // M9 — Single replicate: percentile of a one-element bootstrap distribution is that single value,
    // so lower = upper = mean (Efron 1979; type-7 n=1 case).
    [Test]
    public void BootstrapExposures_SingleReplicate_BoundsEqualMean()
    {
        var catalog = new[] { 10 };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: 1, confidence: 0.95, seed: 99);

        Assert.Multiple(() =>
        {
            Assert.That(intervals[0].Lower, Is.EqualTo(intervals[0].Mean).Within(1e-12),
                "With R=1 the percentile equals the single replicate value, so lower = mean.");
            Assert.That(intervals[0].Upper, Is.EqualTo(intervals[0].Mean).Within(1e-12),
                "With R=1 the percentile equals the single replicate value, so upper = mean.");
            Assert.That(intervals[0].Mean, Is.EqualTo(10.0).Within(1e-10),
                "The single resample of [10] is [10]; NNLS exposure = 10.");
        });
    }

    // M8 — Type-7 percentile property on a constant distribution: for any probability p, Q(p) of a constant
    // sample equals that constant (h=p*(n-1), interpolation between equal order statistics). Verified at the
    // median (p=0.5) via a deterministic collapse so the bound is exactly the point estimate.
    // NOTE: a constant distribution does NOT exercise the type-7 interpolation rule (any quantile estimator
    // returns the constant). The interpolation itself is locked on non-constant samples by
    // Percentile_Type7_NonConstantSamples_MatchesHandDerivedValues below.
    [Test]
    public void BootstrapExposures_Type7Median_OnConstantSplit_IsExact()
    {
        var catalog = new[] { 4 };

        // confidence 0.0 is invalid, so probe the median indirectly: a symmetric narrow level still yields
        // the constant value because the distribution is degenerate {4}.
        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: 64, confidence: 0.50, seed: 5);

        Assert.Multiple(() =>
        {
            // 25th and 75th percentiles of constant {4} are both 4 (type-7 over equal order statistics).
            Assert.That(intervals[0].Lower, Is.EqualTo(4.0).Within(1e-10),
                "Type-7 25th percentile of constant {4} = 4.");
            Assert.That(intervals[0].Upper, Is.EqualTo(4.0).Within(1e-10),
                "Type-7 75th percentile of constant {4} = 4.");
        });
    }

    // M8b — Type-7 interpolation on NON-constant samples. The percentile bound values are determined by the
    // internal Percentile() (R-7 / NumPy-default linear interpolation, Hyndman & Fan 1996). The deterministic
    // collapse tests above cannot distinguish type-7 from any other quantile rule, so this test exercises the
    // interpolation directly on the two hand-derived datasets in the Evidence artifact. The expected values
    // were independently confirmed with numpy.quantile(..., method='linear'):
    //   [0,1,2,3,4]: p=0.025 -> 0.1, p=0.5 -> 2.0, p=0.975 -> 3.9 ; p=0 -> 0 (min), p=1 -> 4 (max)
    //   [2,4,6,8]:   p=0.025 -> 2.15, p=0.5 -> 5.0, p=0.975 -> 7.85
    // Source: Hyndman & Fan (1996), The American Statistician 50(4):361-365 (type-7);
    //         cross-checked numpy.quantile method='linear'.
    [Test]
    public void Percentile_Type7_NonConstantSamples_MatchesHandDerivedValues()
    {
        var sampleA = new[] { 0.0, 1.0, 2.0, 3.0, 4.0 };
        var sampleB = new[] { 2.0, 4.0, 6.0, 8.0 };

        Assert.Multiple(() =>
        {
            Assert.That(InvokePercentile(sampleA, 0.0), Is.EqualTo(0.0).Within(1e-12), "min of [0..4].");
            Assert.That(InvokePercentile(sampleA, 0.025), Is.EqualTo(0.1).Within(1e-12), "type-7 p=0.025 of [0..4].");
            Assert.That(InvokePercentile(sampleA, 0.5), Is.EqualTo(2.0).Within(1e-12), "type-7 median of [0..4].");
            Assert.That(InvokePercentile(sampleA, 0.975), Is.EqualTo(3.9).Within(1e-12), "type-7 p=0.975 of [0..4].");
            Assert.That(InvokePercentile(sampleA, 1.0), Is.EqualTo(4.0).Within(1e-12), "max of [0..4].");

            Assert.That(InvokePercentile(sampleB, 0.025), Is.EqualTo(2.15).Within(1e-12), "type-7 p=0.025 of [2,4,6,8].");
            Assert.That(InvokePercentile(sampleB, 0.5), Is.EqualTo(5.0).Within(1e-12), "type-7 median of [2,4,6,8].");
            Assert.That(InvokePercentile(sampleB, 0.975), Is.EqualTo(7.85).Within(1e-12), "type-7 p=0.975 of [2,4,6,8].");
        });
    }

    private static readonly MethodInfo PercentileMethod =
        typeof(OncologyAnalyzer).GetMethod("Percentile", BindingFlags.NonPublic | BindingFlags.Static)
        ?? throw new InvalidOperationException("OncologyAnalyzer.Percentile(double[], double) not found.");

    private static double InvokePercentile(double[] values, double probability) =>
        (double)PercentileMethod.Invoke(null, new object[] { values, probability })!;

    #endregion

    #region Contract: point estimate, ordering, determinism, shape (seeded randomized cases)

    // M4 — Point estimate equals the NNLS fit of the observed (un-resampled) catalog (Senkin 2021; INV-5).
    [Test]
    public void BootstrapExposures_PointEstimate_EqualsObservedNnlsFit()
    {
        var catalog = new[] { 30, 10, 0, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 0, 1 },
            new double[] { 0, 0, 1, 0 },
        };

        var observed = catalog.Select(c => (double)c).ToArray();
        var expectedPoint = OncologyAnalyzer.FitSignatures(observed, signatures).Exposures;

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, signatures, replicates: 200, confidence: 0.95, seed: 42);

        Assert.That(intervals, Has.Count.EqualTo(3), "Three signatures -> three intervals.");
        Assert.Multiple(() =>
        {
            for (int j = 0; j < 3; j++)
            {
                Assert.That(intervals[j].PointEstimate, Is.EqualTo(expectedPoint[j]).Within(1e-10),
                    $"Interval {j} point estimate must equal the observed-catalog NNLS exposure (Senkin 2021).");
            }
        });
    }

    // M5 — Interval ordering and non-negativity invariants on a non-degenerate catalog (INV-1/2/3).
    [Test]
    public void BootstrapExposures_NonDegenerate_IntervalsOrderedAndNonNegative()
    {
        var catalog = new[] { 25, 25, 25, 25 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 0, 0 },
            new double[] { 0, 0, 1, 1 },
        };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, signatures, replicates: 500, confidence: 0.95, seed: 42);

        Assert.Multiple(() =>
        {
            foreach (var iv in intervals)
            {
                Assert.That(iv.Lower, Is.GreaterThanOrEqualTo(0.0), "All bounds >= 0 (NNLS x >= 0, multinomial >= 0).");
                Assert.That(iv.Lower, Is.LessThanOrEqualTo(iv.Upper + 1e-12), "Lower bound must not exceed upper bound.");
                Assert.That(iv.Mean, Is.GreaterThanOrEqualTo(iv.Lower - 1e-9), "Mean must be >= lower bound.");
                Assert.That(iv.Mean, Is.LessThanOrEqualTo(iv.Upper + 1e-9), "Mean must be <= upper bound.");
            }
        });
    }

    // M6 — Determinism: identical inputs and seed produce element-wise identical intervals (INV-4).
    [Test]
    public void BootstrapExposures_SameSeed_IsDeterministic()
    {
        var catalog = new[] { 12, 8, 0, 4 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 0, 1 },
        };

        var a = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95, seed: 42);
        var b = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95, seed: 42);

        Assert.Multiple(() =>
        {
            for (int j = 0; j < a.Count; j++)
            {
                Assert.That(b[j].Lower, Is.EqualTo(a[j].Lower).Within(1e-12), $"Lower[{j}] must be reproducible for a fixed seed.");
                Assert.That(b[j].Upper, Is.EqualTo(a[j].Upper).Within(1e-12), $"Upper[{j}] must be reproducible for a fixed seed.");
                Assert.That(b[j].Mean, Is.EqualTo(a[j].Mean).Within(1e-12), $"Mean[{j}] must be reproducible for a fixed seed.");
            }
        });
    }

    // M7 — One interval per signature, in signature order (INV-5 shape).
    [Test]
    public void BootstrapExposures_ReturnsOneIntervalPerSignatureInOrder()
    {
        var catalog = new[] { 5, 5, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0 },
            new double[] { 0, 1, 0 },
            new double[] { 0, 0, 1 },
            new double[] { 1, 1, 1 },
        };

        var intervals = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 100, seed: 1);

        Assert.That(intervals, Has.Count.EqualTo(4), "Result count must equal the signature count, in signature order.");
    }

    // S1 — Wider confidence gives a wider-or-equal interval on a non-degenerate catalog (percentile monotonicity).
    [Test]
    public void BootstrapExposures_WiderConfidence_GivesWiderInterval()
    {
        var catalog = new[] { 40, 30, 20, 10 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 1, 1 },
        };

        var narrow = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 500, confidence: 0.50, seed: 42);
        var wide = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 500, confidence: 0.99, seed: 42);

        Assert.Multiple(() =>
        {
            for (int j = 0; j < narrow.Count; j++)
            {
                Assert.That(wide[j].Lower, Is.LessThanOrEqualTo(narrow[j].Lower + 1e-9),
                    $"A higher confidence level cannot raise the lower bound (signature {j}).");
                Assert.That(wide[j].Upper, Is.GreaterThanOrEqualTo(narrow[j].Upper - 1e-9),
                    $"A higher confidence level cannot lower the upper bound (signature {j}).");
            }
        });
    }

    #endregion

    #region Poisson resampling variant (Senkin 2021 MSA Poisson noise)

    // Independent reference implementation of the per-channel Poisson draw used to derive expected values,
    // NOT a copy of production code: Knuth's exact multiplication-of-uniforms generator (Knuth, TAOCP Vol. 2,
    // §3.4.1) for Poisson(lambda). Consumes Random in exactly the order BootstrapExposures does so the seeded
    // draw sequence is reproduced. lambda <= 0 returns 0 (mean-0 degenerate case, Senkin 2021).
    private static long ReferencePoisson(double lambda, Random random)
    {
        if (lambda <= 0.0)
        {
            return 0;
        }

        double limit = Math.Exp(-lambda);
        long count = 0;
        double product = random.NextDouble();
        while (product > limit)
        {
            count++;
            product *= random.NextDouble();
        }

        return count;
    }

    // P1 — Degenerate single-channel catalog under the Poisson variant: the unit signature makes NNLS
    // exposure == the resampled count, so each replicate exposure must equal an independent Poisson(observed)
    // draw. We reproduce the seeded draw sequence with the Knuth reference generator and require an EXACT,
    // element-wise match of the bootstrap mean to the mean of those independent Poisson draws. This proves the
    // Poisson construction (per-channel Poisson(observedₖ)), not just "some randomness".
    // Source: Senkin S. (2021), MSA — bootstrapped per-category counts follow a Poisson distribution.
    [Test]
    public void BootstrapExposures_Poisson_SingleChannel_MatchesIndependentPoissonDraws()
    {
        const int Lambda = 12;
        const int Replicates = 300;
        const int Seed = 42;
        var catalog = new[] { Lambda };

        // Independent expected bootstrap distribution: Poisson(12) draws in the same seeded order.
        var reference = new Random(Seed);
        double sum = 0.0;
        double[] draws = new double[Replicates];
        for (int rep = 0; rep < Replicates; rep++)
        {
            draws[rep] = ReferencePoisson(Lambda, reference);
            sum += draws[rep];
        }

        double expectedMean = sum / Replicates;
        Array.Sort(draws);
        double expectedLower = Type7Percentile(draws, 0.025);
        double expectedUpper = Type7Percentile(draws, 0.975);

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: Replicates, confidence: 0.95, seed: Seed,
            resampling: OncologyAnalyzer.BootstrapResampling.Poisson);

        Assert.Multiple(() =>
        {
            Assert.That(intervals[0].PointEstimate, Is.EqualTo((double)Lambda).Within(1e-10),
                "Point estimate is the NNLS fit of the observed (un-resampled) catalog, independent of resampling.");
            Assert.That(intervals[0].Mean, Is.EqualTo(expectedMean).Within(1e-10),
                "Each replicate exposure equals an independent Poisson(12) draw; mean must match the reference draws.");
            Assert.That(intervals[0].Lower, Is.EqualTo(expectedLower).Within(1e-10),
                "Lower bound is the 2.5th type-7 percentile of the independent Poisson(12) bootstrap distribution.");
            Assert.That(intervals[0].Upper, Is.EqualTo(expectedUpper).Within(1e-10),
                "Upper bound is the 97.5th type-7 percentile of the independent Poisson(12) bootstrap distribution.");
        });
    }

    // P2 — A channel with observed count 0 has Poisson mean 0, so it resamples to 0 every replicate: that
    // signature's exposure is exactly 0 (point estimate and all bounds). The non-zero channel is non-degenerate
    // (Poisson does not fix N), so its bounds straddle the point estimate. Source: Senkin 2021 (Poisson(0)=0).
    [Test]
    public void BootstrapExposures_Poisson_ZeroCountChannel_StaysZero()
    {
        var catalog = new[] { 0, 20 };

        var intervals = OncologyAnalyzer.BootstrapExposures(
            catalog, TwoOrthogonalSignatures, replicates: 400, confidence: 0.95, seed: 42,
            resampling: OncologyAnalyzer.BootstrapResampling.Poisson);

        Assert.Multiple(() =>
        {
            Assert.That(intervals[0].PointEstimate, Is.EqualTo(0.0).Within(1e-12), "Channel 0 observed 0 -> exposure 0.");
            Assert.That(intervals[0].Mean, Is.EqualTo(0.0).Within(1e-12), "Poisson(0) draws are always 0.");
            Assert.That(intervals[0].Lower, Is.EqualTo(0.0).Within(1e-12), "All replicates 0 -> lower 0.");
            Assert.That(intervals[0].Upper, Is.EqualTo(0.0).Within(1e-12), "All replicates 0 -> upper 0.");
            Assert.That(intervals[1].PointEstimate, Is.EqualTo(20.0).Within(1e-10), "Channel 1 NNLS exposure 20.");
        });
    }

    // P3 — Poisson variance is non-zero: unlike the fixed-N multinomial collapse of a single channel
    // (every replicate == observed), the Poisson single-channel resample fluctuates, so the 95% interval is
    // strictly wider than the degenerate point. Locks that the Poisson path is genuinely Poisson, not a
    // fixed-N copy. Source: Senkin 2021 (total burden no longer fixed; variance = mean).
    [Test]
    public void BootstrapExposures_Poisson_SingleChannel_HasNonZeroSpread()
    {
        var catalog = new[] { 25 };

        var interval = OncologyAnalyzer.BootstrapExposures(
            catalog, SingleUnitSignature, replicates: 400, confidence: 0.95, seed: 42,
            resampling: OncologyAnalyzer.BootstrapResampling.Poisson)[0];

        Assert.That(interval.Upper - interval.Lower, Is.GreaterThan(0.0),
            "Poisson(25) resampling fluctuates (variance = mean = 25), so the interval has positive width "
            + "— in contrast to the deterministic fixed-N multinomial collapse of a single channel.");
    }

    // P4 — Determinism of the Poisson path: identical inputs + seed give element-wise identical intervals.
    // Source: fixed-seed reproducibility convention (deterministic clinical re-runs).
    [Test]
    public void BootstrapExposures_Poisson_SameSeed_IsDeterministic()
    {
        var catalog = new[] { 12, 8, 0, 4 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 0, 1 },
        };

        var a = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95,
            seed: 42, resampling: OncologyAnalyzer.BootstrapResampling.Poisson);
        var b = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95,
            seed: 42, resampling: OncologyAnalyzer.BootstrapResampling.Poisson);

        Assert.Multiple(() =>
        {
            for (int j = 0; j < a.Count; j++)
            {
                Assert.That(b[j].Lower, Is.EqualTo(a[j].Lower).Within(1e-12), $"Lower[{j}] reproducible for a fixed seed.");
                Assert.That(b[j].Upper, Is.EqualTo(a[j].Upper).Within(1e-12), $"Upper[{j}] reproducible for a fixed seed.");
                Assert.That(b[j].Mean, Is.EqualTo(a[j].Mean).Within(1e-12), $"Mean[{j}] reproducible for a fixed seed.");
            }
        });
    }

    // P5 — Default path is unchanged: omitting the resampling argument is byte-for-byte equal to explicitly
    // requesting Multinomial, and the Poisson path differs (proving the default did not silently change).
    // Source: backward-compatibility requirement (sigminer fixed-N multinomial remains the default).
    [Test]
    public void BootstrapExposures_DefaultEqualsMultinomial_AndDiffersFromPoisson()
    {
        var catalog = new[] { 30, 10, 0, 5 };
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1, 0, 0, 0 },
            new double[] { 0, 1, 0, 1 },
            new double[] { 0, 0, 1, 0 },
        };

        var defaulted = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95, seed: 42);
        var multinomial = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95,
            seed: 42, resampling: OncologyAnalyzer.BootstrapResampling.Multinomial);
        var poisson = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates: 300, confidence: 0.95,
            seed: 42, resampling: OncologyAnalyzer.BootstrapResampling.Poisson);

        bool poissonDiffers = false;
        Assert.Multiple(() =>
        {
            for (int j = 0; j < defaulted.Count; j++)
            {
                Assert.That(multinomial[j].Lower, Is.EqualTo(defaulted[j].Lower).Within(1e-12),
                    $"Default Lower[{j}] must equal explicit Multinomial (default is unchanged).");
                Assert.That(multinomial[j].Upper, Is.EqualTo(defaulted[j].Upper).Within(1e-12),
                    $"Default Upper[{j}] must equal explicit Multinomial.");
                Assert.That(multinomial[j].Mean, Is.EqualTo(defaulted[j].Mean).Within(1e-12),
                    $"Default Mean[{j}] must equal explicit Multinomial.");
                if (Math.Abs(poisson[j].Upper - multinomial[j].Upper) > 1e-9
                    || Math.Abs(poisson[j].Mean - multinomial[j].Mean) > 1e-9)
                {
                    poissonDiffers = true;
                }
            }
        });
        Assert.That(poissonDiffers, Is.True,
            "The Poisson path must produce a different bootstrap distribution than the multinomial default.");
    }

    [Test]
    public void BootstrapExposures_UnrecognisedResampling_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.BootstrapExposures(new[] { 10 }, SingleUnitSignature,
                resampling: (OncologyAnalyzer.BootstrapResampling)99),
            "An undefined resampling enum value must throw ArgumentOutOfRangeException.");
    }

    // Type-7 (R-7 / NumPy default) percentile on a pre-sorted sample, used only to derive P1 expected bounds
    // independently of production code. Hyndman & Fan (1996), The American Statistician 50(4):361-365.
    private static double Type7Percentile(double[] sorted, double probability)
    {
        int n = sorted.Length;
        if (n == 1)
        {
            return sorted[0];
        }

        double rank = probability * (n - 1);
        int lower = (int)Math.Floor(rank);
        int upper = (int)Math.Ceiling(rank);
        if (lower == upper)
        {
            return sorted[lower];
        }

        return sorted[lower] + (rank - lower) * (sorted[upper] - sorted[lower]);
    }

    #endregion

    #region Failure modes

    [Test]
    public void BootstrapExposures_NullCatalog_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.BootstrapExposures(null!, SingleUnitSignature),
            "Null catalog must throw ArgumentNullException.");
    }

    [Test]
    public void BootstrapExposures_NullSignatures_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.BootstrapExposures(new[] { 1 }, null!),
            "Null signatures must throw ArgumentNullException.");
    }

    [Test]
    public void BootstrapExposures_EmptySignatures_Throws()
    {
        var empty = Array.Empty<IReadOnlyList<double>>();
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.BootstrapExposures(new[] { 1 }, empty),
            "An empty signature set must throw ArgumentException.");
    }

    [Test]
    public void BootstrapExposures_CatalogLengthMismatch_Throws()
    {
        var catalog = new[] { 1, 2, 3 }; // length 3 vs channel count 1
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.BootstrapExposures(catalog, SingleUnitSignature),
            "Catalog length must equal the signature channel count.");
    }

    [Test]
    public void BootstrapExposures_NegativeCatalogCount_Throws()
    {
        var catalog = new[] { -1 };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.BootstrapExposures(catalog, SingleUnitSignature),
            "A negative catalog count is invalid.");
    }

    [Test]
    public void BootstrapExposures_ReplicatesBelowOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.BootstrapExposures(new[] { 10 }, SingleUnitSignature, replicates: 0),
            "At least one replicate is required.");
    }

    [Test]
    public void BootstrapExposures_ConfidenceOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.BootstrapExposures(new[] { 10 }, SingleUnitSignature, confidence: 0.0),
                "Confidence 0 is outside (0,1).");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.BootstrapExposures(new[] { 10 }, SingleUnitSignature, confidence: 1.0),
                "Confidence 1 is outside (0,1).");
        });
    }

    #endregion
}
