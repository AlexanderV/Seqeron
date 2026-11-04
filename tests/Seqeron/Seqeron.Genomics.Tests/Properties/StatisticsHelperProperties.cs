using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the canonical statistics helpers added in the 2026-09 review (B24 ONCO-HETERO-001):
/// <see cref="StatisticsHelper.Median"/> (R <c>median.default</c> / <c>numpy.median</c>) and
/// <see cref="StatisticsHelper.ShannonIndex"/> (Shannon 1948, natural log; = <c>scipy.stats.entropy(counts)</c>).
/// Oracles are independent restatements of the definitions.
/// </summary>
[TestFixture]
[Category("Property")]
public class StatisticsHelperProperties
{
    private static Arbitrary<(double[] Values, int Seed)> ValuesArbitrary() =>
        (from n in Gen.Choose(1, 40)
         from raw in Gen.Choose(-100_000, 100_000).ArrayOf(n)
         from seed in Gen.Choose(0, int.MaxValue)
         select (raw.Select(v => v / 997.0).ToArray(), seed)).ToArbitrary();

    /// <summary>
    /// Definition oracle (R <c>median.default</c>): the central order statistic for odd n, the mean of the two central
    /// order statistics for even n; the result lies in [min, max] and the input is not mutated.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Median_EqualsCentralOrderStatistic_WithinRange_InputUntouched()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
        {
            double[] copy = (double[])t.Values.Clone();
            double median = StatisticsHelper.Median(t.Values);

            double[] sorted = t.Values.OrderBy(v => v).ToArray();
            int n = sorted.Length;
            double oracle = n % 2 == 1 ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;

            return (median == oracle && median >= sorted[0] && median <= sorted[^1] && copy.SequenceEqual(t.Values))
                .Label($"median {median} vs oracle {oracle}");
        });
    }

    /// <summary>The median is a function of the multiset only: any permutation gives a bit-identical value.</summary>
    [FsCheck.NUnit.Property]
    public Property Median_IsPermutationInvariant()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
        {
            var rng = new Random(t.Seed);
            double[] shuffled = t.Values.OrderBy(_ => rng.Next()).ToArray();
            return (StatisticsHelper.Median(shuffled) == StatisticsHelper.Median(t.Values)).Label("permutation changed the median");
        });
    }

    /// <summary>Median(−x) = −Median(x) exactly (order statistics reverse; negation is exact in IEEE arithmetic).</summary>
    [FsCheck.NUnit.Property]
    public Property Median_IsOddUnderNegation()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
            (StatisticsHelper.Median(t.Values.Select(v => -v).ToArray()) == -StatisticsHelper.Median(t.Values))
                .Label("Median(−x) ≠ −Median(x)"));
    }

    /// <summary>A NaN anywhere makes the median NaN (R: <c>median(c(1, NaN))</c> is NA); empty input is rejected.</summary>
    [FsCheck.NUnit.Property]
    public Property Median_PropagatesNaN()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
        {
            var rng = new Random(t.Seed);
            var withNaN = t.Values.ToList();
            withNaN.Insert(rng.Next(withNaN.Count + 1), double.NaN);
            return double.IsNaN(StatisticsHelper.Median(withNaN)).Label("NaN did not propagate");
        });
    }

    [Test]
    public void Median_EmptyOrNull_Throws()
    {
        Assert.Throws<ArgumentException>(() => StatisticsHelper.Median(Array.Empty<double>()));
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.Median(null!));
    }

    private static Arbitrary<(int[] Counts, int Seed, int Scale)> CountsArbitrary() =>
        (from n in Gen.Choose(1, 12)
         from counts in Gen.Choose(0, 50).ArrayOf(n)
         from seed in Gen.Choose(0, int.MaxValue)
         from scale in Gen.Choose(2, 7)
         select (counts, seed, scale)).ToArbitrary();

    /// <summary>
    /// Definition oracle H = −Σ pᵢ ln pᵢ over the non-zero classes, and the bounds 0 ≤ H ≤ ln(#non-zero classes)
    /// (maximum for equal counts).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ShannonIndex_MatchesDefinition_AndIsBoundedByLnRichness()
    {
        return Prop.ForAll(CountsArbitrary(), t =>
        {
            long total = t.Counts.Sum();
            if (total == 0)
            {
                return true.Label("zero total is rejected (separate test)");
            }

            double h = StatisticsHelper.ShannonIndex(t.Counts);
            double oracle = -t.Counts.Where(c => c > 0).Sum(c => (double)c / total * Math.Log((double)c / total));
            int richness = t.Counts.Count(c => c > 0);
            return (Math.Abs(h - oracle) <= 1e-12 && h >= 0.0 && h <= Math.Log(richness) + 1e-12)
                .Label($"H={h}, oracle={oracle}, ln S={Math.Log(richness)}");
        });
    }

    /// <summary>
    /// H depends only on the class proportions: permuting classes, inserting empty classes and multiplying every count by
    /// the same integer leave it unchanged (inserting zeros is exact; the other two up to rounding).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ShannonIndex_DependsOnlyOnProportions()
    {
        return Prop.ForAll(CountsArbitrary(), t =>
        {
            if (t.Counts.Sum() == 0)
            {
                return true.Label("zero total");
            }

            var rng = new Random(t.Seed);
            double h = StatisticsHelper.ShannonIndex(t.Counts);
            double permuted = StatisticsHelper.ShannonIndex(t.Counts.OrderBy(_ => rng.Next()).ToArray());
            var padded = t.Counts.ToList();
            padded.Insert(rng.Next(padded.Count + 1), 0);
            double withZero = StatisticsHelper.ShannonIndex(padded);
            double scaled = StatisticsHelper.ShannonIndex(t.Counts.Select(c => c * t.Scale).ToArray());

            return (Math.Abs(permuted - h) <= 1e-12 && withZero == h && Math.Abs(scaled - h) <= 1e-12)
                .Label($"H={h}, permuted={permuted}, withZero={withZero}, scaled={scaled}");
        });
    }

    /// <summary>k equal non-zero counts give H = ln k (maximum evenness).</summary>
    [FsCheck.NUnit.Property]
    public Property ShannonIndex_EqualCounts_IsLnK()
    {
        var arb = (from k in Gen.Choose(1, 30)
                   from c in Gen.Choose(1, 1000)
                   select (k, c)).ToArbitrary();
        return Prop.ForAll(arb, t =>
            (Math.Abs(StatisticsHelper.ShannonIndex(Enumerable.Repeat(t.c, t.k).ToArray()) - Math.Log(t.k)) <= 1e-12)
                .Label($"k={t.k}, c={t.c}"));
    }

    [Test]
    public void ShannonIndex_InvalidCounts_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.ShannonIndex(null!));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndex(new[] { 0, 0 }));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndex(Array.Empty<int>()));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndex(new[] { 3, -1 }));
    }

    // =========================================================================
    // FIN-B24 heavy tier: StatisticsHelper members added by the finisher (F39 Student t / incomplete beta, F33/F66
    // R mean/var/bw.nrd0, F57 Shannon of weights). Oracles are definitional identities, never the production code.
    // =========================================================================

    private static Arbitrary<(double[] Weights, int Seed, int ScaleExp)> WeightsArbitrary() =>
        (from n in Gen.Choose(1, 15)
         from raw in Gen.Choose(0, 10_000).ArrayOf(n)
         from seed in Gen.Choose(0, int.MaxValue)
         from scaleExp in Gen.Choose(-6, 6)
         where raw.Any(v => v > 0)
         select (raw.Select(v => v / 1013.0).ToArray(), seed, scaleExp)).ToArbitrary();

    /// <summary>
    /// F57: H′ of real weights depends only on the proportions (permutation exact up to rounding; scaling every weight by
    /// 2^j is exact in IEEE arithmetic, so bit-identical), satisfies 0 ≤ H′ ≤ ln(#non-zero), and on integer-valued weights
    /// is bit-identical to <see cref="StatisticsHelper.ShannonIndex"/> (shared accumulation core).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ShannonIndexOfWeights_ProportionInvariant_Bounded_AgreesWithCountsForm()
    {
        return Prop.ForAll(WeightsArbitrary(), t =>
        {
            var rng = new Random(t.Seed);
            double h = StatisticsHelper.ShannonIndexOfWeights(t.Weights);
            double permuted = StatisticsHelper.ShannonIndexOfWeights(t.Weights.OrderBy(_ => rng.Next()).ToArray());
            double scale = Math.Pow(2, t.ScaleExp);
            double scaled = StatisticsHelper.ShannonIndexOfWeights(t.Weights.Select(w => w * scale).ToArray());
            int richness = t.Weights.Count(w => w > 0);
            int[] counts = t.Weights.Select(w => (int)Math.Round(w * 1013.0)).ToArray();
            double fromWeights = StatisticsHelper.ShannonIndexOfWeights(counts.Select(c => (double)c).ToArray());
            double fromCounts = StatisticsHelper.ShannonIndex(counts);
            return (Math.Abs(permuted - h) <= 1e-12 && scaled == h && h >= 0.0 && h <= Math.Log(richness) + 1e-12
                    && fromWeights == fromCounts)
                .Label($"H={h}, permuted={permuted}, scaled={scaled}, lnS={Math.Log(richness)}, w={fromWeights}, c={fromCounts}");
        });
    }

    private static Arbitrary<(double X, double A, double B)> BetaArbitrary() =>
        (from xi in Gen.Choose(0, 1000)
         from ai in Gen.Choose(1, 4000)
         from bi in Gen.Choose(1, 4000)
         select (xi / 1000.0, ai / 40.0, bi / 40.0)).ToArbitrary(); // a, b ∈ [0.025, 100]

    /// <summary>
    /// F39: I_x(a, b) ∈ [0, 1], the reflection identity I_x(a, b) + I_{1−x}(b, a) = 1 (Abramowitz &amp; Stegun 26.5.2),
    /// the boundary values I_0 = 0, I_1 = 1, and monotonicity in x (a CDF in x).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RegularizedIncompleteBeta_Reflection_Bounds_Monotone()
    {
        return Prop.ForAll(BetaArbitrary(), t =>
        {
            double i = StatisticsHelper.RegularizedIncompleteBeta(t.X, t.A, t.B);
            double reflected = StatisticsHelper.RegularizedIncompleteBeta(1.0 - t.X, t.B, t.A);
            double next = StatisticsHelper.RegularizedIncompleteBeta(Math.Min(1.0, t.X + 0.001), t.A, t.B);
            bool bounds = StatisticsHelper.RegularizedIncompleteBeta(0.0, t.A, t.B) == 0.0
                          && StatisticsHelper.RegularizedIncompleteBeta(1.0, t.A, t.B) == 1.0;
            return (i >= 0.0 && i <= 1.0 && Math.Abs(i + reflected - 1.0) <= 1e-12 && next >= i - 1e-15 && bounds)
                .Label($"I={i}, I_reflected={reflected}, I(x+0.001)={next}");
        });
    }

    /// <summary>
    /// F39: Student t CDF — symmetry F(−t) = 1 − F(t), F(0) = ½, F(±∞) = 0/1, monotone in t, and the closed forms for
    /// ν = 1 (Cauchy: ½ + atan(t)/π) and ν = 2 (½ + t/(2√(2 + t²))).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property StudentTCdf_Symmetry_Limits_Monotone_ClosedForms()
    {
        var arb = (from ti in Gen.Choose(-20_000, 20_000)
                   from nui in Gen.Choose(1, 2000)
                   select (T: ti / 1000.0, Nu: nui / 10.0)).ToArbitrary();
        return Prop.ForAll(arb, g =>
        {
            double f = StatisticsHelper.StudentTCdf(g.T, g.Nu);
            double fNeg = StatisticsHelper.StudentTCdf(-g.T, g.Nu);
            double fNext = StatisticsHelper.StudentTCdf(g.T + 0.01, g.Nu);
            double cauchy = StatisticsHelper.StudentTCdf(g.T, 1.0);
            double nu2 = StatisticsHelper.StudentTCdf(g.T, 2.0);
            bool limits = StatisticsHelper.StudentTCdf(0.0, g.Nu) == 0.5
                          && StatisticsHelper.StudentTCdf(double.NegativeInfinity, g.Nu) == 0.0
                          && StatisticsHelper.StudentTCdf(double.PositiveInfinity, g.Nu) == 1.0;
            return (f >= 0.0 && f <= 1.0 && Math.Abs(f + fNeg - 1.0) <= 1e-13 && fNext >= f && limits
                    && Math.Abs(cauchy - (0.5 + Math.Atan(g.T) / Math.PI)) <= 1e-13
                    && Math.Abs(nu2 - (0.5 + g.T / (2.0 * Math.Sqrt(2.0 + g.T * g.T)))) <= 1e-13)
                .Label($"F={f}, F(-t)={fNeg}, F(t+0.01)={fNext}, cauchy={cauchy}, nu2={nu2}");
        });
    }

    private static Arbitrary<(double[] X, double[] Y, double Shift, int Seed)> TwoSamplesArbitrary() =>
        (from nx in Gen.Choose(2, 25)
         from ny in Gen.Choose(2, 25)
         from xs in Gen.Choose(-4000, 4000).ArrayOf(nx)
         from ys in Gen.Choose(-4000, 4000).ArrayOf(ny)
         from shift in Gen.Choose(-64, 64)
         from seed in Gen.Choose(0, int.MaxValue)
         where xs.Distinct().Count() > 1 && ys.Distinct().Count() > 1
         select (xs.Select(v => v / 128.0).ToArray(), ys.Select(v => v / 128.0).ToArray(), (double)shift, seed))
        .ToArbitrary();

    /// <summary>
    /// F39: Welch and one-sample t-test p-values lie in [0, 1]; Welch is symmetric in its samples, invariant to a common
    /// location shift (t and the Satterthwaite ν depend only on x̄ − ȳ and the variances) and to a common positive scale;
    /// the one-sample test is invariant to shifting the data and μ together. Dyadic inputs keep the shifts exact, so the
    /// tolerance covers only the accumulation order of the mean. Scaling is not bit-exact: ν uses stderr⁴ through
    /// <c>pow</c> (as R's <c>stderr^4</c> → <c>R_pow</c> → libm <c>pow</c>), which is not exactly homogeneous (seen: 1.1e−15).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TTests_PValueRange_Symmetry_LocationAndScaleInvariance()
    {
        return Prop.ForAll(TwoSamplesArbitrary(), t =>
        {
            double p = StatisticsHelper.WelchTTestPValue(t.X, t.Y);
            double pSwap = StatisticsHelper.WelchTTestPValue(t.Y, t.X);
            double pShift = StatisticsHelper.WelchTTestPValue(t.X.Select(v => v + t.Shift).ToArray(), t.Y.Select(v => v + t.Shift).ToArray());
            double pScale = StatisticsHelper.WelchTTestPValue(t.X.Select(v => v * 4.0).ToArray(), t.Y.Select(v => v * 4.0).ToArray());
            double mu = t.Y[0];
            double q = StatisticsHelper.OneSampleTTestPValue(t.X, mu);
            double qShift = StatisticsHelper.OneSampleTTestPValue(t.X.Select(v => v + t.Shift).ToArray(), mu + t.Shift);
            return (p >= 0.0 && p <= 1.0 && q >= 0.0 && q <= 1.0 && pSwap == p && Math.Abs(pScale - p) <= 1e-12
                    && Math.Abs(pShift - p) <= 1e-9 && Math.Abs(qShift - q) <= 1e-9)
                .Label($"p={p}, swap={pSwap}, shift={pShift}, scale={pScale}, q={q}, qShift={qShift}");
        });
    }

    /// <summary>
    /// F39: R's <c>t.test</c> error cases are NaN — fewer than two observations and "data are essentially constant".
    /// </summary>
    [Test]
    public void TTests_DegenerateInputs_AreNaN()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.OneSampleTTestPValue(new[] { 1.0 }, 0.0), Is.NaN);
            Assert.That(StatisticsHelper.OneSampleTTestPValue(new[] { 2.5, 2.5, 2.5 }, 0.0), Is.NaN);
            Assert.That(StatisticsHelper.WelchTTestPValue(new[] { 1.0 }, new[] { 1.0, 2.0 }), Is.NaN);
            Assert.That(StatisticsHelper.WelchTTestPValue(new[] { 3.0, 3.0 }, new[] { 3.0, 3.0, 3.0 }), Is.NaN);
        });
    }

    /// <summary>
    /// F33/F66: R <c>mean</c>/<c>var</c> (long-double <see cref="StatisticsHelper.ExtendedPrecisionMean"/> /
    /// <see cref="StatisticsHelper.ExtendedPrecisionVariance"/> and the double <see cref="StatisticsHelper.RMean"/> /
    /// <see cref="StatisticsHelper.RSampleVariance"/>): mean ∈ [min, max]; permutation invariance and shift equivariance
    /// (mean(x + c) = mean(x) + c, var(x + c) = var(x)) within rounding; the long-double and double forms agree to a few
    /// ulps; var ≥ 0; on dyadic integers (exact sums) the mean is the exact rational rounded once.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RMeanVariance_PermutationShift_AgreeAcrossPrecisions()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
        {
            if (t.Values.Length < 2)
            {
                return true.Label("variance needs n ≥ 2");
            }

            var rng = new Random(t.Seed);
            double[] perm = t.Values.OrderBy(_ => rng.Next()).ToArray();
            const double c = 37.25;
            double[] shifted = t.Values.Select(v => v + c).ToArray();
            double scale = t.Values.Max(Math.Abs) + c;

            double m = StatisticsHelper.ExtendedPrecisionMean(t.Values);
            double v = StatisticsHelper.ExtendedPrecisionVariance(t.Values);
            double mPerm = StatisticsHelper.ExtendedPrecisionMean(perm);
            double vPerm = StatisticsHelper.ExtendedPrecisionVariance(perm);
            double mShift = StatisticsHelper.ExtendedPrecisionMean(shifted);
            double vShift = StatisticsHelper.ExtendedPrecisionVariance(shifted);
            double rm = StatisticsHelper.RMean(t.Values);
            double rv = StatisticsHelper.RSampleVariance(t.Values);

            double tolM = 1e-14 * scale;
            double tolV = 1e-12 * (v + scale * scale * 1e-3);
            return (m >= t.Values.Min() && m <= t.Values.Max() && v >= 0.0
                    && Math.Abs(mPerm - m) <= tolM && Math.Abs(vPerm - v) <= tolV
                    && Math.Abs(mShift - (m + c)) <= tolM && Math.Abs(vShift - v) <= tolV
                    && Math.Abs(rm - m) <= tolM && Math.Abs(rv - v) <= tolV)
                .Label($"m={m}, mPerm={mPerm}, mShift={mShift}, rm={rm}; v={v}, vPerm={vPerm}, vShift={vShift}, rv={rv}");
        });
    }

    /// <summary>On integers the sum Σx is exact in long double, so <see cref="StatisticsHelper.ExtendedPrecisionMean"/>
    /// (R <c>mean</c>: long-double Σx/n plus the residual, rounded once to double) is within one ulp of the correctly
    /// rounded Σx/n (the long-double intermediate permits a double-rounding tie, never more).</summary>
    [FsCheck.NUnit.Property]
    public Property ExtendedPrecisionMean_IntegerInput_WithinOneUlpOfExact()
    {
        var arb = (from n in Gen.Choose(1, 40)
                   from xs in Gen.Choose(-1_000_000, 1_000_000).ArrayOf(n)
                   select xs).ToArbitrary();
        return Prop.ForAll(arb, xs =>
        {
            double m = StatisticsHelper.ExtendedPrecisionMean(xs.Select(x => (double)x).ToArray());
            double exact = (double)xs.Sum(x => (long)x) / xs.Length;
            double ulp = Math.BitIncrement(Math.Abs(exact)) - Math.Abs(exact);
            return (Math.Abs(m - exact) <= ulp).Label($"mean={m}, exact={exact}");
        });
    }

    /// <summary>
    /// F33: R <c>bw.nrd0</c> is scale-equivariant — bw(2ʲ·x) = 2ʲ·bw(x) bit-identically (power-of-two scaling is exact
    /// in every step: sd, type-7 quantiles, the |x₁| fallback), invariant to a common shift for non-constant data
    /// (within rounding), positive, and permutation-invariant for non-constant data (sd and IQR are symmetric).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property BandwidthNrd0_ScaleEquivariant_ShiftInvariant_Positive()
    {
        return Prop.ForAll(ValuesArbitrary(), t =>
        {
            if (t.Values.Length < 2)
            {
                return true.Label("bw.nrd0 needs n ≥ 2");
            }

            double bw = StatisticsHelper.BandwidthNrd0(t.Values);
            double bwScaled = StatisticsHelper.BandwidthNrd0(t.Values.Select(v => v * 8.0).ToArray());
            if (t.Values.Distinct().Count() == 1)
            {
                return (bw > 0 && bwScaled == (t.Values[0] == 0 ? bw : 8.0 * bw)).Label($"constant: bw={bw}, scaled={bwScaled}");
            }

            var rng = new Random(t.Seed);
            double bwPerm = StatisticsHelper.BandwidthNrd0(t.Values.OrderBy(_ => rng.Next()).ToArray());
            double bwShift = StatisticsHelper.BandwidthNrd0(t.Values.Select(v => v + 11.5).ToArray());
            double tol = 1e-9 * bw;
            return (bw > 0 && bwScaled == 8.0 * bw && Math.Abs(bwPerm - bw) <= tol && Math.Abs(bwShift - bw) <= tol)
                .Label($"bw={bw}, scaled={bwScaled}, perm={bwPerm}, shift={bwShift}");
        });
    }
}
