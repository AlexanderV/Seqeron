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
}
