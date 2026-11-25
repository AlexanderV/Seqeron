namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.5, the long navigation engine over generated sections (default set: n of a few thousand; the same checks at 2^32 + 10 and 2^33 are <see cref="NavigationGeneratedLargeTests"/>). PSI, the RMQ levels, the suffix array, the LCP array
/// and the child cells of <c>a^n</c>, <c>(a^k b)*</c> and <c>u^m</c> of a de Bruijn word are computed by the closed forms of S2.1 (5-byte readers, U64 LCP array); the navigator and the pack run over them exactly as over a container, and every
/// answer is compared with the closed form's own (<see cref="NavChecks"/>).
/// </summary>
[TestFixture]
public class NavigationGeneratedTests
{
    public static IEnumerable<TestCaseData> Shapes()
    {
        foreach (int n in new[] { 2, 3, 10, 33, 300, 3_000 }) yield return Case($"a^{n}", () => new PowerOfA(n));
        foreach ((int n, int k) in new[] { (60, 1), (3_000, 1), (3_000, 2), (2_999, 5), (4_000, 9), (1_000, 40) }) yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(n, k));
        foreach ((int sigma, int order, int m) in new[] { (2, 6, 9), (2, 8, 12), (3, 4, 10), (4, 3, 30) }) yield return Case($"u^{m} de Bruijn {sigma}^{order}", () => new DeBruijnRepeat(sigma, order, m));
    }

    private static TestCaseData Case(string name, Func<ClosedFormText> make) => new TestCaseData(make).SetArgDisplayNames(name);

    [TestCaseSource(nameof(Shapes))]
    public void TheSuffixLink_OfEveryInternalNodeAmongTheFactors_IsTheNodeOfTheFactorWithoutItsFirstSymbol(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        int count = NavChecks.SuffixLinks(rig, rig.Form.Name);
        if (rig.Form.Length >= 10) Assert.That(count, Is.GreaterThan(0), $"{rig.Form.Name}: the factors must reach an internal node");
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheParent_OfEveryInternalNodeAmongTheFactors_IsTheNodeOfThePrefixOfItsParentDepth(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        int count = NavChecks.Parents(rig, rig.Form.Name);
        if (rig.Form.Length >= 10) Assert.That(count, Is.GreaterThan(0));
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheLcpAndTheLca_OfRanks_AreTheRangeMinimumAndTheNodeOfTheCommonPrefix(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        NavChecks.LcpAndLca(rig, rig.Form.Name, randomPairs: 40);
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheChildren_OfEveryInternalNodeAmongTheFactors_PartitionItsIntervalInTheOrderOfTheKeys(Func<ClosedFormText> make)
    {
        ClosedFormText f = make();
        if (!f.HasChildTable) Assert.Ignore("the closed form of u^m has no child table: it is the scan of the engine over the generated LCP (ADR-005 section 3.7)");
        var rig = new NavRig(f);
        int count = NavChecks.Children(rig, f.Name);
        if (f.Length >= 10) Assert.That(count, Is.GreaterThan(0));
    }

    [Test]
    public void TheGeneratedPyramid_EqualsThePyramidBuiltFromTheLcpArray()
    {
        foreach (ClosedFormText f in new ClosedFormText[] { new PowerOfA(700), new PeriodicAkB(1_999, 3), new DeBruijnRepeat(2, 7, 10) })
        {
            var rig = new NavRig(f);
            using WideLcpPyramid built = WideLcpPyramid.Build(rig.Lcp, 32);
            long n = f.Length;
            var rng = new TestSupport.SplitMix64(11);
            for (int i = 0; i < 400; i++)
            {
                long a = (long)(rng.Next() % (ulong)(n + 1)), b = (long)(rng.Next() % (ulong)(n + 1));
                long lo = Math.Min(a, b), hi = Math.Max(a, b);
                Assert.That(rig.Ranges.Min(lo, hi), Is.EqualTo(built.Min(lo, hi)), $"{f.Name} Min({lo}, {hi})");
                long bound = (long)(rng.Next() % 40);
                Assert.That(rig.Ranges.PrevLess(a, bound), Is.EqualTo(built.PrevLess(a, bound)), $"{f.Name} PrevLess({a}, {bound})");
                Assert.That(rig.Ranges.NextLess(a, bound), Is.EqualTo(built.NextLess(a, bound)), $"{f.Name} NextLess({a}, {bound})");
            }
        }
    }
}
