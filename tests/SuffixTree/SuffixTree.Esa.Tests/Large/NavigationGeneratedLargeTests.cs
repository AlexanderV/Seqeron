namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.5, the navigation side at 2^32 + 10 and 2^33 symbols over generated sections (ADR-005 section 3.7 layer 1; category <c>Large</c>, explicit): the same checks as <see cref="NavigationGeneratedTests"/> (the link of the node of a
/// factor is the node of the factor without its first symbol: for <c>a^n</c> the link of <c>a^k</c> is <c>a^(k-1)</c>; the parent; Lcp and Lca of ranks; the children) over 5-byte readers, a U64 LCP array, PSI and the RMQ pyramid at a size
/// where nothing of them exists in memory. A suffix link is a few reads of PSI and of the levels of the pyramid, so each test is a matter of seconds. Written before the engine in the sense that the sparse-shell tests of the same package
/// were: against the pack of S2.4 (32-bit, heap ISA, pyramid scan) they are the refusal by name.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class NavigationGeneratedLargeTests : LargeFixtureBase
{
    private const long TwoTo32 = 1L << 32;
    private const long TwoTo33 = 1L << 33;

    public static IEnumerable<TestCaseData> Shapes()
    {
        foreach (long n in new[] { TwoTo32 + 10, TwoTo33 }) yield return Case($"a^{n}", () => new PowerOfA(n));
        foreach ((long n, long k) in new[] { (TwoTo32 + 10, 1L), (TwoTo32 + 10, 2L), (TwoTo32 + 10, 7L), (TwoTo33 + 3, 1000L) }) yield return Case($"(a^{k} b)* n={n}", () => new PeriodicAkB(n, k));
    }

    private static TestCaseData Case(string name, Func<ClosedFormText> make) => new TestCaseData(make).SetArgDisplayNames(name);

    [TestCaseSource(nameof(Shapes))]
    public void TheSuffixLink_OfTheInternalNodeOfAFactor_IsTheNodeOfTheFactorWithoutItsFirstSymbol(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        Assert.That(NavChecks.SuffixLinks(rig, rig.Form.Name), Is.GreaterThan(5));
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheParent_OfTheInternalNodeOfAFactor_IsTheNodeOfThePrefixOfItsParentDepth(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        Assert.That(NavChecks.Parents(rig, rig.Form.Name), Is.GreaterThan(5));
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheLcpAndTheLca_OfRanksAcrossTheBoundaries_AreTheRangeMinimumAndTheNodeOfTheCommonPrefix(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        NavChecks.LcpAndLca(rig, rig.Form.Name, randomPairs: 60);
    }

    [TestCaseSource(nameof(Shapes))]
    public void TheChildren_OfTheInternalNodeOfAFactor_PartitionItsIntervalInTheOrderOfTheKeys(Func<ClosedFormText> make)
    {
        var rig = new NavRig(make());
        Assert.That(NavChecks.Children(rig, rig.Form.Name), Is.GreaterThan(5));
    }

    /// <summary>a^n: the link of a^k is a^(k-1) for k at the ends, across 2^31 and at 2^32 and above; the parent of a^k is a^(k-1) too.</summary>
    [TestCase(TwoTo32 + 10)]
    [TestCase(TwoTo33)]
    public void ThePowerOfA_HasTheChainOfLinksAndParents_AtEveryK(long n)
    {
        var rig = new NavRig(new PowerOfA(n));
        foreach (long k in new long[] { 1, 2, 3, 100, (1L << 31) - 1, 1L << 31, (1L << 31) + 1, TwoTo32 - 1, TwoTo32, TwoTo32 + 1, n / 2, n - 100, n - 3, n - 2, n - 1 }.Where(k => k >= 1 && k <= n - 1))
        {
            var node = new EsaNode64(k, n, k, k - 1, n - k, true);
            EsaNode64 link = rig.Nav.SuffixLink(node);
            Assert.That((link.Lb, link.Rb, link.Depth, link.ParentDepth), Is.EqualTo(k == 1 ? (0L, n, 0L, 0L) : (k - 1, n, k - 1, Math.Max(0, k - 2))), $"the link of a^{k}");
            EsaNode64 parent = rig.Nav.Parent(node);
            Assert.That((parent.Lb, parent.Rb, parent.Depth), Is.EqualTo(k == 1 ? (0L, n, 0L) : (k - 1, n, k - 1)), $"the parent of a^{k}");
        }
        Assert.That(rig.Pack.PayloadBytes, Is.EqualTo(0), "nothing was built");
    }
}
