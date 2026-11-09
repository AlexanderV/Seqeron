using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Links;

/// <summary>The definition-based check of every kernel on every internal node of one text (shared by the default and the exhaustive fixtures).</summary>
internal static class LinkKernelChecks
{
    private static readonly LinkKernel[] Kernels = Enum.GetValues<LinkKernel>();
    private static readonly EndpointSource[] Sources = Enum.GetValues<EndpointSource>();

    internal sealed class Prepared
    {
        internal required string Text { get; init; }

        internal required EsaIndex Esa { get; init; }

        internal required int[] Sa { get; init; }

        internal required int[] Lcp { get; init; }

        internal required LinkContext[] Contexts { get; init; }
    }

    internal static Prepared Prepare(string text, int[] fanouts)
    {
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        EsaIndex esa = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
        int[] psi = LinkTables.BuildPsi(esa.SuffixArray, isa);
        return new Prepared
        {
            Text = text,
            Esa = esa,
            Sa = sa,
            Lcp = lcp,
            Contexts = fanouts.Select(f => new LinkContext(esa, isa, psi, new LcpRangeHelper(esa.Lcp, f))).ToArray(),
        };
    }

    /// <summary>The definitional link target of the node: the rank interval of the label minus its first symbol.</summary>
    internal static (LinkTarget Target, string Y) Expected(Prepared p, BruteInterval node)
    {
        int n = p.Text.Length;
        if (node.Ell == 1) return (new LinkTarget(0, n, -1, true), "");
        string y = p.Text.Substring(p.Sa[node.Lb] + 1, node.Ell - 1);
        int lo = -1, hi = -1;
        for (int r = 1; r <= n; r++)
        {
            if (!p.Text.AsSpan(p.Sa[r]).StartsWith(y, StringComparison.Ordinal)) continue;
            if (lo < 0) lo = r;
            else if (hi != r - 1) Assert.Fail($"'{p.Text}': the suffixes starting with '{y}' are not contiguous");
            hi = r;
        }
        if (lo < 0) Assert.Fail($"'{p.Text}': the label of node ({node.Lb}, {node.Rb}, {node.Ell}) minus its first symbol does not occur");
        if (p.Text.Length <= 400 && (!p.Esa.TryMatch(y, out int slo, out int shi) || slo != lo || shi != hi))
            Assert.Fail($"'{p.Text}': the two definitional oracles disagree for '{y}'");
        int right = hi + 1 <= n ? p.Lcp[hi + 1] : -1;
        int minInside = int.MaxValue;
        for (int k = lo + 1; k <= hi; k++) minInside = Math.Min(minInside, p.Lcp[k]);
        if (lo < hi && minInside != node.Ell - 1) Assert.Fail($"'{p.Text}': target of ({node.Lb}, {node.Rb}, {node.Ell}) is not an lcp-interval of value {node.Ell - 1}");
        return (new LinkTarget(lo, hi, Math.Max(p.Lcp[lo], right), p.Lcp[lo] > right), y);
    }

    internal static void CheckText(string text, int[] fanouts, string label = "")
    {
        Prepared p = Prepare(text, fanouts);
        foreach (BruteInterval node in BruteLcpTree.Intervals(p.Lcp))
        {
            if (node.Ell < 1) continue;
            var (expected, y) = Expected(p, node);
            var v = new LinkNode(node.Lb, node.Rb, node.Ell, p.Sa[node.Lb]);
            foreach (LinkContext ctx in p.Contexts)
                foreach (LinkKernel kernel in Kernels)
                    foreach (EndpointSource source in Sources)
                    {
                        LinkTarget actual = LinkKernels.Link(kernel, source, ctx, v);
                        if (actual != expected)
                            Assert.Fail($"{label} '{(text.Length > 60 ? text[..60] + "..." : text)}' node ({node.Lb}, {node.Rb}, depth {node.Ell}) Y '{(y.Length > 30 ? y[..30] + "..." : y)}': {kernel}/{source} gave {actual}, expected {expected}");
                    }
        }
    }
}

/// <summary>
/// Suffix-link kernels of prototype P1 (work package 1b.4a). Oracle (definition, ADR-001 section 5 fact F1): for an internal node <c>v</c> with path label <c>aY</c> the link is the
/// node with label <c>Y</c>, i.e. the rank interval of the suffixes that start with <c>Y</c> (the root for <c>|Y| = 0</c>). The interval of <c>Y</c> is taken from two independent
/// sources: a scan of the naive suffix array and the pattern search of <see cref="EsaIndex"/>. Every kernel (scan, child-table descent, RMQ + PSV/NSV, PSV/NSV) with both endpoint
/// sources (ISA, Psi) must return exactly that interval, with the parent depth and the last-child flag a child-table walk needs.
/// </summary>
[TestFixture]
public class LinkKernelTests
{
    private static readonly LinkKernel[] Kernels = Enum.GetValues<LinkKernel>();
    private static readonly EndpointSource[] Sources = Enum.GetValues<EndpointSource>();
    private static readonly int[] SmallFanouts = [2, 16];

    [Test]
    public void ExhaustiveSmallStrings_AllKernelsAllSources()
    {
        foreach (string s in ExhaustiveStrings.All(2, 10)) LinkKernelChecks.CheckText(s, SmallFanouts);
        foreach (string s in ExhaustiveStrings.All(3, 7)) LinkKernelChecks.CheckText(s, SmallFanouts);
        foreach (string s in ExhaustiveStrings.All(4, 5)) LinkKernelChecks.CheckText(s, SmallFanouts);
    }

    [Test]
    public void AdversarialFamilies_AllNodes()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 17, 64, 100, 257, 400 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
                LinkKernelChecks.CheckText(text, SmallFanouts, $"{name} n {n}");
    }

    [Test]
    public void RandomStrings_AllNodes()
    {
        var rng = new SplitMix64(0x1B4A01);
        foreach (int sigma in new[] { 1, 2, 3, 4, 20, 256 })
            for (int rep = 0; rep < 25; rep++)
            {
                int n = 1 + rng.NextInt(sigma >= 20 ? 600 : 400);
                LinkKernelChecks.CheckText(StringFamilies.Random(n, sigma, rng.Next()), SmallFanouts, $"random sigma {sigma} rep {rep}");
            }
    }

    [Test]
    public void NodesWithDepthOne_LinkToTheRoot()
    {
        EsaIndex esa = EsaIndex.Build("abracadabra", ChildLayoutKind.Tagged);
        int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
        var ctx = new LinkContext(esa, isa, LinkTables.BuildPsi(esa.SuffixArray, isa), new LcpRangeHelper(esa.Lcp, 4));
        var set = LinkNodeSet.Enumerate(esa.Lcp);
        int checkedNodes = 0;
        for (int i = 0; i < set.Count; i++)
        {
            if (set.Depth[i] != 1) continue;
            checkedNodes++;
            foreach (LinkKernel kernel in Kernels)
                foreach (EndpointSource source in Sources)
                    Assert.That(LinkKernels.Link(kernel, source, ctx, set.Node(i, esa.SuffixArray)), Is.EqualTo(new LinkTarget(0, 11, -1, true)));
        }
        Assert.That(checkedNodes, Is.GreaterThan(0));
    }

    [Test]
    public void LargerTexts_AllKernelsAgreeWithTheScan_AndASampleWithTheDefinition()
    {
        var cases = new List<(string Name, string Text)>
        {
            ("random-dna-30000", StringFamilies.RandomDna(30000, 21)),
            ("random-protein-20000", StringFamilies.RandomProtein(20000, 22)),
            ("fibonacci-30000", StringFamilies.Fibonacci(30000)),
            ("thue-morse-30000", StringFamilies.ThueMorse(30000)),
            ("a^n-3000", StringFamilies.AllSame(3000)),
            ("runs-L256-20000", StringFamilies.Runs(20000, 256, 23)),
            ("periodic-65-20000", StringFamilies.Periodic(20000, 65, 24)),
            ("akb-100-20000", StringFamilies.AkB(20000, 100)),
            ("debruijn-4-20000", StringFamilies.DeBruijn(4, 8, 20000)),
        };
        foreach (var (name, text) in cases)
        {
            EsaIndex esa = EsaIndex.Build(text, ChildLayoutKind.Tagged);
            int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
            int[] psi = LinkTables.BuildPsi(esa.SuffixArray, isa);
            var ctx = new LinkContext(esa, isa, psi, new LcpRangeHelper(esa.Lcp, 16));
            var set = LinkNodeSet.Enumerate(esa.Lcp);
            int[] sample = set.Sample(400, 5)
                .Concat(LinkNodeSet.Largest(set.Depth, 20)).Concat(LinkNodeSet.Largest(set.NodeDepth, 20)).Concat(LinkNodeSet.Largest(set.Rb.Zip(set.Lb, (r, l) => r - l).ToArray(), 20))
                .ToArray();
            int searched = 0;
            foreach (int i in sample)
            {
                LinkNode v = set.Node(i, esa.SuffixArray);
                LinkTarget oracle = LinkKernels.Link(LinkKernel.Scan, EndpointSource.Isa, ctx, v);
                foreach (LinkKernel kernel in Kernels)
                    foreach (EndpointSource source in Sources)
                        if (LinkKernels.Link(kernel, source, ctx, v) != oracle)
                            Assert.Fail($"{name}: node ({v.Lb}, {v.Rb}, depth {v.Depth}) {kernel}/{source} differs from the scan");
                if (searched++ < 60)
                {
                    // the definition: search the label without its first symbol
                    bool found = v.Depth == 1
                        ? (oracle.Lb == 0 && oracle.Rb == text.Length)
                        : esa.TryMatch(text.AsSpan(v.Pos + 1, v.Depth - 1), out int lo, out int hi) && lo == oracle.Lb && hi == oracle.Rb;
                    if (!found) Assert.Fail($"{name}: node ({v.Lb}, {v.Rb}, depth {v.Depth}): the scan differs from the pattern search");
                }
            }
        }
    }

    [Test]
    public void TheCountingDelegate_ReportsSteps_AndTheNoCountOneNone()
    {
        LinkKernelChecks.Prepared p = LinkKernelChecks.Prepare(StringFamilies.RandomDna(2000, 31), [16]);
        var set = LinkNodeSet.Enumerate(p.Esa.Lcp);
        LinkNode v = set.Node(set.Count / 2, p.Esa.SuffixArray);
        foreach (LinkKernel kernel in Kernels)
            foreach (EndpointSource source in Sources)
            {
                var counting = new CountSteps();
                Assert.That(LinkKernels.Get<CountSteps>(kernel, source)(p.Contexts[0], v, ref counting, out LinkTarget t1), Is.True);
                var none = new NoCount();
                Assert.That(LinkKernels.Get<NoCount>(kernel, source)(p.Contexts[0], v, ref none, out LinkTarget t2), Is.True);
                Assert.That(t2, Is.EqualTo(t1), $"{kernel}/{source}");
                Assert.That(counting.Steps, Is.GreaterThan(0), $"{kernel}/{source}");
                Assert.That(none.Steps, Is.EqualTo(0));
            }
    }

    [Test]
    public void StepBudget_StopsTheUnboundedKernels_AndTheBoundedOnesStayInsideIt()
    {
        const int n = 6000;
        LinkKernelChecks.Prepared p = LinkKernelChecks.Prepare(StringFamilies.AllSame(n), [16]);
        var set = LinkNodeSet.Enumerate(p.Esa.Lcp);
        // the deepest node of a^n has depth n - 1 and its link target is a^(n-2): wide (n) for the scan, node depth n - 2 for the descent
        int deepest = LinkNodeSet.Largest(set.Depth, 1)[0];
        LinkNode v = set.Node(deepest, p.Esa.SuffixArray);
        Assert.That(v.Depth, Is.EqualTo(n - 1));

        foreach (EndpointSource source in Sources)
        {
            Assert.That(LinkKernels.TryLink(LinkKernel.Scan, source, p.Contexts[0], v, 100, out _, out long scanSteps), Is.False, "scan on a^n exceeds a small budget");
            Assert.That(scanSteps, Is.GreaterThan(100));
            Assert.That(LinkKernels.TryLink(LinkKernel.Descent, source, p.Contexts[0], v, 100, out _, out long descentSteps), Is.False, "descent on a^n exceeds a small budget");
            Assert.That(descentSteps, Is.GreaterThan(100));
            Assert.That(LinkKernels.TryLink(LinkKernel.Scan, source, p.Contexts[0], v, 10_000_000, out LinkTarget scanned, out long used), Is.True);
            Assert.That(used, Is.GreaterThan(n / 2), "the scan reads the whole target interval");
            foreach (LinkKernel bounded in new[] { LinkKernel.RmqPsvNsv, LinkKernel.PsvNsv })
            {
                Assert.That(LinkKernels.TryLink(bounded, source, p.Contexts[0], v, 1500, out LinkTarget t, out long steps), Is.True, $"{bounded} stays inside a small budget on a^n");
                Assert.That(t, Is.EqualTo(scanned));
                Assert.That(steps, Is.LessThan(1500));
            }
        }
    }

    [Test]
    public void MissingStructures_AreReported()
    {
        EsaIndex esa = EsaIndex.Build("mississippi", ChildLayoutKind.Tagged);
        int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
        var onlyIsa = new LinkContext(esa, isa, null, new LcpRangeHelper(esa.Lcp, 4));
        var noHelper = new LinkContext(esa, isa, LinkTables.BuildPsi(esa.SuffixArray, isa), null);
        var set = LinkNodeSet.Enumerate(esa.Lcp);
        LinkNode v = set.Node(set.Count - 1, esa.SuffixArray);
        Assert.Throws<InvalidOperationException>(() => LinkKernels.Link(LinkKernel.Scan, EndpointSource.Psi, onlyIsa, v));
        Assert.Throws<InvalidOperationException>(() => LinkKernels.Link(LinkKernel.RmqPsvNsv, EndpointSource.Isa, noHelper, v));
        Assert.DoesNotThrow(() => LinkKernels.Link(LinkKernel.Scan, EndpointSource.Isa, noHelper, v), "the scan and the descent need no helper");
        Assert.DoesNotThrow(() => LinkKernels.Link(LinkKernel.Descent, EndpointSource.Isa, noHelper, v));
        EsaIndex notTagged = EsaIndex.Build("mississippi", ChildLayoutKind.Packed);
        Assert.Throws<NotSupportedException>(() => _ = new LinkContext(notTagged, isa, null, null));
    }
}

/// <summary>Full enumeration of the pre-registered exhaustive sizes (sigma 2 up to 14, sigma 3 up to 9, sigma 4 up to 7). Heavy; explicit, category <c>Exhaustive</c>.</summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.4a)")]
public class LinkKernelExhaustiveTests
{
    private static readonly int[] Fanouts = [4, 16];

    private static void All(int sigma, int maxLength)
    {
        foreach (string s in ExhaustiveStrings.All(sigma, maxLength)) LinkKernelChecks.CheckText(s, Fanouts);
    }

    [Test]
    public void AllBinaryStrings_UpTo14() => All(2, 14);

    [Test]
    public void AllTernaryStrings_UpTo9() => All(3, 9);

    [Test]
    public void AllFourLetterStrings_UpTo7() => All(4, 7);
}
