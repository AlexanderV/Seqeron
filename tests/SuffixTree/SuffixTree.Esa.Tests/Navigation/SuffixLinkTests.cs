using SuffixTree.Esa.Tests.Analysis;
using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>Helpers shared by the suffix-link fixtures: node handles from the definitional tree, the node of a label found by descending with the navigator only.</summary>
internal static class SuffixLinkChecks
{
    /// <summary>The handle of the lcp-interval <c>[lb..rb]</c> of value <paramref name="depth"/>: parent depth and last-child flag from the two LCP values around it, as the definition of the tree gives them.</summary>
    internal static EsaNode HandleOf(EsaIndex index, int lb, int rb, int depth)
    {
        int n = index.Length;
        int left = index.Lcp[lb];
        int right = rb + 1 <= n ? index.Lcp[rb + 1] : -1;
        return new EsaNode(lb, rb, depth, Math.Max(left, right), index.SuffixArray[lb], left > right);
    }

    /// <summary>The node whose path label is exactly <paramref name="label"/>, found by <c>TryGetChild</c> and the edge symbols only; null when the label ends inside an edge or does not occur.</summary>
    internal static EsaNode? NodeOfLabel(EsaNavigator nav, string label)
    {
        EsaNode node = nav.Root;
        int i = 0;
        while (i < label.Length)
        {
            if (!nav.TryGetChild(node, label[i], out EsaNode child)) return null;
            int length = nav.LengthOf(child);
            if (i + length > label.Length) return null;
            for (int o = 0; o < length; o++)
                if (nav.GetEdgeSymbol(child, o) != label[i + o]) return null;
            i += length;
            node = child;
        }
        return node;
    }

    /// <summary>Every node below the root (preorder) with its path label (a leaf: the suffix without the terminator), through <c>GetChildren</c>.</summary>
    internal static List<(EsaNode Node, string Path)> AllNodes(EsaNavigator nav, string text)
    {
        var result = new List<(EsaNode, string)>();
        var stack = new Stack<EsaNode>();
        var kids = new List<EsaNode>();
        stack.Push(nav.Root);
        while (stack.Count > 0)
        {
            EsaNode node = stack.Pop();
            nav.GetChildren(node, kids);
            foreach (EsaNode child in kids)
            {
                result.Add((child, child.IsLeaf ? text.Substring(Math.Min(child.Pos, text.Length)) : text.Substring(child.Pos, child.Depth)));
                stack.Push(child);
            }
        }
        return result;
    }

    /// <summary>Checks the link of every node of the tree through the navigator: internal nodes against the node found by descending the label minus its first symbol, leaves against the leaf of the next suffix, the root against itself; also the prototypes (a) and (c2) on the internal nodes.</summary>
    internal static void CheckNavigatorLinks(string text, int fanout, string context)
    {
        EsaIndex index = EsaIndex.Build(text);
        EsaNavigationPack pack = EsaNavigationPack.Build(index, fanout);
        EsaNavigator nav = index.CreateNavigator();
        Assert.That(nav.GetSuffixLink(nav.Root), Is.EqualTo(nav.Root), $"{context}: the link of the root is the root");
        var nodes = AllNodes(nav, text);
        var leafByPosition = nodes.Where(t => t.Node.IsLeaf).ToDictionary(t => t.Node.Pos, t => t.Node);
        foreach (var (node, path) in nodes)
        {
            EsaNode viaNavigator = nav.GetSuffixLink(node);
            EsaNode viaPack = pack.SuffixLink(node);
            Assert.That(viaPack, Is.EqualTo(viaNavigator), $"{context}: pack (F{fanout}) and navigator agree on '{path}'");
            if (node.IsLeaf)
            {
                EsaNode expectedLeaf = node.Pos >= text.Length ? nav.Root : leafByPosition[node.Pos + 1];
                Assert.That(viaNavigator, Is.EqualTo(expectedLeaf), $"{context}: leaf of '{path}'");
                continue;
            }
            EsaNode? expected = NodeOfLabel(nav, path.Substring(1));
            Assert.That(expected, Is.Not.Null, $"{context}: the label of '{path}' without its first symbol is a node (fact F1)");
            Assert.That(viaNavigator, Is.EqualTo(expected!.Value), $"{context}: link of '{path}'");
            Assert.That(viaNavigator.Depth, Is.EqualTo(node.Depth - 1), context);
        }
    }
}

/// <summary>
/// The production suffix-link kernel (c2) behind <see cref="EsaNavigationPack.SuffixLink"/> and <c>EsaNavigator.GetSuffixLink</c> (work package 1b.4). Oracles: the definition (the link of the
/// node with label <c>aY</c> is the node with label <c>Y</c>: the rank interval of the suffixes that start with <c>Y</c>, from a scan of the naive suffix array and from the pattern search, the
/// root for <c>|Y| = 0</c>), the node found by descending <c>Y</c> with the navigator, and the prototype kernels (a) outward scan and (c2) of work package 1b.4a.
/// </summary>
[TestFixture]
public class SuffixLinkTests
{
    private static readonly int[] Fanouts = [2, 32];

    internal static void CheckKernel(string text, string context)
    {
        LinkKernelChecks.Prepared p = LinkKernelChecks.Prepare(text, [16]);
        var packs = Fanouts.Select(f => EsaNavigationPack.Build(p.Esa, f)).ToArray();
        int n = text.Length;
        foreach (BruteInterval interval in BruteLcpTree.Intervals(p.Lcp))
        {
            if (interval.Ell < 1) continue;
            var (expected, y) = LinkKernelChecks.Expected(p, interval);
            EsaNode node = SuffixLinkChecks.HandleOf(p.Esa, interval.Lb, interval.Rb, interval.Ell);
            var v = new LinkNode(interval.Lb, interval.Rb, interval.Ell, p.Sa[interval.Lb]);
            foreach (EsaNavigationPack pack in packs)
            {
                EsaNode actual = pack.SuffixLink(node);
                if (interval.Ell == 1)
                {
                    if (actual != new EsaNode(0, n, 0, 0, n, true)) Assert.Fail($"{context} '{text}' node ({interval.Lb}, {interval.Rb}) depth 1: expected the root, got {actual}");
                    continue;
                }
                var shape = new LinkTarget(actual.Lb, actual.Rb, actual.ParentDepth, actual.LastChild);
                if (shape != expected || actual.Depth != interval.Ell - 1 || actual.Pos != p.Sa[actual.Lb])
                    Assert.Fail($"{context} '{(text.Length > 60 ? text[..60] + "..." : text)}' node ({interval.Lb}, {interval.Rb}, depth {interval.Ell}) Y '{(y.Length > 30 ? y[..30] + "..." : y)}': got {actual}, expected {expected}");
            }
            // the prototypes agree with the definition (they are the extra oracles of the production kernel)
            foreach (LinkKernel kernel in new[] { LinkKernel.Scan, LinkKernel.PsvNsv })
                if (LinkKernels.Link(kernel, EndpointSource.Isa, p.Contexts[0], v) != expected)
                    Assert.Fail($"{context} '{text}': prototype {kernel} differs from the definition at ({interval.Lb}, {interval.Rb}, {interval.Ell})");
        }
    }

    [Test]
    public void ExhaustiveSmallStrings_AgainstDefinitionAndPrototypes()
    {
        foreach (string s in ExhaustiveStrings.All(2, 10)) CheckKernel(s, "all");
        foreach (string s in ExhaustiveStrings.All(3, 7)) CheckKernel(s, "all");
        foreach (string s in ExhaustiveStrings.All(4, 5)) CheckKernel(s, "all");
    }

    [Test]
    public void AdversarialFamilies_AllNodes()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 17, 64, 100, 257, 400 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
                CheckKernel(text, $"{name} n {n}");
    }

    [Test]
    public void RandomStrings_AllNodes()
    {
        var rng = new SplitMix64(0x1B4B20);
        foreach (int sigma in new[] { 1, 2, 3, 4, 20, 256 })
            for (int rep = 0; rep < 25; rep++)
            {
                int n = 1 + rng.NextInt(sigma >= 20 ? 600 : 400);
                CheckKernel(StringFamilies.Random(n, sigma, rng.Next()), $"random sigma {sigma} rep {rep}");
            }
    }

    [Test]
    public void EdgeTexts_AllNodes()
    {
        foreach (string s in AnalysisCorpus.Edge()) CheckKernel(s, "edge");
    }

    [Test]
    public void NavigatorLinks_AllNodesIncludingLeavesAndRoot_ExhaustiveSmall()
    {
        foreach (string s in ExhaustiveStrings.All(2, 8)) SuffixLinkChecks.CheckNavigatorLinks(s, 32, "all");
        foreach (string s in ExhaustiveStrings.All(3, 5)) SuffixLinkChecks.CheckNavigatorLinks(s, 2, "all");
    }

    [Test]
    public void NavigatorLinks_AdversarialAndRandom()
    {
        foreach (int n in new[] { 1, 2, 7, 33, 64, 150 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 31))
                SuffixLinkChecks.CheckNavigatorLinks(text, n % 2 == 0 ? 2 : 32, $"{name} n {n}");
        var rng = new SplitMix64(0x1B4B21);
        foreach (int sigma in new[] { 2, 4, 20 })
            for (int rep = 0; rep < 15; rep++)
                SuffixLinkChecks.CheckNavigatorLinks(StringFamilies.Random(1 + rng.NextInt(200), sigma, rng.Next()), 32, $"random sigma {sigma}");
    }

    [Test]
    public void DepthOneNodes_AndTheEmptySuffixLeaf_LinkToTheRoot()
    {
        EsaIndex index = EsaIndex.Build("abracadabra");
        EsaNavigator nav = index.CreateNavigator();
        EsaNode root = nav.Root;
        int depthOne = 0;
        foreach (var (node, _) in SuffixLinkChecks.AllNodes(nav, index.Text))
        {
            if (node.Depth == 1 && !node.IsLeaf)
            {
                depthOne++;
                Assert.That(nav.GetSuffixLink(node), Is.EqualTo(root));
            }
            if (node.IsLeaf && node.Pos == index.Length) Assert.That(nav.GetSuffixLink(node), Is.EqualTo(root), "the terminator-only leaf");
        }
        Assert.That(depthOne, Is.GreaterThan(0));
        Assert.That(nav.GetSuffixLink(root), Is.EqualTo(root));
    }

    [Test]
    public void LargerTexts_AgainstTheScanPrototype_AndASampleAgainstThePatternSearch()
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
            EsaIndex esa = EsaIndex.Build(text);
            int[] isa = LinkTables.BuildIsa(esa.SuffixArray);
            var ctx = new LinkContext(esa, isa, null, new LcpRangeHelper(esa.Lcp, 16));
            var set = LinkNodeSet.Enumerate(esa.Lcp);
            int[] sample = set.Sample(400, 5)
                .Concat(LinkNodeSet.Largest(set.Depth, 20)).Concat(LinkNodeSet.Largest(set.NodeDepth, 20)).Concat(LinkNodeSet.Largest(set.Rb.Zip(set.Lb, (r, l) => r - l).ToArray(), 20))
                .ToArray();
            EsaNavigationPack pack = esa.Navigation;
            int searched = 0;
            foreach (int i in sample)
            {
                LinkNode v = set.Node(i, esa.SuffixArray);
                LinkTarget oracle = LinkKernels.Link(LinkKernel.Scan, EndpointSource.Isa, ctx, v);
                EsaNode actual = pack.SuffixLink(SuffixLinkChecks.HandleOf(esa, v.Lb, v.Rb, v.Depth));
                if (v.Depth == 1)
                {
                    if (actual != pack.Root) Assert.Fail($"{name}: depth-1 node ({v.Lb}, {v.Rb}) does not link to the root");
                    continue;
                }
                if (new LinkTarget(actual.Lb, actual.Rb, actual.ParentDepth, actual.LastChild) != oracle)
                    Assert.Fail($"{name}: node ({v.Lb}, {v.Rb}, depth {v.Depth}) got {actual}, the scan says {oracle}");
                if (searched++ < 60 && !(esa.TryMatch(text.AsSpan(v.Pos + 1, v.Depth - 1), out int lo, out int hi) && lo == actual.Lb && hi == actual.Rb))
                    Assert.Fail($"{name}: node ({v.Lb}, {v.Rb}, depth {v.Depth}): the link differs from the pattern search");
            }
        }
    }

    [Test]
    public void AllNodesOfRealFibonacciAndRuns_LinkToTheNodeOfTheLabel()
    {
        // the unbounded cases of the outward scan: every node of a Fibonacci word and of a^n through the production kernel, checked by descending the label
        foreach (string text in new[] { StringFamilies.Fibonacci(1500), StringFamilies.ThueMorse(1200), StringFamilies.AllSame(500), StringFamilies.PeriodicAaB(900, 40) })
            SuffixLinkChecks.CheckNavigatorLinks(text, 32, $"n {text.Length}");
    }
}

/// <summary>Full enumeration of the pre-registered exhaustive sizes for the production link kernel (sigma 2 up to 14, sigma 3 up to 9, sigma 4 up to 7), every internal node against the definition and the prototypes. Heavy; explicit, category <c>Exhaustive</c>.</summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.4)")]
public class SuffixLinkExhaustiveTests
{
    private static void All(int sigma, int maxLength)
    {
        foreach (string s in ExhaustiveStrings.All(sigma, maxLength)) SuffixLinkTests.CheckKernel(s, $"all sigma {sigma}");
    }

    [Test]
    public void AllBinaryStrings_UpTo14() => All(2, 14);

    [Test]
    public void AllTernaryStrings_UpTo9() => All(3, 9);

    [Test]
    public void AllFourLetterStrings_UpTo7() => All(4, 7);
}
