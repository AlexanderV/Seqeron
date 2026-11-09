using SuffixTree.Esa.Tests.Analysis;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// The node model and the child navigation of <see cref="EsaNavigator"/> (work package 1b.4, ADR-001 section 4): the ESA seen through <see cref="ISuffixTreeNavigator{TNode}"/> must be node for
/// node the terminated Ukkonen tree. Oracles: a hand-built tree ("banana"), the definition (path labels, occurrences by brute force) and the Ukkonen tree through its public <c>Traverse</c>
/// (matched record by record: path label, leaf flag, leaf count, child count, branch key). Every present child key and some absent ones go through <c>TryGetChild</c>.
/// </summary>
[TestFixture]
public class EsaNavigatorTests
{
    private static readonly ChildLayoutKind[] Layouts = Enum.GetValues<ChildLayoutKind>();

    [Test]
    public void Contract_IsAValueTypeNavigatorOfEsaNode()
    {
        Assert.That(typeof(ISuffixTreeNavigator<EsaNode>).IsAssignableFrom(typeof(EsaNavigator)), "EsaNavigator : ISuffixTreeNavigator<EsaNode>");
        Assert.That(typeof(EsaNavigator).IsValueType, "a struct, so that the JIT specialises the generic algorithms");
        Assert.That(typeof(EsaNode).IsValueType, "the handle is a value");
        Assert.That(EsaNode.Null.IsNull, Is.True);
        Assert.That(new EsaNode(0, 0, 1, 0, 3, false).IsNull, Is.False);
        Assert.That(new EsaNode(4, 4, 7, 2, 3, true).IsLeaf);
        Assert.That(new EsaNode(4, 9, 7, 2, 3, true).IsLeaf, Is.False);
        Assert.That(new EsaNode(4, 9, 7, 2, 3, true).EdgeLength, Is.EqualTo(5));
        Assert.That(new EsaNode(4, 9, 7, 2, 3, true).RankCount, Is.EqualTo(6));
    }

    [Test]
    public void Contract_TextIsTheIndexedString()
    {
        EsaIndex index = EsaIndex.Build("banana");
        EsaNavigator nav = index.CreateNavigator();
        Assert.That(nav.Text.Length, Is.EqualTo(6));
        Assert.That(nav.Text.Substring(0, 6), Is.EqualTo("banana"));
        Assert.That(nav.Text[3], Is.EqualTo('a'));
    }

    [Test]
    public void Banana_HandBuiltTree()
    {
        // SA = [6 5 3 1 0 4 2], LCP = [0 0 1 3 0 0 2]; interior intervals: root [0..6], "a" [1..3], "ana" [2..3], "na" [5..6]
        EsaIndex index = EsaIndex.Build("banana");
        EsaNavigator nav = index.CreateNavigator();
        EsaNode root = nav.Root;
        Assert.That(root, Is.EqualTo(new EsaNode(0, 6, 0, 0, 6, true)));

        var level0 = new List<EsaNode>();
        nav.GetChildren(root, level0);
        Assert.That(level0.OrderBy(c => nav.GetEdgeSymbol(c, 0)).ToList(), Is.EqualTo(new[]
        {
            new EsaNode(0, 0, 1, 0, 6, false),   // the empty suffix: the terminator child
            new EsaNode(1, 3, 1, 0, 5, false),   // "a"
            new EsaNode(4, 4, 7, 0, 0, false),   // "banana$"
            new EsaNode(5, 6, 2, 0, 4, true),    // "na"
        }));

        var a = new EsaNode(1, 3, 1, 0, 5, false);
        var underA = new List<EsaNode>();
        nav.GetChildren(a, underA);
        Assert.That(underA.OrderBy(c => nav.GetEdgeSymbol(c, 0)).ToList(), Is.EqualTo(new[]
        {
            new EsaNode(1, 1, 2, 1, 5, false),   // "a$"
            new EsaNode(2, 3, 3, 1, 3, true),    // "ana"
        }));
        var ana = new EsaNode(2, 3, 3, 1, 3, true);
        var underAna = new List<EsaNode>();
        nav.GetChildren(ana, underAna);
        Assert.That(underAna.OrderBy(c => nav.GetEdgeSymbol(c, 0)).ToList(), Is.EqualTo(new[]
        {
            new EsaNode(2, 2, 4, 3, 3, false),   // "ana$"
            new EsaNode(3, 3, 6, 3, 1, true),    // "anana$"
        }));
        var na = new EsaNode(5, 6, 2, 0, 4, true);
        var underNa = new List<EsaNode>();
        nav.GetChildren(na, underNa);
        Assert.That(underNa.OrderBy(c => nav.GetEdgeSymbol(c, 0)).ToList(), Is.EqualTo(new[]
        {
            new EsaNode(5, 5, 3, 2, 4, false),   // "na$"
            new EsaNode(6, 6, 5, 2, 2, true),    // "nana$"
        }));

        Assert.That(nav.TryGetChild(root, 'a', out EsaNode viaKey) && viaKey == a);
        Assert.That(nav.TryGetChild(root, -1, out EsaNode term) && term == new EsaNode(0, 0, 1, 0, 6, false), "the terminator child of the root");
        Assert.That(nav.TryGetChild(a, -1, out EsaNode aTerm) && aTerm == new EsaNode(1, 1, 2, 1, 5, false), "the suffix 'a' is the terminator child of the node 'a'");
        Assert.That(nav.TryGetChild(a, 'n', out EsaNode an) && an == ana);
        Assert.That(nav.TryGetChild(a, 'b', out _), Is.False);
        Assert.That(nav.TryGetChild(new EsaNode(4, 4, 7, 0, 0, false), 'b', out _), Is.False, "a leaf has no children");
        Assert.That(nav.LengthOf(ana), Is.EqualTo(2));
        Assert.That(nav.GetEdgeSymbol(ana, 0), Is.EqualTo('n'));
        Assert.That(nav.GetEdgeSymbol(ana, 1), Is.EqualTo('a'));
        Assert.That(nav.GetEdgeSymbol(new EsaNode(6, 6, 5, 2, 2, true), 2), Is.EqualTo(-1), "the last symbol of a leaf edge is the terminator");
        Assert.That(nav.FindAnyLeafPosition(ana, 1), Is.AnyOf(1, 3));
        var leaves = new List<int>();
        nav.CollectLeaves(ana, 1, leaves);
        Assert.That(leaves.Order().ToList(), Is.EqualTo(new[] { 1, 3 }));
        leaves.Clear();
        nav.CollectLeaves(root, 0, leaves);
        Assert.That(leaves.Order().ToList(), Is.EqualTo(new[] { 0, 1, 2, 3, 4, 5 }), "the empty suffix is not a position");
        Assert.That(nav.FindAnyLeafPosition(new EsaNode(0, 0, 1, 0, 6, false), 0), Is.EqualTo(-1), "the terminator leaf of the root has no position");
    }

    [Test]
    public void EmptyAndSingleCharacterTexts()
    {
        NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(""), true, "empty");
        NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build("a"), true, "a");
        NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build("\0"), true, "NUL");
        NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build("￿"), true, "U+FFFF");
    }

    [Test]
    public void ExhaustiveSmallStrings_WalkEqualsUkkonen()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 10), (3, 7), (4, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(s), true, $"'{s}'");
    }

    [Test]
    public void AllChildLayouts_WalkEqualsUkkonen()
    {
        foreach (ChildLayoutKind layout in Layouts)
        {
            foreach (string s in ExhaustiveStrings.All(2, 8))
                NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(s, layout), true, $"{layout} '{s}'");
            foreach (var (name, text) in StringFamilies.Adversarial(60, 5))
                NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(text, layout), true, $"{layout} {name}");
        }
    }

    [Test]
    public void AdversarialFamilies_WalkEqualsUkkonen()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 17, 64, 100, 257 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 21))
                NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(text), n <= 64, $"{name} n {n}");
    }

    [Test]
    public void RandomStrings_WalkEqualsUkkonen()
    {
        var rng = new SplitMix64(0x1B4B01);
        foreach (int sigma in new[] { 1, 2, 3, 4, 20, 256 })
            for (int rep = 0; rep < 25; rep++)
            {
                int n = rng.NextInt(300);
                string s = StringFamilies.Random(n, sigma, rng.Next());
                NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(s), n <= 64, $"random sigma {sigma} n {n} #{rep}");
            }
    }

    [Test]
    public void EdgeTexts_WalkEqualsUkkonen()
    {
        foreach (string s in AnalysisCorpus.Edge())
            NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(s), true, $"edge '{s}'");
    }

    [Test]
    public void LargerTexts_WalkEqualsUkkonen()
    {
        foreach (var (name, text) in new[]
        {
            ("fibonacci", StringFamilies.Fibonacci(4000)),
            ("thue-morse", StringFamilies.ThueMorse(3000)),
            ("a^2500", StringFamilies.AllSame(2500)),
            ("dna", StringFamilies.RandomDna(12000, 7)),
            ("periodic", StringFamilies.Periodic(5000, 37, 9)),
            ("a^(p-1)b", StringFamilies.PeriodicAaB(3000, 50)),
        })
            NavigatorWalk.AssertEqualsUkkonen(EsaIndex.Build(text), false, name);
    }

    [Test]
    public void GetChildren_ReplacesTheList_AndDoesNotAllocateNodesForLeaves()
    {
        EsaIndex index = EsaIndex.Build("abracadabra");
        EsaNavigator nav = index.CreateNavigator();
        var list = new List<EsaNode> { EsaNode.Null, EsaNode.Null };
        nav.GetChildren(nav.Root, list);
        Assert.That(list.Any(nav.IsNull), Is.False, "the previous content is gone");
        Assert.That(list.Count, Is.EqualTo(6), "root of abracadabra: $, a, b, c, d, r");
        nav.GetChildren(list.First(c => c.IsLeaf), list);
        Assert.That(list, Is.Empty, "a leaf has no children");
    }

    [Test]
    public void Navigators_OfOneIndexAgree_AndAreIndependent()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(2000, 3));
        EsaNavigator a = index.CreateNavigator();
        EsaNavigator b = index.CreateNavigator();
        Assert.That(a.Root, Is.EqualTo(b.Root));
        List<TreeRecord> first = NavigatorWalk.Walk(index, false, "x");
        List<TreeRecord> second = NavigatorWalk.Walk(index, false, "y");
        Assert.That(second, Is.EqualTo(first));
    }
}
