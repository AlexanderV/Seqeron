using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// The oracle of the index contract is itself checked, before it is trusted against the engine: its tree of "banana" is the hand-derived tree, and its nodes are exactly the lcp-intervals of the
/// definitional interval tree of <see cref="BruteLcpTree"/> plus the leaves (one per suffix, the empty one included), for every string of the small alphabets.
/// </summary>
[TestFixture]
public class IndexContractOracleTests
{
    private static string Show(IEnumerable<int[]> strings) => string.Join(" ", strings.Select(IndexContractOracle.Render));

    [Test]
    public void Banana_NodesInPreorderAndTheirLoci()
    {
        var o = new IndexContractOracle("banana");
        Assert.That(o.Sa, Is.EqualTo(new[] { 6, 5, 3, 1, 0, 4, 2 }));
        Assert.That(o.Lcp, Is.EqualTo(new[] { 0, 0, 1, 3, 0, 0, 2 }));
        Assert.That(Show(o.Nodes()), Is.EqualTo("\"\" \"$\" \"a\" \"a$\" \"ana\" \"ana$\" \"anana$\" \"banana$\" \"na\" \"na$\" \"nana$\""));
        Assert.That(o.LocusOf([]), Is.EqualTo(new IndexLocus(0, 6, 0, 0)));
        Assert.That(o.LocusOf(IndexContractOracle.Symbols("a")), Is.EqualTo(new IndexLocus(1, 3, 1, 1)));
        Assert.That(o.LocusOf(IndexContractOracle.Symbols("an")), Is.EqualTo(new IndexLocus(2, 3, 2, 3)));
        Assert.That(o.LocusOf(IndexContractOracle.Symbols("ana")), Is.EqualTo(new IndexLocus(2, 3, 3, 3)));
        Assert.That(o.LocusOf(IndexContractOracle.Symbols("ban")), Is.EqualTo(new IndexLocus(4, 4, 3, 7)));
        Assert.That(o.LocusOf(IndexContractOracle.WithTerminator(IndexContractOracle.Symbols("na"))), Is.EqualTo(new IndexLocus(5, 5, 3, 3)));
        Assert.That(o.LocusOf([-1]), Is.EqualTo(new IndexLocus(0, 0, 1, 1)));
        Assert.That(o.Count(IndexContractOracle.Symbols("a")), Is.EqualTo(3));
        Assert.That(o.Count([-1]), Is.EqualTo(0));
        Assert.That(o.Count([]), Is.EqualTo(6));
        Assert.That(o.Positions(IndexContractOracle.Symbols("ana")).OrderBy(x => x), Is.EqualTo(new[] { 1, 3 }));
        Assert.That(IndexContractOracle.Render(o.Parent(IndexContractOracle.Symbols("an"))!), Is.EqualTo("\"a\""));
        Assert.That(IndexContractOracle.Render(o.Parent(IndexContractOracle.Symbols("ana"))!), Is.EqualTo("\"a\""));
        Assert.That(IndexContractOracle.Render(o.Parent(IndexContractOracle.Symbols("ban"))!), Is.EqualTo("\"\""));
        Assert.That(o.Parent([]), Is.Null);
        Assert.That(IndexContractOracle.Render(IndexContractOracle.SuffixLink(IndexContractOracle.Symbols("ana"))!), Is.EqualTo("\"na\""));
        Assert.That(o.ExtendLeft('b', IndexContractOracle.Symbols("an")), Is.EqualTo(IndexContractOracle.Symbols("ban")));
        Assert.That(o.ExtendLeft('a', IndexContractOracle.Symbols("an")), Is.Null);
        Assert.That(o.LcpOfRanks(2, 3), Is.EqualTo(3));
        Assert.That(o.LcpOfRanks(4, 4), Is.EqualTo(6));
        Assert.That(IndexContractOracle.Render(o.LcaString(2, 3)), Is.EqualTo("\"ana\""));
        Assert.That(IndexContractOracle.Render(o.LcaString(2, 2)), Is.EqualTo("\"ana$\""));
        Assert.That(IndexContractOracle.Render(o.LcaString(4, 5)), Is.EqualTo("\"\""));
    }

    [Test]
    public void EmptyText_OnlyTheRoot()
    {
        var o = new IndexContractOracle("");
        Assert.That(Show(o.Nodes()), Is.EqualTo("\"\""));
        Assert.That(o.LocusOf([]), Is.EqualTo(default(IndexLocus)));
        Assert.That(o.Count([]), Is.EqualTo(0));
        Assert.That(o.Children([]), Is.Empty);
        Assert.That(o.LcaString(0, 0), Is.Empty);
        Assert.That(() => o.LocusOf([-1]), Throws.TypeOf<InvalidOperationException>());
    }

    [Test]
    public void Nodes_AreTheLcpIntervalsOfTheDefinitionalTreePlusOneLeafPerSuffix()
    {
        foreach (string text in ExhaustiveStrings.All(2, 8).Concat(ExhaustiveStrings.All(3, 5)).Concat(ExhaustiveStrings.All(4, 4)))
        {
            var o = new IndexContractOracle(text);
            var nodes = o.Nodes();
            if (text.Length == 0)
            {
                Assert.That(nodes, Has.Count.EqualTo(1));
                continue;
            }
            var intervals = BruteLcpTree.Intervals(o.Lcp);
            Assert.That(nodes, Has.Count.EqualTo(intervals.Count + text.Length + 1), $"'{text}' node count 1 + I + (n + 1)");
            var internalLoci = nodes.Where(w => w.Length == 0 || w[^1] >= 0).Select(o.LocusOf).ToHashSet();
            var expected = intervals.Select(i => new IndexLocus(i.Lb, i.Rb, i.Ell, i.Ell)).ToHashSet();
            Assert.That(internalLoci, Is.EquivalentTo(expected), $"'{text}' internal nodes are the lcp-intervals");
            var leaves = nodes.Where(w => w.Length > 0 && w[^1] < 0).ToList();
            Assert.That(leaves, Has.Count.EqualTo(text.Length + 1), $"'{text}' one leaf per suffix");
            for (int r = 0; r <= text.Length; r++)
                Assert.That(leaves.Select(IndexContractOracle.Key), Does.Contain(IndexContractOracle.Key(o.LeafString(r))), $"'{text}' leaf of rank {r}");
            foreach (int[] w in nodes.Where(w => w.Length > 0 && w[^1] >= 0))
                Assert.That(o.IsExplicit(w[1..]), Is.True, $"'{text}' the suffix link target of {IndexContractOracle.Render(w)} is a node");
        }
    }
}
