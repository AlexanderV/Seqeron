using SuffixTree.Esa.Tests.Contract;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The result of a locus-returning member may be the very variable that holds the locus, on the public surface (work package 1b.12, finding R01 of the review of 1b.10): the calls here are made on the concrete
/// <see cref="EsaSuffixTree"/> (the facade, every child-table layout), not through the interfaces, and they are the calls of the review (<c>tree.TryExtendRight(l, c, out l)</c> in the loop that advances a locus,
/// <c>tree.GetChildren(kids[1], kids)</c>). The member-by-member and the random-walk suites of <see cref="IndexAliasingTests"/> and <see cref="IndexAliasingPropertyTests"/> run on the facade, the index and the container reader too.
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.EsaFixtures))]
public class FacadeAliasingTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    private EsaSuffixTree Esa(string text) => (EsaSuffixTree)Tree(text);

    [Test]
    public void TheResultMayAliasTheLocus()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        IndexLocus apart = tree.Root, alias = tree.Root;
        foreach (char c in "abr")
        {
            Assert.That(tree.TryExtendRight(apart, c, out IndexLocus next), Is.True);
            apart = next;
            Assert.That(tree.TryExtendRight(alias, c, out alias), Is.True);          // failed at 'b': false
            Assert.That(alias, Is.EqualTo(apart));
        }
        IndexLocus left = apart;
        Assert.That(tree.TryExtendLeft(apart, 'd', out IndexLocus leftApart), Is.True);
        Assert.That(tree.TryExtendLeft(left, 'd', out left), Is.True);               // threw ArgumentException
        Assert.That(left, Is.EqualTo(leftApart));
    }

    [Test]
    public void TryExtendRight_SymbolOverload_AliasedCallEqualsTheSeparateCall()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        foreach (string pattern in new[] { "abr", "abra", "bra", "cad", "dabra", "ra", "a", "abracadabra" })
        {
            IndexLocus apart = tree.Root, alias = tree.Root;
            foreach (char c in pattern)
            {
                bool okApart = tree.TryExtendRight(apart, c, out IndexLocus next);
                apart = next;
                bool okAlias = tree.TryExtendRight(alias, c, out alias);
                Assert.That(okAlias, Is.EqualTo(okApart), $"'{pattern}' at '{c}'");
                Assert.That(alias, Is.EqualTo(apart), $"'{pattern}' at '{c}'");
                Assert.That(okAlias, Is.True, $"'{pattern}' occurs");
            }
            Assert.That(tree.Count(alias), Is.EqualTo(tree.CountOccurrences(pattern)), pattern);
        }
        IndexLocus absent = tree.Root;
        Assert.That(tree.TryExtendRight(absent, 'a', out absent), Is.True);
        Assert.That(tree.TryExtendRight(absent, 'a', out absent), Is.False, "aa does not occur");
        Assert.That(absent, Is.EqualTo(default(IndexLocus)), "a failed call leaves the default locus in the one variable");
    }

    [Test]
    public void TryExtendRight_SpanOverload_AliasedCallEqualsTheSeparateCall()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        IndexLocus abr = default;
        Assert.That(tree.TryExtendRight(tree.Root, "abr".AsSpan(), out abr), Is.True);
        foreach (string pattern in new[] { "a", "acad", "abra", "", "z", "ax" })
        {
            IndexLocus alias = abr;
            bool okApart = tree.TryExtendRight(abr, pattern.AsSpan(), out IndexLocus apart);
            bool okAlias = tree.TryExtendRight(alias, pattern.AsSpan(), out alias);
            Assert.That(okAlias, Is.EqualTo(okApart), $"'abr' + '{pattern}'");
            Assert.That(alias, Is.EqualTo(apart), $"'abr' + '{pattern}'");
            if (okAlias) Assert.That(alias.Length, Is.EqualTo(3 + pattern.Length), $"'abr' + '{pattern}': the depth of the extension");
        }
        IndexLocus s2 = abr;
        Assert.That(tree.TryExtendRight(s2, "a".AsSpan(), out s2), Is.True);
        Assert.That(s2.Length, Is.EqualTo(4), "abra: the depth was 1 when the result was cleared before it was read");
        Assert.That(tree.Count(s2), Is.EqualTo(2));
    }

    [Test]
    public void TryExtendLeft_AliasedCallEqualsTheSeparateCall()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        foreach (string pattern in new[] { "abr", "a", "ra", "bra", "cad", "" })
        {
            Assert.That(tree.TryExtendRight(tree.Root, pattern.AsSpan(), out IndexLocus locus), Is.True, pattern);
            foreach (char c in "abcdrz")
            {
                bool okApart = tree.TryExtendLeft(locus, c, out IndexLocus apart);
                IndexLocus alias = locus;
                bool okAlias = tree.TryExtendLeft(alias, c, out alias);
                Assert.That(okAlias, Is.EqualTo(okApart), $"'{c}' + '{pattern}'");
                Assert.That(alias, Is.EqualTo(apart), $"'{c}' + '{pattern}'");
                Assert.That(okAlias, Is.EqualTo(tree.CountOccurrences(c + pattern) > 0), $"'{c}' + '{pattern}' occurs");
            }
        }
    }

    [Test]
    public void TryParent_TryChild_TrySuffixLink_AliasedCallEqualsTheSeparateCall()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        Assert.That(tree.TryExtendRight(tree.Root, "abra".AsSpan(), out IndexLocus abra), Is.True);
        Assert.That(tree.IsExplicitNode(abra), Is.True, "abra is a branching node: it occurs twice, followed by $ and by c");

        IndexLocus parent = abra;
        bool okParent = tree.TryParent(abra, out IndexLocus parentApart);
        Assert.That(tree.TryParent(parent, out parent), Is.EqualTo(okParent));
        Assert.That(parent, Is.EqualTo(parentApart));

        IndexLocus child = abra;
        bool okChild = tree.TryChild(abra, 'c', out IndexLocus childApart);
        Assert.That(okChild, Is.True);
        Assert.That(tree.TryChild(child, 'c', out child), Is.EqualTo(okChild));
        Assert.That(child, Is.EqualTo(childApart));

        IndexLocus missing = abra;
        Assert.That(tree.TryChild(missing, 'z', out missing), Is.False);
        Assert.That(missing, Is.EqualTo(default(IndexLocus)));

        IndexLocus link = abra;
        bool okLink = tree.TrySuffixLink(abra, out IndexLocus linkApart);
        Assert.That(okLink, Is.True);
        Assert.That(tree.TrySuffixLink(link, out link), Is.EqualTo(okLink));
        Assert.That(link, Is.EqualTo(linkApart));

        IndexLocus root = tree.Root;
        Assert.That(tree.TrySuffixLink(root, out root), Is.False, "the root has no link");
        Assert.That(root, Is.EqualTo(default(IndexLocus)));
        root = tree.Root;
        Assert.That(tree.TryParent(root, out root), Is.False, "the root has no parent");
        Assert.That(root, Is.EqualTo(default(IndexLocus)));
    }

    [Test]
    public void GetChildren_ChildrenWrittenOverTheArrayThatHoldsTheNode()
    {
        EsaSuffixTree tree = Esa("abracadabra");
        var o = new IndexContractOracle("abracadabra");
        var kids = new IndexLocus[8];
        int count = tree.GetChildren(tree.Root, kids);
        Assert.That(count, Is.EqualTo(o.Children([]).Count), "the terminator leaf, a, b, c, d, r");
        IndexLocus second = kids[1];
        var expected = new IndexLocus[8];
        int expectedCount = tree.GetChildren(second, expected);
        Assert.That(expectedCount, Is.GreaterThan(2), "the node of 'a' has the children $, b, c, d");
        Assert.That(tree.GetChildren(kids[1], kids), Is.EqualTo(expectedCount));
        Assert.That(kids.Take(expectedCount), Is.EqualTo(expected.Take(expectedCount)));
        // and an array that is too small: nothing is written, the node stays where it is
        var small = new IndexLocus[2];
        small[1] = second;
        Assert.That(tree.GetChildren(small[1], small), Is.EqualTo(expectedCount));
        Assert.That(small[0], Is.EqualTo(default(IndexLocus)));
        Assert.That(small[1], Is.EqualTo(second));
    }

    [Test]
    public void EveryLocusOfAText_ThroughTheFacade_IsTheSameWhetherTheResultAliasesOrNot()
    {
        // every member on every locus of "mississippi" through the concrete class: the aliased answer is the answer of two variables
        const string text = "mississippi";
        EsaSuffixTree tree = Esa(text);
        var o = new IndexContractOracle(text);
        foreach (int[] w in AliasCatalog.Loci(o))
        {
            IndexLocus locus = o.LocusOf(w);
            foreach (int c in new[] { 'i', 's', 'p', 'm', 'z' })
            {
                IndexLocus a = locus;
                bool ok = tree.TryExtendRight(locus, c, out IndexLocus b);
                Assert.That(tree.TryExtendRight(a, c, out a), Is.EqualTo(ok), $"{IndexContractOracle.Render(w)} right {(char)c}");
                Assert.That(a, Is.EqualTo(b), $"{IndexContractOracle.Render(w)} right {(char)c}");
                a = locus;
                ok = tree.TryExtendRight(locus, new[] { (char)c, (char)c }, out b);
                Assert.That(tree.TryExtendRight(a, new[] { (char)c, (char)c }, out a), Is.EqualTo(ok), $"{IndexContractOracle.Render(w)} right pair {(char)c}");
                Assert.That(a, Is.EqualTo(b), $"{IndexContractOracle.Render(w)} right pair {(char)c}");
                a = locus;
                ok = tree.TryExtendLeft(locus, c, out b);
                Assert.That(tree.TryExtendLeft(a, c, out a), Is.EqualTo(ok), $"{IndexContractOracle.Render(w)} left {(char)c}");
                Assert.That(a, Is.EqualTo(b), $"{IndexContractOracle.Render(w)} left {(char)c}");
            }
            IndexLocus p = locus;
            bool okParent = tree.TryParent(locus, out IndexLocus parent);
            Assert.That(tree.TryParent(p, out p), Is.EqualTo(okParent), $"{IndexContractOracle.Render(w)} parent");
            Assert.That(p, Is.EqualTo(parent), $"{IndexContractOracle.Render(w)} parent");
            if (!o.IsExplicit(w)) continue;
            foreach (int key in new[] { -1, 'i', 's', 'p', 'm' })
            {
                IndexLocus k = locus;
                bool okChild = tree.TryChild(locus, key, out IndexLocus child);
                Assert.That(tree.TryChild(k, key, out k), Is.EqualTo(okChild), $"{IndexContractOracle.Render(w)} child {key}");
                Assert.That(k, Is.EqualTo(child), $"{IndexContractOracle.Render(w)} child {key}");
            }
            IndexLocus l = locus;
            bool okLink = tree.TrySuffixLink(locus, out IndexLocus link);
            Assert.That(tree.TrySuffixLink(l, out l), Is.EqualTo(okLink), $"{IndexContractOracle.Render(w)} link");
            Assert.That(l, Is.EqualTo(link), $"{IndexContractOracle.Render(w)} link");
        }
    }
}
