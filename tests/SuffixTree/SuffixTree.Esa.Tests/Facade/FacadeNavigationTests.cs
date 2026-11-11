using SuffixTree.Esa;
using SuffixTree.Esa.Tests.Contract;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The navigation and suffix-array surfaces (<see cref="ISuffixTreeNavigation"/>, <see cref="ISuffixArrayAccess"/>) of the facade, through the facade object: the contract checker of
/// work package 1b.5 (every operation on every kind of locus against an oracle derived from the definitions) runs on <see cref="EsaSuffixTree"/> itself, and the int surface and the
/// long surface of the one object must tell the same story. The conformance suite of work package 1b.0 runs on the facade through <c>IndexEngineFactories.All</c> ("esa-facade").
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.EsaFixtures))]
public class FacadeNavigationTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    private static readonly string[] Texts =
    [
        "", "a", "aa", "ab", "aba", "banana", "mississippi", "abracadabra", "aaaaaaaa", "abababab", "abaababaabaab", "xyzxyz", "\uD800a\uD800", "\0￿\0", "GATTACAGATTACA",
    ];

    private (ISuffixTreeNavigation Nav, ISuffixArrayAccess Sa, ISuffixTree Tree) Facade(string text)
    {
        ISuffixTree tree = Tree(text);
        return ((ISuffixTreeNavigation)tree, (ISuffixArrayAccess)tree, tree);
    }

    [Test]
    public void TheFacade_ImplementsTheThreeContractFamilies()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(tree, Is.InstanceOf<ISuffixTreeNavigation>());
            Assert.That(tree, Is.InstanceOf<ISuffixArrayAccess>());
            Assert.That(tree, Is.InstanceOf<ISuffixIndex>());
            Assert.That(tree, Is.InstanceOf<IDisposable>());
            Assert.That(tree, Is.InstanceOf<EsaSuffixTree>());
        });
    }

    [Test]
    public void ContractChecker_RunsOnTheFacade_ForEveryLocusOfTheSmallTexts()
    {
        foreach (string text in Texts)
        {
            var (nav, sa, _) = Facade(text);
            Assert.That(nav, Is.SameAs(sa));
            EsaContractChecker.CheckAll((EsaSuffixTree)nav, new IndexContractOracle(text), $"[{Backend.Name}] '{text}'", CheckOptions.Thorough);
        }
    }

    [Test]
    public void ContractChecker_RunsOnTheFacade_ForAdversarialFamilies()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(60, 5))
        {
            var (nav, _, _) = Facade(text);
            EsaContractChecker.CheckAll((EsaSuffixTree)nav, new IndexContractOracle(text), $"[{Backend.Name}] {name}", CheckOptions.Sampled);
        }
    }

    [Test]
    public void EngineNameAndLengths()
    {
        foreach (string text in Texts)
        {
            var (nav, _, tree) = Facade(text);
            Assert.Multiple(() =>
            {
                Assert.That(nav.Engine, Is.EqualTo("esa"));
                Assert.That(nav.TextLength, Is.EqualTo(text.Length));
                Assert.That(nav.Text.Length, Is.EqualTo(text.Length));
                Assert.That(tree.Text.Length, Is.EqualTo(text.Length));
                Assert.That(nav.Text.SymbolAt(0), Is.EqualTo(text.Length == 0 ? -1 : text[0]));
                Assert.That(nav.Text.SymbolAt(text.Length), Is.EqualTo(-1), "the terminator");
            });
        }
    }

    [Test]
    public void TextProperty_ISuffixTreeSearch_And_ISuffixIndex_AreTwoViewsOfTheSameText()
    {
        var (nav, _, tree) = Facade("abracadabra");
        var buffer = new char[11];
        Assert.Multiple(() =>
        {
            Assert.That(tree.Text.ToString(), Is.EqualTo("abracadabra"));
            Assert.That(nav.Text.CopyTo(0, buffer), Is.EqualTo(11));
            Assert.That(new string(buffer), Is.EqualTo("abracadabra"));
            for (int i = 0; i < 11; i++) Assert.That(nav.Text.SymbolAt(i), Is.EqualTo(tree.Text[i]), $"position {i}");
        });
    }

    [Test]
    public void IntSurface_And_LongSurface_AgreeOnCountsAndPositions()
    {
        foreach (string text in Texts)
        {
            var (nav, _, tree) = Facade(text);
            Assert.That(nav.Count(nav.Root), Is.EqualTo(text.Length), $"Count(Root) of '{text}'");
            Assert.That(tree.CountOccurrences(""), Is.EqualTo(nav.Count(nav.Root)));
            var positions = new long[Math.Max(1, text.Length)];
            foreach (string pattern in FacadeChecks.Patterns(text, 4).Where(p => p.Length > 0))
            {
                IndexLocus locus = nav.Root;
                bool found = nav.TryExtendRight(locus, pattern.AsSpan(), out locus);
                int expectedCount = tree.CountOccurrences(pattern);
                Assert.That(found, Is.EqualTo(expectedCount > 0), $"'{text}' / '{pattern}'");
                if (!found) continue;
                Assert.That(nav.Count(locus), Is.EqualTo(expectedCount), $"'{text}' / '{pattern}': Count");
                var located = new long[expectedCount];
                Assert.That(nav.Locate(locus, 0, located), Is.EqualTo(expectedCount));
                Assert.That(located.Order().Select(p => (int)p), Is.EqualTo(tree.FindAllOccurrences(pattern).Order()), $"'{text}' / '{pattern}': Locate");
            }
        }
    }

    [Test]
    public void SuffixArrayAccess_IsTheSortedSuffixes_AndRankOfIsItsInverse()
    {
        foreach (string text in Texts)
        {
            var (_, sa, tree) = Facade(text);
            string[] sorted = AnalysisReference.SortedSuffixes(text);
            for (int rank = 1; rank <= text.Length; rank++)
            {
                long position = sa.SuffixAt(rank);
                Assert.That(text.Substring((int)position), Is.EqualTo(sorted[rank - 1]), $"'{text}' rank {rank}");
                Assert.That(sa.RankOf(position), Is.EqualTo(rank), $"'{text}' rank of {position}");
            }
            Assert.That(sa.SuffixAt(0), Is.EqualTo(text.Length), "rank 0 is the empty suffix");
            Assert.That(tree.GetAllSuffixes(), Is.EqualTo(sorted));
        }
    }

    [Test]
    public void CapabilityMatrix_FollowsThePack_ThroughTheFacade()
    {
        var (nav, _, tree) = Facade("mississippi");
        IndexOperation[] needPack = [IndexOperation.Parent, IndexOperation.SuffixLink, IndexOperation.InverseSuffixArrayAccess, IndexOperation.Lcp, IndexOperation.Lca];
        foreach (IndexOperation op in needPack) Assert.That(nav.Support(op), Is.EqualTo(IndexSupport.Degraded), $"{op} before the pack");
        Assert.That(nav.Support(IndexOperation.Count), Is.EqualTo(IndexSupport.Native));
        Assert.That(nav.Support(IndexOperation.ExtendLeft), Is.EqualTo(IndexSupport.Degraded));
        Assert.That(nav.SupportNote(IndexOperation.SuffixLink), Is.Not.Null.And.Contain("O(n)"));
        // the search and analysis members of the int surface that do not need links leave the pack alone
        _ = tree.CountOccurrences("ssi");
        _ = tree.LongestRepeatedSubstring();
        _ = tree.CountDistinctSubstrings();
        _ = tree.FindAllLongestRepeatedSubstrings();
        _ = tree.Traverse_Count();
        foreach (IndexOperation op in needPack) Assert.That(nav.Support(op), Is.EqualTo(IndexSupport.Degraded), $"{op} after core operations");
        // a matching member needs links, so it builds the pack once; the matrix then reads Native
        _ = tree.FindMaximalExactMatches("ssippix", 2);
        foreach (IndexOperation op in needPack) Assert.That(nav.Support(op), Is.EqualTo(IndexSupport.Native), $"{op} after the pack exists");
        Assert.That(((EsaSuffixTree)tree).Index.NavigationBuildCount, Is.EqualTo(1));
    }

    [Test]
    public void PrintTree_BuildsThePackForTheLinkMarks_OnlyOnce()
    {
        EsaSuffixTree tree = (EsaSuffixTree)Tree("abracadabra");
        _ = tree.PrintTree();
        _ = tree.PrintTree();
        Assert.That(tree.Index.NavigationBuildCount, Is.EqualTo(1));
    }

    [Test]
    public void InvalidLoci_AreRejected_ThroughTheFacade()
    {
        var (nav, sa, _) = Facade("banana");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => nav.Count(new IndexLocus(5, 2, 1, 1)));
            Assert.Throws<ArgumentException>(() => nav.IsExplicitNode(new IndexLocus(0, 99, 1, 1)));
            Assert.Throws<ArgumentOutOfRangeException>(() => sa.SuffixAt(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => sa.SuffixAt(7));
            Assert.Throws<ArgumentOutOfRangeException>(() => sa.RankOf(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => sa.LcpAt(99));
            Assert.Throws<ArgumentOutOfRangeException>(() => nav.Support((IndexOperation)99));
        });
    }
}

/// <summary>Small helpers used by the navigation tests.</summary>
internal static class TraverseCountExtensions
{
    /// <summary>Runs <c>Traverse</c> with a counting visitor and returns the number of visited nodes.</summary>
    internal static int Traverse_Count(this ISuffixTree tree)
    {
        var collector = new FacadeChecks.EventCollector();
        tree.Traverse(collector);
        return collector.Events.Count;
    }
}
