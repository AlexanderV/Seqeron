using System.Reflection;
using SuffixTree.Esa;
using SuffixTree.Indexing;
using SuffixTree.Persistent;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Lifetime of the facade (ADR-001 section 7: disposable like <see cref="PersistentSuffixTree"/>): <c>Dispose</c> is idempotent, every member of a disposed tree throws
/// <see cref="ObjectDisposedException"/> (the property getters, <c>EnumerateSuffixes</c> eagerly, the argument checks after the disposal check), the navigation pack is released, and
/// a tree that a library routine builds and disposes (<see cref="SuffixTreeAlgorithms.FindLongestCommonSubstrings"/>) is unusable afterwards. The persistent tree is in the fixture list as
/// the reference behaviour; the in-memory Ukkonen tree is not disposable and not part of it.
/// </summary>
[TestFixtureSource(typeof(FacadeLifetimeTests), nameof(Fixtures))]
public class FacadeLifetimeTests(SuffixTreeBackend backend)
{
    public static IEnumerable<TestFixtureData> Fixtures()
    {
        foreach (SuffixTreeBackend b in new[] { SuffixTreeBackends.PersistentHeap, SuffixTreeBackends.Esa, SuffixTreeBackends.EsaPacked, SuffixTreeBackends.EsaReference3 })
            yield return new TestFixtureData(b).SetArgDisplayNames(b.Name);
    }

    private ISuffixTree Disposed(string text = "banana")
    {
        ISuffixTree tree = backend.FromString(text);
        ((IDisposable)tree).Dispose();
        return tree;
    }

    [Test]
    public void TheTree_IsDisposable()
    {
        ISuffixTree tree = backend.FromString("banana");
        Assert.That(tree, Is.InstanceOf<IDisposable>());
        ((IDisposable)tree).Dispose();
    }

    [Test]
    public void Dispose_IsIdempotent()
    {
        ISuffixTree tree = backend.FromString("banana");
        var d = (IDisposable)tree;
        Assert.DoesNotThrow(d.Dispose);
        Assert.DoesNotThrow(d.Dispose);
        Assert.DoesNotThrow(d.Dispose);
    }

    [Test]
    public void Dispose_OfAnEmptyTextTree_Works()
    {
        ISuffixTree tree = backend.FromString("");
        Assert.DoesNotThrow(((IDisposable)tree).Dispose);
        Assert.Throws<ObjectDisposedException>(() => tree.Contains(""));
    }

    [Test]
    public void Dispose_FromSeveralThreads_DisposesOnce_WithoutException()
    {
        for (int round = 0; round < 20; round++)
        {
            ISuffixTree tree = backend.FromString("mississippi");
            var d = (IDisposable)tree;
            var errors = new System.Collections.Concurrent.ConcurrentBag<Exception>();
            Parallel.For(0, 8, _ =>
            {
                try { d.Dispose(); }
                catch (Exception ex) { errors.Add(ex); }
            });
            Assert.That(errors, Is.Empty);
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("i"));
        }
    }

    [Test]
    public void UseAfterDispose_EverySearchAnalysisAndDiagnosticsMember_ThrowsObjectDisposedException()
    {
        ISuffixTree tree = Disposed();
        Assert.Multiple(() =>
        {
            Assert.Throws<ObjectDisposedException>(() => _ = tree.Text);
            Assert.Throws<ObjectDisposedException>(() => _ = tree.NodeCount);
            Assert.Throws<ObjectDisposedException>(() => _ = tree.LeafCount);
            Assert.Throws<ObjectDisposedException>(() => _ = tree.MaxDepth);
            Assert.Throws<ObjectDisposedException>(() => _ = tree.IsEmpty);
            Assert.Throws<ObjectDisposedException>(() => tree.Contains("a"));
            Assert.Throws<ObjectDisposedException>(() => tree.Contains("a".AsSpan()));
            Assert.Throws<ObjectDisposedException>(() => tree.Contains(""));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllOccurrences("a"));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllOccurrences("a".AsSpan()));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllOccurrences(""));
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("a"));
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("a".AsSpan()));
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences(""));
            Assert.Throws<ObjectDisposedException>(() => tree.LongestRepeatedSubstring());
            Assert.Throws<ObjectDisposedException>(() => tree.LongestRepeatedSubstringMemory());
            Assert.Throws<ObjectDisposedException>(() => tree.GetAllSuffixes());
            Assert.Throws<ObjectDisposedException>(() => tree.EnumerateSuffixes());
            Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstring("x"));
            Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstring("x".AsSpan()));
            Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstringInfo("x"));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllLongestCommonSubstrings("x"));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllDistinctLongestCommonSubstrings("x"));
            Assert.Throws<ObjectDisposedException>(() => tree.FindExactMatchAnchors("a", 1));
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalExactMatches("a", 1));
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalUniqueMatches("a", 1));
            Assert.Throws<ObjectDisposedException>(() => tree.CountDistinctSubstrings());
            Assert.Throws<ObjectDisposedException>(() => tree.CountDistinctSubstringsByLength(2));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllLongestRepeatedSubstrings());
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalRepeatedPairs(1));
            Assert.Throws<ObjectDisposedException>(() => tree.PrintTree());
            Assert.Throws<ObjectDisposedException>(() => tree.Traverse(new FacadeChecks.EventCollector()));
        });
    }

    [Test]
    public void UseAfterDispose_TheDisposalCheckComesBeforeTheArgumentCheck()
    {
        ISuffixTree tree = Disposed();
        Assert.Multiple(() =>
        {
            Assert.Throws<ObjectDisposedException>(() => tree.Contains((string)null!));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllOccurrences((string)null!));
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences((string)null!));
            Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstring((string)null!));
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalExactMatches(null!, 0));
            Assert.Throws<ObjectDisposedException>(() => tree.CountDistinctSubstringsByLength(0));
            Assert.Throws<ObjectDisposedException>(() => tree.Traverse(null!));
        });
    }

    [Test]
    public void UseAfterDispose_NamesTheTreeInTheException()
    {
        ISuffixTree tree = Disposed();
        var ex = Assert.Throws<ObjectDisposedException>(() => tree.Contains("a"))!;
        Assert.That(ex.ObjectName, Is.EqualTo(tree.GetType().FullName).Or.EqualTo(tree.GetType().Name));
    }

    [Test]
    public void UseAfterDispose_NavigationAndSuffixArrayAccess_ThrowObjectDisposedException()
    {
        ISuffixTree tree = backend.FromString("banana");
        if (tree is not ISuffixTreeNavigation nav || tree is not ISuffixArrayAccess sa)
        {
            Assert.Ignore("the backend has no navigation / suffix-array surface");
            return;
        }
        IndexLocus root = nav.Root;
        bool built = nav.TryChild(root, 'b', out IndexLocus child);
        Assert.That(built, Is.True);
        ((IDisposable)tree).Dispose();
        Assert.Multiple(() =>
        {
            Assert.Throws<ObjectDisposedException>(() => _ = nav.Engine);
            Assert.Throws<ObjectDisposedException>(() => _ = nav.TextLength);
            Assert.Throws<ObjectDisposedException>(() => _ = ((ISuffixIndex)tree).Text);
            Assert.Throws<ObjectDisposedException>(() => _ = nav.Root);
            Assert.Throws<ObjectDisposedException>(() => nav.Support(IndexOperation.Count));
            Assert.Throws<ObjectDisposedException>(() => nav.SupportNote(IndexOperation.Count));
            Assert.Throws<ObjectDisposedException>(() => nav.Count(root));
            Assert.Throws<ObjectDisposedException>(() => nav.TryExtendRight(root, 'a', out _));
            Assert.Throws<ObjectDisposedException>(() => nav.TryExtendRight(root, "an".AsSpan(), out _));
            Assert.Throws<ObjectDisposedException>(() => nav.TryExtendLeft(child, 'a', out _));
            Assert.Throws<ObjectDisposedException>(() => nav.Locate(root, 0, new long[3]));
            Assert.Throws<ObjectDisposedException>(() => nav.IsExplicitNode(root));
            Assert.Throws<ObjectDisposedException>(() => nav.TryParent(child, out _));
            Assert.Throws<ObjectDisposedException>(() => nav.TryChild(root, 'b', out _));
            Assert.Throws<ObjectDisposedException>(() => nav.GetChildren(root, new IndexLocus[8]));
            Assert.Throws<ObjectDisposedException>(() => nav.TrySuffixLink(child, out _));
            Assert.Throws<ObjectDisposedException>(() => nav.EdgeSymbol(child, 0));
            Assert.Throws<ObjectDisposedException>(() => sa.SuffixAt(1));
            Assert.Throws<ObjectDisposedException>(() => sa.RankOf(1));
            Assert.Throws<ObjectDisposedException>(() => sa.LcpAt(1));
            Assert.Throws<ObjectDisposedException>(() => sa.Lcp(1, 2));
            Assert.Throws<ObjectDisposedException>(() => sa.Lca(1, 2));
        });
    }

    [Test]
    public void EnumerateSuffixes_ThrowsAtTheCallNotAtTheFirstMoveNext()
    {
        ISuffixTree tree = Disposed();
        IEnumerable<string>? held = null;
        Assert.Throws<ObjectDisposedException>(() => held = tree.EnumerateSuffixes());
        Assert.That(held, Is.Null);
    }

    [Test]
    public void FindLongestCommonSubstrings_DisposesTheTreesItBuilds()
    {
        var created = new List<ISuffixTree>();
        ISuffixTree Build(string text)
        {
            ISuffixTree tree = backend.FromString(text);
            created.Add(tree);
            return tree;
        }

        var result = SuffixTreeAlgorithms.FindLongestCommonSubstrings(["ACGTAC", "GTACGT", "TACGTA"], 3, Build);
        Assert.That(result, Is.Not.Empty);
        Assert.That(created, Has.Count.GreaterThanOrEqualTo(1));
        foreach (ISuffixTree tree in created)
            Assert.Throws<ObjectDisposedException>(() => tree.Contains("A"));
        created.Clear();
        _ = SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(["abab", "baba"], Build);
        Assert.That(created, Has.Count.GreaterThanOrEqualTo(1));
        foreach (ISuffixTree tree in created)
            Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("a"));
    }

    [Test]
    public void UsingBlock_DisposesAtTheEnd()
    {
        ISuffixTree tree = backend.FromString("abab");
        using (var d = (IDisposable)tree)
        {
            Assert.That(tree.CountOccurrences("ab"), Is.EqualTo(2));
            Assert.That(d, Is.SameAs(tree));
        }
        Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("ab"));
    }

    [Test]
    public void ResultsHandedOutBeforeDispose_StayValid()
    {
        ISuffixTree tree = backend.FromString("abracadabra");
        IReadOnlyList<int> positions = tree.FindAllOccurrences("abra");
        string lrs = tree.LongestRepeatedSubstring();
        IReadOnlyList<string> suffixes = tree.GetAllSuffixes();
        ((IDisposable)tree).Dispose();
        Assert.Multiple(() =>
        {
            Assert.That(positions.Order(), Is.EqualTo(new[] { 0, 7 }));
            Assert.That(lrs, Is.EqualTo("abra"));
            Assert.That(suffixes, Has.Count.EqualTo(11));
        });
    }
}

/// <summary>What <c>Dispose</c> releases in the ESA facade: the pins of the navigation pack, if it was built; a pack that was never built is not built by the disposal.</summary>
[TestFixture]
public class FacadeEsaDisposeTests
{
    private static bool PackHoldsHandle(EsaIndex index)
    {
        EsaNavigationPack pack = index.Navigation;
        FieldInfo field = typeof(EsaNavigationPack).GetField("_sa", BindingFlags.Instance | BindingFlags.NonPublic)!;
        return ((PinnedArray<int>)field.GetValue(pack)!).HoldsHandle;
    }

    [Test]
    public void Dispose_WithoutAPack_DoesNotBuildOne()
    {
        EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa("banana");
        Assert.That(tree.Index.NavigationBuilt, Is.False);
        tree.Dispose();
        Assert.That(tree.Index.NavigationBuilt, Is.False, "the disposal must not pay for the pack");
        Assert.That(tree.Index.NavigationBuildCount, Is.EqualTo(0));
    }

    [Test]
    public void Dispose_AfterAPackWasBuilt_ReleasesThePinOfTheSuffixArray()
    {
        EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa("mississippi");
        _ = tree.FindMaximalExactMatches("ssippix", 2);     // builds the pack
        Assert.That(tree.Index.NavigationBuilt, Is.True);
        Assert.That(PackHoldsHandle(tree.Index), Is.True, "the pack pins the suffix array while it is in use");
        tree.Dispose();
        Assert.That(PackHoldsHandle(tree.Index), Is.False, "Dispose releases the GC handle of the pack");
    }

    [Test]
    public void Dispose_AfterTheNavigationSurfaceBuiltThePack_ReleasesIt()
    {
        EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa("abracadabra");
        _ = tree.RankOf(3);
        Assert.That(tree.Index.NavigationBuilt, Is.True);
        tree.Dispose();
        Assert.That(PackHoldsHandle(tree.Index), Is.False);
    }

    [Test]
    public void TheFacade_DoesNotDisposeTheCallersText()
    {
        string text = "banana";
        EsaSuffixTree tree = EsaSuffixTreeFactory.CreateEsa(text);
        tree.Dispose();
        Assert.That(text, Is.EqualTo("banana"));
    }
}
