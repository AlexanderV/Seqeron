using SuffixTree.Esa.Tests.Analysis;
using System.Runtime;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// The navigation pack (ISA plus the LCP pyramid; work package 1b.4, ADR-001 R-SL): correctness of ISA, its size, and above all the lazy life cycle: built on first use, once, thread-safe,
/// never by the core operations (search, count, locate, Traverse, the native analysis of 1b.3), reported as <see cref="IndexSupport.Degraded"/> until it exists.
/// </summary>
[TestFixture]
public class NavigationPackTests
{
    [Test]
    public void Isa_IsTheInverseOfTheSuffixArray()
    {
        foreach (string text in AnalysisCorpus.Edge().Concat(ExhaustiveStrings.All(2, 7)).Concat(ExhaustiveStrings.All(3, 4)))
        {
            EsaIndex index = EsaIndex.Build(text);
            EsaNavigationPack pack = EsaNavigationPack.Build(index);
            int n = text.Length;
            Assert.That(pack.IsaArray.Length, Is.EqualTo(n + 1), $"'{text}'");
            for (int rank = 0; rank <= n; rank++)
            {
                if (pack.IsaArray[index.SuffixArray[rank]] != rank) Assert.Fail($"'{text}': ISA[SA[{rank}]] != {rank}");
                if (pack.Isa(index.SuffixArray[rank]) != rank) Assert.Fail($"'{text}': Isa(SA[{rank}]) != {rank}");
            }
            Assert.That(pack.Isa(n), Is.EqualTo(0), $"'{text}': the empty suffix has rank 0");
            Assert.That(pack.Index, Is.SameAs(index));
            Assert.That(pack.Root, Is.EqualTo(new EsaNode(0, n, 0, 0, n, true)));
            Assert.That(pack.Pyramid.Length, Is.EqualTo(n + 1));
        }
    }

    [Test]
    public void Isa_LargeRandomText_IsAPermutation()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(50_000, 4));
        EsaNavigationPack pack = index.Navigation;
        var seen = new bool[index.Length + 1];
        for (int p = 0; p <= index.Length; p++)
        {
            int r = pack.Isa(p);
            Assert.That(seen[r], Is.False);
            seen[r] = true;
            Assert.That(index.SuffixArray[r], Is.EqualTo(p));
        }
    }

    [Test]
    public void Bytes_FourPerRankPlusAboutAThirtyThirdOfThat()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(200_000, 5));
        EsaNavigationPack pack = index.Navigation;
        Assert.That(pack.PayloadBytes, Is.EqualTo((4L * (index.Length + 1)) + pack.Pyramid.LevelBytes));
        double perChar = pack.PayloadBytes / (double)index.Length;
        Assert.That(perChar, Is.InRange(4.1, 4.15), "ISA 4.0 + pyramid about 0.13 bytes per character at fanout 32");
        Assert.That(pack.Pyramid.Fanout, Is.EqualTo(32));
    }

    [Test]
    public void Life_NotBuiltByCoreOperations()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(3000, 6));
        Assert.That(index.NavigationBuilt, Is.False);
        Assert.That(index.NavigationBuildCount, Is.EqualTo(0));
        Assert.That(index.SuffixLinkSupport, Is.EqualTo(IndexSupport.Degraded), "the first suffix-link user pays O(n)");

        string pattern = index.Text.Substring(100, 12);
        Assert.That(index.Contains(pattern));
        Assert.That(index.CountOccurrences(pattern), Is.GreaterThanOrEqualTo(1));
        Assert.That(index.FindAllOccurrences(pattern), Is.Not.Empty);
        index.TryMatch(pattern, out int lo, out int hi);
        Assert.That(index.Locate(lo, hi, 0, new int[4]), Is.GreaterThan(0));
        index.Traverse(new TraverseRecorder());
        _ = index.LongestRepeatedSubstring();
        _ = index.FindAllLongestRepeatedSubstrings();
        _ = index.CountDistinctSubstrings();
        _ = index.CountDistinctSubstringsByLength(10);
        _ = index.GetAllSuffixes();
        _ = index.EnumerateSuffixes().Count();
        Assert.That(index.NavigationBuilt, Is.False, "the core operations never touch the pack");
        Assert.That(index.SuffixLinkSupport, Is.EqualTo(IndexSupport.Degraded));
    }

    [Test]
    public void Life_BuiltOnFirstUse_Once()
    {
        EsaIndex index = EsaIndex.Build("mississippi");
        EsaNavigationPack first = index.Navigation;
        Assert.That(index.NavigationBuilt);
        Assert.That(index.NavigationBuildCount, Is.EqualTo(1));
        Assert.That(index.SuffixLinkSupport, Is.EqualTo(IndexSupport.Native));
        Assert.That(index.Navigation, Is.SameAs(first));
        EsaNavigator nav = index.CreateNavigator();
        _ = nav.GetSuffixLink(nav.Root);
        Assert.That(index.NavigationBuildCount, Is.EqualTo(1), "creating a navigator reuses the pack");

        EsaIndex other = EsaIndex.Build("mississippi");
        _ = other.CreateNavigator();
        Assert.That(other.NavigationBuilt, Is.True, "CreateNavigator builds the pack when it is missing");
        Assert.That(other.NavigationBuildCount, Is.EqualTo(1));
    }

    [Test]
    public void Life_ConcurrentFirstUse_BuildsOnce()
    {
        for (int round = 0; round < 12; round++)
        {
            EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(2000 + round, (ulong)round));
            const int Threads = 8;
            var start = new Barrier(Threads);
            var packs = new EsaNavigationPack[Threads];
            var threads = Enumerable.Range(0, Threads).Select(t => new Thread(() =>
            {
                start.SignalAndWait();
                packs[t] = index.Navigation;
            })).ToArray();
            foreach (Thread thread in threads) thread.Start();
            foreach (Thread thread in threads) thread.Join();
            Assert.That(index.NavigationBuildCount, Is.EqualTo(1), $"round {round}");
            Assert.That(packs.All(p => ReferenceEquals(p, packs[0])), Is.True, $"round {round}: one instance");
        }
    }

    [Test]
    public void Life_ConcurrentReaders_AgreeWithTheSequentialWalk()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.Fibonacci(1200));
        var expected = NavigatorWalk.Walk(index, false, "sequential");
        var results = new List<TreeRecord>[6];
        Parallel.For(0, results.Length, i => results[i] = NavigatorWalk.Walk(index, false, $"thread {i}"));
        foreach (var r in results) Assert.That(r, Is.EqualTo(expected));
    }

    [Test]
    public void Pins_SurviveCompactingCollections()
    {
        EsaIndex index = EsaIndex.Build(StringFamilies.RandomDna(1500, 8));
        EsaNavigationPack pack = index.Navigation;
        var nodes = SuffixLinkChecks.AllNodes(index.CreateNavigator(), index.Text).Where(t => !t.Node.IsLeaf).Select(t => t.Node).Take(300).ToList();
        EsaNode[] before = nodes.Select(node => pack.SuffixLink(node)).ToArray();
        var garbage = new List<byte[]>();
        for (int round = 0; round < 4; round++)
        {
            for (int i = 0; i < 200; i++) garbage.Add(new byte[1000 + i]);
            garbage.RemoveAll(_ => true);
            GCSettings.LargeObjectHeapCompactionMode = GCLargeObjectHeapCompactionMode.CompactOnce;
            GC.Collect(2, GCCollectionMode.Forced, blocking: true, compacting: true);
        }
        Assert.That(nodes.Select(node => pack.SuffixLink(node)).ToArray(), Is.EqualTo(before));
    }

    [Test]
    public void Contract_ArgumentChecks_AndTheEmptyText()
    {
        Assert.Throws<ArgumentNullException>(() => EsaNavigationPack.Build(null!));
        EsaIndex index = EsaIndex.Build("abc");
        foreach (int bad in new[] { 0, 1, 3, 33 })
            Assert.Throws<ArgumentOutOfRangeException>(() => EsaNavigationPack.Build(index, bad), $"fanout {bad}");
        EsaIndex empty = EsaIndex.Build("");
        EsaNavigationPack pack = empty.Navigation;
        Assert.That(pack.IsaArray, Is.EqualTo(new[] { 0 }));
        Assert.That(pack.Root, Is.EqualTo(new EsaNode(0, 0, 0, 0, 0, true)));
        EsaNavigator nav = empty.CreateNavigator();
        Assert.That(nav.GetSuffixLink(nav.Root), Is.EqualTo(nav.Root));
        Assert.DoesNotThrow(pack.Dispose);
        Assert.DoesNotThrow(pack.Dispose);
    }

    [Test]
    public void SuffixLink_OfTheNullNode_Throws()
    {
        EsaNavigationPack pack = EsaIndex.Build("abc").Navigation;
        Assert.Throws<ArgumentException>(() => pack.SuffixLink(EsaNode.Null));
        EsaNavigator nav = EsaIndex.Build("abc").CreateNavigator();
        Assert.Throws<ArgumentException>(() => nav.GetSuffixLink(nav.NullNode));
    }

    [Test]
    public void Borders_RejectsEndpointsOutsideTheArrayOrOutOfOrder()
    {
        using LcpPyramid pyramid = LcpPyramid.Build(CompactLcp.FromDense(new[] { 0, 0, 1, 2, 0 }), 2);
        Assert.Throws<ArgumentOutOfRangeException>(() => pyramid.Borders(3, 2, 1, out _, out _, out _, out _), "b before a");
        Assert.Throws<ArgumentOutOfRangeException>(() => pyramid.Borders(5, 5, 1, out _, out _, out _, out _), "a at the length");
        Assert.Throws<ArgumentOutOfRangeException>(() => pyramid.Borders(1, 5, 1, out _, out _, out _, out _), "b at the length");
        Assert.Throws<ArgumentOutOfRangeException>(() => pyramid.Borders(-1, 2, 1, out _, out _, out _, out _), "a negative");
        pyramid.Borders(2, 3, 1, out int left, out int right, out int leftValue, out int rightValue);
        Assert.That((left, right, leftValue, rightValue), Is.EqualTo((1, 4, 0, 0)), "ranks 2..3 hold LCP >= 1; the borders are LCP[1] = 0 and LCP[4] = 0");
        pyramid.Borders(2, 4, 1, out left, out right, out leftValue, out rightValue);
        Assert.That((left, right, leftValue, rightValue), Is.EqualTo((1, 5, 0, -1)), "no border to the right: the length and -1");
    }
}
