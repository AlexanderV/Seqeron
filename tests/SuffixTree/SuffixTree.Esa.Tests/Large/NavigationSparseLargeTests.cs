using SuffixTree.Esa.External;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.5, the navigation side over a container at the boundaries of the stored width (ADR-005 section 3.7 layer 2; category <c>Large</c>, explicit): the sparse shell of <c>a^n</c> with PSI and the persisted RMQ pyramid (the real
/// ends of every section, the closed-form levels of the pyramid, and a window around the rank 2^31 or n / 2; the middle is a hole), opened through <see cref="EsaIndex.FromFile"/> and read through the long contract at
/// n = 2^31 - 2, 2^31 - 1, 2^31, 2^31 + 10 and 2^32 - 4 (the last of the 4-byte width) in both access modes. <c>a^n</c> has n - 1 internal nodes <c>a^k</c> = <c>[k..n]</c> in a chain: the link and the parent of <c>a^k</c> are <c>a^(k-1)</c>
/// (the root for k = 1), the Lcp of the ranks i and j is min(i, j), their Lca is <c>a^min(i, j)</c>, and the children of <c>a^k</c> are the leaf of the suffix <c>a^k</c> and <c>a^(k+1)</c>. Written before the engine: against the engine of S2.4
/// (a pack that scans all LCP entries to build its pyramid and holds a 32-bit inverse suffix array) the navigation is red; the 32-bit surface that cannot hold the numbers refuses by name.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class NavigationSparseLargeTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;

    internal static readonly long[] Boundaries = [TwoTo31 - 2, TwoTo31 - 1, TwoTo31, TwoTo31 + 10, TwoTo32 - 4];

    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    public static IEnumerable<object[]> Cases() => Boundaries.SelectMany(n => Modes.Select(m => new object[] { n, m }));

    private static EsaOpenOptions At(EsaAccessMode access) => new() { Validation = EsaValidationLevel.Header, Access = access };

    /// <summary>The rank around which the shell has real entries besides its ends: 2^31 where the text reaches it, else the middle.</summary>
    private static long WindowOf(long n) => n >= TwoTo31 + 100 ? TwoTo31 : n / 2;

    /// <summary>The shell of <c>a^n</c> with PSI and RMQ: real ends of 4096 ranks, the window of <see cref="WindowOf"/> (40 ranks on either side), the closed-form levels around all of them.</summary>
    private string Shell(long n)
    {
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f, psi: true, rmq: true);
        SparseLayout layout = SparseContainer.Plan(spec);
        long window = WindowOf(n);
        var patches = new List<SparsePatch>(SparseContainer.EndPatches(f, layout, 4096));
        patches.AddRange(SparseContainer.WindowPatches(f, layout, window - 40, window + 40));
        patches.AddRange(SparseContainer.RmqPatches(f, layout, (0, 4095), (n - 4095, n), (window - 40, window + 40)));
        return WriteContainer("nav", spec, patches);
    }

    /// <summary>The ranks whose nodes the tests look at: the start, the end and the window.</summary>
    private static long[] Ks(long n)
    {
        long w = WindowOf(n);
        return [.. new long[] { 1, 2, 3, 40, 4000, w - 20, w - 1, w, w + 1, w + 20, n - 4000, n - 40, n - 3, n - 2, n - 1 }.Where(k => k >= 1 && k <= n - 1).Distinct().Order()];
    }

    private static IndexLocus Ak(long n, long k) => new(k, n, k, k);

    [TestCaseSource(nameof(Cases))]
    public void TheSuffixLink_OfTheNodeAk_IsTheNodeAkMinus1_AcrossTheBoundary_WithNoHeap(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaIndex index = EsaIndex.FromFile(file);
        var idx = new EsaSuffixIndex(index);
        foreach (long k in Ks(n))
        {
            Assert.That(idx.TrySuffixLink(Ak(n, k), out IndexLocus link), Is.True, $"a^{k}");
            Assert.That(link, Is.EqualTo(k == 1 ? idx.Root : Ak(n, k - 1)), $"the link of a^{k}");
        }
        Assert.That(index.Navigation.PayloadBytes, Is.EqualTo(0), "PSI and RMQ are read where they lie: no inverse suffix array, no pyramid in memory");
        Assert.That(idx.TrySuffixLink(idx.Root, out _), Is.False);
        Assert.That(idx.TrySuffixLink(new IndexLocus(n, n, n + 1, n + 1), out _), Is.False, "a leaf has no link in the contract");
    }

    [TestCaseSource(nameof(Cases))]
    public void TheParent_OfTheNodeAk_IsTheNodeAkMinus1_AcrossTheBoundary(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        foreach (long k in Ks(n))
        {
            Assert.That(idx.TryParent(Ak(n, k), out IndexLocus parent), Is.True, $"a^{k}");
            Assert.That(parent, Is.EqualTo(k == 1 ? idx.Root : Ak(n, k - 1)), $"the parent of a^{k}");
        }
        Assert.That(idx.TryParent(idx.Root, out _), Is.False);
    }

    [TestCaseSource(nameof(Cases))]
    public void TheLcpAndTheLca_OfRanks_AcrossTheBoundary_AreTheMinimumOfTheTwo(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        long[] ranks = [.. Ks(n).Append(n).Distinct().Order()];
        foreach (long i in ranks)
        {
            foreach (long j in ranks.Where(x => x >= i))
            {
                long m = Math.Min(i, j);
                Assert.That(idx.Lcp(i, j), Is.EqualTo(m), $"Lcp({i}, {j})");
                Assert.That(idx.Lcp(j, i), Is.EqualTo(m), $"Lcp({j}, {i})");
                IndexLocus lca = idx.Lca(i, j);
                if (i == j) Assert.That(lca, Is.EqualTo(new IndexLocus(i, i, i + 1, i + 1)), $"Lca({i}, {i}): the leaf, its depth counts the terminator");
                else Assert.That(lca, Is.EqualTo(m == 0 ? idx.Root : Ak(n, m)), $"Lca({i}, {j})");
            }
        }
    }

    [TestCaseSource(nameof(Cases))]
    public void TheChildren_OfTheNodeAk_AreTheLeafOfItsSuffixAndTheNextNode_AcrossTheBoundary(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        var children = new IndexLocus[4];
        foreach (long k in Ks(n))
        {
            Assert.That(idx.GetChildren(Ak(n, k), children), Is.EqualTo(2), $"a^{k}");
            Assert.That(children[0], Is.EqualTo(new IndexLocus(k, k, k + 1, k + 1)), $"the leaf of the suffix a^{k}");
            Assert.That(children[1], Is.EqualTo(k + 1 <= n - 1 ? Ak(n, k + 1) : new IndexLocus(n, n, n + 1, n + 1)), $"the second child of a^{k}");
            Assert.That(idx.TryChild(Ak(n, k), 'a', out IndexLocus next), Is.True);
            Assert.That(next, Is.EqualTo(children[1]));
            Assert.That(idx.TryChild(Ak(n, k), -1, out IndexLocus leaf), Is.True, "the terminator is the first key");
            Assert.That(leaf, Is.EqualTo(children[0]));
        }
    }

    [TestCaseSource(nameof(Cases))]
    public void TheInverseSuffixArray_AboveThePackHeapLimit_IsUnsupportedByName_NotAnOutOfMemoryError(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => idx.RankOf(5))!;
        Assert.That(e.Message, Does.Contain("RankOf").And.Contain(EsaNavigationPack.PackHeapLimit.ToString(System.Globalization.CultureInfo.InvariantCulture)).And.Contain("ISA"));
        NotSupportedException e2 = Assert.Throws<NotSupportedException>(() => idx.TryExtendLeft(Ak(n, 3), 'a', out _))!;
        Assert.That(e2.Message, Does.Contain("TryExtendLeft").And.Contain(EsaNavigationPack.PackHeapLimit.ToString(System.Globalization.CultureInfo.InvariantCulture)));
        Assert.That(idx.TryExtendLeft(idx.Root, 'a', out IndexLocus a), Is.True, "the empty string extends to a symbol on the right: no inverse suffix array needed");
        Assert.That(a, Is.EqualTo(Ak(n, 1)));
        Assert.That(idx.Lcp(7, 9), Is.EqualTo(7), "Lcp and Lca need no inverse suffix array");
    }
}
