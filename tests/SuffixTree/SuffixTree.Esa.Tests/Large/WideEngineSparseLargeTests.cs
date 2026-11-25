using SuffixTree.Esa.External;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.4, the operation matrix of the index over a container at the boundaries of the stored width (ADR-005 section 3.7 layer 2; category <c>Large</c>, explicit): the sparse shell of <c>a^n</c> (the real ends of every section
/// and the pages that a test patches, the middle a hole) opened through <see cref="EsaIndex.FromFile"/> and read through the long contract (<see cref="EsaSuffixIndex"/>: <see cref="ISuffixIndex"/>, <see cref="ISuffixArrayAccess"/>) at
/// n = 2^31 - 2 (the last text of the 32-bit views of S1), 2^31 - 1 (the first of the wide ones: the S1 constructor said <see cref="NotSupportedException"/>), 2^31, 2^31 + 10 and 2^32 - 4 (the last of the 4-byte width) in both access
/// modes. Written before the long engine (it was red against the S1 engine, which cannot make an index over a text of 2^31 - 1 symbols or more); an off-by-one at a boundary shows as a wrong rank, a wrong count or the wrong exception.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class WideEngineSparseLargeTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;

    /// <summary>The texts at the edges: the last of the 32-bit views, the first of the wide ones, a power of two, a text just above it, the last of the 4-byte width.</summary>
    internal static readonly long[] Boundaries = [TwoTo31 - 2, TwoTo31 - 1, TwoTo31, TwoTo31 + 10, TwoTo32 - 4];

    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    public static IEnumerable<object[]> Cases() => Boundaries.SelectMany(n => Modes.Select(m => new object[] { n, m }));

    private static EsaOpenOptions At(EsaAccessMode access) => new() { Validation = EsaValidationLevel.Header, Access = access };

    /// <summary>The shell of <c>a^n</c> with the real ends and, when <paramref name="window"/> is at least 0, the ranks <c>window - 40 .. window + 40</c> as well.</summary>
    private string Shell(long n, long window = -1)
    {
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        SparseLayout layout = SparseContainer.Plan(spec);
        var patches = new List<SparsePatch>(SparseContainer.EndPatches(f, layout, 4096));
        // a window that the real ends (the first and the last 4096 ranks) already cover needs no patch of its own (patches must not overlap)
        if (window >= 0 && window - 40 >= 4096 && window + 40 <= n - 4096) patches.AddRange(SparseContainer.WindowPatches(f, layout, window - 40, window + 40));
        return WriteContainer("wide", spec, patches);
    }

    private static string A(int k) => new('a', k);

    [TestCaseSource(nameof(Cases))]
    public void TheIndexOverTheContainer_IsMade_AndAnswersContains_AtTheEndsOfTheText(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(index.Contains(A(1)), Is.True);
        Assert.That(index.Contains(A(4)), Is.True, "a^4 occurs n - 3 times");
        Assert.That(index.Contains(A(3000)), Is.True);
        Assert.That(index.Contains("b"), Is.False);
        Assert.That(index.Contains("ab"), Is.False);
        Assert.That(index.Contains("aaab"), Is.False);
    }

    [TestCaseSource(nameof(Cases))]
    public void TheStatisticsOfTheHeader_AreTheNumbersOfTheText_At64Bits(long n, EsaAccessMode access)
    {
        ClosedFormText f = new PowerOfA(n);
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaStatistics st = file.Statistics;
        Assert.That(st.Length, Is.EqualTo(n));
        Assert.That(st.NodeCount, Is.EqualTo(f.NodeCount));
        Assert.That(st.InternalCount, Is.EqualTo(f.InternalCount));
        Assert.That(st.LeafCount, Is.EqualTo(n));
        Assert.That(st.MaxDepth, Is.EqualTo(n));
        Assert.That(st.LrsDepth, Is.EqualTo(n - 1), "LRS_DEPTH of a^n");
        Assert.That(st.LrsRank, Is.EqualTo(f.LrsRank));
    }

    [TestCaseSource(nameof(Cases))]
    public void TheContractIndex_ExtendsRightAndCounts_WithRanksAndCountsAbove2To31(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        Assert.That(idx.TextLength, Is.EqualTo(n));
        Assert.That(idx.Root, Is.EqualTo(new IndexLocus(0, n, 0, 0)));
        Assert.That(idx.Count(idx.Root), Is.EqualTo(n));
        IndexLocus locus = idx.Root;
        for (int k = 1; k <= 40; k++)
        {
            Assert.That(idx.TryExtendRight(locus, 'a', out IndexLocus next), Is.True, $"a^{k}");
            locus = next;
            Assert.That(locus, Is.EqualTo(new IndexLocus(k, n, k, k)), $"a^{k} is the node [{k}..{n}]");
            Assert.That(idx.Count(locus), Is.EqualTo(n - k + 1), $"count of a^{k}");
            Assert.That(idx.IsExplicitNode(locus), Is.True);
        }
        Assert.That(idx.TryExtendRight(idx.Root, A(3000).AsSpan(), out IndexLocus far), Is.True);
        Assert.That(far, Is.EqualTo(new IndexLocus(3000, n, 3000, 3000)));
        Assert.That(idx.Count(far), Is.EqualTo(n - 2999));
        Assert.That(idx.TryExtendRight(locus, 'b', out _), Is.False);
        Assert.That(idx.TryExtendRight(idx.Root, "aab".AsSpan(), out _), Is.False);
    }

    [TestCaseSource(nameof(Cases))]
    public void TheContractIndex_ReadsTheSuffixArrayAndTheLcpArray_AtTheEnds(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        foreach (long r in new long[] { 0, 1, 2, 3, 4000, n - 4000, n - 2, n - 1, n })
        {
            Assert.That(idx.SuffixAt(r), Is.EqualTo(n - r), $"SA[{r}]");
            Assert.That(idx.LcpAt(r), Is.EqualTo(r <= 1 ? 0 : r - 1), $"LCP[{r}]");
        }
        Assert.That(() => idx.SuffixAt(n + 1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(idx.Text.Length, Is.EqualTo(n));
        Assert.That(idx.Text.SymbolAt(0), Is.EqualTo((int)'a'));
        Assert.That(idx.Text.SymbolAt(n - 1), Is.EqualTo((int)'a'));
        Assert.That(idx.Text.SymbolAt(n), Is.EqualTo(-1));
    }

    [TestCaseSource(nameof(Cases))]
    public void TheContractIndex_ListsTheChildrenOfTheNodesAtTheStart_AsTheClosedFormSays(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        var children = new IndexLocus[4];
        Assert.That(idx.GetChildren(idx.Root, children), Is.EqualTo(2));
        Assert.That(children[0], Is.EqualTo(new IndexLocus(0, 0, 1, 1)), "the leaf of the empty suffix");
        Assert.That(children[1], Is.EqualTo(new IndexLocus(1, n, 1, 1)), "a");
        idx.TryExtendRight(idx.Root, A(4).AsSpan(), out IndexLocus a4);
        Assert.That(idx.GetChildren(a4, children), Is.EqualTo(2));
        Assert.That(children[0], Is.EqualTo(new IndexLocus(4, 4, 5, 5)), "the leaf of the suffix a^4 (its terminator is the fifth symbol)");
        Assert.That(children[1], Is.EqualTo(new IndexLocus(5, n, 5, 5)));
        Assert.That(idx.TryChild(a4, 'a', out IndexLocus child), Is.True);
        Assert.That(child, Is.EqualTo(new IndexLocus(5, n, 5, 5)));
        Assert.That(idx.TryChild(a4, 'b', out _), Is.False);
    }

    private sealed class CountingVisitor : IWideTreeVisitor
    {
        internal int Events { get; private set; }

        public void VisitNode(long startIndex, long endIndex, long leafCount, long childCount, long depth) => Events++;

        public void EnterBranch(int key) => Events++;

        public void ExitBranch() => Events++;
    }

    [TestCaseSource(nameof(Cases))]
    public void TheSearchSideOfTheIndex_CountsAndRanksAt64Bits_AndTheThirtyTwoBitShapesRefuseByName(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(index.LongCount("aaaa"), Is.EqualTo(n - 3));
        Assert.That(index.LongCount(A(3000)), Is.EqualTo(n - 2999));
        Assert.That(index.LongCount("b"), Is.EqualTo(0));
        Assert.That(index.LongCount(string.Empty), Is.EqualTo(n));
        Assert.That(index.TryMatch("aaaa".AsSpan(), out long lo, out long hi), Is.True);
        Assert.That((lo, hi), Is.EqualTo((4L, n)));
        Assert.That(index.TryMatch("ab".AsSpan(), out lo, out hi), Is.False);
        Assert.That((lo, hi), Is.EqualTo((0L, -1L)));
        if (n > int.MaxValue)
        {
            NotSupportedException e = Assert.Throws<NotSupportedException>(() => index.CountOccurrences("aaaa"))!;
            Assert.That(e.Message, Does.Contain("LongCount"));
            Assert.Throws<NotSupportedException>(() => index.TryMatch("aaaa".AsSpan(), out int _, out int _));
            Assert.Throws<NotSupportedException>(() => index.CountOccurrences(string.Empty));
            Assert.That(index.CountOccurrences("b"), Is.EqualTo(0), "an absent pattern has no count to narrow");
            Assert.That(index.TryMatch("b".AsSpan(), out int l32, out int h32), Is.False);
            Assert.That((l32, h32), Is.EqualTo((0, -1)));
        }
        else
        {
            Assert.That((long)index.CountOccurrences("aaaa"), Is.EqualTo(n - 3));
        }
    }

    [TestCaseSource(nameof(Cases))]
    public void TheLongestRepeat_IsReportedAt64Bits_AndTheStringOfAnOverlongRepeatRefusesByName(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(index.LrsLength, Is.EqualTo(n - 1), "a^n: the longest repeat is a^(n - 1)");
        Assert.That(index.LrsPosition, Is.EqualTo(1), "SA[n - 1] = 1: the rank of the header is n - 1");
        if (n - 1 > int.MaxValue)
        {
            NotSupportedException e = Assert.Throws<NotSupportedException>(() => index.LongestRepeatedSubstring())!;
            Assert.That(e.Message, Does.Contain("LrsLength"));
            Assert.Throws<NotSupportedException>(() => index.LongestRepeatedSubstringMemory());
        }
    }

    [TestCaseSource(nameof(Cases))]
    public void TheTraversalOfATextWhoseTableOfOccurrencesCannotExist_IsRefusedByName_BeforeTheFirstEvent(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(access));
        EsaIndex index = EsaIndex.FromFile(file);
        var visitor = new CountingVisitor();
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => index.Traverse(visitor))!;
        Assert.That(e.Message, Does.Contain("traversal"));
        Assert.That(visitor.Events, Is.Zero, "the refusal comes before any event");
        if (n > int.MaxValue) Assert.Throws<NotSupportedException>(() => index.GetAllSuffixes());
    }

    /// <summary>The page of <c>Locate</c> that straddles the rank 2^31 (and the pages at the ends): the positions are <c>n - rank</c>, the first of the page is at a rank below 2^31 and the last above it.</summary>
    [TestCase(TwoTo31 + 10, EsaAccessMode.Checked)]
    [TestCase(TwoTo31 + 10, EsaAccessMode.Trusted)]
    [TestCase(TwoTo32 - 4, EsaAccessMode.Checked)]
    [TestCase(TwoTo32 - 4, EsaAccessMode.Trusted)]
    public void Locate_PagesAcrossTwoToThe31_AndReadsPositionsAbove2To31(long n, EsaAccessMode access)
    {
        using SesaFile file = SesaFile.Open(Shell(n, window: TwoTo31), At(access));
        var idx = new EsaSuffixIndex(EsaIndex.FromFile(file));
        idx.TryExtendRight(idx.Root, A(4).AsSpan(), out IndexLocus a4); // the ranks 4 .. n
        long skip = TwoTo31 - 4 - 8; // the page starts at the rank 2^31 - 8 and ends at 2^31 + 7
        var page = new long[16];
        Assert.That(idx.Locate(a4, skip, page), Is.EqualTo(16));
        for (int i = 0; i < page.Length; i++) Assert.That(page[i], Is.EqualTo(n - (TwoTo31 - 8 + i)), $"element {i}: the suffix at the rank {TwoTo31 - 8 + i}");
        // the first page and the last page of the node
        var head = new long[3];
        Assert.That(idx.Locate(a4, 0, head), Is.EqualTo(3));
        Assert.That(head, Is.EqualTo(new long[] { n - 4, n - 5, n - 6 }));
        var tail = new long[10];
        Assert.That(idx.Locate(a4, n - 4 - 2, tail), Is.EqualTo(3), "the last three ranks of the node: n - 2, n - 1, n");
        Assert.That(tail.Take(3), Is.EqualTo(new long[] { 2, 1, 0 }));
        Assert.That(idx.Locate(a4, n - 3, tail), Is.EqualTo(0), "skip beyond the interval");
        Assert.That(() => idx.Locate(a4, -1, tail), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [TestCase(TwoTo31 - 2)]
    [TestCase(TwoTo31 - 1)]
    [TestCase(TwoTo31)]
    [TestCase(TwoTo32 - 4)]
    public void TheIndexOverTheContainer_AtTheBoundary_RefusesTheInt32SurfaceByName_OrServesIt(long n)
    {
        using SesaFile file = SesaFile.Open(Shell(n), At(EsaAccessMode.Checked));
        Assert.That(file.IsNarrow, Is.EqualTo(n < int.MaxValue), "the 32-bit loop of S1 serves a text of fewer than 2^31 - 1 symbols and no other");
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(new EsaSuffixIndex(index).TextLength, Is.EqualTo(n));
        if (n > int.MaxValue)
        {
            NotSupportedException e = Assert.Throws<NotSupportedException>(() => _ = index.Length)!;
            Assert.That(e.Message, Does.Contain("EsaIndex.Length").And.Contain("TextLength"), "the bridge names the member and the 64-bit one");
        }
        else
        {
            Assert.That((long)index.Length, Is.EqualTo(n));
        }
    }
}
