using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gap that the audit of work package 1b.9c showed in the comparer of the leaf sort (<c>External/Sort/SuffixOrder.cs</c>, mutant C38 of <c>analysis/mutation/1b7d_sort_suffixorder_listed.json</c>).
/// </summary>
[TestFixture]
public class SaSortStatsGap1b9cTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    [Test]
    public void ALeafOfEqualKeys_IsSortedByTheTiedComparison_AndTheSortCountsEveryOneOfThem()
    {
        // C38 (TiedComparer is the full comparison): a text of one symbol is one leaf (no partition, whose splitters use the tied comparison directly) in which every key ties, so the comparer of the leaf sort is the only
        // source of tied comparisons, and every comparison that it makes is one
        SaCase c = SaCase.Of(StringFamilies.AllSame(200));
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader text = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        SortOptions options = SaTiny.Sort(_dir.Sub("scratch"), 1, 1) with { LeafCapacity = 1000 };
        var sink = new MemorySaSink(c.N + 1);
        SortStats stats = BucketSorter.Run(text, options, sink, ranks);
        Assert.That(sink.Data, Is.EqualTo(c.Sa));
        Assert.That((stats.Leaves, stats.Partitions), Is.EqualTo((1L, 0L)), "the root is the only leaf");
        Assert.That(stats.Work.TiedCompares, Is.GreaterThan(0), "every key of the leaf ties");
        Assert.That(stats.Work.TiedCompares, Is.EqualTo(stats.Work.Compares), "the leaf sort compares tied suffixes only");
    }
}
