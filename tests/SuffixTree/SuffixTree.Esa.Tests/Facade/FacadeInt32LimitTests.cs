using System.Collections;
using SuffixTree.Esa;
using SuffixTree.Persistent;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The int-limit guards of the facade (work package 1b.6; ADR-001 sections 7 and 8; precedent: the Int32LimitGuardTests of the Ukkonen and persistent suites): an oversized text is
/// refused with a clear exception before any content is read, the node count above <see cref="int.MaxValue"/> is refused by the facade and never wraps, and none of it allocates
/// gigabytes: the oversized texts are fakes that report a length and refuse to serve content, the node count is a synthetic statistics value.
/// </summary>
[TestFixture]
public class FacadeInt32LimitTests
{
    private const int LargestAcceptedLength = 1_073_741_823;
    private const int SmallestRejectedLength = 1_073_741_824;

    /// <summary>Reports a configurable length and throws if anybody touches the content.</summary>
    private sealed class HugeLengthTextSource(int length) : ITextSource
    {
        public int Reads { get; private set; }

        public int Length { get; } = length;

        public char this[int index] => throw Refuse();

        public string Substring(int start, int length) => throw Refuse();

        public ReadOnlySpan<char> Slice(int start, int length) => throw Refuse();

        public IEnumerator<char> GetEnumerator() => throw Refuse();

        IEnumerator IEnumerable.GetEnumerator() => throw Refuse();

        private NotSupportedException Refuse()
        {
            Reads++;
            return new NotSupportedException("content must not be read");
        }
    }

    // ------------------------------------------------------------ text length, through the factory

    [Test]
    public void Create_TextOneBeyondTheLimit_ThrowsArgumentException_BeforeReadingContent()
    {
        var source = new HugeLengthTextSource(SmallestRejectedLength);
        var ex = Assert.Throws<ArgumentException>(() => EsaSuffixTreeFactory.Create(source))!;
        Assert.Multiple(() =>
        {
            Assert.That(ex.ParamName, Is.EqualTo("text"));
            Assert.That(ex.Message, Does.Contain("1073741824").And.Contain("1073741823"));
            Assert.That(source.Reads, Is.EqualTo(0), "the length is checked before the content is touched");
        });
    }

    [Test]
    public void Create_TextOfIntMaxValue_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() => EsaSuffixTreeFactory.Create(new HugeLengthTextSource(int.MaxValue)));
        Assert.Throws<ArgumentException>(() => EsaSuffixTreeFactory.CreateEsa(new HugeLengthTextSource(int.MaxValue)));
    }

    [Test]
    public void Create_TextExactlyAtTheLimit_PassesTheLengthGuard_AndAllocatesNoGigabytes()
    {
        var source = new HugeLengthTextSource(LargestAcceptedLength);
        long before = GC.GetAllocatedBytesForCurrentThread();
        Assert.Throws<NotSupportedException>(() => EsaSuffixTreeFactory.Create(source), "the fake refuses to serve content, which is the first thing a legal length reaches");
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.Multiple(() =>
        {
            Assert.That(source.Reads, Is.GreaterThan(0), "a legal length proceeds to the content");
            Assert.That(allocated, Is.LessThan(64L * 1024 * 1024), "the guard and the first read must not allocate the gigabytes of a full buffer");
        });
    }

    [Test]
    public void MaxTextLength_IsTheLargestLengthWhoseWorstCaseNodeCountFitsInt32()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SuffixTreeLimits.MaxTextLength, Is.EqualTo(LargestAcceptedLength));
            Assert.That(2L * SuffixTreeLimits.MaxTextLength + 1, Is.LessThanOrEqualTo(int.MaxValue));
            Assert.That(2L * (SuffixTreeLimits.MaxTextLength + 1L) + 1, Is.GreaterThan(int.MaxValue));
        });
    }

    [Test]
    public void TheLimitException_IsTheOneOfTheOtherBackends()
    {
        // Ukkonen: ArgumentException (ParamName "value"); v6: ArgumentException (ParamName "text"); the facade: ArgumentException (ParamName "text")
        var uk = Assert.Throws<ArgumentException>(() => global::SuffixTree.SuffixTree.Build(new HugeLengthTextSource(SmallestRejectedLength)))!;
        var pe = Assert.Throws<ArgumentException>(() => PersistentSuffixTreeFactory.Create(new HugeLengthTextSource(SmallestRejectedLength)))!;
        var esa = Assert.Throws<ArgumentException>(() => EsaSuffixTreeFactory.Create(new HugeLengthTextSource(SmallestRejectedLength)))!;
        Assert.Multiple(() =>
        {
            Assert.That(esa.GetType(), Is.EqualTo(uk.GetType()));
            Assert.That(esa.GetType(), Is.EqualTo(pe.GetType()));
            Assert.That(pe.ParamName, Is.EqualTo("text"));
            Assert.That(esa.ParamName, Is.EqualTo(pe.ParamName));
            Assert.That(esa.Message, Does.Contain("1073741824"));
        });
    }

    // ------------------------------------------------------------ the guard logic of ADR-001 section 8 (container-opened texts), without any allocation

    [Test]
    public void EnsureTextLength_AcceptsUpToIntMaxValueMinusOne_AndRefusesBeyond()
    {
        Assert.DoesNotThrow(() => EsaSuffixTree.EnsureTextLength(0));
        Assert.DoesNotThrow(() => EsaSuffixTree.EnsureTextLength(1));
        Assert.DoesNotThrow(() => EsaSuffixTree.EnsureTextLength(int.MaxValue - 1L));
        var ex = Assert.Throws<NotSupportedException>(() => EsaSuffixTree.EnsureTextLength(int.MaxValue))!;
        Assert.That(ex.Message, Does.Contain("2147483647").And.Contain("long"));
        Assert.Throws<NotSupportedException>(() => EsaSuffixTree.EnsureTextLength(int.MaxValue + 1L));
        Assert.Throws<NotSupportedException>(() => EsaSuffixTree.EnsureTextLength(long.MaxValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => EsaSuffixTree.EnsureTextLength(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => EsaSuffixTree.EnsureTextLength(long.MinValue));
    }

    [Test]
    public void CheckedNodeCount_PassesUpToIntMaxValue_AndRefusesAboveWithAClearException()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EsaSuffixTree.CheckedNodeCount(1), Is.EqualTo(1));
            Assert.That(EsaSuffixTree.CheckedNodeCount(int.MaxValue), Is.EqualTo(int.MaxValue));
            var ex = Assert.Throws<NotSupportedException>(() => EsaSuffixTree.CheckedNodeCount(int.MaxValue + 1L))!;
            Assert.That(ex.Message, Does.Contain("2147483648").And.Contain("2147483647").And.Contain("64-bit"));
            Assert.Throws<NotSupportedException>(() => EsaSuffixTree.CheckedNodeCount(long.MaxValue));
            Assert.Throws<NotSupportedException>(() => EsaSuffixTree.CheckedNodeCount(4_000_000_000L));
        });
    }

    [Test]
    public void NodeCount_AboveInt32_IsRefusedByTheFacade_AndEverythingElseKeepsWorking()
    {
        EsaIndex index = EsaIndex.Build("banana");
        EsaStatistics real = index.Statistics;
        EsaStatistics big = real with { NodeCount = int.MaxValue + 1L };
        using var refused = new EsaSuffixTree(index, big);
        using var justFits = new EsaSuffixTree(index, real with { NodeCount = int.MaxValue });
        using var normal = new EsaSuffixTree(index);
        Assert.Multiple(() =>
        {
            var ex = Assert.Throws<NotSupportedException>(() => _ = refused.NodeCount)!;
            Assert.That(ex.Message, Does.Contain("2147483648"));
            Assert.That(justFits.NodeCount, Is.EqualTo(int.MaxValue));
            Assert.That(normal.NodeCount, Is.EqualTo(11));
            // the rest of the surface does not need the node count
            Assert.That(refused.LeafCount, Is.EqualTo(6));
            Assert.That(refused.MaxDepth, Is.EqualTo(6));
            Assert.That(refused.CountOccurrences("an"), Is.EqualTo(2));
            Assert.That(refused.LongestRepeatedSubstring(), Is.EqualTo("ana"));
        });
    }

    [Test]
    public void TheNodeCountOfTheLargestAcceptedText_FitsInt32_ByTheFormulaOfTheStatistics()
    {
        // worst case: n leaves, the terminator-only leaf, at most n - 1 further internal nodes and the root, i.e. 2n + 1
        long worst = (2L * SuffixTreeLimits.MaxTextLength) + 1;
        Assert.That(EsaSuffixTree.CheckedNodeCount(worst), Is.EqualTo((int)worst));
    }

    // ------------------------------------------------------------ argument checks of the factory

    [Test]
    public void Create_NullText_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => EsaSuffixTreeFactory.Create((string)null!));
        Assert.Throws<ArgumentNullException>(() => EsaSuffixTreeFactory.Create((ITextSource)null!));
        Assert.Throws<ArgumentNullException>(() => EsaSuffixTreeFactory.CreateEsa((string)null!));
        Assert.Throws<ArgumentNullException>(() => EsaSuffixTreeFactory.CreateEsa((ITextSource)null!));
        Assert.Throws<ArgumentNullException>(() => EsaSuffixTree.Build(null!));
    }
}
