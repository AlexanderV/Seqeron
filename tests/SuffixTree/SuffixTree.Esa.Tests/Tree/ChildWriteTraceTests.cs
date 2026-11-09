using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// The write pattern of the child-table scan (ADR-004 section 3.2: child-table writes must be sequential in one direction or LIFO-local so that a
/// streaming builder can write into a mapped file). <see cref="ChildWriteTrace"/> replays the scan events as the writes of each layout.
/// </summary>
[TestFixture]
public class ChildWriteTraceTests
{
    private static ChildWriteTrace Trace(string text, ChildLayoutKind layout, int pageSlots = 1024, int cachePages = 64)
    {
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        return ChildWriteTrace.Run(layout, text.Length, CompactLcp.FromDense(lcp), pageSlots, cachePages);
    }

    [Test]
    public void PackedLayouts_NeverWriteTheSameSlotTwice()
    {
        // the single-array layouts share one slot per index between the next, up and down fields: a collision would overwrite a field
        var rng = new SplitMix64(0x1B03);
        for (int i = 0; i < 30; i++)
        {
            int n = 1 + rng.NextInt(400);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 4 == 0))
                foreach (ChildLayoutKind layout in new[] { ChildLayoutKind.Packed, ChildLayoutKind.Tagged })
                    Assert.That(Trace(text, layout).DuplicateSlots, Is.EqualTo(0), $"{layout} {name}");
        }
        foreach (string s in ExhaustiveStrings.All(2, 9).Concat(ExhaustiveStrings.All(3, 5)))
            foreach (ChildLayoutKind layout in new[] { ChildLayoutKind.Packed, ChildLayoutKind.Tagged })
                Assert.That(Trace(s, layout).DuplicateSlots, Is.EqualTo(0), $"{layout} '{s}'");
    }

    [Test]
    public void ReferenceLayout_NeverWritesTheSameSlotOfOneArrayTwice()
    {
        foreach (string s in ExhaustiveStrings.All(2, 9).Concat(ExhaustiveStrings.All(3, 5)))
            Assert.That(Trace(s, ChildLayoutKind.Reference3).DuplicateSlots, Is.EqualTo(0), $"'{s}'");
    }

    [Test]
    public void EveryWriteIsAtTheCursorOrBehindIt_AndBehindWritesAreBounded()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(300, 7))
            foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
            {
                ChildWriteTrace t = Trace(text, layout);
                Assert.That(t.Writes, Is.EqualTo(t.AtCursor + t.Behind), $"{layout} {name}");
                Assert.That(t.MaxBackDistance, Is.LessThanOrEqualTo(text.Length + 1), $"{layout} {name}");
                Assert.That(t.Behind == 0, Is.EqualTo(t.MaxBackDistance == 0), $"{layout} {name}");
            }
    }

    [Test]
    public void AllSameText_IsOneBackwardSweepAtTheEnd()
    {
        // a^n: the intervals close from the innermost to the root at the very end, their left boundaries descend by one: every write after the
        // scan reached the end is the previous write minus one (a single sequential backward sweep)
        ChildWriteTrace t = Trace(new string('a', 5000), ChildLayoutKind.Packed, pageSlots: 1024, cachePages: 4);
        Assert.That(t.Behind, Is.GreaterThan(4000));
        Assert.That(t.BackwardRuns, Is.LessThanOrEqualTo(3));
        Assert.That(t.PageMisses, Is.LessThanOrEqualTo(2 * (5001 / 1024 + 2)));
    }

    [Test]
    public void RandomText_WritesStayInASmallWindow()
    {
        ChildWriteTrace t = Trace(StringFamilies.RandomDna(60_000, 3), ChildLayoutKind.Packed, pageSlots: 1024, cachePages: 8);
        Assert.That(t.Writes, Is.GreaterThan(30_000));
        // a few huge intervals close far behind the cursor; everything else is local (a window of 8 pages = 8,192 slots)
        Assert.That(t.PageMisses, Is.LessThan(t.Writes / 20));
    }

    [Test]
    public void CountsAreConsistentWithTheScanEvents()
    {
        string text = "mississippi";
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        var log = new ScanLog();
        var stream = new ArrayLcpStream(lcp);
        var sink = new LoggingSink(log);
        LcpIntervalScan.Run(text.Length, ref stream, ref sink);
        // reference layout: one write per l-index link plus one per interval; packed: the links to the next index are implied and not written
        ChildWriteTrace reference = Trace(text, ChildLayoutKind.Reference3);
        ChildWriteTrace packed = Trace(text, ChildLayoutKind.Packed);
        ChildWriteTrace tagged = Trace(text, ChildLayoutKind.Tagged);
        int links = log.LIndices.Count(x => x.Previous >= 0);
        int implied = log.LIndices.Count(x => x.Previous >= 0 && x.K == x.Previous + 1);
        Assert.That(reference.Writes, Is.EqualTo(links + log.Closes.Count));
        Assert.That(tagged.Writes, Is.EqualTo(links + log.Closes.Count));
        Assert.That(packed.Writes, Is.EqualTo(links - implied + log.Closes.Count));
    }
}
