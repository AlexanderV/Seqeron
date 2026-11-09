using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Tree;

/// <summary>
/// The stack scan over an LCP stream (work package 1b.1): the events it emits are the lcp-interval tree of the definition. The scan is the part of
/// the child-table builder that the streaming builder of ADR-004 reuses, so it is tested on its own, without any layout.
/// </summary>
[TestFixture]
public class LcpIntervalScanTests
{
    private static ScanLog Scan(string text, out int[] lcp, out int stackPeak)
    {
        int[] sa = NaiveSuffixArray.Build(text);
        lcp = NaiveLcp.Build(text, sa);
        var log = new ScanLog();
        var stream = new ArrayLcpStream(lcp);
        var sink = new LoggingSink(log);
        stackPeak = LcpIntervalScan.Run(text.Length, ref stream, ref sink);
        return log;
    }

    private static void AssertMatchesDefinition(string text, string context)
    {
        ScanLog log = Scan(text, out int[] lcp, out int stackPeak);
        List<BruteInterval> expected = BruteLcpTree.Intervals(lcp);

        var expectedCloses = expected
            .Select(x => (x.Ell, x.Lb, x.Rb, First: x.LIndices[0], x.LastChild))
            .OrderBy(x => x.Lb).ThenByDescending(x => x.Rb).ToList();
        var actualCloses = log.Closes.OrderBy(x => x.Lb).ThenByDescending(x => x.Rb).ToList();
        Assert.That(actualCloses, Is.EqualTo(expectedCloses), context);

        var expectedLIndices = new List<(int, int)>();
        foreach (BruteInterval x in expected)
        {
            int previous = -1;
            foreach (int k in x.LIndices)
            {
                expectedLIndices.Add((previous, k));
                previous = k;
            }
        }
        Assert.That(log.LIndices.OrderBy(x => x.K).ToList(), Is.EqualTo(expectedLIndices.OrderBy(x => x.Item2).ToList()), context);

        if (text.Length > 0)
        {
            // intervals close in post-order: a child before its parent, the root last, every rb is the scan position - 1 (non-decreasing)
            Assert.That(log.Closes[^1].Ell, Is.EqualTo(0), context);
            Assert.That(log.Closes[^1].Lb, Is.EqualTo(0), context);
            Assert.That(log.Closes.Select(x => x.Rb), Is.Ordered, context);
            // the stack holds the intervals that are open at a scan position: pushed at their first l-index, popped after their right boundary
            Assert.That(stackPeak, Is.EqualTo(MaxOpen(expected, text.Length)), context);
        }
        else
        {
            Assert.That(log.Closes, Is.Empty, context);
            Assert.That(log.LIndices, Is.Empty, context);
            Assert.That(stackPeak, Is.EqualTo(0), context);
        }
    }

    private static int MaxOpen(List<BruteInterval> intervals, int n)
    {
        int best = 0;
        for (int i = 1; i <= n; i++)
            best = Math.Max(best, intervals.Count(x => x.LIndices[0] <= i && i <= x.Rb));
        return best;
    }

    [Test]
    public void Banana_EventsMatchTheDefinition()
    {
        AssertMatchesDefinition("banana", "banana");
        AssertMatchesDefinition("mississippi", "mississippi");
    }

    [Test]
    public void EmptyAndTinyTexts()
    {
        foreach (string s in new[] { "", "a", "aa", "ab", "aaa", "aba" })
            AssertMatchesDefinition(s, $"'{s}'");
    }

    [Test]
    public void BananaClosesInPostOrderWithTheExpectedIntervals()
    {
        // SA of banana$ = 6 5 3 1 0 4 2, LCP = 0 0 1 3 0 0 2 : intervals [1..3] ell 1 (a), [2..3] ell 3 (ana), [5..6] ell 2 (na), root [0..6]
        ScanLog log = Scan("banana", out _, out _);
        Assert.That(log.Closes, Is.EqualTo(new[]
        {
            (3, 2, 3, 3, true),    // ana: [2..3], first l-index 3, last child of "a"
            (1, 1, 3, 2, false),   // a: [1..3], first l-index 2, not the last child of the root
            (2, 5, 6, 6, true),    // na: [5..6], last child of the root
            (0, 0, 6, 1, true),    // root
        }));
        Assert.That(log.LIndices, Is.EqualTo(new[] { (-1, 1), (-1, 2), (-1, 3), (1, 4), (4, 5), (-1, 6) }));
    }

    [Test]
    public void AllSmallStrings_EventsMatchTheDefinition()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 9), (3, 6), (4, 4) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                AssertMatchesDefinition(s, $"'{s}'");
    }

    [Test]
    public void AdversarialAndRandomTexts_EventsMatchTheDefinition()
    {
        var rng = new SplitMix64(0x1B01);
        for (int i = 0; i < 40; i++)
        {
            int n = rng.NextInt(250);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 4 == 0))
                AssertMatchesDefinition(text, $"{name} n {n}");
        }
    }

    [Test]
    public void StackPeakIsTheNestingDepth_ChainForAllSameAndLogarithmicForRandom()
    {
        // a^n: LCP = 0,0,1,2,...: every entry opens an interval, nothing closes before the end: n intervals nested + the root -> height n
        ScanLog log = Scan(new string('a', 40), out _, out int peak);
        Assert.That(peak, Is.EqualTo(40));
        Assert.That(log.Closes, Has.Count.EqualTo(40));
        Scan(StringFamilies.RandomDna(500, 5), out _, out int randomPeak);
        Assert.That(randomPeak, Is.LessThan(40));
    }

    [Test]
    public void ArrayAndCompactStreamsDeliverTheSameSequence()
    {
        string text = StringFamilies.Fibonacci(300);
        int[] sa = NaiveSuffixArray.Build(text);
        int[] lcp = NaiveLcp.Build(text, sa);
        CompactLcp compact = CompactLcp.FromDense(lcp);
        var a = new ArrayLcpStream(lcp);
        var c = new CompactLcpStream(compact);
        for (int r = 1; r <= text.Length; r++)
            Assert.That(c.Next(), Is.EqualTo(a.Next()), $"rank {r}");
    }

    [Test]
    public void CorruptStreams_AreRejected()
    {
        // LCP[1] must be 0 (the predecessor of rank 1 is the empty suffix) and no value may be negative: the scan uses the values as stack comparisons
        foreach (int[] lcp in new[] { new[] { 0, 1, 0, 0 }, new[] { 0, 0, -1, 0 }, new[] { 0, 0, 1, -5 } })
        {
            var stream = new ArrayLcpStream(lcp);
            var sink = new LoggingSink(new ScanLog());
            Assert.That(() => LcpIntervalScan.Run(lcp.Length - 1, ref stream, ref sink), Throws.TypeOf<ArgumentException>(), string.Join(",", lcp));
        }
        var good = new ArrayLcpStream([0, 0, 1, 0]);
        var goodSink = new LoggingSink(new ScanLog());
        Assert.That(LcpIntervalScan.Run(3, ref good, ref goodSink), Is.GreaterThan(0));
        var none = new ArrayLcpStream([0]);
        var noneSink = new LoggingSink(new ScanLog());
        Assert.That(() => LcpIntervalScan.Run(-1, ref none, ref noneSink), Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void DeepStacks_GrowPastTheInitialCapacity()
    {
        // a^n needs a stack of n entries; the initial capacity is 64 and doubles
        foreach (int n in new[] { 63, 64, 65, 127, 128, 129, 1000 })
        {
            ScanLog log = Scan(new string('a', n), out _, out int peak);
            Assert.That(peak, Is.EqualTo(n), $"n {n}");
            Assert.That(log.Closes, Has.Count.EqualTo(n), $"n {n}");
        }
    }
}
