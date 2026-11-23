using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// GENOMIC-REPEAT-001 (audit R14-4b, B09 F292): the repeat finders compute the size of their output from the suffix tree (leaf counts x label lengths) before they materialise anything and refuse a larger one with the
/// typed <see cref="RepeatOutputLimitException"/> (an <see cref="ArgumentException"/>): poly-A of length 20 000 took 30.5 s / 1.92 GB for 19 999 repeats with about 200 M occurrences, 50 000 about 12 GB (reachable from
/// MCP find_repeats). The estimate is locked to the real output on seeded sequences; <see cref="RepeatOutputLimits"/> chooses the limit.
/// </summary>
[TestFixture]
[Category("GENOMIC-REPEAT-001")]
public class GenomicAnalyzer_RepeatOutputLimits_Tests
{
    private static readonly RepeatOutputLimits Tiny = new(1);

    private static (long Repeats, long Occurrences, long Characters) Totals(IEnumerable<RepeatInfo> repeats)
    {
        long r = 0, o = 0, c = 0;
        foreach (var x in repeats)
        {
            r++;
            o += x.Positions.Count;
            c += x.Sequence.Length;
        }

        return (r, o, c);
    }

    private static string Random(Random rnd, int n, int alphabet) => new([.. Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(alphabet)])]);

    [Test]
    public void TheEstimate_EqualsTheRealOutput_ForEveryMode()
    {
        var rnd = new Random(14008);
        for (int t = 0; t < 60; t++)
        {
            string s = Random(rnd, rnd.Next(2, 120), new[] { 1, 2, 4 }[rnd.Next(3)]);
            int min = rnd.Next(0, 5);
            foreach (bool maximal in new[] { false, true })
            {
                var dna = new DnaSequence(s);
                var real = Totals(maximal ? GenomicAnalyzer.FindMaximalRepeats(dna, min, RepeatOutputLimits.Unlimited) : GenomicAnalyzer.FindRepeats(dna, min, RepeatOutputLimits.Unlimited));
                if (real.Repeats == 0)
                {
                    Assert.DoesNotThrow(() => (maximal ? GenomicAnalyzer.FindMaximalRepeats(dna, min, Tiny) : GenomicAnalyzer.FindRepeats(dna, min, Tiny)).ToList(), "nothing to materialise");
                    continue;
                }

                var ex = Assert.Throws<RepeatOutputLimitException>(() => (maximal ? GenomicAnalyzer.FindMaximalRepeats(dna, min, Tiny) : GenomicAnalyzer.FindRepeats(dna, min, Tiny)).ToList(), $"{s} {min} {maximal}");
                Assert.That((ex!.Repeats, ex.Occurrences, ex.Characters), Is.EqualTo(real), $"{s} {min} {maximal}");
                Assert.That(ex.EstimatedBytes, Is.EqualTo(2 * real.Characters + 4 * real.Occurrences + 96 * real.Repeats));
            }
        }
    }

    [Test]
    public void PolyA_20000_IsRefusedByDefault_BeforeAnythingIsMaterialised()
    {
        var dna = new DnaSequence(new string('A', 20_000));
        long before = GC.GetAllocatedBytesForCurrentThread();
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var ex = Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(dna, 1).ToList());
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(ex!.Repeats, Is.EqualTo(19_999));
        Assert.That(ex.Occurrences, Is.EqualTo(20_000L * 20_001 / 2 - 1), "A^k (k = 1 .. n - 1) occurs n - k + 1 times: about 200 M occurrences");
        Assert.That(ex.EstimatedBytes, Is.GreaterThan(RepeatOutputLimits.DefaultMaxBytes));
        Assert.That(ex.Message, Does.Contain("exceeds the limit of 512 MB").And.Contain("minLength"));
        Assert.That(allocated, Is.LessThan(300L * 1024 * 1024), $"{allocated / 1048576} MB allocated (1.92 GB before)");
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(20), "30.5 s before");
        // the start of the text is its own left context, so every A^k of poly-A is a maximal repeat as well
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindMaximalRepeats(dna, 1).ToList());
        // a random sequence of the same length is far below the limit
        var rnd = new Random(14009);
        var random = new DnaSequence(Random(rnd, 20_000, 4));
        Assert.DoesNotThrow(() => GenomicAnalyzer.FindRepeats(random, 1).Count());
    }

    [Test]
    public void TheLimit_IsTheMemoryEstimate_AndMaxOccurrencesIsSeparate()
    {
        var dna = new DnaSequence(new string('A', 300));
        var ex = Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(dna, 1, Tiny).ToList());
        long estimate = ex!.EstimatedBytes;
        Assert.DoesNotThrow(() => GenomicAnalyzer.FindRepeats(dna, 1, new RepeatOutputLimits(estimate)).ToList(), "exactly the estimate passes");
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(dna, 1, new RepeatOutputLimits(estimate - 1)).ToList());
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(dna, 1, new RepeatOutputLimits(long.MaxValue, ex.Occurrences - 1)).ToList());
        Assert.DoesNotThrow(() => GenomicAnalyzer.FindRepeats(dna, 1, new RepeatOutputLimits(long.MaxValue, ex.Occurrences)).ToList());
        Assert.That(ex.Message, Does.Contain("1 MB").Or.Contain("0 MB"));
        Assert.That(ex.Limits, Is.EqualTo(Tiny));
        Assert.That(RepeatOutputLimits.Default.MaxBytes, Is.EqualTo(1L << 29));
        Assert.That(RepeatOutputLimits.Unlimited.MaxBytes, Is.EqualTo(long.MaxValue));
    }

    [Test]
    public void StringAndBothStrandOverloads_AreGuardedToo_AndUnlimitedReturnsTheSameAsTheDefault()
    {
        string s = new('A', 400);
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(s, 1, Tiny).ToList());
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(s, 1, "ACGT", Tiny).ToList());
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindMaximalRepeats("ACGTACGTACGTAAAAAAAAAAAAAAAA", 1, Tiny).ToList());
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeatsBothStrands(s, 1, limits: Tiny));
        Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeatsBothStrands(new DnaSequence(s), 1, false, default, Tiny));
        var small = new DnaSequence("ACGTACGTAAACCCGGGTTTAAACCC");
        Assert.That(Totals(GenomicAnalyzer.FindRepeats(small, 2, RepeatOutputLimits.Unlimited)), Is.EqualTo(Totals(GenomicAnalyzer.FindRepeats(small, 2))));
        Assert.That(GenomicAnalyzer.FindRepeatsBothStrands(small, 2, true, default, RepeatOutputLimits.Unlimited).Count, Is.EqualTo(GenomicAnalyzer.FindRepeatsBothStrands(small, 2, true).Count));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindRepeats(small, 2, (RepeatOutputLimits)null!));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindMaximalRepeats(s, 2, (RepeatOutputLimits)null!));
    }

    [Test]
    public void TheException_IsAnArgumentException_WithTheNumbers()
    {
        var ex = Assert.Throws<RepeatOutputLimitException>(() => GenomicAnalyzer.FindRepeats(new DnaSequence("AAAAAAAAAA"), 1, Tiny).ToList());
        Assert.That(ex, Is.InstanceOf<ArgumentException>());
        Assert.That((ex!.Repeats, ex.Occurrences, ex.Characters), Is.EqualTo((9L, 54L, 45L)));      // A^k, k = 1 .. 9, n - k + 1 occurrences each
    }
}
