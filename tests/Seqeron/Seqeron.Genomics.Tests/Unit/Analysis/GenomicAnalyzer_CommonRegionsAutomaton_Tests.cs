using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// GENOMIC-COMMON-001 (audit R14-4a, B09 F291): the right-maximal common regions are computed from the matching statistics of the suffix automaton of the reverse of the first sequence
/// (O(n + m), no suffix tree) instead of the suffix-tree MEM stream; <see cref="GenomicAnalyzer.FindCommonRegionSpans(DnaSequence, DnaSequence, int, StrandSelection, CancellationToken)"/> returns the
/// same regions as coordinates only. Measured before (poly-A n against n): 2 500 1.45 s, 5 000 5.97 s, 10 000 24.3 s, 20 000 99.6 s / 1.5 GB; (AC)n 10 000 12.5 s.
/// The reference is an independent brute force: ms(i) by naive extension with string.Contains, distinct strings in order of their first start in the query, the lowest occurrence by string.IndexOf.
/// </summary>
[TestFixture]
[Category("GENOMIC-COMMON-001")]
public class GenomicAnalyzer_CommonRegionsAutomaton_Tests
{
    private static List<(string Sequence, int First, int Second)> BruteForce(string text, string query, int minLength)
    {
        var list = new List<(string, int, int)>();
        var seen = new HashSet<string>();
        minLength = Math.Max(1, minLength);
        for (int i = 0; i < query.Length; i++)
        {
            int len = 0;
            while (i + len < query.Length && text.Contains(query.AsSpan(i, len + 1), StringComparison.Ordinal))
                len++;
            if (len < minLength)
                continue;
            string match = query.Substring(i, len);
            if (seen.Add(match))
                list.Add((match, text.IndexOf(match, StringComparison.Ordinal), i));
        }

        return list;
    }

    private static string Random(Random rnd, int n, int alphabet)
        => new([.. Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(alphabet)])]);

    [Test]
    public void DocumentedExample_ACGTACGT_TTACGTGG_Min3_IsTACGT_ACGT_CGT()
    {
        var r = GenomicAnalyzer.FindCommonRegions(new DnaSequence("ACGTACGT"), new DnaSequence("TTACGTGG"), 3).Select(x => (x.Sequence, x.PositionInFirst, x.PositionInSecond)).ToList();
        Assert.That(r, Is.EqualTo(new[] { ("TACGT", 3, 1), ("ACGT", 0, 2), ("CGT", 1, 3) }));
    }

    [Test]
    public void FindCommonRegions_EqualsBruteForce_OnRandomAndLowComplexityPairs()
    {
        var rnd = new Random(14005);
        for (int t = 0; t < 400; t++)
        {
            int alphabet = new[] { 1, 2, 2, 4 }[rnd.Next(4)];
            string a = Random(rnd, rnd.Next(0, 40), alphabet);
            string b = Random(rnd, rnd.Next(0, 40), alphabet);
            int min = rnd.Next(-1, 6);
            if (a.Length == 0 || b.Length == 0)
                continue;
            var expected = BruteForce(a, b, min);
            var got = GenomicAnalyzer.FindCommonRegions(new DnaSequence(a), new DnaSequence(b), min).Select(x => (x.Sequence, x.PositionInFirst, x.PositionInSecond)).ToList();
            Assert.That(got, Is.EqualTo(expected), $"{a} / {b} / {min}");
            var spans = GenomicAnalyzer.FindCommonRegionSpans(new DnaSequence(a), new DnaSequence(b), min);
            Assert.That(spans.Select(s => (a.Substring(s.PositionInFirst, s.Length), s.PositionInFirst, s.PositionInSecond)).ToList(), Is.EqualTo(expected), "the coordinates-only result is the same list");
            Assert.That(spans.Any(s => s.IsReverseComplement), Is.False);
            var strings = GenomicAnalyzer.FindCommonRegions(a, b, min).Select(x => (x.Sequence, x.PositionInFirst, x.PositionInSecond)).ToList();
            Assert.That(strings, Is.EqualTo(expected), "string overload");
        }
    }

    [Test]
    public void ReverseStrand_IsTheForwardSearchOfTheReverseComplement_WithTheMummerPositionConvention()
    {
        var rnd = new Random(14006);
        for (int t = 0; t < 120; t++)
        {
            string a = Random(rnd, rnd.Next(1, 35), rnd.Next(2) == 0 ? 2 : 4);
            string b = Random(rnd, rnd.Next(1, 35), rnd.Next(2) == 0 ? 2 : 4);
            string rc = DnaSequence.GetReverseComplementString(b);
            var expected = BruteForce(a, rc, 2).Select(x => (x.Sequence, x.First, b.Length - x.Second - x.Sequence.Length)).ToList();
            var regions = GenomicAnalyzer.FindCommonRegions(new DnaSequence(a), new DnaSequence(b), 2, CommonRegionMode.RightMaximal, StrandSelection.ReverseComplement);
            Assert.That(regions.Select(x => (x.Sequence, x.PositionInFirst, x.PositionInSecond)).ToList(), Is.EqualTo(expected), $"{a} / {b}");
            Assert.That(regions.All(x => x.IsReverseComplement), Is.True);
            var spans = GenomicAnalyzer.FindCommonRegionSpans(new DnaSequence(a), new DnaSequence(b), 2, StrandSelection.ReverseComplement);
            Assert.That(spans.Select(s => (a.Substring(s.PositionInFirst, s.Length), s.PositionInFirst, s.PositionInSecond)).ToList(), Is.EqualTo(expected));
            // both strands: forward first, then the reverse-complement ones
            var both = GenomicAnalyzer.FindCommonRegionSpans(new DnaSequence(a), new DnaSequence(b), 2, StrandSelection.Forward | StrandSelection.ReverseComplement);
            Assert.That(both.Count, Is.EqualTo(BruteForce(a, b, 2).Count + expected.Count));
        }
    }

    [Test]
    public void StringOverload_NRunsAndLowerCase_NeverMatch_AndEndRegions()
    {
        // only A/C/G/T match, an N ends a region (MUMmer -n); lower case is folded
        var spans = GenomicAnalyzer.FindCommonRegionSpans("acgtNNacgt", "ACGTACGTNN", 3);
        var text = "ACGTNNACGT";
        var query = "ACGTACGTNN";
        var expected = BruteForce(text.Replace('N', '#'), query.Replace('N', '$'), 3);
        Assert.That(spans.Select(s => (text.Substring(s.PositionInFirst, s.Length), s.PositionInFirst, s.PositionInSecond)).ToList(), Is.EqualTo(expected));
        Assert.That(GenomicAnalyzer.FindCommonRegionSpans("", "ACGT", 1), Is.Empty);
        Assert.That(GenomicAnalyzer.FindCommonRegionSpans("ACGT", "", 1), Is.Empty);
    }

    [Test]
    public void PolyA_And_Dinucleotide_RunInLinearTime_AndBuildNoSuffixTree()
    {
        foreach (string unit in new[] { "A", "AC" })
        {
            var a = new DnaSequence(string.Concat(Enumerable.Repeat(unit, 20_000 / unit.Length)));
            var b = new DnaSequence(string.Concat(Enumerable.Repeat(unit, 20_000 / unit.Length)));
            long before = GC.GetAllocatedBytesForCurrentThread();
            var sw = System.Diagnostics.Stopwatch.StartNew();
            var spans = GenomicAnalyzer.FindCommonRegionSpans(a, b, 1);
            sw.Stop();
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
            Assert.That(a.HasSuffixTree, Is.False, "no suffix tree");
            Assert.That(allocated, Is.LessThan(60L * 1024 * 1024), $"{unit}: {allocated / 1048576} MB (the suffix-tree engine took 1.5 GB at this size)");
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(10), $"{unit}: {sw.Elapsed.TotalSeconds:F1} s (99.6 s before)");
            // poly-A: one region per length 1 .. n; the first occurrence of A^k is 0
            if (unit == "A")
            {
                Assert.That(spans.Count, Is.EqualTo(20_000));
                Assert.That(spans.All(s => s.PositionInFirst == 0), Is.True);
                Assert.That(spans.Select(s => s.Length).OrderBy(x => x), Is.EqualTo(Enumerable.Range(1, 20_000)));
                Assert.That(spans[0], Is.EqualTo(new CommonRegionSpan(0, 0, 20_000)));
            }
        }
    }

    [Test]
    public void PolyA_StringRegions_AreTheSameAsTheSpans()
    {
        var a = new DnaSequence(new string('A', 800));
        var regions = GenomicAnalyzer.FindCommonRegions(a, a, 1).ToList();
        var spans = GenomicAnalyzer.FindCommonRegionSpans(a, a, 1);
        Assert.That(regions.Count, Is.EqualTo(800));
        Assert.That(regions.Select(r => (r.Length, r.PositionInFirst, r.PositionInSecond)), Is.EqualTo(spans.Select(s => (s.Length, s.PositionInFirst, s.PositionInSecond))));
    }

    [Test]
    public void Cancellation_IsObserved()
    {
        var a = new DnaSequence(new string('A', 50_000));
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => GenomicAnalyzer.FindCommonRegionSpans(a, a, 1, cancellationToken: cts.Token));
        Assert.Throws<OperationCanceledException>(() => GenomicAnalyzer.FindCommonRegions(a, a, 1, cts.Token).ToList());
    }

    [Test]
    public void NullArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindCommonRegionSpans((DnaSequence)null!, new DnaSequence("A"), 1));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindCommonRegionSpans("A", (string)null!, 1));
    }
}
