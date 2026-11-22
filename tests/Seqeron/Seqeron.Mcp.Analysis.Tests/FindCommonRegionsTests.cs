using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_common_regions</c> MCP tool.
/// Expected results derived from GenomicAnalyzer.FindCommonRegions: for each start position in
/// sequence2 the single longest substring (>= minLength) also occurring in sequence1, keyed by
/// distinct match. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindCommonRegionsTests
{
    [Test]
    public void FindCommonRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindCommonRegions("ACGT", "ACGT", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindCommonRegions("", "ACGT", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindCommonRegions("ACGT", null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindCommonRegions("AUGC", "ACGT", 2)); // not DNA
    }

    [Test]
    public void FindCommonRegions_Binding_InvokesSuccessfully()
    {
        // seq1="ACGT", seq2="ACGT", minLength 2. Per seq2 start:
        //   i=0 longest "ACGT" (len4) at seq1 pos 0; i=1 "CGT" (len3) at pos1; i=2 "GT" (len2) at pos2.
        var items = AnalysisTools.FindCommonRegions("ACGT", "ACGT", 2).Items;
        Assert.That(items, Has.Length.EqualTo(3));
        Assert.Multiple(() =>
        {
            Assert.That(items[0].Sequence, Is.EqualTo("ACGT"));
            Assert.That(items[0].PositionInFirst, Is.EqualTo(0));
            Assert.That(items[0].PositionInSecond, Is.EqualTo(0));
            Assert.That(items[0].Length, Is.EqualTo(4));

            Assert.That(items[1].Sequence, Is.EqualTo("CGT"));
            Assert.That(items[1].PositionInSecond, Is.EqualTo(1));

            Assert.That(items[2].Sequence, Is.EqualTo("GT"));
            Assert.That(items[2].PositionInSecond, Is.EqualTo(2));
        });
    }

    // B09 F15: positionInFirst is the FIRST occurrence in sequence1 (Python str.find oracle), not tree order.
    [Test]
    public void FindCommonRegions_PositionInFirst_IsFirstOccurrence()
    {
        var items = AnalysisTools.FindCommonRegions("AGGAAGAAGGAGGCCGATTACT", "CGTGGAAGGAGGCCGA", 2).Items;
        Assert.That(items, Has.Length.EqualTo(13));
        Assert.Multiple(() =>
        {
            Assert.That(items[^1].Sequence, Is.EqualTo("GA"));
            Assert.That(items[^1].PositionInFirst, Is.EqualTo(2));
            Assert.That(items[^1].PositionInSecond, Is.EqualTo(14));
        });
    }

    [Test]
    public void FindCommonRegions_IncludeSequenceFalse_IsTheSameListAsCoordinatesOnly()
    {
        // R14-4a: coordinates only; the strings are the substrings of the first sequence at positionInFirst
        const string reference = "AGGAAGAAGGAGGCCGATTACT";
        const string query = "CGTGGAAGGAGGCCGA";
        foreach (string strands in new[] { "Forward", "ReverseComplement", "Both" })
        {
            var full = AnalysisTools.FindCommonRegions(reference, query, 2, strands: strands).Items;
            var spans = AnalysisTools.FindCommonRegions(reference, query, 2, strands: strands, includeSequence: false).Items;
            Assert.That(spans.Select(s => s.Sequence), Is.All.Null, strands);
            Assert.That(spans.Select(s => (s.PositionInFirst, s.PositionInSecond, s.Length, s.IsReverseComplement)),
                Is.EqualTo(full.Select(s => (s.PositionInFirst, s.PositionInSecond, s.Length, s.IsReverseComplement))), strands);
            Assert.That(full.Select(s => s.Sequence), Is.All.Not.Null);
        }

        // the other modes materialise and then omit the string
        var mem = AnalysisTools.FindCommonRegions(reference, query, 4, "MaximalExactMatches", includeSequence: false).Items;
        Assert.That(mem, Is.Not.Empty);
        Assert.That(mem.Select(s => s.Sequence), Is.All.Null);
    }

    [Test]
    public void FindCommonRegions_PolyA20000_CoordinatesOnly_IsFast()
    {
        string polyA = new('A', 20_000);
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var items = AnalysisTools.FindCommonRegions(polyA, polyA, 1, includeSequence: false).Items;
        Assert.That(items, Has.Length.EqualTo(20_000));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(20), "99.6 s with the suffix-tree engine (R14-4a)");
    }
}
