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
}
