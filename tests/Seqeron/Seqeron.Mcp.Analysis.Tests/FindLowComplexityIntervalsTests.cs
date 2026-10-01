using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_low_complexity_intervals</c> MCP tool. Expected = NCBI dustmasker -outfmt interval / lh3/sdust
/// output, locked in SequenceComplexity_CalculateDustScore_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindLowComplexityIntervalsTests
{
    private const string LinkerSeq =
        "ACGTGCATGCAAAAAAAAAAAAAAAAGCTAGCATCGACTGCAGCACACACACACACACAGATCGATCGTACGGTGCATGACAAAAAAAAAAAAACT";

    [Test]
    public void FindLowComplexityIntervals_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindLowComplexityIntervals(LinkerSeq));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityIntervals(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityIntervals(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLowComplexityIntervals("XYZ"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLowComplexityIntervals(LinkerSeq, linker: 33));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLowComplexityIntervals(LinkerSeq, windowSize: 2));
    }

    [Test]
    public void FindLowComplexityIntervals_Binding_DefaultsEqualSdust()
    {
        // lh3/sdust -w 64 -t 20: x 10 26 / x 43 59 / x 81 94.
        var items = AnalysisTools.FindLowComplexityIntervals(LinkerSeq).Items;
        Assert.That(items.Select(i => (i.Start, i.End, i.Length)), Is.EqualTo(new[] { (10, 26, 16), (43, 59, 16), (81, 94, 13) }));
    }

    [TestCase(18, new[] { 10, 59, 81, 94 })]
    [TestCase(32, new[] { 10, 94 })]
    public void FindLowComplexityIntervals_Linker_MatchesDustmasker(int linker, int[] flat)
    {
        var expected = Enumerable.Range(0, flat.Length / 2).Select(i => (flat[2 * i], flat[2 * i + 1])).ToArray();
        var items = AnalysisTools.FindLowComplexityIntervals(LinkerSeq, linker: linker).Items;
        Assert.That(items.Select(i => (i.Start, i.End)), Is.EqualTo(expected));
    }

    [Test]
    public void FindLowComplexityIntervals_NSplitsTheScan_EqualsPerPieceSdust()
    {
        var items = AnalysisTools.FindLowComplexityIntervals("ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA").Items;
        Assert.That(items.Select(i => (i.Start, i.End)), Is.EqualTo(new[] { (6, 18), (23, 38) }));
    }
}
