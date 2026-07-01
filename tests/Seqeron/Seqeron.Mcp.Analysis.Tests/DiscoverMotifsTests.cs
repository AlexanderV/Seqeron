using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>discover_motifs</c> MCP tool.
/// Expected values computed by hand from MotifFinder.DiscoverMotifs:
/// enrichment = observed / ((N-k+1)/4^k). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class DiscoverMotifsTests
{
    [Test]
    public void DiscoverMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DiscoverMotifs("ATATATAT", 2, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DiscoverMotifs("", 2, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DiscoverMotifs(null!, 2, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DiscoverMotifs("AUGC", 2, 2)); // not DNA
    }

    [Test]
    public void DiscoverMotifs_Binding_InvokesSuccessfully()
    {
        // "ATATATAT" (N=8), k=2, minCount=2: AT x4 at [0,2,4,6]; TA x3 at [1,3,5].
        // windowCount = 8-2+1 = 7; expected = 7/16 = 0.4375.
        //   AT enrichment = 4/0.4375 = 9.142857...; TA enrichment = 3/0.4375 = 6.857142...
        var items = AnalysisTools.DiscoverMotifs("ATATATAT", 2, 2).Items;
        var at = items.Single(m => m.Sequence == "AT");
        var ta = items.Single(m => m.Sequence == "TA");
        Assert.Multiple(() =>
        {
            Assert.That(at.Count, Is.EqualTo(4));
            Assert.That(at.Positions, Is.EqualTo(new[] { 0, 2, 4, 6 }));
            Assert.That(at.Enrichment, Is.EqualTo(4.0 / 0.4375).Within(1e-9));

            Assert.That(ta.Count, Is.EqualTo(3));
            Assert.That(ta.Positions, Is.EqualTo(new[] { 1, 3, 5 }));
            Assert.That(ta.Enrichment, Is.EqualTo(3.0 / 0.4375).Within(1e-9));
        });
    }
}
