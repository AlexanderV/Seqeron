using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>reversal_distance</c> MCP tool.
/// Expected values from the breakpoint lower bound d = ceil(breakpoints/2)
/// (Bafna & Pevzner 1998), NOT the wrapper output.
/// </summary>
[TestFixture]
public class ReversalDistanceTests
{
    [Test]
    public void ReversalDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.ReversalDistance(new[] { 1, 2, 3 }, new[] { 3, 2, 1 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ReversalDistance(null!, new[] { 1 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ReversalDistance(new[] { 1 }, null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ReversalDistance(new[] { 1, 2 }, new[] { 1 }));
    }

    [Test]
    public void ReversalDistance_Binding_InvokesSuccessfully()
    {
        // Identical permutations -> 0 breakpoints -> distance 0.
        var same = AnalysisTools.ReversalDistance(new[] { 1, 2, 3, 4 }, new[] { 1, 2, 3, 4 }).Distance;
        Assert.That(same, Is.EqualTo(0));

        // Full reversal -> 2 breakpoints -> ceil(2/2) = 1 reversal.
        var reversed = AnalysisTools.ReversalDistance(new[] { 1, 2, 3, 4 }, new[] { 4, 3, 2, 1 }).Distance;
        Assert.That(reversed, Is.EqualTo(1));
    }
}
