using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>can_pair</c> MCP tool.
/// Expected results follow RNA-PAIR-001 / RnaSecondaryStructure_CanPair_Tests:
/// A-U and G-C canonical, G-U wobble; everything else false.
/// </summary>
[TestFixture]
public class CanPairTests
{
    [Test]
    public void CanPair_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CanPair("A", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CanPair("", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CanPair("A", null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CanPair("AU", "U"));
    }

    [Test]
    public void CanPair_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Watson-Crick.
            Assert.That(AnalysisTools.CanPair("A", "U").Result, Is.True);
            Assert.That(AnalysisTools.CanPair("U", "A").Result, Is.True);
            Assert.That(AnalysisTools.CanPair("G", "C").Result, Is.True);
            Assert.That(AnalysisTools.CanPair("C", "G").Result, Is.True);
            // Wobble.
            Assert.That(AnalysisTools.CanPair("G", "U").Result, Is.True);
            Assert.That(AnalysisTools.CanPair("U", "G").Result, Is.True);
            // Non-pairs.
            Assert.That(AnalysisTools.CanPair("A", "A").Result, Is.False);
            Assert.That(AnalysisTools.CanPair("A", "G").Result, Is.False);
            Assert.That(AnalysisTools.CanPair("C", "U").Result, Is.False);
            // Case-insensitive.
            Assert.That(AnalysisTools.CanPair("g", "u").Result, Is.True);
        });
    }
}
