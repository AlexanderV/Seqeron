using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>is_disorder_promoting</c> MCP tool.
/// Expected values taken from Dunker's disorder-promoting set {A,R,G,Q,S,P,E,K}
/// (DisorderPredictor.DisorderPromotingSet). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class IsDisorderPromotingTests
{
    [Test]
    public void IsDisorderPromoting_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.IsDisorderPromoting("A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.IsDisorderPromoting(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.IsDisorderPromoting(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.IsDisorderPromoting("AR"));
    }

    [Test]
    public void IsDisorderPromoting_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Disorder-promoting {A,R,G,Q,S,P,E,K}
            Assert.That(AnalysisTools.IsDisorderPromoting("A").Result, Is.True);
            Assert.That(AnalysisTools.IsDisorderPromoting("P").Result, Is.True);
            Assert.That(AnalysisTools.IsDisorderPromoting("K").Result, Is.True);
            Assert.That(AnalysisTools.IsDisorderPromoting("k").Result, Is.True); // case-insensitive
            // Order-promoting / ambiguous are NOT promoting.
            Assert.That(AnalysisTools.IsDisorderPromoting("W").Result, Is.False);
            Assert.That(AnalysisTools.IsDisorderPromoting("C").Result, Is.False);
            Assert.That(AnalysisTools.IsDisorderPromoting("H").Result, Is.False); // ambiguous
        });
    }
}
