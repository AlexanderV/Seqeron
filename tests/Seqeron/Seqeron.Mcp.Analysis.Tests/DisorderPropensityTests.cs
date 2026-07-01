using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>disorder_propensity</c> MCP tool.
/// Expected values taken directly from the TOP-IDP scale (Campen 2008, Table 2) embedded in
/// DisorderPredictor.DisorderPropensity. NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class DisorderPropensityTests
{
    [Test]
    public void DisorderPropensity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DisorderPropensity("P"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DisorderPropensity(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DisorderPropensity(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DisorderPropensity("PR"));
    }

    [Test]
    public void DisorderPropensity_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.DisorderPropensity("P").Propensity, Is.EqualTo(0.987).Within(1e-12));  // max
            Assert.That(AnalysisTools.DisorderPropensity("W").Propensity, Is.EqualTo(-0.884).Within(1e-12)); // min
            Assert.That(AnalysisTools.DisorderPropensity("A").Propensity, Is.EqualTo(0.060).Within(1e-12));
            Assert.That(AnalysisTools.DisorderPropensity("E").Propensity, Is.EqualTo(0.736).Within(1e-12));
            // Case-insensitive, unknown -> 0.
            Assert.That(AnalysisTools.DisorderPropensity("p").Propensity, Is.EqualTo(0.987).Within(1e-12));
            Assert.That(AnalysisTools.DisorderPropensity("Z").Propensity, Is.EqualTo(0.0).Within(1e-12));
        });
    }
}
