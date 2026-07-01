using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_chou_fasman</c> MCP tool.
/// Expected values are the published Chou-Fasman (1978) propensities for A (1.42,0.83,0.66)
/// and E (1.51,0.37,0.74), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictChouFasmanTests
{
    [Test]
    public void PredictChouFasman_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictChouFasman("AE", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictChouFasman("", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictChouFasman(null!, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictChouFasman("AE", 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictChouFasman("AE", -1));
    }

    [Test]
    public void PredictChouFasman_Binding_InvokesSuccessfully()
    {
        // Single Alanine window -> its published (Pa, Pb, Pt).
        var a = AnalysisTools.PredictChouFasman("A", 1).Items;
        Assert.Multiple(() =>
        {
            Assert.That(a, Has.Length.EqualTo(1));
            Assert.That(a[0].Helix, Is.EqualTo(1.42).Within(1e-9));
            Assert.That(a[0].Sheet, Is.EqualTo(0.83).Within(1e-9));
            Assert.That(a[0].Turn, Is.EqualTo(0.66).Within(1e-9));
        });

        // "AE" window 1 -> A then E.
        var ae = AnalysisTools.PredictChouFasman("AE", 1).Items;
        Assert.Multiple(() =>
        {
            Assert.That(ae, Has.Length.EqualTo(2));
            Assert.That(ae[1].Helix, Is.EqualTo(1.51).Within(1e-9));
            Assert.That(ae[1].Sheet, Is.EqualTo(0.37).Within(1e-9));
            Assert.That(ae[1].Turn, Is.EqualTo(0.74).Within(1e-9));
        });

        // Window larger than sequence -> empty.
        var empty = AnalysisTools.PredictChouFasman("A", 5).Items;
        Assert.That(empty, Is.Empty);
    }
}
