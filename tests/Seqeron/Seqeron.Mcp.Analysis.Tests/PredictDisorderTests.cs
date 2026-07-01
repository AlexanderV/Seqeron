using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_disorder</c> MCP tool.
/// Expected values from DisorderPredictor's own unit test
/// (DisorderPredictor_DisorderedRegion_Tests: 30xP -> one region [0,29] fully disordered;
/// 30xW -> no regions), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictDisorderTests
{
    [Test]
    public void PredictDisorder_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictDisorder(new string('P', 30)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictDisorder(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictDisorder(null!));
    }

    [Test]
    public void PredictDisorder_Binding_InvokesSuccessfully()
    {
        // 30xP: most disorder-promoting (normalized 1.0) -> one full-length region.
        var disordered = AnalysisTools.PredictDisorder(new string('P', 30), minRegionLength: 5);
        Assert.Multiple(() =>
        {
            Assert.That(disordered.DisorderedRegions, Has.Length.EqualTo(1));
            Assert.That(disordered.DisorderedRegions[0].Start, Is.EqualTo(0));
            Assert.That(disordered.DisorderedRegions[0].End, Is.EqualTo(29));
            Assert.That(disordered.OverallDisorderContent, Is.EqualTo(1.0).Within(1e-9));
        });

        // 30xW: most order-promoting (normalized 0.0) -> no disordered regions.
        var ordered = AnalysisTools.PredictDisorder(new string('W', 30), minRegionLength: 5);
        Assert.Multiple(() =>
        {
            Assert.That(ordered.DisorderedRegions, Is.Empty);
            Assert.That(ordered.OverallDisorderContent, Is.EqualTo(0.0).Within(1e-9));
        });
    }
}
