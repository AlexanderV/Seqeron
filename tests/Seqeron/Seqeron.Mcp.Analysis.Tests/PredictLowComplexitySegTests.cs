using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_low_complexity_seg</c> MCP tool.
/// Expected values from DisorderPredictor's own SEG unit test
/// (DisorderPredictor_LowComplexity_Tests: 26xQ -> region [0,25] "Q-rich";
/// all-20-AA-twice -> none), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictLowComplexitySegTests
{
    [Test]
    public void PredictLowComplexitySeg_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictLowComplexitySeg(new string('Q', 26)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictLowComplexitySeg(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictLowComplexitySeg(null!));
    }

    [Test]
    public void PredictLowComplexitySeg_Binding_InvokesSuccessfully()
    {
        // 26xQ -> one Q-rich low-complexity region spanning [0,25].
        var polyQ = AnalysisTools.PredictLowComplexitySeg(new string('Q', 26)).Items;
        Assert.Multiple(() =>
        {
            Assert.That(polyQ, Has.Length.EqualTo(1));
            Assert.That(polyQ[0].Start, Is.EqualTo(0));
            Assert.That(polyQ[0].End, Is.EqualTo(25));
            Assert.That(polyQ[0].Type, Is.EqualTo("Q-rich"));
        });

        // All 20 amino acids twice -> maximal complexity -> no region.
        var complex = AnalysisTools.PredictLowComplexitySeg("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY").Items;
        Assert.That(complex, Is.Empty);
    }
}
