using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_morfs</c> MCP tool.
/// Expected values from DisorderPredictor's own MoRF unit test
/// (DisorderPredictor_MoRF_Tests: P25-L30-P25 -> one MoRF [29,50] score 0.275934;
/// 40xL -> none), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictMorfsTests
{
    [Test]
    public void PredictMorfs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictMorfs(new string('P', 25) + new string('L', 30) + new string('P', 25)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictMorfs(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictMorfs(null!));
    }

    [Test]
    public void PredictMorfs_Binding_InvokesSuccessfully()
    {
        // A hydrophobic L30 dip within a disordered P flank -> one MoRF at [29,50].
        string dip = new string('P', 25) + new string('L', 30) + new string('P', 25);
        var morfs = AnalysisTools.PredictMorfs(dip).Items;
        Assert.Multiple(() =>
        {
            Assert.That(morfs, Has.Length.EqualTo(1));
            Assert.That(morfs[0].Start, Is.EqualTo(29));
            Assert.That(morfs[0].End, Is.EqualTo(50));
            Assert.That(morfs[0].Score, Is.EqualTo(0.275934).Within(1e-6));
        });

        // A fully ordered poly-L has no IDR and therefore no MoRF.
        var none = AnalysisTools.PredictMorfs(new string('L', 40)).Items;
        Assert.That(none, Is.Empty);
    }
}
