using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_transmembrane_helices</c> MCP tool.
/// Expected values from ProteinMotifFinder's own unit test
/// (ProteinMotifFinder_PredictTransmembraneHelices_Tests, D10-L20-D10 -> segment 5..34,
/// peak KD(L)=3.8), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictTransmembraneHelicesTests
{
    [Test]
    public void PredictTransmembraneHelices_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictTransmembraneHelices("LLLLLLLLLLLLLLLLLLLL"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictTransmembraneHelices(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictTransmembraneHelices(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictTransmembraneHelices("LLLL", 0));
    }

    [Test]
    public void PredictTransmembraneHelices_Binding_InvokesSuccessfully()
    {
        // D10 - L20 - D10, window 19, threshold 1.6 -> one segment 5..34, peak 3.8.
        string sequence = new string('D', 10) + new string('L', 20) + new string('D', 10);
        var segments = AnalysisTools.PredictTransmembraneHelices(sequence).Items;
        Assert.Multiple(() =>
        {
            Assert.That(segments, Has.Length.EqualTo(1));
            Assert.That(segments[0].Start, Is.EqualTo(5));
            Assert.That(segments[0].End, Is.EqualTo(34));
            Assert.That(segments[0].Score, Is.EqualTo(3.8).Within(1e-10));
        });
    }
}
