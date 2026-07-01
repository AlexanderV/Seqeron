using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>hairpin_loop_energy</c> MCP tool.
/// Expected values from RnaSecondaryStructure's NNDB hairpin unit test
/// (RnaSecondaryStructure_HairpinEnergy_Tests: AAAAAA/A-U -> 4.6; GAAAG/A-U -> 4.1),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class HairpinLoopEnergyTests
{
    [Test]
    public void HairpinLoopEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.HairpinLoopEnergy("AAAAAA", "A", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.HairpinLoopEnergy("AAAAAA", "", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.HairpinLoopEnergy("AAAAAA", "AA", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.HairpinLoopEnergy("AAAAAA", "A", null!));
    }

    [Test]
    public void HairpinLoopEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // 6-nt loop AAAAAA closed by A-U -> init(5.4) + tm(AAAU=-0.8) = 4.6.
            Assert.That(AnalysisTools.HairpinLoopEnergy("AAAAAA", "A", "U").Energy, Is.EqualTo(4.6).Within(0.01));
            // 3-nt loop AAC closed by C-G -> special triloop 6.8.
            Assert.That(AnalysisTools.HairpinLoopEnergy("AAC", "C", "G").Energy, Is.EqualTo(6.8).Within(0.01));
        });
    }
}
