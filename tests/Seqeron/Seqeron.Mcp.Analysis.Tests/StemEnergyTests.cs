using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>stem_energy</c> MCP tool.
/// Expected values from RnaSecondaryStructure's NNDB stacking unit test
/// (single pair -> 0; GG/CC two-pair stack -> -3.26), NOT the wrapper output.
/// </summary>
[TestFixture]
public class StemEnergyTests
{
    private static BasePairItem[] TwoPairsGgCc() => new[]
    {
        new BasePairItem(0, 5, 'G', 'C', "WatsonCrick"),
        new BasePairItem(1, 4, 'G', 'C', "WatsonCrick")
    };

    [Test]
    public void StemEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.StemEnergy("GGAACC", TwoPairsGgCc()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.StemEnergy("", TwoPairsGgCc()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.StemEnergy("GGAACC", null!));
    }

    [Test]
    public void StemEnergy_Binding_InvokesSuccessfully()
    {
        // Single base pair -> no stacking -> 0.
        var single = AnalysisTools.StemEnergy("GAAAC", new[]
        {
            new BasePairItem(0, 4, 'G', 'C', "WatsonCrick")
        }).Energy;
        Assert.That(single, Is.EqualTo(0.0).Within(0.01));

        // GG/CC stacking of two consecutive WC pairs -> -3.26 (no terminal penalty).
        var stack = AnalysisTools.StemEnergy("dummy", TwoPairsGgCc()).Energy;
        Assert.That(stack, Is.EqualTo(-3.26).Within(0.01));
    }
}
