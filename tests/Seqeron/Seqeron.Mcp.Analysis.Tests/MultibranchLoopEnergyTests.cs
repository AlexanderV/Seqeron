using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>multibranch_loop_energy</c> MCP tool.
/// Expected values from RnaSecondaryStructure's NNDB affine-model unit test
/// (3-way junction, 6 unpaired -> 9.18), NOT the wrapper output.
/// </summary>
[TestFixture]
public class MultibranchLoopEnergyTests
{
    [Test]
    public void MultibranchLoopEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MultibranchLoopEnergy(3, 6));
    }

    [Test]
    public void MultibranchLoopEnergy_Binding_InvokesSuccessfully()
    {
        // a + b*(unpaired/helices) + c*helices = 9.25 + 0.91*2 + (-0.63)*3 = 9.18.
        var e = AnalysisTools.MultibranchLoopEnergy(3, 6).Energy;
        Assert.That(e, Is.EqualTo(9.18).Within(0.01));

        // With strain and stacking: adds -2.0 stacking + 3.14 strain.
        var strained = AnalysisTools.MultibranchLoopEnergy(3, 1, hasStrain: true, stackingEnergy: -2.0).Energy;
        double expected = System.Math.Round(9.25 + 0.91 * (1.0 / 3.0) + (-0.63) * 3 + (-2.0) + 3.14, 2);
        Assert.That(e, Is.Not.EqualTo(strained));
        Assert.That(strained, Is.EqualTo(expected).Within(0.01));
    }
}
