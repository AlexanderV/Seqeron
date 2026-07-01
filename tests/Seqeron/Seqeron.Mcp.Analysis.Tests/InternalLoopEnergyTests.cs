using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>internal_loop_energy</c> MCP tool.
/// Expected value from RnaSecondaryStructure's NNDB 1x1 int11 unit test
/// (CG/CG with G-G mismatch -> -2.2), NOT the wrapper output.
/// </summary>
[TestFixture]
public class InternalLoopEnergyTests
{
    [Test]
    public void InternalLoopEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.InternalLoopEnergy(1, 1, "C", "G", "C", "G", "G", "G", "G", "G"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.InternalLoopEnergy(1, 1, "", "G", "C", "G", "G", "G", "G", "G"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.InternalLoopEnergy(1, 1, "CC", "G", "C", "G", "G", "G", "G", "G"));
    }

    [Test]
    public void InternalLoopEnergy_Binding_InvokesSuccessfully()
    {
        // 1x1 loop CG/CG with G-G mismatch -> -2.2 (int11 table).
        var e = AnalysisTools.InternalLoopEnergy(1, 1, "C", "G", "C", "G", "G", "G", "G", "G").Energy;
        Assert.That(e, Is.EqualTo(-2.2).Within(0.01));
    }
}
