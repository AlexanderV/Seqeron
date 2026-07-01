using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>minimum_free_energy</c> MCP tool.
/// Expected value from RnaSecondaryStructure's own Turner-2004 manual calculation
/// (RnaSecondaryStructureTests: "GGGAAACCC" -> -1.12 kcal/mol), NOT the wrapper output.
/// </summary>
[TestFixture]
public class MinimumFreeEnergyTests
{
    [Test]
    public void MinimumFreeEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MinimumFreeEnergy("GGGAAACCC"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MinimumFreeEnergy(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MinimumFreeEnergy(null!));
    }

    [Test]
    public void MinimumFreeEnergy_Binding_InvokesSuccessfully()
    {
        // GGGAAACCC: 3 GC pairs (stack -6.52) + 3-nt hairpin init (5.4) = -1.12.
        var mfe = AnalysisTools.MinimumFreeEnergy("GGGAAACCC").Mfe;
        Assert.That(mfe, Is.EqualTo(-1.12).Within(0.01));

        // Poly-A cannot pair -> MFE 0.
        var flat = AnalysisTools.MinimumFreeEnergy("AAAAAAAA").Mfe;
        Assert.That(flat, Is.EqualTo(0.0).Within(1e-9));
    }
}
