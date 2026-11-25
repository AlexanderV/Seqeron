using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_structure_energy</c> MCP tool. Expected values: ViennaRNA 2.7.2
/// <c>eval_structure</c> (dangles = 2, NNDB-precision Turner 2004 parameter file), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaStructureEnergyTests
{
    [Test]
    public void RnaStructureEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaStructureEnergy("GGGAAACCC", "(((...)))"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaStructureEnergy("", "."));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaStructureEnergy("GGGAAACCC", null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaStructureEnergy("GGGAAACCC", "((....))"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaStructureEnergy("GGGAAACCC", "[[[...]]]"));
        // RNA-HAIRPIN-001: hairpin loops < 3 nt are prohibited (NNDB; ViennaRNA eval_structure = INF).
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaStructureEnergy("GGGAACCC", "(((..)))"));
    }

    [Test]
    public void RnaStructureEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.RnaStructureEnergy("GGGAAACCC", "(((...)))").Energy, Is.EqualTo(-1.12).Within(1e-9));
            Assert.That(AnalysisTools.RnaStructureEnergy("GAGGGAAAACCCAGGGAAAACCCAC", "(.(((....))).(((....))).)").Energy,
                Is.EqualTo(-1.54).Within(1e-9));
        });
    }
}
