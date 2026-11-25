using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_mfe_structure</c> MCP tool (RNA-MFE-001). Expected values: ViennaRNA 2.7.2 <c>RNA.fold</c>
/// (dangles = 2, NNDB-precision Turner 2004 parameter file), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaMfeStructureTests
{
    [Test]
    public void RnaMfeStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaMfeStructure("GGGAAAACCC"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaMfeStructure(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaMfeStructure(null!));
    }

    [Test]
    public void RnaMfeStructure_Binding_InvokesSuccessfully()
    {
        var hairpin = AnalysisTools.RnaMfeStructure("gggaaaaccc");
        var twoHairpins = AnalysisTools.RnaMfeStructure("GAGGGAAAACCCAGGGAAAACCCAC");
        var minLoop5 = AnalysisTools.RnaMfeStructure("GGGGAAAACCCC", 5);
        Assert.Multiple(() =>
        {
            Assert.That(hairpin.Sequence, Is.EqualTo("GGGAAAACCC"));
            Assert.That(hairpin.DotBracket, Is.EqualTo("(((....)))"));
            Assert.That(hairpin.BasePairs, Has.Length.EqualTo(3));
            Assert.That(hairpin.MinimumFreeEnergy, Is.EqualTo(-2.02).Within(1e-9));
            Assert.That(twoHairpins.DotBracket, Is.EqualTo("..(((....))).(((....)))..")); // RNA.fold
            Assert.That(twoHairpins.MinimumFreeEnergy, Is.EqualTo(-7.04).Within(1e-9));
            Assert.That(twoHairpins.MinimumFreeEnergy,
                Is.EqualTo(AnalysisTools.MinimumFreeEnergy("GAGGGAAAACCCAGGGAAAACCCAC").Mfe).Within(1e-9));
            Assert.That(minLoop5.DotBracket, Is.EqualTo("(((.....)))."));                 // md.min_loop_size = 5
            Assert.That(minLoop5.MinimumFreeEnergy, Is.EqualTo(-4.12).Within(1e-9));
        });
    }
}
