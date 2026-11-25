using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_pseudoknot_structure</c> MCP tool (RNA-PKPREDICT-001 / RNA-PKRECURSIVE-001). Expected values: independent
/// ViennaRNA-primitive implementation of the pknotsRG grammar (docs/Evidence/RNA-PKPREDICT-001-Evidence.md), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaPseudoknotStructureTests
{
    [Test]
    public void RnaPseudoknotStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaPseudoknotStructure("GGGGAACCCCAACCCCAAGGGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaPseudoknotStructure(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaPseudoknotStructure(null!));
    }

    [Test]
    public void RnaPseudoknotStructure_Binding_InvokesSuccessfully()
    {
        var knot = AnalysisTools.RnaPseudoknotStructure("GGGGAACCCCAACCCCAAGGGG");
        var twin = AnalysisTools.RnaPseudoknotStructure("GGGGAACCCCAACCCCAAGGGGAAAGGGGAACCCCAACCCCAAGGGG");
        var twinSingle = AnalysisTools.RnaPseudoknotStructure("GGGGAACCCCAACCCCAAGGGGAAAGGGGAACCCCAACCCCAAGGGG", allowNestedKnots: false);
        var hairpin = AnalysisTools.RnaPseudoknotStructure("gggaaaaccc");
        Assert.Multiple(() =>
        {
            Assert.That(knot.DotBracket, Is.EqualTo("((((..[[[[..))))..]]]]"));
            Assert.That(knot.FreeEnergy, Is.EqualTo(-13.56).Within(1e-9));
            Assert.That(knot.HasPseudoknot, Is.True);
            Assert.That(knot.BasePairs, Has.Length.EqualTo(8));
            Assert.That(twin.FreeEnergy, Is.EqualTo(-27.52).Within(1e-9));
            Assert.That(twinSingle.FreeEnergy, Is.EqualTo(-24.46).Within(1e-9));
            Assert.That(hairpin.HasPseudoknot, Is.False);
            Assert.That(hairpin.DotBracket, Is.EqualTo("(((....)))"));
            Assert.That(hairpin.FreeEnergy, Is.EqualTo(AnalysisTools.RnaMfeStructure("gggaaaaccc").MinimumFreeEnergy).Within(1e-9));
        });
    }
}
