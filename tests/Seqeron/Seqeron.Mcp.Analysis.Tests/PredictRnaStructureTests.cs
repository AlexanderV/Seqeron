using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_rna_structure</c> MCP tool.
/// Expected values from RnaSecondaryStructure's own tests (GGGAAAACCC hairpin ->
/// dot-bracket "(((....)))"; empty sequence -> empty structure), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictRnaStructureTests
{
    [Test]
    public void PredictRnaStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictRnaStructure("GGGAAAACCC"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictRnaStructure(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictRnaStructure(null!));
    }

    [Test]
    public void PredictRnaStructure_Binding_InvokesSuccessfully()
    {
        // GGGAAAACCC folds into a single 3bp hairpin -> "(((....)))".
        var s = AnalysisTools.PredictRnaStructure("GGGAAAACCC", 3, 4, 4);
        Assert.Multiple(() =>
        {
            Assert.That(s.Sequence, Is.EqualTo("GGGAAAACCC"));
            Assert.That(s.DotBracket, Is.EqualTo("(((....)))"));
            Assert.That(s.BasePairs, Is.Not.Empty);
            Assert.That(s.StemLoops, Has.Length.EqualTo(1));
            Assert.That(s.MinimumFreeEnergy, Is.EqualTo(-2.02).Within(1e-9)); // RNA.eval_structure (= RNA.fold)
        });
        // RNA-MFE-001: a destabilising stem-loop (+1.91 in context) is not selected; reported energy is the
        // free energy of the returned structure (RNA.eval_structure), not a sum of isolated stem-loops.
        var open = AnalysisTools.PredictRnaStructure("ACAAAUACACAGGUUAAU");
        var two = AnalysisTools.PredictRnaStructure("GGGAAACCCAGGGAAACCC");
        Assert.Multiple(() =>
        {
            Assert.That(open.DotBracket, Is.EqualTo(".................."));
            Assert.That(open.MinimumFreeEnergy, Is.EqualTo(0.0));
            Assert.That(two.DotBracket, Is.EqualTo("(((...))).(((...)))"));
            Assert.That(two.MinimumFreeEnergy, Is.EqualTo(-4.14).Within(1e-9));
        });
    }
}
