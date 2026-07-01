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
        });
    }
}
