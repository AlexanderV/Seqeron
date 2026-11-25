using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_stem_loops</c> MCP tool.
/// Expected values from RnaSecondaryStructure's own unit test
/// (RnaSecondaryStructureTests.FindStemLoops_SimpleHairpin: "GGGAAAACCC" ->
/// 3bp stem, 4nt AAAA loop, dot-bracket "(((....)))"), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindStemLoopsTests
{
    [Test]
    public void FindStemLoops_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindStemLoops("GGGAAAACCC", 3, 4, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindStemLoops("", 3, 4, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindStemLoops(null!, 3, 4, 4));
    }

    [Test]
    public void FindStemLoops_Binding_InvokesSuccessfully()
    {
        // GGGAAAACCC -> one hairpin: 3bp stem GGG/CCC, 4nt loop AAAA.
        var loops = AnalysisTools.FindStemLoops("GGGAAAACCC", 3, 4, 4).Items;
        Assert.Multiple(() =>
        {
            Assert.That(loops, Has.Length.EqualTo(1));
            var sl = loops[0];
            Assert.That(sl.Start, Is.EqualTo(0));
            Assert.That(sl.End, Is.EqualTo(9));
            Assert.That(sl.Stem.Length, Is.EqualTo(3));
            Assert.That(sl.Stem.Start5Prime, Is.EqualTo(0));
            Assert.That(sl.Stem.End3Prime, Is.EqualTo(9));
            Assert.That(sl.Loop.Sequence, Is.EqualTo("AAAA"));
            Assert.That(sl.DotBracketNotation, Is.EqualTo("(((....)))"));
        });
    }

    // RNA-STEMLOOP-001 (B12): opt-in stems with bulges. GGGACAGGAAACUGCCC: ViennaRNA 2.7.2 RNA.fold MFE
    // (NNDB-precision parameters) = (((.(((....)))))) −7.23 kcal/mol; perfect-stem mode finds only (((....))).
    [Test]
    public void FindStemLoops_MaxInternalLoopSize_DelegatesToUnbranchedStems()
    {
        var perfect = AnalysisTools.FindStemLoops("GGGACAGGAAACUGCCC", 3, 4, 4).Items;
        var bulged = AnalysisTools.FindStemLoops("GGGACAGGAAACUGCCC", 3, 4, 4, true, 1).Items;
        Assert.Multiple(() =>
        {
            Assert.That(perfect.Single().DotBracketNotation, Is.EqualTo("(((....)))"));
            Assert.That(bulged.Single().DotBracketNotation, Is.EqualTo("(((.(((....))))))"));
            Assert.That(bulged.Single().TotalFreeEnergy, Is.EqualTo(-7.23).Within(0.005));
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindStemLoops("GGGAAAACCC", 3, 4, 4, true, -1));
    }
}
