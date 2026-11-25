using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ComputeBranchLengthScoreTests
{
    // Kit targetscan_70_output.BL_PCT.txt (CDC2L6 let-7 rows): species list on Tree.bin_01 -> Branch length score 1.6379.
    private static readonly string[] KitSpecies = { "9031", "9606", "9615", "10090", "10116" };

    [Test]
    public void ComputeBranchLengthScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies, targetScanBin: 1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeBranchLengthScore(null!, targetScanBin: 1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies), "no tree");
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies, "(A:1,B:2);", 1), "both trees");
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies, targetScanBin: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies, targetScanBin: 11));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeBranchLengthScore(KitSpecies, "(A:1,B:"), "bad Newick");
    }

    [Test]
    public void ComputeBranchLengthScore_Binding_KitRow_Bin1_1_6379()
    {
        var r = AnnotationTools.ComputeBranchLengthScore(KitSpecies, targetScanBin: 1);
        var lib = MiRnaAnalyzer.ComputeBranchLengthScore(MiRnaAnalyzer.GetTargetScanBinTree(1), KitSpecies);
        Assert.That(r.BranchLengthScore, Is.EqualTo(1.6379));
        Assert.That(r.BranchLengthScore, Is.EqualTo(lib));
    }

    [Test]
    public void ComputeBranchLengthScore_Newick_SumsConnectingSubtree()
    {
        // ((A:1,B:2):3,C:4): A+B connect through edges 1+2 = 3 (the parent edge 3 is not on the path).
        Assert.That(AnnotationTools.ComputeBranchLengthScore(new[] { "A", "B" }, "((A:1,B:2):3,C:4);").BranchLengthScore, Is.EqualTo(3.0));
        Assert.That(AnnotationTools.ComputeBranchLengthScore(new[] { "A", "C" }, "((A:1,B:2):3,C:4);").BranchLengthScore, Is.EqualTo(8.0));
        Assert.That(AnnotationTools.ComputeBranchLengthScore(new[] { "A" }, "((A:1,B:2):3,C:4);").BranchLengthScore, Is.EqualTo(0.0), "one species");
        Assert.That(AnnotationTools.ComputeBranchLengthScore(Array.Empty<string>(), "((A:1,B:2):3,C:4);").BranchLengthScore, Is.EqualTo(0.0), "empty list");
        Assert.That(AnnotationTools.ComputeBranchLengthScore(new[] { "A", "Z" }, "((A:1,B:2):3,C:4);").BranchLengthScore, Is.EqualTo(0.0), "non-leaf ignored");
    }

    [Test]
    public void ComputeBranchLengthScore_WhitespaceNewick_PassedThroughToParser_NoWrapperNormalisation()
    {
        // The wrapper does not normalise whitespace: whatever PhylogeneticAnalyzer.ParseNewick does with it is what the tool does.
        const string multi = "((A:1,\nB:2):3,C:4);";
        double expected = MiRnaAnalyzer.ComputeBranchLengthScore(
            Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(multi), new[] { "A", "B" });
        Assert.That(AnnotationTools.ComputeBranchLengthScore(new[] { "A", "B" }, multi).BranchLengthScore, Is.EqualTo(expected));
    }
}
