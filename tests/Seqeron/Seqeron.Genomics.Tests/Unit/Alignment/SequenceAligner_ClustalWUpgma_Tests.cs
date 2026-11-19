using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14, B09 F92): ClustalW 2.1 <c>-CLUSTERING=UPGMA</c> (<c>UPGMAAlgorithm</c>, <c>Node</c>, <c>RootedGuideTree</c>, <c>TreeInterface::getWeightsAndStepsFromDistMatUPGMA</c>).
/// The installed Debian <c>clustalw</c> segmentation-faults on every UPGMA run, so there is no row oracle yet (hand-off package <c>handoff/A2-14-clustalw-upgma</c>, built from the unmodified source
/// by the user); these tests lock the algorithm to hand-derived values of the source's definitions (average-linkage joins, join heights = distance / 2, the weights of
/// <c>RootedGuideTree::doWeightCalc</c> in which a node carries the height of its PARENT's join) and to the program's two-sequence behaviour read from the source.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWUpgma_Tests
{
    [Test]
    public void ThreeTaxa_JoinsByAverageLinkage_StepsAndHeights()
    {
        // d(1,2) = 0.2, d(1,3) = d(2,3) = 0.6: (1,2) joins at height 0.1, then {1,2} + 3 at the mean distance 0.6 -> height 0.3
        var d = new double[3, 3];
        d[1, 0] = d[0, 1] = 0.2;
        d[2, 0] = d[0, 2] = 0.6;
        d[2, 1] = d[1, 2] = 0.6;
        var t = SequenceAligner.ClustalWUpgmaTree(d, 3);
        Assert.Multiple(() =>
        {
            Assert.That(t.Steps, Has.Count.EqualTo(2));
            Assert.That(t.Steps[0], Is.EqualTo(new[] { 0, 1, 2, 0 }), "step 1: sequences 1 and 2 (group 1 = the lower index, group 2 = the joined node)");
            Assert.That(t.Steps[1], Is.EqualTo(new[] { 0, 1, 1, 2 }), "step 2: the cluster {1, 2} against sequence 3");
            Assert.That(t.Root.Height, Is.EqualTo(0.3).Within(1e-12));
            Assert.That(t.Root.Size, Is.EqualTo(3));
            Assert.That(t.Root.Left!.Height, Is.EqualTo(0.3).Within(1e-12), "a child carries the height of the join that absorbed it (Node::merge)");
        });
    }

    [Test]
    public void ThreeTaxa_SequenceWeights_FollowDoWeightCalc()
    {
        var d = new double[3, 3];
        d[1, 0] = d[0, 1] = 0.2;
        d[2, 0] = d[0, 2] = 0.6;
        d[2, 1] = d[1, 2] = 0.6;
        var t = SequenceAligner.ClustalWUpgmaTree(d, 3);
        var w = SequenceAligner.ClustalWUpgmaSeqWeights(t.Root, 3);
        // leaf 3: root 0.3 / 3 = 0.1; leaves 1, 2: 0.1 + 0.3 / 2 = 0.25 (the child {1, 2} carries the root's height) -> (int)(w x 100) = 25, 25, 10 (single precision) -> x 1000 / 60 = 416, 416, 166
        Assert.That(w, Is.EqualTo(new[] { 416, 416, 166 }));
    }

    [Test]
    public void NoTree_GivesEqualWeights_AsTheProgramWithTwoSequences()
    {
        Assert.That(SequenceAligner.ClustalWUpgmaSeqWeights(null, 2), Is.EqualTo(new[] { 500, 500 }), "all weights 0 -> 1 -> 1000 / 2");
    }

    [Test]
    public void MultipleAlignClustalW_Upgma_KeepsEverySequenceIntact_ForDnaAndProtein()
    {
        string[] dna = ["ACGTACGTAAGCT", "ACGTACGTCAGCT", "ACGTTCGTCAGGT", "TCGTTCGTCAGGT", "ACGTACGTCAGCTAA"];
        var res = SequenceAligner.MultipleAlignClustalW(dna, new ClustalWOptions { Clustering = ClustalWClustering.Upgma });
        Assert.That(res.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(dna));
        Assert.That(res.AlignedSequences.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));

        string[] prot = ["MKVLAAGIVGLLAQ", "MKVLAGIVGLLAQP", "MRVLAAGIIGLLQP", "KVLAAGIVALLAQ"];
        var pres = SequenceAligner.MultipleAlignClustalW(prot, new ClustalWProteinOptions { Clustering = ClustalWClustering.Upgma });
        Assert.That(pres.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(prot));
        Assert.That(pres.AlignedSequences.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
    }

    [Test]
    public void Default_IsNeighborJoining_AndDiffersFromUpgmaOnlyThroughTheTree()
    {
        string[] dna = ["ACGTACGTAAGCT", "ACGTACGTCAGCT", "ACGTTCGTCAGGT", "TCGTTCGTCAGGT", "ACGTACGTCAGCTAA"];
        Assert.That(new ClustalWOptions().Clustering, Is.EqualTo(ClustalWClustering.NeighborJoining));
        var nj = SequenceAligner.MultipleAlignClustalW(dna);
        Assert.That(nj.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(dna));
    }
}
