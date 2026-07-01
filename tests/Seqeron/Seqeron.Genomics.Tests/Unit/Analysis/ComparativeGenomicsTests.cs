using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomicsTests
{
    // NOTE: FindSyntenicBlocks tests moved to the canonical unit file
    // ComparativeGenomics_FindSyntenicBlocks_Tests.cs (COMPGEN-SYNTENY-001).

    // FindOrthologs tests moved to the canonical file for COMPGEN-ORTHO-001:
    // ComparativeGenomics_FindOrthologs_Tests.cs (RBH + FindParalogs, evidence-based).

    // FindReciprocalBestHits tests moved to the canonical file for COMPGEN-RBH-001:
    // ComparativeGenomics_FindReciprocalBestHits_Tests.cs (evidence-based RBH).

    // DetectRearrangements tests moved to ComparativeGenomics_DetectRearrangements_Tests.cs
    // (COMPGEN-REARR-001), which tests the corrected breakpoint-based behavior with exact
    // evidence-derived counts. The prior permissive heuristic tests were removed there.

    // CompareGenomes tests moved to the canonical unit file
    // ComparativeGenomics_CompareGenomes_Tests.cs (COMPGEN-COMPARE-001): evidence-based
    // core/dispensable partition and syntenic-gene fraction with exact expected values.

    #region CalculateReversalDistance Tests

    [Test]
    public void CalculateReversalDistance_IdenticalPermutations_ReturnsZero()
    {
        var perm1 = new List<int> { 1, 2, 3, 4, 5 };
        var perm2 = new List<int> { 1, 2, 3, 4, 5 };

        int distance = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(distance, Is.EqualTo(0));
    }

    [Test]
    public void CalculateReversalDistance_ReversedPermutation_ReturnsPositive()
    {
        var perm1 = new List<int> { 1, 2, 3, 4, 5 };
        var perm2 = new List<int> { 5, 4, 3, 2, 1 };

        int distance = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(distance, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateReversalDistance_SingleElement_ReturnsZero()
    {
        var perm1 = new List<int> { 1 };
        var perm2 = new List<int> { 1 };

        int distance = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(distance, Is.EqualTo(0));
    }

    [Test]
    public void CalculateReversalDistance_DifferentLengths_ThrowsException()
    {
        var perm1 = new List<int> { 1, 2, 3 };
        var perm2 = new List<int> { 1, 2 };

        Assert.Throws<ArgumentException>(() =>
            ComparativeGenomics.CalculateReversalDistance(perm1, perm2));
    }

    [Test]
    public void CalculateReversalDistance_PartialReversal_ReturnsExpectedRange()
    {
        var perm1 = new List<int> { 1, 2, 3, 4, 5 };
        var perm2 = new List<int> { 1, 4, 3, 2, 5 }; // Middle reversed

        int distance = ComparativeGenomics.CalculateReversalDistance(perm1, perm2);

        Assert.That(distance, Is.GreaterThanOrEqualTo(1));
    }

    #endregion

    // FindConservedClusters tests moved to the canonical evidence-based fixture
    // ComparativeGenomics_FindConservedClusters_Tests.cs (COMPGEN-CLUSTER-001).

    // CalculateANI tests are consolidated into the canonical fixture
    // ComparativeGenomics_CalculateANI_Tests.cs (COMPGEN-ANI-001). The previous permissive
    // tests here (GreaterThan/LessThan ranges) were superseded by exact evidence-based cases.

    #region GenerateDotPlot Tests

    [Test]
    public void GenerateDotPlot_IdenticalSequences_ReturnsDiagonal()
    {
        string sequence = "ATGCATGCATGCATGCATGC";

        var points = ComparativeGenomics.GenerateDotPlot(sequence, sequence, wordSize: 5).ToList();

        Assert.That(points, Is.Not.Empty);
        // Check for diagonal points
        Assert.That(points.Any(p => p.x == p.y));
    }

    [Test]
    public void GenerateDotPlot_NoMatch_ReturnsEmpty()
    {
        string seq1 = "AAAAAAAAAAAAAAAAAAAAAA";
        string seq2 = "TTTTTTTTTTTTTTTTTTTTTT";

        var points = ComparativeGenomics.GenerateDotPlot(seq1, seq2, wordSize: 5).ToList();

        Assert.That(points, Is.Empty);
    }

    [Test]
    public void GenerateDotPlot_RepeatedSequence_ReturnsMultipleHits()
    {
        string seq1 = "ATGCATGCATGCATGC";
        string seq2 = "ATGCATGCATGCATGC";

        var points = ComparativeGenomics.GenerateDotPlot(seq1, seq2, wordSize: 4, stepSize: 1).ToList();

        Assert.That(points.Count, Is.GreaterThan(5));
    }

    [Test]
    public void GenerateDotPlot_EmptySequence_ReturnsEmpty()
    {
        string seq1 = "";
        string seq2 = "ATGCATGC";

        var points = ComparativeGenomics.GenerateDotPlot(seq1, seq2).ToList();

        Assert.That(points, Is.Empty);
    }

    [Test]
    public void GenerateDotPlot_InvertedRepeat_DetectsAntiDiagonal()
    {
        string seq1 = "ATGCATGCATGC";
        string seq2 = "GCATGCATGCAT"; // Shifted version

        var points = ComparativeGenomics.GenerateDotPlot(seq1, seq2, wordSize: 4).ToList();

        Assert.That(points, Is.Not.Empty);
    }

    #endregion

    #region Edge Cases

    [Test]
    public void SyntenicBlock_RecordProperties_Work()
    {
        var block = new ComparativeGenomics.SyntenicBlock(
            Genome1Id: "genome1",
            Start1: 0,
            End1: 1000,
            Genome2Id: "genome2",
            Start2: 500,
            End2: 1500,
            IsInverted: false,
            GeneCount: 10,
            Identity: 0.95);

        Assert.That(block.Genome1Id, Is.EqualTo("genome1"));
        Assert.That(block.GeneCount, Is.EqualTo(10));
        Assert.That(block.Identity, Is.EqualTo(0.95));
    }

    [Test]
    public void OrthologPair_RecordProperties_Work()
    {
        var pair = new ComparativeGenomics.OrthologPair(
            Gene1Id: "gene1",
            Gene2Id: "geneA",
            Identity: 0.85,
            Coverage: 0.90,
            AlignmentLength: 300);

        Assert.That(pair.Gene1Id, Is.EqualTo("gene1"));
        Assert.That(pair.Identity, Is.EqualTo(0.85));
    }

    [Test]
    public void RearrangementEvent_RecordProperties_Work()
    {
        var rearrangement = new ComparativeGenomics.RearrangementEvent(
            Type: ComparativeGenomics.RearrangementType.Inversion,
            GenomeId: "genome1",
            Position: 1000,
            Length: 500,
            TargetPosition: "2000");

        Assert.That(rearrangement.Type, Is.EqualTo(ComparativeGenomics.RearrangementType.Inversion));
        Assert.That(rearrangement.Length, Is.EqualTo(500));
    }

    [Test]
    public void Gene_RecordProperties_Work()
    {
        var gene = new ComparativeGenomics.Gene(
            Id: "gene1",
            GenomeId: "genome1",
            Start: 100,
            End: 500,
            Strand: '+',
            Sequence: "ATGC");

        Assert.That(gene.Id, Is.EqualTo("gene1"));
        Assert.That(gene.Strand, Is.EqualTo('+'));
    }

    [Test]
    public void RearrangementType_AllValuesExist()
    {
        var types = Enum.GetValues<ComparativeGenomics.RearrangementType>();

        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Inversion));
        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Translocation));
        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Deletion));
        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Insertion));
        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Duplication));
        Assert.That(types, Contains.Item(ComparativeGenomics.RearrangementType.Transposition));
    }

    #endregion
}
