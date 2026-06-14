using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ComparativeGenomicsTests
{
    #region Test Data

    private static List<ComparativeGenomics.Gene> CreateTestGenome1()
    {
        return new List<ComparativeGenomics.Gene>
        {
            new("gene1", "genome1", 0, 100, '+', "ATGCATGCATGC"),
            new("gene2", "genome1", 150, 250, '+', "GCTAGCTAGCTA"),
            new("gene3", "genome1", 300, 400, '+', "TATATATATAT"),
            new("gene4", "genome1", 450, 550, '+', "CGCGCGCGCG"),
            new("gene5", "genome1", 600, 700, '+', "AAAATTTTCCCC"),
        };
    }

    private static List<ComparativeGenomics.Gene> CreateTestGenome2()
    {
        return new List<ComparativeGenomics.Gene>
        {
            new("geneA", "genome2", 0, 100, '+', "ATGCATGCATGC"), // Ortholog of gene1
            new("geneB", "genome2", 150, 250, '+', "GCTAGCTAGCTA"), // Ortholog of gene2
            new("geneC", "genome2", 300, 400, '+', "TATATATATAT"), // Ortholog of gene3
            new("geneD", "genome2", 450, 550, '+', "CGCGCGCGCG"), // Ortholog of gene4
            new("geneE", "genome2", 600, 700, '+', "GGGGTTTTAAAA"), // Different from gene5
        };
    }

    private static Dictionary<string, string> CreateOrthologMap()
    {
        return new Dictionary<string, string>
        {
            { "gene1", "geneA" },
            { "gene2", "geneB" },
            { "gene3", "geneC" },
            { "gene4", "geneD" },
        };
    }

    #endregion

    // NOTE: FindSyntenicBlocks tests moved to the canonical unit file
    // ComparativeGenomics_FindSyntenicBlocks_Tests.cs (COMPGEN-SYNTENY-001).

    // FindOrthologs tests moved to the canonical file for COMPGEN-ORTHO-001:
    // ComparativeGenomics_FindOrthologs_Tests.cs (RBH + FindParalogs, evidence-based).

    // FindReciprocalBestHits tests moved to the canonical file for COMPGEN-RBH-001:
    // ComparativeGenomics_FindReciprocalBestHits_Tests.cs (evidence-based RBH).

    // DetectRearrangements tests moved to ComparativeGenomics_DetectRearrangements_Tests.cs
    // (COMPGEN-REARR-001), which tests the corrected breakpoint-based behavior with exact
    // evidence-derived counts. The prior permissive heuristic tests were removed there.

    #region CompareGenomes Tests

    [Test]
    public void CompareGenomes_SimilarGenomes_ReturnsComprehensiveResult()
    {
        var genome1 = CreateTestGenome1();
        var genome2 = CreateTestGenome2();

        var result = ComparativeGenomics.CompareGenomes(genome1, genome2, minOrthologIdentity: 0.3);

        Assert.That(result.ConservedGenes, Is.GreaterThan(0));
        Assert.That(result.Orthologs.Count, Is.GreaterThan(0));
        Assert.That(result.OverallSynteny, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void CompareGenomes_EmptyGenomes_HandlesGracefully()
    {
        var genome1 = new List<ComparativeGenomics.Gene>();
        var genome2 = new List<ComparativeGenomics.Gene>();

        var result = ComparativeGenomics.CompareGenomes(genome1, genome2);

        Assert.That(result.ConservedGenes, Is.EqualTo(0));
        Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(0));
        Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(0));
    }

    [Test]
    public void CompareGenomes_CompletelyDifferent_ReturnsNoConservation()
    {
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            new("gene1", "genome1", 0, 100, '+', "AAAAAAAAAAAAAAAAAAAAAAAAA"),
            new("gene2", "genome1", 150, 250, '+', "AAAAAAAAAAAAAAAAAAAAAAAA"),
        };
        var genome2 = new List<ComparativeGenomics.Gene>
        {
            new("geneA", "genome2", 0, 100, '+', "TTTTTTTTTTTTTTTTTTTTTTTTT"),
            new("geneB", "genome2", 150, 250, '+', "TTTTTTTTTTTTTTTTTTTTTTTT"),
        };

        var result = ComparativeGenomics.CompareGenomes(genome1, genome2);

        Assert.That(result.ConservedGenes, Is.EqualTo(0));
        Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(2));
        Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(2));
    }

    #endregion

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

    #region FindConservedClusters Tests

    [Test]
    public void FindConservedClusters_ConservedCluster_ReturnsCluster()
    {
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            new("gene1", "g1", 0, 100, '+'),
            new("gene2", "g1", 150, 250, '+'),
            new("gene3", "g1", 300, 400, '+'),
        };
        var genome2 = new List<ComparativeGenomics.Gene>
        {
            new("geneA", "g2", 0, 100, '+'),
            new("geneB", "g2", 150, 250, '+'),
            new("geneC", "g2", 300, 400, '+'),
        };

        var orthologGroups = new Dictionary<string, string>
        {
            { "gene1", "group1" }, { "geneA", "group1" },
            { "gene2", "group2" }, { "geneB", "group2" },
            { "gene3", "group3" }, { "geneC", "group3" },
        };

        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { genome1, genome2 };
        var clusters = ComparativeGenomics.FindConservedClusters(genomes, orthologGroups, minClusterSize: 3).ToList();

        Assert.That(clusters, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(clusters[0], Has.Count.GreaterThanOrEqualTo(3));
    }

    [Test]
    public void FindConservedClusters_SingleGenome_ReturnsEmpty()
    {
        var genome1 = CreateTestGenome1();
        var orthologGroups = new Dictionary<string, string> { { "gene1", "group1" } };
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>> { genome1 };

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, orthologGroups).ToList();

        Assert.That(clusters, Is.Empty);
    }

    #endregion

    #region CalculateANI Tests

    [Test]
    public void CalculateANI_IdenticalSequences_ReturnsHigh()
    {
        string genome1 = string.Concat(Enumerable.Repeat("ATGCATGCATGC", 200)); // 2400 bp
        string genome2 = genome1;

        double ani = ComparativeGenomics.CalculateANI(genome1, genome2, fragmentSize: 100);

        Assert.That(ani, Is.GreaterThan(0.9));
    }

    [Test]
    public void CalculateANI_CompletelyDifferent_ReturnsLow()
    {
        string genome1 = string.Concat(Enumerable.Repeat("AAAAAAAAAA", 200));
        string genome2 = string.Concat(Enumerable.Repeat("TTTTTTTTTT", 200));

        double ani = ComparativeGenomics.CalculateANI(genome1, genome2, fragmentSize: 100);

        Assert.That(ani, Is.LessThan(0.3));
    }

    [Test]
    public void CalculateANI_EmptySequence_ReturnsZero()
    {
        string genome1 = "";
        string genome2 = "ATGCATGCATGC";

        double ani = ComparativeGenomics.CalculateANI(genome1, genome2);

        Assert.That(ani, Is.EqualTo(0));
    }

    [Test]
    public void CalculateANI_NullSequence_ReturnsZero()
    {
        string genome1 = null!;
        string genome2 = "ATGCATGCATGC";

        double ani = ComparativeGenomics.CalculateANI(genome1, genome2);

        Assert.That(ani, Is.EqualTo(0));
    }

    [Test]
    public void CalculateANI_PartialSimilarity_ReturnsIntermediate()
    {
        // Create sequences where only a portion matches
        // Use completely different patterns to ensure distinct fragments
        string pattern1 = "ATGCATGCATGCATGCATGCATGCATGCATGC"; // 32bp repeating pattern
        string pattern2 = "TTTTAAAATTTTAAAATTTTAAAATTTTAAAA"; // Completely different pattern

        // Build genomes with distinct halves (2000bp each for reliable fragment comparison)
        string sharedPart = string.Concat(Enumerable.Repeat(pattern1, 63)); // ~2000bp
        string diffPart1 = string.Concat(Enumerable.Repeat(pattern2, 63));
        string diffPart2 = string.Concat(Enumerable.Repeat("CCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGG", 63));

        string genome1 = sharedPart + diffPart1;
        string genome2 = sharedPart + diffPart2;

        double ani = ComparativeGenomics.CalculateANI(genome1, genome2, fragmentSize: 200, minFragmentIdentity: 0.5);

        // Should find matches in the shared part but not in the different parts
        Assert.That(ani, Is.GreaterThan(0.0)); // Some fragments match
        Assert.That(ani, Is.LessThanOrEqualTo(1.0));
    }

    #endregion

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
