using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ChromosomeAnalyzerTests
{
    #region Karyotype Analysis Tests

    [Test]
    public void AnalyzeKaryotype_WithNormalDiploidSet_ReturnsCorrectCounts()
    {
        // For diploid, each autosome should appear twice (with copy suffix)
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false),
            ("chr1_2", 248956422, false),
            ("chr2_1", 242193529, false),
            ("chr2_2", 242193529, false),
            ("chrX", 156040895, true),
            ("chrY", 57227415, true)
        };

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        Assert.That(result.TotalChromosomes, Is.EqualTo(6));
        Assert.That(result.AutosomeCount, Is.EqualTo(4));
        Assert.That(result.SexChromosomes.Count, Is.EqualTo(2));
        Assert.That(result.HasAneuploidy, Is.False);
    }

    [Test]
    public void AnalyzeKaryotype_WithTrisomy_DetectsAneuploidy()
    {
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr21_1", 46709983, false),
            ("chr21_2", 46709983, false),
            ("chr21_3", 46709983, false) // Trisomy 21
        };

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        Assert.That(result.HasAneuploidy, Is.True);
        Assert.That(result.Abnormalities, Has.Some.Contains("Trisomy"));
    }

    [Test]
    public void AnalyzeKaryotype_WithMonosomy_DetectsAneuploidy()
    {
        var chromosomes = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            ("chr1_1", 248956422, false) // Monosomy - only one copy
        };

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes);

        Assert.That(result.HasAneuploidy, Is.True);
        Assert.That(result.Abnormalities, Has.Some.Contains("Monosomy"));
    }

    [Test]
    public void AnalyzeKaryotype_EmptyInput_ReturnsEmptyKaryotype()
    {
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(
            Enumerable.Empty<(string, long, bool)>());

        Assert.That(result.TotalChromosomes, Is.EqualTo(0));
        Assert.That(result.HasAneuploidy, Is.False);
    }

    [Test]
    public void DetectPloidy_WithDiploidDepth_ReturnsPloidy2()
    {
        var depths = Enumerable.Repeat(1.0, 100);

        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        Assert.That(ploidy, Is.EqualTo(2));
        Assert.That(confidence, Is.GreaterThan(0.9));
    }

    [Test]
    public void DetectPloidy_WithTetraploidDepth_ReturnsPloidy4()
    {
        var depths = Enumerable.Repeat(2.0, 100);

        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(depths);

        Assert.That(ploidy, Is.EqualTo(4));
        Assert.That(confidence, Is.GreaterThan(0.9));
    }

    [Test]
    public void DetectPloidy_EmptyInput_ReturnsDefault()
    {
        var (ploidy, confidence) = ChromosomeAnalyzer.DetectPloidy(
            Enumerable.Empty<double>());

        Assert.That(ploidy, Is.EqualTo(2));
        Assert.That(confidence, Is.EqualTo(0));
    }

    #endregion

    // NOTE: Telomere Analysis Tests moved to ChromosomeAnalyzer_Telomere_Tests.cs
    // as part of CHROM-TELO-001 test unit consolidation

    // NOTE: Centromere Analysis Tests moved to ChromosomeAnalyzer_Centromere_Tests.cs
    // as part of CHROM-CENT-001 test unit consolidation

    #region Cytogenetic Band Tests

    [Test]
    public void PredictGBands_GeneratesBands()
    {
        // Create sequence with varying GC content
        string atRichRegion = new string('A', 5000000);
        string gcRichRegion = new string('G', 5000000);
        string sequence = atRichRegion + gcRichRegion;

        var bands = ChromosomeAnalyzer.PredictGBands(
            "chr1", sequence, bandSize: 5000000).ToList();

        Assert.That(bands.Count, Is.GreaterThan(0));
        Assert.That(bands.All(b => b.Chromosome == "chr1"), Is.True);

        // First band should be dark (AT-rich)
        Assert.That(bands[0].Stain, Does.Contain("gpos"));

        // Second band should be light (GC-rich)
        Assert.That(bands[1].Stain, Is.EqualTo("gneg"));
    }

    [Test]
    public void PredictGBands_EmptySequence_YieldsNoBands()
    {
        var bands = ChromosomeAnalyzer.PredictGBands("chr1", "").ToList();

        Assert.That(bands, Is.Empty);
    }

    [Test]
    public void PredictGBands_IncludesGcContentAndGeneDensity()
    {
        string sequence = new string('G', 10000000);

        var bands = ChromosomeAnalyzer.PredictGBands(
            "chr1", sequence, bandSize: 5000000).ToList();

        Assert.That(bands.All(b => b.GcContent > 0.9), Is.True);
        Assert.That(bands.All(b => b.GeneDensity > 0), Is.True);
    }

    [Test]
    public void FindHeterochromatinRegions_EmptySequence_YieldsNoRegions()
    {
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions("")
            .ToList();

        Assert.That(regions, Is.Empty);
    }

    #endregion

    #region Synteny Analysis Tests

    [Test]
    public void FindSyntenyBlocks_WithCollinearGenes_FindsBlocks()
    {
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
            ("chr1", 7000, 8000, "gene4", "chrA", 7000, 8000, "geneD"),
        };

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        Assert.That(blocks.Count, Is.GreaterThanOrEqualTo(1));
        Assert.That(blocks[0].GeneCount, Is.GreaterThanOrEqualTo(3));
        Assert.That(blocks[0].Strand, Is.EqualTo('+'));
    }

    [Test]
    public void FindSyntenyBlocks_WithInvertedGenes_DetectsInversion()
    {
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 8000, 9000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 6000, 7000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 4000, 5000, "geneC"),
            ("chr1", 7000, 8000, "gene4", "chrA", 2000, 3000, "geneD"),
        };

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        Assert.That(blocks.Any(b => b.Strand == '-'), Is.True);
    }

    [Test]
    public void FindSyntenyBlocks_TooFewGenes_ReturnsEmpty()
    {
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
        };

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3).ToList();

        Assert.That(blocks, Is.Empty);
    }

    [Test]
    public void DetectRearrangements_WithTranslocation_DetectsTranslocation()
    {
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrB", 1000, 40000, '+', 8, 0.93)
        };

        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        Assert.That(rearrangements.Any(r => r.Type == "Translocation"), Is.True);
    }

    [Test]
    public void DetectRearrangements_WithInversion_DetectsInversion()
    {
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrA", 60000, 100000, '-', 8, 0.93)
        };

        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        Assert.That(rearrangements.Any(r => r.Type == "Inversion"), Is.True);
    }

    #endregion

    #region Aneuploidy Detection Tests

    [Test]
    public void DetectAneuploidy_WithNormalDepth_ReturnsCopyNumber2()
    {
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chr1", i * 1000000, 30.0));

        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: 30.0, binSize: 1000000).ToList();

        Assert.That(results.All(r => r.CopyNumber == 2), Is.True);
    }

    [Test]
    public void DetectAneuploidy_WithTrisomy_ReturnsCopyNumber3()
    {
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chr21", i * 1000000, 45.0)); // 1.5x normal depth

        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: 30.0, binSize: 1000000).ToList();

        Assert.That(results.All(r => r.CopyNumber == 3), Is.True);
    }

    [Test]
    public void DetectAneuploidy_WithMonosomy_ReturnsCopyNumber1()
    {
        var depthData = Enumerable.Range(0, 100)
            .Select(i => ("chrX", i * 1000000, 15.0)); // 0.5x normal depth

        var results = ChromosomeAnalyzer.DetectAneuploidy(
            depthData, medianDepth: 30.0, binSize: 1000000).ToList();

        Assert.That(results.All(r => r.CopyNumber == 1), Is.True);
    }

    [Test]
    public void DetectAneuploidy_EmptyInput_ReturnsEmpty()
    {
        var results = ChromosomeAnalyzer.DetectAneuploidy(
            Enumerable.Empty<(string, int, double)>(), 30.0).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_WithTrisomy_IdentifiesTrisomy()
    {
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr21", i * 1000000, (i + 1) * 1000000 - 1, 3, 0.58, 0.9));

        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        Assert.That(results.Count, Is.EqualTo(1));
        Assert.That(results[0].Chromosome, Is.EqualTo("chr21"));
        Assert.That(results[0].CopyNumber, Is.EqualTo(3));
        Assert.That(results[0].Type, Is.EqualTo("Trisomy"));
    }

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_WithNormalChromosome_ReturnsEmpty()
    {
        var states = Enumerable.Range(0, 100)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState(
                "chr1", i * 1000000, (i + 1) * 1000000 - 1, 2, 0, 0.95));

        var results = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        Assert.That(results, Is.Empty);
    }

    #endregion

    #region Utility Method Tests

    [Test]
    public void CalculateArmRatio_Metacentric_ReturnsNearOne()
    {
        double ratio = ChromosomeAnalyzer.CalculateArmRatio(
            centromerePosition: 50000000, chromosomeLength: 100000000);

        Assert.That(ratio, Is.InRange(0.9, 1.1));
    }

    [Test]
    public void CalculateArmRatio_Acrocentric_ReturnsLow()
    {
        double ratio = ChromosomeAnalyzer.CalculateArmRatio(
            centromerePosition: 10000000, chromosomeLength: 100000000);

        Assert.That(ratio, Is.LessThan(0.15));
    }

    [Test]
    public void CalculateArmRatio_InvalidInput_ReturnsZero()
    {
        Assert.That(ChromosomeAnalyzer.CalculateArmRatio(0, 100), Is.EqualTo(0));
        Assert.That(ChromosomeAnalyzer.CalculateArmRatio(50, 0), Is.EqualTo(0));
    }

    [Test]
    [TestCase(1.0, "Metacentric")]
    [TestCase(0.7, "Submetacentric")]
    [TestCase(0.3, "Acrocentric")]
    [TestCase(0.1, "Telocentric")]
    [TestCase(1.5, "Submetacentric")]
    [TestCase(3.0, "Acrocentric")]
    [TestCase(10.0, "Telocentric")]
    public void ClassifyChromosomeByArmRatio_ClassifiesCorrectly(
        double armRatio, string expectedType)
    {
        string result = ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(armRatio);

        Assert.That(result, Is.EqualTo(expectedType));
    }

    [Test]
    public void EstimateCellDivisionsFromTelomereLength_CalculatesCorrectly()
    {
        double divisions = ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength(
            currentLength: 10000, birthLength: 15000, lossPerDivision: 50);

        Assert.That(divisions, Is.EqualTo(100));
    }

    [Test]
    public void EstimateCellDivisionsFromTelomereLength_LongerThanBirth_ReturnsZero()
    {
        double divisions = ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength(
            currentLength: 20000, birthLength: 15000, lossPerDivision: 50);

        Assert.That(divisions, Is.EqualTo(0));
    }

    [Test]
    public void EstimateCellDivisionsFromTelomereLength_ZeroLoss_ReturnsZero()
    {
        double divisions = ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength(
            currentLength: 10000, birthLength: 15000, lossPerDivision: 0);

        Assert.That(divisions, Is.EqualTo(0));
    }

    #endregion

    #region Constants Tests

    [Test]
    public void HumanTelomereRepeat_IsCorrect()
    {
        Assert.That(ChromosomeAnalyzer.HumanTelomereRepeat, Is.EqualTo("TTAGGG"));
    }

    [Test]
    public void AlphaSatelliteConsensus_IsNotEmpty()
    {
        Assert.That(ChromosomeAnalyzer.AlphaSatelliteConsensus, Is.Not.Empty);
        Assert.That(ChromosomeAnalyzer.AlphaSatelliteConsensus.Length, Is.GreaterThan(50));
    }

    #endregion

    #region Helper Methods

    private static string CreateRandomSequence(System.Random random, int length)
    {
        var bases = new char[] { 'A', 'C', 'G', 'T' };
        var sequence = new char[length];

        for (int i = 0; i < length; i++)
        {
            sequence[i] = bases[random.Next(4)];
        }

        return new string(sequence);
    }

    #endregion
}
