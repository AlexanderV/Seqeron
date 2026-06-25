using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ChromosomeAnalyzerTests
{
    // NOTE: Karyotype Analysis Tests moved to ChromosomeAnalyzer_Karyotype_Tests.cs
    // as part of CHROM-KARYO-001 test unit consolidation

    // NOTE: Telomere Analysis Tests moved to ChromosomeAnalyzer_Telomere_Tests.cs
    // as part of CHROM-TELO-001 test unit consolidation

    // NOTE: Centromere Analysis Tests moved to ChromosomeAnalyzer_Centromere_Tests.cs
    // as part of CHROM-CENT-001 test unit consolidation

    // NOTE: Aneuploidy Detection Tests moved to ChromosomeAnalyzer_Aneuploidy_Tests.cs
    // as part of CHROM-ANEU-001 test unit consolidation

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

    // NOTE: Synteny Analysis Tests moved to ChromosomeAnalyzer_Synteny_Tests.cs
    // as part of CHROM-SYNT-001 test unit consolidation

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

    // Levan, Fredga & Sandberg (1964): r = long/short ≥ 1 → metacentric ≤1.7,
    // submetacentric ≤3.0, subtelocentric <7.0, acrocentric ≥7.0. The method normalises
    // either p/q or q/p input to r ≥ 1. Expected values hand-derived from those thresholds.
    [Test]
    [TestCase(1.0, "Metacentric")]      // r = 1.0
    [TestCase(0.7, "Metacentric")]      // q/p = 1.43 < 1.7
    [TestCase(0.5, "Submetacentric")]   // q/p = 2.0
    [TestCase(0.3, "Subtelocentric")]   // q/p = 3.33
    [TestCase(0.1, "Acrocentric")]      // q/p = 10.0
    [TestCase(1.5, "Metacentric")]      // r = 1.5 < 1.7
    [TestCase(3.0, "Submetacentric")]   // r = 3.0 boundary
    [TestCase(10.0, "Acrocentric")]     // r = 10.0
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

    // AlphaSatelliteConsensus tests consolidated into ChromosomeAnalyzer_Centromere_Tests.cs
    // as part of CHROM-CENT-001 test unit

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
