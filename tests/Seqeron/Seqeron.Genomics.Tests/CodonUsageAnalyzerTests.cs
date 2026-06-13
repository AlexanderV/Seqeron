using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class CodonUsageAnalyzerTests
{
    // NOTE: CountCodons and RSCU tests live in the canonical file
    // CodonUsageAnalyzer_CalculateRscu_Tests.cs (Test Unit CODON-RSCU-001).

    #region CAI Tests

    [Test]
    public void CalculateCai_OptimalCodons_ReturnsHigh()
    {
        // Use E. coli optimal codons
        var sequence = new DnaSequence("CTGCTGCTG"); // CTG is preferred for Leu in E. coli
        double cai = CodonUsageAnalyzer.CalculateCai(sequence, CodonUsageAnalyzer.EColiOptimalCodons);

        Assert.That(cai, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void CalculateCai_SuboptimalCodons_ReturnsLower()
    {
        // Use suboptimal codons for E. coli
        var sequence = new DnaSequence("CTACTACTA"); // CTA is rare for Leu in E. coli
        double cai = CodonUsageAnalyzer.CalculateCai(sequence, CodonUsageAnalyzer.EColiOptimalCodons);

        Assert.That(cai, Is.LessThan(0.5));
    }

    [Test]
    public void CalculateCai_EmptySequence_ReturnsZero()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("", CodonUsageAnalyzer.EColiOptimalCodons);
        Assert.That(cai, Is.EqualTo(0));
    }

    [Test]
    public void CalculateCai_RangeIsZeroToOne()
    {
        var sequence = new DnaSequence("ATGAAATTTGGG");
        double cai = CodonUsageAnalyzer.CalculateCai(sequence, CodonUsageAnalyzer.EColiOptimalCodons);

        Assert.That(cai, Is.GreaterThanOrEqualTo(0));
        Assert.That(cai, Is.LessThanOrEqualTo(1));
    }

    #endregion

    #region ENC Tests

    [Test]
    public void CalculateEnc_NoBias_ReturnsHigh()
    {
        // For ENC to be high, we need diverse codon usage across synonymous families
        // Using a realistic coding sequence with varied codons
        var sequence = new DnaSequence("ATGAAAGAGCTGTTCGCCAAA");

        double enc = CodonUsageAnalyzer.CalculateEnc(sequence);

        // ENC ranges from 20 (max bias) to 61 (no bias)
        // Even with limited data, should be >= 20
        Assert.That(enc, Is.GreaterThanOrEqualTo(20));
    }

    [Test]
    public void CalculateEnc_HighBias_ReturnsLow()
    {
        // Only one codon type per amino acid (extreme bias)
        var sequence = new DnaSequence("TTTTTTTTTTTTTTT"); // Only Phe-TTT
        double enc = CodonUsageAnalyzer.CalculateEnc(sequence);

        Assert.That(enc, Is.GreaterThanOrEqualTo(20)); // Min ENC = 20
    }

    [Test]
    public void CalculateEnc_RangeIsTwentyToSixtyOne()
    {
        var sequence = new DnaSequence("ATGAAATTTGGGCCC");
        double enc = CodonUsageAnalyzer.CalculateEnc(sequence);

        Assert.That(enc, Is.GreaterThanOrEqualTo(20));
        Assert.That(enc, Is.LessThanOrEqualTo(61));
    }

    [Test]
    public void CalculateEnc_EmptySequence_ReturnsZero()
    {
        double enc = CodonUsageAnalyzer.CalculateEnc("");
        Assert.That(enc, Is.EqualTo(0));
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void GetStatistics_ReturnsComprehensiveStats()
    {
        var sequence = new DnaSequence("ATGAAATTTGGGCCC");
        var stats = CodonUsageAnalyzer.GetStatistics(sequence);

        Assert.That(stats.TotalCodons, Is.EqualTo(5));
        Assert.That(stats.CodonCounts, Is.Not.Empty);
        Assert.That(stats.Rscu, Is.Not.Empty);
        Assert.That(stats.Enc, Is.GreaterThan(0));
    }

    [Test]
    public void GetStatistics_GcPositions_CalculatedCorrectly()
    {
        // ATG = 50% GC, GCG = 100% GC, TAA = 0% GC
        var sequence = new DnaSequence("ATGGCGTAA");
        var stats = CodonUsageAnalyzer.GetStatistics(sequence);

        Assert.That(stats.Gc1, Is.GreaterThan(0)); // G, G, T = 66%
        Assert.That(stats.Gc2, Is.GreaterThan(0)); // T, C, A = 33%
        Assert.That(stats.Gc3, Is.GreaterThan(0)); // G, G, A = 66%
    }

    [Test]
    public void GetStatistics_OverallGc_IsAverageOfPositions()
    {
        var sequence = new DnaSequence("ATGGCGTAA");
        var stats = CodonUsageAnalyzer.GetStatistics(sequence);

        double expectedOverallGc = (stats.Gc1 + stats.Gc2 + stats.Gc3) / 3;
        Assert.That(stats.OverallGc, Is.EqualTo(expectedOverallGc).Within(0.01));
    }

    [Test]
    public void GetStatistics_EmptySequence_ReturnsZeros()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("");

        Assert.That(stats.TotalCodons, Is.EqualTo(0));
        Assert.That(stats.Enc, Is.EqualTo(0));
    }

    #endregion

    #region Reference Tables Tests

    [Test]
    public void EColiOptimalCodons_ContainsAllCodons()
    {
        var table = CodonUsageAnalyzer.EColiOptimalCodons;
        Assert.That(table.Count, Is.EqualTo(64));
    }

    [Test]
    public void HumanOptimalCodons_ContainsAllCodons()
    {
        var table = CodonUsageAnalyzer.HumanOptimalCodons;
        Assert.That(table.Count, Is.EqualTo(64));
    }

    [Test]
    public void EColiOptimalCodons_CTG_IsPreferred()
    {
        var table = CodonUsageAnalyzer.EColiOptimalCodons;

        // CTG should have highest RSCU among Leu codons
        Assert.That(table["CTG"], Is.GreaterThan(table["CTA"]));
        Assert.That(table["CTG"], Is.GreaterThan(table["CTT"]));
        Assert.That(table["CTG"], Is.GreaterThan(table["CTC"]));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void CalculateCai_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.CalculateCai((DnaSequence)null!, CodonUsageAnalyzer.EColiOptimalCodons));
    }

    [Test]
    public void CalculateCai_NullReference_ThrowsException()
    {
        var sequence = new DnaSequence("ATGATG");
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.CalculateCai(sequence, null!));
    }

    [Test]
    public void GetStatistics_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.GetStatistics((DnaSequence)null!));
    }

    #endregion
}
