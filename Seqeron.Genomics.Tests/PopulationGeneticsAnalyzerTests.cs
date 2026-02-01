using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class PopulationGeneticsAnalyzerTests
{
    // Note: Allele Frequency tests (POP-FREQ-001) have been moved to
    // PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs

    // Note: Diversity Statistics tests (POP-DIV-001) have been moved to
    // PopulationGeneticsAnalyzer_Diversity_Tests.cs

    // Note: Hardy-Weinberg tests (POP-HW-001) have been moved to
    // PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs

    // Note: F-Statistics tests (POP-FST-001) have been moved to
    // PopulationGeneticsAnalyzer_FStatistics_Tests.cs

    #region Linkage Disequilibrium Tests

    [Test]
    public void CalculateLD_PerfectLD_ReturnsHighValues()
    {
        // Perfect LD: genotypes always match
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 0), (1, 1), (1, 1), (2, 2), (2, 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        Assert.That(ld.RSquared, Is.GreaterThan(0.4)); // High but not perfect due to estimation
    }

    [Test]
    public void CalculateLD_NoLD_ReturnsLowValues()
    {
        // No LD: random pairing
        var genotypes = new List<(int, int)>
        {
            (0, 2), (2, 0), (1, 1), (0, 1), (2, 1), (1, 0)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        Assert.That(ld.RSquared, Is.LessThan(0.5));
    }

    [Test]
    public void CalculateLD_RecordsDistance()
    {
        var genotypes = new List<(int, int)> { (0, 0), (1, 1) };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 5000);

        Assert.That(ld.Distance, Is.EqualTo(5000));
        Assert.That(ld.Variant1, Is.EqualTo("V1"));
        Assert.That(ld.Variant2, Is.EqualTo("V2"));
    }

    [Test]
    public void FindHaplotypeBlocks_HighLD_CreatesBlock()
    {
        // Create variants in strong LD
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 200, genotypes),
            ("V3", 300, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        Assert.That(blocks, Has.Count.GreaterThanOrEqualTo(1));
    }

    [Test]
    public void FindHaplotypeBlocks_LowLD_NoBlock()
    {
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, new List<int> { 0, 1, 2, 0, 1, 2 }),
            ("V2", 200, new List<int> { 2, 1, 0, 2, 1, 0 })
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.9).ToList();

        Assert.That(blocks, Is.Empty);
    }

    #endregion

    #region Selection Tests

    [Test]
    public void CalculateIHS_BalancedEHH_ReturnsNearZero()
    {
        var ehh0 = new List<double> { 1.0, 0.8, 0.5, 0.2, 0.1 };
        var ehh1 = new List<double> { 1.0, 0.8, 0.5, 0.2, 0.1 };
        var positions = new List<int> { 0, 1000, 2000, 3000, 4000 };

        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions);

        Assert.That(Math.Abs(ihs), Is.LessThan(0.1));
    }

    [Test]
    public void CalculateIHS_ExtendedDerived_ReturnsPositive()
    {
        var ehh0 = new List<double> { 1.0, 0.5, 0.1, 0.05, 0.01 }; // Rapid decay
        var ehh1 = new List<double> { 1.0, 0.9, 0.8, 0.7, 0.6 };   // Extended
        var positions = new List<int> { 0, 1000, 2000, 3000, 4000 };

        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions);

        Assert.That(ihs, Is.GreaterThan(0));
    }

    [Test]
    public void ScanForSelection_NegativeTajimaD_DetectsSignal()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, -2.5, 0.1, 0.5) // Significant Tajima's D
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(
            regions, tajimaDThreshold: -2.0).ToList();

        Assert.That(signals, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(signals.Any(s => s.TestType == "TajimasD"), Is.True);
    }

    [Test]
    public void ScanForSelection_HighFst_DetectsSignal()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, 0, 0.5, 0) // High Fst
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(
            regions, fstThreshold: 0.25).ToList();

        Assert.That(signals.Any(s => s.TestType == "Fst"), Is.True);
    }

    [Test]
    public void ScanForSelection_NoSignificantSignals_ReturnsEmpty()
    {
        var regions = new List<(string, int, int, double, double, double)>
        {
            ("Region1", 0, 10000, 0, 0.1, 0.5) // All neutral
        };

        var signals = PopulationGeneticsAnalyzer.ScanForSelection(regions).ToList();

        Assert.That(signals, Is.Empty);
    }

    #endregion

    #region Ancestry Analysis Tests

    [Test]
    public void EstimateAncestry_SinglePopulation_Returns100Percent()
    {
        var individuals = new List<(string, IReadOnlyList<int>)>
        {
            ("IND1", new List<int> { 2, 2, 2, 2, 2 })
        };

        var refPops = new List<(string, IReadOnlyList<double>)>
        {
            ("POP1", new List<double> { 1.0, 1.0, 1.0, 1.0, 1.0 })
        };

        var ancestry = PopulationGeneticsAnalyzer.EstimateAncestry(
            individuals, refPops, maxIterations: 10).ToList();

        Assert.That(ancestry, Has.Count.EqualTo(1));
        Assert.That(ancestry[0].Proportions["POP1"], Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void EstimateAncestry_TwoPopulations_SumsToOne()
    {
        var individuals = new List<(string, IReadOnlyList<int>)>
        {
            ("IND1", new List<int> { 1, 1, 1, 1, 1 })
        };

        var refPops = new List<(string, IReadOnlyList<double>)>
        {
            ("POP1", new List<double> { 0.9, 0.9, 0.9, 0.9, 0.9 }),
            ("POP2", new List<double> { 0.1, 0.1, 0.1, 0.1, 0.1 })
        };

        var ancestry = PopulationGeneticsAnalyzer.EstimateAncestry(
            individuals, refPops).ToList();

        double sum = ancestry[0].Proportions.Values.Sum();
        Assert.That(sum, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void EstimateAncestry_EmptyInput_ReturnsEmpty()
    {
        var ancestry = PopulationGeneticsAnalyzer.EstimateAncestry(
            new List<(string, IReadOnlyList<int>)>(),
            new List<(string, IReadOnlyList<double>)>()).ToList();

        Assert.That(ancestry, Is.Empty);
    }

    #endregion

    #region Inbreeding Tests

    [Test]
    public void CalculateInbreedingFromROH_NoROH_ReturnsZero()
    {
        var roh = new List<(int, int)>();

        double f = PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(roh, 300_000_000);

        Assert.That(f, Is.EqualTo(0));
    }

    [Test]
    public void CalculateInbreedingFromROH_WithROH_CalculatesCorrectly()
    {
        var roh = new List<(int, int)>
        {
            (0, 10_000_000),
            (50_000_000, 60_000_000)
        };

        double f = PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(roh, 100_000_000);

        Assert.That(f, Is.EqualTo(0.2).Within(0.001)); // 20M / 100M
    }

    [Test]
    public void FindROH_LongHomozygousRun_DetectsROH()
    {
        // Create 100 homozygous SNPs spanning 2Mb
        var genotypes = Enumerable.Range(0, 100)
            .Select(i => (Position: i * 20000, Genotype: 0))
            .ToList();

        var roh = PopulationGeneticsAnalyzer.FindROH(
            genotypes,
            minSnps: 50,
            minLength: 1_000_000,
            maxHeterozygotes: 1).ToList();

        Assert.That(roh, Has.Count.EqualTo(1));
    }

    [Test]
    public void FindROH_TooManyHeterozygotes_NoROH()
    {
        var genotypes = Enumerable.Range(0, 100)
            .Select(i => (Position: i * 20000, Genotype: i % 5 == 0 ? 1 : 0))
            .ToList();

        var roh = PopulationGeneticsAnalyzer.FindROH(
            genotypes,
            maxHeterozygotes: 1).ToList();

        Assert.That(roh, Is.Empty);
    }

    [Test]
    public void FindROH_ShortRun_NotDetected()
    {
        var genotypes = Enumerable.Range(0, 20)
            .Select(i => (Position: i * 1000, Genotype: 0))
            .ToList();

        var roh = PopulationGeneticsAnalyzer.FindROH(
            genotypes,
            minSnps: 50).ToList();

        Assert.That(roh, Is.Empty);
    }

    #endregion

    #region Edge Cases Tests

    // Note: Diversity edge cases (CalculateNucleotideDiversity, CalculateWattersonTheta,
    // CalculateTajimasD) have been moved to PopulationGeneticsAnalyzer_Diversity_Tests.cs

    // Note: CalculateFst edge case tests (POP-FST-001) have been moved to
    // PopulationGeneticsAnalyzer_FStatistics_Tests.cs

    [Test]
    public void CalculateLD_EmptyGenotypes_ReturnsZeroLD()
    {
        var ld = PopulationGeneticsAnalyzer.CalculateLD(
            "V1", "V2", new List<(int, int)>(), 1000);

        Assert.That(ld.RSquared, Is.EqualTo(0));
        Assert.That(ld.DPrime, Is.EqualTo(0));
    }

    [Test]
    public void CalculateIHS_InsufficientData_ReturnsZero()
    {
        double ihs = PopulationGeneticsAnalyzer.CalculateIHS(
            new List<double> { 1.0 },
            new List<double> { 1.0 },
            new List<int> { 0 });

        Assert.That(ihs, Is.EqualTo(0));
    }

    [Test]
    public void FindHaplotypeBlocks_SingleVariant_NoBlocks()
    {
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, new List<int> { 0, 1, 2 })
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants).ToList();

        Assert.That(blocks, Is.Empty);
    }

    #endregion
}
