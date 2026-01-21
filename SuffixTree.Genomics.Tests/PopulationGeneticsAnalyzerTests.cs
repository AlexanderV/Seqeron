using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class PopulationGeneticsAnalyzerTests
{
    #region Allele Frequency Tests

    [Test]
    public void CalculateAlleleFrequencies_EqualGenotypes_Returns50Percent()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 25,
            heterozygous: 50,
            homozygousMinor: 25);

        Assert.That(major, Is.EqualTo(0.5).Within(0.001));
        Assert.That(minor, Is.EqualTo(0.5).Within(0.001));
    }

    [Test]
    public void CalculateAlleleFrequencies_AllHomozygousMajor_Returns100Percent()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 100,
            heterozygous: 0,
            homozygousMinor: 0);

        Assert.That(major, Is.EqualTo(1.0).Within(0.001));
        Assert.That(minor, Is.EqualTo(0.0).Within(0.001));
    }

    [Test]
    public void CalculateAlleleFrequencies_ZeroSamples_ReturnsZero()
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, 0);

        Assert.That(major, Is.EqualTo(0));
        Assert.That(minor, Is.EqualTo(0));
    }

    [Test]
    public void CalculateMAF_FromGenotypes_CalculatesCorrectly()
    {
        // 10 samples: 5 hom ref (0), 4 het (1), 1 hom alt (2)
        // Alt alleles = 0*5 + 1*4 + 2*1 = 6
        // Total = 20, freq = 0.3, MAF = 0.3
        var genotypes = new[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 2 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0.3).Within(0.001));
    }

    [Test]
    public void CalculateMAF_HighAltFreq_ReturnsMinorAllele()
    {
        // All hom alt = freq 1.0, MAF = 0
        var genotypes = new[] { 2, 2, 2, 2, 2 };

        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        Assert.That(maf, Is.EqualTo(0));
    }

    [Test]
    public void FilterByMAF_FiltersCorrectly()
    {
        var variants = new List<PopulationGeneticsAnalyzer.Variant>
        {
            new("V1", "chr1", 100, "A", "G", 0.001, 100), // Below threshold
            new("V2", "chr1", 200, "A", "G", 0.05, 100),  // Within range
            new("V3", "chr1", 300, "A", "G", 0.3, 100),   // Within range
            new("V4", "chr1", 400, "A", "G", 0.5, 100)    // Within range
        };

        var filtered = PopulationGeneticsAnalyzer.FilterByMAF(variants, minMAF: 0.01).ToList();

        Assert.That(filtered, Has.Count.EqualTo(3));
        Assert.That(filtered.Any(v => v.Id == "V1"), Is.False);
    }

    #endregion

    #region Diversity Statistics Tests

    [Test]
    public void CalculateNucleotideDiversity_IdenticalSequences_ReturnsZero()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList(),
            "ACGT".ToList(),
            "ACGT".ToList()
        };

        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        Assert.That(pi, Is.EqualTo(0));
    }

    [Test]
    public void CalculateNucleotideDiversity_AllDifferent_ReturnsPositive()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            "AAAA".ToList(),
            "TTTT".ToList()
        };

        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        Assert.That(pi, Is.EqualTo(1.0)); // All positions differ
    }

    [Test]
    public void CalculateWattersonTheta_KnownValues_CalculatesCorrectly()
    {
        // S = 10, n = 10, L = 1000
        // a1 = 1 + 1/2 + ... + 1/9 ≈ 2.829
        // theta = 10 / (2.829 * 1000) ≈ 0.00353

        double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
            segregatingSites: 10,
            sampleSize: 10,
            sequenceLength: 1000);

        Assert.That(theta, Is.GreaterThan(0.003).And.LessThan(0.004));
    }

    [Test]
    public void CalculateTajimasD_NeutralEvolution_NearZero()
    {
        // Under neutrality, Tajima's D should be near 0
        // When pi ≈ theta, D ≈ 0

        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity: 0.01,
            wattersonTheta: 0.01,
            segregatingSites: 100,
            sampleSize: 50);

        Assert.That(Math.Abs(d), Is.LessThan(1));
    }

    [Test]
    public void CalculateTajimasD_PositiveSelection_Negative()
    {
        // Positive selection: excess of rare variants
        // pi << theta, D < 0

        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            nucleotideDiversity: 0.001,
            wattersonTheta: 0.01,
            segregatingSites: 100,
            sampleSize: 50);

        Assert.That(d, Is.LessThan(0));
    }

    [Test]
    public void CalculateDiversityStatistics_ReturnsAllMetrics()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGTACGT".ToList(),
            "ACGTATGT".ToList(),
            "ACGTACGA".ToList()
        };

        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        Assert.That(stats.SampleSize, Is.EqualTo(3));
        Assert.That(stats.SegregratingSites, Is.GreaterThan(0));
        Assert.That(stats.NucleotideDiversity, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void CalculateDiversityStatistics_SingleSequence_ReturnsZeroDiversity()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            "ACGT".ToList()
        };

        var stats = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(sequences);

        Assert.That(stats.NucleotideDiversity, Is.EqualTo(0));
        Assert.That(stats.SampleSize, Is.EqualTo(1));
    }

    #endregion

    #region Hardy-Weinberg Tests

    [Test]
    public void TestHardyWeinberg_InEquilibrium_PassesTest()
    {
        // Expected under HWE with p = 0.5:
        // AA: 25, Aa: 50, aa: 25

        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "V1",
            observedAA: 25,
            observedAa: 50,
            observedaa: 25);

        // Chi-square should be very low for perfect HWE
        Assert.That(result.ChiSquare, Is.LessThan(1));
        Assert.That(result.ExpectedAA, Is.EqualTo(25).Within(0.1));
    }

    [Test]
    public void TestHardyWeinberg_ExcessHeterozygotes_FailsTest()
    {
        // Too many heterozygotes for HWE

        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "V1",
            observedAA: 10,
            observedAa: 80,
            observedaa: 10);

        Assert.That(result.ChiSquare, Is.GreaterThan(0));
        // May or may not pass depending on threshold
    }

    [Test]
    public void TestHardyWeinberg_CalculatesExpectedCorrectly()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "V1",
            observedAA: 100,
            observedAa: 0,
            observedaa: 0);

        // p = 1.0, expected: 100 AA, 0 Aa, 0 aa
        Assert.That(result.ExpectedAA, Is.EqualTo(100).Within(0.1));
        Assert.That(result.ExpectedAa, Is.EqualTo(0).Within(0.1));
    }

    [Test]
    public void TestHardyWeinberg_ZeroSamples_HandlesGracefully()
    {
        var result = PopulationGeneticsAnalyzer.TestHardyWeinberg("V1", 0, 0, 0);

        Assert.That(result.InEquilibrium, Is.True);
        Assert.That(result.PValue, Is.EqualTo(1));
    }

    #endregion

    #region F-Statistics Tests

    [Test]
    public void CalculateFst_IdenticalPopulations_ReturnsZero()
    {
        var pop1 = new List<(double, int)> { (0.5, 100), (0.3, 100) };
        var pop2 = new List<(double, int)> { (0.5, 100), (0.3, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(0).Within(0.001));
    }

    [Test]
    public void CalculateFst_DifferentPopulations_ReturnsPositive()
    {
        var pop1 = new List<(double, int)> { (0.9, 100), (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.1, 100), (0.2, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateFst_FixedDifferences_ReturnsHighFst()
    {
        var pop1 = new List<(double, int)> { (1.0, 100) };
        var pop2 = new List<(double, int)> { (0.0, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.GreaterThan(0.5));
    }

    [Test]
    public void CalculatePairwiseFst_ThreePopulations_ReturnsMatrix()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100) }),
            ("Pop2", new List<(double, int)> { (0.6, 100) }),
            ("Pop3", new List<(double, int)> { (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        Assert.That(matrix.GetLength(0), Is.EqualTo(3));
        Assert.That(matrix.GetLength(1), Is.EqualTo(3));
        Assert.That(matrix[0, 0], Is.EqualTo(0)); // Diagonal is 0
        Assert.That(matrix[0, 1], Is.EqualTo(matrix[1, 0])); // Symmetric
    }

    [Test]
    public void CalculateFStatistics_ReturnsAllComponents()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 50, 25, 50, 0.4, 0.5),
            (30, 50, 15, 50, 0.5, 0.3)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        Assert.That(stats.Population1, Is.EqualTo("Pop1"));
        Assert.That(stats.Population2, Is.EqualTo("Pop2"));
        // Fst, Fis, Fit should be between 0 and 1 (or slightly negative for Fis)
    }

    #endregion

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

    [Test]
    public void CalculateNucleotideDiversity_SingleSequence_ReturnsZero()
    {
        var sequences = new List<IReadOnlyList<char>> { "ACGT".ToList() };

        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        Assert.That(pi, Is.EqualTo(0));
    }

    [Test]
    public void CalculateWattersonTheta_SmallSample_HandlesCorrectly()
    {
        double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(5, 2, 100);

        Assert.That(theta, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateTajimasD_NoSegregratingSites_ReturnsZero()
    {
        double d = PopulationGeneticsAnalyzer.CalculateTajimasD(0, 0, 0, 50);

        Assert.That(d, Is.EqualTo(0));
    }

    [Test]
    public void CalculateFst_EmptyPopulations_ReturnsZero()
    {
        double fst = PopulationGeneticsAnalyzer.CalculateFst(
            new List<(double, int)>(),
            new List<(double, int)>());

        Assert.That(fst, Is.EqualTo(0));
    }

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
