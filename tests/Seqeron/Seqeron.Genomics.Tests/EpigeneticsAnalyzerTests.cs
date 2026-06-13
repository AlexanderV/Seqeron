using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class EpigeneticsAnalyzerTests
{
    #region Methylation Site Detection Tests

    [Test]
    public void FindMethylationSites_IdentifiesAllContexts()
    {
        // CpG at 0, CHG at 3 (CAG), CHH at 6 (CAA)
        string sequence = "CGACAGCAA";

        var sites = EpigeneticsAnalyzer.FindMethylationSites(sequence).ToList();

        var cpg = sites.FirstOrDefault(s => s.Type == EpigeneticsAnalyzer.MethylationType.CpG);
        var chg = sites.FirstOrDefault(s => s.Type == EpigeneticsAnalyzer.MethylationType.CHG);
        var chh = sites.FirstOrDefault(s => s.Type == EpigeneticsAnalyzer.MethylationType.CHH);

        Assert.That(cpg.Position, Is.EqualTo(0));
        Assert.That(chg.Position, Is.EqualTo(3));
        Assert.That(chh.Position, Is.EqualTo(6));
    }

    #endregion

    #region Bisulfite Conversion Tests

    [Test]
    public void SimulateBisulfiteConversion_UnmethylatedCytosines_ConvertsToThymine()
    {
        string sequence = "ACGT";

        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence);

        Assert.That(converted, Is.EqualTo("ATGT"));
    }

    [Test]
    public void SimulateBisulfiteConversion_MethylatedCytosines_Protected()
    {
        string sequence = "ACGT";
        var methylated = new HashSet<int> { 1 }; // C at position 1 is methylated

        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, methylated);

        Assert.That(converted, Is.EqualTo("ACGT")); // C is protected
    }

    [Test]
    public void SimulateBisulfiteConversion_MultipleCytosines_SelectiveConversion()
    {
        string sequence = "CCCC";
        var methylated = new HashSet<int> { 0, 2 }; // Positions 0 and 2 methylated

        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, methylated);

        Assert.That(converted, Is.EqualTo("CTCT"));
    }

    [Test]
    public void SimulateBisulfiteConversion_EmptySequence_ReturnsEmpty()
    {
        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion("");

        Assert.That(converted, Is.Empty);
    }

    [Test]
    public void SimulateBisulfiteConversion_PreservesCase()
    {
        string sequence = "AcGt";
        var methylated = new HashSet<int>();

        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, methylated);

        Assert.That(converted, Is.EqualTo("AtGt"));
    }

    #endregion

    #region Methylation from Bisulfite Tests

    [Test]
    public void CalculateMethylationFromBisulfite_FullyMethylated_Returns100Percent()
    {
        string reference = "ACGTCGACGT"; // CpG at positions 2 and 5
        var reads = new List<(string, int)>
        {
            ("ACGTCGACGT", 0), // C preserved = methylated
            ("ACGTCGACGT", 0),
            ("ACGTCGACGT", 0)
        };

        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, reads).ToList();

        Assert.That(sites, Has.Count.GreaterThan(0));
        Assert.That(sites.All(s => s.MethylationLevel == 1.0), Is.True);
    }

    [Test]
    public void CalculateMethylationFromBisulfite_UnmethylatedCpG_ReturnsZero()
    {
        string reference = "ACGTCG"; // CpG at positions 2 and 4
        var reads = new List<(string, int)>
        {
            ("ATGTTG", 0), // C converted to T = unmethylated
            ("ATGTTG", 0),
            ("ATGTTG", 0)
        };

        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, reads).ToList();

        Assert.That(sites, Has.Count.GreaterThan(0));
        Assert.That(sites.All(s => s.MethylationLevel == 0.0), Is.True);
    }

    [Test]
    public void CalculateMethylationFromBisulfite_PartialMethylation_ReturnsCorrectLevel()
    {
        string reference = "ACGT"; // CpG at position 2
        var reads = new List<(string, int)>
        {
            ("ACGT", 0), // Methylated
            ("ATGT", 0), // Unmethylated
            ("ACGT", 0), // Methylated
            ("ATGT", 0)  // Unmethylated
        };

        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, reads).ToList();

        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].MethylationLevel, Is.EqualTo(0.5));
    }

    #endregion

    #region Methylation Profile Tests

    [Test]
    public void GenerateMethylationProfile_MixedSites_CorrectAverages()
    {
        var sites = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.8, 10),
            new(10, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.6, 10),
            new(5, EpigeneticsAnalyzer.MethylationType.CHG, "CAG", 0.3, 10),
            new(15, EpigeneticsAnalyzer.MethylationType.CHH, "CAA", 0.1, 10)
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.That(profile.CpGMethylation, Is.EqualTo(0.7).Within(0.001));
        Assert.That(profile.CHGMethylation, Is.EqualTo(0.3).Within(0.001));
        Assert.That(profile.CHHMethylation, Is.EqualTo(0.1).Within(0.001));
        Assert.That(profile.TotalCpGSites, Is.EqualTo(2));
        Assert.That(profile.MethylatedCpGSites, Is.EqualTo(2)); // Both >= 0.5
    }

    [Test]
    public void GenerateMethylationProfile_NoSites_ReturnsZeros()
    {
        var sites = new List<EpigeneticsAnalyzer.MethylationSite>();

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.That(profile.GlobalMethylation, Is.EqualTo(0));
        Assert.That(profile.TotalCpGSites, Is.EqualTo(0));
    }

    [Test]
    public void GenerateMethylationProfile_TracksPositions()
    {
        var sites = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(50, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10)
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.That(profile.MethylationByPosition, Has.Count.EqualTo(2));
        Assert.That(profile.MethylationByPosition[0].Position, Is.EqualTo(50)); // Sorted
        Assert.That(profile.MethylationByPosition[1].Position, Is.EqualTo(100));
    }

    #endregion

    #region Differentially Methylated Region Tests

    [Test]
    public void FindDMRs_SignificantDifference_DetectsRegion()
    {
        var sample1 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10),
            new(300, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10)
        };

        var sample2 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(300, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10)
        };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(
            sample1, sample2,
            windowSize: 1000,
            minDifference: 0.25,
            minCpGCount: 3).ToList();

        Assert.That(dmrs, Has.Count.EqualTo(1));
        Assert.That(dmrs[0].Annotation, Is.EqualTo("Hypermethylated"));
        Assert.That(dmrs[0].MeanDifference, Is.EqualTo(0.8).Within(0.01));
    }

    [Test]
    public void FindDMRs_NoDifference_ReturnsEmpty()
    {
        var sample1 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10)
        };

        var sample2 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.5, 10)
        };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(sample1, sample2, minDifference: 0.25).ToList();

        Assert.That(dmrs, Is.Empty);
    }

    [Test]
    public void FindDMRs_Hypomethylated_AnnotatesCorrectly()
    {
        var sample1 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.9, 10)
        };

        var sample2 = new List<EpigeneticsAnalyzer.MethylationSite>
        {
            new(0, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10),
            new(100, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10),
            new(200, EpigeneticsAnalyzer.MethylationType.CpG, "CGA", 0.1, 10)
        };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(sample1, sample2).ToList();

        Assert.That(dmrs, Has.Count.EqualTo(1));
        Assert.That(dmrs[0].Annotation, Is.EqualTo("Hypomethylated"));
    }

    #endregion

    // NOTE: Chromatin State Prediction and Chromatin Accessibility tests were moved to
    // the canonical evidence-based fixture EpigeneticsAnalyzer_ChromatinState_Tests.cs
    // (EPIGEN-CHROM-001).

    #region Imprinting Analysis Tests

    [Test]
    public void PredictImprintedGenes_MaternallyExpressed_Detected()
    {
        var genes = new List<(string, int, int, double, double)>
        {
            ("GENE1", 0, 1000, 0.9, 0.1) // Maternal methylated, paternal unmethylated
        };

        var imprinted = EpigeneticsAnalyzer.PredictImprintedGenes(genes, minDifference: 0.4).ToList();

        Assert.That(imprinted, Has.Count.EqualTo(1));
        Assert.That(imprinted[0].ParentalOrigin, Is.EqualTo("Maternal"));
        Assert.That(imprinted[0].HasDMR, Is.True);
    }

    [Test]
    public void PredictImprintedGenes_PaternallyExpressed_Detected()
    {
        var genes = new List<(string, int, int, double, double)>
        {
            ("GENE2", 0, 1000, 0.1, 0.9) // Paternal methylated
        };

        var imprinted = EpigeneticsAnalyzer.PredictImprintedGenes(genes, minDifference: 0.4).ToList();

        Assert.That(imprinted, Has.Count.EqualTo(1));
        Assert.That(imprinted[0].ParentalOrigin, Is.EqualTo("Paternal"));
    }

    [Test]
    public void PredictImprintedGenes_NoDifference_NotImprinted()
    {
        var genes = new List<(string, int, int, double, double)>
        {
            ("GENE3", 0, 1000, 0.5, 0.5)
        };

        var imprinted = EpigeneticsAnalyzer.PredictImprintedGenes(genes, minDifference: 0.4).ToList();

        Assert.That(imprinted, Is.Empty);
    }

    [Test]
    public void PredictImprintedGenes_CalculatesImprintingScore()
    {
        var genes = new List<(string, int, int, double, double)>
        {
            ("GENE1", 0, 1000, 1.0, 0.0) // Maximum difference
        };

        var imprinted = EpigeneticsAnalyzer.PredictImprintedGenes(genes).ToList();

        Assert.That(imprinted[0].ImprintingScore, Is.GreaterThan(0.9));
    }

    #endregion

    #region Epigenetic Age Tests

    [Test]
    public void CalculateEpigeneticAge_WithDefaultCoefficients_ReturnsAge()
    {
        var methylation = new Dictionary<string, double>
        {
            { "cg00000029", 0.5 },
            { "cg00000165", 0.3 },
            { "cg00000236", 0.7 }
        };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);

        Assert.That(age, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void CalculateEpigeneticAge_EmptyInput_ReturnsZero()
    {
        var methylation = new Dictionary<string, double>();

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);

        Assert.That(age, Is.EqualTo(0));
    }

    [Test]
    public void CalculateEpigeneticAge_CustomCoefficients_Works()
    {
        var methylation = new Dictionary<string, double>
        {
            { "custom_cpg1", 0.5 },
            { "custom_cpg2", 0.5 }
        };

        var coefficients = new Dictionary<string, double>
        {
            { "custom_cpg1", 1.0 },
            { "custom_cpg2", 1.0 }
        };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation, coefficients);

        // exp(1.0) - 1 ≈ 1.718
        Assert.That(age, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateEpigeneticAge_UnknownCpGs_Ignored()
    {
        var methylation = new Dictionary<string, double>
        {
            { "unknown_cpg", 1.0 }
        };

        double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylation);

        Assert.That(age, Is.EqualTo(0)); // No matching coefficients
    }

    #endregion

    #region Edge Cases Tests

    [Test]
    public void FindMethylationSites_EmptySequence_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindMethylationSites("").ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindDMRs_EmptySamples_ReturnsEmpty()
    {
        var dmrs = EpigeneticsAnalyzer.FindDMRs(
            new List<EpigeneticsAnalyzer.MethylationSite>(),
            new List<EpigeneticsAnalyzer.MethylationSite>()).ToList();

        Assert.That(dmrs, Is.Empty);
    }

    #endregion
}
