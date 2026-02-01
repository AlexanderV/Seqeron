using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class EpigeneticsAnalyzerTests
{
    #region CpG Site Detection Tests

    [Test]
    public void FindCpGSites_SimpleCpG_ReturnsCorrectPositions()
    {
        string sequence = "ACGTCGACG";
        //                 012345678
        // CpG at positions: 1 (CG at pos 1-2), 4 (CG at pos 4-5), 7 (CG at pos 7-8)

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Contains.Item(1)); // First CG at position 1
        Assert.That(sites, Contains.Item(4)); // Second CG at position 4
        Assert.That(sites, Contains.Item(7)); // Third CG at position 7
    }

    [Test]
    public void FindCpGSites_NoCpG_ReturnsEmpty()
    {
        string sequence = "AATTAATT";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindCpGSites_NullSequence_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites(null!).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindCpGSites_LowercaseSequence_Works()
    {
        string sequence = "acgtcg";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(2));
    }

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

    #region CpG Island Analysis Tests

    [Test]
    public void CalculateCpGObservedExpected_HighCpGDensity_ReturnsHighRatio()
    {
        // CpG-rich sequence
        string sequence = "CGCGCGCGCGCGCGCGCGCG";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.GreaterThan(0.6));
    }

    [Test]
    public void CalculateCpGObservedExpected_LowCpGDensity_ReturnsLowRatio()
    {
        // AT-rich sequence with rare CpG
        string sequence = "AATTAATTAATTAATTAATT";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.EqualTo(0));
    }

    [Test]
    public void CalculateCpGObservedExpected_NullSequence_ReturnsZero()
    {
        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(null!);

        Assert.That(ratio, Is.EqualTo(0));
    }

    [Test]
    public void FindCpGIslands_CpGRichRegion_DetectsIsland()
    {
        // Create a CpG island: high GC, high CpG O/E ratio
        string cpgIsland = string.Concat(Enumerable.Repeat("CGCG", 100)); // 400bp CpG-rich

        var islands = EpigeneticsAnalyzer.FindCpGIslands(
            cpgIsland,
            minLength: 100,
            minGc: 0.5,
            minCpGRatio: 0.6).ToList();

        Assert.That(islands, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(islands[0].GcContent, Is.EqualTo(1.0));
    }

    [Test]
    public void FindCpGIslands_NoCpGRichRegion_ReturnsEmpty()
    {
        // AT-rich sequence
        string atRich = string.Concat(Enumerable.Repeat("AATT", 100));

        var islands = EpigeneticsAnalyzer.FindCpGIslands(atRich).ToList();

        Assert.That(islands, Is.Empty);
    }

    [Test]
    public void FindCpGIslands_ShortSequence_ReturnsEmpty()
    {
        string shortSeq = "CGCGCG";

        var islands = EpigeneticsAnalyzer.FindCpGIslands(shortSeq, minLength: 200).ToList();

        Assert.That(islands, Is.Empty);
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

    #region Chromatin State Prediction Tests

    [Test]
    public void PredictChromatinState_ActivePromoter_Detected()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.8,
            h3k4me1: 0.1,
            h3k27ac: 0.8,
            h3k36me3: 0.1,
            h3k27me3: 0.0,
            h3k9me3: 0.0);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.ActivePromoter));
    }

    [Test]
    public void PredictChromatinState_ActiveEnhancer_Detected()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.1,
            h3k4me1: 0.8,
            h3k27ac: 0.8,
            h3k36me3: 0.1,
            h3k27me3: 0.0,
            h3k9me3: 0.0);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.ActiveEnhancer));
    }

    [Test]
    public void PredictChromatinState_Transcribed_Detected()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.1,
            h3k4me1: 0.1,
            h3k27ac: 0.1,
            h3k36me3: 0.8,
            h3k27me3: 0.0,
            h3k9me3: 0.0);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.Transcribed));
    }

    [Test]
    public void PredictChromatinState_Repressed_Detected()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.0,
            h3k4me1: 0.0,
            h3k27ac: 0.0,
            h3k36me3: 0.0,
            h3k27me3: 0.8,
            h3k9me3: 0.0);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.Repressed));
    }

    [Test]
    public void PredictChromatinState_Heterochromatin_Detected()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.0,
            h3k4me1: 0.0,
            h3k27ac: 0.0,
            h3k36me3: 0.0,
            h3k27me3: 0.0,
            h3k9me3: 0.9);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.Heterochromatin));
    }

    [Test]
    public void PredictChromatinState_LowSignal_Default()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.1,
            h3k4me1: 0.1,
            h3k27ac: 0.1,
            h3k36me3: 0.1,
            h3k27me3: 0.1,
            h3k9me3: 0.1);

        Assert.That(state, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.LowSignal));
    }

    [Test]
    public void AnnotateHistoneModifications_ReturnsAnnotations()
    {
        var modifications = new List<(int, int, string, double)>
        {
            (0, 1000, "H3K4me3", 0.8),
            (2000, 3000, "H3K27me3", 0.7)
        };

        var annotations = EpigeneticsAnalyzer.AnnotateHistoneModifications(modifications).ToList();

        Assert.That(annotations, Has.Count.EqualTo(2));
        Assert.That(annotations[0].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.ActivePromoter));
        Assert.That(annotations[1].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.Repressed));
    }

    #endregion

    #region Chromatin Accessibility Tests

    [Test]
    public void FindAccessibleRegions_HighSignalRegion_Detected()
    {
        var signal = Enumerable.Range(0, 200)
            .Select(i => (i * 10, 0.7))
            .ToList();

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal,
            threshold: 0.5,
            minWidth: 100).ToList();

        Assert.That(regions, Has.Count.GreaterThan(0));
        Assert.That(regions[0].AccessibilityScore, Is.EqualTo(0.7));
    }

    [Test]
    public void FindAccessibleRegions_LowSignal_NoRegions()
    {
        var signal = Enumerable.Range(0, 200)
            .Select(i => (i * 10, 0.2))
            .ToList();

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal,
            threshold: 0.5).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void FindAccessibleRegions_ClassifiesPeakStrength()
    {
        var signal = new List<(int, double)>
        {
            (0, 0.9), (10, 0.9), (20, 0.9), (30, 0.9), (40, 0.9),
            (50, 0.9), (60, 0.9), (70, 0.9), (80, 0.9), (90, 0.9),
            (100, 0.9), (110, 0.9)
        };

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal,
            threshold: 0.5,
            minWidth: 100).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].PeakType, Is.EqualTo("Strong"));
    }

    [Test]
    public void FindAccessibleRegions_EmptySignal_ReturnsEmpty()
    {
        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            new List<(int, double)>()).ToList();

        Assert.That(regions, Is.Empty);
    }

    #endregion

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
    public void FindCpGSites_CpGAtEnd_Detected()
    {
        string sequence = "AACG";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0], Is.EqualTo(2));
    }

    [Test]
    public void FindCpGSites_SingleNucleotide_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites("C").ToList();

        Assert.That(sites, Is.Empty);
    }

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

    [Test]
    public void AnnotateHistoneModifications_LowSignal_ReturnsLowSignalState()
    {
        var modifications = new List<(int, int, string, double)>
        {
            (0, 1000, "H3K4me3", 0.1) // Below threshold
        };

        var annotations = EpigeneticsAnalyzer.AnnotateHistoneModifications(modifications).ToList();

        Assert.That(annotations[0].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.LowSignal));
    }

    [Test]
    public void FindAccessibleRegions_GapInSignal_CreatesSeparateRegions()
    {
        var signal = new List<(int, double)>();

        // First region
        for (int i = 0; i <= 150; i += 10)
            signal.Add((i, 0.8));

        // Gap > maxGap (150)
        for (int i = 300; i <= 450; i += 10)
            signal.Add((i, 0.8));

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal,
            threshold: 0.5,
            minWidth: 100,
            maxGap: 50).ToList();

        Assert.That(regions, Has.Count.EqualTo(2));
    }

    #endregion
}
