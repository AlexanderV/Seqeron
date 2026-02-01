using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class TranscriptomeAnalyzerTests
{
    #region TPM and FPKM Tests

    [Test]
    public void CalculateTPM_EqualGenes_EqualTPM()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 1000, 1000),
            ("GENE2", 1000, 1000)
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        Assert.That(results, Has.Count.EqualTo(2));
        Assert.That(results[0].TPM, Is.EqualTo(500000).Within(1));
        Assert.That(results[1].TPM, Is.EqualTo(500000).Within(1));
    }

    [Test]
    public void CalculateTPM_DifferentLengths_AdjustsForLength()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 1000, 1000),  // Rate = 1.0
            ("GENE2", 1000, 2000)   // Rate = 0.5
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        // Gene1 should have higher TPM (shorter gene with same counts)
        Assert.That(results[0].TPM, Is.GreaterThan(results[1].TPM));
    }

    [Test]
    public void CalculateTPM_SumToMillion()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 500, 1000),
            ("GENE2", 1000, 1500),
            ("GENE3", 2000, 2000)
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        double totalTPM = results.Sum(r => r.TPM);
        Assert.That(totalTPM, Is.EqualTo(1_000_000).Within(1));
    }

    [Test]
    public void CalculateTPM_EmptyInput_ReturnsEmpty()
    {
        var results = TranscriptomeAnalyzer.CalculateTPM(
            new List<(string, double, int)>()).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void CalculateTPM_ZeroReads_ReturnsZeroTPM()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 0, 1000),
            ("GENE2", 0, 1000)
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        Assert.That(results.All(r => r.TPM == 0), Is.True);
    }

    [Test]
    public void CalculateTPM_IncludesFPKM()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 1000, 1000)
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        Assert.That(results[0].FPKM, Is.GreaterThan(0));
    }

    #endregion

    #region Quantile Normalization Tests

    [Test]
    public void QuantileNormalize_TwoSamples_EqualDistribution()
    {
        var samples = new List<List<double>>
        {
            new() { 1, 2, 3 },
            new() { 4, 5, 6 }
        };

        var normalized = TranscriptomeAnalyzer.QuantileNormalize(samples).ToList();

        Assert.That(normalized, Has.Count.EqualTo(2));
        // After normalization, both samples should have same distribution
        var sorted1 = normalized[0].OrderBy(x => x).ToList();
        var sorted2 = normalized[1].OrderBy(x => x).ToList();

        for (int i = 0; i < sorted1.Count; i++)
        {
            Assert.That(sorted1[i], Is.EqualTo(sorted2[i]).Within(0.001));
        }
    }

    [Test]
    public void QuantileNormalize_EmptyInput_ReturnsEmpty()
    {
        var normalized = TranscriptomeAnalyzer.QuantileNormalize(
            new List<List<double>>()).ToList();

        Assert.That(normalized, Is.Empty);
    }

    [Test]
    public void Log2Transform_PositiveValues_TransformsCorrectly()
    {
        var values = new List<double> { 1, 2, 4, 8 };

        var transformed = TranscriptomeAnalyzer.Log2Transform(values, pseudocount: 0).ToList();

        Assert.That(transformed[0], Is.EqualTo(0).Within(0.001));
        Assert.That(transformed[1], Is.EqualTo(1).Within(0.001));
        Assert.That(transformed[2], Is.EqualTo(2).Within(0.001));
        Assert.That(transformed[3], Is.EqualTo(3).Within(0.001));
    }

    [Test]
    public void Log2Transform_WithPseudocount_HandleZeros()
    {
        var values = new List<double> { 0, 1, 3 };

        var transformed = TranscriptomeAnalyzer.Log2Transform(values, pseudocount: 1).ToList();

        Assert.That(transformed[0], Is.EqualTo(0).Within(0.001)); // log2(0+1) = 0
        Assert.That(transformed[1], Is.EqualTo(1).Within(0.001)); // log2(1+1) = 1
        Assert.That(transformed[2], Is.EqualTo(2).Within(0.001)); // log2(3+1) = 2
    }

    #endregion

    #region Differential Expression Tests

    [Test]
    public void AnalyzeDifferentialExpression_SignificantChange_DetectsDE()
    {
        var data = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>
        {
            ("GENE1",
             new List<double> { 10, 12, 11, 13 },
             new List<double> { 100, 110, 105, 115 })
        };

        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(
            data, foldChangeThreshold: 1.0, pValueThreshold: 0.05).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Log2FoldChange, Is.GreaterThan(2)); // ~3x fold change
        Assert.That(results[0].Regulation, Is.EqualTo("Upregulated"));
    }

    [Test]
    public void AnalyzeDifferentialExpression_NoChange_NotSignificant()
    {
        var data = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>
        {
            ("GENE1",
             new List<double> { 100, 102, 98, 101 },
             new List<double> { 101, 99, 100, 102 })
        };

        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(
            data, foldChangeThreshold: 1.0).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].IsSignificant, Is.False);
    }

    [Test]
    public void AnalyzeDifferentialExpression_Downregulated_AnnotatesCorrectly()
    {
        var data = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>
        {
            ("GENE1",
             new List<double> { 100, 110, 105 },
             new List<double> { 10, 12, 11 })
        };

        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(data).ToList();

        Assert.That(results[0].Log2FoldChange, Is.LessThan(0));
        Assert.That(results[0].Regulation, Is.EqualTo("Downregulated"));
    }

    [Test]
    public void AnalyzeDifferentialExpression_MultipleGenes_AppliesBHCorrection()
    {
        var data = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>();

        for (int i = 0; i < 100; i++)
        {
            data.Add(($"GENE{i}",
                new List<double> { 100 + i, 102 + i, 98 + i },
                new List<double> { 101 + i, 99 + i, 100 + i }));
        }

        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(data).ToList();

        // Adjusted p-values should be >= raw p-values
        Assert.That(results.All(r => r.AdjustedPValue >= r.PValue), Is.True);
    }

    [Test]
    public void AnalyzeDifferentialExpression_EmptyInput_ReturnsEmpty()
    {
        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(
            new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>()).ToList();

        Assert.That(results, Is.Empty);
    }

    #endregion

    #region Gene Set Enrichment Tests

    [Test]
    public void PerformOverRepresentationAnalysis_EnrichedPathway_DetectsEnrichment()
    {
        var deGenes = new HashSet<string> { "A", "B", "C", "D", "E" };
        var pathways = new List<(string, string, IReadOnlySet<string>)>
        {
            ("P1", "Pathway1", new HashSet<string> { "A", "B", "C", "X", "Y" })
        };

        var results = TranscriptomeAnalyzer.PerformOverRepresentationAnalysis(
            deGenes, pathways, backgroundGeneCount: 1000).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].OverlappingGenes, Is.EqualTo(3));
        Assert.That(results[0].EnrichmentScore, Is.GreaterThan(1)); // Enriched
    }

    [Test]
    public void PerformOverRepresentationAnalysis_NoOverlap_NotIncluded()
    {
        var deGenes = new HashSet<string> { "A", "B", "C" };
        var pathways = new List<(string, string, IReadOnlySet<string>)>
        {
            ("P1", "Pathway1", new HashSet<string> { "X", "Y", "Z" })
        };

        var results = TranscriptomeAnalyzer.PerformOverRepresentationAnalysis(
            deGenes, pathways, backgroundGeneCount: 1000).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void CalculateEnrichmentScore_GenesAtTop_PositiveScore()
    {
        var rankedGenes = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H" };
        var geneSet = new HashSet<string> { "A", "B" }; // At top of list

        double score = TranscriptomeAnalyzer.CalculateEnrichmentScore(rankedGenes, geneSet);

        Assert.That(score, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateEnrichmentScore_GenesAtBottom_NegativeScore()
    {
        var rankedGenes = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H" };
        var geneSet = new HashSet<string> { "G", "H" }; // At bottom of list

        double score = TranscriptomeAnalyzer.CalculateEnrichmentScore(rankedGenes, geneSet);

        Assert.That(score, Is.LessThan(0));
    }

    [Test]
    public void CalculateEnrichmentScore_EmptyGeneSet_ReturnsZero()
    {
        var rankedGenes = new List<string> { "A", "B", "C" };

        double score = TranscriptomeAnalyzer.CalculateEnrichmentScore(
            rankedGenes, new HashSet<string>());

        Assert.That(score, Is.EqualTo(0));
    }

    #endregion

    #region Alternative Splicing Tests

    [Test]
    public void FindSkippedExonEvents_CalculatesPSI()
    {
        var exonData = new List<(string, int, int, double, double)>
        {
            ("GENE1", 100, 200, 80, 20) // 80% inclusion
        };

        var events = TranscriptomeAnalyzer.FindSkippedExonEvents(exonData).ToList();

        Assert.That(events, Has.Count.EqualTo(1));
        Assert.That(events[0].InclusionLevel, Is.EqualTo(0.8).Within(0.001));
        Assert.That(events[0].EventType, Is.EqualTo("SkippedExon"));
    }

    [Test]
    public void DetectDifferentialSplicing_SignificantChange_Detected()
    {
        var splicingData = new List<(string, int, int, double, double)>
        {
            ("GENE1", 100, 200, 0.3, 0.7) // DeltaPSI = 0.4
        };

        var events = TranscriptomeAnalyzer.DetectDifferentialSplicing(
            splicingData, deltaPsiThreshold: 0.1).ToList();

        Assert.That(events, Has.Count.EqualTo(1));
        Assert.That(events[0].DeltaPSI, Is.EqualTo(0.4).Within(0.001));
        Assert.That(events[0].EventType, Is.EqualTo("IncreasedInclusion"));
    }

    [Test]
    public void DetectDifferentialSplicing_SmallChange_NotDetected()
    {
        var splicingData = new List<(string, int, int, double, double)>
        {
            ("GENE1", 100, 200, 0.5, 0.52)
        };

        var events = TranscriptomeAnalyzer.DetectDifferentialSplicing(
            splicingData, deltaPsiThreshold: 0.1).ToList();

        Assert.That(events, Is.Empty);
    }

    [Test]
    public void DetectDifferentialSplicing_DecreasedInclusion_AnnotatesCorrectly()
    {
        var splicingData = new List<(string, int, int, double, double)>
        {
            ("GENE1", 100, 200, 0.8, 0.3)
        };

        var events = TranscriptomeAnalyzer.DetectDifferentialSplicing(splicingData).ToList();

        Assert.That(events[0].EventType, Is.EqualTo("IncreasedSkipping"));
    }

    #endregion

    #region Transcript Isoform Tests

    [Test]
    public void FindDominantIsoforms_IdentifiesHighestExpression()
    {
        var isoforms = new List<TranscriptomeAnalyzer.TranscriptIsoform>
        {
            new("TX1", "GENE1", 1000, 5, 100, true, new List<(int, int)>()),
            new("TX2", "GENE1", 800, 4, 50, true, new List<(int, int)>()),
            new("TX3", "GENE1", 1200, 6, 25, true, new List<(int, int)>())
        };

        var results = TranscriptomeAnalyzer.FindDominantIsoforms(isoforms).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].DominantIsoform.TranscriptId, Is.EqualTo("TX1"));
        Assert.That(results[0].DominanceRatio, Is.EqualTo(100.0 / 175.0).Within(0.01));
    }

    [Test]
    public void FindDominantIsoforms_MultipleGenes_GroupsCorrectly()
    {
        var isoforms = new List<TranscriptomeAnalyzer.TranscriptIsoform>
        {
            new("TX1", "GENE1", 1000, 5, 100, true, new List<(int, int)>()),
            new("TX2", "GENE2", 800, 4, 200, true, new List<(int, int)>())
        };

        var results = TranscriptomeAnalyzer.FindDominantIsoforms(isoforms).ToList();

        Assert.That(results, Has.Count.EqualTo(2));
    }

    [Test]
    public void DetectIsoformSwitching_SignificantSwitch_Detected()
    {
        var isoformData = new List<(TranscriptomeAnalyzer.TranscriptIsoform, double, double)>
        {
            (new TranscriptomeAnalyzer.TranscriptIsoform(
                "TX1", "GENE1", 1000, 5, 0, true, new List<(int, int)>()), 80, 20),
            (new TranscriptomeAnalyzer.TranscriptIsoform(
                "TX2", "GENE1", 800, 4, 0, true, new List<(int, int)>()), 20, 80)
        };

        var switches = TranscriptomeAnalyzer.DetectIsoformSwitching(
            isoformData, switchThreshold: 0.3).ToList();

        Assert.That(switches, Has.Count.EqualTo(1));
        Assert.That(switches[0].GeneId, Is.EqualTo("GENE1"));
    }

    #endregion

    #region Co-Expression Tests

    [Test]
    public void CalculatePearsonCorrelation_PerfectPositive_ReturnsOne()
    {
        var exp1 = new List<double> { 1, 2, 3, 4, 5 };
        var exp2 = new List<double> { 2, 4, 6, 8, 10 };

        double corr = TranscriptomeAnalyzer.CalculatePearsonCorrelation(exp1, exp2);

        Assert.That(corr, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculatePearsonCorrelation_PerfectNegative_ReturnsNegativeOne()
    {
        var exp1 = new List<double> { 1, 2, 3, 4, 5 };
        var exp2 = new List<double> { 5, 4, 3, 2, 1 };

        double corr = TranscriptomeAnalyzer.CalculatePearsonCorrelation(exp1, exp2);

        Assert.That(corr, Is.EqualTo(-1.0).Within(0.001));
    }

    [Test]
    public void CalculatePearsonCorrelation_NoCorrelation_ReturnsNearZero()
    {
        var exp1 = new List<double> { 1, 2, 1, 2, 1 };
        var exp2 = new List<double> { 1, 1, 2, 2, 1 };

        double corr = TranscriptomeAnalyzer.CalculatePearsonCorrelation(exp1, exp2);

        Assert.That(Math.Abs(corr), Is.LessThan(0.5));
    }

    [Test]
    public void CalculatePearsonCorrelation_DifferentLengths_ReturnsZero()
    {
        var exp1 = new List<double> { 1, 2, 3 };
        var exp2 = new List<double> { 1, 2 };

        double corr = TranscriptomeAnalyzer.CalculatePearsonCorrelation(exp1, exp2);

        Assert.That(corr, Is.EqualTo(0));
    }

    [Test]
    public void BuildCoExpressionNetwork_HighCorrelation_IncludesEdge()
    {
        var genes = new List<(string, IReadOnlyList<double>)>
        {
            ("GENE1", new List<double> { 1, 2, 3, 4, 5 }),
            ("GENE2", new List<double> { 2, 4, 6, 8, 10 }), // Perfect correlation with GENE1
            ("GENE3", new List<double> { 5, 4, 3, 2, 1 })   // Negative correlation
        };

        var edges = TranscriptomeAnalyzer.BuildCoExpressionNetwork(
            genes, correlationThreshold: 0.7).ToList();

        Assert.That(edges.Any(e =>
            (e.Gene1 == "GENE1" && e.Gene2 == "GENE2") ||
            (e.Gene1 == "GENE2" && e.Gene2 == "GENE1")), Is.True);
    }

    [Test]
    public void ClusterGenesByExpression_ReturnsRequestedClusters()
    {
        var genes = new List<(string, IReadOnlyList<double>)>();

        for (int i = 0; i < 20; i++)
        {
            genes.Add(($"GENE{i}", new List<double> { i, i + 1, i + 2 }));
        }

        var clusters = TranscriptomeAnalyzer.ClusterGenesByExpression(
            genes, numClusters: 4).ToList();

        Assert.That(clusters, Has.Count.EqualTo(4));
        Assert.That(clusters.Sum(c => c.Genes.Count), Is.EqualTo(20));
    }

    [Test]
    public void ClusterGenesByExpression_EmptyInput_ReturnsEmpty()
    {
        var clusters = TranscriptomeAnalyzer.ClusterGenesByExpression(
            new List<(string, IReadOnlyList<double>)>()).ToList();

        Assert.That(clusters, Is.Empty);
    }

    #endregion

    #region Quality Control Tests

    [Test]
    public void CalculateQualityMetrics_TypicalValues_CalculatesCorrectly()
    {
        double totalReads = 100_000_000;
        double mappedReads = 85_000_000;
        double exonicReads = 70_000_000;
        double rrnaReads = 5_000_000;
        var geneCounts = Enumerable.Range(0, 20000).Select(i => i < 15000 ? (double)i : 0);

        var metrics = TranscriptomeAnalyzer.CalculateQualityMetrics(
            totalReads, mappedReads, exonicReads, rrnaReads, geneCounts);

        Assert.That(metrics.MappingRate, Is.EqualTo(0.85).Within(0.001));
        Assert.That(metrics.ExonicRate, Is.EqualTo(70.0 / 85).Within(0.001));
        Assert.That(metrics.RRNARate, Is.EqualTo(5.0 / 85).Within(0.001));
        Assert.That(metrics.DetectedGenes, Is.EqualTo(14999)); // i=1 to 14999 have counts > 0
    }

    [Test]
    public void CalculateQualityMetrics_ZeroReads_HandlesGracefully()
    {
        var metrics = TranscriptomeAnalyzer.CalculateQualityMetrics(
            0, 0, 0, 0, new List<double>());

        Assert.That(metrics.MappingRate, Is.EqualTo(0));
        Assert.That(metrics.DetectedGenes, Is.EqualTo(0));
    }

    [Test]
    public void PerformPCA_MultipleSamples_ReturnsCoordinates()
    {
        var samples = new List<(string, IReadOnlyList<double>)>
        {
            ("S1", Enumerable.Range(0, 100).Select(i => (double)i).ToList()),
            ("S2", Enumerable.Range(0, 100).Select(i => (double)(i + 10)).ToList()),
            ("S3", Enumerable.Range(0, 100).Select(i => (double)(i * 2)).ToList())
        };

        var pca = TranscriptomeAnalyzer.PerformPCA(samples).ToList();

        Assert.That(pca, Has.Count.EqualTo(3));
        Assert.That(pca.All(p => p.SampleId != null), Is.True);
    }

    [Test]
    public void PerformPCA_SingleSample_ReturnsZeroCoordinates()
    {
        var samples = new List<(string, IReadOnlyList<double>)>
        {
            ("S1", new List<double> { 1, 2, 3 })
        };

        var pca = TranscriptomeAnalyzer.PerformPCA(samples).ToList();

        Assert.That(pca, Has.Count.EqualTo(1));
        Assert.That(pca[0].PC1, Is.EqualTo(0));
        Assert.That(pca[0].PC2, Is.EqualTo(0));
    }

    #endregion

    #region Edge Cases Tests

    [Test]
    public void CalculateTPM_SingleGene_FullMillion()
    {
        var genes = new List<(string, double, int)>
        {
            ("GENE1", 1000, 1000)
        };

        var results = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();

        Assert.That(results[0].TPM, Is.EqualTo(1_000_000));
    }

    [Test]
    public void AnalyzeDifferentialExpression_EmptyGroup_Skipped()
    {
        var data = new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>
        {
            ("GENE1", new List<double>(), new List<double> { 1, 2, 3 })
        };

        var results = TranscriptomeAnalyzer.AnalyzeDifferentialExpression(data).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void BuildCoExpressionNetwork_SingleGene_NoEdges()
    {
        var genes = new List<(string, IReadOnlyList<double>)>
        {
            ("GENE1", new List<double> { 1, 2, 3 })
        };

        var edges = TranscriptomeAnalyzer.BuildCoExpressionNetwork(genes).ToList();

        Assert.That(edges, Is.Empty);
    }

    [Test]
    public void PerformOverRepresentationAnalysis_EmptyDE_ReturnsEmpty()
    {
        var deGenes = new HashSet<string>();
        var pathways = new List<(string, string, IReadOnlySet<string>)>
        {
            ("P1", "Pathway1", new HashSet<string> { "A", "B" })
        };

        var results = TranscriptomeAnalyzer.PerformOverRepresentationAnalysis(
            deGenes, pathways, 1000).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void FindSkippedExonEvents_ZeroReads_Skipped()
    {
        var exonData = new List<(string, int, int, double, double)>
        {
            ("GENE1", 100, 200, 0, 0)
        };

        var events = TranscriptomeAnalyzer.FindSkippedExonEvents(exonData).ToList();

        Assert.That(events, Is.Empty);
    }

    #endregion
}
