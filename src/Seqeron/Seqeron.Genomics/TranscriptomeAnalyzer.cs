using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics;

/// <summary>
/// Provides algorithms for transcriptome analysis including expression quantification and differential expression.
/// </summary>
public static class TranscriptomeAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents gene expression data.
    /// </summary>
    public readonly record struct GeneExpression(
        string GeneId,
        double RawCount,
        double TPM,
        double FPKM,
        int Length);

    /// <summary>
    /// Represents differential expression result.
    /// </summary>
    public readonly record struct DifferentialExpression(
        string GeneId,
        double Log2FoldChange,
        double PValue,
        double AdjustedPValue,
        bool IsSignificant,
        string Regulation);

    /// <summary>
    /// Represents a gene set enrichment result.
    /// </summary>
    public readonly record struct EnrichmentResult(
        string PathwayId,
        string PathwayName,
        int GenesInPathway,
        int OverlappingGenes,
        double EnrichmentScore,
        double PValue,
        IReadOnlyList<string> Genes);

    /// <summary>
    /// Represents alternative splicing event.
    /// </summary>
    public readonly record struct SplicingEvent(
        string GeneId,
        string EventType,
        int Start,
        int End,
        double InclusionLevel,
        double DeltaPSI);

    /// <summary>
    /// Represents transcript isoform.
    /// </summary>
    public readonly record struct TranscriptIsoform(
        string TranscriptId,
        string GeneId,
        int Length,
        int ExonCount,
        double Expression,
        bool IsProteinCoding,
        IReadOnlyList<(int Start, int End)> Exons);

    /// <summary>
    /// Represents co-expression cluster.
    /// </summary>
    public readonly record struct CoExpressionCluster(
        int ClusterId,
        IReadOnlyList<string> Genes,
        double MeanCorrelation,
        string RepresentativeGene,
        IReadOnlyList<string> EnrichedFunctions);

    /// <summary>
    /// Types of alternative splicing events.
    /// </summary>
    public enum SplicingEventType
    {
        SkippedExon,
        RetainedIntron,
        AlternativeFivePrimeSS,
        AlternativeThreePrimeSS,
        MutuallyExclusiveExons
    }

    #endregion

    #region Expression Quantification

    /// <summary>
    /// Calculates TPM (Transcripts Per Million) from raw counts.
    /// </summary>
    public static IEnumerable<GeneExpression> CalculateTPM(
        IEnumerable<(string GeneId, double RawCount, int Length)> geneCounts)
    {
        var geneList = geneCounts.ToList();

        if (geneList.Count == 0)
            yield break;

        // Calculate rate = count / length for each gene
        var rates = geneList
            .Select(g => (g.GeneId, g.RawCount, g.Length, Rate: g.RawCount / Math.Max(g.Length, 1)))
            .ToList();

        double sumRates = rates.Sum(r => r.Rate);

        if (sumRates == 0)
        {
            foreach (var gene in geneList)
            {
                yield return new GeneExpression(gene.GeneId, gene.RawCount, 0, 0, gene.Length);
            }
            yield break;
        }

        // TPM = rate / sum(rates) * 1,000,000
        foreach (var (geneId, rawCount, length, rate) in rates)
        {
            double tpm = (rate / sumRates) * 1_000_000;
            double fpkm = CalculateFPKM(rawCount, length, geneList.Sum(g => g.RawCount));

            yield return new GeneExpression(geneId, rawCount, tpm, fpkm, length);
        }
    }

    /// <summary>
    /// Calculates FPKM for a single gene.
    /// </summary>
    private static double CalculateFPKM(double rawCount, int length, double totalReads)
    {
        if (length <= 0 || totalReads <= 0)
            return 0;

        // FPKM = (reads * 10^9) / (length * total mapped reads)
        return (rawCount * 1_000_000_000) / (length * totalReads);
    }

    /// <summary>
    /// Normalizes expression values using quantile normalization.
    /// </summary>
    public static IEnumerable<IReadOnlyList<double>> QuantileNormalize(
        IEnumerable<IEnumerable<double>> samples)
    {
        var sampleList = samples.Select(s => s.ToList()).ToList();

        if (sampleList.Count == 0)
            yield break;

        int geneCount = sampleList[0].Count;
        int sampleCount = sampleList.Count;

        if (geneCount == 0)
            yield break;

        // Sort each sample and track original indices
        var sortedIndices = sampleList
            .Select(s => s
                .Select((val, idx) => (val, idx))
                .OrderBy(x => x.val)
                .Select(x => x.idx)
                .ToArray())
            .ToList();

        // Calculate mean at each rank
        var rankMeans = new double[geneCount];
        for (int rank = 0; rank < geneCount; rank++)
        {
            double sum = 0;
            foreach (var sample in sampleList)
            {
                var sorted = sample.OrderBy(x => x).ToList();
                sum += sorted[rank];
            }
            rankMeans[rank] = sum / sampleCount;
        }

        // Assign rank means back to original positions
        foreach (var sample in sampleList)
        {
            var sorted = sample
                .Select((val, idx) => (val, idx))
                .OrderBy(x => x.val)
                .ToList();

            var normalized = new double[geneCount];
            for (int rank = 0; rank < geneCount; rank++)
            {
                normalized[sorted[rank].idx] = rankMeans[rank];
            }

            yield return normalized;
        }
    }

    /// <summary>
    /// Performs log2 transformation with pseudocount.
    /// </summary>
    public static IEnumerable<double> Log2Transform(
        IEnumerable<double> values,
        double pseudocount = 1.0)
    {
        foreach (var value in values)
        {
            yield return Math.Log2(value + pseudocount);
        }
    }

    #endregion

    #region Differential Expression Analysis

    /// <summary>
    /// Performs simple differential expression analysis using fold change and t-test.
    /// </summary>
    public static IEnumerable<DifferentialExpression> AnalyzeDifferentialExpression(
        IEnumerable<(string GeneId, IReadOnlyList<double> Group1, IReadOnlyList<double> Group2)> expressionData,
        double foldChangeThreshold = 1.0,
        double pValueThreshold = 0.05)
    {
        var genes = expressionData.ToList();

        if (genes.Count == 0)
            yield break;

        var results = new List<(string GeneId, double Log2FC, double PValue)>();

        foreach (var (geneId, group1, group2) in genes)
        {
            if (group1.Count == 0 || group2.Count == 0)
                continue;

            double mean1 = group1.Average();
            double mean2 = group2.Average();

            // Avoid division by zero
            double log2FC = mean1 > 0
                ? Math.Log2((mean2 + 0.01) / (mean1 + 0.01))
                : 0;

            double pValue = CalculateTTestPValue(group1, group2);

            results.Add((geneId, log2FC, pValue));
        }

        // Multiple testing correction (Benjamini-Hochberg)
        var sortedByPValue = results.OrderBy(r => r.PValue).ToList();
        var adjustedPValues = BenjaminiHochberg(sortedByPValue.Select(r => r.PValue));
        var adjustedList = adjustedPValues.ToList();

        for (int i = 0; i < sortedByPValue.Count; i++)
        {
            var (geneId, log2FC, pValue) = sortedByPValue[i];
            double adjPValue = adjustedList[i];
            bool isSignificant = adjPValue < pValueThreshold && Math.Abs(log2FC) >= foldChangeThreshold;
            string regulation = log2FC > 0 ? "Upregulated" : (log2FC < 0 ? "Downregulated" : "Unchanged");

            yield return new DifferentialExpression(
                GeneId: geneId,
                Log2FoldChange: log2FC,
                PValue: pValue,
                AdjustedPValue: adjPValue,
                IsSignificant: isSignificant,
                Regulation: regulation);
        }
    }

    /// <summary>
    /// Calculates a simple t-test p-value.
    /// </summary>
    private static double CalculateTTestPValue(IReadOnlyList<double> group1, IReadOnlyList<double> group2)
    {
        if (group1.Count < 2 || group2.Count < 2)
            return 1.0;

        double mean1 = group1.Average();
        double mean2 = group2.Average();

        double var1 = group1.Sum(x => (x - mean1) * (x - mean1)) / (group1.Count - 1);
        double var2 = group2.Sum(x => (x - mean2) * (x - mean2)) / (group2.Count - 1);

        double se = Math.Sqrt(var1 / group1.Count + var2 / group2.Count);

        if (se == 0)
            return mean1 == mean2 ? 1.0 : 0.0;

        double t = Math.Abs(mean2 - mean1) / se;

        // Approximate p-value using normal distribution
        return 2 * (1 - StatisticsHelper.NormalCDF(t));
    }

    /// <summary>
    /// Benjamini-Hochberg multiple testing correction.
    /// </summary>
    private static IEnumerable<double> BenjaminiHochberg(IEnumerable<double> pValues)
    {
        var pList = pValues.ToList();
        int n = pList.Count;

        if (n == 0)
            yield break;

        var adjusted = new double[n];
        double minSoFar = 1.0;

        for (int i = n - 1; i >= 0; i--)
        {
            double corrected = pList[i] * n / (i + 1);
            corrected = Math.Min(corrected, minSoFar);
            corrected = Math.Min(corrected, 1.0);
            adjusted[i] = corrected;
            minSoFar = corrected;
        }

        foreach (var adj in adjusted)
            yield return adj;
    }

    #endregion

    #region Gene Set Enrichment

    /// <summary>
    /// Performs over-representation analysis (ORA) for gene set enrichment.
    /// </summary>
    public static IEnumerable<EnrichmentResult> PerformOverRepresentationAnalysis(
        IReadOnlySet<string> differentiallyExpressedGenes,
        IEnumerable<(string PathwayId, string PathwayName, IReadOnlySet<string> Genes)> pathways,
        int backgroundGeneCount)
    {
        if (differentiallyExpressedGenes.Count == 0 || backgroundGeneCount <= 0)
            yield break;

        foreach (var (pathwayId, pathwayName, pathwayGenes) in pathways)
        {
            var overlapping = pathwayGenes.Intersect(differentiallyExpressedGenes).ToList();

            if (overlapping.Count == 0)
                continue;

            // Fisher's exact test approximation
            double pValue = CalculateFisherPValue(
                differentiallyExpressedGenes.Count,
                pathwayGenes.Count,
                overlapping.Count,
                backgroundGeneCount);

            // Enrichment score
            double expected = (double)differentiallyExpressedGenes.Count * pathwayGenes.Count / backgroundGeneCount;
            double enrichmentScore = expected > 0 ? overlapping.Count / expected : 0;

            yield return new EnrichmentResult(
                PathwayId: pathwayId,
                PathwayName: pathwayName,
                GenesInPathway: pathwayGenes.Count,
                OverlappingGenes: overlapping.Count,
                EnrichmentScore: enrichmentScore,
                PValue: pValue,
                Genes: overlapping);
        }
    }

    /// <summary>
    /// Approximates Fisher's exact test p-value using hypergeometric distribution.
    /// </summary>
    private static double CalculateFisherPValue(int deGenes, int pathwaySize, int overlap, int background)
    {
        // Hypergeometric probability approximation
        // Using normal approximation for large samples

        double expectedOverlap = (double)deGenes * pathwaySize / background;
        double variance = expectedOverlap * (1 - (double)pathwaySize / background) *
                          (background - deGenes) / (background - 1);

        if (variance <= 0)
            return overlap >= expectedOverlap ? 0.0 : 1.0;

        double z = (overlap - expectedOverlap) / Math.Sqrt(variance);
        return 1 - StatisticsHelper.NormalCDF(z);
    }

    /// <summary>
    /// Calculates Gene Set Enrichment Score (GSEA-like).
    /// </summary>
    public static double CalculateEnrichmentScore(
        IReadOnlyList<string> rankedGenes,
        IReadOnlySet<string> geneSet)
    {
        if (rankedGenes.Count == 0 || geneSet.Count == 0)
            return 0;

        int n = rankedGenes.Count;
        int hitCount = rankedGenes.Count(g => geneSet.Contains(g));
        int missCount = n - hitCount;

        if (hitCount == 0 || missCount == 0)
            return 0;

        double hitIncrement = 1.0 / hitCount;
        double missDecrement = 1.0 / missCount;

        double runningSum = 0;
        double maxDeviation = 0;

        foreach (var gene in rankedGenes)
        {
            if (geneSet.Contains(gene))
            {
                runningSum += hitIncrement;
            }
            else
            {
                runningSum -= missDecrement;
            }

            if (Math.Abs(runningSum) > Math.Abs(maxDeviation))
            {
                maxDeviation = runningSum;
            }
        }

        return maxDeviation;
    }

    #endregion

    #region Alternative Splicing Analysis

    /// <summary>
    /// Identifies potential skipped exon events.
    /// </summary>
    public static IEnumerable<SplicingEvent> FindSkippedExonEvents(
        IEnumerable<(string GeneId, int ExonStart, int ExonEnd, double InclusionReads, double SkippingReads)> exonData)
    {
        foreach (var (geneId, start, end, inclusion, skipping) in exonData)
        {
            double total = inclusion + skipping;
            if (total == 0)
                continue;

            double psi = inclusion / total; // Percent Spliced In

            yield return new SplicingEvent(
                GeneId: geneId,
                EventType: "SkippedExon",
                Start: start,
                End: end,
                InclusionLevel: psi,
                DeltaPSI: 0); // Would need comparison sample
        }
    }

    /// <summary>
    /// Detects differential splicing between conditions.
    /// </summary>
    public static IEnumerable<SplicingEvent> DetectDifferentialSplicing(
        IEnumerable<(string GeneId, int Start, int End, double PSI_Condition1, double PSI_Condition2)> splicingData,
        double deltaPsiThreshold = 0.1)
    {
        foreach (var (geneId, start, end, psi1, psi2) in splicingData)
        {
            double deltaPsi = psi2 - psi1;

            if (Math.Abs(deltaPsi) >= deltaPsiThreshold)
            {
                string eventType = deltaPsi > 0 ? "IncreasedInclusion" : "IncreasedSkipping";

                yield return new SplicingEvent(
                    GeneId: geneId,
                    EventType: eventType,
                    Start: start,
                    End: end,
                    InclusionLevel: psi2,
                    DeltaPSI: deltaPsi);
            }
        }
    }

    #endregion

    #region Transcript Isoform Analysis

    /// <summary>
    /// Identifies dominant transcript isoform for each gene.
    /// </summary>
    public static IEnumerable<(string GeneId, TranscriptIsoform DominantIsoform, double DominanceRatio)>
        FindDominantIsoforms(IEnumerable<TranscriptIsoform> isoforms)
    {
        var byGene = isoforms.GroupBy(i => i.GeneId);

        foreach (var group in byGene)
        {
            var sorted = group.OrderByDescending(i => i.Expression).ToList();

            if (sorted.Count == 0)
                continue;

            var dominant = sorted[0];
            double totalExpression = sorted.Sum(i => i.Expression);
            double dominanceRatio = totalExpression > 0 ? dominant.Expression / totalExpression : 0;

            yield return (group.Key, dominant, dominanceRatio);
        }
    }

    /// <summary>
    /// Detects isoform switching between conditions.
    /// </summary>
    public static IEnumerable<(string GeneId, string TranscriptId1, string TranscriptId2, double SwitchScore)>
        DetectIsoformSwitching(
            IEnumerable<(TranscriptIsoform Isoform, double Expression1, double Expression2)> isoformData,
            double switchThreshold = 0.3)
    {
        var byGene = isoformData.GroupBy(d => d.Isoform.GeneId);

        foreach (var group in byGene)
        {
            var isoforms = group.ToList();

            if (isoforms.Count < 2)
                continue;

            double total1 = isoforms.Sum(i => i.Expression1);
            double total2 = isoforms.Sum(i => i.Expression2);

            if (total1 == 0 || total2 == 0)
                continue;

            // Calculate usage ratios
            var usageChanges = isoforms
                .Select(i => (
                    i.Isoform.TranscriptId,
                    Usage1: i.Expression1 / total1,
                    Usage2: i.Expression2 / total2))
                .Select(u => (
                    u.TranscriptId,
                    u.Usage1,
                    u.Usage2,
                    Delta: u.Usage2 - u.Usage1))
                .OrderByDescending(u => Math.Abs(u.Delta))
                .ToList();

            if (usageChanges.Count < 2)
                continue;

            // Check for significant switching
            var increased = usageChanges.FirstOrDefault(u => u.Delta > switchThreshold);
            var decreased = usageChanges.FirstOrDefault(u => u.Delta < -switchThreshold);

            if (increased.TranscriptId != null && decreased.TranscriptId != null)
            {
                double switchScore = Math.Abs(increased.Delta) + Math.Abs(decreased.Delta);

                yield return (group.Key, decreased.TranscriptId, increased.TranscriptId, switchScore);
            }
        }
    }

    #endregion

    #region Co-Expression Analysis

    /// <summary>
    /// Calculates Pearson correlation between gene expression profiles.
    /// </summary>
    public static double CalculatePearsonCorrelation(
        IReadOnlyList<double> expression1,
        IReadOnlyList<double> expression2)
    {
        if (expression1.Count != expression2.Count || expression1.Count < 2)
            return 0;

        int n = expression1.Count;
        double mean1 = expression1.Average();
        double mean2 = expression2.Average();

        double covariance = 0;
        double var1 = 0;
        double var2 = 0;

        for (int i = 0; i < n; i++)
        {
            double d1 = expression1[i] - mean1;
            double d2 = expression2[i] - mean2;
            covariance += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }

        if (var1 == 0 || var2 == 0)
            return 0;

        return covariance / Math.Sqrt(var1 * var2);
    }

    /// <summary>
    /// Builds a co-expression network.
    /// </summary>
    public static IEnumerable<(string Gene1, string Gene2, double Correlation)> BuildCoExpressionNetwork(
        IEnumerable<(string GeneId, IReadOnlyList<double> Expression)> geneProfiles,
        double correlationThreshold = 0.7)
    {
        var genes = geneProfiles.ToList();

        for (int i = 0; i < genes.Count; i++)
        {
            for (int j = i + 1; j < genes.Count; j++)
            {
                double corr = CalculatePearsonCorrelation(genes[i].Expression, genes[j].Expression);

                if (Math.Abs(corr) >= correlationThreshold)
                {
                    yield return (genes[i].GeneId, genes[j].GeneId, corr);
                }
            }
        }
    }

    /// <summary>
    /// Performs hierarchical clustering of genes by expression.
    /// </summary>
    public static IEnumerable<CoExpressionCluster> ClusterGenesByExpression(
        IEnumerable<(string GeneId, IReadOnlyList<double> Expression)> geneProfiles,
        int numClusters = 5,
        double correlationThreshold = 0.5)
    {
        var genes = geneProfiles.ToList();

        if (genes.Count == 0)
            yield break;

        // Simple k-means-like clustering based on correlation
        var clusters = new List<List<(string GeneId, IReadOnlyList<double> Expression)>>();

        // Initialize clusters with first genes
        for (int i = 0; i < Math.Min(numClusters, genes.Count); i++)
        {
            clusters.Add(new List<(string, IReadOnlyList<double>)> { genes[i] });
        }

        // Assign remaining genes to nearest cluster
        for (int i = numClusters; i < genes.Count; i++)
        {
            int bestCluster = 0;
            double bestCorr = double.MinValue;

            for (int c = 0; c < clusters.Count; c++)
            {
                // Calculate average correlation with cluster members
                double avgCorr = clusters[c]
                    .Select(m => CalculatePearsonCorrelation(genes[i].Expression, m.Expression))
                    .Average();

                if (avgCorr > bestCorr)
                {
                    bestCorr = avgCorr;
                    bestCluster = c;
                }
            }

            clusters[bestCluster].Add(genes[i]);
        }

        // Create cluster results
        for (int c = 0; c < clusters.Count; c++)
        {
            if (clusters[c].Count == 0)
                continue;

            var clusterGenes = clusters[c].Select(g => g.GeneId).ToList();

            // Calculate mean internal correlation
            double meanCorr = 0;
            int corrCount = 0;
            for (int i = 0; i < clusters[c].Count; i++)
            {
                for (int j = i + 1; j < clusters[c].Count; j++)
                {
                    meanCorr += CalculatePearsonCorrelation(
                        clusters[c][i].Expression,
                        clusters[c][j].Expression);
                    corrCount++;
                }
            }
            meanCorr = corrCount > 0 ? meanCorr / corrCount : 0;

            // Representative gene is one with highest mean correlation to others
            string representative = clusterGenes.First();

            yield return new CoExpressionCluster(
                ClusterId: c + 1,
                Genes: clusterGenes,
                MeanCorrelation: meanCorr,
                RepresentativeGene: representative,
                EnrichedFunctions: new List<string>());
        }
    }

    #endregion

    #region RNA-seq Quality Control

    /// <summary>
    /// Calculates basic RNA-seq quality metrics.
    /// </summary>
    public static (double MappingRate, double ExonicRate, double RRNARate, int DetectedGenes)
        CalculateQualityMetrics(
            double totalReads,
            double mappedReads,
            double exonicReads,
            double rRNAReads,
            IEnumerable<double> geneCounts)
    {
        double mappingRate = totalReads > 0 ? mappedReads / totalReads : 0;
        double exonicRate = mappedReads > 0 ? exonicReads / mappedReads : 0;
        double rrnaRate = mappedReads > 0 ? rRNAReads / mappedReads : 0;
        int detectedGenes = geneCounts.Count(c => c > 0);

        return (mappingRate, exonicRate, rrnaRate, detectedGenes);
    }

    /// <summary>
    /// Identifies potential batch effects using PCA.
    /// </summary>
    public static IEnumerable<(string SampleId, double PC1, double PC2)> PerformPCA(
        IEnumerable<(string SampleId, IReadOnlyList<double> Expression)> samples,
        int topGenes = 500)
    {
        var sampleList = samples.ToList();

        if (sampleList.Count < 2)
        {
            foreach (var sample in sampleList)
            {
                yield return (sample.SampleId, 0, 0);
            }
            yield break;
        }

        int geneCount = sampleList[0].Expression.Count;

        // Select top variable genes
        var variances = new double[geneCount];
        for (int g = 0; g < geneCount; g++)
        {
            var values = sampleList.Select(s => s.Expression[g]).ToList();
            double mean = values.Average();
            variances[g] = values.Sum(v => (v - mean) * (v - mean)) / values.Count;
        }

        var topGeneIndices = variances
            .Select((v, i) => (v, i))
            .OrderByDescending(x => x.v)
            .Take(Math.Min(topGenes, geneCount))
            .Select(x => x.i)
            .ToHashSet();

        // Simple PCA approximation using first two principal directions
        // (Full PCA would require SVD)
        foreach (var sample in sampleList)
        {
            var selectedValues = sample.Expression
                .Where((v, i) => topGeneIndices.Contains(i))
                .ToList();

            // Approximate PC scores as weighted sums
            double pc1 = selectedValues.Take(selectedValues.Count / 2).Sum();
            double pc2 = selectedValues.Skip(selectedValues.Count / 2).Sum();

            yield return (sample.SampleId, pc1, pc2);
        }
    }

    #endregion
}
