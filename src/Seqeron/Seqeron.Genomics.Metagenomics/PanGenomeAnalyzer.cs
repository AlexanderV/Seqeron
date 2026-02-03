using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Metagenomics;

/// <summary>
/// Provides algorithms for pan-genome analysis.
/// </summary>
public static class PanGenomeAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a pan-genome analysis result.
    /// </summary>
    public readonly record struct PanGenomeResult(
        IReadOnlyList<string> CoreGenes,
        IReadOnlyList<string> AccessoryGenes,
        IReadOnlyList<string> UniqueGenes,
        IReadOnlyDictionary<string, IReadOnlyList<string>> GenomeToGenes,
        PanGenomeStatistics Statistics);

    /// <summary>
    /// Statistics about the pan-genome.
    /// </summary>
    public readonly record struct PanGenomeStatistics(
        int TotalGenomes,
        int TotalGenes,
        int CoreGeneCount,
        int AccessoryGeneCount,
        int UniqueGeneCount,
        double CoreFraction,
        double GenomeFluidity,
        PanGenomeType Type);

    /// <summary>
    /// Type of pan-genome (open vs closed).
    /// </summary>
    public enum PanGenomeType
    {
        Open,     // New genomes add new genes
        Closed    // Gene content is largely conserved
    }

    /// <summary>
    /// Represents a gene cluster (ortholog group).
    /// </summary>
    public readonly record struct GeneCluster(
        string ClusterId,
        IReadOnlyList<string> GeneIds,
        IReadOnlyList<string> GenomeIds,
        int GenomeCount,
        double AverageIdentity,
        string ConsensusSequence);

    /// <summary>
    /// Represents a gene presence/absence matrix entry.
    /// </summary>
    public readonly record struct GenePresenceRow(
        string GenomeId,
        IReadOnlyDictionary<string, bool> GenePresence,
        int TotalGenes,
        int PresentGenes);

    /// <summary>
    /// Result of heaps law fitting for pan-genome size prediction.
    /// </summary>
    public readonly record struct HeapsLawFit(
        double K,
        double Gamma,
        double RSquared,
        Func<int, double> PredictPanGenomeSize);

    #endregion

    #region Pan-Genome Construction

    /// <summary>
    /// Constructs a pan-genome from multiple genomes.
    /// </summary>
    public static PanGenomeResult ConstructPanGenome(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double identityThreshold = 0.9,
        double coreFraction = 0.99)
    {
        if (genomes == null || genomes.Count == 0)
        {
            return new PanGenomeResult(
                new List<string>(),
                new List<string>(),
                new List<string>(),
                new Dictionary<string, IReadOnlyList<string>>(),
                new PanGenomeStatistics(0, 0, 0, 0, 0, 0, 0, PanGenomeType.Closed));
        }

        // Cluster genes into ortholog groups
        var clusters = ClusterGenes(genomes, identityThreshold).ToList();

        int totalGenomes = genomes.Count;
        int coreThreshold = (int)(totalGenomes * coreFraction);

        var coreGenes = new List<string>();
        var accessoryGenes = new List<string>();
        var uniqueGenes = new List<string>();

        foreach (var cluster in clusters)
        {
            if (cluster.GenomeCount >= coreThreshold)
            {
                coreGenes.Add(cluster.ClusterId);
            }
            else if (cluster.GenomeCount == 1)
            {
                uniqueGenes.Add(cluster.ClusterId);
            }
            else
            {
                accessoryGenes.Add(cluster.ClusterId);
            }
        }

        var genomeToGenes = genomes.ToDictionary(
            g => g.Key,
            g => (IReadOnlyList<string>)g.Value.Select(gene => gene.GeneId).ToList());

        int totalGenes = clusters.Count;
        double coreFrac = totalGenes > 0 ? (double)coreGenes.Count / totalGenes : 0;

        // Calculate genome fluidity
        double fluidity = CalculateGenomeFluidity(genomes, clusters);

        // Determine if pan-genome is open or closed
        var type = DeterminePanGenomeType(genomes, clusters);

        var stats = new PanGenomeStatistics(
            TotalGenomes: totalGenomes,
            TotalGenes: totalGenes,
            CoreGeneCount: coreGenes.Count,
            AccessoryGeneCount: accessoryGenes.Count,
            UniqueGeneCount: uniqueGenes.Count,
            CoreFraction: coreFrac,
            GenomeFluidity: fluidity,
            Type: type);

        return new PanGenomeResult(coreGenes, accessoryGenes, uniqueGenes, genomeToGenes, stats);
    }

    /// <summary>
    /// Clusters genes into ortholog groups based on sequence similarity.
    /// </summary>
    public static IEnumerable<GeneCluster> ClusterGenes(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double identityThreshold = 0.9)
    {
        var allGenes = new List<(string GenomeId, string GeneId, string Sequence)>();

        foreach (var (genomeId, genes) in genomes)
        {
            foreach (var (geneId, sequence) in genes)
            {
                allGenes.Add((genomeId, geneId, sequence));
            }
        }

        if (allGenes.Count == 0)
            yield break;

        var assigned = new HashSet<int>();
        int clusterId = 1;

        for (int i = 0; i < allGenes.Count; i++)
        {
            if (assigned.Contains(i))
                continue;

            var clusterMembers = new List<int> { i };
            assigned.Add(i);

            for (int j = i + 1; j < allGenes.Count; j++)
            {
                if (assigned.Contains(j))
                    continue;

                double identity = CalculateSequenceIdentity(allGenes[i].Sequence, allGenes[j].Sequence);
                if (identity >= identityThreshold)
                {
                    clusterMembers.Add(j);
                    assigned.Add(j);
                }
            }

            var geneIds = clusterMembers.Select(m => allGenes[m].GeneId).ToList();
            var genomeIds = clusterMembers.Select(m => allGenes[m].GenomeId).Distinct().ToList();

            // Calculate average identity within cluster
            double avgIdentity = 1.0;
            if (clusterMembers.Count > 1)
            {
                var identities = new List<double>();
                for (int a = 0; a < clusterMembers.Count; a++)
                {
                    for (int b = a + 1; b < clusterMembers.Count; b++)
                    {
                        identities.Add(CalculateSequenceIdentity(
                            allGenes[clusterMembers[a]].Sequence,
                            allGenes[clusterMembers[b]].Sequence));
                    }
                }
                avgIdentity = identities.Count > 0 ? identities.Average() : 1.0;
            }

            // Get consensus (just use first sequence for simplicity)
            string consensus = allGenes[clusterMembers[0]].Sequence;

            yield return new GeneCluster(
                ClusterId: $"cluster_{clusterId++}",
                GeneIds: geneIds,
                GenomeIds: genomeIds,
                GenomeCount: genomeIds.Count,
                AverageIdentity: avgIdentity,
                ConsensusSequence: consensus);
        }
    }

    private static double CalculateSequenceIdentity(string seq1, string seq2)
    {
        if (string.IsNullOrEmpty(seq1) || string.IsNullOrEmpty(seq2))
            return 0;

        int minLen = Math.Min(seq1.Length, seq2.Length);
        int maxLen = Math.Max(seq1.Length, seq2.Length);

        if (maxLen == 0)
            return 1;

        // Quick k-mer based similarity
        int k = 7;
        if (minLen < k)
            return seq1 == seq2 ? 1 : 0;

        var kmers1 = new HashSet<string>();
        var kmers2 = new HashSet<string>();

        for (int i = 0; i <= seq1.Length - k; i++)
            kmers1.Add(seq1.Substring(i, k));

        for (int i = 0; i <= seq2.Length - k; i++)
            kmers2.Add(seq2.Substring(i, k));

        int shared = kmers1.Intersect(kmers2).Count();
        int total = kmers1.Union(kmers2).Count();

        return total > 0 ? (double)shared / total : 0;
    }

    #endregion

    #region Gene Presence/Absence Matrix

    /// <summary>
    /// Creates a gene presence/absence matrix.
    /// </summary>
    public static IEnumerable<GenePresenceRow> CreatePresenceAbsenceMatrix(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IEnumerable<GeneCluster> clusters)
    {
        var clusterList = clusters.ToList();
        var clusterToGenes = new Dictionary<string, HashSet<string>>();

        foreach (var cluster in clusterList)
        {
            clusterToGenes[cluster.ClusterId] = new HashSet<string>(cluster.GeneIds);
        }

        foreach (var (genomeId, genes) in genomes)
        {
            var geneSet = new HashSet<string>(genes.Select(g => g.GeneId));
            var presence = new Dictionary<string, bool>();

            foreach (var cluster in clusterList)
            {
                bool present = cluster.GeneIds.Any(g => geneSet.Contains(g));
                presence[cluster.ClusterId] = present;
            }

            int presentCount = presence.Values.Count(v => v);

            yield return new GenePresenceRow(
                GenomeId: genomeId,
                GenePresence: presence,
                TotalGenes: clusterList.Count,
                PresentGenes: presentCount);
        }
    }

    #endregion

    #region Pan-Genome Statistics

    /// <summary>
    /// Calculates genome fluidity (dissimilarity between genome pairs).
    /// </summary>
    private static double CalculateGenomeFluidity(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IEnumerable<GeneCluster> clusters)
    {
        var genomeList = genomes.Keys.ToList();
        if (genomeList.Count < 2)
            return 0;

        var clusterList = clusters.ToList();

        // Build genome to cluster map
        var genomeToClusterIds = new Dictionary<string, HashSet<string>>();
        foreach (var genomeId in genomeList)
        {
            genomeToClusterIds[genomeId] = new HashSet<string>();
        }

        foreach (var cluster in clusterList)
        {
            foreach (var genomeId in cluster.GenomeIds)
            {
                if (genomeToClusterIds.ContainsKey(genomeId))
                {
                    genomeToClusterIds[genomeId].Add(cluster.ClusterId);
                }
            }
        }

        double totalFluidity = 0;
        int pairCount = 0;

        for (int i = 0; i < genomeList.Count; i++)
        {
            for (int j = i + 1; j < genomeList.Count; j++)
            {
                var set1 = genomeToClusterIds[genomeList[i]];
                var set2 = genomeToClusterIds[genomeList[j]];

                int unique = set1.Except(set2).Count() + set2.Except(set1).Count();
                int total = set1.Count + set2.Count;

                if (total > 0)
                {
                    totalFluidity += (double)unique / total;
                    pairCount++;
                }
            }
        }

        return pairCount > 0 ? totalFluidity / pairCount : 0;
    }

    /// <summary>
    /// Determines if pan-genome is open or closed based on gene accumulation.
    /// </summary>
    private static PanGenomeType DeterminePanGenomeType(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IEnumerable<GeneCluster> clusters)
    {
        var clusterList = clusters.ToList();

        // Simple heuristic: if unique genes > 10% of total, likely open
        int uniqueGenes = clusterList.Count(c => c.GenomeCount == 1);
        int totalGenes = clusterList.Count;

        double uniqueFraction = totalGenes > 0 ? (double)uniqueGenes / totalGenes : 0;

        return uniqueFraction > 0.1 ? PanGenomeType.Open : PanGenomeType.Closed;
    }

    /// <summary>
    /// Fits Heaps' law to predict pan-genome size growth.
    /// </summary>
    public static HeapsLawFit FitHeapsLaw(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double identityThreshold = 0.9,
        int permutations = 10)
    {
        if (genomes.Count < 3)
        {
            return new HeapsLawFit(0, 0, 0, n => 0);
        }

        var genomeList = genomes.Keys.ToList();
        var dataPoints = new List<(int N, double PanSize)>();
        var random = new Random(42);

        for (int perm = 0; perm < permutations; perm++)
        {
            var shuffled = genomeList.OrderBy(_ => random.Next()).ToList();
            var accumulatedGenes = new HashSet<string>();

            for (int i = 0; i < shuffled.Count; i++)
            {
                var genomeGenes = genomes[shuffled[i]];

                foreach (var (geneId, sequence) in genomeGenes)
                {
                    bool isNew = true;
                    foreach (var existingGene in accumulatedGenes.ToList())
                    {
                        // Simplified: just use gene ID matching
                        if (geneId == existingGene)
                        {
                            isNew = false;
                            break;
                        }
                    }

                    if (isNew)
                        accumulatedGenes.Add(geneId);
                }

                dataPoints.Add((i + 1, accumulatedGenes.Count));
            }
        }

        // Average pan-genome size at each N
        var avgSizes = dataPoints
            .GroupBy(p => p.N)
            .ToDictionary(g => g.Key, g => g.Average(p => p.PanSize));

        // Fit Heaps' law: P(N) = K * N^gamma
        // Using log-linear regression: log(P) = log(K) + gamma * log(N)
        var logN = avgSizes.Keys.Select(n => Math.Log(n)).ToList();
        var logP = avgSizes.Values.Select(p => Math.Log(Math.Max(p, 1))).ToList();

        var (slope, intercept, rSquared) = LinearRegression(logN, logP);

        double gamma = slope;
        double k = Math.Exp(intercept);

        Func<int, double> predictor = n => k * Math.Pow(n, gamma);

        return new HeapsLawFit(k, gamma, rSquared, predictor);
    }

    private static (double Slope, double Intercept, double RSquared) LinearRegression(
        List<double> x, List<double> y)
    {
        if (x.Count != y.Count || x.Count < 2)
            return (0, 0, 0);

        double n = x.Count;
        double sumX = x.Sum();
        double sumY = y.Sum();
        double sumXY = x.Zip(y, (a, b) => a * b).Sum();
        double sumX2 = x.Sum(a => a * a);

        double slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
        double intercept = (sumY - slope * sumX) / n;

        // Calculate R-squared
        double meanY = sumY / n;
        double ssTotal = y.Sum(yi => (yi - meanY) * (yi - meanY));
        double ssResidual = x.Zip(y, (xi, yi) => yi - (slope * xi + intercept))
            .Sum(residual => residual * residual);

        double rSquared = ssTotal > 0 ? 1 - ssResidual / ssTotal : 0;

        return (slope, intercept, rSquared);
    }

    #endregion

    #region Core Genome Analysis

    /// <summary>
    /// Extracts core genes with optional filtering.
    /// </summary>
    public static IEnumerable<GeneCluster> GetCoreGeneClusters(
        IEnumerable<GeneCluster> clusters,
        int totalGenomes,
        double threshold = 0.99)
    {
        int minGenomes = (int)(totalGenomes * threshold);
        return clusters.Where(c => c.GenomeCount >= minGenomes);
    }

    /// <summary>
    /// Calculates the core genome alignment (concatenated core genes).
    /// </summary>
    public static string CreateCoreGenomeAlignment(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IEnumerable<GeneCluster> coreClusters,
        string genomeId)
    {
        if (!genomes.TryGetValue(genomeId, out var genes))
            return "";

        var geneDict = genes.ToDictionary(g => g.GeneId, g => g.Sequence);
        var sb = new StringBuilder();

        foreach (var cluster in coreClusters)
        {
            var matchingGene = cluster.GeneIds.FirstOrDefault(g => geneDict.ContainsKey(g));
            if (matchingGene != null)
            {
                sb.Append(geneDict[matchingGene]);
            }
        }

        return sb.ToString();
    }

    #endregion

    #region Accessory Genome Analysis

    /// <summary>
    /// Analyzes accessory genome patterns.
    /// </summary>
    public static IEnumerable<(string ClusterId, IReadOnlyList<string> GenomesWithGene, double Frequency)>
        AnalyzeAccessoryGenes(
            IEnumerable<GeneCluster> clusters,
            int totalGenomes)
    {
        return clusters
            .Where(c => c.GenomeCount > 1 && c.GenomeCount < totalGenomes)
            .Select(c => (
                c.ClusterId,
                c.GenomeIds,
                (double)c.GenomeCount / totalGenomes));
    }

    /// <summary>
    /// Finds genes unique to specific genomes.
    /// </summary>
    public static IEnumerable<(string GenomeId, IReadOnlyList<string> UniqueGeneIds)>
        FindGenomeSpecificGenes(
            IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
            IEnumerable<GeneCluster> clusters)
    {
        var uniqueClusters = clusters.Where(c => c.GenomeCount == 1).ToList();

        var genomeToUnique = new Dictionary<string, List<string>>();
        foreach (var genomeId in genomes.Keys)
        {
            genomeToUnique[genomeId] = new List<string>();
        }

        foreach (var cluster in uniqueClusters)
        {
            var genomeId = cluster.GenomeIds.First();
            if (genomeToUnique.ContainsKey(genomeId))
            {
                genomeToUnique[genomeId].Add(cluster.ClusterId);
            }
        }

        return genomeToUnique
            .Where(kv => kv.Value.Count > 0)
            .Select(kv => (kv.Key, (IReadOnlyList<string>)kv.Value));
    }

    #endregion

    #region Phylogenetic Marker Selection

    /// <summary>
    /// Selects informative markers from core genes for phylogenetic analysis.
    /// </summary>
    public static IEnumerable<GeneCluster> SelectPhylogeneticMarkers(
        IEnumerable<GeneCluster> coreClusters,
        int maxMarkers = 100,
        double minIdentity = 0.7,
        double maxIdentity = 0.99)
    {
        return coreClusters
            .Where(c => c.AverageIdentity >= minIdentity && c.AverageIdentity <= maxIdentity)
            .OrderByDescending(c => c.ConsensusSequence.Length)
            .ThenBy(c => c.AverageIdentity)
            .Take(maxMarkers);
    }

    #endregion
}
