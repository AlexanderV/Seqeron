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
    /// Result of fitting Heaps' law to the pan-genome new-gene-discovery curve
    /// (Tettelin et al., 2008; micropan <c>heaps()</c>). The fitted model is
    /// <c>n(N) = Intercept · N^(-Alpha)</c>, the number of <em>new</em> gene clusters
    /// contributed by the N-th genome.
    /// </summary>
    /// <param name="Intercept">
    /// The Heaps' law intercept K (micropan "Intercept" parameter), constrained to
    /// [0, 10000] per the micropan <c>heaps()</c> reference implementation.
    /// </param>
    /// <param name="Alpha">
    /// The Heaps' law decay exponent α (micropan "alpha" parameter), constrained to
    /// [0, 2]. Per Tettelin et al. (2008) the pan-genome is OPEN when α &lt; 1 and
    /// CLOSED when α &gt; 1.
    /// </param>
    /// <param name="IsOpen">
    /// <c>true</c> when α &lt; 1 (open pan-genome), otherwise <c>false</c>
    /// (Tettelin et al., 2008).
    /// </param>
    /// <param name="PredictNewGenes">
    /// Predictor n(N) = Intercept · N^(-Alpha): the expected number of new gene
    /// clusters added by the N-th genome.
    /// </param>
    public readonly record struct HeapsLawFit(
        double Intercept,
        double Alpha,
        bool IsOpen,
        Func<int, double> PredictNewGenes);

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

    // CD-HIT default sequence-identity clustering threshold (-c option, default 0.9):
    // "sequence identity threshold, default 0.9" (Li & Godzik 2006; CD-HIT User's Guide).
    private const double DefaultIdentityThreshold = 0.9;

    /// <summary>
    /// Clusters genes into ortholog (homolog) groups using the CD-HIT greedy incremental
    /// clustering model (Li &amp; Godzik, 2006): sequences are sorted from longest to
    /// shortest; the longest sequence becomes the representative of the first cluster, and
    /// each remaining sequence joins the first existing representative whose global
    /// sequence identity meets <paramref name="identityThreshold"/>, otherwise it becomes
    /// the representative of a new cluster. Global sequence identity is the number of
    /// identical residues in the (ungapped) alignment divided by the length of the shorter
    /// sequence (CD-HIT "-G 1" default).
    /// </summary>
    /// <param name="genomes">Genome id → list of (gene id, sequence) entries.</param>
    /// <param name="identityThreshold">
    /// Global sequence-identity cutoff in [0,1]; a gene joins a cluster only when its
    /// identity to that cluster's representative is &gt;= this value (CD-HIT -c, default 0.9).
    /// </param>
    public static IEnumerable<GeneCluster> ClusterGenes(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double identityThreshold = DefaultIdentityThreshold)
    {
        if (genomes == null)
            yield break;

        var allGenes = new List<(string GenomeId, string GeneId, string Sequence)>();
        foreach (var (genomeId, genes) in genomes)
        {
            if (genes == null)
                continue;
            foreach (var (geneId, sequence) in genes)
                allGenes.Add((genomeId, geneId, sequence));
        }

        if (allGenes.Count == 0)
            yield break;

        // CD-HIT: "sorts the input sequences from long to short, and processes them
        // sequentially from the longest to the shortest." Ties are broken by original
        // order via a stable sort so the result is deterministic.
        var order = Enumerable.Range(0, allGenes.Count)
            .OrderByDescending(i => allGenes[i].Sequence?.Length ?? 0)
            .ToList();

        // Each cluster keeps its representative's index plus its accumulated members.
        var representatives = new List<int>();
        var members = new List<List<int>>();

        foreach (int idx in order)
        {
            string seq = allGenes[idx].Sequence ?? string.Empty;

            int joined = -1;
            // "each query sequence ... is compared to the representative sequences found
            // before it, and is classified as redundant ... if it is similar to one of the
            // existing representative sequences." First match wins (greedy).
            for (int c = 0; c < representatives.Count; c++)
            {
                double identity = CalculateSequenceIdentity(seq, allGenes[representatives[c]].Sequence ?? string.Empty);
                if (identity >= identityThreshold)
                {
                    joined = c;
                    break;
                }
            }

            if (joined >= 0)
            {
                members[joined].Add(idx);
            }
            else
            {
                representatives.Add(idx);
                members.Add(new List<int> { idx });
            }
        }

        int clusterId = 1;
        foreach (var clusterMembers in members)
        {
            var geneIds = clusterMembers.Select(m => allGenes[m].GeneId).ToList();
            var genomeIds = clusterMembers.Select(m => allGenes[m].GenomeId).Distinct().ToList();

            // Average pairwise global identity within the cluster (1.0 for singletons,
            // which are 100% identical to themselves).
            double avgIdentity = 1.0;
            if (clusterMembers.Count > 1)
            {
                double sum = 0;
                int pairs = 0;
                for (int a = 0; a < clusterMembers.Count; a++)
                {
                    for (int b = a + 1; b < clusterMembers.Count; b++)
                    {
                        sum += CalculateSequenceIdentity(
                            allGenes[clusterMembers[a]].Sequence ?? string.Empty,
                            allGenes[clusterMembers[b]].Sequence ?? string.Empty);
                        pairs++;
                    }
                }
                avgIdentity = pairs > 0 ? sum / pairs : 1.0;
            }

            // The cluster representative (longest member, per CD-HIT) is its consensus.
            string consensus = allGenes[clusterMembers[0]].Sequence ?? string.Empty;

            yield return new GeneCluster(
                ClusterId: $"cluster_{clusterId++}",
                GeneIds: geneIds,
                GenomeIds: genomeIds,
                GenomeCount: genomeIds.Count,
                AverageIdentity: avgIdentity,
                ConsensusSequence: consensus);
        }
    }

    /// <summary>
    /// Computes CD-HIT global sequence identity between two sequences: the number of
    /// identical residues in the ungapped positional alignment divided by the length of
    /// the shorter sequence (Li &amp; Godzik, 2006; CD-HIT "-G 1" default). Two empty
    /// sequences are defined as identical (1.0); one empty and one non-empty are 0.0.
    /// </summary>
    private static double CalculateSequenceIdentity(string seq1, string seq2)
    {
        seq1 ??= string.Empty;
        seq2 ??= string.Empty;

        int len1 = seq1.Length;
        int len2 = seq2.Length;

        // Two empty sequences: identical by convention (0 differences over 0 positions).
        if (len1 == 0 && len2 == 0)
            return 1.0;

        // One empty, one not: no shared residues -> 0% identity.
        if (len1 == 0 || len2 == 0)
            return 0.0;

        int shorter = Math.Min(len1, len2);

        // Ungapped alignment: compare residues position-by-position over the shared prefix.
        // Positions beyond the shorter length count as non-identical (no residue to match).
        int identical = 0;
        for (int i = 0; i < shorter; i++)
        {
            if (seq1[i] == seq2[i])
                identical++;
        }

        // Global identity denominator is the full length of the shorter sequence.
        return (double)identical / shorter;
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

    // Heaps' law openness criterion: the number of NEW gene clusters added by the
    // k-th genome follows a power law n_new(k) = K * k^(-alpha). The pan-genome is
    // OPEN when alpha < 1 (new genes keep accumulating without bound) and CLOSED when
    // alpha > 1. Per Tettelin et al. (2008) Curr Opin Microbiol 11:472, and the
    // micropan heaps() reference implementation: "If alpha<1.0 the pan-genome is open,
    // if alpha>1.0 it is closed."
    private const double HeapsOpennessThreshold = 1.0;

    // The decay exponent is only meaningful once several genomes have accumulated; the
    // new-gene curve is degenerate below this many genomes (Tettelin 2008; micropan).
    private const int MinGenomesForOpennessFit = 3;

    /// <summary>
    /// Determines whether the pan-genome is open or closed using the Heaps' law decay
    /// exponent of newly observed gene clusters per added genome (Tettelin et al., 2008).
    /// Open when the decay exponent alpha &lt; 1, closed when alpha &gt; 1.
    /// </summary>
    private static PanGenomeType DeterminePanGenomeType(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IEnumerable<GeneCluster> clusters)
    {
        var clusterList = clusters.ToList();

        // Below the minimum count the decay exponent cannot be estimated; the conservative
        // default is Closed (no evidence of unbounded growth).
        if (genomes.Count < MinGenomesForOpennessFit || clusterList.Count == 0)
            return PanGenomeType.Closed;

        double alpha = EstimateHeapsDecayExponent(genomes, clusterList);

        // alpha == threshold is the boundary; treat the non-open boundary as Closed.
        return alpha < HeapsOpennessThreshold ? PanGenomeType.Open : PanGenomeType.Closed;
    }

    /// <summary>
    /// Estimates the Heaps' law decay exponent alpha of the new-gene-cluster curve
    /// n_new(k) = K * k^(-alpha) by log-log least-squares regression over the cumulative
    /// gene-cluster accumulation as genomes are added in dictionary order.
    /// </summary>
    private static double EstimateHeapsDecayExponent(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        IReadOnlyList<GeneCluster> clusterList)
    {
        // Map each genome to the set of cluster IDs it contains.
        var genomeToClusters = new Dictionary<string, HashSet<string>>();
        foreach (var genomeId in genomes.Keys)
            genomeToClusters[genomeId] = new HashSet<string>();

        foreach (var cluster in clusterList)
        {
            foreach (var genomeId in cluster.GenomeIds)
            {
                if (genomeToClusters.TryGetValue(genomeId, out var set))
                    set.Add(cluster.ClusterId);
            }
        }

        // Accumulate genomes in order; record the count of NEW clusters contributed at
        // each step k = 2, 3, ... (the first genome has no "new" baseline to compare).
        var accumulated = new HashSet<string>();
        var logK = new List<double>();
        var logNew = new List<double>();

        int k = 0;
        foreach (var genomeId in genomes.Keys)
        {
            k++;
            int before = accumulated.Count;
            accumulated.UnionWith(genomeToClusters[genomeId]);
            int newClusters = accumulated.Count - before;

            // Only positions k >= 2 are part of the decay curve; log requires positive
            // new-cluster counts, so zero-novelty steps are mapped to the minimum (1) to
            // keep the curve defined (cumulative curve flattening drives alpha upward).
            if (k >= 2)
            {
                logK.Add(Math.Log(k));
                logNew.Add(Math.Log(Math.Max(newClusters, 1)));
            }
        }

        if (logK.Count < 2)
            return HeapsOpennessThreshold; // not enough points to fit -> boundary (Closed)

        // log(n_new) = log(K) - alpha * log(k); slope of the regression is -alpha.
        var (slope, _, _) = LinearRegression(logK, logNew);
        return -slope;
    }

    // micropan heaps() classification rule (Tettelin et al., 2008): the pan-genome is
    // OPEN when the decay exponent alpha < 1 and CLOSED when alpha > 1.
    private const double HeapsAlphaOpenThreshold = 1.0;

    // Lower/upper bounds on the fitted parameters, copied verbatim from the micropan
    // heaps() reference implementation:
    //   optim(p0, objectFun, ..., lower = c(0, 0), upper = c(10000, 2))
    // i.e. Intercept K in [0, 10000] and decay exponent alpha in [0, 2].
    private const double HeapsInterceptLowerBound = 0.0;
    private const double HeapsInterceptUpperBound = 10000.0;
    private const double HeapsAlphaLowerBound = 0.0;
    private const double HeapsAlphaUpperBound = 2.0;

    // Heaps' law is fit to new-gene counts at genome indices N = 2, 3, ... ; the first
    // genome has no predecessor against which "new" can be defined (micropan: x = 2:ng).
    private const int HeapsFirstFittedIndex = 2;

    // micropan default number of random genome-order permutations (heaps() n.perm
    // default = 100; "The default value of 100 is certainly a minimum.").
    private const int DefaultHeapsPermutations = 100;

    // Deterministic RNG seed so the permutation-averaged fit is reproducible.
    private const int HeapsPermutationSeed = 42;

    /// <summary>
    /// Fits Heaps' law to the pan-genome new-gene-discovery curve from a gene
    /// presence/absence matrix, following the micropan <c>heaps()</c> reference
    /// implementation (Tettelin et al., 2008). For each of <paramref name="permutations"/>
    /// random genome orderings, the number of new gene clusters contributed by the N-th
    /// genome is recorded; the model <c>n(N) = K · N^(-alpha)</c> is then fit to the pooled
    /// (N, new-gene-count) points for N = 2..G by minimizing
    /// <c>J = sqrt(Σ (y − K·x^(−alpha))²) / |x|</c> over K ∈ [0,10000], alpha ∈ [0,2].
    /// </summary>
    /// <param name="matrix">
    /// Gene presence/absence rows (one per genome), e.g. from
    /// <see cref="CreatePresenceAbsenceMatrix"/>. A cluster is present in a genome when its
    /// entry in <see cref="GenePresenceRow.GenePresence"/> is <c>true</c>.
    /// </param>
    /// <param name="permutations">Number of random genome-order permutations to pool over
    /// (micropan n.perm; default 100).</param>
    /// <returns>The fitted Intercept (K), decay exponent (alpha), open/closed flag, and a
    /// predictor for new genes at the N-th genome.</returns>
    public static HeapsLawFit FitHeapsLaw(
        IEnumerable<GenePresenceRow> matrix,
        int permutations = DefaultHeapsPermutations)
    {
        var rows = matrix?.ToList() ?? new List<GenePresenceRow>();
        return FitHeapsLawFromOrderings(BuildPermutedNewGeneCurves(rows, permutations));
    }

    /// <summary>
    /// Convenience overload that clusters <paramref name="genomes"/> into ortholog groups,
    /// builds the presence/absence matrix, and fits Heaps' law (micropan <c>heaps()</c>;
    /// Tettelin et al., 2008). See <see cref="FitHeapsLaw(IEnumerable{GenePresenceRow}, int)"/>.
    /// </summary>
    public static HeapsLawFit FitHeapsLaw(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes,
        double identityThreshold = DefaultIdentityThreshold,
        int permutations = DefaultHeapsPermutations)
    {
        if (genomes == null || genomes.Count == 0)
            return EmptyHeapsFit();

        var clusters = ClusterGenes(genomes, identityThreshold).ToList();
        var matrix = CreatePresenceAbsenceMatrix(genomes, clusters);
        return FitHeapsLaw(matrix, permutations);
    }

    /// <summary>
    /// Builds, for each random genome ordering, the per-step count of newly observed gene
    /// clusters. A cluster is "new" at ordering position i (1-based, i &gt;= 2) when it is
    /// present at position i but absent in all positions 1..i-1 — micropan's
    /// <c>rowSums((cm == 1)[2:ng,] &amp; (cm == 0)[1:(ng-1),])</c> over the binary matrix.
    /// </summary>
    private static List<int[]> BuildPermutedNewGeneCurves(
        IReadOnlyList<GenePresenceRow> rows, int permutations)
    {
        var curves = new List<int[]>();
        int genomeCount = rows.Count;
        if (genomeCount < HeapsFirstFittedIndex)
            return curves;

        // Stable cluster ordering across all rows (the matrix columns).
        var clusterIds = rows
            .SelectMany(r => r.GenePresence.Keys)
            .Distinct()
            .ToList();
        if (clusterIds.Count == 0)
            return curves;

        // Binary presence vector per genome (micropan: pan.matrix[pan.matrix > 0] <- 1).
        var presence = rows
            .Select(r => clusterIds
                .Select(c => r.GenePresence.TryGetValue(c, out bool p) && p)
                .ToArray())
            .ToArray();

        int permCount = Math.Max(1, permutations);
        var random = new Random(HeapsPermutationSeed);

        for (int perm = 0; perm < permCount; perm++)
        {
            int[] order = Enumerable.Range(0, genomeCount).ToArray();
            // First permutation uses the natural input order so a single-permutation,
            // fixed-order fit is exactly reproducible; the rest are Fisher-Yates shuffles.
            if (perm > 0)
            {
                for (int i = genomeCount - 1; i > 0; i--)
                {
                    int j = random.Next(i + 1);
                    (order[i], order[j]) = (order[j], order[i]);
                }
            }

            var seen = new bool[clusterIds.Count];
            // Seed "seen" with the first genome; new counts are recorded from position 2.
            for (int c = 0; c < clusterIds.Count; c++)
                seen[c] = presence[order[0]][c];

            var newCounts = new int[genomeCount - 1]; // for positions 2..G
            for (int pos = 1; pos < genomeCount; pos++)
            {
                int added = 0;
                var pv = presence[order[pos]];
                for (int c = 0; c < clusterIds.Count; c++)
                {
                    if (pv[c] && !seen[c])
                    {
                        seen[c] = true;
                        added++;
                    }
                }
                newCounts[pos - 1] = added;
            }
            curves.Add(newCounts);
        }

        return curves;
    }

    /// <summary>
    /// Fits n(N) = K · N^(-alpha) to the pooled (N, new-gene-count) points and returns the
    /// micropan-style result. Pooling and the objective follow micropan <c>heaps()</c>:
    /// x = rep(2:ng, n.perm), y = stacked new-gene counts.
    /// </summary>
    private static HeapsLawFit FitHeapsLawFromOrderings(List<int[]> curves)
    {
        if (curves.Count == 0)
            return EmptyHeapsFit();

        var x = new List<double>();
        var y = new List<double>();
        foreach (var curve in curves)
        {
            for (int k = 0; k < curve.Length; k++)
            {
                x.Add(HeapsFirstFittedIndex + k); // genome index N = 2, 3, ...
                y.Add(curve[k]);
            }
        }

        if (x.Count == 0)
            return EmptyHeapsFit();

        // micropan start values: p0 = c(mean(y at x==2), 1).
        double meanAtFirst = y.Where((_, i) => x[i] == HeapsFirstFittedIndex).DefaultIfEmpty(0).Average();
        double k0 = Math.Clamp(meanAtFirst, HeapsInterceptLowerBound, HeapsInterceptUpperBound);

        var (intercept, alpha) = MinimizeHeapsObjective(x, y, k0, HeapsAlphaOpenThreshold);

        bool isOpen = alpha < HeapsAlphaOpenThreshold;
        Func<int, double> predictor = n => intercept * Math.Pow(n, -alpha);
        return new HeapsLawFit(intercept, alpha, isOpen, predictor);
    }

    private static HeapsLawFit EmptyHeapsFit()
    {
        // Below 2 genomes (or no clusters) the new-gene curve is undefined: report a
        // degenerate, closed fit with a zero predictor rather than inventing parameters.
        return new HeapsLawFit(0, 0, IsOpen: false, n => 0);
    }

    /// <summary>
    /// Heaps' law least-squares objective from micropan heaps():
    /// <c>J(K, alpha) = sqrt(Σ (y − K·x^(−alpha))²) / |x|</c>.
    /// </summary>
    private static double HeapsObjective(IReadOnlyList<double> x, IReadOnlyList<double> y,
        double k, double alpha)
    {
        double ss = 0;
        for (int i = 0; i < x.Count; i++)
        {
            double yHat = k * Math.Pow(x[i], -alpha);
            double r = y[i] - yHat;
            ss += r * r;
        }
        return Math.Sqrt(ss) / x.Count;
    }

    /// <summary>
    /// Bounded minimization of the Heaps' law objective over K ∈ [0,10000], alpha ∈ [0,2]
    /// (the micropan optim() bounds). Deterministic coordinate descent with geometric
    /// step refinement from the micropan start point; for data lying exactly on a power
    /// curve within the bounds it recovers the analytic (K, alpha) to machine precision.
    /// </summary>
    private static (double K, double Alpha) MinimizeHeapsObjective(
        IReadOnlyList<double> x, IReadOnlyList<double> y, double k0, double alpha0)
    {
        double k = Math.Clamp(k0, HeapsInterceptLowerBound, HeapsInterceptUpperBound);
        double alpha = Math.Clamp(alpha0, HeapsAlphaLowerBound, HeapsAlphaUpperBound);
        double best = HeapsObjective(x, y, k, alpha);

        double kStep = Math.Max(k, 1.0);
        double alphaStep = HeapsAlphaUpperBound;
        const double minStep = 1e-12;
        const int maxIterations = 200_000;

        int iterations = 0;
        while ((kStep > minStep || alphaStep > minStep) && iterations < maxIterations)
        {
            bool improved = false;
            foreach (double dk in new[] { kStep, -kStep })
            {
                double nk = Math.Clamp(k + dk, HeapsInterceptLowerBound, HeapsInterceptUpperBound);
                double j = HeapsObjective(x, y, nk, alpha);
                if (j < best) { best = j; k = nk; improved = true; }
            }
            foreach (double da in new[] { alphaStep, -alphaStep })
            {
                double na = Math.Clamp(alpha + da, HeapsAlphaLowerBound, HeapsAlphaUpperBound);
                double j = HeapsObjective(x, y, k, na);
                if (j < best) { best = j; alpha = na; improved = true; }
            }
            if (!improved)
            {
                kStep *= 0.5;
                alphaStep *= 0.5;
            }
            iterations++;
        }

        return (k, alpha);
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
