using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Tumor Phylogeny Reconstruction (ONCO-PHYLO-001)

    /// <summary>
    /// Cancer cell fraction (CCF) of the synthetic root (normal / germline) clone in every sample: it is, by
    /// definition, present in 100% of cells. Source: Popic et al. (2015), <i>Genome Biology</i> 16:91 — the
    /// lineage tree is a spanning tree rooted at the population that contains all observed clones.
    /// </summary>
    private const double RootCcf = 1.0;

    /// <summary>
    /// Default noise margin ε for the lineage-precedence (Eq. 2) and sum-rule (Eq. 5) inequalities. The cited
    /// sources relax both rules by a configurable ε (Popic et al. 2015 ϵ; Zheng et al. 2022 ε₁=0.1, ε₂=0.2). Because
    /// this unit consumes already-clustered CCF point estimates (clustering and its noise model are ONCO-CCF-001),
    /// the default is the strict ε = 0; callers may pass a positive tolerance to reproduce the source defaults.
    /// </summary>
    public const double DefaultPhylogenyTolerance = 0.0;

    /// <summary>
    /// One CCF cluster (a candidate clone/subclone) used as input to phylogeny reconstruction. Each cluster carries
    /// its cancer cell fraction in each sequenced sample. CCF clustering itself is out of scope (ONCO-CCF-001).
    /// </summary>
    /// <param name="Id">Caller-assigned cluster identifier (e.g. a subclone label). Used for deterministic tie-breaking.</param>
    /// <param name="CcfPerSample">Cancer cell fraction in each sample, each value in [0, 1]; all clusters must share length.</param>
    public readonly record struct CcfCluster(int Id, IReadOnlyList<double> CcfPerSample);

    /// <summary>A single parent → child edge of the reconstructed clonal tree.</summary>
    /// <param name="ParentId">Id of the ancestral cluster, or <see cref="ClonalPhylogeny.RootId"/> for the normal root.</param>
    /// <param name="ChildId">Id of the descendant cluster.</param>
    public readonly record struct ClonalEdge(int ParentId, int ChildId);

    /// <summary>
    /// The reconstructed rooted clonal tree: a synthetic normal root (<see cref="RootId"/>, CCF = 1 in every sample)
    /// plus one node per input CCF cluster, connected by parent→child <see cref="Edges"/>.
    /// </summary>
    /// <param name="RootId">Identifier of the synthetic normal/germline root node.</param>
    /// <param name="Clusters">Input clusters keyed by id, in input order.</param>
    /// <param name="Edges">Tree edges (parent→child), one per non-root cluster.</param>
    /// <param name="SampleCount">Number of samples per cluster.</param>
    public readonly record struct ClonalPhylogeny(
        int RootId,
        IReadOnlyList<CcfCluster> Clusters,
        IReadOnlyList<ClonalEdge> Edges,
        int SampleCount)
    {
        /// <summary>Returns the parent id of <paramref name="clusterId"/>, or null if it is the root or absent.</summary>
        public int? ParentOf(int clusterId)
        {
            foreach (ClonalEdge e in Edges)
            {
                if (e.ChildId == clusterId)
                {
                    return e.ParentId;
                }
            }

            return null;
        }

        /// <summary>Returns the ids of the direct children of <paramref name="clusterId"/>, in input order.</summary>
        public IReadOnlyList<int> ChildrenOf(int clusterId)
        {
            var children = new List<int>();
            foreach (ClonalEdge e in Edges)
            {
                if (e.ParentId == clusterId)
                {
                    children.Add(e.ChildId);
                }
            }

            return children;
        }
    }

    /// <summary>
    /// Reconstructs a rooted clonal (tumor) phylogeny from per-sample CCF clusters, applying the two lineage
    /// constraints from the multi-sample perfect-phylogeny model:
    /// <list type="number">
    /// <item><description><b>Lineage precedence (ancestor ≥ descendant), Eq. 2:</b> an edge u→v is admissible only if,
    /// for every sample i, <c>u.CCF[i] ≥ v.CCF[i] − ε</c> and (presence) <c>u.CCF[i] = 0 ⇒ v.CCF[i] = 0</c>. Source:
    /// Popic et al. (2015), <i>Genome Biology</i> 16:91, Eq. 2; Zheng et al. (2022) PICTograph, <i>Bioinformatics</i>
    /// 38(15):3677–3683 — "the CCF of any mutation cannot exceed the CCF of its ancestor".</description></item>
    /// <item><description><b>Sum rule, Eq. 5:</b> for every node u and every sample i, the children CCFs may not exceed
    /// the parent: <c>Σ_children v.CCF[i] ≤ u.CCF[i] + ε</c>. Source: Popic et al. (2015) Eq. 5; Zheng et al. (2022) —
    /// "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants".</description></item>
    /// </list>
    /// The constraints leave a set of valid trees; to return a single deterministic tree this method attaches each
    /// cluster (processed in descending order of total CCF) to its <i>deepest valid ancestor</i> — the admissible
    /// parent with the smallest total CCF whose remaining per-sample sum-rule budget still admits the child — with
    /// ties broken by ascending cluster id (Evidence Assumption 1).
    /// </summary>
    /// <param name="clusters">CCF clusters to place; each cluster's <see cref="CcfCluster.CcfPerSample"/> must have the same length.</param>
    /// <param name="tolerance">Noise margin ε for both inequalities; default <see cref="DefaultPhylogenyTolerance"/> (0).</param>
    /// <returns>The reconstructed <see cref="ClonalPhylogeny"/> rooted at a synthetic normal node.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="clusters"/> or any cluster's CCF list is null.</exception>
    /// <exception cref="ArgumentException">CCF lists differ in length, are empty, or contain NaN / out-of-[0,1] values; or two clusters share an id.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tolerance"/> is negative or NaN.</exception>
    public static ClonalPhylogeny ReconstructPhylogeny(
        IReadOnlyList<CcfCluster> clusters,
        double tolerance = DefaultPhylogenyTolerance)
    {
        ArgumentNullException.ThrowIfNull(clusters);
        if (double.IsNaN(tolerance) || tolerance < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tolerance), tolerance, "Phylogeny tolerance ε must be a non-negative number.");
        }

        int rootId = RootIdFor(clusters);
        if (clusters.Count == 0)
        {
            // Empty cohort: tree is the root alone, no clusters, no edges. One synthetic sample of CCF 1.
            return new ClonalPhylogeny(rootId, Array.Empty<CcfCluster>(), Array.Empty<ClonalEdge>(), 1);
        }

        int sampleCount = ValidateAndGetSampleCount(clusters);

        // Synthetic root: present in 100% of cells in every sample, so it is a valid ancestor of every cluster and
        // its sum-rule budget is the full RootCcf. The root participates in placement as node id = rootId.
        double[] rootCcf = new double[sampleCount];
        Array.Fill(rootCcf, RootCcf);

        // Per-node remaining sum-rule budget = node CCF minus the sum of CCFs of children already attached (per sample).
        var remainingBudget = new Dictionary<int, double[]>(clusters.Count + 1);
        remainingBudget[rootId] = (double[])rootCcf.Clone();
        var ccfById = new Dictionary<int, double[]>(clusters.Count + 1) { [rootId] = rootCcf };
        foreach (CcfCluster c in clusters)
        {
            double[] ccf = c.CcfPerSample.ToArray();
            ccfById[c.Id] = ccf;
            remainingBudget[c.Id] = (double[])ccf.Clone();
        }

        // Process clusters most-clonal-first (descending total CCF) so that ancestors are placed before descendants;
        // ties broken by ascending id for determinism (Evidence Assumption 1).
        var ordered = clusters
            .OrderByDescending(c => TotalCcf(ccfById[c.Id]))
            .ThenBy(c => c.Id)
            .ToList();

        var edges = new List<ClonalEdge>(clusters.Count);
        foreach (CcfCluster child in ordered)
        {
            double[] childCcf = ccfById[child.Id];

            // Candidate parents: the root plus every already-placed cluster. Choose the deepest valid ancestor =
            // the candidate with the smallest total CCF that (a) satisfies lineage precedence and (b) still has
            // per-sample budget for this child. Ties broken by ascending id.
            int bestParent = rootId;
            double bestParentTotal = double.PositiveInfinity;
            bool found = false;
            foreach (int candidateId in EnumerateCandidates(rootId, edges))
            {
                if (candidateId == child.Id)
                {
                    continue;
                }

                double[] parentCcf = ccfById[candidateId];
                if (!SatisfiesLineagePrecedence(parentCcf, childCcf, tolerance))
                {
                    continue;
                }

                if (!FitsSumRule(remainingBudget[candidateId], childCcf, tolerance))
                {
                    continue;
                }

                double candidateTotal = TotalCcf(parentCcf);
                if (!found
                    || candidateTotal < bestParentTotal
                    || (candidateTotal == bestParentTotal && candidateId < bestParent))
                {
                    bestParent = candidateId;
                    bestParentTotal = candidateTotal;
                    found = true;
                }
            }

            // The root always satisfies lineage precedence (CCF=1 ≥ any child) and, per the sum rule, can only be
            // exhausted if children already consumed its full budget; in that degenerate case we still attach to the
            // root because every cluster must have a parent (Popic 2015: spanning tree). 'found' is therefore the
            // generic path; the explicit root fallback preserves the spanning-tree invariant.
            int parentId = found ? bestParent : rootId;
            edges.Add(new ClonalEdge(parentId, child.Id));

            // Debit the chosen parent's per-sample budget by the child's CCF.
            double[] budget = remainingBudget[parentId];
            for (int i = 0; i < sampleCount; i++)
            {
                budget[i] -= childCcf[i];
            }
        }

        var clusterList = clusters.ToArray();
        return new ClonalPhylogeny(rootId, clusterList, edges, sampleCount);
    }

    /// <summary>
    /// Returns the ids of clusters on the <b>trunk</b> of the phylogeny — the clonal mutations shared by every
    /// tumor cell. The trunk is the path from the root down to (but excluding) the first branch point: a maximal
    /// chain of single-child nodes starting at the root's unique child. Source: Popic et al. (2015) — the trunk
    /// holds the mutations of the common predecessor present across all samples.
    /// </summary>
    /// <param name="phylogeny">A phylogeny produced by <see cref="ReconstructPhylogeny"/>.</param>
    /// <returns>Trunk cluster ids, ordered from the root downward; empty if the tree has no clusters.</returns>
    public static IReadOnlyList<int> IdentifyTrunkMutations(ClonalPhylogeny phylogeny)
    {
        var trunk = new List<int>();
        if (phylogeny.Clusters is null || phylogeny.Clusters.Count == 0)
        {
            return trunk;
        }

        // Walk down from the root while each node has exactly one child (no branching yet) and the root itself has
        // exactly one child. The first node with ≠1 children is the branch point; nodes below it are subclonal.
        int current = phylogeny.RootId;
        while (true)
        {
            IReadOnlyList<int> children = phylogeny.ChildrenOf(current);
            if (children.Count != 1)
            {
                break;
            }

            int only = children[0];
            trunk.Add(only);
            current = only;
        }

        return trunk;
    }

    /// <summary>
    /// Returns the ids of <b>branch</b> (subclonal) clusters — every input cluster that is not on the trunk.
    /// Source: Popic et al. (2015) — mutations off the common-predecessor trunk are subclonal lineage branches.
    /// </summary>
    /// <param name="phylogeny">A phylogeny produced by <see cref="ReconstructPhylogeny"/>.</param>
    /// <returns>Branch cluster ids in input order.</returns>
    public static IReadOnlyList<int> IdentifyBranchMutations(ClonalPhylogeny phylogeny)
    {
        var branches = new List<int>();
        if (phylogeny.Clusters is null || phylogeny.Clusters.Count == 0)
        {
            return branches;
        }

        var trunk = new HashSet<int>(IdentifyTrunkMutations(phylogeny));
        foreach (CcfCluster c in phylogeny.Clusters)
        {
            if (!trunk.Contains(c.Id))
            {
                branches.Add(c.Id);
            }
        }

        return branches;
    }

    /// <summary>
    /// Lineage precedence (Popic 2015 Eq. 2): for every sample i, <c>parent ≥ child − ε</c> and a parent absent in
    /// a sample (CCF = 0) cannot have a child present there (presence pattern, constraint (1)).
    /// </summary>
    private static bool SatisfiesLineagePrecedence(double[] parentCcf, double[] childCcf, double tolerance)
    {
        for (int i = 0; i < parentCcf.Length; i++)
        {
            if (parentCcf[i] < childCcf[i] - tolerance)
            {
                return false;
            }

            // Presence: if the parent is absent in sample i, the child must also be absent there.
            if (parentCcf[i] <= 0.0 && childCcf[i] > 0.0)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>
    /// Sum rule (Popic 2015 Eq. 5): the child's CCF fits a parent only if, in every sample, the parent's remaining
    /// budget (parent CCF minus already-attached children) is ≥ child CCF − ε.
    /// </summary>
    private static bool FitsSumRule(double[] remainingParentBudget, double[] childCcf, double tolerance)
    {
        for (int i = 0; i < remainingParentBudget.Length; i++)
        {
            if (remainingParentBudget[i] < childCcf[i] - tolerance)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>Sum of a cluster's CCF over all samples (proxy for clonal dominance / processing order).</summary>
    private static double TotalCcf(double[] ccf)
    {
        double total = 0.0;
        foreach (double v in ccf)
        {
            total += v;
        }

        return total;
    }

    /// <summary>Candidate parents = the root plus every cluster already attached to the tree (in attachment order).</summary>
    private static IEnumerable<int> EnumerateCandidates(int rootId, List<ClonalEdge> edges)
    {
        yield return rootId;
        foreach (ClonalEdge e in edges)
        {
            yield return e.ChildId;
        }
    }

    /// <summary>Chooses a synthetic root id distinct from every cluster id (one less than the minimum, or -1).</summary>
    private static int RootIdFor(IReadOnlyList<CcfCluster> clusters)
    {
        int minId = int.MaxValue;
        foreach (CcfCluster c in clusters)
        {
            if (c.Id < minId)
            {
                minId = c.Id;
            }
        }

        // -1 is conventional for "no real cluster"; if a caller used -1, step below the minimum to stay unique.
        return clusters.Count == 0 ? -1 : Math.Min(-1, minId - 1);
    }

    /// <summary>Validates every cluster's CCF list and returns the common sample count.</summary>
    private static int ValidateAndGetSampleCount(IReadOnlyList<CcfCluster> clusters)
    {
        int sampleCount = -1;
        var seenIds = new HashSet<int>(clusters.Count);
        foreach (CcfCluster c in clusters)
        {
            if (c.CcfPerSample is null)
            {
                throw new ArgumentNullException(nameof(clusters), $"Cluster {c.Id} has a null CCF list.");
            }

            if (c.CcfPerSample.Count == 0)
            {
                throw new ArgumentException($"Cluster {c.Id} has an empty CCF list; at least one sample is required.", nameof(clusters));
            }

            if (sampleCount < 0)
            {
                sampleCount = c.CcfPerSample.Count;
            }
            else if (c.CcfPerSample.Count != sampleCount)
            {
                throw new ArgumentException(
                    $"All clusters must have the same number of samples; cluster {c.Id} has {c.CcfPerSample.Count}, expected {sampleCount}.",
                    nameof(clusters));
            }

            if (!seenIds.Add(c.Id))
            {
                throw new ArgumentException($"Duplicate cluster id {c.Id}; cluster ids must be unique.", nameof(clusters));
            }

            foreach (double v in c.CcfPerSample)
            {
                if (double.IsNaN(v) || v < 0.0 || v > 1.0)
                {
                    throw new ArgumentException(
                        $"Cancer cell fraction must be in [0, 1]; cluster {c.Id} has {v}.", nameof(clusters));
                }
            }
        }

        return sampleCount;
    }

    #endregion


    #region Tumor Heterogeneity Analysis (ONCO-HETERO-001)

    /// <summary>
    /// MAD consistency-scaling constant 1.4826 = 1/Φ⁻¹(3/4): scales the raw median absolute deviation so that for
    /// a normally distributed variable the expected MAD equals its standard deviation. Source: Mroz &amp; Rocco
    /// (2013), <i>Oral Oncology</i> 49(3):211–215 / Mroz et al. (2015), <i>PLOS Medicine</i> 12(2):e1001786 —
    /// "The median [absolute deviation] is then multiplied by a factor of 1.4826, so that the expected MAD of a
    /// normally distributed variable is equal to its SD"; maftools <c>mathScore.R</c> uses the same 1.4826.
    /// </summary>
    private const double MadConsistencyConstant = 1.4826;

    /// <summary>
    /// Percentage-scaling factor in the MATH score. Source: Mroz &amp; Rocco (2013) / Mroz et al. (2015):
    /// "MATH = 100 × MAD/median"; maftools <c>mathScore.R</c>: <c>pat.math = pat.mad * 1.4826 / median(vaf)</c>
    /// with <c>pat.mad = median(abs.med.dev) * 100</c>.
    /// </summary>
    private const double MathPercentScale = 100.0;

    /// <summary>
    /// Result of a tumour intratumour-heterogeneity (ITH) analysis over a set of somatic mutations.
    /// </summary>
    /// <param name="MathScore">Mutant-Allele Tumour Heterogeneity (MATH) score = 100·1.4826·MAD(VAF)/median(VAF),
    /// computed over the mutant-allele (variant) fractions (Mroz &amp; Rocco 2013).</param>
    /// <param name="ShannonDiversity">Shannon diversity index H = −Σ pᵢ·ln(pᵢ) over the clone fractions pᵢ
    /// (fraction of mutations assigned to each CCF cluster), using the natural logarithm (Shannon 1948).</param>
    /// <param name="SubcloneCount">Number of distinct clones/subclones = number of non-empty CCF clusters.</param>
    /// <param name="SubclonalFraction">Fraction of mutations whose CCF is below the clonal threshold (CCF &lt; 0.95,
    /// Landau et al. 2013) and are therefore subclonal.</param>
    public readonly record struct HeterogeneityResult(
        double MathScore,
        double ShannonDiversity,
        int SubcloneCount,
        double SubclonalFraction);

    /// <summary>
    /// Computes the intratumour-heterogeneity (ITH) score as the Mutant-Allele Tumour Heterogeneity (MATH) score
    /// over a distribution of mutant-allele (variant) fractions:
    /// <c>MATH = 100 · 1.4826 · median(|fᵢ − median(f)|) / median(f)</c>, the ratio of the width (scaled median
    /// absolute deviation) to the centre (median) of the VAF distribution, expressed as a percentage. A wider,
    /// more dispersed VAF distribution (more genetic heterogeneity) gives a higher MATH. Source: Mroz &amp; Rocco
    /// (2013), <i>Oral Oncology</i> 49(3):211–215; Mroz et al. (2015), <i>PLOS Medicine</i> 12(2):e1001786
    /// ("MATH = 100 × MAD/median"); maftools <c>mathScore.R</c>.
    /// </summary>
    /// <param name="ccfDistribution">Mutant-allele fractions (VAFs) of the tumour's somatic mutations, each in
    /// [0, 1]; the median must be strictly positive (MATH divides by it).</param>
    /// <returns>The MATH score (≥ 0; 0 when every value equals the median).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfDistribution"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, a value is non-finite or outside [0, 1], or the
    /// median is 0 (the MATH ratio is undefined).</exception>
    public static double CalculateITH(IReadOnlyList<double> ccfDistribution)
    {
        ArgumentNullException.ThrowIfNull(ccfDistribution);

        int n = ccfDistribution.Count;
        if (n == 0)
        {
            throw new ArgumentException("At least one allele fraction is required.", nameof(ccfDistribution));
        }

        double[] values = new double[n];
        for (int i = 0; i < n; i++)
        {
            double v = ccfDistribution[i];
            if (double.IsNaN(v) || double.IsInfinity(v) || v < 0.0 || v > 1.0)
            {
                throw new ArgumentException(
                    $"Allele fraction must be a finite value in [0, 1]; got {v} at index {i}.", nameof(ccfDistribution));
            }

            values[i] = v;
        }

        double median = Median(values);
        if (median == 0.0)
        {
            throw new ArgumentException(
                "The median allele fraction is 0; the MATH score divides by the median and is undefined.",
                nameof(ccfDistribution));
        }

        // Raw MAD = median of absolute deviations from the median; scale by 1.4826 for normal consistency.
        double[] absDeviations = new double[n];
        for (int i = 0; i < n; i++)
        {
            absDeviations[i] = Math.Abs(values[i] - median);
        }

        double rawMad = Median(absDeviations);
        double scaledMad = MadConsistencyConstant * rawMad;
        return MathPercentScale * scaledMad / median;
    }

    /// <summary>
    /// Counts the number of distinct clones/subclones in a tumour as the number of non-empty CCF clusters produced
    /// by <see cref="ClusterCcfValues"/> (ONCO-CCF-001). Each cluster centroid is a clonal population; the count is
    /// the tumour's clonal richness. Source: Mroz et al. (2015) treat genetic heterogeneity as the number/spread of
    /// subpopulations; Liu et al. (2017), <i>BMC Genomics</i> 18:457 (PMC5468233) define richness as "the number of
    /// clones present" when computing Shannon-based ITH scores.
    /// </summary>
    /// <param name="ccfClusters">A CCF clustering (its <see cref="CcfClustering.Assignments"/> determine which
    /// clusters actually contain at least one mutation).</param>
    /// <returns>The number of clusters that contain at least one assigned mutation (≥ 1).</returns>
    /// <exception cref="ArgumentException"><paramref name="ccfClusters"/> has no centroids or no assignments.</exception>
    public static int InferSubclones(CcfClustering ccfClusters)
    {
        IReadOnlyList<int> assignments = ccfClusters.Assignments;
        if (assignments is null || assignments.Count == 0 || ccfClusters.Centroids is null || ccfClusters.Centroids.Count == 0)
        {
            throw new ArgumentException("The CCF clustering must contain at least one cluster and one assignment.", nameof(ccfClusters));
        }

        var occupied = new HashSet<int>();
        foreach (int label in assignments)
        {
            occupied.Add(label);
        }

        return occupied.Count;
    }

    /// <summary>
    /// Performs a tumour intratumour-heterogeneity (ITH) analysis from per-mutation variant allele fractions and
    /// cancer cell fractions, returning four standard ITH metrics: the MATH score over the VAFs (Mroz &amp; Rocco
    /// 2013), the Shannon diversity index H = −Σ pᵢ·ln(pᵢ) over the clone fractions (Shannon 1948), the number of
    /// subclones (CCF clusters), and the fraction of subclonal mutations (CCF &lt; 0.95, Landau et al. 2013). CCF
    /// values are clustered with <see cref="ClusterCcfValues"/> (ONCO-CCF-001) into <paramref name="clusterCount"/>
    /// clones; the clone fractions pᵢ are the proportions of mutations assigned to each cluster.
    /// </summary>
    /// <param name="variantAlleleFractions">Per-mutation VAFs in [0, 1] (used for the MATH score).</param>
    /// <param name="ccfValues">Per-mutation cancer cell fractions in [0, 1], aligned with
    /// <paramref name="variantAlleleFractions"/> (same count and order).</param>
    /// <param name="clusterCount">Number of clones k to cluster the CCF values into, in [1, count].</param>
    /// <returns>The MATH score, Shannon diversity, subclone count, and subclonal fraction.</returns>
    /// <exception cref="ArgumentNullException">either list is null.</exception>
    /// <exception cref="ArgumentException">the lists are empty or of different length, or values are out of range.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count].</exception>
    public static HeterogeneityResult AnalyzeHeterogeneity(
        IReadOnlyList<double> variantAlleleFractions,
        IReadOnlyList<double> ccfValues,
        int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(variantAlleleFractions);
        ArgumentNullException.ThrowIfNull(ccfValues);

        if (variantAlleleFractions.Count == 0)
        {
            throw new ArgumentException("At least one mutation is required.", nameof(variantAlleleFractions));
        }

        if (variantAlleleFractions.Count != ccfValues.Count)
        {
            throw new ArgumentException(
                $"VAF and CCF lists must have the same length; got {variantAlleleFractions.Count} and {ccfValues.Count}.",
                nameof(ccfValues));
        }

        double math = CalculateITH(variantAlleleFractions);

        // Cluster CCF values into clones (ONCO-CCF-001); validates ccfValues and clusterCount.
        CcfClustering clustering = ClusterCcfValues(ccfValues, clusterCount);
        int subcloneCount = InferSubclones(clustering);

        // Clone fractions pᵢ = proportion of mutations in each occupied cluster; Shannon H = −Σ pᵢ·ln pᵢ.
        int n = ccfValues.Count;
        var clusterSizes = new Dictionary<int, int>();
        foreach (int label in clustering.Assignments)
        {
            clusterSizes[label] = clusterSizes.TryGetValue(label, out int existing) ? existing + 1 : 1;
        }

        double shannon = 0.0;
        foreach (int size in clusterSizes.Values)
        {
            double p = (double)size / n;
            shannon -= p * Math.Log(p);
        }

        // Subclonal mutations: CCF strictly below the clonal threshold (Landau et al. 2013, reused from ONCO-CLONAL-001).
        int subclonal = 0;
        for (int i = 0; i < n; i++)
        {
            if (ccfValues[i] < ClonalCcfThreshold)
            {
                subclonal++;
            }
        }

        double subclonalFraction = (double)subclonal / n;
        return new HeterogeneityResult(math, shannon, subcloneCount, subclonalFraction);
    }

    /// <summary>
    /// Returns the median of the supplied values. For an even count the median is the arithmetic mean of the two
    /// central order statistics; for an odd count it is the central value (standard definition; matches R's
    /// <c>median</c> used by maftools <c>mathScore.R</c>). Does not mutate the input.
    /// </summary>
    private static double Median(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        int mid = n / 2;
        return (n % 2 == 1) ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2.0;
    }

    #endregion

}
