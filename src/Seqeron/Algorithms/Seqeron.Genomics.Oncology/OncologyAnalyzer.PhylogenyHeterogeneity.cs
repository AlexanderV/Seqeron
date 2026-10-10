namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Tumor Phylogeny Reconstruction (ONCO-PHYLO-001)

    /// <summary>
    /// Cancer cell fraction (CCF) of the synthetic root (normal / germline) node in every sample. LICHeE's cell-prevalence
    /// mode (<c>-cp</c>) sets the root value <c>VAF_MAX = 1.0</c> (LICHeE <c>LineageEngine</c> / <c>PHYNode.getAAF</c>);
    /// Popic et al. (2015), <i>Genome Biology</i> 16:91.
    /// </summary>
    private const double RootCcf = 1.0;

    /// <summary>
    /// Default noise margin ε for the lineage-precedence (Eq. 2) and sum-rule (Eq. 5) inequalities — LICHeE's default
    /// <c>-e</c> = <c>Parameters.VAF_ERROR_MARGIN</c> = 0.1, which also applies in cell-prevalence mode (<c>-cp</c> only
    /// sets <c>VAF_MAX</c> = <c>MAX_ALLOWED_VAF</c> = 1; <c>LineageEngine</c> option parsing, github.com/viq854/lichee).
    /// Pass 0 for the strict inequalities.
    /// </summary>
    public const double DefaultPhylogenyTolerance = 0.1;

    /// <summary>
    /// Maximum number of valid spanning trees enumerated before the search stops — LICHeE
    /// <c>Parameters.MAX_NUM_TREES</c> = 100 000 (github.com/viq854/lichee).
    /// </summary>
    public const int MaxPhylogenyTreesEnumerated = 100_000;

    /// <summary>
    /// Maximum number of recursive <c>grow</c> calls of the Gabow–Myers spanning-tree enumeration — LICHeE
    /// <c>Parameters.MAX_NUM_GROW_CALLS</c> = 10⁸. When the budget is exhausted the enumeration stops and the best tree
    /// found so far is returned. LICHeE's <c>grow</c> returns from the frame that hit the cap and its callers keep looping
    /// over their remaining stack edges, but a recursion is only entered below the cap, so no tree is added after it:
    /// the tree list — and hence the ranking — equals the trees completed within the first 10⁸ calls, as here.
    /// </summary>
    public const int MaxPhylogenyGrowCalls = 100_000_000;

    /// <summary>
    /// One CCF cluster (a candidate clone/subclone) used as input to phylogeny reconstruction. Each cluster carries
    /// its cancer cell fraction in each sequenced sample. CCF clustering itself is out of scope (ONCO-CCF-001).
    /// </summary>
    /// <param name="Id">Caller-assigned cluster identifier (e.g. a subclone label).</param>
    /// <param name="CcfPerSample">Cancer cell fraction in each sample, each value in [0, 1]; all clusters must share length.
    /// A value of exactly 0 means "absent from that sample" (LICHeE presence profile).</param>
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
    /// <param name="Edges">Tree edges (parent→child), one per cluster, in input order of the child.</param>
    /// <param name="SampleCount">Number of samples per cluster.</param>
    public readonly record struct ClonalPhylogeny(
        int RootId,
        IReadOnlyList<CcfCluster> Clusters,
        IReadOnlyList<ClonalEdge> Edges,
        int SampleCount)
    {
        /// <summary>Noise margin ε the tree was reconstructed with (also used by <see cref="IdentifyTrunkMutations"/>).</summary>
        public double Tolerance { get; init; }

        /// <summary>
        /// LICHeE error score of the returned tree: √(Σ_nodes Σ_samples max(0, Σ_children CCF − node CCF)²)
        /// (<c>PHYTree.computeErrorScore</c>). 0 when every sum rule holds without using the margin ε.
        /// </summary>
        public double ErrorScore { get; init; }

        /// <summary>Number of valid spanning trees enumerated (LICHeE "Found N valid tree(s)"), capped at <see cref="MaxPhylogenyTreesEnumerated"/>.</summary>
        public int ValidTreeCount { get; init; }

        /// <summary>
        /// True when the default constraint network admitted no valid tree and the complete network
        /// (LICHeE <c>-c</c> / <c>ALL_EDGES</c> fallback) was used.
        /// </summary>
        public bool UsedCompleteNetwork { get; init; }

        /// <summary>
        /// Ids of the non-robust clusters LICHeE's <c>fixNetwork</c> dropped (in removal order) before a valid tree was
        /// found; they are absent from <see cref="Clusters"/> and <see cref="Edges"/>. Empty when nothing was removed
        /// (always for <see cref="ReconstructPhylogeny"/>).
        /// </summary>
        public IReadOnlyList<int> RemovedClusterIds { get; init; } = Array.Empty<int>();

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

        /// <summary>Returns the ids of the direct children of <paramref name="clusterId"/>, in edge order.</summary>
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
    /// Reconstructs a rooted clonal (tumor) phylogeny from per-sample CCF clusters with the LICHeE algorithm
    /// (Popic et al. 2015, <i>Genome Biology</i> 16:91; reference code github.com/viq854/lichee, run in cell-prevalence
    /// mode <c>-cp</c> on pre-computed clusters):
    /// <list type="number">
    /// <item><description><b>Constraint network</b> (<c>PHYNetwork</c>): nodes are the clusters, levelled by the number of
    /// samples they are present in (CCF &gt; 0); the root sits above every level. An edge u→v is admissible iff for every
    /// sample i <c>u.CCF[i] ≥ v.CCF[i] − ε</c> and <c>u.CCF[i] = 0 ⇒ v.CCF[i] = 0</c> (Eq. 2); when both directions are
    /// admissible the one with the smaller one-sided CCF excess is kept. Edges are added within a presence profile,
    /// between each level and the next non-empty lower level, and nodes left without a parent are linked to the closest
    /// higher level that admits them, else to the root.</description></item>
    /// <item><description><b>Tree search</b>: every spanning tree of the network rooted at the normal node that satisfies
    /// the sum rule <c>Σ_children v.CCF[i] ≤ u.CCF[i] + ε</c> (Eq. 5) is enumerated (Gabow &amp; Myers 1978, as in
    /// LICHeE <c>grow</c>).</description></item>
    /// <item><description><b>Ranking</b>: trees are ranked by the LICHeE error score (stable order); the top-ranking
    /// tree is returned. If the default network admits no valid tree, the complete network (<c>ALL_EDGES</c>) is
    /// searched, as LICHeE does.</description></item>
    /// </list>
    /// Presence profiles are grouped in order of first appearance in <paramref name="clusters"/> (LICHeE iterates a
    /// <c>HashMap</c>), which fixes the enumeration order and therefore the tie-break among equal-score trees.
    /// </summary>
    /// <param name="clusters">CCF clusters to place; each cluster's <see cref="CcfCluster.CcfPerSample"/> must have the same length.</param>
    /// <param name="tolerance">Noise margin ε for both inequalities; default <see cref="DefaultPhylogenyTolerance"/> (0.1, LICHeE <c>-e</c>).</param>
    /// <returns>The top-ranking <see cref="ClonalPhylogeny"/> rooted at a synthetic normal node.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="clusters"/> or any cluster's CCF list is null.</exception>
    /// <exception cref="ArgumentException">CCF lists differ in length, are empty, or contain NaN / out-of-[0,1] values; or two clusters share an id.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tolerance"/> is negative or NaN.</exception>
    /// <exception cref="InvalidOperationException">No spanning tree satisfies the sum rule (LICHeE finds no valid tree); use <see cref="TryReconstructPhylogeny"/>.</exception>
    public static ClonalPhylogeny ReconstructPhylogeny(
        IReadOnlyList<CcfCluster> clusters,
        double tolerance = DefaultPhylogenyTolerance)
    {
        if (!TryReconstructPhylogeny(clusters, out ClonalPhylogeny phylogeny, tolerance))
        {
            throw new InvalidOperationException(
                "No lineage tree satisfies the sum rule (LICHeE Eq. 5) for these CCF clusters at the given tolerance; "
                + "increase the tolerance or re-cluster the CCFs.");
        }

        return phylogeny;
    }

    /// <summary>
    /// Non-throwing variant of <see cref="ReconstructPhylogeny"/>: returns false (and a default phylogeny) when no
    /// spanning tree of the LICHeE constraint network satisfies the sum rule, also after the complete-network fallback.
    /// Argument validation still throws as in <see cref="ReconstructPhylogeny"/>.
    /// </summary>
    public static bool TryReconstructPhylogeny(
        IReadOnlyList<CcfCluster> clusters,
        out ClonalPhylogeny phylogeny,
        double tolerance = DefaultPhylogenyTolerance)
    {
        ArgumentNullException.ThrowIfNull(clusters);
        ValidatePhylogenyTolerance(tolerance);

        int rootId = RootIdFor(clusters);
        if (clusters.Count == 0)
        {
            // Empty cohort: tree is the root alone, no clusters, no edges. One synthetic sample of CCF 1.
            phylogeny = new ClonalPhylogeny(rootId, Array.Empty<CcfCluster>(), Array.Empty<ClonalEdge>(), 1)
            {
                Tolerance = tolerance,
                ValidTreeCount = 1,
            };
            return true;
        }

        int sampleCount = ValidateAndGetSampleCount(clusters);
        return TryReconstructLichee(clusters, sampleCount, tolerance, standardErrors: null, memberCounts: null, robust: null, out phylogeny);
    }

    /// <summary>
    /// Test hook: <see cref="TryReconstructPhylogeny"/> with the grow-call budget <paramref name="maxGrowCalls"/> instead of
    /// <see cref="MaxPhylogenyGrowCalls"/> (LICHeE <c>Parameters.MAX_NUM_GROW_CALLS</c>), to lock the post-cap result.
    /// </summary>
    internal static bool TryReconstructPhylogeny(
        IReadOnlyList<CcfCluster> clusters,
        out ClonalPhylogeny phylogeny,
        double tolerance,
        int maxGrowCalls)
    {
        ArgumentNullException.ThrowIfNull(clusters);
        ValidatePhylogenyTolerance(tolerance);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxGrowCalls, 1);
        if (clusters.Count == 0)
        {
            return TryReconstructPhylogeny(clusters, out phylogeny, tolerance);
        }

        int sampleCount = ValidateAndGetSampleCount(clusters);
        return TryReconstructLichee(
            clusters, sampleCount, tolerance, standardErrors: null, memberCounts: null, robust: null, out phylogeny, maxGrowCalls);
    }

    /// <summary>
    /// LICHeE <c>LineageEngine.buildLineage</c> steps 4–6 on validated clusters: constraint network, tree search,
    /// network adjustment, ranking. <paramref name="standardErrors"/> = per-cluster <c>1.96·sd/√n</c> (null = static ε).
    /// When no tree exists: with <paramref name="memberCounts"/>/<paramref name="robust"/> supplied, LICHeE's
    /// <c>fixNetwork</c> loop drops the smallest non-robust cluster (first in node-id order on ties) and searches again
    /// until a tree is found or no non-robust cluster is left; then the complete network (<c>ALL_EDGES</c>) is
    /// searched on the remaining clusters.
    /// </summary>
    private static bool TryReconstructLichee(
        IReadOnlyList<CcfCluster> clusters,
        int sampleCount,
        double tolerance,
        double[][]? standardErrors,
        int[]? memberCounts,
        bool[]? robust,
        out ClonalPhylogeny phylogeny,
        int maxGrowCalls = MaxPhylogenyGrowCalls)
    {
        int rootId = RootIdFor(clusters);
        var nodeOrder = new List<int>(LicheeNodeOrder(clusters, sampleCount));
        var removed = new List<int>();
        bool usedComplete = false;
        LicheeSearchResult result = new LicheeNetwork(
            clusters, sampleCount, tolerance, completeNetwork: false, nodeOrder, standardErrors, maxGrowCalls).Search();
        if (result.TreeCount == 0 && memberCounts is not null && robust is not null)
        {
            // fixNetwork: iterate nodes in id order, keep the first non-robust cluster of strictly smallest size.
            while (result.TreeCount == 0)
            {
                int toRemove = -1;
                foreach (int c in nodeOrder)
                {
                    if (!robust[c] && (toRemove < 0 || memberCounts[c] < memberCounts[toRemove]))
                    {
                        toRemove = c;
                    }
                }

                if (toRemove < 0)
                {
                    break; // no node removed (delta = 0): the rebuilt network is unchanged.
                }

                nodeOrder.Remove(toRemove);
                removed.Add(clusters[toRemove].Id);
                result = new LicheeNetwork(
                    clusters, sampleCount, tolerance, completeNetwork: false, nodeOrder, standardErrors, maxGrowCalls).Search();
            }
        }

        if (result.TreeCount == 0)
        {
            usedComplete = true;
            result = new LicheeNetwork(
                clusters, sampleCount, tolerance, completeNetwork: true, nodeOrder, standardErrors, maxGrowCalls).Search();
        }

        if (result.TreeCount == 0)
        {
            phylogeny = default;
            return false;
        }

        var kept = new List<CcfCluster>(clusters.Count - removed.Count);
        var edges = new List<ClonalEdge>(clusters.Count - removed.Count);
        for (int c = 0; c < clusters.Count; c++)
        {
            int parentCluster = result.ParentClusterIndex[c];
            if (parentCluster == int.MinValue)
            {
                continue; // removed by fixNetwork
            }

            kept.Add(clusters[c]);
            edges.Add(new ClonalEdge(parentCluster < 0 ? rootId : clusters[parentCluster].Id, clusters[c].Id));
        }

        phylogeny = new ClonalPhylogeny(rootId, kept, edges, sampleCount)
        {
            Tolerance = tolerance,
            ErrorScore = result.ErrorScore,
            ValidTreeCount = result.TreeCount,
            UsedCompleteNetwork = usedComplete,
            RemovedClusterIds = removed,
        };
        return true;
    }

    /// <summary>
    /// LICHeE's per-cluster standard-error multiplier: the edge margin uses <c>1.96·sd/√n</c> per cluster and sample
    /// (<c>PHYNetwork.getAAFErrorMargin</c>, github.com/viq854/lichee) — the two-sided 95 % normal quantile.
    /// </summary>
    public const double LicheeStandardErrorZ = 1.96;

    /// <summary>
    /// LICHeE <c>Parameters.MIN_ROBUST_CLUSTER_SUPPORT</c> = 2: a cluster is robust iff at least this many of its member
    /// mutations are robust (<c>SNVGroup.setSubPopulations</c>); <c>PHYNetwork.fixNetwork</c> drops non-robust clusters.
    /// </summary>
    public const int LicheeMinRobustClusterSupport = 2;

    /// <summary>
    /// A CCF cluster with the dispersion summary LICHeE keeps per cluster (<c>AAFClusterer.Cluster</c>): centroid
    /// (mean member CCF per sample), per-sample standard deviation of the member CCFs and member count. Used by
    /// <see cref="ReconstructPhylogenyFromClusterSummaries"/> for LICHeE's per-cluster edge error margins.
    /// </summary>
    /// <param name="Id">Caller-assigned cluster identifier.</param>
    /// <param name="CcfPerSample">Centroid CCF per sample, each in [0, 1] (0 = absent, LICHeE presence profile).</param>
    /// <param name="StdDevPerSample">Standard deviation of the member CCFs per sample (LICHeE: population SD, divisor
    /// n — <c>Cluster.recomputeCentroidAndStdDev</c> / EM <c>setStdDev</c>), finite and ≥ 0; ignored (treated as 0)
    /// where the centroid is 0, as LICHeE <c>PHYNode.getStdDev</c> returns 0 for samples outside the profile.</param>
    /// <param name="MemberCount">Number of member mutations n ≥ 1 (<c>getMembership().size()</c>).</param>
    public readonly record struct CcfClusterSummary(
        int Id,
        IReadOnlyList<double> CcfPerSample,
        IReadOnlyList<double> StdDevPerSample,
        int MemberCount)
    {
        /// <summary>
        /// Number of robust member mutations (LICHeE <c>SNVEntry.isRobust</c>: no sample VAF in the ambiguous band
        /// [<c>MAX_VAF_ABSENT</c>, <c>MIN_VAF_PRESENT</c>), empty at the defaults 0.005/0.005). Null (default) = every
        /// member is robust (= <see cref="MemberCount"/>). The cluster is robust iff this is ≥
        /// <see cref="LicheeMinRobustClusterSupport"/>; 0 reproduces LICHeE's <c>--clustersFile</c> input, whose clusters
        /// are never marked robust.
        /// </summary>
        public int? RobustMemberCount { get; init; }

        /// <summary>
        /// Builds the summary from member-level CCFs exactly as LICHeE does
        /// (<c>Cluster.recomputeCentroidAndStdDev</c>): centroid = Σ member / n, SD = √(Σ (member − centroid)² / n)
        /// per sample (population SD, <see cref="StatisticsHelper.PopulationVariance"/>), n = member count.
        /// </summary>
        /// <param name="id">Cluster identifier.</param>
        /// <param name="memberCcfs">One CCF vector per member mutation (all of equal length ≥ 1, values in [0, 1]).</param>
        /// <exception cref="ArgumentNullException"><paramref name="memberCcfs"/> or a member vector is null.</exception>
        /// <exception cref="ArgumentException">No members, empty or ragged vectors.</exception>
        public static CcfClusterSummary FromMembers(int id, IReadOnlyList<IReadOnlyList<double>> memberCcfs)
        {
            ArgumentNullException.ThrowIfNull(memberCcfs);
            if (memberCcfs.Count == 0)
            {
                throw new ArgumentException($"Cluster {id} has no members.", nameof(memberCcfs));
            }

            int samples = -1;
            foreach (IReadOnlyList<double> m in memberCcfs)
            {
                if (m is null)
                {
                    throw new ArgumentNullException(nameof(memberCcfs), $"Cluster {id} has a null member CCF vector.");
                }

                if (samples < 0)
                {
                    samples = m.Count;
                }
                else if (m.Count != samples)
                {
                    throw new ArgumentException($"Cluster {id}: all member CCF vectors must have the same length.", nameof(memberCcfs));
                }
            }

            if (samples == 0)
            {
                throw new ArgumentException($"Cluster {id}: member CCF vectors are empty.", nameof(memberCcfs));
            }

            int n = memberCcfs.Count;
            var centroid = new double[samples];
            var sd = new double[samples];
            var column = new double[n];
            for (int s = 0; s < samples; s++)
            {
                double sum = 0.0;
                for (int m = 0; m < n; m++)
                {
                    column[m] = memberCcfs[m][s];
                    sum += column[m];
                }

                centroid[s] = sum / n;
                sd[s] = Math.Sqrt(StatisticsHelper.PopulationVariance(column));
            }

            return new CcfClusterSummary(id, centroid, sd, n);
        }
    }

    /// <summary>
    /// <see cref="ReconstructPhylogeny"/> with LICHeE's per-cluster edge error margins
    /// (<c>PHYNetwork.getAAFErrorMargin</c>): the lineage-precedence test (Eq. 2) of an edge u→v in sample i uses
    /// <c>max(ε, se_u,i + se_v,i)</c> with <c>se = 1.96·sd/√n</c> for a cluster (<see cref="LicheeStandardErrorZ"/>)
    /// and <c>se = ε</c> for the root, instead of the static ε. The sum rule (Eq. 5, <c>PHYTree.checkConstraint</c>)
    /// keeps the static ε, as in LICHeE. With every SD = 0 and every cluster robust the result is identical to
    /// <see cref="ReconstructPhylogeny"/>.
    /// <para>
    /// Network adjustment (LICHeE <c>LineageEngine.buildLineage</c> step 5, <c>PHYNetwork.fixNetwork</c>): when no tree
    /// exists, the smallest non-robust cluster (fewer than <see cref="LicheeMinRobustClusterSupport"/> robust members,
    /// see <see cref="CcfClusterSummary.RobustMemberCount"/>; ties → first in network node order) is dropped and the
    /// search repeated until a tree is found or no non-robust cluster is left; then the complete network
    /// (<c>ALL_EDGES</c>) is searched on the remaining clusters. Dropped ids are reported in
    /// <see cref="ClonalPhylogeny.RemovedClusterIds"/>. LICHeE always runs this step; pass
    /// <paramref name="removeNonRobustClusters"/> = false to keep every cluster.
    /// </para>
    /// </summary>
    /// <param name="clusters">Cluster summaries (centroid, SD, member count), e.g. from <see cref="CcfClusterSummary.FromMembers"/>.</param>
    /// <param name="tolerance">Static margin ε (LICHeE <c>-e</c>); default <see cref="DefaultPhylogenyTolerance"/>.</param>
    /// <param name="removeNonRobustClusters">Run LICHeE's <c>fixNetwork</c> cluster removal (default true, as LICHeE).</param>
    /// <returns>The top-ranking phylogeny; <see cref="ClonalPhylogeny.Clusters"/> holds the centroids of the kept clusters.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="clusters"/>, a CCF or an SD list is null.</exception>
    /// <exception cref="ArgumentException">As <see cref="ReconstructPhylogeny"/>; SD lists of the wrong length or with
    /// negative / non-finite values; member count &lt; 1; robust member count outside [0, member count].</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tolerance"/> is negative or NaN.</exception>
    /// <exception cref="InvalidOperationException">No valid lineage tree exists.</exception>
    public static ClonalPhylogeny ReconstructPhylogenyFromClusterSummaries(
        IReadOnlyList<CcfClusterSummary> clusters,
        double tolerance = DefaultPhylogenyTolerance,
        bool removeNonRobustClusters = true)
    {
        if (!TryReconstructPhylogenyFromClusterSummaries(clusters, out ClonalPhylogeny phylogeny, tolerance, removeNonRobustClusters))
        {
            throw new InvalidOperationException(
                "No lineage tree satisfies the sum rule (LICHeE Eq. 5) for these CCF clusters at the given tolerance; "
                + "increase the tolerance or re-cluster the CCFs.");
        }

        return phylogeny;
    }

    /// <summary>
    /// Non-throwing variant of <see cref="ReconstructPhylogenyFromClusterSummaries"/>: false (and a default phylogeny)
    /// when no valid lineage tree exists, also after <c>fixNetwork</c> and the complete-network fallback. Argument
    /// validation still throws.
    /// </summary>
    public static bool TryReconstructPhylogenyFromClusterSummaries(
        IReadOnlyList<CcfClusterSummary> clusters,
        out ClonalPhylogeny phylogeny,
        double tolerance = DefaultPhylogenyTolerance,
        bool removeNonRobustClusters = true)
    {
        ArgumentNullException.ThrowIfNull(clusters);
        ValidatePhylogenyTolerance(tolerance);
        var centroids = new CcfCluster[clusters.Count];
        for (int c = 0; c < clusters.Count; c++)
        {
            centroids[c] = new CcfCluster(clusters[c].Id, clusters[c].CcfPerSample);
        }

        if (clusters.Count == 0)
        {
            return TryReconstructPhylogeny(centroids, out phylogeny, tolerance);
        }

        int sampleCount = ValidateAndGetSampleCount(centroids);
        var standardErrors = new double[clusters.Count][];
        var memberCounts = new int[clusters.Count];
        var robust = new bool[clusters.Count];
        for (int c = 0; c < clusters.Count; c++)
        {
            CcfClusterSummary summary = clusters[c];
            if (summary.StdDevPerSample is null)
            {
                throw new ArgumentNullException(nameof(clusters), $"Cluster {summary.Id} has a null SD list.");
            }

            if (summary.StdDevPerSample.Count != sampleCount)
            {
                throw new ArgumentException(
                    $"Cluster {summary.Id} has {summary.StdDevPerSample.Count} SD values, expected {sampleCount}.", nameof(clusters));
            }

            if (summary.MemberCount < 1)
            {
                throw new ArgumentException($"Cluster {summary.Id} must have at least one member.", nameof(clusters));
            }

            int robustMembers = summary.RobustMemberCount ?? summary.MemberCount;
            if (robustMembers < 0 || robustMembers > summary.MemberCount)
            {
                throw new ArgumentException(
                    $"Cluster {summary.Id}: robust member count {robustMembers} must be in [0, {summary.MemberCount}].", nameof(clusters));
            }

            memberCounts[c] = summary.MemberCount;
            robust[c] = robustMembers >= LicheeMinRobustClusterSupport;
            double rootN = Math.Sqrt((double)summary.MemberCount);
            standardErrors[c] = new double[sampleCount];
            for (int s = 0; s < sampleCount; s++)
            {
                double sd = summary.StdDevPerSample[s];
                if (!double.IsFinite(sd) || sd < 0.0)
                {
                    throw new ArgumentException(
                        $"Standard deviation must be finite and non-negative; cluster {summary.Id} has {sd}.", nameof(clusters));
                }

                // PHYNode.getStdDev: 0 for a sample outside the presence profile.
                double effectiveSd = summary.CcfPerSample[s] > 0.0 ? sd : 0.0;
                standardErrors[c][s] = LicheeStandardErrorZ * effectiveSd / rootN;
            }
        }

        return removeNonRobustClusters
            ? TryReconstructLichee(centroids, sampleCount, tolerance, standardErrors, memberCounts, robust, out phylogeny)
            : TryReconstructLichee(centroids, sampleCount, tolerance, standardErrors, memberCounts: null, robust: null, out phylogeny);
    }

    /// <summary>
    /// Returns the ids of clusters on the <b>trunk</b> of the phylogeny — the truncal (clonal) mutations present in
    /// every cancer cell of every sample. A trunk alteration "must be present in all cells of the tumour" (Werner et
    /// al. 2017, <i>Sci. Rep.</i> 7:44991), i.e. CCF = 1 in every sample; in the tree these are the nodes on the path
    /// from the root whose CCF equals the root's (≥ 1 − ε, ε = <see cref="ClonalPhylogeny.Tolerance"/>) in every sample.
    /// A single-child descendant with CCF &lt; 1 is subclonal (the parent keeps a residual population without it).
    /// </summary>
    /// <param name="phylogeny">A phylogeny produced by <see cref="ReconstructPhylogeny"/>.</param>
    /// <returns>Trunk cluster ids, ordered from the root downward; empty if no cluster is clonal in every sample.</returns>
    public static IReadOnlyList<int> IdentifyTrunkMutations(ClonalPhylogeny phylogeny)
    {
        var trunk = new List<int>();
        if (phylogeny.Clusters is null || phylogeny.Clusters.Count == 0)
        {
            return trunk;
        }

        var ccfById = new Dictionary<int, IReadOnlyList<double>>(phylogeny.Clusters.Count);
        foreach (CcfCluster c in phylogeny.Clusters)
        {
            ccfById[c.Id] = c.CcfPerSample;
        }

        double clonalThreshold = RootCcf - phylogeny.Tolerance;
        int current = phylogeny.RootId;
        while (true)
        {
            int next = 0;
            int clonalChildren = 0;
            foreach (int child in phylogeny.ChildrenOf(current))
            {
                if (ccfById.TryGetValue(child, out IReadOnlyList<double>? ccf) && IsClonalInEverySample(ccf, clonalThreshold))
                {
                    next = child;
                    clonalChildren++;
                }
            }

            // Exactly one clonal child continues the trunk; none (or an ambiguous pair, only possible for ε ≥ 0.5)
            // ends it.
            if (clonalChildren != 1)
            {
                break;
            }

            trunk.Add(next);
            current = next;
        }

        return trunk;
    }

    /// <summary>
    /// Returns the ids of <b>branch</b> (subclonal) clusters — every input cluster that is not on the trunk.
    /// Source: Werner et al. (2017) — alterations off the trunk are present in only a subset of tumour cells.
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

    private static bool IsClonalInEverySample(IReadOnlyList<double> ccf, double threshold)
    {
        foreach (double v in ccf)
        {
            if (v < threshold)
            {
                return false;
            }
        }

        return true;
    }

    private readonly record struct LicheeSearchResult(int TreeCount, double ErrorScore, int[] ParentClusterIndex);

    /// <summary>
    /// Port of LICHeE's <c>PHYNetwork</c> (constraint network construction, <c>checkAndAddEdge</c>,
    /// <c>getLineageTrees</c>/<c>grow</c> Gabow–Myers enumeration with the <c>PHYTree.checkConstraint</c> sum rule) and
    /// <c>PHYTree.computeErrorScore</c> ranking. Node 0 is the root; nodes 1..n are the clusters grouped by presence
    /// profile (first-appearance order), within a profile in input order — the LICHeE node-id order.
    /// </summary>
    private sealed class LicheeNetwork
    {
        private readonly int _samples;
        private readonly double _eps;
        private readonly int _nodeCount;
        private readonly int _clusterTotal;
        private readonly double[][] _ccf;
        private readonly double[][]? _se;
        private readonly int[] _level;
        private readonly int[] _clusterIndex;
        private readonly SortedDictionary<int, List<int>> _levels = new();
        private readonly List<int>[] _net;
        private readonly int _maxGrowCalls;

        // Search state (LICHeE grow): f stack, working tree t (= L after the first complete tree).
        private readonly List<(int From, int To)> _f = new();
        private readonly List<int> _treeNodes = new();
        private readonly Dictionary<int, List<int>> _treeEdges = new();
        private bool _haveL;
        private int _growCalls;
        private int _treeCount;
        private bool _stop;
        private double _bestError = double.PositiveInfinity;
        private int[]? _bestParent;

        /// <param name="clusters">All input clusters (indexed by input position).</param>
        /// <param name="samples">Sample count.</param>
        /// <param name="eps">Static margin ε (<c>VAF_ERROR_MARGIN</c>).</param>
        /// <param name="completeNetwork">LICHeE <c>ALL_EDGES</c>.</param>
        /// <param name="nodeOrder">Input indices of the clusters in the network, in LICHeE node-id order (grouped by
        /// presence profile, see <see cref="LicheeNodeOrder"/>); clusters not listed are absent (removed by fixNetwork).</param>
        /// <param name="standardErrors">Per input cluster and sample, LICHeE's <c>1.96·sd/√n</c> (null = static ε,
        /// i.e. every cluster's standard error is 0).</param>
        /// <param name="maxGrowCalls">Grow-call budget (<see cref="MaxPhylogenyGrowCalls"/>; smaller only from the test hook).</param>
        public LicheeNetwork(
            IReadOnlyList<CcfCluster> clusters,
            int samples,
            double eps,
            bool completeNetwork,
            IReadOnlyList<int> nodeOrder,
            double[][]? standardErrors,
            int maxGrowCalls)
        {
            _maxGrowCalls = maxGrowCalls;
            _samples = samples;
            _eps = eps;
            _nodeCount = nodeOrder.Count + 1;
            _clusterTotal = clusters.Count;
            _ccf = new double[_nodeCount][];
            _se = standardErrors is null ? null : new double[_nodeCount][];
            _level = new int[_nodeCount];
            _clusterIndex = new int[_nodeCount];
            _net = new List<int>[_nodeCount];
            for (int i = 0; i < _nodeCount; i++)
            {
                _net[i] = new List<int>();
            }

            // Root (LICHeE: level numSamples + 1, AAF = VAF_MAX = 1 in -cp mode).
            _ccf[0] = new double[samples];
            Array.Fill(_ccf[0], RootCcf);
            _clusterIndex[0] = -1;
            AddNode(0, samples + 1);

            // Clusters grouped by presence profile (nodeOrder is already grouped); node ids follow nodeOrder.
            int first = 1;
            string? groupKey = null;
            for (int k = 0; k < nodeOrder.Count; k++)
            {
                int c = nodeOrder[k];
                int id = k + 1;
                string key = PresenceProfile(clusters[c].CcfPerSample, samples);
                if (groupKey is not null && key != groupKey)
                {
                    AddIntraGroupEdges(first, id);
                    first = id;
                }

                groupKey = key;
                _ccf[id] = clusters[c].CcfPerSample.ToArray();
                if (_se is not null)
                {
                    _se[id] = standardErrors![c];
                }

                _clusterIndex[id] = c;
                int present = 0;
                foreach (double v in _ccf[id])
                {
                    if (v > 0.0)
                    {
                        present++;
                    }
                }

                AddNode(id, present);
            }

            if (nodeOrder.Count > 0)
            {
                AddIntraGroupEdges(first, _nodeCount);
            }

            // Inter-level edges: each level to the next non-empty lower level.
            for (int i = samples + 1; i > 0; i--)
            {
                if (!_levels.TryGetValue(i, out List<int>? fromLevel))
                {
                    continue;
                }

                int j = i - 1;
                _levels.TryGetValue(j, out List<int>? toLevel);
                while (toLevel is null && j > 0)
                {
                    j--;
                    _levels.TryGetValue(j, out toLevel);
                }

                if (toLevel is null)
                {
                    continue;
                }

                foreach (int n1 in fromLevel)
                {
                    foreach (int n2 in toLevel)
                    {
                        CheckAndAddEdge(n1, n2);
                    }
                }
            }

            if (completeNetwork)
            {
                // addAllHiddenEdges: every higher level to every lower level ≥ 1.
                for (int i = samples + 1; i > 0; i--)
                {
                    if (!_levels.TryGetValue(i, out List<int>? fromLevel))
                    {
                        continue;
                    }

                    for (int j = i - 1; j >= 1; j--)
                    {
                        if (!_levels.TryGetValue(j, out List<int>? toLevel))
                        {
                            continue;
                        }

                        foreach (int n1 in fromLevel)
                        {
                            foreach (int n2 in toLevel)
                            {
                                CheckAndAddEdge(n1, n2);
                            }
                        }
                    }
                }
            }

            // Nodes with no incoming edge: connect to a valid node in the closest higher level, else to the root.
            var hasParent = new bool[_nodeCount];
            for (int u = 0; u < _nodeCount; u++)
            {
                foreach (int v in _net[u])
                {
                    hasParent[v] = true;
                }
            }

            for (int n = 1; n < _nodeCount; n++)
            {
                if (hasParent[n])
                {
                    continue;
                }

                bool found = false;
                for (int j = _level[n] + 2; j <= samples + 1 && !found; j++)
                {
                    if (!_levels.TryGetValue(j, out List<int>? fromLevel))
                    {
                        continue;
                    }

                    foreach (int n2 in fromLevel)
                    {
                        if (CheckAndAddEdge(n2, n) == 0)
                        {
                            found = true;
                            break;
                        }
                    }
                }

                if (!found)
                {
                    AddNetEdge(0, n);
                }
            }
        }

        public LicheeSearchResult Search()
        {
            // getLineageTrees: t = {root}; f = all (root, v).
            _treeNodes.Add(0);
            if (_net[0].Count == 0)
            {
                return new LicheeSearchResult(0, double.NaN, Array.Empty<int>());
            }

            foreach (int v in _net[0])
            {
                _f.Add((0, v));
            }

            Grow();
            return new LicheeSearchResult(_treeCount, _bestError, _bestParent ?? Array.Empty<int>());
        }

        private void AddNode(int id, int level)
        {
            _level[id] = level;
            if (!_levels.TryGetValue(level, out List<int>? list))
            {
                list = new List<int>();
                _levels[level] = list;
            }

            list.Add(id);
        }

        /// <summary>Edges between each group's sub-population nodes (node ids [first, end)).</summary>
        private void AddIntraGroupEdges(int first, int end)
        {
            for (int i = first; i < end; i++)
            {
                for (int j = i + 1; j < end; j++)
                {
                    CheckAndAddEdge(i, j);
                }
            }
        }

        /// <summary>
        /// LICHeE <c>PHYNetwork.getAAFErrorMargin(from, to, i)</c>: <c>max(ε, se_from + se_to)</c> with
        /// <c>se = 1.96·sd/√n</c> for a cluster and <c>se = ε</c> for the root; the static margin ε when no per-cluster
        /// dispersion was supplied.
        /// </summary>
        private double ErrorMargin(int from, int to, int i)
        {
            if (_se is null)
            {
                return _eps;
            }

            double parentStdError = from == 0 ? _eps : _se[from]![i];
            double childStdError = to == 0 ? _eps : _se[to]![i];
            double standardError = parentStdError + childStdError;
            return standardError > _eps ? standardError : _eps;
        }

        private void AddNetEdge(int from, int to)
        {
            if (!_net[from].Contains(to))
            {
                _net[from].Add(to);
            }
        }

        /// <summary>LICHeE <c>checkAndAddEdge</c>: 0 = added n1→n2, 1 = added n2→n1, −1 = none.</summary>
        private int CheckAndAddEdge(int n1, int n2)
        {
            double[] a1 = _ccf[n1];
            double[] a2 = _ccf[n2];
            int comp12 = 0;
            int comp21 = 0;
            double err12 = 0.0;
            double err21 = 0.0;
            for (int i = 0; i < _samples; i++)
            {
                if (a1[i] == 0.0 && a2[i] != 0.0)
                {
                    break;
                }

                comp12 += a1[i] >= a2[i] - ErrorMargin(n1, n2, i) ? 1 : 0;
                if (a1[i] < a2[i])
                {
                    err12 += a2[i] - a1[i];
                }
            }

            for (int i = 0; i < _samples; i++)
            {
                if (a2[i] == 0.0 && a1[i] != 0.0)
                {
                    break;
                }

                comp21 += a2[i] >= a1[i] - ErrorMargin(n2, n1, i) ? 1 : 0;
                if (a2[i] < a1[i])
                {
                    err21 += a1[i] - a2[i];
                }
            }

            if (comp12 == _samples)
            {
                if (comp21 == _samples && !(err12 < err21))
                {
                    AddNetEdge(n2, n1);
                    return 1;
                }

                AddNetEdge(n1, n2);
                return 0;
            }

            if (comp21 == _samples)
            {
                AddNetEdge(n2, n1);
                return 1;
            }

            return -1;
        }

        // ---- Gabow & Myers (1978) spanning-tree enumeration, LICHeE PHYNetwork.grow ----

        private void Grow()
        {
            _growCalls++;
            if (_treeNodes.Count == _nodeCount)
            {
                _haveL = true;
                _treeCount++;
                double error = ComputeErrorScore();
                if (error < _bestError)
                {
                    // Stable sort by error score: the first tree reaching the minimum ranks first.
                    _bestError = error;
                    _bestParent = CurrentParents();
                }

                if (_treeCount == MaxPhylogenyTreesEnumerated)
                {
                    _stop = true;
                }

                return;
            }

            var ff = new List<(int From, int To)>();
            bool b = false;
            while (!b && _f.Count > 0)
            {
                (int From, int To) e = _f[^1];
                _f.RemoveAt(_f.Count - 1);
                int v = e.To;
                TreeAddNode(v);
                TreeAddEdge(e.From, v);

                if (CheckConstraint(e.From))
                {
                    var edgesAdded = new List<(int From, int To)>();
                    foreach (int w in _net[v])
                    {
                        if (!_treeNodes.Contains(w))
                        {
                            _f.Add((v, w));
                            edgesAdded.Add((v, w));
                        }
                    }

                    var edgesRemoved = new List<(int From, int To)>();
                    foreach ((int From, int To) wv in _f)
                    {
                        if (_treeNodes.Contains(wv.From) && wv.To == v)
                        {
                            edgesRemoved.Add(wv);
                        }
                    }

                    _f.RemoveAll(edgesRemoved.Contains);

                    // LICHeE returns from this frame only; its callers keep looping but can never recurse again
                    // (the counter stays ≥ the cap), so no further tree is found — stopping here is equivalent.
                    if (_growCalls >= _maxGrowCalls)
                    {
                        _stop = true;
                        return;
                    }

                    Grow();
                    if (_stop)
                    {
                        return;
                    }

                    _f.RemoveAll(edgesAdded.Contains);
                    _f.AddRange(edgesRemoved);
                }

                TreeRemoveEdge(e.From, e.To);
                _net[e.From].Remove(e.To);
                ff.Add(e);

                // Bridge test: stop when every remaining network edge into v comes from a descendant of v in L.
                b = true;
                for (int w = 0; w < _nodeCount && b; w++)
                {
                    foreach (int n in _net[w])
                    {
                        if (n == v && (!_haveL || !IsDescendant(v, w)))
                        {
                            b = false;
                            break;
                        }
                    }
                }
            }

            for (int i = ff.Count - 1; i >= 0; i--)
            {
                _f.Add(ff[i]);
                AddNetEdge(ff[i].From, ff[i].To);
            }
        }

        private void TreeAddNode(int n)
        {
            if (!_treeNodes.Contains(n))
            {
                _treeNodes.Add(n);
            }
        }

        private void TreeAddEdge(int from, int to)
        {
            if (!_treeEdges.TryGetValue(from, out List<int>? children))
            {
                children = new List<int>();
                _treeEdges[from] = children;
            }

            if (!children.Contains(to))
            {
                children.Add(to);
            }
        }

        private void TreeRemoveEdge(int from, int to)
        {
            if (_treeEdges.TryGetValue(from, out List<int>? children))
            {
                children.Remove(to);
            }

            foreach (List<int> list in _treeEdges.Values)
            {
                if (list.Contains(to))
                {
                    return;
                }
            }

            _treeNodes.Remove(to);
        }

        private bool IsDescendant(int v, int w)
        {
            if (!_treeEdges.TryGetValue(v, out List<int>? start))
            {
                return false;
            }

            var queue = new Queue<int>(start);
            while (queue.Count > 0)
            {
                int n = queue.Dequeue();
                if (n == w)
                {
                    return true;
                }

                if (_treeEdges.TryGetValue(n, out List<int>? next))
                {
                    foreach (int c in next)
                    {
                        queue.Enqueue(c);
                    }
                }
            }

            return false;
        }

        /// <summary>LICHeE <c>PHYTree.checkConstraint</c>: sum rule (Eq. 5) at node n with margin ε.</summary>
        private bool CheckConstraint(int n)
        {
            if (!_treeEdges.TryGetValue(n, out List<int>? children))
            {
                return true;
            }

            for (int i = 0; i < _samples; i++)
            {
                double sum = 0.0;
                foreach (int c in children)
                {
                    sum += _ccf[c][i];
                }

                if (sum > _ccf[n][i] + _eps)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>LICHeE <c>PHYTree.computeErrorScore</c> (nodes in descending id order).</summary>
        private double ComputeErrorScore()
        {
            double err = 0.0;
            for (int n = _nodeCount - 1; n >= 0; n--)
            {
                if (!_treeEdges.TryGetValue(n, out List<int>? children))
                {
                    continue;
                }

                for (int i = 0; i < _samples; i++)
                {
                    double sum = 0.0;
                    foreach (int c in children)
                    {
                        sum += _ccf[c][i];
                    }

                    if (sum > _ccf[n][i])
                    {
                        double d = sum - _ccf[n][i];
                        err += d * d;
                    }
                }
            }

            return Math.Sqrt(err);
        }

        /// <summary>Parent of every input cluster (by input index) in the current complete tree; −1 = root.</summary>
        private int[] CurrentParents()
        {
            var parent = new int[_clusterTotal];
            Array.Fill(parent, int.MinValue);
            foreach (KeyValuePair<int, List<int>> kv in _treeEdges)
            {
                foreach (int child in kv.Value)
                {
                    parent[_clusterIndex[child]] = _clusterIndex[kv.Key];
                }
            }

            return parent;
        }
    }

    private static void ValidatePhylogenyTolerance(double tolerance)
    {
        if (double.IsNaN(tolerance) || tolerance < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tolerance), tolerance, "Phylogeny tolerance ε must be a non-negative number.");
        }
    }

    /// <summary>LICHeE presence profile of a cluster: '1' where CCF &gt; 0, else '0'.</summary>
    private static string PresenceProfile(IReadOnlyList<double> ccf, int samples)
    {
        var profile = new char[samples];
        for (int s = 0; s < samples; s++)
        {
            profile[s] = ccf[s] > 0.0 ? '1' : '0';
        }

        return new string(profile);
    }

    /// <summary>
    /// LICHeE node-id order of the input clusters: presence-profile groups in order of first appearance, within a
    /// group in input order (node 0 is the root).
    /// </summary>
    private static int[] LicheeNodeOrder(IReadOnlyList<CcfCluster> clusters, int samples)
    {
        var groups = new List<List<int>>();
        var groupByProfile = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        for (int c = 0; c < clusters.Count; c++)
        {
            string key = PresenceProfile(clusters[c].CcfPerSample, samples);
            if (!groupByProfile.TryGetValue(key, out List<int>? members))
            {
                members = new List<int>();
                groupByProfile[key] = members;
                groups.Add(members);
            }

            members.Add(c);
        }

        var order = new int[clusters.Count];
        int k = 0;
        foreach (List<int> group in groups)
        {
            foreach (int c in group)
            {
                order[k++] = c;
            }
        }

        return order;
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
    /// <param name="SubclonalFraction">Fraction of mutations whose CCF does not exceed the clonal threshold (CCF ≤ 0.95, i.e. not clonal under
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

        double median = StatisticsHelper.Median(values);
        if (median == 0.0)
        {
            throw new ArgumentException(
                "The median allele fraction is 0; the MATH score divides by the median and is undefined.",
                nameof(ccfDistribution));
        }

        // Raw MAD = median of absolute deviations from the median (shared helper, also behind R mad in ASPCF).
        // Same operation order as maftools mathScore.R (pat.mad = median(abs.med.dev) * 100;
        // pat.math = pat.mad * 1.4826 / median(vaf)) so the result is bit-identical to the reference.
        double rawMad = RawMedianAbsoluteDeviation(values, median);
        double percentMad = rawMad * MathPercentScale;
        return percentMad * MadConsistencyConstant / median;
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
    /// <exception cref="ArgumentException"><paramref name="ccfClusters"/> has no centroids or no assignments, or an
    /// assignment label lies outside [0, centroid count).</exception>
    public static int InferSubclones(CcfClustering ccfClusters)
    {
        IReadOnlyList<int> assignments = ccfClusters.Assignments;
        if (assignments is null || assignments.Count == 0 || ccfClusters.Centroids is null || ccfClusters.Centroids.Count == 0)
        {
            throw new ArgumentException("The CCF clustering must contain at least one cluster and one assignment.", nameof(ccfClusters));
        }

        int clusterCount = ccfClusters.Centroids.Count;
        var occupied = new HashSet<int>();
        for (int i = 0; i < assignments.Count; i++)
        {
            int label = assignments[i];
            if (label < 0 || label >= clusterCount)
            {
                throw new ArgumentException(
                    $"Assignment {i} refers to cluster {label}, outside the {clusterCount} centroid(s) [0, {clusterCount - 1}].",
                    nameof(ccfClusters));
            }

            occupied.Add(label);
        }

        return occupied.Count;
    }

    /// <summary>
    /// Performs a tumour intratumour-heterogeneity (ITH) analysis from per-mutation variant allele fractions and
    /// cancer cell fractions, returning four standard ITH metrics: the MATH score over the VAFs (Mroz &amp; Rocco
    /// 2013), the Shannon diversity index H = −Σ pᵢ·ln(pᵢ) over the clone fractions (Shannon 1948), the number of
    /// subclones (CCF clusters), and the fraction of subclonal mutations (CCF ≤ 0.95 — Landau et al. 2013: clonal if CCF &gt; 0.95, "subclonal otherwise"). CCF
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

        // Clone fractions pᵢ = proportion of mutations in each occupied cluster; Shannon H = −Σ pᵢ·ln pᵢ
        // (canonical StatisticsHelper.ShannonIndex; empty clusters contribute 0).
        int n = ccfValues.Count;
        var clusterSizes = new int[clustering.Centroids.Count];
        foreach (int label in clustering.Assignments)
        {
            clusterSizes[label]++;
        }

        double shannon = StatisticsHelper.ShannonIndex(clusterSizes);

        // Subclonal mutations = those not clonal under the canonical Landau et al. (2013) rule
        // (IdentifyClonalMutations: clonal ⇔ CCF > 0.95, "subclonal otherwise"), so CCF = 0.95 is subclonal.
        int subclonal = n - IdentifyClonalMutations(ccfValues).Count;

        double subclonalFraction = (double)subclonal / n;
        return new HeterogeneityResult(math, shannon, subcloneCount, subclonalFraction);
    }

    #endregion

}
