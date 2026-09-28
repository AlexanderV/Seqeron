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
    /// Default noise margin ε for the lineage-precedence (Eq. 2) and sum-rule (Eq. 5) inequalities. LICHeE's default
    /// (<c>-e</c>, <c>Parameters.VAF_ERROR_MARGIN</c>) is 0.1 and PICTograph uses ε₁=0.1, ε₂=0.2; this unit consumes
    /// already-clustered CCF point estimates, so the library default is the strict ε = 0 and callers pass a positive
    /// tolerance to reproduce the source defaults.
    /// </summary>
    public const double DefaultPhylogenyTolerance = 0.0;

    /// <summary>
    /// Maximum number of valid spanning trees enumerated before the search stops — LICHeE
    /// <c>Parameters.MAX_NUM_TREES</c> = 100 000 (github.com/viq854/lichee).
    /// </summary>
    public const int MaxPhylogenyTreesEnumerated = 100_000;

    /// <summary>
    /// Maximum number of recursive <c>grow</c> calls of the Gabow–Myers spanning-tree enumeration — LICHeE
    /// <c>Parameters.MAX_NUM_GROW_CALLS</c> = 10⁸. When the budget is exhausted the enumeration stops and the best tree
    /// found so far is returned.
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
    /// <param name="tolerance">Noise margin ε for both inequalities; default <see cref="DefaultPhylogenyTolerance"/> (0).</param>
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
        if (double.IsNaN(tolerance) || tolerance < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tolerance), tolerance, "Phylogeny tolerance ε must be a non-negative number.");
        }

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

        bool usedComplete = false;
        LicheeSearchResult result = new LicheeNetwork(clusters, sampleCount, tolerance, completeNetwork: false).Search();
        if (result.TreeCount == 0)
        {
            usedComplete = true;
            result = new LicheeNetwork(clusters, sampleCount, tolerance, completeNetwork: true).Search();
        }

        if (result.TreeCount == 0)
        {
            phylogeny = default;
            return false;
        }

        var edges = new ClonalEdge[clusters.Count];
        for (int c = 0; c < clusters.Count; c++)
        {
            int parentCluster = result.ParentClusterIndex[c];
            edges[c] = new ClonalEdge(parentCluster < 0 ? rootId : clusters[parentCluster].Id, clusters[c].Id);
        }

        phylogeny = new ClonalPhylogeny(rootId, clusters.ToArray(), edges, sampleCount)
        {
            Tolerance = tolerance,
            ErrorScore = result.ErrorScore,
            ValidTreeCount = result.TreeCount,
            UsedCompleteNetwork = usedComplete,
        };
        return true;
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
        private readonly double[][] _ccf;
        private readonly int[] _level;
        private readonly int[] _clusterIndex;
        private readonly SortedDictionary<int, List<int>> _levels = new();
        private readonly List<int>[] _net;

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

        public LicheeNetwork(IReadOnlyList<CcfCluster> clusters, int samples, double eps, bool completeNetwork)
        {
            _samples = samples;
            _eps = eps;
            _nodeCount = clusters.Count + 1;
            _ccf = new double[_nodeCount][];
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

            // Group clusters by presence profile in order of first appearance.
            var groups = new List<List<int>>();
            var groupByProfile = new Dictionary<string, List<int>>(StringComparer.Ordinal);
            for (int c = 0; c < clusters.Count; c++)
            {
                var profile = new char[samples];
                for (int s = 0; s < samples; s++)
                {
                    profile[s] = clusters[c].CcfPerSample[s] > 0.0 ? '1' : '0';
                }

                string key = new(profile);
                if (!groupByProfile.TryGetValue(key, out List<int>? members))
                {
                    members = new List<int>();
                    groupByProfile[key] = members;
                    groups.Add(members);
                }

                members.Add(c);
            }

            int nextId = 1;
            foreach (List<int> group in groups)
            {
                int first = nextId;
                foreach (int c in group)
                {
                    int id = nextId++;
                    _ccf[id] = clusters[c].CcfPerSample.ToArray();
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

                // Edges between each group's sub-population nodes.
                for (int i = first; i < nextId; i++)
                {
                    for (int j = i + 1; j < nextId; j++)
                    {
                        CheckAndAddEdge(i, j);
                    }
                }
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

                comp12 += a1[i] >= a2[i] - _eps ? 1 : 0;
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

                comp21 += a2[i] >= a1[i] - _eps ? 1 : 0;
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

                    if (_growCalls >= MaxPhylogenyGrowCalls)
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
            var parent = new int[_nodeCount - 1];
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
