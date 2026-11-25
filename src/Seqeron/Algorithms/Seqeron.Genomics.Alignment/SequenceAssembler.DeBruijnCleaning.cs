using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAssembler
{
    /// <summary>
    /// Velvet graph cleaning on the unitig graph of a k-mer set (Zerbino &amp; Birney 2008, Genome Res 18:821;
    /// Velvet 1.2.10 <c>velvetg</c> <c>run2.c</c>). Order of the Velvet pipeline, minus the bubble
    /// removal (<c>correctGraph</c> / Tour Bus, not implemented here):
    /// <list type="number">
    /// <item><c>clipTips_pg</c> (<c>concatenatedPreGraph.c</c>) — tips of fewer than 2k k-mers whose arc is a
    /// minority at the junction;</item>
    /// <item><c>clipTipsHard</c> (<c>correctedGraph.c</c>, the end of <c>correctGraph</c>) — dead-end chains of
    /// fewer than 2k k-mers, regardless of multiplicity;</item>
    /// <item><c>removeLowCoverageNodesAndDenounceDubiousReads</c> (<c>graphStats.c</c>, <c>-cov_cutoff</c>) —
    /// two passes (the second after re-concatenation), a node goes when the <i>integer</i> quotient
    /// total k-mer coverage / length (C <c>Coordinate</c> division) is below the cut-off;
    /// <c>-cov_cutoff auto</c> = <c>estimated_cov / 2</c>, the length-weighted median node coverage of the
    /// graph after the tip stages (<c>estimated_cov</c>, <c>graphStats.c</c>);</item>
    /// <item><c>clipTipsHard</c> again.</item>
    /// </list>
    /// Nodes are the unitigs (maximal non-branching paths, length = number of k-mers); an arc A → B exists when
    /// A's last k-mer overlaps B's first k-mer by k-1 characters. Arc multiplicity = number of (k+1)-mer
    /// occurrences in the reads (both strands when double-stranded), node coverage = number of k-mer
    /// occurrences. A removed node takes its reverse-complement twin with it.
    /// </summary>
    private static HashSet<string> CleanDeBruijnKmers(
        HashSet<string> kmers,
        int k,
        bool doubleStranded,
        Dictionary<string, int> kmerCounts,
        Dictionary<string, int> arcCounts,
        bool clipTips,
        long tipCutoff,
        double coverageCutoff,
        bool autoCoverageCutoff)
    {
        var graph = new UnitigGraph(kmers, k, doubleStranded, kmerCounts, arcCounts);
        if (clipTips)
        {
            graph.ClipTipsMultiplicity(tipCutoff);
            graph.ClipTipsHard(tipCutoff);
            if (graph.AnyRemoved)
                graph = Rebuild(graph);
        }

        double cutoff = coverageCutoff;
        if (cutoff <= 0 && autoCoverageCutoff)
            cutoff = graph.EstimatedCoverage() / 2;

        if (cutoff > 0)
        {
            graph.RemoveLowCoverage(cutoff);
            graph = Rebuild(graph);
            graph.RemoveLowCoverage(cutoff);
            if (clipTips)
            {
                graph = Rebuild(graph);
                graph.ClipTipsHard(tipCutoff);
            }
        }

        return graph.RemainingKmers();

        // Velvet concatenateGraph: re-merge the nodes that became 1-in-1-out.
        UnitigGraph Rebuild(UnitigGraph g) =>
            g.AnyRemoved ? new UnitigGraph(g.RemainingKmers(), k, doubleStranded, kmerCounts, arcCounts) : g;
    }

    /// <summary>Mutable unitig graph with the Velvet tip / coverage operations.</summary>
    private sealed class UnitigGraph
    {
        private readonly int _k;
        private readonly bool _doubleStranded;
        private readonly Dictionary<string, int> _kmerCounts;
        private readonly Dictionary<string, int> _arcCounts;
        private readonly string[] _spelled;
        private readonly int[] _length;
        private readonly long[] _coverage;
        private readonly int[] _twin;
        private readonly string[] _startVertex;
        private readonly string[] _endVertex;
        private readonly Dictionary<string, List<int>> _startingAt = new(StringComparer.Ordinal);
        private readonly Dictionary<string, List<int>> _endingAt = new(StringComparer.Ordinal);
        private readonly bool[] _removed;

        public UnitigGraph(
            HashSet<string> kmers, int k, bool doubleStranded,
            Dictionary<string, int> kmerCounts, Dictionary<string, int> arcCounts)
        {
            _k = k;
            _doubleStranded = doubleStranded;
            _kmerCounts = kmerCounts;
            _arcCounts = arcCounts;

            // Unitigs in ordinal order of their spelling: the (deterministic) node numbering.
            var kmerGraph = BuildKmerGraph(kmers);
            IEnumerable<IReadOnlyList<string>> paths = FindMaximalNonBranchingPaths(kmerGraph);
            if (doubleStranded)
                paths = SplitAtHairpins(paths, kmerGraph);
            var nodes = paths
                .Select(p => (Path: p, Spelled: SpellPath(p)))
                .OrderBy(x => x.Spelled, StringComparer.Ordinal)
                .ToArray();

            int n = nodes.Length;
            _spelled = new string[n];
            _length = new int[n];
            _coverage = new long[n];
            _twin = new int[n];
            _startVertex = new string[n];
            _endVertex = new string[n];
            _removed = new bool[n];

            var nodeOfKmer = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int i = 0; i < n; i++)
            {
                _spelled[i] = nodes[i].Spelled;
                _length[i] = _spelled[i].Length - k + 1;
                _startVertex[i] = nodes[i].Path[0];
                _endVertex[i] = nodes[i].Path[^1];
                Index(_startingAt, _startVertex[i], i);
                Index(_endingAt, _endVertex[i], i);
                foreach (ReadOnlySpan<char> kmer in _spelled[i].AsSpan().EnumerateKmers(k))
                {
                    string key = kmer.ToString();
                    nodeOfKmer[key] = i;
                    _coverage[i] += KmerCoverage(key);
                }
            }

            for (int i = 0; i < n; i++)
            {
                _twin[i] = !doubleStranded
                    ? -1
                    : nodeOfKmer[ReverseComplement(_spelled[i][..k])];
            }
        }

        /// <summary>
        /// Velvet never concatenates across a hairpin arc (<c>concatenatePreGraph_pg</c>: <c>isLoop_pg</c> is true
        /// for an arc between a node and its own twin), so a unitig that is its own reverse complement — it passes
        /// through a self-complementary (k-1)-mer, possible for odd k — is two twin nodes joined by a hairpin arc.
        /// Paths are cut at their internal self-complementary vertices; an isolated cycle through such a vertex is
        /// first rotated to start at it.
        /// </summary>
        private static List<IReadOnlyList<string>> SplitAtHairpins(
            IEnumerable<IReadOnlyList<string>> paths, Dictionary<string, List<string>> graph)
        {
            var inDegree = new Dictionary<string, int>(StringComparer.Ordinal);
            foreach (var successors in graph.Values)
            {
                foreach (string w in successors)
                    inDegree[w] = inDegree.GetValueOrDefault(w) + 1;
            }

            static bool IsSelfComplementary(string vertex) => ReverseComplement(vertex) == vertex;

            var result = new List<IReadOnlyList<string>>();
            foreach (var path in paths)
            {
                IReadOnlyList<string> q = path;
                bool isolatedCycle = path[0] == path[^1] && inDegree[path[0]] == 1 && graph[path[0]].Count == 1;
                if (isolatedCycle)
                {
                    int start = -1;
                    for (int j = 0; j < path.Count - 1 && start < 0; j++)
                    {
                        if (IsSelfComplementary(path[j]))
                            start = j;
                    }

                    if (start > 0)
                    {
                        var ring = path.Skip(start).Take(path.Count - 1 - start).Concat(path.Take(start)).ToList();
                        ring.Add(ring[0]);
                        q = ring;
                    }
                }

                int from = 0;
                for (int j = 1; j < q.Count - 1; j++)
                {
                    if (!IsSelfComplementary(q[j]))
                        continue;
                    result.Add(q.Skip(from).Take(j - from + 1).ToList());
                    from = j;
                }

                result.Add(from == 0 ? q : q.Skip(from).ToList());
            }

            return result;
        }

        public bool AnyRemoved => _removed.Contains(true);

        private static void Index(Dictionary<string, List<int>> map, string vertex, int node)
        {
            if (!map.TryGetValue(vertex, out var list))
                map[vertex] = list = new List<int>();
            list.Add(node);
        }

        /// <summary>Velvet node coverage unit: occurrences of the k-mer in the reads, both strands for a twin node.</summary>
        private long KmerCoverage(string kmer)
        {
            long c = _kmerCounts.GetValueOrDefault(kmer);
            if (_doubleStranded)
            {
                string rc = ReverseComplement(kmer);
                if (rc != kmer)
                    c += _kmerCounts.GetValueOrDefault(rc);
            }

            return c;
        }

        /// <summary>Multiplicity of the arc <paramref name="pred"/> → <paramref name="succ"/> = occurrences of its (k+1)-mer.</summary>
        private long ArcMultiplicity(int pred, int succ)
        {
            string w = string.Concat(_spelled[pred].AsSpan(_spelled[pred].Length - _k), _spelled[succ].AsSpan(_k - 1, 1));
            long m = _arcCounts.GetValueOrDefault(w);
            if (_doubleStranded)
                m += _arcCounts.GetValueOrDefault(ReverseComplement(w));
            return m;
        }

        /// <summary>Live neighbours of <paramref name="node"/> through its end (forward) or start (backward).</summary>
        private IEnumerable<int> Arcs(int node, bool forward)
        {
            var map = forward ? _startingAt : _endingAt;
            string vertex = forward ? _endVertex[node] : _startVertex[node];
            if (!map.TryGetValue(vertex, out var list))
                yield break;
            foreach (int other in list)
            {
                if (!_removed[other])
                    yield return other;
            }
        }

        private int ArcCount(int node, bool forward) => Arcs(node, forward).Count();

        /// <summary>Velvet numbers a node pair once and tries both orientations; the larger strand is the same node.</summary>
        private bool IsSecondOfPair(int node) => _twin[node] >= 0 && _twin[node] < node;

        private void Remove(int node)
        {
            _removed[node] = true;
            if (_twin[node] >= 0)
                _removed[_twin[node]] = true;
        }

        /// <summary>
        /// Velvet <c>isEligibleTip</c> (<c>concatenatedPreGraph.c</c>) for node <paramref name="x"/> whose
        /// <paramref name="deadForward"/> side (end or start) has no arc: walk from the dead end through simple
        /// nodes to the first junction; a tip when the chain is shorter than the cut-off and it ends in an
        /// isolated snippet, or at a junction where the chain's arc has multiplicity 1 or is a strict minority.
        /// </summary>
        private bool IsEligibleTip(int x, bool deadForward, long cutoff)
        {
            if (ArcCount(x, deadForward) != 0)
                return false;

            bool walkForward = !deadForward;
            long total = 0;
            int cur = x, last = -1, steps = 0;
            while (cur >= 0 && ArcCount(cur, deadForward) < 2 && ArcCount(cur, walkForward) < 2)
            {
                if (++steps > _length.Length)
                    return false;
                total += _length[cur];
                last = cur;
                cur = Arcs(cur, walkForward).DefaultIfEmpty(-1).First();
            }

            if (total >= cutoff)
                return false;
            if (cur < 0)
                return true;                                  // isolated snippet
            if (ArcCount(cur, deadForward) < 2)
                return false;                                 // joined tips

            long mult = Multiplicity(last, cur);
            if (mult == 1)
                return true;                                  // unique event
            long max = Arcs(cur, deadForward).Max(p => Multiplicity(p, cur));
            return max > mult;                                // minority

            // Arc between a walk-previous node and a walk-next node, in the graph's forward orientation.
            long Multiplicity(int walkPrevious, int walkNext) =>
                walkForward ? ArcMultiplicity(walkPrevious, walkNext) : ArcMultiplicity(walkNext, walkPrevious);
        }

        /// <summary>Velvet <c>getTipLength</c> (<c>correctedGraph.c</c>), for a node with no arc on its dead side.</summary>
        private long HardTipLength(int x, bool deadForward)
        {
            bool walkForward = !deadForward;
            if (ArcCount(x, walkForward) > 1)
                return _length[x];

            long length = 0;
            int cur = x, steps = 0;
            while (cur >= 0 && ArcCount(cur, deadForward) < 2 && ArcCount(cur, walkForward) < 2)
            {
                if (++steps > _length.Length)
                    break;
                length += _length[cur];
                cur = Arcs(cur, walkForward).DefaultIfEmpty(-1).First();
            }

            return length;
        }

        /// <summary>Velvet <c>clipTips_pg</c>: repeat until no tip is removed.</summary>
        public void ClipTipsMultiplicity(long cutoff)
        {
            bool modified = true;
            while (modified)
            {
                modified = false;
                for (int i = 0; i < _removed.Length; i++)
                {
                    if (_removed[i] || IsSecondOfPair(i))
                        continue;
                    if (IsEligibleTip(i, true, cutoff) || IsEligibleTip(i, false, cutoff))
                    {
                        Remove(i);
                        modified = true;
                    }
                }
            }
        }

        /// <summary>Velvet <c>clipTipsHard</c>: a dead-end chain shorter than the cut-off goes, whatever its multiplicity.</summary>
        public void ClipTipsHard(long cutoff)
        {
            bool modified = true;
            while (modified)
            {
                modified = false;
                for (int i = 0; i < _removed.Length; i++)
                {
                    if (_removed[i] || IsSecondOfPair(i))
                        continue;
                    if ((ArcCount(i, true) == 0 && HardTipLength(i, true) < cutoff)
                        || (ArcCount(i, false) == 0 && HardTipLength(i, false) < cutoff))
                    {
                        Remove(i);
                        modified = true;
                    }
                }
            }
        }

        /// <summary>
        /// One pass of Velvet <c>removeLowCoverageNodesAndDenounceDubiousReads</c>: remove every node whose
        /// <c>getTotalCoverage / getNodeLength</c> — an integer quotient in C — is below the cut-off.
        /// </summary>
        public void RemoveLowCoverage(double cutoff)
        {
            for (int i = 0; i < _removed.Length; i++)
            {
                if (!_removed[i] && _coverage[i] / _length[i] < cutoff)
                    Remove(i);
            }
        }

        /// <summary>
        /// Velvet <c>estimated_cov</c>: the coverage (total / length, floating point) of the node at which the
        /// cumulative length of the coverage-sorted nodes first reaches half the total length (integer
        /// halving); each strand pair counts once. 0 for an empty graph.
        /// </summary>
        public double EstimatedCoverage()
        {
            var live = Enumerable.Range(0, _removed.Length)
                .Where(i => !_removed[i] && (_twin[i] < 0 || i <= _twin[i]))
                .OrderBy(i => (double)_coverage[i] / _length[i])
                .ToList();
            if (live.Count == 0)
                return 0;

            long half = live.Sum(i => (long)_length[i]) / 2;
            long sum = 0;
            foreach (int i in live)
            {
                sum += _length[i];
                if (sum >= half)
                    return (double)_coverage[i] / _length[i];
            }

            return -1;
        }

        public HashSet<string> RemainingKmers()
        {
            var result = new HashSet<string>(StringComparer.Ordinal);
            for (int i = 0; i < _removed.Length; i++)
            {
                if (_removed[i])
                    continue;
                foreach (ReadOnlySpan<char> kmer in _spelled[i].AsSpan().EnumerateKmers(_k))
                    result.Add(kmer.ToString());
            }

            return result;
        }
    }
}
