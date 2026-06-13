using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Provides genome assembly algorithms for constructing contigs from sequence reads.
/// Supports overlap-layout-consensus and de Bruijn graph approaches.
/// </summary>
public static class SequenceAssembler
{
    /// <summary>
    /// Result of sequence assembly.
    /// </summary>
    public readonly record struct AssemblyResult(
        IReadOnlyList<string> Contigs,
        int TotalReads,
        int AssembledReads,
        double N50,
        int LongestContig,
        int TotalLength);

    /// <summary>
    /// Represents an overlap between two reads.
    /// </summary>
    public readonly record struct Overlap(
        int ReadIndex1,
        int ReadIndex2,
        int OverlapLength,
        int Position1,
        int Position2);

    /// <summary>
    /// Assembly parameters.
    /// </summary>
    public readonly record struct AssemblyParameters(
        int MinOverlap = 20,
        double MinIdentity = 0.9,
        int KmerSize = 31,
        int MinContigLength = 100);

    /// <summary>
    /// Assembles reads using the Overlap-Layout-Consensus (OLC) paradigm:
    /// (1) Overlap — build the overlap graph by finding the longest suffix-prefix
    /// overlap (length ≥ <see cref="AssemblyParameters.MinOverlap"/>) between every
    /// ordered pair of reads; (2) Layout — greedily chain reads by their best
    /// (longest) overlap into contigs; (3) Consensus — emit the merged superstring
    /// of each chain. Exact OLC layout (a Hamiltonian path through the overlap graph)
    /// is NP-complete, so a greedy heuristic is used; it reconstructs unambiguous
    /// (non-repeat) tilings but, like all OLC heuristics, may split or not optimally
    /// resolve repeats longer than the read length.
    /// </summary>
    /// <remarks>
    /// References: Compeau, Pevzner &amp; Tesler (2011), Nat Biotechnol 29:987–991,
    /// DOI 10.1038/nbt.2023 (overlap graph, Hamiltonian-path layout, NP-completeness);
    /// Langmead, "Overlap Layout Consensus assembly" (JHU lecture notes), p.4–5, 25, 28.
    /// </remarks>
    public static AssemblyResult AssembleOLC(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
    {
        var param = parameters ?? new AssemblyParameters();

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        // Step 1: Find overlaps
        var overlaps = FindAllOverlaps(reads, param.MinOverlap, param.MinIdentity);

        // Step 2: Build overlap graph and find layout
        var contigs = BuildContigsFromOverlaps(reads, overlaps, param);

        // Filter by minimum length
        contigs = contigs.Where(c => c.Length >= param.MinContigLength).ToList();

        // Calculate statistics
        var stats = CalculateStats(contigs, reads.Count);

        return stats;
    }

    /// <summary>
    /// Assembles reads using the de Bruijn graph (DBG) paradigm. The reads are decomposed
    /// into k-mers; each k-mer is a directed edge from its (k-1)-prefix to its (k-1)-suffix
    /// in a directed multigraph whose nodes are (k-1)-mers (see <see cref="BuildDeBruijnGraph"/>).
    /// Each weakly-connected component is then spelled out by an Eulerian walk — a walk that
    /// traverses every edge exactly once — and the resulting superstring is emitted as one
    /// contig. Reconstruction is exact when a component's Eulerian walk is unique (no repeated
    /// (k-1)-mer); when a (k-1)-mer repeats, several Eulerian walks exist and only one is the
    /// true genome, so the emitted contig may not be the original sequence.
    /// </summary>
    /// <remarks>
    /// References: Langmead B., "De Bruijn Graph assembly" (JHU lecture notes), p.5-11
    /// (nodes = (k-1)-mers, edges = k-mers), p.15-19 (Eulerian walk spells the genome,
    /// <c>path[0] + concat(last char of each subsequent node)</c>); Jones &amp; Pevzner (2004),
    /// <i>An Introduction to Bioinformatics Algorithms</i>, MIT Press, Theorems 8.1/8.2
    /// (a connected graph has an Eulerian path iff it has at most two semi-balanced nodes,
    /// all others balanced); Compeau, Pevzner &amp; Tesler (2011), Nat Biotechnol 29:987-991,
    /// DOI 10.1038/nbt.2023.
    /// </remarks>
    public static AssemblyResult AssembleDeBruijn(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
    {
        var param = parameters ?? new AssemblyParameters();

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        int k = param.KmerSize;

        // Build the (k-1)-mer de Bruijn multigraph (edges = input k-mers).
        var graph = BuildDeBruijnGraph(reads, k);

        // One Eulerian walk per weakly-connected component -> one contig per component.
        var contigs = ReconstructContigs(graph);

        // Filter by minimum contig length.
        contigs = contigs.Where(c => c.Length >= param.MinContigLength).ToList();

        return CalculateStats(contigs, reads.Count);
    }

    /// <summary>
    /// Builds the directed overlap graph for a set of reads: for every ordered pair
    /// (i, j), i ≠ j, reports the longest suffix-of-read[i] / prefix-of-read[j] overlap
    /// whose length is ≥ <paramref name="minOverlap"/> and whose identity is
    /// ≥ <paramref name="minIdentity"/>. Each returned <see cref="Overlap"/> is an edge
    /// i → j with weight <see cref="Overlap.OverlapLength"/>. Self-overlaps are excluded.
    /// </summary>
    /// <remarks>
    /// Overlap definition and "report only the longest suffix/prefix match" per
    /// Langmead, "Overlap Layout Consensus assembly" (JHU lecture notes), p.5, p.10;
    /// overlap-graph edge semantics per Compeau, Pevzner &amp; Tesler (2011),
    /// Nat Biotechnol 29:987–991, DOI 10.1038/nbt.2023. All-pairs scan is O(N²); chosen
    /// over the O(N+a) suffix-tree approach because <paramref name="minIdentity"/> permits
    /// mismatches, which the exact-match suffix tree cannot model (Langmead OLC p.16).
    /// </remarks>
    public static IReadOnlyList<Overlap> FindAllOverlaps(
        IReadOnlyList<string> reads,
        int minOverlap = 20,
        double minIdentity = 0.9)
    {
        var overlaps = new List<Overlap>();

        for (int i = 0; i < reads.Count; i++)
        {
            for (int j = 0; j < reads.Count; j++)
            {
                if (i == j) continue;

                var overlap = FindOverlap(reads[i], reads[j], minOverlap, minIdentity);
                if (overlap.HasValue)
                {
                    overlaps.Add(new Overlap(i, j, overlap.Value.length,
                        overlap.Value.pos1, overlap.Value.pos2));
                }
            }
        }

        return overlaps;
    }

    /// <summary>
    /// Finds all overlaps between reads with cancellation support.
    /// Useful for large read sets.
    /// </summary>
    /// <param name="reads">Collection of sequence reads.</param>
    /// <param name="minOverlap">Minimum overlap length.</param>
    /// <param name="minIdentity">Minimum identity threshold (0.0-1.0).</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <param name="progress">Optional progress reporter (0.0 to 1.0).</param>
    /// <returns>List of overlaps found.</returns>
    public static IReadOnlyList<Overlap> FindAllOverlaps(
        IReadOnlyList<string> reads,
        int minOverlap,
        double minIdentity,
        CancellationToken cancellationToken,
        IProgress<double>? progress = null)
    {
        var overlaps = new List<Overlap>();
        int total = reads.Count * reads.Count;
        int processed = 0;

        for (int i = 0; i < reads.Count; i++)
        {
            for (int j = 0; j < reads.Count; j++)
            {
                if (processed % 100 == 0)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    progress?.Report((double)processed / total);
                }
                processed++;

                if (i == j) continue;

                var overlap = FindOverlap(reads[i], reads[j], minOverlap, minIdentity);
                if (overlap.HasValue)
                {
                    overlaps.Add(new Overlap(
                        i, j, overlap.Value.length, overlap.Value.pos1, overlap.Value.pos2));
                }
            }
        }

        progress?.Report(1.0);
        return overlaps;
    }

    /// <summary>
    /// Finds the longest suffix-of-<paramref name="seq1"/> / prefix-of-<paramref name="seq2"/>
    /// overlap whose length is ≥ <paramref name="minOverlap"/> and whose identity is
    /// ≥ <paramref name="minIdentity"/>; returns the overlap length and the 0-based start
    /// positions (<c>pos1</c> in seq1, <c>pos2</c> = 0 in seq2), or <c>null</c> if none.
    /// Only the single longest qualifying overlap is reported (Langmead OLC p.5, p.10).
    /// </summary>
    public static (int length, int pos1, int pos2)? FindOverlap(
        string seq1, string seq2,
        int minOverlap = 20,
        double minIdentity = 0.9)
    {
        // Check if suffix of seq1 overlaps with prefix of seq2
        int maxPossible = Math.Min(seq1.Length, seq2.Length);

        for (int overlapLen = maxPossible; overlapLen >= minOverlap; overlapLen--)
        {
            string suffix = seq1.Substring(seq1.Length - overlapLen);
            string prefix = seq2.Substring(0, overlapLen);

            double identity = CalculateIdentity(suffix, prefix);
            if (identity >= minIdentity)
            {
                return (overlapLen, seq1.Length - overlapLen, 0);
            }
        }

        return null;
    }

    /// <summary>
    /// Calculates sequence identity between two strings of equal length.
    /// </summary>
    public static double CalculateIdentity(string seq1, string seq2)
    {
        if (seq1.Length != seq2.Length) return 0;
        if (seq1.Length == 0) return 1;

        int matches = 0;
        for (int i = 0; i < seq1.Length; i++)
        {
            if (char.ToUpperInvariant(seq1[i]) == char.ToUpperInvariant(seq2[i]))
                matches++;
        }

        return (double)matches / seq1.Length;
    }

    private static List<string> BuildContigsFromOverlaps(
        IReadOnlyList<string> reads,
        IReadOnlyList<Overlap> overlaps,
        AssemblyParameters param)
    {
        var contigs = new List<string>();
        var used = new HashSet<int>();

        // Build adjacency: for each read, best successor
        var bestSuccessor = new Dictionary<int, (int next, int overlap)>();
        var hasPredecessor = new HashSet<int>();

        foreach (var ov in overlaps.OrderByDescending(o => o.OverlapLength))
        {
            int r1 = ov.ReadIndex1;
            int r2 = ov.ReadIndex2;

            if (!bestSuccessor.ContainsKey(r1))
            {
                bestSuccessor[r1] = (r2, ov.OverlapLength);
                hasPredecessor.Add(r2);
            }
        }

        // Find starting reads (no predecessor in best-overlap chain)
        var starters = new List<int>();
        for (int i = 0; i < reads.Count; i++)
        {
            if (!hasPredecessor.Contains(i))
                starters.Add(i);
        }

        // Build contigs from each starter
        foreach (int start in starters)
        {
            if (used.Contains(start)) continue;

            var sb = new StringBuilder();
            int current = start;

            while (current != -1 && !used.Contains(current))
            {
                used.Add(current);
                string read = reads[current];

                if (sb.Length == 0)
                {
                    sb.Append(read);
                }
                else if (bestSuccessor.TryGetValue(current, out var info))
                {
                    // Already added previous, now extend
                }

                if (bestSuccessor.TryGetValue(current, out var next))
                {
                    int overlap = next.overlap;
                    string nextRead = reads[next.next];
                    if (!used.Contains(next.next))
                    {
                        sb.Append(nextRead.Substring(overlap));
                    }
                    current = next.next;
                }
                else
                {
                    current = -1;
                }
            }

            if (sb.Length > 0)
                contigs.Add(sb.ToString());
        }

        // Add unused reads as singleton contigs
        for (int i = 0; i < reads.Count; i++)
        {
            if (!used.Contains(i) && reads[i].Length >= param.MinContigLength)
            {
                contigs.Add(reads[i]);
            }
        }

        return contigs;
    }

    /// <summary>
    /// Builds the de Bruijn multigraph for a set of reads at k-mer length
    /// <paramref name="k"/>. Each input k-mer becomes one directed edge from its left
    /// (k-1)-mer (prefix, <c>kmer[0..k-1]</c>) to its right (k-1)-mer (suffix,
    /// <c>kmer[1..k]</c>); the graph nodes are the distinct (k-1)-mers. The result is the
    /// out-adjacency multimap: <c>graph[u]</c> lists the suffix node for every k-mer whose
    /// prefix is <c>u</c>, so a k-mer that occurs more than once yields a repeated
    /// (multi-)edge. Reads shorter than <paramref name="k"/> contribute no k-mers and are
    /// silently skipped.
    /// </summary>
    /// <param name="reads">Input reads (k-mers are read left-to-right, no reverse-complement).</param>
    /// <param name="k">k-mer length; must be ≥ 2 so that the (k-1)-mer nodes are non-empty.</param>
    /// <returns>Out-adjacency multimap from each (k-1)-mer node to its successor nodes.</returns>
    /// <exception cref="ArgumentOutOfRangeException">If <paramref name="k"/> &lt; 2.</exception>
    /// <remarks>
    /// Node/edge definition and the prefix/suffix split per Langmead B., "De Bruijn Graph
    /// assembly" (JHU lecture notes), p.5-9, p.16; multiedges per the same notes p.8.
    /// </remarks>
    public static Dictionary<string, List<string>> BuildDeBruijnGraph(
        IReadOnlyList<string> reads, int k)
    {
        // k-mers split into two (k-1)-mers, so k must be at least 2 (Langmead DBG p.5).
        const int MinKmerLength = 2;
        if (k < MinKmerLength)
            throw new ArgumentOutOfRangeException(nameof(k), k,
                $"k-mer length must be at least {MinKmerLength} so that (k-1)-mer nodes are non-empty.");

        var graph = new Dictionary<string, List<string>>();

        if (reads == null)
            return graph;

        foreach (string read in reads)
        {
            if (read == null || read.Length < k)
                continue; // A read shorter than k contributes no k-mers (Langmead DBG p.16).

            // Chop the read into k-mers: i in [0, read.Length - k].
            for (int i = 0; i <= read.Length - k; i++)
            {
                string prefix = read.Substring(i, k - 1);     // left (k-1)-mer
                string suffix = read.Substring(i + 1, k - 1); // right (k-1)-mer

                if (!graph.TryGetValue(prefix, out var neighbors))
                {
                    neighbors = new List<string>();
                    graph[prefix] = neighbors;
                }

                neighbors.Add(suffix);
            }
        }

        return graph;
    }

    /// <summary>
    /// Reconstructs contigs from a de Bruijn multigraph by computing one Eulerian walk per
    /// weakly-connected component and spelling each walk out into a superstring. Each
    /// component contributes one contig: the first node is emitted in full and the last
    /// character of every subsequent node is appended (Langmead DBG p.18). A walk starts at a
    /// semi-balanced source node (out-degree − in-degree = 1) when one exists, otherwise at any
    /// node of the component (an Eulerian cycle); this realizes the Eulerian-path existence
    /// condition of Jones &amp; Pevzner Theorems 8.1/8.2.
    /// </summary>
    private static List<string> ReconstructContigs(Dictionary<string, List<string>> graph)
    {
        var contigs = new List<string>();
        if (graph.Count == 0)
            return contigs;

        // Mutable copy of the out-adjacency (edges are consumed by the walk).
        var adjacency = new Dictionary<string, List<string>>();
        var inDegree = new Dictionary<string, int>();
        var outDegree = new Dictionary<string, int>();
        var allNodes = new HashSet<string>();

        foreach (var kvp in graph)
        {
            adjacency[kvp.Key] = new List<string>(kvp.Value);
            allNodes.Add(kvp.Key);
            outDegree[kvp.Key] = outDegree.GetValueOrDefault(kvp.Key, 0) + kvp.Value.Count;
            foreach (string target in kvp.Value)
            {
                allNodes.Add(target);
                inDegree[target] = inDegree.GetValueOrDefault(target, 0) + 1;
            }
        }

        // Undirected adjacency for connectivity (weakly-connected components).
        var undirected = BuildUndirectedAdjacency(graph, allNodes);

        var visitedNodes = new HashSet<string>();
        // Deterministic component order: by lexicographically smallest node.
        foreach (string node in allNodes.OrderBy(n => n, StringComparer.Ordinal))
        {
            if (visitedNodes.Contains(node))
                continue;

            var component = CollectComponent(node, undirected, visitedNodes);
            string start = ChooseEulerStart(component, inDegree, outDegree);
            string contig = SpellEulerianWalk(adjacency, start);
            if (contig.Length > 0)
                contigs.Add(contig);
        }

        return contigs;
    }

    private static Dictionary<string, List<string>> BuildUndirectedAdjacency(
        Dictionary<string, List<string>> graph, HashSet<string> allNodes)
    {
        var undirected = new Dictionary<string, List<string>>();
        foreach (string n in allNodes)
            undirected[n] = new List<string>();

        foreach (var kvp in graph)
        {
            foreach (string target in kvp.Value)
            {
                undirected[kvp.Key].Add(target);
                undirected[target].Add(kvp.Key);
            }
        }

        return undirected;
    }

    private static List<string> CollectComponent(
        string start,
        Dictionary<string, List<string>> undirected,
        HashSet<string> visitedNodes)
    {
        var component = new List<string>();
        var stack = new Stack<string>();
        stack.Push(start);
        visitedNodes.Add(start);

        while (stack.Count > 0)
        {
            string current = stack.Pop();
            component.Add(current);
            foreach (string neighbor in undirected[current])
            {
                if (visitedNodes.Add(neighbor))
                    stack.Push(neighbor);
            }
        }

        return component;
    }

    /// <summary>
    /// Picks the Eulerian-walk start node for one component: the unique semi-balanced node
    /// with one more outgoing than incoming edge (Langmead DBG p.15; J&amp;P Theorem 8.2),
    /// or — for a balanced component (Eulerian cycle) or any imperfect graph — the
    /// lexicographically smallest node, for determinism.
    /// </summary>
    private static string ChooseEulerStart(
        List<string> component,
        Dictionary<string, int> inDegree,
        Dictionary<string, int> outDegree)
    {
        string fallback = component[0];
        foreach (string node in component)
        {
            if (StringComparer.Ordinal.Compare(node, fallback) < 0)
                fallback = node;

            int outD = outDegree.GetValueOrDefault(node, 0);
            int inD = inDegree.GetValueOrDefault(node, 0);
            if (outD - inD == 1)
                return node; // semi-balanced source: start of the Eulerian path
        }

        return fallback;
    }

    /// <summary>
    /// Hierholzer's algorithm: traverses each out-edge of the component exactly once starting
    /// from <paramref name="start"/>, then spells the resulting node walk into a superstring
    /// (<c>walk[0] + concat(last char of walk[i] for i &gt; 0)</c>). Runs in O(|E|) time
    /// (Langmead DBG p.17).
    /// </summary>
    private static string SpellEulerianWalk(
        Dictionary<string, List<string>> adjacency, string start)
    {
        // Index of next unused out-edge per node (edges consumed in stored order).
        var edgeCursor = new Dictionary<string, int>();
        var stack = new Stack<string>();
        var walk = new List<string>();
        stack.Push(start);

        while (stack.Count > 0)
        {
            string node = stack.Peek();
            int cursor = edgeCursor.GetValueOrDefault(node, 0);

            if (adjacency.TryGetValue(node, out var neighbors) && cursor < neighbors.Count)
            {
                edgeCursor[node] = cursor + 1;
                stack.Push(neighbors[cursor]);
            }
            else
            {
                walk.Add(node);
                stack.Pop();
            }
        }

        // Hierholzer produces the Eulerian walk in reverse order.
        walk.Reverse();

        if (walk.Count == 0)
            return string.Empty;

        var sb = new StringBuilder(walk[0]);
        for (int i = 1; i < walk.Count; i++)
            sb.Append(walk[i][^1]); // append last character of each subsequent (k-1)-mer node

        return sb.ToString();
    }

    /// <summary>
    /// Calculates N50 and other assembly statistics.
    /// </summary>
    public static AssemblyResult CalculateStats(
        IReadOnlyList<string> contigs, int totalReads)
    {
        if (contigs.Count == 0)
            return new AssemblyResult(contigs, totalReads, 0, 0, 0, 0);

        var sortedLengths = contigs.Select(c => c.Length).OrderDescending().ToList();
        int totalLength = sortedLengths.Sum();
        int halfLength = totalLength / 2;

        // Calculate N50
        int cumulative = 0;
        double n50 = 0;
        foreach (int len in sortedLengths)
        {
            cumulative += len;
            if (cumulative >= halfLength)
            {
                n50 = len;
                break;
            }
        }

        return new AssemblyResult(
            contigs,
            totalReads,
            totalReads, // Simplified: assume all reads used
            n50,
            sortedLengths.First(),
            totalLength);
    }

    /// <summary>
    /// Merges two contigs whose suffix/prefix overlap length is already known, collapsing the
    /// shared region so it appears once: the merged superstring is
    /// <paramref name="contig1"/> followed by <paramref name="contig2"/> with its length-<c>l</c>
    /// prefix removed (<c>l = overlapLength</c>). For a valid overlap the result length is
    /// <c>|contig1| + |contig2| − overlapLength</c>.
    /// </summary>
    /// <remarks>
    /// An overlap is a length-<c>l</c> suffix of <paramref name="contig1"/> that matches a
    /// length-<c>l</c> prefix of <paramref name="contig2"/>; merging keeps a single copy of that
    /// region (Langmead, "Assembly &amp; shortest common superstring", JHU notes — overlap
    /// definition and the greedy merge trace, e.g. <c>BAA</c> + <c>AAB</c> at overlap 2 → <c>BAAB</c>).
    /// A non-positive overlap, or one larger than the shorter contig, is not a usable suffix/prefix
    /// overlap, so the contigs are simply concatenated (Langmead: "without requirement of 'shortest',
    /// it's easy: just concatenate them"; a valid overlap is bounded by <c>min(|x|,|y|)</c>).
    /// The supplied <paramref name="overlapLength"/> is trusted; computing/verifying it is the
    /// responsibility of <see cref="FindOverlap"/> / <c>FindAllOverlaps</c>.
    /// </remarks>
    /// <param name="contig1">The left contig; its suffix overlaps the prefix of <paramref name="contig2"/>.</param>
    /// <param name="contig2">The right contig; its length-<paramref name="overlapLength"/> prefix is the shared region.</param>
    /// <param name="overlapLength">Length of the suffix(contig1)/prefix(contig2) overlap to collapse.</param>
    /// <returns>The merged superstring.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contig1"/> or <paramref name="contig2"/> is null.</exception>
    public static string MergeContigs(string contig1, string contig2, int overlapLength)
    {
        ArgumentNullException.ThrowIfNull(contig1);
        ArgumentNullException.ThrowIfNull(contig2);

        // Overlap length 0 means "no overlap → concatenate" and a valid overlap cannot exceed the
        // shorter contig (it is both a suffix of contig1 and a prefix of contig2): Langmead SCS notes.
        const int NoOverlap = 0;
        if (overlapLength <= NoOverlap || overlapLength > Math.Min(contig1.Length, contig2.Length))
            return contig1 + contig2;

        return contig1 + contig2.Substring(overlapLength);
    }

    /// <summary>
    /// Default gap length, in characters, used when a paired-end link reports a non-positive
    /// (zero, negative, or unknown) gap estimate. Negative estimates indicate the contigs should
    /// overlap, but resolving that overlap is out of scope here, so a placeholder gap is emitted
    /// instead. 100 is the GenBank/EMBL/DDBJ standard length for a gap of unknown size
    /// (NCBI AGP Specification v2.1, §"Gaps": "use ... 100 as the gap size, since 100 is the
    /// GenBank/EMBL/DDBJ standard for gaps of unknown size"; gap lengths must be positive).
    /// </summary>
    private const int UnknownGapLength = 100;

    /// <summary>
    /// Joins ordered contigs into scaffolds using paired-end link information. Following the
    /// scaffold construction of Jackman et al. (ABySS 2.0, Genome Research 2017), the contigs
    /// along a link path are concatenated, "interspersed with gaps represented by a run of the
    /// character <c>N</c>, whose length corresponds to the estimate of the distance between those
    /// two contigs". Each link <c>(contig1, contig2, gapSize)</c> places <paramref name="contigs"/>
    /// [contig2] immediately after [contig1] with a gap of <c>gapSize</c> copies of
    /// <paramref name="gapCharacter"/> between them.
    /// </summary>
    /// <remarks>
    /// Gap length: a positive <c>gapSize</c> emits exactly that many gap characters (Jackman et al.
    /// 2017). A non-positive estimate (<c>gapSize ≤ 0</c>) is treated as a gap of unknown size: a
    /// negative estimate indicates the contigs should overlap (Sahlin et al. 2012, "Improved gap
    /// size estimation for scaffolding algorithms"; Jackman et al. 2017: "It is possible that the
    /// distance estimate is negative, indicating that the two contigs should in fact overlap"), but
    /// overlap resolution is out of scope here, so the GenBank/EMBL/DDBJ standard unknown-gap length
    /// of <see cref="UnknownGapLength"/> characters is emitted (NCBI AGP Specification v2.1).
    /// Each contig is placed into at most one scaffold; a link to an already-placed contig is
    /// skipped (a contig cannot appear twice). Contigs not reached by any followed link start their
    /// own single-contig scaffold, in ascending index order.
    /// </remarks>
    /// <param name="contigs">The contigs to scaffold, indexed by position; links reference these indices.</param>
    /// <param name="links">Paired-end links: each <c>(contig1, contig2, gapSize)</c> orders contig2 after contig1 with the given gap estimate (in characters). Indices out of range are ignored.</param>
    /// <param name="gapCharacter">The character used to fill gaps between contigs (default <c>'N'</c>).</param>
    /// <returns>The assembled scaffolds, one string per scaffold.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contigs"/> or <paramref name="links"/> is null.</exception>
    public static IReadOnlyList<string> Scaffold(
        IReadOnlyList<string> contigs,
        IReadOnlyList<(int contig1, int contig2, int gapSize)> links,
        char gapCharacter = 'N')
    {
        ArgumentNullException.ThrowIfNull(contigs);
        ArgumentNullException.ThrowIfNull(links);

        if (contigs.Count == 0) return Array.Empty<string>();

        // Index the first usable forward link out of each contig (links to in-range, distinct
        // contigs), preserving input order so the first declared link wins on ties.
        var linkMap = new Dictionary<int, List<(int contig1, int contig2, int gapSize)>>();
        foreach (var link in links)
        {
            if (link.contig1 < 0 || link.contig1 >= contigs.Count) continue;
            if (link.contig2 < 0 || link.contig2 >= contigs.Count) continue;
            if (link.contig2 == link.contig1) continue;

            if (!linkMap.TryGetValue(link.contig1, out var bucket))
            {
                bucket = new List<(int, int, int)>();
                linkMap[link.contig1] = bucket;
            }
            bucket.Add(link);
        }

        var scaffolds = new List<string>();
        var used = new HashSet<int>();

        for (int i = 0; i < contigs.Count; i++)
        {
            if (used.Contains(i)) continue;

            var sb = new StringBuilder(contigs[i]);
            used.Add(i);

            // Follow the link path from the current contig, appending each unplaced successor
            // separated by a run of gap characters (Jackman et al. 2017 scaffold construction).
            int current = i;
            while (linkMap.TryGetValue(current, out var nextLinks))
            {
                int nextIndex = -1;
                int gapSize = 0;
                foreach (var link in nextLinks)
                {
                    if (used.Contains(link.contig2)) continue;
                    nextIndex = link.contig2;
                    gapSize = link.gapSize;
                    break;
                }

                if (nextIndex < 0) break;

                int gapLength = gapSize > 0 ? gapSize : UnknownGapLength;
                sb.Append(gapCharacter, gapLength);
                sb.Append(contigs[nextIndex]);
                used.Add(nextIndex);
                current = nextIndex;
            }

            scaffolds.Add(sb.ToString());
        }

        return scaffolds;
    }

    /// <summary>
    /// Default minimum number of matching characters required to place a read against the
    /// reference. This is a usability default for read placement, not a biological constant;
    /// the depth-counting arithmetic (per-base depth = number of placed reads spanning the
    /// position) is independent of its value (Metagenomics Wiki, "SAMtools: get breadth of
    /// coverage"; Cook, "Calculate Depth and Breadth of Coverage From a bam File").
    /// </summary>
    private const int DefaultMinOverlap = 20;

    /// <summary>
    /// Calculates per-base sequencing coverage (read depth) along the reference. Each read is
    /// placed at its best ungapped match position (requiring at least <paramref name="minOverlap"/>
    /// matching characters); a placed read of length L starting at position p increments the depth
    /// of every reference position in the half-open interval [p, min(p + L, reference.Length)).
    /// Per-base depth is the number of reads covering that position (Metagenomics Wiki, "SAMtools:
    /// get breadth of coverage"; Cook, "Calculate Depth and Breadth of Coverage From a bam File").
    /// </summary>
    /// <param name="reference">Reference sequence; depth is reported per position.</param>
    /// <param name="reads">Reads to map against the reference.</param>
    /// <param name="minOverlap">Minimum matching characters required to place a read.</param>
    /// <returns>
    /// An array of length <c>reference.Length</c> whose i-th element is the number of placed reads
    /// covering reference position i. Reads that fail to place contribute nothing; portions of a
    /// read extending past the reference end are clipped and do not contribute.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> or <paramref name="reads"/> is null.</exception>
    public static int[] CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = DefaultMinOverlap)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(reads);

        var coverage = new int[reference.Length];

        foreach (string read in reads)
        {
            // Place the read at its best ungapped match, then increment the depth of every
            // reference position the read spans (clipped at the reference end).
            int pos = FindBestAlignment(reference, read, minOverlap);
            if (pos >= 0)
            {
                for (int i = pos; i < pos + read.Length && i < reference.Length; i++)
                {
                    coverage[i]++;
                }
            }
        }

        return coverage;
    }

    private static int FindBestAlignment(string reference, string read, int minOverlap)
    {
        int bestPos = -1;
        int bestScore = minOverlap - 1;

        for (int pos = 0; pos <= reference.Length - read.Length; pos++)
        {
            int matches = 0;
            for (int i = 0; i < read.Length; i++)
            {
                if (char.ToUpperInvariant(reference[pos + i]) ==
                    char.ToUpperInvariant(read[i]))
                    matches++;
            }

            if (matches > bestScore)
            {
                bestScore = matches;
                bestPos = pos;
            }
        }

        return bestPos;
    }

    /// <summary>
    /// Computes the consensus sequence from multiple aligned reads.
    /// </summary>
    public static string ComputeConsensus(IReadOnlyList<string> alignedReads)
    {
        if (alignedReads.Count == 0) return "";

        int length = alignedReads[0].Length;
        var sb = new StringBuilder();

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            foreach (string read in alignedReads)
            {
                if (pos < read.Length)
                {
                    char c = char.ToUpperInvariant(read[pos]);
                    if (c != '-' && c != 'N')
                    {
                        counts[c] = counts.GetValueOrDefault(c, 0) + 1;
                    }
                }
            }

            if (counts.Count > 0)
            {
                char consensus = counts.MaxBy(kvp => kvp.Value).Key;
                sb.Append(consensus);
            }
            else
            {
                sb.Append('N');
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Quality trims reads based on quality scores.
    /// </summary>
    public static IReadOnlyList<string> QualityTrimReads(
        IReadOnlyList<(string sequence, string quality)> reads,
        int minQuality = 20,
        int minLength = 50)
    {
        var trimmed = new List<string>();

        foreach (var (sequence, quality) in reads)
        {
            int start = 0;
            int end = sequence.Length;

            // Trim from start
            while (start < end && (quality[start] - 33) < minQuality)
                start++;

            // Trim from end
            while (end > start && (quality[end - 1] - 33) < minQuality)
                end--;

            if (end - start >= minLength)
            {
                trimmed.Add(sequence.Substring(start, end - start));
            }
        }

        return trimmed;
    }

    /// <summary>
    /// Error corrects reads using k-mer frequency analysis.
    /// </summary>
    public static IReadOnlyList<string> ErrorCorrectReads(
        IReadOnlyList<string> reads,
        int kmerSize = 21,
        int minKmerFrequency = 3)
    {
        // Build k-mer frequency table
        var kmerCounts = new Dictionary<string, int>();
        foreach (string read in reads)
        {
            for (int i = 0; i <= read.Length - kmerSize; i++)
            {
                string kmer = read.Substring(i, kmerSize).ToUpperInvariant();
                kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer, 0) + 1;
            }
        }

        var corrected = new List<string>();

        foreach (string read in reads)
        {
            var sb = new StringBuilder(read);

            for (int i = 0; i <= read.Length - kmerSize; i++)
            {
                string kmer = read.Substring(i, kmerSize).ToUpperInvariant();

                if (kmerCounts.GetValueOrDefault(kmer, 0) < minKmerFrequency)
                {
                    // Try to correct by substituting middle base
                    int midPos = i + kmerSize / 2;
                    char original = char.ToUpperInvariant(sb[midPos]);

                    foreach (char replacement in new[] { 'A', 'C', 'G', 'T' })
                    {
                        if (replacement == original) continue;

                        sb[midPos] = replacement;
                        string newKmer = sb.ToString().Substring(i, kmerSize).ToUpperInvariant();

                        if (kmerCounts.GetValueOrDefault(newKmer, 0) >= minKmerFrequency)
                        {
                            break; // Keep correction
                        }
                        else
                        {
                            sb[midPos] = original; // Revert
                        }
                    }
                }
            }

            corrected.Add(sb.ToString());
        }

        return corrected;
    }
}

