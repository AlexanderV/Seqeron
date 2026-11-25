using System.Text;

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
    /// <param name="Contigs">The reported contigs.</param>
    /// <param name="TotalReads">Number of input reads.</param>
    /// <param name="AssembledReads">Reads that contribute to at least one reported contig (Velvet
    /// <c>usedReads</c>, "using X/Y reads"); see <see cref="CalculateStats(IReadOnlyList{string}, int, int)"/>.</param>
    /// <param name="N50">Contig N50 (QUAST <c>N50()</c>; 0 for no contigs).</param>
    /// <param name="LongestContig">Length of the longest contig.</param>
    /// <param name="TotalLength">Sum of the contig lengths (64-bit).</param>
    public readonly record struct AssemblyResult(
        IReadOnlyList<string> Contigs,
        int TotalReads,
        int AssembledReads,
        double N50,
        int LongestContig,
        long TotalLength);

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
    /// The documented default assembly parameters (MinOverlap 20, MinIdentity 0.9, KmerSize 31,
    /// MinContigLength 100), built through the primary constructor so the declared defaults apply.
    /// </summary>
    public static readonly AssemblyParameters DefaultParameters = new(
        MinOverlap: 20, MinIdentity: 0.9, KmerSize: 31, MinContigLength: 100);

    /// <summary>
    /// Assembles reads using the Overlap-Layout-Consensus (OLC) paradigm:
    /// (1) Overlap — reads contained in another read (a full-length window of the container
    /// with identity ≥ <see cref="AssemblyParameters.MinIdentity"/>) are set aside, and the
    /// directed overlap graph of the remaining (substring-free) reads is built from the longest
    /// suffix-prefix overlap (length ≥ <see cref="AssemblyParameters.MinOverlap"/>) of every
    /// ordered pair; (2) Layout — the GREEDY algorithm: edges are scanned in non-increasing
    /// overlap order and an edge i → j is taken iff i has no successor yet, j has no predecessor
    /// yet and the edge does not close a cycle; every resulting path is one contig, and each
    /// contained read is placed at its offset inside its container; (3) Consensus — each layout
    /// column is resolved by majority vote (<see cref="ComputeConsensus"/>, threshold 0.5,
    /// ties → 'N'). Exact OLC layout (a Hamiltonian path through the overlap graph) is
    /// NP-complete, so GREEDY is a heuristic: it reconstructs unambiguous tilings exactly but
    /// may split or mis-resolve repeats longer than the reads.
    /// </summary>
    /// <remarks>
    /// References: Compeau, Pevzner &amp; Tesler (2011), Nat Biotechnol 29:987–991,
    /// DOI 10.1038/nbt.2023 (overlap graph, Hamiltonian-path layout, NP-completeness);
    /// Blum, Jiang, Li, Tromp &amp; Yannakakis (1994), J ACM 41:630–647 (GREEDY: merge the pair
    /// with maximum overlap; equivalently take overlap-graph edges in non-increasing order that
    /// share no head/tail with a taken edge and close no cycle); Myers (2005), Bioinformatics
    /// 21 Suppl 2:ii79 (contained reads are removed from the overlap/string graph);
    /// Langmead, "Overlap Layout Consensus assembly" (JHU lecture notes), p.4–5, 25, 28
    /// (consensus = per-column majority vote).
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <see cref="AssemblyParameters.MinOverlap"/> &lt; 1, or <see cref="AssemblyParameters.MinIdentity"/>
    /// outside [0, 1].
    /// </exception>
    public static AssemblyResult AssembleOLC(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
    {
        // NB: `new AssemblyParameters()` on a record struct is zero-initialisation (MinOverlap 0,
        // MinIdentity 0, MinContigLength 0), NOT the declared defaults — use DefaultParameters.
        var param = parameters ?? DefaultParameters;
        ValidateOverlapThresholds(param.MinOverlap, param.MinIdentity);

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        // Overlap + Layout + Consensus (greedy layout over the substring-free read set).
        // Every read is placed in exactly one layout (as a path member or a contained read).
        var layouts = BuildContigsFromOverlaps(reads, param);

        // Filter by minimum length; a read is assembled when its contig is reported
        // (Velvet usedReads: reads on nodes of at least the minimum contig length).
        var kept = layouts.Where(c => c.Contig.Length >= param.MinContigLength).ToList();
        var contigs = kept.Select(c => c.Contig).ToList();
        int assembledReads = kept.Sum(c => c.ReadCount);

        return CalculateStats(contigs, reads.Count, assembledReads);
    }

    /// <summary>
    /// Assembles reads using the de Bruijn graph (DBG) paradigm and returns its contigs — the
    /// strings spelled by the <b>maximal non-branching paths</b> (unitigs) of the graph.
    /// (1) The reads are decomposed into their <b>distinct</b> k-mers (case-insensitive,
    /// upper-cased; reads shorter than k contribute none), so read redundancy (coverage) does not
    /// add edges; (2) each distinct k-mer is a directed edge from its (k-1)-prefix to its
    /// (k-1)-suffix node; (3) every maximal non-branching path — a path whose internal nodes are
    /// 1-in-1-out and whose end nodes are not — and every isolated cycle is spelled out
    /// (<c>path[0]</c> + last character of each subsequent node) as one contig
    /// (<see cref="FindMaximalNonBranchingPaths"/>); (4) contigs shorter than
    /// <see cref="AssemblyParameters.MinContigLength"/> are discarded. Every distinct k-mer lies in
    /// exactly one contig; repeats longer than k-2 split the genome at the branch nodes, so a
    /// repeat-free genome is reconstructed as one contig while repeats yield several contigs, each
    /// of which must occur in any genome with this k-mer set. Only the forward strand is used;
    /// see the overload with <c>includeReverseComplement</c> for double-stranded reads.
    /// </summary>
    /// <remarks>
    /// References: Compeau &amp; Pevzner, <i>Bioinformatics Algorithms</i> ch.3 (Rosalind BA3K
    /// "contigs correspond to strings spelled by maximal non-branching paths in the de Bruijn
    /// graph"; BA3M MaximalNonBranchingPaths, incl. isolated cycles); Zerbino &amp; Birney (2008),
    /// Genome Res 18:821–829, DOI 10.1101/gr.074492.107 ("the fundamental data structure ... is
    /// based on k-mers, not reads, thus high redundancy is naturally handled by the graph without
    /// affecting the number of nodes"; Simplification: a node with one outgoing arc to a node with
    /// one ingoing arc is merged); Compeau, Pevzner &amp; Tesler (2011), Nat Biotechnol
    /// 29:987–991, DOI 10.1038/nbt.2023; Langmead B., "De Bruijn Graph assembly" (JHU notes),
    /// p.5–9 (nodes/edges), p.18 (spelling), p.24–28 (gaps, coverage and errors make the graph
    /// non-Eulerian, so a single Eulerian walk is not a practical contig definition).
    /// When <paramref name="parameters"/> is null, <see cref="DefaultParameters"/> (k = 31,
    /// MinContigLength 100) apply.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><see cref="AssemblyParameters.KmerSize"/> &lt; 2.</exception>
    public static AssemblyResult AssembleDeBruijn(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
        => AssembleDeBruijn(reads, parameters, includeReverseComplement: false);

    /// <summary>
    /// De Bruijn graph assembly as in <see cref="AssembleDeBruijn(IReadOnlyList{string}, AssemblyParameters?)"/>,
    /// optionally treating the reads as double-stranded DNA. With
    /// <paramref name="includeReverseComplement"/> = true every k-mer is entered together with its
    /// reverse complement (Velvet: "each k-mer is recorded simultaneously to its reverse complement
    /// ... this ensures that overlaps between reads from opposite strands are taken into account"),
    /// so reads from either strand join into one contig. The resulting graph is
    /// reverse-complement symmetric — each unitig's reverse complement is also a unitig — so each
    /// such pair is reported once, as the ordinally smaller of the two strings. Reverse
    /// complements use the canonical IUPAC complement (<see cref="SequenceExtensions.TryGetReverseComplement"/>).
    /// An odd k guarantees no k-mer is its own reverse complement (Velvet requires odd k);
    /// an even k is accepted.
    /// </summary>
    /// <remarks>Reference: Zerbino &amp; Birney (2008), Genome Res 18:821, "Construction" and "Structure" (twin nodes).</remarks>
    /// <exception cref="ArgumentOutOfRangeException"><see cref="AssemblyParameters.KmerSize"/> &lt; 2.</exception>
    public static AssemblyResult AssembleDeBruijn(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters,
        bool includeReverseComplement)
    {
        // NB: `new AssemblyParameters()` on a record struct is zero-initialisation (KmerSize 0,
        // MinContigLength 0), NOT the declared defaults — use DefaultParameters.
        var param = parameters ?? DefaultParameters;
        int k = param.KmerSize;
        ValidateDeBruijnK(k);

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        // (1) Distinct k-mers through the canonical counter (case-insensitive, upper-case keys).
        var kmers = new HashSet<string>(StringComparer.Ordinal);
        foreach (string read in reads)
        {
            if (read == null || read.Length < k)
                continue; // a read shorter than k contributes no k-mers (Langmead DBG p.16)

            foreach (string kmer in read.AsSpan().CountKmersSpan(k).Keys)
            {
                kmers.Add(kmer);
                if (includeReverseComplement)
                    kmers.Add(ReverseComplement(kmer));
            }
        }

        // (2) One edge per distinct k-mer; nodes and edges in ordinal order for determinism.
        var graph = new Dictionary<string, List<string>>(StringComparer.Ordinal);
        foreach (string kmer in kmers.Order(StringComparer.Ordinal))
        {
            string prefix = kmer[..^1];
            if (!graph.TryGetValue(prefix, out var successors))
                graph[prefix] = successors = new List<string>();
            successors.Add(kmer[1..]);
        }

        // (3) Contigs = spelled maximal non-branching paths and isolated cycles.
        var paths = FindMaximalNonBranchingPaths(graph);
        var contigs = includeReverseComplement
            ? SpellStrandPairsOnce(paths)
            : paths.Select(SpellPath).ToList();

        // (4) Filter by minimum contig length.
        contigs = contigs.Where(c => c.Length >= param.MinContigLength).ToList();

        // (5) A read is assembled when at least one of its k-mers (or, double-stranded, their
        // reverse complement) lies in a reported contig — Velvet usedReads (graphStats.c).
        return CalculateStats(contigs, reads.Count, CountReadsInContigs(reads, contigs, k, includeReverseComplement));
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
        ValidateOverlapThresholds(minOverlap, minIdentity);
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
        ValidateOverlapThresholds(minOverlap, minIdentity);
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
    /// Identity is <c>1 − Hamming/L</c> over the length-<c>L</c> window (case-insensitive).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minOverlap"/> &lt; 1 (an overlap of length 0 is no overlap), or
    /// <paramref name="minIdentity"/> outside [0, 1].
    /// </exception>
    public static (int length, int pos1, int pos2)? FindOverlap(
        string seq1, string seq2,
        int minOverlap = 20,
        double minIdentity = 0.9)
    {
        ValidateOverlapThresholds(minOverlap, minIdentity);

        // Check if suffix of seq1 overlaps with prefix of seq2
        int maxPossible = Math.Min(seq1.Length, seq2.Length);

        for (int overlapLen = maxPossible; overlapLen >= minOverlap; overlapLen--)
        {
            var suffix = seq1.AsSpan(seq1.Length - overlapLen);
            var prefix = seq2.AsSpan(0, overlapLen);

            if (WindowIdentity(suffix, prefix) >= minIdentity)
            {
                return (overlapLen, seq1.Length - overlapLen, 0);
            }
        }

        return null;
    }

    /// <summary>
    /// Calculates the fractional identity of two equal-length strings,
    /// <c>1 − HammingDistance / length</c> (case-insensitive). Returns 0 for unequal lengths
    /// and 1 for two empty strings.
    /// </summary>
    public static double CalculateIdentity(string seq1, string seq2)
    {
        if (seq1.Length != seq2.Length) return 0;
        return WindowIdentity(seq1, seq2);
    }

    /// <summary>
    /// Identity of two equal-length windows via the canonical Hamming distance
    /// (<see cref="SequenceExtensions.HammingDistance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>, PAT-APPROX-001).
    /// </summary>
    private static double WindowIdentity(ReadOnlySpan<char> a, ReadOnlySpan<char> b)
    {
        if (a.Length == 0) return 1;
        return (double)(a.Length - a.HammingDistance(b)) / a.Length; // matches / L
    }

    private static void ValidateOverlapThresholds(int minOverlap, double minIdentity)
    {
        if (minOverlap < 1)
            throw new ArgumentOutOfRangeException(nameof(minOverlap), minOverlap,
                "Minimum overlap length must be at least 1.");
        if (double.IsNaN(minIdentity) || minIdentity < 0.0 || minIdentity > 1.0)
            throw new ArgumentOutOfRangeException(nameof(minIdentity), minIdentity,
                "Minimum identity must be in [0, 1].");
    }

    /// <summary>
    /// Finds, for read <paramref name="j"/>, the read that contains it: a read that dominates it
    /// (strictly longer, or equal length and lower index — so identical copies keep the first)
    /// and has a length-|j| window with identity ≥ <paramref name="minIdentity"/>. Among several
    /// containers the longest (then lowest-index) wins; within it the leftmost best window.
    /// Reads shorter than <paramref name="minOverlap"/> are never treated as contained (an
    /// overlap shorter than the threshold is not an overlap).
    /// </summary>
    private static (int container, int offset)? FindContainer(
        IReadOnlyList<string> reads, int j, int minOverlap, double minIdentity)
    {
        string read = reads[j];
        if (read.Length < minOverlap) return null;

        (int container, int offset)? best = null;
        for (int i = 0; i < reads.Count; i++)
        {
            if (i == j) continue;
            string host = reads[i];
            bool dominates = host.Length > read.Length || (host.Length == read.Length && i < j);
            if (!dominates) continue;
            if (best.HasValue)
            {
                string cur = reads[best.Value.container];
                if (host.Length <= cur.Length) continue; // keep the longest (ties: lowest index)
            }

            int bestPos = -1;
            double bestId = -1;
            for (int p = 0; p + read.Length <= host.Length; p++)
            {
                double id = WindowIdentity(host.AsSpan(p, read.Length), read);
                if (id >= minIdentity && id > bestId)
                {
                    bestId = id;
                    bestPos = p;
                    if (id == 1.0) break;
                }
            }

            if (bestPos >= 0)
                best = (i, bestPos);
        }

        return best;
    }

    /// <summary>
    /// Overlap → Layout → Consensus for <see cref="AssembleOLC"/>: removes contained reads,
    /// lays out the substring-free reads with the GREEDY edge-selection rule (non-increasing
    /// overlap; out(i) = 0, in(j) = 0, no cycle), places contained reads inside their containers,
    /// and computes each contig as the majority-vote consensus of its ungapped layout columns.
    /// </summary>
    private static List<(string Contig, int ReadCount)> BuildContigsFromOverlaps(
        IReadOnlyList<string> reads,
        AssemblyParameters param)
    {
        int n = reads.Count;

        // --- Containment (Myers 2005): contained reads are removed from the overlap graph.
        var containedIn = new (int container, int offset)?[n];
        for (int j = 0; j < n; j++)
            containedIn[j] = FindContainer(reads, j, param.MinOverlap, param.MinIdentity);

        var isNode = new bool[n];
        var nodeIndex = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (!containedIn[i].HasValue)
            {
                isNode[i] = true;
                nodeIndex.Add(i);
            }
        }

        // --- Overlap graph over the substring-free reads (FindAllOverlaps on the sub-list).
        var nodeReads = nodeIndex.Select(i => reads[i]).ToList();
        var edges = FindAllOverlaps(nodeReads, param.MinOverlap, param.MinIdentity)
            .Select(o => new Overlap(nodeIndex[o.ReadIndex1], nodeIndex[o.ReadIndex2],
                o.OverlapLength, o.Position1, o.Position2))
            .OrderByDescending(o => o.OverlapLength)
            .ThenBy(o => o.ReadIndex1)
            .ThenBy(o => o.ReadIndex2)
            .ToList();

        // --- Layout: GREEDY (Blum et al. 1994) with union-find cycle check.
        var successor = new int[n];
        var overlapToSuccessor = new int[n];
        var hasPredecessor = new bool[n];
        var component = new int[n];
        for (int i = 0; i < n; i++)
        {
            successor[i] = -1;
            component[i] = i;
        }

        int Find(int x)
        {
            while (component[x] != x)
            {
                component[x] = component[component[x]];
                x = component[x];
            }
            return x;
        }

        foreach (var e in edges)
        {
            int a = e.ReadIndex1, b = e.ReadIndex2;
            if (successor[a] != -1 || hasPredecessor[b]) continue;
            int ra = Find(a), rb = Find(b);
            if (ra == rb) continue; // would close a cycle
            successor[a] = b;
            overlapToSuccessor[a] = e.OverlapLength;
            hasPredecessor[b] = true;
            component[rb] = ra;
        }

        // Resolve each contained read to its root (non-contained) container and cumulative offset.
        (int root, int offset) ResolveRoot(int j)
        {
            int offset = 0;
            while (containedIn[j].HasValue)
            {
                offset += containedIn[j]!.Value.offset;
                j = containedIn[j]!.Value.container;
            }
            return (j, offset);
        }

        var placedByRoot = new Dictionary<int, List<(int read, int offset)>>();
        for (int j = 0; j < n; j++)
        {
            if (isNode[j]) continue;
            var (root, off) = ResolveRoot(j);
            if (!placedByRoot.TryGetValue(root, out var list))
                placedByRoot[root] = list = new List<(int, int)>();
            list.Add((j, off));
        }

        // --- Consensus: one contig per greedy path (every path starts at a node with no predecessor).
        var contigs = new List<(string Contig, int ReadCount)>();
        foreach (int start in nodeIndex)
        {
            if (hasPredecessor[start]) continue;

            var rows = new List<(int offset, string read)>();
            int current = start, offsetInContig = 0;
            while (current != -1)
            {
                rows.Add((offsetInContig, reads[current]));
                if (placedByRoot.TryGetValue(current, out var contained))
                {
                    foreach (var (read, off) in contained)
                        rows.Add((offsetInContig + off, reads[read]));
                }

                int next = successor[current];
                if (next != -1)
                    offsetInContig += reads[current].Length - overlapToSuccessor[current];
                current = next;
            }

            if (rows.Count == 1)
            {
                // A one-read layout: its consensus is the read itself (returned verbatim).
                contigs.Add((rows[0].read, 1));
                continue;
            }

            // Ungapped layout: pad each read with gap symbols to its column offset and take the
            // column-wise majority vote (Biopython dumb_consensus rule; ties -> 'N', upper-cased).
            var aligned = rows.Select(r => new string(GapDash, r.offset) + r.read).ToList();
            contigs.Add((ComputeConsensus(aligned), rows.Count));
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
        ValidateDeBruijnK(k);

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
    /// Finds all maximal non-branching paths of a directed (multi)graph given as an
    /// out-adjacency list (Compeau &amp; Pevzner, <i>Bioinformatics Algorithms</i> ch.3; Rosalind
    /// BA3M MaximalNonBranchingPaths). A node is 1-in-1-out when its in-degree and out-degree
    /// are both 1. For every node that is not 1-in-1-out and has out-degree &gt; 0, each of its
    /// outgoing edges (v, w) starts a path that is extended through w while w is 1-in-1-out;
    /// afterwards every isolated cycle (a component made only of 1-in-1-out nodes) is added as a
    /// path that starts and ends at its ordinally smallest node. Parallel edges are separate
    /// paths. Every edge lies on exactly one returned path.
    /// </summary>
    /// <param name="graph">Out-adjacency list (e.g. from <see cref="BuildDeBruijnGraph"/>); nodes
    /// that appear only as successors are included.</param>
    /// <returns>The paths as node sequences (each of at least two nodes), the non-1-in-1-out start
    /// nodes taken in ordinal order and their edges in stored order, followed by the isolated
    /// cycles.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="graph"/> is null.</exception>
    public static IReadOnlyList<IReadOnlyList<string>> FindMaximalNonBranchingPaths(
        IReadOnlyDictionary<string, List<string>> graph)
    {
        ArgumentNullException.ThrowIfNull(graph);

        var inDegree = new Dictionary<string, int>(StringComparer.Ordinal);
        var nodes = new SortedSet<string>(StringComparer.Ordinal);
        foreach (var (node, successors) in graph)
        {
            nodes.Add(node);
            foreach (string w in successors)
            {
                nodes.Add(w);
                inDegree[w] = inDegree.GetValueOrDefault(w) + 1;
            }
        }

        int OutDegree(string v) => graph.TryGetValue(v, out var s) ? s.Count : 0;
        bool IsOneInOneOut(string v) => inDegree.GetValueOrDefault(v) == 1 && OutDegree(v) == 1;

        var paths = new List<IReadOnlyList<string>>();
        // A 1-in-1-out node has exactly one out-edge, so "its edge was used" identifies it.
        var usedOneInOneOut = new HashSet<string>(StringComparer.Ordinal);

        foreach (string v in nodes)
        {
            if (IsOneInOneOut(v) || OutDegree(v) == 0)
                continue;

            foreach (string first in graph[v])
            {
                var path = new List<string> { v, first };
                string w = first;
                while (IsOneInOneOut(w))
                {
                    usedOneInOneOut.Add(w);
                    w = graph[w][0];
                    path.Add(w);
                }
                paths.Add(path);
            }
        }

        // Isolated cycles: 1-in-1-out nodes not reached from any branching / end node.
        foreach (string v in nodes)
        {
            if (!IsOneInOneOut(v) || usedOneInOneOut.Contains(v))
                continue;

            var cycle = new List<string> { v };
            string w = v;
            do
            {
                usedOneInOneOut.Add(w);
                w = graph[w][0];
                cycle.Add(w);
            }
            while (w != v);
            paths.Add(cycle);
        }

        return paths;
    }

    /// <summary>
    /// Spells a node path: the first (k-1)-mer in full, then the last character of every
    /// subsequent node (Langmead DBG p.18; Compeau &amp; Pevzner PathToGenome).
    /// </summary>
    private static string SpellPath(IReadOnlyList<string> path)
    {
        var sb = new StringBuilder(path[0], path[0].Length + path.Count - 1);
        for (int i = 1; i < path.Count; i++)
            sb.Append(path[i][^1]);
        return sb.ToString();
    }

    /// <summary>
    /// For a reverse-complement-symmetric graph, keeps one contig of each {unitig, reverse
    /// complement} pair: a path is skipped when its first k-mer already belongs to an emitted
    /// path or to that path's reverse complement; an emitted contig is the ordinally smaller of
    /// the spelled path and its reverse complement.
    /// </summary>
    private static List<string> SpellStrandPairsOnce(IReadOnlyList<IReadOnlyList<string>> paths)
    {
        var covered = new HashSet<string>(StringComparer.Ordinal);
        var contigs = new List<string>();
        foreach (var path in paths)
        {
            string contig = SpellPath(path);
            int k = path[0].Length + 1;
            if (covered.Contains(contig[..k]))
                continue;

            for (int i = 0; i + k <= contig.Length; i++)
            {
                string kmer = contig.Substring(i, k);
                covered.Add(kmer);
                covered.Add(ReverseComplement(kmer));
            }

            string reverse = ReverseComplement(contig);
            contigs.Add(StringComparer.Ordinal.Compare(reverse, contig) < 0 ? reverse : contig);
        }

        return contigs;
    }

    /// <summary>
    /// Number of reads with at least one k-mer in the reported contigs (either strand when
    /// <paramref name="bothStrands"/>); Velvet <c>usedReads</c> (Zerbino &amp; Birney 2008,
    /// <c>src/graphStats.c</c>: a read is used when it has a marker on a node of at least the
    /// minimum contig length, the twin node included).
    /// </summary>
    private static int CountReadsInContigs(
        IReadOnlyList<string> reads, IReadOnlyList<string> contigs, int k, bool bothStrands)
    {
        var contigKmers = new HashSet<string>(StringComparer.Ordinal);
        foreach (string contig in contigs)
        {
            foreach (string kmer in contig.AsSpan().CountKmersSpan(k).Keys)
            {
                contigKmers.Add(kmer);
                if (bothStrands)
                    contigKmers.Add(ReverseComplement(kmer));
            }
        }

        int used = 0;
        foreach (string read in reads)
        {
            if (read == null || read.Length < k)
                continue;
            if (read.AsSpan().CountKmersSpan(k).Keys.Any(contigKmers.Contains))
                used++;
        }

        return used;
    }

    /// <summary>Reverse complement through the canonical Core helper.</summary>
    private static string ReverseComplement(string sequence) =>
        string.Create(sequence.Length, sequence,
            static (dest, src) => src.AsSpan().TryGetReverseComplement(dest));

    /// <summary>k-mers split into two (k-1)-mer nodes, so k must be at least 2 (Langmead DBG p.5).</summary>
    private static void ValidateDeBruijnK(int k)
    {
        const int MinKmerLength = 2;
        if (k < MinKmerLength)
            throw new ArgumentOutOfRangeException(nameof(k), k,
                $"k-mer length must be at least {MinKmerLength} so that (k-1)-mer nodes are non-empty.");
    }

    /// <summary>
    /// Calculates the assembly statistics of a precomputed contig set: contig N50 (through the
    /// canonical <see cref="CalculateNx(IReadOnlyList{int}, long, int)"/>), longest contig and total
    /// length. <paramref name="totalReads"/> is the number of reads the contigs were built from and
    /// is reported as <see cref="AssemblyResult.TotalReads"/> and, when there is at least one contig,
    /// as <see cref="AssemblyResult.AssembledReads"/>: this overload has no read placement, so every
    /// supplied read is taken as used (no contig → 0 used reads). Use
    /// <see cref="CalculateStats(IReadOnlyList{string}, int, int)"/> when the number of reads that
    /// lie in the contigs is known (the assemblers do this).
    /// </summary>
    /// <remarks>
    /// N50 is QUAST's <c>N50()</c> (<c>quast_libs/N50.py</c>): the length of the shortest contig in the
    /// fewest largest contigs whose combined length is at least 50% of the total (Miller, Koren &amp;
    /// Sutton 2010, Genomics 95:315, §1.2). Lengths are summed in 64-bit, so totals above
    /// <see cref="int.MaxValue"/> are exact.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="contigs"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="totalReads"/> is negative.</exception>
    public static AssemblyResult CalculateStats(
        IReadOnlyList<string> contigs, int totalReads)
        => CalculateStats(contigs, totalReads, contigs is { Count: > 0 } ? totalReads : 0);

    /// <summary>
    /// Calculates the assembly statistics of a contig set (as
    /// <see cref="CalculateStats(IReadOnlyList{string}, int)"/>) with an explicit count of
    /// <paramref name="assembledReads"/> — the reads that contribute to at least one reported
    /// contig (Velvet <c>usedReads</c>: "using X/Y reads", reads with a marker on a node of at
    /// least the minimum contig length; Zerbino &amp; Birney 2008, <c>src/graphStats.c</c>).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="contigs"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="totalReads"/> is negative, or <paramref name="assembledReads"/> is outside [0, totalReads].
    /// </exception>
    public static AssemblyResult CalculateStats(
        IReadOnlyList<string> contigs, int totalReads, int assembledReads)
    {
        ArgumentNullException.ThrowIfNull(contigs);
        ArgumentOutOfRangeException.ThrowIfNegative(totalReads);
        if (assembledReads < 0 || assembledReads > totalReads)
            throw new ArgumentOutOfRangeException(nameof(assembledReads), assembledReads,
                "The number of assembled reads must lie in [0, totalReads].");

        if (contigs.Count == 0)
            return new AssemblyResult(contigs, totalReads, assembledReads, 0, 0, 0);

        var sortedLengths = contigs.Select(c => c?.Length ?? 0).OrderDescending().ToArray();
        long totalLength = 0;
        foreach (int length in sortedLengths)
            totalLength += length;

        var (n50, _, _) = CalculateNx(sortedLengths, totalLength, N50Percent);

        return new AssemblyResult(
            contigs,
            totalReads,
            assembledReads,
            n50,
            sortedLengths[0],
            totalLength);
    }

    // N50 = Nx at x = 50 (Miller, Koren & Sutton 2010, §1.2).
    private const int N50Percent = 50;

    /// <summary>
    /// Canonical Nx/Lx (and NGx/LGx) core shared by every assembly-statistics method
    /// (<see cref="CalculateStats(IReadOnlyList{string}, int, int)"/> and
    /// <c>GenomeAssemblyAnalyzer.CalculateNx</c>). Walking the lengths largest-first, Nx is the
    /// length of the first sequence at which the cumulative length reaches at least
    /// <paramref name="percentage"/>% of <paramref name="referenceLength"/>, and Lx is the number
    /// of sequences taken. With <paramref name="referenceLength"/> = the sum of the lengths this is
    /// Nx/Lx; with an (estimated) genome size it is NGx/LGx.
    /// </summary>
    /// <remarks>
    /// Port of QUAST <c>NG50_and_LG50</c> (<c>quast_libs/N50.py</c>): <c>limit = ref·(100 − x)/100</c>,
    /// stop at the first <c>ref − cumulative ≤ limit</c>, evaluated here exactly as
    /// <c>cumulative·100 ≥ ref·x</c> in 64-bit integers (no floating-point cutoff, no integer
    /// halving). Definition: Miller, Koren &amp; Sutton (2010), Genomics 95(6):315-327, §1.2
    /// ("... whose combined length represents at least 50% of the assembly"). Where QUAST returns
    /// <c>None</c> — an empty list, or a reference length the lengths never reach (NGx) — this
    /// returns <c>(0, 0, cumulative)</c>. A non-empty all-zero list gives <c>(0, 1, 0)</c> as in QUAST.
    /// Examples (QUAST): {4,3,2} → N50 3, L50 2; {2,1,1,1} → N50 1, L50 2; {80,70,50,40,30,20} →
    /// N50 70, L50 2, N90 30, L90 5; same lengths with reference 400 → NG50 50, LG50 3.
    /// </remarks>
    /// <param name="lengthsDescending">Non-negative sequence lengths in non-increasing order.</param>
    /// <param name="referenceLength">Assembly size (sum of the lengths) for Nx, or genome size for NGx; ≥ 0.</param>
    /// <param name="percentage">x in [0, 100].</param>
    /// <returns>(Nx, Lx, cumulative length of the first Lx sequences).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="lengthsDescending"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="percentage"/> outside [0, 100] (QUAST asserts the same range) or a negative <paramref name="referenceLength"/>.
    /// </exception>
    /// <exception cref="ArgumentException">A length is negative or the lengths are not in non-increasing order.</exception>
    public static (int Nx, int Lx, long CumulativeLength) CalculateNx(
        IReadOnlyList<int> lengthsDescending, long referenceLength, int percentage)
    {
        ArgumentNullException.ThrowIfNull(lengthsDescending);
        if (percentage is < 0 or > 100)
            throw new ArgumentOutOfRangeException(nameof(percentage), percentage,
                "The Nx percentage must lie in [0, 100].");
        ArgumentOutOfRangeException.ThrowIfNegative(referenceLength);

        long cumulative = 0;
        int count = 0;
        int previous = int.MaxValue;
        (int, int, long)? result = null;
        foreach (int length in lengthsDescending)
        {
            if (length < 0)
                throw new ArgumentException("Sequence lengths must be non-negative.", nameof(lengthsDescending));
            if (length > previous)
                throw new ArgumentException("Sequence lengths must be sorted in non-increasing order.", nameof(lengthsDescending));
            previous = length;

            if (result.HasValue)
                continue; // keep validating the remaining lengths

            cumulative += length;
            count++;
            // QUAST: ref - cumulative <= ref * (100 - x) / 100  <=>  cumulative * 100 >= ref * x.
            if (cumulative * 100 >= referenceLength * percentage)
                result = (length, count, cumulative);
        }

        return result ?? (0, 0, cumulative);
    }

    /// <summary>
    /// Merges two contigs whose suffix/prefix overlap length is already known, collapsing the
    /// shared region so it appears once: the merged superstring is
    /// <paramref name="contig1"/> followed by <paramref name="contig2"/> with its length-<c>l</c>
    /// prefix removed (<c>l = overlapLength</c>), i.e. Langmead's <c>sup += y[olen:]</c>.
    /// The result length is <c>|contig1| + |contig2| − overlapLength</c>.
    /// </summary>
    /// <remarks>
    /// An overlap is a length-<c>l</c> suffix of <paramref name="contig1"/> that matches a
    /// length-<c>l</c> prefix of <paramref name="contig2"/>; merging keeps a single copy of that
    /// region (Langmead, "Assembly &amp; shortest common superstring", JHU notes — overlap
    /// definition and the greedy merge trace, e.g. <c>BAA</c> + <c>AAB</c> at overlap 2 → <c>BAAB</c>).
    /// <c>overlapLength = 0</c> is plain concatenation ("just concatenate them"). A valid overlap is
    /// bounded by <c>min(|x|,|y|)</c>, so a negative or larger length is rejected.
    /// The supplied <paramref name="overlapLength"/> is trusted, not re-verified (computing it is
    /// the job of <see cref="FindOverlap"/> / <c>FindAllOverlaps</c>); in the overlap the bases of
    /// <paramref name="contig1"/> are kept. To resolve mismatches/<c>N</c>s in a non-exact
    /// overlap use <see cref="MergeContigs(string, string, int, bool)"/>.
    /// </remarks>
    /// <param name="contig1">The left contig; its suffix overlaps the prefix of <paramref name="contig2"/>.</param>
    /// <param name="contig2">The right contig; its length-<paramref name="overlapLength"/> prefix is the shared region.</param>
    /// <param name="overlapLength">Length of the suffix(contig1)/prefix(contig2) overlap to collapse, in [0, min(|contig1|, |contig2|)].</param>
    /// <returns>The merged superstring.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contig1"/> or <paramref name="contig2"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="overlapLength"/> is negative or exceeds the shorter contig.</exception>
    public static string MergeContigs(string contig1, string contig2, int overlapLength)
        => MergeContigs(contig1, contig2, overlapLength, resolveOverlapConflicts: false);

    /// <summary>
    /// Merges two contigs over a known ungapped suffix/prefix overlap, optionally resolving
    /// disagreements inside the overlap with the conflict rule of EMBOSS <c>merger</c>
    /// (emboss/merger.c, <c>merger_Merge</c>/<c>merger_bestquality</c>/<c>merger_quality</c>).
    /// </summary>
    /// <remarks>
    /// With <paramref name="resolveOverlapConflicts"/> = <c>false</c> this is the textbook collapse
    /// <c>contig1 + contig2[overlapLength..]</c> (contig1's bases kept in the overlap). With
    /// <c>true</c>, each overlap column (case-insensitive) is resolved as EMBOSS <c>merger</c> does
    /// on an alignment column: agreeing bases are kept; an <c>N</c> in contig1 takes contig2's base,
    /// an <c>N</c> in contig2 takes contig1's base; otherwise the base of the contig with the higher
    /// local "quality" wins — quality = points/(2w+1) over a window of w positions forward
    /// (including the column) and w−1 backward in the whole contig, with ACGTU = 2 points,
    /// the IUPAC ambiguity codes MRWSYKVHDB = 1, anything else (N, off-end) = 0; windows
    /// w = 5, then 10, then 20 break ties, and a remaining tie keeps contig1. This favours the
    /// interior of a read over its end, where sequencing errors accumulate. Outside the overlap the
    /// two layouts coincide. The overlap is taken as given (ungapped); EMBOSS first finds it with
    /// a global alignment, which here is <see cref="FindOverlap"/>'s job.
    /// </remarks>
    /// <param name="contig1">The left contig; its suffix overlaps the prefix of <paramref name="contig2"/>.</param>
    /// <param name="contig2">The right contig; its length-<paramref name="overlapLength"/> prefix is the shared region.</param>
    /// <param name="overlapLength">Overlap length, in [0, min(|contig1|, |contig2|)].</param>
    /// <param name="resolveOverlapConflicts">Resolve overlap mismatches / <c>N</c>s by the EMBOSS <c>merger</c> rule.</param>
    /// <returns>The merged sequence, of length <c>|contig1| + |contig2| − overlapLength</c>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contig1"/> or <paramref name="contig2"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="overlapLength"/> is negative or exceeds the shorter contig.</exception>
    public static string MergeContigs(string contig1, string contig2, int overlapLength, bool resolveOverlapConflicts)
    {
        ArgumentNullException.ThrowIfNull(contig1);
        ArgumentNullException.ThrowIfNull(contig2);

        // A length-l overlap is both a suffix of contig1 and a prefix of contig2 (Langmead SCS notes),
        // so 0 ≤ l ≤ min(|contig1|, |contig2|); 0 = concatenation.
        int maxOverlap = Math.Min(contig1.Length, contig2.Length);
        if (overlapLength < 0 || overlapLength > maxOverlap)
            throw new ArgumentOutOfRangeException(nameof(overlapLength), overlapLength,
                $"Overlap length must be in [0, {maxOverlap}] (the shorter contig's length).");

        if (!resolveOverlapConflicts)
            return string.Concat(contig1, contig2.AsSpan(overlapLength));

        int start1 = contig1.Length - overlapLength;
        var merged = new StringBuilder(contig1.Length + contig2.Length - overlapLength);
        merged.Append(contig1, 0, start1);
        for (int i = 0; i < overlapLength; i++)
        {
            char p = contig1[start1 + i];
            char q = contig2[i];
            if (char.ToUpperInvariant(p) == 'N') merged.Append(q);
            else if (char.ToUpperInvariant(q) == 'N') merged.Append(p);
            else if (char.ToUpperInvariant(p) != char.ToUpperInvariant(q))
                merged.Append(MergerFirstHasBestQuality(contig1, contig2, start1 + i, i) ? p : q);
            else merged.Append(p);
        }

        merged.Append(contig2, overlapLength, contig2.Length - overlapLength);
        return merged.ToString();
    }

    /// <summary>EMBOSS merger_bestquality: windows 5, 10, 20; ties keep the first sequence.</summary>
    private static bool MergerFirstHasBestQuality(string a, string b, int apos, int bpos)
    {
        float qa = 0, qb = 0;
        foreach (int window in MergerQualityWindows)
        {
            qa = MergerWindowQuality(a, apos, window);
            qb = MergerWindowQuality(b, bpos, window);
            if (qa != qb) break;
        }

        return qa >= qb;
    }

    private static readonly int[] MergerQualityWindows = { 5, 10, 20 };

    /// <summary>
    /// EMBOSS merger_quality: points in [pos, pos+window) and (pos−window, pos) divided by
    /// 2·window+1; ACGTU = 2, IUPAC ambiguity codes MRWSYKVHDB = 1, other symbols and off-end = 0.
    /// </summary>
    private static float MergerWindowQuality(string seq, int pos, int window)
    {
        int value = 0;
        for (int i = pos; i < pos + window && i < seq.Length; i++)
            value += MergerBasePoints(seq[i]);
        for (int j = pos - 1; j > pos - window && j >= 0; j--)
            value += MergerBasePoints(seq[j]);
        return (float)(value / (window * 2.0 + 1.0));
    }

    private static int MergerBasePoints(char c) => char.ToUpperInvariant(c) switch
    {
        'A' or 'C' or 'G' or 'T' or 'U' => 2,
        'M' or 'R' or 'W' or 'S' or 'Y' or 'K' or 'V' or 'H' or 'D' or 'B' => 1,
        _ => 0,
    };

    /// <summary>
    /// Gap length, in characters, emitted when a link's distance estimate does not give a
    /// representable positive gap: a zero estimate, or a negative estimate (= expected overlap)
    /// whose overlap is not found in the sequences. 100 is the GenBank/EMBL/DDBJ standard length
    /// for a gap of unknown size (NCBI AGP Specification v2.1, §"Gaps": "For negative gaps, or gaps
    /// of unknown size, use U as the component_type and 100 as the gap size"; "Negative gaps and
    /// gap lines with zero length are not valid").
    /// </summary>
    private const int UnknownGapLength = 100;

    /// <summary>
    /// A scaffolding link (a mate-edge of the contig-mate-pair graph, Huson, Reinert &amp; Myers 2002):
    /// in the scaffold, <see cref="Contig1"/> read in its own orientation (reverse complement when
    /// <see cref="Contig1Reversed"/>) is followed, after <see cref="GapSize"/> characters, by
    /// <see cref="Contig2"/> read in its orientation (reverse complement when
    /// <see cref="Contig2Reversed"/>) — the ABySS distance-estimate notation <c>A+ B-,d</c>, and the
    /// AGP <c>+</c>/<c>-</c> component orientation. The same link read from the other strand is
    /// <c>Contig2(¬o2) → Contig1(¬o1)</c>, so a link joins one end of each contig.
    /// </summary>
    /// <param name="Contig1">Index of the left contig.</param>
    /// <param name="Contig2">Index of the right contig.</param>
    /// <param name="GapSize">Estimated distance between the two contigs; negative = expected overlap of −GapSize bases.</param>
    /// <param name="Contig1Reversed">The left contig is placed reverse-complemented (<c>-</c>).</param>
    /// <param name="Contig2Reversed">The right contig is placed reverse-complemented (<c>-</c>).</param>
    /// <param name="Weight">Link support (number of mate pairs bundled into this edge, Huson et al. 2002 §2.2); at least 1.</param>
    public readonly record struct ScaffoldLink(
        int Contig1,
        int Contig2,
        int GapSize,
        bool Contig1Reversed = false,
        bool Contig2Reversed = false,
        int Weight = 1);

    /// <summary>
    /// Joins contigs into scaffolds from forward-strand paired-end links
    /// <c>(contig1, contig2, gapSize)</c> = "contig2 follows contig1 after gapSize characters".
    /// Each tuple is a <see cref="ScaffoldLink"/> with both contigs forward and weight 1; see
    /// <see cref="Scaffold(IReadOnlyList{string}, IReadOnlyList{ScaffoldLink}, char)"/> for the
    /// greedy rule, gap/overlap handling and output order.
    /// </summary>
    /// <param name="contigs">The contigs to scaffold, indexed by position; links reference these indices.</param>
    /// <param name="links">Links <c>(contig1, contig2, gapSize)</c>; out-of-range and self links are ignored.</param>
    /// <param name="gapCharacter">The character used to fill gaps between contigs (default <c>'N'</c>).</param>
    /// <returns>The scaffolds, one string per scaffold, ordered by their lowest contig index.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contigs"/> or <paramref name="links"/> is null.</exception>
    public static IReadOnlyList<string> Scaffold(
        IReadOnlyList<string> contigs,
        IReadOnlyList<(int contig1, int contig2, int gapSize)> links,
        char gapCharacter = 'N')
    {
        ArgumentNullException.ThrowIfNull(contigs);
        ArgumentNullException.ThrowIfNull(links);

        var oriented = new ScaffoldLink[links.Count];
        for (int i = 0; i < links.Count; i++)
            oriented[i] = new ScaffoldLink(links[i].contig1, links[i].contig2, links[i].gapSize);
        return Scaffold(contigs, oriented, gapCharacter);
    }

    /// <summary>
    /// Orders and orients contigs into scaffolds with the greedy path-merging rule of Huson, Reinert
    /// &amp; Myers (2002, J. ACM 49:603, Algorithm 5.1; the Bambus "most links first" greedy, Pop et al.
    /// 2004), then spells each scaffold as the oriented contigs "interspersed with gaps represented by
    /// a run of the character N, whose length corresponds to the estimate of the distance"
    /// (Jackman et al. 2017, ABySS 2.0).
    /// </summary>
    /// <remarks>
    /// <para><b>Graph and greedy selection.</b> Each contig has two end nodes (5′, 3′); a link joins the
    /// exit end of <c>Contig1</c> (3′ if forward, 5′ if reversed) to the entry end of <c>Contig2</c>
    /// (5′ if forward, 3′ if reversed). Links are taken in descending <see cref="ScaffoldLink.Weight"/>
    /// (ties: input order). A link is selected iff neither end node already carries a selected link
    /// (every node is adjacent to at most two selected edges: its contig edge and one link) and the
    /// two contigs lie on different selected paths (a link inside one path would close a cycle).
    /// Rejected links are the conflicts: a second successor/predecessor at an end (branching), or a
    /// cycle-closing link. Huson et al. would additionally try to interleave ("zipper") two paths
    /// when the link hits a path interior, which needs per-link standard deviations to judge
    /// happiness; links here carry none, so such a link is not used.</para>
    /// <para><b>Orientation.</b> Each selected path is spelled from one terminal contig, every contig in
    /// the orientation the links imply (reverse complement via the canonical Core helper). Of the two
    /// readings (a path and its reverse complement) the one with more contigs forward is reported
    /// (SSPACE reverses a scaffold that is wholly reverse-complemented); on a tie, the one with the
    /// lowest-index contig forward. Forward-only links therefore spell forward scaffolds.</para>
    /// <para><b>Gaps.</b> A positive estimate <c>d</c> emits <c>d</c> gap characters (Jackman et al.
    /// 2017). A negative estimate means the contigs overlap by <c>−d</c> bases (Huson et al. 2002
    /// §2.1; Jackman et al. 2017): as in ABySS <c>MergeContigs</c>, if the last <c>−d</c> bases placed
    /// (since the previous gap) and the first <c>−d</c> bases of the next contig agree — equal
    /// ignoring case, or one of them <c>N</c> — the overlap is collapsed with
    /// <see cref="MergeContigs(string, string, int, bool)"/> (an <c>N</c> takes the other base).
    /// Otherwise (no such overlap, overlap longer than a contig, or a zero estimate) the gap is of
    /// unknown size and <see cref="UnknownGapLength"/> (100) gap characters are emitted (NCBI AGP
    /// v2.1). ABySS falls back to a gapped overlap alignment and IUPAC-subset consensus before
    /// giving up; those are not attempted here.</para>
    /// <para><b>Output.</b> Every contig is placed in exactly one scaffold; contigs without a selected
    /// link are singleton scaffolds. Scaffolds are listed by their lowest contig index.</para>
    /// </remarks>
    /// <param name="contigs">The contigs to scaffold, indexed by position; links reference these indices.</param>
    /// <param name="links">Oriented, weighted links; out-of-range and self links are ignored.</param>
    /// <param name="gapCharacter">The character used to fill gaps between contigs (default <c>'N'</c>).</param>
    /// <returns>The scaffolds, one string per scaffold, ordered by their lowest contig index.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="contigs"/> or <paramref name="links"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A link has <see cref="ScaffoldLink.Weight"/> &lt; 1.</exception>
    public static IReadOnlyList<string> Scaffold(
        IReadOnlyList<string> contigs,
        IReadOnlyList<ScaffoldLink> links,
        char gapCharacter = 'N')
    {
        ArgumentNullException.ThrowIfNull(contigs);
        ArgumentNullException.ThrowIfNull(links);
        for (int i = 0; i < links.Count; i++)
        {
            if (links[i].Weight < 1)
                throw new ArgumentOutOfRangeException(nameof(links), links[i].Weight,
                    $"Link {i} has weight {links[i].Weight}; a link is supported by at least one mate pair.");
        }

        int n = contigs.Count;
        if (n == 0) return Array.Empty<string>();

        // End node e = 2·contig + (0 = 5′/left, 1 = 3′/right). partner[e] = joined end node or −1.
        var partner = new int[2 * n];
        var joinGap = new int[2 * n];
        Array.Fill(partner, -1);
        var parent = new int[n];
        for (int i = 0; i < n; i++) parent[i] = i;

        int Find(int x)
        {
            while (parent[x] != x) x = parent[x] = parent[parent[x]];
            return x;
        }

        // Huson et al. Algorithm 5.1: mate-edges in descending weight; stable on input order.
        var order = Enumerable.Range(0, links.Count)
            .OrderByDescending(i => links[i].Weight)
            .ThenBy(i => i);
        foreach (int li in order)
        {
            var link = links[li];
            if ((uint)link.Contig1 >= (uint)n || (uint)link.Contig2 >= (uint)n || link.Contig1 == link.Contig2)
                continue;

            int exit = 2 * link.Contig1 + (link.Contig1Reversed ? 0 : 1);
            int entry = 2 * link.Contig2 + (link.Contig2Reversed ? 1 : 0);
            if (partner[exit] >= 0 || partner[entry] >= 0) continue;   // end already joined (branch)
            int r1 = Find(link.Contig1), r2 = Find(link.Contig2);
            if (r1 == r2) continue;                                      // same path (cycle)

            partner[exit] = entry;
            partner[entry] = exit;
            joinGap[exit] = joinGap[entry] = link.GapSize;
            parent[r1] = r2;
        }

        var placed = new bool[n];
        var scaffolds = new List<string>();
        for (int lowest = 0; lowest < n; lowest++)
        {
            if (placed[lowest]) continue;

            // Walk to one terminal of lowest's path, then read the path from there.
            var path = ReadScaffoldPath(TerminalOf(lowest, partner), partner);
            int forward = path.Count(p => !p.Reversed);
            bool flip = 2 * forward < path.Count
                || (2 * forward == path.Count && path.First(p => p.Contig == lowest).Reversed);
            if (flip) path = ReverseScaffoldPath(path);

            foreach (var p in path) placed[p.Contig] = true;
            scaffolds.Add(SpellScaffold(path, contigs, joinGap, gapCharacter));
        }

        return scaffolds;
    }

    private readonly record struct PlacedContig(int Contig, bool Reversed, int EntryEnd);

    /// <summary>Follows the path from <paramref name="contig"/> to a contig with a free end.</summary>
    private static int TerminalOf(int contig, int[] partner)
    {
        int current = contig;
        int outEnd = 2 * contig;                       // leave through the 5′ end first
        while (partner[outEnd] >= 0)
        {
            int entered = partner[outEnd];
            current = entered / 2;
            outEnd = entered ^ 1;                      // the other end of the contig just entered
            if (current == contig) break;              // defensive: paths are acyclic by construction
        }

        return current;
    }

    /// <summary>Reads a path from a terminal contig, oriented so its free end comes first.</summary>
    private static List<PlacedContig> ReadScaffoldPath(int terminal, int[] partner)
    {
        var path = new List<PlacedContig>();
        int entryEnd = partner[2 * terminal] < 0 ? 2 * terminal : 2 * terminal + 1;
        while (true)
        {
            int contig = entryEnd / 2;
            path.Add(new PlacedContig(contig, Reversed: (entryEnd & 1) == 1, entryEnd));
            int exitEnd = entryEnd ^ 1;
            if (partner[exitEnd] < 0) break;
            entryEnd = partner[exitEnd];
        }

        return path;
    }

    /// <summary>The same path read from the other strand: order reversed, every orientation flipped.</summary>
    private static List<PlacedContig> ReverseScaffoldPath(List<PlacedContig> path)
    {
        var reversed = new List<PlacedContig>(path.Count);
        for (int i = path.Count - 1; i >= 0; i--)
            reversed.Add(new PlacedContig(path[i].Contig, !path[i].Reversed, path[i].EntryEnd ^ 1));
        return reversed;
    }

    /// <summary>
    /// Spells a path: positive gap → that many gap characters; negative gap with an agreeing overlap
    /// → collapsed (ABySS MergeContigs); otherwise the AGP unknown-size gap of 100.
    /// </summary>
    private static string SpellScaffold(
        List<PlacedContig> path, IReadOnlyList<string> contigs, int[] joinGap, char gapCharacter)
    {
        var sb = new StringBuilder();
        int sinceGap = 0;                              // bases placed after the last gap run
        for (int i = 0; i < path.Count; i++)
        {
            string seq = contigs[path[i].Contig] ?? string.Empty;
            if (path[i].Reversed) seq = ReverseComplement(seq);

            if (i == 0)
            {
                sb.Append(seq);
                sinceGap = seq.Length;
                continue;
            }

            int gap = joinGap[path[i].EntryEnd];
            if (gap > 0)
            {
                sb.Append(gapCharacter, gap);
                sb.Append(seq);
                sinceGap = seq.Length;
            }
            else if (gap < 0 && -(long)gap <= Math.Min(sinceGap, seq.Length)
                     && OverlapAgrees(sb, seq, -gap))
            {
                int overlap = -gap;
                string tail = sb.ToString(sb.Length - overlap, overlap);
                sb.Length -= overlap;
                sb.Append(MergeContigs(tail, seq, overlap, resolveOverlapConflicts: true));
                sinceGap += seq.Length - overlap;
            }
            else
            {
                sb.Append(gapCharacter, UnknownGapLength);
                sb.Append(seq);
                sinceGap = seq.Length;
            }
        }

        return sb.ToString();
    }

    /// <summary>ABySS createConsensus agreement: equal ignoring case, or either base is N.</summary>
    private static bool OverlapAgrees(StringBuilder placed, string next, int overlap)
    {
        int start = placed.Length - overlap;
        for (int i = 0; i < overlap; i++)
        {
            char a = char.ToUpperInvariant(placed[start + i]);
            char b = char.ToUpperInvariant(next[i]);
            if (a != b && a != 'N' && b != 'N') return false;
        }

        return true;
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
    /// Calculates per-base sequencing coverage (read depth) along the reference, forward strand
    /// only. Equivalent to
    /// <see cref="CalculateCoverage(string, IReadOnlyList{string}, int, bool, bool)"/> with
    /// <c>includeReverseComplement: false, excludeMultiMapped: false</c>.
    /// </summary>
    /// <param name="reference">Reference sequence; depth is reported per position.</param>
    /// <param name="reads">Reads to map against the reference.</param>
    /// <param name="minOverlap">Minimum matching characters required to place a read (≥ 1).</param>
    /// <returns>Per-base depth array of length <c>reference.Length</c>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> or <paramref name="reads"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="reads"/> contains a null element.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minOverlap"/> is less than 1.</exception>
    public static int[] CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = DefaultMinOverlap) =>
        CalculateCoverage(reference, reads, minOverlap, includeReverseComplement: false, excludeMultiMapped: false);

    /// <summary>
    /// Calculates per-base sequencing coverage (read depth) along the reference.
    /// <para>
    /// <b>Placement.</b> Each read is aligned ungapped and end-to-end (the whole read must lie
    /// inside the reference, no soft clipping — Bowtie 2's default <c>--end-to-end</c> mode) at
    /// the window of minimum Hamming distance, found by the canonical
    /// <see cref="ApproximateMatcher.FindBestMatch(string, string)"/> (case-insensitive, leftmost on
    /// ties). With <paramref name="includeReverseComplement"/> the reverse complement of the read
    /// is also aligned (the reverse/Crick strand; Bowtie 2 and BWA search both strands by default,
    /// <c>--nofw/--norc</c> disable one); on an equal score the forward placement is kept. A read is
    /// placed only when its best placement has at least <paramref name="minOverlap"/> matching
    /// characters (length − Hamming distance); a read longer than the reference cannot align.
    /// Each read is placed at most once (one primary alignment per read, as reported by Bowtie 2
    /// in its default mode).
    /// </para>
    /// <para>
    /// <b>Depth.</b> A read placed at position p with length L increments every reference
    /// position in [p, p + L): per-base depth is the number of reads whose aligned span covers
    /// the position (<c>samtools depth</c>; Metagenomics Wiki, "SAMtools: get breadth of
    /// coverage"; Cook, "Calculate Depth and Breadth of Coverage From a bam File"). Placements
    /// are ungapped, so there are no deletions whose counting would depend on <c>samtools depth -J</c>.
    /// </para>
    /// <para>
    /// <b>Multi-mapping reads.</b> By default a read whose best score is reached at several
    /// loci counts once, at the forward-first, leftmost one (Bowtie 2 / BWA report one of the
    /// equally good hits, chosen pseudo-randomly). With <paramref name="excludeMultiMapped"/> such
    /// reads are skipped: Bowtie 2 assigns MAPQ ≤ 1 when the best and second-best alignment scores
    /// are equal (<c>unique.h</c>, <c>BowtieMapq2</c>), so this reproduces
    /// <c>samtools depth -Q 2</c>. Loci are distinct reference spans, so a reverse-palindromic read
    /// matching both strands at the same position is one locus.
    /// </para>
    /// </summary>
    /// <param name="reference">Reference sequence; depth is reported per position.</param>
    /// <param name="reads">Reads to map against the reference.</param>
    /// <param name="minOverlap">Minimum matching characters required to place a read (≥ 1).</param>
    /// <param name="includeReverseComplement">Also align each read's reverse complement (both strands).</param>
    /// <param name="excludeMultiMapped">Skip reads whose best score is attained at more than one locus.</param>
    /// <returns>
    /// An array of length <c>reference.Length</c> whose i-th element is the number of placed reads
    /// covering reference position i. Reads that fail to place contribute nothing.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> or <paramref name="reads"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="reads"/> contains a null element.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minOverlap"/> is less than 1.</exception>
    public static int[] CalculateCoverage(
        string reference,
        IReadOnlyList<string> reads,
        int minOverlap,
        bool includeReverseComplement,
        bool excludeMultiMapped = false)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(reads);
        if (minOverlap < 1)
            throw new ArgumentOutOfRangeException(nameof(minOverlap), minOverlap,
                "At least one matching character is required to place a read.");

        var coverage = new int[reference.Length];

        foreach (string read in reads)
        {
            if (read is null)
                throw new ArgumentException("Reads must not contain null elements.", nameof(reads));

            int pos = PlaceRead(reference, read, minOverlap, includeReverseComplement, excludeMultiMapped);
            if (pos < 0) continue;
            for (int i = pos; i < pos + read.Length; i++)
                coverage[i]++;
        }

        return coverage;
    }

    /// <summary>
    /// Start of the read's primary end-to-end ungapped placement, or -1 when it does not place
    /// (fewer than <paramref name="minOverlap"/> matches, or multi-mapped and excluded).
    /// </summary>
    private static int PlaceRead(
        string reference, string read, int minOverlap, bool bothStrands, bool excludeMultiMapped)
    {
        if (ApproximateMatcher.FindBestMatch(reference, read) is not { } forward)
            return -1; // empty read / reference, or read longer than reference

        string? reverse = bothStrands ? ReverseComplement(read) : null;
        var best = forward;
        if (reverse is not null)
        {
            var rc = ApproximateMatcher.FindBestMatch(reference, reverse)!.Value;
            if (rc.Distance < forward.Distance) best = rc; // equal score → keep forward strand
        }

        if (read.Length - best.Distance < minOverlap) return -1;

        if (excludeMultiMapped)
        {
            var loci = new HashSet<int>();
            foreach (var hit in ApproximateMatcher.FindWithMismatches(reference, read, best.Distance))
                loci.Add(hit.Position);
            if (reverse is not null)
            {
                foreach (var hit in ApproximateMatcher.FindWithMismatches(reference, reverse, best.Distance))
                    loci.Add(hit.Position);
            }

            if (loci.Count > 1) return -1;
        }

        return best.Position;
    }

    // Default plurality cut-off for emitting a committed residue. EMBOSS `cons` defines the
    // default plurality as a simple majority (half the total weight); we use 0.5 so the default
    // is a true majority vote. Biopython's `dumb_consensus` documented default (0.7) is reachable
    // by passing threshold: 0.7. Source: EMBOSS cons (plurality cut-off); Biopython dumb_consensus.
    private const double DefaultConsensusThreshold = 0.5;

    // Default IUPAC "any base" symbol emitted when a column has no committed residue. Biopython's
    // dumb_consensus default ambiguous symbol is 'X' (protein); for DNA/RNA the IUPAC degenerate
    // code for "any base" is 'N'. Source: Wikipedia "Consensus sequence" (IUPAC N = any base).
    private const char DefaultAmbiguousSymbol = 'N';

    // Gap symbols excluded from the per-column tally per Biopython dumb_consensus
    // (residue counted only when it is neither '-' nor '.').
    private const char GapDash = '-';
    private const char GapDot = '.';

    /// <summary>
    /// Computes the consensus sequence from a set of aligned reads using the column-wise
    /// majority/threshold rule of Biopython's <c>dumb_consensus</c>: for each alignment column,
    /// non-gap residues ('-' and '.' are skipped) are tallied; a residue is emitted only when a
    /// single residue holds the strict maximum count AND its frequency among the non-gap residues
    /// is at least <paramref name="threshold"/>. Ties for the maximum count, sub-threshold majorities,
    /// and all-gap columns emit <paramref name="ambiguous"/>. The consensus length equals the longest
    /// input read (the full alignment length); reads shorter than a column contribute nothing there.
    /// </summary>
    /// <remarks>
    /// Reference: Biopython <c>Bio.Align.AlignInfo.SummaryInfo.dumb_consensus</c> (v1.79),
    /// https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py
    /// (decision rule, gap skipping, tie→ambiguous, alignment-length consensus);
    /// EMBOSS <c>cons</c> (plurality cut-off below which there is no consensus).
    /// </remarks>
    /// <param name="alignedReads">Pre-aligned reads (equal columns; ragged lengths allowed).</param>
    /// <param name="threshold">
    /// Minimum frequency (max count / non-gap count in the column) required to commit a residue.
    /// Default <see cref="DefaultConsensusThreshold"/> (simple majority); pass 0.7 to reproduce
    /// Biopython's documented default.
    /// </param>
    /// <param name="ambiguous">Symbol emitted when a column has no committed residue.</param>
    /// <returns>The consensus string, length equal to the longest read.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedReads"/> is null.</exception>
    public static string ComputeConsensus(
        IReadOnlyList<string> alignedReads,
        double threshold = DefaultConsensusThreshold,
        char ambiguous = DefaultAmbiguousSymbol)
    {
        ArgumentNullException.ThrowIfNull(alignedReads);

        if (alignedReads.Count == 0) return "";

        // Consensus spans the full alignment length = the longest read (Biopython: con_len).
        int length = 0;
        foreach (string read in alignedReads)
        {
            if (read.Length > length) length = read.Length;
        }

        var sb = new StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            int numAtoms = 0; // non-gap residues contributing to this column (Biopython: num_atoms)

            foreach (string read in alignedReads)
            {
                if (pos >= read.Length) continue;

                char c = char.ToUpperInvariant(read[pos]);
                if (c == GapDash || c == GapDot) continue;

                counts[c] = counts.GetValueOrDefault(c, 0) + 1;
                numAtoms++;
            }

            // Determine the residue(s) with the maximum count (Biopython: max_atoms / max_size).
            int maxSize = 0;
            int maxCount = 0; // how many residues share the maximum (>1 means a tie)
            char maxResidue = ambiguous;

            foreach (KeyValuePair<char, int> kvp in counts)
            {
                if (kvp.Value > maxSize)
                {
                    maxSize = kvp.Value;
                    maxCount = 1;
                    maxResidue = kvp.Key;
                }
                else if (kvp.Value == maxSize)
                {
                    maxCount++;
                }
            }

            // Emit the residue only when exactly one residue holds the max AND it meets the
            // threshold among non-gap residues; otherwise emit the ambiguous symbol. The
            // numAtoms > 0 guard mirrors Biopython's len(max_atoms)==1 short-circuit, which
            // prevents division by zero on all-gap columns.
            bool committed = maxCount == 1
                             && numAtoms > 0
                             && (double)maxSize / numAtoms >= threshold;

            sb.Append(committed ? maxResidue : ambiguous);
        }

        return sb.ToString();
    }

    // Sanger/Phred+33 ASCII offset: quality char Q decodes to Phred value (Q - 33).
    // Cock et al. (2010) NAR 38(6):1767–1771; BWA bwa_trim_read decodes as (qual - 33).
    private const int PhredAsciiOffset = 33;

    /// <summary>
    /// Quality-trims reads using the BWA / cutadapt running-sum algorithm, then drops
    /// reads shorter than <paramref name="minLength"/>.
    /// </summary>
    /// <remarks>
    /// For each read the running-sum method is applied to the 3' end and then the 5' end:
    /// subtract the quality cutoff from every Phred score, compute the partial sum from
    /// each index to that end, and cut at the index where the partial sum is minimal.
    /// This is the algorithm used by BWA (<c>bwa_trim_read</c>) and cutadapt; it removes
    /// low-quality bases from the ends while allowing some high-quality bases among the
    /// bad ones. Qualities are decoded as Phred+33 (Sanger encoding).
    /// A cutoff &lt; 1 disables trimming (BWA <c>trim_qual &lt; 1</c> guard).
    /// </remarks>
    /// <param name="reads">Reads as (sequence, Phred+33 quality string) pairs; the two strings must be equal length per read.</param>
    /// <param name="minQuality">Quality cutoff subtracted from each Phred score. Values &lt; 1 disable trimming.</param>
    /// <param name="minLength">Reads whose trimmed length is below this are dropped.</param>
    /// <returns>The surviving trimmed sequences, in input order.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="reads"/> is null.</exception>
    public static IReadOnlyList<string> QualityTrimReads(
        IReadOnlyList<(string sequence, string quality)> reads,
        int minQuality = 20,
        int minLength = 50)
    {
        ArgumentNullException.ThrowIfNull(reads);

        var trimmed = new List<string>(reads.Count);

        foreach (var (sequence, quality) in reads)
        {
            int n = sequence.Length;
            int start = 0;
            int end = n;

            // BWA: a cutoff below 1 disables trimming (bwa_trim_read returns 0).
            if (minQuality >= 1)
            {
                // cutadapt quality_trim_index: the 3' and 5' passes are run independently
                // over the FULL read (not on each other's surviving window), then the read
                // is dropped entirely if the two windows cross (start >= stop).
                // 3' end: cut at the index minimizing the partial sum of (q - cutoff) from
                // that index to the 3' end — keep [start, end).
                end = TrimEnd(quality, n, minQuality);
                // 5' end: mirror the procedure from the 5' end over the full read.
                start = TrimStart(quality, n, minQuality);

                // cutadapt: if start >= stop the good-quality segment is empty.
                if (start >= end)
                {
                    start = 0;
                    end = 0;
                }
            }

            if (end - start >= minLength)
            {
                trimmed.Add(sequence.Substring(start, end - start));
            }
        }

        return trimmed;
    }

    /// <summary>
    /// Running-sum 3'-end cut over <c>quality[0..length)</c>: returns the new exclusive upper
    /// bound (the index of the minimal partial sum of (Phred − cutoff) accumulated from the 3'
    /// end). Scanning stops once the partial sum becomes positive — the cutadapt/BWA
    /// <c>s &lt; 0</c> early break (sign-flipped here), which is what lets a few good-quality
    /// bases survive inside a low-quality tail.
    /// </summary>
    private static int TrimEnd(string quality, int length, int cutoff)
    {
        int sum = 0;
        int min = 0;
        int cut = length; // default: keep everything (minimum at the end position)

        for (int i = length - 1; i >= 0; i--)
        {
            sum += (quality[i] - PhredAsciiOffset) - cutoff;
            if (sum > 0)
            {
                break; // cutadapt bwa_trim_read: stop when the accumulated (cutoff - q) < 0.
            }

            if (sum < min)
            {
                min = sum;
                cut = i;
            }
        }

        return cut;
    }

    /// <summary>
    /// Running-sum 5'-end cut over <c>quality[0..length)</c>: returns the new inclusive lower
    /// bound (one past the index of the minimal partial sum of (Phred − cutoff) accumulated from
    /// the 5' end). Scanning stops once the partial sum becomes positive — the cutadapt/BWA
    /// <c>s &lt; 0</c> early break (sign-flipped here).
    /// </summary>
    private static int TrimStart(string quality, int length, int cutoff)
    {
        int sum = 0;
        int min = 0;
        int cut = 0; // default: keep everything (minimum at the start position)

        for (int i = 0; i < length; i++)
        {
            sum += (quality[i] - PhredAsciiOffset) - cutoff;
            if (sum > 0)
            {
                break; // cutadapt bwa_trim_read: stop when the accumulated (cutoff - q) < 0.
            }

            if (sum < min)
            {
                min = sum;
                cut = i + 1;
            }
        }

        return cut;
    }

    // Canonical DNA alphabet over which substitution candidates are enumerated, in fixed
    // A,C,G,T order so the correction is deterministic. Quake/Musket correct only the four
    // standard nucleotides. Source: Kelley et al. (2010) Genome Biol 11:R116 (nucleotide edits);
    // Liu et al. (2013) Musket, Bioinformatics 29(3):308-315 (substitution model).
    private static readonly char[] DnaBases = { 'A', 'C', 'G', 'T' };

    // Default k-mer length. Quake/Musket leave k user-selectable; the de Bruijn default in this
    // class is 31, but error correction uses smaller k so that error k-mers stay rare and solid
    // k-mers are well sampled at typical read depth. We follow Quake's recommended order of
    // magnitude (k chosen so the genome is well covered); k is exposed as a parameter so callers
    // set it from their data. Source: Kelley et al. (2010) Genome Biol 11:R116 (k selection).
    private const int DefaultErrorCorrectKmerSize = 15;

    // Default solidity (trusted) coverage cut-off: a k-mer is "trusted/solid" when its
    // multiplicity is at least this value, "untrusted/weak" otherwise. Quake and Musket choose
    // this cut-off automatically from the coverage histogram valley; we expose it as a parameter
    // (callers pass the data-driven cut-off) and default to 2 so a k-mer must be seen more than
    // once to be trusted. Source: Liu et al. (2013) Musket (multiplicity > cut-off => trusted);
    // Kelley et al. (2010) Quake (high-coverage k-mers trusted, low-coverage untrusted).
    private const int DefaultMinKmerFrequency = 2;

    /// <summary>
    /// Corrects substitution errors in reads using the k-mer spectrum (two-sided) method of
    /// Musket / Quake. A k-mer is <em>trusted</em> (solid) when its multiplicity across all reads
    /// is at least <paramref name="minKmerFrequency"/> and <em>untrusted</em> (weak) otherwise.
    /// For every read position covered only by untrusted k-mers, the method looks for a single
    /// alternative base that makes <em>every</em> k-mer covering that position trusted; the base is
    /// changed only when that alternative is <em>unique</em> — if zero or more than one alternative
    /// works, the base is left unchanged (ambiguity is not resolved). The k-mer frequency table is
    /// computed once from the original reads and is not updated during correction.
    /// </summary>
    /// <remarks>
    /// Algorithm: Liu, Schmidt &amp; Maskell (2013), "Musket: a multistage k-mer spectrum-based
    /// error corrector for Illumina sequence data", Bioinformatics 29(3):308-315,
    /// DOI 10.1093/bioinformatics/bts690 — two-sided correction: "a base covered by any trusted
    /// k-mer is trusted"; find a unique alternative base making the k-mers covering position i
    /// trusted; "if more than one alternative is found ... the base will keep unchanged as a
    /// result of ambiguity". Kelley, Schatz &amp; Salzberg (2010), "Quake", Genome Biol 11:R116,
    /// DOI 10.1186/gb-2010-11-11-r116 — trusted = high-coverage k-mers, untrusted = low-coverage;
    /// corrections are single-nucleotide edits over a region of untrusted k-mers.
    /// Corrections are restricted to single-base substitutions (no insertions/deletions), matching
    /// the substitution model of both tools.
    /// </remarks>
    /// <param name="reads">Reads to correct. Compared case-insensitively (upper-cased internally).</param>
    /// <param name="kmerSize">k-mer length k (k ≥ 1). Reads shorter than k contribute no k-mers and are returned unchanged.</param>
    /// <param name="minKmerFrequency">Trusted cut-off: k-mers with multiplicity ≥ this value are trusted.</param>
    /// <returns>The corrected reads, upper-cased, in input order; length and count are preserved.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="reads"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="kmerSize"/> &lt; 1.</exception>
    public static IReadOnlyList<string> ErrorCorrectReads(
        IReadOnlyList<string> reads,
        int kmerSize = DefaultErrorCorrectKmerSize,
        int minKmerFrequency = DefaultMinKmerFrequency)
    {
        ArgumentNullException.ThrowIfNull(reads);
        if (kmerSize < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(kmerSize), kmerSize, "k-mer size must be at least 1.");
        }

        // Stage 1: build the k-mer spectrum (multiplicity table) from the original reads once.
        var kmerCounts = BuildKmerSpectrum(reads, kmerSize);

        var corrected = new List<string>(reads.Count);

        foreach (string read in reads)
        {
            corrected.Add(CorrectRead(read, kmerSize, minKmerFrequency, kmerCounts));
        }

        return corrected;
    }

    /// <summary>
    /// Builds the k-mer multiplicity table (k-mer spectrum) over all reads, counting overlapping
    /// k-mers case-insensitively. Reads shorter than k contribute nothing.
    /// </summary>
    private static Dictionary<string, int> BuildKmerSpectrum(IReadOnlyList<string> reads, int kmerSize)
    {
        var kmerCounts = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (string read in reads)
        {
            if (read is null)
            {
                continue;
            }

            for (int i = 0; i + kmerSize <= read.Length; i++)
            {
                string kmer = read.Substring(i, kmerSize).ToUpperInvariant();
                kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer, 0) + 1;
            }
        }

        return kmerCounts;
    }

    /// <summary>
    /// Applies the Musket two-sided substitution rule to one read against the global k-mer
    /// spectrum: for each position covered only by untrusted k-mers, substitute the unique base
    /// (if any) that makes every k-mer covering that position trusted.
    /// </summary>
    private static string CorrectRead(
        string read,
        int kmerSize,
        int minKmerFrequency,
        Dictionary<string, int> kmerCounts)
    {
        // Reads with no k-mers (shorter than k) cannot be evaluated; return upper-cased copy.
        if (read.Length < kmerSize)
        {
            return read.ToUpperInvariant();
        }

        var sb = new StringBuilder(read.ToUpperInvariant());

        for (int pos = 0; pos < sb.Length; pos++)
        {
            if (IsPositionTrusted(sb, pos, kmerSize, minKmerFrequency, kmerCounts))
            {
                continue; // Base is covered by a trusted k-mer -> leave it unchanged.
            }

            char original = sb[pos];
            char unique = '\0';
            int candidateCount = 0;

            foreach (char candidate in DnaBases)
            {
                if (candidate == original)
                {
                    continue;
                }

                sb[pos] = candidate;
                if (AllCoveringKmersTrusted(sb, pos, kmerSize, minKmerFrequency, kmerCounts))
                {
                    candidateCount++;
                    unique = candidate;
                }
            }

            // Apply only an unambiguous correction; otherwise restore the original base.
            sb[pos] = candidateCount == 1 ? unique : original;
        }

        return sb.ToString();
    }

    /// <summary>
    /// True when position <paramref name="pos"/> is covered by at least one trusted k-mer
    /// (multiplicity ≥ <paramref name="minKmerFrequency"/>) — the Musket "trusted base" rule.
    /// </summary>
    private static bool IsPositionTrusted(
        StringBuilder sb,
        int pos,
        int kmerSize,
        int minKmerFrequency,
        Dictionary<string, int> kmerCounts)
    {
        foreach (int start in CoveringKmerStarts(sb.Length, pos, kmerSize))
        {
            if (kmerCounts.GetValueOrDefault(KmerAt(sb, start, kmerSize), 0) >= minKmerFrequency)
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// True when EVERY k-mer covering position <paramref name="pos"/> is trusted. Used to test a
    /// candidate substitution: a correction is valid only if it makes all covering k-mers solid.
    /// </summary>
    private static bool AllCoveringKmersTrusted(
        StringBuilder sb,
        int pos,
        int kmerSize,
        int minKmerFrequency,
        Dictionary<string, int> kmerCounts)
    {
        foreach (int start in CoveringKmerStarts(sb.Length, pos, kmerSize))
        {
            if (kmerCounts.GetValueOrDefault(KmerAt(sb, start, kmerSize), 0) < minKmerFrequency)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>
    /// Enumerates the start indices of all length-k windows that cover read position
    /// <paramref name="pos"/> within a read of length <paramref name="length"/>.
    /// </summary>
    private static IEnumerable<int> CoveringKmerStarts(int length, int pos, int kmerSize)
    {
        int first = Math.Max(0, pos - kmerSize + 1);
        int last = Math.Min(pos, length - kmerSize);
        for (int start = first; start <= last; start++)
        {
            yield return start;
        }
    }

    /// <summary>Reads the length-k window starting at <paramref name="start"/> from the buffer.</summary>
    private static string KmerAt(StringBuilder sb, int start, int kmerSize)
    {
        char[] buffer = new char[kmerSize];
        sb.CopyTo(start, buffer, 0, kmerSize);
        return new string(buffer);
    }
}

