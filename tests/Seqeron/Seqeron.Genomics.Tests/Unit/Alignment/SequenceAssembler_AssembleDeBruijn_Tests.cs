// ASSEMBLY-DBG-001 — De Bruijn graph assembly
// Evidence: docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-DBG-001.md
// Source: Compeau & Pevzner, Bioinformatics Algorithms ch.3 (Rosalind BA3K/BA3M:
//         contigs = maximal non-branching paths, incl. isolated cycles); Zerbino & Birney (2008),
//         Genome Res 18:821 (Velvet: distinct k-mers, chain simplification, reverse-complement
//         twins); Langmead B., "De Bruijn Graph assembly" (JHU lecture notes), p.5-28;
//         Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987-991, DOI 10.1038/nbt.2023.

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_AssembleDeBruijn_Tests
{
    // Langmead DBG p.19/22 example sentence.
    private const string ToEvery = "to_every_thing_turn_turn_turn_there_is_a_season";

    // Flattens the (k-1)-mer out-adjacency multigraph into an ordered (from, to) edge list.
    private static List<(string From, string To)> Edges(Dictionary<string, List<string>> g) =>
        g.SelectMany(kvp => kvp.Value.Select(to => (From: kvp.Key, To: to)))
         .OrderBy(e => e.From, StringComparer.Ordinal)
         .ThenBy(e => e.To, StringComparer.Ordinal)
         .ToList();

    #region BuildDeBruijnGraph

    // M1 — k=3 graph of "AAABBBA": nodes {AA,AB,BB,BA}, edges AA->AA, AA->AB, AB->BB, BB->BB, BB->BA.
    // Source: Langmead DBG p.5-11.
    [Test]
    public void BuildDeBruijnGraph_AaabbbaK3_ProducesExpectedNodesAndEdges()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AAABBBA" }, 3);

        var nodes = new HashSet<string>(graph.Keys);
        foreach (var succ in graph.Values)
            foreach (var s in succ) nodes.Add(s);

        var expectedEdges = new List<(string, string)>
        {
            ("AA", "AA"), ("AA", "AB"), ("AB", "BB"), ("BB", "BA"), ("BB", "BB"),
        };

        Assert.Multiple(() =>
        {
            Assert.That(nodes.OrderBy(n => n, StringComparer.Ordinal),
                Is.EqualTo(new[] { "AA", "AB", "BA", "BB" }),
                "The (k-1)-mer node set of AAABBBA at k=3 is exactly {AA, AB, BB, BA} (Langmead DBG p.6).");
            Assert.That(Edges(graph), Is.EqualTo(expectedEdges),
                "INV-01: each edge runs from a k-mer's (k-1)-prefix to its (k-1)-suffix (Langmead DBG p.7).");
        });
    }

    // M2 — A repeated k-mer produces a multiedge: AAABBBBA at k=3 has two BB->BB edges.
    // Source: Langmead DBG p.8.
    [Test]
    public void BuildDeBruijnGraph_RepeatedKmer_ProducesMultiedge()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AAABBBBA" }, 3);

        int bbToBb = graph.TryGetValue("BB", out var succ)
            ? succ.Count(s => s == "BB")
            : 0;

        Assert.Multiple(() =>
        {
            Assert.That(bbToBb, Is.EqualTo(2),
                "AAABBBBA contributes the 3-mer BBB twice, so the node BB has two BB->BB edges (multiedge, Langmead DBG p.8).");
            Assert.That(Edges(graph).Count, Is.EqualTo(6),
                "INV-02: AAABBBBA (length 8) yields 8-3+1 = 6 k-mers, hence 6 directed edges.");
        });
    }

    // M3 — INV-02: total #edges equals total #k-mers chopped (single read).
    // Source: Langmead DBG p.7.
    [Test]
    public void BuildDeBruijnGraph_EdgeCount_EqualsNumberOfKmers()
    {
        const string read = "AAABBBA"; // length 7, k=3 -> 5 k-mers
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { read }, 3);

        Assert.That(Edges(graph).Count, Is.EqualTo(read.Length - 3 + 1),
            "INV-02: the number of directed edges equals the number of k-mers (7-3+1 = 5).");
    }

    // C2 — INV-02 property across multiple reads: edge count equals the summed k-mer count.
    // Source: Langmead DBG p.7.
    [Test]
    public void BuildDeBruijnGraph_MultiRead_EdgeCountEqualsKmerCount()
    {
        var reads = new[] { "ATGGCGT", "GGCGTGCA", "ACGTACGT" };
        const int k = 4;
        int expectedKmers = reads.Sum(r => r.Length - k + 1);

        var graph = SequenceAssembler.BuildDeBruijnGraph(reads, k);

        Assert.That(Edges(graph).Count, Is.EqualTo(expectedKmers),
            "INV-02: total edges across reads equals the total number of k-mers chopped (multigraph).");
    }

    // S1 — A read shorter than k contributes no k-mers (INV-03). Source: Langmead DBG p.16.
    [Test]
    public void BuildDeBruijnGraph_ReadShorterThanK_ProducesEmptyGraph()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AC" }, 5);

        Assert.That(graph.Count, Is.EqualTo(0),
            "INV-03: a read of length 2 < k=5 yields no k-mers, so the graph is empty.");
    }

    // S2 — Empty read list yields an empty graph. Trivial identity.
    [Test]
    public void BuildDeBruijnGraph_EmptyReads_ProducesEmptyGraph()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(Array.Empty<string>(), 3);

        Assert.That(graph.Count, Is.EqualTo(0), "No reads -> no nodes or edges.");
    }

    // C1 — A repeated (k-1)-mer creates a branch node (out-degree >= 2), the structural cause
    // of multiple Eulerian walks / unresolvable repeats. Source: Langmead DBG p.21-22.
    [Test]
    public void BuildDeBruijnGraph_ToEveryK3_HasBranchNode()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { ToEvery }, 3);

        Assert.That(graph.Values.Any(succ => succ.Count >= 2), Is.True,
            "A repeated (k-1)-mer at k=3 produces a node with out-degree >= 2 (branching => multiple Eulerian walks, Langmead DBG p.21).");
    }

    // Edge-case — k < 2 is rejected (nodes would be empty). Implementation contract.
    [Test]
    public void BuildDeBruijnGraph_KLessThanTwo_Throws()
    {
        Assert.That(() => SequenceAssembler.BuildDeBruijnGraph(new[] { "AAAA" }, 1),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(),
            "k must be >= 2 so that (k-1)-mer nodes are non-empty.");
    }

    // Edge-case — null read list yields an empty graph (no exception). Implementation contract.
    [Test]
    public void BuildDeBruijnGraph_NullReads_ProducesEmptyGraph()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(null!, 3);

        Assert.That(graph.Count, Is.EqualTo(0), "Null reads are treated as empty -> no nodes or edges.");
    }

    #endregion

    #region AssembleDeBruijn — contigs = maximal non-branching paths of the distinct-k-mer graph

    // All expected contig sets below were computed by two independent Python oracles (scratch
    // dbg_oracle.py, recorded in the Evidence): (A) the Compeau & Pevzner BA3M
    // MaximalNonBranchingPaths pseudocode over DeBruijn(distinct k-mers), and (B) Velvet's
    // Simplification rule (merge node A -> B when A has one outgoing and B one ingoing arc) over
    // the k-mer node graph. A and B agree on every case.

    private static SequenceAssembler.AssemblyParameters K(int k, int minContigLength = 1) =>
        new(KmerSize: k, MinContigLength: minContigLength);

    private static string Spell(IReadOnlyList<string> path) =>
        path[0] + string.Concat(path.Skip(1).Select(n => n[^1]));

    // M4 — AAABBBA, k=3: AA and BB carry self-loops, so neither is 1-in-1-out and the graph
    // splits into the unitigs AAA, AABB, BBA, BBB (BA3M; Velvet merge rule).
    [Test]
    public void AssembleDeBruijn_Aaabbba_EmitsMaximalNonBranchingPaths()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "AAABBBA" }, K(3));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "AAA", "AABB", "BBA", "BBB" }),
            "Contigs are the spelled maximal non-branching paths, starting nodes in ordinal order (BA3K/BA3M).");
    }

    // M5 — a_long_long_long_time, k=5: the repeated 4-mer nodes "_lon"/"ong_" branch, so the
    // repeat is not traversed three times; the contigs are the four unitigs (upper-cased).
    [Test]
    public void AssembleDeBruijn_ALongRepeat_SplitsAtRepeatBranches()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "a_long_long_long_time" }, K(5));

        Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal),
            Is.EqualTo(new[] { "A_LON", "ONG_LON", "ONG_TIME", "_LONG_" }),
            "A repeat longer than k-2 creates branch nodes; contigs stop there (BA3K: contigs must be present in any assembly).");
    }

    // M6 — to_every…, k=4: the "turn" repeat yields six unitigs.
    [Test]
    public void AssembleDeBruijn_ToEveryK4_SplitsAtTurnRepeat()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { ToEvery }, K(4));

        Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal),
            Is.EqualTo(new[] { "N_TH", "N_TU", "TO_EVERY_TH", "_THERE_IS_A_SEASON", "_THING_TU", "_TURN_T" }));
    }

    // M7 — A repeat-free DNA string is one non-branching path, reconstructed verbatim.
    [Test]
    public void AssembleDeBruijn_DnaUniqueWalk_ReconstructsInput()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ATGGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // F-DBG-1 — Overlapping reads share k-mers; coverage must not add edges (Velvet: "high
    // redundancy is naturally handled by the graph without affecting the number of nodes").
    // Before the fix the shared 4-mers GGCG/GCGT became parallel edges and the Eulerian walk
    // spelled a chimera.
    [Test]
    public void AssembleDeBruijn_OverlappingReads_AssembleIntoGenome()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ATGGCGT", "GGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // F-DBG-1 — Two copies of the same read (coverage 2) assemble to that read; before the fix the
    // doubled edges produced "ACGGTT", which contains 3-mers (CGG, GGT, GTT) absent from the input.
    [Test]
    public void AssembleDeBruijn_DuplicateReads_DoNotDoubleTheGenome()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ACGT", "ACGT" }, K(3));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGT" }));
    }

    // F-DBG-1 — MCP example reads AAABBB / AABBBC at k=3: former output "AAABBBBBBC" (BBB walked
    // three times); the distinct-k-mer unitigs are AAA, AABB, BBB, BBC.
    [Test]
    public void AssembleDeBruijn_SharedKmersAcrossReads_CountedOnce()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "AAABBB", "AABBBC" }, K(3));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "AAA", "AABB", "BBB", "BBC" }));
    }

    // Isolated cycle — every node is 1-in-1-out; BA3M emits the cycle (from its ordinally
    // smallest node ACG back to ACG): ACGTTACG.
    [Test]
    public void AssembleDeBruijn_IsolatedCycle_IsEmitted()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ACGTTACG" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTACG" }),
            "BA3M: 'MaximalNonBranchingPaths finds all isolated cycles in the graph'.");
    }

    // Cycle attached to a branch: ACGTACGTACGT + ACGTACGTTTTT at k=5 -> ACGTACGT (the loop from the
    // branching node ACGT back to itself), ACGTTTT and the TTTT self-loop TTTTT.
    [Test]
    public void AssembleDeBruijn_LoopThroughBranchNode_EmittedAsPath()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ACGTACGTACGT", "ACGTACGTTTTT" }, K(5));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTACGT", "ACGTTTT", "TTTTT" }));
    }

    // INV — every distinct input k-mer lies in exactly one contig, and no contig contains a k-mer
    // that is absent from the input (the unitigs partition the k-mer set; Chikhi et al. 2016).
    [Test]
    public void AssembleDeBruijn_ToEveryK3_ContigsPartitionTheKmerSet()
    {
        const int k = 3;
        string upper = ToEvery.ToUpperInvariant();
        var input = Enumerable.Range(0, upper.Length - k + 1).Select(i => upper.Substring(i, k))
            .Distinct().OrderBy(s => s, StringComparer.Ordinal).ToList();

        var contigs = SequenceAssembler.AssembleDeBruijn(new[] { ToEvery }, K(k)).Contigs;
        var emitted = contigs.SelectMany(c => Enumerable.Range(0, c.Length - k + 1).Select(i => c.Substring(i, k)))
            .OrderBy(s => s, StringComparer.Ordinal).ToList();

        Assert.That(emitted, Is.EqualTo(input), "each distinct k-mer is spelled exactly once over all contigs");
    }

    // Case — k-mers are counted case-insensitively by the canonical counter (CountKmersSpan);
    // contigs are upper-case.
    [Test]
    public void AssembleDeBruijn_MixedCaseReads_JoinCaseInsensitively()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "atggcGT", "GGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // Reads shorter than k contribute nothing but are still counted in TotalReads.
    [Test]
    public void AssembleDeBruijn_ReadShorterThanK_IsSkipped()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "AC", "ATGGCGTGCA" }, K(4));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
            Assert.That(result.TotalReads, Is.EqualTo(2));
        });
    }

    // Double-stranded option — the second read is the reverse complement of the genome's 3' half
    // (genome GATTACAGGTCCTTA). Forward-only gives two contigs; with reverse complements the
    // strands join (Velvet: each k-mer recorded with its reverse complement) and the pair
    // {genome, revcomp} is reported once as the ordinally smaller string.
    [Test]
    public void AssembleDeBruijn_ReverseComplementOption_JoinsOppositeStrandReads()
    {
        var reads = new[] { "GATTACAGGT", "TAAGGACCTG" };

        var forward = SequenceAssembler.AssembleDeBruijn(reads, K(5));
        var both = SequenceAssembler.AssembleDeBruijn(reads, K(5), includeReverseComplement: true);

        Assert.Multiple(() =>
        {
            Assert.That(forward.Contigs, Is.EqualTo(new[] { "GATTACAGGT", "TAAGGACCTG" }));
            Assert.That(both.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
        });
    }

    // Double-stranded option on a disconnected pair: one contig per strand pair, canonical orientation.
    [Test]
    public void AssembleDeBruijn_ReverseComplementOption_ReportsEachStrandPairOnce()
    {
        var both = SequenceAssembler.AssembleDeBruijn(new[] { "GATTACAGGTCCTTA" }, K(5), includeReverseComplement: true);
        var fromReverse = SequenceAssembler.AssembleDeBruijn(new[] { "TAAGGACCTGTAATC" }, K(5), includeReverseComplement: true);

        Assert.Multiple(() =>
        {
            Assert.That(both.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
            Assert.That(fromReverse.Contigs, Is.EqualTo(both.Contigs), "either strand yields the same canonical contig");
        });
    }

    // S3 — INV-05: every input k-mer is a substring of the reconstructed contig.
    [Test]
    public void AssembleDeBruijn_DnaUniqueWalk_EveryKmerIsSubstringOfContig()
    {
        const string genome = "ATGGCGTGCA";
        const int k = 4;
        var kmers = Enumerable.Range(0, genome.Length - k + 1).Select(i => genome.Substring(i, k)).ToList();

        string contig = SequenceAssembler.AssembleDeBruijn(new[] { genome }, K(k)).Contigs.Single();

        Assert.That(kmers.All(km => contig.Contains(km)), Is.True);
    }

    // S4 — INV-06: result statistics are consistent with the emitted contigs.
    [Test]
    public void AssembleDeBruijn_ALongRepeat_StatisticsAreConsistent()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "a_long_long_long_time" }, K(5));

        Assert.Multiple(() =>
        {
            Assert.That(result.TotalLength, Is.EqualTo(5 + 7 + 8 + 6), "A_LON + ONG_LON + ONG_TIME + _LONG_");
            Assert.That(result.LongestContig, Is.EqualTo(8));
            Assert.That(result.TotalReads, Is.EqualTo(1));
        });
    }

    // Two reads sharing no (k-1)-mer -> two disconnected simple paths -> two contigs.
    [Test]
    public void AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ATGGCGTGCA", "GATTACAGGTC" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA", "GATTACAGGTC" }));
    }

    // MinContigLength filter — the length-10 contig is dropped at MinContigLength = 11.
    [Test]
    public void AssembleDeBruijn_MinContigLength_FiltersShortContigs()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "ATGGCGTGCA", "GATTACAGGTC" }, K(4, 11));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "GATTACAGGTC" }));
    }

    // F-DBG-2 — null parameters use DefaultParameters (k = 31, MinContigLength 100). Before the fix
    // `new AssemblyParameters()` zero-initialised KmerSize and the call threw. Two overlapping reads
    // (0-80, 50-120) of a 120-bp genome with unique 30-mers assemble to the genome; a 99-bp
    // contig is filtered out.
    [Test]
    public void AssembleDeBruijn_NullParameters_UseDocumentedDefaults()
    {
        const string genome =
            "CGCGAGAGCCAACGTTAAGTCTGCCAACATTTTTCATCCTATCGCAAGTCGACTAGTAGTGTTCGATAGCCCTTTCCTAGTTGCGCCACGGCGAACAGCTGGGACGCACGGCGTATTTTC";

        var assembled = SequenceAssembler.AssembleDeBruijn(new[] { genome[..80], genome[50..] });
        var tooShort = SequenceAssembler.AssembleDeBruijn(new[] { genome[..99] });

        Assert.Multiple(() =>
        {
            Assert.That(assembled.Contigs, Is.EqualTo(new[] { genome }));
            Assert.That(tooShort.Contigs, Is.Empty, "a 99-bp contig is below the default MinContigLength 100");
        });
    }

    // k < 2 is rejected up front, even for an empty read set.
    [TestCase(1)]
    [TestCase(0)]
    [TestCase(-3)]
    public void AssembleDeBruijn_KmerSizeBelowTwo_Throws(int k)
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => SequenceAssembler.AssembleDeBruijn(new[] { "ACGT" }, K(k)),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.AssembleDeBruijn(Array.Empty<string>(), K(k)),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        });
    }

    // M8 — Empty read set returns an empty AssemblyResult. ASSUMPTION-2.
    [Test]
    public void AssembleDeBruijn_EmptyReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new List<string>());

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(0));
            Assert.That(result.TotalReads, Is.EqualTo(0));
            Assert.That(result.TotalLength, Is.EqualTo(0));
        });
    }

    // M9 — Null read set is handled like empty (no exception). Contract (mirrors OLC).
    [Test]
    public void AssembleDeBruijn_NullReads_ReturnsEmptyResult()
    {
        Assert.That(SequenceAssembler.AssembleDeBruijn(null!).Contigs.Count, Is.EqualTo(0));
    }

    // Determinism and read-order independence: the graph is built from the ordinally sorted
    // distinct k-mer set, so permuting the reads cannot change the contigs.
    [Test]
    public void AssembleDeBruijn_ReadOrder_DoesNotChangeContigs()
    {
        var reads = new[] { "GATCGAAACCC", "AAACCCGATCC", "GATCCTTTAAA" };

        var first = SequenceAssembler.AssembleDeBruijn(reads, K(4)).Contigs;
        var reversed = SequenceAssembler.AssembleDeBruijn(reads.Reverse().ToList(), K(4)).Contigs;

        Assert.That(reversed, Is.EqualTo(first));
    }

    #endregion

    #region FindMaximalNonBranchingPaths

    // Rosalind BA3M sample (Compeau & Pevzner ch.3): 1->2, 2->3, 3->4,5, 6->7, 7->6 ->
    // 1->2->3, 3->4, 3->5 and the isolated cycle 6->7->6 (the sample prints the same cycle as 7->6->7).
    [Test]
    public void FindMaximalNonBranchingPaths_Ba3mSample()
    {
        var graph = new Dictionary<string, List<string>>
        {
            ["1"] = new() { "2" }, ["2"] = new() { "3" }, ["3"] = new() { "4", "5" },
            ["6"] = new() { "7" }, ["7"] = new() { "6" },
        };

        var paths = SequenceAssembler.FindMaximalNonBranchingPaths(graph);

        Assert.That(paths.Select(p => string.Join("->", p)),
            Is.EqualTo(new[] { "1->2->3", "3->4", "3->5", "6->7->6" }));
    }

    // Rosalind BA3K sample (Contig Generation): k-mers ATG ATG TGT TGG CAT GGA GAT AGA ->
    // contigs AGA ATG ATG CAT GAT TGGA TGT. The textbook multigraph keeps the repeated ATG as two
    // parallel edges, each its own path.
    [Test]
    public void FindMaximalNonBranchingPaths_Ba3kSample_OverKmerMultigraph()
    {
        var patterns = new[] { "ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA" };
        var graph = SequenceAssembler.BuildDeBruijnGraph(patterns, 3);

        var contigs = SequenceAssembler.FindMaximalNonBranchingPaths(graph).Select(Spell)
            .OrderBy(c => c, StringComparer.Ordinal);

        Assert.That(contigs, Is.EqualTo(new[] { "AGA", "ATG", "ATG", "CAT", "GAT", "TGGA", "TGT" }));
    }

    [Test]
    public void FindMaximalNonBranchingPaths_EmptyGraph_NoPaths()
    {
        Assert.That(SequenceAssembler.FindMaximalNonBranchingPaths(new Dictionary<string, List<string>>()), Is.Empty);
    }

    [Test]
    public void FindMaximalNonBranchingPaths_NullGraph_Throws()
    {
        Assert.That(() => SequenceAssembler.FindMaximalNonBranchingPaths(null!),
            NUnit.Framework.Throws.ArgumentNullException);
    }

    #endregion
}
