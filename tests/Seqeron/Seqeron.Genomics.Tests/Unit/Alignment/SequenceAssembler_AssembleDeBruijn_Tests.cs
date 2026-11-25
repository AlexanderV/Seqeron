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

    // The oracle values of the BA3K/BA3M tests below are for the single-strand textbook graph; AssembleDeBruijn is
    // double-stranded by default (Velvet, A2-4), so they select includeReverseComplement: false explicitly. The
    // double-stranded default has its own region (Double-stranded default, Velvet read rule).
    private static SequenceAssembler.AssemblyResult SingleStrand(
        IReadOnlyList<string> reads, SequenceAssembler.AssemblyParameters parameters) =>
        SequenceAssembler.AssembleDeBruijn(reads, parameters, includeReverseComplement: false, clipTips: false);

    // The unitig oracle values below are for the uncleaned graph (Velvet's tip clipping, A2-2, would drop their
    // toy-sized contigs, which are shorter than 2k k-mers); cleaning has its own region (Velvet graph cleaning).
    private static SequenceAssembler.AssemblyResult Unclipped(
        IReadOnlyList<string> reads, SequenceAssembler.AssemblyParameters? parameters = null,
        bool includeReverseComplement = true) =>
        SequenceAssembler.AssembleDeBruijn(reads, parameters, includeReverseComplement, clipTips: false);

    // Contigs of an arbitrary-alphabet text (Langmead / Rosalind letters are not DNA, and AssembleDeBruijn reads every
    // non-ACGT symbol as A, as Velvet does): the BA3K/BA3M construction over the distinct upper-cased k-mers.
    private static List<string> TextUnitigs(string text, int k) => TextUnitigs(new[] { text }, k);

    private static List<string> TextUnitigs(IEnumerable<string> texts, int k)
    {
        var kmers = texts.Select(t => t.ToUpperInvariant())
            .SelectMany(upper => Enumerable.Range(0, upper.Length - k + 1).Select(i => upper.Substring(i, k)))
            .Distinct().OrderBy(x => x, StringComparer.Ordinal).ToList();
        var graph = SequenceAssembler.BuildDeBruijnGraph(kmers, k);
        return SequenceAssembler.FindMaximalNonBranchingPaths(graph).Select(Spell)
            .OrderBy(c => c, StringComparer.Ordinal).ToList();
    }

    private static string Spell(IReadOnlyList<string> path) =>
        path[0] + string.Concat(path.Skip(1).Select(n => n[^1]));

    // M4 — AAABBBA, k=3: AA and BB carry self-loops, so neither is 1-in-1-out and the graph
    // splits into the unitigs AAA, AABB, BBA, BBB (BA3M; Velvet merge rule).
    [Test]
    public void AssembleDeBruijn_Aaabbba_EmitsMaximalNonBranchingPaths()
    {
        Assert.Multiple(() =>
        {
            Assert.That(TextUnitigs("AAABBBA", 3), Is.EqualTo(new[] { "AAA", "AABB", "BBA", "BBB" }),
                "Contigs are the spelled maximal non-branching paths of the textbook text graph (BA3K/BA3M).");

            // The same topology in DNA (B -> C): the assembler itself, single strand.
            var result = SingleStrand(new[] { "AAACCCA" }, K(3));
            Assert.That(result.Contigs, Is.EqualTo(new[] { "AAA", "AACC", "CCA", "CCC" }),
                "Contigs are the spelled maximal non-branching paths, starting nodes in ordinal order (BA3K/BA3M).");
        });
    }

    // M5 — a_long_long_long_time, k=5: the repeated 4-mer nodes "_lon"/"ong_" branch, so the
    // repeat is not traversed three times; the contigs are the four unitigs (upper-cased).
    [Test]
    public void AssembleDeBruijn_ALongRepeat_SplitsAtRepeatBranches()
    {
        Assert.That(TextUnitigs("a_long_long_long_time", 5),
            Is.EqualTo(new[] { "A_LON", "ONG_LON", "ONG_TIME", "_LONG_" }),
            "A repeat longer than k-2 creates branch nodes; contigs stop there (BA3K: contigs must be present in any assembly).");
    }

    // M6 — to_every…, k=4: the "turn" repeat yields six unitigs.
    [Test]
    public void AssembleDeBruijn_ToEveryK4_SplitsAtTurnRepeat()
    {
        Assert.That(TextUnitigs(ToEvery, 4),
            Is.EqualTo(new[] { "N_TH", "N_TU", "TO_EVERY_TH", "_THERE_IS_A_SEASON", "_THING_TU", "_TURN_T" }));
    }

    // M7 — A repeat-free DNA string is one non-branching path, reconstructed verbatim.
    [Test]
    public void AssembleDeBruijn_DnaUniqueWalk_ReconstructsInput()
    {
        var result = SingleStrand(new[] { "ATGGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // F-DBG-1 — Overlapping reads share k-mers; coverage must not add edges (Velvet: "high
    // redundancy is naturally handled by the graph without affecting the number of nodes").
    // Before the fix the shared 4-mers GGCG/GCGT became parallel edges and the Eulerian walk
    // spelled a chimera.
    [Test]
    public void AssembleDeBruijn_OverlappingReads_AssembleIntoGenome()
    {
        var result = SingleStrand(new[] { "ATGGCGT", "GGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // F-DBG-1 — Two copies of the same read (coverage 2) assemble to that read; before the fix the
    // doubled edges produced "ACGGTT", which contains 3-mers (CGG, GGT, GTT) absent from the input.
    [Test]
    public void AssembleDeBruijn_DuplicateReads_DoNotDoubleTheGenome()
    {
        var result = SingleStrand(new[] { "ACGT", "ACGT" }, K(3));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGT" }));
    }

    // F-DBG-1 — MCP example reads AAABBB / AABBBC at k=3: former output "AAABBBBBBC" (BBB walked
    // three times); the distinct-k-mer unitigs are AAA, AABB, BBB, BBC.
    [Test]
    public void AssembleDeBruijn_SharedKmersAcrossReads_CountedOnce()
    {
        Assert.Multiple(() =>
        {
            Assert.That(TextUnitigs(new[] { "AAABBB", "AABBBC" }, 3), Is.EqualTo(new[] { "AAA", "AABB", "BBB", "BBC" }));
            var result = SingleStrand(new[] { "AAACCC", "AACCCG" }, K(3));
            Assert.That(result.Contigs, Is.EqualTo(new[] { "AAA", "AACC", "CCC", "CCG" }));
        });
    }

    // Isolated cycle — every node is 1-in-1-out; BA3M emits the cycle (from its ordinally
    // smallest node ACG back to ACG): ACGTTACG.
    [Test]
    public void AssembleDeBruijn_IsolatedCycle_IsEmitted()
    {
        var result = SingleStrand(new[] { "ACGTTACG" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTACG" }),
            "BA3M: 'MaximalNonBranchingPaths finds all isolated cycles in the graph'.");
    }

    // Cycle attached to a branch: ACGTACGTACGT + ACGTACGTTTTT at k=5 -> ACGTACGT (the loop from the
    // branching node ACGT back to itself), ACGTTTT and the TTTT self-loop TTTTT.
    [Test]
    public void AssembleDeBruijn_LoopThroughBranchNode_EmittedAsPath()
    {
        var result = SingleStrand(new[] { "ACGTACGTACGT", "ACGTACGTTTTT" }, K(5));

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

        var contigs = TextUnitigs(ToEvery, k);
        var emitted = contigs.SelectMany(c => Enumerable.Range(0, c.Length - k + 1).Select(i => c.Substring(i, k)))
            .OrderBy(s => s, StringComparer.Ordinal).ToList();

        Assert.That(emitted, Is.EqualTo(input), "each distinct k-mer is spelled exactly once over all contigs");
    }

    // Case — k-mers are counted case-insensitively by the canonical counter (CountKmersSpan);
    // contigs are upper-case.
    [Test]
    public void AssembleDeBruijn_MixedCaseReads_JoinCaseInsensitively()
    {
        var result = SingleStrand(new[] { "atggcGT", "GGCGTGCA" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA" }));
    }

    // Reads shorter than k contribute nothing but are still counted in TotalReads.
    [Test]
    public void AssembleDeBruijn_ReadShorterThanK_IsSkipped()
    {
        var result = SingleStrand(new[] { "AC", "ATGGCGTGCA" }, K(4));

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

        var forward = SingleStrand(reads, K(5));
        var both = Unclipped(reads, K(5), includeReverseComplement: true);
        var byDefault = Unclipped(reads, K(5));

        Assert.Multiple(() =>
        {
            Assert.That(forward.Contigs, Is.EqualTo(new[] { "GATTACAGGT", "TAAGGACCTG" }));
            Assert.That(both.Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
            Assert.That(byDefault.Contigs, Is.EqualTo(both.Contigs), "double-stranded is the default (Velvet)");
        });
    }

    // Double-stranded option on a disconnected pair: one contig per strand pair, canonical orientation.
    [Test]
    public void AssembleDeBruijn_ReverseComplementOption_ReportsEachStrandPairOnce()
    {
        var both = Unclipped(new[] { "GATTACAGGTCCTTA" }, K(5), includeReverseComplement: true);
        var fromReverse = Unclipped(new[] { "TAAGGACCTGTAATC" }, K(5), includeReverseComplement: true);

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

        string contig = SingleStrand(new[] { genome }, K(k)).Contigs.Single();

        Assert.That(kmers.All(km => contig.Contains(km)), Is.True);
    }

    // S4 — INV-06: result statistics are consistent with the emitted contigs.
    [Test]
    public void AssembleDeBruijn_ALongRepeat_StatisticsAreConsistent()
    {
        // DNA repeat example (oracle: contigs ACGTACGT, ACGTTTT, TTTTT at k = 5; see the loop test above).
        var result = SingleStrand(new[] { "ACGTACGTACGT", "ACGTACGTTTTT" }, K(5));

        Assert.Multiple(() =>
        {
            Assert.That(result.TotalLength, Is.EqualTo(8 + 7 + 5), "ACGTACGT + ACGTTTT + TTTTT");
            Assert.That(result.LongestContig, Is.EqualTo(8));
            Assert.That(result.TotalReads, Is.EqualTo(2));
            Assert.That(result.AssembledReads, Is.EqualTo(2));
        });
    }

    // Two reads sharing no (k-1)-mer -> two disconnected simple paths -> two contigs.
    [Test]
    public void AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent()
    {
        var result = SingleStrand(new[] { "ATGGCGTGCA", "GATTACAGGTC" }, K(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA", "GATTACAGGTC" }));
    }

    // MinContigLength filter — the length-10 contig is dropped at MinContigLength = 11.
    [Test]
    public void AssembleDeBruijn_MinContigLength_FiltersShortContigs()
    {
        var result = SingleStrand(new[] { "ATGGCGTGCA", "GATTACAGGTC" }, K(4, 11));

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
            Assert.That(() => SingleStrand(new[] { "ACGT" }, K(k)),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => SingleStrand(Array.Empty<string>(), K(k)),
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

        var first = SingleStrand(reads, K(4)).Contigs;
        var reversed = SingleStrand(reads.Reverse().ToList(), K(4)).Contigs;

        Assert.That(reversed, Is.EqualTo(first));
    }

    #endregion

    #region Double-stranded default (Velvet) and the Velvet read rule (non-ACGT -> A)

    // Expected values: tools/review/oracles/B10/wp4_oracle.py (dbg_expected): the folded distinct k-mer set closed
    // under reverse complement, unitigs by two independent constructions (Rosalind BA3M pseudocode and Velvet's
    // node-centric chain merge, identical on every case), each {U, revcomp(U)} pair reported once as the
    // ordinally smaller string. 7 seeds x 1500 random cases (even/odd k, both strands, N/IUPAC/U/lower-case noise)
    // matched the C# output exactly (Evidence R1-WP4).

    // F-A2-4 — the default of both overloads is double-stranded: GATTACAGGT + the reverse complement of the
    // genome's 3' half join without any option.
    [Test]
    public void AssembleDeBruijn_Default_IsDoubleStranded()
    {
        var reads = new[] { "GATTACAGGT", "TAAGGACCTG" };

        Assert.Multiple(() =>
        {
            Assert.That(Unclipped(reads, K(5)).Contigs, Is.EqualTo(new[] { "GATTACAGGTCCTTA" }));
            Assert.That(Unclipped(reads, K(5), includeReverseComplement: false).Contigs,
                Is.EqualTo(new[] { "GATTACAGGT", "TAAGGACCTG" }));
        });
    }

    // Strand symmetry of the double-stranded graph (Velvet twin nodes): reverse-complementing every read cannot
    // change the contig set. ATGGCGTGCA contains the palindromic 4-mer TGCA, so at k = 5 its strand continues
    // into its reverse complement: the single contig ATGGCGTGCACGCCAT equals its own reverse complement.
    [Test]
    public void AssembleDeBruijn_Default_ReverseComplementedReadsGiveTheSameContigs()
    {
        var forward = Unclipped(new[] { "ATGGCGTGCA" }, K(5)).Contigs;
        var reverse = Unclipped(new[] { "TGCACGCCAT" }, K(5)).Contigs;
        var overlapping = Unclipped(new[] { "ATGGCGT", "GGCGTGCA" }, K(5)).Contigs;

        Assert.Multiple(() =>
        {
            Assert.That(forward, Is.EqualTo(new[] { "ATGGCGTGCACGCCAT" }));
            Assert.That(reverse, Is.EqualTo(forward));
            Assert.That(overlapping, Is.EqualTo(forward));
        });
    }

    // Palindrome-free double-stranded example with a repeat (oracle): the two loops of
    // ACGTACGTACGT/ACGTACGTTTTT join with their reverse complements into AAAAA, AAAACGT, ACGTACGT.
    [Test]
    public void AssembleDeBruijn_Default_RepeatReadsJoinWithReverseComplementUnitigs()
    {
        var result = Unclipped(new[] { "ACGTACGTACGT", "ACGTACGTTTTT" }, K(5));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "AAAAA", "AAAACGT", "ACGTACGT" }));
    }

    // Velvet read loader: velvetifySequence maps every symbol other than A/C/G/T (either case) to N and
    // writeNucleotide stores N as Adenine, so a k-mer with N, an IUPAC code or U is the k-mer with an A there.
    [Test]
    public void AssembleDeBruijn_NonAcgtSymbols_AreReadAsA()
    {
        var withN = Unclipped(new[] { "ACGTNGGTCA" }, K(4), includeReverseComplement: false);
        var withIupacAndU = Unclipped(new[] { "ACGTRYGGTCAU" }, K(5), includeReverseComplement: false);
        var withOtherSymbols = Unclipped(new[] { "AAGCTTCCGANNGGATC" }, K(5), includeReverseComplement: false);
        var bothStrands = Unclipped(new[] { "ACGTNGGTCA" }, K(5));

        Assert.Multiple(() =>
        {
            Assert.That(withN.Contigs, Is.EqualTo(new[] { "ACGTAGGTCA" }), "N -> A");
            Assert.That(withIupacAndU.Contigs, Is.EqualTo(new[] { "ACGTAAGGTCAA" }), "R, Y and U -> A");
            Assert.That(withOtherSymbols.Contigs, Is.EqualTo(new[] { "AAGCTTCCGAAAGGATC" }));
            Assert.That(bothStrands.Contigs, Is.EqualTo(new[] { "TGACCTACGTAGGTCA" }),
                "the folded read is ACGTAGGTCA; its reverse-complement k-mers join it into a 16-mer");
        });
    }

    // A k-mer containing N counts as the k-mer with A, so reads that differ only by N/n/lower-case/IUPAC at the
    // same position are the same read for the graph (read-order and case independent).
    [Test]
    public void AssembleDeBruijn_NonAcgtAndLowerCase_ShareKmersWithTheAFoldedRead()
    {
        var mixed = Unclipped(new[] { "acgtnGGTCA", "ACGTAGGTCA" }, K(5), includeReverseComplement: false);
        var plain = Unclipped(new[] { "ACGTAGGTCA" }, K(5), includeReverseComplement: false);

        Assert.Multiple(() =>
        {
            Assert.That(mixed.Contigs, Is.EqualTo(new[] { "ACGTAGGTCA" }));
            Assert.That(mixed.Contigs, Is.EqualTo(plain.Contigs));
            Assert.That(mixed.AssembledReads, Is.EqualTo(2));
        });
    }

    // The non-ACGT rule also applies to non-ASCII letters (U+017F is not an 'S' that folds: it is not A/C/G/T).
    [Test]
    public void AssembleDeBruijn_NonAsciiSymbol_IsReadAsA()
    {
        var result = Unclipped(new[] { "ACGT\u017fGGTCA" }, K(4), includeReverseComplement: false);

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTAGGTCA" }));
    }

    // AssembledReads (Velvet usedReads): a read counts when one of its k-mers (folded; on either strand by
    // default) lies in a reported contig. Oracle: GATTACAGGTC + ATGGCGTGCAATC at k = 5 -> contigs ATGGCGTGCA,
    // GACCTGTAATC, GATTGCA (both strands); both reads used.
    [Test]
    public void AssembleDeBruijn_Default_AssembledReadsCountsEitherStrand()
    {
        var result = Unclipped(new[] { "GATTACAGGTC", "ATGGCGTGCAATC" }, K(5));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs, Is.EqualTo(new[] { "ATGGCGTGCA", "GACCTGTAATC", "GATTGCA" }));
            Assert.That(result.AssembledReads, Is.EqualTo(2));
            Assert.That(result.TotalReads, Is.EqualTo(2));
        });
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

    #region Velvet graph cleaning (A2-1 coverage cut-off, A2-2 tip clipping)

    // Sources: Velvet 1.2.10 (dzerbino/velvet) src/concatenatedPreGraph.c isEligibleTip/clipTips_pg, src/correctedGraph.c
    // getTipLength/clipTipsHard, src/graphStats.c removeLowCoverageNodesAndDenounceDubiousReads/estimated_cov, src/run2.c
    // (pipeline order); Zerbino & Birney (2008). Oracle: tools/review/oracles/B10/wp5_oracle.py (literal transcription,
    // node-centric unitigs, bidirected twins); the values below are the oracle's, the C# side agrees on 1200 random cases
    // (3 seeds x 400, see ASSEMBLY-DBG-001-Evidence.md "R1-WP5").

    // 50-nt repeat-free genome (no shared 6-mer between it and its reverse complement) read at k = 7.
    private const string TipGenome = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCA";

    // Two reads cover the genome at 2x; the third carries a substitution at its last base (TipGenome[..29] + 'C'),
    // i.e. one erroneous k-mer hanging off the genome (a tip of length 1 < 2k).
    private static readonly string[] TipReads =
    {
        "GCTAAAGACAATTACATAACATACACGTCAGCACGA", "GCTAAAGACAATTACATAACATACACGTCAGCACGA",
        "CATAACATACACGTCAGCACGAAACTTGTTGGCCCA", "CATAACATACACGTCAGCACGAAACTTGTTGGCCCA",
        "GCTAAAGACAATTACATAACATACACGTCC",
    };

    [Test]
    public void AssembleDeBruijn_TipClipping_IsOnByDefault_AndRemovesTheErroneousTip()
    {
        var clipped = SequenceAssembler.AssembleDeBruijn(TipReads, K(7));
        var unclipped = SequenceAssembler.AssembleDeBruijn(TipReads, K(7), includeReverseComplement: true, clipTips: false);

        Assert.Multiple(() =>
        {
            Assert.That(clipped.Contigs, Is.EqualTo(new[] { TipGenome }), "Velvet clipTips: the erroneous k-mer is clipped, the genome is one contig.");
            Assert.That(unclipped.Contigs, Is.EqualTo(new[]
            {
                "CACGTCAGCACGAAACTTGTTGGCCCA", "CACGTCC", "GACGTGTATGTTATGTAATTGTCTTTAGC",
            }), "Without cleaning the tip forks the unitig graph (BA3M contigs).");
            Assert.That(clipped.AssembledReads, Is.EqualTo(5), "the erroneous read still has genome k-mers in the contig (Velvet usedReads)");
        });
    }

    [Test]
    public void AssembleDeBruijn_TipClipping_SingleStrandModel()
    {
        // Same genome, forward strand only: the oracle (rc = false) gives the genome as one contig.
        var clipped = SequenceAssembler.AssembleDeBruijn(TipReads, K(7), includeReverseComplement: false);
        Assert.That(clipped.Contigs, Is.EqualTo(new[] { TipGenome }));
    }

    [Test]
    public void AssembleDeBruijn_MinTipLength_IsInKmers_AndALengthEqualToTheCutoffIsNotATip()
    {
        // Velvet: a chain is a tip when its length (in k-mers) is < 2k; the erroneous tail is 1 k-mer.
        var atLength = SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, clipTips: true, minTipLength: 1);
        var above = SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, clipTips: true, minTipLength: 2);
        var zeroIsTwoK = SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, clipTips: true, minTipLength: 0);

        Assert.Multiple(() =>
        {
            Assert.That(atLength.Contigs.Count, Is.EqualTo(3), "1 < 1 is false: the tip stays");
            Assert.That(above.Contigs, Is.EqualTo(new[] { TipGenome }));
            Assert.That(zeroIsTwoK.Contigs, Is.EqualTo(new[] { TipGenome }));
        });
    }

    [Test]
    public void AssembleDeBruijn_IsolatedSnippetShorterThanTwoK_IsRemoved_UnlessItsLengthReachesTheCutoff()
    {
        // Velvet isEligibleTip: "If isolated snippet" (the dead-end walk reaches no junction) -> removed when < 2k k-mers.
        // The 12-nt read TTTCCTCATGCA shares no 6-mer with the genome or its reverse complement: 6 k-mers at k = 7.
        string[] reads = { TipGenome, TipGenome, "TTTCCTCATGCA" };
        var byDefault = SequenceAssembler.AssembleDeBruijn(reads, K(7));
        var cutoff7 = SequenceAssembler.AssembleDeBruijn(reads, K(7), true, minTipLength: 7);
        var cutoff6 = SequenceAssembler.AssembleDeBruijn(reads, K(7), true, minTipLength: 6);
        var unclipped = SequenceAssembler.AssembleDeBruijn(reads, K(7), true, clipTips: false);

        Assert.Multiple(() =>
        {
            Assert.That(byDefault.Contigs, Is.EqualTo(new[] { TipGenome }));
            Assert.That(cutoff7.Contigs, Is.EqualTo(new[] { TipGenome }), "6 < 7");
            Assert.That(cutoff6.Contigs, Is.EqualTo(new[] { TipGenome, "TGCATGAGGAAA" }), "6 < 6 is false");
            Assert.That(unclipped.Contigs, Is.EqualTo(new[] { TipGenome, "TGCATGAGGAAA" }));
            Assert.That(byDefault.AssembledReads, Is.EqualTo(2));
        });
    }

    [Test]
    public void AssembleDeBruijn_TipClipping_CleanReadsAreUnchanged()
    {
        // Error-free reads of a repeat-free genome have no tips: clipping changes nothing.
        string[] reads = { TipGenome[..36], TipGenome[14..] };
        Assert.That(SequenceAssembler.AssembleDeBruijn(reads, K(7)).Contigs,
            Is.EqualTo(SequenceAssembler.AssembleDeBruijn(reads, K(7), true, clipTips: false).Contigs));
    }

    // Velvet: getTotalCoverage / getNodeLength is an integer division of two int64 values compared with the double cut-off.
    // Single strand, k = 3: unitig ACGT = k-mers ACG (4 reads) + CGT (3): coverage 7, length 2 -> quotient 3.
    [TestCase(3.5, new[] { "GCATTAGC", "GGC" })]
    [TestCase(3.0, new[] { "ACGT", "GCATTAGC", "GGC" })]
    [TestCase(0.0, new[] { "ACGT", "GCATTAGC", "GGC" })]
    public void AssembleDeBruijn_CoverageCutoff_UsesTheIntegerQuotient(double cutoff, string[] expected)
    {
        string[] reads = { "ACGT", "ACGT", "ACGT", "ACG", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC", "GGCATTAGC" };
        var result = SequenceAssembler.AssembleDeBruijn(reads, K(3), includeReverseComplement: false, clipTips: false, coverageCutoff: cutoff);
        Assert.That(result.Contigs, Is.EqualTo(expected));
    }

    [Test]
    public void AssembleDeBruijn_CoverageCutoff_CountsBothStrandsOfADoubleStrandedRead()
    {
        // Velvet node coverage = k-mer occurrences of either strand (a node and its twin share it): the 3 occurrences
        // (2 forward + 1 reverse complement) of every k-mer give quotient 3.
        string[] reads = { TipGenome, TipGenome, "TGGGCCAACAAGTTTCGTGCTGACGTGTATGTTATGTAATTGTCTTTAGC" };
        string rc = "TGGGCCAACAAGTTTCGTGCTGACGTGTATGTTATGTAATTGTCTTTAGC";
        Assert.That(rc, Is.EqualTo(ReverseComplementOf(TipGenome)));
        var at3 = SequenceAssembler.AssembleDeBruijn(reads, K(7), true, clipTips: false, coverageCutoff: 3);
        var at35 = SequenceAssembler.AssembleDeBruijn(reads, K(7), true, clipTips: false, coverageCutoff: 3.5);
        Assert.Multiple(() =>
        {
            Assert.That(at3.Contigs, Is.EqualTo(new[] { TipGenome }));
            Assert.That(at35.Contigs, Is.Empty);
        });
    }

    private static string ReverseComplementOf(string s) =>
        new(s.Reverse().Select(c => c switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', _ => c }).ToArray());

    [Test]
    public void AssembleDeBruijn_AutoCoverageCutoff_IsHalfTheLengthWeightedMedianCoverage()
    {
        // Velvet estimated_cov: unitigs GAATTCAGGTCCAC (11 k-mers, coverage 3 each), ACATG and GCTGG (2 k-mers, coverage 1):
        // total length 15, half = 7: sorted by coverage B(2) C(4) A(15 >= 7) -> median 3.0 -> cut-off 1.5. An unweighted
        // median over the three nodes is 1 (cut-off 0.5) and would keep B and C.
        string[] reads = { "GAATTCAGGTCCAC", "GAATTCAGGTCCAC", "GAATTCAGGTCCAC", "ACATG", "GCTGG" };
        var auto = SequenceAssembler.AssembleDeBruijn(reads, K(4), false, clipTips: false, autoCoverageCutoff: true);
        var explicitCutoff = SequenceAssembler.AssembleDeBruijn(reads, K(4), false, clipTips: false, coverageCutoff: 1.5);
        var off = SequenceAssembler.AssembleDeBruijn(reads, K(4), false, clipTips: false);
        var explicitWins = SequenceAssembler.AssembleDeBruijn(reads, K(4), false, clipTips: false, coverageCutoff: 0.5, autoCoverageCutoff: true);

        Assert.Multiple(() =>
        {
            Assert.That(auto.Contigs, Is.EqualTo(new[] { "GAATTCAGGTCCAC" }));
            Assert.That(explicitCutoff.Contigs, Is.EqualTo(auto.Contigs));
            Assert.That(off.Contigs, Is.EqualTo(new[] { "ACATG", "GAATTCAGGTCCAC", "GCTGG" }));
            Assert.That(explicitWins.Contigs, Is.EqualTo(off.Contigs), "an explicit cut-off > 0 overrides auto (run2.c: coverageCutoff < 0 && estimateCutoff)");
            Assert.That(auto.AssembledReads, Is.EqualTo(3));
        });
    }

    [Test]
    public void AssembleDeBruijn_CoverageCutoff_RemovesTheLowCoverageBranchThenTheJoinedUnitigsMerge()
    {
        // The error read forks the genome with a 7-k-mer dead-end branch (ACCTGAG after position 24). With clipping off only
        // the coverage rule can remove it: genome coverage 2 per k-mer (read x2), branch coverage 1; cut-off 1.5 removes the
        // branch unitig and the surviving halves concatenate into the genome (Velvet re-concatenates, second pass).
        string errorRead = TipGenome[..24] + "ACCTGAG";
        string[] reads = { TipGenome, TipGenome, errorRead };
        var withTail = SequenceAssembler.AssembleDeBruijn(reads, K(7), false, clipTips: false);
        var cut = SequenceAssembler.AssembleDeBruijn(reads, K(7), false, clipTips: false, coverageCutoff: 1.5);

        Assert.Multiple(() =>
        {
            Assert.That(withTail.Contigs, Is.EqualTo(new[]
            {
                "ATACACCTGAG", "ATACACGTCAGCACGAAACTTGTTGGCCCA", "GCTAAAGACAATTACATAACATACAC",
            }));
            Assert.That(cut.Contigs, Is.EqualTo(new[] { TipGenome }));
        });
    }

    [Test]
    public void AssembleDeBruijn_InvalidCleaningParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, minTipLength: -1),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, coverageCutoff: -1),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.AssembleDeBruijn(TipReads, K(7), true, coverageCutoff: double.NaN),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void AssembleDeBruijn_Cleaning_SelfComplementaryUnitig_IsTwoTwinNodesJoinedByAHairpinArc()
    {
        // GCTTGACATCGATGTCAAGC is its own reverse complement (k = 7): its unitig passes the self-complementary 6-mer vertex
        // ATCGAT; Velvet never concatenates across a hairpin arc (concatenatePreGraph_pg skips
        // isLoop_pg), so it is two twin nodes of 7 k-mers each; as an isolated snippet the walk sums 7 + 7 = 14 k-mers
        // (Velvet isEligibleTip), so it is a tip below a cut-off of 15 and not at 14 (and 14 is the default 2k).
        string[] reads = { "TTGACATCGATGTCAAGC", "TTGACATCGATGTCAAGC" };
        string[] palindrome = { "GCTTGACATCGATGTCAAGC" };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.AssembleDeBruijn(reads, K(7)).Contigs, Is.EqualTo(palindrome), "wp5_oracle.py");
            Assert.That(SequenceAssembler.AssembleDeBruijn(reads, K(7), true, clipTips: false).Contigs, Is.EqualTo(palindrome));
            Assert.That(SequenceAssembler.AssembleDeBruijn(reads, K(7), true, minTipLength: 14).Contigs, Is.EqualTo(palindrome));
            Assert.That(SequenceAssembler.AssembleDeBruijn(reads, K(7), true, minTipLength: 15).Contigs, Is.Empty);
        });
    }

    // Pinned random cases of wp5_oracle.py (seeds 1-3): reads, k, double strand, min contig length, clipTips, cut-off, auto,
    // expected contigs (as reported: strand pairs once, ordinally smaller), expected assembled reads.
    private static readonly object[][] OracleCases =
    {
        new object[] { new[] { "CTACAGCTCGCTAC", "AGCTCTACAGCTCG", "GCTCTACAGCTCGC", "TCGCTACTGCAGAC", "GCTCGCTACTGCAG", "CTCGCTACTGCAGA", "GTAGAGCTTGGGAT", "CAACGATCCCAAGC", "CTACAGCTCGCTAC" }, 5, true, 1, true, 4.5d, true, new[] { "GTAGCGAGCTCGCTAC" }, 8 },
        new object[] { new[] { "ACTTTACA", "CGGTCTTTGGG", "TAGAG", "AAGTAC", "GAATCC", "AGGACTATTCCATCG", "CTGGGACG" }, 3, false, 3, true, 2d, true, new[] { "GGACTT", "GGG", "TTT" }, 4 },
        new object[] { new[] { "GTGT", "TAAATCTTCTATG", "GTCGACTT", "GCAGAATCTCCTAAG", "GAGTATCG", "AGCACACCGTCCAG", "GTTTCCA" }, 3, true, 3, true, 3d, false, new[] { "GAAGA", "GATC", "TCGA" }, 5 },
        new object[] { new[] { "AATTCTTGTCAGACGCC", "AGTTCTTGTCAGACGCC", "TTATCGGCGTCTGACAA", "TCAGACGCCGATAAGCC", "TGTCAGACGCCGATCAG", "CTGGCTTATCGGCGTCT", "TCTGGCTTATCAGCGTC", "CAGACGCCGATAAGCCA", "TTGTCAGACGCCGATAA", "TGTCAGACGCGGATAAG" }, 9, false, 1, true, 1.5d, false, new[] { "AGTTCTTGTCAGACGCCGATAAGCCA" }, 7 },
        new object[] { new[] { "TTTCGCTAGCCCCCACT", "CCTCACTGCAGACCATC", "TGGTCTGCAGTGGGGGC", "TTCGCTAGCCCCCACTG", "CAGTGGGGGCTAGCGAA", "CCCCCGCTGCAGACCAT", "AGCCCCCACTGCAGACC", "GCCCCCACTGCAGACCA" }, 5, true, 1, true, 0d, true, new[] { "CCCCACTGCAGTGGGG", "CCCCC", "GGGGCTAGCCCC" }, 8 },
        new object[] { new[] { "CGGAATTGGCCA", "AACCGACG", "TGTTAGTTAGAGT", "ATGTTTCCGCT", "GACAACATG", "TCGCAGAAAATCGG", "GACGGATGTGCGAGTA", "CCATGGAAGTTTTAG" }, 7, true, 1, true, 0d, false, new[] { "CTAAAACTTCCATGGAAGTTTTAG" }, 1 },
        new object[] { new[] { "GGTGCGGCTCGCCCGG", "GCCTCTTACT", "TCCGGGCCATCGGGT", "GTCTG", "CGCGTTGG" }, 5, false, 5, true, 0d, false, new[] { "CGGGCCATCGGG", "GGTGCGGCTCGCCCGGG" }, 2 },
        new object[] { new[] { "CACATTCAGGTAAAGC", "TGATCCTTGCACAAAG", "AGCTGTACCTGAATGT", "TTTGTGCACGGATCAG", "CTGATCCTTGCACAAA", "TTCAGGTAAAGCTTTG", "AAGCTTTGTGCAAGGA", "TTCAGGTAAAGCTTTG" }, 7, true, 7, false, 2.5d, true, new[] { "AAAGCTTT" }, 3 },
        new object[] { new[] { "AGCACGTCCCAATCC", "ACGTGCTGGCAACGT", "TACGTTGCCAGCACG", "TGCCAGCACGTCCCA", "CGTCCCAATCCGTGG", "AATTACGTTGCCAGC", "ATTACGTTGCCAGCA", "CCAATTACGTTGCCA", "CGTGCTGGCAACGTA", "GATTACGTTGCCAGC" }, 7, false, 1, false, 1.5d, true, new[] { "ATTACGTTGCCAGCACGTCCCAATCCGTGG" }, 8 },
        new object[] { new[] { "GTTAAT", "ATAAGC", "TTATCCGGGATGGAGAG", "CCAAAAG", "GTCGGACCGT", "ACACTAAAGCTGGCAA", "TGGAACTCGTAAGGAG", "TCTAATTCTGCCACT" }, 5, true, 5, false, 2d, false, new[] { "TGCCA" }, 2 },
        new object[] { new[] { "CATGTGAAAATGATCACTG", "AATGATCACTGTACCTATT", "ATCGAATAGGTACAGTGAT", "CATGTGAAAATGATCACTG", "AATAGGTACAGTGATCATT", "CATGTGAAAATGATCACTG", "AAAATGATCACTGTCCCTA", "ATCACTGTACCTATTCGAT" }, 5, false, 10, false, 2d, false, new[] { "GTGAAAATGAT", "TGATCACTGTAC" }, 7 },
        new object[] { new[] { "CGAACATCTGGTAACGAT", "CCACGGGTCCCAGCCATC", "TAATATCTGATCCTCGA", "AAAGGGATATAAAGA" }, 7, true, 14, false, 0d, true, new[] { "AAAGGGATATAAAGA", "ATCGTTACCAGATGTTCG", "CCACGGGTCCCAGCCATC", "TAATATCTGATCCTCGA" }, 4 },
        new object[] { new[] { "TTGTCTAATGATGA" }, 7, false, 7, false, 0d, true, new[] { "TTGTCTAATGATGA" }, 1 },
        new object[] { new[] { "GGTAAAGGGTTCTCCATTGGG", "GTAAAGGGTTCTCCATTGGGT", "GTAACCTTACCCAATGGAGAA", "CCAATGGAGAACCCTTTACCG", "GAACCCCTTACCGACACGTCC", "TAAAGGGTTCTCCATTGGGTA", "TAACGTTACCCAATGGAGAAC" }, 9, true, 1, false, 0d, false, new[] { "CGGTAAAGGGTTCTCCATTGGGTAA", "GAACCCCTTACCGACACGTCC", "GTAACCTTACCCAA", "TTGGGTAACGTTACCCAA" }, 7 },
        new object[] { new[] { "ATTTTGAGGTCACGCT", "GGGGTGACGGGCACGC", "ATCTCCGC", "GTTAGGCAGCGGTGCTC", "TGGAGATGGTGCCTGAGTC" }, 7, false, 14, false, 0d, false, new[] { "ATTTTGAGGTCACGCT", "GGGGTGACGGGCACGC", "GTTAGGCAGCGGTGCTC", "TGGAGATGGTGCCTGAGTC" }, 4 },
    };

    [TestCaseSource(nameof(OracleCases))]
    public void AssembleDeBruijn_Cleaning_MatchesVelvetOracle(
        string[] reads, int k, bool doubleStrand, int minContig, bool clipTips, double cutoff, bool auto,
        string[] expected, int expectedAssembled)
    {
        var result = SequenceAssembler.AssembleDeBruijn(reads, K(k, minContig), doubleStrand, clipTips, 0, cutoff, auto);

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal), Is.EqualTo(expected.OrderBy(c => c, StringComparer.Ordinal)));
            Assert.That(result.AssembledReads, Is.EqualTo(expectedAssembled));
        });
    }

    #endregion
}
