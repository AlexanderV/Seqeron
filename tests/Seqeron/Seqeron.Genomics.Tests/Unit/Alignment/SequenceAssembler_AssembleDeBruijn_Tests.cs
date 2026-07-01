// ASSEMBLY-DBG-001 — De Bruijn graph assembly
// Evidence: docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-DBG-001.md
// Source: Langmead B., "De Bruijn Graph assembly" (JHU lecture notes), p.5-22;
//         Jones NC, Pevzner PA (2004), An Introduction to Bioinformatics Algorithms,
//         MIT Press, Theorems 8.1/8.2; Compeau, Pevzner & Tesler (2011),
//         Nat Biotechnol 29:987-991, DOI 10.1038/nbt.2023.

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_AssembleDeBruijn_Tests
{
    // Published exact-reconstruction inputs (unique Eulerian walk).
    // Source: Langmead DBG p.11 (AAABBBA), p.18 (a_long...), p.19/22 (to_every...).
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

    #region AssembleDeBruijn

    // M4 — k=3 single read "AAABBBA" reconstructs one contig "AAABBBA".
    // Source: Langmead DBG p.11 (Eulerian walk AA->AA->AB->BB->BB->BA).
    [Test]
    public void AssembleDeBruijn_Aaabbba_ReconstructsSingleContig()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { "AAABBBA" },
            new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1),
                "A connected, two-semi-balanced-node graph yields a single Eulerian-walk contig.");
            Assert.That(result.Contigs[0], Is.EqualTo("AAABBBA"),
                "INV-04: spelling the unique Eulerian walk reconstructs the source genome (Langmead DBG p.11).");
        });
    }

    // M5 — k=5 reconstructs "a_long_long_long_time" exactly. Source: Langmead DBG p.18.
    [Test]
    public void AssembleDeBruijn_ALongRepeat_ReconstructsInput()
    {
        const string genome = "a_long_long_long_time";

        var result = SequenceAssembler.AssembleDeBruijn(new[] { genome },
            new SequenceAssembler.AssemblyParameters(KmerSize: 5, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1), "The graph is connected and Eulerian -> one contig.");
            Assert.That(result.Contigs[0], Is.EqualTo(genome),
                "INV-04: the printed Eulerian walk for k=5 spells a_long_long_long_time (Langmead DBG p.18).");
        });
    }

    // M6 — k=4 reconstructs the to_every... sentence exactly. Source: Langmead DBG p.19/22.
    [Test]
    public void AssembleDeBruijn_ToEveryK4_ReconstructsInput()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new[] { ToEvery },
            new SequenceAssembler.AssemblyParameters(KmerSize: 4, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1), "Connected Eulerian graph -> one contig.");
            Assert.That(result.Contigs[0], Is.EqualTo(ToEvery),
                "INV-04: at k=4 the Eulerian walk reconstructs the input exactly (Langmead DBG p.22).");
        });
    }

    // M7 — A repeat-free DNA string with a unique Eulerian path reconstructs verbatim.
    // Source: construction + spelling rule (Langmead DBG p.6-7, p.18).
    [Test]
    public void AssembleDeBruijn_DnaUniqueWalk_ReconstructsInput()
    {
        const string genome = "ATGGCGTGCA"; // all 4-mers distinct, no 3-mer repeats

        var result = SequenceAssembler.AssembleDeBruijn(new[] { genome },
            new SequenceAssembler.AssemblyParameters(KmerSize: 4, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1), "Unique Eulerian path -> single contig.");
            Assert.That(result.Contigs[0], Is.EqualTo(genome),
                "INV-04: a repeat-free k-mer set reconstructs the source DNA exactly.");
        });
    }

    // S3 — INV-05: every input k-mer is a substring of the reconstructed contig.
    // Source: Langmead DBG p.7, p.18.
    [Test]
    public void AssembleDeBruijn_DnaUniqueWalk_EveryKmerIsSubstringOfContig()
    {
        const string genome = "ATGGCGTGCA";
        const int k = 4;
        var kmers = Enumerable.Range(0, genome.Length - k + 1)
            .Select(i => genome.Substring(i, k)).ToList();

        var result = SequenceAssembler.AssembleDeBruijn(new[] { genome },
            new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));
        string contig = result.Contigs.Single();

        Assert.That(kmers.All(km => contig.Contains(km)), Is.True,
            "INV-05: each input k-mer (an edge traversed by the walk) appears in the spelled contig.");
    }

    // S4 — INV-06: result statistics are consistent with the emitted contigs.
    [Test]
    public void AssembleDeBruijn_ALongRepeat_StatisticsAreConsistent()
    {
        const string genome = "a_long_long_long_time";

        var result = SequenceAssembler.AssembleDeBruijn(new[] { genome },
            new SequenceAssembler.AssemblyParameters(KmerSize: 5, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.TotalLength, Is.EqualTo(result.Contigs.Sum(c => c.Length)),
                "INV-06: TotalLength equals the sum of contig lengths.");
            Assert.That(result.LongestContig, Is.EqualTo(result.Contigs.Max(c => c.Length)),
                "INV-06: LongestContig equals the maximum contig length.");
            Assert.That(result.TotalReads, Is.EqualTo(1), "TotalReads equals the input read count.");
        });
    }

    // Disconnected graph — two reads that share no (k-1)-mer form two weakly-connected
    // components, each individually Eulerian, so assembly yields one contig per component,
    // each reconstructing its source read exactly. Source: Langmead DBG p.24-25
    // ("Connected components are individually Eulerian, overall graph is not"); algorithm
    // step "one Eulerian walk per weakly-connected component".
    // ATGGCGTGCA and GATTACAGGTC each have all-distinct 3-mers (unique walk) and share no
    // 3-mer node; expected contig set re-derived independently with a Hierholzer reference
    // (see report ASSEMBLY-DBG-001.md), not from this implementation's output.
    [Test]
    public void AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent()
    {
        var reads = new[] { "ATGGCGTGCA", "GATTACAGGTC" };

        var result = SequenceAssembler.AssembleDeBruijn(reads,
            new SequenceAssembler.AssemblyParameters(KmerSize: 4, MinContigLength: 1));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(2),
                "Two disjoint reads => two weakly-connected components => two contigs (Langmead DBG p.24).");
            Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal),
                Is.EqualTo(new[] { "ATGGCGTGCA", "GATTACAGGTC" }),
                "INV-04: each component's unique Eulerian walk reconstructs its source read exactly.");
        });
    }

    // MinContigLength filter — contigs shorter than MinContigLength are discarded (contract §3.1).
    // With the two disjoint reads above (lengths 10 and 11), MinContigLength=11 keeps only the
    // 11-mer contig GATTACAGGTC and drops ATGGCGTGCA.
    [Test]
    public void AssembleDeBruijn_MinContigLength_FiltersShortContigs()
    {
        var reads = new[] { "ATGGCGTGCA", "GATTACAGGTC" }; // lengths 10 and 11

        var result = SequenceAssembler.AssembleDeBruijn(reads,
            new SequenceAssembler.AssemblyParameters(KmerSize: 4, MinContigLength: 11));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1),
                "The length-10 contig is below MinContigLength=11 and is discarded (contract §3.1).");
            Assert.That(result.Contigs.Single(), Is.EqualTo("GATTACAGGTC"),
                "Only the length-11 contig survives the MinContigLength filter.");
        });
    }

    // M8 — Empty read set returns an empty AssemblyResult. ASSUMPTION-2.
    [Test]
    public void AssembleDeBruijn_EmptyReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new List<string>());

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(0), "No reads -> no contigs.");
            Assert.That(result.TotalReads, Is.EqualTo(0), "No reads -> TotalReads 0.");
            Assert.That(result.TotalLength, Is.EqualTo(0), "No reads -> TotalLength 0.");
        });
    }

    // M9 — Null read set is handled like empty (no exception). Contract (mirrors OLC).
    [Test]
    public void AssembleDeBruijn_NullReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(null!);

        Assert.That(result.Contigs.Count, Is.EqualTo(0),
            "Null input is treated as empty, returning no contigs.");
    }

    // Determinism — repeated assembly of the same reads yields identical contigs (order-independent).
    [Test]
    public void AssembleDeBruijn_RepeatedRuns_AreDeterministic()
    {
        var reads = new[] { "ATGGCGTGCA" };
        var p = new SequenceAssembler.AssemblyParameters(KmerSize: 4, MinContigLength: 1);

        var first = SequenceAssembler.AssembleDeBruijn(reads, p).Contigs;
        var second = SequenceAssembler.AssembleDeBruijn(reads, p).Contigs;

        Assert.That(second, Is.EqualTo(first),
            "Reconstruction is deterministic: the same input always yields the same contigs.");
    }

    #endregion
}
