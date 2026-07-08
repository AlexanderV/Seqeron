using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — De Bruijn Graph (DBG) assembly
/// (ASSEMBLY-DBG-001): the k-mer de Bruijn graph builder
/// <see cref="SequenceAssembler.BuildDeBruijnGraph(IReadOnlyList{string}, int)"/> and
/// the Eulerian-walk assembler
/// <see cref="SequenceAssembler.AssembleDeBruijn(IReadOnlyList{string}, SequenceAssembler.AssemblyParameters?)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// the code NEVER fails in an undisciplined way: no hang / infinite loop (the
/// Hierholzer Eulerian traversal must ALWAYS terminate, even on the cyclic multigraph
/// produced by repeated / all-identical reads), no negative-length substring or
/// IndexOutOfRange (notably when k &gt; read length, where the chop bound
/// i ∈ [0, |r|−k] is empty), no state corruption, and no nondeterministic contig
/// output. Every input must resolve to EITHER a well-formed, theory-correct result OR
/// a *documented, intentional* validation exception (ArgumentOutOfRangeException for
/// k &lt; 2 — contract §3.3, §6.1; null/empty reads → empty result, §6.1). A raw
/// runtime exception, a hang, a contig containing characters not present in the input,
/// a contig shorter than k−1, or an order-dependent contig set is a bug, not a passing
/// test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-DBG-001 — De Bruijn Graph Assembly
/// Checklist: docs/checklists/03_FUZZING.md, row 143.
/// Algorithm doc: docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row:
///          – k &gt; READ LENGTH: every read shorter than k → the chop bound
///            i ∈ [0, |r|−k] is EMPTY (INV-03), so NO k-mers, an empty graph, and an
///            empty contig list; NO negative-length Substring / IndexOutOfRange,
///            NO hang (§6.1 "Read length &lt; k", INV-03).
///          – SINGLE READ (length ≥ k, no repeated (k-1)-mer): the graph is a simple
///            path → exactly ONE contig that reconstructs the read verbatim
///            (INV-04, §4.1 spelling), no crash.
///          – ALL-IDENTICAL READS: many copies of the same read collapse onto the SAME
///            (k-1)-mer nodes (the k-mer space is shared) → a DETERMINISTIC contig set,
///            independent of read count and order, with NO infinite loop / no quadratic
///            hang on the multigraph (Hierholzer is O(|E|), §4.3); guarded with
///            [CancelAfter].
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (De_Bruijn_Graph_Assembly.md §2, §3, §6)
/// ───────────────────────────────────────────────────────────────────────────
/// BuildDeBruijnGraph: every length-k substring r[i..i+k] of a read is a directed
/// edge from its left (k-1)-mer r[i..i+k-1] (prefix) to its right (k-1)-mer
/// r[i+1..i+k] (suffix); nodes are the distinct (k-1)-mers; the graph is a MULTIGRAPH
/// (a repeated k-mer → a repeated edge, INV-02). Total #edges = total #k-mers chopped
/// (INV-02). A read of length &lt; k contributes ZERO edges (INV-03). k &lt; 2 throws
/// ArgumentOutOfRangeException (§3.3, §6.1); null/empty reads → empty graph (§3.3).
///   SequenceAssembler.BuildDeBruijnGraph(IReadOnlyList&lt;string&gt; reads, int k)
///   → Dictionary&lt;string, List&lt;string&gt;&gt;  (out-adjacency multimap)
///
/// AssembleDeBruijn: builds the DBG, walks each weakly-connected component with one
/// Eulerian walk (Hierholzer, semi-balanced source else lex-smallest node), spells
/// p₀ + last-char(p_i) into one contig per component, then drops contigs shorter than
/// MinContigLength. When a component's Eulerian walk is unique the contig equals the
/// source genome (INV-04); every input k-mer is a substring of some emitted contig
/// (INV-05). null/empty reads → empty result (§6.1). The choice of component order
/// and start node is the lexicographically smallest node → DETERMINISTIC, order- and
/// count-independent output (§4.2, §5.2).
///   SequenceAssembler.AssembleDeBruijn(
///       IReadOnlyList&lt;string&gt; reads, AssemblyParameters? parameters = null)
///   → AssemblyResult (Contigs, TotalReads, N50, LongestContig, TotalLength)
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyDeBruijnFuzzTests
{
    // Canonical DNA alphabet for fuzzed reads (k-mers are taken verbatim, case-sensitive — §3.3).
    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A fuzzed read of the given length over the DNA alphabet.</summary>
    private static string RandomRead(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaAlphabet[rng.Next(DnaAlphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle for the EDGE COUNT of <see cref="SequenceAssembler.BuildDeBruijnGraph"/>:
    /// the total number of edges equals the total number of k-mers chopped from all reads,
    /// Σ max(0, |r| − k + 1) (INV-02; INV-03 makes reads shorter than k contribute 0).
    /// </summary>
    private static int ExpectedEdgeCount(IReadOnlyList<string> reads, int k)
        => reads.Where(r => r != null).Sum(r => Math.Max(0, r.Length - k + 1));

    /// <summary>Total edges actually stored in the out-adjacency multimap.</summary>
    private static int EdgeCount(Dictionary<string, List<string>> graph)
        => graph.Values.Sum(list => list.Count);

    /// <summary>
    /// Asserts every emitted contig is WELL-FORMED per the documented contract,
    /// regardless of the (possibly degenerate) input:
    ///   • length ≥ k−1 — a spelled walk emits at least its first (k-1)-mer node (§4.1);
    ///   • composed only of characters present in the supplied reads (contigs are
    ///     spelled from (k-1)-mer nodes, which are read substrings — §2.2);
    ///   • the recorded statistics are consistent: TotalLength = Σ contig lengths,
    ///     LongestContig = max contig length (INV-06).
    /// </summary>
    private static void AssertWellFormed(
        SequenceAssembler.AssemblyResult result, IReadOnlyList<string> reads, int k)
    {
        var alphabet = reads.Where(r => r != null).SelectMany(r => r).ToHashSet();

        foreach (string contig in result.Contigs)
        {
            contig.Length.Should().BeGreaterThanOrEqualTo(k - 1,
                "a spelled Eulerian walk emits at least its first (k-1)-mer node (§4.1)");
            contig.All(ch => alphabet.Contains(ch)).Should().BeTrue(
                "a contig is spelled from (k-1)-mer nodes that are read substrings (§2.2)");
        }

        result.TotalLength.Should().Be(result.Contigs.Sum(c => c.Length),
            "INV-06: TotalLength = Σ contig lengths");
        if (result.Contigs.Count > 0)
        {
            result.LongestContig.Should().Be(result.Contigs.Max(c => c.Length),
                "INV-06: LongestContig = max contig length");
        }
    }

    #endregion

    #region ASSEMBLY-DBG-001 — De Bruijn Graph Assembly (BE: k>read length, single read, all-identical reads)

    #region Positive sanity — documented Eulerian reconstruction

    // Documented worked example (§7.1): AAABBBA, k=3 → single contig "AAABBBA".
    [Test]
    public void AssembleDeBruijn_WorkedExample_ReconstructsGenome()
    {
        var reads = new[] { "AAABBBA" };

        var result = SequenceAssembler.AssembleDeBruijn(
            reads, new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));

        result.Contigs.Should().ContainSingle("the graph is connected and Eulerian → one contig (§7.1)");
        result.Contigs[0].Should().Be("AAABBBA", "the unique Eulerian walk spells the source genome (INV-04, §7.1)");
        AssertWellFormed(result, reads, 3);
    }

    // A single read spanning a known genome whose (k-1)-mer walk is UNIQUE (no repeated (k-1)-mer
    // anywhere) reconstructs the genome verbatim as one contig, and every input k-mer is a
    // contiguous substring of it (INV-04, INV-05). This is the canonical exact-reconstruction case:
    // each k-mer occurs exactly once, so the Eulerian walk is unique. (Note: tiling the SAME genome
    // with multiple OVERLAPPING reads instead makes each shared k-mer a parallel multiedge — many
    // Eulerian walks then exist and the contig is no longer the genome, per ASSUMPTION-1 / §6.2;
    // that documented multigraph behavior is exercised by the all-identical-reads tests below.)
    [Test]
    public void AssembleDeBruijn_UniqueWalkGenome_ReconstructsExactlyWithKmerCoverage()
    {
        const string genome = "GACAAATCATGGTTAACATCCCACCTGAGC"; // (k-1)=4-mer walk is unique at k=5
        const int k = 5;
        var reads = new[] { genome };

        var result = SequenceAssembler.AssembleDeBruijn(
            reads, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

        result.Contigs.Should().ContainSingle("a unique Eulerian walk yields one connected contig");
        result.Contigs[0].Should().Be(genome, "the unique Eulerian walk reconstructs the genome (INV-04)");

        // INV-05: each input k-mer is a contiguous substring of the (exactly reconstructed) contig.
        for (int i = 0; i + k <= genome.Length; i++)
            result.Contigs[0].Should().Contain(genome.Substring(i, k), "INV-05: every k-mer is covered");

        AssertWellFormed(result, reads, k);
    }

    // BuildDeBruijnGraph: documented node/edge construction ((k-1)-mer nodes, k-mer edges, multigraph).
    [Test]
    public void BuildDeBruijnGraph_WorkedExample_ProducesDocumentedNodesAndEdges()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AAABBBA" }, 3);

        // Edges: AA→AA, AA→AB, AB→BB, BB→BB, BB→BA  (5 k-mers → 5 edges, INV-02).
        EdgeCount(graph).Should().Be(5, "INV-02: #edges = #k-mers chopped");
        graph["AA"].Should().BeEquivalentTo("AA", "AB");
        graph["AB"].Should().BeEquivalentTo("BB");
        graph["BB"].Should().BeEquivalentTo("BB", "BA");
    }

    #endregion

    #region BE — Boundary: k > read length (no k-mers; empty graph; no negative substring / IndexOutOfRange)

    // k strictly greater than the read length → chop bound i ∈ [0, |r|−k] is empty → empty graph.
    [Test]
    public void BuildDeBruijnGraph_KGreaterThanReadLength_EmptyGraph_NoCrash()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AC" }, 5);

        graph.Should().BeEmpty("a read shorter than k contributes no k-mers (INV-03), no negative-length Substring");
    }

    // k > read length for the assembler → no edges → no contigs, empty result, no IndexOutOfRange.
    [Test]
    public void AssembleDeBruijn_KGreaterThanReadLength_EmptyResult()
    {
        var reads = new[] { "ACGTA" };

        var result = SequenceAssembler.AssembleDeBruijn(
            reads, new SequenceAssembler.AssemblyParameters(KmerSize: 10, MinContigLength: 1));

        result.Contigs.Should().BeEmpty("no read reaches length k → no k-mers → no contigs (INV-03)");
        result.TotalReads.Should().Be(1, "TotalReads still reflects the input count");
        AssertWellFormed(result, reads, 10);
    }

    // k exactly one larger than EVERY read (the immediate boundary) → still empty, no off-by-one crash.
    [Test]
    public void BuildDeBruijnGraph_KOneMoreThanLongestRead_EmptyGraph()
    {
        var reads = new[] { "ACGT", "AC", "ACGTA" }; // longest = 5
        var graph = SequenceAssembler.BuildDeBruijnGraph(reads, 6);

        graph.Should().BeEmpty("k = longestRead + 1 → no read can be chopped (boundary of INV-03)");
    }

    // Fuzz: for ANY random reads, choosing k strictly larger than every read never throws and yields
    // an empty graph and an empty assembly (the k > read-length boundary, swept over many shapes).
    [Test]
    public void AssembleDeBruijn_KAboveAllReadLengths_AlwaysEmpty_NeverThrows()
    {
        var rng = new Random(143_001);
        for (int trial = 0; trial < 600; trial++)
        {
            int count = rng.Next(1, 6);
            var reads = Enumerable.Range(0, count).Select(_ => RandomRead(rng, rng.Next(0, 12))).ToList();
            int maxLen = reads.Max(r => r.Length);
            // Strictly above every read AND ≥ 2 (the k < 2 boundary throws and is tested separately).
            int k = Math.Max(2, maxLen + rng.Next(1, 5));

            var graph = SequenceAssembler.BuildDeBruijnGraph(reads, k);
            graph.Should().BeEmpty("every read is shorter than k → no k-mers (INV-03)");

            var result = SequenceAssembler.AssembleDeBruijn(
                reads, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));
            result.Contigs.Should().BeEmpty("no k-mers → no contigs");
            AssertWellFormed(result, reads, k);
        }
    }

    // k < 2 is the documented validation boundary: ArgumentOutOfRangeException (NOT silent / not a crash).
    [Test]
    public void BuildDeBruijnGraph_KBelowTwo_ThrowsArgumentOutOfRange()
    {
        Action act = () => SequenceAssembler.BuildDeBruijnGraph(new[] { "ACGT" }, 1);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "k must be ≥ 2 so (k-1)-mer nodes are non-empty (§3.3, §6.1)");
    }

    #endregion

    #region BE — Boundary: single read (simple path ⇒ one contig reconstructing the read)

    // A single read with a UNIQUE (k-1)-mer walk → one contig equal to the read (INV-04).
    [Test]
    public void AssembleDeBruijn_SingleReadUniqueWalk_ReconstructsRead()
    {
        const string read = "ACGTACTGACATG"; // distinct 3-mers ⇒ unique Eulerian path at k=4
        const int k = 4;

        var result = SequenceAssembler.AssembleDeBruijn(
            new[] { read }, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

        result.Contigs.Should().ContainSingle("a single read is one connected component → one contig");
        result.Contigs[0].Should().Be(read, "a unique Eulerian walk reconstructs the read verbatim (INV-04)");
        AssertWellFormed(result, new[] { read }, k);
    }

    // A single read of length EXACTLY k → exactly one k-mer → a single edge → contig == the read.
    [Test]
    public void AssembleDeBruijn_SingleReadLengthEqualsK_ReconstructsRead()
    {
        const string read = "ACGT";
        const int k = 4;

        var result = SequenceAssembler.AssembleDeBruijn(
            new[] { read }, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

        result.Contigs.Should().ContainSingle("one k-mer → one edge → one contig");
        result.Contigs[0].Should().Be(read, "the lone edge spells the whole read (INV-05)");
    }

    // Fuzz: a single random read with no repeated (k-1)-mer reconstructs itself; every input
    // k-mer is a substring of the contig (INV-05); never throws, always well-formed.
    [Test]
    public void AssembleDeBruijn_SingleRandomRead_ReconstructsOrCoversKmers()
    {
        var rng = new Random(143_002);
        int reconstructed = 0;
        for (int trial = 0; trial < 800; trial++)
        {
            int k = rng.Next(2, 8);
            int len = rng.Next(k, k + 20); // length ≥ k so at least one k-mer exists
            string read = RandomRead(rng, len);

            var result = SequenceAssembler.AssembleDeBruijn(
                new[] { read }, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

            AssertWellFormed(result, new[] { read }, k);

            // INV-05: every input k-mer must appear in SOME emitted contig.
            for (int i = 0; i + k <= read.Length; i++)
            {
                string kmer = read.Substring(i, k);
                result.Contigs.Any(c => c.Contains(kmer)).Should().BeTrue(
                    "INV-05: every input k-mer is a substring of an emitted contig");
            }

            // When the (k-1)-mers along the read are all distinct, the walk is unique → exact read.
            var nodes = Enumerable.Range(0, read.Length - (k - 1) + 1)
                                  .Select(i => read.Substring(i, k - 1)).ToList();
            if (nodes.Distinct().Count() == nodes.Count)
            {
                result.Contigs.Should().ContainSingle("unique (k-1)-mer walk → one contig");
                result.Contigs[0].Should().Be(read, "unique Eulerian walk reconstructs the read (INV-04)");
                reconstructed++;
            }
        }

        reconstructed.Should().BeGreaterThan(0, "the fuzz corpus must include unique-walk reads (sanity on coverage)");
    }

    #endregion

    #region BE — Boundary: all-identical reads (collapsed graph; deterministic; no hang / no infinite loop)

    // Many copies of the SAME read collapse onto the SAME (k-1)-mer nodes (the node set is shared),
    // so the graph is ONE weakly-connected component regardless of copy count. Edge multiplicity
    // scales linearly with the copy count (INV-02 multigraph). The assembly must be DETERMINISTIC
    // (a fixed input → byte-identical output) and must TERMINATE on the resulting cyclic-looking
    // multigraph (Hierholzer is O(|E|), §4.3). NOTE: with ≥ 2 copies each k-mer is a PARALLEL
    // multiedge, so several Eulerian walks exist and the spelled contig is NOT the original read
    // (ASSUMPTION-1, §6.2) — only that it is well-formed and deterministic is guaranteed.
    [Test]
    [CancelAfter(20_000)]
    public void AssembleDeBruijn_AllIdenticalReads_OneComponentDeterministicTerminates()
    {
        const string read = "ACGTACTGAC"; // distinct (k-1)-mers at k=4 (1 copy reconstructs exactly)
        const int k = 4;
        var prm = new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1);

        for (int copies = 1; copies <= 50; copies++)
        {
            var reads = Enumerable.Repeat(read, copies).ToList();

            var graph = SequenceAssembler.BuildDeBruijnGraph(reads, k);
            EdgeCount(graph).Should().Be(copies * (read.Length - k + 1),
                "INV-02: identical reads add repeated (multi-)edges, edge count = copies × #k-mers");

            var result = SequenceAssembler.AssembleDeBruijn(reads, prm);

            result.Contigs.Should().ContainSingle("identical reads share all nodes → one connected component");
            AssertWellFormed(result, reads, k);

            // A single copy is a SIMPLE path → the read is reconstructed exactly (INV-04).
            if (copies == 1)
                result.Contigs[0].Should().Be(read, "1 copy is a unique walk → exact reconstruction (INV-04)");

            // Determinism: re-running on the identical input is byte-identical.
            SequenceAssembler.AssembleDeBruijn(reads, prm).Contigs
                .Should().Equal(result.Contigs, "a fixed input produces byte-identical contigs (§5.2)");
        }
    }

    // Determinism across read ORDER: shuffling identical (and merely duplicated) reads must not
    // change the emitted contig set (lex-smallest start node, §4.2, §5.2).
    [Test]
    [CancelAfter(20_000)]
    public void AssembleDeBruijn_DuplicatedReads_OrderIndependentContigs()
    {
        var rng = new Random(143_003);
        for (int trial = 0; trial < 300; trial++)
        {
            int k = rng.Next(3, 6);
            string read = RandomRead(rng, rng.Next(k, k + 12));
            int copies = rng.Next(2, 8);
            var reads = Enumerable.Repeat(read, copies).ToList();

            var baseline = SequenceAssembler.AssembleDeBruijn(
                reads, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1)).Contigs;

            for (int shuffle = 0; shuffle < 3; shuffle++)
            {
                var permuted = reads.OrderBy(_ => rng.Next()).ToList();
                var contigs = SequenceAssembler.AssembleDeBruijn(
                    permuted, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1)).Contigs;

                contigs.Should().BeEquivalentTo(baseline,
                    "the contig set is deterministic and order-independent (§4.2, §5.2)");
            }
        }
    }

    // An all-identical read whose (k-1)-mers REPEAT internally (a cyclic / repeat-rich read) is the
    // worst case for the Eulerian traversal — it MUST terminate (no infinite loop on the cycle).
    [Test]
    [CancelAfter(20_000)]
    public void AssembleDeBruijn_IdenticalRepeatRichReads_TerminateAndCoverKmers()
    {
        // Tandem repeat: (k-1)-mer "AAA" etc. recur → multiple parallel edges and a cycle.
        const string read = "ACACACACACACAC"; // 2-periodic ⇒ repeated (k-1)-mers at k=4
        const int k = 4;
        var reads = Enumerable.Repeat(read, 30).ToList();

        var prm = new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1);
        var result = SequenceAssembler.AssembleDeBruijn(reads, prm);

        AssertWellFormed(result, reads, k);
        result.Contigs.Should().NotBeEmpty("a repeat-rich read still yields at least one Eulerian contig");
        // The cyclic, multi-parallel-edge graph must not hang and must be deterministic for a fixed input.
        SequenceAssembler.AssembleDeBruijn(reads, prm).Contigs
            .Should().Equal(result.Contigs, "a fixed repeat-rich input is deterministic (§5.2), no hang on the cycle");
    }

    #endregion

    #region BE — Boundary: empty / null reads (documented empty result, no crash)

    [Test]
    public void AssembleDeBruijn_EmptyReads_EmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(
            Array.Empty<string>(), new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));

        result.Contigs.Should().BeEmpty("empty reads → empty assembly (§6.1)");
        result.TotalReads.Should().Be(0);
    }

    [Test]
    public void AssembleDeBruijn_NullReads_EmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(
            null!, new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));

        result.Contigs.Should().BeEmpty("null reads → empty assembly (§6.1, ASSUMPTION-2)");
    }

    [Test]
    public void BuildDeBruijnGraph_NullAndEmptyReads_EmptyGraph()
    {
        SequenceAssembler.BuildDeBruijnGraph(null!, 3).Should().BeEmpty("null reads → empty graph (§3.3)");
        SequenceAssembler.BuildDeBruijnGraph(Array.Empty<string>(), 3).Should().BeEmpty("empty reads → empty graph");
    }

    // Reads that are all empty strings (length 0 < k) contribute nothing — boundary of INV-03.
    [Test]
    public void AssembleDeBruijn_AllEmptyStringReads_EmptyResult()
    {
        var reads = new[] { "", "", "" };

        var result = SequenceAssembler.AssembleDeBruijn(
            reads, new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1));

        result.Contigs.Should().BeEmpty("zero-length reads are shorter than k → no k-mers (INV-03)");
        result.TotalReads.Should().Be(3, "TotalReads still counts the (empty) inputs");
    }

    #endregion

    #region BE — Broad fuzz: arbitrary reads × k never crash, always well-formed, always terminate

    [Test]
    [CancelAfter(30_000)]
    public void AssembleDeBruijn_RandomReadsAndK_NeverThrows_WellFormed_Terminates()
    {
        var rng = new Random(143_004);
        for (int trial = 0; trial < 1500; trial++)
        {
            int k = rng.Next(2, 9);
            int count = rng.Next(0, 8);
            // Read lengths straddle k (some shorter → skipped, some longer → contribute k-mers).
            var reads = Enumerable.Range(0, count)
                                  .Select(_ => RandomRead(rng, rng.Next(0, k + 14)))
                                  .ToList();

            var graph = SequenceAssembler.BuildDeBruijnGraph(reads, k);
            EdgeCount(graph).Should().Be(ExpectedEdgeCount(reads, k),
                "INV-02: total edges = total k-mers chopped (reads < k contribute 0, INV-03)");

            var result = SequenceAssembler.AssembleDeBruijn(
                reads, new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

            AssertWellFormed(result, reads, k);
            result.TotalReads.Should().Be(count, "TotalReads reflects the input count");
        }
    }

    // Determinism over arbitrary reads: repeated invocation on the SAME input is byte-identical
    // (§5.2 — the component order and start node are the lex-smallest node, so a fixed input always
    // produces the same contigs). NOTE: across a read-order PERMUTATION the contig set may differ
    // when the graph has a repeated (k-1)-mer (multiple Eulerian walks; edges are consumed in
    // insertion order, which depends on read order — §5.2, ASSUMPTION-1); so cross-permutation
    // invariance is asserted only for inputs that share all nodes by construction (identical reads,
    // covered by AssembleDeBruijn_DuplicatedReads_OrderIndependentContigs).
    [Test]
    [CancelAfter(30_000)]
    public void AssembleDeBruijn_ArbitraryReads_FixedInputDeterministic()
    {
        var rng = new Random(143_005);
        for (int trial = 0; trial < 400; trial++)
        {
            int k = rng.Next(2, 7);
            int count = rng.Next(1, 7);
            var reads = Enumerable.Range(0, count)
                                  .Select(_ => RandomRead(rng, rng.Next(0, k + 10)))
                                  .ToList();
            var prm = new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1);

            var a = SequenceAssembler.AssembleDeBruijn(reads, prm).Contigs;
            var b = SequenceAssembler.AssembleDeBruijn(reads, prm).Contigs;
            a.Should().Equal(b, "the same input produces byte-identical contigs (deterministic, §5.2)");
        }
    }

    #endregion

    #endregion
}
