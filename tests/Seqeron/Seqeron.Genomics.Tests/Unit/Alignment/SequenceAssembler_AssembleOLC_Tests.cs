// ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly
// Evidence: docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-OLC-001.md
// Source: Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987-991, DOI 10.1038/nbt.2023;
//         Langmead B., "Overlap Layout Consensus assembly" / "Assembly & Shortest Common
//         Superstring" (JHU lecture notes).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_AssembleOLC_Tests
{
    // The 6 distinct 6-mers of GTACGTACGAT, in genome order.
    // Source: Langmead, "Assembly & Shortest Common Superstring", p.24-25.
    private static readonly string[] GtacgtacgatSixMers =
        { "GTACGT", "TACGTA", "ACGTAC", "CGTACG", "GTACGA", "TACGAT" };

    // The 12 directed overlap-graph edges (from-6mer, to-6mer, overlap length) for the
    // GTACGTACGAT 6-mers at minOverlap 4, derived from the longest suffix-prefix definition
    // and matching the edge weights (4, 5) drawn in Langmead SCS p.24-25.
    private static readonly (string From, string To, int Len)[] GtacgtacgatEdges =
    {
        ("GTACGT", "TACGTA", 5), ("GTACGT", "ACGTAC", 4),
        ("TACGTA", "ACGTAC", 5), ("TACGTA", "CGTACG", 4),
        ("ACGTAC", "GTACGT", 4), ("ACGTAC", "CGTACG", 5), ("ACGTAC", "GTACGA", 4),
        ("CGTACG", "GTACGT", 5), ("CGTACG", "TACGTA", 4), ("CGTACG", "GTACGA", 5), ("CGTACG", "TACGAT", 4),
        ("GTACGA", "TACGAT", 5),
    };

    #region FindAllOverlaps

    // M1 — FindAllOverlaps on the GTACGTACGAT 6-mers (minOverlap 4, identity 1.0) must produce
    // exactly the 12 directed edges with the published overlap lengths. Source: Langmead SCS p.24-25.
    [Test]
    public void FindAllOverlaps_GtacgtacgatSixMers_ReturnsExactTwelveEdgeGraph()
    {
        var reads = GtacgtacgatSixMers.ToList();

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap: 4, minIdentity: 1.0);

        // Map each edge to (fromString, toString, length) for evidence comparison.
        var actual = overlaps
            .Select(o => (From: reads[o.ReadIndex1], To: reads[o.ReadIndex2], Len: o.OverlapLength))
            .OrderBy(e => e.From).ThenBy(e => e.To).ThenBy(e => e.Len)
            .ToList();

        var expected = GtacgtacgatEdges
            .OrderBy(e => e.From).ThenBy(e => e.To).ThenBy(e => e.Len)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(actual.Count, Is.EqualTo(12),
                "GTACGTACGAT 6-mers at minOverlap 4 form exactly 12 directed overlap edges (Langmead SCS p.24-25).");
            Assert.That(actual, Is.EqualTo(expected),
                "Each directed edge and its overlap length must match the published overlap graph exactly.");
        });
    }

    // M2 — The overlap graph must never contain a self-edge (INV-01). Source: Compeau/Pevzner; Langmead OLC p.5.
    [Test]
    public void FindAllOverlaps_SixMers_ContainsNoSelfOverlap()
    {
        var reads = GtacgtacgatSixMers.ToList();

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap: 4, minIdentity: 1.0);

        Assert.That(overlaps.All(o => o.ReadIndex1 != o.ReadIndex2), Is.True,
            "INV-01: a read is never overlapped against itself; the overlap graph has no self-edges.");
    }

    // S3 — Reads sharing only a 3-base suffix-prefix overlap must yield no edges at minOverlap 4 (INV-02).
    [Test]
    public void FindAllOverlaps_OverlapBelowThreshold_ReturnsNoEdges()
    {
        // "ACGTACGT" suffix "TAC..." vs "CGTAAAAA": longest suffix-prefix match is 3 ("CGT"), below 4.
        var reads = new List<string> { "ACGTACGT", "CGTAAAAA" };

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap: 4, minIdentity: 1.0);

        Assert.That(overlaps.Count, Is.EqualTo(0),
            "INV-02: an overlap shorter than minOverlap is not an edge.");
    }

    // C2 — Overlap detection is case-insensitive; lowercase reads give the same edge set as uppercase.
    [Test]
    public void FindAllOverlaps_LowercaseReads_SameEdgesAsUppercase()
    {
        var upper = GtacgtacgatSixMers.ToList();
        var lower = GtacgtacgatSixMers.Select(s => s.ToLowerInvariant()).ToList();

        var upperEdges = SequenceAssembler.FindAllOverlaps(upper, minOverlap: 4, minIdentity: 1.0)
            .Select(o => (o.ReadIndex1, o.ReadIndex2, o.OverlapLength)).OrderBy(e => e).ToList();
        var lowerEdges = SequenceAssembler.FindAllOverlaps(lower, minOverlap: 4, minIdentity: 1.0)
            .Select(o => (o.ReadIndex1, o.ReadIndex2, o.OverlapLength)).OrderBy(e => e).ToList();

        Assert.That(lowerEdges, Is.EqualTo(upperEdges),
            "Identity is computed case-insensitively, so case does not change the overlap graph.");
    }

    // S5 — The cancellable FindAllOverlaps overload (CancellationToken.None) returns the same edges.
    [Test]
    public void FindAllOverlaps_CancellableOverload_SameResultAsBasic()
    {
        var reads = GtacgtacgatSixMers.ToList();

        var basic = SequenceAssembler.FindAllOverlaps(reads, 4, 1.0)
            .Select(o => (o.ReadIndex1, o.ReadIndex2, o.OverlapLength)).OrderBy(e => e).ToList();
        var cancellable = SequenceAssembler.FindAllOverlaps(reads, 4, 1.0, CancellationToken.None)
            .Select(o => (o.ReadIndex1, o.ReadIndex2, o.OverlapLength)).OrderBy(e => e).ToList();

        Assert.That(cancellable, Is.EqualTo(basic),
            "The cancellable overload must compute the identical overlap graph (delegation smoke test).");
    }

    // S6 — The cancellable overload throws OperationCanceledException for an already-cancelled token.
    [Test]
    public void FindAllOverlaps_CancellableOverload_AlreadyCancelled_Throws()
    {
        var reads = GtacgtacgatSixMers.ToList();
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.That(() => SequenceAssembler.FindAllOverlaps(reads, 4, 1.0, cts.Token),
            NUnit.Framework.Throws.InstanceOf<System.OperationCanceledException>(),
            "An already-cancelled token must abort overlap detection.");
    }

    #endregion

    #region FindOverlap

    // M3 — FindOverlap reports the single longest suffix-prefix overlap and its 0-based positions.
    // Source: Langmead OLC p.5 (CTCTAGGCC / TAGGCCCTC share the length-6 suffix-prefix "TAGGCC").
    [Test]
    public void FindOverlap_LongestSuffixPrefix_ReturnsLengthAndPositions()
    {
        var overlap = SequenceAssembler.FindOverlap("CTCTAGGCC", "TAGGCCCTC", minOverlap: 3, minIdentity: 1.0);

        Assert.That(overlap, Is.Not.Null, "A length-6 suffix-prefix overlap exists and must be found.");
        Assert.Multiple(() =>
        {
            Assert.That(overlap!.Value.length, Is.EqualTo(6),
                "Longest suffix of 'CTCTAGGCC' equal to a prefix of 'TAGGCCCTC' is 'TAGGCC' (length 6).");
            Assert.That(overlap.Value.pos1, Is.EqualTo(3),
                "pos1 is the 0-based start of the overlapping suffix in seq1 (9 - 6 = 3).");
            Assert.That(overlap.Value.pos2, Is.EqualTo(0),
                "pos2 is always 0: the overlap is a prefix of seq2.");
        });
    }

    // S1 — Identity threshold gates overlap acceptance: 7/8 = 0.875 accepted at 0.85, rejected at 0.95.
    // Source: Langmead OLC p.11-15 (approximate overlap via the identity fraction).
    [Test]
    public void FindOverlap_OneMismatchInEight_RespectsIdentityThreshold()
    {
        // seq1 suffix "ACGTACGT" vs seq2 prefix "ACGTACCT": 1 mismatch -> 7/8 = 0.875 identity.
        string seq1 = "ACGTACGTACGT";
        string seq2 = "ACGTACCTAAAA";

        var accepted = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 8, minIdentity: 0.85);
        var rejected = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 8, minIdentity: 0.95);

        Assert.Multiple(() =>
        {
            Assert.That(accepted, Is.Not.Null,
                "0.875 identity >= 0.85 threshold: the overlap is accepted.");
            Assert.That(accepted!.Value.length, Is.EqualTo(8),
                "The accepted overlap spans the full 8-base window.");
            Assert.That(rejected, Is.Null,
                "0.875 identity < 0.95 threshold: the overlap is rejected.");
        });
    }

    // S2 — An overlap exactly equal to minOverlap is accepted; one base shorter is rejected (INV-02).
    [Test]
    public void FindOverlap_MinOverlapBoundary_AcceptsAtThresholdRejectsBelow()
    {
        // Longest suffix-prefix match between "AAAACGTT" and "CGTTGGGG" is "CGTT" (length 4).
        string seq1 = "AAAACGTT";
        string seq2 = "CGTTGGGG";

        var atThreshold = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 4, minIdentity: 1.0);
        var aboveThreshold = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 5, minIdentity: 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(atThreshold, Is.Not.Null, "Overlap length 4 == minOverlap 4 is accepted.");
            Assert.That(atThreshold!.Value.length, Is.EqualTo(4), "The qualifying overlap is exactly 4.");
            Assert.That(aboveThreshold, Is.Null, "No suffix-prefix overlap of length >= 5 exists, so none is reported.");
        });
    }

    #endregion

    #region AssembleOLC

    // M4 — An unambiguous 5-overlap tiling reconstructs a single contig that is the superstring of
    // all reads. Source: Langmead OLC p.5 + chain consensus; INV-04.
    [Test]
    public void AssembleOLC_UnambiguousChain_ProducesSingleSuperstringContig()
    {
        var reads = new List<string> { "AAAAACCCCC", "CCCCCGGGGG", "GGGGGTTTTT" };

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, MinContigLength: 10));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1),
                "Three reads forming one unambiguous overlap chain collapse to a single contig.");
            Assert.That(result.Contigs[0], Is.EqualTo("AAAAACCCCCGGGGGTTTTT"),
                "Merging along the chain (A + B[overlap:]) yields the 20-base superstring of all reads.");
            Assert.That(result.TotalReads, Is.EqualTo(3), "TotalReads equals the input read count.");
            Assert.That(result.LongestContig, Is.EqualTo(20), "The single contig is 20 bases long.");
        });
    }

    // M5 — Three non-overlapping reads form an edgeless graph; each read is its own contig (INV-05).
    // Source: Compeau/Pevzner overlap-graph definition (no edge below threshold).
    [Test]
    public void AssembleOLC_NoOverlaps_ReturnsSingletonContigs()
    {
        var reads = new List<string> { "AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG" };

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, MinContigLength: 5));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(3),
                "INV-05: with no above-threshold overlap, each read is its own singleton contig.");
            Assert.That(result.Contigs.OrderBy(c => c), Is.EqualTo(reads.OrderBy(c => c)),
                "Each singleton contig equals an input read verbatim.");
            Assert.That(result.TotalLength, Is.EqualTo(30), "Total length is the sum of the three 10-base reads.");
        });
    }

    // M5b — MinContigLength discards short contigs: a singleton read below the threshold is
    // dropped while one at/above it is kept (algorithm doc §3.2: "Contigs shorter than this are
    // discarded"). Two non-overlapping reads (5 and 10 bases) at MinContigLength 8 -> only the
    // 10-base read survives.
    [Test]
    public void AssembleOLC_MinContigLength_DiscardsShortContigs()
    {
        var reads = new List<string> { "AAAAA", "CCCCCCCCCC" };

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, MinContigLength: 8));

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(1),
                "The 5-base read is shorter than MinContigLength 8 and must be discarded.");
            Assert.That(result.Contigs[0], Is.EqualTo("CCCCCCCCCC"),
                "Only the 10-base read meets the MinContigLength 8 threshold.");
            Assert.That(result.TotalReads, Is.EqualTo(2),
                "TotalReads reflects the input read count regardless of contig filtering.");
            Assert.That(result.TotalLength, Is.EqualTo(10),
                "TotalLength sums only the surviving contigs (the 5-base read was discarded).");
        });
    }

    // M6 — Empty read set returns an empty AssemblyResult (trivial identity; ASSUMPTION-2).
    [Test]
    public void AssembleOLC_EmptyReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleOLC(new List<string>());

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs.Count, Is.EqualTo(0), "No reads -> no contigs.");
            Assert.That(result.TotalReads, Is.EqualTo(0), "No reads -> TotalReads 0.");
            Assert.That(result.TotalLength, Is.EqualTo(0), "No reads -> TotalLength 0.");
        });
    }

    // M6b — Null read set is handled like empty (contract; no exception).
    [Test]
    public void AssembleOLC_NullReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleOLC(null!);

        Assert.That(result.Contigs.Count, Is.EqualTo(0), "Null input is treated as empty, returning no contigs.");
    }

    // S4 — INV-04 property: for the unambiguous chain, the contig length lies within
    // [longest read length, sum of read lengths]. Source: superstring property (Langmead SCS p.26).
    [Test]
    public void AssembleOLC_UnambiguousChain_ContigLengthWithinBounds()
    {
        var reads = new List<string> { "AAAAACCCCC", "CCCCCGGGGG", "GGGGGTTTTT" };
        int sumLen = reads.Sum(r => r.Length);   // 30
        int longest = reads.Max(r => r.Length);  // 10

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, MinContigLength: 10));

        Assert.Multiple(() =>
        {
            foreach (var contig in result.Contigs)
            {
                Assert.That(contig.Length, Is.GreaterThanOrEqualTo(longest),
                    "INV-04: a merged contig is at least as long as the longest single read.");
                Assert.That(contig.Length, Is.LessThanOrEqualTo(sumLen),
                    "INV-04: a merged contig is no longer than the concatenation of its reads.");
            }
            // And every input read appears as a substring of some emitted contig.
            Assert.That(reads.All(r => result.Contigs.Any(c => c.Contains(r))), Is.True,
                "INV-04: each input read is a substring of an emitted contig.");
        });
    }

    // C1 — Repeat limitation (ASM-02): reads with an internal repeat are not collapsed below the
    // longest read length; the lower bound of INV-04 still holds. Source: Langmead SCS p.58-62.
    [Test]
    public void AssembleOLC_RepeatContainingReads_DoesNotCollapseBelowLongestRead()
    {
        // Reads tiling "ACGACGACGTTT" where "ACG" repeats; the period (3) < read length.
        var reads = new List<string> { "ACGACGACG", "GACGTTT", "ACGACGTTT" };
        int longest = reads.Max(r => r.Length);

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 4, MinIdentity: 1.0, MinContigLength: 1));

        Assert.That(result.Contigs.All(c => c.Length >= longest), Is.True,
            "ASM-02/INV-04: even with an internal repeat the assembler never emits a contig shorter than the longest read.");
    }

    #endregion

    #region 2026-09 review — GREEDY layout, containment, consensus, threshold validation

    // Reference values below were produced by the Langmead greedy_scs code (Coursera
    // "Algorithms for DNA Sequencing", mirrored in kywertheim/Greedy_shortest_common_superstring
    // main.py; list-of-contigs variant) and an independent edge-GREEDY + Biopython 1.79
    // dumb_consensus re-implementation (scratch ref_olc.py); both agree on every exact case.

    private static SequenceAssembler.AssemblyParameters Exact(int minOverlap) =>
        new(MinOverlap: minOverlap, MinIdentity: 1.0, MinContigLength: 1);

    // F1 — GREEDY must not close a cycle. The GTACGTACGAT 6-mers' best successors form the cycle
    // GTACGT→TACGTA→ACGTAC→CGTACG→GTACGT; GREEDY (Blum et al. 1994) rejects the cycle-closing
    // edge and takes CGTACG→GTACGA instead, spelling the whole genome. greedy_scs(k=4) → ["GTACGTACGAT"].
    [Test]
    public void AssembleOLC_GtacgtacgatSixMers_GreedyRejectsCycle_ReconstructsGenome()
    {
        var result = SequenceAssembler.AssembleOLC(GtacgtacgatSixMers.ToList(), Exact(4));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "GTACGTACGAT" }),
            "GREEDY layout of the 6-mers (l = 4) reconstructs GTACGTACGAT (Langmead greedy_scs reference).");
    }

    // F1 — An edge whose head already has a predecessor is skipped, so the tail can take its
    // next-best successor: TTTTCATGCA→CATGCAAAAA (6) wins; GGGGGATGCA→CATGCAAAAA (5) is rejected and
    // GGGGGATGCA→TGCACCCCCC (4) is taken. greedy_scs(k=4) → 2 contigs.
    [Test]
    public void AssembleOLC_TakenHead_TailUsesNextBestSuccessor()
    {
        var reads = new List<string> { "TTTTCATGCA", "GGGGGATGCA", "CATGCAAAAA", "TGCACCCCCC" };

        var result = SequenceAssembler.AssembleOLC(reads, Exact(4));

        Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal),
            Is.EqualTo(new[] { "GGGGGATGCACCCCCC", "TTTTCATGCAAAAA" }),
            "GREEDY takes the 4-overlap GGGGGATGCA→TGCACCCCCC once CATGCAAAAA already has a predecessor.");
    }

    // F1/F2 — Duplicate reads (mutual full-length overlaps) must not create a 2-cycle that drops them
    // out of the layout. greedy_scs(k=5) → ["AAAAACCCCCGGGGG"].
    [Test]
    public void AssembleOLC_DuplicateReads_AssembleIntoOneContig()
    {
        var reads = new List<string> { "AAAAACCCCC", "AAAAACCCCC", "CCCCCGGGGG" };

        var result = SequenceAssembler.AssembleOLC(reads, Exact(5));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "AAAAACCCCCGGGGG" }),
            "An identical copy is contained in its twin and is removed before layout.");
    }

    // F2 — A read contained in another read is removed from the overlap graph (Myers 2005) and does
    // not become a spurious extra contig. greedy_scs on the substring-free set (k=3) → 1 contig.
    [Test]
    public void AssembleOLC_ContainedRead_DoesNotProduceExtraContig()
    {
        var reads = new List<string> { "AAAAACCCCC", "CCCCCGGGGG", "GGGGGTTTTT", "CCCGG" };

        var result = SequenceAssembler.AssembleOLC(reads, Exact(3));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "AAAAACCCCCGGGGGTTTTT" }),
            "CCCGG is a substring of CCCCCGGGGG; it is placed inside it, not emitted separately.");
    }

    // F3 — Consensus is the per-column majority vote over the layout (Langmead OLC p.28), not the
    // first read's bases. ACGTTGCTAC has an error (T) at column 7; the two other reads covering that
    // column carry A. Edges at l=5, identity 0.8: TGCAACGGAT→GCAACGGATT 9, ACGTTGCTAC→TGCAACGGAT 6
    // (5/6). Reference (edge-GREEDY + dumb_consensus 0.5): ACGTTGCAACGGATT.
    [Test]
    public void AssembleOLC_MismatchInOverlap_ResolvedByMajorityVote()
    {
        var reads = new List<string> { "ACGTTGCTAC", "TGCAACGGAT", "GCAACGGATT" };

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 0.8, MinContigLength: 1));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTGCAACGGATT" }),
            "Column 7 is T,A,A → majority A.");
    }

    // F3 — A 1:1 disagreement has no majority: dumb_consensus emits the ambiguity symbol ('N').
    [Test]
    public void AssembleOLC_TwoReadTieInOverlap_EmitsAmbiguitySymbol()
    {
        var reads = new List<string> { "ACGTTGCTAC", "TGCAACGGAT" };

        var result = SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 0.8, MinContigLength: 1));

        Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTGCNACGGAT" }),
            "Column 7 is T vs A (tie) → 'N' (Biopython dumb_consensus rule).");
    }

    // F4 — A zero-length "overlap" is no overlap (overlap graph edges need length ≥ l ≥ 1);
    // identity is a fraction in [0, 1].
    [TestCase(0, 0.9)]
    [TestCase(-1, 0.9)]
    [TestCase(5, -0.1)]
    [TestCase(5, 1.1)]
    [TestCase(5, double.NaN)]
    public void OverlapThresholds_OutOfRange_Throw(int minOverlap, double minIdentity)
    {
        var reads = new List<string> { "AAAAACCCCC", "CCCCCGGGGG" };

        Assert.Multiple(() =>
        {
            Assert.That(() => SequenceAssembler.FindOverlap(reads[0], reads[1], minOverlap, minIdentity),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.FindAllOverlaps(reads, minOverlap, minIdentity),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.FindAllOverlaps(reads, minOverlap, minIdentity, CancellationToken.None),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.AssembleOLC(reads,
                    new SequenceAssembler.AssemblyParameters(MinOverlap: minOverlap, MinIdentity: minIdentity)),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        });
    }

    // F5 — Omitting parameters must apply the documented defaults (MinOverlap 20, MinIdentity 0.9,
    // MinContigLength 100). Previously `new AssemblyParameters()` zero-initialised the record struct
    // (MinOverlap 0) and every pair was chained by a 0-length "overlap". Two 120-base reads sharing
    // only a 10-base suffix-prefix overlap (< 20) must therefore stay separate.
    [Test]
    public void AssembleOLC_NoParameters_UsesDocumentedDefaults()
    {
        string a = new string('A', 110) + "CGTCGTCGTC";
        string b = "CGTCGTCGTC" + new string('T', 110);

        var result = SequenceAssembler.AssembleOLC(new List<string> { a, b });

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.DefaultParameters,
                Is.EqualTo(new SequenceAssembler.AssemblyParameters(20, 0.9, 31, 100)));
            Assert.That(result.Contigs.OrderBy(c => c, StringComparer.Ordinal), Is.EqualTo(new[] { a, b }),
                "A 10-base overlap is below the default MinOverlap 20: no merge.");
        });
    }

    #endregion
}
