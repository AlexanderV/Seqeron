// ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly
// Evidence: docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-OLC-001.md
// Source: Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987-991, DOI 10.1038/nbt.2023;
//         Langmead B., "Overlap Layout Consensus assembly" / "Assembly & Shortest Common
//         Superstring" (JHU lecture notes).

using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Seqeron.Genomics.Tests;

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
}
