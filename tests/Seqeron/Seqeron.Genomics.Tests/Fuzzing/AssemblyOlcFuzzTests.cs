using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Assembly area — Overlap-Layout-Consensus assembly
/// (ASSEMBLY-OLC-001), the full OLC pipeline
/// <see cref="SequenceAssembler.AssembleOLC(IReadOnlyList{string}, SequenceAssembler.AssemblyParameters?)"/>
/// (overlap detection → greedy best-successor layout → chain-merge consensus).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to the unit and
/// asserts the code NEVER fails in an undisciplined way: no hang / infinite loop
/// in the greedy chain walk, no IndexOutOfRange / negative-length
/// <c>Substring(overlap)</c> in the chain merge (especially when
/// <c>MinOverlap</c> &gt; every read length so NO edge can exist), no false merge
/// of reads whose suffix/prefix overlap is below threshold, no dropped or
/// duplicated read, and no non-deterministic contig set across repeated runs.
/// Every input must resolve to a well-defined, theory-correct
/// <see cref="SequenceAssembler.AssemblyResult"/>: an empty result for null/empty
/// reads (no exception), and otherwise a set of contigs each of which is a
/// superstring of its constituent reads (INV-04) drawn from the input alphabet.
/// A raw runtime exception, a hang, a wrong-length contig, a phantom merge below
/// the overlap threshold, or an order-dependent / non-deterministic contig set is
/// a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly
/// Checklist: docs/checklists/03_FUZZING.md, row 145.
/// Algorithm doc: docs/algorithms/Assembly/Overlap_Layout_Consensus.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row (minOverlap&gt;read length, single read, no overlaps):
///          – minOverlap &gt; read length: the threshold <c>l</c> exceeds EVERY read
///            length, so <see cref="SequenceAssembler.FindOverlap"/>'s loop
///            (<c>overlapLen = min(|A|,|B|) .. l</c>) never runs → the overlap graph
///            is EDGELESS → each read is its own singleton contig (INV-05), with NO
///            IndexOutOfRange / negative <c>Substring</c> in the chain merge.
///          – single read: one read → one contig equal to that read (the trivial
///            superstring, INV-04), no crash, no self-overlap (INV-01).
///          – no overlaps: reads sharing no suffix/prefix match of length ≥ <c>l</c>
///            are each returned as a SEPARATE contig (INV-05) — no false merge below
///            the threshold, no crash, deterministic.
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Overlap_Layout_Consensus.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   INV-01 The overlap graph has no self-edge (a read is not overlapped with itself).
///   INV-02 Every reported overlap length is ≥ MinOverlap and ≤ min(|A|,|B|).
///   INV-03 The reported overlap for an ordered pair is the single longest qualifying
///          suffix-prefix match.
///   INV-04 An assembled contig is a superstring of its constituent reads; its length
///          is ≤ Σ read lengths and ≥ the longest single read.
///   INV-05 When the overlap graph is edgeless, each read becomes its own contig.
/// Null / empty <c>reads</c> → empty AssemblyResult, no exception (§3.3, §6.1).
/// Contigs shorter than <c>MinContigLength</c> are discarded (§3.2, §4.1) — the BE
/// tests below set MinContigLength small so the singleton boundary is observable.
///   SequenceAssembler.AssembleOLC(IReadOnlyList&lt;string&gt; reads,
///                                 AssemblyParameters? parameters) → AssemblyResult
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyOlcFuzzTests
{
    #region Helpers

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };
    // A wider alphabet to exercise the "alphabet-unconstrained, case-insensitive identity"
    // contract (§3.3): lowercase, N, gap/degenerate chars and non-biological symbols.
    private static readonly char[] WideAlphabet =
        { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'N', '-', '.', 'x', '7' };

    private static string RandomString(Random rng, int length, char[] alphabet)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Asserts an <see cref="SequenceAssembler.AssemblyResult"/> is WELL-FORMED per the
    /// documented contract, independent of how the layout chained the reads:
    ///   • every contig has length ≥ <paramref name="minContigLength"/> (the min-length
    ///     filter, §3.2);
    ///   • every contig is composed only of characters that occur in the input reads
    ///     (the assembler invents no symbols — it only concatenates read substrings);
    ///   • the reported statistics are self-consistent: TotalReads = read count,
    ///     TotalLength = Σ |contig|, LongestContig = max |contig| (or 0 when empty),
    ///     and N50 lies within the contig length range (INV/§3.2).
    /// This is the shared "valid result" oracle: it never re-implements the layout, only
    /// the universal post-conditions every OLC result must satisfy.
    /// </summary>
    private static void AssertWellFormed(
        SequenceAssembler.AssemblyResult result, IReadOnlyList<string> reads, int minContigLength)
    {
        result.TotalReads.Should().Be(reads.Count, "TotalReads echoes the input read count (§3.2)");

        var alphabet = new HashSet<char>();
        foreach (string read in reads)
            foreach (char c in read)
                alphabet.Add(c);

        foreach (string contig in result.Contigs)
        {
            contig.Length.Should().BeGreaterThanOrEqualTo(minContigLength,
                "every emitted contig survives the MinContigLength filter (§3.2)");
            foreach (char c in contig)
                alphabet.Should().Contain(c,
                    "the assembler only concatenates read substrings; it invents no character");
        }

        int sumLen = result.Contigs.Sum(c => c.Length);
        result.TotalLength.Should().Be(sumLen, "TotalLength = Σ contig lengths (§3.2)");

        if (result.Contigs.Count == 0)
        {
            result.LongestContig.Should().Be(0, "no contig ⇒ longest length 0");
        }
        else
        {
            int longest = result.Contigs.Max(c => c.Length);
            result.LongestContig.Should().Be(longest, "LongestContig = max contig length (§3.2)");
            result.N50.Should().BeInRange(result.Contigs.Min(c => c.Length), longest,
                "N50 is a contig length in the distribution (§3.2)");
        }
    }

    private static SequenceAssembler.AssemblyResult Run(
        IReadOnlyList<string> reads, int minOverlap, double minIdentity, int minContigLength)
        => SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(
                MinOverlap: minOverlap, MinIdentity: minIdentity, MinContigLength: minContigLength));

    #endregion

    #region ASSEMBLY-OLC-001 — OLC assembly (BE: minOverlap>read length, single read, no overlaps)

    #region Positive sanity — documented overlap detection, layout and consensus

    // Doc §7.1 worked example: three reads forming one unambiguous 5-overlap tiling collapse
    // into a single contig that is the superstring of all reads (INV-04).
    [Test]
    public void AssembleOLC_DocWorkedExample_UnambiguousChain_SingleSuperstringContig()
    {
        var reads = new List<string> { "AAAAACCCCC", "CCCCCGGGGG", "GGGGGTTTTT" };

        var result = Run(reads, minOverlap: 5, minIdentity: 1.0, minContigLength: 10);

        result.Contigs.Should().ContainSingle("an unambiguous overlap chain collapses to one contig (INV-04)");
        result.Contigs[0].Should().Be("AAAAACCCCCGGGGGTTTTT",
            "doc §7.1: merging A + B[overlap:] along the chain yields the 20-base superstring");
        AssertWellFormed(result, reads, 10);
    }

    // Reads with NO sufficient pairwise overlap stay separate — each is its own singleton contig
    // (INV-05); proves merges are driven by real overlaps, not produced spuriously.
    [Test]
    public void AssembleOLC_NoSufficientOverlap_ReadsStaySeparate()
    {
        var reads = new List<string> { "AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG" };

        var result = Run(reads, minOverlap: 5, minIdentity: 1.0, minContigLength: 5);

        result.Contigs.Should().HaveCount(3, "INV-05: no above-threshold overlap ⇒ three singleton contigs");
        result.Contigs.Should().BeEquivalentTo(reads, "each non-overlapping read passes through verbatim");
        AssertWellFormed(result, reads, 5);
    }

    // A contig is always a superstring of its reads, so it reconstructs the source genome the reads
    // were tiled from. Tile a known genome into overlapping reads and recover it.
    [Test]
    public void AssembleOLC_TiledGenome_ReconstructsSourceAsSuperstring()
    {
        const string genome = "ACGTACGTGGGGCCCCAATTTTGGCCAATTGGCCAA"; // 36 bases
        // 8 reads of length 12 stepping by 3 → consecutive reads share a 9-base suffix/prefix overlap.
        var reads = new List<string>();
        for (int start = 0; start + 12 <= genome.Length; start += 3)
            reads.Add(genome.Substring(start, 12));

        var result = Run(reads, minOverlap: 6, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().ContainSingle("the overlapping tiling forms one unambiguous chain");
        result.Contigs[0].Should().Be(genome, "the merged superstring of the chain is the source genome (INV-04)");
        AssertWellFormed(result, reads, 1);
    }

    #endregion

    #region BE — Boundary: minOverlap > read length (edgeless graph, no negative Substring, INV-05)

    // The threshold l exceeds EVERY read length, so FindOverlap's loop never runs → no edge →
    // each read is its own contig. The chain-merge must NOT attempt a negative/oversized
    // Substring(overlap). Reads kept long enough (and minContigLength small) so they survive.
    [Test]
    public void AssembleOLC_MinOverlapExceedsEveryReadLength_EachReadOwnContig_NoCrash()
    {
        var reads = new List<string> { "ACGTACGTAC", "ACGTACGTAC", "TTTTTTGGGG" }; // each len 10

        // minOverlap 50 ≫ 10 = every read length ⇒ impossible to satisfy ⇒ edgeless graph.
        Action act = () => Run(reads, minOverlap: 50, minIdentity: 1.0, minContigLength: 1);
        act.Should().NotThrow("minOverlap > read length ⇒ no overlap loop iteration, no Substring overflow");

        var result = Run(reads, minOverlap: 50, minIdentity: 1.0, minContigLength: 1);
        result.Contigs.Should().HaveCount(3, "INV-05: edgeless graph ⇒ every read is a singleton contig");
        result.Contigs.Should().BeEquivalentTo(reads, "no read is merged, dropped or altered");
        AssertWellFormed(result, reads, 1);
    }

    // Even reads that DO share a long identical suffix/prefix must NOT merge when minOverlap exceeds
    // their length — the threshold dominates; proves no false merge above the read length.
    [Test]
    public void AssembleOLC_MinOverlapAboveLength_IdenticalReadsStillNotMerged()
    {
        var reads = new List<string> { "ACGTACGT", "ACGTACGT" }; // identical, len 8

        var result = Run(reads, minOverlap: 9, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().HaveCount(2,
            "minOverlap (9) > read length (8) ⇒ no edge even between identical reads (INV-02/INV-05)");
        result.Contigs.Should().AllBe("ACGTACGT", "both reads pass through unchanged");
        AssertWellFormed(result, reads, 1);
    }

    // int.MaxValue minOverlap is the extreme boundary: still edgeless, still no overflow.
    [Test]
    public void AssembleOLC_MinOverlapIntMaxValue_EdgelessNoOverflow()
    {
        var reads = new List<string> { "AAAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAAA" };

        Action act = () => Run(reads, minOverlap: int.MaxValue, minIdentity: 1.0, minContigLength: 1);
        act.Should().NotThrow("a colossal minOverlap is simply unsatisfiable, not an overflow (BE: MaxInt)");

        var result = Run(reads, minOverlap: int.MaxValue, minIdentity: 1.0, minContigLength: 1);
        result.Contigs.Should().HaveCount(3, "INV-05: nothing can overlap ⇒ three singletons");
        AssertWellFormed(result, reads, 1);
    }

    // Fuzz: random reads, each strictly shorter than minOverlap ⇒ ALWAYS edgeless ⇒ contig set
    // equals the input reads (after the min-length filter), never throws, never merges.
    [Test]
    [CancelAfter(30_000)]
    public void AssembleOLC_RandomReadsBelowMinOverlap_AlwaysSingletons_NeverThrows()
    {
        var rng = new Random(145_001);
        for (int trial = 0; trial < 300; trial++)
        {
            int n = rng.Next(1, 8);
            var reads = new List<string>(n);
            int maxLen = 0;
            for (int i = 0; i < n; i++)
            {
                string r = RandomString(rng, rng.Next(1, 16), WideAlphabet);
                reads.Add(r);
                maxLen = Math.Max(maxLen, r.Length);
            }

            // minOverlap strictly above the longest read ⇒ no pair can reach the threshold.
            int minOverlap = maxLen + rng.Next(1, 50);

            SequenceAssembler.AssemblyResult result = default;
            Action act = () => result = Run(reads, minOverlap, minIdentity: 1.0, minContigLength: 1);
            act.Should().NotThrow($"minOverlap {minOverlap} > every read length ⇒ no crash (seed trial {trial})");

            result.Contigs.Should().BeEquivalentTo(reads,
                "edgeless graph ⇒ each read is returned verbatim as its own contig (INV-05)");
            AssertWellFormed(result, reads, 1);
        }
    }

    #endregion

    #region BE — Boundary: single read (one contig equal to that read, no self-overlap)

    // One read → one contig equal to that read (the trivial superstring, INV-04); no self-edge
    // (INV-01), no crash.
    [Test]
    public void AssembleOLC_SingleRead_OneContigEqualToRead()
    {
        var reads = new List<string> { "ACGTACGTACGTACGT" };

        var result = Run(reads, minOverlap: 4, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().ContainSingle("a single read assembles to a single contig (INV-01: no self-overlap)");
        result.Contigs[0].Should().Be("ACGTACGTACGTACGT", "the lone contig is exactly the read (trivial superstring, INV-04)");
        result.TotalReads.Should().Be(1);
        AssertWellFormed(result, reads, 1);
    }

    // A single read that is internally self-repetitive must NOT self-overlap into a longer/shorter
    // contig — INV-01 forbids a self-edge regardless of minOverlap.
    [Test]
    public void AssembleOLC_SingleRepetitiveRead_NoSelfOverlap()
    {
        var reads = new List<string> { "AAAAAAAAAAAAAAAAAAAA" }; // 20 A's, would "self-overlap" if a self-edge existed

        var result = Run(reads, minOverlap: 4, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().ContainSingle();
        result.Contigs[0].Should().Be("AAAAAAAAAAAAAAAAAAAA", "INV-01: a read is never overlapped with itself");
        AssertWellFormed(result, reads, 1);
    }

    // Fuzz: a single random read of any length / alphabet ALWAYS yields exactly that read as the
    // sole contig (when it survives the min-length filter); never throws.
    [Test]
    [CancelAfter(30_000)]
    public void AssembleOLC_SingleRandomRead_AlwaysReturnsItself()
    {
        var rng = new Random(145_002);
        for (int trial = 0; trial < 500; trial++)
        {
            string read = RandomString(rng, rng.Next(1, 60), WideAlphabet);
            var reads = new List<string> { read };
            int minOverlap = rng.Next(1, 40);

            SequenceAssembler.AssemblyResult result = default;
            Action act = () => result = Run(reads, minOverlap, minIdentity: 1.0, minContigLength: 1);
            act.Should().NotThrow($"a single read never crashes (seed trial {trial})");

            result.Contigs.Should().ContainSingle("one read ⇒ one contig (INV-01, INV-04)");
            result.Contigs[0].Should().Be(read, "the lone contig is exactly the read");
            AssertWellFormed(result, reads, 1);
        }
    }

    #endregion

    #region BE — Boundary: no overlaps (reads share no sufficient overlap, no false merge)

    // Reads with no qualifying suffix/prefix overlap stay separate; count and content preserved.
    [Test]
    public void AssembleOLC_PairwiseDisjointReads_NoFalseMerge()
    {
        // No suffix of any read equals a prefix of another at length ≥ 5.
        var reads = new List<string> { "AAAAACCCCC", "GGGGGTTTTT", "TTACCGGTAC" };

        var result = Run(reads, minOverlap: 5, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().HaveCount(3, "INV-05: no false merge below the overlap threshold");
        result.Contigs.Should().BeEquivalentTo(reads, "each read is returned as its own contig");
        AssertWellFormed(result, reads, 1);
    }

    // A near-miss overlap (one base below the threshold) must NOT trigger a merge — the threshold is
    // strict (INV-02: reported overlap ≥ minOverlap).
    [Test]
    public void AssembleOLC_OverlapOneBelowThreshold_NotMerged()
    {
        // "...CCCC" suffix == "CCCC..." prefix is a 4-base overlap; with minOverlap 5 it is rejected.
        var reads = new List<string> { "AAAAAACCCC", "CCCCGGGGGG" };

        var result = Run(reads, minOverlap: 5, minIdentity: 1.0, minContigLength: 1);

        result.Contigs.Should().HaveCount(2,
            "a 4-base overlap is below minOverlap 5 ⇒ no merge (INV-02 threshold is strict)");
        AssertWellFormed(result, reads, 1);

        // Lowering the threshold to 4 admits the same overlap ⇒ the two reads now merge into one
        // contig: proves the separation above was caused by the threshold, not a missing overlap.
        var merged = Run(reads, minOverlap: 4, minIdentity: 1.0, minContigLength: 1);
        merged.Contigs.Should().ContainSingle("at minOverlap 4 the 4-base overlap is accepted ⇒ one merged contig");
        merged.Contigs[0].Should().Be("AAAAAACCCCGGGGGG", "the shared 'CCCC' appears once in the superstring (INV-04)");
        AssertWellFormed(merged, reads, 1);
    }

    // Fuzz: many random reads with a HIGH threshold so overlaps are essentially impossible ⇒ the
    // contig set is exactly the surviving reads; no merge, no crash, no hang in the chain walk.
    [Test]
    [CancelAfter(60_000)]
    public void AssembleOLC_RandomReadsHighThreshold_NoFalseMerge_NeverHangs()
    {
        var rng = new Random(145_003);
        for (int trial = 0; trial < 200; trial++)
        {
            int n = rng.Next(1, 10);
            var reads = new List<string>(n);
            for (int i = 0; i < n; i++)
                reads.Add(RandomString(rng, rng.Next(1, 30), DnaAlphabet));

            // A threshold larger than the longest read makes every overlap impossible.
            int minOverlap = reads.Max(r => r.Length) + 1;

            SequenceAssembler.AssemblyResult result = default;
            Action act = () => result = Run(reads, minOverlap, minIdentity: 1.0, minContigLength: 1);
            act.Should().NotThrow($"no overlap is possible ⇒ no crash and the greedy walk terminates (trial {trial})");

            result.Contigs.Should().BeEquivalentTo(reads, "no false merge ⇒ contig set equals the input reads");
            AssertWellFormed(result, reads, 1);
        }
    }

    #endregion

    #region BE — Degenerate inputs and determinism

    // Null / empty reads → empty AssemblyResult, no exception (§3.3, §6.1).
    [Test]
    public void AssembleOLC_NullOrEmptyReads_EmptyResultNoThrow()
    {
        var fromNull = SequenceAssembler.AssembleOLC(null!);
        fromNull.Contigs.Should().BeEmpty("null reads ⇒ empty result (§6.1)");
        fromNull.TotalReads.Should().Be(0);
        fromNull.TotalLength.Should().Be(0);

        var fromEmpty = SequenceAssembler.AssembleOLC(new List<string>());
        fromEmpty.Contigs.Should().BeEmpty("empty reads ⇒ empty result (§6.1)");
        fromEmpty.TotalReads.Should().Be(0);
    }

    // Empty-string reads: they cannot overlap (length 0) and are filtered out by the min-length
    // filter; this must not crash (no negative Substring, no division surprise in stats).
    [Test]
    public void AssembleOLC_EmptyStringReads_NoCrash()
    {
        var reads = new List<string> { "", "", "" };

        Action act = () => Run(reads, minOverlap: 1, minIdentity: 1.0, minContigLength: 1);
        act.Should().NotThrow("empty reads are degenerate but must not crash the pipeline");

        var result = Run(reads, minOverlap: 1, minIdentity: 1.0, minContigLength: 1);
        result.Contigs.Should().BeEmpty("zero-length contigs are dropped by the MinContigLength filter");
        result.TotalReads.Should().Be(3, "TotalReads still echoes the input count");
        AssertWellFormed(result, reads, 1);
    }

    // Determinism: assembling the SAME read set repeatedly yields the IDENTICAL contig sequence
    // (same contigs, same order) — no non-deterministic layout/ordering. Covers all three BE
    // boundaries (single, no-overlap, above-threshold) plus a real merge.
    [Test]
    [CancelAfter(30_000)]
    public void AssembleOLC_RepeatedRuns_DeterministicContigSet()
    {
        var rng = new Random(145_004);
        for (int trial = 0; trial < 100; trial++)
        {
            int n = rng.Next(1, 8);
            var reads = new List<string>(n);
            for (int i = 0; i < n; i++)
                reads.Add(RandomString(rng, rng.Next(1, 25), DnaAlphabet));

            int minOverlap = rng.Next(1, 12);

            var first = Run(reads, minOverlap, minIdentity: 1.0, minContigLength: 1);
            for (int repeat = 0; repeat < 5; repeat++)
            {
                var again = Run(reads, minOverlap, minIdentity: 1.0, minContigLength: 1);
                again.Contigs.Should().Equal(first.Contigs,
                    "AssembleOLC is a deterministic pure function of (reads, parameters)");
            }
            AssertWellFormed(first, reads, 1);
        }
    }

    // Broad fuzz across the whole parameter space (low/high minOverlap, exact and approximate
    // identity, mixed read lengths incl. empties) — the pipeline NEVER throws, NEVER hangs, and
    // ALWAYS returns a well-formed result whose contigs are read-derived superstrings.
    [Test]
    [CancelAfter(120_000)]
    public void AssembleOLC_BroadFuzz_NeverThrows_WellFormed()
    {
        var rng = new Random(145_005);
        for (int trial = 0; trial < 400; trial++)
        {
            int n = rng.Next(0, 12);
            var reads = new List<string>(n);
            for (int i = 0; i < n; i++)
                reads.Add(RandomString(rng, rng.Next(0, 30), WideAlphabet));

            int minOverlap = rng.Next(6) switch
            {
                0 => 1,                                   // permissive
                1 => rng.Next(1, 10),                     // mid
                2 => reads.Count == 0 ? 5 : reads.Max(r => r.Length) + 1, // above every read
                3 => int.MaxValue,                        // BE: MaxInt
                4 => rng.Next(1, 40),                     // arbitrary
                _ => 20,                                  // default
            };
            double minIdentity = rng.Next(3) switch { 0 => 1.0, 1 => 0.9, _ => 0.8 };
            int minContigLength = rng.Next(3); // 0, 1, 2

            SequenceAssembler.AssemblyResult result = default;
            Action act = () => result = Run(reads, minOverlap, minIdentity, minContigLength);
            act.Should().NotThrow($"AssembleOLC never throws on fuzzed input (trial {trial})");

            AssertWellFormed(result, reads, minContigLength);
        }
    }

    #endregion

    #endregion
}
