using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Analysis area — Open Reading Frame detection (GENOMIC-ORF-001), the six-frame
/// canonical ORF enumerator
/// <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)"/>: for a DNA sequence it
/// enumerates every span that begins at a start codon (ATG) and ends at the first in-frame stop
/// codon (TAA, TAG, TGA) with no internal in-frame stop, scanned in all six reading frames (three on
/// the forward strand, three on the reverse complement), filtered by a minimum nucleotide length.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary sequences to the unit and asserts the code NEVER fails in an
/// undisciplined way: no crash / IndexOutOfRange when a frame ends mid-codon, when an ATG is opened
/// near the strand end with NO downstream in-frame stop (the canonical "read past the end" trap — the
/// inner codon scan must stop at <c>seq.Length - 3</c>, never index past it), when the sequence is
/// empty or shorter than one codon, or when minLength is a boundary value (0, 1, negative, int.Max
/// large). It also asserts the code NEVER produces nonsense: never an ORF whose Sequence does not
/// start with ATG, never one whose Sequence does not end in TAA/TAG/TGA, never a Length not divisible
/// by 3, never a Length below minLength, never an out-of-bounds Position, never a Frame outside {1,2,3},
/// never a fabricated ORF where the strand does not spell one, never a MISSED or DUPLICATED nested ORF.
/// The only declared validation exception is a null sequence → ArgumentNullException (§3.3, §6.1); an
/// empty subject and a no-ATG / no-stop subject are all accepted and yield an empty (no-crash) result.
/// A raw runtime exception, a hang, a malformed ORF, a wrong coordinate/frame, a missed or duplicated
/// nested ORF, or a non-deterministic output is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-ORF-001 — six-frame ORF detection (ATG → first in-frame stop, ATG-only / standard code)
/// Checklist: docs/checklists/03_FUZZING.md, row 177 (BE — "no ATG, no stop, nested ORFs").
/// Algorithm doc: docs/algorithms/Analysis/Open_Reading_Frame_Detection.md
/// Distinct from row 28 (ANNOT-ORF-001, GenomeAnnotator.FindOrfs — the genetic-code-parameterized
/// Annotation-area ORF finder): THIS unit is the fixed ATG-only / standard-genetic-code six-frame
/// scan on GenomicAnalyzer (no alternative start codons, no non-standard genetic codes — §5.3, §6.2).
///
/// Fuzz strategy exercised for THIS unit (BE = Boundary Exploitation — граничні значення: 0, -1,
/// MaxInt, empty — docs/checklists/03_FUZZING.md §Description), mapped to the row's three targets:
///   • NO ATG — a subject containing NO start codon in any frame on either strand (e.g. an all-stop
///     or A/C/T-only sequence with no "ATG" substring) → EMPTY result, NO crash, NO fabricated ORF
///     (§6.1 "no ATG → empty"). The empty-string and sub-codon subjects are the EMPTY boundary.
///   • NO STOP — an ATG opened with NO downstream in-frame stop before the strand end → that reading
///     is NOT a complete ORF and is NOT reported; the inner codon scan must NOT read a codon past
///     <c>seq.Length - 3</c> (the IndexOutOfRange trap). The opposite strand may still produce ORFs;
///     each must remain well-formed (§6.1 "ATG with no in-frame stop → not reported").
///   • NESTED ORFs — two (or more) in-frame ATGs sharing the SAME downstream stop (outer ATG..stop
///     containing an inner in-frame ATG): ALL are reported, each ending at that shared stop, with
///     correct 0-based Position and decreasing length (canonical Rosalind semantics — §6.1 "nested
///     ATGs sharing one stop → all reported"; worked example §7.1).
///
/// Note on Malformed Content / Injection: the subject is a <see cref="DnaSequence"/> (uppercased and
/// validated to {A,C,G,T} at construction, so out-of-domain residues / null bytes / unicode cannot
/// reach the scan); thus this is a pure boundary (BE) row over the subject SHAPE (empty / sub-codon /
/// no-ATG / no-stop / nested) and the minLength filter boundary, exactly as the row says.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Open_Reading_Frame_Detection.md §2.4 INV-01..05, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   INV-01 every reported ORF Sequence begins with ATG;
///   INV-02 every reported ORF Sequence ends with TAA/TAG/TGA (the first in-frame stop), no internal stop;
///   INV-03 Length is divisible by 3 (whole codons, start through stop inclusive);
///   INV-04 Length ≥ minLength (length filter applied before yielding);
///   INV-05 Frame ∈ {1,2,3}; IsReverseComplement selects the strand; Position is the 0-based start
///          offset within the scanned strand (forward = original; reverse = reverse complement).
///   Null sequence → ArgumentNullException. No ATG / no in-frame stop / sub-codon / empty → no ORF.
///   GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int) → IEnumerable&lt;OrfInfo&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicOrfFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };
    private const int CodonLength = 3;
    private static readonly string[] StopCodons = { "TAA", "TAG", "TGA" };

    #region Helpers

    /// <summary>A random ACGT string of the given length (length 0 ⇒ empty string).</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    private static bool IsStop(string codon) => StopCodons.Contains(codon);

    /// <summary>
    /// Independent oracle for the single-strand ORF spans, built directly from the documented model
    /// (§2.2, §4.1) — NOT from the unit. For each frame offset f∈{0,1,2}, scan codons; at each ATG
    /// find the FIRST in-frame stop downstream and emit the ATG→stop span (start..i+3, inclusive of
    /// the stop) if its length ≥ minLength. Every ATG is considered independently, so nested ATGs
    /// sharing a stop are all emitted. An ATG with no downstream in-frame stop emits nothing. The
    /// inner scan never reads past <c>seq.Length - 3</c>. Returns (position, frame 1..3, sequence).
    /// </summary>
    private static List<(int Position, int Frame, string Seq)> OracleStrand(string seq, int minLength)
    {
        var orfs = new List<(int, int, string)>();
        for (int frame = 0; frame < CodonLength; frame++)
        {
            for (int start = frame; start <= seq.Length - CodonLength; start += CodonLength)
            {
                if (string.CompareOrdinal(seq, start, "ATG", 0, CodonLength) != 0)
                    continue;
                for (int i = start; i <= seq.Length - CodonLength; i += CodonLength)
                {
                    if (!IsStop(seq.Substring(i, CodonLength)))
                        continue;
                    int length = i + CodonLength - start;
                    if (length >= minLength)
                        orfs.Add((start, frame + 1, seq.Substring(start, length)));
                    break; // first in-frame stop terminates this reading
                }
            }
        }
        return orfs;
    }

    /// <summary>
    /// Asserts a single ORF is WELL-FORMED against the documented invariants and the strand it was
    /// scanned from: starts with ATG (INV-01), ends with a stop codon (INV-02), Length % 3 == 0 and
    /// Length == Sequence.Length (INV-03), Length ≥ minLength (INV-04), Frame ∈ {1,2,3} (INV-05),
    /// Position is an in-bounds 0-based start in the correct frame, the strand SUBSTRING at Position
    /// exactly equals the reported Sequence (no fabricated/shifted span), and there is NO internal
    /// in-frame stop strictly before the terminal stop (first-in-frame-stop rule, INV-02).
    /// </summary>
    private static void AssertOrfWellFormed(OrfInfo orf, string forward, int minLength)
    {
        string strand = orf.IsReverseComplement
            ? new DnaSequence(forward).ReverseComplement().Sequence
            : forward;

        orf.Frame.Should().BeInRange(1, 3, "Frame ∈ {1,2,3} (INV-05)");
        orf.Position.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based offset (INV-05)");
        (orf.Position % CodonLength).Should().Be(orf.Frame - 1,
            "Position lies on the reported frame's codon grid (INV-05)");

        orf.Sequence.Should().NotBeNull();
        orf.Length.Should().Be(orf.Sequence.Length, "Length is the span length (INV-03)");
        (orf.Length % CodonLength).Should().Be(0, "Length is divisible by 3 (INV-03)");
        orf.Length.Should().BeGreaterThanOrEqualTo(minLength, "Length ≥ minLength (INV-04)");
        orf.CodonCount.Should().Be(orf.Length / CodonLength, "CodonCount is Length/3");

        (orf.Position + orf.Length).Should().BeLessThanOrEqualTo(strand.Length,
            "the span is in-bounds on the scanned strand (no read past the end)");
        strand.Substring(orf.Position, orf.Length).Should().Be(orf.Sequence,
            "the scanned strand spells exactly the reported ORF at Position (no fabricated/shifted span)");

        orf.Sequence.Substring(0, CodonLength).Should().Be("ATG", "ORF starts at ATG (INV-01)");
        IsStop(orf.Sequence.Substring(orf.Length - CodonLength, CodonLength)).Should()
            .BeTrue("ORF ends at an in-frame stop TAA/TAG/TGA (INV-02)");
        for (int i = 0; i < orf.Length - CodonLength; i += CodonLength)
            IsStop(orf.Sequence.Substring(i, CodonLength)).Should()
                .BeFalse("the terminal stop is the FIRST in-frame stop (no internal stop, INV-02)");
    }

    /// <summary>
    /// Cross-checks the full six-frame result against the independent oracle (forward strand + reverse
    /// complement strand), and asserts every yielded ORF is well-formed. The reported (Position, Frame,
    /// IsReverseComplement, Sequence) tuple set must EXACTLY equal the oracle's — no missed ORF, no
    /// fabricated ORF, no duplicate. Idempotent: two enumerations must agree (determinism).
    /// </summary>
    private static List<OrfInfo> AssertWellFormed(string forward, int minLength)
    {
        var dna = new DnaSequence(forward);
        var result = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength).ToList();

        foreach (var orf in result)
            AssertOrfWellFormed(orf, forward, minLength);

        var actual = result
            .Select(o => (o.Position, o.Frame, o.IsReverseComplement, o.Sequence))
            .ToHashSet();
        actual.Should().HaveCount(result.Count, "no duplicate ORF is yielded");

        string revComp = dna.ReverseComplement().Sequence;
        var expected = OracleStrand(dna.Sequence, minLength)
            .Select(o => (o.Position, o.Frame, false, o.Seq))
            .Concat(OracleStrand(revComp, minLength).Select(o => (o.Position, o.Frame, true, o.Seq)))
            .ToHashSet();

        actual.Should().BeEquivalentTo(expected,
            "the six-frame ORF set must exactly match the spec oracle (none missed, none fabricated)");

        // Determinism: re-enumerate and compare.
        var again = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength)
            .Select(o => (o.Position, o.Frame, o.IsReverseComplement, o.Sequence)).ToHashSet();
        again.Should().BeEquivalentTo(actual, "enumeration is deterministic");

        return result;
    }

    #endregion

    #region GENOMIC-ORF-001 — Open Reading Frame detection (BE: no ATG, no stop, nested)

    // ─── Positive sanity: a known ORF is found at the documented coordinates/frame ──────────────

    [Test]
    public void PositiveSanity_KnownForwardOrf_FoundAtDocumentedCoordinates()
    {
        // Worked example (doc §7.1): "ATGAAAAAATAA" with minLength 1 ⇒ one forward ORF, frame 1,
        // position 0, the whole span "ATGAAAAAATAA" (protein candidate "MKK").
        var dna = new DnaSequence("ATGAAAAAATAA");
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();

        var fwd = orfs.Where(o => !o.IsReverseComplement).ToList();
        fwd.Should().ContainSingle("the only forward ORF is the documented ATG..TAA span");
        fwd[0].Sequence.Should().Be("ATGAAAAAATAA");
        fwd[0].Position.Should().Be(0);
        fwd[0].Frame.Should().Be(1);
        fwd[0].IsReverseComplement.Should().BeFalse();
        fwd[0].Length.Should().Be(12);
        fwd[0].CodonCount.Should().Be(4);

        AssertWellFormed("ATGAAAAAATAA", minLength: 1);
    }

    [Test]
    public void NullSequence_Throws()
    {
        Action act = () => GenomicAnalyzer.FindOpenReadingFrames(null!, minLength: 1).ToList();
        act.Should().Throw<ArgumentNullException>("null sequence is the only declared validation error (§3.3)");
    }

    // ─── NO ATG: subjects with no start codon in any frame on either strand ─────────────────────

    [Test]
    public void NoAtg_Empty_NoCrash()
    {
        AssertWellFormed(string.Empty, minLength: 1).Should().BeEmpty("empty subject yields no ORF");
        // Empty + default minLength (100) likewise.
        GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(string.Empty)).Should().BeEmpty();
    }

    [Test]
    public void NoAtg_SubCodon_NoCrash()
    {
        foreach (var s in new[] { "A", "AT", "TG", "G" })
            AssertWellFormed(s, minLength: 1).Should().BeEmpty($"'{s}' is shorter than/cannot host an ORF");
    }

    [Test]
    public void NoAtg_AllStopCodons_NoOrf()
    {
        // No "ATG" substring on EITHER strand: only stop codons; reverse complement of stops also lacks ATG.
        var orfs = AssertWellFormed("TAATAGTGATAA", minLength: 1);
        orfs.Should().BeEmpty("no ATG anywhere ⇒ no ORF, none fabricated (§6.1)");
    }

    [Test]
    public void NoAtg_Random_NoStartCodonAnywhere_NoOrf()
    {
        var rng = new Random(177_001);
        int hits = 0;
        for (int t = 0; t < 400; t++)
        {
            // Build an A/C/T-only string: "ATG" requires a 'G', so removing G from the alphabet
            // guarantees no ATG on the forward strand. The reverse complement of an {A,C,T} string is
            // an {A,G,T} string (no C) — which CAN contain ATG — so we only assert on the forward strand.
            var sb = new StringBuilder();
            foreach (char c in RandomDna(rng, rng.Next(0, 40)))
                sb.Append(c == 'G' ? 'A' : c);
            string fwd = sb.ToString();

            var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), minLength: 1).ToList();
            orfs.Where(o => !o.IsReverseComplement).Should()
                .BeEmpty("no ATG on a G-free forward strand ⇒ no forward ORF");
            if (orfs.Any(o => o.IsReverseComplement)) hits++;
        }
        hits.Should().BeGreaterThan(0, "reverse-complement ORFs should still be reachable (sanity)");
    }

    // ─── NO STOP: an ATG opened with no downstream in-frame stop (read-past-end trap) ───────────

    [Test]
    public void NoStop_AtgWithoutInFrameStop_NotReported_NoIndexOutOfRange()
    {
        // ATG followed by codons that are never a stop, running to the very end with no in-frame stop.
        // "ATG" + "AAA"*k leaves the reading open ⇒ not a complete ORF ⇒ not reported (§6.1).
        foreach (int k in new[] { 0, 1, 3, 10 })
        {
            string fwd = "ATG" + string.Concat(Enumerable.Repeat("AAA", k));
            var orfs = AssertWellFormed(fwd, minLength: 1);
            orfs.Where(o => !o.IsReverseComplement).Should()
                .BeEmpty("an ATG with no downstream in-frame stop is incomplete, not reported (§6.1)");
        }
    }

    [Test]
    public void NoStop_AtgAtVeryEnd_FrameEndsMidCodon_NoCrash()
    {
        // ATG flush at the strand end, then a trailing partial codon: the inner scan must stop at
        // seq.Length-3 and never read the partial tail (IndexOutOfRange trap).
        foreach (var fwd in new[] { "AAATG", "CCCATGA", "GGGATGAA", "ATGA", "ATGAA" })
        {
            Action act = () => AssertWellFormed(fwd, minLength: 1);
            act.Should().NotThrow($"'{fwd}' must not read a codon past the end");
        }
    }

    [Test]
    public void NoStop_RandomNoStopOnForward_NoForwardOrf()
    {
        var rng = new Random(177_002);
        for (int t = 0; t < 400; t++)
        {
            // Forward strand of {A,C,G} only: A/C/G codons can spell ATG (start) but NEVER a stop
            // (every stop codon TAA/TAG/TGA contains a T as a non-first base... actually TAA/TAG/TGA
            // all contain T; an alphabet without T cannot spell any stop). So no forward ORF can close.
            var sb = new StringBuilder();
            foreach (char c in RandomDna(rng, rng.Next(0, 45)))
                sb.Append(c == 'T' ? 'C' : c);
            string fwd = sb.ToString();

            var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), minLength: 1).ToList();
            orfs.Where(o => !o.IsReverseComplement).Should()
                .BeEmpty("no stop codon possible on a T-free strand ⇒ no closed forward ORF");
        }
    }

    // ─── NESTED ORFs: in-frame ATGs sharing a stop ⇒ all reported ───────────────────────────────

    [Test]
    public void Nested_TwoInFrameAtgsShareStop_BothReported()
    {
        // Outer ATG at 0, inner in-frame ATG at 6 (codon boundary), shared stop "TAA" at the end.
        // Layout codons: ATG AAA ATG AAA TAA  (positions 0..14). Both ATGs are in frame 1 and reach
        // the SAME stop ⇒ BOTH reported (canonical Rosalind nesting, §6.1).
        string fwd = "ATGAAAATGAAATAA"; // 15 nt, 5 codons
        var orfs = AssertWellFormed(fwd, minLength: 1);
        var fwdF1 = orfs.Where(o => !o.IsReverseComplement && o.Frame == 1)
            .OrderBy(o => o.Position).ToList();

        fwdF1.Should().HaveCount(2, "both nested in-frame ATGs sharing the stop are reported");
        fwdF1[0].Position.Should().Be(0);
        fwdF1[0].Sequence.Should().Be("ATGAAAATGAAATAA");
        fwdF1[1].Position.Should().Be(6);
        fwdF1[1].Sequence.Should().Be("ATGAAATAA");
        // Outer is strictly longer than the inner; both end at the same stop.
        fwdF1[0].Length.Should().BeGreaterThan(fwdF1[1].Length);
        fwdF1[0].Sequence.Should().EndWith("TAA");
        fwdF1[1].Sequence.Should().EndWith("TAA");
    }

    [Test]
    public void Nested_ThreeAtgsShareStop_AllThreeReported_DecreasingLength()
    {
        // Three in-frame ATGs (0, 3, 6) sharing the stop at the end ⇒ three ORFs of decreasing length.
        string fwd = "ATGATGATGTAA"; // ATG ATG ATG TAA
        var orfs = AssertWellFormed(fwd, minLength: 1);
        var f1 = orfs.Where(o => !o.IsReverseComplement && o.Frame == 1)
            .OrderBy(o => o.Position).ToList();

        f1.Select(o => o.Position).Should().Equal(0, 3, 6);
        f1.Select(o => o.Sequence).Should().Equal("ATGATGATGTAA", "ATGATGTAA", "ATGTAA");
        f1.Select(o => o.Length).Should().BeInDescendingOrder();
    }

    [Test]
    public void Nested_MinLengthFiltersInnerButKeepsOuter()
    {
        // Same nested layout; raising minLength above the inner ORF's length drops only the inner one.
        // Inner ORF "ATGTAA" has length 6; outer "ATGATGATGTAA" has length 12. minLength 7 keeps outer
        // two (12, 9) and the inner (6) is filtered out.
        string fwd = "ATGATGATGTAA";
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), minLength: 7).ToList();
        foreach (var o in orfs) AssertOrfWellFormed(o, fwd, minLength: 7);

        var f1 = orfs.Where(o => !o.IsReverseComplement && o.Frame == 1)
            .Select(o => o.Length).OrderByDescending(x => x).ToList();
        f1.Should().Equal(new[] { 12, 9 }, "the 6-nt inner ORF is filtered by minLength 7 (INV-04)");
    }

    [Test]
    public void InternalStop_TerminatesAtFirstStop_NotSecond()
    {
        // ATG ... TAA ... TAG: the first in-frame stop (TAA) terminates the ORF; the second stop is
        // NOT part of the span (off-by-one / first-stop boundary check).
        string fwd = "ATGAAATAAGGGTAG"; // ATG AAA TAA GGG TAG
        var orfs = AssertWellFormed(fwd, minLength: 1);
        var f1 = orfs.Where(o => !o.IsReverseComplement && o.Frame == 1).ToList();
        f1.Should().ContainSingle();
        f1[0].Sequence.Should().Be("ATGAAATAA", "ORF ends at the FIRST in-frame stop, not the later one");
    }

    // ─── minLength boundary values (BE: 0, 1, negative, large) ──────────────────────────────────

    [Test]
    public void MinLength_BoundaryValues_NoCrash_RespectFilter()
    {
        string fwd = "ATGAAATAA"; // one forward ORF, length 9
        foreach (int min in new[] { int.MinValue, -1, 0, 1, 9 })
        {
            var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), min).ToList();
            foreach (var o in orfs) AssertOrfWellFormed(o, fwd, min);
            orfs.Where(o => !o.IsReverseComplement && o.Frame == 1).Should()
                .ContainSingle("the length-9 ORF passes any minLength ≤ 9");
        }

        // minLength larger than every span ⇒ all filtered out, no crash.
        GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), 10).ToList()
            .Should().BeEmpty("minLength 10 exceeds the only ORF's length 9 (INV-04)");
        GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(fwd), int.MaxValue).ToList()
            .Should().BeEmpty("int.MaxValue minLength filters everything, no crash");
    }

    // ─── Broad random fuzz: well-formedness + oracle equivalence over many shapes ────────────────

    [Test]
    [CancelAfter(60_000)]
    public void RandomSequences_WellFormed_AndMatchOracle()
    {
        var rng = new Random(177_003);
        for (int t = 0; t < 1500; t++)
        {
            string fwd = RandomDna(rng, rng.Next(0, 90));
            int bucket = t % 5;
            int min;
            if (bucket == 0) min = 1;
            else if (bucket == 1) min = 3;
            else if (bucket == 2) min = 6;
            else if (bucket == 3) min = 100; // the production default
            else min = rng.Next(-5, 30);
            AssertWellFormed(fwd, min);
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void RandomSeededOrfs_AlwaysFound_AndWellFormed()
    {
        // Embed a guaranteed in-frame ORF (ATG..TAA) at a random codon-aligned offset so the result is
        // non-empty often, then verify full well-formedness and oracle equivalence.
        var rng = new Random(177_004);
        for (int t = 0; t < 800; t++)
        {
            int prefixCodons = rng.Next(0, 6);
            int bodyCodons = rng.Next(1, 8);
            var sb = new StringBuilder();
            for (int i = 0; i < prefixCodons; i++) sb.Append("CCC"); // never ATG, never stop
            sb.Append("ATG");
            for (int i = 0; i < bodyCodons; i++) sb.Append("AAA");    // never stop
            sb.Append(StopCodons[rng.Next(StopCodons.Length)]);
            string fwd = sb.ToString();

            var orfs = AssertWellFormed(fwd, minLength: 1);
            orfs.Where(o => !o.IsReverseComplement && o.Frame == 1 && o.Position == prefixCodons * 3)
                .Should().ContainSingle("the embedded forward ORF is found at its codon-aligned offset");
        }
    }

    #endregion
}
