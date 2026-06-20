using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Annotation area — open reading frame (ORF) detection
/// (ANNOT-ORF-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, malformed and boundary inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// out-of-bounds indexing (IndexOutOfRangeException / ArgumentOutOfRangeException
/// leaking from internal Substring), no DivideByZero, no NullReferenceException,
/// and no nonsense output (an ORF whose span is not codon-aligned, or that runs off
/// the end of the sequence, or that is reported when the biological contract says it
/// must not be). Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a documented, intentional validation outcome. A raw runtime exception,
/// a hang, or a blow-up of bogus ORFs on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-ORF-001 — ORF finding
/// Checklist: docs/checklists/03_FUZZING.md, row 28.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a sequence with NO start codon, an all-start-no-stop
///          sequence, and minLength = 0.
///   • MC = Malformed Content — non-DNA characters (N, digits, IUPAC ambiguity
///          codes, garbage) embedded in the scanned sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ORF-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// An ORF is a contiguous stretch of codons beginning with a start codon and ending
/// with an in-frame stop codon, with no intervening in-frame stop. The canonical
/// Annotation-area detector is
///   GenomeAnnotator.FindOrfs(string dnaSequence,
///                            int minLength = 100,
///                            bool searchBothStrands = true,
///                            bool requireStartCodon = true)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 91–120).
/// This is the surface the checklist row probes (it is the canonical annotation ORF
/// finder named in ORF_Detection.md §1, §5.1). NOTE there is a distinct Analysis-area
/// surface GenomicAnalyzer.FindOpenReadingFrames (GENOMIC-ORF-001, row 177) with
/// different semantics (ATG-only start, minLength compared to nucleotide span); it is
/// NOT the unit under fuzz here.
///
/// Fixed codon sets (ORF_Detection.md §4.2; GenomeAnnotator.cs StartCodons/StopCodons):
///   • Start codons: ATG, GTG, TTG (common prokaryotic starts).
///   • Stop codons : TAA, TAG, TGA (standard).
/// The scan walks each of the three forward frames (plus three reverse-complement
/// frames when searchBothStrands = true) in codon-sized steps; it accumulates pending
/// start positions and, on reaching an in-frame stop, emits every pending start whose
/// stop-delimited amino-acid length L_aa = (t − s)/3 is ≥ minLength
/// (GenomeAnnotator.cs lines 139–213; ORF_Detection.md §2.2, §4.1).
///
/// Documented parameter / validation contract (ORF_Detection.md §3.3, §6.1):
///   • dnaSequence null or empty → NO ORFs (the method `yield break`s on
///     string.IsNullOrEmpty; GenomeAnnotator.cs lines 97–98). NOT an exception — this
///     is the documented contract for the canonical surface, so the fuzz tests pin the
///     empty/null result, NOT a throw.
///   • The canonical finder performs NO numeric validation of minLength and NO
///     alphabet pre-validation. minLength is simply compared with `>=` against the
///     amino-acid length, and each codon is upper-cased then tested for exact
///     membership in the fixed A/C/G/T codon sets. The fuzz tests below therefore
///     assert the THEORY-CORRECT consequences of that contract, not an exception:
///       – No start codon with requireStartCodon = true → no pending start is ever
///         recorded → empty result (ORF_Detection.md §6.1 row 3).
///       – ALL-START-NO-STOP (e.g. "ATGATGATG…" with no in-frame stop) with the
///         DEFAULT requireStartCodon = true → NOT reported: a terminating stop is
///         REQUIRED before emission, so a start that never reaches a stop yields
///         nothing (ORF_Detection.md §6.1 row 4, INV-02). This is THE key biological
///         boundary for this unit — a partial/unterminated ORF is suppressed by
///         default. Only when requireStartCodon = false does the finder additionally
///         emit trailing frame segments that run to the (codon-aligned) sequence end
///         (GenomeAnnotator.cs lines 183–212). The fuzz tests pin BOTH halves of that
///         contract so the default-suppression cannot silently drift.
///       – minLength = 0 is degenerate: `orfAaLength >= 0` is ALWAYS true, so a
///         zero-amino-acid ORF (a start codon immediately followed by an in-frame
///         stop, e.g. "ATGTAA", with L_aa = 0) IS reported. The canonical contract
///         does not reject minLength = 0 (no validation gate, ORF_Detection.md §3.3);
///         the disciplined outcome here is that the degenerate threshold produces a
///         WELL-FORMED trivial ORF (codon-aligned, in-bounds, correct start/stop),
///         not a crash and not a blow-up of garbage. The fuzz test pins that
///         well-formedness rather than asserting an exception the contract never
///         promises.
///       – Non-DNA characters (N, digits, IUPAC ambiguity codes, punctuation) are
///         NOT pre-validated and NOT rejected: a codon containing such a character
///         simply fails exact start/stop membership and is treated as an ordinary
///         non-boundary codon (ORF_Detection.md §3.3, §6.1 rows 6–7). The fuzz tests
///         pin that malformed content never crashes (no out-of-range Substring on the
///         3-char codon window), never invents a spurious start/stop from a non-ACGT
///         triplet, and that an ORF surrounding internal garbage can still close at a
///         genuine downstream stop.
///
/// Documented invariants pinned on every positive result (ORF_Detection.md §2.4):
///   INV-01 every reported ORF begins with ATG/GTG/TTG (requireStartCodon = true);
///   INV-02 every reported ORF ends with TAA/TAG/TGA (requireStartCodon = true);
///   INV-03 every reported ORF span (End − Start) is divisible by 3;
///   INV-04 coordinates satisfy 0 ≤ Start < End ≤ dnaSequence.Length.
///
/// FindOrfs is a LAZY iterator (`yield`); even the null/empty short-circuit sits
/// inside the iterator body, so it only fires on enumeration. Every test therefore
/// forces enumeration (`.ToList()`) so the documented behavior actually surfaces and
/// any hang would manifest as a non-terminating materialization. Deterministic only:
/// any random fuzz input is generated from a locally fixed `new Random(seed)`.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AnnotationFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-ORF-001 — ORF finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-ORF-001 — ORF finding

    #region BE — Boundary: no start codon (no ATG/GTG/TTG) in the sequence

    /// <summary>
    /// BE: a sequence containing NO start codon at all. With the default
    /// requireStartCodon = true, no pending ORF start is ever recorded, so even a
    /// sequence riddled with stop codons emits NOTHING — no crash, no spurious ORF
    /// (ORF_Detection.md §6.1 row 3, INV-01). We use "CCC...TAA..." built entirely from
    /// codons that are neither starts nor (for the body) stops, plus a stop, and scan
    /// both strands so a start cannot sneak in on the reverse complement either. The
    /// body codon "CCC" reverse-complements to "GGG", still not a start; the stop "TAA"
    /// reverse-complements within frame to "TTA"/"TAA" variants — none are starts. The
    /// disciplined outcome is the empty result.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NoStartCodon_ReturnsEmptyAndDoesNotThrow()
    {
        // No ATG/GTG/TTG anywhere; only C/A bodies and a TAA stop.
        const string seq = "CCCCCCCCCCCCCCCCCCTAA";

        var act = () => GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true).ToList();
        act.Should().NotThrow("a sequence with no start codon never records a pending ORF; it simply yields nothing");

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true).ToList();
        orfs.Should().BeEmpty(
            "with requireStartCodon = true and no ATG/GTG/TTG present, no ORF can be emitted on either strand");
    }

    #endregion

    #region BE — Boundary: all start, no stop (the key biological boundary)

    /// <summary>
    /// BE: "ATGATGATG…" — a start codon tiled repeatedly with NO in-frame stop. This is
    /// THE key biological boundary for ORF finding: a start that never reaches a stop is
    /// an UNTERMINATED ORF. Under the DEFAULT requireStartCodon = true the canonical
    /// finder REQUIRES a terminating stop before emission, so the unterminated frame is
    /// suppressed — empty result, no crash, no off-the-end Substring
    /// (ORF_Detection.md §6.1 row 4, INV-02; GenomeAnnotator.cs lines 158–187 emit only
    /// on a stop). We restrict to the forward strand so a reverse-complement stop cannot
    /// terminate it: revcomp("ATG") = "CAT", giving "CATCATCAT…", which has no stop
    /// either, so even both-strand search would stay empty — but pinning forward-only
    /// keeps the assertion about the all-start frame unambiguous.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_AllStartNoStop_RequireStartCodon_ReturnsEmpty()
    {
        // 30 nt of pure 'ATG' — in frame +1 it is ten start codons and never a stop.
        string allStart = string.Concat(Enumerable.Repeat("ATG", 10));

        var act = () => GenomeAnnotator.FindOrfs(allStart, minLength: 1, searchBothStrands: false).ToList();
        act.Should().NotThrow("an unterminated start-only frame is scanned to its end without indexing past it");

        var orfs = GenomeAnnotator.FindOrfs(allStart, minLength: 1, searchBothStrands: false).ToList();
        orfs.Should().BeEmpty(
            "INV-02 / §6.1 row 4: with requireStartCodon = true a stop is REQUIRED, so a start that never reaches a stop is not reported");
    }

    /// <summary>
    /// BE (other half of the contract): the SAME all-start-no-stop input WITH
    /// requireStartCodon = false. Here the canonical finder additionally emits trailing
    /// frame segments that run to the codon-aligned sequence end even without a stop
    /// (GenomeAnnotator.cs lines 183–212). Pinning this half proves the default
    /// (requireStartCodon = true) suppression above is a deliberate switch, not an
    /// accident, and that the partial-ORF path is itself disciplined: the reported span
    /// is codon-aligned (INV-03) and never runs past the sequence end (INV-04). For
    /// "ATGATGATG" (9 nt, frame +1) the trailing segment from the first pending start
    /// runs to the aligned end, length divisible by 3.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_AllStartNoStop_NoRequireStartCodon_EmitsCodonAlignedInBoundsPartialOrf()
    {
        string allStart = string.Concat(Enumerable.Repeat("ATG", 3)); // "ATGATGATG", 9 nt

        var orfs = GenomeAnnotator.FindOrfs(allStart, minLength: 0, searchBothStrands: false, requireStartCodon: false).ToList();

        orfs.Should().NotBeEmpty(
            "with requireStartCodon = false the canonical finder emits the trailing unterminated frame segment");
        orfs.Should().OnlyContain(o =>
            (o.End - o.Start) % 3 == 0 &&                       // INV-03: codon-aligned span
            o.Start >= 0 && o.Start < o.End && o.End <= allStart.Length, // INV-04: in bounds
            "even the partial-ORF path stays codon-aligned and within the sequence; it never runs off the end");
    }

    #endregion

    #region BE — Boundary: minLength = 0 (degenerate threshold)

    /// <summary>
    /// BE: minLength = 0 is the degenerate length floor. `orfAaLength >= 0` is ALWAYS
    /// true, so a zero-amino-acid ORF — a start codon immediately followed by an in-frame
    /// stop, e.g. "ATGTAA" with L_aa = (3 − 0)/3 − 1 = 0 — IS reported. The canonical
    /// contract has NO validation gate on minLength (ORF_Detection.md §3.3), so the
    /// disciplined outcome is NOT an exception but a WELL-FORMED trivial ORF: it must
    /// begin with a start, end with a stop, span exactly 6 nt (start + stop), be
    /// codon-aligned, and stay in bounds. We pin that the degenerate threshold yields a
    /// clean trivial result, never a crash and never a blow-up of garbage.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_MinLengthZero_EmitsWellFormedTrivialOrf()
    {
        const string seq = "ATGTAA"; // start immediately followed by stop → 0-aa ORF

        var act = () => GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false).ToList();
        act.Should().NotThrow("minLength = 0 is not rejected by the canonical finder; it admits the 0-aa ORF without crashing");

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false).ToList();

        var trivial = orfs.Should().ContainSingle(
            "exactly one frame (+1) contains the immediate start→stop pair").Subject;
        trivial.Sequence.Should().Be("ATGTAA", "the 0-aa ORF spans the start codon plus the terminal stop");
        trivial.Start.Should().Be(0);
        trivial.End.Should().Be(6);
        ((trivial.End - trivial.Start) % 3).Should().Be(0, "INV-03: the span is codon-aligned");
        trivial.End.Should().BeLessThanOrEqualTo(seq.Length, "INV-04: the ORF stays within the sequence bounds");
    }

    /// <summary>
    /// BE: minLength = 0 must NOT turn into a non-terminating or out-of-range scan on a
    /// longer mixed sequence — the threshold relaxation must not derail the codon walk.
    /// On a fixed-seed random 600-nt sequence the finder must complete promptly and every
    /// reported ORF (however short) must still satisfy all four invariants. This is the
    /// "degenerate parameter does not corrupt the scan" guard for minLength = 0.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindOrfs_MinLengthZero_RandomSequence_AllResultsWellFormed()
    {
        string seq = RandomDna(600, seed: 28_001);

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: true).ToList();

        orfs.Should().OnlyContain(o =>
            (o.End - o.Start) % 3 == 0 &&                           // INV-03
            o.Start >= 0 && o.Start < o.End && o.End <= seq.Length, // INV-04
            "even at the degenerate minLength = 0 every ORF on random input is codon-aligned and in bounds");
    }

    #endregion

    #region MC — Malformed content: non-DNA characters

    /// <summary>
    /// MC: non-DNA characters (N, IUPAC ambiguity codes, digits, punctuation) embedded
    /// in the sequence must NOT crash. The canonical finder does not pre-validate the
    /// alphabet (ORF_Detection.md §3.3): each 3-char codon window is upper-cased and
    /// tested for exact membership in the fixed A/C/G/T codon sets, so a triplet
    /// containing a non-ACGT character simply fails start/stop matching and is treated as
    /// an ordinary non-boundary codon — no out-of-range Substring, no spurious start or
    /// stop invented from garbage. Here the whole sequence is non-DNA filler with no
    /// valid start, so the result is empty on both strands.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NonDnaCharacters_DoesNotCrashAndEmitsNoSpuriousOrf()
    {
        // Length divisible by 3 (12) but full of non-ACGT garbage and ambiguity codes;
        // no triplet here is a start codon, so nothing can be emitted.
        const string garbage = "NNNRYK123!?WSM";

        var act = () => GenomeAnnotator.FindOrfs(garbage, minLength: 0, searchBothStrands: true).ToList();
        act.Should().NotThrow("non-DNA codons fail set membership but never index past the 3-char window or throw");

        var orfs = GenomeAnnotator.FindOrfs(garbage, minLength: 0, searchBothStrands: true).ToList();
        orfs.Should().BeEmpty(
            "no non-ACGT triplet is a start codon, so malformed filler invents no spurious ORF");
    }

    /// <summary>
    /// MC: an ORF whose INTERNAL coding codons contain non-DNA characters must still
    /// close at a genuine downstream stop. Only start/stop boundaries are matched against
    /// the fixed codon sets, so an internal 'NNN' codon is just a non-boundary codon and
    /// the ORF continues until a recognized stop (ORF_Detection.md §6.1 row 7). We embed
    /// "ATG" + "NNN" + "TAA" in frame +1: the ORF must open at the ATG, ride through the
    /// ambiguous internal codon, and close at the TAA — codon-aligned and in bounds.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NonDnaInternalCodon_OrfStillClosesAtGenuineStop()
    {
        const string seq = "ATGNNNTAA"; // start, ambiguous internal codon, stop — frame +1

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false, requireStartCodon: true).ToList();

        var orf = orfs.Should().ContainSingle(o => o.Start == 0).Subject;
        orf.Sequence.Should().Be("ATGNNNTAA",
            "the ORF opens at ATG, rides through the ambiguous internal NNN codon, and closes at the genuine TAA stop");
        orf.End.Should().Be(9);
        ((orf.End - orf.Start) % 3).Should().Be(0, "INV-03: the span is codon-aligned despite the internal garbage");
        orf.End.Should().BeLessThanOrEqualTo(seq.Length, "INV-04: the ORF stays within the sequence bounds");
    }

    #endregion

    #region Positive sanity — a clean ORF is found with correct boundaries

    /// <summary>
    /// Positive sanity: alongside the degenerate/malformed probes, a textbook ORF —
    /// "ATG" + a run of coding codons + "TAA" — must be detected with the CORRECT start,
    /// stop, frame, codon-aligned length and translated protein, so the boundary
    /// hardening never silently breaks the core function. We build
    /// ATG (GCT)×6 TAA in frame +1: start at 0, stop closing the 8-codon ORF, span 24 nt
    /// (divisible by 3), 6 internal amino acids (L_aa = 6 ≥ minLength). Pinned per
    /// INV-01..INV-04 (ORF_Detection.md §2.4): starts ATG, ends TAA, span %3 == 0, in
    /// bounds; ProteinSequence includes the terminal '*' stop symbol (ORF_Detection.md
    /// §3.2, §5.2).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_CleanOrf_DetectedWithCorrectStartStopFrameAndLength()
    {
        // ATG + (GCT)×6 [Ala] + TAA  → 1 + 6 + 1 = 8 codons, 24 nt, L_aa = 6.
        string coding = string.Concat(Enumerable.Repeat("GCT", 6));
        string seq = "ATG" + coding + "TAA";

        var orf = GenomeAnnotator.FindOrfs(seq, minLength: 6, searchBothStrands: false, requireStartCodon: true)
            .Should().ContainSingle("frame +1 holds exactly one clean ATG…TAA ORF at minLength = 6").Subject;

        orf.Start.Should().Be(0, "the ORF opens at the leading ATG");
        orf.End.Should().Be(seq.Length, "the ORF closes at the terminal TAA, spanning the whole sequence");
        orf.Frame.Should().Be(1, "an ORF starting at index 0 is in reading frame +1");
        orf.IsReverseComplement.Should().BeFalse("this ORF is on the forward strand");
        orf.Sequence.Should().StartWith("ATG", "INV-01: the ORF begins with a start codon");
        orf.Sequence.Should().EndWith("TAA", "INV-02: the ORF ends with a stop codon");
        ((orf.End - orf.Start) % 3).Should().Be(0, "INV-03: the ORF span is divisible by 3");
        (orf.Start >= 0 && orf.Start < orf.End && orf.End <= seq.Length).Should().BeTrue(
            "INV-04: 0 ≤ Start < End ≤ sequence length");
        orf.ProteinSequence.Should().EndWith("*",
            "the canonical translation includes the terminal '*' stop symbol for a stop-terminated ORF");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — every reported ORF begins with a start codon,
    /// ends with a stop codon, is codon-aligned and in bounds — so the degenerate-boundary
    /// and malformed-content probes do not corrupt the scan on ordinary input. Pinned per
    /// INV-01..INV-04 (ORF_Detection.md §2.4) with the default requireStartCodon = true.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindOrfs_RandomSequence_ProducesOnlyWellFormedOrfs()
    {
        var starts = new[] { "ATG", "GTG", "TTG" };
        var stops = new[] { "TAA", "TAG", "TGA" };
        string seq = RandomDna(1500, seed: 28_011);

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true, requireStartCodon: true).ToList();

        orfs.Should().OnlyContain(o =>
            o.Sequence.Length >= 6 &&                                       // start + stop minimum
            starts.Contains(o.Sequence.Substring(0, 3)) &&                  // INV-01
            stops.Contains(o.Sequence.Substring(o.Sequence.Length - 3, 3)) && // INV-02
            (o.End - o.Start) % 3 == 0 &&                                   // INV-03
            o.Start >= 0 && o.Start < o.End && o.End <= seq.Length,         // INV-04
            "every ORF on random input begins with a start, ends with a stop, is codon-aligned, and stays in bounds");
    }

    #endregion

    #endregion
}
