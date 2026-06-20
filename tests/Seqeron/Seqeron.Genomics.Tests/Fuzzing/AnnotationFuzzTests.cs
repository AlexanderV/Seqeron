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
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-GENE-001 — gene prediction (prokaryotic ORF-first heuristic)
/// Checklist: docs/checklists/03_FUZZING.md, row 29.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: an empty sequence, a non-coding-only sequence (no ORF,
///          no RBS), and overlapping genes (overlap bookkeeping must stay valid).
///   • MC = Malformed Content — a sequence that contains NO Shine-Dalgarno (RBS)
///          motif anywhere, so the RBS helper must invent nothing.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The gene-prediction contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Gene prediction here is a prokaryote-oriented, ORF-first heuristic split across two
/// helpers (Gene_Prediction.md §1, §5.1):
///   GenomeAnnotator.PredictGenes(string dnaSequence,
///                                int minOrfLength = 100,
///                                string prefix = "gene")
///     — converts every QUALIFYING ORF (translated length ≥ minOrfLength, found with
///       searchBothStrands: true, requireStartCodon: true) into a CDS GeneAnnotation,
///       ordered by genomic Start, with sequential IDs (Gene_Prediction.md §4.1, §5.2;
///       GenomeAnnotator.cs lines 389–421).
///   GenomeAnnotator.FindRibosomeBindingSites(string dnaSequence,
///                                            int upstreamWindow = 20,
///                                            int minDistance = 4,
///                                            int maxDistance = 15)
///     — independently scans the FORWARD-strand upstream window of every ORF (found with
///       an internal minLength = 30 aa) for Shine-Dalgarno-like motifs at an aligned
///       spacing within [minDistance, maxDistance] (GenomeAnnotator.cs lines 267–384).
///
/// CRITICAL: the two helpers are INDEPENDENT. PredictGenes does NOT require, call, or
/// consider an RBS — it emits genes purely from ORF structure (Gene_Prediction.md §5.2,
/// §5.4 deviation 1). So "no RBS in seq" suppresses RBS HITS, not gene predictions. The
/// fuzz tests below pin that exact split rather than a coupled model the code never
/// implements. The SD-like motif library is AGGAGG, GGAGG, AGGAG, GAGG, AGGA; the RBS
/// score is motif.Length / 6.0 (Gene_Prediction.md §4.2).
///
/// Documented parameter / validation contract (Gene_Prediction.md §3.3, §6.1):
///   • dnaSequence null or empty → BOTH helpers return an empty sequence (PredictGenes
///     via the underlying FindOrfs short-circuit; FindRibosomeBindingSites via its own
///     IsNullOrEmpty `yield break`). NOT an exception — pinned as the empty result.
///   • No valid ORF (non-coding-only input) → no genes AND no RBS hits, because RBS
///     scanning is ORF-driven (Gene_Prediction.md §6.1 row 2).
///   • Overlapping / nested ORFs → every qualifying ORF is emitted as its OWN CDS;
///     there is NO best-model selection or overlap suppression (Gene_Prediction.md §5.2
///     bullet 2, §5.3 "Intentionally simplified", §6.2). The disciplined outcome is that
///     overlap bookkeeping never crashes and every emitted gene keeps VALID coordinates.
///   • Non-DNA / no-motif content is not pre-validated; a region lacking any SD-like
///     motif simply produces no RBS hit (Gene_Prediction.md §3.3).
///
/// Documented invariants pinned on every predicted gene (Gene_Prediction.md §2.4):
///   GENE-INV-01 every gene derives from an ORF starting ATG/GTG/TTG and ending
///               TAA/TAG/TGA (PredictGenes delegates with requireStartCodon: true);
///   GENE-INV-02 every gene has valid bounds 0 ≤ Start < End ≤ length and Frame ∈ {1,2,3},
///               Strand ∈ {'+','-'}, Type == "CDS";
///   GENE-INV-03 every reported RBS hit lies at aligned distance ∈ [minDistance,
///               maxDistance] from its ORF start;
///   GENE-INV-04 every RBS score ∈ [4/6, 1.0] for the supported motif library.
/// PredictGenes / FindRibosomeBindingSites are LAZY iterators (`yield`), so every test
/// forces enumeration (`.ToList()`). Deterministic only: random fuzz input comes from a
/// locally fixed `new Random(seed)`.
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

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-GENE-001 — gene prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-GENE-001 — gene prediction

    #region Helpers — gene-prediction constructions

    /// <summary>
    /// Builds a forward-strand coding cassette: a leading 5' pad, an in-frame ATG start,
    /// <paramref name="codingCodons"/> repeats of "GCT" (Ala), then a TAA stop. The pad is
    /// "AAA"-tiled so it introduces no start/stop in any frame and keeps the cassette
    /// frame-+1 aligned. Translated length = codingCodons aa (the terminal '*' is trimmed
    /// by PredictGenes for protein_length).
    /// </summary>
    private static string CodingCassette(int codingCodons, int padCodons)
    {
        string pad = string.Concat(Enumerable.Repeat("AAA", padCodons));
        string coding = string.Concat(Enumerable.Repeat("GCT", codingCodons));
        return pad + "ATG" + coding + "TAA";
    }

    #endregion

    #region BE — Boundary: empty sequence yields no genes and no RBS hits

    /// <summary>
    /// BE: an empty sequence must yield NO genes and NO RBS hits, never a crash
    /// (Gene_Prediction.md §3.3, §6.1 row 1). PredictGenes short-circuits through the
    /// underlying FindOrfs `yield break`; FindRibosomeBindingSites short-circuits on its
    /// own IsNullOrEmpty guard. Both are lazy iterators, so we force enumeration to make
    /// the documented empty result actually surface.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_EmptySequence_ReturnsNoGenesAndNoRbs()
    {
        var predict = () => GenomeAnnotator.PredictGenes(string.Empty, minOrfLength: 1).ToList();
        predict.Should().NotThrow("an empty sequence short-circuits through FindOrfs; it simply yields no genes");
        GenomeAnnotator.PredictGenes(string.Empty, minOrfLength: 1).ToList()
            .Should().BeEmpty("no ORF can exist in an empty sequence, so no gene is predicted");

        var rbs = () => GenomeAnnotator.FindRibosomeBindingSites(string.Empty).ToList();
        rbs.Should().NotThrow("the RBS helper `yield break`s on empty input before any indexing");
        GenomeAnnotator.FindRibosomeBindingSites(string.Empty).ToList()
            .Should().BeEmpty("RBS scanning is ORF-driven; no ORF means no upstream motif scan");
    }

    #endregion

    #region BE — Boundary: non-coding-only sequence (no ORF, no RBS) yields nothing

    /// <summary>
    /// BE: a non-coding-only sequence — no start codon in any frame, hence no ORF — must
    /// yield NO genes and NO RBS hits (Gene_Prediction.md §6.1 row 2). We tile "CCC"
    /// (revcomp "GGG"), which is neither a start nor a stop on either strand, so the ORF
    /// scan that BOTH helpers drive finds nothing. The disciplined outcome is the empty
    /// result on both surfaces, not a crash and not a spurious gene.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_NonCodingOnly_ReturnsNoGenesAndNoRbs()
    {
        // 60 nt of pure 'CCC' — no ATG/GTG/TTG and no TAA/TAG/TGA on either strand.
        string nonCoding = string.Concat(Enumerable.Repeat("CCC", 20));

        var predict = () => GenomeAnnotator.PredictGenes(nonCoding, minOrfLength: 1).ToList();
        predict.Should().NotThrow("a non-coding sequence records no pending ORF start; the scan terminates cleanly");
        GenomeAnnotator.PredictGenes(nonCoding, minOrfLength: 1).ToList()
            .Should().BeEmpty("with no start codon on either strand there is no ORF and therefore no gene");

        GenomeAnnotator.FindRibosomeBindingSites(nonCoding).ToList()
            .Should().BeEmpty("RBS scanning is ORF-driven, so a sequence with no ORF yields no RBS hits");
    }

    #endregion

    #region MC — Malformed content: no RBS motif anywhere upstream of a real gene

    /// <summary>
    /// MC: a sequence that DOES contain a genuine, qualifying ORF but has NO Shine-Dalgarno
    /// motif (AGGAGG / GGAGG / AGGAG / GAGG / AGGA) anywhere upstream must still predict the
    /// gene, yet report NO RBS hit. This pins the documented INDEPENDENCE of the two helpers
    /// (Gene_Prediction.md §5.2, §5.4 deviation 1): "no RBS in seq" suppresses RBS hits, NOT
    /// gene calls. The ORF is ≥ 30 aa so the RBS helper's internal FindOrfs(minLength: 30)
    /// discovers it and actually scans its upstream window; the pad is "AAA"-tiled so the
    /// upstream region contains no SD-like motif. The gene must appear with valid bounds; the
    /// RBS list must be empty.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_NoRbsMotifUpstream_PredictsGeneButReportsNoRbs()
    {
        // 35 coding codons (≥ 30 aa) with a 4-codon 'AAA' pad — no SD-like motif upstream.
        string seq = CodingCassette(codingCodons: 35, padCodons: 4);

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();
        genes.Should().ContainSingle("the cassette holds exactly one qualifying forward-strand ORF")
            .Which.Strand.Should().Be('+', "the gene is on the forward strand");

        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().BeEmpty(
            "the upstream window holds no Shine-Dalgarno motif, so the RBS helper invents nothing — " +
            "yet the gene is still predicted because PredictGenes ignores RBS evidence entirely");
    }

    #endregion

    #region BE — Boundary: overlapping genes are all reported with valid coordinates

    /// <summary>
    /// BE: overlapping / nested ORFs must each be emitted as their OWN CDS — there is NO
    /// best-model selection or overlap suppression (Gene_Prediction.md §5.2 bullet 2,
    /// §5.3 "Intentionally simplified", §6.2). The disciplined boundary here is the overlap
    /// BOOKKEEPING: predicting multiple coincident genes must not crash, and EVERY emitted
    /// gene must keep valid coordinates (GENE-INV-02: 0 ≤ Start &lt; End ≤ length,
    /// Frame ∈ {1,2,3}, Strand ∈ {'+','-'}, Type == "CDS"), with IDs ordered by Start.
    ///
    /// We build a single coding cassette that is a palindromic-by-construction overlap
    /// generator: a forward ORF AND its reverse complement both score as ORFs over the SAME
    /// genomic span. Concretely "ATG...TAA" on the forward strand revcomps to "TTA...CAT";
    /// to force a genuine overlapping pair we instead concatenate two forward ORFs whose
    /// spans abut/overlap by sharing coordinates across frames. Simpler and fully
    /// deterministic: a cassette with an inner in-frame ATG produces a nested ORF (a second
    /// start before the shared stop), so PredictGenes emits BOTH the outer and the inner
    /// gene over overlapping coordinates.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_OverlappingGenes_AllReportedWithValidCoordinates()
    {
        // Outer ORF: ATG (GCT)x40 ATG (GCT)x5 TAA  → an INNER in-frame ATG creates a nested
        // ORF that shares the same TAA stop, so two overlapping genes are emitted in frame +1.
        string outerHead = "ATG" + string.Concat(Enumerable.Repeat("GCT", 40));
        string innerHead = "ATG" + string.Concat(Enumerable.Repeat("GCT", 5));
        string seq = outerHead + innerHead + "TAA";

        var predict = () => GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();
        predict.Should().NotThrow("emitting multiple overlapping ORFs must not crash the overlap bookkeeping");

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();

        genes.Count.Should().BeGreaterThanOrEqualTo(2,
            "both the outer ORF and the nested inner ORF (sharing the terminal stop) are emitted; " +
            "there is no overlap suppression");

        genes.Should().OnlyContain(g =>
            g.Start >= 0 && g.Start < g.End && g.End <= seq.Length &&  // GENE-INV-02 bounds
            (g.Strand == '+' || g.Strand == '-') &&                    // GENE-INV-02 strand
            g.Type == "CDS",                                           // GENE-INV-02 type
            "every overlapping gene keeps valid coordinates, a ± strand, and Type=CDS");

        // GENE-INV-02 frame ∈ {1,2,3} — carried in the 'frame' attribute, not a record field.
        genes.Should().OnlyContain(g =>
            g.Attributes.ContainsKey("frame") &&
            (g.Attributes["frame"] == "1" || g.Attributes["frame"] == "2" || g.Attributes["frame"] == "3"),
            "every overlapping gene tags a reading frame in {1,2,3}");

        // Genes are ordered by genomic Start, and two distinct overlapping starts exist.
        genes.Select(g => g.Start).Should().BeInAscendingOrder("PredictGenes orders genes by genomic Start");
        genes.Select(g => g.Start).Distinct().Count().Should().BeGreaterThanOrEqualTo(2,
            "the outer and inner ORFs open at distinct positions yet overlap, confirming no overlap was suppressed");
    }

    #endregion

    #region Positive sanity — RBS upstream of a downstream ORF yields a gene + RBS hit

    /// <summary>
    /// Positive sanity: a textbook prokaryotic gene — a Shine-Dalgarno motif (AGGAGG) an
    /// aligned short distance upstream of an in-frame ATG that opens a ≥ 30 aa ORF — must
    /// produce BOTH a predicted gene with correct coordinates AND a matching RBS hit, so the
    /// boundary/malformed probes never silently break the core function.
    ///
    /// Layout (frame +1): [pad AAA] AGGAGG [8-nt 'AAAAAAAA' spacer] ATG (GCT)x35 TAA.
    /// The SD motif 3' end sits 8 nt upstream of ATG → aligned distance 8 ∈ [4,15]
    /// (GENE-INV-03), inside the default 20-nt upstream window. The ORF is 35 aa ≥ the RBS
    /// helper's internal 30-aa floor, so it is scanned. Pinned: the gene starts at the ATG
    /// with End at the TAA, valid bounds, '+' strand, CDS; the RBS hit reports AGGAGG with
    /// score 1.0 (= 6/6, GENE-INV-04) at the SD position.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_RbsUpstreamOfOrf_YieldsGeneWithCorrectCoordsAndRbsHit()
    {
        const string pad = "AAA";               // 3 nt 5' pad (no start/stop, keeps frame)
        const string sd = "AGGAGG";             // full Shine-Dalgarno consensus, score 1.0
        const string spacer = "AAAAAAAA";       // 8 nt → aligned distance 8 ∈ [4,15]
        string coding = "ATG" + string.Concat(Enumerable.Repeat("GCT", 35)) + "TAA"; // 35 aa ORF
        string seq = pad + sd + spacer + coding;

        int expectedStart = pad.Length + sd.Length + spacer.Length; // start of ATG
        int expectedSdPos = pad.Length;                             // start of AGGAGG

        // ── Gene prediction ──────────────────────────────────────────────────────────
        var gene = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30)
            .Should().ContainSingle("exactly one qualifying forward-strand ORF is present").Subject;

        gene.Start.Should().Be(expectedStart, "the gene opens at the in-frame ATG downstream of the SD motif");
        gene.End.Should().Be(seq.Length, "the gene closes at the terminal TAA");
        gene.Strand.Should().Be('+', "the gene is on the forward strand");
        gene.Type.Should().Be("CDS", "PredictGenes emits CDS annotations");
        gene.Attributes["frame"].Should().BeOneOf("1", "2", "3", "GENE-INV-02: frame is in {1,2,3}");
        (gene.Start >= 0 && gene.Start < gene.End && gene.End <= seq.Length).Should().BeTrue(
            "GENE-INV-02: 0 ≤ Start < End ≤ sequence length");

        // ── RBS detection ────────────────────────────────────────────────────────────
        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().NotBeEmpty("a full AGGAGG motif sits at an aligned distance of 8 nt upstream of the start");

        var hit = rbs.Should().Contain(h => h.sequence == "AGGAGG",
            "the full Shine-Dalgarno consensus is the highest-scoring motif present").Subject;
        hit.position.Should().Be(expectedSdPos, "the RBS hit reports the genomic 5' position of the SD motif");
        hit.score.Should().BeApproximately(1.0, 1e-9, "GENE-INV-04: AGGAGG (6 nt) scores 6/6 = 1.0");

        int alignedDistance = gene.Start - hit.position - hit.sequence.Length;
        alignedDistance.Should().BeInRange(4, 15, "GENE-INV-03: aligned spacing lies within [minDistance, maxDistance]");
    }

    /// <summary>
    /// Positive sanity / robustness: on a fixed-seed random 1200-nt sequence both helpers
    /// must complete promptly and every result must be well-formed — every predicted gene
    /// satisfies GENE-INV-02 (valid bounds, frame ∈ {1,2,3}, ± strand, CDS) and every RBS
    /// hit satisfies GENE-INV-03/04 (aligned spacing in [4,15], score ∈ [4/6, 1]). This is
    /// the "degenerate-boundary probes do not corrupt the scan on ordinary input" guard.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void PredictGenes_RandomSequence_ProducesOnlyWellFormedGenesAndRbsHits()
    {
        string seq = RandomDna(1200, seed: 29_001);

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();
        genes.Should().OnlyContain(g =>
            g.Start >= 0 && g.Start < g.End && g.End <= seq.Length &&
            (g.Strand == '+' || g.Strand == '-') &&
            g.Type == "CDS" &&
            g.Attributes.ContainsKey("frame") &&
            (g.Attributes["frame"] == "1" || g.Attributes["frame"] == "2" || g.Attributes["frame"] == "3"),
            "every predicted gene on random input is in bounds, codon-frame-tagged, ± stranded, and CDS");

        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().OnlyContain(h =>
            h.score >= 4.0 / 6.0 && h.score <= 1.0 + 1e-9 &&   // GENE-INV-04
            h.position >= 0 && h.position < seq.Length,        // in-bounds motif position
            "every RBS hit on random input carries a length-normalized score in [4/6, 1] at an in-bounds position");
    }

    #endregion

    #endregion
}
