using System;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Alignment area — pairwise GLOBAL alignment (Needleman-Wunsch)
/// and pairwise LOCAL alignment (Smith-Waterman).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds boundary and malformed inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no infinite loop, no state
/// corruption, no IndexOutOfRange while indexing the DP matrix, no OutOfMemory on
/// long inputs, and no *unhandled* runtime exception leaking from the dynamic
/// program or traceback. Every input must yield EITHER a well-defined,
/// theory-correct alignment, OR a *documented, intentional* validation exception
/// (here: ArgumentNullException on a null DnaSequence). A raw runtime exception, a
/// structurally inconsistent alignment, or a hang on boundary input is a bug, not
/// a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ALIGN-GLOBAL-001 — global alignment (Needleman-Wunsch) (Alignment)
/// Checklist: docs/checklists/03_FUZZING.md, row 35.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — empty vs empty, empty vs non-empty, the
///          single-character vs single-character extreme, and very long sequences
///          (the O(mn) DP must still complete in bounded time/memory). These are
///          the size boundaries of the dynamic-programming matrix.
///   • MC = Malformed Content — non-DNA characters on the lenient raw-string
///          surface (matched by plain character comparison, never a crash) and the
///          strict typed surface (the DnaSequence ctor REJECTS them up front), plus
///          the documented null-DnaSequence validation gate.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The global-alignment (Needleman-Wunsch) contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Global alignment aligns BOTH sequences end-to-end, inserting gaps until the two
/// aligned strings occupy one common coordinate system. The recurrence (linear gap
/// penalty d = ScoringMatrix.GapExtend; GapOpen is NOT used — the repo implements a
/// LINEAR, not affine, gap model):
///   F(i,j) = max( F(i-1,j-1) + S(Aᵢ,Bⱼ),  F(i-1,j) + d,  F(i,j-1) + d )
/// with boundary F(i,0) = d·i, F(0,j) = d·j; traceback starts at F(m,n).
///   — docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md §2.2, §2.4,
///     §3, §5, §6, §7. Sources: Needleman-Wunsch algorithm; Sequence alignment
///     (Wikipedia).
///
/// KEY STRUCTURAL INVARIANTS asserted across the fuzz fodder (the contract that
/// must never silently break):
///   • INV-01 — the two aligned output strings have EQUAL length (gaps pad both
///     into one coordinate system).
///   • INV-02 — removing the gap character '-' from each aligned output reproduces
///     the (upper-cased) input sequence exactly (traceback emits original
///     characters or gaps; it never rewrites a base).
///   • INV-03 — under the linear model the reported Score equals the per-column
///     sum of match/mismatch/gap contributions of the returned alignment.
///   • INV-04 — traceback covers the FULL length of both inputs (end-to-end).
///   — Global_Alignment_Needleman_Wunsch.md §2.4 (INV-01..INV-04).
///
/// DOCUMENTED SURFACE-SPECIFIC BOUNDARY behaviour (pinned so the strict/lenient and
/// typed/string splits cannot drift — §3.3, §5.2, §6.1):
///   • GlobalAlign(string, string, …): if EITHER raw string is null or empty it
///     returns AlignmentResult.Empty (early exit), WITHOUT DnaSequence alphabet
///     validation. It upper-cases the inputs but does NOT enforce the DNA alphabet,
///     so non-DNA characters flow through and are aligned by plain comparison.
///   • GlobalAlign(DnaSequence, DnaSequence) — NON-cancellation typed overload:
///     null input → ArgumentNullException; empty DnaSequence is a LEGAL input that
///     is aligned through the normal traceback path (empty vs empty → trivial
///     empty alignment, score 0; empty vs length-n → an all-gap alignment of score
///     n·d). End coordinates are sequence.Length − 1 (so −1 for an empty input).
///   • GlobalAlign(DnaSequence, DnaSequence, scoring, CancellationToken, progress)
///     — cancellation-aware typed overload: it delegates to the raw-string
///     cancellation path, so an empty typed input returns AlignmentResult.Empty.
///   • DnaSequence(string) — the strict typed surface — upper-cases and REJECTS any
///     non-A/C/G/T character with ArgumentException, so a non-DNA sequence can never
///     reach the typed aligner at all (the MC rejection happens at construction).
///
/// Default scoring is SequenceAligner.SimpleDna: Match=+1, Mismatch=−1,
/// GapExtend=−1 (so the gap penalty per position d = −1). All inputs here are ASCII
/// DNA / boundary strings; randomness (where used) is from a locally fixed-seed
/// Random so the fuzz fodder is fully reproducible and adding tests cannot perturb
/// any other fixture.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ALIGN-LOCAL-001 — local alignment (Smith-Waterman) (Alignment)
/// Checklist: docs/checklists/03_FUZZING.md, row 36.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the boundary fodder that strains the
///          zero-floored DP: empty seqs (either side), identical seqs, completely
///          DIFFERENT seqs (no shared base), and the 1-char vs 1-char extreme.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The local-alignment (Smith-Waterman) contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Local alignment returns the best-scoring SUBsequence pair, not an end-to-end
/// alignment. The recurrence adds a ZERO FLOOR so a negative-running region is
/// dropped and the alignment can restart (linear gap penalty d = ScoringMatrix.GapExtend):
///   H(i,j) = max( 0,  H(i-1,j-1)+S(Aᵢ,Bⱼ),  H(i-1,j)+d,  H(i,j-1)+d )
/// with H(i,0)=H(0,j)=0; traceback starts at the matrix MAXIMUM cell and stops at
/// the first zero-valued cell. — docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md
/// §2.2, §2.4, §4.1, §4.2. Sources: Smith-Waterman algorithm; Sequence alignment.
///
/// KEY LOCAL invariants asserted across the fuzz fodder (Local_Alignment_Smith_Waterman.md §2.4):
///   • L-INV-01 — every reported local-alignment Score is ≥ 0 (the KEY zero-floor
///     property; a negative-scoring region is dropped, never propagated).
///   • L-INV-02 — the two aligned output strings have EQUAL length.
///   • L-INV-03 — stripping gaps from each aligned output yields a CONTIGUOUS
///     SUBSTRING of the corresponding (upper-cased) input (local = subsequence pair,
///     not the full input and never a reordered/rewritten base).
///
/// DOCUMENTED LOCAL boundary behaviour (pinned so the strict/lenient and typed/string
/// splits cannot drift — Local_Alignment_Smith_Waterman.md §3.3, §5.2, §6.1):
///   • LocalAlign(string, string, …): a null OR empty raw string short-circuits to
///     AlignmentResult.Empty WITHOUT DnaSequence validation; otherwise the inputs are
///     upper-cased and aligned by plain character comparison.
///   • LocalAlign(DnaSequence, DnaSequence, …): null input → ArgumentNullException;
///     an EMPTY DnaSequence is a LEGAL input that flows through the core DP (no
///     short-circuit), yielding a Local result with empty aligned strings and all
///     coordinate fields = −1 (the zero-size traceback endpoints).
///   • When NO positive-scoring region exists (empty input, or completely different
///     seqs under negative mismatch/gap), the matrix max stays 0, traceback never
///     runs, Score = 0 and both aligned strings are empty (the zero floor).
///   • Identical sequences → the full diagonal dominates: aligned == upper-cased
///     input on BOTH sides, gap-free, Score = length × Match.
///   • AlignmentResult.Empty carries AlignmentType.Global (shared sentinel), so the
///     typed empty-DnaSequence path is the way to observe a tagged-Local empty result.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AlignmentFuzzTests
{
    #region Helpers

    /// <summary>
    /// Deterministic, LOCALLY-scoped RNG — seed fixed so generated fuzz inputs are
    /// reproducible. Kept private to this fixture and never shared/mutated across
    /// fixtures, so adding tests here cannot shift any other fixture's randomness.
    /// </summary>
    private static readonly Random Rng = new(20260620);

    /// <summary>Generates a random valid DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the four global-alignment structural invariants against the raw
    /// (upper-cased) inputs the aligner actually consumed: equal aligned length
    /// (INV-01), gap-stripped aligned strings reproduce the inputs (INV-02, INV-04),
    /// and the reported Score equals the per-column linear-model recomputation
    /// (INV-03). Used by every non-empty-result fuzz case.
    /// </summary>
    private static void AssertGlobalInvariants(
        AlignmentResult result, string input1, string input2, ScoringMatrix scoring)
    {
        string expected1 = input1.ToUpperInvariant();
        string expected2 = input2.ToUpperInvariant();

        result.AlignmentType.Should().Be(AlignmentType.Global,
            "a global alignment is tagged Global");

        // INV-01: the two aligned strings share one coordinate system.
        result.AlignedSequence1.Length.Should().Be(result.AlignedSequence2.Length,
            "INV-01: global alignment pads both sequences to equal length");

        // INV-02 / INV-04: removing gaps recovers the FULL original inputs.
        new string(result.AlignedSequence1.Where(c => c != '-').ToArray())
            .Should().Be(expected1,
                "INV-02/INV-04: gap-stripped aligned sequence 1 must reproduce the full input (traceback never rewrites a base)");
        new string(result.AlignedSequence2.Where(c => c != '-').ToArray())
            .Should().Be(expected2,
                "INV-02/INV-04: gap-stripped aligned sequence 2 must reproduce the full input (end-to-end)");

        // INV-03: the reported score equals the per-column linear-model sum.
        int recomputed = 0;
        for (int i = 0; i < result.AlignedSequence1.Length; i++)
        {
            char a = result.AlignedSequence1[i];
            char b = result.AlignedSequence2[i];
            if (a == '-' || b == '-')
                recomputed += scoring.GapExtend;
            else
                recomputed += a == b ? scoring.Match : scoring.Mismatch;
        }
        result.Score.Should().Be(recomputed,
            "INV-03: the reported score equals the sum of per-column match/mismatch/gap contributions");
    }

    /// <summary>
    /// Asserts the three LOCAL (Smith-Waterman) structural invariants against the raw
    /// (upper-cased) inputs the aligner actually consumed: Score ≥ 0 (L-INV-01, the
    /// zero floor), equal aligned length (L-INV-02), and each gap-stripped aligned
    /// segment is a CONTIGUOUS substring of the corresponding input (L-INV-03 — the
    /// local-subsequence property). Used by every local fuzz case.
    /// — Local_Alignment_Smith_Waterman.md §2.4.
    /// </summary>
    private static void AssertLocalInvariants(
        AlignmentResult result, string input1, string input2)
    {
        string expected1 = input1.ToUpperInvariant();
        string expected2 = input2.ToUpperInvariant();

        // L-INV-01: the zero floor guarantees a non-negative local score.
        result.Score.Should().BeGreaterThanOrEqualTo(0,
            "L-INV-01: the Smith-Waterman zero floor makes every local-alignment score ≥ 0");

        // L-INV-02: traceback emits one column per step, padding gaps on both sides.
        result.AlignedSequence1.Length.Should().Be(result.AlignedSequence2.Length,
            "L-INV-02: the two aligned local segments share one coordinate system (equal length)");

        // L-INV-03: gap-stripped aligned segments are CONTIGUOUS substrings of the
        // inputs. An EMPTY segment (zero-score / dropped region) is vacuously a
        // substring of any input, so we only assert containment for non-empty segments.
        string seg1 = new string(result.AlignedSequence1.Where(c => c != '-').ToArray());
        string seg2 = new string(result.AlignedSequence2.Where(c => c != '-').ToArray());
        if (seg1.Length > 0)
            expected1.Should().Contain(seg1,
                "L-INV-03: the gap-stripped aligned segment 1 is a contiguous substring of input 1 (local subsequence)");
        if (seg2.Length > 0)
            expected2.Should().Contain(seg2,
                "L-INV-03: the gap-stripped aligned segment 2 is a contiguous substring of input 2 (local subsequence)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ALIGN-GLOBAL-001 — global alignment (Needleman-Wunsch) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ALIGN-GLOBAL-001 — global alignment (Needleman-Wunsch)

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty vs empty
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty vs empty

    /// <summary>
    /// BE: empty vs empty is the lowest size boundary. On the STRICT typed
    /// (non-cancellation) overload an empty DnaSequence is a legal input aligned
    /// through the normal traceback path: the trivial alignment is two EMPTY aligned
    /// strings, score 0, type Global — never an IndexOutOfRange on the degenerate
    /// 1×1 DP matrix (Global_Alignment_Needleman_Wunsch.md §3.3, §6.1). INV-01..04
    /// hold trivially (both gap-stripped outputs are "").
    /// </summary>
    [Test]
    public void GlobalAlign_TypedEmptyVsEmpty_IsTrivialZeroScoreAlignment()
    {
        var result = SequenceAligner.GlobalAlign(new DnaSequence(""), new DnaSequence(""));

        result.AlignedSequence1.Should().BeEmpty("aligning two empty sequences yields empty aligned strings");
        result.AlignedSequence2.Should().BeEmpty("aligning two empty sequences yields empty aligned strings");
        result.Score.Should().Be(0, "there are no columns, so the global-alignment score is 0");
        AssertGlobalInvariants(result, "", "", SequenceAligner.SimpleDna);
    }

    /// <summary>
    /// BE: the lenient RAW-STRING surface short-circuits on null/empty input and
    /// returns AlignmentResult.Empty (Global_Alignment_Needleman_Wunsch.md §6.1) —
    /// the documented strict/lenient split versus the typed overload above (which
    /// traces an empty result). Pinned for empty-vs-empty AND the null variants so
    /// the early exit can never drift into a NullReferenceException.
    /// </summary>
    [Test]
    public void GlobalAlign_RawStringEmptyOrNull_ReturnsEmptyResult()
    {
        SequenceAligner.GlobalAlign("", "").Should().Be(AlignmentResult.Empty,
            "the raw-string overload exits early to AlignmentResult.Empty on empty input");
        SequenceAligner.GlobalAlign("ACGT", "").Should().Be(AlignmentResult.Empty,
            "an empty second raw string short-circuits before the DP runs");
        SequenceAligner.GlobalAlign("", "ACGT").Should().Be(AlignmentResult.Empty,
            "an empty first raw string short-circuits before the DP runs");

        var nullFirst = () => SequenceAligner.GlobalAlign((string)null!, "ACGT");
        var nullSecond = () => SequenceAligner.GlobalAlign("ACGT", (string)null!);
        nullFirst.Should().NotThrow("a null raw string is guarded as empty, never dereferenced");
        nullSecond.Should().NotThrow("a null raw string is guarded as empty, never dereferenced");
        nullFirst().Should().Be(AlignmentResult.Empty, "a null raw string yields the empty result");
        nullSecond().Should().Be(AlignmentResult.Empty, "a null raw string yields the empty result");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty vs non-empty sequence
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty vs sequence of length n (all-gap alignment)

    /// <summary>
    /// BE: empty vs a length-n sequence must produce the ALL-GAP alignment — n gap
    /// columns aligning the n bases of the non-empty side against '-', a Score of
    /// n·d (here d = GapExtend = −1, so −n), and aligned length exactly n
    /// (Global_Alignment_Needleman_Wunsch.md §2.2 boundary F(0,j)=d·j, §6.1). Pinned
    /// in BOTH orientations on the typed traceback surface. INV-01..04 confirm the
    /// gap-stripped non-empty side reproduces the input and the empty side is all
    /// gaps. No crash on the degenerate single-row / single-column DP.
    /// </summary>
    [TestCase("ACGT", TestName = "GlobalAlign_EmptyVsSeq_AllGap_Len4")]
    [TestCase("A", TestName = "GlobalAlign_EmptyVsSeq_AllGap_SingleBase")]
    [TestCase("ACGTACGTACGT", TestName = "GlobalAlign_EmptyVsSeq_AllGap_Len12")]
    public void GlobalAlign_TypedEmptyVsSequence_IsAllGapAlignment(string seq)
    {
        var scoring = SequenceAligner.SimpleDna;

        // empty (seq1) vs non-empty (seq2)
        var r1 = SequenceAligner.GlobalAlign(new DnaSequence(""), new DnaSequence(seq));
        r1.AlignedSequence1.Should().Be(new string('-', seq.Length),
            "the empty side aligns as n gap characters");
        r1.AlignedSequence2.Should().Be(seq.ToUpperInvariant(),
            "the non-empty side aligns as itself across the n gap columns");
        r1.Score.Should().Be(seq.Length * scoring.GapExtend,
            "the all-gap alignment scores n × gap penalty (boundary F(0,j) = d·j)");
        r1.AlignedSequence1.Length.Should().Be(seq.Length,
            "the aligned length equals n for an empty-vs-length-n alignment");
        AssertGlobalInvariants(r1, "", seq, scoring);

        // non-empty (seq1) vs empty (seq2) — symmetric.
        var r2 = SequenceAligner.GlobalAlign(new DnaSequence(seq), new DnaSequence(""));
        r2.AlignedSequence2.Should().Be(new string('-', seq.Length),
            "the empty side aligns as n gap characters in the mirror orientation");
        r2.AlignedSequence1.Should().Be(seq.ToUpperInvariant(),
            "the non-empty side aligns as itself in the mirror orientation");
        r2.Score.Should().Be(seq.Length * scoring.GapExtend,
            "the mirror all-gap alignment scores n × gap penalty (boundary F(i,0) = d·i)");
        AssertGlobalInvariants(r2, seq, "", scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: single char vs single char
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: single char vs single char (match / mismatch)

    /// <summary>
    /// BE: the minimal non-empty alignment — one base vs one base. The optimal
    /// global alignment is the single diagonal column (one match or one mismatch),
    /// NOT a two-gap detour, because under SimpleDna a mismatch (−1) still beats two
    /// gaps (−2): so the aligned length is exactly 1 and the Score is Match for
    /// equal bases or Mismatch for unequal bases
    /// (Global_Alignment_Needleman_Wunsch.md §2.2 recurrence, §7). Pinned for a
    /// match and for several mismatching pairs. INV-01..04 hold.
    /// </summary>
    [TestCase("A", "A", 1, TestName = "GlobalAlign_SingleChar_Match_A")]
    [TestCase("G", "G", 1, TestName = "GlobalAlign_SingleChar_Match_G")]
    [TestCase("A", "C", -1, TestName = "GlobalAlign_SingleChar_Mismatch_AC")]
    [TestCase("G", "T", -1, TestName = "GlobalAlign_SingleChar_Mismatch_GT")]
    public void GlobalAlign_TypedSingleCharVsSingleChar_IsOneColumn(string a, string b, int expectedScore)
    {
        var scoring = SequenceAligner.SimpleDna;
        var result = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.AlignedSequence1.Length.Should().Be(1,
            "the optimal one-base global alignment is a single diagonal column, not a two-gap detour");
        result.AlignedSequence1.Should().NotContain("-", "the optimal single-column alignment has no gap");
        result.AlignedSequence2.Should().NotContain("-", "the optimal single-column alignment has no gap");
        result.Score.Should().Be(expectedScore,
            "a single column scores Match for equal bases and Mismatch for unequal bases");
        AssertGlobalInvariants(result, a, b, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: non-DNA characters
    // ───────────────────────────────────────────────────────────────────

    #region MC — Malformed content: non-DNA characters (lenient string vs strict typed)

    /// <summary>
    /// MC: the LENIENT raw-string surface does NOT enforce the DNA alphabet beyond
    /// upper-casing (Global_Alignment_Needleman_Wunsch.md §3.3, §6.2). Non-DNA
    /// characters (digits, punctuation, IUPAC ambiguity codes, lowercase) must flow
    /// through and be aligned by plain character comparison — never a crash, never a
    /// KeyNotFound — and STILL satisfy every structural invariant: equal aligned
    /// length, gap-stripped outputs reproduce the upper-cased inputs, and the score
    /// equals the per-column recomputation. Two identical non-DNA strings still
    /// align gap-free (every column an equal-character match). (The gap symbol '-'
    /// itself is intentionally excluded from the fuzz fodder: it would be
    /// indistinguishable from an inserted gap when invariants strip gaps.)
    /// </summary>
    [TestCase("AC9T", "ACGT", TestName = "GlobalAlign_RawNonDna_Digit")]
    [TestCase("AC*GT", "ACGT", TestName = "GlobalAlign_RawNonDna_Punctuation")]
    [TestCase("ACRYT", "ACGTT", TestName = "GlobalAlign_RawNonDna_IupacCodes")]
    [TestCase("acgt", "ACGT", TestName = "GlobalAlign_RawNonDna_Lowercase_FoldsToUpper")]
    [TestCase("XYZ@!", "XYZ@!", TestName = "GlobalAlign_RawNonDna_IdenticalGarbage")]
    public void GlobalAlign_RawStringNonDna_AlignsWithoutCrashAndHoldsInvariants(string a, string b)
    {
        var scoring = SequenceAligner.SimpleDna;

        var act = () => SequenceAligner.GlobalAlign(a, b, scoring);
        act.Should().NotThrow(
            "the raw-string surface does not enforce the DNA alphabet; non-DNA input is aligned, not rejected");

        var result = act();
        AssertGlobalInvariants(result, a, b, scoring);

        // Lowercase folds to the same alignment as its uppercase form (case-insensitive).
        if (a == "acgt")
            result.Score.Should().Be(SequenceAligner.GlobalAlign("ACGT", b, scoring).Score,
                "the raw-string surface upper-cases input, so case must not change the score");

        // Two identical strings (even garbage) align gap-free with an all-match score.
        if (a == b)
        {
            result.AlignedSequence1.Should().Be(a.ToUpperInvariant(), "identical inputs align without gaps");
            result.Score.Should().Be(a.Length * scoring.Match,
                "every column of an identical-vs-identical alignment is a match");
        }
    }

    /// <summary>
    /// MC: the STRICT typed surface rejects non-DNA up front — the DnaSequence ctor
    /// validates the A/C/G/T alphabet and throws ArgumentException, so malformed
    /// content can never reach the typed aligner (Global_Alignment_Needleman_Wunsch.md
    /// §3.3). This is the documented validation gate, NOT a downstream KeyNotFound or
    /// IndexOutOfRange. Pinned against the lenient string surface above so the
    /// strict/lenient split is explicit and cannot silently converge.
    /// </summary>
    [TestCase("AC9T", TestName = "TypedNonDna_Digit_Rejected")]
    [TestCase("ACNGT", TestName = "TypedNonDna_AmbiguityCode_Rejected")]
    [TestCase("AC GT", TestName = "TypedNonDna_Space_Rejected")]
    [TestCase("ACαGT", TestName = "TypedNonDna_Unicode_Rejected")]
    public void GlobalAlign_TypedNonDna_RejectedAtConstruction(string invalid)
    {
        var act = () => new DnaSequence(invalid);

        act.Should().Throw<ArgumentException>(
            "the strict typed surface validates the DNA alphabet at construction, so non-DNA never reaches the aligner");
    }

    /// <summary>
    /// MC: a null DnaSequence on the typed overload is the documented validation
    /// gate — it must throw ArgumentNullException (the ThrowIfNull guard,
    /// Global_Alignment_Needleman_Wunsch.md §3.3, §6.1), NOT a NullReferenceException
    /// from a downstream dereference. Pinned in both argument positions.
    /// </summary>
    [Test]
    public void GlobalAlign_TypedNullSequence_ThrowsArgumentNullException()
    {
        var seq = new DnaSequence("ACGT");

        var nullFirst = () => SequenceAligner.GlobalAlign(null!, seq);
        var nullSecond = () => SequenceAligner.GlobalAlign(seq, null!);

        nullFirst.Should().Throw<ArgumentNullException>(
            "a null first DnaSequence hits the documented null gate, not a raw dereference");
        nullSecond.Should().Throw<ArgumentNullException>(
            "a null second DnaSequence hits the documented null gate, not a raw dereference");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: very long sequences (O(mn) DP must complete bounded)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: very long sequences (no OOM, no hang)

    /// <summary>
    /// BE: very long sequences stress the O(mn) DP matrix and traceback. With a few
    /// hundred bases per side the matrix is tens of thousands of cells — modest, but
    /// enough that an off-by-one or unbounded loop would surface. The alignment MUST
    /// complete within the CancelAfter budget (no hang, no OutOfMemory) and STILL
    /// satisfy every structural invariant on the random (fixed-seed) inputs
    /// (Global_Alignment_Needleman_Wunsch.md §4.3 complexity, §6.2). A sequence
    /// aligned against ITSELF must additionally be gap-free with the all-match score.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void GlobalAlign_VeryLongSequences_CompletesAndHoldsInvariants()
    {
        var scoring = SequenceAligner.SimpleDna;
        string a = RandomDna(400);
        string b = RandomDna(400);

        var result = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b), scoring);
        AssertGlobalInvariants(result, a, b, scoring);

        // A long sequence vs itself: the optimal alignment is gap-free, all matches.
        var self = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(a), scoring);
        self.AlignedSequence1.Should().Be(a, "a sequence aligned to itself is gap-free");
        self.AlignedSequence2.Should().Be(a, "a sequence aligned to itself is gap-free");
        self.Score.Should().Be(a.Length * scoring.Match,
            "self-alignment scores Match per base over the full length");
        AssertGlobalInvariants(self, a, a, scoring);
    }

    /// <summary>
    /// BE: the cancellation-aware raw-string overload must honour a pre-cancelled
    /// token promptly on a large input — the periodic ThrowIfCancellationRequested
    /// in the matrix fill must surface as OperationCanceledException rather than
    /// running the full O(mn) DP to completion or hanging
    /// (Global_Alignment_Needleman_Wunsch.md §5.2). CancelAfter bounds the test so a
    /// regression that ignored the token would still fail rather than wedge the run.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void GlobalAlign_CancellationAware_AlreadyCancelled_ThrowsPromptly()
    {
        string a = RandomDna(600);
        string b = RandomDna(600);
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        var act = () => SequenceAligner.GlobalAlign(a, b, SequenceAligner.SimpleDna, cts.Token);

        act.Should().Throw<OperationCanceledException>(
            "the cancellation-aware overload must observe an already-cancelled token, not run to completion");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity: a known pair aligns with the expected score/structure
    // ───────────────────────────────────────────────────────────────────

    #region Positive sanity — the documented GCATGCG vs GATTACA example

    /// <summary>
    /// Positive sanity: the canonical Needleman-Wunsch worked example
    /// GCATGCG vs GATTACA under Match=+1, Mismatch=−1, GapExtend=−1 has optimal
    /// score 0 (Global_Alignment_Needleman_Wunsch.md §7.1; Wikipedia). The exact
    /// aligned strings may vary among optimal solutions, so we assert the score and
    /// the structural invariants rather than a fixed gap layout: equal aligned
    /// length, gap-stripped outputs reproduce the inputs, score equals the
    /// per-column recomputation. This proves the fuzz harness asserts a real,
    /// theory-correct alignment — not merely "did not crash".
    /// </summary>
    [Test]
    public void GlobalAlign_DocExample_GcatgcgVsGattaca_HasScoreZeroAndValidStructure()
    {
        var scoring = SequenceAligner.SimpleDna; // Match=+1, Mismatch=-1, GapExtend=-1
        const string a = "GCATGCG";
        const string b = "GATTACA";

        var result = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.Score.Should().Be(0,
            "the canonical GCATGCG vs GATTACA Needleman-Wunsch example has optimal score 0");
        AssertGlobalInvariants(result, a, b, scoring);
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ALIGN-LOCAL-001 — local alignment (Smith-Waterman) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ALIGN-LOCAL-001 — local alignment (Smith-Waterman)

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty sequences (either side)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty sequences (no positive-scoring region)

    /// <summary>
    /// BE: an empty sequence on EITHER side has no positive-scoring region, so the
    /// zero-floored matrix maximum stays 0, traceback never runs, and the result is a
    /// tagged-Local empty alignment — Score 0, both aligned strings empty, every
    /// coordinate field −1 from the zero-size traceback endpoints
    /// (Local_Alignment_Smith_Waterman.md §3.3, §5.2, §6.1). Never an IndexOutOfRange
    /// on the degenerate single-row / single-column DP. Pinned in both orientations
    /// (empty-vs-empty and empty-vs-length-n). L-INV-01..03 hold trivially.
    /// </summary>
    [TestCase("", "", TestName = "LocalAlign_EmptyVsEmpty")]
    [TestCase("", "ACGT", TestName = "LocalAlign_EmptyVsSeq")]
    [TestCase("ACGT", "", TestName = "LocalAlign_SeqVsEmpty")]
    [TestCase("", "A", TestName = "LocalAlign_EmptyVsSingleBase")]
    public void LocalAlign_TypedEmptySequence_IsTaggedLocalEmptyZeroScore(string a, string b)
    {
        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b));

        result.AlignmentType.Should().Be(AlignmentType.Local,
            "the typed overload flows an empty input through the core DP, tagging the result Local");
        result.Score.Should().Be(0, "an empty side has no positive-scoring region, so the local score is 0");
        result.AlignedSequence1.Should().BeEmpty("no positive-scoring region means an empty local alignment");
        result.AlignedSequence2.Should().BeEmpty("no positive-scoring region means an empty local alignment");
        result.StartPosition1.Should().Be(-1, "the zero-size traceback yields −1 coordinates (§5.2)");
        result.StartPosition2.Should().Be(-1, "the zero-size traceback yields −1 coordinates (§5.2)");
        result.EndPosition1.Should().Be(-1, "the zero-size traceback yields −1 coordinates (§5.2)");
        result.EndPosition2.Should().Be(-1, "the zero-size traceback yields −1 coordinates (§5.2)");
        AssertLocalInvariants(result, a, b);
    }

    /// <summary>
    /// BE: the lenient RAW-STRING surface short-circuits on null OR empty input and
    /// returns AlignmentResult.Empty (Local_Alignment_Smith_Waterman.md §3.3, §6.1) —
    /// the documented strict/lenient split versus the typed overload above (which
    /// traces a tagged-Local empty result). Pinned for empty AND null variants so the
    /// early exit can never drift into a NullReferenceException.
    /// </summary>
    [Test]
    public void LocalAlign_RawStringEmptyOrNull_ReturnsEmptyResult()
    {
        SequenceAligner.LocalAlign("", "").Should().Be(AlignmentResult.Empty,
            "the raw-string overload exits early to AlignmentResult.Empty on empty input");
        SequenceAligner.LocalAlign("ACGT", "").Should().Be(AlignmentResult.Empty,
            "an empty second raw string short-circuits before the DP runs");
        SequenceAligner.LocalAlign("", "ACGT").Should().Be(AlignmentResult.Empty,
            "an empty first raw string short-circuits before the DP runs");

        var nullFirst = () => SequenceAligner.LocalAlign((string)null!, "ACGT");
        var nullSecond = () => SequenceAligner.LocalAlign("ACGT", (string)null!);
        nullFirst.Should().NotThrow("a null raw string is guarded as empty, never dereferenced");
        nullSecond.Should().NotThrow("a null raw string is guarded as empty, never dereferenced");
        nullFirst().Should().Be(AlignmentResult.Empty, "a null raw string yields the empty result");
        nullSecond().Should().Be(AlignmentResult.Empty, "a null raw string yields the empty result");
    }

    /// <summary>
    /// BE: a null DnaSequence on the typed overload is the documented validation gate —
    /// it must throw ArgumentNullException (the ThrowIfNull guard,
    /// Local_Alignment_Smith_Waterman.md §3.3), NOT a NullReferenceException from a
    /// downstream dereference. Pinned in both argument positions.
    /// </summary>
    [Test]
    public void LocalAlign_TypedNullSequence_ThrowsArgumentNullException()
    {
        var seq = new DnaSequence("ACGT");

        var nullFirst = () => SequenceAligner.LocalAlign(null!, seq);
        var nullSecond = () => SequenceAligner.LocalAlign(seq, null!);

        nullFirst.Should().Throw<ArgumentNullException>(
            "a null first DnaSequence hits the documented null gate, not a raw dereference");
        nullSecond.Should().Throw<ArgumentNullException>(
            "a null second DnaSequence hits the documented null gate, not a raw dereference");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: identical sequences (full diagonal optimum)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: identical sequences (whole sequence is the optimum)

    /// <summary>
    /// BE: identical sequences are the case where the local optimum spans the WHOLE
    /// sequence — the full diagonal dominates under positive-match / negative-gap
    /// scoring, so both aligned strings equal the (upper-cased) input, the alignment
    /// is gap-free, and the Score is length × Match
    /// (Local_Alignment_Smith_Waterman.md §6.1 "Identical sequences"). The 0-based
    /// inclusive coordinates therefore span [0, length−1] on both sides. L-INV-01..03
    /// hold (the gap-stripped segment is the whole input, trivially a substring).
    /// </summary>
    [TestCase("A", TestName = "LocalAlign_Identical_SingleBase")]
    [TestCase("ACGT", TestName = "LocalAlign_Identical_Len4")]
    [TestCase("AAAA", TestName = "LocalAlign_Identical_Homopolymer")]
    [TestCase("GATTACAGATTACA", TestName = "LocalAlign_Identical_Len14")]
    public void LocalAlign_TypedIdenticalSequences_AlignWholeSequenceGapFree(string seq)
    {
        var scoring = SequenceAligner.SimpleDna;
        string upper = seq.ToUpperInvariant();

        var result = SequenceAligner.LocalAlign(new DnaSequence(seq), new DnaSequence(seq), scoring);

        result.AlignedSequence1.Should().Be(upper, "identical sequences align as the whole sequence on side 1");
        result.AlignedSequence2.Should().Be(upper, "identical sequences align as the whole sequence on side 2");
        result.AlignedSequence1.Should().NotContain("-", "the full-diagonal optimum is gap-free");
        result.Score.Should().Be(seq.Length * scoring.Match,
            "the whole-sequence diagonal scores Match per base");
        result.StartPosition1.Should().Be(0, "the optimal local region starts at index 0 on side 1");
        result.StartPosition2.Should().Be(0, "the optimal local region starts at index 0 on side 2");
        result.EndPosition1.Should().Be(seq.Length - 1, "the optimal local region ends at the last base on side 1");
        result.EndPosition2.Should().Be(seq.Length - 1, "the optimal local region ends at the last base on side 2");
        AssertLocalInvariants(result, seq, seq);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: completely different sequences (no shared base)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: completely different sequences (zero floor wins)

    /// <summary>
    /// BE: sequences sharing NO base, under negative mismatch and gap penalties, have
    /// no positive-scoring region — every candidate is clamped by the zero floor, the
    /// matrix maximum stays 0, traceback never runs, and the result is Score 0 with an
    /// empty alignment (Local_Alignment_Smith_Waterman.md §6.1 "Completely dissimilar
    /// sequences"). This is the KEY local-vs-global distinction: global would force a
    /// negative end-to-end alignment, but local drops it entirely. L-INV-01..03 hold.
    /// </summary>
    [TestCase("AAAA", "GGGG", TestName = "LocalAlign_Different_AsVsGs")]
    [TestCase("AAAA", "TTTT", TestName = "LocalAlign_Different_AsVsTs")]
    [TestCase("ACAC", "GTGT", TestName = "LocalAlign_Different_DisjointAlphabets")]
    [TestCase("AAAAAAAA", "CCCCGGGG", TestName = "LocalAlign_Different_NoSharedBase")]
    public void LocalAlign_TypedCompletelyDifferent_IsZeroScoreEmptyAlignment(string a, string b)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.AlignmentType.Should().Be(AlignmentType.Local, "the typed overload tags the result Local");
        result.Score.Should().Be(0, "no shared base means no positive-scoring region — the zero floor wins");
        result.AlignedSequence1.Should().BeEmpty("a zero-score local alignment is empty");
        result.AlignedSequence2.Should().BeEmpty("a zero-score local alignment is empty");
        AssertLocalInvariants(result, a, b);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: single char vs single char (match → positive, mismatch → 0)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: single char vs single char

    /// <summary>
    /// BE: the minimal non-empty local alignment — one base vs one base. A MATCH gives
    /// a single positive-scoring cell: Score = Match, the aligned strings are that one
    /// base on each side, coordinates [0,0]. A MISMATCH is clamped by the zero floor:
    /// Score 0, empty alignment (Local_Alignment_Smith_Waterman.md §2.2 zero floor,
    /// §6.1). This pins the defining local behaviour at the smallest scale. L-INV-01..03 hold.
    /// </summary>
    [TestCase("A", "A", TestName = "LocalAlign_SingleChar_Match_A")]
    [TestCase("G", "G", TestName = "LocalAlign_SingleChar_Match_G")]
    public void LocalAlign_TypedSingleCharMatch_IsPositiveOneColumn(string a, string b)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.Score.Should().Be(scoring.Match, "a single matching base scores exactly Match");
        result.AlignedSequence1.Should().Be(a.ToUpperInvariant(), "the one matching base aligns on side 1");
        result.AlignedSequence2.Should().Be(b.ToUpperInvariant(), "the one matching base aligns on side 2");
        result.StartPosition1.Should().Be(0, "the single matching base is at index 0 on side 1");
        result.EndPosition1.Should().Be(0, "the single matching base is at index 0 on side 1");
        result.StartPosition2.Should().Be(0, "the single matching base is at index 0 on side 2");
        result.EndPosition2.Should().Be(0, "the single matching base is at index 0 on side 2");
        AssertLocalInvariants(result, a, b);
    }

    /// <summary>
    /// BE: a 1-char vs 1-char MISMATCH under negative mismatch scoring is clamped to 0
    /// by the zero floor, so the local alignment is empty with Score 0 — never a
    /// negative score, never a forced mismatch column (Local_Alignment_Smith_Waterman.md
    /// §2.2, §6.1). The local-vs-global contrast at the smallest scale: global keeps
    /// the −1 mismatch column, local drops it. L-INV-01..03 hold.
    /// </summary>
    [TestCase("A", "C", TestName = "LocalAlign_SingleChar_Mismatch_AC")]
    [TestCase("G", "T", TestName = "LocalAlign_SingleChar_Mismatch_GT")]
    public void LocalAlign_TypedSingleCharMismatch_IsZeroScoreEmpty(string a, string b)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.Score.Should().Be(0, "a single mismatching base is clamped to 0 by the zero floor");
        result.AlignedSequence1.Should().BeEmpty("a zero-score local alignment is empty");
        result.AlignedSequence2.Should().BeEmpty("a zero-score local alignment is empty");
        AssertLocalInvariants(result, a, b);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity: a known shared substring aligns locally with the expected score
    // ───────────────────────────────────────────────────────────────────

    #region Positive sanity — shared substring + canonical Wikipedia example

    /// <summary>
    /// Positive sanity: two sequences that share a clear common substring align
    /// LOCALLY on exactly that substring. Here the flanks are disjoint and the shared
    /// core is "ACGTACGT": under SimpleDna the optimal local alignment is that core,
    /// gap-free, Score = 8 × Match = 8, with both aligned segments equal to the shared
    /// substring and the coordinates pinpointing it inside each input. This proves the
    /// harness asserts a real, theory-correct LOCAL alignment — the best-scoring
    /// subsequence pair — not merely "did not crash". L-INV-01..03 hold.
    /// </summary>
    [Test]
    public void LocalAlign_SharedSubstring_AlignsOnTheCommonCore()
    {
        var scoring = SequenceAligner.SimpleDna; // Match=+1, Mismatch=-1, GapExtend=-1
        const string core = "ACGTACGT";
        const string a = "TTTTT" + core + "TTTTT"; // core embedded between T-flanks
        const string b = "GGG" + core + "GGGGG";    // same core, disjoint G-flanks

        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.AlignmentType.Should().Be(AlignmentType.Local, "a successful local alignment is tagged Local");
        result.Score.Should().Be(core.Length * scoring.Match,
            "the optimal local alignment is the shared core, scoring Match per base");
        new string(result.AlignedSequence1.Where(c => c != '-').ToArray())
            .Should().Be(core, "the gap-stripped aligned segment 1 is exactly the shared core");
        new string(result.AlignedSequence2.Where(c => c != '-').ToArray())
            .Should().Be(core, "the gap-stripped aligned segment 2 is exactly the shared core");
        result.AlignedSequence1.Should().NotContain("-", "the shared-core local alignment is gap-free");
        AssertLocalInvariants(result, a, b);
    }

    /// <summary>
    /// Positive sanity: the canonical Smith-Waterman Wikipedia example — TGTTACGG vs
    /// GGTTGACTA under Match=+3, Mismatch=−3, GapExtend=−2 — has optimal local
    /// alignment GTT-AC / GTTGAC with Score 13
    /// (Local_Alignment_Smith_Waterman.md §6.1, §7.1; validated independently in the
    /// unit-test suite). Pinned here as a cross-check that the fuzz harness agrees with
    /// the documented worked example, including a single internal gap column.
    /// </summary>
    [Test]
    public void LocalAlign_WikipediaExample_HasScore13AndExpectedAlignment()
    {
        var scoring = new ScoringMatrix(Match: 3, Mismatch: -3, GapOpen: -2, GapExtend: -2);
        const string a = "TGTTACGG";
        const string b = "GGTTGACTA";

        var result = SequenceAligner.LocalAlign(new DnaSequence(a), new DnaSequence(b), scoring);

        result.Score.Should().Be(13, "the canonical Smith-Waterman example has optimal local score 13");
        result.AlignedSequence1.Should().Be("GTT-AC", "the documented worked-example alignment for side 1");
        result.AlignedSequence2.Should().Be("GTTGAC", "the documented worked-example alignment for side 2");
        AssertLocalInvariants(result, a, b);
    }

    #endregion

    #endregion
}
