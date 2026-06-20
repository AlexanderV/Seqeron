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
/// Fuzz tests for the Alignment area — pairwise GLOBAL alignment (Needleman-Wunsch).
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
}
