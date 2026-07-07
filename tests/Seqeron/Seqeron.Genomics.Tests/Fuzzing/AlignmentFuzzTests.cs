namespace Seqeron.Genomics.Tests.Fuzzing;

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
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ALIGN-SEMI-001 — semi-global (fitting / query-in-reference) alignment (Alignment)
/// Checklist: docs/checklists/03_FUZZING.md, row 37.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the boundary fodder for the fitting DP: empty
///          seqs (query and/or reference), NO overlap (disjoint alphabets), COMPLETE
///          overlap (identical, or the query embedded as a substring of the
///          reference), and the 1-char vs 1-char extreme. These probe the degenerate
///          single-row / single-column DP and the trailing-suffix traceback branch.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The semi-global (fitting / query-in-reference) contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The repository implements ONE member of the semi-global family: the FITTING
/// (query-in-reference) variant, NOT overlap and NOT all-ends-free. The query
/// (sequence1) is aligned GLOBALLY end-to-end; only the REFERENCE (sequence2) ends are
/// FREE. The recurrence is the standard linear-gap Needleman-Wunsch (d = GapExtend; the
/// GapOpen field is ignored — linear, not affine) with a FITTING initialization:
///   F(0,j) = 0            ← free LEADING reference gaps (any reference prefix is free),
///   F(i,0) = d·i          ← the query first column IS penalized (the full query must be aligned),
///   score = max_j F(m,j)  ← free TRAILING reference gaps (the score is the LAST-ROW maximum).
/// There is NO Smith-Waterman zero floor, so the fitting score MAY be negative.
///   — docs/algorithms/Alignment/Semi_Global_Alignment.md §2.2, §2.4, §4.1, §5.2.
///     Sources: Sequence alignment; Needleman-Wunsch; Rosalind SIMS/SMGB; Brudno 2003 (glocal).
///
/// KEY SEMI-GLOBAL invariants asserted across the fuzz fodder (Semi_Global_Alignment.md §2.4):
///   • S-INV-01 — the two aligned output strings have EQUAL length (S §2.4 INV-02).
///   • S-INV-02 — gap-stripping AlignedSequence1 reproduces the FULL (upper-cased) query:
///     the query is aligned globally and never drops bases (S §2.4 INV-01).
///   • S-INV-03 — gap-stripping AlignedSequence2 reproduces the FULL (upper-cased)
///     reference: the implementation preserves the unmatched reference prefix and suffix
///     explicitly, so the entire reference is represented (S §5.2). [This is the invariant
///     that surfaced the trailing-suffix bug fixed in this campaign — see below.]
///   • S-INV-04 — the reported Score equals the per-column linear-model sum over the
///     FITTED BODY only, i.e. excluding the leading and trailing maximal runs of FREE
///     reference end-gap columns (columns where AlignedSequence1 == '-'). This is the
///     fitting score = max_j F(m,j) recomputed independently (S §2.4 INV-03, §2.2).
///   • S-INV-05 — the result is tagged AlignmentType.SemiGlobal; StartPosition* = 0,
///     EndPosition1 = query.Length−1, EndPosition2 = reference.Length−1 (S §3.2).
///
/// SOURCE BUG FOUND & FIXED during this campaign (SequenceAligner.Traceback, semi-global
/// trailing-suffix branch): the unmatched trailing reference suffix seq2[j..n−1] was
/// appended in FORWARD order and then reversed together with the rest of chars2, so the
/// suffix came out REVERSED — AlignedSequence2 became an anagram of the reference (e.g.
/// query "ATGC" vs reference "ATGCTTAGG" produced "ATGCGGATT"; an empty query vs "ACGT"
/// produced "TGCA"). The Score was unaffected (it is the last-row maximum), but S-INV-03
/// was violated: the alignment misrepresented the reference content. The bug was masked
/// in earlier hand-checked cases because the suffixes there were homopolymers/palindromic.
/// FIX: append the trailing suffix in REVERSE order so the final reverse restores forward
/// reference order. The cases below (distinct, non-palindromic suffixes) pin the fix.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ALIGN-MULTI-001 — multiple sequence alignment (Alignment)
/// Checklist: docs/checklists/03_FUZZING.md, row 38 (the LAST Alignment-area unit).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the size boundaries of the COLLECTION input to the
///          MSA: an empty list, a single sequence, exactly two sequences, 100+
///          sequences (the heuristic must complete in bounded time — kept SHORT), and
///          a collection of length-1 sequences. These probe the degenerate trivial
///          paths, the center-selection + anchor reconciliation heuristic at scale, and
///          the zero-/one-length sequence corners of the gap-merge/padding stage.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The multiple-sequence-alignment (center-star) contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// `SequenceAligner.MultipleAlign(IEnumerable&lt;DnaSequence&gt;, ScoringMatrix?)` aligns
/// N validated DNA sequences into rows of EQUAL length by inserting gap characters,
/// using a simplified CENTER-STAR heuristic: select one center by 4-mer cosine
/// similarity, align every other sequence to it with the anchor-based pairwise aligner,
/// reconcile the independent gap patterns into one coordinate system, pad all rows to
/// equal length, build a majority-vote consensus, and report a sum-of-pairs (SP) score.
/// Returns a MultipleAlignmentResult(AlignedSequences[], Consensus, TotalScore).
///   — docs/algorithms/Alignment/Multiple_Sequence_Alignment.md §1, §2.2, §4.1, §5.2.
///     Sources: Multiple sequence alignment; Clustal; Consensus sequence (Wikipedia).
///
/// KEY STRUCTURAL INVARIANTS asserted across the fuzz fodder (the MSA contract that must
/// never silently break — Multiple_Sequence_Alignment.md §2.4):
///   • M-INV-01 — ALL output rows have the SAME length (the defining MSA representation;
///     gaps pad every row into one coordinate system). This is the KEY invariant.
///   • M-INV-02 — removing the gap character '-' from each output row reproduces the
///     corresponding (upper-cased) input sequence exactly (the algorithm only inserts
///     gaps; it never rewrites or reorders a base — content + order preserved).
///   • M-INV-03 — no output column consists ENTIRELY of gaps (a structural MSA validity
///     condition: an all-gap column carries no information and is excluded).
///   • M-INV-04 — Consensus.Length equals the aligned row length (the consensus is built
///     one column at a time across the final rows).
///   • M-INV-COUNT — the number of output rows equals the number of input sequences
///     (count preservation; no sequence is dropped or duplicated).
///
/// DOCUMENTED BOUNDARY behaviour (pinned so the trivial-case fast paths cannot drift —
/// Multiple_Sequence_Alignment.md §3.3, §6.1):
///   • Null collection → ArgumentNullException (the ThrowIfNull guard), NOT a downstream
///     NullReferenceException while enumerating.
///   • Empty collection → MultipleAlignmentResult.Empty (empty rows, empty consensus,
///     TotalScore 0) — an explicit early return, never an IndexOutOfRange on Max() over
///     an empty row set.
///   • Single sequence → that one row UNCHANGED as both the sole aligned row and the
///     consensus, with TotalScore = 0 (no pairs to score) — no gaps are introduced.
///   • Two sequences → a well-formed 2-row MSA consistent with the pairwise path; equal
///     row lengths, gap-stripping recovers both inputs.
///   • 100+ sequences (kept SHORT, length 5–20) → the center-star heuristic COMPLETES
///     within a CancelAfter budget (no hang, no OutOfMemory on the O(k²·L) center
///     selection + parallel reconciliation) and still satisfies every M-INV.
///   • Length-1 sequences / empty sequences inside the collection → handled without
///     rejecting the whole alignment; rows stay equal length and gap-strip back to their
///     (possibly empty) inputs.
/// All inputs here are validated DNA (the strict DnaSequence surface); randomness (where
/// used) is from the fixture's locally fixed-seed Rng so the fuzz fodder is reproducible
/// and adding tests cannot perturb any other fixture.
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

    /// <summary>
    /// Asserts the five semi-global (FITTING / query-in-reference) structural invariants
    /// against the raw (upper-cased) query/reference the aligner actually consumed:
    /// equal aligned length (S-INV-01); gap-stripped AlignedSequence1 reproduces the FULL
    /// query (S-INV-02, the global-query property); gap-stripped AlignedSequence2
    /// reproduces the FULL reference (S-INV-03, the prefix/suffix preservation property);
    /// the reported Score equals the per-column linear-model sum over the FITTED BODY
    /// only — excluding the leading and trailing maximal runs of FREE reference end-gap
    /// columns (S-INV-04, fitting score = max_j F(m,j)); and the result is tagged
    /// SemiGlobal with the documented coordinate bounds (S-INV-05).
    /// — docs/algorithms/Alignment/Semi_Global_Alignment.md §2.2, §2.4, §3.2, §5.2.
    /// </summary>
    private static void AssertSemiGlobalInvariants(
        AlignmentResult result, string query, string reference, ScoringMatrix scoring)
    {
        string expectedQuery = query.ToUpperInvariant();
        string expectedRef = reference.ToUpperInvariant();
        string a1 = result.AlignedSequence1;
        string a2 = result.AlignedSequence2;

        // S-INV-05: tag + documented coordinate bounds (full-input, not fitted interval).
        result.AlignmentType.Should().Be(AlignmentType.SemiGlobal,
            "S-INV-05: a semi-global (fitting) alignment is tagged SemiGlobal");
        result.StartPosition1.Should().Be(0, "S-INV-05: the implementation returns StartPosition1 = 0");
        result.StartPosition2.Should().Be(0, "S-INV-05: the implementation returns StartPosition2 = 0");
        result.EndPosition1.Should().Be(query.Length - 1,
            "S-INV-05: EndPosition1 is query.Length − 1 (−1 for an empty query)");
        result.EndPosition2.Should().Be(reference.Length - 1,
            "S-INV-05: EndPosition2 is reference.Length − 1 (−1 for an empty reference)");

        // S-INV-01: the two aligned strings share one coordinate system.
        a1.Length.Should().Be(a2.Length,
            "S-INV-01: semi-global alignment pads query and reference to equal length");

        // S-INV-02: the query is aligned globally — gap-stripping reproduces it in full.
        new string(a1.Where(c => c != '-').ToArray())
            .Should().Be(expectedQuery,
                "S-INV-02: gap-stripped AlignedSequence1 reproduces the FULL query (the fitting variant aligns the whole query)");

        // S-INV-03: the full reference is preserved (prefix + body + suffix). This is the
        // invariant the trailing-suffix bug violated — the reference must NOT be reordered.
        new string(a2.Where(c => c != '-').ToArray())
            .Should().Be(expectedRef,
                "S-INV-03: gap-stripped AlignedSequence2 reproduces the FULL reference (unmatched prefix/suffix preserved, never reversed)");

        // S-INV-04: recompute the fitting score over the FITTED BODY only. The leading and
        // trailing maximal runs of columns where the QUERY is a gap are the FREE reference
        // end-gaps (unpenalized); everything in between is scored by the linear model.
        int lead = 0;
        while (lead < a1.Length && a1[lead] == '-') lead++;
        int trail = a1.Length;
        while (trail > lead && a1[trail - 1] == '-') trail--;

        int recomputed = 0;
        for (int k = lead; k < trail; k++)
        {
            char qc = a1[k];
            char rc = a2[k];
            if (qc == '-' || rc == '-')
                recomputed += scoring.GapExtend;
            else
                recomputed += qc == rc ? scoring.Match : scoring.Mismatch;
        }
        result.Score.Should().Be(recomputed,
            "S-INV-04: the fitting score equals the per-column linear-model sum over the fitted body, with leading/trailing reference end-gaps free");
    }

    /// <summary>
    /// Asserts the multiple-sequence-alignment structural invariants against the raw
    /// (upper-cased) input sequences the aligner actually consumed: count preservation
    /// (one output row per input, M-INV-COUNT); ALL rows share one length (M-INV-01, the
    /// KEY MSA invariant); each gap-stripped row reproduces the corresponding input
    /// (M-INV-02, content + order preserved); no column is entirely gaps (M-INV-03); and
    /// Consensus.Length equals the aligned row length (M-INV-04). Used by every non-empty
    /// MSA fuzz case. — docs/algorithms/Alignment/Multiple_Sequence_Alignment.md §2.4.
    /// </summary>
    private static void AssertMsaInvariants(MultipleAlignmentResult result, string[] inputs)
    {
        string[] rows = result.AlignedSequences;

        // M-INV-COUNT: one output row per input, no sequence dropped or duplicated.
        rows.Length.Should().Be(inputs.Length,
            "M-INV-COUNT: the MSA emits exactly one aligned row per input sequence");

        // M-INV-01 (KEY): every aligned row has the SAME length — one coordinate system.
        // (Vacuous for an empty MSA, which has no rows.)
        if (rows.Length > 0)
            rows.Select(r => r.Length).Distinct().Should().HaveCount(1,
                "M-INV-01: all aligned rows of an MSA share one common length (the defining representation)");
        int length = rows.Length == 0 ? 0 : rows[0].Length;

        // M-INV-02: gap-stripping each row reproduces its (upper-cased) input exactly.
        for (int i = 0; i < rows.Length; i++)
        {
            string expected = inputs[i].ToUpperInvariant();
            new string(rows[i].Where(c => c != '-').ToArray())
                .Should().Be(expected,
                    "M-INV-02: gap-stripping aligned row {0} reproduces input {0} (gaps inserted, bases never rewritten/reordered)", i);
        }

        // M-INV-03: no output column consists entirely of gaps (MSA validity condition).
        for (int col = 0; col < length; col++)
        {
            bool allGaps = rows.All(r => col >= r.Length || r[col] == '-');
            allGaps.Should().BeFalse(
                "M-INV-03: column {0} must not consist entirely of gaps", col);
        }

        // M-INV-04: the consensus is built one column per aligned row position.
        result.Consensus.Length.Should().Be(length,
            "M-INV-04: Consensus.Length equals the aligned row length");
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

    // ═══════════════════════════════════════════════════════════════════
    //  ALIGN-SEMI-001 — semi-global (fitting / query-in-reference) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ALIGN-SEMI-001 — semi-global (overlap) alignment

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty sequences (query and/or reference)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty sequences

    /// <summary>
    /// BE: empty query vs empty reference is the lowest size boundary. The 1×1 DP
    /// degenerates to F(0,0)=0, the last-row maximum is 0, and traceback emits nothing —
    /// a trivial empty alignment, Score 0, tagged SemiGlobal, never an IndexOutOfRange
    /// on the degenerate matrix (Semi_Global_Alignment.md §2.2, §6.1). Coordinates are
    /// −1 on both ends (Length−1 = −1). S-INV-01..05 hold trivially.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_EmptyVsEmpty_IsTrivialZeroScoreAlignment()
    {
        var scoring = SequenceAligner.SimpleDna;
        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(""), new DnaSequence(""), scoring);

        result.AlignedSequence1.Should().BeEmpty("an empty query/reference pair yields empty aligned strings");
        result.AlignedSequence2.Should().BeEmpty("an empty query/reference pair yields empty aligned strings");
        result.Score.Should().Be(0, "there are no columns to score, so the fitting score is 0");
        AssertSemiGlobalInvariants(result, "", "", scoring);
    }

    /// <summary>
    /// BE: an EMPTY QUERY against a length-n reference. The query is aligned globally
    /// (it contributes nothing), and the ENTIRE reference is a free end-gap overhang
    /// under the fitting initialization (F(0,j)=0, free leading reference gaps), so the
    /// fitting score is 0 and the reference is fully preserved against an all-gap query
    /// (Semi_Global_Alignment.md §2.2, §5.2). This is the case that originally exposed
    /// the trailing-suffix REVERSAL bug — an empty query vs "ACGT" used to return the
    /// anagram "TGCA"; S-INV-03 now pins the reference order. No crash on the single-row
    /// DP. S-INV-01..05 hold.
    /// </summary>
    [TestCase("ACGT", TestName = "SemiGlobalAlign_EmptyQuery_Len4")]
    [TestCase("A", TestName = "SemiGlobalAlign_EmptyQuery_SingleBase")]
    [TestCase("GACATGCTTAG", TestName = "SemiGlobalAlign_EmptyQuery_DistinctSuffix")]
    public void SemiGlobalAlign_EmptyQueryVsReference_IsZeroScoreFreeOverhang(string reference)
    {
        var scoring = SequenceAligner.SimpleDna;
        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(""), new DnaSequence(reference), scoring);

        result.Score.Should().Be(0,
            "the whole reference is a free end-gap overhang for an empty query, so the fitting score is 0");
        result.AlignedSequence1.Should().Be(new string('-', reference.Length),
            "the empty query aligns as n gap characters across the reference");
        new string(result.AlignedSequence2.Where(c => c != '-').ToArray())
            .Should().Be(reference.ToUpperInvariant(),
                "the reference is preserved in forward order against the all-gap query (no suffix reversal)");
        AssertSemiGlobalInvariants(result, "", reference, scoring);
    }

    /// <summary>
    /// BE: a length-m query against an EMPTY reference. There is no reference to fit
    /// into, so the query first column is penalized in full (F(i,0)=d·i, the query is
    /// always aligned globally) — the fitting score is m·d (here −m) with the query
    /// aligned against m gaps and NO free end-gaps available (the reference is empty)
    /// (Semi_Global_Alignment.md §2.2, §6.1 "Query longer than the reference"). S-INV-01..05 hold.
    /// </summary>
    [TestCase("ACGT", TestName = "SemiGlobalAlign_EmptyReference_Len4")]
    [TestCase("A", TestName = "SemiGlobalAlign_EmptyReference_SingleBase")]
    public void SemiGlobalAlign_QueryVsEmptyReference_IsAllGapPenalizedQuery(string query)
    {
        var scoring = SequenceAligner.SimpleDna;
        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(""), scoring);

        result.Score.Should().Be(query.Length * scoring.GapExtend,
            "with no reference the full query is penalized at d per base (F(i,0) = d·i)");
        result.AlignedSequence1.Should().Be(query.ToUpperInvariant(),
            "the query is aligned in full against the empty reference");
        result.AlignedSequence2.Should().Be(new string('-', query.Length),
            "the empty reference aligns as m gap characters");
        AssertSemiGlobalInvariants(result, query, "", scoring);
    }

    /// <summary>
    /// BE: a null DnaSequence on the typed overload is the documented validation gate —
    /// it must throw ArgumentNullException (the ThrowIfNull guard,
    /// Semi_Global_Alignment.md §3.3, §6.1), NOT a NullReferenceException from a
    /// downstream dereference. Pinned in both argument positions.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_TypedNullSequence_ThrowsArgumentNullException()
    {
        var seq = new DnaSequence("ACGT");

        var nullFirst = () => SequenceAligner.SemiGlobalAlign(null!, seq);
        var nullSecond = () => SequenceAligner.SemiGlobalAlign(seq, null!);

        nullFirst.Should().Throw<ArgumentNullException>(
            "a null query hits the documented null gate, not a raw dereference");
        nullSecond.Should().Throw<ArgumentNullException>(
            "a null reference hits the documented null gate, not a raw dereference");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: no overlap (disjoint sequences, no shared base)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: no overlap (disjoint alphabets)

    /// <summary>
    /// BE: a query and reference that share NO base. The fitting DP still aligns the
    /// full query globally; the best placement carries the query as m gaps in the
    /// reference (or all mismatches) and leaves the rest of the reference as free
    /// end-gaps — there is NO zero floor, so the score is NEGATIVE (S-INV-04), never
    /// clamped to 0 the way Smith-Waterman would (Semi_Global_Alignment.md §2.2, §2.4
    /// INV-04). The key check: no crash, the full query and full reference are both
    /// preserved, and the negative fitting score is the body-only recomputation.
    /// S-INV-01..05 hold.
    /// </summary>
    [TestCase("AAAA", "GGGG", TestName = "SemiGlobalAlign_NoOverlap_AsVsGs")]
    [TestCase("ACAC", "GTGT", TestName = "SemiGlobalAlign_NoOverlap_DisjointAlphabets")]
    [TestCase("AAA", "GGGGGG", TestName = "SemiGlobalAlign_NoOverlap_ShortVsLong")]
    [TestCase("A", "GGGG", TestName = "SemiGlobalAlign_NoOverlap_SingleQueryNoMatch")]
    public void SemiGlobalAlign_NoOverlap_IsNegativeScoreFullContentPreserved(string query, string reference)
    {
        var scoring = SequenceAligner.SimpleDna;

        var act = () => SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);
        act.Should().NotThrow("disjoint sequences must not crash the fitting DP");

        var result = act();
        result.Score.Should().BeLessThan(0,
            "with no shared base and no zero floor the fitting score is negative (S-INV-04, §2.4 INV-04)");
        AssertSemiGlobalInvariants(result, query, reference, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: complete overlap (identical, or query embedded in reference)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: complete overlap (identical / query is a substring of reference)

    /// <summary>
    /// BE: IDENTICAL query and reference — the maximal complete-overlap case. The full
    /// diagonal dominates, both aligned strings equal the (upper-cased) input gap-free,
    /// and the fitting score is length × Match (every base matched, no free end-gaps
    /// needed because the reference IS the query) (Semi_Global_Alignment.md §6.1). The
    /// 0-based coordinate bounds span the whole input. S-INV-01..05 hold.
    /// </summary>
    [TestCase("A", TestName = "SemiGlobalAlign_Identical_SingleBase")]
    [TestCase("ACGT", TestName = "SemiGlobalAlign_Identical_Len4")]
    [TestCase("GATTACAGATTACA", TestName = "SemiGlobalAlign_Identical_Len14")]
    public void SemiGlobalAlign_IdenticalSequences_FullMatchGapFree(string seq)
    {
        var scoring = SequenceAligner.SimpleDna;
        string upper = seq.ToUpperInvariant();

        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(seq), new DnaSequence(seq), scoring);

        result.AlignedSequence1.Should().Be(upper, "identical sequences align as the whole sequence on the query side");
        result.AlignedSequence2.Should().Be(upper, "identical sequences align as the whole sequence on the reference side");
        result.AlignedSequence1.Should().NotContain("-", "the full-diagonal complete overlap is gap-free");
        result.Score.Should().Be(seq.Length * scoring.Match,
            "an identical query/reference scores Match per base over the full length");
        AssertSemiGlobalInvariants(result, seq, seq, scoring);
    }

    /// <summary>
    /// BE: COMPLETE overlap where the query is a SUBSTRING of a longer reference (the
    /// query "fits" inside the reference). The fitting score is the embedded match score
    /// query.Length × Match; the unmatched reference PREFIX and SUFFIX are free end-gaps
    /// aligned against query gaps; the full reference (prefix + matched core + suffix) is
    /// preserved IN ORDER (Semi_Global_Alignment.md §6.1 "Query embedded exactly inside
    /// the reference"). The distinct, non-palindromic flanks here would expose the
    /// trailing-suffix reversal bug if it regressed (S-INV-03). S-INV-01..05 hold.
    /// </summary>
    [TestCase("ATGC", "GACATGCTTAG", TestName = "SemiGlobalAlign_Embedded_DistinctFlanks")]
    [TestCase("ATGC", "ATGCTTAGG", TestName = "SemiGlobalAlign_Embedded_PrefixMatch_DistinctSuffix")]
    [TestCase("ATGC", "GGATTATGC", TestName = "SemiGlobalAlign_Embedded_SuffixMatch_DistinctPrefix")]
    [TestCase("ACGTACGT", "TTTTTACGTACGTGGGGG", TestName = "SemiGlobalAlign_Embedded_LongCore")]
    public void SemiGlobalAlign_QueryEmbeddedInReference_FitsOnTheMatchedCore(string query, string reference)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

        result.Score.Should().Be(query.Length * scoring.Match,
            "an exactly embedded query fits with Match per query base (free reference flanks)");
        new string(result.AlignedSequence1.Where(c => c != '-').ToArray())
            .Should().Be(query.ToUpperInvariant(), "the gap-stripped query side is the full query");
        new string(result.AlignedSequence2.Where(c => c != '-').ToArray())
            .Should().Be(reference.ToUpperInvariant(),
                "the full reference is preserved in forward order (prefix + core + suffix, no reversal)");
        AssertSemiGlobalInvariants(result, query, reference, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: single char vs single char
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: single char vs single char (match / mismatch)

    /// <summary>
    /// BE: the minimal non-empty fitting alignment — one query base vs one reference
    /// base. A MATCH is a single diagonal column, Score = Match, gap-free, coordinates
    /// [0,0]. A MISMATCH has no zero floor and no free interior, so the query is forced
    /// against the reference: under SimpleDna the optimal one-base fitting is the single
    /// mismatch/gap arrangement scoring Mismatch (−1), NOT clamped to 0 like local
    /// alignment (Semi_Global_Alignment.md §2.2, §2.4 INV-04). S-INV-01..05 (and the
    /// body-only S-INV-04 recompute) hold for both.
    /// </summary>
    [TestCase("A", "A", 1, TestName = "SemiGlobalAlign_SingleChar_Match_A")]
    [TestCase("G", "G", 1, TestName = "SemiGlobalAlign_SingleChar_Match_G")]
    [TestCase("A", "C", -1, TestName = "SemiGlobalAlign_SingleChar_Mismatch_AC")]
    [TestCase("G", "T", -1, TestName = "SemiGlobalAlign_SingleChar_Mismatch_GT")]
    public void SemiGlobalAlign_SingleCharVsSingleChar_IsDefined(string query, string reference, int expectedScore)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

        result.Score.Should().Be(expectedScore,
            "a single matching base scores Match; a single mismatching base scores Mismatch (no zero floor)");
        AssertSemiGlobalInvariants(result, query, reference, scoring);
    }

    /// <summary>
    /// BE: a single query base that occurs INSIDE a longer reference is the smallest
    /// fitting case — the one base matches at its position and the rest of the reference
    /// is split into free leading/trailing end-gaps, so Score = Match (1) regardless of
    /// the reference length (Semi_Global_Alignment.md §2.2, §6.1). S-INV-01..05 hold.
    /// </summary>
    [TestCase("A", "GGGAGGG", TestName = "SemiGlobalAlign_SingleQueryInReference_Interior")]
    [TestCase("A", "AGGGGGG", TestName = "SemiGlobalAlign_SingleQueryInReference_Prefix")]
    [TestCase("A", "GGGGGGA", TestName = "SemiGlobalAlign_SingleQueryInReference_Suffix")]
    public void SemiGlobalAlign_SingleQueryBaseInReference_FitsAtMatch(string query, string reference)
    {
        var scoring = SequenceAligner.SimpleDna;

        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

        result.Score.Should().Be(scoring.Match,
            "a single query base present in the reference fits at Match, the flanks being free end-gaps");
        AssertSemiGlobalInvariants(result, query, reference, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity: a short query fits inside a long reference on its matched core
    // ───────────────────────────────────────────────────────────────────

    #region Positive sanity — short query fits inside a long reference

    /// <summary>
    /// Positive sanity: the documented worked example — query "ATGC" against reference
    /// "AAAATGCAAA" — fits exactly at reference positions 4–7 with Score 4 (4 matches ×
    /// +1), the unmatched A-flanks being free end-gaps (Semi_Global_Alignment.md §6.1;
    /// validated independently in SequenceAligner_SemiGlobalAlign_Tests). This proves the
    /// fuzz harness asserts a real, theory-correct FITTING alignment — the query placed at
    /// its best-scoring position with free reference overhangs — not merely "did not
    /// crash". S-INV-01..05 hold.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_DocExample_AtgcInAaaatgcaaa_FitsWithScore4()
    {
        var scoring = SequenceAligner.SimpleDna; // Match=+1, Mismatch=-1, GapExtend=-1
        const string query = "ATGC";
        const string reference = "AAAATGCAAA";

        var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

        result.Score.Should().Be(4,
            "the query ATGC fits perfectly inside AAAATGCAAA, scoring 4 matches with free flanks");
        new string(result.AlignedSequence1.Where(c => c != '-').ToArray())
            .Should().Be(query, "the gap-stripped query side is exactly the query");
        new string(result.AlignedSequence2.Where(c => c != '-').ToArray())
            .Should().Be(reference, "the gap-stripped reference side is the full reference in order");
        AssertSemiGlobalInvariants(result, query, reference, scoring);
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ALIGN-MULTI-001 — multiple sequence alignment (center-star) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ALIGN-MULTI-001 — multiple sequence alignment

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty list
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty list (trivial Empty result)

    /// <summary>
    /// BE: an empty collection is the lowest size boundary of the MSA input. It must
    /// return MultipleAlignmentResult.Empty via the explicit early return — empty rows,
    /// empty consensus, TotalScore 0 — never an IndexOutOfRange from a Max() over an empty
    /// row set or a NullReference while enumerating
    /// (Multiple_Sequence_Alignment.md §3.3, §6.1). M-INV-* hold vacuously (zero rows).
    /// </summary>
    [Test]
    public void MultipleAlign_EmptyList_ReturnsEmptyResult()
    {
        var act = () => SequenceAligner.MultipleAlign(Array.Empty<DnaSequence>());

        act.Should().NotThrow("an empty collection is a defined boundary, not a crash");
        var result = act();

        result.Should().Be(MultipleAlignmentResult.Empty,
            "an empty collection returns the Empty MSA sentinel (explicit early return)");
        result.AlignedSequences.Should().BeEmpty("there are no sequences to align");
        result.Consensus.Should().BeEmpty("an empty MSA has no consensus");
        result.TotalScore.Should().Be(0, "an empty MSA has no columns to score");
        AssertMsaInvariants(result, Array.Empty<string>());
    }

    /// <summary>
    /// BE: a null collection is the documented validation gate — it must throw
    /// ArgumentNullException (the ThrowIfNull guard, Multiple_Sequence_Alignment.md §6.1),
    /// NOT a NullReferenceException from a downstream enumeration.
    /// </summary>
    [Test]
    public void MultipleAlign_NullCollection_ThrowsArgumentNullException()
    {
        var act = () => SequenceAligner.MultipleAlign(null!);

        act.Should().Throw<ArgumentNullException>(
            "a null collection hits the documented null gate, not a raw dereference");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: single sequence (trivial unchanged row)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: single sequence (trivially that one row)

    /// <summary>
    /// BE: a single-sequence collection is the trivial MSA — the algorithm short-circuits
    /// before any center-selection/anchor work and returns that one sequence UNCHANGED as
    /// both the sole aligned row and the consensus, with TotalScore 0 (no pairs to score),
    /// no gaps introduced (Multiple_Sequence_Alignment.md §3.3, §6.1). Pinned across a
    /// short, a long, and a homopolymer single input. M-INV-* hold.
    /// </summary>
    [TestCase("A", TestName = "MultipleAlign_Single_SingleBase")]
    [TestCase("ATGCATGC", TestName = "MultipleAlign_Single_Len8")]
    [TestCase("AAAAAAAA", TestName = "MultipleAlign_Single_Homopolymer")]
    public void MultipleAlign_SingleSequence_IsTriviallyThatRowUnchanged(string seq)
    {
        string upper = seq.ToUpperInvariant();

        var result = SequenceAligner.MultipleAlign(new[] { new DnaSequence(seq) });

        result.AlignedSequences.Should().HaveCount(1, "a single input yields a single aligned row");
        result.AlignedSequences[0].Should().Be(upper, "the single row is the input unchanged (no gaps needed)");
        result.AlignedSequences[0].Should().NotContain("-", "a single sequence needs no gaps");
        result.Consensus.Should().Be(upper, "the consensus of one sequence is that sequence");
        result.TotalScore.Should().Be(0, "a single sequence has no pairs to score");
        AssertMsaInvariants(result, new[] { seq });
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: 2 sequences (minimal non-trivial MSA)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: two sequences (minimal non-trivial MSA)

    /// <summary>
    /// BE: exactly two sequences is the minimal non-trivial MSA — equal row lengths,
    /// gap-stripping recovers both inputs, count preserved, no all-gap column. Pinned for
    /// identical inputs (gap-free, both rows equal), different-length inputs (the shorter
    /// padded with gaps), and partially overlapping inputs
    /// (Multiple_Sequence_Alignment.md §6.1). M-INV-* hold across all three.
    /// </summary>
    [TestCase("ATGC", "ATGC", TestName = "MultipleAlign_Two_Identical")]
    [TestCase("ATGCATGC", "ATGC", TestName = "MultipleAlign_Two_DifferentLengths")]
    [TestCase("ATGCATGC", "GCATGCAT", TestName = "MultipleAlign_Two_PartiallyOverlapping")]
    [TestCase("AAAA", "TTTT", TestName = "MultipleAlign_Two_NoSharedBase")]
    public void MultipleAlign_TwoSequences_IsWellFormedMsa(string a, string b)
    {
        var result = SequenceAligner.MultipleAlign(new[] { new DnaSequence(a), new DnaSequence(b) });

        result.AlignedSequences.Should().HaveCount(2, "two inputs yield two aligned rows");
        AssertMsaInvariants(result, new[] { a, b });

        // Identical inputs align to identical, gap-free rows.
        if (a == b)
        {
            result.AlignedSequences[0].Should().Be(result.AlignedSequences[1],
                "identical inputs produce identical aligned rows");
            result.AlignedSequences[0].Should().NotContain("-",
                "identical inputs need no gaps");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: 100+ sequences (heuristic must complete bounded)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: 100+ short sequences (no hang, no OOM)

    /// <summary>
    /// BE: 100+ sequences stress the O(k²·L) center selection and the parallel anchor
    /// reconciliation. Kept SHORT (length 5–20) so the heuristic completes well within the
    /// CancelAfter budget (no hang, no OutOfMemory) and still satisfies every M-INV on the
    /// random fixed-seed inputs (Multiple_Sequence_Alignment.md §4.3 complexity, §5.2).
    /// Exercised on diverse random sequences AND on a 100-row all-identical set (the case
    /// where the reconciliation must keep every row gap-free and equal).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void MultipleAlign_ManyShortSequences_CompletesAndHoldsInvariants()
    {
        // 120 diverse short sequences (length 5–20) from the fixture's fixed-seed Rng.
        var inputs = Enumerable.Range(0, 120)
            .Select(_ => RandomDna(Rng.Next(5, 21)))
            .ToArray();
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();

        var act = () => SequenceAligner.MultipleAlign(seqs);
        act.Should().NotThrow("100+ short sequences must align without crashing");

        var result = act();
        result.AlignedSequences.Should().HaveCount(120, "count is preserved at scale");
        AssertMsaInvariants(result, inputs);

        // 100 identical short sequences: every row equals the input, gap-free.
        const string same = "ACGTACGT";
        var identical = Enumerable.Range(0, 100).Select(_ => new DnaSequence(same)).ToArray();
        var identicalResult = SequenceAligner.MultipleAlign(identical);
        identicalResult.AlignedSequences.Should().HaveCount(100, "count preserved for the identical set");
        identicalResult.AlignedSequences.Should().OnlyContain(r => r == same,
            "100 identical inputs align to 100 identical, gap-free rows");
        AssertMsaInvariants(identicalResult, Enumerable.Repeat(same, 100).ToArray());
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: sequences of length 1 (and empty sequences in the set)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: length-1 sequences and empty sequences in the set

    /// <summary>
    /// BE: a collection of length-1 sequences is the degenerate per-sequence size corner.
    /// The MSA must still produce equal-length rows that gap-strip back to each single
    /// base, with count preserved and no all-gap column
    /// (Multiple_Sequence_Alignment.md §6.1). Pinned for all-equal single bases (the
    /// trivial single-column MSA) and for a mix of distinct single bases.
    /// </summary>
    private static readonly string[][] LengthOneCases =
    {
        new[] { "A", "A", "A" },
        new[] { "A", "C", "G", "T" },
        new[] { "A", "A", "C", "C", "G" },
    };

    [TestCaseSource(nameof(LengthOneCases))]
    public void MultipleAlign_LengthOneSequences_AreWellFormedMsa(string[] inputs)
    {
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();

        var act = () => SequenceAligner.MultipleAlign(seqs);
        act.Should().NotThrow("length-1 sequences are a defined boundary, not a crash");

        var result = act();
        result.AlignedSequences.Should().HaveCount(inputs.Length, "count is preserved");
        AssertMsaInvariants(result, inputs);
    }

    /// <summary>
    /// BE: empty sequences inside the collection (DnaSequence allows the empty string) must
    /// be handled WITHOUT rejecting the whole alignment — equal row lengths, the empty
    /// rows gap-strip to "", the non-empty rows recover their inputs, and no column is all
    /// gaps (Multiple_Sequence_Alignment.md §3.3, §6.1). Pinned with one and several empty
    /// sequences mixed among non-empty ones, and the all-empty extreme.
    /// </summary>
    private static readonly string[][] EmptyInCollectionCases =
    {
        new[] { "ATGC", "", "ATGC" },
        new[] { "", "ATGC", "", "GCAT" },
        new[] { "", "", "" },
    };

    [TestCaseSource(nameof(EmptyInCollectionCases))]
    public void MultipleAlign_EmptySequencesInCollection_HandledGracefully(string[] inputs)
    {
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();

        var act = () => SequenceAligner.MultipleAlign(seqs);
        act.Should().NotThrow("empty sequences inside the collection must not crash the MSA");

        var result = act();
        result.AlignedSequences.Should().HaveCount(inputs.Length, "count is preserved");

        // All-empty inputs collapse to zero-length rows; M-INV-* still hold (no columns).
        AssertMsaInvariants(result, inputs);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity: a few related sequences align to equal-length, recoverable rows
    // ───────────────────────────────────────────────────────────────────

    #region Positive sanity — related sequences align to recoverable equal-length rows

    /// <summary>
    /// Positive sanity: three related sequences that differ only by a single internal
    /// deletion align into one common coordinate system. The shorter sequence is padded
    /// with a gap, every row has equal length, gap-stripping each row recovers its exact
    /// input, the consensus length matches the rows, and the SP TotalScore is the
    /// independently recomputed sum-of-pairs over the final columns. This proves the fuzz
    /// harness asserts a real, theory-correct MSA — equal-length rows that gap-strip back
    /// to the inputs with a consistent score — not merely "did not crash".
    /// — Multiple_Sequence_Alignment.md §2.2 (SP score), §2.4.
    /// </summary>
    [Test]
    public void MultipleAlign_RelatedSequences_AlignToRecoverableEqualLengthRows()
    {
        var scoring = SequenceAligner.SimpleDna;
        var inputs = new[] { "ATGCATGC", "ATGCATGC", "ATGATGC" }; // third drops the 4th base
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();

        var result = SequenceAligner.MultipleAlign(seqs, scoring);

        AssertMsaInvariants(result, inputs);

        // At least one row must carry a gap (the shorter sequence is padded into place).
        result.AlignedSequences.Any(r => r.Contains('-'))
            .Should().BeTrue("the shorter related sequence is padded with a gap into the common layout");

        // The public SP TotalScore equals an independent column-based recomputation
        // (match/mismatch from the matrix, gap-vs-base = GapExtend, gap-vs-gap = 0).
        string[] rows = result.AlignedSequences;
        int width = rows[0].Length;
        int recomputed = 0;
        for (int col = 0; col < width; col++)
            for (int i = 0; i < rows.Length; i++)
                for (int j = i + 1; j < rows.Length; j++)
                {
                    char x = rows[i][col];
                    char y = rows[j][col];
                    if (x == '-' && y == '-')
                        recomputed += 0;
                    else if (x == '-' || y == '-')
                        recomputed += scoring.GapExtend;
                    else
                        recomputed += x == y ? scoring.Match : scoring.Mismatch;
                }

        result.TotalScore.Should().Be(recomputed,
            "TotalScore is the sum-of-pairs score recomputed over the final alignment columns");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ALIGN-STATS-001 — pairwise alignment statistics : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ALIGN-STATS-001 — pairwise alignment statistics (Alignment)
    // Checklist: docs/checklists/03_FUZZING.md, row 226.
    // Fuzz strategy exercised for THIS unit:
    //   • BE = Boundary Exploitation — the degenerate column-classification corners of
    //          the statistics pass: IDENTICAL rows (100% identity, 0 gaps), NO OVERLAP
    //          / completely different rows (0% identity, all mismatch), an ALL-GAP
    //          alignment (every column a gap — Identity 0, Gaps 100; the denominator L
    //          is the gap-INCLUSIVE alignment length so there is no divide-by-zero/NaN
    //          on a positive-length all-gap alignment), and the EMPTY alignment (the
    //          documented Empty fast path — denominator undefined, returns
    //          AlignmentStatistics.Empty). — docs/checklists/03_FUZZING.md §Description
    //          (strategy code BE = 0, -1, MaxInt, empty).
    //
    // The statistics contract under test (docs/algorithms/Alignment/Alignment_Statistics.md):
    //   Each of the L equal-length columns is classified (§2.2): GAP if either row has
    //   '-', else IDENTICAL if the two characters are equal, else MISMATCH. Let M =
    //   identical, X = mismatch, G = gap; then M + X + G = L (INV-01). With the
    //   gap-INCLUSIVE denominator L (§2.1, EMBOSS needle convention):
    //       Identity%   = M / L × 100
    //       Similarity% = (M + Sim⁺) / L × 100   (Sim⁺ = mismatch columns whose
    //                                              substitution score is POSITIVE, §2.2)
    //       Gaps%       = G / L × 100
    //   For every DNA model on SequenceAligner (Mismatch < 0) no mismatch is similar, so
    //   Similarity == Identity (§5.2); a caller-supplied Mismatch > 0 model makes
    //   Similarity exceed Identity (§2.5, §5.2). Empty alignment → AlignmentStatistics.Empty
    //   (§3.3, §6.1); null alignment → ArgumentNullException (§3.3).
    //
    // KEY invariants asserted across the fuzz fodder (Alignment_Statistics.md §2.4):
    //   • ST-INV-01 — Matches + Mismatches + Gaps = AlignmentLength (the three classes
    //     partition every column, §2.4 INV-01).
    //   • ST-INV-02 — Identity% ≤ Similarity% ≤ 100% (the similar set ⊇ the identical
    //     set; §2.4 INV-02). Combined with Identity ≥ 0 this pins every percentage into
    //     [0,100] — the fuzz bar's "no out-of-range / NaN / Infinity" guarantee.
    //   • ST-INV-03 — all three percentages are finite (never NaN/Infinity) on any
    //     POSITIVE-length alignment, including the all-gap alignment whose denominator L
    //     is the gap-inclusive length (NOT the zero ungapped-column count).
    // — fuzz bar additionally per docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing": no
    //   crash/hang/corruption; a documented validation throw (null) is acceptable.
    //
    // Expected counts/percentages below are derived INDEPENDENTLY from the §2.2 column
    // model on HAND-BUILT AlignmentResult rows (the statistics pass operates on the two
    // gapped rows directly — §1, §3.1 — so we construct the rows by hand rather than
    // echo an aligner's output), never read back off the code's own arrays.
    // ───────────────────────────────────────────────────────────────────────────

    #region ALIGN-STATS-001 — pairwise alignment statistics

    #region Helpers (ALIGN-STATS-001)

    /// <summary>
    /// Builds an <see cref="AlignmentResult"/> directly from two equal-length gapped
    /// rows so the statistics pass (which consumes the rows alone — Alignment_Statistics.md
    /// §1, §3.1) can be exercised on hand-chosen boundary alignments whose expected
    /// counts are derived independently from the §2.2 column model.
    /// </summary>
    private static AlignmentResult Rows(string row1, string row2) =>
        new(row1, row2, 0, AlignmentType.Global, 0, 0,
            row1.Length - 1, row2.Length - 1);

    /// <summary>
    /// Re-derives Matches/Mismatches/Gaps independently from the §2.2 column model and
    /// asserts the statistics structural invariants: ST-INV-01 (M+X+G=L), the
    /// gap-inclusive percentage formulas (§2.2/§3.2), ST-INV-02 (Identity ≤ Similarity ≤
    /// 100 with Identity ≥ 0), and ST-INV-03 (all percentages finite — no NaN/Infinity).
    /// For SimpleDna (Mismatch &lt; 0) Similarity must equal Identity (§5.2).
    /// </summary>
    private static void AssertStatsInvariants(
        AlignmentStatistics stats, string row1, string row2, ScoringMatrix scoring)
    {
        int len = row1.Length;
        int m = 0, x = 0, g = 0;
        for (int i = 0; i < len; i++)
        {
            char a = row1[i];
            char b = row2[i];
            if (a == '-' || b == '-') g++;
            else if (a == b) m++;
            else x++;
        }
        int simPlus = scoring.Mismatch > 0 ? x : 0; // §2.2: mismatch column is similar iff its score is positive

        // The reported counts must match the independent column classification.
        stats.AlignmentLength.Should().Be(len, "L is the gap-inclusive number of columns (§2.1, §3.2)");
        stats.Matches.Should().Be(m, "Matches = identical columns (§2.2)");
        stats.Mismatches.Should().Be(x, "Mismatches = non-gap non-identical columns (§2.2)");
        stats.Gaps.Should().Be(g, "Gaps = columns with a gap on either row (§2.2)");

        // ST-INV-01: the three classes partition every column.
        (stats.Matches + stats.Mismatches + stats.Gaps).Should().Be(stats.AlignmentLength,
            "ST-INV-01: Matches + Mismatches + Gaps = AlignmentLength (§2.4 INV-01)");

        // Gap-inclusive percentage formulas (§2.2, §3.2), recomputed from the counts.
        double expIdentity = (double)m / len * 100;
        double expSimilarity = (double)(m + simPlus) / len * 100;
        double expGap = (double)g / len * 100;
        stats.Identity.Should().BeApproximately(expIdentity, 1e-9, "Identity% = M/L × 100 (§2.2)");
        stats.Similarity.Should().BeApproximately(expSimilarity, 1e-9, "Similarity% = (M + Sim⁺)/L × 100 (§2.2)");
        stats.GapPercent.Should().BeApproximately(expGap, 1e-9, "Gaps% = G/L × 100 (§2.2)");

        // ST-INV-03: every percentage is finite (no divide-by-zero/NaN on positive L).
        double.IsNaN(stats.Identity).Should().BeFalse("ST-INV-03: Identity% is finite, never NaN");
        double.IsNaN(stats.Similarity).Should().BeFalse("ST-INV-03: Similarity% is finite, never NaN");
        double.IsNaN(stats.GapPercent).Should().BeFalse("ST-INV-03: Gaps% is finite, never NaN");
        double.IsInfinity(stats.Identity).Should().BeFalse("ST-INV-03: Identity% is finite, never Infinity");
        double.IsInfinity(stats.Similarity).Should().BeFalse("ST-INV-03: Similarity% is finite, never Infinity");
        double.IsInfinity(stats.GapPercent).Should().BeFalse("ST-INV-03: Gaps% is finite, never Infinity");

        // ST-INV-02: Identity ≤ Similarity ≤ 100, and Identity ≥ 0 ⇒ every % ∈ [0,100].
        stats.Identity.Should().BeGreaterThanOrEqualTo(0).And.BeLessThanOrEqualTo(100,
            "ST-INV-02: Identity% ∈ [0,100]");
        stats.Similarity.Should().BeGreaterThanOrEqualTo(stats.Identity).And.BeLessThanOrEqualTo(100,
            "ST-INV-02: Identity% ≤ Similarity% ≤ 100 (the similar set ⊇ the identical set, §2.4 INV-02)");
        stats.GapPercent.Should().BeGreaterThanOrEqualTo(0).And.BeLessThanOrEqualTo(100,
            "ST-INV-02: Gaps% ∈ [0,100]");

        // §5.2: for a negative-mismatch DNA model no mismatch is similar ⇒ Similarity == Identity.
        if (scoring.Mismatch < 0)
            stats.Similarity.Should().BeApproximately(stats.Identity, 1e-9,
                "§5.2: with Mismatch < 0 no mismatch is similar, so Similarity equals Identity");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: identical alignment (100% identity, 0 gaps)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: identical rows (Identity = Similarity = 100%, Gaps 0%)

    /// <summary>
    /// BE: two IDENTICAL gapped rows are the perfect-identity boundary — every column is
    /// identical, none is a gap or mismatch, so M = L, X = G = 0 and Identity =
    /// Similarity = 100%, Gaps = 0% (Alignment_Statistics.md §2.2, §6.1 "Perfect
    /// identity"). Expected counts derived independently from the column model; pinned
    /// for several lengths including the single-column extreme.
    /// </summary>
    [TestCase("A", TestName = "Stats_Identical_SingleColumn")]
    [TestCase("ACGT", TestName = "Stats_Identical_Len4")]
    [TestCase("AAAA", TestName = "Stats_Identical_Homopolymer")]
    [TestCase("GATTACAGATTACA", TestName = "Stats_Identical_Len14")]
    public void Stats_IdenticalRows_AreHundredPercentIdentityNoGaps(string seq)
    {
        var scoring = SequenceAligner.SimpleDna;
        var stats = SequenceAligner.CalculateStatistics(Rows(seq, seq), scoring);

        stats.Matches.Should().Be(seq.Length, "every column of identical rows is a match");
        stats.Mismatches.Should().Be(0, "identical rows have no mismatch column");
        stats.Gaps.Should().Be(0, "identical rows have no gap column");
        stats.AlignmentLength.Should().Be(seq.Length, "L equals the row length");
        stats.Identity.Should().BeApproximately(100.0, 1e-9, "M = L ⇒ Identity = 100%");
        stats.Similarity.Should().BeApproximately(100.0, 1e-9, "all-identical ⇒ Similarity = 100%");
        stats.GapPercent.Should().BeApproximately(0.0, 1e-9, "no gap column ⇒ Gaps = 0%");
        AssertStatsInvariants(stats, seq, seq, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: no overlap (completely different — 0% identity)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: no overlap / completely different (0% identity, all mismatch)

    /// <summary>
    /// BE: two equal-length rows that share NO base in any column (no gaps) are the
    /// no-overlap / completely-different boundary — every column is a mismatch, so
    /// X = L, M = G = 0 and Identity = 0%, Gaps = 0%; under SimpleDna (Mismatch &lt; 0)
    /// no mismatch is similar so Similarity = 0% too (Alignment_Statistics.md §2.2,
    /// §5.2). Expected counts derived independently from the column model.
    /// </summary>
    [TestCase("AAAA", "GGGG", TestName = "Stats_NoOverlap_AsVsGs")]
    [TestCase("ACAC", "GTGT", TestName = "Stats_NoOverlap_DisjointAlphabets")]
    [TestCase("A", "C", TestName = "Stats_NoOverlap_SingleColumn")]
    public void Stats_NoOverlapRows_AreZeroPercentIdentity(string row1, string row2)
    {
        var scoring = SequenceAligner.SimpleDna;
        var stats = SequenceAligner.CalculateStatistics(Rows(row1, row2), scoring);

        stats.Matches.Should().Be(0, "completely different rows share no column ⇒ no match");
        stats.Mismatches.Should().Be(row1.Length, "every column is a mismatch");
        stats.Gaps.Should().Be(0, "no-overlap fodder carries no gap column");
        stats.Identity.Should().BeApproximately(0.0, 1e-9, "M = 0 ⇒ Identity = 0%");
        stats.Similarity.Should().BeApproximately(0.0, 1e-9, "SimpleDna Mismatch < 0 ⇒ no similar mismatch ⇒ 0%");
        stats.GapPercent.Should().BeApproximately(0.0, 1e-9, "no gap column ⇒ Gaps = 0%");
        AssertStatsInvariants(stats, row1, row2, scoring);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: all-gap alignment (denominator must not divide by zero)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: all-gap alignment (Identity 0%, Gaps 100%, no divide-by-zero)

    /// <summary>
    /// BE: an ALL-GAP alignment — one row entirely gaps against a row of bases (an
    /// empty-vs-length-n alignment, the canonical all-gap output) — has G = L gap
    /// columns, M = X = 0, so Identity = 0% and Gaps = 100% (Alignment_Statistics.md
    /// §6.1 "All-gap column run"). The KEY fuzz point: the denominator L is the
    /// gap-INCLUSIVE alignment length (§2.1), so even though there are ZERO ungapped
    /// columns the percentages are computed over a POSITIVE L and must NOT divide by
    /// zero / produce NaN / Infinity (ST-INV-03). Pinned in both row orientations.
    /// </summary>
    [TestCase("ACGT", TestName = "Stats_AllGap_Len4")]
    [TestCase("A", TestName = "Stats_AllGap_SingleColumn")]
    [TestCase("ACGTACGTACGT", TestName = "Stats_AllGap_Len12")]
    public void Stats_AllGapAlignment_IsZeroIdentityHundredGapsNoNaN(string bases)
    {
        var scoring = SequenceAligner.SimpleDna;
        string gapRow = new string('-', bases.Length);

        // gaps on row1 vs bases on row2, and the mirror orientation.
        foreach (var (r1, r2) in new[] { (gapRow, bases), (bases, gapRow) })
        {
            var stats = SequenceAligner.CalculateStatistics(Rows(r1, r2), scoring);

            stats.Matches.Should().Be(0, "an all-gap alignment has no identical column");
            stats.Mismatches.Should().Be(0, "an all-gap alignment has no mismatch column");
            stats.Gaps.Should().Be(bases.Length, "every column is a gap column");
            stats.AlignmentLength.Should().Be(bases.Length, "L is the gap-inclusive length (§2.1)");
            stats.Identity.Should().BeApproximately(0.0, 1e-9, "M = 0 over positive L ⇒ Identity = 0%");
            stats.GapPercent.Should().BeApproximately(100.0, 1e-9, "G = L ⇒ Gaps = 100%");
            double.IsNaN(stats.Identity).Should().BeFalse("no divide-by-zero: L is gap-inclusive and positive");
            double.IsNaN(stats.GapPercent).Should().BeFalse("no divide-by-zero on the all-gap denominator");
            AssertStatsInvariants(stats, r1, r2, scoring);
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty / null alignment (documented fast path + guard)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Boundary: empty alignment (Empty fast path) + null guard

    /// <summary>
    /// BE: an EMPTY alignment (AlignedSequence1 null or empty) has an UNDEFINED
    /// denominator, so the documented fast path returns AlignmentStatistics.Empty —
    /// all-zero counts and percentages, never a divide-by-zero NaN/Infinity
    /// (Alignment_Statistics.md §3.3, §6.1 "Empty alignment"). Pinned for the
    /// empty-string row, the AlignmentResult.Empty sentinel, and a null first row.
    /// </summary>
    [Test]
    public void Stats_EmptyAlignment_ReturnsEmptyStatistics()
    {
        var fromEmptyRows = SequenceAligner.CalculateStatistics(Rows("", ""));
        fromEmptyRows.Should().Be(AlignmentStatistics.Empty,
            "an empty alignment has an undefined denominator and returns AlignmentStatistics.Empty (§6.1)");

        var fromSentinel = SequenceAligner.CalculateStatistics(AlignmentResult.Empty);
        fromSentinel.Should().Be(AlignmentStatistics.Empty,
            "the AlignmentResult.Empty sentinel maps to AlignmentStatistics.Empty");

        // The Empty statistics carry no NaN/Infinity (all fields are 0).
        double.IsNaN(fromEmptyRows.Identity).Should().BeFalse("Empty statistics are all-zero, never NaN");
        double.IsNaN(fromEmptyRows.GapPercent).Should().BeFalse("Empty statistics are all-zero, never NaN");
        fromEmptyRows.AlignmentLength.Should().Be(0, "an empty alignment has length 0");
    }

    /// <summary>
    /// BE/validation: a NULL alignment is the documented validation gate — it must throw
    /// ArgumentNullException (the ThrowIfNull guard, Alignment_Statistics.md §3.3, §6.1
    /// "Null alignment"), NOT a NullReferenceException from a downstream dereference.
    /// </summary>
    [Test]
    public void Stats_NullAlignment_ThrowsArgumentNullException()
    {
        var act = () => SequenceAligner.CalculateStatistics(null!);

        act.Should().Throw<ArgumentNullException>(
            "a null alignment hits the documented null gate, not a raw dereference (§3.3)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity: a hand-checkable mixed alignment + positive-mismatch model
    // ───────────────────────────────────────────────────────────────────

    #region Positive sanity — hand-checked mixed alignment and the similar-mismatch rule

    /// <summary>
    /// Positive sanity: a hand-built mixed alignment with a known column breakdown pins
    /// the exact counts and gap-inclusive percentages, proving the harness asserts a
    /// real, theory-correct statistic — not merely "did not crash". Rows
    /// "AC-GTA" / "ACCGAA" over L = 6 columns classify (§2.2) as:
    ///   col0 A|A match, col1 C|C match, col2 -|C gap, col3 G|G match, col4 T|A mismatch,
    ///   col5 A|A match ⇒ M = 4, X = 1, G = 1, L = 6.
    /// So Identity = 4/6 × 100 = 66.666…%, Gaps = 1/6 × 100 = 16.666…%, and under
    /// SimpleDna (Mismatch &lt; 0) Similarity = Identity (§5.2). Values derived
    /// independently from the doc's column model, not read off the implementation.
    /// </summary>
    [Test]
    public void Stats_HandCheckedMixedAlignment_HasExactCountsAndPercentages()
    {
        var scoring = SequenceAligner.SimpleDna;
        const string row1 = "AC-GTA";
        const string row2 = "ACCGAA";

        var stats = SequenceAligner.CalculateStatistics(Rows(row1, row2), scoring);

        stats.Matches.Should().Be(4, "hand count: cols 0,1,3,5 are identical");
        stats.Mismatches.Should().Be(1, "hand count: col 4 (T vs A) is the only mismatch");
        stats.Gaps.Should().Be(1, "hand count: col 2 (- vs C) is the only gap");
        stats.AlignmentLength.Should().Be(6, "L = 6 columns");
        stats.Identity.Should().BeApproximately(4.0 / 6.0 * 100, 1e-9, "Identity = 4/6 × 100 (§2.2)");
        stats.GapPercent.Should().BeApproximately(1.0 / 6.0 * 100, 1e-9, "Gaps = 1/6 × 100 (§2.2)");
        stats.Similarity.Should().BeApproximately(stats.Identity, 1e-9, "SimpleDna ⇒ Similarity == Identity (§5.2)");
        AssertStatsInvariants(stats, row1, row2, scoring);
    }

    /// <summary>
    /// Positive sanity: the §2.5/§5.2 distinction — a caller-supplied model with
    /// Mismatch &gt; 0 makes a mismatch column count as SIMILAR, so Similarity EXCEEDS
    /// Identity (the similar set strictly ⊇ the identical set). On the same
    /// "AC-GTA"/"ACCGAA" alignment (M = 4, X = 1, G = 1, L = 6) a positive-mismatch
    /// model gives Similarity = (4 + 1)/6 × 100 > Identity = 4/6 × 100, while Identity
    /// and Gaps are model-independent. This pins the documented parameterised-similarity
    /// behaviour the campaign's §5.4 fix corrected (away from the non-gap-fraction rule).
    /// </summary>
    [Test]
    public void Stats_PositiveMismatchModel_SimilarityExceedsIdentity()
    {
        // A model with Mismatch > 0: any non-identical column scores positively ⇒ similar.
        var positiveMismatch = new ScoringMatrix(Match: 2, Mismatch: 1, GapOpen: -2, GapExtend: -1);
        const string row1 = "AC-GTA";
        const string row2 = "ACCGAA";

        var stats = SequenceAligner.CalculateStatistics(Rows(row1, row2), positiveMismatch);

        stats.Identity.Should().BeApproximately(4.0 / 6.0 * 100, 1e-9, "Identity is model-independent: 4/6 × 100");
        stats.Similarity.Should().BeApproximately(5.0 / 6.0 * 100, 1e-9,
            "Mismatch > 0 ⇒ the 1 mismatch column is similar ⇒ (4+1)/6 × 100 (§2.2, §2.5)");
        stats.Similarity.Should().BeGreaterThan(stats.Identity,
            "a positive-mismatch model makes Similarity exceed Identity (§5.2)");
        AssertStatsInvariants(stats, row1, row2, positiveMismatch);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Randomized boundary sweep (no crash/hang/NaN; contract always holds)
    // ───────────────────────────────────────────────────────────────────

    #region BE — Randomized boundary sweep (invariants hold on arbitrary gapped rows)

    /// <summary>
    /// BE: a randomized sweep over arbitrary equal-length gapped rows (random mix of
    /// A/C/G/T and the gap symbol on both sides, including the all-gap and gap-free
    /// extremes that fall out by chance) must NEVER crash/hang/produce NaN/Infinity and
    /// must ALWAYS satisfy the statistics contract: M + X + G = L, every percentage in
    /// [0,100], Identity ≤ Similarity ≤ 100, and Similarity == Identity under SimpleDna
    /// (Alignment_Statistics.md §2.4). Locally fixed-seed RNG keeps the fodder
    /// reproducible; CancelAfter bounds the run.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Stats_RandomizedGappedRows_AlwaysHoldStatisticsContract()
    {
        var scoring = SequenceAligner.SimpleDna;
        const string alphabet = "ACGT-"; // includes the gap symbol so gap columns arise

        for (int trial = 0; trial < 500; trial++)
        {
            int len = Rng.Next(1, 40);
            var r1 = new char[len];
            var r2 = new char[len];
            for (int i = 0; i < len; i++)
            {
                r1[i] = alphabet[Rng.Next(alphabet.Length)];
                r2[i] = alphabet[Rng.Next(alphabet.Length)];
            }
            string row1 = new string(r1);
            string row2 = new string(r2);

            AlignmentStatistics stats = default;
            var act = () => stats = SequenceAligner.CalculateStatistics(Rows(row1, row2), scoring);
            act.Should().NotThrow("the statistics pass never crashes on arbitrary equal-length gapped rows");

            AssertStatsInvariants(stats, row1, row2, scoring);
        }
    }

    #endregion

    #endregion
}
