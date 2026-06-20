using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Matching area — exact pattern matching.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds boundary and malformed inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no infinite loop, no state
/// corruption, and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, ArgumentOutOfRangeException leaking from indexing).
/// Every input must result in EITHER a well-defined, theory-correct result, OR a
/// *documented, intentional* validation exception (here: ArgumentNullException on
/// a null pattern on the core API). A raw runtime exception, a wrong-result
/// corruption, or a hang on boundary input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-EXACT-001 — exact pattern matching (Matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 8.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — pattern longer than the text, the empty
///          pattern, a single-character pattern, the empty sequence, and the
///          pattern == sequence boundary (the whole text as a pattern).
///   • MC = Malformed Content — patterns that cannot occur (longer than the
///          text, absent symbols) must yield the defined "no match" result, not
///          a crash; the empty/degenerate pattern must hit the documented
///          empty-pattern contract rather than an indexing fault.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The exact-pattern-matching contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Exact pattern matching finds EVERY occurrence of a pattern P in a text T. A
/// position i is reported iff T[i..i+m-1] == P (INV-01). Occurrences may overlap
/// and are all reported. Positions are 0-based.
///   — docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md §2.2, §2.4
///     (INV-01 every reported i satisfies T[i..i+m-1]=P; INV-02 count == list
///     length; INV-03 Contains ⇔ at least one occurrence), §3.3, §6.1 (edge
///     cases). Sources: Gusfield (1997), Ukkonen (1995), Rosalind SUBS.
///
/// PAT-EXACT-001 has THREE documented surfaces with DIFFERENT, intentional
/// contracts at the boundary. Fuzzing pins all three, and the boundary between
/// them, so none can silently drift:
///
/// (1) The CORE suffix-tree API — SuffixTree.FindAllOccurrences /
///     CountOccurrences / Contains (string and ReadOnlySpan&lt;char&gt; overloads,
///     src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs). Built from a
///     DnaSequence's cached tree (DnaSequence.SuffixTree, DnaSequence.cs line 48).
///     Documented contract (Exact_Pattern_Search.md §3.3, §5.2, §6.1):
///       • null string pattern        → ArgumentNullException (explicit guard,
///         SuffixTreeSearchContracts.EnsureNotNull) — the validation gate, NOT a
///         NullReferenceException;
///       • empty pattern              → all valid start positions [0..n-1];
///         CountOccurrences == text length; Contains == true. (For an empty TEXT
///         this collapses to the empty list / 0 / true — n == 0.)
///       • pattern longer than text   → no full traversal possible → empty list,
///         count 0, Contains false;
///       • single-char pattern        → every occurrence of that base, overlap-
///         aware;
///       • pattern == text            → exactly one occurrence at position 0.
///     INV-02 (count == list length) and INV-03 (Contains ⇔ ≥1 occurrence) must
///     hold across ALL boundary inputs, so no surface can disagree with another.
///
/// (2) The WRAPPER MotifFinder.FindExactMotif(DnaSequence, string)
///     (MotifFinder.cs lines 24–37). By documented design it:
///       • returns an EMPTY result for a null OR empty motif (the wrapper-level
///         IsNullOrEmpty guard, yield break — it does NOT inherit the core's
///         "empty pattern ⇒ all positions" contract);
///       • upper-cases the motif before searching (case-insensitive DNA);
///       • yields positions SORTED ascending.
///     The strict/lenient boundary: the same empty pattern the CORE turns into
///     [0..n-1], the wrapper turns into the empty result. We pin both so the
///     boundary is explicit and cannot drift.
///
/// (3) The WRAPPER GenomicAnalyzer.FindMotif(DnaSequence, string)
///     (GenomicAnalyzer.cs lines 136–143). Same empty-motif guard (empty result
///     for null/empty motif) and upper-casing, but returns the suffix-tree list
///     DIRECTLY (DFS order, not re-sorted). We pin its boundary results too.
///
/// All inputs here are ASCII DNA / boundary strings; randomness (where used) is
/// from a locally fixed-seed Random so the fuzz fodder is fully reproducible and
/// adding tests cannot perturb any other fixture.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MatchingFuzzTests
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

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-EXACT-001 — exact pattern matching : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PAT-EXACT-001 — exact pattern matching

    // ───────────────────────────────────────────────────────────────────
    //  Core suffix-tree API: FindAllOccurrences / CountOccurrences / Contains
    // ───────────────────────────────────────────────────────────────────

    #region BE/MC — Boundary: pattern longer than the sequence

    /// <summary>
    /// BE/MC: a pattern strictly longer than the text can never occur — no full
    /// edge traversal is possible (Exact_Pattern_Search.md §6.1). The core API
    /// must return the DEFINED "no match" result on every surface — an empty
    /// position list, a count of 0, and Contains == false — never an
    /// IndexOutOfRange from over-reading the text. INV-02/INV-03 must agree.
    /// </summary>
    [TestCase("ACGT", "ACGTA", TestName = "FindOccurrences_PatternLongerByOne_NoMatch")]
    [TestCase("ACGT", "ACGTACGTACGT", TestName = "FindOccurrences_PatternMuchLonger_NoMatch")]
    [TestCase("A", "AA", TestName = "FindOccurrences_PatternLongerThanSingleBase_NoMatch")]
    public void FindOccurrences_PatternLongerThanSequence_IsNoMatch(string text, string pattern)
    {
        var tree = new DnaSequence(text).SuffixTree;

        IReadOnlyList<int> positions = tree.FindAllOccurrences(pattern);
        int count = tree.CountOccurrences(pattern);
        bool contains = tree.Contains(pattern);

        positions.Should().BeEmpty(
            "a pattern longer than the text cannot occur (no full traversal is possible)");
        count.Should().Be(0, "the occurrence count must equal the empty list's length (INV-02)");
        contains.Should().BeFalse("Contains is true iff there is at least one occurrence (INV-03)");
    }

    #endregion

    #region BE — Boundary: empty pattern (core "all positions" contract)

    /// <summary>
    /// BE: the empty pattern is the lower size boundary. The CORE suffix-tree API
    /// defines it as matching at EVERY position: FindAllOccurrences returns
    /// [0..n-1], CountOccurrences returns the text length, and Contains is true
    /// (Exact_Pattern_Search.md §5.2, §6.1; SuffixTree.Search.cs lines 84–85,
    /// 99, 26). No indexing, no division, no exception. INV-02/INV-03 must agree.
    /// </summary>
    [TestCase("ACGT", TestName = "FindOccurrences_EmptyPattern_AllPositions_Len4")]
    [TestCase("A", TestName = "FindOccurrences_EmptyPattern_AllPositions_Len1")]
    [TestCase("AAAAAAAA", TestName = "FindOccurrences_EmptyPattern_AllPositions_Homopolymer")]
    public void FindOccurrences_EmptyPattern_MatchesEveryPosition(string text)
    {
        var tree = new DnaSequence(text).SuffixTree;

        IReadOnlyList<int> positions = tree.FindAllOccurrences(string.Empty);
        int count = tree.CountOccurrences(string.Empty);
        bool contains = tree.Contains(string.Empty);

        positions.Should().Equal(Enumerable.Range(0, text.Length),
            "the empty pattern matches at every start position [0..n-1] on the core API");
        count.Should().Be(text.Length,
            "CountOccurrences of the empty pattern is the text length (INV-02)");
        contains.Should().BeTrue("the empty pattern is contained in any text (INV-03)");
    }

    #endregion

    #region BE — Boundary: empty sequence (empty text)

    /// <summary>
    /// BE: the empty SEQUENCE is the other lower boundary. Over an empty text the
    /// suffix tree has no positions, so: a non-empty pattern → no match (empty
    /// list / 0 / false); and the empty pattern collapses [0..n-1] to the EMPTY
    /// list with n == 0 (CountOccurrences == 0), while Contains stays true (the
    /// empty pattern is trivially present). No crash, no hang on the empty index.
    /// </summary>
    [Test]
    public void FindOccurrences_EmptySequence_HasNoNonEmptyMatchesAndNoCrash()
    {
        var tree = new DnaSequence(string.Empty).SuffixTree;

        var act = () =>
        {
            // A real pattern cannot occur in an empty text.
            tree.FindAllOccurrences("A").Should().BeEmpty(
                "an empty text contains no occurrence of any non-empty pattern");
            tree.CountOccurrences("A").Should().Be(0);
            tree.Contains("A").Should().BeFalse();

            // The empty pattern over an empty text: [0..n-1] with n == 0 is empty.
            tree.FindAllOccurrences(string.Empty).Should().BeEmpty(
                "[0..n-1] over an empty text (n == 0) is the empty position list");
            tree.CountOccurrences(string.Empty).Should().Be(0,
                "the empty text has length 0, so the empty-pattern count is 0 (INV-02)");
            tree.Contains(string.Empty).Should().BeTrue(
                "the empty pattern is trivially contained even in the empty text (INV-03)");
        };

        act.Should().NotThrow("the empty sequence is a defined boundary index, not an error");
    }

    #endregion

    #region BE — Boundary: single-character pattern (all occurrences, overlap-aware)

    /// <summary>
    /// BE: a one-character pattern is the minimal non-empty pattern. It must report
    /// EVERY position of that base, overlap-aware. Verified against a brute-force
    /// scan over a fixed text and over a fixed-seed random text, so the result is
    /// neither under- nor over-counted at the size-1 boundary. A base ABSENT from
    /// the text (MC) yields the defined empty result. INV-02/INV-03 must agree.
    /// </summary>
    [TestCase("ACGTACGT", 'A', TestName = "FindOccurrences_SingleChar_A_AllPositions")]
    [TestCase("AAAA", 'A', TestName = "FindOccurrences_SingleChar_Homopolymer_EveryPosition")]
    [TestCase("ACGT", 'T', TestName = "FindOccurrences_SingleChar_T_OnePosition")]
    [TestCase("AAAA", 'C', TestName = "FindOccurrences_SingleChar_AbsentBase_NoMatch")]
    public void FindOccurrences_SingleCharPattern_ReportsAllOccurrences(string text, char baseChar)
    {
        var tree = new DnaSequence(text).SuffixTree;
        string pattern = baseChar.ToString();

        var expected = Enumerable.Range(0, text.Length)
            .Where(i => text[i] == baseChar)
            .ToList();

        IReadOnlyList<int> positions = tree.FindAllOccurrences(pattern);
        int count = tree.CountOccurrences(pattern);
        bool contains = tree.Contains(pattern);

        positions.OrderBy(p => p).Should().Equal(expected,
            "a single-base pattern occurs at exactly the positions holding that base");
        count.Should().Be(expected.Count, "the count must equal the number of occurrences (INV-02)");
        contains.Should().Be(expected.Count > 0,
            "Contains is true iff the base appears at least once (INV-03)");
    }

    /// <summary>
    /// BE: single-character pattern over a fixed-seed random text — a stronger,
    /// reproducible fuzz of the size-1 boundary against a brute-force oracle for
    /// every one of the four bases. Pins that overlap-aware single-base matching
    /// stays correct on arbitrary (random) composition without crashing.
    /// </summary>
    [Test]
    public void FindOccurrences_SingleCharPattern_RandomText_MatchesBruteForce()
    {
        string text = RandomDna(500);
        var tree = new DnaSequence(text).SuffixTree;

        foreach (char baseChar in "ACGT")
        {
            var expected = Enumerable.Range(0, text.Length)
                .Where(i => text[i] == baseChar)
                .ToList();

            tree.FindAllOccurrences(baseChar.ToString()).OrderBy(p => p).Should().Equal(expected,
                $"every '{baseChar}' position must be reported, overlap-aware, in a random text");
            tree.CountOccurrences(baseChar.ToString()).Should().Be(expected.Count,
                "count must equal the brute-force occurrence count (INV-02)");
        }
    }

    #endregion

    #region BE — Boundary: pattern == sequence (whole text as pattern)

    /// <summary>
    /// BE: the pattern == text boundary. The whole text occurs exactly ONCE, at
    /// position 0 — never more (the suffix that starts at 0 is the only one with
    /// the full length) and never an off-by-one extra. CountOccurrences == 1,
    /// Contains == true. Verified over several lengths including the length-1 and
    /// homopolymer extremes. INV-02/INV-03 must agree.
    /// </summary>
    [TestCase("A", TestName = "FindOccurrences_PatternEqualsSequence_SingleBase")]
    [TestCase("ACGT", TestName = "FindOccurrences_PatternEqualsSequence_Distinct")]
    [TestCase("AAAA", TestName = "FindOccurrences_PatternEqualsSequence_Homopolymer")]
    [TestCase("GATATATGCATATACTT", TestName = "FindOccurrences_PatternEqualsSequence_DocExample")]
    public void FindOccurrences_PatternEqualsSequence_IsExactlyOneMatchAtZero(string text)
    {
        var tree = new DnaSequence(text).SuffixTree;

        IReadOnlyList<int> positions = tree.FindAllOccurrences(text);
        int count = tree.CountOccurrences(text);
        bool contains = tree.Contains(text);

        positions.Should().Equal(new[] { 0 },
            "the whole text occurs exactly once, at position 0");
        count.Should().Be(1, "pattern == text occurs exactly once (INV-02)");
        contains.Should().BeTrue("the whole text is trivially contained in itself (INV-03)");
    }

    #endregion

    #region INJ/MC — Injection: null pattern on the core API throws documented exception

    /// <summary>
    /// MC/INJ: a null string pattern on the core string API is the explicit
    /// validation gate — it must throw the *documented, intentional*
    /// ArgumentNullException (SuffixTreeSearchContracts.EnsureNotNull;
    /// Exact_Pattern_Search.md §3.1, §6.1), NOT a NullReferenceException and not a
    /// silent empty result. Pinned on all three core string entry points.
    /// </summary>
    [Test]
    public void CoreApi_NullStringPattern_ThrowsArgumentNullException()
    {
        var tree = new DnaSequence("ACGT").SuffixTree;

        var findAct = () => tree.FindAllOccurrences((string)null!);
        var countAct = () => tree.CountOccurrences((string)null!);
        var containsAct = () => tree.Contains((string)null!);

        findAct.Should().Throw<ArgumentNullException>(
            "a null pattern is rejected by the documented core validation gate, not dereferenced");
        countAct.Should().Throw<ArgumentNullException>(
            "a null pattern is rejected by the documented core validation gate, not dereferenced");
        containsAct.Should().Throw<ArgumentNullException>(
            "a null pattern is rejected by the documented core validation gate, not dereferenced");
    }

    #endregion

    #region BE/OVF — Boundary: extremely long text and pattern do not hang

    /// <summary>
    /// BE/OVF: an extremely long valid text searched for a long pattern must
    /// complete without hang, infinite loop, or IndexOutOfRange, and the result
    /// must be theory-correct. A known-construction text ("AC" repeated) pins the
    /// exact overlap-aware occurrence set for the "ACAC" pattern, and the
    /// whole-text pattern still yields exactly one match at scale.
    /// </summary>
    [Test]
    public void FindOccurrences_ExtremelyLong_StaysCorrectAndDoesNotHang()
    {
        const int repeats = 100_000;                 // 200k-base text
        string text = string.Concat(Enumerable.Repeat("AC", repeats));
        var tree = new DnaSequence(text).SuffixTree;

        // "ACAC" occurs starting at every even index 0,2,4,... except the last
        // window cannot start past length-4. Overlap-aware count = repeats - 1.
        int count = tree.CountOccurrences("ACAC");
        count.Should().Be(repeats - 1,
            "'ACAC' in (AC)^n starts at every even index that leaves room for 4 chars");
        tree.FindAllOccurrences("ACAC").Count.Should().Be(count,
            "the listed occurrences must equal the count at scale (INV-02)");

        // The whole text is still exactly one match at position 0.
        tree.CountOccurrences(text).Should().Be(1,
            "pattern == text is one match even for a very long sequence");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Wrapper MotifFinder.FindExactMotif — empty-motif guard, sorted output
    // ───────────────────────────────────────────────────────────────────

    #region BE — wrapper FindExactMotif: empty/null motif → empty (NOT all positions)

    /// <summary>
    /// BE: the WRAPPER MotifFinder.FindExactMotif does NOT inherit the core's
    /// "empty pattern ⇒ all positions" contract — by documented design a null OR
    /// empty motif yields an EMPTY result (the wrapper-level IsNullOrEmpty guard,
    /// MotifFinder.cs line 27). This pins the strict/lenient boundary: the same
    /// empty pattern the core turns into [0..n-1], the wrapper turns into nothing.
    /// </summary>
    [Test]
    public void FindExactMotif_EmptyOrNullMotif_IsEmptyResult()
    {
        var seq = new DnaSequence("ACGTACGT");

        MotifFinder.FindExactMotif(seq, string.Empty).Should().BeEmpty(
            "the wrapper guards the empty motif and yields nothing, unlike the core 'all positions'");
        MotifFinder.FindExactMotif(seq, null!).Should().BeEmpty(
            "the wrapper guards a null motif and yields nothing, never a NullReferenceException");
    }

    #endregion

    #region BE/MC — wrapper FindExactMotif: boundary motifs, sorted, case-folded

    /// <summary>
    /// BE/MC: the wrapper must honor the same occurrence theory as the core for
    /// real motifs, but emit positions SORTED ascending (MotifFinder.cs line 33),
    /// and it upper-cases the motif (case-insensitive DNA). Pinned at the
    /// boundaries: a pattern longer than the text → empty; the whole text → exactly
    /// {0}; a 1-char motif → all positions; a lowercase motif → identical to its
    /// uppercase form.
    /// </summary>
    [Test]
    public void FindExactMotif_BoundaryMotifs_AreSortedAndCaseFolded()
    {
        var seq = new DnaSequence("ACGTACGT");

        // Pattern longer than the text → no match (MC boundary).
        MotifFinder.FindExactMotif(seq, "ACGTACGTA").Should().BeEmpty(
            "a motif longer than the sequence cannot occur");

        // Whole text as the motif → exactly one match at 0.
        MotifFinder.FindExactMotif(seq, "ACGTACGT").Should().Equal(new[] { 0 },
            "the whole sequence as a motif occurs exactly once, at position 0");

        // Single-char motif → every occurrence, SORTED.
        MotifFinder.FindExactMotif(seq, "A").Should().Equal(new[] { 0, 4 },
            "the wrapper yields all 'A' positions sorted ascending");
        MotifFinder.FindExactMotif(seq, "ACGT").Should().Equal(new[] { 0, 4 },
            "overlapping/repeated motif occurrences are all reported, sorted");

        // Lowercase motif folds to the same result (case-insensitive DNA).
        MotifFinder.FindExactMotif(seq, "acgt").Should().Equal(
            MotifFinder.FindExactMotif(seq, "ACGT"),
            "the motif is upper-cased before searching, so case must not change the result");
    }

    /// <summary>
    /// BE: the wrapper over an EMPTY sequence with a real motif must yield the
    /// empty result — no occurrence, no crash on the empty index.
    /// </summary>
    [Test]
    public void FindExactMotif_EmptySequence_RealMotif_IsEmpty()
    {
        var seq = new DnaSequence(string.Empty);

        MotifFinder.FindExactMotif(seq, "ACGT").Should().BeEmpty(
            "no motif occurs in an empty sequence");
    }

    /// <summary>
    /// MC: a null DnaSequence is rejected by the wrapper's explicit guard with the
    /// documented ArgumentNullException (MotifFinder.cs line 26), never a
    /// NullReferenceException. Because FindExactMotif is a lazy iterator, the guard
    /// only runs on enumeration — so the assertion materializes the result.
    /// </summary>
    [Test]
    public void FindExactMotif_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => MotifFinder.FindExactMotif(null!, "ACGT").ToList();

        act.Should().Throw<ArgumentNullException>(
            "a null sequence is rejected by the documented guard, not dereferenced");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Wrapper GenomicAnalyzer.FindMotif — empty-motif guard, direct list
    // ───────────────────────────────────────────────────────────────────

    #region BE/MC — wrapper FindMotif: boundary motifs

    /// <summary>
    /// BE/MC: GenomicAnalyzer.FindMotif shares the empty-motif guard (empty result
    /// for a null OR empty motif, GenomicAnalyzer.cs lines 138–139) and upper-cases
    /// the motif, but returns the suffix-tree list directly. Pinned at the
    /// boundaries: empty/null motif → empty; pattern longer than text → empty;
    /// whole text → exactly {0}; 1-char motif → all positions; lowercase folds.
    /// Positions are compared order-insensitively since this surface does not
    /// re-sort.
    /// </summary>
    [Test]
    public void FindMotif_BoundaryMotifs_AreDefined()
    {
        var seq = new DnaSequence("ACGTACGT");

        // Empty-motif guard (does NOT inherit core "all positions").
        GenomicAnalyzer.FindMotif(seq, string.Empty).Should().BeEmpty(
            "the wrapper guards the empty motif and yields nothing");
        GenomicAnalyzer.FindMotif(seq, null!).Should().BeEmpty(
            "the wrapper guards a null motif and yields nothing, never a NullReferenceException");

        // Pattern longer than the text → no match.
        GenomicAnalyzer.FindMotif(seq, "ACGTACGTA").Should().BeEmpty(
            "a motif longer than the sequence cannot occur");

        // Whole text → exactly one match at 0.
        GenomicAnalyzer.FindMotif(seq, "ACGTACGT").Should().Equal(new[] { 0 },
            "the whole sequence as a motif occurs exactly once, at position 0");

        // Single-char and repeated motifs → all occurrences (order-insensitive).
        GenomicAnalyzer.FindMotif(seq, "A").Should().BeEquivalentTo(new[] { 0, 4 },
            "all 'A' positions are reported");
        GenomicAnalyzer.FindMotif(seq, "ACGT").Should().BeEquivalentTo(new[] { 0, 4 },
            "all occurrences of the repeated motif are reported");

        // Lowercase folds to the same occurrence set.
        GenomicAnalyzer.FindMotif(seq, "acgt").Should().BeEquivalentTo(
            GenomicAnalyzer.FindMotif(seq, "ACGT"),
            "the motif is upper-cased before searching; case must not change the result");
    }

    /// <summary>
    /// BE: FindMotif over an EMPTY sequence with a real motif yields the empty
    /// result — no occurrence, no crash on the empty index.
    /// </summary>
    [Test]
    public void FindMotif_EmptySequence_RealMotif_IsEmpty()
    {
        var seq = new DnaSequence(string.Empty);

        GenomicAnalyzer.FindMotif(seq, "ACGT").Should().BeEmpty(
            "no motif occurs in an empty sequence");
    }

    #endregion

    #endregion
}
