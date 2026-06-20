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
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-APPROX-001 — approximate (Hamming) pattern matching (Matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 9.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the empty sequence and the empty pattern; a
///          single-base pattern; maxMismatches = 0 (reduces to EXACT matching);
///          maxMismatches ≥ pattern.Length (every equal-length window qualifies);
///          a negative maxMismatches (the documented validation gate).
///   • MC = Malformed Content — unequal-length comparisons on the direct
///          HammingDistance metric (rejected by contract); a pattern longer than
///          the sequence (no equal-length window exists → no match); and non-DNA
///          / control / non-ASCII characters on the string surface (handled by
///          plain case-folded comparison, never a crash).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The approximate-(Hamming)-matching contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PAT-APPROX-001 is the SUBSTITUTIONS-ONLY (Hamming) approximate-match family on
/// ApproximateMatcher: FindWithMismatches / CountApproximateOccurrences /
/// FindBestMatch and the direct HammingDistance metric. (The maxEdits /
/// Levenshtein family — FindWithEdits — is PAT-APPROX-002, row 10, and is NOT
/// exercised here.) Documented contract
/// (docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md §2.2, §3.1,
/// §3.3, §6.1; sources: Hamming 1950, Rosalind HAMM, Navarro 2001, Gusfield 1997):
///   • HammingDistance(s1, s2) — equal-length substitution count:
///       – null s1 or s2          → ArgumentNullException (the validation gate,
///         NOT a NullReferenceException);
///       – UNEQUAL lengths        → ArgumentException (Hamming is defined only on
///         equal-length strings — the doc rejects, it does not silently truncate);
///       – equal length           → d_H ≥ 0 (INV-01), symmetric d_H(s,t)=d_H(t,s)
///         (INV-03), case-insensitive (uppercased before comparison, §5.2);
///   • FindWithMismatches(string sequence, string pattern, int maxMismatches):
///       – empty/null sequence OR pattern → no matches (explicit source guard,
///         yield break — NOT the core exact-match "all positions" contract);
///       – pattern LONGER than sequence   → no matches (no equal-length window);
///       – maxMismatches < 0              → ArgumentOutOfRangeException;
///       – maxMismatches = 0              → reduces to EXACT matching (only
///         zero-mismatch windows qualify);
///       – maxMismatches ≥ pattern.Length → EVERY equal-length window qualifies
///         (every achievable mismatch count is within threshold) → all start
///         positions [0 .. n−m] reported;
///       – non-DNA / control / non-ASCII chars on the string surface → matched by
///         plain case-folded character comparison, never a crash (the string path
///         does NOT enforce the DNA alphabet — only the typed DnaSequence ctor
///         does);
///   • CountApproximateOccurrences == FindWithMismatches(...).Count() on every
///     boundary input (the count surface must agree with the listing surface);
///   • the DnaSequence wrapper overloads add NO null guard → a null DnaSequence
///     throws (NullReferenceException) because they dereference sequence.Sequence
///     (§3.3, §6.1) — pinned so the strict-vs-lenient surface split cannot drift.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-APPROX-002 — approximate (edit-distance / Levenshtein) matching (Matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 10.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the empty sequence and the empty pattern; a
///          single-base pattern; maxEdits = 0 (reduces to EXACT, equal-length-only
///          matching); maxEdits ≥ sequence length (every window is reachable, so
///          every start position matches); a negative maxEdits (the documented
///          validation gate, at −1 and int.MinValue).
///   • MC = Malformed Content — non-DNA / control / non-ASCII characters on the
///          string surface (matched by the case-folded variable-window scan, never
///          a crash); a null DnaSequence on the typed wrapper (documented throw).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The approximate-(edit-distance)-matching contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PAT-APPROX-002 is the LEVENSHTEIN (insertions + deletions + substitutions)
/// approximate-match family on ApproximateMatcher: FindWithEdits(...) and the
/// direct EditDistance(...) metric. (The substitutions-only Hamming family —
/// FindWithMismatches / HammingDistance — is PAT-APPROX-001, row 9, and is NOT
/// re-exercised here.) Unlike Hamming, edit distance permits UNEQUAL pattern and
/// window lengths, so FindWithEdits scans variable-length windows whose length
/// ranges over [pat.Length − maxEdits .. pat.Length + maxEdits] and may emit
/// MULTIPLE results at the SAME start position (one per qualifying window length).
/// Documented contract (docs/algorithms/Pattern_Matching/Edit_Distance.md §2.2,
/// §2.4, §3.1, §3.3, §6.1; sources: Levenshtein 1966, Wagner &amp; Fischer 1974,
/// Navarro 2001, Berger et al. 2021):
///   • EditDistance(s1, s2) — Wagner-Fischer Levenshtein:
///       – null s1 or s2          → ArgumentNullException (the validation gate,
///         NOT a NullReferenceException; ApproximateMatcher.cs lines 197–198);
///       – d = 0 iff identical (INV-01); d ≥ |len(a)−len(b)| (INV-02);
///         d ≤ max(len(a),len(b)) (INV-03); "" vs "abc" → 3 (§6.1);
///       – CASE-SENSITIVE: the core metric compares characters as-is (§5.2), so
///         it does NOT uppercase (unlike the search surface);
///   • FindWithEdits(string sequence, string pattern, int maxEdits):
///       – empty/null sequence OR pattern → no matches (explicit IsNullOrEmpty
///         guard, yield break — NOT the exact-match core's "all positions");
///       – maxEdits < 0              → ArgumentOutOfRangeException (lazy guard,
///         fires on enumeration; ApproximateMatcher.cs lines 124–125);
///       – maxEdits = 0              → minLen == maxLen == pattern length, so only
///         equal-length windows are checked and only zero-distance ones qualify →
///         reduces EXACTLY to exact matching;
///       – maxEdits ≥ sequence length → every pattern is within that many edits of
///         every reachable window, so every start position [0 .. n−minLen] yields
///         at least one match (no crash on the widened window range);
///       – non-DNA / control / non-ASCII chars on the string surface → matched by
///         the case-folded variable-window scan, never a crash (the string path
///         does NOT enforce the DNA alphabet — only the typed DnaSequence ctor does);
///       – the scan is verified against a brute-force oracle that REPLAYS the exact
///         (start i, window length len) double-loop and EditDistance call, so the
///         duplicate-position semantics cannot drift;
///   • the DnaSequence wrapper overload adds NO null guard → a null DnaSequence
///     throws (NullReferenceException) because it dereferences sequence.Sequence
///     (§3.1, §3.3, §6.1) — pinned so the strict-vs-lenient surface split is explicit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-IUPAC-001 — IUPAC (degenerate / ambiguity) pattern matching (Matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 11.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — motifs containing characters that are NOT one of
///          the 15 standard IUPAC DNA codes: amino-acid-only protein letters
///          (E/F/I/L/P/Q/Z, and J/O/U/X which are not nucleotide codes), and a
///          mixed valid-then-invalid motif so validation cannot be skipped because
///          the first window already failed to match.
///   • INJ = Injection — digits (0–9) inside the pattern, punctuation/gap symbols,
///          control / null bytes, and non-ASCII Unicode (e.g. 'α', 'ñ', full-width
///          digits) — every one must be REJECTED deterministically, never crash
///          with KeyNotFound (the dictionary lookup) or IndexOutOfRange.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The IUPAC-degenerate-matching contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PAT-IUPAC-001 expands ambiguity codes (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T,
/// M=A/C, B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G, N=A/C/G/T) over the DNA alphabet and
/// reports every length-m window whose every base lies in the per-position code
/// set. Documented contract (docs/algorithms/Pattern_Matching/
/// IUPAC_Degenerate_Matching.md §2.2, §3.1, §3.3, §6.1; INV-01/INV-02/INV-03;
/// sources: IUPAC-IUB 1970, NC-IUB 1984):
///   • the matching SURFACE is MotifFinder.FindDegenerateMotif (MotifFinder.cs
///     lines 90–119), a LAZY iterator over a typed DnaSequence whose sequence is
///     therefore always valid uppercase A/C/G/T — the fuzz target is the MOTIF:
///       – null DnaSequence            → ArgumentNullException (explicit guard,
///         line 92) — but it is lazy, so it fires only on enumeration;
///       – null OR empty motif         → no matches (IsNullOrEmpty guard, line 93,
///         yield break — NOT a crash);
///       – the motif is UPPER-CASED (ToUpperInvariant, line 96) BEFORE validation,
///         so lowercase IUPAC codes are accepted and fold to the same result;
///       – an INVALID IUPAC code in the motif (an amino-acid-only protein letter,
///         a digit, punctuation, a control byte, or non-ASCII Unicode) →
///         ArgumentException naming "motif" (ValidateIupacPattern, lines 70–82) —
///         a documented, intentional validation exception, NEVER a KeyNotFound
///         from the IupacCodes dictionary lookup or an IndexOutOfRange. Because the
///         method is a lazy iterator, validation fires on ENUMERATION, so the fuzz
///         assertions materialize the result;
///   • the DECISION-TABLE primitive IupacHelper.MatchesIupac(nucleotide, code)
///     (IupacHelper.cs lines 16–37) is the per-character oracle:
///       – a valid code returns the membership truth (R matches A and G, not C/T;
///         N matches all four bases) — the positive theory contract;
///       – an INVALID code throws ArgumentOutOfRangeException (the switch default,
///         lines 33–35), NOT KeyNotFound — a DIFFERENT documented exception type
///         from the MotifFinder surface, pinned so the two surfaces cannot drift.
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

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-APPROX-001 — approximate (Hamming) pattern matching : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PAT-APPROX-001 — approximate (Hamming) pattern matching

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty strings (empty sequence / empty pattern)
    // ───────────────────────────────────────────────────────────────────

    #region BE — empty strings: empty sequence / pattern → no matches, no crash

    /// <summary>
    /// BE: the empty sequence and the empty pattern are the lower size boundaries.
    /// FindWithMismatches guards both with an explicit IsNullOrEmpty check and
    /// yields nothing — it does NOT inherit the exact-match core's
    /// "empty pattern ⇒ all positions" contract (Approximate_Matching_Hamming.md
    /// §3.1, §6.1). Every empty combination — over every maxMismatches, including
    /// 0 and a large value — must be the defined empty result, never an
    /// IndexOutOfRange on the empty index. CountApproximateOccurrences must agree.
    /// </summary>
    [TestCase("", "ACGT", 0, TestName = "FindWithMismatches_EmptySequence_NoMatch")]
    [TestCase("", "ACGT", 5, TestName = "FindWithMismatches_EmptySequence_LargeMaxDist_NoMatch")]
    [TestCase("ACGTACGT", "", 0, TestName = "FindWithMismatches_EmptyPattern_NoMatch")]
    [TestCase("ACGTACGT", "", 5, TestName = "FindWithMismatches_EmptyPattern_LargeMaxDist_NoMatch")]
    [TestCase("", "", 0, TestName = "FindWithMismatches_BothEmpty_NoMatch")]
    public void FindWithMismatches_EmptyStrings_AreNoMatch(string sequence, string pattern, int maxMismatches)
    {
        IEnumerable<ApproximateMatchResult> results = ApproximateMatcher.FindWithMismatches(
            sequence, pattern, maxMismatches);

        results.Should().BeEmpty(
            "an empty sequence or pattern is guarded by the source and yields no approximate matches");
        ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches).Should().Be(0,
            "the count surface must agree with the empty listing surface on empty input");
    }

    /// <summary>
    /// BE: null sequence/pattern on the string surface is treated the same as
    /// empty by the IsNullOrEmpty guard — it yields nothing rather than throwing a
    /// NullReferenceException (Approximate_Matching_Hamming.md §3.1). FindBestMatch
    /// likewise returns null on any empty/null input. Pinned so the lenient string
    /// surface cannot drift into a raw dereference.
    /// </summary>
    [Test]
    public void FindWithMismatches_NullStrings_AreNoMatchNotCrash()
    {
        var findNullSeq = () => ApproximateMatcher.FindWithMismatches((string)null!, "ACGT", 1).ToList();
        var findNullPat = () => ApproximateMatcher.FindWithMismatches("ACGT", null!, 1).ToList();

        findNullSeq.Should().NotThrow("a null sequence is guarded as empty, never dereferenced");
        findNullPat.Should().NotThrow("a null pattern is guarded as empty, never dereferenced");

        findNullSeq().Should().BeEmpty("a null sequence yields no matches");
        findNullPat().Should().BeEmpty("a null pattern yields no matches");

        ApproximateMatcher.FindBestMatch(null!, "ACGT").Should().BeNull(
            "FindBestMatch returns null on a null sequence");
        ApproximateMatcher.FindBestMatch("ACGT", null!).Should().BeNull(
            "FindBestMatch returns null on a null pattern");
        ApproximateMatcher.FindBestMatch("", "").Should().BeNull(
            "FindBestMatch returns null when both inputs are empty");
    }

    /// <summary>
    /// BE/MC: the typed DnaSequence overload adds NO null guard and dereferences
    /// sequence.Sequence, so a null DnaSequence throws — the documented
    /// strict-surface behavior (Approximate_Matching_Hamming.md §3.3, §6.1). Pinned
    /// against the lenient string surface (which guards null as empty above) so the
    /// strict/lenient split is explicit and cannot silently converge.
    /// </summary>
    [Test]
    public void FindWithMismatches_NullDnaSequence_Throws()
    {
        var act = () => ApproximateMatcher.FindWithMismatches((DnaSequence)null!, "ACGT", 1).ToList();

        act.Should().Throw<NullReferenceException>(
            "the typed wrapper dereferences sequence.Sequence with no null guard (documented strict surface)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: unequal lengths
    // ───────────────────────────────────────────────────────────────────

    #region BE/MC — unequal lengths: HammingDistance rejects; search uses equal windows

    /// <summary>
    /// MC: the direct Hamming metric is defined ONLY on equal-length strings, so an
    /// unequal-length comparison must throw the documented ArgumentException —
    /// never silently truncate to the shorter string and never IndexOutOfRange
    /// (Approximate_Matching_Hamming.md §2.2, §3.3, §6.2). Pinned in both
    /// orientations (longer-first and shorter-first) and including the empty-vs-
    /// non-empty boundary, which is itself an unequal-length case.
    /// </summary>
    [TestCase("ACGT", "ACG", TestName = "HammingDistance_LongerVsShorter_Throws")]
    [TestCase("ACG", "ACGT", TestName = "HammingDistance_ShorterVsLonger_Throws")]
    [TestCase("ACGT", "", TestName = "HammingDistance_NonEmptyVsEmpty_Throws")]
    [TestCase("", "A", TestName = "HammingDistance_EmptyVsNonEmpty_Throws")]
    public void HammingDistance_UnequalLengths_ThrowsArgumentException(string s1, string s2)
    {
        var act = () => ApproximateMatcher.HammingDistance(s1, s2);

        act.Should().Throw<ArgumentException>(
            "Hamming distance is defined only on equal-length strings; unequal lengths are rejected")
            .Which.Should().NotBeOfType<ArgumentNullException>(
                "the unequal-length rejection is an ArgumentException, distinct from the null gate");
    }

    /// <summary>
    /// MC/INJ: null strings on the direct Hamming metric hit the explicit null gate
    /// (ArgumentNullException), NOT a NullReferenceException
    /// (Approximate_Matching_Hamming.md §3.3). The null check precedes the length
    /// check, so null-vs-shorter still reports the null gate.
    /// </summary>
    [Test]
    public void HammingDistance_NullArguments_ThrowArgumentNullException()
    {
        var nullFirst = () => ApproximateMatcher.HammingDistance(null!, "ACGT");
        var nullSecond = () => ApproximateMatcher.HammingDistance("ACGT", null!);
        var nullVsDifferentLength = () => ApproximateMatcher.HammingDistance(null!, "AC");

        nullFirst.Should().Throw<ArgumentNullException>("a null first string hits the null gate");
        nullSecond.Should().Throw<ArgumentNullException>("a null second string hits the null gate");
        nullVsDifferentLength.Should().Throw<ArgumentNullException>(
            "the null check precedes the length check, so null wins over the length mismatch");
    }

    /// <summary>
    /// MC: an UNEQUAL-LENGTH pattern vs sequence in the sliding-window search is a
    /// well-defined no-op when the pattern is longer than the sequence — no
    /// equal-length window exists, so the result is empty and never an
    /// IndexOutOfRange from over-reading the window (Approximate_Matching_Hamming.md
    /// §6.1). A maxMismatches large enough to "forgive" the length gap does NOT
    /// rescue it: Hamming compares only equal-length windows, so length still wins.
    /// </summary>
    [TestCase("ACGT", "ACGTA", 0, TestName = "FindWithMismatches_PatternLongerByOne_NoMatch")]
    [TestCase("ACGT", "ACGTACGT", 10, TestName = "FindWithMismatches_PatternMuchLonger_LargeMaxDist_NoMatch")]
    [TestCase("A", "AC", 5, TestName = "FindWithMismatches_PatternLongerThanSingleBase_NoMatch")]
    public void FindWithMismatches_PatternLongerThanSequence_IsNoMatch(
        string sequence, string pattern, int maxMismatches)
    {
        ApproximateMatcher.FindWithMismatches(sequence, pattern, maxMismatches).Should().BeEmpty(
            "no equal-length window exists when the pattern is longer than the sequence");
        ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches).Should().Be(0,
            "the count surface must agree with the empty listing surface");
        ApproximateMatcher.FindBestMatch(sequence, pattern).Should().BeNull(
            "FindBestMatch returns null when the pattern is longer than the sequence");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: maxDist > len  (and the maxDist = 0 / negative boundaries)
    // ───────────────────────────────────────────────────────────────────

    #region BE — maxMismatches ≥ pattern length, = 0, and negative

    /// <summary>
    /// BE: maxMismatches ≥ pattern.Length is the upper threshold boundary — every
    /// achievable mismatch count is within tolerance, so EVERY equal-length window
    /// qualifies and ALL start positions [0 .. n−m] are reported
    /// (Approximate_Matching_Hamming.md §6.1). Verified against a brute-force
    /// oracle for a threshold equal to the length, well above it, and at int.MaxValue
    /// (no overflow / no IndexOutOfRange). CountApproximateOccurrences must agree.
    /// </summary>
    [TestCase("ACGTACGT", "TTTT", 4, TestName = "FindWithMismatches_MaxDistEqualsLen_AllWindows")]
    [TestCase("ACGTACGT", "TTTT", 100, TestName = "FindWithMismatches_MaxDistAboveLen_AllWindows")]
    [TestCase("ACGTACGT", "TTTT", int.MaxValue, TestName = "FindWithMismatches_MaxDistMaxInt_AllWindows")]
    [TestCase("AAAA", "C", 1, TestName = "FindWithMismatches_SingleBase_MaxDistEqualsLen_AllWindows")]
    public void FindWithMismatches_MaxDistAtLeastPatternLength_MatchesEveryWindow(
        string sequence, string pattern, int maxMismatches)
    {
        int windowCount = sequence.Length - pattern.Length + 1;

        IReadOnlyList<ApproximateMatchResult> results =
            ApproximateMatcher.FindWithMismatches(sequence, pattern, maxMismatches).ToList();

        results.Select(r => r.Position).Should().Equal(Enumerable.Range(0, windowCount),
            "with maxMismatches ≥ pattern length every equal-length window qualifies, in left-to-right order");
        results.Should().OnlyContain(r => r.Distance <= pattern.Length,
            "no reported window can exceed the pattern length in mismatches (INV-01 bound)");
        ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches).Should().Be(windowCount,
            "the count surface must equal the number of equal-length windows (all qualify)");
    }

    /// <summary>
    /// BE: maxMismatches = 0 is the lower threshold boundary — it reduces EXACTLY
    /// to exact matching (only zero-mismatch windows qualify,
    /// Approximate_Matching_Hamming.md §6.1). The reported zero-distance positions
    /// must equal the brute-force set of exact occurrences over a fixed-seed random
    /// text, and every reported result must carry distance 0 with no mismatch
    /// positions. This pins the "approx(d=0) == exact" contract at the boundary.
    /// </summary>
    [Test]
    public void FindWithMismatches_MaxDistZero_ReducesToExactMatch()
    {
        string sequence = RandomDna(400);
        const int k = 5;
        string pattern = sequence.Substring(123, k); // a pattern guaranteed to occur

        var expected = Enumerable.Range(0, sequence.Length - k + 1)
            .Where(i => sequence.Substring(i, k) == pattern)
            .ToList();

        IReadOnlyList<ApproximateMatchResult> results =
            ApproximateMatcher.FindWithMismatches(sequence, pattern, 0).ToList();

        results.Select(r => r.Position).Should().Equal(expected,
            "maxMismatches = 0 reports exactly the exact-match positions (approx(d=0) == exact)");
        results.Should().OnlyContain(r => r.Distance == 0 && r.MismatchPositions.Count == 0,
            "an exact (zero-mismatch) match carries distance 0 and no mismatch positions");
        expected.Should().NotBeEmpty("the planted pattern must occur at least once for a meaningful boundary check");
    }

    /// <summary>
    /// BE: a NEGATIVE maxMismatches is the documented validation gate — it must
    /// throw ArgumentOutOfRangeException (Approximate_Matching_Hamming.md §3.1,
    /// §3.3). Because FindWithMismatches is a lazy iterator, the guard fires on
    /// enumeration, so the assertion materializes the result. Pinned at −1 and
    /// int.MinValue. The same guard governs CountApproximateOccurrences, which
    /// enumerates internally.
    /// </summary>
    [TestCase(-1, TestName = "FindWithMismatches_NegativeMaxDist_MinusOne_Throws")]
    [TestCase(int.MinValue, TestName = "FindWithMismatches_NegativeMaxDist_MinInt_Throws")]
    public void FindWithMismatches_NegativeMaxDist_ThrowsArgumentOutOfRange(int maxMismatches)
    {
        var findAct = () => ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGT", maxMismatches).ToList();
        var countAct = () => ApproximateMatcher.CountApproximateOccurrences("ACGTACGT", "ACGT", maxMismatches);

        findAct.Should().Throw<ArgumentOutOfRangeException>(
            "a negative mismatch budget is rejected by the documented validation gate");
        countAct.Should().Throw<ArgumentOutOfRangeException>(
            "the count surface inherits the same negative-budget validation gate");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: non-DNA characters
    // ───────────────────────────────────────────────────────────────────

    #region MC/INJ — non-DNA / control / non-ASCII characters on the string surface

    /// <summary>
    /// MC: the STRING surface does not enforce the DNA alphabet (only the typed
    /// DnaSequence ctor does), so non-DNA symbols — IUPAC ambiguity codes, digits,
    /// punctuation, gap characters — are matched by plain case-folded character
    /// comparison and must never crash (Approximate_Matching_Hamming.md §5.2).
    /// A pattern of non-DNA chars present verbatim in the text is an EXACT
    /// (distance-0) match at d = 0; a one-symbol difference is a distance-1 match
    /// reported once the budget allows it. Verified against a brute-force oracle.
    /// </summary>
    [Test]
    public void FindWithMismatches_NonDnaCharacters_AreMatchedByPlainComparison()
    {
        const string sequence = "NNN-ACGT-NNN?1234";
        const string pattern = "-NNN"; // occurs verbatim starting at index 8 ("-NNN" before "?1234")

        // d = 0 must find the verbatim occurrence (and only it).
        var exact = ApproximateMatcher.FindWithMismatches(sequence, pattern, 0).ToList();
        exact.Select(r => r.Position).Should().Equal(new[] { 8 },
            "a non-DNA pattern present verbatim is an exact match at its true position");
        exact.Single().Distance.Should().Be(0, "a verbatim non-DNA window has zero Hamming distance");

        // A brute-force oracle over a tolerant budget: every equal-length window
        // whose Hamming distance ≤ 2 must be reported, with the correct distance.
        const int budget = 2;
        var oracle = Enumerable.Range(0, sequence.Length - pattern.Length + 1)
            .Select(i => (Pos: i, Dist: ApproximateMatcher.HammingDistance(
                sequence.Substring(i, pattern.Length), pattern)))
            .Where(x => x.Dist <= budget)
            .ToList();

        var actual = ApproximateMatcher.FindWithMismatches(sequence, pattern, budget).ToList();
        actual.Select(r => (r.Position, r.Distance)).Should().Equal(oracle.Select(x => (x.Pos, x.Dist)),
            "non-DNA windows must be matched/counted exactly as the case-folded Hamming oracle dictates");
    }

    /// <summary>
    /// MC/INJ: control characters, null bytes, and non-ASCII Unicode in the string
    /// surface must be handled as ordinary code points by the case-folded
    /// comparison — no crash, no hang. A pattern equal to the whole text is a single
    /// exact match at position 0; an absent pattern of the same length is no match.
    /// This fuzzes the INJ boundary (null byte / unicode) on the lenient surface.
    /// </summary>
    [Test]
    public void FindWithMismatches_ControlAndUnicodeCharacters_DoNotCrash()
    {
        string sequence = "AC\0GTé\tNNαβ";

        var act = () =>
        {
            // Whole-text pattern → exactly one exact match at 0 on any code points.
            var whole = ApproximateMatcher.FindWithMismatches(sequence, sequence, 0).ToList();
            whole.Select(r => r.Position).Should().Equal(new[] { 0 },
                "the whole text matches itself exactly once at position 0, control/unicode included");
            whole.Single().Distance.Should().Be(0);

            // An absent same-length pattern (all 'Z') → no match at d = 0.
            string absent = new string('Z', sequence.Length);
            ApproximateMatcher.FindWithMismatches(sequence, absent, 0).Should().BeEmpty(
                "a same-length pattern that occurs nowhere yields no exact match");

            // The direct metric over equal-length unicode strings is still symmetric.
            int d = ApproximateMatcher.HammingDistance(sequence, absent);
            d.Should().Be(ApproximateMatcher.HammingDistance(absent, sequence),
                "Hamming distance stays symmetric over control/unicode code points (INV-03)");
        };

        act.Should().NotThrow("control bytes and unicode are ordinary code points to the string surface");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-APPROX-002 — approximate (edit-distance) pattern matching : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PAT-APPROX-002 — approximate (edit-distance) pattern matching

    #region Helper — brute-force oracle replaying the variable-window scan

    /// <summary>
    /// Brute-force oracle for FindWithEdits: replays the EXACT (start i, window
    /// length len) double-loop of the source (ApproximateMatcher.cs lines 134–154)
    /// and the same case-folding, calling the public EditDistance(...) metric on
    /// each candidate window. Returns the ordered (Position, MatchedSequence,
    /// Distance) tuples a correct implementation must yield — including DUPLICATE
    /// positions, since several window lengths can qualify at one start.
    /// </summary>
    private static List<(int Position, string Window, int Distance)> EditScanOracle(
        string sequence, string pattern, int maxEdits)
    {
        var result = new List<(int, string, int)>();
        string seq = sequence.ToUpperInvariant();
        string pat = pattern.ToUpperInvariant();

        int minLen = Math.Max(1, pat.Length - maxEdits);
        int maxLen = pat.Length + maxEdits;

        for (int i = 0; i <= seq.Length - minLen; i++)
        {
            for (int len = minLen; len <= maxLen && i + len <= seq.Length; len++)
            {
                string window = seq.Substring(i, len);
                int distance = ApproximateMatcher.EditDistance(pat, window);
                if (distance <= maxEdits)
                    result.Add((i, window, distance));
            }
        }

        return result;
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: empty strings (empty sequence / empty pattern)
    // ───────────────────────────────────────────────────────────────────

    #region BE — empty strings: empty sequence / pattern → no matches, no crash

    /// <summary>
    /// BE: the empty sequence and the empty pattern are the lower size boundaries.
    /// FindWithEdits guards both with an explicit IsNullOrEmpty check and yields
    /// nothing — it does NOT inherit the exact-match core's "empty pattern ⇒ all
    /// positions" contract (Edit_Distance.md §3.3, §6.1). Every empty combination —
    /// over every maxEdits, including 0 and a large value — must be the defined
    /// empty result, never an IndexOutOfRange from the widened window range. The
    /// guard precedes the negative-maxEdits check, so an empty input with maxEdits
    /// large is still simply empty (the negative case is pinned separately).
    /// </summary>
    [TestCase("", "ACGT", 0, TestName = "FindWithEdits_EmptySequence_NoMatch")]
    [TestCase("", "ACGT", 5, TestName = "FindWithEdits_EmptySequence_LargeMaxEdits_NoMatch")]
    [TestCase("ACGTACGT", "", 0, TestName = "FindWithEdits_EmptyPattern_NoMatch")]
    [TestCase("ACGTACGT", "", 5, TestName = "FindWithEdits_EmptyPattern_LargeMaxEdits_NoMatch")]
    [TestCase("", "", 0, TestName = "FindWithEdits_BothEmpty_NoMatch")]
    public void FindWithEdits_EmptyStrings_AreNoMatch(string sequence, string pattern, int maxEdits)
    {
        ApproximateMatcher.FindWithEdits(sequence, pattern, maxEdits).Should().BeEmpty(
            "an empty sequence or pattern is guarded by the source and yields no edit-distance matches");
    }

    /// <summary>
    /// BE/MC: null sequence/pattern on the string surface is treated the same as
    /// empty by the IsNullOrEmpty guard — it yields nothing rather than throwing a
    /// NullReferenceException (Edit_Distance.md §3.1, §3.3). Because FindWithEdits
    /// is a lazy iterator, the guard runs on enumeration, so the assertions
    /// materialize the result. Pinned so the lenient string surface cannot drift
    /// into a raw dereference.
    /// </summary>
    [Test]
    public void FindWithEdits_NullStrings_AreNoMatchNotCrash()
    {
        var findNullSeq = () => ApproximateMatcher.FindWithEdits((string)null!, "ACGT", 1).ToList();
        var findNullPat = () => ApproximateMatcher.FindWithEdits("ACGT", null!, 1).ToList();

        findNullSeq.Should().NotThrow("a null sequence is guarded as empty, never dereferenced");
        findNullPat.Should().NotThrow("a null pattern is guarded as empty, never dereferenced");

        findNullSeq().Should().BeEmpty("a null sequence yields no matches");
        findNullPat().Should().BeEmpty("a null pattern yields no matches");
    }

    /// <summary>
    /// MC/INJ: the direct EditDistance metric has an explicit null gate
    /// (ArgumentNullException), NOT a NullReferenceException (Edit_Distance.md §3.1,
    /// §3.3; ApproximateMatcher.cs lines 197–198). Pinned in both orientations,
    /// including null-vs-non-empty.
    /// </summary>
    [Test]
    public void EditDistance_NullArguments_ThrowArgumentNullException()
    {
        var nullFirst = () => ApproximateMatcher.EditDistance(null!, "ACGT");
        var nullSecond = () => ApproximateMatcher.EditDistance("ACGT", null!);
        var bothNull = () => ApproximateMatcher.EditDistance(null!, null!);

        nullFirst.Should().Throw<ArgumentNullException>("a null first string hits the null gate");
        nullSecond.Should().Throw<ArgumentNullException>("a null second string hits the null gate");
        bothNull.Should().Throw<ArgumentNullException>("two null strings still hit the null gate");
    }

    /// <summary>
    /// BE: the direct EditDistance metric over the empty string is the canonical
    /// length-difference boundary — "" vs "abc" is 3 insertions/deletions, in both
    /// orientations, and "" vs "" is 0 (Edit_Distance.md §6.1, INV-01/INV-02). The
    /// empty string here is a metric input, distinct from the FindWithEdits empty
    /// guard above. No crash on the zero-length DP row.
    /// </summary>
    [TestCase("", "ABC", 3, TestName = "EditDistance_EmptyVsThree_IsThree")]
    [TestCase("ABC", "", 3, TestName = "EditDistance_ThreeVsEmpty_IsThree")]
    [TestCase("", "", 0, TestName = "EditDistance_BothEmpty_IsZero")]
    [TestCase("ACGT", "ACGT", 0, TestName = "EditDistance_Identical_IsZero")]
    public void EditDistance_EmptyStringBoundary_IsLengthDifference(string s1, string s2, int expected)
    {
        ApproximateMatcher.EditDistance(s1, s2).Should().Be(expected,
            "edit distance over an empty operand is the other string's length (INV-02 length bound)");
    }

    /// <summary>
    /// MC: a null DnaSequence on the typed wrapper adds NO null guard and
    /// dereferences sequence.Sequence, so it throws — the documented strict-surface
    /// behavior (Edit_Distance.md §3.1, §3.3, §6.1). Pinned against the lenient
    /// string surface (which guards null as empty above) so the split is explicit.
    /// Because FindWithEdits is lazy, the dereference happens on enumeration.
    /// </summary>
    [Test]
    public void FindWithEdits_NullDnaSequence_Throws()
    {
        var act = () => ApproximateMatcher.FindWithEdits((DnaSequence)null!, "ACGT", 1).ToList();

        act.Should().Throw<NullReferenceException>(
            "the typed wrapper dereferences sequence.Sequence with no null guard (documented strict surface)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: maxEdits negative (the documented validation gate)
    // ───────────────────────────────────────────────────────────────────

    #region BE — maxEdits negative → ArgumentOutOfRangeException; maxEdits = 0 → exact

    /// <summary>
    /// BE: a NEGATIVE maxEdits is the documented validation gate — it must throw
    /// ArgumentOutOfRangeException (Edit_Distance.md §3.3, §6.1;
    /// ApproximateMatcher.cs lines 124–125). Because FindWithEdits is a lazy
    /// iterator AND the IsNullOrEmpty guard precedes the negative check, the
    /// assertion uses NON-empty inputs and materializes the result so the guard
    /// fires. Pinned at −1 and int.MinValue (no overflow on the extreme value).
    /// </summary>
    [TestCase(-1, TestName = "FindWithEdits_NegativeMaxEdits_MinusOne_Throws")]
    [TestCase(int.MinValue, TestName = "FindWithEdits_NegativeMaxEdits_MinInt_Throws")]
    public void FindWithEdits_NegativeMaxEdits_ThrowsArgumentOutOfRange(int maxEdits)
    {
        var act = () => ApproximateMatcher.FindWithEdits("ACGTACGT", "ACGT", maxEdits).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a negative edit budget is rejected by the documented validation gate")
            .Which.ParamName.Should().Be("maxEdits", "the gate names the offending parameter");
    }

    /// <summary>
    /// BE: maxEdits = 0 is the lower threshold boundary. With maxEdits = 0 the
    /// window-length range collapses to exactly the pattern length (minLen ==
    /// maxLen == pat.Length), so only equal-length windows are scanned and only
    /// zero-distance ones qualify — FindWithEdits reduces EXACTLY to exact matching
    /// (Edit_Distance.md §3.3). The reported positions must equal the brute-force
    /// set of exact occurrences over a fixed-seed random text, each at distance 0,
    /// with no duplicate positions (a single window length is checked per start).
    /// </summary>
    [Test]
    public void FindWithEdits_MaxEditsZero_ReducesToExactMatch()
    {
        string sequence = RandomDna(400);
        const int k = 6;
        string pattern = sequence.Substring(150, k); // a pattern guaranteed to occur

        var expected = Enumerable.Range(0, sequence.Length - k + 1)
            .Where(i => sequence.Substring(i, k) == pattern)
            .ToList();

        IReadOnlyList<ApproximateMatchResult> results =
            ApproximateMatcher.FindWithEdits(sequence, pattern, 0).ToList();

        results.Select(r => r.Position).Should().Equal(expected,
            "maxEdits = 0 reports exactly the exact-match positions (approx(d=0) == exact)");
        results.Should().OnlyContain(r => r.Distance == 0,
            "an exact (zero-edit) match carries distance 0");
        results.Select(r => r.Position).Should().OnlyHaveUniqueItems(
            "with maxEdits = 0 only the pattern-length window is scanned, so no start repeats");
        expected.Should().NotBeEmpty("the planted pattern must occur at least once for a meaningful boundary check");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: maxEdits > seq length  (every start position matches)
    // ───────────────────────────────────────────────────────────────────

    #region BE — maxEdits ≥ sequence length → every window reachable

    /// <summary>
    /// BE: maxEdits ≥ sequence length is the upper threshold boundary — every
    /// pattern is within that many edits of every reachable window, so EVERY start
    /// position [0 .. n−minLen] yields at least one match and the widened window
    /// range (pat.Length + maxEdits, far past the text) does NOT over-read or crash
    /// (Edit_Distance.md §2.4 INV-03, §6.1). Verified against the variable-window
    /// oracle for the exact (position, window, distance) multiset, then the "every
    /// reachable start matches" property is checked at a finite, non-overflowing
    /// budget. (The int.MaxValue case is fuzzed separately for no-crash / no-hang —
    /// the source's pat.Length + maxEdits window arithmetic overflows there, which
    /// is a DEFINED empty result, not a corruption.)
    /// </summary>
    [TestCase("ACGTACGT", "GGG", 8, TestName = "FindWithEdits_MaxEditsEqualsLen_AllReachableStarts")]
    [TestCase("ACGTACGT", "GGG", 100, TestName = "FindWithEdits_MaxEditsAboveLen_AllReachableStarts")]
    [TestCase("AAAA", "C", 4, TestName = "FindWithEdits_SingleBasePattern_MaxEditsLen_AllStarts")]
    public void FindWithEdits_MaxEditsAtLeastSequenceLength_MatchesEveryReachableStart(
        string sequence, string pattern, int maxEdits)
    {
        var oracle = EditScanOracle(sequence, pattern, maxEdits);

        IReadOnlyList<ApproximateMatchResult> results =
            ApproximateMatcher.FindWithEdits(sequence, pattern, maxEdits).ToList();

        results.Select(r => (r.Position, r.MatchedSequence, r.Distance))
            .Should().Equal(oracle,
                "the variable-window scan must yield exactly the oracle's (position, window, distance) sequence");

        // Every reachable start position is covered when the budget swamps the length.
        // minLen collapses to 1 once maxEdits ≥ pat.Length, so starts span [0 .. n-1].
        int minLen = Math.Max(1, pattern.Length - maxEdits);
        var expectedStarts = Enumerable.Range(0, sequence.Length - minLen + 1).ToList();
        results.Select(r => r.Position).Distinct().OrderBy(p => p).Should().Equal(expectedStarts,
            "with an edit budget at least the sequence length every reachable start position matches");

        results.Should().OnlyContain(r => r.Distance <= Math.Max(pattern.Length, r.MatchedSequence.Length),
            "no reported edit distance can exceed max(pattern, window) length (INV-03)");
    }

    /// <summary>
    /// BE/OVF: maxEdits = int.MaxValue is the extreme upper boundary. The source
    /// computes the window range as pat.Length ± maxEdits, which overflows int at
    /// this value — the fuzz contract is that this must NOT crash, hang, or corrupt:
    /// enumeration completes and the result is whatever the (overflowed but
    /// deterministic) window arithmetic yields, mirrored exactly by the oracle that
    /// replays the same arithmetic. Pinned so the overflow path stays a defined,
    /// reproducible outcome rather than an IndexOutOfRange or infinite loop.
    /// </summary>
    [Test]
    public void FindWithEdits_MaxEditsMaxInt_DoesNotCrashAndMatchesOracle()
    {
        const string sequence = "ACGTACGT";
        const string pattern = "GGG";

        var act = () =>
            ApproximateMatcher.FindWithEdits(sequence, pattern, int.MaxValue)
                .Select(r => (r.Position, r.MatchedSequence, r.Distance))
                .ToList();

        act.Should().NotThrow("an int.MaxValue edit budget must not crash or hang on the widened window range");
        act().Should().Equal(EditScanOracle(sequence, pattern, int.MaxValue),
            "the int.MaxValue scan yields the same deterministic result as the overflow-faithful oracle");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Fuzz target: non-DNA characters
    // ───────────────────────────────────────────────────────────────────

    #region MC/INJ — non-DNA / control / non-ASCII characters on the string surface

    /// <summary>
    /// MC: the STRING surface does not enforce the DNA alphabet (only the typed
    /// DnaSequence ctor does), so non-DNA symbols — IUPAC ambiguity codes, gap
    /// characters, digits, punctuation — flow through the case-folded variable-
    /// window scan and must never crash (Edit_Distance.md §5.2). A non-DNA pattern
    /// present verbatim is an exact (distance-0) edit match; an off-by-one-edit
    /// neighbour is reported once the budget allows. Verified against the
    /// variable-window oracle so the duplicate-position semantics are pinned too.
    /// </summary>
    [Test]
    public void FindWithEdits_NonDnaCharacters_AreMatchedByVariableWindowScan()
    {
        const string sequence = "NNN-ACGT-NNN?1234";
        const string pattern = "-NNN"; // occurs verbatim starting at index 8

        // d = 0: the verbatim occurrence is found at distance 0.
        var exact = ApproximateMatcher.FindWithEdits(sequence, pattern, 0).ToList();
        exact.Should().Contain(r => r.Position == 8 && r.Distance == 0,
            "a non-DNA pattern present verbatim is an exact (zero-edit) match at its true position");
        exact.Select(r => (r.Position, r.MatchedSequence, r.Distance))
            .Should().Equal(EditScanOracle(sequence, pattern, 0),
                "the d=0 scan over non-DNA symbols must equal the variable-window oracle");

        // A tolerant budget over non-DNA symbols must match the oracle exactly,
        // including duplicate start positions from multiple qualifying window lengths.
        const int budget = 2;
        ApproximateMatcher.FindWithEdits(sequence, pattern, budget)
            .Select(r => (r.Position, r.MatchedSequence, r.Distance))
            .Should().Equal(EditScanOracle(sequence, pattern, budget),
                "non-DNA windows must be matched exactly as the variable-window edit oracle dictates");
    }

    /// <summary>
    /// MC/INJ: control characters, null bytes, and non-ASCII Unicode in the string
    /// surface must flow through the variable-window scan as ordinary code points —
    /// no crash, no hang, no over-read. A pattern equal to the whole text is an
    /// exact match at position 0; an absent same-length pattern still yields no
    /// zero-edit match. This fuzzes the INJ boundary (null byte / unicode) on the
    /// lenient surface and pins the whole run against the oracle.
    /// </summary>
    [Test]
    public void FindWithEdits_ControlAndUnicodeCharacters_DoNotCrash()
    {
        string sequence = "AC\0GT\tNNZZ"; // ASCII-folding-safe code points incl. null byte

        var act = () =>
        {
            // Whole-text pattern → exactly one exact (zero-edit) match at 0 at d = 0.
            var whole = ApproximateMatcher.FindWithEdits(sequence, sequence, 0).ToList();
            whole.Select(r => (r.Position, r.MatchedSequence, r.Distance))
                .Should().Equal(EditScanOracle(sequence, sequence, 0),
                    "the whole-text edit scan over control/unicode must equal the oracle");
            whole.Should().ContainSingle(r => r.Position == 0 && r.Distance == 0,
                "the whole text matches itself exactly once at position 0, control bytes included");

            // An absent same-length pattern (all 'Z' over the folded length) → no exact match.
            string absent = new string('Z', sequence.Length + 3); // longer than text → unreachable at d=0
            ApproximateMatcher.FindWithEdits(sequence, absent, 0).Should().BeEmpty(
                "a pattern longer than the text cannot match at zero edits (minLen exceeds the text)");

            // A tolerant budget still tracks the oracle exactly over the odd code points.
            ApproximateMatcher.FindWithEdits(sequence, "NN", 1)
                .Select(r => (r.Position, r.MatchedSequence, r.Distance))
                .Should().Equal(EditScanOracle(sequence, "NN", 1),
                    "a tolerant edit scan over control/unicode must match the variable-window oracle");
        };

        act.Should().NotThrow("control bytes and unicode are ordinary code points to the edit-distance scan");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-IUPAC-001 — IUPAC (degenerate) pattern matching : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PAT-IUPAC-001 — IUPAC (degenerate) pattern matching

    #region Sanity — valid ambiguity codes match exactly their allowed base set

    /// <summary>
    /// Positive theory sanity (NOT a fuzz target, the anchor the fuzz tests defend):
    /// a valid ambiguity code matches EXACTLY the bases in its degenerate set and no
    /// others — R matches A and G but not C/T, Y matches C and T, N matches all four
    /// — on BOTH surfaces (IupacHelper.MatchesIupac and the MotifFinder scan over a
    /// single-base sequence). IUPAC_Degenerate_Matching.md §2.2 (code table),
    /// INV-01. This pins the contract the malformed-input tests must not weaken.
    /// </summary>
    [Test]
    public void Iupac_ValidAmbiguityCodes_MatchExactlyTheirBaseSet()
    {
        // Per-character oracle: MatchesIupac membership for every (code, base) pair.
        var codeSets = new Dictionary<char, string>
        {
            ['A'] = "A", ['C'] = "C", ['G'] = "G", ['T'] = "T",
            ['R'] = "AG", ['Y'] = "CT", ['S'] = "GC", ['W'] = "AT",
            ['K'] = "GT", ['M'] = "AC", ['B'] = "CGT", ['D'] = "AGT",
            ['H'] = "ACT", ['V'] = "ACG", ['N'] = "ACGT",
        };

        foreach (var (code, allowed) in codeSets)
        {
            foreach (char baseChar in "ACGT")
            {
                bool expected = allowed.Contains(baseChar);

                IupacHelper.MatchesIupac(baseChar, code).Should().Be(expected,
                    $"'{code}' must match base '{baseChar}' iff '{baseChar}' is in its degenerate set");

                // The MotifFinder scan over a single-base sequence must agree.
                var matches = MotifFinder.FindDegenerateMotif(new DnaSequence(baseChar.ToString()), code.ToString())
                    .ToList();
                (matches.Count == 1).Should().Be(expected,
                    $"the scan of motif '{code}' over '{baseChar}' must match iff '{baseChar}' is in the set");
            }
        }

        // The canonical R-vs-{A,G} / Y-vs-{C,T} / N-vs-all spot checks, explicitly.
        IupacHelper.MatchesIupac('A', 'R').Should().BeTrue("R = puRine includes A");
        IupacHelper.MatchesIupac('G', 'R').Should().BeTrue("R = puRine includes G");
        IupacHelper.MatchesIupac('C', 'R').Should().BeFalse("R = puRine excludes C");
        IupacHelper.MatchesIupac('T', 'R').Should().BeFalse("R = puRine excludes T");
        foreach (char b in "ACGT")
            IupacHelper.MatchesIupac(b, 'N').Should().BeTrue("N matches every base");
    }

    /// <summary>
    /// Positive sanity: a real degenerate motif over a known sequence reports
    /// EXACTLY the theory-correct windows (every base in its per-position code set),
    /// and a lowercase motif folds to the identical result (ToUpperInvariant before
    /// validation). IUPAC_Degenerate_Matching.md §2.2, §6.1, INV-01.
    /// </summary>
    [Test]
    public void Iupac_DegenerateMotif_MatchesTheoryCorrectWindowsAndFoldsCase()
    {
        var seq = new DnaSequence("AGCTAACTGT");
        const string motif = "RY"; // R=A/G then Y=C/T

        var positions = MotifFinder.FindDegenerateMotif(seq, motif).Select(m => m.Position).ToList();

        // Brute-force oracle: window [i,i+1] where seq[i]∈{A,G} and seq[i+1]∈{C,T}.
        string s = seq.Sequence;
        var expected = Enumerable.Range(0, s.Length - 1)
            .Where(i => "AG".Contains(s[i]) && "CT".Contains(s[i + 1]))
            .ToList();

        positions.Should().Equal(expected, "an 'RY' motif matches exactly the puRine-pYrimidine windows");
        expected.Should().NotBeEmpty("the chosen sequence must contain at least one RY window for a real check");

        MotifFinder.FindDegenerateMotif(seq, "ry").Select(m => m.Position).Should().Equal(positions,
            "a lowercase motif is uppercased before validation, so case must not change the result");
    }

    #endregion

    #region MC — invalid IUPAC codes (protein-only letters) → ArgumentException

    /// <summary>
    /// MC: an amino-acid-only protein letter that is NOT a standard IUPAC nucleotide
    /// code (E, F, I, L, P, Q, Z — and J/O/U/X which are not nucleotide codes) must
    /// be rejected by the documented validation gate with an ArgumentException naming
    /// "motif" (ValidateIupacPattern, MotifFinder.cs lines 70–82;
    /// IUPAC_Degenerate_Matching.md §3.3, §6.1) — NEVER a KeyNotFoundException from
    /// the IupacCodes dictionary lookup. Because the scan is a lazy iterator, the
    /// guard fires on enumeration, so the assertion materializes the result.
    /// </summary>
    [TestCase("Z", TestName = "Iupac_InvalidCode_Z_Throws")]
    [TestCase("E", TestName = "Iupac_InvalidCode_E_Throws")]
    [TestCase("J", TestName = "Iupac_InvalidCode_J_Throws")]
    [TestCase("O", TestName = "Iupac_InvalidCode_O_Throws")]
    [TestCase("U", TestName = "Iupac_InvalidCode_U_Throws")]
    [TestCase("X", TestName = "Iupac_InvalidCode_X_Throws")]
    [TestCase("ACGTZ", TestName = "Iupac_InvalidCode_ProteinAfterValidPrefix_Throws")]
    public void Iupac_InvalidProteinLetterMotif_ThrowsArgumentException(string motif)
    {
        var seq = new DnaSequence("ACGTACGTACGT");

        var act = () => MotifFinder.FindDegenerateMotif(seq, motif).ToList();

        act.Should().Throw<ArgumentException>(
                "a non-IUPAC protein letter is rejected by the documented validation gate, not looked up blindly")
            .And.ParamName.Should().Be("motif", "the validation gate names the offending parameter");
    }

    /// <summary>
    /// MC: a motif whose FIRST window already fails to match still validates the
    /// WHOLE pattern — the invalid code at the end is NOT skipped just because the
    /// match short-circuits earlier. Validation precedes the scan (MotifFinder.cs
    /// line 97), so 'ACGTZ' over a non-matching sequence still throws rather than
    /// silently yielding no matches and hiding the malformed code.
    /// </summary>
    [Test]
    public void Iupac_InvalidCode_AfterNonMatchingPrefix_StillThrows()
    {
        var seq = new DnaSequence("TTTTTTTT"); // the 'A...' prefix never matches

        var act = () => MotifFinder.FindDegenerateMotif(seq, "ACGTZ").ToList();

        act.Should().Throw<ArgumentException>(
            "validation runs over the whole motif before scanning, so a trailing invalid code is never skipped");
    }

    /// <summary>
    /// MC: the decision-table primitive IupacHelper.MatchesIupac rejects an invalid
    /// code with ArgumentOutOfRangeException (the switch default, IupacHelper.cs
    /// lines 33–35) — a DIFFERENT documented exception type from the MotifFinder
    /// surface's ArgumentException, NEVER a KeyNotFoundException. Pinned so the two
    /// surfaces' rejection contracts cannot silently converge or degrade to a crash.
    /// </summary>
    [TestCase('Z', TestName = "MatchesIupac_InvalidCode_Z_ThrowsOutOfRange")]
    [TestCase('E', TestName = "MatchesIupac_InvalidCode_E_ThrowsOutOfRange")]
    [TestCase('5', TestName = "MatchesIupac_InvalidCode_Digit_ThrowsOutOfRange")]
    [TestCase('α', TestName = "MatchesIupac_InvalidCode_Unicode_ThrowsOutOfRange")]
    public void MatchesIupac_InvalidCode_ThrowsArgumentOutOfRange(char code)
    {
        var act = () => IupacHelper.MatchesIupac('A', code);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "the decision-table default rejects unknown codes, never KeyNotFound")
            .And.ParamName.Should().Be("iupacCode", "the gate names the offending code parameter");
    }

    #endregion

    #region INJ — digits in the pattern → ArgumentException

    /// <summary>
    /// INJ: digits (0–9) injected into the motif are not IUPAC codes and must be
    /// rejected with an ArgumentException — never a KeyNotFound or a silent match.
    /// Pinned for a bare digit, a digit embedded after a valid prefix, and an
    /// all-digit motif. Lazy iterator → assertion materializes the result.
    /// </summary>
    [TestCase("1", TestName = "Iupac_Digit_Bare_Throws")]
    [TestCase("0", TestName = "Iupac_Digit_Zero_Throws")]
    [TestCase("AC1GT", TestName = "Iupac_Digit_Embedded_Throws")]
    [TestCase("1234", TestName = "Iupac_Digit_AllDigits_Throws")]
    public void Iupac_DigitsInMotif_ThrowArgumentException(string motif)
    {
        var seq = new DnaSequence("ACGTACGTACGT");

        var act = () => MotifFinder.FindDegenerateMotif(seq, motif).ToList();

        act.Should().Throw<ArgumentException>(
                "a digit is not a valid IUPAC code and is rejected by the documented validation gate")
            .And.ParamName.Should().Be("motif");
    }

    #endregion

    #region INJ — Unicode, control bytes, punctuation in the pattern → ArgumentException

    /// <summary>
    /// INJ: non-ASCII Unicode, a null byte, a tab, punctuation, and a gap character
    /// injected into the motif are all NON-IUPAC and must be rejected with an
    /// ArgumentException — never a KeyNotFound from the dictionary lookup, never an
    /// IndexOutOfRange, never a hang. Greek/accented letters survive ToUpperInvariant
    /// as non-ASCII code points and still hit the default validation branch. Pinned
    /// across a spread of injection code points; the lazy scan materializes on ToList.
    /// </summary>
    [TestCase("α", TestName = "Iupac_Unicode_Greek_Throws")]
    [TestCase("ñ", TestName = "Iupac_Unicode_Accented_Throws")]
    [TestCase("AC\0GT", TestName = "Iupac_NullByte_Embedded_Throws")]
    [TestCase("AC\tGT", TestName = "Iupac_Tab_Embedded_Throws")]
    [TestCase("A-C", TestName = "Iupac_GapChar_Throws")]
    [TestCase("A.C", TestName = "Iupac_Punctuation_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "Iupac_FullWidthLetters_Throws")]
    public void Iupac_UnicodeAndControlAndPunctuationMotif_ThrowArgumentException(string motif)
    {
        var seq = new DnaSequence("ACGTACGTACGT");

        var act = () => MotifFinder.FindDegenerateMotif(seq, motif).ToList();

        act.Should().Throw<ArgumentException>(
                "Unicode, control bytes, and punctuation are not IUPAC codes and are rejected deterministically")
            .And.ParamName.Should().Be("motif");
    }

    /// <summary>
    /// INJ/MC: a battery of malformed motifs must NEVER fault with an undisciplined
    /// runtime exception (KeyNotFoundException from the IupacCodes lookup,
    /// IndexOutOfRangeException, NullReferenceException) and must NEVER hang — each
    /// must resolve to EITHER the documented ArgumentException OR a defined empty/
    /// match result. This is the core fuzzing guarantee for the malformed-pattern
    /// surface (docs/ADVANCED_TESTING_CHECKLIST.md §8). The empty motif is the one
    /// "malformed" input that is a DEFINED no-match (the IsNullOrEmpty guard), so it
    /// is allowed through without throwing; every other entry must be rejected.
    /// </summary>
    [Test]
    public void Iupac_MalformedMotifBattery_NeverFaultsUndisciplined()
    {
        var seq = new DnaSequence("ACGTACGTACGTACGT");

        var alwaysRejected = new[]
        {
            "Z", "EE", "ACGTZ", "1", "AC1GT", "α", "ñ", "AC\0GT", "AC\tGT",
            "A-C", "A.C", "ＡＣＧＴ", "RYZ", "NNNNNNQ", "?", "*", "𝓐",
        };

        foreach (string motif in alwaysRejected)
        {
            var act = () => MotifFinder.FindDegenerateMotif(seq, motif).ToList();
            act.Should().Throw<ArgumentException>(
                $"malformed motif '{motif}' must be rejected by the documented gate, never an undisciplined fault");
        }

        // The empty motif is the one DEFINED no-match (guarded), not an error.
        var emptyAct = () => MotifFinder.FindDegenerateMotif(seq, string.Empty).ToList();
        emptyAct.Should().NotThrow("the empty motif is a defined no-match, not a malformed-input error");
        MotifFinder.FindDegenerateMotif(seq, string.Empty).Should().BeEmpty(
            "the empty motif yields no matches by the explicit guard");
    }

    #endregion

    #endregion
}
