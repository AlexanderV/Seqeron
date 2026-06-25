using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Matching area — PAT-APPROX-003, "Best Match and Frequency
/// Analysis": the THIRD approximate (Hamming) pattern-matching unit on
/// ApproximateMatcher, distinct from PAT-APPROX-001 (FindWithMismatches /
/// HammingDistance, row 9) and PAT-APPROX-002 (FindWithEdits / EditDistance,
/// row 10). Checklist: docs/checklists/03_FUZZING.md, row 174.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds boundary and malformed inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no infinite loop, no state
/// corruption, and no *unhandled* runtime exception (IndexOutOfRangeException,
/// a negative-length Substring, NullReferenceException, DivideByZero). Every
/// input must yield EITHER a well-defined, theory-correct result, OR a
/// *documented, intentional* validation exception (here:
/// ArgumentOutOfRangeException on a negative d / non-positive k). A raw runtime
/// exception, a wrong-result corruption, or a hang on boundary input is a bug,
/// not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-APPROX-003 — approximate (Hamming) best-match + frequency (Matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 174.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the three documented boundaries from row 174:
///       – pattern > text — a pattern STRICTLY LONGER than the text. No
///         equal-length window exists, so FindBestMatch → null and
///         CountApproximateOccurrences → 0, even when the mismatch budget is huge.
///         The internal `for i in 0 .. seq.Length - pat.Length` must NOT run with
///         a negative bound, and `seq.Substring(i, pat.Length)` must NOT throw
///         IndexOutOfRange / negative-length. FindFrequentKmersWithMismatches with
///         k > sequence length must likewise enumerate no window (and therefore
///         hit the documented empty-tally path).
///       – empty — empty sequence and/or empty pattern. FindBestMatch and
///         CountApproximateOccurrences guard with IsNullOrEmpty → null / 0 (NOT a
///         crash, NOT "matches everywhere"). FindFrequentKmersWithMismatches on an
///         empty sequence yield-breaks BEFORE the `counts.Values.Max()` call, so it
///         never throws InvalidOperationException ("Sequence contains no elements")
///         on the empty tally — pinned as the degenerate guard.
///       – exact present — a pattern that occurs EXACTLY (Hamming distance 0) is a
///         special case of an approximate occurrence and MUST be found at the
///         correct 0-based position(s) even when the allowed mismatch budget is
///         strictly > 0. FindBestMatch must report distance 0 / IsExact at the
///         leftmost exact position; CountApproximateOccurrences must count it.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The PAT-APPROX-003 contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PAT-APPROX-003 is the substitutions-only (Hamming) best-match + frequency
/// family on ApproximateMatcher:
///   • FindBestMatch(string sequence, string pattern)
///   • CountApproximateOccurrences(string sequence, string pattern, int maxMismatches)
///   • FindFrequentKmersWithMismatches(string sequence, int k, int d)
/// Documented contract (docs/algorithms/Pattern_Matching/
/// Frequent_Words_With_Mismatches.md §2.2, §2.4, §3.1, §3.2, §3.3, §6.1; INV-01..
/// INV-05; sources: Compeau &amp; Pevzner ch.1, ROSALIND BA1H/BA1I/BA1N):
///   • Approximate occurrence: position i counts iff
///     HammingDistance(P, T[i..i+m]) ≤ d (INV-02, BA1H). Count_d(T,P) is the number
///     of such positions (§2.2). Count_0 = exact count; Count_d ≥ exact count
///     (INV-01).
///   • FindBestMatch returns the leftmost minimum-Hamming-distance equal-length
///     window (INV-04 distance = min over windows; INV-05 leftmost tie-break),
///     `ApproximateMatchResult?`; null when either input is empty OR the pattern is
///     longer than the sequence (§3.2, §6.1). IsExact iff that minimum is 0.
///   • FindFrequentKmersWithMismatches returns ALL k-mers over {A,C,G,T}
///     maximizing Count_d (every tie, INV-03); the result may include k-mers that
///     never occur exactly (counting is over the d-neighborhood / Hamming ball,
///     §2.4). Output order of tied maxima is unspecified — callers compare as a
///     SET (§4.2).
///   • Inputs are upper-cased (case-insensitive). Empty sequence/pattern, or a
///     pattern longer than the sequence, yields the empty result / 0 / null.
///     maxMismatches &lt; 0 and d &lt; 0 → ArgumentOutOfRangeException; k ≤ 0 →
///     ArgumentOutOfRangeException (§3.1, §3.3, §6.1).
///
/// A well-formed FindBestMatch result is checked by <see cref="AssertWellFormedBestMatch"/>:
/// the reported Position is in-bounds, the MatchedSequence is exactly the
/// pat.Length window at that position, the Distance equals the recomputed Hamming
/// distance of that window (so the reported distance can never drift from the
/// actual window), the MismatchPositions are the exact differing indices, and
/// IsExact ⇔ Distance == 0.
///
/// All inputs here are ASCII DNA / boundary strings; randomness (where used) is
/// from a locally fixed-seed Random so the fuzz fodder is fully reproducible and
/// adding tests here cannot perturb any other fixture.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PatternApprox003FuzzTests
{
    #region Helpers

    /// <summary>
    /// Deterministic, LOCALLY-scoped RNG — seed fixed so generated fuzz inputs are
    /// reproducible. Private to this fixture and never shared/mutated across
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
    /// Independent Hamming-distance oracle over equal-length, already-uppercased
    /// strings — recomputed here, NOT borrowed from the unit under test, so the
    /// well-formed-result check cannot echo the implementation's own arithmetic.
    /// </summary>
    private static int HammingOracle(string a, string b)
    {
        a.Length.Should().Be(b.Length, "Hamming distance is defined only on equal-length strings");
        int d = 0;
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) d++;
        return d;
    }

    /// <summary>
    /// Asserts a FindBestMatch result is well-formed against the documented
    /// contract: in-bounds position, MatchedSequence == the pat.Length window at
    /// that position (case-folded), reported Distance == recomputed Hamming
    /// distance of that window, MismatchPositions == the exact differing indices,
    /// and IsExact ⇔ Distance == 0 (INV-04). All inputs are compared upper-cased
    /// because the unit upper-cases internally (§3.3).
    /// </summary>
    private static void AssertWellFormedBestMatch(string sequence, string pattern, ApproximateMatchResult r)
    {
        var seq = sequence.ToUpperInvariant();
        var pat = pattern.ToUpperInvariant();

        r.Position.Should().BeGreaterThanOrEqualTo(0, "a match position is never negative");
        (r.Position + pat.Length).Should().BeLessThanOrEqualTo(seq.Length,
            "the reported window must lie entirely inside the sequence (no over-read)");

        string window = seq.Substring(r.Position, pat.Length);
        r.MatchedSequence.Should().Be(window,
            "MatchedSequence must be exactly the pat-length window at the reported position");

        int expectedDistance = HammingOracle(pat, window);
        r.Distance.Should().Be(expectedDistance,
            "the reported distance must equal the recomputed Hamming distance of the window (INV-04)");
        r.Distance.Should().BeGreaterThanOrEqualTo(0, "Hamming distance is a non-negative count");

        var expectedMismatchPositions = Enumerable.Range(0, pat.Length).Where(j => window[j] != pat[j]).ToList();
        r.MismatchPositions.Should().Equal(expectedMismatchPositions,
            "MismatchPositions must be the exact zero-based differing indices");

        r.IsExact.Should().Be(r.Distance == 0, "IsExact ⇔ Distance == 0 (INV-04)");
        r.MismatchType.Should().Be(MismatchType.Substitution, "Hamming-based matching is substitutions only");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-APPROX-003 — best match + frequency : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region BE — pattern > text: no equal-length window → null / 0 / empty, no over-read

    /// <summary>
    /// BE — pattern STRICTLY LONGER than the text. The internal scan bounds the
    /// start index by `seq.Length - pat.Length`, which is negative here, so the
    /// loop must not run and `seq.Substring(i, pat.Length)` must never be reached
    /// with an out-of-range / negative length. FindBestMatch → null,
    /// CountApproximateOccurrences → 0 — even with a huge mismatch budget, because
    /// no equal-length window exists at all (§3.2, §6.1, INV-04).
    /// </summary>
    [TestCase("ACGT", "ACGTA", 0, TestName = "PatternLongerByOne_ZeroBudget")]
    [TestCase("ACGT", "ACGTA", 100, TestName = "PatternLongerByOne_HugeBudget")]
    [TestCase("ACGT", "ACGTACGT", int.MaxValue, TestName = "PatternMuchLonger_MaxIntBudget")]
    [TestCase("A", "AC", 5, TestName = "PatternLongerThanSingleBase")]
    [TestCase("A", "AAAAAAAA", 3, TestName = "OneBaseText_LongPattern")]
    public void PatternLongerThanText_BestMatchNull_CountZero(
        string sequence, string pattern, int maxMismatches)
    {
        ApproximateMatcher.FindBestMatch(sequence, pattern).Should().BeNull(
            "no equal-length window exists when the pattern is longer than the text");

        ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches).Should().Be(0,
            "Count_d is 0 when the pattern is longer than the text, regardless of the budget");
    }

    /// <summary>
    /// BE — FindFrequentKmersWithMismatches with k STRICTLY GREATER than the
    /// sequence length: no length-k window can be cut, so the per-window loop
    /// `for i in 0 .. seq.Length - k` never executes. The method must therefore
    /// return an EMPTY tally WITHOUT reaching `counts.Values.Max()` on an empty
    /// dictionary (which would throw InvalidOperationException). Pinned as the
    /// degenerate empty-window path (no crash, no over-read).
    /// </summary>
    [TestCase("ACGT", 5, 0, TestName = "KOneAboveLength_ZeroD")]
    [TestCase("ACGT", 8, 2, TestName = "KFarAboveLength_NonZeroD")]
    [TestCase("A", 2, 1, TestName = "OneBaseSeq_KTwo")]
    public void FrequentKmers_KLongerThanSequence_IsEmptyNoCrash(string sequence, int k, int d)
    {
        Func<List<(string Kmer, int Count)>> act = () =>
            ApproximateMatcher.FindFrequentKmersWithMismatches(sequence, k, d).ToList();

        act.Should().NotThrow("no length-k window can be cut, so the scan must not run, not crash");
        act().Should().BeEmpty("there are no k-mers to tally when k exceeds the sequence length");
    }

    #endregion

    #region BE — empty: empty sequence / pattern → null / 0 / empty, no crash / DivideByZero

    /// <summary>
    /// BE — empty sequence and/or empty pattern on FindBestMatch and
    /// CountApproximateOccurrences. The IsNullOrEmpty guard yields null / 0
    /// (NOT "matches everywhere", NOT a crash). The empty case is special: an
    /// empty pattern is NOT treated as matching every position here — the unit
    /// returns the degenerate empty result (§3.3, §6.1).
    /// </summary>
    [TestCase("", "ACGT", 0, TestName = "EmptySequence_ZeroBudget")]
    [TestCase("", "ACGT", 5, TestName = "EmptySequence_LargeBudget")]
    [TestCase("ACGTACGT", "", 0, TestName = "EmptyPattern_ZeroBudget")]
    [TestCase("ACGTACGT", "", 5, TestName = "EmptyPattern_LargeBudget")]
    [TestCase("", "", 0, TestName = "BothEmpty")]
    [TestCase("", "", 99, TestName = "BothEmpty_LargeBudget")]
    public void EmptyInputs_BestMatchNull_CountZero(string sequence, string pattern, int maxMismatches)
    {
        ApproximateMatcher.FindBestMatch(sequence, pattern).Should().BeNull(
            "FindBestMatch returns null when either input is empty (§3.2)");

        ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches).Should().Be(0,
            "Count_d is 0 when either input is empty (§6.1)");
    }

    /// <summary>
    /// BE — empty sequence on FindFrequentKmersWithMismatches must yield-break on
    /// the IsNullOrEmpty guard BEFORE `counts.Values.Max()`. If the guard were
    /// missing, Max() over the empty tally would throw InvalidOperationException
    /// ("Sequence contains no elements") — the degenerate-empty crash this pins.
    /// </summary>
    [TestCase(4, 1, TestName = "EmptySeq_K4_D1")]
    [TestCase(1, 0, TestName = "EmptySeq_K1_D0")]
    [TestCase(12, 3, TestName = "EmptySeq_K12_D3")]
    public void FrequentKmers_EmptySequence_IsEmptyNoCrash(int k, int d)
    {
        Func<List<(string Kmer, int Count)>> act = () =>
            ApproximateMatcher.FindFrequentKmersWithMismatches("", k, d).ToList();

        act.Should().NotThrow("the empty-sequence guard must fire before Max() over an empty tally");
        act().Should().BeEmpty("there are no k-mers in an empty sequence");
    }

    /// <summary>
    /// BE/MC — invalid k and d on FindFrequentKmersWithMismatches throw the
    /// documented ArgumentOutOfRangeException (k ≤ 0; d &lt; 0), never an
    /// IndexOutOfRange / negative-length crash. Because the method is a lazy
    /// iterator the guard fires on enumeration, so the action materializes the
    /// result. The k ≤ 0 guard takes precedence over the empty-sequence
    /// yield-break only on a non-empty sequence (§3.1, §3.3).
    /// </summary>
    [TestCase(0, 1, TestName = "FrequentKmers_KZero_Throws")]
    [TestCase(-1, 1, TestName = "FrequentKmers_KNegative_Throws")]
    [TestCase(int.MinValue, 0, TestName = "FrequentKmers_KMinInt_Throws")]
    [TestCase(4, -1, TestName = "FrequentKmers_DNegative_Throws")]
    [TestCase(4, int.MinValue, TestName = "FrequentKmers_DMinInt_Throws")]
    public void FrequentKmers_InvalidKorD_ThrowsArgumentOutOfRange(int k, int d)
    {
        Func<List<(string Kmer, int Count)>> act = () =>
            ApproximateMatcher.FindFrequentKmersWithMismatches("ACGTACGTACGT", k, d).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
            "k ≤ 0 and d < 0 are the documented validation gates (§3.1)");
    }

    /// <summary>
    /// BE — negative maxMismatches on CountApproximateOccurrences throws the
    /// documented ArgumentOutOfRangeException (it delegates to FindWithMismatches,
    /// whose guard fires on enumeration), never a negative-length crash (§3.3).
    /// </summary>
    [TestCase(-1, TestName = "Count_NegativeBudget_MinusOne_Throws")]
    [TestCase(int.MinValue, TestName = "Count_NegativeBudget_MinInt_Throws")]
    public void Count_NegativeMaxMismatches_ThrowsArgumentOutOfRange(int maxMismatches)
    {
        Action act = () => ApproximateMatcher.CountApproximateOccurrences("ACGTACGT", "ACGT", maxMismatches);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "maxMismatches < 0 is the documented validation gate");
    }

    #endregion

    #region BE — exact present: an exact occurrence is found even when budget > 0

    /// <summary>
    /// BE — exact present: a pattern occurring EXACTLY (Hamming distance 0) is a
    /// special case of an approximate occurrence and MUST be found at the correct
    /// 0-based position even when the budget is strictly &gt; 0. FindBestMatch must
    /// short-circuit to the leftmost exact window with Distance 0 / IsExact, and
    /// CountApproximateOccurrences with a budget &gt; 0 must include it (INV-01,
    /// INV-04, §6.1). An exact match must never be MISSED because the budget is
    /// generous.
    /// </summary>
    // sequence, pattern, expectedFirstExactPosition
    [TestCase("ACGTACGT", "CGTA", 1, TestName = "ExactInterior")]
    [TestCase("GGGGACGTGGGG", "ACGT", 4, TestName = "ExactAfterPrefix")]
    [TestCase("ACGTGGGG", "ACGT", 0, TestName = "ExactAtStart")]
    [TestCase("GGGGACGT", "ACGT", 4, TestName = "ExactAtEnd")]
    [TestCase("ACGT", "ACGT", 0, TestName = "ExactWholeText")]
    [TestCase("AAAAAAAA", "AAAA", 0, TestName = "ExactLeftmostAmongMany")]
    public void ExactOccurrence_FoundEvenWithGenerousBudget(
        string sequence, string pattern, int expectedFirstExactPosition)
    {
        // FindBestMatch: exact occurrence ⇒ best distance 0 at the leftmost exact position.
        var best = ApproximateMatcher.FindBestMatch(sequence, pattern);
        best.Should().NotBeNull("an exact occurrence exists");
        var r = best!.Value;
        AssertWellFormedBestMatch(sequence, pattern, r);
        r.Distance.Should().Be(0, "an exact occurrence has Hamming distance 0 (INV-04)");
        r.IsExact.Should().BeTrue("distance 0 means IsExact");
        r.Position.Should().Be(expectedFirstExactPosition,
            "FindBestMatch returns the LEFTMOST exact (distance-0) window (INV-05)");

        // CountApproximateOccurrences with a generous budget still counts the exact occurrence.
        int exactCount = ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, 0);
        exactCount.Should().BeGreaterThanOrEqualTo(1, "the exact occurrence is counted at d=0");

        for (int budget = 1; budget <= pattern.Length; budget++)
        {
            int approxCount = ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, budget);
            approxCount.Should().BeGreaterThanOrEqualTo(exactCount,
                $"Count_d is monotone non-decreasing in d (INV-01); budget {budget} must include all exact occurrences");
        }
    }

    /// <summary>
    /// BE — exact present, positive sanity over the THREE surfaces at once on a
    /// known case. With pattern "GATG" in the BA1I sample text, the exact
    /// occurrences are found by FindBestMatch (distance 0 at the leftmost), counted
    /// by Count_0, and "GATG" is among the frequent 4-mers with 1 mismatch
    /// (BA1I sample: {GATG, ATGC, ATGT}, each Count_1 = 5). Pins that the exact
    /// case and the frequency case agree on the same text.
    /// </summary>
    [Test]
    public void ExactPresent_AcrossAllThreeSurfaces_BA1ISample()
    {
        const string text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";

        // "GATG" occurs exactly once (at position 16, 0-based) — find it at d=0.
        var best = ApproximateMatcher.FindBestMatch(text, "GATG");
        best.Should().NotBeNull();
        AssertWellFormedBestMatch(text, "GATG", best!.Value);
        best!.Value.IsExact.Should().BeTrue("GATG occurs exactly in the BA1I text");
        best.Value.Position.Should().Be(16, "the single exact GATG window starts at 0-based position 16");

        ApproximateMatcher.CountApproximateOccurrences(text, "GATG", 0).Should().Be(1,
            "GATG occurs exactly once at d=0");

        // BA1I sample: most frequent 4-mers with 1 mismatch are {GATG, ATGC, ATGT}, each count 5.
        var freq = ApproximateMatcher.FindFrequentKmersWithMismatches(text, k: 4, d: 1).ToList();
        var maxCount = freq.Select(f => f.Count).Distinct().Single();
        maxCount.Should().Be(5, "BA1I sample: the maximum Count_1 is 5");
        freq.Select(f => f.Kmer).Should().BeEquivalentTo(new[] { "GATG", "ATGC", "ATGT" },
            "BA1I sample returns all three tied frequent 4-mers (compared as a set, §4.2)");
    }

    #endregion

    #region Positive sanity — known approximate (1-mismatch) occurrence + BA1H/Count_1

    /// <summary>
    /// Positive sanity — a pattern with a known ONE-mismatch occurrence (and no
    /// exact one) is the best match at the correct position with distance 1. Text
    /// "TTTACGTTT" vs pattern "ACGA": the window "ACGT" at position 3 differs only
    /// at the last base (T vs A), so the minimum Hamming distance is 1 there. With
    /// budget 0 there is NO match (Count_0 = 0); with budget ≥ 1 it is counted.
    /// </summary>
    [Test]
    public void ApproximateOccurrence_OneMismatch_FoundAtCorrectPosition()
    {
        const string text = "TTTACGTTT";
        const string pattern = "ACGA";

        var best = ApproximateMatcher.FindBestMatch(text, pattern);
        best.Should().NotBeNull();
        var r = best!.Value;
        AssertWellFormedBestMatch(text, pattern, r);
        r.Distance.Should().Be(1, "the closest window 'ACGT' at position 3 differs from 'ACGA' at one base");
        r.IsExact.Should().BeFalse("there is no exact occurrence of ACGA");
        r.Position.Should().Be(3, "the minimum-distance window starts at 0-based position 3");
        r.MismatchPositions.Should().Equal(new[] { 3 }, "the single mismatch is at index 3 of the window");

        ApproximateMatcher.CountApproximateOccurrences(text, pattern, 0).Should().Be(0,
            "no exact occurrence ⇒ Count_0 = 0");
        ApproximateMatcher.CountApproximateOccurrences(text, pattern, 1).Should().Be(1,
            "the single 1-mismatch window is counted once at d=1");
    }

    /// <summary>
    /// Positive sanity — the BA1H worked example: Pattern=ATTCTGGA in a 99-nt text
    /// with d=3 has Count_3 = 5 (occurrence positions 6, 7, 26, 27, 78). Pins
    /// CountApproximateOccurrences against a sourced primary-source value, and
    /// confirms the leftmost best match is the closest of those windows.
    /// </summary>
    [Test]
    public void BA1H_CountApproximateOccurrences_MatchesSourcedValue()
    {
        const string text =
            "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC";
        const string pattern = "ATTCTGGA";

        ApproximateMatcher.CountApproximateOccurrences(text, pattern, 3).Should().Be(5,
            "BA1H sample: Count_3(text, ATTCTGGA) = 5 (positions 6 7 26 27 78)");

        // The best (minimum-distance) window must itself be a valid <= 3-mismatch window.
        var best = ApproximateMatcher.FindBestMatch(text, pattern);
        best.Should().NotBeNull();
        AssertWellFormedBestMatch(text, pattern, best!.Value);
        best!.Value.Distance.Should().BeLessThanOrEqualTo(3,
            "the minimum-distance window is one of the d=3 approximate occurrences");
    }

    #endregion

    #region BE/OVF — randomized boundary scan: well-formed best match, no hang, count agreement

    /// <summary>
    /// BE/OVF — randomized fuzz over the three boundaries: random pattern length
    /// relative to text (shorter, equal, longer), random budgets including 0 and
    /// huge. For every input EITHER FindBestMatch is null (when empty or pattern >
    /// text) OR it is a well-formed result whose recomputed distance is the true
    /// minimum over all windows (INV-04). CountApproximateOccurrences must agree
    /// with an independent count of windows within budget. Bounded sizes +
    /// CancelAfter guard against any hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void RandomizedBoundaries_BestMatchWellFormed_CountAgrees()
    {
        for (int trial = 0; trial < 400; trial++)
        {
            int textLen = Rng.Next(0, 24);
            string text = RandomDna(textLen);

            // Pattern length intentionally spans shorter / equal / longer than the text.
            int patLen = Rng.Next(0, textLen + 4);
            string pattern = RandomDna(patLen);

            int budget = Rng.Next(0, 6) == 0 ? int.MaxValue : Rng.Next(0, patLen + 2);

            var best = ApproximateMatcher.FindBestMatch(text, pattern);

            bool noWindow = text.Length == 0 || pattern.Length == 0 || pattern.Length > text.Length;
            if (noWindow)
            {
                best.Should().BeNull(
                    $"no equal-length window for text='{text}' pattern='{pattern}'");
                ApproximateMatcher.CountApproximateOccurrences(text, pattern, budget).Should().Be(0,
                    "Count_d is 0 when no equal-length window exists");
                continue;
            }

            best.Should().NotBeNull($"a window exists for text='{text}' pattern='{pattern}'");
            AssertWellFormedBestMatch(text, pattern, best!.Value);

            // Independent minimum-Hamming over all windows: the reported distance must equal it,
            // and the reported position must be the LEFTMOST achieving that minimum (INV-04/INV-05).
            int n = text.Length, m = pattern.Length;
            string up = text.ToUpperInvariant();
            string pu = pattern.ToUpperInvariant();
            int trueMin = int.MaxValue;
            int leftmostMinPos = -1;
            int within = 0;
            for (int i = 0; i <= n - m; i++)
            {
                int dist = HammingOracle(pu, up.Substring(i, m));
                if (dist < trueMin) { trueMin = dist; leftmostMinPos = i; }
                if (dist <= budget) within++;
            }

            best!.Value.Distance.Should().Be(trueMin,
                "FindBestMatch distance must be the true minimum Hamming distance over all windows");
            best.Value.Position.Should().Be(leftmostMinPos,
                "FindBestMatch must report the leftmost window achieving the minimum (INV-05)");

            ApproximateMatcher.CountApproximateOccurrences(text, pattern, budget).Should().Be(within,
                "Count_d must equal the independent count of windows within the budget");
        }
    }

    /// <summary>
    /// OVF — large inputs must not hang: a long random text scanned by all three
    /// surfaces completes well under the cancel budget, and Count_d is monotone
    /// non-decreasing in d (INV-01). Guards against an accidental quadratic blow-up
    /// or infinite loop on sizeable boundary input.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void LargeInput_DoesNotHang_AndCountIsMonotone()
    {
        string text = RandomDna(4000);
        string pattern = RandomDna(8);

        var best = ApproximateMatcher.FindBestMatch(text, pattern);
        best.Should().NotBeNull("an 8-mer pattern fits in a 4000-nt text");
        AssertWellFormedBestMatch(text, pattern, best!.Value);

        int prev = -1;
        for (int d = 0; d <= pattern.Length; d++)
        {
            int c = ApproximateMatcher.CountApproximateOccurrences(text, pattern, d);
            c.Should().BeGreaterThanOrEqualTo(prev, "Count_d is monotone non-decreasing in d (INV-01)");
            prev = c;
        }

        // FindFrequentKmersWithMismatches on a sizeable text with small k,d must also complete.
        Func<List<(string Kmer, int Count)>> freq = () =>
            ApproximateMatcher.FindFrequentKmersWithMismatches(text, k: 6, d: 1).ToList();
        freq.Should().NotThrow();
        freq().Should().NotBeEmpty("a non-empty text has at least one frequent k-mer");
    }

    #endregion
}
