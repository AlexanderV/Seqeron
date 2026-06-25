using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the K-mer area — UNIQUE K-MERS (KMER-UNIQUE-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no DivideByZero (an empty count table has no division), no IndexOutOfRange
/// from an internal Substring, and no *unhandled* runtime exception. Every input
/// must resolve to EITHER a well-defined, theory-correct unique-k-mer set, OR a
/// *documented, intentional* validation exception
/// (<see cref="ArgumentOutOfRangeException"/> for k ≤ 0 on a non-empty sequence).
/// A raw runtime exception, a hang, a duplicated entry, a wrongly-included repeated
/// k-mer, or a missing singleton on a boundary input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-UNIQUE-001 — unique (frequency-1) k-mers
/// Checklist: docs/checklists/03_FUZZING.md, row 162.
/// Algorithm doc: docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: ALL-IDENTICAL (a homopolymer), ALL-DISTINCT (a
///          maximally-diverse sequence where every window differs), and the EMPTY
///          sequence. We also pin the off-by-one neighbours k = L (a single,
///          trivially-unique window), k = L + 1 (no window), and the k ≤ 0 floor.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The unique-k-mer contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The API entry under test is
///   KmerAnalyzer.FindUniqueKmers(string sequence, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines
///    253–257), delegating counting to CountKmers (lines 20–42).
///
/// "Unique" here means OCCURRENCE COUNT EXACTLY 1 (singletons) — NOT the distinct
/// set. The doc is explicit (Unique_And_MinCount_Kmers.md §2.2):
///   • Unique k-mers:  { P : c(P) = 1 } — "Unique k-mers are those that appear
///                     only once" [BioInfoLogics 2018].
///   • distinct  = each different string once;  unique ⊆ distinct (INV-02).
/// For a sequence of length L there are T = L − k + 1 overlapping, step-1 windows
/// (§2.1). Filter predicate is c(P) == 1 (§4.1, INV-01).
///
/// Documented invariants (Unique_And_MinCount_Kmers.md §2.4):
///   INV-01 every returned k-mer has c(P) = 1.
///   INV-02 unique ⊆ distinct ⇒ |unique| ≤ |distinct| ≤ L − k + 1.
///
/// Documented edge cases (§3.3, §6.1):
///   • null / empty sequence → empty result (no k-mers exist; the null/empty guard
///     in CountKmers wins BEFORE k is validated).
///   • k > L                 → empty result (L − k + 1 ≤ 0).
///   • k ≤ 0 (non-empty seq) → ArgumentOutOfRangeException(nameof(k)).
///   • Homopolymer (AAAAA,k=3) → NO unique k-mers ("AAA" has c = 3 > 1).
///
/// The three checklist BE targets map to these documented behaviours:
///   • all-identical (homopolymer "AAA…", L > k) → the single distinct k-mer occurs
///       L − k + 1 > 1 times ⇒ NOT a singleton ⇒ EMPTY unique set. The lone
///       exception is L − k + 1 == 1 (k == L): one window, trivially count-1.
///   • all-distinct (every window different) → EVERY window is a singleton ⇒ the
///       unique set equals the full distinct set, |unique| = L − k + 1.
///   • empty sequence → EMPTY unique set, no crash / DivideByZero.
/// A positive-sanity test pins the documented worked example
/// (ATCGATCAC, k=3, Unique_And_MinCount_Kmers.md §2.2/§7.1): 7 total, 6 distinct,
/// 5 unique = {TCG, CGA, GAT, TCA, CAC}; ATC (c = 2) is distinct but excluded.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerUniqueFuzzTests
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

    /// <summary>
    /// Asserts the result of <see cref="KmerAnalyzer.FindUniqueKmers"/> is a
    /// WELL-FORMED unique set for the given sequence and k, judged against an
    /// INDEPENDENT reference count over the overlapping step-1 windows:
    ///   • every returned k-mer has length exactly k (INV — well-formed members);
    ///   • the OUTPUT has no duplicates (it is a set of distinct strings);
    ///   • every returned k-mer occurs EXACTLY ONCE in the sequence (INV-01:
    ///     c(P) = 1) — re-derived here, not trusted from the SUT;
    ///   • the returned set equals { P : c(P) = 1 } exactly (no missing singleton,
    ///     no included repeat);
    ///   • |unique| ≤ |distinct| ≤ L − k + 1 (INV-02).
    /// This is the load-bearing structural oracle reused across the fuzz cases
    /// (Unique_And_MinCount_Kmers.md §2.4 INV-01/INV-02).
    /// </summary>
    private static void AssertWellFormed(IEnumerable<string> result, string sequence, int k)
    {
        var unique = result.ToList();

        // Output is a SET of distinct strings — no duplicate entries.
        unique.Should().OnlyHaveUniqueItems("the unique-k-mer output is a set of distinct strings");

        // Every member is a length-k substring.
        unique.Should().OnlyContain(km => km.Length == k, "every k-mer has length k");

        // Independent reference: count overlapping step-1 windows of the UPPER-CASED
        // sequence (CountKmers upper-cases internally) and derive the singleton set.
        var counts = new Dictionary<string, int>();
        if (!string.IsNullOrEmpty(sequence) && k > 0 && k <= sequence.Length)
        {
            string seq = sequence.ToUpperInvariant();
            for (int i = 0; i <= seq.Length - k; i++)
            {
                string km = seq.Substring(i, k);
                counts[km] = counts.TryGetValue(km, out int c) ? c + 1 : 1;
            }
        }

        int totalWindows = (string.IsNullOrEmpty(sequence) || k <= 0 || k > sequence.Length)
            ? 0
            : sequence.Length - k + 1;
        var expectedUnique = counts.Where(kvp => kvp.Value == 1).Select(kvp => kvp.Key).ToHashSet();

        // INV-01: the set is EXACTLY the singletons — no missing one, no repeat included.
        unique.Should().BeEquivalentTo(expectedUnique,
            "FindUniqueKmers returns exactly { P : c(P) = 1 } (INV-01)");

        // INV-01 re-derived per member: each returned k-mer truly occurs once.
        foreach (string km in unique)
            counts[km].Should().Be(1, $"'{km}' is reported unique ⇒ it must occur exactly once");

        // INV-02: unique ⊆ distinct ⇒ |unique| ≤ |distinct| ≤ L − k + 1.
        unique.Count.Should().BeLessThanOrEqualTo(counts.Count, "INV-02: |unique| ≤ |distinct|");
        unique.Count.Should().BeLessThanOrEqualTo(Math.Max(totalWindows, 0),
            "INV-02: |unique| ≤ L − k + 1");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-UNIQUE-001 — unique (frequency-1) k-mers : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-UNIQUE-001 — unique k-mers

    #region BE — Boundary: empty / null sequence

    /// <summary>
    /// BE / headline target: the EMPTY sequence yields an EMPTY unique set — no
    /// k-mers exist, no crash, no DivideByZero (the count table is empty so nothing
    /// is divided). Unique_And_MinCount_Kmers.md §6.1.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_EmptySequence_ReturnsEmpty()
    {
        var unique = KmerAnalyzer.FindUniqueKmers("", 3).ToList();

        unique.Should().BeEmpty("empty sequence has no k-mers (§6.1)");
        AssertWellFormed(unique, "", 3);
    }

    /// <summary>
    /// BE: a null sequence is short-circuited to the empty result before k is
    /// validated (CountKmers null/empty guard). No NullReferenceException.
    /// Unique_And_MinCount_Kmers.md §3.3.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_NullSequence_ReturnsEmpty()
    {
        var unique = KmerAnalyzer.FindUniqueKmers(null!, 4).ToList();

        unique.Should().BeEmpty("null sequence has no k-mers");
    }

    /// <summary>
    /// BE: the empty/null guard wins even when k is itself degenerate (k ≤ 0):
    /// CountKmers returns the empty map for null/empty BEFORE validating k, so
    /// FindUniqueKmers returns empty rather than throwing.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_EmptyOrNull_WithNonPositiveK_StillEmpty()
    {
        foreach (int k in new[] { 0, -1, int.MinValue })
        {
            KmerAnalyzer.FindUniqueKmers("", k).Should().BeEmpty(
                $"empty input short-circuits before k={k} is validated");
            KmerAnalyzer.FindUniqueKmers(null!, k).Should().BeEmpty(
                $"null input short-circuits before k={k} is validated");
        }
    }

    #endregion

    #region BE — Boundary: all-identical (homopolymer)

    /// <summary>
    /// BE / headline target: an ALL-IDENTICAL sequence ("AAAA…", L > k) has exactly
    /// ONE distinct k-mer "AA…A" that occurs L − k + 1 > 1 times ⇒ it is NOT a
    /// singleton ⇒ the unique set is EMPTY. Off-by-one trap: this only holds while
    /// L − k + 1 > 1; the k == L case (one window) is pinned separately below.
    /// Unique_And_MinCount_Kmers.md §6.1 (homopolymer → no unique k-mers).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_Homopolymer_MoreThanOneWindow_ReturnsEmpty()
    {
        foreach (int len in new[] { 2, 3, 5, 10, 64 })
        {
            foreach (int k in new[] { 1, 2, len - 1 })
            {
                if (k <= 0 || k >= len) continue; // ensure L − k + 1 ≥ 2 (the repeated case)
                var unique = KmerAnalyzer.FindUniqueKmers(new string('A', len), k).ToList();

                unique.Should().BeEmpty(
                    $"homopolymer of length {len}, k={k}: the single k-mer occurs {len - k + 1}× (> 1) ⇒ not unique");
                AssertWellFormed(unique, new string('A', len), k);
            }
        }
    }

    /// <summary>
    /// BE off-by-one: a homopolymer with k == L has EXACTLY ONE window, so the lone
    /// k-mer is trivially count-1 ⇒ the unique set is { "AA…A" } (size 1), NOT empty.
    /// This is the boundary where the "homopolymer ⇒ no unique" rule flips, and a
    /// naive "homopolymers never have unique k-mers" shortcut would be wrong here.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_Homopolymer_KEqualsLength_SingleTriviallyUnique()
    {
        foreach (int len in new[] { 1, 2, 5, 10 })
        {
            string seq = new string('A', len); // k = len ⇒ exactly one window
            var unique = KmerAnalyzer.FindUniqueKmers(seq, len).ToList();

            unique.Should().ContainSingle().Which.Should().Be(seq,
                $"k == L == {len}: one window ⇒ that single k-mer is trivially unique");
            AssertWellFormed(unique, seq, len);
        }
    }

    #endregion

    #region BE — Boundary: all-distinct (every window a singleton)

    /// <summary>
    /// BE / headline target: an ALL-DISTINCT sequence — every overlapping window is
    /// a different string ⇒ EVERY k-mer is a singleton ⇒ unique = the full distinct
    /// set, |unique| = L − k + 1. We use a de-Bruijn-flavoured construction where no
    /// k-mer repeats, then assert the count and that the unique set covers all windows.
    /// Unique_And_MinCount_Kmers.md §2.4 (INV-02 upper bound attained).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_AllDistinct_EveryWindowIsUnique()
    {
        // "ACGTACGTACGT…" is NOT all-distinct for small k; instead build a sequence
        // whose every length-k window differs. The lexicographic ramp
        // A,C,G,T then 2-char, etc. is awkward; simplest robust generator:
        // take a long de-Bruijn-like string and pick (seq, k) pairs that we VERIFY
        // are all-distinct via the reference count.
        foreach ((string seq, int k) in AllDistinctCases())
        {
            int windows = seq.Length - k + 1;
            var unique = KmerAnalyzer.FindUniqueKmers(seq, k).ToList();

            unique.Should().HaveCount(windows,
                $"every one of the {windows} windows of '{seq}' (k={k}) is distinct ⇒ all are unique");
            AssertWellFormed(unique, seq, k);

            // The unique set must cover EVERY window position.
            var allWindows = Enumerable.Range(0, windows)
                .Select(i => seq.ToUpperInvariant().Substring(i, k))
                .ToHashSet();
            unique.Should().BeEquivalentTo(allWindows,
                "all-distinct ⇒ unique set = the set of every window");
        }
    }

    /// <summary>Yields (sequence, k) pairs verified all-distinct (no repeated window).</summary>
    private static IEnumerable<(string seq, int k)> AllDistinctCases()
    {
        // de Bruijn B(4,2) over ACGT yields every 2-mer exactly once when read as
        // a length-17 cyclic-but-linearised string; we use a known all-distinct
        // hand string plus a couple of generated ones, each filtered for distinctness.
        string[] seqs = { "ACGTAGT", "ACGTCAGT", "AACGTCGA", "ACAGCTGT", "GTAGAGCTGT" };
        foreach (string s in seqs)
            foreach (int k in new[] { 3, 4 })
                if (k <= s.Length && AllWindowsDistinct(s, k))
                    yield return (s, k);

        // Generated DNA at large k relative to L is overwhelmingly all-distinct;
        // filter to be certain.
        for (int seed = 0; seed < 40; seed++)
        {
            string s = RandomDna(12, seed);
            int k = 6;
            if (AllWindowsDistinct(s, k))
                yield return (s, k);
        }
    }

    private static bool AllWindowsDistinct(string s, int k)
    {
        string seq = s.ToUpperInvariant();
        var seen = new HashSet<string>();
        for (int i = 0; i <= seq.Length - k; i++)
            if (!seen.Add(seq.Substring(i, k)))
                return false;
        return true;
    }

    #endregion

    #region BE — Boundary: k > L, k = L (window-count edge)

    /// <summary>
    /// BE: k strictly greater than the sequence length yields an EMPTY unique set —
    /// L − k + 1 ≤ 0 means no windows. No negative-length Substring, no crash.
    /// Unique_And_MinCount_Kmers.md §6.1 (k > L → empty).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_KGreaterThanLength_ReturnsEmpty()
    {
        const string seq = "ACGT"; // L = 4
        foreach (int k in new[] { 5, 6, 100, int.MaxValue })
        {
            var unique = KmerAnalyzer.FindUniqueKmers(seq, k).ToList();
            unique.Should().BeEmpty($"k={k} > L=4 ⇒ no windows ⇒ no unique k-mers");
            AssertWellFormed(unique, seq, k);
        }
    }

    /// <summary>
    /// BE off-by-one: k = L gives EXACTLY one window ⇒ the whole sequence is a single
    /// trivially-unique k-mer; k = L + 1 gives no window ⇒ empty. Pins the
    /// one-window / no-window boundary precisely on a non-homopolymer.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindUniqueKmers_KEqualsLengthBoundary_OneWindowVsNone()
    {
        const string seq = "ACGTA"; // L = 5

        var atLen = KmerAnalyzer.FindUniqueKmers(seq, seq.Length).ToList(); // one window
        atLen.Should().ContainSingle().Which.Should().Be("ACGTA",
            "k == L ⇒ one window, trivially unique");
        AssertWellFormed(atLen, seq, seq.Length);

        var pastLen = KmerAnalyzer.FindUniqueKmers(seq, seq.Length + 1).ToList();
        pastLen.Should().BeEmpty("k = L + 1 ⇒ no window");
    }

    #endregion

    #region BE — Boundary: k ≤ 0 floor (documented validation exception)

    /// <summary>
    /// BE: k ≤ 0 with NON-EMPTY input throws ArgumentOutOfRangeException(nameof(k)),
    /// surfaced unchanged from CountKmers — a 0-length k-mer is meaningless. This is
    /// the documented, intentional validation exception (Unique_And_MinCount_Kmers.md
    /// §3.3). Materialise the enumerable (ToList) to force the deferred call.
    /// </summary>
    [Test]
    public void FindUniqueKmers_NonPositiveK_NonEmptySequence_Throws()
    {
        foreach (int k in new[] { 0, -1, -100, int.MinValue })
        {
            Action act = () => KmerAnalyzer.FindUniqueKmers("ACGTACGT", k).ToList();
            act.Should().Throw<ArgumentOutOfRangeException>()
               .WithParameterName("k", $"k={k} ≤ 0 is rejected by the documented contract");
        }
    }

    #endregion

    #region BE/RB — Random fuzz: structural invariants hold for arbitrary inputs

    /// <summary>
    /// Random-input sweep over assorted lengths and k values (including k &gt; len
    /// and k = len boundaries): EVERY result must satisfy the well-formed oracle —
    /// a duplicate-free set of length-k strings, each occurring exactly once, equal
    /// to the independently-counted singleton set, with |unique| ≤ L − k + 1. No
    /// undisciplined failure on any boundary k. (k ≤ 0 on a non-empty sequence is
    /// the throwing contract, pinned separately, so it is excluded here.)
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindUniqueKmers_RandomInputs_AlwaysWellFormed()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 60);
            string seq = RandomDna(len, seed);
            int k = rng.Next(1, len + 5); // valid k, k = L, and k > L; never k ≤ 0

            var unique = KmerAnalyzer.FindUniqueKmers(seq, k).ToList();
            AssertWellFormed(unique, seq, k);
        }
    }

    #endregion

    #region POSITIVE sanity — documented worked example and a mixed spectrum

    /// <summary>
    /// POSITIVE: the documented worked example (ATCGATCAC, k=3,
    /// Unique_And_MinCount_Kmers.md §2.2/§7.1). 7 windows: ATC, TCG, CGA, GAT, ATC,
    /// TCA, CAC. ATC occurs twice (distinct but NOT unique); the other 5 are
    /// singletons ⇒ unique = {TCG, CGA, GAT, TCA, CAC}. Hand-computed — pins the
    /// business contract (singleton, not distinct), not just green.
    /// </summary>
    [Test]
    public void FindUniqueKmers_WorkedExample_MatchesDocumentedSet()
    {
        var unique = KmerAnalyzer.FindUniqueKmers("ATCGATCAC", 3).ToHashSet();

        unique.Should().BeEquivalentTo(new[] { "TCG", "CGA", "GAT", "TCA", "CAC" },
            "documented §2.2: 5 unique; ATC (c=2) is distinct but excluded");
        unique.Should().NotContain("ATC", "ATC occurs twice ⇒ distinct but NOT unique");
        unique.Should().HaveCount(5);
        AssertWellFormed(unique, "ATCGATCAC", 3);
    }

    /// <summary>
    /// POSITIVE: a hand-crafted MIXED spectrum that exercises the singleton/repeat
    /// boundary directly. "ATATATCG" (L=8), k=2 ⇒ windows AT,TA,AT,TA,AT,TC,CG.
    /// Counts: AT=3, TA=2, TC=1, CG=1 ⇒ unique = {TC, CG} (the count-1 set), NOT the
    /// distinct set {AT,TA,TC,CG}. Confirms repeats (AT, TA) are EXCLUDED.
    /// </summary>
    [Test]
    public void FindUniqueKmers_MixedSpectrum_ReturnsOnlySingletons()
    {
        var unique = KmerAnalyzer.FindUniqueKmers("ATATATCG", 2).ToHashSet();

        unique.Should().BeEquivalentTo(new[] { "TC", "CG" },
            "only TC and CG occur exactly once; AT (×3) and TA (×2) are repeats");
        unique.Should().NotContain("AT").And.NotContain("TA");
        AssertWellFormed(unique, "ATATATCG", 2);
    }

    /// <summary>
    /// POSITIVE: case-insensitivity — lower-case input yields the same unique set as
    /// its upper-cased form (input is upper-cased internally), and the returned
    /// k-mers are upper-case. Unique_And_MinCount_Kmers.md §3.3.
    /// </summary>
    [Test]
    public void FindUniqueKmers_LowerCase_IdenticalToUpperCase()
    {
        var lower = KmerAnalyzer.FindUniqueKmers("atcgatcac", 3).ToHashSet();
        var upper = KmerAnalyzer.FindUniqueKmers("ATCGATCAC", 3).ToHashSet();

        lower.Should().BeEquivalentTo(upper, "input is upper-cased internally");
        lower.Should().OnlyContain(km => km == km.ToUpperInvariant(), "returned k-mers are upper-cased");
    }

    #endregion

    #endregion
}
