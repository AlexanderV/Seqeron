using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the K-mer area — K-MER STATISTICS (KMER-STATS-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (DivideByZero on total/distinct, negative TotalKmers when
/// k > L, NaN/Infinity entropy from log(0) or an empty distribution,
/// IndexOutOfRange from an internal Substring). Every input must resolve to
/// EITHER a well-defined, theory-correct <see cref="KmerStatistics"/>, OR a
/// *documented, intentional* validation exception
/// (<see cref="ArgumentOutOfRangeException"/> for k ≤ 0). A raw runtime exception,
/// a hang, a NaN entropy, or a negative/wrong statistic on a boundary input is a
/// bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-STATS-001 — k-mer composition statistics
/// Checklist: docs/checklists/03_FUZZING.md, row 161.
/// Algorithm doc: docs/algorithms/K-mer/K-mer_Statistics.md.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the EMPTY sequence, k &gt; len, and a HOMOPOLYMER
///          ("AAAA…"). We also pin the off-by-one neighbours k = L (one window)
///          and k = L + 1 (no window), plus the k ≤ 0 floor.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The statistics contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The API entry under test is
///   KmerAnalyzer.AnalyzeKmers(string sequence, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines
///    528–560), returning the KmerStatistics record (lines 572–578).
/// It builds the k-mer count table via CountKmers, then aggregates. For the count
/// table mult(α) of each distinct k-mer α, with T = Σ mult(α) = L − k + 1 windows
/// (K-mer_Statistics.md §2.2):
///   • TotalKmers  = T = L − k + 1            (INV-01, INV-02)
///   • UniqueKmers = D = #distinct k-mers      (INV-03)
///   • MaxCount    = max_α mult(α)
///   • MinCount    = min_α mult(α)
///   • AverageCount= T / D, Math.Round(·, 2)   (INV-04)
///   • Entropy     = −Σ p(α) log₂ p(α), p(α)=mult(α)/T, in bits (INV-05)
/// Bounds (K-mer_Statistics.md §2.4):
///   INV-04 MinCount ≤ AverageCount ≤ MaxCount
///   INV-05 0 ≤ Entropy ≤ log₂(D); Entropy = 0 iff D = 1 (one distinct k-mer);
///          Entropy = log₂(D) iff all multiplicities equal (uniform).
///   INV-06 empty seq or k > L ⇒ ALL fields 0 (L − k + 1 ≤ 0 ⇒ no k-mers).
///
/// The three checklist BE targets map to these documented behaviours:
///   • empty seq    → KmerStatistics(0,0,0,0,0,0). CountKmers short-circuits
///                    null/empty to the empty map BEFORE k validation, and
///                    AnalyzeKmers returns the all-zero record (no DivideByZero on
///                    the T/D mean, no NaN entropy). K-mer_Statistics.md §6.1.
///   • k &gt; len      → all-zero record. T = L − k + 1 ≤ 0 is GUARDED to 0 (never
///                    a negative TotalKmers); the empty count table triggers the
///                    all-zero early return. K-mer_Statistics.md §2.4 INV-06.
///   • homopolymer  → ONE distinct k-mer ("AAA…"): D = 1, TotalKmers = L − k + 1,
///                    Max = Min = Average = T, Entropy = 0 (one-component
///                    distribution). K-mer_Statistics.md §6.1, INV-05.
/// A positive-sanity test pins the documented worked example
/// (GTAGAGCTGT, K-mer_Statistics.md §7.1): k=1 ⇒ Total=10, Distinct=4, Max=4,
/// Min=1, Avg=2.5, Entropy ≈ 1.846439; k=3 ⇒ Total=Distinct=8, Max=Min=1,
/// Avg=1.0, Entropy = 3.0 (= log₂8, all windows distinct ⇒ Distinct = Total).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerStatsFuzzTests
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
    /// Asserts the result is a WELL-FORMED <see cref="KmerStatistics"/> for the
    /// given sequence and k: every field finite and non-negative, the documented
    /// total/distinct relationships, the AverageCount bracket
    /// (Min ≤ Avg ≤ Max, Avg ≈ T/D rounded), and the entropy bounds
    /// 0 ≤ Entropy ≤ log₂(D). This is the load-bearing structural oracle reused
    /// across the fuzz cases (K-mer_Statistics.md §2.4 INV-01…INV-06).
    /// </summary>
    private static void AssertWellFormed(KmerStatistics s, string sequence, int k)
    {
        int len = sequence?.Length ?? 0;
        int expectedTotal = len - k + 1;
        if (expectedTotal < 0) expectedTotal = 0; // INV-06: k > L is guarded to 0.

        // No field is NaN / Infinity — entropy and average must be real numbers.
        double.IsFinite(s.AverageCount).Should().BeTrue(
            "AverageCount must be a finite real number, never NaN/Infinity");
        double.IsFinite(s.Entropy).Should().BeTrue(
            "Entropy must be a finite real number — never NaN from log(0) or an empty distribution");

        // No field is negative — including TotalKmers when k > L (must be guarded to 0).
        s.TotalKmers.Should().BeGreaterThanOrEqualTo(0, "TotalKmers = max(L − k + 1, 0); never negative");
        s.UniqueKmers.Should().BeGreaterThanOrEqualTo(0);
        s.MaxCount.Should().BeGreaterThanOrEqualTo(0);
        s.MinCount.Should().BeGreaterThanOrEqualTo(0);
        s.AverageCount.Should().BeGreaterThanOrEqualTo(0);
        s.Entropy.Should().BeGreaterThanOrEqualTo(0, "INV-05: Shannon entropy is non-negative");

        // INV-01/INV-02: total = L − k + 1 = Σ mult(α).
        s.TotalKmers.Should().Be(expectedTotal, "INV-01: TotalKmers = max(L − k + 1, 0)");

        if (s.TotalKmers == 0)
        {
            // INV-06: all-zero record when no k-mers exist.
            s.UniqueKmers.Should().Be(0);
            s.MaxCount.Should().Be(0);
            s.MinCount.Should().Be(0);
            s.AverageCount.Should().Be(0);
            s.Entropy.Should().Be(0);
            return;
        }

        // INV-03: distinct count is between 1 and total (each window ≥ 1 distinct,
        // at most T distinct when all windows differ).
        s.UniqueKmers.Should().BeInRange(1, s.TotalKmers,
            "INV-03: 1 ≤ distinct ≤ total");

        // INV-04: MinCount ≤ AverageCount ≤ MaxCount and Average ≈ T/D (rounded 2dp).
        s.MinCount.Should().BeLessThanOrEqualTo(s.MaxCount, "min ≤ max");
        s.MaxCount.Should().BeLessThanOrEqualTo(s.TotalKmers, "a multiplicity cannot exceed the window count");
        s.MinCount.Should().BeGreaterThanOrEqualTo(1, "every distinct k-mer occurs ≥ once");
        double exactMean = (double)s.TotalKmers / s.UniqueKmers;
        // Reported mean is the exact ratio rounded to 2 decimals ⇒ |reported − exact| ≤ 0.005.
        s.AverageCount.Should().BeApproximately(exactMean, 0.0051,
            "INV-04: AverageCount = TotalKmers / UniqueKmers, rounded to 2 decimals");
        s.AverageCount.Should().BeInRange(s.MinCount - 0.0051, s.MaxCount + 0.0051,
            "INV-04: Min ≤ Average ≤ Max");

        // INV-05: 0 ≤ Entropy ≤ log₂(D).
        double maxEntropy = Math.Log2(s.UniqueKmers);
        s.Entropy.Should().BeLessThanOrEqualTo(maxEntropy + 1e-9,
            "INV-05: Entropy ≤ log₂(distinct)");

        // INV-05 boundary: D = 1 ⇒ Entropy = 0 (one-component distribution).
        if (s.UniqueKmers == 1)
            s.Entropy.Should().BeApproximately(0.0, 1e-12,
                "INV-05: a single distinct k-mer has zero entropy");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-STATS-001 — k-mer composition statistics : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-STATS-001 — k-mer statistics

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the EMPTY sequence yields the documented all-zero record — no
    /// DivideByZero on the T/D mean (D = 0), no NaN entropy. K-mer_Statistics.md
    /// §6.1, INV-06.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_EmptySequence_ReturnsAllZero()
    {
        var stats = KmerAnalyzer.AnalyzeKmers("", 3);

        stats.Should().Be(new KmerStatistics(0, 0, 0, 0, 0, 0),
            "empty sequence has no k-mers ⇒ all fields 0 (INV-06)");
        AssertWellFormed(stats, "", 3);
    }

    /// <summary>
    /// BE: null sequence is short-circuited to the all-zero record before k is
    /// validated (CountKmers null/empty guard). No NullReference, no DivideByZero.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_NullSequence_ReturnsAllZero()
    {
        var stats = KmerAnalyzer.AnalyzeKmers(null!, 4);

        stats.Should().Be(new KmerStatistics(0, 0, 0, 0, 0, 0),
            "null sequence has no k-mers ⇒ all fields 0");
    }

    /// <summary>
    /// BE: the empty-sequence guard wins even when k is itself degenerate
    /// (k ≤ 0): CountKmers returns the empty map for null/empty BEFORE validating
    /// k, so AnalyzeKmers returns the all-zero record rather than throwing.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_EmptySequence_WithNonPositiveK_StillAllZero()
    {
        foreach (int k in new[] { 0, -1, int.MinValue })
        {
            var stats = KmerAnalyzer.AnalyzeKmers("", k);
            stats.Should().Be(new KmerStatistics(0, 0, 0, 0, 0, 0),
                $"empty input short-circuits before k={k} is validated");
        }
    }

    #endregion

    #region BE — Boundary: k > len (no windows, guarded total)

    /// <summary>
    /// BE: k strictly greater than the sequence length yields the all-zero record.
    /// L − k + 1 ≤ 0 must be GUARDED to 0 — TotalKmers is never negative, and the
    /// empty count table gives an all-zero, NaN-free result. K-mer_Statistics.md
    /// §6.1 INV-06.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_KGreaterThanLength_ReturnsAllZeroNoNegativeTotal()
    {
        const string seq = "ACGT"; // L = 4

        foreach (int k in new[] { 5, 6, 100, int.MaxValue })
        {
            var stats = KmerAnalyzer.AnalyzeKmers(seq, k);

            stats.TotalKmers.Should().Be(0,
                $"k={k} > L=4 ⇒ no windows; L − k + 1 ≤ 0 guarded to 0, never negative");
            stats.Should().Be(new KmerStatistics(0, 0, 0, 0, 0, 0));
            AssertWellFormed(stats, seq, k);
        }
    }

    /// <summary>
    /// BE off-by-one: k = L + 1 gives no window (all-zero); k = L gives EXACTLY one
    /// window ⇒ Total = Distinct = 1, Max = Min = Avg = 1, Entropy = 0. Pins the
    /// no-window / one-window boundary precisely.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_KEqualsLengthBoundary_OneWindowVsNone()
    {
        const string seq = "ACGTA"; // L = 5

        var atLen = KmerAnalyzer.AnalyzeKmers(seq, seq.Length);       // k = 5 ⇒ 1 window
        atLen.TotalKmers.Should().Be(1);
        atLen.UniqueKmers.Should().Be(1);
        atLen.MaxCount.Should().Be(1);
        atLen.MinCount.Should().Be(1);
        atLen.AverageCount.Should().Be(1.0);
        atLen.Entropy.Should().BeApproximately(0.0, 1e-12, "a single window ⇒ one distinct k-mer ⇒ entropy 0");
        AssertWellFormed(atLen, seq, seq.Length);

        var pastLen = KmerAnalyzer.AnalyzeKmers(seq, seq.Length + 1); // k = 6 ⇒ 0 windows
        pastLen.Should().Be(new KmerStatistics(0, 0, 0, 0, 0, 0));
    }

    #endregion

    #region BE — Boundary: homopolymer (single distinct k-mer)

    /// <summary>
    /// BE / headline target: a HOMOPOLYMER ("AAAA…") has exactly ONE distinct
    /// k-mer "AA…A". For L bases and k ≤ L: Distinct = 1, Total = L − k + 1,
    /// Max = Min = Average = Total (the single k-mer's multiplicity), Entropy = 0
    /// (one-component distribution). K-mer_Statistics.md §6.1, INV-05.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_Homopolymer_OneDistinctKmerEntropyZero()
    {
        foreach (int len in new[] { 1, 2, 5, 10, 64 })
        {
            string seq = new string('A', len);
            foreach (int k in new[] { 1, 2, len })
            {
                if (k > len) continue;
                var stats = KmerAnalyzer.AnalyzeKmers(seq, k);
                int total = len - k + 1;

                stats.UniqueKmers.Should().Be(1,
                    $"homopolymer of length {len} has a single distinct {k}-mer");
                stats.TotalKmers.Should().Be(total, "Total = L − k + 1");
                stats.MaxCount.Should().Be(total, "the one k-mer occurs in every window");
                stats.MinCount.Should().Be(total, "Max == Min for a single distinct k-mer");
                stats.AverageCount.Should().Be(total, "Average = Total / 1 = Total");
                stats.Entropy.Should().BeApproximately(0.0, 1e-12,
                    "INV-05: one distinct k-mer ⇒ entropy 0");
                AssertWellFormed(stats, seq, k);
            }
        }
    }

    /// <summary>
    /// BE: a homopolymer over a NON-DNA character forms ordinary k-mers (no
    /// alphabet restriction; K-mer_Statistics.md §6.2). Still exactly one distinct
    /// k-mer and zero entropy — confirms the stats are alphabet-agnostic.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void AnalyzeKmers_NonDnaHomopolymer_StillSingleDistinctZeroEntropy()
    {
        string seq = new string('N', 12);
        var stats = KmerAnalyzer.AnalyzeKmers(seq, 4);

        stats.UniqueKmers.Should().Be(1, "no alphabet restriction — 'N' forms ordinary k-mers");
        stats.TotalKmers.Should().Be(9); // 12 − 4 + 1
        stats.Entropy.Should().BeApproximately(0.0, 1e-12);
        AssertWellFormed(stats, seq, 4);
    }

    #endregion

    #region BE — Boundary: k ≤ 0 floor (documented validation exception)

    /// <summary>
    /// BE: k ≤ 0 with NON-EMPTY input throws ArgumentOutOfRangeException(nameof(k)),
    /// surfaced unchanged from CountKmers — a 0-length k-mer is meaningless. This is
    /// the documented, intentional validation exception (K-mer_Statistics.md §3.3).
    /// </summary>
    [Test]
    public void AnalyzeKmers_NonPositiveK_NonEmptySequence_Throws()
    {
        foreach (int k in new[] { 0, -1, -100, int.MinValue })
        {
            Action act = () => KmerAnalyzer.AnalyzeKmers("ACGTACGT", k);
            act.Should().Throw<ArgumentOutOfRangeException>()
               .WithParameterName("k", $"k={k} ≤ 0 is rejected by the documented contract");
        }
    }

    #endregion

    #region BE/RB — Random fuzz: structural invariants hold for arbitrary inputs

    /// <summary>
    /// Random-input sweep over assorted lengths and k values (including k &gt; len
    /// and k = len boundaries): EVERY result must satisfy the well-formed oracle —
    /// finite non-negative fields, Total = max(L−k+1, 0), Min ≤ Avg ≤ Max,
    /// 0 ≤ Entropy ≤ log₂(distinct), no crash, no NaN. No undisciplined failure
    /// on any boundary k.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void AnalyzeKmers_RandomInputs_AlwaysWellFormed()
    {
        for (int seed = 0; seed < 250; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 60);
            string seq = RandomDna(len, seed);
            // Mix valid k, k = L, and k > L (never k ≤ 0 here — that path is the
            // throwing contract pinned separately).
            int k = rng.Next(1, len + 5);

            var stats = KmerAnalyzer.AnalyzeKmers(seq, k);
            AssertWellFormed(stats, seq, k);
        }
    }

    #endregion

    #region POSITIVE sanity — documented worked example and uniform-case bounds

    /// <summary>
    /// POSITIVE: the documented worked example (GTAGAGCTGT, K-mer_Statistics.md
    /// §7.1). k=1 ⇒ G=4,T=3,A=2,C=1 ⇒ Total=10, Distinct=4, Max=4, Min=1,
    /// Avg=2.5, Entropy ≈ 1.846439 bits. Hand-computed — pins the business contract,
    /// not just green.
    /// </summary>
    [Test]
    public void AnalyzeKmers_WorkedExample_K1_MatchesDocumentedStats()
    {
        var stats = KmerAnalyzer.AnalyzeKmers("GTAGAGCTGT", 1);

        stats.TotalKmers.Should().Be(10);
        stats.UniqueKmers.Should().Be(4);
        stats.MaxCount.Should().Be(4);
        stats.MinCount.Should().Be(1);
        stats.AverageCount.Should().Be(2.5);

        // E = −(0.4·log₂0.4 + 0.3·log₂0.3 + 0.2·log₂0.2 + 0.1·log₂0.1)
        double[] p = { 0.4, 0.3, 0.2, 0.1 };
        double expected = -p.Sum(x => x * Math.Log2(x));
        stats.Entropy.Should().BeApproximately(expected, 1e-9);
        stats.Entropy.Should().BeApproximately(1.846439, 1e-5, "documented value §7.1");
    }

    /// <summary>
    /// POSITIVE: the documented worked example, k=3 — a MAXIMALLY-DIVERSE case where
    /// all 8 windows are distinct ⇒ Distinct = Total = 8, Max = Min = 1, Avg = 1.0,
    /// Entropy = log₂8 = 3.0 bits (uniform distribution, INV-05 upper bound).
    /// K-mer_Statistics.md §7.1.
    /// </summary>
    [Test]
    public void AnalyzeKmers_WorkedExample_K3_AllWindowsDistinct()
    {
        var stats = KmerAnalyzer.AnalyzeKmers("GTAGAGCTGT", 3);

        stats.TotalKmers.Should().Be(8);
        stats.UniqueKmers.Should().Be(8, "all windows distinct ⇒ Distinct = Total");
        stats.MaxCount.Should().Be(1);
        stats.MinCount.Should().Be(1);
        stats.AverageCount.Should().Be(1.0);
        stats.Entropy.Should().BeApproximately(3.0, 1e-9, "uniform over 8 ⇒ log₂8 = 3 bits (INV-05 max)");
        stats.Entropy.Should().BeApproximately(Math.Log2(stats.UniqueKmers), 1e-9);
    }

    /// <summary>
    /// POSITIVE: case-insensitivity — lower-case input yields identical stats to its
    /// upper-cased form (input is upper-cased internally). K-mer_Statistics.md §6.1.
    /// </summary>
    [Test]
    public void AnalyzeKmers_LowerCase_IdenticalToUpperCase()
    {
        var lower = KmerAnalyzer.AnalyzeKmers("gtagagctgt", 2);
        var upper = KmerAnalyzer.AnalyzeKmers("GTAGAGCTGT", 2);

        lower.Should().Be(upper, "input is upper-cased internally ⇒ identical stats");
    }

    /// <summary>
    /// POSITIVE: a repeat-laden sequence pins the AverageCount = Total/Distinct and
    /// the entropy upper bound. "ATATAT" (L=6) with k=2 ⇒ k-mers AT×3, TA×2 ⇒
    /// Total=5, Distinct=2, Max=3, Min=2, Avg=2.5, Entropy = −(0.6·log₂0.6 +
    /// 0.4·log₂0.4) ≈ 0.970951 bits ≤ log₂2 = 1.
    /// </summary>
    [Test]
    public void AnalyzeKmers_RepeatLaden_MatchesHandComputedStats()
    {
        var stats = KmerAnalyzer.AnalyzeKmers("ATATAT", 2);

        stats.TotalKmers.Should().Be(5);   // 6 − 2 + 1
        stats.UniqueKmers.Should().Be(2);  // AT, TA
        stats.MaxCount.Should().Be(3);     // AT at 0,2,4
        stats.MinCount.Should().Be(2);     // TA at 1,3
        stats.AverageCount.Should().Be(2.5); // 5 / 2

        double expected = -(0.6 * Math.Log2(0.6) + 0.4 * Math.Log2(0.4));
        stats.Entropy.Should().BeApproximately(expected, 1e-9);
        stats.Entropy.Should().BeLessThan(Math.Log2(2) + 1e-9, "INV-05: ≤ log₂(distinct)");
    }

    #endregion

    #endregion
}
