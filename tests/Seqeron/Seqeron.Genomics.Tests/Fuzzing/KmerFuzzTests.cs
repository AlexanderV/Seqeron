using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the K-mer area — overlapping sliding-window k-mer counting
/// (KMER-COUNT-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (IndexOutOfRangeException, ArgumentOutOfRangeException
/// leaking from internal indexing, OutOfMemoryException). Every input must resolve
/// to EITHER a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception (ArgumentException / ArgumentOutOfRangeException). A raw
/// runtime exception, a hang, or a miscount on a boundary k value is a bug, not a
/// passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-COUNT-001 — k-mer counting
/// Checklist: docs/checklists/03_FUZZING.md, row 32.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate k / sequence-size boundaries
///          called out in the checklist row: k = 0, k &gt; seqLen, k = seqLen, the
///          empty sequence, and k = 1.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The k-mer-counting contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// K-mer counting tallies every overlapping length-k substring of a sequence. For
/// a sequence of length L it slides a window of length k across every valid start
/// position and increments a dictionary keyed by the (uppercased) window:
///   Count(w) = Σ_{i=0..L-k} 𝟙(sequence[i..i+k-1] = w)
/// (K-mer_Counting.md §2.2). The API entry under test is
///   KmerAnalyzer.CountKmers(string sequence, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 20–42),
/// the "canonical string-based counting" surface (K-mer_Counting.md §4.2). This is
/// the surface the metamorphic checklist also names for this row
/// (02_METAMORPHIC_TESTING.md row 32: "INV: total k-mer instances = seqLen − k + 1").
/// The DnaSequence overload delegates to this string path (KmerAnalyzer.cs lines
/// 105–108), so the empty-DnaSequence boundary is also probed here.
///
/// THE KEY INVARIANT (K-mer_Counting.md §2.4 INV-01): whenever k ≤ L the SUM of all
/// k-mer counts equals exactly the number of windows, L − k + 1 — the sliding window
/// emits exactly one k-mer per valid start position. Every positive-result test below
/// pins this window-count/sum invariant; it is the single load-bearing correctness
/// check that distinguishes a correct count from a miscount.
///
/// Documented parameter contract (K-mer_Counting.md §3.1, §3.3, §6.1; KmerAnalyzer.cs
/// lines 22–29):
///   • Null or empty sequence → returns the EMPTY dictionary, NOT an exception
///     (explicit `string.IsNullOrEmpty` guard, KmerAnalyzer.cs lines 22–23). This
///     guard runs BEFORE k is validated, so the empty/null short-circuit wins even
///     when k is itself degenerate.
///   • k ≤ 0 with NON-EMPTY input → ArgumentOutOfRangeException (nameof(k),
///     "K must be positive.", KmerAnalyzer.cs lines 25–26). This is THE k = 0
///     boundary the checklist row probes: a 0-length k-mer is biologically
///     meaningless (there is no "length-0 substring" to tally), so a *literal* k = 0
///     would make the window loop bound `i ≤ L − 0` run L+1 times emitting the empty
///     string "" at every position — a nonsense degenerate tally. The contract
///     REJECTS it outright rather than defining it. We pin that k = 0 (and any k &lt; 0)
///     throws on non-empty input — VERIFIED against the doc (§6.1: "k ≤ 0 with
///     non-empty string input → Throws ArgumentOutOfRangeException").
///   • k &gt; sequence.Length → returns the EMPTY dictionary (KmerAnalyzer.cs lines
///     28–29): no length-k window fits, so there are zero windows (L − k + 1 ≤ 0) and
///     the result is empty — never an out-of-range Substring, never a crash
///     (§6.1: "k > sequence.Length → Returns an empty dictionary; no valid windows
///     exist").
/// The implementation uppercases the input before keying (INV-03, case-insensitive),
/// and the raw-string path does NOT filter the alphabet — non-ACGT symbols are kept
/// as literal k-mer characters (§5.2, §6.1). These tests exercise only the boundary
/// k / size targets of THIS fuzz row; alphabet/injection targets are out of scope here.
///
/// The five checklist targets map to these documented behaviours:
///   • k = 0        → ArgumentOutOfRangeException on non-empty input (rejected, not defined).
///   • k &gt; seqLen   → empty dictionary (no windows), no crash.
///   • k = seqLen   → exactly ONE window (the whole sequence), count 1, sum = L − k + 1 = 1.
///   • empty seq    → empty dictionary (and empty seq + degenerate k still empty, not a throw).
///   • k = 1        → counts == single-base composition; sum = L (one window per base).
/// A positive-sanity test pins the window-count/sum invariant on a known sequence
/// (2-mers of "AAAA" → {AA:3}, sum = 4 − 2 + 1 = 3).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerFuzzTests
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
    //  KMER-COUNT-001 — k-mer counting : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-COUNT-001 — k-mer counting

    #region BE — Boundary: k = 0

    /// <summary>
    /// BE: k = 0 is the degenerate floor and the KEY meaningless-length target. A
    /// length-0 k-mer has no biological meaning — there is no "length-0 substring" to
    /// tally, and a literal k = 0 would run the window loop `i ≤ L − 0` for every
    /// position emitting the empty string "" each time, a nonsense degenerate count.
    /// The documented contract REJECTS k ≤ 0 on non-empty input with
    /// ArgumentOutOfRangeException (KmerAnalyzer.cs lines 25–26; K-mer_Counting.md
    /// §6.1). We pin that k = 0 throws and carries the documented "k" parameter name,
    /// so the floor cannot silently drift into defining a 0-length k-mer.
    /// </summary>
    [Test]
    public void CountKmers_KZero_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.CountKmers("ACGTACGT", 0);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; the contract rejects k <= 0 on non-empty input rather than defining it")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: a negative k is below the floor too and must be rejected the same way —
    /// pinning that the rejection boundary is exactly k ≤ 0, not merely k == 0, so a
    /// negative length can never slip through into the window loop.
    /// </summary>
    [Test]
    public void CountKmers_NegativeK_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.CountKmers("ACGTACGT", -3);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0 on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: k > sequence length

    /// <summary>
    /// BE: k far larger than the sequence length must NOT crash. No length-k window
    /// fits, so the window count L − k + 1 is ≤ 0 and the loop bound `i ≤ L − k` is
    /// negative — the loop never runs, no Substring is taken past the end, and the
    /// result is the empty dictionary (KmerAnalyzer.cs lines 28–29; K-mer_Counting.md
    /// §6.1). We pin no-throw AND emptiness so an oversized k can never index past the
    /// end nor invent a k-mer longer than the sequence.
    /// </summary>
    [Test]
    public void CountKmers_KGreaterThanSequenceLength_IsEmptyAndDoesNotThrow()
    {
        var act = () => KmerAnalyzer.CountKmers("ACGT", 1000);
        act.Should().NotThrow(
            "k > L makes the window count L − k + 1 negative; the loop never runs, so there is nothing to index past the end");

        KmerAnalyzer.CountKmers("ACGT", 1000).Should().BeEmpty(
            "no length-1000 window fits a 4-base sequence; the result is empty, not a crash");

        // k = L + 1 is the exact off-by-one boundary above the sequence length.
        KmerAnalyzer.CountKmers("ACGT", 5).Should().BeEmpty(
            "k = L + 1 is one past the last fitting window; still empty, never an out-of-range Substring");
    }

    #endregion

    #region BE — Boundary: k = sequence length

    /// <summary>
    /// BE: k = sequence length is the upper boundary where exactly ONE window fits:
    /// the whole sequence. The window count L − k + 1 = 1, so the single k-mer (the
    /// entire uppercased sequence) is reported with count 1 and the SUM of all counts
    /// is exactly 1 (INV-01, K-mer_Counting.md §2.4). This pins both the off-by-one
    /// upper edge (one window, not zero, not two) and the window-count/sum invariant
    /// at that edge.
    /// </summary>
    [Test]
    public void CountKmers_KEqualsSequenceLength_IsSingleWholeSequenceKmer()
    {
        const string seq = "ACGTA";
        var counts = KmerAnalyzer.CountKmers(seq, seq.Length);

        counts.Should().ContainSingle("k = L admits exactly one window — the whole sequence");
        counts.Should().ContainKey("ACGTA").WhoseValue.Should().Be(1, "the single whole-sequence k-mer occurs once");
        counts.Values.Sum().Should().Be(seq.Length - seq.Length + 1,
            "INV-01: the sum of all counts equals the window count L − k + 1 = 1 at k = L");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The string overload's
    /// `string.IsNullOrEmpty` guard short-circuits to the empty dictionary BEFORE k is
    /// validated (KmerAnalyzer.cs lines 22–23), so empty/null input NEVER throws — even
    /// when k is itself degenerate (k = 0). We pin that empty, null, and the
    /// empty-DnaSequence surface (which delegates to the string path) all return the
    /// empty dictionary with no exception, and that empty input wins over a degenerate
    /// k rather than throwing (K-mer_Counting.md §6.1: "Empty or null string input
    /// with k ≤ 0 → Returns an empty dictionary").
    /// </summary>
    [Test]
    public void CountKmers_EmptyOrNullSequence_IsEmptyAndDoesNotThrow()
    {
        var emptyAct = () => KmerAnalyzer.CountKmers(string.Empty, 3);
        var nullAct = () => KmerAnalyzer.CountKmers((string)null!, 3);
        var emptyDnaAct = () => KmerAnalyzer.CountKmers(new DnaSequence(string.Empty), 3);
        var emptyDegenerateKAct = () => KmerAnalyzer.CountKmers(string.Empty, 0);

        emptyAct.Should().NotThrow("an empty sequence has no windows; the explicit guard returns an empty dictionary");
        nullAct.Should().NotThrow("null input is treated as empty, not as an error, by the explicit guard");
        emptyDnaAct.Should().NotThrow("the DnaSequence overload delegates to the string path's empty guard");
        emptyDegenerateKAct.Should().NotThrow(
            "the empty/null guard runs BEFORE k is validated, so empty input wins even with a degenerate k = 0");

        KmerAnalyzer.CountKmers(string.Empty, 3).Should().BeEmpty();
        KmerAnalyzer.CountKmers((string)null!, 3).Should().BeEmpty();
        KmerAnalyzer.CountKmers(new DnaSequence(string.Empty), 3).Should().BeEmpty();
        KmerAnalyzer.CountKmers(string.Empty, 0).Should().BeEmpty(
            "empty input short-circuits to an empty dictionary before the k <= 0 check can throw");
    }

    #endregion

    #region BE — Boundary: k = 1

    /// <summary>
    /// BE: k = 1 is the smallest legal k-mer length and reduces counting to the
    /// single-base composition — one window per base, so the SUM of all counts equals
    /// L (= L − 1 + 1) and each key is a single uppercased base whose count is exactly
    /// how many times that base occurs (INV-01, K-mer_Counting.md §2.4). We pin the
    /// per-base counts AND the sum so the smallest-k edge cannot silently miscount the
    /// composition.
    /// </summary>
    [Test]
    public void CountKmers_KOne_EqualsSingleBaseComposition()
    {
        const string seq = "AACGTA"; // A:3, C:1, G:1, T:1 — length 6
        var counts = KmerAnalyzer.CountKmers(seq, 1);

        counts.Should().HaveCount(4, "four distinct single bases appear in 'AACGTA'");
        counts["A"].Should().Be(3, "'A' occurs three times");
        counts["C"].Should().Be(1);
        counts["G"].Should().Be(1);
        counts["T"].Should().Be(1);
        counts.Values.Sum().Should().Be(seq.Length,
            "INV-01: at k = 1 the sum of all counts equals L − 1 + 1 = L (one window per base)");
    }

    #endregion

    #region Positive sanity — the window-count/sum invariant on a known sequence

    /// <summary>
    /// Positive sanity: the textbook homopolymer count from the algorithm doc —
    /// 2-mers of "AAAA" — must yield the single key "AA" with count L − k + 1 = 4 − 2 +
    /// 1 = 3 (K-mer_Counting.md §6.1: "Homopolymer such as AAAA, k = 2 → one key with
    /// count L − k + 1"). This pins the KEY window-count/sum invariant (INV-01) on a
    /// known sequence so the boundary hardening never comes at the cost of the core
    /// count silently breaking.
    /// </summary>
    [Test]
    public void CountKmers_HomopolymerTwoMers_HasSingleKeyWithWindowCount()
    {
        var counts = KmerAnalyzer.CountKmers("AAAA", 2);

        counts.Should().ContainSingle("every length-2 window of 'AAAA' is the identical k-mer 'AA'");
        counts.Should().ContainKey("AA").WhoseValue.Should().Be(3, "L − k + 1 = 4 − 2 + 1 = 3 overlapping 'AA' windows");
        counts.Values.Sum().Should().Be(3, "INV-01: the sum of all counts equals the window count L − k + 1 = 3");
    }

    /// <summary>
    /// Positive sanity: a heterogeneous sequence with a known multi-key distribution.
    /// "ACGTACGT" (L = 8), k = 3 → 8 − 3 + 1 = 6 windows: ACG, CGT, GTA, TAC, ACG, CGT.
    /// 'ACG' and 'CGT' each occur twice; the other four occur once. We pin every count
    /// explicitly AND the sum-equals-window-count invariant (INV-01).
    /// </summary>
    [Test]
    public void CountKmers_HeterogeneousThreeMers_MatchKnownDistributionAndWindowCount()
    {
        const string seq = "ACGTACGT";
        const int k = 3;
        var counts = KmerAnalyzer.CountKmers(seq, k);

        counts["ACG"].Should().Be(2, "'ACG' starts at positions 0 and 4");
        counts["CGT"].Should().Be(2, "'CGT' starts at positions 1 and 5");
        counts["GTA"].Should().Be(1);
        counts["TAC"].Should().Be(1);
        counts.Should().HaveCount(4, "four distinct 3-mers appear in 'ACGTACGT'");
        counts.Values.Sum().Should().Be(seq.Length - k + 1,
            "INV-01: the sum of all counts equals the window count L − k + 1 = 6");
    }

    /// <summary>
    /// Positive sanity: case-insensitivity (INV-03). A lowercase sequence must produce
    /// the same uppercased keys and the same counts as its uppercase form, because the
    /// implementation normalizes to uppercase before keying (KmerAnalyzer.cs line 31).
    /// </summary>
    [Test]
    public void CountKmers_LowercaseInput_IsUppercasedBeforeCounting()
    {
        var lower = KmerAnalyzer.CountKmers("acgtacgt", 3);
        var upper = KmerAnalyzer.CountKmers("ACGTACGT", 3);

        lower.Should().BeEquivalentTo(upper,
            "INV-03: input is uppercased before keying, so case does not change the counts");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// satisfy the KEY window-count invariant for several k values — the sum of all
    /// counts must equal L − k + 1 whenever k ≤ L (INV-01, K-mer_Counting.md §2.4),
    /// regardless of the random content. [CancelAfter] guards against any hang on the
    /// largest k scanned.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void CountKmers_RandomSequence_SumEqualsWindowCountForEveryK()
    {
        const int length = 2000;
        string seq = RandomDna(length, seed: 32_001);

        foreach (int k in new[] { 1, 2, 3, 5, 8, 13 })
        {
            var counts = KmerAnalyzer.CountKmers(seq, k);

            counts.Values.Sum().Should().Be(length - k + 1,
                $"INV-01: at k = {k} the sum of all counts equals the window count L − k + 1");
            counts.Keys.Should().OnlyContain(key => key.Length == k,
                $"every emitted key is a length-{k} window — no shorter/longer fragments");
        }
    }

    #endregion

    #endregion
}
