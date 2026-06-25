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
/// (KMER-COUNT-001), normalized k-mer frequencies (KMER-FREQ-001), and the
/// over-represented / frequent-k-mer filter (KMER-FIND-001).
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
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-FREQ-001 — k-mer frequencies (normalized counts)
/// Checklist: docs/checklists/03_FUZZING.md, row 33.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate k / sequence-size boundaries
///          called out in the checklist row: k = 0, k &gt; seqLen, the empty
///          sequence, and the single-character sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The k-mer-frequency contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A k-mer FREQUENCY is a normalized count — the count of a k-mer divided by the
/// total number of k-mer windows — i.e. a probability over the observed k-mer
/// distribution: fᵢ = cᵢ / Σⱼ cⱼ (K-mer_Frequency_Analysis.md §2.2). The API entry
/// under test is
///   KmerAnalyzer.GetKmerFrequencies(string sequence, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 177–189),
/// which counts via CountKmers, sums the counts into `total`, and divides each count
/// by `total`.
///
/// THE KEY INVARIANT (K-mer_Frequency_Analysis.md §2.4 INV-01): whenever at least
/// one k-mer window exists, the SUM of all returned frequencies is exactly 1.0 (each
/// count is divided by the same total). Every positive-result test below pins this
/// sum-to-1 invariant alongside non-negativity; it is the single load-bearing
/// correctness check that distinguishes a true probability distribution from a
/// miscount or a mis-normalization.
///
/// THE KEY FUZZ CONCERN — the 0-window division guard. The frequency is count/total.
/// At zero windows (empty sequence, or k &gt; sequence length) the count dictionary is
/// empty, so `total == 0`. A naïve cᵢ/total would be a 0/0 → NaN or a
/// DivideByZeroException. The implementation guards this explicitly: `if (total == 0)
/// return new Dictionary<string, double>()` (KmerAnalyzer.cs lines 182–183), returning
/// the EMPTY distribution — never NaN, never a DivideByZeroException, never a crash
/// (K-mer_Frequency_Analysis.md §3.3, §6.1: "Empty sequence / k &gt; sequence.Length →
/// Empty frequency map"). These tests pin that guard at every 0-window boundary.
///
/// Documented parameter contract (K-mer_Frequency_Analysis.md §3.1, §3.3, §6.1):
///   • k = 0 (and k &lt; 0) on non-empty input → ArgumentOutOfRangeException, surfaced
///     unchanged from the underlying CountKmers k ≤ 0 rejection (KmerAnalyzer.cs
///     lines 25–26): a 0-length k-mer is meaningless, so frequencies of it are too.
///   • k &gt; sequence.Length → empty frequency map (0 windows → total == 0 → guard).
///   • empty / null sequence → empty frequency map (0 windows → total == 0 → guard);
///     the CountKmers empty/null guard runs BEFORE k is validated, so empty input
///     wins even when k is itself degenerate (k = 0 on empty → still empty, no throw).
///   • single-character sequence with k = 1 → exactly one window → one k-mer with
///     frequency 1.0 (a distribution with a single outcome; §6.1 "Single possible
///     k-mer → Frequency 1.0").
/// The implementation uppercases before keying (it counts via CountKmers, INV-03),
/// so frequencies are case-insensitive. These tests exercise only the boundary
/// k / size targets of THIS fuzz row.
///
/// The four checklist targets map to these documented behaviours:
///   • k = 0          → ArgumentOutOfRangeException on non-empty input (rejected).
///   • k &gt; seqLen     → empty frequency map; no DivideByZero, no NaN.
///   • empty seq      → empty frequency map; no DivideByZero, no NaN (degenerate k too).
///   • single char    → one k-mer with frequency 1.0 (k = 1).
/// A positive-sanity test pins the sum-to-1 invariant AND the exact frequency ratios
/// on a known sequence ("ACGTACGT", k = 3 → {ACG:2/6, CGT:2/6, GTA:1/6, TAC:1/6}).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-FIND-001 — finding frequent / over-represented k-mers
/// Checklist: docs/checklists/03_FUZZING.md, row 34 (the LAST K-mer unit in the
/// first block).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate k / threshold / sequence-size
///          boundaries called out in the checklist row: k = 0, minFreq = 0,
///          minFreq &gt; the maximum achievable count, and the empty sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The frequent-k-mer contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// "Finding frequent k-mers" is the recurrent / over-represented filter: the
/// k-mers whose overlapping occurrence count is AT LEAST a threshold t, i.e. the
/// solution to the Frequent Words / Count(Text, Pattern) ≥ t problem (Compeau &amp;
/// Pevzner, <i>Bioinformatics Algorithms</i>; Rosalind BA1B; K-mer_Search.md §2).
/// The API entry under test is
///   KmerAnalyzer.FindKmersWithMinCount(string sequence, int k, int minCount)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 274–282).
/// It counts via CountKmers, keeps the (k-mer, count) pairs whose count ≥ minCount,
/// and orders them by count DESCENDING (KmerAnalyzer.cs lines 277–281; the doc
/// comment there names minCount the "Inclusive minimum occurrence count threshold").
///
/// minFreq SEMANTICS — minFreq is a raw COUNT, NOT a frequency fraction. The
/// parameter is `int minCount`, an inclusive occurrence-count threshold (count ≥
/// minCount), NOT a normalized 0.0–1.0 frequency. So "minFreq = 0" means a count
/// threshold of 0, under which EVERY observed k-mer qualifies (every observed count
/// is ≥ 1 ≥ 0); it is NOT a 0% fraction. This is the load-bearing fact the whole
/// row turns on, and the tests pin it explicitly.
///
/// THE KEY CONTRACT (soundness + completeness, K-mer_Search.md §2.2 / §2.4 model):
/// the result is EXACTLY {(w, c) : c = Count(w) and c ≥ minCount}. Soundness — every
/// returned k-mer truly meets the threshold (Count ≥ minCount); completeness — no
/// above-threshold k-mer is omitted. Every positive-result test below cross-checks
/// the returned set against an INDEPENDENT CountKmers map: the returned (k-mer,count)
/// pairs must equal exactly the count-map entries with value ≥ minCount, the reported
/// counts must match the map, and the counts must be non-increasing (the documented
/// descending order). This is the single load-bearing correctness check.
///
/// DEFERRED-ENUMERATION NOTE: FindKmersWithMinCount is a LINQ pipeline
/// (Where/OrderByDescending) over CountKmers, so the k ≤ 0 rejection from CountKmers
/// only surfaces when the result is enumerated. The throw-asserting tests therefore
/// force enumeration (ToList) inside the asserted delegate — the same pattern the
/// existing KMER-UNIQUE-001 suite uses — so the documented exception is actually
/// observed, not silently deferred away.
///
/// Documented parameter contract (KmerAnalyzer.cs lines 264–273; K-mer_Search.md
/// §3.3, §6.1):
///   • k = 0 (and k &lt; 0) on non-empty input → ArgumentOutOfRangeException
///     (nameof(k)), surfaced unchanged from CountKmers (KmerAnalyzer.cs lines 25–26):
///     a 0-length k-mer is meaningless, so "frequent 0-mers" are too. (§6.1 "k ≤ 0 →
///     Throws for count-based searches".)
///   • minCount = 0 → every observed k-mer qualifies — every overlapping count is
///     ≥ 1 ≥ 0, so the filter c ≥ 0 admits the entire distinct-k-mer set. The doc
///     states "With minCount ≤ 1 every distinct k-mer qualifies" (KmerAnalyzer.cs
///     lines 270–271); minCount = 0 is below that floor and behaves identically.
///     This is the count-vs-fraction probe: a 0 COUNT threshold returns ALL k-mers,
///     not "0% of them".
///   • minCount &gt; the maximum achievable count (e.g. &gt; n − k + 1, the count an
///     all-identical run would reach) → NO k-mer qualifies → EMPTY result, never a
///     crash. The filter c ≥ minCount selects nothing.
///   • empty / null sequence → CountKmers short-circuits to the empty count map
///     before k is validated (KmerAnalyzer.cs lines 22–23), so the filter sees zero
///     k-mers and returns EMPTY — no k-mers exist (§6.1 "Empty sequence → empty
///     result"). Empty input wins even when k is itself degenerate (k = 0 on empty →
///     still empty, NOT a throw), because the empty guard runs before the k check.
///   • k &gt; sequence.Length → empty count map → empty result (no windows).
/// The implementation uppercases before keying (it counts via CountKmers), so the
/// filter is case-insensitive. These tests exercise only the BE targets of THIS row.
///
/// The four checklist targets map to these documented behaviours:
///   • k = 0              → ArgumentOutOfRangeException on non-empty input (rejected).
///   • minFreq = 0        → every observed k-mer qualifies (count ≥ 0 admits all);
///                          NOT a 0% fraction — proves minFreq is a COUNT.
///   • minFreq &gt; possible → empty result; the c ≥ minCount filter selects nothing.
///   • empty seq          → empty result; no k-mers (degenerate k on empty too).
/// A positive-sanity test pins soundness + completeness on a known sequence
/// ("ACGTACGT", k = 4, minCount = 2 → only {ACGT:2} qualifies; CGTA/GTAC/TACG are
/// below threshold and excluded), cross-checked against an independent CountKmers map.
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

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-FREQ-001 — k-mer frequencies (normalized counts) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-FREQ-001 — k-mer frequencies

    #region BE — Boundary: k = 0

    /// <summary>
    /// BE: k = 0 is the degenerate floor and a meaningless k-mer length — there is no
    /// "length-0 substring" whose frequency could be defined. GetKmerFrequencies counts
    /// via CountKmers, which rejects k ≤ 0 on non-empty input with
    /// ArgumentOutOfRangeException (KmerAnalyzer.cs lines 25–26); the frequency layer
    /// adds no k-handling of its own, so that rejection surfaces unchanged
    /// (K-mer_Frequency_Analysis.md §3.3: "If k ≤ 0, the underlying counting routine
    /// throws ArgumentOutOfRangeException"). We pin that k = 0 throws and carries the
    /// documented "k" parameter name, so a 0-length k-mer can never reach the
    /// count/total division.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_KZero_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.GetKmerFrequencies("ACGTACGT", 0);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; the frequency layer surfaces the underlying k <= 0 rejection on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: a negative k is below the floor too and must be rejected the same way —
    /// pinning that the rejection boundary is exactly k ≤ 0 (surfaced from CountKmers),
    /// not merely k == 0, so a negative length can never slip into the normalization.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_NegativeK_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.GetKmerFrequencies("ACGTACGT", -3);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0 on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: k > sequence length (the 0-window division guard)

    /// <summary>
    /// BE / KEY: k far larger than the sequence length yields ZERO windows — the count
    /// dictionary is empty, so total = Σ counts = 0. A frequency is count/total, so a
    /// naïve normalization here would be 0/0 → NaN or a DivideByZeroException. The
    /// explicit `if (total == 0) return new Dictionary<string,double>()` guard
    /// (KmerAnalyzer.cs lines 182–183) instead returns the EMPTY distribution. We pin
    /// no-throw AND emptiness so an oversized k can never divide by zero, never emit a
    /// NaN frequency, and never index past the end (K-mer_Frequency_Analysis.md §6.1:
    /// "k > sequence.Length → Empty frequency map").
    /// </summary>
    [Test]
    public void GetKmerFrequencies_KGreaterThanSequenceLength_IsEmptyNoDivideByZero()
    {
        var act = () => KmerAnalyzer.GetKmerFrequencies("ACGT", 1000);
        act.Should().NotThrow(
            "k > L gives zero windows so total == 0; the guard returns an empty map instead of dividing by zero");

        var freqs = KmerAnalyzer.GetKmerFrequencies("ACGT", 1000);
        freqs.Should().BeEmpty("no length-1000 window fits a 4-base sequence; the distribution is empty, not a crash");
        freqs.Values.Should().NotContain(double.NaN, "the 0-window guard avoids 0/0; no NaN frequency is ever produced");

        // k = L + 1 is the exact off-by-one boundary above the sequence length.
        KmerAnalyzer.GetKmerFrequencies("ACGT", 5).Should().BeEmpty(
            "k = L + 1 is one past the last fitting window; still empty, still no division by zero");
    }

    #endregion

    #region BE — Boundary: empty sequence (the 0-window division guard)

    /// <summary>
    /// BE / KEY: the empty sequence is the lower size boundary and the other 0-window
    /// source. CountKmers short-circuits empty/null input to the empty dictionary
    /// (KmerAnalyzer.cs lines 22–23), so total == 0 and the frequency layer's guard
    /// returns the empty distribution — no division by zero, no NaN, no throw. Because
    /// the CountKmers empty/null guard runs BEFORE k is validated, empty input wins even
    /// when k is itself degenerate (k = 0 on empty → still an empty map, not a throw).
    /// We pin that empty, null, and empty-with-degenerate-k all return an empty,
    /// NaN-free distribution (K-mer_Frequency_Analysis.md §6.1: "Empty sequence → Empty
    /// frequency map").
    /// </summary>
    [Test]
    public void GetKmerFrequencies_EmptyOrNullSequence_IsEmptyNoDivideByZero()
    {
        var emptyAct = () => KmerAnalyzer.GetKmerFrequencies(string.Empty, 3);
        var nullAct = () => KmerAnalyzer.GetKmerFrequencies(null!, 3);
        var emptyDegenerateKAct = () => KmerAnalyzer.GetKmerFrequencies(string.Empty, 0);

        emptyAct.Should().NotThrow(
            "an empty sequence has no windows so total == 0; the guard returns an empty map, not a 0/0 division");
        nullAct.Should().NotThrow("null input is treated as empty by the underlying count guard, not as an error");
        emptyDegenerateKAct.Should().NotThrow(
            "the empty/null guard runs BEFORE k is validated, so empty input wins even with a degenerate k = 0");

        var emptyFreqs = KmerAnalyzer.GetKmerFrequencies(string.Empty, 3);
        emptyFreqs.Should().BeEmpty();
        emptyFreqs.Values.Should().NotContain(double.NaN, "the 0-window guard avoids 0/0; no NaN frequency is produced");
        KmerAnalyzer.GetKmerFrequencies(null!, 3).Should().BeEmpty();
        KmerAnalyzer.GetKmerFrequencies(string.Empty, 0).Should().BeEmpty(
            "empty input short-circuits to an empty map before the k <= 0 check can throw and before any division");
    }

    #endregion

    #region BE — Boundary: single-character sequence

    /// <summary>
    /// BE: a single-character sequence with k = 1 is the minimal non-empty case — it
    /// admits exactly ONE window, so the distribution has a single outcome whose
    /// frequency is 1.0 (count 1 / total 1). This is the §6.1 "Single possible k-mer →
    /// Frequency 1.0" edge: the sum-to-1 invariant (INV-01) holds with a single mass
    /// concentrated on one k-mer, and no other key appears. We also confirm the k = 1
    /// frequency of a longer single-base run still concentrates on that one base.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_SingleCharSequence_IsSingleKmerFrequencyOne()
    {
        var single = KmerAnalyzer.GetKmerFrequencies("A", 1);

        single.Should().ContainSingle("a 1-base sequence at k = 1 admits exactly one window");
        single.Should().ContainKey("A").WhoseValue.Should().Be(1.0,
            "the single observed k-mer carries all of the probability mass (count 1 / total 1)");
        single.Values.Sum().Should().BeApproximately(1.0, 1e-9,
            "INV-01: with one outcome the sole frequency is 1.0, so the distribution sums to 1");

        // A single-character RUN (homopolymer) at k = 1 still concentrates on one base.
        var run = KmerAnalyzer.GetKmerFrequencies("AAAA", 1);
        run.Should().ContainSingle("'AAAA' has a single distinct 1-mer");
        run["A"].Should().BeApproximately(1.0, 1e-9, "all four windows are 'A', so its frequency is 4/4 = 1.0");
    }

    #endregion

    #region Positive sanity — the sum-to-1 invariant and exact ratios on known sequences

    /// <summary>
    /// Positive sanity: the textbook homopolymer frequency from the algorithm doc —
    /// 2-mers of "AAAA" — must yield the single key "AA" with frequency 1.0 (all three
    /// windows are the identical k-mer, so 3/3 = 1.0), and the distribution sums to 1
    /// (K-mer_Frequency_Analysis.md §6.1: "Homopolymer such as AAAA, k = 2 →
    /// Frequency {\"AA\": 1.0}"). This pins INV-01 on a known sequence so the boundary
    /// hardening never comes at the cost of the core normalization silently breaking.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_HomopolymerTwoMers_SingleKeyFrequencyOne()
    {
        var freqs = KmerAnalyzer.GetKmerFrequencies("AAAA", 2);

        freqs.Should().ContainSingle("every length-2 window of 'AAAA' is the identical k-mer 'AA'");
        freqs.Should().ContainKey("AA").WhoseValue.Should().BeApproximately(1.0, 1e-9,
            "all L − k + 1 = 3 windows are 'AA', so its frequency is 3/3 = 1.0");
        freqs.Values.Sum().Should().BeApproximately(1.0, 1e-9, "INV-01: the frequencies of a non-empty distribution sum to 1");
    }

    /// <summary>
    /// Positive sanity: a heterogeneous sequence with a known multi-key distribution.
    /// "ACGTACGT" (L = 8), k = 3 → 8 − 3 + 1 = 6 windows: ACG, CGT, GTA, TAC, ACG, CGT.
    /// 'ACG' and 'CGT' occur twice (frequency 2/6 each); 'GTA' and 'TAC' once (1/6 each).
    /// We pin every frequency explicitly, non-negativity, AND the KEY sum-to-1 invariant
    /// (INV-01) — the load-bearing check that the counts were normalized into a true
    /// probability distribution.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_HeterogeneousThreeMers_MatchKnownRatiosAndSumToOne()
    {
        const string seq = "ACGTACGT";
        const int k = 3;
        var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);

        const double total = 6.0; // L − k + 1
        freqs.Should().HaveCount(4, "four distinct 3-mers appear in 'ACGTACGT'");
        freqs["ACG"].Should().BeApproximately(2.0 / total, 1e-9, "'ACG' starts at positions 0 and 4 → 2/6");
        freqs["CGT"].Should().BeApproximately(2.0 / total, 1e-9, "'CGT' starts at positions 1 and 5 → 2/6");
        freqs["GTA"].Should().BeApproximately(1.0 / total, 1e-9);
        freqs["TAC"].Should().BeApproximately(1.0 / total, 1e-9);
        freqs.Values.Should().OnlyContain(f => f > 0.0, "every observed k-mer has a strictly positive frequency");
        freqs.Values.Sum().Should().BeApproximately(1.0, 1e-9,
            "INV-01: the sum of all returned frequencies is exactly 1.0 whenever a window exists");
    }

    /// <summary>
    /// Positive sanity: case-insensitivity (INV-03). A lowercase sequence must produce
    /// the same uppercased keys and the same frequencies as its uppercase form, because
    /// GetKmerFrequencies counts via CountKmers, which uppercases before keying.
    /// </summary>
    [Test]
    public void GetKmerFrequencies_LowercaseInput_IsUppercasedBeforeNormalizing()
    {
        var lower = KmerAnalyzer.GetKmerFrequencies("acgtacgt", 3);
        var upper = KmerAnalyzer.GetKmerFrequencies("ACGTACGT", 3);

        lower.Should().BeEquivalentTo(upper,
            "INV-03: input is uppercased before counting, so case does not change the frequencies");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// satisfy the KEY frequency invariant for several k values — the sum of all
    /// frequencies must be 1.0 (INV-01) and every frequency must be non-negative and
    /// NaN-free, regardless of the random content. [CancelAfter] guards against any hang
    /// on the largest k scanned.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void GetKmerFrequencies_RandomSequence_SumsToOneForEveryK()
    {
        const int length = 2000;
        string seq = RandomDna(length, seed: 33_001);

        foreach (int k in new[] { 1, 2, 3, 5, 8, 13 })
        {
            var freqs = KmerAnalyzer.GetKmerFrequencies(seq, k);

            freqs.Values.Should().OnlyContain(f => f >= 0.0 && !double.IsNaN(f),
                $"at k = {k} every frequency is a valid non-negative probability");
            freqs.Values.Sum().Should().BeApproximately(1.0, 1e-9,
                $"INV-01: at k = {k} the sum of all frequencies is exactly 1.0");
            freqs.Keys.Should().OnlyContain(key => key.Length == k,
                $"every emitted key is a length-{k} window — no shorter/longer fragments");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-FIND-001 — finding frequent k-mers (count ≥ minCount) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-FIND-001 — finding frequent k-mers

    #region BE — Boundary: k = 0

    /// <summary>
    /// BE: k = 0 is the degenerate floor and a meaningless k-mer length — there are no
    /// "length-0 k-mers" whose frequency could be over-represented. FindKmersWithMinCount
    /// counts via CountKmers, which rejects k ≤ 0 on non-empty input with
    /// ArgumentOutOfRangeException (KmerAnalyzer.cs lines 25–26; K-mer_Search.md §6.1
    /// "k ≤ 0 → Throws for count-based searches"). Because the method is a deferred LINQ
    /// pipeline, we force enumeration (ToList) inside the asserted delegate so the throw
    /// actually surfaces. We pin that k = 0 throws and carries the documented "k"
    /// parameter name, so a 0-length k-mer can never reach the count filter.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_KZero_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.FindKmersWithMinCount("ACGTACGT", 0, 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; the frequent-k-mer filter surfaces the underlying k <= 0 rejection on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: a negative k is below the floor too and must be rejected the same way —
    /// pinning that the rejection boundary is exactly k ≤ 0 (surfaced from CountKmers),
    /// not merely k == 0, so a negative length can never slip into the count filter.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_NegativeK_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.FindKmersWithMinCount("ACGTACGT", -3, 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0 on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: minFreq = 0 (the count-vs-fraction probe)

    /// <summary>
    /// BE / KEY: minFreq = 0 is the degenerate threshold floor and the count-vs-fraction
    /// probe. minFreq is a raw COUNT (the `int minCount` "inclusive minimum occurrence
    /// count threshold", KmerAnalyzer.cs lines 270–271), NOT a 0.0–1.0 frequency fraction.
    /// So a threshold of 0 means "count ≥ 0", which EVERY observed k-mer satisfies (every
    /// overlapping count is ≥ 1 ≥ 0) — the filter must return the ENTIRE distinct-k-mer
    /// set, NOT "0% of them" and NOT an empty set. We pin completeness against an
    /// independent CountKmers map: the returned (k-mer, count) pairs must equal exactly
    /// the full count map, with matching counts and the documented descending order.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_MinFreqZero_ReturnsEveryObservedKmer()
    {
        const string seq = "ACGTACGT";
        const int k = 3;
        var counts = KmerAnalyzer.CountKmers(seq, k); // independent oracle

        var result = KmerAnalyzer.FindKmersWithMinCount(seq, k, 0).ToList();

        result.Select(p => p.Kmer).Should().BeEquivalentTo(counts.Keys,
            "minCount = 0 is a COUNT threshold of zero (not a 0% fraction); every observed k-mer has count >= 1 >= 0, so ALL distinct k-mers qualify");
        result.Should().OnlyContain(p => p.Count == counts[p.Kmer],
            "soundness: each reported count must equal the independent CountKmers value for that k-mer");
        result.Select(p => p.Count).Should().BeInDescendingOrder(
            "the documented contract orders the qualifying k-mers by count descending");
    }

    /// <summary>
    /// BE: minFreq = 1 is the off-by-one neighbour just above the 0 floor and the
    /// documented "minCount ≤ 1 → every distinct k-mer qualifies" edge (KmerAnalyzer.cs
    /// lines 270–271): since every observed count is ≥ 1, the c ≥ 1 filter STILL admits
    /// the whole distinct-k-mer set, identically to minCount = 0. Pinning both 0 and 1
    /// fixes that the "all-qualify" region is exactly minCount ≤ 1.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_MinFreqOne_AlsoReturnsEveryObservedKmer()
    {
        const string seq = "ACGTACGT";
        const int k = 3;
        var counts = KmerAnalyzer.CountKmers(seq, k);

        var atZero = KmerAnalyzer.FindKmersWithMinCount(seq, k, 0).Select(p => p.Kmer);
        var atOne = KmerAnalyzer.FindKmersWithMinCount(seq, k, 1).Select(p => p.Kmer);

        atOne.Should().BeEquivalentTo(counts.Keys,
            "every observed count is >= 1, so minCount = 1 still admits the whole distinct-k-mer set");
        atOne.Should().BeEquivalentTo(atZero,
            "minCount = 0 and minCount = 1 select the identical set — the 'all qualify' region is exactly minCount <= 1");
    }

    #endregion

    #region BE — Boundary: minFreq > the maximum achievable count

    /// <summary>
    /// BE: a threshold ABOVE the maximum achievable count must select nothing without
    /// crashing. The largest count any k-mer can reach is the window count n − k + 1
    /// (attained only by an all-identical run); a threshold strictly greater than that is
    /// unsatisfiable, so the c ≥ minCount filter yields the EMPTY result. We probe two
    /// flavours: a threshold above the actual max observed count in a heterogeneous
    /// sequence (max 2, threshold 3 → empty), and an absurdly large threshold above the
    /// absolute ceiling (threshold &gt; n − k + 1 → empty). Neither may throw.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_MinFreqAboveMaxPossible_IsEmptyAndDoesNotThrow()
    {
        const string seq = "ACGTACGT"; // k = 4: max count is 2 (ACGT); n − k + 1 = 5 windows.

        var aboveObservedMax = () => KmerAnalyzer.FindKmersWithMinCount(seq, 4, 3).ToList();
        aboveObservedMax.Should().NotThrow(
            "a threshold above the maximum observed count is simply unsatisfiable; the filter selects nothing, it does not crash");
        aboveObservedMax().Should().BeEmpty(
            "no 4-mer occurs 3 times in 'ACGTACGT' (max is 2), so minCount = 3 selects nothing");

        // Absurd threshold far above the absolute ceiling n − k + 1 = 5.
        KmerAnalyzer.FindKmersWithMinCount(seq, 4, 1_000_000).ToList().Should().BeEmpty(
            "no k-mer can occur 1,000,000 times in an 8-base sequence (ceiling n − k + 1 = 5); the result is empty, never a crash");

        // Even a homopolymer's single k-mer maxes out at the window count; one above is empty.
        const string run = "AAAAAA"; // k = 2: AA occurs n − k + 1 = 5 times — the ceiling.
        KmerAnalyzer.FindKmersWithMinCount(run, 2, 5).Should().ContainSingle(
            "AA reaches exactly the ceiling count 5, so minCount = 5 still includes it");
        KmerAnalyzer.FindKmersWithMinCount(run, 2, 6).ToList().Should().BeEmpty(
            "minCount = 6 is one above the ceiling count 5; no k-mer can reach it, so the result is empty");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary and has no k-mers at all, so the
    /// frequent-k-mer filter has nothing to select and returns EMPTY. CountKmers
    /// short-circuits empty/null input to the empty map BEFORE k is validated
    /// (KmerAnalyzer.cs lines 22–23), so empty input NEVER throws — even when k is itself
    /// degenerate (k = 0 on empty → still empty, NOT a throw) and even with a degenerate
    /// minFreq. We pin that empty, null, k &gt; length, and empty-with-degenerate-k all
    /// return the empty result with no exception (K-mer_Search.md §6.1: "Empty sequence →
    /// Returns an empty result").
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_EmptyOrNullSequence_IsEmptyAndDoesNotThrow()
    {
        var emptyAct = () => KmerAnalyzer.FindKmersWithMinCount(string.Empty, 3, 2).ToList();
        var nullAct = () => KmerAnalyzer.FindKmersWithMinCount(null!, 3, 2).ToList();
        var emptyDegenerateKAct = () => KmerAnalyzer.FindKmersWithMinCount(string.Empty, 0, 2).ToList();
        var emptyDegenerateMinAct = () => KmerAnalyzer.FindKmersWithMinCount(string.Empty, 3, 0).ToList();

        emptyAct.Should().NotThrow("an empty sequence has no k-mers; the count short-circuit yields an empty filter result");
        nullAct.Should().NotThrow("null input is treated as empty by the underlying count guard, not as an error");
        emptyDegenerateKAct.Should().NotThrow(
            "the empty/null guard runs BEFORE k is validated, so empty input wins even with a degenerate k = 0");
        emptyDegenerateMinAct.Should().NotThrow("an empty sequence with any minCount still yields an empty result");

        KmerAnalyzer.FindKmersWithMinCount(string.Empty, 3, 2).Should().BeEmpty();
        KmerAnalyzer.FindKmersWithMinCount(null!, 3, 2).Should().BeEmpty();
        KmerAnalyzer.FindKmersWithMinCount(string.Empty, 0, 2).Should().BeEmpty(
            "empty input short-circuits to an empty result before the k <= 0 check can throw");

        // k > sequence length is the other 0-window source: no windows fit, so empty.
        KmerAnalyzer.FindKmersWithMinCount("ACG", 5, 1).Should().BeEmpty(
            "k = 5 > length 3 gives no valid windows, so no frequent k-mer exists");
    }

    #endregion

    #region Positive sanity — soundness + completeness on a known sequence

    /// <summary>
    /// Positive sanity: the worked recurrent-filter example — 'ACGTACGT', k = 4,
    /// minCount = 2. The six 4-mers are ACGT (pos 0, 4), CGTA, GTAC, TACG; only ACGT
    /// reaches count 2, so the frequent set is exactly {ACGT:2} and the three below-
    /// threshold k-mers are excluded. This pins BOTH halves of the contract: soundness
    /// (every returned k-mer truly has count ≥ minCount) and completeness (no above-
    /// threshold k-mer omitted), cross-checked against an independent CountKmers map so
    /// the boundary hardening never comes at the cost of the core filter silently breaking.
    /// </summary>
    [Test]
    public void FindKmersWithMinCount_KnownRecurrent_IsSoundAndComplete()
    {
        const string seq = "ACGTACGT";
        const int k = 4;
        const int minCount = 2;
        var counts = KmerAnalyzer.CountKmers(seq, k); // independent oracle

        var result = KmerAnalyzer.FindKmersWithMinCount(seq, k, minCount).ToList();

        // Soundness: every returned k-mer actually meets the threshold.
        result.Should().OnlyContain(p => p.Count >= minCount && p.Count == counts[p.Kmer],
            "soundness: each returned k-mer truly occurs >= minCount times, with the count matching the independent map");

        // Completeness: no above-threshold k-mer is omitted.
        var expected = counts.Where(kvp => kvp.Value >= minCount).Select(kvp => kvp.Key).ToHashSet();
        result.Select(p => p.Kmer).Should().BeEquivalentTo(expected,
            "completeness: the returned set equals exactly the count-map entries with count >= minCount");

        // The concrete worked answer: only ACGT qualifies, others excluded.
        result.Should().ContainSingle("only 'ACGT' occurs >= 2 times in 'ACGTACGT' (k = 4)");
        result[0].Kmer.Should().Be("ACGT");
        result[0].Count.Should().Be(2, "'ACGT' starts at positions 0 and 4");
        result.Select(p => p.Kmer).Should().NotContain(new[] { "CGTA", "GTAC", "TACG" },
            "the count-1 4-mers are below the threshold and must be excluded");

        // Documented descending order.
        result.Select(p => p.Count).Should().BeInDescendingOrder(
            "the contract orders qualifying k-mers by count descending");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// satisfy soundness + completeness + descending order for several k and minCount
    /// values — the returned set must always equal exactly {(w, c) : c = Count(w),
    /// c ≥ minCount} cross-checked against an independent CountKmers map, regardless of
    /// the random content. [CancelAfter] guards against any hang on the largest k scanned.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindKmersWithMinCount_RandomSequence_MatchesCountMapFilterForEveryThreshold()
    {
        const int length = 2000;
        string seq = RandomDna(length, seed: 34_001);

        foreach (int k in new[] { 1, 3, 5, 8 })
        {
            var counts = KmerAnalyzer.CountKmers(seq, k);
            int maxCount = counts.Values.Max();

            foreach (int minCount in new[] { 0, 1, 2, maxCount, maxCount + 1 })
            {
                var result = KmerAnalyzer.FindKmersWithMinCount(seq, k, minCount).ToList();

                var expected = counts.Where(kvp => kvp.Value >= minCount)
                    .Select(kvp => kvp.Key).ToHashSet();
                result.Select(p => p.Kmer).Should().BeEquivalentTo(expected,
                    $"at k = {k}, minCount = {minCount} the result equals exactly the count-map entries with count >= minCount");
                result.Should().OnlyContain(p => p.Count == counts[p.Kmer],
                    $"at k = {k}, minCount = {minCount} every reported count matches the independent CountKmers value");
                result.Select(p => p.Count).Should().BeInDescendingOrder(
                    $"at k = {k}, minCount = {minCount} the qualifying k-mers are ordered by count descending");
            }

            // A threshold strictly above the max achievable count is always empty.
            KmerAnalyzer.FindKmersWithMinCount(seq, k, maxCount + 1).ToList().Should().BeEmpty(
                $"at k = {k} no k-mer reaches count {maxCount + 1} (max observed is {maxCount}); the result is empty");
        }
    }

    #endregion

    #endregion
}
