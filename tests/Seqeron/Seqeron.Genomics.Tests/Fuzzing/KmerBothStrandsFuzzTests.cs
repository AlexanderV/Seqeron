using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the K-mer area — BOTH-STRANDS / strand-aware k-mer counting
/// (KMER-BOTH-001).
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
/// Unit: KMER-BOTH-001 — both-strands k-mer counting
/// Checklist: docs/checklists/03_FUZZING.md, row 157.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a PALINDROMIC k-mer (a k-mer equal to its own reverse
///          complement), k &gt; sequence length, and the empty sequence. We also pin
///          the off-by-one neighbours k = L + 1 and k = L, plus the k ≤ 0 floor.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The both-strands contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Double-stranded DNA carries information on both strands, so a strand-aware k-mer
/// profile sums each k-mer's forward count with that of its complementary reading.
/// Counting a k-mer w on the reverse-complement strand (read 5'→3') equals counting
/// its reverse complement RC(w) on the forward strand, so the method yields:
///   count[w] = forward[w] + forward[RC(w)]
/// — the "balance" operation of kPAL ("adding the values of each k-mer to its
/// reverse complement", Anvar et al., 2014), reflecting the generalized second
/// Chargaff rule / inversion symmetry (Shporer et al., 2016). The API entry under
/// test is
///   KmerAnalyzer.CountKmersBothStrands(string sequence, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 476–490),
/// with the DnaSequence overload delegating to it (lines 501–502). The method counts
/// the forward sequence and the reverse-complement sequence INDEPENDENTLY via
/// CountKmers, then sums the two count maps key-by-key (KmerAnalyzer.cs lines
/// 478–489; K-mer_Counting.md §5.2: "counts forward and reverse-complement sequences
/// independently before summing counts").
///
/// IMPORTANT — this is NOT canonical counting. Canonical k-mer counting (Marçais &amp;
/// Kingsford, 2011) COLLAPSES the pair {w, RC(w)} onto a single lexicographically-
/// smaller key. This method does the opposite: it RETAINS a key for every observed
/// k-mer; w and RC(w) carry EQUAL counts and remain separate dictionary entries
/// unless they are the same string (K-mer_Counting.md §5.3 "Intentionally simplified:
/// Both-strand counting sums forward and reverse-complement counts rather than
/// canonicalizing… forward and reverse-complement words remain separate dictionary
/// entries unless they are identical strings"; §5.3 "Not implemented: Canonical
/// reverse-complement collapsing"). The tests assert the DOCUMENTED both-strands
/// contract, not a canonical-collapse rule.
///
/// THE KEY INVARIANT (KmerAnalyzer.cs lines 466): the total over ALL k-mer counts is
/// exactly 2·(L − k + 1) — the forward strand contributes L − k + 1 windows and the
/// reverse-complement strand (same length) contributes another L − k + 1. Every
/// positive-result test below pins this total-instance invariant; it is the single
/// load-bearing correctness check that distinguishes a correct strand-aware tally
/// from a miscount or a wrongly-collapsed canonical count.
///
/// THE PALINDROME RULE (the headline BE target). A reverse-complement PALINDROME is a
/// k-mer with w == RC(w) (e.g. "AT", "GATC", "ACGT" — necessarily even-length over
/// the unambiguous A/C/G/T alphabet). For such a k-mer, forward[w] and forward[RC(w)]
/// are the SAME forward entry, so the documented sum count[w] = forward[w] +
/// forward[RC(w)] equals 2·forward[w]: a palindrome legitimately gets DOUBLE its
/// forward count BECAUSE it genuinely occurs on both strands at every forward
/// position. This is the both-strands semantics — distinct from canonical counting,
/// where a palindrome would be its own canonical representative and carry only its
/// forward count. The tests pin count[palindrome] == 2·forward[palindrome] exactly,
/// and confirm a palindrome does NOT spawn a second separate revcomp key (it IS its
/// own revcomp), so the total-instance invariant 2·(L − k + 1) still holds.
///
/// Documented parameter contract (KmerAnalyzer.cs lines 468–475; K-mer_Counting.md
/// §3.1, §3.3, §6.1):
///   • k ≤ 0 with NON-EMPTY input → ArgumentOutOfRangeException (nameof(k)),
///     surfaced unchanged from the underlying CountKmers (KmerAnalyzer.cs lines
///     25–26): a 0-length k-mer is meaningless on either strand.
///   • k &gt; sequence.Length → EMPTY dictionary (no windows fit on EITHER strand:
///     the reverse complement has the same length, so it has no windows either);
///     never an out-of-range Substring, never a crash (K-mer_Counting.md §6.1).
///   • empty / null sequence → EMPTY dictionary, NOT an exception. The forward
///     CountKmers short-circuits empty/null to the empty map before k is validated
///     (KmerAnalyzer.cs lines 22–23), and the reverse complement of null/empty is
///     itself empty (DnaSequence.GetReverseComplementString returns the input for
///     null/empty, lines 151–152; the method coalesces null to string.Empty at line
///     479) → the combined map is empty. Empty input wins even with a degenerate k.
/// The implementation uppercases before keying (CountKmers, INV-03), and the reverse
/// complement is built per-base via GetComplementBase (case-insensitive), so the
/// whole surface is case-insensitive. These tests exercise only the BE targets of
/// THIS row.
///
/// The three checklist targets map to these documented behaviours:
///   • palindromic k-mer → count[w] = 2·forward[w] (occurs on both strands); the
///     palindrome is its own revcomp, so no extra key, and the 2·(L−k+1) total holds.
///   • k &gt; len           → empty dictionary; no windows on either strand, no crash.
///   • empty seq          → empty dictionary; revcomp of empty is empty, no crash.
/// A positive-sanity test pins the both-strands rule on a known sequence with BOTH a
/// palindromic k-mer (doubled) and a non-palindromic k-mer (its count includes its
/// revcomp's forward occurrences), cross-checked against the documented
/// count[w] = forward[w] + forward[RC(w)] formula and the 2·(L−k+1) total.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerBothStrandsFuzzTests
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
    /// Asserts the result is a WELL-FORMED both-strands count map for the given
    /// sequence and k: every count is a positive integer, every key is a length-k
    /// uppercased k-mer over the observed forward∪revcomp alphabet, each key's count
    /// equals the DOCUMENTED forward[w] + forward[RC(w)], and the grand total equals
    /// 2·(L − k + 1). This is the load-bearing structural oracle reused across the
    /// fuzz cases (KmerAnalyzer.cs lines 459, 466; K-mer_Counting.md §5.2/§5.3).
    /// </summary>
    private static void AssertWellFormedBothStrands(
        IReadOnlyDictionary<string, int> result, string sequence, int k)
    {
        // Independent oracles: forward counts and the forward counts of the revcomp seq.
        var forward = KmerAnalyzer.CountKmers(sequence, k);
        var revComp = KmerAnalyzer.CountKmers(
            DnaSequence.GetReverseComplementString(sequence ?? string.Empty), k);

        // Every count is a strictly positive integer; no negative/zero phantom entries.
        result.Values.Should().OnlyContain(v => v > 0,
            "every observed k-mer has a strictly positive both-strands count");

        // Every key is a length-k uppercased k-mer.
        result.Keys.Should().OnlyContain(key => key.Length == k,
            "every emitted key is a length-k window — no shorter/longer fragments");
        result.Keys.Should().OnlyContain(key => key.ToUpperInvariant() == key,
            "INV-03: keys are uppercased before being emitted");

        // The documented per-key rule: count[w] = forward[w] + forward[RC(w)],
        // realized here as forward[w] + revCompForward[w] (the two summed maps).
        foreach (var (kmer, count) in result)
        {
            int expected = forward.GetValueOrDefault(kmer, 0) + revComp.GetValueOrDefault(kmer, 0);
            count.Should().Be(expected,
                $"count['{kmer}'] must equal forward + reverse-complement strand occurrences (documented sum)");
        }

        // The total-instance invariant: 2·(L − k + 1) when k ≤ L, else 0.
        int windows = (sequence?.Length ?? 0) - k + 1;
        int expectedTotal = windows > 0 ? 2 * windows : 0;
        result.Values.Sum().Should().Be(expectedTotal,
            "the total over all k-mers is 2·(L − k + 1): each strand contributes one window per valid start position");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-BOTH-001 — both-strands k-mer counting : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-BOTH-001 — both-strands k-mer counting

    #region BE — Boundary: palindromic k-mer (the headline both-strands rule)

    /// <summary>
    /// BE / KEY: a reverse-complement PALINDROME (w == RC(w)) is the headline target.
    /// "AT" is its own reverse complement, so in "AT" (L = 2, k = 2) the forward count
    /// of "AT" is 1 and the documented sum count[AT] = forward[AT] + forward[RC(AT)] =
    /// forward[AT] + forward[AT] = 2. The palindrome legitimately gets DOUBLE its
    /// forward count because it genuinely occurs on both strands at the same position —
    /// this is the both-strands semantics, NOT a double-counting bug. Crucially, because
    /// "AT" IS its own revcomp, it does NOT spawn a second separate key, and the
    /// total-instance invariant 2·(L − k + 1) = 2·1 = 2 still holds exactly (the whole
    /// total sits on the single palindromic key). We pin the doubled count, the single
    /// key, and the total.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_PalindromeAtKEqualsLength_IsDoubledForwardCount()
    {
        const string seq = "AT"; // RC("AT") == "AT": a reverse-complement palindrome.
        const int k = 2;

        var result = KmerAnalyzer.CountKmersBothStrands(seq, k);

        result.Should().ContainSingle("'AT' is its own reverse complement, so it is the only key — no separate revcomp entry");
        result.Should().ContainKey("AT").WhoseValue.Should().Be(2,
            "count[AT] = forward[AT] + forward[RC(AT)] = 1 + 1 = 2: the palindrome occurs on BOTH strands");
        result.Values.Sum().Should().Be(2,
            "total = 2·(L − k + 1) = 2·1 = 2, all concentrated on the single palindromic key");

        AssertWellFormedBothStrands(result, seq, k);
    }

    /// <summary>
    /// BE: a longer palindrome inside a longer sequence still doubles per the rule.
    /// "GATC" is a reverse-complement palindrome (RC("GATC") == "GATC"). In a sequence
    /// where "GATC" appears once on the forward strand, its both-strands count is
    /// 2·forward = 2, while NO separate "GATC"-revcomp key is created (it is its own
    /// revcomp). We compare directly against forward counts to prove count[palindrome]
    /// == 2·forward[palindrome], and against the well-formed oracle for the rest.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_PalindromeInsideSequence_IsExactlyTwiceForward()
    {
        const string seq = "GGATCC"; // contains the palindrome "GATC" once at position 1.
        const int k = 4;

        var forward = KmerAnalyzer.CountKmers(seq, k);
        var both = KmerAnalyzer.CountKmersBothStrands(seq, k);

        forward.Should().ContainKey("GATC").WhoseValue.Should().Be(1, "'GATC' occurs once on the forward strand");
        DnaSequence.GetReverseComplementString("GATC").Should().Be("GATC", "'GATC' is a reverse-complement palindrome");

        both.Should().ContainKey("GATC").WhoseValue.Should().Be(2 * forward["GATC"],
            "count[GATC] = forward[GATC] + forward[RC(GATC)] = 2·forward[GATC] because GATC is its own revcomp");

        AssertWellFormedBothStrands(both, seq, k);
    }

    /// <summary>
    /// BE: an ALL-palindrome homopolymer-like run — every length-2 window of an
    /// alternating-palindrome sequence is a palindrome. "ATAT" has windows AT, TA, AT;
    /// both "AT" and "TA" are reverse-complement palindromes (RC("AT")=="AT",
    /// RC("TA")=="TA"). Each therefore doubles, and the total stays 2·(L − k + 1) = 2·3
    /// = 6. This pins that MULTIPLE palindromic keys each double independently without
    /// over- or under-counting and without spawning spurious revcomp keys.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_AllPalindromeWindows_EachDoublesIndependently()
    {
        const string seq = "ATAT"; // 2-mer windows: AT, TA, AT — all palindromes.
        const int k = 2;

        var both = KmerAnalyzer.CountKmersBothStrands(seq, k);

        DnaSequence.GetReverseComplementString("AT").Should().Be("AT");
        DnaSequence.GetReverseComplementString("TA").Should().Be("TA");

        both.Should().HaveCount(2, "only the two palindromic 2-mers AT and TA appear — no separate revcomp keys");
        both["AT"].Should().Be(4, "forward[AT] = 2 (positions 0, 2), doubled → 4");
        both["TA"].Should().Be(2, "forward[TA] = 1 (position 1), doubled → 2");
        both.Values.Sum().Should().Be(6, "total = 2·(L − k + 1) = 2·3 = 6");

        AssertWellFormedBothStrands(both, seq, k);
    }

    #endregion

    #region BE — Boundary: k > sequence length

    /// <summary>
    /// BE: k far larger than the sequence length must NOT crash. No length-k window fits
    /// on the forward strand, and the reverse complement has the SAME length so no window
    /// fits there either — both underlying CountKmers calls return the empty map, so the
    /// combined result is empty (KmerAnalyzer.cs lines 28–29 via both branches;
    /// K-mer_Counting.md §6.1). We pin no-throw AND emptiness at a far-oversized k and at
    /// the exact off-by-one boundary k = L + 1, so an oversized k can never index past the
    /// end of either strand nor invent a k-mer longer than the sequence.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_KGreaterThanSequenceLength_IsEmptyAndDoesNotThrow()
    {
        var act = () => KmerAnalyzer.CountKmersBothStrands("ACGT", 1000);
        act.Should().NotThrow(
            "k > L makes the window count L − k + 1 negative on BOTH strands; neither loop runs, so nothing is indexed past the end");

        KmerAnalyzer.CountKmersBothStrands("ACGT", 1000).Should().BeEmpty(
            "no length-1000 window fits a 4-base sequence (or its 4-base revcomp); the result is empty, not a crash");

        // k = L + 1 is the exact off-by-one boundary above the sequence length.
        KmerAnalyzer.CountKmersBothStrands("ACGT", 5).Should().BeEmpty(
            "k = L + 1 is one past the last fitting window on either strand; still empty, never an out-of-range Substring");
    }

    /// <summary>
    /// BE: k = sequence length is the upper boundary where exactly ONE window fits per
    /// strand. For the non-palindromic whole sequence "ACGTAA" (RC = "TTACGT"), the
    /// forward whole-sequence k-mer and the revcomp whole-sequence k-mer are DIFFERENT
    /// strings, so the result has TWO keys each with count 1, and the total is
    /// 2·(L − k + 1) = 2·1 = 2. This pins the off-by-one upper edge (one window per
    /// strand) and the both-strands total at that edge.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_KEqualsLengthNonPalindrome_IsForwardAndRevCompKmers()
    {
        const string seq = "ACGTAA"; // RC = "TTACGT", different from the forward whole sequence.
        var rc = DnaSequence.GetReverseComplementString(seq);
        rc.Should().Be("TTACGT", "the whole sequence is not a reverse-complement palindrome");

        var both = KmerAnalyzer.CountKmersBothStrands(seq, seq.Length);

        both.Should().HaveCount(2, "k = L admits one window per strand; the two whole-strand k-mers differ, so two keys");
        both[seq].Should().Be(1, "the forward whole-sequence k-mer occurs once");
        both[rc].Should().Be(1, "the reverse-complement whole-sequence k-mer occurs once");
        both.Values.Sum().Should().Be(2, "total = 2·(L − k + 1) = 2·1 = 2");

        AssertWellFormedBothStrands(both, seq, seq.Length);
    }

    #endregion

    #region BE — Boundary: k = 0 / negative k

    /// <summary>
    /// BE: k = 0 is the degenerate floor — a 0-length k-mer is meaningless on either
    /// strand. CountKmersBothStrands counts via CountKmers, which rejects k ≤ 0 on
    /// non-empty input with ArgumentOutOfRangeException (KmerAnalyzer.cs lines 25–26);
    /// the forward count is evaluated first, so that rejection surfaces unchanged. We pin
    /// that k = 0 throws and carries the documented "k" parameter name, so a 0-length
    /// k-mer can never reach the strand-summing logic.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_KZero_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.CountKmersBothStrands("ACGTACGT", 0);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; both-strands counting surfaces the underlying k <= 0 rejection on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: a negative k is below the floor too and must be rejected the same way —
    /// pinning that the rejection boundary is exactly k ≤ 0, not merely k == 0, so a
    /// negative length can never slip into either strand's window loop.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_NegativeK_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.CountKmersBothStrands("ACGTACGT", -3);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0 on non-empty input")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The forward CountKmers
    /// short-circuits empty/null input to the empty dictionary BEFORE k is validated
    /// (KmerAnalyzer.cs lines 22–23), and the reverse complement of null/empty is itself
    /// empty (the method coalesces null to string.Empty at line 479;
    /// GetReverseComplementString returns the input for null/empty, DnaSequence.cs lines
    /// 151–152) — so the combined map is empty with no exception, even when k is itself
    /// degenerate (k = 0). We pin that empty, null, and the empty-DnaSequence surface all
    /// return the empty dictionary, and that empty input wins over a degenerate k rather
    /// than throwing (K-mer_Counting.md §6.1).
    /// </summary>
    [Test]
    public void CountKmersBothStrands_EmptyOrNullSequence_IsEmptyAndDoesNotThrow()
    {
        var emptyAct = () => KmerAnalyzer.CountKmersBothStrands(string.Empty, 3);
        var nullAct = () => KmerAnalyzer.CountKmersBothStrands((string)null!, 3);
        var emptyDnaAct = () => KmerAnalyzer.CountKmersBothStrands(new DnaSequence(string.Empty), 3);
        var emptyDegenerateKAct = () => KmerAnalyzer.CountKmersBothStrands(string.Empty, 0);

        emptyAct.Should().NotThrow("an empty sequence has no windows on either strand; the result is the empty dictionary");
        nullAct.Should().NotThrow("null input is treated as empty (forward guard + null-coalesced revcomp), not as an error");
        emptyDnaAct.Should().NotThrow("the DnaSequence overload delegates to the string path's empty handling");
        emptyDegenerateKAct.Should().NotThrow(
            "the forward empty/null guard runs BEFORE k is validated, so empty input wins even with a degenerate k = 0");

        KmerAnalyzer.CountKmersBothStrands(string.Empty, 3).Should().BeEmpty();
        KmerAnalyzer.CountKmersBothStrands((string)null!, 3).Should().BeEmpty();
        KmerAnalyzer.CountKmersBothStrands(new DnaSequence(string.Empty), 3).Should().BeEmpty();
        KmerAnalyzer.CountKmersBothStrands(string.Empty, 0).Should().BeEmpty(
            "empty input short-circuits to an empty dictionary before the k <= 0 check can throw");
    }

    #endregion

    #region Positive sanity — the both-strands rule on known sequences

    /// <summary>
    /// Positive sanity / KEY: a known sequence with BOTH a palindromic and a
    /// non-palindromic k-mer, hand-computed from count[w] = forward[w] + forward[RC(w)].
    /// "AATT" (L = 4), k = 2. Forward 2-mers: AA, AT, TT.
    ///   forward = { AA:1, AT:1, TT:1 }.
    ///   RC("AATT") = "AATT" (the whole sequence is a palindrome), so revcomp 2-mers are
    ///   the SAME multiset { AA:1, AT:1, TT:1 }.
    /// Documented sums:
    ///   • "AT" is a PALINDROME (RC("AT")=="AT"): count[AT] = forward[AT] + forward[AT]
    ///     = 1 + 1 = 2 — doubled.
    ///   • "AA" is NON-palindromic (RC("AA")=="TT"): count[AA] = forward[AA] +
    ///     forward[RC(AA)=TT] = 1 + 1 = 2 — its count INCLUDES the forward occurrence of
    ///     its reverse complement "TT". Likewise count[TT] = forward[TT] + forward[AA]
    ///     = 1 + 1 = 2. So AA and TT remain SEPARATE keys carrying EQUAL counts (the
    ///     inversion-symmetry property), NOT collapsed onto one canonical key.
    /// Total = 2·(L − k + 1) = 2·3 = 6. We pin every count, the palindrome-vs-non-
    /// palindrome distinction, the separate-but-equal AA/TT keys, and the total.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_KnownPalindromeAndNonPalindrome_MatchHandComputedSums()
    {
        const string seq = "AATT";
        const int k = 2;

        var both = KmerAnalyzer.CountKmersBothStrands(seq, k);

        // Palindromic 2-mer: doubled forward count.
        DnaSequence.GetReverseComplementString("AT").Should().Be("AT", "'AT' is a reverse-complement palindrome");
        both["AT"].Should().Be(2, "count[AT] = forward[AT] + forward[RC(AT)=AT] = 1 + 1 = 2 (palindrome doubles)");

        // Non-palindromic pair AA / TT: separate keys, each count = forward + revcomp-of-itself.
        DnaSequence.GetReverseComplementString("AA").Should().Be("TT", "'AA' is NOT a palindrome; its revcomp is 'TT'");
        both["AA"].Should().Be(2, "count[AA] = forward[AA] + forward[RC(AA)=TT] = 1 + 1 = 2");
        both["TT"].Should().Be(2, "count[TT] = forward[TT] + forward[RC(TT)=AA] = 1 + 1 = 2 (equal to AA — inversion symmetry)");
        both.Should().ContainKeys("AA", "TT", "AT");
        both.Should().HaveCount(3, "the result retains a key per observed k-mer; AA and TT are SEPARATE (not canonicalized)");

        both.Values.Sum().Should().Be(6, "total = 2·(L − k + 1) = 2·3 = 6 (each strand contributes 3 windows)");

        AssertWellFormedBothStrands(both, seq, k);
    }

    /// <summary>
    /// Positive sanity: the inversion-symmetry property w / RC(w) carry EQUAL counts
    /// (Shporer et al., 2016; KmerAnalyzer.cs lines 461–466). For any sequence and any
    /// k, count[w] must equal count[RC(w)] whenever both are present — because
    /// count[w] = forward[w] + forward[RC(w)] and count[RC(w)] = forward[RC(w)] +
    /// forward[w] are the SAME sum. We pin this on a heterogeneous sequence.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_HeterogeneousSequence_KmerAndRevCompCarryEqualCounts()
    {
        const string seq = "ACGTACGAA";
        const int k = 3;

        var both = KmerAnalyzer.CountKmersBothStrands(seq, k);

        foreach (var kmer in both.Keys)
        {
            string rc = DnaSequence.GetReverseComplementString(kmer);
            both.Should().ContainKey(rc,
                $"both-strands counting retains a key for the reverse complement of every observed k-mer ('{kmer}' → '{rc}')");
            both[rc].Should().Be(both[kmer],
                $"inversion symmetry: count['{kmer}'] == count['{rc}'] (both equal forward['{kmer}'] + forward['{rc}'])");
        }

        AssertWellFormedBothStrands(both, seq, k);
    }

    /// <summary>
    /// Positive sanity: case-insensitivity. A lowercase sequence must produce the same
    /// uppercased keys and the same both-strands counts as its uppercase form, because
    /// the forward count uppercases (INV-03) and the reverse complement is built via the
    /// case-insensitive GetComplementBase (SequenceExtensions.cs lines 138–157).
    /// </summary>
    [Test]
    public void CountKmersBothStrands_LowercaseInput_IsUppercasedAndComplementInsensitive()
    {
        var lower = KmerAnalyzer.CountKmersBothStrands("acgtacgt", 3);
        var upper = KmerAnalyzer.CountKmersBothStrands("ACGTACGT", 3);

        lower.Should().BeEquivalentTo(upper,
            "case does not change the strand-aware counts: forward uppercases and complement is case-insensitive");
    }

    /// <summary>
    /// Positive sanity: the DnaSequence overload must agree with the string overload it
    /// delegates to (KmerAnalyzer.cs lines 501–502), so the strand-aware contract is
    /// identical on both surfaces.
    /// </summary>
    [Test]
    public void CountKmersBothStrands_DnaSequenceOverload_AgreesWithStringOverload()
    {
        const string seq = "GGATCCACGT";
        const int k = 4;

        var viaString = KmerAnalyzer.CountKmersBothStrands(seq, k);
        var viaDna = KmerAnalyzer.CountKmersBothStrands(new DnaSequence(seq), k);

        viaDna.Should().BeEquivalentTo(viaString,
            "the DnaSequence overload delegates to the string overload, so the counts must be identical");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// satisfy the KEY both-strands invariants for several k values — the total equals
    /// 2·(L − k + 1), every count equals the documented forward + revcomp sum, and
    /// w / RC(w) carry equal counts — regardless of the random content. [CancelAfter]
    /// guards against any hang on the largest k scanned.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void CountKmersBothStrands_RandomSequence_SatisfiesBothStrandsInvariantsForEveryK()
    {
        const int length = 2000;
        string seq = RandomDna(length, seed: 157_001);

        foreach (int k in new[] { 1, 2, 3, 5, 8, 13 })
        {
            var both = KmerAnalyzer.CountKmersBothStrands(seq, k);

            // Total-instance invariant + per-key documented sum + structural checks.
            AssertWellFormedBothStrands(both, seq, k);

            // Inversion symmetry across every observed key.
            foreach (var kmer in both.Keys)
            {
                string rc = DnaSequence.GetReverseComplementString(kmer);
                both.Should().ContainKey(rc, $"at k = {k} the revcomp of every key is also a key");
                both[rc].Should().Be(both[kmer],
                    $"at k = {k} inversion symmetry holds: count['{kmer}'] == count['{rc}']");
            }
        }
    }

    #endregion

    #endregion
}
