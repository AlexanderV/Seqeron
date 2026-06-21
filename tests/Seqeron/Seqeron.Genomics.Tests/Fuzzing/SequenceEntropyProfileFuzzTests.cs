using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area sliding-window Shannon-entropy-profile unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no infinite loop, no
/// state corruption, and no *unhandled* runtime exception (IndexOutOfRange,
/// NullReference, DivideByZero, OverflowException, ArgumentOutOfRange, …) and no
/// NaN / Infinity leaking into a profile value. Every input must result in EITHER
/// a well-defined, theory-correct profile OR a *documented, intentional* outcome
/// (here: an empty profile). A raw runtime exception, a hang, or a NaN on garbage
/// input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-ENTROPY-PROFILE-001 — Shannon entropy profile (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 232.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — windowSize &gt; length, windowSize = 0, empty /
///          null sequence, single char, plus windowSize = length, very long input,
///          and a randomized boundary sweep.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes; BE = 0, -1,
///   MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The entropy-profile contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateEntropyProfile(string, int windowSize = 50,
///            int stepSize = 1) → IEnumerable&lt;double&gt; (lazy / deferred).
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)
///
/// Documented behaviour (docs/algorithms/Statistics/Entropy_Profile.md,
/// Test Unit ID SEQ-ENTROPY-PROFILE-001):
///   • §2.2 / §4.1: for each window of width W slid by step over the sequence,
///     emit H = −Σᵢ pᵢ log₂ pᵢ (bits), pᵢ = symbol count / counted symbols. Terms
///     with pᵢ = 0 contribute 0. The per-window kernel case-folds (ToUpperInvariant)
///     and counts ONLY char.IsLetter symbols; there is NO T↔U normalization.
///   • §2.4 INV-01: every profile value H ≥ 0.
///   • §2.4 INV-02 / §6.1: every value H ≤ log₂ k (k = distinct symbols in the
///     window) ≤ 2 bits for the 4-letter DNA alphabet.
///   • §2.4 INV-03 / §6.1: a homopolymer window yields H = 0.0.
///   • §2.4 INV-04 / §6.1: a window with all symbols equally frequent yields
///     H = log₂ k (a uniform 4-base DNA window = 2.0 bits).
///   • §2.4 INV-05: number of windows = ⌊(n − W)/step⌋ + 1 when W ≤ n, else 0;
///     windows start at offsets 0, step, 2·step, … ≤ n − W.
///   • §3.3 / §6.1: null / empty sequence, OR windowSize &gt; length → EMPTY profile
///     (no exception). Indexing of window offsets is 0-based.
///   • §6.1: windowSize == length → exactly one window (single value).
///
/// The Shannon formula H = −Σ p·log₂ p is re-implemented independently in the
/// <see cref="Shannon"/> oracle below — expected values are derived from the spec,
/// never read off the implementation's own output.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceEntropyProfileFuzzTests
{
    #region SEQ-ENTROPY-PROFILE-001 — Helpers and independent Shannon oracle

    private const double Tolerance = 1e-9;

    /// <summary>The DNA alphabet — its max entropy is log₂ 4 = 2 bits (INV-02/§6.1).</summary>
    private const string Dna = "ACGT";

    /// <summary>The 2-bit DNA ceiling, log₂ 4 (Entropy_Profile.md §2.2 / INV-02).</summary>
    private static readonly double Log2OfFour = Math.Log2(4);

    /// <summary>
    /// Independent oracle for the per-window Shannon entropy kernel: re-derives
    /// H = −Σ pᵢ log₂ pᵢ from the SPEC (§2.2) — case-fold with ToUpperInvariant and
    /// count ONLY letters, exactly as §5.2 documents, with the 0·log 0 ≡ 0
    /// convention. Written from the formula, NOT copied from the implementation.
    /// </summary>
    private static double Shannon(string window)
    {
        var counts = new Dictionary<char, int>();
        int total = 0;
        foreach (char ch in window.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
            {
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
                total++;
            }
        }
        if (total == 0)
            return 0.0;

        double h = 0.0;
        foreach (int c in counts.Values)
        {
            double p = (double)c / total;
            h -= p * Math.Log2(p); // p > 0 here, so no 0·log 0 term
        }
        return h;
    }

    /// <summary>
    /// Independent oracle for the WHOLE profile: the documented window-offset rule
    /// (§4.1 / INV-05). Yields one Shannon value per offset i = 0, step, 2·step, …
    /// while i ≤ length − W, when W ≤ length; otherwise nothing. Derived from the
    /// spec, not from the implementation's loop.
    /// </summary>
    private static List<double> OracleProfile(string seq, int windowSize, int stepSize)
    {
        var result = new List<double>();
        if (string.IsNullOrEmpty(seq) || windowSize > seq.Length)
            return result;
        for (int i = 0; i <= seq.Length - windowSize; i += stepSize)
            result.Add(Shannon(seq.Substring(i, windowSize)));
        return result;
    }

    /// <summary>
    /// The documented window count: 0 for the null/empty-sequence guard (§3.3, which
    /// precedes the window logic — so even W = 0 on an empty sequence yields nothing),
    /// 0 when W &gt; n, else ⌊(n − W)/step⌋ + 1 (Entropy_Profile.md §2.4 INV-05).
    /// Computed here from the spec directly, independent of the implementation.
    /// </summary>
    private static int ExpectedWindowCount(int n, int windowSize, int stepSize)
    {
        if (n == 0)
            return 0; // empty/null guard wins before any window is produced (§3.3)
        return windowSize > n ? 0 : (n - windowSize) / stepSize + 1;
    }

    /// <summary>
    /// Universal well-formedness contract for ANY DNA-alphabet input: every value is
    /// finite (never NaN/Infinity), ≥ 0 (INV-01) and ≤ log₂ 4 = 2 bits (INV-02/§6.1
    /// — the DNA ceiling). Use the DNA ceiling only when the window symbols are over
    /// {A,C,G,T}; the general bound log₂ k is asserted via the oracle elsewhere.
    /// </summary>
    private static void AssertDnaWellFormed(IReadOnlyList<double> profile)
    {
        profile.Should().NotBeNull("the method must always return a sequence, never null");
        foreach (double h in profile)
        {
            double.IsFinite(h).Should().BeTrue($"entropy {h} must be finite — never NaN/Infinity");
            h.Should().BeGreaterThanOrEqualTo(-Tolerance, "INV-01: every profile value H ≥ 0");
            h.Should().BeLessThanOrEqualTo(Log2OfFour + Tolerance,
                "INV-02/§6.1: a DNA window's entropy never exceeds log₂ 4 = 2 bits");
        }
    }

    /// <summary>Random DNA string over {A,C,G,T} of the given length.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = Dna[rng.Next(Dna.Length)];
        return new string(chars);
    }

    /// <summary>Random BMP code points (control chars, null byte, lone surrogate
    /// halves, unicode letters/digits) — random-byte fuzz fodder for the profiler.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-ENTROPY-PROFILE-001 — entropy profile : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-checkable worked example (§7.1)

    /// <summary>
    /// Positive baseline (not a boundary): the doc's worked example must be reproduced
    /// EXACTLY. "AAATGC" with window 4, step 1 has ⌊(6−4)/1⌋+1 = 3 windows:
    ///   AAAT → −(3/4·log₂(3/4) + 1/4·log₂(1/4)) = 0.8112781244591328
    ///   AATG → −(1/2·log₂(1/2) + 2·(1/4·log₂(1/4))) = 1.5
    ///   ATGC → log₂ 4 = 2.0
    /// The three Shannon values are computed here independently from the formula
    /// (literals re-derived by hand) and re-checked against the kernel oracle — this
    /// pins the BUSINESS contract (Shannon H in bits per window), not mere
    /// absence of crashes. — Entropy_Profile.md §7.1.
    /// </summary>
    [Test]
    public void EntropyProfile_WorkedExample_MatchesHandComputedBits()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AAATGC", windowSize: 4, stepSize: 1)
            .ToArray();

        profile.Should().HaveCount(3, "⌊(6−4)/1⌋+1 = 3 windows (INV-05)");

        // Hand-derived literals from §7.1.
        profile[0].Should().BeApproximately(0.8112781244591328, Tolerance, "AAAT → 0.8112781… bits");
        profile[1].Should().BeApproximately(1.5, Tolerance, "AATG → 1.5 bits");
        profile[2].Should().BeApproximately(2.0, Tolerance, "ATGC → log₂ 4 = 2.0 bits");

        // Cross-check the hand literals against the independent Shannon oracle.
        profile[0].Should().BeApproximately(Shannon("AAAT"), Tolerance);
        profile[1].Should().BeApproximately(Shannon("AATG"), Tolerance);
        profile[2].Should().BeApproximately(Shannon("ATGC"), Tolerance);

        AssertDnaWellFormed(profile);
    }

    /// <summary>
    /// Positive: INV-03 (homopolymer → 0) and INV-04 (uniform 4-base → log₂ 4 = 2)
    /// as exact numeric anchors. A homopolymer "AAAAAAAA" (window 4, step 1) yields
    /// five windows ALL equal to 0.0; a tiled "ACGTACGT" yields windows that are each
    /// a permutation of {A,C,G,T} ⇒ all 2.0 bits. — Entropy_Profile.md §6.1 / INV-03 / INV-04.
    /// </summary>
    [Test]
    public void EntropyProfile_HomopolymerIsZero_UniformIsTwoBits()
    {
        var homo = SequenceStatistics.CalculateEntropyProfile("AAAAAAAA", windowSize: 4, stepSize: 1)
            .ToArray();
        homo.Should().HaveCount(5, "⌊(8−4)/1⌋+1 = 5 windows");
        homo.Should().OnlyContain(h => h == 0.0, "INV-03: a homopolymer window is 0.0 bits");

        var uniform = SequenceStatistics.CalculateEntropyProfile("ACGTACGT", windowSize: 4, stepSize: 1)
            .ToArray();
        uniform.Should().HaveCount(5);
        uniform.Should().OnlyContain(h => Math.Abs(h - 2.0) <= Tolerance,
            "INV-04: every window is a permutation of {A,C,G,T} ⇒ log₂ 4 = 2.0 bits");

        AssertDnaWellFormed(homo);
        AssertDnaWellFormed(uniform);
    }

    /// <summary>
    /// Positive: the step parameter must select offsets 0, step, 2·step, … and the
    /// window count must follow INV-05 exactly. "AAATGC" with window 4, step 2 picks
    /// offsets {0, 2}: ⌊(6−4)/2⌋+1 = 2 windows → AAAT (0.8112781…) and TGC?… here
    /// offset 2 = "ATGC" (2.0). Confirms step is honoured (not always 1).
    /// — Entropy_Profile.md §4.1 / INV-05.
    /// </summary>
    [Test]
    public void EntropyProfile_StepGreaterThanOne_SelectsStridedOffsets()
    {
        var profile = SequenceStatistics.CalculateEntropyProfile("AAATGC", windowSize: 4, stepSize: 2)
            .ToArray();

        profile.Should().HaveCount(2, "offsets {0, 2}: ⌊(6−4)/2⌋+1 = 2 windows");
        profile[0].Should().BeApproximately(Shannon("AAAT"), Tolerance, "offset 0 → AAAT");
        profile[1].Should().BeApproximately(Shannon("ATGC"), Tolerance, "offset 2 → ATGC = 2.0 bits");
        AssertDnaWellFormed(profile);
    }

    #endregion

    #region BE — Boundary: windowSize > length (no full window ⇒ empty profile)

    /// <summary>
    /// BE: when the window is wider than the sequence, no full window exists, so the
    /// documented result is the EMPTY profile (no exception, no negative-length
    /// Substring). Covered just over the edge (W = n+1) and far past it (W = 1000),
    /// across short DNA inputs. — Entropy_Profile.md §3.3 / §6.1 / INV-05 (W &gt; n ⇒ 0).
    /// </summary>
    [TestCase("ACGT", 5)]
    [TestCase("ACGT", 1000)]
    [TestCase("A", 2)]
    [TestCase("ACGTACGTAC", 11)]
    public void EntropyProfile_WindowWiderThanSequence_IsEmpty(string seq, int windowSize)
    {
        var act = () => SequenceStatistics.CalculateEntropyProfile(seq, windowSize, stepSize: 1).ToArray();

        act.Should().NotThrow("windowSize > length is a defined boundary, not an error");
        act().Should().BeEmpty("no full window exists when W > length (INV-05)");
    }

    #endregion

    #region BE — Boundary: windowSize == length (exactly one window)

    /// <summary>
    /// BE: when the window exactly equals the sequence length, INV-05 gives
    /// ⌊(n−n)/step⌋+1 = 1 — a single window spanning the whole sequence, whose value
    /// is the Shannon entropy of the entire input. "ACGT" (window 4) ⇒ one value,
    /// 2.0 bits; "AAAA" ⇒ one value, 0.0 bits. — Entropy_Profile.md §6.1 (W == length).
    /// </summary>
    [Test]
    public void EntropyProfile_WindowEqualsLength_IsSingleWholeSequenceValue()
    {
        var uniform = SequenceStatistics.CalculateEntropyProfile("ACGT", windowSize: 4).ToArray();
        uniform.Should().HaveCount(1, "INV-05: ⌊(4−4)/1⌋+1 = 1 window");
        uniform[0].Should().BeApproximately(2.0, Tolerance, "whole 'ACGT' = log₂ 4 = 2.0 bits");

        var homo = SequenceStatistics.CalculateEntropyProfile("AAAA", windowSize: 4).ToArray();
        homo.Should().HaveCount(1);
        homo[0].Should().BeApproximately(0.0, Tolerance, "whole 'AAAA' homopolymer = 0.0 bits");

        AssertDnaWellFormed(uniform);
        AssertDnaWellFormed(homo);
    }

    #endregion

    #region BE — Boundary: windowSize == 0 (degenerate parameter, defined result)

    /// <summary>
    /// BE: windowSize = 0 is the degenerate zero-window-width boundary. It is NOT
    /// guarded by the W &gt; length check (0 ≤ length), so by INV-05 the profile has
    /// ⌊(n−0)/step⌋+1 = n+1 windows, each a length-0 substring whose Shannon entropy
    /// (no counted symbols, total = 0) is the documented 0.0. The key fuzz guarantee:
    /// NO divide-by-zero (the kernel returns 0 on a zero total) and NO infinite loop
    /// (step ≥ 1 advances the offset to termination at i = n). We pin the exact
    /// shape: n+1 zero values. — Entropy_Profile.md §2.2 (0·log 0 ≡ 0; total = 0 ⇒ 0)
    /// and §4.1 (offsets while i ≤ n − W = n) / INV-05.
    /// </summary>
    [TestCase("ACGT", 5)]   // n = 4 → 5 windows
    [TestCase("A", 2)]      // n = 1 → 2 windows
    [TestCase("ACGTACGTAC", 11)]
    [CancelAfter(5000)]
    public void EntropyProfile_WindowZero_YieldsZeroBitsPerOffset_NoCrashNoHang(string seq, int expectedCount)
    {
        IReadOnlyList<double> profile = null!;
        var act = () => profile = SequenceStatistics.CalculateEntropyProfile(seq, windowSize: 0, stepSize: 1)
            .ToArray();

        act.Should().NotThrow("windowSize = 0 must not divide by zero or throw");
        profile.Should().HaveCount(expectedCount,
            "INV-05: ⌊(n−0)/1⌋+1 = n+1 zero-width windows");
        profile.Should().OnlyContain(h => h == 0.0,
            "a zero-width window has no counted symbols ⇒ entropy 0.0 (0·log 0 ≡ 0)");
        AssertDnaWellFormed(profile);
    }

    /// <summary>
    /// BE: windowSize = 0 on the EMPTY string is guarded by the null/empty
    /// short-circuit (§3.3) ⇒ empty profile, even though 0 ≤ 0. Confirms the empty
    /// guard precedes the zero-window logic. — Entropy_Profile.md §3.3.
    /// </summary>
    [Test]
    public void EntropyProfile_WindowZero_OnEmptySequence_IsEmpty()
    {
        var act = () => SequenceStatistics.CalculateEntropyProfile(string.Empty, windowSize: 0).ToArray();

        act.Should().NotThrow();
        act().Should().BeEmpty("empty sequence is guarded before any window is produced (§3.3)");
    }

    #endregion

    #region BE — Boundary: empty / null sequence (defined ⇒ empty profile)

    /// <summary>
    /// BE: the empty string is the lower size boundary — null/empty yields the EMPTY
    /// profile (no exception). Asserted across several window sizes (incl. the default
    /// 50). — Entropy_Profile.md §3.3 / §6.1 (null/empty → empty profile).
    /// </summary>
    [TestCase(1)]
    [TestCase(4)]
    [TestCase(50)]
    public void EntropyProfile_EmptySequence_IsEmpty(int windowSize)
    {
        var act = () => SequenceStatistics.CalculateEntropyProfile(string.Empty, windowSize).ToArray();

        act.Should().NotThrow("the empty string is a defined boundary, not an error");
        act().Should().BeEmpty("no symbols ⇒ no window ⇒ empty profile");
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — empty
    /// profile, no NullReferenceException. — Entropy_Profile.md §3.3 (null → empty).
    /// </summary>
    [Test]
    public void EntropyProfile_NullSequence_IsEmpty_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateEntropyProfile(null!, windowSize: 4).ToArray();

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().BeEmpty();
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a single-character sequence. With the default window (50 &gt; 1) no full
    /// window exists ⇒ empty profile. With window 1 there is exactly one window whose
    /// single symbol gives p = 1, H = −1·log₂ 1 = 0.0 (homopolymer of length 1,
    /// INV-03). With window 1 and a 2-char homopolymer "AA" there are two windows,
    /// both 0.0. — Entropy_Profile.md §6.1 (W &gt; length → empty; INV-03).
    /// </summary>
    [Test]
    public void EntropyProfile_SingleChar_DefaultWindowEmpty_Window1IsZero()
    {
        // Default window 50 > length 1 → empty.
        SequenceStatistics.CalculateEntropyProfile("A").ToArray()
            .Should().BeEmpty("default window 50 > single-char length ⇒ empty (INV-05)");

        // Window 1 over a single char → one window, 0.0 bits.
        var single = SequenceStatistics.CalculateEntropyProfile("A", windowSize: 1).ToArray();
        single.Should().HaveCount(1, "⌊(1−1)/1⌋+1 = 1 window");
        single[0].Should().BeApproximately(0.0, Tolerance, "one symbol ⇒ p=1 ⇒ H=0 (INV-03)");

        // Window 1 over "AA" → two windows, both 0.0.
        var two = SequenceStatistics.CalculateEntropyProfile("AA", windowSize: 1).ToArray();
        two.Should().HaveCount(2);
        two.Should().OnlyContain(h => h == 0.0, "each length-1 window is a single symbol ⇒ 0 bits");

        AssertDnaWellFormed(single);
        AssertDnaWellFormed(two);
    }

    #endregion

    #region BE — non-letter symbols ignored, no T↔U normalization

    /// <summary>
    /// BE/MC: the per-window kernel counts ONLY char.IsLetter symbols (§5.2). A window
    /// of pure non-letters (digits, punctuation, the null byte) has total = 0 ⇒ 0.0
    /// bits, with no divide-by-zero. And there is NO T↔U normalization (§3.3): a window
    /// "ATUU" has THREE distinct symbols {A,T,U} (counts 1,1,2) ⇒ H = 1.5 bits, NOT the
    /// 2 distinct symbols that a T≡U merge would give. — Entropy_Profile.md §3.3 / §5.2.
    /// </summary>
    [Test]
    public void EntropyProfile_NonLettersIgnored_AndNoTUEquivalence()
    {
        // Pure non-letter window → 0.0, no crash.
        var junk = SequenceStatistics.CalculateEntropyProfile("12-+", windowSize: 4).ToArray();
        junk.Should().HaveCount(1);
        junk[0].Should().BeApproximately(0.0, Tolerance, "no counted letters ⇒ total 0 ⇒ 0 bits");

        // U and T are distinct symbols (no normalization): ATUU = {A:1, T:1, U:2}.
        var withU = SequenceStatistics.CalculateEntropyProfile("ATUU", windowSize: 4).ToArray();
        withU.Should().HaveCount(1);
        withU[0].Should().BeApproximately(Shannon("ATUU"), Tolerance,
            "U is its own symbol; H over {A,T,U} = 1.5 bits (no T↔U merge)");
        withU[0].Should().BeApproximately(1.5, Tolerance);
    }

    #endregion

    #region BE — Boundary: very long sequence (O(n·W), no overflow / hang)

    /// <summary>
    /// BE: a long sequence must be processed without overflow, hang or NaN. A
    /// 200 000-base "ACGT" tiling with window 4, step 4 gives 50 000 windows, each a
    /// permutation of {A,C,G,T} ⇒ all 2.0 bits. Bounds runtime under [CancelAfter] and
    /// pins INV-05 at scale. — Entropy_Profile.md §4.3 (O(n·W)) / INV-04 / INV-05.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void EntropyProfile_VeryLongTiledSequence_AllTwoBits()
    {
        const int tiles = 50_000;
        string seq = string.Concat(Enumerable.Repeat("ACGT", tiles)); // 200 000 bases

        var profile = SequenceStatistics.CalculateEntropyProfile(seq, windowSize: 4, stepSize: 4)
            .ToArray();

        profile.Should().HaveCount(tiles, "non-overlapping ACGT tiles: 50 000 windows");
        profile.Should().OnlyContain(h => Math.Abs(h - 2.0) <= Tolerance,
            "every tiled window is a {A,C,G,T} permutation ⇒ 2.0 bits");
        AssertDnaWellFormed(profile);
    }

    #endregion

    #region BE / RB — Randomized boundary sweep: never crash, match oracle

    /// <summary>
    /// BE: a randomized sweep over the documented boundary space — random DNA length
    /// (incl. 0 and 1), random windowSize that straddles 0, &lt; len, == len and &gt; len,
    /// and random step ≥ 1. The profile must EXACTLY match the independent
    /// offset-rule + Shannon oracle (count = INV-05; values = §2.2) and be DNA-well-formed,
    /// with no crash, hang, NaN or Infinity. Locally seeded Random.
    /// — Entropy_Profile.md §2.2 / §4.1 / INV-01..INV-05.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void EntropyProfile_RandomDnaBoundarySweep_MatchesOracle()
    {
        var rng = new Random(2026_06_21);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 40);
            string seq = RandomDna(rng, len);
            // windowSize straddles every boundary: 0, < len, == len, > len.
            int windowSize = rng.Next(0, len + 4);
            int stepSize = rng.Next(1, 5);

            IReadOnlyList<double> profile = null!;
            var act = () => profile = SequenceStatistics
                .CalculateEntropyProfile(seq, windowSize, stepSize).ToArray();

            act.Should().NotThrow(
                $"DNA len {len}, window {windowSize}, step {stepSize} must never crash");

            var oracle = OracleProfile(seq, windowSize, stepSize);
            profile.Should().HaveCount(ExpectedWindowCount(len, windowSize, stepSize),
                $"INV-05 window count (len {len}, W {windowSize}, step {stepSize})");
            profile.Should().HaveCount(oracle.Count);
            for (int k = 0; k < oracle.Count; k++)
                profile[k].Should().BeApproximately(oracle[k], Tolerance,
                    $"window {k} entropy (len {len}, W {windowSize}, step {stepSize})");

            AssertDnaWellFormed(profile);
        }
    }

    /// <summary>
    /// BE/RB: arbitrary BMP garbage (control chars, null byte, lone surrogate halves,
    /// unicode letters/digits/symbols) with random window/step must NEVER throw and
    /// must ALWAYS match the spec oracle (letters counted, non-letters ignored; §5.2).
    /// Each value must be finite and ≥ 0; the upper bound is the GENERAL log₂ k from
    /// the oracle (not the DNA 2-bit ceiling, since unicode windows may have k > 4).
    /// Core fuzz guarantee: no IndexOutOfRange / DivideByZero / NaN / hang on garbage.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void EntropyProfile_RandomGarbage_NeverThrows_MatchesOracle()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 60);
            string seq = RandomBmpChars(rng, len);
            int windowSize = rng.Next(0, len + 4);
            int stepSize = rng.Next(1, 5);

            IReadOnlyList<double> profile = null!;
            var act = () => profile = SequenceStatistics
                .CalculateEntropyProfile(seq, windowSize, stepSize).ToArray();

            act.Should().NotThrow(
                $"garbage len {len}, window {windowSize}, step {stepSize} must never crash");

            var oracle = OracleProfile(seq, windowSize, stepSize);
            profile.Should().HaveCount(oracle.Count,
                $"INV-05 count on garbage (len {len}, W {windowSize}, step {stepSize})");
            for (int k = 0; k < oracle.Count; k++)
            {
                double.IsFinite(profile[k]).Should().BeTrue("entropy must be finite on garbage");
                profile[k].Should().BeGreaterThanOrEqualTo(-Tolerance, "INV-01: H ≥ 0");
                profile[k].Should().BeApproximately(oracle[k], Tolerance,
                    $"window {k} matches spec oracle (letters only)");
            }
        }
    }

    #endregion
}
