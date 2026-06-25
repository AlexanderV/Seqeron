using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area sliding-window GC-content-profile unit.
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
/// Unit: SEQ-GC-PROFILE-001 — GC-content profile (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 234 (the FINAL row).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — windowSize &gt; length, windowSize = 0, empty /
///          null sequence, single char, plus windowSize = length, very long input,
///          and a randomized boundary sweep.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes; BE = 0, -1,
///   MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GC-content-profile contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateGcContentProfile(string, int windowSize = 100,
///            int stepSize = 1) → IEnumerable&lt;double&gt; (lazy / deferred).
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)
///
/// Documented behaviour (docs/algorithms/Statistics/GC_Content_Profile.md,
/// Test Unit ID SEQ-GC-PROFILE-001):
///   • §2.2 / §4.1: for each window of width W slid by step over the sequence,
///     emit GC% = (G + C) / (A + T + G + C) × 100. The per-window kernel
///     case-folds (ToUpperInvariant). The denominator counts ONLY the standard
///     bases A/T/U/G/C; ambiguous symbols (e.g. N) and all non-standard symbols
///     are EXCLUDED from the denominator (§2.4 INV-02).
///   • §2.4 INV-01: every profile value lies in [0, 100].
///   • §2.4 INV-02: denominator counts only A/T/U/G/C; N etc. excluded.
///   • §2.4 INV-03: window count = ⌊(n − W)/step⌋ + 1 for W ≤ n, else 0; offsets
///     are 0, step, 2·step, … ≤ n − W.
///   • §2.4 INV-04: U is a non-GC base equivalent to T.
///   • §2.4 INV-05 / §6.1: a window with no A/T/U/G/C base yields 0 (zero-division
///     convention) — NO divide-by-zero.
///   • §3.3 / §6.1: null / empty sequence, OR windowSize &gt; length → EMPTY profile
///     (no exception).
///   • §6.1: windowSize == length → exactly one window (whole-sequence GC%).
///   • §6.1: lowercase input == uppercase (case-folded counting).
///
/// The GC% formula GC% = (G+C)/(A+T+U+G+C) × 100 is re-implemented independently in
/// the <see cref="GcPercent"/> oracle below — expected values are derived from the
/// spec, never read off the implementation's own output.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceGcContentProfileFuzzTests
{
    #region SEQ-GC-PROFILE-001 — Helpers and independent GC% oracle

    private const double Tolerance = 1e-9;

    /// <summary>The DNA alphabet — standard bases used for random fixtures.</summary>
    private const string Dna = "ACGT";

    /// <summary>
    /// Independent oracle for the per-window GC% kernel: re-derives
    /// GC% = (G + C) / (A + T + U + G + C) × 100 from the SPEC (§2.2 / §2.4).
    /// Case-folds with ToUpperInvariant (§3.3); counts G/C in the numerator and
    /// A/T/U/G/C in the denominator; every other symbol (N, junk, non-letters) is
    /// EXCLUDED from BOTH (INV-02). A window with no standard base ⇒ 0 (INV-05,
    /// zero-division convention). Written from the formula, NOT copied from the
    /// implementation.
    /// </summary>
    private static double GcPercent(string window)
    {
        int gc = 0;
        int total = 0;
        foreach (char ch in window.ToUpperInvariant())
        {
            switch (ch)
            {
                case 'G':
                case 'C':
                    gc++;
                    total++;
                    break;
                case 'A':
                case 'T':
                case 'U':
                    total++;
                    break;
                // every other symbol is excluded from the denominator (INV-02).
            }
        }
        return total > 0 ? (double)gc / total * 100.0 : 0.0;
    }

    /// <summary>
    /// Independent oracle for the WHOLE profile: the documented window-offset rule
    /// (§4.1 / INV-03). Yields one GC% value per offset i = 0, step, 2·step, … while
    /// i ≤ length − W, when W ≤ length; otherwise nothing. Derived from the spec, not
    /// from the implementation's loop.
    /// </summary>
    private static List<double> OracleProfile(string seq, int windowSize, int stepSize)
    {
        var result = new List<double>();
        if (string.IsNullOrEmpty(seq) || windowSize > seq.Length)
            return result;
        for (int i = 0; i <= seq.Length - windowSize; i += stepSize)
            result.Add(GcPercent(seq.Substring(i, windowSize)));
        return result;
    }

    /// <summary>
    /// The documented window count: 0 for the null/empty-sequence guard (§3.3, which
    /// precedes the window logic — so even W = 0 on an empty sequence yields nothing),
    /// 0 when W &gt; n, else ⌊(n − W)/step⌋ + 1 (GC_Content_Profile.md §2.4 INV-03).
    /// Computed here from the spec directly, independent of the implementation.
    /// </summary>
    private static int ExpectedWindowCount(int n, int windowSize, int stepSize)
    {
        if (n == 0)
            return 0; // empty/null guard wins before any window is produced (§3.3)
        return windowSize > n ? 0 : (n - windowSize) / stepSize + 1;
    }

    /// <summary>
    /// Universal well-formedness contract for ANY input: every value is finite (never
    /// NaN/Infinity) and lies in [0, 100] (INV-01) — GC content as a percentage can
    /// never exceed 100 nor drop below 0, whatever the input alphabet.
    /// </summary>
    private static void AssertWellFormed(IReadOnlyList<double> profile)
    {
        profile.Should().NotBeNull("the method must always return a sequence, never null");
        foreach (double v in profile)
        {
            double.IsFinite(v).Should().BeTrue($"GC% {v} must be finite — never NaN/Infinity");
            v.Should().BeGreaterThanOrEqualTo(-Tolerance, "INV-01: every profile value ≥ 0");
            v.Should().BeLessThanOrEqualTo(100.0 + Tolerance, "INV-01: every profile value ≤ 100");
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
    //  SEQ-GC-PROFILE-001 — GC-content profile : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-checkable worked example (§7.1)

    /// <summary>
    /// Positive baseline (not a boundary): the doc's worked example must be reproduced
    /// EXACTLY. "GGGAAATGCC" with window 4, step 3 has ⌊(10−4)/3⌋+1 = 3 windows at
    /// offsets {0, 3, 6}:
    ///   GGGA → G+C=3, total=4 → 3/4×100 = 75.0
    ///   AAAT → G+C=0, total=4 → 0.0
    ///   TGCC → G+C=3, total=4 → 75.0
    /// The three GC% values are computed here independently from the formula (literals
    /// re-derived by hand) and re-checked against the kernel oracle — this pins the
    /// BUSINESS contract (GC% per window), not mere absence of crashes.
    /// — GC_Content_Profile.md §7.1.
    /// </summary>
    [Test]
    public void GcProfile_WorkedExample_MatchesHandComputedPercents()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGAAATGCC", windowSize: 4, stepSize: 3)
            .ToArray();

        profile.Should().HaveCount(3, "⌊(10−4)/3⌋+1 = 3 windows at offsets {0,3,6} (INV-03)");

        // Hand-derived literals from §7.1.
        profile[0].Should().BeApproximately(75.0, Tolerance, "GGGA → 3/4×100 = 75.0");
        profile[1].Should().BeApproximately(0.0, Tolerance, "AAAT → 0/4×100 = 0.0");
        profile[2].Should().BeApproximately(75.0, Tolerance, "TGCC → 3/4×100 = 75.0");

        // Cross-check the hand literals against the independent GC% oracle.
        profile[0].Should().BeApproximately(GcPercent("GGGA"), Tolerance);
        profile[1].Should().BeApproximately(GcPercent("AAAT"), Tolerance);
        profile[2].Should().BeApproximately(GcPercent("TGCC"), Tolerance);

        AssertWellFormed(profile);
    }

    /// <summary>
    /// Positive: the boundary GC% anchors. An all-GC window "GCGC" yields 100.0; an
    /// all-AT window "ATAT" yields 0.0; a balanced "ACGT" yields 50.0 (2 of 4). The
    /// step parameter selects offsets 0, step, 2·step, … — "GCGCATAT" with window 4,
    /// step 4 picks offsets {0,4} → GCGC (100.0), ATAT (0.0). Confirms step is honoured
    /// (not always 1) and pins the [0,100] extremes. — GC_Content_Profile.md §2.2 /
    /// INV-01 / §4.1.
    /// </summary>
    [Test]
    public void GcProfile_ExtremesAndStep_AreHonoured()
    {
        var allGc = SequenceStatistics.CalculateGcContentProfile("GCGC", windowSize: 4).ToArray();
        allGc.Should().HaveCount(1);
        allGc[0].Should().BeApproximately(100.0, Tolerance, "all G/C ⇒ 100.0");

        var allAt = SequenceStatistics.CalculateGcContentProfile("ATAT", windowSize: 4).ToArray();
        allAt.Should().HaveCount(1);
        allAt[0].Should().BeApproximately(0.0, Tolerance, "all A/T ⇒ 0.0");

        var balanced = SequenceStatistics.CalculateGcContentProfile("ACGT", windowSize: 4).ToArray();
        balanced[0].Should().BeApproximately(50.0, Tolerance, "2 of 4 are G/C ⇒ 50.0");

        var strided = SequenceStatistics.CalculateGcContentProfile("GCGCATAT", windowSize: 4, stepSize: 4)
            .ToArray();
        strided.Should().HaveCount(2, "offsets {0,4}: ⌊(8−4)/4⌋+1 = 2 windows");
        strided[0].Should().BeApproximately(100.0, Tolerance, "offset 0 → GCGC = 100.0");
        strided[1].Should().BeApproximately(0.0, Tolerance, "offset 4 → ATAT = 0.0");

        AssertWellFormed(allGc);
        AssertWellFormed(allAt);
        AssertWellFormed(strided);
    }

    /// <summary>
    /// Positive: INV-02 (N excluded from denominator) and INV-04 (U ≡ T, a non-GC
    /// base) as exact numeric anchors, mirroring the Biopython `gc_fraction("ACTGN")
    /// = 0.50` example in §2.2 [2]. "ACTGN" (window 5) ⇒ G+C = 2 (C,G), denominator =
    /// 4 (A,C,T,G; N excluded) ⇒ 2/4×100 = 50.0, NOT 2/5×100 = 40.0. A window "GGUU"
    /// counts U as a non-GC base like T ⇒ G+C = 2, total = 4 ⇒ 50.0.
    /// — GC_Content_Profile.md §2.2 / INV-02 / INV-04.
    /// </summary>
    [Test]
    public void GcProfile_AmbiguousNExcluded_AndUisNonGc()
    {
        var withN = SequenceStatistics.CalculateGcContentProfile("ACTGN", windowSize: 5).ToArray();
        withN.Should().HaveCount(1, "⌊(5−5)/1⌋+1 = 1 window");
        withN[0].Should().BeApproximately(50.0, Tolerance,
            "N excluded from denominator: 2 G+C / 4 standard ×100 = 50.0 (Biopython gc_fraction 0.50)");
        withN[0].Should().BeApproximately(GcPercent("ACTGN"), Tolerance);

        var withU = SequenceStatistics.CalculateGcContentProfile("GGUU", windowSize: 4).ToArray();
        withU.Should().HaveCount(1);
        withU[0].Should().BeApproximately(50.0, Tolerance,
            "U is a non-GC base ≡ T: 2 G+C / 4 standard ×100 = 50.0 (INV-04)");
        withU[0].Should().BeApproximately(GcPercent("GGUU"), Tolerance);
    }

    /// <summary>
    /// Positive: lowercase input is case-folded before counting (§3.3 / §6.1) and must
    /// yield the same profile as the uppercase equivalent. "gcat" with window 2 step 1
    /// equals "GCAT". — GC_Content_Profile.md §3.3 / §6.1 (lowercase == uppercase).
    /// </summary>
    [Test]
    public void GcProfile_LowercaseEqualsUppercase()
    {
        var lower = SequenceStatistics.CalculateGcContentProfile("gcat", windowSize: 2, stepSize: 1)
            .ToArray();
        var upper = SequenceStatistics.CalculateGcContentProfile("GCAT", windowSize: 2, stepSize: 1)
            .ToArray();

        lower.Should().Equal(upper, "case-folded counting: lowercase == uppercase (§3.3)");
        // GC, CA, AT → 100.0, 50.0, 0.0
        lower.Should().HaveCount(3);
        lower[0].Should().BeApproximately(100.0, Tolerance, "gc → 100.0");
        lower[1].Should().BeApproximately(50.0, Tolerance, "ca → 50.0");
        lower[2].Should().BeApproximately(0.0, Tolerance, "at → 0.0");
        AssertWellFormed(lower);
    }

    #endregion

    #region BE — Boundary: windowSize > length (no full window ⇒ empty profile)

    /// <summary>
    /// BE: when the window is wider than the sequence, no full window exists, so the
    /// documented result is the EMPTY profile (no exception, no negative-length
    /// Substring). Covered just over the edge (W = n+1) and far past it (W = 1000),
    /// across short DNA inputs. — GC_Content_Profile.md §3.3 / §6.1 / INV-03 (W &gt; n ⇒ 0).
    /// </summary>
    [TestCase("ACGT", 5)]
    [TestCase("ACGT", 1000)]
    [TestCase("A", 2)]
    [TestCase("GCGCATATGC", 11)]
    public void GcProfile_WindowWiderThanSequence_IsEmpty(string seq, int windowSize)
    {
        var act = () => SequenceStatistics.CalculateGcContentProfile(seq, windowSize, stepSize: 1).ToArray();

        act.Should().NotThrow("windowSize > length is a defined boundary, not an error");
        act().Should().BeEmpty("no full window exists when W > length (INV-03)");
    }

    #endregion

    #region BE — Boundary: windowSize == length (exactly one window)

    /// <summary>
    /// BE: when the window exactly equals the sequence length, INV-03 gives
    /// ⌊(n−n)/step⌋+1 = 1 — a single window spanning the whole sequence, whose value
    /// is the GC% of the entire input. "GCGC" (window 4) ⇒ one value, 100.0; "ATAT" ⇒
    /// one value, 0.0; "ACGT" ⇒ 50.0. — GC_Content_Profile.md §6.1 (W == length).
    /// </summary>
    [Test]
    public void GcProfile_WindowEqualsLength_IsSingleWholeSequenceValue()
    {
        SequenceStatistics.CalculateGcContentProfile("GCGC", windowSize: 4).ToArray()
            .Should().ContainSingle().Which.Should().BeApproximately(100.0, Tolerance,
                "INV-03: one window = whole 'GCGC' = 100.0");

        SequenceStatistics.CalculateGcContentProfile("ATAT", windowSize: 4).ToArray()
            .Should().ContainSingle().Which.Should().BeApproximately(0.0, Tolerance,
                "whole 'ATAT' = 0.0");

        SequenceStatistics.CalculateGcContentProfile("ACGT", windowSize: 4).ToArray()
            .Should().ContainSingle().Which.Should().BeApproximately(50.0, Tolerance,
                "whole 'ACGT' = 50.0");
    }

    #endregion

    #region BE — Boundary: windowSize == 0 (degenerate parameter, defined result)

    /// <summary>
    /// BE: windowSize = 0 is the degenerate zero-window-width boundary. It is NOT
    /// guarded by the W &gt; length check (0 ≤ length), so by INV-03 the profile has
    /// ⌊(n−0)/step⌋+1 = n+1 windows, each a length-0 substring whose GC% (no counted
    /// bases, total = 0) is the documented 0.0 (INV-05, zero-division convention). The
    /// key fuzz guarantee: NO divide-by-zero (the kernel returns 0 on a zero total) and
    /// NO infinite loop (step ≥ 1 advances the offset to termination at i = n). We pin
    /// the exact shape: n+1 zero values. — GC_Content_Profile.md §4.1 (offsets while
    /// i ≤ n − W = n) / INV-03 / INV-05.
    /// </summary>
    [TestCase("ACGT", 5)]   // n = 4 → 5 windows
    [TestCase("A", 2)]      // n = 1 → 2 windows
    [TestCase("GCGCATATGC", 11)]
    [CancelAfter(5000)]
    public void GcProfile_WindowZero_YieldsZeroPercentPerOffset_NoCrashNoHang(string seq, int expectedCount)
    {
        IReadOnlyList<double> profile = null!;
        var act = () => profile = SequenceStatistics.CalculateGcContentProfile(seq, windowSize: 0, stepSize: 1)
            .ToArray();

        act.Should().NotThrow("windowSize = 0 must not divide by zero or throw");
        profile.Should().HaveCount(expectedCount, "INV-03: ⌊(n−0)/1⌋+1 = n+1 zero-width windows");
        profile.Should().OnlyContain(v => v == 0.0,
            "a zero-width window has no counted bases ⇒ total 0 ⇒ GC% 0.0 (INV-05)");
        AssertWellFormed(profile);
    }

    /// <summary>
    /// BE: windowSize = 0 on the EMPTY string is guarded by the null/empty
    /// short-circuit (§3.3) ⇒ empty profile, even though 0 ≤ 0. Confirms the empty
    /// guard precedes the zero-window logic. — GC_Content_Profile.md §3.3.
    /// </summary>
    [Test]
    public void GcProfile_WindowZero_OnEmptySequence_IsEmpty()
    {
        var act = () => SequenceStatistics.CalculateGcContentProfile(string.Empty, windowSize: 0).ToArray();

        act.Should().NotThrow();
        act().Should().BeEmpty("empty sequence is guarded before any window is produced (§3.3)");
    }

    #endregion

    #region BE — Boundary: empty / null sequence (defined ⇒ empty profile)

    /// <summary>
    /// BE: the empty string is the lower size boundary — null/empty yields the EMPTY
    /// profile (no exception). Asserted across several window sizes (incl. the default
    /// 100). — GC_Content_Profile.md §3.3 / §6.1 (null/empty → empty profile).
    /// </summary>
    [TestCase(1)]
    [TestCase(4)]
    [TestCase(100)]
    public void GcProfile_EmptySequence_IsEmpty(int windowSize)
    {
        var act = () => SequenceStatistics.CalculateGcContentProfile(string.Empty, windowSize).ToArray();

        act.Should().NotThrow("the empty string is a defined boundary, not an error");
        act().Should().BeEmpty("no bases ⇒ no window ⇒ empty profile");
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — empty
    /// profile, no NullReferenceException. — GC_Content_Profile.md §3.3 (null → empty).
    /// </summary>
    [Test]
    public void GcProfile_NullSequence_IsEmpty_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateGcContentProfile(null!, windowSize: 4).ToArray();

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().BeEmpty();
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a single-character sequence. With the default window (100 &gt; 1) no full
    /// window exists ⇒ empty profile. With window 1 there is exactly one window whose
    /// single base gives the GC% of that base: "G" ⇒ 100.0, "A" ⇒ 0.0. With window 1
    /// and a 2-char "GA" there are two windows: 100.0, 0.0.
    /// — GC_Content_Profile.md §6.1 (W &gt; length → empty; INV-01).
    /// </summary>
    [Test]
    public void GcProfile_SingleChar_DefaultWindowEmpty_Window1IsBaseGc()
    {
        // Default window 100 > length 1 → empty.
        SequenceStatistics.CalculateGcContentProfile("G").ToArray()
            .Should().BeEmpty("default window 100 > single-char length ⇒ empty (INV-03)");

        // Window 1 over a single G → one window, 100.0.
        SequenceStatistics.CalculateGcContentProfile("G", windowSize: 1).ToArray()
            .Should().ContainSingle().Which.Should().BeApproximately(100.0, Tolerance,
                "one G ⇒ 1/1×100 = 100.0");

        // Window 1 over a single A → one window, 0.0.
        SequenceStatistics.CalculateGcContentProfile("A", windowSize: 1).ToArray()
            .Should().ContainSingle().Which.Should().BeApproximately(0.0, Tolerance,
                "one A ⇒ 0/1×100 = 0.0");

        // Window 1 over "GA" → two windows, 100.0 then 0.0.
        var two = SequenceStatistics.CalculateGcContentProfile("GA", windowSize: 1).ToArray();
        two.Should().HaveCount(2);
        two[0].Should().BeApproximately(100.0, Tolerance, "G ⇒ 100.0");
        two[1].Should().BeApproximately(0.0, Tolerance, "A ⇒ 0.0");
        AssertWellFormed(two);
    }

    /// <summary>
    /// BE: a single NON-standard character with window 1 ⇒ one window whose only base
    /// is excluded from the denominator (total = 0) ⇒ GC% 0.0 (INV-05, zero-division
    /// convention), with NO divide-by-zero. — GC_Content_Profile.md §6.1 (all-N
    /// window → 0) / INV-05.
    /// </summary>
    [Test]
    public void GcProfile_SingleNonStandardChar_IsZero_NoDivideByZero()
    {
        var n = SequenceStatistics.CalculateGcContentProfile("N", windowSize: 1).ToArray();
        n.Should().ContainSingle().Which.Should().BeApproximately(0.0, Tolerance,
            "N excluded from denominator ⇒ total 0 ⇒ 0.0 (INV-05)");
        AssertWellFormed(n);
    }

    #endregion

    #region BE — non-standard symbols excluded from the denominator

    /// <summary>
    /// BE/MC: the per-window kernel counts ONLY A/T/U/G/C in the denominator (§2.4
    /// INV-02). A window of pure non-standard symbols (digits, punctuation, the null
    /// byte, N) has total = 0 ⇒ 0.0, with no divide-by-zero. A mixed window "GCNN"
    /// (window 4) has G+C = 2 over a denominator of 2 (the two N excluded) ⇒ 100.0,
    /// NOT 50.0 — proving N is removed from the denominator rather than treated as a
    /// non-GC base. — GC_Content_Profile.md §2.2 / INV-02.
    /// </summary>
    [Test]
    public void GcProfile_NonStandardExcludedFromDenominator()
    {
        // Pure non-standard window → 0.0, no crash.
        var junk = SequenceStatistics.CalculateGcContentProfile("12-+", windowSize: 4).ToArray();
        junk.Should().ContainSingle().Which.Should().BeApproximately(0.0, Tolerance,
            "no counted bases ⇒ total 0 ⇒ 0.0");

        // N removed from denominator (not counted as a non-GC base).
        var mixed = SequenceStatistics.CalculateGcContentProfile("GCNN", windowSize: 4).ToArray();
        mixed.Should().ContainSingle().Which.Should().BeApproximately(100.0, Tolerance,
            "GCNN: 2 G+C / 2 standard ×100 = 100.0 (N excluded, INV-02)");
        mixed[0].Should().BeApproximately(GcPercent("GCNN"), Tolerance);
    }

    #endregion

    #region BE — Boundary: very long sequence (O(n·W), no overflow / hang)

    /// <summary>
    /// BE: a long sequence must be processed without overflow, hang or NaN. A
    /// 200 000-base "GCAT" tiling with window 4, step 4 gives 50 000 windows, each
    /// "GCAT" ⇒ 2 G+C of 4 ⇒ 50.0. Bounds runtime under [CancelAfter] and pins INV-03
    /// at scale. — GC_Content_Profile.md §4.3 (O(W·windowSize)) / INV-01 / INV-03.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GcProfile_VeryLongTiledSequence_AllFiftyPercent()
    {
        const int tiles = 50_000;
        string seq = string.Concat(Enumerable.Repeat("GCAT", tiles)); // 200 000 bases

        var profile = SequenceStatistics.CalculateGcContentProfile(seq, windowSize: 4, stepSize: 4)
            .ToArray();

        profile.Should().HaveCount(tiles, "non-overlapping GCAT tiles: 50 000 windows");
        profile.Should().OnlyContain(v => Math.Abs(v - 50.0) <= Tolerance,
            "every tiled window 'GCAT' ⇒ 2/4×100 = 50.0");
        AssertWellFormed(profile);
    }

    #endregion

    #region BE / RB — Randomized boundary sweep: never crash, match oracle

    /// <summary>
    /// BE: a randomized sweep over the documented boundary space — random DNA length
    /// (incl. 0 and 1), random windowSize that straddles 0, &lt; len, == len and &gt; len,
    /// and random step ≥ 1. The profile must EXACTLY match the independent offset-rule
    /// + GC% oracle (count = INV-03; values = §2.2) and be well-formed in [0,100], with
    /// no crash, hang, NaN or Infinity. Locally seeded Random.
    /// — GC_Content_Profile.md §2.2 / §4.1 / INV-01..INV-05.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void GcProfile_RandomDnaBoundarySweep_MatchesOracle()
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
                .CalculateGcContentProfile(seq, windowSize, stepSize).ToArray();

            act.Should().NotThrow(
                $"DNA len {len}, window {windowSize}, step {stepSize} must never crash");

            var oracle = OracleProfile(seq, windowSize, stepSize);
            profile.Should().HaveCount(ExpectedWindowCount(len, windowSize, stepSize),
                $"INV-03 window count (len {len}, W {windowSize}, step {stepSize})");
            profile.Should().HaveCount(oracle.Count);
            for (int k = 0; k < oracle.Count; k++)
                profile[k].Should().BeApproximately(oracle[k], Tolerance,
                    $"window {k} GC% (len {len}, W {windowSize}, step {stepSize})");

            AssertWellFormed(profile);
        }
    }

    /// <summary>
    /// BE/RB: arbitrary BMP garbage (control chars, null byte, lone surrogate halves,
    /// unicode letters/digits/symbols) with random window/step must NEVER throw and
    /// must ALWAYS match the spec oracle (only A/T/U/G/C counted in the denominator,
    /// G/C in the numerator; everything else excluded; §2.4). Each value must be finite
    /// and within [0,100]. Core fuzz guarantee: no IndexOutOfRange / DivideByZero / NaN
    /// / hang on garbage. — GC_Content_Profile.md §2.4 / INV-01 / INV-02.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void GcProfile_RandomGarbage_NeverThrows_MatchesOracle()
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
                .CalculateGcContentProfile(seq, windowSize, stepSize).ToArray();

            act.Should().NotThrow(
                $"garbage len {len}, window {windowSize}, step {stepSize} must never crash");

            var oracle = OracleProfile(seq, windowSize, stepSize);
            profile.Should().HaveCount(oracle.Count,
                $"INV-03 count on garbage (len {len}, W {windowSize}, step {stepSize})");
            for (int k = 0; k < oracle.Count; k++)
            {
                double.IsFinite(profile[k]).Should().BeTrue("GC% must be finite on garbage");
                profile[k].Should().BeInRange(-Tolerance, 100.0 + Tolerance, "INV-01: 0 ≤ GC% ≤ 100");
                profile[k].Should().BeApproximately(oracle[k], Tolerance,
                    $"window {k} matches spec oracle (A/T/U/G/C only)");
            }
        }
    }

    #endregion
}
