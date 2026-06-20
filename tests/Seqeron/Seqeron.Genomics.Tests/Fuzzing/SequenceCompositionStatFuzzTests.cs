using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area sequence-composition unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception. A raw runtime exception or a
/// hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPOSITION-001 — Sequence composition (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 121.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — empty, single base, all-N, lowercase,
///          non-ACGT (degenerate IUPAC / junk).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The composition contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateNucleotideComposition(string)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 49–94), returning a NucleotideComposition record struct.
///
/// NOTE: SEQ-COMPOSITION-001 is the Statistics-area composition unit; it is a
/// CONSOLIDATED duplicate of SEQ-STATS-001 and shares the same method/doc. It is
/// distinct from row 2 SEQ-COMP-001 (DNA complement, CompositionFuzzTests.cs).
///   — docs/algorithms/Sequence_Composition/Sequence_Composition.md (Test Unit ID
///     SEQ-COMPOSITION-001, "Consolidation note"); doc §2.2, §3.3, §6.1;
///     docs/algorithms/Sequence_Composition/Sequence_Composition_Statistics.md
///     (SEQ-STATS-001) §2.2, §2.4 (INV-01..INV-04), §3.3.
///
/// Documented behaviour (Sequence_Composition.md §3.3 / §6.1, §2.2):
///   • null / empty input → all-zero NucleotideComposition (no exception, no
///     DivideByZero); GcContent = AtContent = GcSkew = AtSkew = 0.
///   • Counting is case-insensitive: input is upper-cased (ToUpperInvariant)
///     before classification, so lowercase round-trips to the identical result.
///   • Standard alphabet is {A,T,G,C,U}. 'N' → CountN; any other non-canonical
///     character (degenerate IUPAC S/W/R/Y/…, digits, whitespace, unicode) →
///     CountOther. N and Other do NOT contribute to the GC/AT totals.
///   • total = A+T+G+C+U.  GcContent = (G+C)/total, AtContent = (A+T+U)/total,
///     each 0 when total = 0 (zero-denominator guard).
///   • GcSkew = (G−C)/(G+C), 0 when G+C = 0; AtSkew = (A−T)/(A+T), 0 when A+T = 0.
///
/// Invariants pinned (Sequence_Composition.md §2.4):
///   • INV-01: 0 ≤ GcContent ≤ 1, 0 ≤ AtContent ≤ 1.
///   • INV-02: −1 ≤ GcSkew ≤ 1, −1 ≤ AtSkew ≤ 1.
///   • INV-03: CountA+CountT+CountG+CountC+CountU+CountN+CountOther = Length
///             (the counts partition every character).
///   • Fractions-sum: when total &gt; 0, GcContent + AtContent = 1 exactly
///             (GC numerator + AT numerator = G+C + A+T+U = total).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceCompositionStatFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The canonical alphabet that contributes to the GC/AT totals.</summary>
    private const string CanonicalBases = "ATGCU";

    /// <summary>Generates a random string of arbitrary BMP code points (0x0000–0xFFFF),
    /// spanning control characters, the null byte, lone surrogate halves, unicode
    /// letters and digits — random-byte fuzz fodder for the classifier.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal well-formedness contract that must hold for ANY input:
    /// non-negative counts, the count partition (INV-03), bounded fractions and skews
    /// (INV-01/INV-02), the zero-denominator guards, and the fractions-sum identity.
    /// </summary>
    private static void AssertWellFormed(
        SequenceStatistics.NucleotideComposition c, int expectedLength)
    {
        c.Length.Should().Be(expectedLength, "Length must equal the input character count");

        c.CountA.Should().BeGreaterThanOrEqualTo(0);
        c.CountT.Should().BeGreaterThanOrEqualTo(0);
        c.CountG.Should().BeGreaterThanOrEqualTo(0);
        c.CountC.Should().BeGreaterThanOrEqualTo(0);
        c.CountU.Should().BeGreaterThanOrEqualTo(0);
        c.CountN.Should().BeGreaterThanOrEqualTo(0);
        c.CountOther.Should().BeGreaterThanOrEqualTo(0);

        // INV-03: counts partition every character.
        (c.CountA + c.CountT + c.CountG + c.CountC + c.CountU + c.CountN + c.CountOther)
            .Should().Be(c.Length, "INV-03: every character lands in exactly one bucket");

        // INV-01: content fractions bounded.
        c.GcContent.Should().BeInRange(0.0, 1.0, "INV-01");
        c.AtContent.Should().BeInRange(0.0, 1.0, "INV-01");

        // INV-02: skews bounded.
        c.GcSkew.Should().BeInRange(-1.0, 1.0, "INV-02");
        c.AtSkew.Should().BeInRange(-1.0, 1.0, "INV-02");

        int total = c.CountA + c.CountT + c.CountG + c.CountC + c.CountU;
        if (total == 0)
        {
            // Zero-denominator guard: no DivideByZero, all defined as 0.
            c.GcContent.Should().Be(0.0, "no canonical base ⇒ GcContent guarded to 0");
            c.AtContent.Should().Be(0.0, "no canonical base ⇒ AtContent guarded to 0");
        }
        else
        {
            // Fractions-sum identity: (G+C)/total + (A+T+U)/total = 1.
            (c.GcContent + c.AtContent).Should().BeApproximately(1.0, Tolerance,
                "GC numerator + AT numerator = G+C + A+T+U = total");
        }

        if (c.CountG + c.CountC == 0)
            c.GcSkew.Should().Be(0.0, "INV-04: GcSkew guarded to 0 when G+C = 0");
        if (c.CountA + c.CountT == 0)
            c.AtSkew.Should().Be(0.0, "INV-04: AtSkew guarded to 0 when A+T = 0");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMPOSITION-001 — sequence composition : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact result

    /// <summary>
    /// Positive baseline (not a boundary): a known mixed sequence must yield the
    /// documented per-base counts and fractions EXACTLY. "AATTGGGC" → A2 T2 G3 C1,
    /// total 8, GcContent = 4/8 = 0.5, AtContent = 4/8 = 0.5, GcSkew = (3−1)/4 = 0.5,
    /// AtSkew = (2−2)/4 = 0. Confirms the suite asserts the BUSINESS contract.
    /// — Sequence_Composition.md §7.1 worked example (GGGC: GcSkew (3−1)/4 = 0.5).
    /// </summary>
    [Test]
    public void Composition_KnownSequence_MatchesHandComputedCountsAndFractions()
    {
        var c = SequenceStatistics.CalculateNucleotideComposition("AATTGGGC");

        c.Length.Should().Be(8);
        c.CountA.Should().Be(2);
        c.CountT.Should().Be(2);
        c.CountG.Should().Be(3);
        c.CountC.Should().Be(1);
        c.CountU.Should().Be(0);
        c.CountN.Should().Be(0);
        c.CountOther.Should().Be(0);

        c.GcContent.Should().BeApproximately(0.5, Tolerance);
        c.AtContent.Should().BeApproximately(0.5, Tolerance);
        c.GcSkew.Should().BeApproximately(0.5, Tolerance, "(G−C)/(G+C) = (3−1)/4");
        c.AtSkew.Should().BeApproximately(0.0, Tolerance, "(A−T)/(A+T) = (2−2)/4");

        // Fractions over the counted alphabet sum to 1.
        (c.GcContent + c.AtContent).Should().BeApproximately(1.0, Tolerance);

        AssertWellFormed(c, 8);
    }

    #endregion

    #region BE — Boundary: empty / null

    /// <summary>
    /// BE: the empty string is the lower size boundary. Documented as an all-zero
    /// composition — NO DivideByZero on length 0, no exception.
    /// — Sequence_Composition.md §6.1 (empty/null → all-zero composition).
    /// </summary>
    [Test]
    public void Composition_EmptyString_IsAllZeroAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateNucleotideComposition(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        var c = act();
        AssertWellFormed(c, 0);
        c.GcContent.Should().Be(0.0);
        c.AtContent.Should().Be(0.0);
        c.GcSkew.Should().Be(0.0);
        c.AtSkew.Should().Be(0.0);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 51) — all-zero composition, no NullReference.
    /// — Sequence_Composition.md §3.3 (null/empty allowed; all-zero, no exception).
    /// </summary>
    [Test]
    public void Composition_Null_IsAllZeroAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateNucleotideComposition(null!);

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        AssertWellFormed(act(), 0);
    }

    #endregion

    #region BE — Boundary: single base

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. The composition is
    /// the binary extreme — the single base's count is 1, all others 0; GcContent
    /// is 1 for G/C and 0 for A/T; AtContent is the complement. U is the RNA base:
    /// counted in CountU, AT-content (A+T+U) numerator, GcContent 0.
    /// — Sequence_Composition.md §2.2 (per-base counts; GC/AT content formulas).
    /// </summary>
    [TestCase('G', 1.0)]
    [TestCase('C', 1.0)]
    [TestCase('A', 0.0)]
    [TestCase('T', 0.0)]
    [TestCase('U', 0.0)]
    public void Composition_SingleCanonicalBase_IsBinaryExtreme(char baseChar, double expectedGc)
    {
        var c = SequenceStatistics.CalculateNucleotideComposition(baseChar.ToString());

        c.Length.Should().Be(1);
        c.GcContent.Should().BeApproximately(expectedGc, Tolerance);
        c.AtContent.Should().BeApproximately(1.0 - expectedGc, Tolerance,
            "with one canonical base GcContent + AtContent = 1");

        // The matching counter is exactly 1 and is the only non-zero one.
        int matched = baseChar switch
        {
            'G' => c.CountG, 'C' => c.CountC, 'A' => c.CountA,
            'T' => c.CountT, 'U' => c.CountU, _ => -1
        };
        matched.Should().Be(1, $"the single '{baseChar}' increments only its own counter");

        AssertWellFormed(c, 1);
    }

    /// <summary>
    /// BE: a single 'N' — minimal all-N input. CountN = 1, no canonical base, so
    /// total = 0 and every fraction/skew is the zero-denominator-guarded 0. No crash.
    /// </summary>
    [Test]
    public void Composition_SingleN_IsCountedAsNWithZeroFractions()
    {
        var c = SequenceStatistics.CalculateNucleotideComposition("N");

        c.Length.Should().Be(1);
        c.CountN.Should().Be(1);
        c.CountOther.Should().Be(0);
        c.GcContent.Should().Be(0.0);
        c.AtContent.Should().Be(0.0);
        AssertWellFormed(c, 1);
    }

    #endregion

    #region BE — Boundary: all-N

    /// <summary>
    /// BE: every base is 'N'. N is the documented IUPAC "any base" symbol routed to
    /// CountN; it is EXCLUDED from the GC/AT totals, so total = 0 and the
    /// zero-denominator guards keep GcContent/AtContent/GcSkew/AtSkew at 0 with NO
    /// DivideByZero. CountN must equal Length.
    /// — Sequence_Composition.md §6.1 (N/degenerate counted as N/Other; excluded
    ///   from GC/AT totals).
    /// </summary>
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(50)]
    [TestCase(1000)]
    public void Composition_AllN_CountedSeparately_NoDivideByZero(int length)
    {
        string seq = new string('N', length);

        var c = SequenceStatistics.CalculateNucleotideComposition(seq);

        c.CountN.Should().Be(length, "every character is N");
        c.CountA.Should().Be(0);
        c.CountT.Should().Be(0);
        c.CountG.Should().Be(0);
        c.CountC.Should().Be(0);
        c.CountU.Should().Be(0);
        c.CountOther.Should().Be(0);
        c.GcContent.Should().Be(0.0, "no canonical base ⇒ guarded to 0");
        c.AtContent.Should().Be(0.0);
        c.GcSkew.Should().Be(0.0);
        c.AtSkew.Should().Be(0.0);

        AssertWellFormed(c, length);
    }

    /// <summary>
    /// BE: lowercase 'n' must be classified identically to 'N' (ToUpperInvariant) —
    /// guards against a case-sensitivity bug in the N branch.
    /// </summary>
    [Test]
    public void Composition_LowercaseN_CountedAsN()
    {
        var c = SequenceStatistics.CalculateNucleotideComposition("nnnn");

        c.CountN.Should().Be(4, "lowercase n is upper-cased before counting");
        c.CountOther.Should().Be(0);
        AssertWellFormed(c, 4);
    }

    #endregion

    #region BE — Boundary: lowercase / mixed case

    /// <summary>
    /// BE: counting is case-insensitive. A lowercase sequence must produce the
    /// IDENTICAL composition to its upper-case form (ToUpperInvariant, line 58) —
    /// guards against a case-sensitivity bug that would route a/c/g/t to CountOther.
    /// — Sequence_Composition.md §6.1 (lowercase/mixed case → same as uppercase).
    /// </summary>
    [TestCase("acgtu")]
    [TestCase("AcGtU")]
    [TestCase("aattgggc")]
    public void Composition_LowercaseOrMixedCase_EqualsUppercase(string seq)
    {
        var lower = SequenceStatistics.CalculateNucleotideComposition(seq);
        var upper = SequenceStatistics.CalculateNucleotideComposition(seq.ToUpperInvariant());

        lower.Should().Be(upper, "ToUpperInvariant makes counting case-insensitive");

        // And nothing leaked into CountOther — every standard base was recognized.
        lower.CountOther.Should().Be(0, "lowercase standard bases must not be 'Other'");
        AssertWellFormed(lower, seq.Length);
    }

    /// <summary>
    /// BE: a single lowercase 'g' is a GC base (GcContent 1.0), proving lowercase
    /// does not silently fall through to the Other bucket (which would give 0).
    /// </summary>
    [Test]
    public void Composition_SingleLowercaseGcBase_IsGcContentOne()
    {
        var c = SequenceStatistics.CalculateNucleotideComposition("g");

        c.CountG.Should().Be(1);
        c.GcContent.Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(c, 1);
    }

    #endregion

    #region BE — Boundary: non-ACGT (degenerate IUPAC / junk)

    /// <summary>
    /// BE: degenerate IUPAC ambiguity letters (S, W, R, Y, K, M, B, D, H, V) are NOT
    /// in the canonical alphabet; they are routed to CountOther and EXCLUDED from the
    /// GC/AT totals (intentional simplification, doc §5.3). No crash, fractions stay
    /// well-defined. Here total = 0, so all fractions are the guarded 0.
    /// — Sequence_Composition.md §5.3 / §6.1 (degenerate → CountOther, not GC/AT).
    /// </summary>
    [Test]
    public void Composition_DegenerateIupacOnly_RoutedToOther_NoCanonicalContribution()
    {
        const string degenerate = "SWRYKMBDHVswrykmbdhv";

        var c = SequenceStatistics.CalculateNucleotideComposition(degenerate);

        c.CountOther.Should().Be(degenerate.Length, "no degenerate code is canonical");
        c.CountA.Should().Be(0);
        c.CountT.Should().Be(0);
        c.CountG.Should().Be(0);
        c.CountC.Should().Be(0);
        c.CountU.Should().Be(0);
        c.CountN.Should().Be(0, "N has its own bucket; other degenerate codes are 'Other'");
        c.GcContent.Should().Be(0.0);
        c.AtContent.Should().Be(0.0);
        AssertWellFormed(c, degenerate.Length);
    }

    /// <summary>
    /// BE: junk symbols (digits, whitespace, punctuation, the null byte) go to
    /// CountOther and never affect the canonical fractions. The canonical bases
    /// interleaved with junk still produce the documented fractions over the
    /// canonical-only total. "A1C\t2G\0" → A1 C1 G1 (total 3), Other = digits/ws/nul.
    /// </summary>
    [Test]
    public void Composition_CanonicalMixedWithJunk_FractionsOverCanonicalTotalOnly()
    {
        const string seq = "A1C\t2G\0";

        var c = SequenceStatistics.CalculateNucleotideComposition(seq);

        c.Length.Should().Be(seq.Length);
        c.CountA.Should().Be(1);
        c.CountC.Should().Be(1);
        c.CountG.Should().Be(1);
        c.CountT.Should().Be(0);
        c.CountU.Should().Be(0);
        c.CountN.Should().Be(0);
        c.CountOther.Should().Be(4, "'1', '\\t', '2', '\\0' are all Other");

        // total = A+C+G = 3; GcContent = (G+C)/3 = 2/3; AtContent = A/3 = 1/3.
        c.GcContent.Should().BeApproximately(2.0 / 3.0, Tolerance);
        c.AtContent.Should().BeApproximately(1.0 / 3.0, Tolerance);
        (c.GcContent + c.AtContent).Should().BeApproximately(1.0, Tolerance);

        AssertWellFormed(c, seq.Length);
    }

    #endregion

    #region BE — Random / RB fuzz: never throws, always well-formed

    /// <summary>
    /// BE/RB: a large batch of arbitrary BMP strings (control chars, null byte, lone
    /// surrogate halves, unicode letters/digits) must NEVER throw and must ALWAYS
    /// produce a well-formed composition satisfying every invariant. This is the core
    /// fuzz guarantee: no DivideByZero, no IndexOutOfRange, no overflow on garbage.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Composition_RandomGarbageStrings_NeverThrow_AlwaysWellFormed()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            SequenceStatistics.NucleotideComposition c = default;
            var act = () => c = SequenceStatistics.CalculateNucleotideComposition(input);

            act.Should().NotThrow($"garbage input (len {len}) must never crash composition");
            AssertWellFormed(c, len);
        }
    }

    /// <summary>
    /// BE: a randomly built canonical-only sequence must have its canonical counters
    /// match an independent re-count, total = Length, and GcContent + AtContent = 1.
    /// Cross-checks the counting logic against a simple oracle over many shapes.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Composition_RandomCanonicalSequences_CountsMatchOracle()
    {
        var rng = new Random(424242);

        for (int iteration = 0; iteration < 1000; iteration++)
        {
            int len = rng.Next(1, 300);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = CanonicalBases[rng.Next(CanonicalBases.Length)];
            string seq = new string(chars);

            var c = SequenceStatistics.CalculateNucleotideComposition(seq);

            c.CountA.Should().Be(seq.Count(ch => ch == 'A'));
            c.CountT.Should().Be(seq.Count(ch => ch == 'T'));
            c.CountG.Should().Be(seq.Count(ch => ch == 'G'));
            c.CountC.Should().Be(seq.Count(ch => ch == 'C'));
            c.CountU.Should().Be(seq.Count(ch => ch == 'U'));
            c.CountN.Should().Be(0);
            c.CountOther.Should().Be(0);

            // All-canonical ⇒ total = Length ⇒ fractions sum to 1.
            (c.GcContent + c.AtContent).Should().BeApproximately(1.0, Tolerance);
            AssertWellFormed(c, len);
        }
    }

    #endregion
}
