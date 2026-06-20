using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area simple melting-temperature (Tm) estimate.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (KeyNotFoundException on a non-ACGT base,
/// DivideByZero on an empty / all-junk sequence, NaN / ±Infinity result, …).
/// Every input must yield EITHER a well-defined, theory-correct Tm, OR a
/// *documented, intentional* outcome (the 0 sentinel for empty / no recognized
/// base). A raw runtime exception, a NaN Tm, or a hang on garbage input is a bug,
/// not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-TM-001 — Melting Temperature (Wallace / Marmur-Doty) (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 130.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation —
///       · empty (L=0 ⇒ the documented 0 sentinel, no DivideByZero on N);
///       · single base (L=1 ⇒ Wallace short-oligo rule, the smallest non-empty
///         input, no crash);
///       · all-AT (GC-content lower boundary ⇒ documented LOW Tm = 2·(A+T));
///       · all-GC (GC-content upper boundary ⇒ documented HIGH Tm = 4·(G+C));
///       · non-ACGT (junk / N / lowercase ⇒ documented case-folding + skip of
///         unrecognized bases, no KeyNotFound, Tm well-defined over A/C/G/T).
///     The formula-switch length threshold (Wallace for N &lt; 14, Marmur-Doty
///     otherwise) is itself a boundary: the N = 13 / N = 14 inclusive/exclusive
///     edge is pinned explicitly. DivideByZero on an all-junk Marmur-Doty input
///     (total = 0) and a case bug are the target failure modes.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE = граничні значення).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The melting-temperature contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry (SequenceStatistics.cs, src/.../Seqeron.Genomics.Analysis):
///   • SequenceStatistics.CalculateMeltingTemperature(string, bool useWallaceRule = true)
///         → double Tm (°C)
///
/// Documented behaviour (Melting_Temperature.md, Test Unit ID SEQ-TM-001):
///   • §2.2 / INV-01 (Wallace rule, short oligos, applied when useWallaceRule and
///         N &lt; 14):            Tm = 2·(A + T) + 4·(G + C)
///   • §2.2 / INV-02 (Marmur-Doty GC formula, otherwise):
///                                Tm = 64.9 + 41·(GC − 16.4) / N
///         where N is the recognized base count (A+T+G+C); Marmur-Doty returns 0
///         when that count is 0.
///   • §4.1 step 2: the formula switch is `useWallaceRule && N_length &lt; 14`,
///         where N_length is `dnaSequence.Length` (the threshold constant is
///         ThermoConstants.WallaceMaxLength = 14 — Wallace for length &lt; 14,
///         Marmur-Doty for length ≥ 14).
///   • §3.3 / §6.1: null/empty ⇒ 0; no exception for guarded input.
///   • §6.1: input is upper-cased internally (case-insensitive); unrecognized
///         characters are simply not counted as A/C/G/T (so they contribute 0 to
///         both Wallace counts and the Marmur-Doty N) — no KeyNotFound, no throw.
///   • §6.1 (all-A·T): A·T pairs are least stable ⇒ low Tm; G·C ⇒ high Tm.
///
/// The Wallace/Marmur-Doty constants pinned below are an INDEPENDENT oracle copy
/// of the documented formulae (Melting_Temperature.md §2.2; ThermoConstants), so
/// the tests assert the real BUSINESS formula, not merely non-throwing.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceMeltingTemperatureFuzzTests
{
    #region Helpers — independent oracle (copy of Melting_Temperature.md §2.2)

    private const double Tolerance = 1e-9;

    // Wallace contributions and threshold (Melting_Temperature.md §2.2 INV-01, §4.1).
    private const int WallaceAt = 2;            // 2 per A·T
    private const int WallaceGc = 4;            // 4 per G·C
    private const int WallaceMaxLength = 14;    // Wallace for length < 14

    // Marmur-Doty constants (Melting_Temperature.md §2.2 INV-02).
    private const double MarmurBase = 64.9;
    private const double MarmurGcCoeff = 41.0;
    private const double MarmurGcOffset = 16.4;

    /// <summary>Counts only the recognized A/C/G/T bases (case-insensitive); everything
    /// else (N, U, gaps, unicode, …) is excluded, mirroring CalculateNucleotideComposition.</summary>
    private static (int at, int gc) CountAcgt(string seq)
    {
        int at = 0, gc = 0;
        foreach (char raw in seq)
        {
            switch (char.ToUpperInvariant(raw))
            {
                case 'A': case 'T': at++; break;
                case 'G': case 'C': gc++; break;
                default: break; // U, N, gaps, junk — not recognized, contributes 0
            }
        }
        return (at, gc);
    }

    /// <summary>Independent oracle mirroring CalculateMeltingTemperature exactly:
    /// 0 for null/empty; Wallace when useWallaceRule and length &lt; 14; otherwise
    /// Marmur-Doty over the recognized base count (0 when no A/C/G/T present).</summary>
    private static double Oracle(string seq, bool useWallaceRule = true)
    {
        if (string.IsNullOrEmpty(seq)) return 0;

        var (at, gc) = CountAcgt(seq);

        if (useWallaceRule && seq.Length < WallaceMaxLength)
            return WallaceAt * at + WallaceGc * gc;

        int total = at + gc;
        if (total == 0) return 0;
        return MarmurBase + MarmurGcCoeff * (gc - MarmurGcOffset) / total;
    }

    /// <summary>Universal well-formedness: Tm is always a finite number (never NaN /
    /// ±Infinity) for ANY input — the core §8 fuzzing guarantee.</summary>
    private static void AssertFinite(double tm) =>
        double.IsFinite(tm).Should().BeTrue("Tm must be finite — no NaN/±Infinity (§8 fuzzing contract)");

    private static void ShouldMatchOracle(double tm, string seq, bool wallace = true)
    {
        tm.Should().BeApproximately(Oracle(seq, wallace), Tolerance,
            $"Tm for '{seq}' (wallace={wallace}) follows the documented formula (§2.2)");
        AssertFinite(tm);
    }

    /// <summary>Random A/C/G/T string of the given length (a well-formed oligo).</summary>
    private static string RandomDna(Random rng, int length)
    {
        const string bases = "ACGT";
        var c = new char[length];
        for (int i = 0; i < length; i++) c[i] = bases[rng.Next(bases.Length)];
        return new string(c);
    }

    /// <summary>Random arbitrary BMP code points (control chars, null byte, lone surrogate
    /// halves, unicode) — fuzz fodder, none of which are recognized as A/C/G/T.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++) c[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(c);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-TM-001 — Wallace / Marmur-Doty melting temperature : fuzz (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact Wallace / Marmur-Doty Tm

    /// <summary>
    /// Positive baseline (not a boundary): a known short oligo reproduces the documented Wallace
    /// Tm EXACTLY = 2·(A+T) + 4·(G+C). "ACGT" has 2 A·T and 2 G·C ⇒ 2·2 + 4·2 = 12.0 °C.
    /// Confirms the suite asserts the real BUSINESS formula, not merely non-throwing.
    /// — Melting_Temperature.md §2.2 INV-01.
    /// </summary>
    [Test]
    public void Tm_KnownShortOligo_MatchesWallaceExactly()
    {
        SequenceStatistics.CalculateMeltingTemperature("ACGT")
            .Should().BeApproximately(12.0, Tolerance, "2·(2 A·T) + 4·(2 G·C) = 12.0 (INV-01)");

        // A 12-mer (< 14) stays on the Wallace branch: AAAAAAGGGGGG = 6 A·T + 6 G·C ⇒ 2·6 + 4·6 = 36.
        SequenceStatistics.CalculateMeltingTemperature("AAAAAAGGGGGG")
            .Should().BeApproximately(36.0, Tolerance, "6·2 + 6·4 = 36.0 (INV-01)");
    }

    /// <summary>
    /// Positive baseline: a longer sequence (length ≥ 14) yields the documented Marmur-Doty
    /// GC%-based Tm EXACTLY = 64.9 + 41·(GC − 16.4)/N. A 20-mer with 10 G·C:
    /// 64.9 + 41·(10 − 16.4)/20 = 64.9 − 13.12 = 51.78 °C.
    /// — Melting_Temperature.md §2.2 INV-02.
    /// </summary>
    [Test]
    public void Tm_LongSequence_MatchesMarmurDotyExactly()
    {
        // 10 A·T + 10 G·C, length 20.
        string seq = "ATATATATATGCGCGCGCGC";
        seq.Length.Should().Be(20);

        double expected = MarmurBase + MarmurGcCoeff * (10 - MarmurGcOffset) / 20; // 51.78
        SequenceStatistics.CalculateMeltingTemperature(seq)
            .Should().BeApproximately(expected, Tolerance, "64.9 + 41·(10−16.4)/20 = 51.78 (INV-02)");
        SequenceStatistics.CalculateMeltingTemperature(seq).Should().BeApproximately(51.78, 1e-9);
    }

    /// <summary>
    /// Positive baseline (case-insensitivity): lowercase input yields the IDENTICAL Tm as
    /// upper-case — input is ToUpperInvariant'd before counting. — Melting_Temperature.md §6.1.
    /// </summary>
    [Test]
    public void Tm_LowercaseInput_EqualsUppercase()
    {
        SequenceStatistics.CalculateMeltingTemperature("gcatgcat")
            .Should().Be(SequenceStatistics.CalculateMeltingTemperature("GCATGCAT"),
                "Tm is case-insensitive (§6.1)");

        SequenceStatistics.CalculateMeltingTemperature("acgtacgtacgtacgt")  // 16-mer, Marmur branch
            .Should().Be(SequenceStatistics.CalculateMeltingTemperature("ACGTACGTACGTACGT"));
    }

    #endregion

    #region BE — formula-switch boundary: Wallace (N<14) vs Marmur-Doty (N≥14)

    /// <summary>
    /// BE: the length threshold is the formula-switch boundary. ThermoConstants.WallaceMaxLength
    /// is 14, and the switch is `length &lt; 14` (Wallace) vs `length ≥ 14` (Marmur-Doty). The
    /// edge is asserted on BOTH sides with the SAME composition (all-A), so only the formula —
    /// not the counts — changes:
    ///   • length 13 (all-A) ⇒ Wallace ⇒ 2·13 = 26.0 °C;
    ///   • length 14 (all-A) ⇒ Marmur-Doty ⇒ 64.9 + 41·(0 − 16.4)/14 = 16.8929 °C.
    /// A regression that moved the threshold (≤ 14 / off-by-one) would flip one of these.
    /// — Melting_Temperature.md §4.1, ThermoConstants.WallaceMaxLength.
    /// </summary>
    [Test]
    public void Tm_FormulaSwitch_AtLength14_Boundary()
    {
        ThermoConstants.WallaceMaxLength.Should().Be(14, "documented Wallace threshold (§4.1)");

        string len13 = new string('A', 13);
        string len14 = new string('A', 14);

        SequenceStatistics.CalculateMeltingTemperature(len13)
            .Should().BeApproximately(2.0 * 13, Tolerance, "length 13 < 14 ⇒ Wallace 2·(A+T) = 26.0");

        double marmur14 = MarmurBase + MarmurGcCoeff * (0 - MarmurGcOffset) / 14; // 16.8928...
        SequenceStatistics.CalculateMeltingTemperature(len14)
            .Should().BeApproximately(marmur14, Tolerance, "length 14 ≥ 14 ⇒ Marmur-Doty (formula switch)");

        // Both branches agree with the oracle and differ (the switch really fired).
        ShouldMatchOracle(SequenceStatistics.CalculateMeltingTemperature(len13), len13);
        ShouldMatchOracle(SequenceStatistics.CalculateMeltingTemperature(len14), len14);
        SequenceStatistics.CalculateMeltingTemperature(len13)
            .Should().NotBe(SequenceStatistics.CalculateMeltingTemperature(len14));
    }

    /// <summary>
    /// BE: with useWallaceRule = false the Wallace branch is suppressed even below the threshold,
    /// so a short oligo takes the Marmur-Doty formula. "ACGT" (length 4, 2 G·C) ⇒
    /// 64.9 + 41·(2 − 16.4)/4 = −82.7 °C — a documented (and physically meaningless for so short
    /// an oligo, §2.3 ASM-03) value that must still be finite and oracle-exact, NOT the Wallace 12.
    /// — Melting_Temperature.md §3.1 (useWallaceRule), §4.1.
    /// </summary>
    [Test]
    public void Tm_UseWallaceRuleFalse_ForcesMarmurDoty_BelowThreshold()
    {
        double wallace = SequenceStatistics.CalculateMeltingTemperature("ACGT", useWallaceRule: true);
        double marmur = SequenceStatistics.CalculateMeltingTemperature("ACGT", useWallaceRule: false);

        wallace.Should().BeApproximately(12.0, Tolerance, "Wallace branch");
        marmur.Should().BeApproximately(MarmurBase + MarmurGcCoeff * (2 - MarmurGcOffset) / 4, Tolerance,
            "useWallaceRule=false forces Marmur-Doty even at length 4");
        marmur.Should().NotBe(wallace);
        ShouldMatchOracle(marmur, "ACGT", wallace: false);
    }

    #endregion

    #region BE — Boundary: empty / null (length 0 ⇒ 0 sentinel, no DivideByZero)

    /// <summary>
    /// BE: the empty string is the L=0 lower size boundary — N = 0 would divide by zero in
    /// Marmur-Doty, so the documented result is the 0 sentinel reached via the IsNullOrEmpty guard
    /// BEFORE any division. Asserted for both Wallace and Marmur branches.
    /// — Melting_Temperature.md §3.3 / §6.1.
    /// </summary>
    [Test]
    public void Tm_EmptyString_IsZero_NoDivideByZero()
    {
        foreach (bool wallace in new[] { true, false })
        {
            var act = () => SequenceStatistics.CalculateMeltingTemperature(string.Empty, wallace);
            act.Should().NotThrow("L=0 is guarded before any N-division (§3.3)");
            act().Should().Be(0.0, "empty ⇒ documented 0 sentinel");
            AssertFinite(act());
        }
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — 0, no
    /// NullReferenceException. — Melting_Temperature.md §3.3.
    /// </summary>
    [Test]
    public void Tm_Null_IsZero_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateMeltingTemperature(null!);
        act.Should().NotThrow("null is 'no sequence', not an error");
        act().Should().Be(0.0);
    }

    #endregion

    #region BE — Boundary: single base (length 1 ⇒ Wallace short-oligo rule)

    /// <summary>
    /// BE: a single base is the smallest non-empty input (length 1 &lt; 14 ⇒ Wallace). Each
    /// recognized base contributes its Wallace weight: A/T ⇒ 2.0 °C, G/C ⇒ 4.0 °C — exactly,
    /// no crash, upper or lower case. — Melting_Temperature.md §2.2 INV-01.
    /// </summary>
    [TestCase("A", 2.0)]
    [TestCase("T", 2.0)]
    [TestCase("G", 4.0)]
    [TestCase("C", 4.0)]
    [TestCase("a", 2.0)]
    [TestCase("g", 4.0)]
    public void Tm_SingleBase_UsesWallaceShortRule(string seq, double expected)
    {
        var act = () => SequenceStatistics.CalculateMeltingTemperature(seq);
        act.Should().NotThrow($"length-1 '{seq}' is the smallest oligo, Wallace branch");
        act().Should().BeApproximately(expected, Tolerance, $"Wallace single-base '{seq}' = {expected}");
        ShouldMatchOracle(act(), seq);
    }

    /// <summary>
    /// BE: a single UNRECOGNIZED character is the length-1 boundary where the only base is also
    /// out of the alphabet. Wallace over zero recognized bases ⇒ 2·0 + 4·0 = 0.0, no KeyNotFound,
    /// no crash. — Melting_Temperature.md §6.1.
    /// </summary>
    [TestCase("N")]
    [TestCase("U")]
    [TestCase("-")]
    [TestCase("*")]
    [TestCase("7")]
    public void Tm_SingleUnrecognizedChar_IsZero_NoThrow(string seq)
    {
        var act = () => SequenceStatistics.CalculateMeltingTemperature(seq);
        act.Should().NotThrow($"'{seq}' is unrecognized, contributes 0 — no KeyNotFound");
        act().Should().Be(0.0, "no A/C/G/T ⇒ Wallace 0");
    }

    #endregion

    #region BE — Boundary: all-AT (GC-content lower bound ⇒ LOW Tm)

    /// <summary>
    /// BE: an all-A·T oligo is the GC-content lower boundary — the least-stable extreme ⇒ low Tm.
    /// Short forms take Wallace = 2·(A+T) with NO G·C term:
    ///   • "ATAT"  ⇒ 2·4 = 8.0;  "AAAA" ⇒ 8.0;  "TTTTTTTT" (8) ⇒ 16.0.
    /// All match the oracle and are strictly below the equal-length all-G·C Tm (cross-boundary
    /// contrast asserted separately). — Melting_Temperature.md §6.1 (A·T least stable), §2.2.
    /// </summary>
    [Test]
    public void Tm_AllAT_IsLow_MatchesWallaceATOnly()
    {
        SequenceStatistics.CalculateMeltingTemperature("ATAT")
            .Should().BeApproximately(8.0, Tolerance, "2·(4 A·T) + 4·0 = 8.0");
        SequenceStatistics.CalculateMeltingTemperature("AAAA")
            .Should().BeApproximately(8.0, Tolerance, "2·4 = 8.0");

        foreach (string seq in new[] { "AT", "TA", "AAAAAAAA", "ATATATAT", "TATATA", "AATTAATT" })
            ShouldMatchOracle(SequenceStatistics.CalculateMeltingTemperature(seq), seq);

        // A long all-A·T sequence (≥14) takes Marmur-Doty with GC = 0 ⇒ 64.9 + 41·(−16.4)/N.
        string longAt = new string('A', 30);
        SequenceStatistics.CalculateMeltingTemperature(longAt)
            .Should().BeApproximately(MarmurBase + MarmurGcCoeff * (0 - MarmurGcOffset) / 30, Tolerance,
                "all-A·T Marmur-Doty, GC term = 0");
    }

    #endregion

    #region BE — Boundary: all-GC (GC-content upper bound ⇒ HIGH Tm)

    /// <summary>
    /// BE: an all-G·C oligo is the GC-content upper boundary — the most-stable extreme ⇒ high Tm.
    /// Short forms take Wallace = 4·(G+C) with NO A·T term:
    ///   • "GCGC" ⇒ 4·4 = 16.0;  "GGGG" ⇒ 16.0;  "CCCCCCCC" (8) ⇒ 32.0.
    /// — Melting_Temperature.md §2.2 INV-01.
    /// </summary>
    [Test]
    public void Tm_AllGC_IsHigh_MatchesWallaceGCOnly()
    {
        SequenceStatistics.CalculateMeltingTemperature("GCGC")
            .Should().BeApproximately(16.0, Tolerance, "4·(4 G·C) = 16.0");
        SequenceStatistics.CalculateMeltingTemperature("GGGG")
            .Should().BeApproximately(16.0, Tolerance, "4·4 = 16.0");

        foreach (string seq in new[] { "GC", "CG", "GGGGGGGG", "GCGCGCGC", "CGCGCG", "GGCCGGCC" })
            ShouldMatchOracle(SequenceStatistics.CalculateMeltingTemperature(seq), seq);
    }

    /// <summary>
    /// BE (the headline physical contrast): for an EQUAL-length oligo, all-G·C melts HIGHER than
    /// all-A·T — strictly higher Tm in BOTH formula regimes. This is the cross-boundary invariant
    /// the two GC-content extremes exist to expose; a sign flip here would mean the Wallace 2/4
    /// weights or the Marmur GC term were swapped. — Melting_Temperature.md §2.1, §6.1.
    /// </summary>
    [Test]
    public void Tm_AllGC_HigherThan_AllAT_OfEqualLength()
    {
        foreach (int n in new[] { 1, 4, 8, 13, 14, 20, 50 })
        {
            string gc = new string('G', n);
            string at = new string('A', n);

            SequenceStatistics.CalculateMeltingTemperature(gc)
                .Should().BeGreaterThan(SequenceStatistics.CalculateMeltingTemperature(at),
                    $"all-G·C (n={n}) melts higher than all-A·T of equal length");
        }
    }

    #endregion

    #region BE / MC — non-ACGT: skipped, never KeyNotFound, Tm over recognized bases

    /// <summary>
    /// BE/contract: unrecognized characters (N, U, gaps, junk) are not counted as A/C/G/T and add
    /// 0 — no KeyNotFound, no throw. Interior junk in a short oligo therefore leaves the Wallace Tm
    /// equal to that of its recognized bases alone:
    ///   • "AC-GT" (length 5 &lt; 14) ⇒ A,C,G,T recognized, '-' skipped ⇒ 2·2 + 4·2 = 12.0;
    ///   • "NNAANN" ⇒ 2 A·T only ⇒ 4.0;  "GCNNGC" ⇒ 4 G·C ⇒ 16.0.
    /// All match the oracle. — Melting_Temperature.md §6.1.
    /// </summary>
    [TestCase("AC-GT", 12.0)]
    [TestCase("NNAANN", 4.0)]
    [TestCase("GCNNGC", 16.0)]
    [TestCase("A C G T", 12.0)]
    public void Tm_NonAcgt_Skipped_NoKeyNotFound(string seq, double expected)
    {
        var act = () => SequenceStatistics.CalculateMeltingTemperature(seq);
        act.Should().NotThrow($"'{seq}' contains non-ACGT — must be skipped, not throw");
        act().Should().BeApproximately(expected, Tolerance,
            $"only recognized A/C/G/T contribute (§6.1) for '{seq}'");
        ShouldMatchOracle(act(), seq);
    }

    /// <summary>
    /// BE: an all-junk sequence long enough (length ≥ 14) to take the Marmur-Doty branch has a
    /// recognized base count of 0 — the documented `total == 0 ⇒ 0` guard prevents a DivideByZero
    /// and yields 0, NOT NaN. This is the critical empty-recognized-set boundary distinct from the
    /// empty-string boundary. — Melting_Temperature.md §2.2 INV-02, ThermoConstants.CalculateMarmurDotyTm.
    /// </summary>
    [TestCase("NNNNNNNNNNNNNN")]   // 14 × N
    [TestCase("--------------------")] // 20 × gap
    [TestCase("UUUUUUUUUUUUUUUU")] // 16 × U (RNA base, not recognized as DNA)
    public void Tm_AllJunkLong_MarmurDoty_NoDivideByZero(string seq)
    {
        seq.Length.Should().BeGreaterThanOrEqualTo(WallaceMaxLength, "long enough for the Marmur branch");

        var act = () => SequenceStatistics.CalculateMeltingTemperature(seq);
        act.Should().NotThrow($"all-junk '{seq}': recognized count 0 ⇒ guarded, no DivideByZero");
        act().Should().Be(0.0, "Marmur-Doty total == 0 ⇒ documented 0 (INV-02 guard)");
        AssertFinite(act());
    }

    #endregion

    #region BE / MC — random garbage and random DNA: never throws, finite, oracle-exact

    /// <summary>
    /// MC/BE: a large batch of arbitrary BMP strings (control chars, the null byte, lone surrogate
    /// halves, unicode, occasionally seeded with real bases) of length 0..200 must NEVER throw and
    /// ALWAYS yield a finite Tm matching the independent oracle, on BOTH the Wallace and Marmur
    /// branches. Core fuzz guarantee: no KeyNotFound on a non-ACGT char, no DivideByZero, no NaN.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Tm_RandomGarbageStrings_NeverThrow_MatchOracle()
    {
        var rng = new Random(130_001);

        for (int iteration = 0; iteration < 5000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);
            bool wallace = rng.Next(2) == 0;

            double tm = 0;
            var act = () => tm = SequenceStatistics.CalculateMeltingTemperature(input, wallace);
            act.Should().NotThrow($"garbage (len {len}) must never crash Tm");

            ShouldMatchOracle(tm, input, wallace);
        }
    }

    /// <summary>
    /// MC/BE: randomly built genuine DNA oligos (A/C/G/T, length 1..300, varied GC content) must
    /// equal the independent oracle EXACTLY across the full random space and across the formula
    /// switch (Wallace below 14, Marmur-Doty at/above), with finite Tm. Confirms the base counting
    /// and the length-threshold branch are correct on both sides of the boundary.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Tm_RandomDnaOligos_MatchOracle()
    {
        var rng = new Random(130_002);

        for (int iteration = 0; iteration < 5000; iteration++)
        {
            int len = rng.Next(1, 300);
            string seq = RandomDna(rng, len);
            bool wallace = rng.Next(2) == 0;

            double tm = SequenceStatistics.CalculateMeltingTemperature(seq, wallace);
            ShouldMatchOracle(tm, seq, wallace);
        }
    }

    /// <summary>
    /// BE: a very long oligo must remain finite and terminate quickly (O(n) base counting) — no
    /// overflow, no hang. The Marmur-Doty closed form is asserted against the oracle up to 1e6.
    /// — Melting_Temperature.md §4.3 (O(n), O(1)).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Tm_VeryLongOligo_IsFinite_MatchesOracle()
    {
        foreach (int n in new[] { 1000, 100_000, 1_000_000 })
        {
            string seq = new string('G', n);
            double tm = SequenceStatistics.CalculateMeltingTemperature(seq);
            AssertFinite(tm);
            tm.Should().BeApproximately(Oracle(seq), 1e-6, $"poly-G n={n} matches Marmur-Doty oracle");
        }
    }

    #endregion
}
