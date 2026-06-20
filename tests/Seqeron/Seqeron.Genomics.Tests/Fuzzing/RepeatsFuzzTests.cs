using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Repeats area — short tandem repeat (microsatellite / STR)
/// detection (REP-STR-001) and general tandem repeat detection (REP-TANDEM-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (IndexOutOfRangeException, ArgumentOutOfRangeException
/// leaking from internal indexing, DivideByZero, OutOfMemory). Every input must
/// resolve to EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* validation exception (ArgumentException / ArgumentNullException /
/// ArgumentOutOfRangeException). A raw runtime exception, a hang, or a blow-up of
/// spurious "repeats" on garbage parameters is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-STR-001 — short tandem repeat (microsatellite / STR) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 13.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate parameter boundaries called out
///          in the checklist row: minRepeats = 0, minRepeats = 1, maxUnitLen >
///          seqLen, and the empty sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The STR-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A microsatellite is a motif U of length m (1 ≤ m ≤ 6 by default) repeated
/// consecutively k times: S[p .. p + k·m) = Uᵏ. Detection enumerates unit lengths
/// minUnitLength..maxUnitLength, skips redundant motifs (motifs that are themselves
/// a repetition of a smaller motif, e.g. "ATAT"), counts consecutive copies, and
/// emits a MicrosatelliteResult only when the count reaches minRepeats.
///   — docs/algorithms/Repeat_Analysis/Microsatellite_Detection.md §2.2, §2.4.
///
/// API entry: RepeatFinder.FindMicrosatellites(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs).
///
/// Documented parameter contract (Microsatellite_Detection.md §3.1, §3.3, §6.1;
/// RepeatFinder.cs lines 25–99) — the SAME validation gate is wired on ALL four
/// overloads (typed DnaSequence, raw string, and both cancellable variants):
///   • minUnitLength &lt; 1            → ArgumentOutOfRangeException;
///   • maxUnitLength &lt; minUnitLength → ArgumentOutOfRangeException;
///   • minRepeats &lt; 2              → ArgumentOutOfRangeException. This is THE
///     boundary the checklist row probes: a *literal* minRepeats = 0 would make
///     every single position a "repeat" (k = 1 ≥ 0) and explode the output, and
///     minRepeats = 1 would report every substring as a trivial 1-copy repeat —
///     so the contract REJECTS both with the documented exception rather than
///     emitting nonsense. The accepted floor is minRepeats = 2
///     (Microsatellite_Detection.md §6.1). We pin that 0 and 1 throw, and that 2
///     is accepted, so the floor cannot silently drift.
///
/// Two documented surfaces with DIFFERENT null/empty handling — fuzzing pins both
/// and the boundary between them (Microsatellite_Detection.md §3.1, §3.3):
///   (1) The TYPED overload — FindMicrosatellites(DnaSequence, …)
///       (RepeatFinder.cs lines 25–37): a null DnaSequence → ArgumentNullException
///       (explicit ThrowIfNull guard, never a NullReferenceException). An empty
///       DnaSequence (the ctor's IsNullOrEmpty short-circuit materializes one) →
///       the empty enumerable: the scan loop bound `i ≤ len − unitLen·minRepeats`
///       is negative, so the loop never runs (Microsatellite_Detection.md §6.1).
///   (2) The RAW-STRING overload — FindMicrosatellites(string, …)
///       (RepeatFinder.cs lines 65–80): null OR empty input → the empty enumerable
///       via `yield break` — NOT an exception. Non-empty input is upper-cased
///       before scanning. Parameter validation (the three OutOfRange checks) runs
///       BEFORE the null/empty short-circuit, so a bad parameter throws even on a
///       null/empty sequence.
///
/// The maxUnitLength &gt; seqLen boundary: when the largest searched unit cannot
/// even fit minRepeats copies inside the sequence, the inner loop bound
/// `i ≤ len − unitLen·minRepeats` is negative for that unit length, so the loop
/// body never executes — no Substring out-of-range, no crash, an empty (or shorter)
/// result. We pin that an oversized maxUnitLength can never index past the end.
///
/// Note these overloads return a LAZY IEnumerable (iterator with `yield`); the
/// OutOfRange parameter validation on the string overloads is itself inside the
/// iterator body, so it only fires on enumeration. Every test therefore forces
/// enumeration (`.ToList()`) so the documented exception actually surfaces and any
/// hang would manifest as a non-terminating materialization.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-TANDEM-001 — general tandem repeat detection
/// Checklist: docs/checklists/03_FUZZING.md, row 14.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: minRepetitions = 0, minUnitLength = 0, a unit length
///          ceiling of 1 (only homopolymer runs count), the empty sequence, and a
///          single-character sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The tandem-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A tandem repeat is a region S[p .. p + k·|U|) = Uᵏ where a unit U of length m
/// occurs consecutively k ≥ 2 times with no gap (Tandem_Repeat_Detection.md §2.1,
/// §2.2). The CANONICAL exact detector is
///   GenomicAnalyzer.FindTandemRepeats(DnaSequence sequence,
///                                     int minUnitLength = 2,
///                                     int minRepetitions = 2)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs).
/// This is the surface the checklist row probes (minReps / minUnitLen are its
/// parameters); it is also the primary surface named in the differential checklist
/// (08_DIFFERENTIAL_TESTING.md row 14: GenomicAnalyzer vs RepeatFinder). The sibling
/// aggregator RepeatFinder.GetTandemRepeatSummary delegates to FindMicrosatellites
/// and is not the parameterised detector under fuzz here.
///
/// Documented parameter contract — BEFORE this unit's work the doc recorded that
/// FindTandemRepeats performed NO validation and that "callers can trigger
/// undefined or exception-driven behavior with nonsensical values" was an *accepted*
/// assumption (Tandem_Repeat_Detection.md §3.3, §5.4). Fuzzing proved that the two
/// degenerate boundaries the checklist row targets were UNDISCIPLINED failures, not
/// merely "exception-driven":
///   • minRepetitions = 0 → the unit-length bound `unitLen ≤ seq.Length / minReps`
///     divides by zero → a raw DivideByZeroException (a crash, not a documented
///     ArgumentException). Per the fuzzing doctrine (ADVANCED §8) a raw runtime
///     exception is a bug.
///   • minUnitLength = 0 → a zero-length unit: `unit = Substring(start, 0) = ""`,
///     and the extension `while (… && Substring(pos, 0) == "")` is ALWAYS true while
///     `pos += 0` never advances → a NON-TERMINATING loop (a hang). A hang is the
///     single worst fuzzing outcome.
/// Both were fixed at the source by adding the same validation gate the sibling
/// repeat finders already use (FindMicrosatellites / FindPalindromes):
///   • sequence == null            → ArgumentNullException (ThrowIfNull);
///   • minUnitLength &lt; 1         → ArgumentOutOfRangeException (a 0-length unit
///     would hang; negative is nonsense);
///   • minRepetitions &lt; 2        → ArgumentOutOfRangeException (k ≥ 2 is the
///     definition of a tandem repeat; 0 would divide-by-zero, 1 would mark every
///     substring a trivial 1-copy "repeat").
/// The validation is hoisted into an eager wrapper (FindTandemRepeatsCore holds the
/// `yield` body) so the exception surfaces at the call, not only on enumeration.
/// The tests below PIN these contracts so the floors cannot silently drift, and
/// pin that the homopolymer-only (unitLen ceiling 1), empty, and single-char
/// boundaries return clean empty/correct results with no crash and no hang.
///
/// Documented invariants pinned on every positive result (Tandem_Repeat_Detection.md
/// §2.4): INV-01 Repetitions ≥ minRepetitions and |Unit| ≥ minUnitLength;
/// INV-02 TotalLength = |Unit| × Repetitions; INV-03 Position + |Unit|×Repetitions
/// ≤ sequence.Length.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RepeatsFuzzTests
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
    //  REP-STR-001 — short tandem repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-STR-001 — short tandem repeat detection

    #region BE — Boundary: minRepeats = 0

    /// <summary>
    /// BE: minRepeats = 0 is the degenerate floor. A literal 0 would make k = 1 ≥ 0
    /// true at every position, so EVERY base would be a "microsatellite" — a
    /// nonsense blow-up of the result set. The documented contract REJECTS it with
    /// ArgumentOutOfRangeException (minRepeats &lt; 2 → throw, RepeatFinder.cs lines
    /// 34/73/95; Microsatellite_Detection.md §3.3, §6.1) on BOTH the typed and the
    /// raw-string surface, BEFORE producing any result. The string overload's check
    /// lives inside the iterator, so we force enumeration to surface it.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsZero_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindMicrosatellites(new DnaSequence("ACACACAC"), 1, 6, 0).ToList();
        var raw = () => RepeatFinder.FindMicrosatellites("ACACACAC", 1, 6, 0).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
            "minRepeats = 0 is below the documented floor of 2; a literal 0 would mark every position a repeat");
        raw.Should().Throw<ArgumentOutOfRangeException>(
            "the raw-string overload enforces the same minRepeats >= 2 floor before scanning");
    }

    #endregion

    #region BE — Boundary: minRepeats = 1

    /// <summary>
    /// BE: minRepeats = 1 is the trivial-repeat boundary. With minRepeats = 1 every
    /// substring is a 1-copy "repeat", which would produce a defined-but-useless
    /// quadratic blow-up. The contract REJECTS it with ArgumentOutOfRangeException
    /// (minRepeats &lt; 2 → throw, RepeatFinder.cs lines 34/73/95;
    /// Microsatellite_Detection.md §6.1) on both surfaces — there is no hang and no
    /// blow-up because the value never reaches the scan. The accepted floor is
    /// exactly 2.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsOne_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindMicrosatellites(new DnaSequence("ACACACAC"), 1, 6, 1).ToList();
        var raw = () => RepeatFinder.FindMicrosatellites("ACACACAC", 1, 6, 1).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
            "minRepeats = 1 is below the documented floor of 2; every substring would be a trivial 1-copy repeat");
        raw.Should().Throw<ArgumentOutOfRangeException>(
            "the raw-string overload enforces the same minRepeats >= 2 floor before scanning");
    }

    /// <summary>
    /// BE: the accepted floor minRepeats = 2 must NOT throw and must behave —
    /// pinning that the rejection boundary is exactly at &lt; 2, not at ≤ 2. A clean
    /// di-nucleotide with two copies ("ATAT" is redundant as a unit, so we use a
    /// non-redundant 2-mer repeated twice: "ATAT" → use "AC" twice = "ACAC") is
    /// reported with count 2.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsTwo_IsAcceptedAndDetects()
    {
        var act = () => RepeatFinder.FindMicrosatellites("ACAC", 2, 2, 2).ToList();
        act.Should().NotThrow("minRepeats = 2 is the accepted floor, not below it");

        var results = RepeatFinder.FindMicrosatellites("ACAC", 2, 2, 2).ToList();

        results.Should().ContainSingle(
            "exactly one non-redundant 2-mer 'AC' repeats twice across 'ACAC'");
        results[0].RepeatUnit.Should().Be("AC");
        results[0].RepeatCount.Should().Be(2);
    }

    #endregion

    #region BE — Boundary: maxUnitLength > sequence length

    /// <summary>
    /// BE: maxUnitLength far larger than the sequence length must NOT crash. For a
    /// unit length that cannot fit minRepeats copies in the sequence, the scan
    /// bound `i ≤ len − unitLen·minRepeats` is negative, so the inner loop never
    /// runs and no Substring is taken past the end (RepeatFinder.cs lines 184–195).
    /// The result is simply the repeats found by the unit lengths that DO fit —
    /// here, the short sequence "AC" has no room for any 2+ copy repeat, so the
    /// result is empty. Pinned on both the typed and raw surfaces.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MaxUnitLengthExceedsSequence_IsEmptyAndDoesNotThrow()
    {
        // 2-char sequence, maxUnitLength = 1000 ≫ length. No unit length 1..1000 can
        // fit ≥ 2 copies in 2 chars except a 1-mer homopolymer, which "AC" is not.
        var typedAct = () => RepeatFinder.FindMicrosatellites(new DnaSequence("AC"), 1, 1000, 2).ToList();
        var rawAct = () => RepeatFinder.FindMicrosatellites("AC", 1, 1000, 2).ToList();

        typedAct.Should().NotThrow(
            "an oversized maxUnitLength makes the scan bound negative for the big units; the loop simply never runs");
        rawAct.Should().NotThrow(
            "the raw-string overload is equally guarded against indexing past the sequence end");

        RepeatFinder.FindMicrosatellites(new DnaSequence("AC"), 1, 1000, 2).Should().BeEmpty(
            "no unit can repeat 2× inside a 2-char heteropolymer; the oversized max yields nothing, not a crash");
        RepeatFinder.FindMicrosatellites("AC", 1, 1000, 2).Should().BeEmpty();
    }

    /// <summary>
    /// BE: maxUnitLength &gt; seqLen on a sequence that DOES contain a short repeat
    /// must still find that short repeat and must not be derailed by the oversized
    /// upper bound. "ACACAC" has a 'AC' di-repeat (3 copies); a maxUnitLength of 100
    /// (far past the 6-char length) leaves that detection intact while the unit
    /// lengths > 2 simply contribute nothing.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MaxUnitLengthExceedsSequence_StillFindsShortRepeat()
    {
        var results = RepeatFinder.FindMicrosatellites("ACACAC", 1, 100, 2).ToList();

        results.Should().Contain(r => r.RepeatUnit == "AC" && r.RepeatCount == 3,
            "the 'AC' di-repeat (3 copies) is found despite maxUnitLength ≫ sequence length");
        results.Should().OnlyContain(r => r.RepeatCount >= 2,
            "INV-01: every reported repeat satisfies RepeatCount >= minRepeats, even at the boundary");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable via `yield break`
    /// (RepeatFinder.cs lines 75–76) — no exception. The typed overload over an
    /// empty DnaSequence (the ctor materializes an empty sequence from "" via its
    /// IsNullOrEmpty short-circuit) yields nothing because the scan bound is
    /// negative. Neither path divides, indexes, or hangs on empty input
    /// (Microsatellite_Detection.md §6.1).
    /// </summary>
    [Test]
    public void FindMicrosatellites_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindMicrosatellites(new DnaSequence(string.Empty), 1, 6, 3).ToList();
        var rawEmptyAct = () => RepeatFinder.FindMicrosatellites(string.Empty, 1, 6, 3).ToList();
        var rawNullAct = () => RepeatFinder.FindMicrosatellites((string)null!, 1, 6, 3).ToList();

        typedAct.Should().NotThrow("an empty sequence has no region long enough to hold a repeat; it yields nothing");
        rawEmptyAct.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNullAct.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RepeatFinder.FindMicrosatellites(new DnaSequence(string.Empty), 1, 6, 3).Should().BeEmpty();
        RepeatFinder.FindMicrosatellites(string.Empty, 1, 6, 3).Should().BeEmpty();
        RepeatFinder.FindMicrosatellites((string)null!, 1, 6, 3).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed
    /// overload guards it with an explicit ArgumentNullException
    /// (ArgumentNullException.ThrowIfNull, RepeatFinder.cs line 31) — never a
    /// NullReferenceException. The exception is raised eagerly at the call (it is
    /// outside the iterator body), so no enumeration is required to surface it.
    /// </summary>
    [Test]
    public void FindMicrosatellites_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindMicrosatellites((DnaSequence)null!, 1, 6, 3);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region Positive sanity — a clear microsatellite is detected correctly

    /// <summary>
    /// Positive sanity: alongside the boundary probes, a textbook microsatellite
    /// must be detected with the CORRECT unit and copy count — fuzzing the degenerate
    /// edges must not come at the cost of the core function silently breaking.
    /// "CACACACA" is a di-nucleotide STR: unit 'CA' repeated 4 times. Pinned per
    /// INV-01..INV-04 (Microsatellite_Detection.md §2.4): RepeatCount = 4,
    /// unit length 2 → Dinucleotide, TotalLength = 2 × 4 = 8.
    /// </summary>
    [Test]
    public void FindMicrosatellites_DinucleotideStr_DetectedWithCorrectUnitAndCount()
    {
        var results = RepeatFinder.FindMicrosatellites("CACACACA", 1, 6, 3).ToList();

        var ca = results.Should().ContainSingle(r => r.RepeatUnit == "CA").Subject;
        ca.RepeatCount.Should().Be(4, "'CA' tiles 'CACACACA' exactly four times");
        ca.RepeatType.Should().Be(RepeatType.Dinucleotide, "a 2-bp unit classifies as a dinucleotide repeat");
        ca.TotalLength.Should().Be(8, "INV-03: TotalLength = unitLength (2) × RepeatCount (4)");
        ca.Position.Should().Be(0, "the repeat tract starts at the first base");
    }

    /// <summary>
    /// Positive sanity: a trinucleotide STR — the clinically important class
    /// (Huntington CAG, Fragile-X CGG, …). "(GAT)5" = "GATGATGATGATGAT" is the
    /// motif 'GAT' repeated 5 times. Pinned: unit 'GAT', count 5, Trinucleotide,
    /// TotalLength = 3 × 5 = 15 (Microsatellite_Detection.md §2.1, §2.4).
    /// </summary>
    [Test]
    public void FindMicrosatellites_TrinucleotideStr_DetectedWithCorrectUnitAndCount()
    {
        var results = RepeatFinder.FindMicrosatellites("GATGATGATGATGAT", 1, 6, 3).ToList();

        var gat = results.Should().ContainSingle(r => r.RepeatUnit == "GAT").Subject;
        gat.RepeatCount.Should().Be(5, "'GAT' tiles '(GAT)5' exactly five times");
        gat.RepeatType.Should().Be(RepeatType.Trinucleotide, "a 3-bp unit classifies as a trinucleotide repeat");
        gat.TotalLength.Should().Be(15, "INV-03: TotalLength = unitLength (3) × RepeatCount (5)");
    }

    /// <summary>
    /// Positive sanity at scale: a long but modest STR tract must complete promptly
    /// (no hang / no quadratic blow-up at the accepted floor) and report a count
    /// consistent with the tract length. 'CA' repeated 500 times → a single
    /// dinucleotide result with RepeatCount 500. Bounds kept modest to guard against
    /// hang while still exercising the extension loop.
    /// </summary>
    [Test]
    public void FindMicrosatellites_LongHomogeneousTract_CompletesAndCountsCorrectly()
    {
        const int copies = 500;
        string tract = string.Concat(Enumerable.Repeat("CA", copies)); // 1000 bases

        var results = RepeatFinder.FindMicrosatellites(tract, 2, 2, 3).ToList();

        var ca = results.Should().ContainSingle(r => r.RepeatUnit == "CA").Subject;
        ca.RepeatCount.Should().Be(copies, "the extension loop counts every consecutive 'CA' copy without overflow or hang");
        ca.TotalLength.Should().Be(copies * 2);
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must produce only
    /// well-formed results and must complete — no spurious negative counts, no
    /// out-of-range positions, no hang. Every result must satisfy the documented
    /// invariants INV-01..INV-03 (Microsatellite_Detection.md §2.4) regardless of
    /// what random content appears.
    /// </summary>
    [Test]
    public void FindMicrosatellites_RandomSequence_ProducesOnlyWellFormedResults()
    {
        string seq = RandomDna(2000, seed: 13_001);

        var results = RepeatFinder.FindMicrosatellites(seq, 1, 6, 3).ToList();

        results.Should().OnlyContain(r =>
            r.RepeatCount >= 3 &&                                   // INV-01
            r.RepeatUnit.Length >= 1 && r.RepeatUnit.Length <= 6 && // INV-02
            r.TotalLength == r.RepeatUnit.Length * r.RepeatCount &&  // INV-03
            r.Position >= 0 && r.Position + r.TotalLength <= seq.Length,
            "every result on random input is well-formed: count ≥ minRepeats, unit in range, total = unit×count, in bounds");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-TANDEM-001 — general tandem repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-TANDEM-001 — tandem repeat detection

    #region BE — Boundary: minRepetitions = 0

    /// <summary>
    /// BE: minRepetitions = 0 is the degenerate floor and was the KEY crash target.
    /// In the pre-fix detector the unit-length bound `unitLen ≤ seq.Length / minReps`
    /// divides by zero → a raw DivideByZeroException (a crash, not a documented
    /// validation error). A tandem repeat is defined as k ≥ 2 copies
    /// (Tandem_Repeat_Detection.md §2.2), so the contract now REJECTS 0 with
    /// ArgumentOutOfRangeException at the call site (eager validation in the wrapper),
    /// never dividing by zero. We force enumeration so a regression to lazy/late
    /// validation would still be caught.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinRepetitionsZero_ThrowsArgumentOutOfRange()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), 2, 0).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "minRepetitions = 0 would divide by zero in the unit-length bound; the contract rejects it")
            .Which.ParamName.Should().Be("minRepetitions");
    }

    #endregion

    #region BE — Boundary: minUnitLength = 0

    /// <summary>
    /// BE: minUnitLength = 0 is the degenerate zero-length-unit boundary and was the
    /// KEY HANG target. A 0-length unit makes `unit = Substring(start, 0) = ""`, and
    /// the extension loop `while (… && Substring(pos, 0) == "")` is always true while
    /// `pos += 0` never advances — a non-terminating loop (the worst fuzzing outcome).
    /// The contract now REJECTS minUnitLength &lt; 1 with ArgumentOutOfRangeException
    /// so the scan can never enter that infinite loop. Pinned eagerly; the guarded
    /// call must return promptly, not hang.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_MinUnitLengthZero_ThrowsArgumentOutOfRange()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), 0, 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a zero-length unit tiles every position infinitely and never advances the scan; the contract rejects it")
            .Which.ParamName.Should().Be("minUnitLength");
    }

    #endregion

    #region BE — Boundary: unit-length ceiling of 1 (homopolymer-only)

    /// <summary>
    /// BE: a unit-length ceiling of 1 — i.e. searching ONLY length-1 units
    /// (minUnitLength = 1) — means only homopolymer runs (AAAA…, the mononucleotide
    /// class) can ever be a tandem repeat; no di-/tri-/longer unit is even considered.
    /// "AAATGC" has a single 'A' run of 3 and nothing else repeats, so exactly one
    /// result — unit "A", 3 copies at position 0 — must be reported, and every result
    /// must be a 1-bp unit. This pins that the smallest legal unit length behaves and
    /// does not over- or under-report (no crash, no spurious multi-bp units).
    /// </summary>
    [Test]
    public void FindTandemRepeats_UnitLengthOne_FindsOnlyHomopolymerRuns()
    {
        var results = GenomicAnalyzer.FindTandemRepeats(new DnaSequence("AAATGC"), minUnitLength: 1, minRepetitions: 2).ToList();

        results.Should().ContainSingle("only the 'A' homopolymer run repeats ≥ 2× in 'AAATGC'");
        results[0].Unit.Should().Be("A");
        results[0].Repetitions.Should().Be(3, "the run 'AAA' is three consecutive 'A' copies");
        results[0].Position.Should().Be(0);
        results.Should().OnlyContain(r => r.Unit.Length == 1,
            "with a unit-length ceiling of 1 only mononucleotide (homopolymer) units are searched");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. With valid thresholds the
    /// scan bounds `unitLen ≤ 0/minReps = 0` and `start ≤ 0 − unitLen·minReps` are
    /// degenerate, so the detector yields nothing — no division, no indexing past the
    /// end, no hang (Tandem_Repeat_Detection.md §6.1). Pinned for the default and a
    /// minimal unit length.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 2, 2).ToList();
        act.Should().NotThrow("an empty sequence has no region long enough to hold a tandem repeat");

        GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 2, 2).Should().BeEmpty();
        GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 1, 2).Should().BeEmpty(
            "even searching 1-bp units, the empty sequence yields no homopolymer run");
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no input". The pre-fix detector
    /// dereferenced `sequence.Sequence` immediately → a raw NullReferenceException.
    /// The contract now guards it with ArgumentNullException (ThrowIfNull), raised
    /// eagerly at the call — never a NullReferenceException
    /// (Tandem_Repeat_Detection.md §3.3, now enforced).
    /// </summary>
    [Test]
    public void FindTandemRepeats_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats((DnaSequence)null!, 2, 2);

        act.Should().Throw<ArgumentNullException>(
            "a null sequence is null-guarded, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region BE — Boundary: single-character sequence

    /// <summary>
    /// BE: a single-character sequence cannot hold a tandem repeat — a tandem needs
    /// ≥ 2 copies (Tandem_Repeat_Detection.md §2.2), and one base is half of even the
    /// shortest possible 1-bp ×2 repeat. The detector must return empty with no crash
    /// and no hang, whether the unit-length floor is the default 2 or the minimal 1.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_SingleCharSequence_IsEmptyAndDoesNotThrow()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 1, 2).ToList();
        act.Should().NotThrow("a single base cannot hold two consecutive copies of any unit");

        GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 1, 2).Should().BeEmpty(
            "one base is too short for even a 1-bp unit repeated twice");
        GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 2, 2).Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — a clear tandem repeat is detected correctly

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, the textbook tandem from the
    /// algorithm doc — "(ATG)3" = "ATGATGATG" — must be detected with the CORRECT unit
    /// and copy count, so the boundary hardening never silently breaks the core
    /// function. Pinned per INV-01..INV-03 (Tandem_Repeat_Detection.md §2.4): unit
    /// "ATG", 3 copies at position 0, TotalLength = 3 × 3 = 9.
    /// </summary>
    [Test]
    public void FindTandemRepeats_ClearTrinucleotideTandem_DetectedWithCorrectUnitAndCount()
    {
        var results = GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), minUnitLength: 3, minRepetitions: 3).ToList();

        var atg = results.Should().ContainSingle(r => r.Unit == "ATG").Subject;
        atg.Repetitions.Should().Be(3, "'ATG' tiles '(ATG)3' exactly three times");
        atg.Position.Should().Be(0, "the tandem block starts at the first base");
        atg.TotalLength.Should().Be(9, "INV-02: TotalLength = |Unit| (3) × Repetitions (3)");
        (atg.Position + atg.TotalLength).Should().BeLessThanOrEqualTo("ATGATGATG".Length,
            "INV-03: the reported tandem stays within the sequence bounds");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — no spurious counts, no out-of-range
    /// positions, no hang — so the degenerate-boundary guards do not corrupt the scan
    /// on ordinary input. Every result must satisfy INV-01..INV-03
    /// (Tandem_Repeat_Detection.md §2.4). Length kept modest because the detector is
    /// O(n²·m); a hang here would trip the timeout.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindTandemRepeats_RandomSequence_ProducesOnlyWellFormedResults()
    {
        const int minUnit = 1;
        const int minReps = 2;
        string seqStr = RandomDna(500, seed: 14_001);
        var seq = new DnaSequence(seqStr);

        var results = GenomicAnalyzer.FindTandemRepeats(seq, minUnit, minReps).ToList();

        results.Should().OnlyContain(r =>
            r.Repetitions >= minReps &&                          // INV-01 (count)
            r.Unit.Length >= minUnit &&                          // INV-01 (unit length)
            r.TotalLength == r.Unit.Length * r.Repetitions &&    // INV-02
            r.Position >= 0 && r.Position + r.TotalLength <= seqStr.Length, // INV-03
            "every result on random input is well-formed: count ≥ minReps, unit ≥ minUnit, total = unit×count, in bounds");
    }

    #endregion

    #endregion
}
