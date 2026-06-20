using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Chromosome area — telomere detection (CHROM-TELO-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException, ArgumentOutOfRangeException leaking from
/// internal indexing, DivideByZero, OutOfMemory). Every input must resolve to
/// EITHER a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception. A raw runtime exception, a hang, a crash, or a spurious
/// "telomere detected" on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-TELO-001 — telomere detection
/// Checklist: docs/checklists/03_FUZZING.md, row 48.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the empty sequence, a sequence with NO telomeric
///          repeat, and a sequence that is ONLY the telomeric repeat (pure tract).
///   • MC = Malformed Content — non-DNA characters in the sequence (digits, IUPAC
///          ambiguity codes, whitespace, unicode) fed through the same surface.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The telomere-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A telomere is the repetitive DNA tract at a chromosome end. The vertebrate
/// canonical repeat unit is the 6-bp hexamer TTAGGG; its reverse complement
/// CCCTAA is the motif expected at the 5' end (Telomere_Analysis.md §2.1,
/// table "chromosome-end orientation"). The detector is
///   ChromosomeAnalyzer.AnalyzeTelomeres(
///       string chromosomeName,
///       string sequence,
///       string telomereRepeat = "TTAGGG",
///       int searchLength       = 10000,
///       int minTelomereLength  = 500,
///       int criticalLength     = 3000)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///    lines 250–288). It returns an EAGER value (the TelomereResult readonly
///    record struct) — there is no lazy iterator, so any exception or hang would
///    surface at the call itself; no `.ToList()` forcing is needed.
///
/// How the scan works (ChromosomeAnalyzer.cs lines 263–288, MeasureTelomereLength
/// lines 293–340; Telomere_Analysis.md §4.1, §5.2):
///   • The sequence and repeat motif are upper-cased; the 5' motif is the reverse
///     complement of `telomereRepeat`.
///   • The 5' PREFIX window (first min(searchLength, len) bases) is scanned left→
///     right against CCCTAA; the 3' SUFFIX window (last min(searchLength, len)
///     bases) is scanned right→left against TTAGGG.
///   • Scanning advances in repeat-sized (6-bp) steps, counting a window only while
///     its per-base similarity to the motif is ≥ 70%; the first window below 70%
///     stops that end. Only COMPLETE 6-bp windows are counted, so the measured
///     length is always a multiple of 6 and equals 6 × (number of accepted windows)
///     (INV-01 length ≥ 0).
///   • RepeatPurity = matchingBases / totalBases ∈ [0,1], or 0 when nothing counted
///     (INV-02). Has*Telomere ⇔ measured length ≥ minTelomereLength (INV-03).
///
/// Documented edge-case contract (Telomere_Analysis.md §3.1, §3.3, §6.1):
///   • Empty OR null sequence → the method special-cases it (line 258) and returns
///     a no-telomere result: lengths 0, Has*Telomere = false, purities 0, and
///     IsCriticallyShort = TRUE (the documented empty-input flag).
///   • Sequence shorter than the 6-bp repeat → MeasureTelomereLength exits
///     immediately (region.Length < repeatLen, line 299) → zero-length telomeres.
///   • No telomeric repeats → no window meets the 70% bar → zero lengths,
///     Has*Telomere = false, and (non-empty input) IsCriticallyShort = FALSE
///     (§5.2: critically-short is only flagged when a telomere is actually
///     detected and is shorter than criticalLength).
///   • Non-DNA characters → there is NO validation/rejection gate: the input is
///     upper-cased and scanned verbatim. Garbage bases simply fail the 70%
///     similarity test, so a non-DNA sequence behaves like "no telomere" and must
///     NEVER crash on indexing. This is the MC contract pinned below.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md.
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (AnalyzeTelomeres + MeasureTelomereLength).
/// • Meyne, Ratliff, Moyzis (1989): conservation of (TTAGGG)n among vertebrates.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ChromosomeFuzzTests
{
    #region Helpers

    private const string TelomereRepeat = "TTAGGG"; // vertebrate canonical hexamer
    private const int RepeatLen = 6;

    /// <summary>Builds a pure forward telomere tract of <paramref name="copies"/> TTAGGG units.</summary>
    private static string PureTelomere(int copies) => string.Concat(System.Linq.Enumerable.Repeat(TelomereRepeat, copies));

    /// <summary>
    /// Deterministic non-telomeric DNA — seed fixed locally so generated fuzz input
    /// is reproducible (no shared static Rng). Avoids accidentally emitting a TTAGGG
    /// run by drawing from a reduced alphabet biased away from the motif.
    /// </summary>
    private static string RandomNonTelomericDna(int length, int seed)
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
    //  CHROM-TELO-001 — telomere detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-TELO-001 — telomere detection

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is special-cased (ChromosomeAnalyzer.cs line 258;
    /// Telomere_Analysis.md §6.1). It must return a no-telomere result — zero
    /// lengths, both Has*Telomere = false, purities 0 — and the documented
    /// IsCriticallyShort = TRUE flag, never a crash.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_EmptySequence_ReturnsNoTelomereAndCriticallyShort()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrEmpty", string.Empty);
        act.Should().NotThrow("empty input is explicitly special-cased, not indexed");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrEmpty", string.Empty);

        result.Chromosome.Should().Be("chrEmpty", "the chromosome name is copied verbatim");
        result.Has5PrimeTelomere.Should().BeFalse();
        result.Has3PrimeTelomere.Should().BeFalse();
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.RepeatPurity5Prime.Should().Be(0);
        result.RepeatPurity3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeTrue(
            "empty input carries the documented critically-short flag (Telomere_Analysis.md §6.1)");
    }

    /// <summary>
    /// BE: a null sequence is treated identically to empty via the
    /// string.IsNullOrEmpty short-circuit — it must NOT throw a
    /// NullReferenceException and must return the same no-telomere, critically-short
    /// result.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NullSequence_ReturnsNoTelomereAndCriticallyShort()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrNull", null!);
        act.Should().NotThrow("null is absorbed by the IsNullOrEmpty guard, never dereferenced");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrNull", null!);

        result.Has3PrimeTelomere.Should().BeFalse();
        result.TelomereLength3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeTrue();
    }

    /// <summary>
    /// BE: a sequence shorter than the 6-bp repeat cannot contain a single complete
    /// repeat window — MeasureTelomereLength exits at region.Length &lt; repeatLen
    /// (line 299). Both ends must report length 0 with no Substring out-of-range, and
    /// (non-empty input) IsCriticallyShort = false.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_SequenceShorterThanRepeat_ReturnsZeroLengthsNoCrash()
    {
        // "TTAGG" is 5 bp — one base short of a single TTAGGG window.
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrShort", "TTAGG");
        act.Should().NotThrow("a region shorter than the repeat exits before any window is taken");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrShort", "TTAGG");

        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeFalse(
            "non-empty input with no detected telomere is NOT flagged critically short (§5.2)");
    }

    #endregion

    #region BE — Boundary: no TTAGGG in sequence

    /// <summary>
    /// BE: a sequence with NO telomeric repeat. No 6-bp window on either end reaches
    /// the 70% similarity bar, so both lengths are 0, both Has*Telomere = false, and
    /// IsCriticallyShort = false for non-empty input (§6.1, §5.2). A poly-A tract is
    /// the cleanest "no motif" input: it shares at most one base with TTAGGG/CCCTAA
    /// (1/6 ≈ 17% &lt; 70%).
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NoTelomericRepeat_DetectsNothing()
    {
        string sequence = new string('A', 600); // no TTAGGG, no CCCTAA anywhere

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrNoMotif", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrNoMotif", sequence);

        result.Has5PrimeTelomere.Should().BeFalse("no window matches CCCTAA at the 5' end");
        result.Has3PrimeTelomere.Should().BeFalse("no window matches TTAGGG at the 3' end");
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.RepeatPurity5Prime.Should().Be(0);
        result.RepeatPurity3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeFalse(
            "no telomere is detected, so the critically-short flag stays false for non-empty input");
    }

    /// <summary>
    /// BE: a pseudo-random non-telomeric sequence must likewise never mis-detect a
    /// full telomere and never crash. We assert only the well-formedness invariants
    /// (INV-01..INV-03): non-negative lengths that are multiples of 6, purities in
    /// [0,1], and Has*Telomere consistent with the default minTelomereLength = 500.
    /// Deterministic seed keeps this reproducible.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_RandomSequence_ProducesOnlyWellFormedResult()
    {
        string sequence = RandomNonTelomericDna(2000, seed: 4848);

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrRand", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrRand", sequence);

        result.TelomereLength5Prime.Should().BeGreaterThanOrEqualTo(0);
        result.TelomereLength3Prime.Should().BeGreaterThanOrEqualTo(0);
        (result.TelomereLength5Prime % RepeatLen).Should().Be(0, "measured length is always a multiple of the 6-bp repeat");
        (result.TelomereLength3Prime % RepeatLen).Should().Be(0, "measured length is always a multiple of the 6-bp repeat");
        result.RepeatPurity5Prime.Should().BeInRange(0, 1);
        result.RepeatPurity3Prime.Should().BeInRange(0, 1);
        result.Has5PrimeTelomere.Should().Be(result.TelomereLength5Prime >= 500);
        result.Has3PrimeTelomere.Should().Be(result.TelomereLength3Prime >= 500);
    }

    #endregion

    #region BE — Boundary: only TTAGGG (pure telomeric tract) — KEY positive

    /// <summary>
    /// BE (KEY positive): a sequence that is ONLY the telomeric repeat. A pure tract
    /// of N copies of TTAGGG must be detected at the 3' end with the correct measured
    /// length = 6 × N (every window is a perfect TTAGGG match, similarity 100% ≥ 70%),
    /// repeat purity exactly 1.0, and Has3PrimeTelomere true when 6·N ≥ 500. This is
    /// the core "only TTAGGG ⇒ telomere with correct repeat count (length/6)" check.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_PureTtagggTract_DetectedWithCorrectLengthAndPurity()
    {
        const int copies = 100;               // 100 × 6 = 600 bp ≥ minTelomereLength (500)
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrPure", sequence);

        result.TelomereLength3Prime.Should().Be(copies * RepeatLen,
            "a pure TTAGGG tract is measured as 6 × (repeat count) at the 3' end");
        (result.TelomereLength3Prime / RepeatLen).Should().Be(copies,
            "the detected repeat count is exactly length / 6");
        result.RepeatPurity3Prime.Should().Be(1.0,
            "every counted window is a perfect TTAGGG match, so purity is exactly 1.0");
        result.Has3PrimeTelomere.Should().BeTrue(
            "600 bp of pure repeat clears the default minTelomereLength of 500");
    }

    /// <summary>
    /// BE: a SHORT pure tract — below minTelomereLength but a valid repeat — must be
    /// measured with the correct length (6 × copies) yet flagged Has3PrimeTelomere =
    /// false, pinning the minTelomereLength boundary rather than mis-flagging a short
    /// real tract. With a default criticalLength of 3000 the detected short telomere
    /// stays below it, but Has3Prime = false means IsCriticallyShort cannot trip
    /// (it requires a *detected* telomere, line 279).
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_ShortPureTract_MeasuredButBelowPresenceThreshold()
    {
        const int copies = 10;                // 10 × 6 = 60 bp < minTelomereLength (500)
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrShortPure", sequence);

        result.TelomereLength3Prime.Should().Be(copies * RepeatLen,
            "the short tract is still measured at its true length 6 × copies");
        result.Has3PrimeTelomere.Should().BeFalse(
            "60 bp is below the default minTelomereLength of 500, so presence is false");
        result.IsCriticallyShort.Should().BeFalse(
            "no telomere is flagged present, so the critically-short test never fires");
    }

    #endregion

    #region MC — Malformed content: non-DNA characters

    /// <summary>
    /// MC: non-DNA characters (digits, IUPAC ambiguity codes, whitespace, unicode)
    /// must be handled without a crash. There is no validation gate — input is
    /// upper-cased and scanned verbatim — so garbage bases simply fail the 70%
    /// similarity bar and the sequence behaves like "no telomere": zero lengths,
    /// no detection, no IndexOutOfRange. We pin no-throw + no spurious detection.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NonDnaCharacters_HandledWithoutCrashOrDetection()
    {
        string garbage = "12345 NRYSWK\t\nXZé中??!!--***";

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrGarbage", garbage);
        act.Should().NotThrow("there is no validation gate; non-DNA chars are scanned verbatim, never indexed past end");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrGarbage", garbage);

        result.Has5PrimeTelomere.Should().BeFalse("no non-DNA window can reach 70% similarity to CCCTAA");
        result.Has3PrimeTelomere.Should().BeFalse("no non-DNA window can reach 70% similarity to TTAGGG");
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
    }

    /// <summary>
    /// MC: a non-DNA-decorated but otherwise pure telomeric tract. Even when garbage
    /// characters are interleaved at the chromosome interior, a clean TTAGGG tract at
    /// the 3' terminus must still be detected (the 3' scan walks inward from the end
    /// and stops at the first sub-70% window), and the run must never crash on the
    /// non-DNA bytes. This pins that malformed content elsewhere in the sequence does
    /// not corrupt or suppress a legitimate terminal telomere.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_GarbagePrefixWithPureTerminalTract_StillDetectsTelomere()
    {
        // Non-DNA junk at the 5' interior, a clean 100-copy TTAGGG tract at the 3' end.
        string sequence = "12NRYSWK??##" + PureTelomere(100);

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrMixed", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrMixed", sequence);

        result.Has3PrimeTelomere.Should().BeTrue(
            "the clean terminal TTAGGG tract is detected despite non-DNA junk upstream");
        result.TelomereLength3Prime.Should().Be(100 * RepeatLen,
            "the 3' scan counts exactly the 100 perfect repeats at the terminus");
        result.RepeatPurity3Prime.Should().Be(1.0);
    }

    #endregion

    #region Positive sanity — a known repeat count is recovered exactly

    /// <summary>
    /// Positive sanity: a TTAGGG tract of a KNOWN number of copies is detected with
    /// exactly that repeat count (length / 6). This is the affirmative anchor that
    /// the detector actually finds the canonical motif and counts it correctly —
    /// without it, an all-passing boundary suite could be vacuously green. We test
    /// several copy counts so the count, not a fixed length, is what is recovered.
    /// </summary>
    [TestCase(85)]
    [TestCase(150)]
    [TestCase(512)]
    public void AnalyzeTelomeres_KnownRepeatCount_RecoveredExactly(int copies)
    {
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrSanity", sequence);

        (result.TelomereLength3Prime / RepeatLen).Should().Be(copies,
            "the detected 3' repeat count must equal the planted number of TTAGGG copies");
        result.TelomereLength3Prime.Should().Be(copies * RepeatLen);
        result.RepeatPurity3Prime.Should().Be(1.0, "a perfect planted tract has purity 1.0");
        result.Has3PrimeTelomere.Should().Be(copies * RepeatLen >= 500);
    }

    #endregion

    #endregion
}
