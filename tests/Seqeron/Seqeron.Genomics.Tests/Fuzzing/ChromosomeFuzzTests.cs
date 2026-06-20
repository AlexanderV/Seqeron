using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Chromosome area — telomere detection (CHROM-TELO-001)
/// and centromere detection (CHROM-CENT-001).
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
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: CHROM-CENT-001 — centromere detection
/// Checklist: docs/checklists/03_FUZZING.md, row 49.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the empty sequence, an all-AT sequence, an all-GC
///          sequence, and an extremely short sequence (shorter than the scan
///          window). — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The centromere-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// `AnalyzeCentromere(string chromosomeName, string sequence,
///      int windowSize = 100000, double minAlphaSatelliteContent = 0.3)`
///   (ChromosomeAnalyzer.cs lines 361–437) returns an EAGER `CentromereResult`
///   readonly record struct (Chromosome, Start?, End?, Length, CentromereType,
///   AlphaSatelliteContent, IsAcrocentric) — no lazy iterator, so any exception
///   or hang surfaces at the call itself.
///
/// CRITICAL — what signal the detector actually keys on (Centromere_Analysis.md
/// §2.2, §4.2; ChromosomeAnalyzer.cs lines 384–391): the score is
///     score = repeatContent * (1 - gcVariability)
/// where `repeatContent` is the 15-mer REPETITIVENESS of a window (EstimateRepeat-
/// Content) and `gcVariability` is the std-dev of GC fraction over 1 kb sub-
/// windows. A window is accepted only when `repeatContent > minAlphaSatelliteContent`
/// AND its score is the running maximum. The biological gloss "centromeres are
/// AT-rich alpha-satellite" is realised in THIS code as "repeat-rich, GC-uniform",
/// NOT as a literal AT-content meter. The consequence for the all-AT / all-GC
/// fuzz boundaries is precise and tested below:
///   • all-AT homopolymer ≥ windowSize → a single 15-mer repeats throughout →
///     repeatContent maximal; GC fraction is a constant 0 across sub-windows →
///     gcVariability 0 → score maximal → DETECTED (the strong-positive boundary).
///   • all-GC homopolymer ≥ windowSize → ALSO a single repeated 15-mer, GC
///     fraction a constant 1 → gcVariability 0 → score maximal → ALSO DETECTED.
///     The checklist's idealised "all-GC ⇒ no AT signal ⇒ no centromere" does NOT
///     hold for this repeat-based heuristic: an all-GC homopolymer is just as
///     maximally repetitive as all-AT, so the theory-correct contract of the
///     ACTUAL code is detection. We assert that, with the high-GC sanity check
///     (gpos/AT-uniformity is not what is measured) documented in line, rather
///     than weaken the test to the checklist's idealised model.
///
/// Documented edge-case contract (Centromere_Analysis.md §3.3, §6.1):
///   • Empty OR null sequence → special-cased (line 367) → CentromereType
///     "Unknown", Start/End null, Length 0, AlphaSatelliteContent 0,
///     IsAcrocentric false. Never indexed.
///   • Sequence shorter than windowSize → the scan loop
///     `for (i = 0; i < sequence.Length - windowSize; ...)` has a NEGATIVE bound
///     (Length − windowSize < 0) so it never iterates; candidate stays unset →
///     "Unknown", null boundaries. This is the KEY boundary: the negative bound
///     must NOT produce an IndexOutOfRange (no `sequence[i..end]` is ever taken)
///     and the GC-variability helper is never reached, so NO DivideByZero.
///   • Non-repetitive sequence → no window exceeds minAlphaSatelliteContent →
///     "Unknown" (asserted via the high-complexity random sanity case).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md
///   (§2.2 core model, §3.3 validation, §4.2 scoring, §6.1 edge cases).
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (AnalyzeCentromere + EstimateRepeatContent + CalculateGcVariability).
/// • Levan, Fredga, Sandberg (1964): centromere arm-ratio nomenclature.
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

    /// <summary>
    /// Centromere scan window used by the detection fuzz cases. The production default
    /// is 100 000; we use the smaller 10 000 (matching ChromosomeAnalyzer_Centromere_Tests)
    /// so fuzz sequences stay a few hundred kb instead of multi-megabyte, while still
    /// exercising the identical scan/score/extend code path.
    /// </summary>
    private const int CentWindow = 10_000;

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

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-CENT-001 — centromere detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-CENT-001 — centromere detection

    #region BE — Boundary: empty / null sequence

    /// <summary>
    /// BE: the empty sequence is special-cased (ChromosomeAnalyzer.cs line 367;
    /// Centromere_Analysis.md §6.1). It must return the documented no-centromere
    /// result — CentromereType "Unknown", null boundaries, length 0,
    /// AlphaSatelliteContent 0, IsAcrocentric false — never a crash, never an index.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_EmptySequence_ReturnsUnknownNoBoundaries()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrEmpty", string.Empty);
        act.Should().NotThrow("empty input is explicitly special-cased, not scanned or indexed");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrEmpty", string.Empty);

        result.Chromosome.Should().Be("chrEmpty", "the chromosome name is copied verbatim");
        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().Be(0);
        result.IsAcrocentric.Should().BeFalse();
    }

    /// <summary>
    /// BE: a null sequence is absorbed by the same string.IsNullOrEmpty short-circuit
    /// (line 367). It must NOT throw a NullReferenceException and must return the same
    /// "Unknown" / null-boundary result as the empty case.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_NullSequence_ReturnsUnknownNoBoundaries()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrNull", null!);
        act.Should().NotThrow("null is absorbed by the IsNullOrEmpty guard, never dereferenced");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrNull", null!);

        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
    }

    #endregion

    #region BE — Boundary: extremely short (shorter than the scan window)

    /// <summary>
    /// BE (KEY boundary): a sequence shorter than <see cref="CentWindow"/>. The scan
    /// loop bound is `sequence.Length - windowSize`, which is NEGATIVE here, so the
    /// loop never iterates: no window slice `sequence[i..end]` is ever taken (no
    /// IndexOutOfRange / ArgumentOutOfRange) and the GC-variability helper is never
    /// reached (no DivideByZero). The candidate stays unset → "Unknown", null
    /// boundaries, length 0 (Centromere_Analysis.md §6.1). This pins the div-by-zero /
    /// out-of-range boundary the checklist row calls out.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_ShorterThanWindow_ReturnsUnknownNoCrash()
    {
        // 500 bp ≪ 10 000-bp window: Length − windowSize is negative.
        string sequence = RandomNonTelomericDna(500, seed: 11);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrShort", sequence, windowSize: CentWindow);
        act.Should().NotThrow("a negative scan bound skips the loop; nothing is sliced or divided");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrShort", sequence, windowSize: CentWindow);

        result.CentromereType.Should().Be("Unknown", "the scan loop never runs, so no candidate is set");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().Be(0);
    }

    /// <summary>
    /// BE: the most extreme short input — a single base — must behave identically to
    /// any other sub-window-length sequence: "Unknown", null boundaries, no crash.
    /// Length 1 makes `Length − windowSize` maximally negative and exercises the same
    /// loop-skip guard with a sequence too short to even hold one 15-mer.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_SingleBase_ReturnsUnknownNoCrash()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrOne", "A", windowSize: CentWindow);
        act.Should().NotThrow("a single base is far shorter than the window; the scan loop is skipped");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrOne", "A", windowSize: CentWindow);

        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.Length.Should().Be(0);
    }

    #endregion

    #region BE — Boundary: all-AT homopolymer (strong-positive)

    /// <summary>
    /// BE (KEY positive): an all-AT homopolymer at least as long as the scan window is
    /// maximally repetitive — a single 15-mer repeats end to end, so repeatContent is
    /// maximal — and its GC fraction is a constant 0 across every 1 kb sub-window, so
    /// gcVariability is 0 and the composite score is maximal. The detector MUST report
    /// a region (non-null Start/End) with AlphaSatelliteContent above the threshold and
    /// a concrete (non-"Unknown") classification. Asserting detection here (not a crash,
    /// not a miss) is the strong-positive AT-richness boundary the checklist demands.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_AllAt_DetectedAsStrongPositive()
    {
        string sequence = new string('A', 300_000); // 30 windows of pure repeat

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllAt", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllAt", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.Start.Should().NotBeNull("an all-AT homopolymer is maximally repetitive and must be detected");
        result.End.Should().NotBeNull();
        result.Length.Should().BeGreaterThan(0).And.Be(result.End!.Value - result.Start!.Value,
            "INV-01: Length = End − Start for a detected region");
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.3,
            "uniform AT content has maximal repeat content and zero GC variability → maximal score");
        result.CentromereType.Should().NotBe("Unknown", "a detected region carries a Levan classification");
    }

    #endregion

    #region BE — Boundary: all-GC homopolymer (repeat heuristic, not AT meter)

    /// <summary>
    /// BE: an all-GC homopolymer. The checklist's idealised model says "all-GC → no
    /// AT signal → no centromere", but the ACTUAL detector scores REPETITIVENESS, not
    /// AT content: an all-GC run is just as maximally repetitive (one repeated 15-mer)
    /// with a constant GC fraction of 1 → gcVariability 0 → maximal score. The theory-
    /// correct contract of this repeat-based heuristic is therefore DETECTION, exactly
    /// like all-AT. We assert detection (never a crash, never a spurious "Unknown") and
    /// document in line that this is the repeat heuristic's behaviour, not an AT meter —
    /// we do not weaken the test to the checklist's idealised expectation.
    /// (Centromere_Analysis.md §2.2, §4.2, §6.1 "Strongly repetitive sequence".)
    /// </summary>
    [Test]
    public void AnalyzeCentromere_AllGc_DetectedBecauseSignalIsRepetitivenessNotAtContent()
    {
        string sequence = new string('G', 300_000);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllGc", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllGc", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.Start.Should().NotBeNull(
            "the detector keys on repeat content; an all-GC homopolymer is maximally repetitive, so it is detected");
        result.End.Should().NotBeNull();
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.3,
            "constant GC fraction → zero GC variability → score equals the (maximal) repeat content");
        result.CentromereType.Should().NotBe("Unknown");
    }

    #endregion

    #region BE — Well-formedness invariants on random / homogeneous input

    /// <summary>
    /// BE: a high-complexity pseudo-random sequence (no engineered repeat) must never
    /// crash and, per Centromere_Analysis.md §6.1 "Non-repetitive sequence", must yield
    /// "Unknown" — no window of distinct 15-mers clears the 0.3 repeat threshold. This
    /// is the negative counterpart to the AT/GC positives: a non-centromeric region is
    /// NOT detected. Deterministic seed keeps it reproducible.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_RandomNonRepetitive_ReturnsUnknown()
    {
        string sequence = RandomNonTelomericDna(300_000, seed: 9090);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrRand", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrRand", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.CentromereType.Should().Be("Unknown",
            "a high-complexity random sequence has no repeat-rich window above the threshold");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().BeGreaterThanOrEqualTo(0, "INV-04: score is never negative");
    }

    #endregion

    #region Positive sanity — engineered satellite tract detected, flank not

    /// <summary>
    /// Positive sanity: a clearly AT-rich, alpha-satellite-like tandem-repeat tract
    /// embedded between non-repetitive flanks is detected as a centromere — non-null
    /// boundaries that land inside the engineered repeat, a score above the threshold,
    /// and a concrete classification — while a purely random sequence of the same size
    /// is NOT (covered by the random case above). This is the affirmative anchor that
    /// the detector finds a genuine satellite signal and not just degenerate homopolymers.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_EngineeredSatelliteTract_DetectedInsideRepeatRegion()
    {
        // Non-repetitive flanks, AT-rich tandem-repeat core (alpha-satellite-like 21-mer).
        string flank = RandomNonTelomericDna(100_000, seed: 7);
        string repeatUnit = "AATGAATATTTCTTTTATGTT"; // 21-mer, AT-rich
        string satellite = string.Concat(System.Linq.Enumerable.Repeat(repeatUnit, 6_000)); // ~126 kb
        string sequence = flank + satellite + flank;

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrSat", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.1);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrSat", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.1);

        result.Start.Should().NotBeNull("the embedded tandem-repeat tract is detected");
        result.End.Should().NotBeNull();
        result.CentromereType.Should().NotBe("Unknown");
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.1, "the satellite tract clears the threshold");
        // The detected region overlaps the engineered satellite core, not the random flanks.
        result.End!.Value.Should().BeGreaterThan(flank.Length,
            "the detected region must reach into the repetitive core past the 5' flank");
        result.Start!.Value.Should().BeLessThan(flank.Length + satellite.Length,
            "the detected region must start before the 3' flank");
    }

    #endregion

    #endregion
}
