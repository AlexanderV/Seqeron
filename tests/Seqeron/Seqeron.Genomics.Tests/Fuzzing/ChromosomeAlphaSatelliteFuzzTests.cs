using System;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Chromosome area — alpha-satellite (alphoid) monomer-level
/// detection (CHROM-ALPHASAT-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang/infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException from monomer/period windowing past the
/// end, DivideByZeroException on a mean/identity computed over zero monomers,
/// NaN in an identity field). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* validation exception. A
/// raw runtime exception, a hang, an identity outside [0,1], or a spurious
/// "alpha-satellite detected" on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-ALPHASAT-001 — alpha-satellite detection
/// Checklist: docs/checklists/03_FUZZING.md, row 257.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a non-satellite (random) sequence, the empty/too-short
///          sequence, and a PERIOD MISMATCH (a tandem array whose monomer period
///          is NOT ~171 bp, so it must fall outside the scanned [166,176] window).
///   • MC = Malformed Content — non-ACGT characters (digits, IUPAC ambiguity
///          codes, whitespace, lowercase, unicode) fed through the same surface.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The alpha-satellite-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Human centromeric alpha satellite is a tandem repeat of an AT-rich ~171 bp
/// "alphoid" monomer; a subset of monomers carries the 17-bp CENP-B box
/// (Willard 1985; Masumoto et al. 1989; review PMC6121732 —
/// docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md §"Alpha-satellite
/// DNA"; Higher_Order_Repeat_Detection.md §2). The monomer-level detector is
///   ChromosomeAnalyzer.DetectAlphaSatellite(string sequence)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///    lines 700–753). It returns an EAGER value — the AlphaSatelliteResult
///    readonly record struct (IsAlphaSatellite, PeriodicityScore, BestPeriod,
///    AtContent, CenpBBoxCount) — so any exception or hang surfaces at the call
///    itself; no `.ToList()` forcing is needed.
///
/// How the scan works (ChromosomeAnalyzer.cs lines 700–752; the sourced
/// thresholds are documented in the same file, lines 665–689):
///   • AlphaSatelliteMonomerLength = 171 bp; the period search half-width
///     (MonomerPeriodTolerance) is 5, so periods 166..176 are scanned. The
///     MINIMUM input length is 171 + 5 + 1 = 177: anything shorter (incl. empty
///     / null) short-circuits to the no-signal result (false, 0, 0, 0, 0).
///   • AtContent = (#A + #T) / (#A+#C+#G+#T), in [0,1]; the A/C/G/T denominator
///     is guarded (`acgtCount > 0 ? … : 0`), so a sequence with NO A/C/G/T base
///     yields AtContent 0 with NO DivideByZero (the MC guarantee).
///   • PeriodicityScore = best over periods 166..176 of the fraction of bases
///     identical to the base `period` positions upstream, in [0,1]; ~1 for a
///     clean tandem array, ~0.25 for random DNA. The comparison count is guarded
///     (`comparisons > 0 ? … : 0`); the inner loop only ever indexes
///     `sequence[i]` / `sequence[i - period]` for `period <= i < length` with
///     `period < length`, so there is no past-the-end / partial-monomer index.
///   • BestPeriod is 0 (none searched) or one of 166..176.
///   • IsAlphaSatellite ⇔ PeriodicityScore ≥ 0.50 (MinPeriodicityScore, the
///     lower bound of the reported 50–70 % intra-array monomer identity) AND
///     AtContent > 0.50 (the AT-rich signature). BOTH gates must fire.
///   • CenpBBoxCount = #forward matches of the 17-bp IUPAC consensus
///     YTTCGTTGGAARCGGGA (Y=C/T, R=A/G); ≥ 0 and ≤ len − 17 + 1.
///
/// Theory-derived invariants pinned below (derived from the doc rule, NOT read
/// back off the code):
///   • POSITIVE control — a 171-bp AT-rich monomer repeated N times is an exact
///     tandem array: PeriodicityScore = 1.0 at BestPeriod = 171, AtContent > 0.5,
///     ⇒ IsAlphaSatellite TRUE, with CenpBBoxCount = (#CENP-B boxes the chosen
///     monomer contains) × N. (BE/positive anchor.)
///   • NEGATIVE control — random DNA scores ~0.25 periodicity (< 0.50) ⇒
///     IsAlphaSatellite FALSE, no crash. (BE non-satellite.)
///   • PERIOD MISMATCH — a tandem array whose monomer period is far from 171
///     (e.g. a 62-bp unit tiled) has high self-similarity only at period 62,
///     which is OUTSIDE the scanned [166,176] window, so the 171-window
///     PeriodicityScore stays low (< 0.50) ⇒ IsAlphaSatellite FALSE.
///   • EMPTY / TOO-SHORT — any input shorter than 177 bp (incl. "" and null)
///     ⇒ the exact no-signal result (false, 0, 0, 0, 0), no DivideByZero.
///   • Output bounds for EVERY input: PeriodicityScore ∈ [0,1] & finite,
///     AtContent ∈ [0,1] & finite, BestPeriod ∈ {0}∪[166,176], CenpBBoxCount ≥ 0.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm docs: docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md;
///   docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md.
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (DetectAlphaSatellite + AlphaSatelliteResult + AlphaSatelliteConsensus +
///    AlphaSatelliteMonomerLength + CountCenpBBoxes/FindCenpBBoxes).
/// • Willard HF 1985; Masumoto et al. 1989 (CENP-B box); McNulty & Sullivan 2018
///   (PMC6121732).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[Category("Fuzzing")]
public class ChromosomeAlphaSatelliteFuzzTests
{
    #region Helpers

    /// <summary>Documented alphoid monomer length (ChromosomeAnalyzer.AlphaSatelliteMonomerLength).</summary>
    private const int MonomerLength = ChromosomeAnalyzer.AlphaSatelliteMonomerLength; // 171

    /// <summary>Period tolerance half-width = 5 (ChromosomeAnalyzer.MonomerPeriodTolerance), so the
    /// minimum analysable length is 171 + 5 + 1 = 177 and the scanned period window is [166,176].</summary>
    private const int MinAnalysableLength = MonomerLength + 5 + 1; // 177
    private const int LowPeriod = MonomerLength - 5;  // 166
    private const int HighPeriod = MonomerLength + 5; // 176

    /// <summary>The published 62-bp AT-rich alpha-satellite consensus fragment
    /// (ChromosomeAnalyzer.AlphaSatelliteConsensus). Tiling it to exactly 171 bp yields a
    /// realistic AT-rich alphoid monomer for the positive control.</summary>
    private static string Monomer171()
    {
        string c = ChromosomeAnalyzer.AlphaSatelliteConsensus;
        var sb = new StringBuilder(MonomerLength);
        while (sb.Length < MonomerLength)
            sb.Append(c);
        return sb.ToString(0, MonomerLength);
    }

    /// <summary>An exact tandem array: <paramref name="copies"/> identical 171-bp monomers.</summary>
    private static string AlphaArray(int copies) => string.Concat(Enumerable.Repeat(Monomer171(), copies));

    /// <summary>Deterministic non-satellite DNA — seed fixed locally (no shared static Rng).</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Asserts the output-contract bounds that must hold for EVERY input.</summary>
    private static void AssertInContract(ChromosomeAnalyzer.AlphaSatelliteResult r)
    {
        double.IsNaN(r.PeriodicityScore).Should().BeFalse("periodicity is a ratio of counts, never NaN");
        double.IsNaN(r.AtContent).Should().BeFalse("AT content is a guarded ratio, never NaN");
        r.PeriodicityScore.Should().BeInRange(0.0, 1.0, "periodicity is matches/comparisons ∈ [0,1]");
        r.AtContent.Should().BeInRange(0.0, 1.0, "AT content is (A+T)/(ACGT) ∈ [0,1]");
        r.CenpBBoxCount.Should().BeGreaterThanOrEqualTo(0, "a count is non-negative");
        // BestPeriod is 0 (nothing searched) or a value inside the scanned tolerance window.
        if (r.BestPeriod != 0)
            r.BestPeriod.Should().BeInRange(LowPeriod, HighPeriod, "best period is within the 171±5 search window");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-ALPHASAT-001 — alpha-satellite detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-ALPHASAT-001 — alpha-satellite detection

    #region Positive control — exact 171-bp tandem array IS alpha-satellite

    /// <summary>
    /// Theory anchor: a 171-bp AT-rich monomer repeated N times is a PERFECT tandem
    /// array. Self-similarity at the monomer period 171 is exact ⇒ PeriodicityScore
    /// = 1.0 and BestPeriod = 171; the AT-rich consensus pushes AtContent > 0.5; both
    /// gates fire ⇒ IsAlphaSatellite TRUE. CENP-B count must equal (boxes in one
    /// monomer)×N (this consensus monomer carries none → 0). This guards against an
    /// implementation that fails to detect a textbook positive.
    /// </summary>
    [TestCase(4)]
    [TestCase(8)]
    [TestCase(20)]
    public void DetectAlphaSatellite_ExactTandemArray_IsDetectedWithPerfectPeriodicity(int copies)
    {
        string array = AlphaArray(copies);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        AssertInContract(result);
        result.PeriodicityScore.Should().Be(1.0, "every base equals the base 171 positions upstream in an exact array");
        result.BestPeriod.Should().Be(MonomerLength, "the best self-similarity is at the true monomer period");
        result.AtContent.Should().BeGreaterThan(0.50, "the alphoid consensus monomer is AT-rich");
        result.IsAlphaSatellite.Should().BeTrue("both signatures (periodicity ≥ 0.50 and AT > 0.50) are met");

        // CENP-B count scales with copy number: (boxes in one monomer) × copies.
        int boxesPerMonomer = ChromosomeAnalyzer.FindCenpBBoxes(Monomer171()).Count;
        result.CenpBBoxCount.Should().Be(boxesPerMonomer * copies,
            "an exact tandem array repeats each monomer's CENP-B boxes once per copy");
    }

    /// <summary>
    /// Positive detection survives MC injection: lowercasing the entire array must
    /// not change the verdict — the detector uppercases internally before scanning.
    /// </summary>
    [Test]
    public void DetectAlphaSatellite_LowercaseTandemArray_StillDetected()
    {
        string array = AlphaArray(6).ToLowerInvariant();

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        AssertInContract(result);
        result.IsAlphaSatellite.Should().BeTrue("the detector uppercases the sequence before scanning");
        result.PeriodicityScore.Should().Be(1.0);
        result.BestPeriod.Should().Be(MonomerLength);
    }

    #endregion

    #region BE — Boundary: non-satellite (random) sequence is NOT detected

    /// <summary>
    /// BE non-satellite: random DNA has ~0.25 periodicity at any period (four-letter
    /// alphabet) — well below the 0.50 bar — so IsAlphaSatellite is FALSE and the
    /// PeriodicityScore stays low. No crash on a long non-repetitive input.
    /// </summary>
    [TestCase(1)]
    [TestCase(7)]
    [TestCase(99)]
    public void DetectAlphaSatellite_RandomDna_IsNotAlphaSatellite(int seed)
    {
        string seq = RandomDna(2000, seed);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.PeriodicityScore.Should().BeLessThan(0.50,
            "random four-letter DNA self-similarity (~0.25) is below the 0.50 periodicity bar");
        result.IsAlphaSatellite.Should().BeFalse("non-repetitive DNA lacks the ~171 bp tandem signature");
    }

    #endregion

    #region BE — Boundary: empty / null / too-short → no-signal result

    /// <summary>
    /// BE empty: the empty (and null) sequence short-circuits to the exact no-signal
    /// result (false, 0, 0, 0, 0). The mean/identity arithmetic is never reached, so
    /// there is NO DivideByZero on zero monomers.
    /// </summary>
    [TestCase("")]
    [TestCase(null)]
    public void DetectAlphaSatellite_EmptyOrNull_ReturnsNoSignalResult(string? seq)
    {
        var act = () => ChromosomeAnalyzer.DetectAlphaSatellite(seq!);
        act.Should().NotThrow("empty/null is short-circuited before any indexing or division");

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq!);

        AssertInContract(result);
        result.IsAlphaSatellite.Should().BeFalse();
        result.PeriodicityScore.Should().Be(0.0);
        result.BestPeriod.Should().Be(0);
        result.AtContent.Should().Be(0.0);
        result.CenpBBoxCount.Should().Be(0);
    }

    /// <summary>
    /// BE too-short boundary: any length below the minimum analysable length
    /// (171 + 5 + 1 = 177) — including the exact length 176 just under the cliff —
    /// returns the no-signal result with no crash; length 177 is the first analysable
    /// length and merely must obey the output contract.
    /// </summary>
    [TestCase(0)]
    [TestCase(1)]
    [TestCase(171)]
    [TestCase(176)]   // exactly one below the analysable cliff
    [TestCase(177)]   // first analysable length
    [TestCase(300)]
    public void DetectAlphaSatellite_LengthAroundMinimum_NeverCrashesAndShortCircuitsBelowCutoff(int length)
    {
        string seq = RandomDna(length, seed: 1234 + length);

        var act = () => ChromosomeAnalyzer.DetectAlphaSatellite(seq);
        act.Should().NotThrow("length boundaries must be handled, never indexed past the end");

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);
        AssertInContract(result);

        if (length < MinAnalysableLength)
        {
            result.IsAlphaSatellite.Should().BeFalse("below the {0}-bp minimum the no-signal result is returned", MinAnalysableLength);
            result.PeriodicityScore.Should().Be(0.0);
            result.BestPeriod.Should().Be(0);
            result.AtContent.Should().Be(0.0);
            result.CenpBBoxCount.Should().Be(0);
        }
    }

    #endregion

    #region BE — Boundary: PERIOD MISMATCH (tandem period ≠ 171)

    /// <summary>
    /// BE period mismatch: a tandem array whose monomer period is far from 171 (the
    /// 62-bp consensus tiled directly) has its perfect self-similarity only at period
    /// 62, which lies OUTSIDE the scanned [166,176] window. The 171-window periodicity
    /// therefore stays low (well under 0.50), so despite being AT-rich and highly
    /// repetitive the sequence is NOT called alpha-satellite — the monomer-PERIOD
    /// signature, not mere repetitiveness, is what the detector keys on.
    /// </summary>
    [TestCase(40)]
    [TestCase(80)]
    public void DetectAlphaSatellite_WrongPeriodTandemArray_IsNotDetected(int unitCopies)
    {
        // Tile the 62-bp consensus directly: period 62, not a multiple/near 171.
        string seq = string.Concat(Enumerable.Repeat(ChromosomeAnalyzer.AlphaSatelliteConsensus, unitCopies));
        seq.Length.Should().BeGreaterThan(MinAnalysableLength);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.AtContent.Should().BeGreaterThan(0.50, "the consensus fragment is AT-rich");
        result.PeriodicityScore.Should().BeLessThan(0.50,
            "a 62-bp period has no strong self-similarity inside the 171±5 scan window");
        result.IsAlphaSatellite.Should().BeFalse("the ~171 bp period gate is not met by a 62-bp-period array");
    }

    /// <summary>
    /// A nearly-correct period — a 168-bp monomer tiled — DOES sit inside the [166,176]
    /// window, so the detector tolerates the small period drift documented in
    /// Higher_Order_Repeat_Detection.md (167–171 bp real monomers). Output must stay
    /// in-contract and the best period must land inside the search window.
    /// </summary>
    [Test]
    public void DetectAlphaSatellite_NearPeriodWithinTolerance_StaysInContract()
    {
        string monomer168 = Monomer171().Substring(0, 168); // 168 ∈ [166,176]
        string seq = string.Concat(Enumerable.Repeat(monomer168, 6));

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.BestPeriod.Should().Be(168, "an exact 168-bp period is the strongest inside the 171±5 window");
        result.PeriodicityScore.Should().BeGreaterThan(0.50, "a 168-bp period is inside the scanned tolerance window");
    }

    #endregion

    #region MC — Malformed Content: non-ACGT characters never crash

    /// <summary>
    /// MC: a sequence of ONLY non-ACGT characters (digits, gaps, IUPAC ambiguity codes,
    /// whitespace, unicode) has an empty A/C/G/T set, so AtContent must be 0 via the
    /// guarded denominator — NO DivideByZero — IsAlphaSatellite false, and no crash even
    /// when the input is long enough to enter the period scan.
    /// </summary>
    [TestCase("N")]
    [TestCase("-")]
    [TestCase("X9?")]
    [TestCase(" \t")]
    [TestCase("ΩβΣ")]
    public void DetectAlphaSatellite_OnlyNonAcgt_AtContentZeroNoDivByZero(string unit)
    {
        string seq = string.Concat(Enumerable.Repeat(unit, 400)); // ≥ 177 bp, all non-ACGT

        var act = () => ChromosomeAnalyzer.DetectAlphaSatellite(seq);
        act.Should().NotThrow("non-ACGT input is scanned verbatim, never validated-then-indexed-out-of-range");

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.AtContent.Should().Be(0.0, "no A/C/G/T base ⇒ guarded denominator yields 0, not NaN/DivByZero");
        result.IsAlphaSatellite.Should().BeFalse("the AT-rich gate cannot be met with zero AT bases");
    }

    /// <summary>
    /// MC: a real alpha-satellite array sprinkled with arbitrary non-ACGT noise must
    /// still be processed without an unhandled exception and stay fully in-contract.
    /// Determinism is per-test (locally seeded Random).
    /// </summary>
    [TestCase(11)]
    [TestCase(202)]
    public void DetectAlphaSatellite_TandemArrayWithInjectedNoise_StaysInContract(int seed)
    {
        const string noise = "NXRYWSKMBDHV-.* 0123\t";
        var rng = new Random(seed);
        var sb = new StringBuilder(AlphaArray(8));
        // Inject ~30 noise characters at random positions.
        for (int k = 0; k < 30; k++)
        {
            int pos = rng.Next(sb.Length + 1);
            sb.Insert(pos, noise[rng.Next(noise.Length)]);
        }

        var act = () => ChromosomeAnalyzer.DetectAlphaSatellite(sb.ToString());
        act.Should().NotThrow("malformed bases inside a real array must not crash the scan");

        AssertInContract(ChromosomeAnalyzer.DetectAlphaSatellite(sb.ToString()));
    }

    #endregion

    #region BE — degenerate homopolymers (period gate vs AT gate)

    /// <summary>
    /// BE degenerate: an all-A homopolymer is trivially periodic at EVERY period
    /// (PeriodicityScore 1.0) and AT content 1.0, so BOTH gates fire and the detector
    /// reports IsAlphaSatellite TRUE. This is a documented heuristic limitation (a
    /// pathological non-biological input), not a crash — the point is that the output
    /// is well-defined and in-contract, never a runtime exception.
    /// </summary>
    [Test]
    public void DetectAlphaSatellite_HomopolymerA_IsWellDefinedAndInContract()
    {
        string seq = new string('A', 600);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.PeriodicityScore.Should().Be(1.0, "a homopolymer is identical at every offset");
        result.AtContent.Should().Be(1.0, "all bases are A");
        result.IsAlphaSatellite.Should().BeTrue("both gates trivially fire for an all-AT homopolymer (heuristic limit)");
    }

    /// <summary>
    /// BE degenerate: an all-G homopolymer is equally periodic (1.0) but has AT content
    /// 0, so the AT-rich gate vetoes detection ⇒ IsAlphaSatellite FALSE. Confirms the
    /// AT gate is a hard AND, not folded away.
    /// </summary>
    [Test]
    public void DetectAlphaSatellite_HomopolymerG_PeriodicButNotAtRich_IsNotDetected()
    {
        string seq = new string('G', 600);

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.PeriodicityScore.Should().Be(1.0, "a homopolymer is identical at every offset");
        result.AtContent.Should().Be(0.0, "no A/T bases");
        result.IsAlphaSatellite.Should().BeFalse("the AT-rich gate vetoes a GC homopolymer");
    }

    #endregion

    #region BE — very long input terminates promptly

    /// <summary>
    /// BE scale: a large exact tandem array (≈ 1.7 Mb) must terminate promptly — the
    /// scan is O(window × length) with a fixed 11-period window — and remain in-contract
    /// with the perfect-array verdict. Bounded by CancelAfter so a hang fails loudly.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void DetectAlphaSatellite_VeryLongArray_TerminatesInContract()
    {
        string seq = AlphaArray(10_000); // 10_000 × 171 ≈ 1.71 Mb

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(seq);

        AssertInContract(result);
        result.IsAlphaSatellite.Should().BeTrue();
        result.PeriodicityScore.Should().Be(1.0);
        result.BestPeriod.Should().Be(MonomerLength);
    }

    #endregion

    #endregion
}
