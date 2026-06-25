using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Epigenetics area — chromatin state annotation
/// (EPIGEN-CHROM-001). The chromatin-state contract under test lives in
/// <see cref="EpigeneticsAnalyzer"/> and is exposed by two public entry points:
///   • <see cref="EpigeneticsAnalyzer.PredictChromatinState"/> — maps a six-mark
///     signal vector (H3K4me3, H3K4me1, H3K27ac, H3K36me3, H3K27me3, H3K9me3) to
///     a single <see cref="ChromatinState"/> via the ChromHMM/Roadmap present-mark
///     signature rules.
///   • <see cref="EpigeneticsAnalyzer.AnnotateHistoneModifications"/> — streams a
///     sequence of (Start, End, Mark, Signal) regions and labels EACH region with
///     the canonical Roadmap state implied by its single histone mark.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (NullReference /
/// DivideByZero / KeyNotFound / IndexOutOfRange). Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a *documented, intentional*
/// outcome. For chromatin-state annotation the headline hazards are:
///   • a NullReferenceException when the region/modification input is empty or the
///     mark string is null;
///   • a DivideByZeroException on an empty input (e.g. averaging over zero marks);
///   • a KeyNotFoundException when a mark name is not one of the recognised marks;
///   • a FALSE active/regulatory state assigned to a region carrying NO present
///     mark (must be the quiescent default LowSignal, never ActivePromoter etc.);
///   • a non-deterministic state for the same mark vector across calls.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-CHROM-001 — chromatin state annotation (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 183.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///     Targets (checklist row 183): "empty, single region, no marks".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// Mapping of the checklist BE targets onto this source API:
///   • "empty"         → AnnotateHistoneModifications over an EMPTY region
///                        sequence → empty result, no crash.
///   • "single region" → AnnotateHistoneModifications over a 1-element sequence →
///                        exactly one labelled region with its documented state.
///   • "no marks"      → PredictChromatinState with every mark below threshold
///                        (and a region whose mark is below threshold) →
///                        LowSignal (Roadmap Quies), never a false active state.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The chromatin-state contract under test (Roadmap 15/18-state model)
/// ───────────────────────────────────────────────────────────────────────────
/// Each mark is binarized to present/absent at presenceThreshold (inclusive
/// `signal >= threshold`; default 0.5). The present-mark SET maps to a state:
///   H3K4me3 + H3K27me3            → BivalentPromoter (TssBiv)   [INV-03]
///   H3K4me1 + H3K27me3            → BivalentEnhancer (EnhBiv)
///   H3K4me3 (± H3K4me1)           → ActivePromoter  (TssA; TSS ranks above Enh)
///   H3K4me1 + H3K27ac            → ActiveEnhancer  (active Enh)
///   H3K4me1, no H3K27ac          → WeakEnhancer    (Enh)
///   H3K36me3                      → Transcribed     (Tx)
///   H3K27me3 (alone)             → Repressed       (ReprPC)
///   H3K9me3 (alone)             → Heterochromatin (Het)
///   none present                  → LowSignal       (Quies)        [INV-02]
///   — docs/algorithms/Epigenetics/Chromatin_State_Prediction.md §2.2, §6.1.
///
/// Invariants pinned here:
///   INV-01 state depends only on the present/absent pattern (equal pattern ⇒
///          equal state) — also gives determinism;
///   INV-02 no mark present ⇒ LowSignal;
///   INV-04 the result is always a DEFINED ChromatinState (total function).
///   — docs/algorithms/Epigenetics/Chromatin_State_Prediction.md §2.4.
///
/// AnnotateHistoneModifications single-mark mapping (source InferStateFromMark):
///   below threshold              → LowSignal
///   H3K4ME3                      → ActivePromoter
///   H3K4ME1                      → WeakEnhancer
///   H3K27AC                      → ActiveEnhancer
///   H3K36ME3                     → Transcribed
///   H3K27ME3                     → Repressed
///   H3K9ME3                      → Heterochromatin
///   H3K9AC                       → ActivePromoter
///   any other (unknown) mark     → LowSignal  (no KeyNotFound)
///   — docs/algorithms/Epigenetics/Chromatin_State_Prediction.md §3.3;
///     EpigeneticsAnalyzer.cs InferStateFromMark.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticChromatinFuzzTests
{
    private const double Threshold = 0.5; // documented default presenceThreshold

    // The complete, documented chromatin-state set (Roadmap 15/18-state model).
    private static readonly HashSet<ChromatinState> DefinedStates =
        new(Enum.GetValues<ChromatinState>());

    // ── Well-formed-result helper ───────────────────────────────────────────
    // Pins INV-04 + the per-region contract: AnnotateHistoneModifications must
    // echo back each input region verbatim and assign EXACTLY ONE state drawn
    // from the documented state set. This stops a test from rubber-stamping a
    // garbage/undefined state green.
    private static void AssertWellFormedAnnotation(
        HistoneModification result,
        (int Start, int End, string Mark, double Signal) input)
    {
        result.Start.Should().Be(input.Start);
        result.End.Should().Be(input.End);
        result.Mark.Should().Be(input.Mark);
        result.Signal.Should().Be(input.Signal);
        DefinedStates.Should().Contain(result.PredictedState,
            "every region must be labelled with a defined ChromatinState (INV-04)");
    }

    #region EPIGEN-CHROM-001 — chromatin state annotation

    // ════════════════════════════════════════════════════════════════════════
    // BE: empty input
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    [CancelAfter(5000)]
    public void AnnotateHistoneModifications_EmptyInput_YieldsEmpty_NoCrash()
    {
        // "empty" target: zero regions in → zero annotations out, no exception.
        var empty = Array.Empty<(int, int, string, double)>();

        Action act = () => AnnotateHistoneModifications(empty).ToList();
        act.Should().NotThrow();

        AnnotateHistoneModifications(empty).Should().BeEmpty();
    }

    [Test]
    [CancelAfter(5000)]
    public void AnnotateHistoneModifications_EmptyEnumerable_YieldsEmpty()
    {
        // A lazily-empty Enumerable (not just an empty array) must also be safe.
        IEnumerable<(int, int, string, double)> none =
            Enumerable.Empty<(int, int, string, double)>();

        AnnotateHistoneModifications(none).ToList().Should().BeEmpty();
    }

    // ════════════════════════════════════════════════════════════════════════
    // BE: single region
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void AnnotateHistoneModifications_SingleRegionPresentMark_GetsCanonicalState()
    {
        // "single region" target: a 1-element input must produce exactly one
        // annotation carrying that mark's canonical Roadmap state.
        var single = new[] { (100, 300, "H3K4me3", 0.9) };

        var result = AnnotateHistoneModifications(single).ToList();

        result.Should().HaveCount(1);
        AssertWellFormedAnnotation(result[0], single[0]);
        result[0].PredictedState.Should().Be(ChromatinState.ActivePromoter,
            "H3K4me3 present is the canonical active-promoter (TssA) mark");
    }

    [Test]
    public void AnnotateHistoneModifications_SingleRegion_AllRecognisedMarks_MapToDocumentedState()
    {
        // One region per recognised mark, each its OWN single-element call, so we
        // exercise the 1-element boundary for every documented branch.
        var expectations = new (string Mark, ChromatinState State)[]
        {
            ("H3K4me3", ChromatinState.ActivePromoter),
            ("H3K4me1", ChromatinState.WeakEnhancer),
            ("H3K27ac", ChromatinState.ActiveEnhancer),
            ("H3K36me3", ChromatinState.Transcribed),
            ("H3K27me3", ChromatinState.Repressed),
            ("H3K9me3", ChromatinState.Heterochromatin),
            ("H3K9ac", ChromatinState.ActivePromoter),
        };

        foreach (var (mark, state) in expectations)
        {
            var single = new[] { (0, 200, mark, 0.8) };
            var result = AnnotateHistoneModifications(single).ToList();

            result.Should().HaveCount(1);
            AssertWellFormedAnnotation(result[0], single[0]);
            result[0].PredictedState.Should().Be(state,
                $"{mark} present must map to its canonical Roadmap state");
        }
    }

    [Test]
    public void AnnotateHistoneModifications_MarkNameCaseInsensitive()
    {
        // Source binarizes mark name via ToUpperInvariant — case must not matter.
        var lower = AnnotateHistoneModifications(new[] { (0, 100, "h3k4me3", 0.9) }).Single();
        var upper = AnnotateHistoneModifications(new[] { (0, 100, "H3K4ME3", 0.9) }).Single();

        lower.PredictedState.Should().Be(ChromatinState.ActivePromoter);
        upper.PredictedState.Should().Be(ChromatinState.ActivePromoter);
    }

    // ════════════════════════════════════════════════════════════════════════
    // BE: no marks  (no present mark / quiescent default)
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void PredictChromatinState_NoMarkPresent_IsLowSignal_NotFalseActive()
    {
        // INV-02: all six marks below threshold → quiescent LowSignal, never a
        // false active/regulatory call.
        var state = PredictChromatinState(
            h3k4me3: 0.0, h3k4me1: 0.0, h3k27ac: 0.0,
            h3k36me3: 0.0, h3k27me3: 0.0, h3k9me3: 0.0);

        state.Should().Be(ChromatinState.LowSignal);
    }

    [Test]
    public void PredictChromatinState_AllJustBelowThreshold_IsLowSignal()
    {
        // Boundary: signals strictly below the inclusive presence call (0.5) are
        // ABSENT → no mark present → LowSignal.
        double justBelow = Threshold - 1e-9;
        var state = PredictChromatinState(
            justBelow, justBelow, justBelow, justBelow, justBelow, justBelow);

        state.Should().Be(ChromatinState.LowSignal);
    }

    [Test]
    public void PredictChromatinState_NegativeSignals_TreatedAsAbsent_IsLowSignal()
    {
        // Docs §6.1: negative/zero signals are treated as absent (below the call),
        // never as a present mark and never a crash.
        Action act = () => PredictChromatinState(-1.0, -5.0, double.MinValue, -0.5, -1e9, -2.0);
        act.Should().NotThrow();

        PredictChromatinState(-1.0, -5.0, double.MinValue, -0.5, -1e9, -2.0)
            .Should().Be(ChromatinState.LowSignal);
    }

    [Test]
    public void AnnotateHistoneModifications_RegionBelowThreshold_IsLowSignal()
    {
        // "no marks" at the region level: a region whose mark signal is below the
        // presence call carries NO present mark → LowSignal, not a false active.
        var single = new[] { (0, 500, "H3K4me3", 0.1) };
        var result = AnnotateHistoneModifications(single).ToList();

        result.Should().HaveCount(1);
        AssertWellFormedAnnotation(result[0], single[0]);
        result[0].PredictedState.Should().Be(ChromatinState.LowSignal,
            "an active mark below the presence call must not yield ActivePromoter");
    }

    [Test]
    public void AnnotateHistoneModifications_UnknownMark_IsLowSignal_NoKeyNotFound()
    {
        // KeyNotFound hazard: an unrecognised mark name (even with strong signal)
        // must map to LowSignal, never throw.
        var single = new[] { (0, 100, "H4K20me3", 0.99) };

        Action act = () => AnnotateHistoneModifications(single).ToList();
        act.Should().NotThrow();

        var result = AnnotateHistoneModifications(single).Single();
        AssertWellFormedAnnotation(result, single[0]);
        result.PredictedState.Should().Be(ChromatinState.LowSignal,
            "an unknown mark is not in the Roadmap signature table → LowSignal");
    }

    [Test]
    public void AnnotateHistoneModifications_EmptyMarkName_IsLowSignal_NoCrash()
    {
        // Degenerate mark string (empty) is just an unrecognised name → LowSignal.
        var single = new[] { (0, 100, "", 0.99) };

        Action act = () => AnnotateHistoneModifications(single).ToList();
        act.Should().NotThrow();

        AnnotateHistoneModifications(single).Single().PredictedState
            .Should().Be(ChromatinState.LowSignal);
    }

    // ════════════════════════════════════════════════════════════════════════
    // POSITIVE sanity: documented mark combinations classify correctly
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void PredictChromatinState_KnownCombinations_MatchDocumentedTable()
    {
        // Each row of the §2.2 mark-set → state table, with present marks at a
        // clearly-above-threshold magnitude and absent marks at 0.
        // H3K4me3 alone → ActivePromoter (TssA).
        PredictChromatinState(0.9, 0, 0, 0, 0, 0)
            .Should().Be(ChromatinState.ActivePromoter);

        // H3K4me1 + H3K27ac → ActiveEnhancer (active Enh).
        PredictChromatinState(0, 0.9, 0.9, 0, 0, 0)
            .Should().Be(ChromatinState.ActiveEnhancer);

        // H3K4me1 without H3K27ac → WeakEnhancer (Enh).
        PredictChromatinState(0, 0.9, 0, 0, 0, 0)
            .Should().Be(ChromatinState.WeakEnhancer);

        // H3K36me3 → Transcribed (Tx).
        PredictChromatinState(0, 0, 0, 0.9, 0, 0)
            .Should().Be(ChromatinState.Transcribed);

        // H3K27me3 alone → Repressed (ReprPC).
        PredictChromatinState(0, 0, 0, 0, 0.9, 0)
            .Should().Be(ChromatinState.Repressed);

        // H3K9me3 alone → Heterochromatin (Het).
        PredictChromatinState(0, 0, 0, 0, 0, 0.9)
            .Should().Be(ChromatinState.Heterochromatin);

        // H3K4me3 + H3K27me3 → BivalentPromoter (TssBiv) — INV-03 override.
        PredictChromatinState(0.9, 0, 0, 0, 0.9, 0)
            .Should().Be(ChromatinState.BivalentPromoter);

        // H3K4me1 + H3K27me3 → BivalentEnhancer (EnhBiv).
        PredictChromatinState(0, 0.9, 0, 0, 0.9, 0)
            .Should().Be(ChromatinState.BivalentEnhancer);
    }

    [Test]
    public void PredictChromatinState_MarkExactlyAtThreshold_CountsAsPresent()
    {
        // Docs §6.1: an inclusive `>=` call — a mark exactly at the threshold is
        // PRESENT.
        PredictChromatinState(Threshold, 0, 0, 0, 0, 0)
            .Should().Be(ChromatinState.ActivePromoter);
    }

    // ════════════════════════════════════════════════════════════════════════
    // Determinism + total-function (INV-01 / INV-04) over fuzzed inputs
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    [CancelAfter(10000)]
    public void PredictChromatinState_RandomVectors_AreDeterministic_AndAlwaysDefined()
    {
        // INV-01 (equal pattern ⇒ equal state, hence deterministic) + INV-04
        // (total function) over a fuzzed cloud of six-mark vectors, including
        // values straddling the threshold and extreme magnitudes.
        var rng = new Random(20260620);

        for (int i = 0; i < 5000; i++)
        {
            double M() => rng.NextDouble() switch
            {
                < 0.15 => 0.0,                          // absent
                < 0.30 => Threshold,                    // exactly at the call
                < 0.45 => -rng.NextDouble() * 10,       // negative
                < 0.60 => rng.NextDouble(),             // [0,1)
                < 0.75 => 1e12 * rng.NextDouble(),      // huge
                _ => rng.NextDouble() < 0.5 ? double.MinValue : double.MaxValue,
            };

            double a = M(), b = M(), c = M(), d = M(), e = M(), f = M();

            var s1 = PredictChromatinState(a, b, c, d, e, f);
            var s2 = PredictChromatinState(a, b, c, d, e, f);

            s2.Should().Be(s1, "same mark vector must yield the same state (INV-01)");
            DefinedStates.Should().Contain(s1, "result must be a defined state (INV-04)");
        }
    }

    [Test]
    [CancelAfter(10000)]
    public void AnnotateHistoneModifications_RandomRegions_WellFormed_AndDeterministic()
    {
        // Fuzz many regions with a mix of known/unknown marks and present/absent
        // signals; every annotation must be well-formed and order/count-preserving.
        var rng = new Random(424242);
        string[] marks =
        {
            "H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3", "H3K27me3",
            "H3K9me3", "H3K9ac", "H4K20me3", "UNKNOWN", "",
        };

        for (int iter = 0; iter < 300; iter++)
        {
            int n = rng.Next(0, 8); // includes the empty (0-region) boundary
            var input = new (int, int, string, double)[n];
            for (int k = 0; k < n; k++)
            {
                int start = rng.Next(0, 1_000_000);
                int end = start + rng.Next(1, 1000);
                string mark = marks[rng.Next(marks.Length)];
                double signal = rng.NextDouble() < 0.5 ? rng.NextDouble() : rng.NextDouble() + 0.5;
                input[k] = (start, end, mark, signal);
            }

            var r1 = AnnotateHistoneModifications(input).ToList();
            var r2 = AnnotateHistoneModifications(input).ToList();

            r1.Should().HaveCount(n, "annotation is one-to-one with input regions");
            for (int k = 0; k < n; k++)
            {
                AssertWellFormedAnnotation(r1[k], input[k]);
                r2[k].PredictedState.Should().Be(r1[k].PredictedState,
                    "annotation must be deterministic");
            }
        }
    }

    #endregion
}
