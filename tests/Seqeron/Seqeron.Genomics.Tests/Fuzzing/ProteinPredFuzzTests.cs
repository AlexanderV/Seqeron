using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.DisorderPredictor;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinPred area — intrinsic disorder PREDICTION
/// (the per-residue TOP-IDP disorder profile, DISORDER-PRED-001), the
/// segmentation of that profile into contiguous disordered REGIONS
/// (DISORDER-REGION-001), and SEG low-complexity region detection
/// (DISORDER-LC-001, <see cref="DisorderPredictor.PredictLowComplexityRegions"/>).
///
/// NOTE — naming: a sibling fixture <c>ProteinLowComplexityFuzzTests</c> covers a
/// DIFFERENT method/unit (<c>ProteinMotifFinder.FindLowComplexityRegions</c>,
/// PROTMOTIF-LC-001, checklist row 165). DISORDER-LC-001 below targets the SEG
/// detector that lives on <see cref="DisorderPredictor"/> and returns
/// (int Start, int End, string Type) — there is NO numeric Complexity field — so it
/// is co-located with the other DisorderPredictor disorder units in this file.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (NaN/±∞ scores), and no
/// *unhandled* runtime exception. For a per-residue propensity scorer that looks
/// each amino acid up in a 20-entry scale and averages a clipped sliding window,
/// the headline hazards are:
///   • a KeyNotFoundException when a NON-STANDARD amino acid (B, Z, J, O, U) or the
///     unknown placeholder 'X' is fed into the propensity table;
///   • a DivideByZeroException / NaN when the sequence is empty or a window contains
///     no recognised residue (an all-X / all-junk window — count == 0);
///   • an IndexOutOfRangeException in the −w/2..+w/2 window extraction when the
///     sequence is shorter than the window (a 1-residue protein).
/// Every input must resolve to EITHER a well-defined, theory-correct result, OR a
/// *documented, intentional* outcome (here: an empty result for empty input; a
/// score of 0.0 for an unrecognised-only window).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: DISORDER-PRED-001 — intrinsic disorder prediction (per-residue profile)
/// Checklist: docs/checklists/03_FUZZING.md, row 80.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length/edge corners that could crash or
///     produce a non-finite score:
///       – empty protein: "" / null → the EMPTY DisorderPredictionResult (no
///         residue predictions, no regions, 0 content, 0 mean score) by an explicit
///         early return, never a DivideByZero on the length-0 `Average()`/division.
///       – 1 AA: a single residue → exactly one ResiduePrediction whose clipped
///         window is the residue itself; the half-window clamp (start = max(0,…),
///         end = min(n,…)) must NOT run off the 1-char string (no IndexOutOfRange).
///       – all-same AA: a homopolymer → every interior window is identical, so all
///         per-residue scores are UNIFORM (the anchor values W→0.0, P→1.0 pin the
///         normalisation endpoints, Campen et al. 2008).
///   • MC = Malformed Content — out-of-table residues:
///       – non-standard AAs (B, Z, J, O, U) and 'X' (unknown): absent from the
///         20-entry TOP-IDP scale. The scorer looks them up with `TryGetValue` and
///         SKIPS misses, so they must be HANDLED deterministically — never a
///         KeyNotFoundException. A window of nothing-but-unknowns scores 0.0
///         (count == 0 branch), never NaN.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC);
///   targets: "Empty protein, 1 AA, all same AA, non-standard AAs, X amino acid".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The disorder-prediction contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Disorder prediction estimates intrinsic disorder from amino-acid composition
/// using the TOP-IDP propensity scale and a sliding window. For residue i the
/// score is the average of NORMALISED TOP-IDP values over the local window W_i:
///     S_i = (1/|W_i|) · Σ_{c ∈ W_i} (p(c) − p_min) / (p_max − p_min),
/// with p_min = −0.884 (Trp) and p_max = 0.987 (Pro); a residue is labelled
/// disordered when S_i ≥ 0.542 (the default TOP-IDP cutoff). Because each term is
/// clamped between the documented min and max, every per-residue score lies in
/// [0, 1] (INV-01). — docs/algorithms/ProteinPred/Disorder_Prediction.md §2.2, §2.4.
///
/// Method under test (src/.../Seqeron.Genomics.Analysis/DisorderPredictor.cs):
///   PredictDisorder(string sequence, int windowSize = 21, double disorderThreshold
///       = 0.542, int minRegionLength = 5) → DisorderPredictionResult
///       { Sequence, ResiduePredictions, DisorderedRegions, OverallDisorderContent,
///         MeanDisorderScore }.
///   Each ResiduePrediction = { Position (0-based), Residue, DisorderScore, IsDisordered }.
///   — Disorder_Prediction.md §3.2, §4.1.
///
/// Documented input handling (Disorder_Prediction.md §3.1, §3.3, §5.4, §6.1):
///   • null / "" → empty DisorderPredictionResult (Sequence "", no predictions/regions,
///     0 content, 0 mean) — explicit short-circuit, never throws, never DivideByZero.
///   • The sequence is uppercased first; lowercase == uppercase.
///   • Characters ABSENT from the TOP-IDP table (non-standard B/Z/J/O/U, X, digits,
///     punctuation) are PRESERVED in ResiduePrediction.Residue but EXCLUDED from the
///     window average (`TryGetValue` skips misses — DEVIATION #1, accepted). A window
///     with no recognised residue scores 0.0 (the count == 0 branch). No KeyNotFound.
///   • Edge windows are clipped to the sequence bounds, not padded, so terminal /
///     single residues are scored from fewer than `windowSize` positions — no
///     IndexOutOfRange.
///
/// Theory-correct invariants asserted (Disorder_Prediction.md §2.4):
///   • INV-01 — every per-residue DisorderScore is in [0, 1] (and FINITE — never NaN/±∞).
///   • INV-02 — the prediction is deterministic for a fixed input (re-running PredictDisorder
///     on the same sequence yields identical scores).
///   • INV-03 — IsDisordered == (DisorderScore ≥ threshold) for the active threshold.
///   • [residue-aligned] — ResiduePredictions has exactly one entry per input residue,
///     Position is the 0-based index, and Residue is the uppercased input character.
///   • [summary-finite] — OverallDisorderContent and MeanDisorderScore are finite and in
///     [0, 1] (mean of in-range scores; content is a fraction).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// PredictDisorder is O(n·w) over a clipped window with O(n) state
/// (Disorder_Prediction.md §4.3); there is no data-dependent branching that could
/// loop. The homopolymer and long-junk targets maximise the per-window work; they
/// are kept modest and [CancelAfter]-guarded so a regression that turned the scan
/// into a hang or super-linear blow-up would FAIL rather than wedge the suite.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: DISORDER-REGION-001 — disorder REGION detection (run segmentation)
/// Checklist: docs/checklists/03_FUZZING.md, row 81.
/// ═══════════════════════════════════════════════════════════════════════════
/// Region detection turns the per-residue disorder CALLS (ResiduePrediction.IsDisordered,
/// itself score ≥ disorderThreshold — DISORDER-PRED-001) into contiguous disordered
/// REGIONS: a single forward scan emits every maximal run of consecutive disordered
/// residues whose run length is ≥ minRegionLength (default 5). Each emitted
/// DisorderedRegion = { Start (0-based incl.), End (0-based incl.), MeanScore (mean of
/// the run's scores), Confidence ∈ [0,1], RegionType }. The run boundaries are derived
/// from the ALREADY-COMPUTED IsDisordered flags; the threshold argument is forwarded but
/// the scan itself only reads the flags (Disordered_Region_Detection.md §5.2).
///   — entry point DisorderPredictor.PredictDisorder(sequence, windowSize=21,
///     disorderThreshold=0.542, minRegionLength=5) → DisorderPredictionResult.DisorderedRegions
///     (Disordered_Region_Detection.md §3.1, §3.2, §4.1).
///
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the threshold/composition corners that probe the
///     run-segmentation boundaries and could crash (IndexOutOfRange merging adjacent
///     residues), hang, or emit a malformed region:
///       – disorderThreshold = 0: a residue is disordered iff score ≥ 0, and every
///         normalised score lies in [0, 1] (INV-01 of DISORDER-PRED-001) → EVERY residue
///         is disordered REGARDLESS of composition → exactly ONE region spanning the whole
///         protein [0, n−1] (for n ≥ minRegionLength), even for an order-promoting (all-W)
///         protein. The boundary is INCLUSIVE (`score >= threshold`).
///       – disorderThreshold = 1: a residue is disordered iff score ≥ 1.0; only the
///         saturating maximum (a window of pure P, the TOP-IDP max p=0.987 → normalised
///         1.0) reaches it. A protein that does NOT saturate to 1.0 everywhere → NO
///         residue is disordered → NO regions (the upper boundary). The saturating
///         exception (all-P) still yields ONE full-span region because P pins exactly 1.0.
///       – all ordered (all scores below threshold, e.g. all-W at the default cutoff) →
///         NO regions (Disordered_Region_Detection.md §6.1).
///       – all disordered (all scores above threshold, e.g. all-P at the default cutoff) →
///         exactly ONE region spanning the whole protein [0, n−1]
///         (Disordered_Region_Detection.md §6.1).
///       – empty: "" / null → the empty result short-circuits before any scan → NO regions,
///         no IndexOutOfRange, no DivideByZero (Disordered_Region_Detection.md §3.3).
///   — docs/checklists/03_FUZZING.md §Description (strategy code BE);
///     targets: "threshold=0, threshold=1, all ordered, all disordered, empty".
///
/// Theory-correct region invariants asserted (Disordered_Region_Detection.md §2.4):
///   • INV-01 — emitted regions are NON-OVERLAPPING and sorted by increasing Start.
///   • INV-02 — every emitted region has (End − Start + 1) ≥ minRegionLength.
///   • [in-bounds] — 0 ≤ Start ≤ End ≤ n−1 (no run-off-the-end IndexOutOfRange when merging
///     adjacent residues into a region).
///   • [run-disordered] — every residue inside an emitted region is IsDisordered, and the
///     residue immediately before Start / after End (when present) is NOT — i.e. the run is
///     a MAXIMAL contiguous disordered stretch.
///   • INV-03 — Confidence ∈ [0, 1]; MeanScore is finite and in [0, 1].
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinPredFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino acids that ARE in the TOP-IDP propensity table.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Default SEG low-complexity parameters (NCBI blast_seg.c: kSegWindow / kSegLocut /
    /// kSegHicut) for DISORDER-LC-001 — Low_Complexity_Region_Detection.md §4.2.</summary>
    private const int SegWindow = 12;   // W (kSegWindow)
    private const double SegK1 = 2.2;   // trigger complexity K1 (kSegLocut), bits/residue
    private const double SegK2 = 2.5;   // extension complexity K2 (kSegHicut), bits/residue

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[rng.Next(StandardAminoAcids.Length)];
        return new string(chars);
    }

    /// <summary>
    /// A maximally diverse mosaic for the SEG low-complexity unit: residue i is the
    /// (i mod 20)-th standard amino acid. Every width-12 window holds 12 DISTINCT residues,
    /// so its Shannon entropy is log₂(12) ≈ 3.585 &gt; the SEG extension cutoff K2 = 2.5 — a
    /// high-complexity sequence with NO low-complexity region
    /// (Low_Complexity_Region_Detection.md §2.4 INV-04, §6.1).
    /// </summary>
    private static string DiverseProtein(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[i % StandardAminoAcids.Length];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every <see cref="DisorderPredictionResult"/>
    /// must satisfy (Disorder_Prediction.md §2.4, §3.2): one ResiduePrediction per residue with a
    /// 0-based Position and the UPPERCASED input character (INV [residue-aligned]); every
    /// DisorderScore FINITE and in [0, 1] (INV-01 — the headline "never NaN, never out of range"
    /// property); IsDisordered == (score ≥ threshold) (INV-03); and finite summary statistics
    /// (OverallDisorderContent and MeanDisorderScore in [0, 1], [summary-finite]).
    /// </summary>
    private static void AssertWellFormedResult(
        DisorderPredictionResult result, string sequence, double threshold = 0.542)
    {
        string normalized = sequence.ToUpperInvariant();

        // Sequence is the uppercased input.
        result.Sequence.Should().Be(normalized, "the result echoes the uppercased input sequence");

        // [residue-aligned] one prediction per residue, in order, on the uppercased char.
        result.ResiduePredictions.Should().HaveCount(normalized.Length,
            "there is exactly one per-residue prediction per input residue");
        for (int i = 0; i < normalized.Length; i++)
        {
            var pred = result.ResiduePredictions[i];
            pred.Position.Should().Be(i, "ResiduePrediction.Position is the 0-based residue index");
            pred.Residue.Should().Be(normalized[i], "ResiduePrediction.Residue is the uppercased input character");

            // INV-01 — the score is finite and normalised to [0, 1].
            double.IsNaN(pred.DisorderScore).Should().BeFalse(
                $"a per-residue score must never be NaN (position {i}, residue '{pred.Residue}')");
            double.IsInfinity(pred.DisorderScore).Should().BeFalse(
                $"a per-residue score must never be infinite (position {i}, residue '{pred.Residue}')");
            pred.DisorderScore.Should().BeInRange(0.0, 1.0,
                $"INV-01: per-residue score is normalised to [0, 1] (position {i})");

            // INV-03 — disorder flag is exactly the cutoff comparison.
            pred.IsDisordered.Should().Be(pred.DisorderScore >= threshold,
                "INV-03: IsDisordered == (DisorderScore >= threshold)");
        }

        // [summary-finite] content + mean score are finite and in [0, 1].
        result.OverallDisorderContent.Should().BeInRange(0.0, 1.0,
            "OverallDisorderContent is a fraction in [0, 1]");
        double.IsNaN(result.MeanDisorderScore).Should().BeFalse("MeanDisorderScore must never be NaN");
        if (normalized.Length > 0)
            result.MeanDisorderScore.Should().BeInRange(0.0, 1.0,
                "MeanDisorderScore is the mean of in-range scores, hence in [0, 1]");
    }

    /// <summary>
    /// Asserts the theory-correct REGION contract every emitted <see cref="DisorderedRegion"/>
    /// must satisfy (DISORDER-REGION-001; Disordered_Region_Detection.md §2.4, §3.2): every
    /// region is in-bounds (0 ≤ Start ≤ End ≤ n−1, so merging adjacent residues never runs off
    /// the end — no IndexOutOfRange); each run length is ≥ minRegionLength (INV-02); regions are
    /// non-overlapping and sorted by increasing Start (INV-01); every residue inside a region is
    /// IsDisordered and the run is MAXIMAL (the residue immediately before/after is not disordered)
    /// ([run-disordered]); and Confidence ∈ [0, 1] with a finite in-range MeanScore (INV-03).
    /// </summary>
    private static void AssertWellFormedRegions(
        DisorderPredictionResult result, int minRegionLength = 5)
    {
        var predictions = result.ResiduePredictions;
        int n = predictions.Count;
        int prevEnd = -1;

        foreach (var region in result.DisorderedRegions)
        {
            // [in-bounds] — coordinates never escape the sequence (no IndexOutOfRange on merge).
            region.Start.Should().BeInRange(0, Math.Max(0, n - 1),
                "a region Start is a valid 0-based residue index");
            region.End.Should().BeInRange(region.Start, Math.Max(0, n - 1),
                "a region End is in-bounds and not before its Start");

            // INV-02 — run length meets the minimum.
            (region.End - region.Start + 1).Should().BeGreaterThanOrEqualTo(minRegionLength,
                "INV-02: an emitted region spans at least minRegionLength residues");

            // INV-01 — non-overlapping and strictly after the previous region.
            region.Start.Should().BeGreaterThan(prevEnd,
                "INV-01: regions are non-overlapping and sorted by increasing Start");
            prevEnd = region.End;

            // [run-disordered] — every residue in the run is disordered; the run is maximal.
            for (int i = region.Start; i <= region.End; i++)
                predictions[i].IsDisordered.Should().BeTrue(
                    $"every residue inside an emitted region is disordered (position {i})");
            if (region.Start > 0)
                predictions[region.Start - 1].IsDisordered.Should().BeFalse(
                    "the residue immediately before a region is ordered (maximal run)");
            if (region.End < n - 1)
                predictions[region.End + 1].IsDisordered.Should().BeFalse(
                    "the residue immediately after a region is ordered (maximal run)");

            // INV-03 — finite, in-range mean and clamped confidence.
            double.IsNaN(region.MeanScore).Should().BeFalse("a region MeanScore must never be NaN");
            region.MeanScore.Should().BeInRange(0.0, 1.0, "a region MeanScore is a mean of in-range scores");
            region.Confidence.Should().BeInRange(0.0, 1.0, "INV-03: region Confidence is clamped to [0, 1]");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  DISORDER-PRED-001 — intrinsic disorder prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region DISORDER-PRED-001 — intrinsic disorder prediction

    #region BE — Empty / null protein: empty result, no DivideByZero

    /// <summary>
    /// Target "Empty protein": "" and null must produce the EMPTY DisorderPredictionResult —
    /// Sequence "", no residue predictions, no regions, OverallDisorderContent 0,
    /// MeanDisorderScore 0 — by the explicit short-circuit, NEVER an exception. This is the
    /// headline no-DivideByZero contract: the length-0 path must not hit the
    /// `disorderedCount / sequence.Length` division or the empty `Average()` that would throw
    /// (Disorder_Prediction.md §3.1, §6.1 "null or empty sequence").
    /// </summary>
    [Test]
    public void PredictDisorder_EmptyOrNull_ReturnsEmptyResult()
    {
        foreach (string? seq in new[] { "", null })
        {
            var act = () => PredictDisorder(seq!);
            var result = act.Should().NotThrow($"empty/null input ('{seq ?? "null"}') must not crash").Subject;

            result.Sequence.Should().BeEmpty("empty/null input yields an empty Sequence");
            result.ResiduePredictions.Should().BeEmpty("empty/null input yields no per-residue predictions");
            result.DisorderedRegions.Should().BeEmpty("empty/null input yields no disordered regions");
            result.OverallDisorderContent.Should().Be(0.0, "empty/null input yields zero disorder content (no DivideByZero)");
            result.MeanDisorderScore.Should().Be(0.0, "empty/null input yields zero mean score (no empty-Average crash)");
        }
    }

    #endregion

    #region BE — Single residue: one score, no IndexOutOfRange on the clipped window

    /// <summary>
    /// Target "1 AA": a single-residue protein must produce exactly ONE ResiduePrediction.
    /// The window clamp (start = max(0, i − w/2), end = min(n, i + w/2 + 1)) collapses to the
    /// whole 1-char string, so the −w/2..+w/2 extraction must NOT run off either end (no
    /// IndexOutOfRange) even though the sequence is far shorter than the default window of 21
    /// (Disorder_Prediction.md §5.2 "clips edge windows to the available sequence bounds").
    /// For a recognised residue the single-residue score equals its normalised TOP-IDP value:
    /// W → 0.0 (the minimum) and P → 1.0 (the maximum), pinning the normalisation endpoints
    /// (Campen et al. 2008; Disorder_Prediction.md §7.1).
    /// </summary>
    [Test]
    public void PredictDisorder_SingleResidue_OneScoreNoIndexOutOfRange()
    {
        // Every standard amino acid, one at a time — never throws, exactly one prediction.
        foreach (char aa in StandardAminoAcids)
        {
            string seq = aa.ToString();
            var act = () => PredictDisorder(seq);
            var result = act.Should().NotThrow($"a 1-residue protein ('{aa}') must not crash on the clipped window").Subject;

            result.ResiduePredictions.Should().ContainSingle("a 1-residue protein has exactly one per-residue prediction");
            AssertWellFormedResult(result, seq);

            // The single-residue score is exactly the residue's normalised TOP-IDP value.
            double expected = (GetDisorderPropensity(aa) - (-0.884)) / (0.987 - (-0.884));
            result.ResiduePredictions[0].DisorderScore.Should().BeApproximately(expected, 1e-9,
                $"a 1-residue window scores the residue's normalised TOP-IDP value ('{aa}')");
        }

        // Pin the documented normalisation endpoints on single residues.
        PredictDisorder("W").ResiduePredictions[0].DisorderScore.Should().BeApproximately(0.0, 1e-9,
            "W is the TOP-IDP minimum → normalised score 0.0 (Campen et al. 2008)");
        PredictDisorder("P").ResiduePredictions[0].DisorderScore.Should().BeApproximately(1.0, 1e-9,
            "P is the TOP-IDP maximum → normalised score 1.0 (Campen et al. 2008)");
    }

    #endregion

    #region BE — All-same AA (homopolymer): uniform scores

    /// <summary>
    /// Target "all same AA": a homopolymer "XXXX…" of a RECOGNISED residue has identical
    /// composition in every interior window, so every interior per-residue score is UNIFORM
    /// (and equals the residue's normalised TOP-IDP value, since the window averages copies of
    /// the same value). Terminal residues see a shorter clipped window but the SAME residue, so
    /// they share the identical value too — the whole profile is flat. This proves the scorer is
    /// composition-driven and that the clipped edge windows do not perturb a uniform profile
    /// (Disorder_Prediction.md §2.2, §5.2). Anchors: all-W → 0.0, all-P → 1.0, all-E → ≈0.866
    /// (above the 0.542 cutoff → disordered) (Disorder_Prediction.md §7.1).
    /// </summary>
    [Test]
    public void PredictDisorder_Homopolymer_ProducesUniformScores()
    {
        foreach (char aa in StandardAminoAcids)
        {
            string seq = new string(aa, 40);
            var act = () => PredictDisorder(seq);
            var result = act.Should().NotThrow($"a homopolymer of '{aa}' must not crash").Subject;

            AssertWellFormedResult(result, seq);

            // Every residue shares the same normalised TOP-IDP value → a flat profile.
            double expected = (GetDisorderPropensity(aa) - (-0.884)) / (0.987 - (-0.884));
            result.ResiduePredictions.Should().OnlyContain(
                p => Math.Abs(p.DisorderScore - expected) < 1e-9,
                $"a homopolymer of '{aa}' produces a uniform per-residue score");
            result.MeanDisorderScore.Should().BeApproximately(expected, 1e-9,
                $"the mean of a flat profile equals the per-residue value ('{aa}')");
        }

        // Documented homopolymer anchors (Disorder_Prediction.md §7.1).
        PredictDisorder(new string('W', 40)).ResiduePredictions.Should().OnlyContain(
            p => Math.Abs(p.DisorderScore - 0.0) < 1e-9, "all-W scores the normalised minimum 0.0");
        PredictDisorder(new string('P', 40)).OverallDisorderContent.Should().Be(1.0,
            "all-P scores 1.0 everywhere → every residue is disordered (content 1.0)");
        var allE = PredictDisorder(new string('E', 40));
        allE.ResiduePredictions[20].DisorderScore.Should().BeApproximately(0.866, 1e-3,
            "all-E scores ≈0.866 (above the 0.542 cutoff)");
        allE.ResiduePredictions[20].IsDisordered.Should().BeTrue("0.866 ≥ 0.542 → disordered");
    }

    #endregion

    #region MC — Non-standard amino acids (B, Z, J, O, U): handled, never KeyNotFound

    /// <summary>
    /// Target "non-standard AAs": the ambiguity / extended IUPAC amino-acid codes that are NOT in
    /// the 20-entry TOP-IDP scale — B (Asx), Z (Glx), J (Leu/Ile), O (pyrrolysine), U
    /// (selenocysteine) — must be HANDLED deterministically, NEVER raising a KeyNotFoundException.
    /// The scorer looks each residue up with `TryGetValue` and SKIPS misses (DEVIATION #1,
    /// accepted), so a window built only from non-standard residues averages over ZERO recognised
    /// residues and scores 0.0 (the count == 0 branch) — finite, never NaN, never a crash
    /// (Disorder_Prediction.md §3.3, §5.4, §6.1 "window containing only unknown residues"). The
    /// residue character is nonetheless PRESERVED in ResiduePrediction.Residue.
    /// </summary>
    [Test]
    public void PredictDisorder_NonStandardAminoAcids_HandledNotKeyNotFound()
    {
        // All-non-standard homopolymers and a mixed non-standard string.
        foreach (string seq in new[]
                 {
                     new string('B', 30), // Asx
                     new string('Z', 30), // Glx
                     new string('J', 30), // Leu/Ile
                     new string('O', 30), // pyrrolysine
                     new string('U', 30), // selenocysteine
                     "BZJOUBZJOUBZJOUBZJOU", // mixed non-standard
                 })
        {
            var act = () => PredictDisorder(seq);
            var result = act.Should().NotThrow(
                $"non-standard amino acids ('{seq[..Math.Min(10, seq.Length)]}…') must be handled, not KeyNotFound").Subject;

            AssertWellFormedResult(result, seq);

            // Every window is built only from out-of-table residues → score 0.0 everywhere.
            result.ResiduePredictions.Should().OnlyContain(p => p.DisorderScore == 0.0,
                "a window with no recognised residue scores 0.0 (count == 0 branch), never NaN");
            result.OverallDisorderContent.Should().Be(0.0,
                "0.0 < 0.542 everywhere → no residue is disordered");

            // The non-standard residue character is preserved verbatim (uppercased).
            string normalized = seq.ToUpperInvariant();
            for (int i = 0; i < normalized.Length; i++)
                result.ResiduePredictions[i].Residue.Should().Be(normalized[i],
                    "an unrecognised residue is preserved in ResiduePrediction.Residue, not dropped");
        }
    }

    #endregion

    #region MC — 'X' unknown amino acid (and other junk): handled, never KeyNotFound

    /// <summary>
    /// Target "X amino acid": 'X' is the standard placeholder for an UNKNOWN/ANY residue and is
    /// absent from the TOP-IDP table. Like the non-standard codes it must be HANDLED via
    /// `TryGetValue`/skip, never KeyNotFound: an all-X protein scores 0.0 at every position
    /// (no recognised residue in any window), finite and in range (Disorder_Prediction.md §3.3,
    /// §6.1). An 'X' EMBEDDED in standard residues must not crash and must not poison its
    /// neighbours: a recognised residue's window still averages only the recognised members it
    /// contains, so X simply reduces the effective window size. We also probe pure junk (digits,
    /// punctuation, whitespace, IUPAC nucleotide ambiguity letters) to confirm out-of-domain
    /// characters are uniformly skipped, never thrown.
    /// </summary>
    [Test]
    public void PredictDisorder_XAndJunk_HandledNotKeyNotFound()
    {
        // (a) All-X → 0.0 everywhere, finite, no KeyNotFound.
        var allX = ((Func<DisorderPredictionResult>)(() => PredictDisorder(new string('X', 30))))
            .Should().NotThrow("an all-X protein must be handled, not KeyNotFound").Subject;
        AssertWellFormedResult(allX, new string('X', 30));
        allX.ResiduePredictions.Should().OnlyContain(p => p.DisorderScore == 0.0,
            "an all-X window has no recognised residue → score 0.0");

        // (b) X embedded among standard residues: no crash; X reduces effective window, never poisons to NaN.
        foreach (string seq in new[]
                 {
                     "XPPPPPPPPPPX",       // X flanking a poly-P core
                     "GXGXGXGXGXGX",       // X interleaved with recognised G
                     "PEKSXXXXSKEP",       // X core flanked by disorder-promoting residues
                 })
        {
            var result = ((Func<DisorderPredictionResult>)(() => PredictDisorder(seq)))
                .Should().NotThrow($"X embedded in '{seq}' must be handled, not crash").Subject;
            AssertWellFormedResult(result, seq);
        }

        // (c) Pure junk — digits, punctuation, whitespace, nucleotide-style ambiguity letters — all skipped.
        foreach (string junk in new[]
                 {
                     "1234567890!@#",      // digits + punctuation
                     "   \t  \n  ",        // whitespace only
                     "BJOUXZ?-*.",         // mixed extended/junk
                 })
        {
            var result = ((Func<DisorderPredictionResult>)(() => PredictDisorder(junk)))
                .Should().NotThrow($"pure junk ('{junk}') must be handled, not crash").Subject;
            AssertWellFormedResult(result, junk);
            result.ResiduePredictions.Should().OnlyContain(p => p.DisorderScore == 0.0,
                "junk windows have no recognised residue → score 0.0, never NaN");
        }
    }

    #endregion

    #region Positive sanity — disordered vs ordered stretch + random determinism

    /// <summary>
    /// Positive sanity: the harness must assert against a predictor that actually DISCRIMINATES
    /// order from disorder, not a no-op. A disorder-promoting low-complexity stretch dominated by
    /// P/E/S/K (the most disorder-promoting residues, Dunker et al. 2001; the top of the TOP-IDP
    /// ranking, Campen et al. 2008) must score HIGHER on average — and be classified disordered —
    /// than a hydrophobic, order-promoting globular stretch dominated by W/F/I/L/V/C/Y. This pins
    /// the documented direction of the scale (W → 0.0 ordered end, P → 1.0 disordered end) on a
    /// realistic pair of sequences, with every score still finite and in [0, 1].
    /// </summary>
    [Test]
    public void PredictDisorder_DisorderedVsOrdered_ScoresHigherAndIsDisordered()
    {
        const string disordered = "PEPSEKSPEKPESKPESKPESKPESKPESKPESKPESKPE"; // P/E/S/K-rich, low complexity
        const string ordered = "WFILVCYWFILVCYWFILVCYWFILVCYWFILVCYWFILVC"; // hydrophobic, order-promoting

        var dis = PredictDisorder(disordered);
        var ord = PredictDisorder(ordered);

        AssertWellFormedResult(dis, disordered);
        AssertWellFormedResult(ord, ordered);

        dis.MeanDisorderScore.Should().BeGreaterThan(ord.MeanDisorderScore,
            "a P/E/S/K-rich low-complexity stretch is more disordered than a hydrophobic globular stretch");
        dis.MeanDisorderScore.Should().BeGreaterThan(0.542,
            "the disorder-promoting stretch scores above the TOP-IDP cutoff on average");
        ord.MeanDisorderScore.Should().BeLessThan(0.542,
            "the order-promoting hydrophobic stretch scores below the TOP-IDP cutoff on average");
        dis.OverallDisorderContent.Should().BeGreaterThan(ord.OverallDisorderContent,
            "more residues cross the disorder cutoff in the disordered stretch");
    }

    /// <summary>
    /// Positive sanity over RANDOM protein sequences: across fixed seeds and lengths the scan must
    /// never crash, hang or emit a non-finite score, and every result must satisfy the full
    /// contract (one in-range finite score per residue, finite summaries). INV-02 determinism is
    /// pinned by re-running PredictDisorder on the same input and requiring identical scores. This
    /// pins finiteness and window-safety on arbitrary sequences, not just hand-built motifs.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictDisorder_RandomProtein_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 7, 31, 137, 2026 })
        {
            foreach (int len in new[] { 1, 5, 21, 60, 250 })
            {
                string seq = RandomProtein(len, seed);

                var act = () => PredictDisorder(seq);
                var result = act.Should().NotThrow($"random protein must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                AssertWellFormedResult(result, seq);

                // INV-02 — deterministic: the same input yields identical scores.
                var again = PredictDisorder(seq);
                again.ResiduePredictions.Select(p => p.DisorderScore)
                    .Should().Equal(result.ResiduePredictions.Select(p => p.DisorderScore),
                        "INV-02: PredictDisorder is deterministic for a fixed input");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  DISORDER-REGION-001 — disorder region detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region DISORDER-REGION-001 — disorder region detection

    #region BE — disorderThreshold = 0: every residue disordered → one full-span region

    /// <summary>
    /// Target "threshold=0": with disorderThreshold = 0 a residue is disordered iff its
    /// normalised score ≥ 0, and EVERY normalised TOP-IDP score lies in [0, 1] (INV-01 of
    /// DISORDER-PRED-001), so EVERY residue is disordered REGARDLESS of composition. The lower
    /// boundary is inclusive (`score >= threshold`), so even the most order-promoting protein
    /// (all-W, which scores exactly 0.0) is fully disordered. The run scan therefore merges every
    /// residue into a SINGLE region spanning the whole protein [0, n−1] — for any sequence whose
    /// length meets minRegionLength — and must never run off the end on that full-width merge
    /// (no IndexOutOfRange). (Disordered_Region_Detection.md §2.2, §5.2.)
    /// </summary>
    [Test]
    public void PredictDisorder_ThresholdZero_OneRegionSpanningWholeProtein()
    {
        // Diverse compositions, including order-promoting all-W and a 1-residue-too-short case.
        foreach (string seq in new[]
                 {
                     new string('W', 30),                 // order-promoting → scores 0.0, still ≥ 0
                     new string('P', 30),                 // disorder-promoting → scores 1.0
                     RandomProtein(60, 41),               // arbitrary mixed protein
                     "WFILVCYWFILVCYWFILVCYWFILVCYWF",     // hydrophobic / ordered stretch
                 })
        {
            var act = () => PredictDisorder(seq, disorderThreshold: 0.0, minRegionLength: 5);
            var result = act.Should().NotThrow(
                $"threshold=0 must not crash on the full-width merge ('{seq[..Math.Min(8, seq.Length)]}…')").Subject;

            AssertWellFormedResult(result, seq, threshold: 0.0);
            AssertWellFormedRegions(result, minRegionLength: 5);

            // Every residue is disordered at threshold 0 → exactly one full-span region.
            result.ResiduePredictions.Should().OnlyContain(p => p.IsDisordered,
                "score ≥ 0 holds for every residue → all disordered at threshold=0");
            result.DisorderedRegions.Should().ContainSingle(
                "threshold=0 merges every residue into one region");
            result.DisorderedRegions[0].Start.Should().Be(0, "the single region starts at residue 0");
            result.DisorderedRegions[0].End.Should().Be(seq.Length - 1,
                "the single region spans to the last residue");
        }

        // A protein shorter than minRegionLength yields NO region even though all residues are disordered.
        var tooShort = PredictDisorder("PEKS", disorderThreshold: 0.0, minRegionLength: 5);
        tooShort.ResiduePredictions.Should().OnlyContain(p => p.IsDisordered,
            "all 4 residues are disordered at threshold=0");
        tooShort.DisorderedRegions.Should().BeEmpty(
            "a 4-residue run is below minRegionLength=5 → no region emitted (INV-02)");
    }

    #endregion

    #region BE — disorderThreshold = 1: non-saturating protein → no regions (upper boundary)

    /// <summary>
    /// Target "threshold=1": with disorderThreshold = 1 a residue is disordered iff its normalised
    /// score ≥ 1.0, which only the SATURATING maximum reaches — a window of pure P (the TOP-IDP
    /// maximum p = 0.987 → normalised 1.0). A protein that does NOT saturate to 1.0 everywhere has
    /// NO disordered residue → NO regions: this is the UPPER boundary of the threshold sweep
    /// (mirroring the lower boundary in the threshold=0 test). The saturating exception (all-P)
    /// still pins exactly 1.0 ≥ 1.0 and is therefore the lone composition that yields a full-span
    /// region even at threshold=1 — asserted explicitly so the boundary is pinned from both sides.
    /// (Disordered_Region_Detection.md §5.2; Campen et al. 2008 TOP-IDP max = P.)
    /// </summary>
    [Test]
    public void PredictDisorder_ThresholdOne_NoRegionsUnlessSaturated()
    {
        // Non-saturating proteins: no window is pure-P → no score reaches 1.0 → no regions.
        foreach (string seq in new[]
                 {
                     new string('W', 30),                 // minimum score 0.0
                     new string('E', 30),                 // ≈0.866 < 1.0 (disorder-promoting but not max)
                     RandomProtein(60, 53),               // arbitrary mixed protein
                     "PEPSEKSPEKPESKPESKPESKPESKPESKPE",   // disorder-rich but P never fills a whole window
                 })
        {
            var act = () => PredictDisorder(seq, disorderThreshold: 1.0, minRegionLength: 5);
            var result = act.Should().NotThrow(
                $"threshold=1 must not crash ('{seq[..Math.Min(8, seq.Length)]}…')").Subject;

            AssertWellFormedResult(result, seq, threshold: 1.0);
            AssertWellFormedRegions(result, minRegionLength: 5);

            result.ResiduePredictions.Should().NotContain(p => p.IsDisordered,
                "no window saturates to 1.0 → no residue clears the threshold=1 boundary");
            result.DisorderedRegions.Should().BeEmpty(
                "threshold=1 leaves no disordered residue in a non-saturating protein → no regions");
        }

        // Saturating exception: all-P scores exactly 1.0 everywhere → still ONE full-span region.
        var allP = PredictDisorder(new string('P', 30), disorderThreshold: 1.0, minRegionLength: 5);
        AssertWellFormedRegions(allP, minRegionLength: 5);
        allP.ResiduePredictions.Should().OnlyContain(p => p.IsDisordered,
            "all-P scores exactly 1.0 ≥ 1.0 → every residue is disordered even at threshold=1");
        allP.DisorderedRegions.Should().ContainSingle(
            "the saturating all-P protein still yields one region at threshold=1");
        allP.DisorderedRegions[0].Start.Should().Be(0);
        allP.DisorderedRegions[0].End.Should().Be(29);
    }

    #endregion

    #region BE — all ordered: every score below threshold → no regions

    /// <summary>
    /// Target "all ordered": a protein whose per-residue scores are ALL below the active disorder
    /// threshold has no disordered residue, so the run scan opens no region → the empty region
    /// list. The canonical all-ordered anchor is all-W (the TOP-IDP minimum, normalised 0.0, far
    /// under the default 0.542 cutoff); a hydrophobic globular stretch (W/F/I/L/V/C/Y) is also
    /// fully ordered. No region must be emitted and nothing may crash
    /// (Disordered_Region_Detection.md §6.1 "all ordered → no regions").
    /// </summary>
    [Test]
    public void PredictDisorder_AllOrdered_NoRegions()
    {
        foreach (string seq in new[]
                 {
                     new string('W', 30),                 // TOP-IDP minimum, scores 0.0
                     new string('F', 30),                 // strongly order-promoting
                     "WFILVCYWFILVCYWFILVCYWFILVCYWFILVC", // hydrophobic globular stretch
                 })
        {
            var act = () => PredictDisorder(seq, minRegionLength: 5);
            var result = act.Should().NotThrow($"an all-ordered protein must not crash ('{seq[..8]}…')").Subject;

            AssertWellFormedResult(result, seq);
            AssertWellFormedRegions(result, minRegionLength: 5);

            result.ResiduePredictions.Should().NotContain(p => p.IsDisordered,
                "every residue scores below the 0.542 cutoff in an order-promoting protein");
            result.DisorderedRegions.Should().BeEmpty(
                "all-ordered input emits no disordered region");
        }
    }

    #endregion

    #region BE — all disordered: every score above threshold → one full-span region

    /// <summary>
    /// Target "all disordered": a protein whose per-residue scores are ALL above the active disorder
    /// threshold is one contiguous disordered run, so the scan merges every residue into a SINGLE
    /// region spanning the whole protein [0, n−1]. The canonical anchor is all-P (the TOP-IDP
    /// maximum, normalised 1.0, well over the default 0.542 cutoff); all-E (≈0.866) is likewise
    /// fully disordered. Exactly one region with the full-span bounds must be emitted, the trailing
    /// run captured by the end-of-loop handler (Disordered_Region_Detection.md §6.1 "fully
    /// disordered → one region spanning the full sequence"; §2.4 INV-04).
    /// </summary>
    [Test]
    public void PredictDisorder_AllDisordered_OneFullSpanRegion()
    {
        foreach (string seq in new[]
                 {
                     new string('P', 30),                 // TOP-IDP maximum, scores 1.0
                     new string('E', 30),                 // ≈0.866, well above cutoff
                     new string('K', 30),                 // strongly disorder-promoting
                 })
        {
            var act = () => PredictDisorder(seq, minRegionLength: 5);
            var result = act.Should().NotThrow($"an all-disordered protein must not crash ('{seq[..8]}…')").Subject;

            AssertWellFormedResult(result, seq);
            AssertWellFormedRegions(result, minRegionLength: 5);

            result.ResiduePredictions.Should().OnlyContain(p => p.IsDisordered,
                "every residue scores above the 0.542 cutoff in a disorder-promoting homopolymer");
            result.DisorderedRegions.Should().ContainSingle(
                "a fully disordered protein is one contiguous region");
            result.DisorderedRegions[0].Start.Should().Be(0, "the region starts at residue 0");
            result.DisorderedRegions[0].End.Should().Be(seq.Length - 1,
                "the region spans to the last residue (trailing run captured)");

            // INV-04 — a homopolymeric region's MeanScore equals the residue's normalised TOP-IDP score.
            double expected = (GetDisorderPropensity(seq[0]) - (-0.884)) / (0.987 - (-0.884));
            result.DisorderedRegions[0].MeanScore.Should().BeApproximately(expected, 1e-9,
                $"INV-04: a homopolymeric region's mean equals the per-residue value ('{seq[0]}')");
        }
    }

    #endregion

    #region BE — empty / null: no regions, no crash

    /// <summary>
    /// Target "empty": "" and null short-circuit in PredictDisorder before any residue scan or
    /// region merge, so the DisorderedRegions list is empty and nothing throws (no IndexOutOfRange
    /// on a zero-length scan, no DivideByZero) (Disordered_Region_Detection.md §3.3). A whitespace /
    /// junk-only string (no recognised residue) likewise produces all-0.0 scores → no disordered
    /// residue → no regions, confirming the region layer is robust to degenerate upstream input.
    /// </summary>
    [Test]
    public void PredictDisorder_EmptyOrNull_NoRegions()
    {
        foreach (string? seq in new[] { "", null })
        {
            var act = () => PredictDisorder(seq!, minRegionLength: 5);
            var result = act.Should().NotThrow($"empty/null input ('{seq ?? "null"}') must not crash the region scan").Subject;

            result.DisorderedRegions.Should().BeEmpty("empty/null input yields no disordered regions");
            AssertWellFormedRegions(result, minRegionLength: 5);
        }

        // Junk-only input: no recognised residue → all 0.0 scores → no disordered residue → no regions.
        var junk = PredictDisorder("12345 \t XJOU?", minRegionLength: 5);
        junk.DisorderedRegions.Should().BeEmpty(
            "junk scores 0.0 everywhere (< 0.542) → no disordered residue → no regions");
        AssertWellFormedRegions(junk, minRegionLength: 5);
    }

    #endregion

    #region Positive sanity — a disordered stretch flanked by ordered regions yields exactly that region

    /// <summary>
    /// Positive sanity: the region detector must actually SEGMENT — isolate a disordered stretch
    /// embedded between ordered flanks, not flag the whole protein or nothing. An ordered hydrophobic
    /// flank (W/F/I/L/V) + a disorder-promoting P/E/K/S core + an ordered hydrophobic flank must
    /// yield EXACTLY ONE disordered region whose Start/End fall WITHIN the core (the window blur at
    /// the boundaries can shift the exact edge by a few residues, but the region must lie strictly
    /// inside the ordered flanks and cover the bulk of the core). This pins real segmentation
    /// behaviour with correct, in-bounds coordinates — the headline "merge adjacent disordered
    /// residues into a region, bounded by the ordered flanks" contract.
    /// </summary>
    [Test]
    public void PredictDisorder_DisorderedStretchFlankedByOrder_YieldsExactlyThatRegion()
    {
        // 25 ordered + 30 disordered + 25 ordered = 80 residues; core spans indices [25, 54].
        const string orderedFlank = "WFILVWFILVWFILVWFILVWFILV";          // 25 hydrophobic residues
        string disorderedCore = new string('P', 30);                      // 30 disorder-promoting residues
        string seq = orderedFlank + disorderedCore + orderedFlank;        // length 80
        const int coreStart = 25;
        const int coreEnd = 54; // inclusive

        var result = PredictDisorder(seq, minRegionLength: 5);

        AssertWellFormedResult(result, seq);
        AssertWellFormedRegions(result, minRegionLength: 5);

        result.DisorderedRegions.Should().ContainSingle(
            "a single disordered core flanked by order yields exactly one region");
        var region = result.DisorderedRegions[0];

        // The region lies strictly inside the ordered flanks (the flanks are NOT disordered).
        region.Start.Should().BeGreaterThanOrEqualTo(coreStart - 2,
            "the region begins at/near the disordered core, inside the leading ordered flank");
        region.End.Should().BeLessThanOrEqualTo(coreEnd + 2,
            "the region ends at/near the disordered core, inside the trailing ordered flank");
        region.Start.Should().BeLessThan(region.End, "the region is a real multi-residue stretch");

        // The core's central residues are unambiguously inside the region (window fully within the P-core).
        region.Start.Should().BeLessThanOrEqualTo(coreStart + 11,
            "the region covers the disorder-saturated centre of the core");
        region.End.Should().BeGreaterThanOrEqualTo(coreEnd - 11,
            "the region covers the disorder-saturated centre of the core");

        region.MeanScore.Should().BeGreaterThan(0.542, "the emitted region is genuinely disordered on average");
        region.RegionType.Should().Be("Proline-rich", "a poly-P core is classified Proline-rich");
    }

    /// <summary>
    /// Positive sanity over RANDOM proteins and a sweep of thresholds / minRegionLengths: across
    /// fixed seeds the region scan must never crash, hang, or emit a malformed region, and the full
    /// region contract (in-bounds, non-overlapping, sorted, maximal, ≥ minRegionLength, finite mean
    /// and clamped confidence) must hold for every parameter combination. Determinism is pinned by
    /// re-running and requiring identical region bounds. This pins segmentation robustness on
    /// arbitrary sequences and threshold corners, not just hand-built motifs.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictDisorder_RandomProtein_RegionsAlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 11, 73, 199, 2026 })
        {
            foreach (int len in new[] { 1, 6, 21, 80, 200 })
            {
                string seq = RandomProtein(len, seed);

                foreach (double threshold in new[] { 0.0, 0.3, 0.542, 0.8, 1.0 })
                {
                    foreach (int minLen in new[] { 1, 5, 20 })
                    {
                        var act = () => PredictDisorder(seq, disorderThreshold: threshold, minRegionLength: minLen);
                        var result = act.Should().NotThrow(
                            $"region scan must not crash (seed {seed}, len {len}, thr {threshold}, minLen {minLen})").Subject;
                        token.ThrowIfCancellationRequested();

                        AssertWellFormedResult(result, seq, threshold);
                        AssertWellFormedRegions(result, minLen);

                        // Deterministic: same input + params yield identical region bounds.
                        var again = PredictDisorder(seq, disorderThreshold: threshold, minRegionLength: minLen);
                        again.DisorderedRegions.Select(r => (r.Start, r.End))
                            .Should().Equal(result.DisorderedRegions.Select(r => (r.Start, r.End)),
                                "region detection is deterministic for a fixed input and parameters");
                    }
                }
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  DISORDER-LC-001 — SEG low-complexity region detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region DISORDER-LC-001 — SEG low-complexity region detection

    // ───────────────────────────────────────────────────────────────────────────
    //  Unit: DISORDER-LC-001 — low-complexity region detection (SEG).
    //  Checklist: docs/checklists/03_FUZZING.md, row 204 (ProteinPred, strategy BE;
    //    targets: "homopolymer, high-complexity, empty").
    //  Method under test (src/.../Seqeron.Genomics.Analysis/DisorderPredictor.cs):
    //    IEnumerable<(int Start,int End,string Type)> PredictLowComplexityRegions(
    //        string sequence, int triggerWindow = 12, double triggerThreshold = 2.2,
    //        double extensionThreshold = 2.5, int minLength = 1)
    //    — Low_Complexity_Region_Detection.md §3.1, §5.1.
    //
    //  The SEG contract under test (Low_Complexity_Region_Detection.md §2.2, §4):
    //    For a width-W window the local compositional complexity is the Shannon entropy of
    //    its residue composition in bits/residue, H = −Σᵢ pᵢ·log₂ pᵢ (pᵢ = nᵢ/W, max
    //    log₂(20) ≈ 4.322). Stage 1 marks every window with H ≤ K1 (trigger); stage 2 grows
    //    each triggered span while the whole growing segment's H ≤ K2; overlapping/adjacent
    //    spans merge; segments shorter than minLength are dropped; each is emitted as
    //    (Start, End inclusive 0-based, Type) where Type is a convenience composition label
    //    (§5.4 — "X-rich" if the dominant residue fraction > 0.5, else "X/Y-rich").
    //
    //  Documented input handling (§3.3, §6.1):
    //    • null sequence → ArgumentNullException (the ONLY documented throw; on enumeration).
    //    • "" or a sequence shorter than triggerWindow → EMPTY (no full trigger window), no throw.
    //    • case-INSENSITIVE: the sequence is upper-cased before counting; only ASCII letters
    //      A–Z contribute to composition counts, but the entropy denominator is the full
    //      window length (non-letter chars are an absence of count, never a NaN).
    //
    //  Theory-correct invariants asserted (§2.4):
    //    • INV-01 — 0 ≤ Start ≤ End < n (in-bounds 0-based inclusive span).
    //    • INV-02 — emitted segments are non-overlapping and non-adjacent, ordered by Start.
    //    • INV-03 — a homopolymer window has H = 0 ≤ K1 (always triggers).
    //    • INV-04 — a window of W distinct residues has H = log₂(W) > K2 (default W=12 ⇒ not flagged).
    //    • INV-05 — every reported segment has length ≥ minLength.
    //
    //  Fuzz bar (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"): degenerate/boundary input
    //  must NEVER crash, hang, corrupt state, or produce a malformed span. The headline hazards
    //  for a windowed Shannon-entropy SEG scanner are: a NullReferenceException on null (must be
    //  the documented ArgumentNullException instead); an IndexOutOfRangeException when the
    //  sequence is shorter than the window (must yield nothing); a NaN from log₂(0) (a count-0
    //  residue must never enter the sum); and an out-of-bounds [Start..End] span. The BE targets
    //  cover the homopolymer (maximally low complexity ⇒ one segment), the high-complexity
    //  mosaic (⇒ no false positive) and the empty / short input corners.
    // ───────────────────────────────────────────────────────────────────────────

    #region Helpers — well-formed SEG segment contract

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted SEG segment must satisfy
    /// against a sequence of length <paramref name="n"/> (Low_Complexity_Region_Detection.md
    /// §2.4, §3.2): INV-01 — in-bounds, ordered 0-based inclusive span 0 ≤ Start ≤ End ≤ n−1;
    /// the segment has positive length and a non-null, non-empty <c>Type</c> label.
    /// </summary>
    private static void AssertWellFormedSegment((int Start, int End, string Type) seg, int n)
    {
        seg.Start.Should().BeInRange(0, n - 1, "INV-01: a segment Start is a valid 0-based residue index");
        seg.End.Should().BeInRange(seg.Start, n - 1, "INV-01: a segment End is in-bounds and not before its Start");
        seg.Type.Should().NotBeNullOrEmpty("every segment carries a composition Type label (§5.4)");
        seg.Type.Should().EndWith("-rich", "the Type label is an 'X-rich' / 'X/Y-rich' composition tag (§5.4)");
    }

    /// <summary>
    /// Asserts INV-02: a list of segments is ordered by Start and non-overlapping / non-adjacent
    /// (the merge step combines any span with <c>start ≤ lastEnd + 1</c>, so distinct emitted
    /// segments are separated by at least one residue) — Low_Complexity_Region_Detection.md §2.4.
    /// </summary>
    private static void AssertSegmentsDisjointOrdered(List<(int Start, int End, string Type)> segs)
    {
        for (int i = 1; i < segs.Count; i++)
            segs[i].Start.Should().BeGreaterThan(segs[i - 1].End + 1,
                "INV-02: emitted segments are non-overlapping AND non-adjacent (merged when start ≤ lastEnd + 1)");
    }

    #endregion

    #region BE — null / empty / shorter-than-window: documented throw vs empty, no crash

    /// <summary>
    /// Targets "empty" and the short-input boundary. A <c>null</c> sequence is the ONLY documented
    /// throw — <c>ArgumentNullException</c> (§3.3) — and because the method is a deferred iterator
    /// the throw must surface on enumeration. An empty string, and ANY sequence strictly shorter
    /// than the trigger window (even a perfect homopolymer that WOULD be low-complexity if long
    /// enough), must yield NO segments — no full trigger window exists — never an
    /// IndexOutOfRangeException running off the absent window (§3.3, §6.1 "sequence shorter than W
    /// → empty"). We probe null, "", every length 1..W−1 of a homopolymer, and the exact W−1
    /// boundary; W exactly is the smallest reportable length and is covered positively below.
    /// </summary>
    [Test]
    public void PredictLcr_NullThrows_EmptyOrShorterThanWindow_NoRegionsNoCrash()
    {
        // null → the single documented ArgumentNullException, surfaced on enumeration.
        var nullAct = () => DisorderPredictor.PredictLowComplexityRegions(null!).ToList();
        nullAct.Should().Throw<ArgumentNullException>("a null sequence is the only documented throw (§3.3)");

        // "" → no window → empty, no crash.
        var emptyAct = () => DisorderPredictor.PredictLowComplexityRegions("").ToList();
        emptyAct.Should().NotThrow("an empty sequence must not crash")
            .Subject.Should().BeEmpty("an empty sequence has no trigger window, so no segments");

        // Every length below the window — even a perfect homopolymer — has no complete window → empty.
        for (int len = 1; len < SegWindow; len++)
        {
            string homo = new string('A', len);
            var act = () => DisorderPredictor.PredictLowComplexityRegions(homo).ToList();
            act.Should().NotThrow($"a length-{len} sequence (< window {SegWindow}) must not crash")
                .Subject.Should().BeEmpty($"a length-{len} sequence is shorter than the window, so no complete window exists");
        }

        // Exactly W−1 of a homopolymer: still one short of a complete window → empty.
        DisorderPredictor.PredictLowComplexityRegions(new string('Q', SegWindow - 1)).Should().BeEmpty(
            "a sequence one residue shorter than the window yields no segments");
    }

    #endregion

    #region BE — homopolymer: the whole tract is one segment, X-rich, H = 0 ≤ K1

    /// <summary>
    /// Target "homopolymer" — the headline POSITIVE outcome (§6.1, §7.1, INV-03). A homopolymer
    /// "AAAA…" of length n ≥ W has H = 0 in EVERY window (a single residue type has p = 1,
    /// −1·log₂1 = 0), so every window triggers (0 ≤ K1) and the entire run merges into ONE segment
    /// spanning the whole tract [0, n−1]; its <c>Type</c> is "X-rich" (dominant fraction 1.0 > 0.5,
    /// §5.4). We assert this for every standard residue and several lengths, then pin the documented
    /// poly-Q worked example: <c>new string('Q', 26)</c> → a single (0, 25, "Q-rich") segment (§7.1).
    /// </summary>
    [Test]
    public void PredictLcr_Homopolymer_SingleSegmentWholeTractXRich()
    {
        foreach (char aa in StandardAminoAcids)
        {
            foreach (int n in new[] { SegWindow, SegWindow + 1, 20, 50 })
            {
                string homo = new string(aa, n);
                var segs = DisorderPredictor.PredictLowComplexityRegions(homo).ToList();

                segs.Should().ContainSingle($"a length-{n} homopolymer of '{aa}' is one merged low-complexity segment");
                var s = segs[0];
                AssertWellFormedSegment(s, n);
                s.Start.Should().Be(0, "INV-03: the homopolymer segment starts at residue 0");
                s.End.Should().Be(n - 1, "the homopolymer segment spans the whole tract [0, n−1]");
                s.Type.Should().Be($"{aa}-rich", "a single-residue tract is dominated > 0.5 by that residue (§5.4)");
            }
        }

        // Documented poly-Q worked example (§7.1): new string('Q', 26) → [ (0, 25, "Q-rich") ].
        var polyQ = DisorderPredictor.PredictLowComplexityRegions(new string('Q', 26)).ToList();
        polyQ.Should().ContainSingle("the §7.1 worked example: a 26-Q tract is one low-complexity segment");
        polyQ[0].Should().Be((0, 25, "Q-rich"), "every 12-window of the Q tract has entropy 0 ≤ 2.2 → (0, 25, \"Q-rich\")");
    }

    /// <summary>
    /// Poly-Q tract embedded in diverse flanks: a 20-Q homopolymer flanked by diverse 8-residue
    /// segments yields ONE segment that covers the Q tract; the surrounding diverse flanks are
    /// high-complexity and must NOT spawn a spurious second LCR (INV-02). The merged segment's
    /// dominant residue is Q (> 0.5 over the merged span), so its Type is "Q-rich" (§5.4).
    /// </summary>
    [Test]
    public void PredictLcr_PolyQTractInDiverseFlanks_ReportsSingleQRichSegment()
    {
        const string flank = "MKLPRDST";                       // 8 diverse residues
        string seq = flank + new string('Q', 20) + flank;       // 36 residues
        int qStart = flank.Length;                              // 8
        int qEnd = qStart + 20 - 1;                             // 27

        var segs = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();

        segs.Should().ContainSingle("a single LCR spans the poly-Q tract amid diverse flanks");
        foreach (var s in segs) AssertWellFormedSegment(s, seq.Length);
        var region = segs[0];
        region.Start.Should().BeLessThanOrEqualTo(qStart, "the segment begins no later than the Q tract start");
        region.End.Should().BeGreaterThanOrEqualTo(qEnd, "the segment extends at least to the Q tract end");
        region.Type.Should().Be("Q-rich", "Q dominates the merged span (> 0.5) → 'Q-rich' (§5.4)");
    }

    #endregion

    #region BE — high-complexity: maximally diverse sequence yields NO segment (no false positive)

    /// <summary>
    /// Target "high-complexity" — the headline NEGATIVE outcome (§6.1 "window of W distinct
    /// residues → not flagged", INV-04). A maximally diverse mosaic where residue i = (i mod 20)-th
    /// amino acid puts 12 DISTINCT residues in every width-12 window, so H = log₂(12) ≈ 3.585 > K2
    /// = 2.5 everywhere: no window can trigger or extend, hence NO segment (no false positive). We
    /// assert emptiness across several lengths, then PROVE the negative is real, not a no-op
    /// scanner: pushing BOTH cutoffs above log₂(12) makes the whole diverse sequence one segment.
    /// </summary>
    [Test]
    public void PredictLcr_HighComplexityDiverseSequence_NoSegments()
    {
        foreach (int n in new[] { SegWindow, 20, 40, 100, 200 })
        {
            string diverse = DiverseProtein(n);
            var segs = DisorderPredictor.PredictLowComplexityRegions(diverse).ToList();
            segs.Should().BeEmpty(
                $"a maximally diverse length-{n} sequence has every window entropy ≈ log₂(12) > K2, so no LCR");
        }

        // Premise: even raising K2 above log₂(12) leaves nothing TRIGGERED (every H ≈ 3.585 > K1 = 2.2).
        string diverse40 = DiverseProtein(40);
        DisorderPredictor.PredictLowComplexityRegions(
                diverse40, triggerThreshold: SegK1, extensionThreshold: 4.0)
            .Should().BeEmpty("even extending over every diverse window, none triggers (H ≈ 3.585 > K1 = 2.2)");

        // With BOTH cutoffs above log₂(12) the whole diverse sequence DOES become one segment —
        // proving the prior emptiness was due to the cutoffs, not a dead scanner.
        var forced = DisorderPredictor.PredictLowComplexityRegions(
            diverse40, triggerThreshold: 4.0, extensionThreshold: 4.0).ToList();
        forced.Should().ContainSingle("with both cutoffs above log₂(12) the diverse sequence is one big 'low-complexity' segment");
        AssertWellFormedSegment(forced[0], diverse40.Length);
        forced[0].Start.Should().Be(0);
        forced[0].End.Should().Be(diverse40.Length - 1, "the forced segment spans the whole sequence");
    }

    #endregion

    #region Positive sanity — documented entropy walk-through and a run found at the correct span

    /// <summary>
    /// Positive sanity reproducing the documented numerical walk-through (§7.1): the window
    /// <c>AAABBBCCCDDD</c> (L=12, four residue types × 3) has pᵢ = 0.25 each and
    /// H = −4·(0.25·log₂0.25) = 2.0 bits, which is ≤ the default K1 = 2.2 ⇒ it TRIGGERS one
    /// segment; with a strict K1 = 0.5, 2.0 > 0.5 ⇒ NO trigger. Its dominant residue fraction is
    /// 0.25 (≤ 0.5), so the Type is the two-residue "A/B-rich" label (§5.4). The 2.0-bit value is
    /// derived independently from the entropy formula, NOT echoed off the code.
    /// </summary>
    [Test]
    public void PredictLcr_FourTypesWindow_TriggersAtDefaultK1NotAtStrictK1()
    {
        const string fourTypes = "AAABBBCCCDDD"; // H = −4·(0.25·log₂0.25) = 2.0 bits, length 12

        // Default K1 = 2.2: 2.0 ≤ 2.2 ⇒ exactly one full-window segment [0, 11], "A/B-rich".
        var triggered = DisorderPredictor.PredictLowComplexityRegions(fourTypes).ToList();
        triggered.Should().ContainSingle("H = 2.0 ≤ K1 = 2.2 ⇒ the four-types window triggers one segment");
        AssertWellFormedSegment(triggered[0], fourTypes.Length);
        triggered[0].Start.Should().Be(0, "the single trigger window spans the whole 12-residue sequence");
        triggered[0].End.Should().Be(fourTypes.Length - 1);
        triggered[0].Type.Should().Be("A/B-rich",
            "no residue exceeds 0.5 (each is 0.25) ⇒ the top-two label A/B-rich (§5.4)");

        // Strict K1 = 0.5: 2.0 > 0.5 ⇒ nothing triggers ⇒ empty.
        DisorderPredictor.PredictLowComplexityRegions(fourTypes, triggerThreshold: 0.5)
            .Should().BeEmpty("with K1 = 0.5, window entropy 2.0 bits > 0.5 ⇒ nothing triggers");
    }

    /// <summary>
    /// Positive sanity: the scanner must FIND a genuinely low-complexity run at the correct span,
    /// not be a no-op. A diverse 30-residue prefix, a 16-A homopolymer, then a diverse 30-residue
    /// suffix yields exactly ONE segment whose span covers the homopolymer run and whose Type is
    /// "A-rich" (A dominates the merged span). The diverse flanks must not spawn extra segments
    /// (INV-02). Coordinates must be in-bounds and ordered (INV-01).
    /// </summary>
    [Test]
    public void PredictLcr_LowComplexityRunInDiverseFlanks_FoundAtCorrectSpan()
    {
        string diversePrefix = DiverseProtein(30);
        string diverseSuffix = DiverseProtein(30);
        const int runLen = 16;
        int runStart = diversePrefix.Length;                 // 30
        int runEnd = runStart + runLen - 1;                  // 45
        string seq = diversePrefix + new string('A', runLen) + diverseSuffix;

        var segs = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();
        segs.Should().ContainSingle("exactly one homopolymer LCR sits amid diverse flanks");
        AssertWellFormedSegment(segs[0], seq.Length);
        segs[0].Start.Should().BeLessThanOrEqualTo(runStart, "the segment begins no later than the homopolymer run");
        segs[0].End.Should().BeGreaterThanOrEqualTo(runEnd, "the segment extends at least to the homopolymer run end");
        segs[0].Type.Should().Be("A-rich", "A dominates the merged low-complexity span (> 0.5) → 'A-rich'");
    }

    #endregion

    #region BE — minLength filter and the no-NaN / determinism sweep over arbitrary input

    /// <summary>
    /// Target the <c>minLength</c> post-filter boundary (INV-05, §6.1 "minLength exceeds segment
    /// length → segment dropped"): a 12-A homopolymer is one length-12 segment; minLength = 12
    /// keeps it, minLength = 13 drops it, minLength = 1 (default) keeps it. A non-positive /
    /// extreme minLength must not crash (BE: 0, −1, int.MaxValue).
    /// </summary>
    [Test]
    public void PredictLcr_MinLengthFilter_DropsSegmentsBelowThreshold()
    {
        string homo12 = new string('A', SegWindow); // one length-12 segment [0, 11]

        DisorderPredictor.PredictLowComplexityRegions(homo12, minLength: SegWindow)
            .Should().ContainSingle("a length-12 segment survives minLength = 12 (length ≥ minLength)");
        DisorderPredictor.PredictLowComplexityRegions(homo12, minLength: SegWindow + 1)
            .Should().BeEmpty("a length-12 segment is dropped by minLength = 13 (INV-05)");
        DisorderPredictor.PredictLowComplexityRegions(homo12, minLength: 1)
            .Should().ContainSingle("the default minLength = 1 keeps every segment");

        // Extreme / non-positive minLength must not crash.
        foreach (int minLen in new[] { 0, -1, -100, int.MaxValue })
        {
            var act = () => DisorderPredictor.PredictLowComplexityRegions(homo12, minLength: minLen).ToList();
            act.Should().NotThrow($"a degenerate minLength ({minLen}) must not crash");
        }
        DisorderPredictor.PredictLowComplexityRegions(homo12, minLength: int.MaxValue)
            .Should().BeEmpty("an int.MaxValue minLength drops every finite-length segment");
    }

    /// <summary>
    /// Headline no-NaN / no-crash / determinism sweep over arbitrary and adversarial composition.
    /// The per-window entropy must NEVER be NaN (a count-0 residue must never enter the −p·log₂ p
    /// sum) and the scan must terminate, on random proteins, biased mixtures and out-of-alphabet
    /// junk (digits, punctuation, X, whitespace, period-2 repeats). Every emitted segment must be
    /// well-formed (in-bounds, ordered, disjoint), and re-running must give the IDENTICAL ordered
    /// result (determinism). Lowercase input must yield the SAME segments as its uppercase form
    /// (case-insensitive, §3.3). [CancelAfter] guards against a regression turning the O(n·W) scan
    /// into a hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictLcr_ArbitraryInput_AlwaysWellFormedDisjointDeterministic(CancellationToken token)
    {
        var inputs = new List<string>();
        foreach (int seed in new[] { 7, 31, 137, 2026 })
            foreach (int len in new[] { 12, 25, 60, 150 })
                inputs.Add(RandomProtein(len, seed));

        inputs.Add(new string('A', 30) + RandomProtein(30, 99));        // homopolymer + random
        inputs.Add("MK1RGD2SP" + new string('S', 20) + "!@#$%^&*()");    // junk + poly-S
        inputs.Add("XXXXXXXXXXXXXXXX");                                   // out-of-alphabet placeholder run
        inputs.Add("ABABABABABABABABABAB");                              // period-2 repeat (entropy 1.0/window)
        inputs.Add("   \t  \n  whitespace and text mixed in here  ");    // whitespace + literals

        foreach (string seq in inputs)
        {
            var act = () => DisorderPredictor.PredictLowComplexityRegions(seq).ToList();
            var segs = act.Should().NotThrow($"arbitrary input must not crash: '{seq[..Math.Min(seq.Length, 20)]}'").Subject;
            token.ThrowIfCancellationRequested();

            foreach (var s in segs) AssertWellFormedSegment(s, seq.Length);
            AssertSegmentsDisjointOrdered(segs);

            // Determinism: the scan is reproducible for a fixed input.
            var again = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();
            again.Should().Equal(segs, "the SEG scan is deterministic for a fixed input");
        }

        // Case-insensitivity: lowercase must give identical segments to uppercase.
        string mixed = "mklprdst" + new string('q', 20) + "MKLPRDST";
        var lower = DisorderPredictor.PredictLowComplexityRegions(mixed).ToList();
        var upper = DisorderPredictor.PredictLowComplexityRegions(mixed.ToUpperInvariant()).ToList();
        lower.Should().Equal(upper, "SEG is case-insensitive: lowercase input yields the same segments as uppercase");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  DISORDER-MORF-001 — MoRF (dip-in-disorder) detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region DISORDER-MORF-001 — MoRF (Molecular Recognition Feature) detection

    // ───────────────────────────────────────────────────────────────────────────
    //  Unit: DISORDER-MORF-001 — MoRF (Molecular Recognition Feature) prediction.
    //  Checklist: docs/checklists/03_FUZZING.md, row 205 (ProteinPred, strategy BE;
    //    targets: "fully ordered, fully disordered, short").
    //  Method under test (src/.../Seqeron.Genomics.Analysis/DisorderPredictor.cs):
    //    IEnumerable<(int Start,int End,double Score)> PredictMoRFs(
    //        string sequence, int minLength = 10, int maxLength = 70)
    //    — MoRF_Prediction.md §3.1, §3.2, §5.1.
    //
    //  The MoRF contract under test (MoRF_Prediction.md §1, §2.2, §4.1):
    //    A MoRF is a short segment of relative ORDER ("dip") embedded WITHIN a longer
    //    intrinsically-disordered region: it folds upon partner binding. Let d(i) ∈ [0,1] be
    //    the per-residue disorder score from PredictDisorder (normalised TOP-IDP). A residue
    //    is ORDERED when d(i) < 0.5 and DISORDERED when d(i) ≥ 0.5 (the 0.5 MoRF threshold,
    //    Cheng/Oldfield PMC2570644 — NOTE this is the literature MoRF cutoff, distinct from the
    //    0.542 TOP-IDP IDR-calling cutoff PredictDisorder uses internally; §5.2). A reported MoRF
    //    is a MAXIMAL interval [s, e] with: (1) d(i) < 0.5 for all i∈[s,e] (the ordered dip);
    //    (2) 10 ≤ (e − s + 1) ≤ 70 (Mohan 2006 length band); (3) d(s−1) ≥ 0.5 AND d(e+1) ≥ 0.5
    //    (flanked by disorder on BOTH immediate sides — embedded, not terminal). Its Score is the
    //    dip depth below the threshold, normalised: Score = clamp₀¹((0.5 − mean_{i∈[s,e]} d(i))/0.5);
    //    because every d(i) < 0.5 in the dip, mean ∈ [0, 0.5) ⇒ Score ∈ (0, 1] and a deeper
    //    (more ordered) dip scores higher (§2.2).
    //
    //  Documented input handling (§3.3, §6.1):
    //    • null / "" → empty result (no throw): a deferred-iterator-free early return. (No
    //      ArgumentNullException — confirmed by C1_NullInput_ReturnsEmpty in the unit suite.)
    //    • case-INSENSITIVE: the sequence is upper-cased internally before scoring.
    //    • residues outside the 20 standard amino acids contribute 0 disorder propensity via
    //      PredictDisorder (never a KeyNotFound / NaN).
    //
    //  Theory-correct invariants asserted (§2.4):
    //    • INV-01 — 0 ≤ Start ≤ End < len(sequence) (in-bounds 0-based inclusive span).
    //    • INV-02 — minLength ≤ (End − Start + 1) ≤ maxLength (length band filter).
    //    • INV-03 — every MoRF is flanked by d ≥ 0.5 on BOTH immediate sides (embedded in disorder),
    //               and d(i) < 0.5 strictly INSIDE the dip (the run is a maximal ordered stretch).
    //    • INV-04 — reported MoRFs are non-overlapping and ordered by Start.
    //    • INV-05 — 0 ≤ Score ≤ 1 (finite, never NaN/±∞), monotone in dip depth.
    //
    //  Fuzz bar (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"): degenerate/boundary input must
    //  NEVER crash, hang, corrupt state, or emit a malformed/non-finite MoRF. The headline hazards
    //  for a dip-in-disorder scanner are: an IndexOutOfRangeException in the flank check d(s−1)/d(e+1)
    //  when a dip touches a terminus (must be SUPPRESSED, not thrown — a terminal dip is simply not
    //  flanked, so not reported); a NaN Score from a length-0 mean; a KeyNotFound/NaN from a
    //  non-standard residue (must contribute 0 propensity); and an out-of-bounds [Start..End] span.
    //  The BE targets cover the FULLY-ORDERED corner (a dip with no surrounding disorder ⇒ no flank
    //  ⇒ no MoRF), the FULLY-DISORDERED corner (no ordered dip exists ⇒ no MoRF) and the SHORT corner
    //  (a sequence — or a dip — shorter than the 10-residue minimum band ⇒ no MoRF).
    // ───────────────────────────────────────────────────────────────────────────

    #region Helpers — well-formed MoRF contract

    /// <summary>
    /// The MoRF order/disorder threshold (Cheng/Oldfield PMC2570644; MoRF_Prediction.md §4.2) —
    /// derived from the PRIMARY source, NOT echoed off the implementation's constant. A residue is
    /// part of an ordered dip iff its disorder score is strictly below this value.
    /// </summary>
    private const double MoRFThreshold = 0.5;

    /// <summary>Mohan et al. (2006) length band lower bound — MoRF_Prediction.md §4.2.</summary>
    private const int MoRFMinLength = 10;

    /// <summary>Mohan et al. (2006) length band upper bound — MoRF_Prediction.md §4.2.</summary>
    private const int MoRFMaxLength = 70;

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted MoRF must satisfy against the
    /// supplied per-residue disorder profile (MoRF_Prediction.md §2.4): INV-01 — in-bounds 0-based
    /// inclusive span 0 ≤ Start ≤ End ≤ n−1; INV-02 — the length lies in [minLength, maxLength];
    /// INV-03 — strictly ordered (d &lt; 0.5) inside the dip, flanked by disorder (d ≥ 0.5) on BOTH
    /// immediate sides; INV-05 — Score is FINITE and in [0, 1] and equals the documented normalised
    /// dip depth (0.5 − mean d)/0.5 recomputed independently from the profile.
    /// </summary>
    private static void AssertWellFormedMoRF(
        (int Start, int End, double Score) morf,
        DisorderPredictionResult profile,
        int minLength, int maxLength)
    {
        var d = profile.ResiduePredictions;
        int n = d.Count;

        // INV-01 — in-bounds inclusive span.
        morf.Start.Should().BeInRange(0, n - 1, "INV-01: a MoRF Start is a valid 0-based residue index");
        morf.End.Should().BeInRange(morf.Start, n - 1, "INV-01: a MoRF End is in-bounds and not before its Start");

        // INV-02 — length lies in the requested band.
        int length = morf.End - morf.Start + 1;
        length.Should().BeInRange(minLength, maxLength,
            "INV-02: a MoRF length lies within [minLength, maxLength]");

        // INV-03 — strictly ordered inside the dip (maximal ordered run).
        double sum = 0;
        for (int i = morf.Start; i <= morf.End; i++)
        {
            d[i].DisorderScore.Should().BeLessThan(MoRFThreshold,
                $"INV-03: every residue inside a MoRF dip is ordered (d < 0.5) (position {i})");
            sum += d[i].DisorderScore;
        }

        // INV-03 — flanked by disorder on BOTH immediate sides (embedded within disorder).
        morf.Start.Should().BeGreaterThan(0, "INV-03: a MoRF cannot touch the N-terminus (no left flank)");
        morf.End.Should().BeLessThan(n - 1, "INV-03: a MoRF cannot touch the C-terminus (no right flank)");
        d[morf.Start - 1].DisorderScore.Should().BeGreaterThanOrEqualTo(MoRFThreshold,
            "INV-03: the residue immediately before a MoRF is disordered (d ≥ 0.5)");
        d[morf.End + 1].DisorderScore.Should().BeGreaterThanOrEqualTo(MoRFThreshold,
            "INV-03: the residue immediately after a MoRF is disordered (d ≥ 0.5)");

        // INV-05 — finite, in [0, 1], and the documented normalised dip depth.
        double.IsNaN(morf.Score).Should().BeFalse("INV-05: a MoRF Score must never be NaN");
        double.IsInfinity(morf.Score).Should().BeFalse("INV-05: a MoRF Score must never be infinite");
        morf.Score.Should().BeInRange(0.0, 1.0, "INV-05: a MoRF Score is normalised to [0, 1]");
        double expected = Math.Max(0.0, Math.Min(1.0, (MoRFThreshold - sum / length) / MoRFThreshold));
        morf.Score.Should().BeApproximately(expected, 1e-9,
            "INV-05: Score = clamp₀¹((0.5 − mean d) / 0.5), recomputed independently from the profile (§2.2)");
    }

    /// <summary>
    /// Asserts INV-04 over a MoRF list: non-overlapping and ordered by increasing Start
    /// (maximal disjoint runs scanned left→right) — MoRF_Prediction.md §2.4.
    /// </summary>
    private static void AssertMoRFsDisjointOrdered(List<(int Start, int End, double Score)> morfs)
    {
        for (int i = 1; i < morfs.Count; i++)
            morfs[i].Start.Should().BeGreaterThan(morfs[i - 1].End,
                "INV-04: reported MoRFs are non-overlapping and ordered by Start");
    }

    #endregion

    #region BE — null / empty / short: empty result, no crash (no flank IndexOutOfRange)

    /// <summary>
    /// Targets "short" and the empty corner (§3.3, §6.1 "null / empty / very short → empty result").
    /// A <c>null</c> and an empty string must return an EMPTY result with NO throw (the early return
    /// — NOT a deferred iterator, so no ArgumentNullException). ANY sequence shorter than the
    /// 10-residue minimum band cannot contain a 10-residue embedded dip, so it must yield no MoRF —
    /// and the dip scan / flank check must never run off the short string (no IndexOutOfRange). We
    /// probe null, "", and every length 0..minLength−1 of an ordered homopolymer (poly-L, the
    /// canonical ordered residue) which WOULD dip if it were long enough but is too short to qualify.
    /// </summary>
    [Test]
    public void PredictMoRFs_NullEmptyOrShorterThanMinLength_EmptyNoCrash()
    {
        // null and "" → empty, no throw (early return, not a deferred iterator).
        foreach (string? seq in new[] { null, "" })
        {
            var act = () => DisorderPredictor.PredictMoRFs(seq!).ToList();
            act.Should().NotThrow($"null/empty input ('{seq ?? "null"}') must not crash")
                .Subject.Should().BeEmpty($"null/empty input yields no MoRFs (§3.3)");
        }

        // Every length below the 10-residue minimum — even a perfect ordered poly-L — has no
        // qualifying dip and must never run off the end in the dip scan / flank check.
        for (int len = 1; len < MoRFMinLength; len++)
        {
            string ordered = new string('L', len);
            var act = () => DisorderPredictor.PredictMoRFs(ordered).ToList();
            act.Should().NotThrow($"a length-{len} sequence (< minLength {MoRFMinLength}) must not crash")
                .Subject.Should().BeEmpty($"a length-{len} sequence cannot contain a {MoRFMinLength}-residue MoRF (§6.1)");
        }
    }

    #endregion

    #region BE — fully ordered: a dip with no surrounding disorder → no MoRF

    /// <summary>
    /// Target "fully ordered" (§6.1 "fully ordered sequence → no MoRFs", INV-03). A protein whose
    /// per-residue disorder is ENTIRELY below 0.5 (an ordered homopolymer such as poly-L, poly-W,
    /// poly-I, or a hydrophobic globular stretch) is one long ordered run that reaches BOTH termini —
    /// it has NO surrounding disorder to embed a dip, so neither flank is disordered and NO MoRF is
    /// reported. The flank check (d(s−1)/d(e+1)) must SUPPRESS the terminal run, never throw an
    /// IndexOutOfRange reaching past position 0 or n−1. We assert emptiness across several ordered
    /// compositions and lengths spanning the 10–70 band, including lengths well inside it.
    /// </summary>
    [Test]
    public void PredictMoRFs_FullyOrdered_NoMoRFs()
    {
        foreach (char aa in new[] { 'L', 'W', 'I', 'F', 'V', 'C' }) // order-promoting residues (d < 0.5)
        {
            foreach (int n in new[] { MoRFMinLength, 20, 40, MoRFMaxLength, 90 })
            {
                string seq = new string(aa, n);

                // Premise: this homopolymer really is fully ordered (every d < 0.5).
                var profile = PredictDisorder(seq);
                profile.ResiduePredictions.Should().OnlyContain(p => p.DisorderScore < MoRFThreshold,
                    $"poly-{aa} (length {n}) is fully ordered — every residue scores below 0.5");

                var act = () => DisorderPredictor.PredictMoRFs(seq).ToList();
                act.Should().NotThrow($"a fully-ordered poly-{aa} (length {n}) must not crash on the flank check")
                    .Subject.Should().BeEmpty(
                        $"a fully-ordered protein has no surrounding disorder to embed a dip → no MoRF (§6.1)");
            }
        }

        // A hydrophobic globular stretch (W/F/I/L/V) is likewise fully ordered → no MoRF.
        DisorderPredictor.PredictMoRFs("WFILVWFILVWFILVWFILVWFILVWFILVWFILV").Should().BeEmpty(
            "a hydrophobic globular stretch is fully ordered → no MoRF (it reaches both termini)");
    }

    #endregion

    #region BE — fully disordered: no ordered dip exists → no MoRF

    /// <summary>
    /// Target "fully disordered" (§6.1 "fully disordered sequence → no MoRFs", INV-03). A protein
    /// whose per-residue disorder is ENTIRELY at or above 0.5 (a disorder-promoting homopolymer such
    /// as poly-P, poly-E, poly-K) contains NO ordered dip (no residue with d &lt; 0.5), so the dip
    /// scan opens no run and NO MoRF is reported. This is the complementary BE corner to "fully
    /// ordered": there, every flank fails; here, no dip ever forms. Nothing may crash.
    /// </summary>
    [Test]
    public void PredictMoRFs_FullyDisordered_NoMoRFs()
    {
        foreach (char aa in new[] { 'P', 'E', 'K', 'S' }) // disorder-promoting residues (d ≥ 0.5)
        {
            foreach (int n in new[] { MoRFMinLength, 20, 40, 90 })
            {
                string seq = new string(aa, n);

                // Premise: this homopolymer really is fully disordered (every d ≥ 0.5).
                var profile = PredictDisorder(seq);
                profile.ResiduePredictions.Should().OnlyContain(p => p.DisorderScore >= MoRFThreshold,
                    $"poly-{aa} (length {n}) is fully disordered — every residue scores at/above 0.5");

                var act = () => DisorderPredictor.PredictMoRFs(seq).ToList();
                act.Should().NotThrow($"a fully-disordered poly-{aa} (length {n}) must not crash")
                    .Subject.Should().BeEmpty(
                        $"a fully-disordered protein contains no ordered dip → no MoRF (§6.1)");
            }
        }
    }

    #endregion

    #region BE — dip outside the 10–70 band: too-short and too-long dips are not reported

    /// <summary>
    /// Targets the LENGTH-band boundaries (§6.1 "dip &lt; 10 or &gt; 70 residues → not reported",
    /// INV-02). A genuine ordered dip EMBEDDED in disorder is reported only when its length is in
    /// [minLength, maxLength]. We embed an ordered poly-L core inside disordered poly-P flanks and
    /// sweep the core length across the default lower boundary (minLength = 10): a 9-residue core is
    /// below the band → no MoRF; a 12-residue core is inside → exactly one MoRF whose length is in
    /// band. We also pin the boundary from the parameter side: with the SAME in-disorder dip, a
    /// minLength set ABOVE the dip length drops it, and a maxLength set BELOW it drops it — proving
    /// the band filter is real and inclusive, not a no-op. The disorder-window smoothing in
    /// PredictDisorder blurs the exact dip edges, so the reported length need not equal the core
    /// length exactly; we assert only that it lies within the requested band (INV-02) and that the
    /// boundary cases flip reported↔empty as the band moves.
    /// </summary>
    [Test]
    public void PredictMoRFs_DipOutsideLengthBand_NotReported()
    {
        const int flank = 25; // disordered poly-P flanks, long enough that the dip stays embedded

        // Below the band: a very short ordered core cannot produce a >= 10-residue reported dip.
        string shortCore = new string('P', flank) + new string('L', 5) + new string('P', flank);
        DisorderPredictor.PredictMoRFs(shortCore).Should().BeEmpty(
            "a 5-residue ordered core is below the 10-residue minimum band → no MoRF (§6.1)");

        // Inside the band: a 30-residue ordered core yields exactly one in-band MoRF.
        string inBand = new string('P', flank) + new string('L', 30) + new string('P', flank);
        var profile = PredictDisorder(inBand);
        var morfs = DisorderPredictor.PredictMoRFs(inBand).ToList();
        morfs.Should().ContainSingle("a 30-residue ordered core embedded in disorder is one in-band MoRF");
        AssertWellFormedMoRF(morfs[0], profile, MoRFMinLength, MoRFMaxLength);
        AssertMoRFsDisjointOrdered(morfs);

        // Parameter side — the SAME embedded dip, band moved around it:
        int dipLen = morfs[0].End - morfs[0].Start + 1;
        // minLength above the dip length drops it.
        DisorderPredictor.PredictMoRFs(inBand, minLength: dipLen + 1, maxLength: MoRFMaxLength).Should().BeEmpty(
            $"a minLength ({dipLen + 1}) above the dip length ({dipLen}) drops the MoRF (INV-02)");
        // maxLength below the dip length drops it.
        DisorderPredictor.PredictMoRFs(inBand, minLength: 1, maxLength: dipLen - 1).Should().BeEmpty(
            $"a maxLength ({dipLen - 1}) below the dip length ({dipLen}) drops the MoRF (INV-02)");
        // A band that contains the dip length keeps it.
        DisorderPredictor.PredictMoRFs(inBand, minLength: dipLen, maxLength: dipLen).Should().ContainSingle(
            $"an exact-fit band [{dipLen}, {dipLen}] keeps the MoRF (inclusive bounds, INV-02)");
    }

    #endregion

    #region Positive sanity — the documented §7.1 worked example, hand-checkable

    /// <summary>
    /// Positive sanity reproducing the DOCUMENTED worked example (MoRF_Prediction.md §7.1): the
    /// sequence <c>P²⁵ L³⁰ P²⁵</c> — 25 disordered prolines, 30 ordered leucines, 25 disordered
    /// prolines — has a smoothed disorder profile that dips below 0.5 over residues [29, 50]
    /// (length 22, inside the 10–70 band) flanked by disorder on both sides; the mean disorder over
    /// the dip ≈ 0.362033, so the reported Score = (0.5 − 0.362033)/0.5 ≈ 0.275934. These coordinates
    /// and the score are stated in the doc INDEPENDENTLY of the code, and re-derived here from the
    /// (0.5 − mean d)/0.5 formula against the actual profile. We also pin score MONOTONICITY (INV-05):
    /// replacing the L core with the more order-promoting I deepens the dip, raising the score
    /// (§7.1: mean ≈ 0.300196, score ≈ 0.399608 &gt; the L-dip score).
    /// </summary>
    [Test]
    public void PredictMoRFs_WorkedExample_OneEmbeddedDipWithDocumentedScore()
    {
        string seqL = new string('P', 25) + new string('L', 30) + new string('P', 25); // §7.1
        var profileL = PredictDisorder(seqL);
        var morfsL = DisorderPredictor.PredictMoRFs(seqL).ToList();

        morfsL.Should().ContainSingle("the §7.1 example: one ordered L-dip embedded in disordered P flanks");
        AssertWellFormedMoRF(morfsL[0], profileL, MoRFMinLength, MoRFMaxLength);
        AssertMoRFsDisjointOrdered(morfsL);

        // Documented coordinates and score (§7.1) — derived independently from the doc, not the code.
        morfsL[0].Start.Should().Be(29, "the smoothed dip begins at residue 29 (§7.1)");
        morfsL[0].End.Should().Be(50, "the smoothed dip ends at residue 50 (§7.1)");
        morfsL[0].Score.Should().BeApproximately(0.275934, 1e-6,
            "Score = (0.5 − mean disorder 0.362033) / 0.5 ≈ 0.275934 (§7.1)");

        // INV-05 monotonicity: a deeper (more ordered) I-dip scores HIGHER than the L-dip.
        string seqI = new string('P', 25) + new string('I', 30) + new string('P', 25);
        var morfsI = DisorderPredictor.PredictMoRFs(seqI).ToList();
        morfsI.Should().ContainSingle("the order-promoting I core also yields one embedded dip");
        morfsI[0].Score.Should().BeApproximately(0.399608, 1e-6,
            "the deeper I-dip score = (0.5 − mean disorder 0.300196) / 0.5 ≈ 0.399608 (§7.1)");
        morfsI[0].Score.Should().BeGreaterThan(morfsL[0].Score,
            "INV-05: a deeper (more ordered) dip scores higher (monotone in dip depth)");
    }

    /// <summary>
    /// Positive sanity: a terminal dip is NOT reported even though it satisfies the length band —
    /// the embedding (flank) requirement is genuinely enforced (§6.1 "dip at sequence terminus →
    /// not reported", INV-03). An ordered poly-L core that runs to the C-terminus (no trailing
    /// disordered flank) and a leading-ordered core (no leading disordered flank) both yield NO MoRF,
    /// while the SAME core embedded between two disordered flanks yields exactly one. This proves the
    /// flank check is the discriminator, not the dip itself, and that the d(e+1)/d(s−1) lookups at the
    /// edges never throw.
    /// </summary>
    [Test]
    public void PredictMoRFs_TerminalDip_NotReportedUnlessFlankedBothSides()
    {
        const int flank = 25;
        string core = new string('L', 30); // an in-band ordered dip on its own

        // Leading-ordered: the dip touches the N-terminus → no left flank → not reported.
        DisorderPredictor.PredictMoRFs(core + new string('P', flank)).Should().BeEmpty(
            "an ordered dip at the N-terminus has no disordered left flank → not reported (§6.1)");

        // Trailing-ordered: the dip touches the C-terminus → no right flank → not reported.
        DisorderPredictor.PredictMoRFs(new string('P', flank) + core).Should().BeEmpty(
            "an ordered dip at the C-terminus has no disordered right flank → not reported (§6.1)");

        // Embedded between two disordered flanks → exactly one MoRF (the flank check is the discriminator).
        string embedded = new string('P', flank) + core + new string('P', flank);
        DisorderPredictor.PredictMoRFs(embedded).Should().ContainSingle(
            "the SAME core embedded between disordered flanks IS reported — the flank check is the discriminator");
    }

    #endregion

    #region BE — randomized boundary sweep: never crash / hang / NaN, contract holds, deterministic

    /// <summary>
    /// Headline no-crash / no-hang / no-NaN sweep over arbitrary and adversarial input and a sweep of
    /// the length-band parameters. Across fixed seeds, lengths and (minLength, maxLength) corners
    /// (including degenerate ones: 0, −1, equal bounds, inverted bounds, and the int.MaxValue extreme)
    /// the scan must terminate, never throw, and every emitted MoRF must satisfy the full contract
    /// (in-bounds, in-band, strictly-ordered-and-flanked, finite in-range Score — AssertWellFormedMoRF)
    /// and be non-overlapping / ordered (INV-04). We also feed out-of-alphabet junk (digits,
    /// punctuation, X, whitespace) — which contributes 0 propensity, never a KeyNotFound/NaN — and
    /// hand-built ordered/disordered mosaics. Re-running must give the IDENTICAL ordered result
    /// (determinism), and lowercase input must yield the SAME MoRFs as its uppercase form (§3.3).
    /// [CancelAfter] guards against a regression turning the O(n·w) scan into a hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictMoRFs_ArbitraryInput_AlwaysWellFormedDeterministic(CancellationToken token)
    {
        var inputs = new List<string>();
        foreach (int seed in new[] { 7, 31, 137, 2026 })
            foreach (int len in new[] { 10, 25, 80, 200 })
                inputs.Add(RandomProtein(len, seed));

        // Hand-built ordered/disordered mosaics that genuinely exercise the dip scan.
        inputs.Add(new string('P', 25) + new string('L', 30) + new string('P', 25));            // one embedded dip
        inputs.Add(new string('P', 25) + new string('L', 30) + new string('P', 30)
                   + new string('L', 30) + new string('P', 25));                                 // two embedded dips
        inputs.Add(new string('L', 40));                                                          // fully ordered
        inputs.Add(new string('P', 40));                                                          // fully disordered
        // Out-of-alphabet junk — 0 propensity, never KeyNotFound / NaN.
        inputs.Add("MK1RGD2SP" + new string('S', 20) + "!@#$%^&*()");
        inputs.Add("XXXXXXXXXXXXXXXXXXXX");
        inputs.Add("   \t  \n  whitespace and text mixed in here  ");

        foreach (string seq in inputs)
        {
            foreach ((int minLen, int maxLen) in new[]
                     {
                         (MoRFMinLength, MoRFMaxLength), // defaults
                         (1, 5), (1, int.MaxValue),      // wide / tiny bands
                         (10, 10),                       // single-length band
                         (0, 0), (-1, -1),               // degenerate non-positive bounds
                         (50, 10),                       // inverted band (max < min) → nothing in band
                     })
            {
                var act = () => DisorderPredictor.PredictMoRFs(seq, minLength: minLen, maxLength: maxLen).ToList();
                var morfs = act.Should().NotThrow(
                    $"arbitrary input must not crash: '{seq[..Math.Min(seq.Length, 16)]}' (band [{minLen}, {maxLen}])").Subject;
                token.ThrowIfCancellationRequested();

                var profile = PredictDisorder(seq);
                foreach (var m in morfs) AssertWellFormedMoRF(m, profile, minLen, maxLen);
                AssertMoRFsDisjointOrdered(morfs);

                // Inverted band can never contain any length → no MoRF.
                if (maxLen < minLen)
                    morfs.Should().BeEmpty("an inverted band (max < min) contains no length → no MoRF");

                // Determinism: the scan is reproducible for a fixed input and band.
                var again = DisorderPredictor.PredictMoRFs(seq, minLength: minLen, maxLength: maxLen).ToList();
                again.Should().Equal(morfs, "PredictMoRFs is deterministic for a fixed input and band");
            }
        }

        // Case-insensitivity (§3.3): lowercase yields identical MoRFs to its uppercase form.
        string mixed = new string('p', 25) + new string('l', 30) + new string('p', 25);
        var lower = DisorderPredictor.PredictMoRFs(mixed).ToList();
        var upper = DisorderPredictor.PredictMoRFs(mixed.ToUpperInvariant()).ToList();
        lower.Should().Equal(upper, "PredictMoRFs is case-insensitive: lowercase yields the same MoRFs as uppercase");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  DISORDER-PROPENSITY-001 — TOP-IDP propensity & Dunker class : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region DISORDER-PROPENSITY-001 — TOP-IDP propensity lookup & Dunker classification

    // ───────────────────────────────────────────────────────────────────────────
    //  Unit: DISORDER-PROPENSITY-001 — per-amino-acid intrinsic-disorder propensity.
    //  Checklist: docs/checklists/03_FUZZING.md, row 206 (ProteinPred, strategies BE, MC;
    //    targets: "empty, non-amino-acid, single residue").
    //  Methods under test (src/.../Seqeron.Genomics.Analysis/DisorderPredictor.cs):
    //    double GetDisorderPropensity(char aminoAcid)
    //    bool   IsDisorderPromoting(char aminoAcid)
    //    IReadOnlyList<char> DisorderPromotingAminoAcids / OrderPromotingAminoAcids / AmbiguousAminoAcids
    //    — Disorder_Propensity.md §3.1, §3.2, §5.1.
    //
    //  The contract under test (Disorder_Propensity.md §1, §2.2, §3):
    //    These are PURE O(1) TABLE LOOKUPS — "no statistics or windowing are involved" (§1).
    //    Two reference tables fully determine behaviour:
    //      • the TOP-IDP per-residue scale (Campen et al. 2008, Table 2) returned by
    //        GetDisorderPropensity — a continuous propensity in [−0.884 (W), 0.987 (P)] for each
    //        of the 20 standard residues, 0.0 for anything outside that set (§2.2, §3.2, §6.1);
    //      • the Dunker et al. (2001) three-way partition of the 20 residues into
    //        disorder-promoting {A,R,G,Q,S,P,E,K} (8), order-promoting {W,C,F,I,Y,V,L,N} (8),
    //        ambiguous {H,M,T,D} (4) — exposed by IsDisorderPromoting and the three sorted
    //        list properties (§2.2, §3.2).
    //
    //  PRIMARY-SOURCE table re-derivation (NOT echoed off the implementation array):
    //    The 20 TOP-IDP values below are transcribed from Campen et al. (2008) "TOP-IDP-Scale"
    //    Protein Pept Lett 15(9):956-963, Table 2, exactly as restated in Disorder_Propensity.md
    //    §2.2 — W=−0.884, F=−0.697, Y=−0.510, I=−0.486, M=−0.397, L=−0.326, V=−0.121, N=0.007,
    //    C=0.020, T=0.059, A=0.060, G=0.166, R=0.180, D=0.192, H=0.303, Q=0.318, S=0.341,
    //    K=0.586, E=0.736, P=0.987. The Dunker sets are transcribed from Dunker et al. (2001)
    //    J Mol Graph Model 19:26-59 as restated in §2.2. Reading "expected" off DisorderPredictor's
    //    own dictionary would be a forbidden code-echo; these constants are the independent oracle.
    //
    //  Documented input handling (§3.3, §6.1):
    //    • Any character NOT among the 20 standard residues — 'X', the ambiguity codes B/Z/J/O/U,
    //      digits, punctuation, whitespace, control / unicode chars, nucleotide letters — returns
    //      propensity 0.0 and IsDisorderPromoting false. NO exception is ever thrown (deviation #1,
    //      accepted): GetValueOrDefault / HashSet.Contains have no failing path (§3.3, §5.4, §6.1).
    //    • Case-INSENSITIVE: the input is upper-cased before lookup, so 'p' == 'P', 'w' == 'W'
    //      (INV-05, §6.1).
    //
    //  Theory-correct invariants asserted (§2.4):
    //    • INV-01 — GetDisorderPropensity returns the EXACT Table 2 value for each of the 20 residues.
    //    • INV-02 — over the 20 residues the value lies in [−0.884, 0.987]; min at W, max at P.
    //    • INV-03 — IsDisorderPromoting(c) ⇔ c ∈ DisorderPromotingAminoAcids.
    //    • INV-04 — the disorder (8) / order (8) / ambiguous (4) sets are pairwise DISJOINT and
    //               together COVER all 20 standard residues.
    //    • INV-05 — both lookups are case-insensitive.
    //
    //  Fuzz bar (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"): degenerate / out-of-alphabet
    //  input must NEVER crash, hang, corrupt state, or produce a non-finite value. The headline
    //  hazards for two constant-table lookups are: a KeyNotFoundException / NullReference when an
    //  out-of-table residue (MC: 'X', B/Z/J/O/U, digits, punctuation, '\0', unicode) is fed to
    //  GetDisorderPropensity or IsDisorderPromoting — both must return the documented neutral
    //  (0.0 / false), never throw; and a NaN / ±∞ ever surfacing from a table that holds only the
    //  20 finite published constants. The BE/MC targets cover the EMPTY-domain corner (the lookup
    //  takes a single char, so "empty" = no residue selected / the empty-set classification view),
    //  every NON-AMINO-ACID character, and the SINGLE-residue lookups that pin the scale element-by
    //  -element against the primary source.
    // ───────────────────────────────────────────────────────────────────────────

    #region Helpers — primary-source oracle (Campen 2008 Table 2 + Dunker 2001 sets)

    /// <summary>
    /// The TOP-IDP scale (Campen et al. 2008 Table 2) transcribed DIRECTLY from the primary source
    /// as restated in Disorder_Propensity.md §2.2 — the INDEPENDENT oracle for INV-01/INV-02. This
    /// is NOT read from <see cref="DisorderPredictor"/>'s own dictionary (that would be a code-echo):
    /// every value here is the published Table 2 number. Min = W = −0.884, max = P = 0.987.
    /// </summary>
    private static readonly IReadOnlyDictionary<char, double> TopIdpExpected = new Dictionary<char, double>
    {
        ['W'] = -0.884, ['F'] = -0.697, ['Y'] = -0.510, ['I'] = -0.486, ['M'] = -0.397,
        ['L'] = -0.326, ['V'] = -0.121, ['N'] = 0.007, ['C'] = 0.020, ['T'] = 0.059,
        ['A'] = 0.060, ['G'] = 0.166, ['R'] = 0.180, ['D'] = 0.192, ['H'] = 0.303,
        ['Q'] = 0.318, ['S'] = 0.341, ['K'] = 0.586, ['E'] = 0.736, ['P'] = 0.987,
    };

    /// <summary>Disorder-promoting set — Dunker et al. (2001), §2.2 — independent oracle.</summary>
    private static readonly char[] DisorderPromotingExpected = { 'A', 'R', 'G', 'Q', 'S', 'P', 'E', 'K' };

    /// <summary>Order-promoting set — Dunker et al. (2001), §2.2 — independent oracle.</summary>
    private static readonly char[] OrderPromotingExpected = { 'W', 'C', 'F', 'I', 'Y', 'V', 'L', 'N' };

    /// <summary>Ambiguous set — Dunker et al. (2001), §2.2 — independent oracle.</summary>
    private static readonly char[] AmbiguousExpected = { 'H', 'M', 'T', 'D' };

    #endregion

    #region Positive sanity — element-by-element TOP-IDP scale vs Campen 2008 Table 2 (INV-01/02)

    /// <summary>
    /// Positive sanity (INV-01, INV-02, §2.2, §7.1): GetDisorderPropensity must return the EXACT
    /// Campen et al. (2008) Table 2 value for every one of the 20 standard residues — checked
    /// element-by-element against the <see cref="TopIdpExpected"/> primary-source oracle (NOT the
    /// code's own array). The §7.1 worked-example anchors are pinned explicitly: P = 0.987 (the
    /// global maximum), W = −0.884 (the global minimum); and INV-02 — every value lies in
    /// [−0.884, 0.987] with the min uniquely at W and the max uniquely at P.
    /// </summary>
    [Test]
    public void GetDisorderPropensity_StandardResidues_MatchCampen2008Table2()
    {
        // INV-01 — element-by-element equality against the primary-source Table 2.
        foreach (char aa in StandardAminoAcids)
        {
            double expected = TopIdpExpected[aa];
            GetDisorderPropensity(aa).Should().BeApproximately(expected, 1e-12,
                $"INV-01: GetDisorderPropensity('{aa}') is the Campen 2008 Table 2 value {expected}");

            // INV-02 — within the documented scale range.
            GetDisorderPropensity(aa).Should().BeInRange(-0.884, 0.987,
                $"INV-02: the propensity of '{aa}' lies within the TOP-IDP scale range");
            double.IsNaN(GetDisorderPropensity(aa)).Should().BeFalse($"a scale value is never NaN ('{aa}')");
            double.IsInfinity(GetDisorderPropensity(aa)).Should().BeFalse($"a scale value is never infinite ('{aa}')");
        }

        // §7.1 / INV-02 anchors: the extrema are W (min) and P (max).
        GetDisorderPropensity('P').Should().BeApproximately(0.987, 1e-12,
            "P is the TOP-IDP global maximum 0.987 (§7.1)");
        GetDisorderPropensity('W').Should().BeApproximately(-0.884, 1e-12,
            "W is the TOP-IDP global minimum −0.884 (§7.1)");

        double maxVal = StandardAminoAcids.Max(GetDisorderPropensity);
        double minVal = StandardAminoAcids.Min(GetDisorderPropensity);
        maxVal.Should().BeApproximately(0.987, 1e-12, "INV-02: the scale maximum is P = 0.987");
        minVal.Should().BeApproximately(-0.884, 1e-12, "INV-02: the scale minimum is W = −0.884");
        StandardAminoAcids.Where(c => GetDisorderPropensity(c) == maxVal).Should().Equal(new[] { 'P' },
            "INV-02: the maximum is attained UNIQUELY at P");
        StandardAminoAcids.Where(c => GetDisorderPropensity(c) == minVal).Should().Equal(new[] { 'W' },
            "INV-02: the minimum is attained UNIQUELY at W");
    }

    #endregion

    #region Positive sanity — Dunker classification sets & IsDisorderPromoting (INV-03/04)

    /// <summary>
    /// Positive sanity (INV-03, INV-04, §2.2, §3.2): the three Dunker (2001) classification sets and
    /// the IsDisorderPromoting predicate. The list properties must equal the documented SORTED sets;
    /// IsDisorderPromoting(c) ⇔ c ∈ DisorderPromotingAminoAcids for every residue (INV-03); and the
    /// three sets must be pairwise DISJOINT and together COVER exactly the 20 standard residues
    /// (INV-04). Sizes are pinned to 8 / 8 / 4. Expected membership is the primary-source oracle
    /// (<see cref="DisorderPromotingExpected"/> etc.), not the code's own backing sets.
    /// </summary>
    [Test]
    public void DunkerClassification_SetsArePartition_AndPredicateAgrees()
    {
        // Property lists equal the documented SORTED Dunker sets (§3.2).
        DisorderPromotingAminoAcids.Should().Equal(DisorderPromotingExpected.OrderBy(c => c),
            "DisorderPromotingAminoAcids is the sorted Dunker disorder set {A,E,G,K,P,Q,R,S}");
        OrderPromotingAminoAcids.Should().Equal(OrderPromotingExpected.OrderBy(c => c),
            "OrderPromotingAminoAcids is the sorted Dunker order set {C,F,I,L,N,V,W,Y}");
        AmbiguousAminoAcids.Should().Equal(AmbiguousExpected.OrderBy(c => c),
            "AmbiguousAminoAcids is the sorted Dunker ambiguous set {D,H,M,T}");

        // Sizes: 8 / 8 / 4 (§2.2).
        DisorderPromotingAminoAcids.Should().HaveCount(8);
        OrderPromotingAminoAcids.Should().HaveCount(8);
        AmbiguousAminoAcids.Should().HaveCount(4);

        // INV-04 — pairwise disjoint and a full cover of the 20 standard residues.
        var disorder = DisorderPromotingAminoAcids.ToHashSet();
        var order = OrderPromotingAminoAcids.ToHashSet();
        var ambiguous = AmbiguousAminoAcids.ToHashSet();
        disorder.Overlaps(order).Should().BeFalse("INV-04: disorder ∩ order = ∅");
        disorder.Overlaps(ambiguous).Should().BeFalse("INV-04: disorder ∩ ambiguous = ∅");
        order.Overlaps(ambiguous).Should().BeFalse("INV-04: order ∩ ambiguous = ∅");
        var union = new HashSet<char>(disorder);
        union.UnionWith(order);
        union.UnionWith(ambiguous);
        union.Should().BeEquivalentTo(StandardAminoAcids.ToHashSet(),
            "INV-04: the three Dunker classes together cover exactly the 20 standard residues");

        // INV-03 — IsDisorderPromoting(c) ⇔ c ∈ DisorderPromotingAminoAcids for every residue.
        foreach (char aa in StandardAminoAcids)
            IsDisorderPromoting(aa).Should().Be(disorder.Contains(aa),
                $"INV-03: IsDisorderPromoting('{aa}') ⇔ '{aa}' ∈ DisorderPromotingAminoAcids");

        // §7.1 worked-example anchors: E disorder-promoting, W not.
        IsDisorderPromoting('E').Should().BeTrue("E is disorder-promoting (§7.1)");
        IsDisorderPromoting('W').Should().BeFalse("W is order-promoting, not disorder-promoting (§7.1)");
    }

    #endregion

    #region INV-05 — case-insensitivity: lowercase == uppercase

    /// <summary>
    /// INV-05 (§2.4, §6.1): both lookups upper-case the input first, so the lowercase form of every
    /// standard residue yields the IDENTICAL propensity and disorder-promoting verdict as its
    /// uppercase form ('p' == 'P', 'w' == 'W'). This pins the documented case folding rather than
    /// echoing the implementation.
    /// </summary>
    [Test]
    public void GetDisorderPropensity_And_IsDisorderPromoting_AreCaseInsensitive()
    {
        foreach (char upper in StandardAminoAcids)
        {
            char lower = char.ToLowerInvariant(upper);
            GetDisorderPropensity(lower).Should().Be(GetDisorderPropensity(upper),
                $"INV-05: GetDisorderPropensity('{lower}') == GetDisorderPropensity('{upper}')");
            IsDisorderPromoting(lower).Should().Be(IsDisorderPromoting(upper),
                $"INV-05: IsDisorderPromoting('{lower}') == IsDisorderPromoting('{upper}')");
        }
    }

    #endregion

    #region MC — non-amino-acid characters: propensity 0.0, not promoting, never throw

    /// <summary>
    /// Target "non-amino-acid" (MC; §3.3, §5.4 deviation #1, §6.1): any character OUTSIDE the 20
    /// standard residues must be handled DETERMINISTICALLY — GetDisorderPropensity returns exactly
    /// 0.0 and IsDisorderPromoting returns false — NEVER a KeyNotFoundException / NullReference and
    /// never NaN/±∞. We sweep the extended IUPAC ambiguity codes (B Asx, Z Glx, J Leu/Ile, O
    /// pyrrolysine, U selenocysteine), the unknown placeholder 'X', digits, punctuation, the gap /
    /// stop symbols, whitespace, the null char '\0', control chars, high unicode, and DNA/RNA
    /// nucleotide letters that are not amino acids in their own right (the only standard residue
    /// letters NOT shared with the nucleotide alphabet are excluded; B/J/O/U/X/Z and the symbols
    /// below are all out-of-table). Every one resolves to the documented neutral, with no throw.
    /// </summary>
    [Test]
    public void GetDisorderPropensity_NonAminoAcidCharacters_ReturnNeutralNeverThrow()
    {
        var nonResidues = new List<char>
        {
            'B', 'Z', 'J', 'O', 'U', 'X',          // extended / unknown amino-acid codes (out of table)
            '0', '1', '9',                          // digits
            '*', '-', '.', '?', '!', '@', '#', ' ', // stop / gap / punctuation / whitespace
            '\0', '\t', '\n', '\r',                 // null + control / whitespace
            'ñ', 'Ω', '日', ' ', '￿',     // unicode / high code points
        };

        foreach (char c in nonResidues)
        {
            string label = c < 0x20 || c == 0x7F ? $"U+{(int)c:X4}" : c.ToString();

            // No throw, exact neutral propensity, finite.
            Func<double> propAct = () => GetDisorderPropensity(c);
            double prop = propAct.Should().NotThrow($"a non-amino-acid char ('{label}') must not crash GetDisorderPropensity").Subject;
            prop.Should().Be(0.0, $"an out-of-table char ('{label}') scores the documented neutral 0.0 (§6.1)");
            double.IsNaN(prop).Should().BeFalse($"a non-amino-acid char ('{label}') never yields NaN");
            double.IsInfinity(prop).Should().BeFalse($"a non-amino-acid char ('{label}') never yields ±∞");

            // No throw, false predicate, and absent from every classification set.
            Func<bool> promAct = () => IsDisorderPromoting(c);
            bool promoting = promAct.Should().NotThrow($"a non-amino-acid char ('{label}') must not crash IsDisorderPromoting").Subject;
            promoting.Should().BeFalse($"an out-of-table char ('{label}') is not disorder-promoting (§6.1)");

            DisorderPromotingAminoAcids.Should().NotContain(c, "an out-of-table char is in no classification set");
            OrderPromotingAminoAcids.Should().NotContain(c, "an out-of-table char is in no classification set");
            AmbiguousAminoAcids.Should().NotContain(c, "an out-of-table char is in no classification set");
        }
    }

    #endregion

    #region BE — empty-domain corner: single residue & empty classification view

    /// <summary>
    /// Target "single residue" / "empty" (BE; §3.1, §6.1): the lookups take ONE char, so the
    /// boundary forms are (a) a SINGLE residue threaded through both lookups — each of the 20
    /// standard residues returns its exact primary-source propensity and the predicate matching its
    /// Dunker class (the single-residue base case of every aggregate that consumes this unit); and
    /// (b) the EMPTY-domain corner — no character selects any entry, modelled by enumerating ALL
    /// char values and asserting that EXACTLY the 20 standard residues are recognised (non-zero or
    /// in-set) while the entire rest of the char domain is the empty / neutral case. This pins both
    /// that nothing outside the 20 leaks a value and that the recognised domain is exactly those 20.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Propensity_SingleResidueExactAndEmptyDomainNeutral(CancellationToken token)
    {
        // (a) Single-residue base case — exact value + matching Dunker verdict for each of the 20.
        foreach (char aa in StandardAminoAcids)
        {
            GetDisorderPropensity(aa).Should().BeApproximately(TopIdpExpected[aa], 1e-12,
                $"a single residue '{aa}' returns its exact Table 2 propensity");
            IsDisorderPromoting(aa).Should().Be(DisorderPromotingExpected.Contains(aa),
                $"a single residue '{aa}' returns its Dunker disorder-promoting verdict");
        }

        // (b) Empty-domain corner — sweep the WHOLE char range; only the 20 standard residues are
        // recognised, EVERY other char is the neutral/empty case (0.0 propensity, not promoting).
        var recognised = new HashSet<char>();
        for (int code = 0; code <= char.MaxValue; code++)
        {
            char c = (char)code;
            double prop = GetDisorderPropensity(c);
            bool promoting = IsDisorderPromoting(c);

            double.IsNaN(prop).Should().BeFalse($"propensity is never NaN (U+{code:X4})");
            double.IsInfinity(prop).Should().BeFalse($"propensity is never ±∞ (U+{code:X4})");

            // A char is "recognised" iff it (case-folded) is one of the 20 standard residues.
            bool isStandard = StandardAminoAcids.Contains(char.ToUpperInvariant(c));
            if (isStandard)
            {
                recognised.Add(char.ToUpperInvariant(c));
            }
            else
            {
                prop.Should().Be(0.0, $"a non-residue char (U+{code:X4}) scores the neutral 0.0 — empty-domain case");
                promoting.Should().BeFalse($"a non-residue char (U+{code:X4}) is not disorder-promoting");
            }

            if ((code & 0x1FFF) == 0)
                token.ThrowIfCancellationRequested();
        }

        recognised.Should().BeEquivalentTo(StandardAminoAcids.ToHashSet(),
            "exactly the 20 standard residues (and their lowercase forms) are recognised; the rest of the char domain is the empty/neutral case");
    }

    #endregion

    #region BE — randomized boundary sweep: never crash / hang / NaN, deterministic

    /// <summary>
    /// Headline no-crash / no-hang / no-NaN sweep over arbitrary and adversarial single-char input
    /// (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"). Across fixed seeds we draw random chars
    /// from the full BMP — a mix of standard residues, extended codes, digits, punctuation, control
    /// chars and unicode — and require that BOTH lookups never throw, never return NaN/±∞, return a
    /// value in the documented scale range OR exactly 0.0 (for out-of-table chars), and agree with
    /// the case-insensitive primary-source oracle. Re-running on the same char gives the IDENTICAL
    /// result (determinism). [CancelAfter] guards against any regression turning an O(1) lookup into
    /// a hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Propensity_RandomChars_NeverCrashNaN_MatchOracleDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 7, 31, 137, 2026 })
        {
            var rng = new Random(seed);
            for (int iter = 0; iter < 5000; iter++)
            {
                char c = (char)rng.Next(0, char.MaxValue + 1);
                char upper = char.ToUpperInvariant(c);

                Func<double> propAct = () => GetDisorderPropensity(c);
                double prop = propAct.Should().NotThrow($"random char U+{(int)c:X4} must not crash propensity lookup").Subject;
                Func<bool> promAct = () => IsDisorderPromoting(c);
                bool promoting = promAct.Should().NotThrow($"random char U+{(int)c:X4} must not crash predicate").Subject;

                double.IsNaN(prop).Should().BeFalse($"random char U+{(int)c:X4} never yields NaN");
                double.IsInfinity(prop).Should().BeFalse($"random char U+{(int)c:X4} never yields ±∞");

                // Oracle agreement (case-insensitive): in-table → exact Table 2 value & Dunker verdict;
                // out-of-table → exactly 0.0 / false.
                if (TopIdpExpected.TryGetValue(upper, out double expectedProp))
                {
                    prop.Should().BeApproximately(expectedProp, 1e-12,
                        $"a recognised residue (folded '{upper}') returns its Table 2 value");
                    prop.Should().BeInRange(-0.884, 0.987, "INV-02: a recognised residue is within the scale range");
                    promoting.Should().Be(DisorderPromotingExpected.Contains(upper),
                        $"a recognised residue (folded '{upper}') matches its Dunker verdict");
                }
                else
                {
                    prop.Should().Be(0.0, $"an out-of-table char (folded '{upper}') scores the neutral 0.0");
                    promoting.Should().BeFalse($"an out-of-table char (folded '{upper}') is not disorder-promoting");
                }

                // Determinism — the same char yields the identical result.
                GetDisorderPropensity(c).Should().Be(prop, "GetDisorderPropensity is deterministic for a fixed char");
                IsDisorderPromoting(c).Should().Be(promoting, "IsDisorderPromoting is deterministic for a fixed char");

                if ((iter & 0x3FF) == 0)
                    token.ThrowIfCancellationRequested();
            }
        }
    }

    #endregion

    #endregion
}
