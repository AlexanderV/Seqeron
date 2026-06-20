using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.DisorderPredictor;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the ProteinPred area — intrinsic disorder PREDICTION
/// (the per-residue TOP-IDP disorder profile, DISORDER-PRED-001) and the
/// segmentation of that profile into contiguous disordered REGIONS
/// (DISORDER-REGION-001).
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
}
