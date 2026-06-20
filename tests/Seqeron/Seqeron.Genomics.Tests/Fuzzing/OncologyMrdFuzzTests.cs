using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology minimal/molecular residual disease (MRD) area — ONCO-MRD-001.
/// The unit under test is tumour-informed MRD detection: tracking a patient-specific panel of somatic
/// markers in a follow-up plasma sample and calling the sample MRD-positive/negative. Implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs as the public entry points of the
/// <c>Minimal/Molecular Residual Disease (ONCO-MRD-001)</c> region:
///   • <see cref="OncologyAnalyzer.DetectMRD(IEnumerable{OncologyAnalyzer.TumorMarker}, int, int, int)"/> —
///       panel-level call: count detected markers D, call positive ⟺ D ≥ τ (default 2); also reports IMAF
///       and the panel-level Poisson detection probability p;
///   • <see cref="OncologyAnalyzer.TrackVariantsOverTime(IEnumerable{IEnumerable{OncologyAnalyzer.TumorMarker}}, int, int, int)"/> —
///       longitudinal per-timepoint MRD with the earliest-positive index;
///   • <see cref="OncologyAnalyzer.IsVariantDetected(OncologyAnalyzer.TumorMarker, int)"/> —
///       per-locus presence (alt reads ≥ minSupportingReads).
/// (ctDNA tumour-fraction entry points — CalculateTumorFraction / CtDnaDetectionProbability — are
///  ONCO-CTDNA-001, row 111 in OncologyCtdnaFuzzTests.cs, and are intentionally OUT OF SCOPE here. This
///  file is scoped to the residual-disease CALL across a tracked panel.)
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts that the code NEVER fails
/// in an undisciplined way: no hang, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / Overflow / NaN). Every input must resolve to EITHER a well-defined, theory-correct MRD
/// call OR a *documented, intentional* outcome (an <see cref="ArgumentNullException"/> for a null panel,
/// an <see cref="ArgumentException"/> for an EMPTY panel, an <see cref="ArgumentOutOfRangeException"/> for a
/// malformed threshold). For the MRD call the headline hazards are:
///   • a DivideByZero / NaN IMAF when the panel has ZERO total depth (Σt = 0) — the contract defines
///     IMAF = 0 when Σt = 0 (§3.3, §6.1, algorithm step 3), so a fraction is still produced, never an
///     exception or NaN;
///   • a DivideByZero / nonsense aggregation over an EMPTY panel — empty is a documented throw, the
///     accumulators never divide by a zero marker count (§6.1);
///   • aggregation OVERFLOW when summing alt / total reads across many high-depth loci — the contract sums
///     into wide accumulators, so Σa / Σt must stay finite and IMAF ∈ [0, 1] (INV-03);
///   • a wrong threshold direction (the ≥ τ rule firing at D = τ−1, or an off-by-one at the boundary) —
///     MRD-positive ⟺ D ≥ τ exactly (INV-01), with 0 ≤ D ≤ m (INV-02).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-MRD-001 — minimal/molecular residual disease detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 112.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 112): "no tracked variants, all-detected, single low-VAF read".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Mapping of the BE targets onto the documented contract:
///   • "no tracked variants" ⇒ an EMPTY patient-specific panel (m = 0). Documented guard: nothing to
///       interrogate ⇒ ArgumentException (§3.3, §6.1). No DivideByZero on the Σ/m aggregation, no NaN IMAF.
///   • "all-detected" ⇒ every tracked marker has supporting mutant reads ⇒ D = m. With m ≥ τ this is the
///       strong positive case (MRD-positive, INV-01); D = m is the documented upper bound (INV-02). No
///       overflow on aggregated alt/total counts; IMAF stays in [0, 1].
///   • "single low-VAF read" ⇒ one marker with a single mutant read at very low VAF (alt = 1 out of a deep
///       total). At default r_min = 1 that locus IS detected (D = 1), but D = 1 &lt; τ = 2 ⇒ MRD-negative
///       (§6.1 "Exactly 1 detected ⇒ negative"); a tiny depth must not DivideByZero and the IMAF must be a
///       finite, near-zero fraction, not nonsense.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// MRD_Detection.md (docs/algorithms/Oncology/MRD_Detection.md):
///   • Per-locus detection: marker i detected ⟺ a_i ≥ r_min (default r_min = 1) (§2.2).
///   • Panel call: D = Σ 1[a_i ≥ r_min]; MRD-positive ⟺ D ≥ τ (default τ = 2) (§2.2, INV-01).
///   • Bounds: 0 ≤ D ≤ m (INV-02); IMAF = Σa / Σt ∈ [0, 1] (INV-03); 0 when Σt = 0.
///   • Panel Poisson p = 1 − e^(−n·IMAF·m) ∈ [0, 1], non-decreasing in m (INV-04).
///   • Longitudinal: order preserved; FirstPositiveIndex = earliest positive timepoint, or −1 (INV-05).
///   • null panel ⇒ ArgumentNullException; empty panel ⇒ ArgumentException; positivityThreshold &lt; 1,
///       minSupportingReads &lt; 1, genomeEquivalents &lt; 0 ⇒ ArgumentOutOfRangeException (§3.3). Negative
///       read counts are clamped to 0 when summing IMAF.
///   • Edge cases (§6.1): exactly 1 detected ⇒ negative; 0 detected ⇒ negative; all total reads 0 ⇒
///       IMAF 0, p 0; empty panel ⇒ ArgumentException; custom τ = 1 ⇒ single detected ⇒ positive.
///   • Worked example (§7.1): panel of 3 markers (3/200, 1/150, 0/180) ⇒ D = 2 ⇒ Positive,
///       IMAF = (3+1+0)/(200+150+180) = 4/530 ≈ 0.0075472.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyMrdFuzzTests
{
    // ── Well-formed-MRD-result assertion helper ──────────────────────────────
    // Pins the documented numeric contract on EVERY accepted MRD result:
    //   • IMAF finite, never NaN/Infinity (no 0/0 on Σt = 0, no DivideByZero) and inside [0, 1] (INV-03);
    //   • DetectionProbability finite and inside [0, 1] (INV-04);
    //   • 0 ≤ D ≤ m (INV-02) and the binary Status agrees with the threshold rule (INV-01).
    // This is what stops a fuzz test from rubber-stamping a NaN IMAF, a probability that escapes [0, 1],
    // a detected count outside [0, m], or a Status that disagrees with D ≥ τ.
    private static void AssertWellFormedMrd(MrdResult r, int positivityThreshold)
    {
        double.IsNaN(r.IntegratedMutantAlleleFraction).Should().BeFalse(
            "IMAF must never be NaN (no 0/0 on Σt = 0, no DivideByZero on tiny depth)");
        double.IsInfinity(r.IntegratedMutantAlleleFraction).Should().BeFalse("IMAF must be finite");
        r.IntegratedMutantAlleleFraction.Should().BeInRange(0.0, 1.0, "IMAF = Σa / Σt ∈ [0, 1] (INV-03)");

        double.IsNaN(r.DetectionProbability).Should().BeFalse("panel Poisson p must never be NaN");
        double.IsInfinity(r.DetectionProbability).Should().BeFalse("panel Poisson p must be finite");
        r.DetectionProbability.Should().BeInRange(0.0, 1.0, "p = 1 − e^(−n·f·m) ∈ [0, 1] (INV-04)");

        r.DetectedVariantCount.Should().BeInRange(0, r.TrackedVariantCount, "0 ≤ D ≤ m (INV-02)");
        bool expectedPositive = r.DetectedVariantCount >= positivityThreshold;
        (r.Status == MrdStatus.Positive).Should().Be(
            expectedPositive, "MRD-positive ⟺ D ≥ τ (INV-01)");
    }

    // A tracked patient-specific marker expressed as plasma read evidence: `alt` mutant supporting reads
    // out of `total` covering reads at one locus.
    private static TumorMarker Marker(int alt, int total)
        => new(
            Chromosome: "chr1",
            Position: 1,
            ReferenceAllele: "A",
            AlternateAllele: "T",
            PlasmaAltReads: alt,
            PlasmaTotalReads: total);

    #region ONCO-MRD-001 — Positive sanity (documented threshold rule / IMAF on hand-built examples)

    [Test]
    public void DetectMRD_DocWorkedExample_TwoDetectedPositiveAndExactImaf()
    {
        // Docs §7.1: panel of 3 markers (3/200, 1/150, 0/180). Loci 1 and 2 have alt ≥ 1 ⇒ D = 2 ⇒
        // Positive (D ≥ τ = 2). IMAF = (3+1+0)/(200+150+180) = 4/530 ≈ 0.0075472. Pins the headline rule.
        var panel = new[]
        {
            Marker(alt: 3, total: 200),
            Marker(alt: 1, total: 150),
            Marker(alt: 0, total: 180),
        };

        MrdResult r = DetectMRD(panel);

        AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        r.Status.Should().Be(MrdStatus.Positive);
        r.DetectedVariantCount.Should().Be(2);
        r.TrackedVariantCount.Should().Be(3);
        r.IntegratedMutantAlleleFraction.Should().BeApproximately(4.0 / 530.0, 1e-12);
    }

    [Test]
    public void DetectMRD_NoMutantSupportAnywhere_NegativeAndZeroImaf()
    {
        // A panel with clear depth but NO mutant support at any locus ⇒ D = 0 ⇒ MRD-negative, IMAF = 0.
        // The documented "no signal" counterpart to the positive call.
        var panel = new[]
        {
            Marker(alt: 0, total: 300),
            Marker(alt: 0, total: 250),
            Marker(alt: 0, total: 400),
        };

        MrdResult r = DetectMRD(panel);

        AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        r.Status.Should().Be(MrdStatus.Negative);
        r.DetectedVariantCount.Should().Be(0);
        r.IntegratedMutantAlleleFraction.Should().Be(0.0);
    }

    [Test]
    public void DetectMRD_ThresholdBoundary_OneBelowNegativeAtTauPositive()
    {
        // INV-01 boundary: with τ = 2, exactly 1 detected ⇒ Negative (§6.1), exactly 2 detected ⇒ Positive.
        // Guards against an off-by-one / wrong threshold direction (>, ≥, <).
        var oneDetected = new[] { Marker(alt: 5, total: 100), Marker(alt: 0, total: 100) };
        var twoDetected = new[] { Marker(alt: 5, total: 100), Marker(alt: 1, total: 100) };

        DetectMRD(oneDetected).Status.Should().Be(MrdStatus.Negative, "D = 1 < τ = 2");
        DetectMRD(twoDetected).Status.Should().Be(MrdStatus.Positive, "D = 2 ≥ τ = 2");
    }

    [Test]
    public void DetectMRD_PanelPoisson_NonDecreasingInPanelSize()
    {
        // INV-04: p = 1 − e^(−n·IMAF·m) is non-decreasing in m. Two panels with the SAME IMAF but more
        // tracked markers ⇒ the larger panel has p ≥ the smaller. Pins the documented Poisson direction.
        var small = new[] { Marker(alt: 5, total: 100), Marker(alt: 5, total: 100) };           // IMAF 0.05, m 2
        var large = new[]
        {
            Marker(alt: 5, total: 100), Marker(alt: 5, total: 100),
            Marker(alt: 5, total: 100), Marker(alt: 5, total: 100),                              // IMAF 0.05, m 4
        };

        MrdResult rs = DetectMRD(small, genomeEquivalents: 1000);
        MrdResult rl = DetectMRD(large, genomeEquivalents: 1000);

        rs.IntegratedMutantAlleleFraction.Should().BeApproximately(0.05, 1e-12);
        rl.IntegratedMutantAlleleFraction.Should().BeApproximately(0.05, 1e-12);
        rl.DetectionProbability.Should().BeGreaterThanOrEqualTo(
            rs.DetectionProbability - 1e-12, "p non-decreasing in m (INV-04)");
    }

    #endregion

    #region ONCO-MRD-001 — BE: no tracked variants (empty panel ⇒ documented throw, no 0/0 aggregation)

    [Test]
    public void DetectMRD_Empty_ThrowsArgumentException_NotNaN()
    {
        // "no tracked variants": an empty patient-specific panel (m = 0) ⇒ nothing to interrogate ⇒
        // documented ArgumentException (§6.1), never a silent NaN/0-division on the Σ/m aggregation.
        Action act = () => DetectMRD(Array.Empty<TumorMarker>());
        act.Should().Throw<ArgumentException>().Which.Should().NotBeOfType<ArgumentNullException>();
    }

    [Test]
    public void DetectMRD_Null_ThrowsArgumentNullException()
    {
        ((Action)(() => DetectMRD(null!))).Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void TrackVariantsOverTime_NullSeriesOrNullPanel_Throw()
    {
        // The longitudinal entry validates the series and delegates each panel to DetectMRD: a null series
        // and a null timepoint panel are both documented throws (§3.3).
        ((Action)(() => TrackVariantsOverTime(null!))).Should().Throw<ArgumentNullException>();

        var seriesWithNullPanel = new IEnumerable<TumorMarker>[] { new[] { Marker(1, 100) }, null! };
        ((Action)(() => TrackVariantsOverTime(seriesWithNullPanel))).Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void TrackVariantsOverTime_EmptyTimepointPanel_ThrowsArgumentException()
    {
        // An empty panel embedded in the series ⇒ ArgumentException via the DetectMRD delegation (§3.3),
        // not a silently-skipped or NaN timepoint.
        var series = new IEnumerable<TumorMarker>[] { new[] { Marker(1, 100) }, Array.Empty<TumorMarker>() };
        ((Action)(() => TrackVariantsOverTime(series)))
            .Should().Throw<ArgumentException>().Which.Should().NotBeOfType<ArgumentNullException>();
    }

    [Test]
    public void DetectMRD_MalformedThresholds_Throw()
    {
        // §3.3 domain limits — each malformed scalar is a documented throw, not a nonsense call.
        var panel = new[] { Marker(1, 100) };
        ((Action)(() => DetectMRD(panel, positivityThreshold: 0)))
            .Should().Throw<ArgumentOutOfRangeException>("positivityThreshold < 1");
        ((Action)(() => DetectMRD(panel, minSupportingReads: 0)))
            .Should().Throw<ArgumentOutOfRangeException>("minSupportingReads < 1");
        ((Action)(() => DetectMRD(panel, genomeEquivalents: -1)))
            .Should().Throw<ArgumentOutOfRangeException>("genomeEquivalents < 0");
    }

    #endregion

    #region ONCO-MRD-001 — BE: all-detected (every marker supported ⇒ D = m, positive, no overflow)

    [Test]
    [CancelAfter(10_000)]
    public void DetectMRD_AllMarkersSupported_AllDetectedAndPositive()
    {
        // "all-detected": every tracked marker has ≥ 1 mutant read ⇒ D = m (the INV-02 upper bound) ⇒
        // with m ≥ τ = 2 the sample is MRD-positive. IMAF stays finite and in [0, 1].
        var rng = new Random(112_001);
        for (int i = 0; i < 300; i++)
        {
            int m = rng.Next(2, 16); // ≥ τ so all-detected ⇒ positive
            var panel = Enumerable.Range(0, m)
                .Select(_ =>
                {
                    int total = rng.Next(2, 5000);
                    int alt = rng.Next(1, total + 1); // ≥ 1 ⇒ detected at r_min = 1
                    return Marker(alt, total);
                })
                .ToArray();

            MrdResult r = DetectMRD(panel);

            AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
            r.DetectedVariantCount.Should().Be(m, "every marker supported ⇒ D = m (INV-02)");
            r.TrackedVariantCount.Should().Be(m);
            r.Status.Should().Be(MrdStatus.Positive, "D = m ≥ τ = 2 ⇒ positive (INV-01)");
        }
    }

    [Test]
    [CancelAfter(10_000)]
    public void DetectMRD_AllDetectedDeepHighDepth_NoOverflowImafInRange()
    {
        // Aggregation OVERFLOW guard: many loci each near int.MaxValue depth. The accumulators must sum
        // into a wide enough type so IMAF = Σa / Σt stays finite and in [0, 1] (no int overflow / no NaN).
        var rng = new Random(112_002);
        for (int i = 0; i < 50; i++)
        {
            int m = rng.Next(8, 16);
            var panel = Enumerable.Range(0, m)
                .Select(_ =>
                {
                    int total = int.MaxValue - rng.Next(0, 1000); // ~2.1e9 each ⇒ Σt overflows a 32-bit int
                    int alt = rng.Next(1, total);                 // some mutant support, < total
                    return Marker(alt, total);
                })
                .ToArray();

            MrdResult r = DetectMRD(panel, genomeEquivalents: int.MaxValue);

            AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
            r.DetectedVariantCount.Should().Be(m, "all supported ⇒ D = m even at extreme depth");
        }
    }

    [Test]
    public void DetectMRD_AllDetectedCustomThresholdOne_SingleMarkerPositive()
    {
        // Custom τ = 1 (§6.1): a single all-detected marker ⇒ D = 1 ≥ τ = 1 ⇒ positive. Exercises the
        // parameterised threshold direction at the m = 1 boundary.
        MrdResult r = DetectMRD(new[] { Marker(alt: 4, total: 200) }, positivityThreshold: 1);

        AssertWellFormedMrd(r, 1);
        r.DetectedVariantCount.Should().Be(1);
        r.Status.Should().Be(MrdStatus.Positive, "D = 1 ≥ τ = 1");
    }

    #endregion

    #region ONCO-MRD-001 — BE: single low-VAF read (one alt read at tiny VAF ⇒ detected but D=1 < τ ⇒ negative)

    [Test]
    public void DetectMRD_SingleLowVafReadOneMarker_DetectedButBelowThresholdNegative()
    {
        // "single low-VAF read": ONE marker with a single mutant read out of deep coverage (VAF ≈ 1/5000).
        // At r_min = 1 the locus IS detected (D = 1), but D = 1 < τ = 2 ⇒ MRD-negative (§6.1). The tiny
        // VAF must produce a finite, near-zero IMAF, not nonsense, and no DivideByZero.
        MrdResult r = DetectMRD(new[] { Marker(alt: 1, total: 5000) });

        AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        r.DetectedVariantCount.Should().Be(1, "alt = 1 ≥ r_min = 1 ⇒ detected");
        r.Status.Should().Be(MrdStatus.Negative, "D = 1 < τ = 2 (§6.1 exactly-1-detected ⇒ negative)");
        r.IntegratedMutantAlleleFraction.Should().BeApproximately(1.0 / 5000.0, 1e-12);
    }

    [Test]
    public void DetectMRD_SingleLowVafReadDepthOne_NoDivideByZero()
    {
        // The smallest possible support at the smallest possible depth: alt = 1, total = 1 (VAF 1.0 at one
        // locus). Detected (D = 1) but below τ = 2 ⇒ negative. The tiny denominator must not DivideByZero,
        // and IMAF = 1/1 = 1.0 stays inside the documented [0, 1] (INV-03).
        MrdResult r = DetectMRD(new[] { Marker(alt: 1, total: 1) });

        AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        r.DetectedVariantCount.Should().Be(1);
        r.Status.Should().Be(MrdStatus.Negative);
        r.IntegratedMutantAlleleFraction.Should().Be(1.0);
    }

    [Test]
    public void DetectMRD_SingleLowVafReadCustomMinSupporting_BelowCutoffNotDetected()
    {
        // With r_min raised to 2, a single low-VAF read (alt = 1) falls BELOW the per-locus cutoff ⇒ NOT
        // detected ⇒ D = 0 ⇒ negative. Confirms the documented configurable r_min direction (§2.2).
        var panel = new[] { Marker(alt: 1, total: 5000), Marker(alt: 1, total: 4000) };

        MrdResult rDefault = DetectMRD(panel);                       // r_min = 1
        MrdResult rStrict = DetectMRD(panel, minSupportingReads: 2); // r_min = 2

        rDefault.DetectedVariantCount.Should().Be(2, "alt = 1 ≥ r_min = 1");
        rStrict.DetectedVariantCount.Should().Be(0, "alt = 1 < r_min = 2 ⇒ neither locus detected");
        rStrict.Status.Should().Be(MrdStatus.Negative);
        AssertWellFormedMrd(rStrict, DefaultMrdPositivityThreshold);
    }

    [Test]
    public void IsVariantDetected_SingleLowVafReadAtCutoffBoundary()
    {
        // Per-locus presence boundary: alt = 1 is detected at r_min = 1, not at r_min = 2.
        IsVariantDetected(Marker(alt: 1, total: 5000)).Should().BeTrue();
        IsVariantDetected(Marker(alt: 1, total: 5000), minSupportingReads: 2).Should().BeFalse();
        IsVariantDetected(Marker(alt: 0, total: 5000)).Should().BeFalse();
        ((Action)(() => IsVariantDetected(Marker(1, 100), minSupportingReads: 0)))
            .Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-MRD-001 — BE: zero / negative read counts (IMAF 0 on Σt = 0, negatives clamped, no NaN)

    [Test]
    public void DetectMRD_AllZeroDepth_ImafZeroAndProbabilityZero()
    {
        // §6.1 "All total reads = 0 ⇒ IMAF = 0, p = 0": a panel where every locus has zero coverage. The
        // Σa / Σt = 0/0 must resolve to the documented 0, never a NaN or DivideByZero. D = 0 ⇒ negative.
        var panel = new[] { Marker(0, 0), Marker(0, 0), Marker(0, 0) };

        MrdResult r = DetectMRD(panel, genomeEquivalents: 10_000);

        AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        r.IntegratedMutantAlleleFraction.Should().Be(0.0);
        r.DetectionProbability.Should().Be(0.0, "p = 1 − e^0 = 0 when IMAF = 0");
        r.DetectedVariantCount.Should().Be(0);
        r.Status.Should().Be(MrdStatus.Negative);
    }

    [Test]
    [CancelAfter(10_000)]
    public void DetectMRD_NegativeReadCounts_ClampedToZeroInImaf_NoNaN()
    {
        // §3.3: negative read counts are clamped to 0 when summing IMAF. Fuzz negative alt/total across a
        // panel and assert the result is still well-formed (IMAF ∈ [0, 1], finite) — never a negative or
        // NaN IMAF. Per-locus detection on a negative alt is simply "not detected" (alt < r_min).
        var rng = new Random(112_003);
        for (int i = 0; i < 200; i++)
        {
            int m = rng.Next(1, 12);
            var panel = Enumerable.Range(0, m)
                .Select(_ => Marker(alt: rng.Next(-50, 50), total: rng.Next(-50, 5000)))
                .ToArray();

            MrdResult r = DetectMRD(panel);

            AssertWellFormedMrd(r, DefaultMrdPositivityThreshold);
        }
    }

    #endregion

    #region ONCO-MRD-001 — Longitudinal tracking (order preserved, earliest-positive index, INV-05)

    [Test]
    public void TrackVariantsOverTime_FirstPositiveIndex_EarliestPositiveTimepoint()
    {
        // INV-05: timepoints are scanned in order; FirstPositiveIndex is the earliest positive (or −1).
        // Build a series negative, negative, positive, positive ⇒ first positive at index 2.
        var negative = new[] { Marker(0, 100), Marker(0, 100) };          // D = 0
        var positive = new[] { Marker(5, 100), Marker(5, 100) };          // D = 2
        var series = new IEnumerable<TumorMarker>[] { negative, negative, positive, positive };

        MrdLongitudinalResult lr = TrackVariantsOverTime(series);

        lr.Timepoints.Should().HaveCount(4);
        lr.Timepoints.Select(t => t.TimepointIndex).Should().Equal(0, 1, 2, 3);
        lr.FirstPositiveIndex.Should().Be(2);
        lr.Timepoints[2].Result.Status.Should().Be(MrdStatus.Positive);
        foreach (MrdTimepoint t in lr.Timepoints)
        {
            AssertWellFormedMrd(t.Result, DefaultMrdPositivityThreshold);
        }
    }

    [Test]
    public void TrackVariantsOverTime_NeverPositive_FirstPositiveIndexMinusOne()
    {
        // A series that never reaches MRD-positivity ⇒ FirstPositiveIndex = −1 (INV-05).
        var negative = new[] { Marker(0, 100), Marker(0, 100) };
        var series = new IEnumerable<TumorMarker>[] { negative, negative, negative };

        MrdLongitudinalResult lr = TrackVariantsOverTime(series);

        lr.FirstPositiveIndex.Should().Be(-1);
    }

    [Test]
    public void TrackVariantsOverTime_EmptySeries_NoTimepointsFirstPositiveMinusOne()
    {
        // An empty SERIES of timepoints (distinct from an empty panel) ⇒ no timepoints, no positive ⇒
        // FirstPositiveIndex = −1. No throw, no NaN.
        MrdLongitudinalResult lr = TrackVariantsOverTime(Array.Empty<IEnumerable<TumorMarker>>());

        lr.Timepoints.Should().BeEmpty();
        lr.FirstPositiveIndex.Should().Be(-1);
    }

    #endregion

    #region ONCO-MRD-001 — Invariant sweep (D ∈ [0, m]; Status ⟺ D ≥ τ; IMAF, p ∈ [0, 1])

    [Test]
    [CancelAfter(20_000)]
    public void DetectMRD_RandomPanels_AlwaysSatisfyDocumentedInvariants()
    {
        // INV-01..INV-04 sweep: for any random panel and threshold, D counts exactly the markers with
        // alt ≥ r_min, Status ⟺ D ≥ τ, IMAF = Σmax(0,alt) / Σmax(0,total) ∈ [0, 1], p ∈ [0, 1]. Never NaN.
        var rng = new Random(112_004);
        for (int i = 0; i < 500; i++)
        {
            int m = rng.Next(1, 16);
            int tau = rng.Next(1, m + 2);
            int rMin = rng.Next(1, 4);
            int n = rng.Next(0, 50_000);

            int[] alts = new int[m];
            int[] totals = new int[m];
            var panel = new TumorMarker[m];
            for (int j = 0; j < m; j++)
            {
                int total = rng.Next(0, 6000);
                int alt = total == 0 ? rng.Next(0, 3) : rng.Next(0, total + 1);
                alts[j] = alt;
                totals[j] = total;
                panel[j] = Marker(alt, total);
            }

            MrdResult r = DetectMRD(panel, positivityThreshold: tau, minSupportingReads: rMin, genomeEquivalents: n);

            AssertWellFormedMrd(r, tau);

            int expectedD = alts.Count(a => a >= rMin);
            long sumAlt = alts.Sum(a => (long)Math.Max(0, a));
            long sumTotal = totals.Sum(t => (long)Math.Max(0, t));
            double expectedImaf = sumTotal == 0 ? 0.0 : (double)sumAlt / sumTotal;

            r.DetectedVariantCount.Should().Be(expectedD, "D = Σ 1[a_i ≥ r_min] (INV-02)");
            r.TrackedVariantCount.Should().Be(m);
            r.IntegratedMutantAlleleFraction.Should().BeApproximately(expectedImaf, 1e-12);
            (r.Status == MrdStatus.Positive).Should().Be(expectedD >= tau, "Status ⟺ D ≥ τ (INV-01)");
        }
    }

    #endregion
}
