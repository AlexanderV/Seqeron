// ONCO-MRD-001 — Minimal/Molecular Residual Disease Detection
// Evidence: docs/Evidence/ONCO-MRD-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MRD-001.md
// Source: Reinert et al. (2019). JAMA Oncol 5(8):1124-1131 (tumour-informed MRD; ≥2 of tracked SNVs => positive).
//         Natera Signatera white paper (2020): 16 tracked SNVs; "at least two" => ctDNA-positive; p = 1 - e^(-nfm).
//         Wan et al. (2020). Sci Transl Med 12(548):eaaz8084 (INVAR integrated mutant allele fraction / IMAF).
//         Quoted positivity rule: PMC9265001 Table 1.

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_DetectMRD_Tests
{
    // Helper: a tracked marker with given plasma alt/total reads. Position/alleles are display-only here.
    private static OncologyAnalyzer.TumorMarker Marker(int altReads, int totalReads) =>
        new("1", 100, "A", "T", altReads, totalReads);

    #region DetectMRD

    // M1 — Reinert 2019 / PMC9265001 Table 1: >=2 detected of the tracked variants => ctDNA-positive.
    [Test]
    public void DetectMRD_TwoOfThreeDetected_PositiveCall()
    {
        // Arrange: 3-marker panel, two with mutant reads (detected), one without.
        var panel = new[] { Marker(5, 200), Marker(3, 200), Marker(0, 200) };

        // Act
        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(2),
                "Two markers have alt reads >= 1, so two are detected.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "2 detected >= threshold 2 => MRD-positive (Reinert 2019; PMC9265001 Table 1).");
        });
    }

    // M2 — exactly 1 detected is BELOW the >=2 rule => MRD-negative (documented corner case).
    [Test]
    public void DetectMRD_OneOfThreeDetected_NegativeCall()
    {
        var panel = new[] { Marker(4, 200), Marker(0, 200), Marker(0, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(1),
                "Only one marker has alt reads => one detected.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "1 detected < threshold 2 => MRD-negative (a single variant is not a positive call).");
        });
    }

    // M3 — no markers detected => MRD-negative.
    [Test]
    public void DetectMRD_ZeroDetected_NegativeCall()
    {
        var panel = new[] { Marker(0, 200), Marker(0, 150), Marker(0, 180) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(0), "No marker has alt reads.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "0 detected < 2 => MRD-negative.");
            Assert.That(result.IntegratedMutantAlleleFraction, Is.EqualTo(0.0).Within(1e-12),
                "IMAF = 0/Sigma(total) = 0 when no alt reads.");
        });
    }

    // M4 — all 3 detected => MRD-positive, detected count = panel size.
    [Test]
    public void DetectMRD_ThreeOfThreeDetected_PositiveCall()
    {
        var panel = new[] { Marker(5, 200), Marker(3, 200), Marker(2, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(3), "All three markers have alt reads.");
            Assert.That(result.TrackedVariantCount, Is.EqualTo(3), "Panel size is 3.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "3 detected >= 2 => MRD-positive.");
        });
    }

    // M5 — Signatera tracks 16 patient-specific SNVs; TrackedVariantCount reports panel size.
    [Test]
    public void DetectMRD_SixteenMarkers_ReportsTrackedCount()
    {
        // 16 markers: first two detected, rest undetected => positive with detected=2, tracked=16.
        var panel = new List<OncologyAnalyzer.TumorMarker> { Marker(3, 1000), Marker(2, 1000) };
        for (int i = 0; i < 14; i++)
        {
            panel.Add(Marker(0, 1000));
        }

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.Multiple(() =>
        {
            Assert.That(result.TrackedVariantCount, Is.EqualTo(16),
                "Signatera panel tracks 16 patient-specific SNVs (white paper).");
            Assert.That(result.DetectedVariantCount, Is.EqualTo(2), "Two markers detected.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "2 of 16 detected => MRD-positive.");
        });
    }

    // M6 — IMAF = depth-weighted (read-pooled) mean VAF across loci (Wan 2020): (3+1+0)/(200+150+180)=4/530.
    [Test]
    public void DetectMRD_ImafWorkedExample_DepthWeightedMeanVaf()
    {
        var panel = new[] { Marker(3, 200), Marker(1, 150), Marker(0, 180) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.Multiple(() =>
        {
            // 4/530 = 0.007547169811320755 (hand-derived: Sigma alt = 4, Sigma total = 530).
            Assert.That(result.IntegratedMutantAlleleFraction, Is.EqualTo(4.0 / 530.0).Within(1e-12),
                "IMAF = Sigma(alt)/Sigma(total) = 4/530 (INVAR integrated mutant allele fraction).");
            Assert.That(result.IntegratedMutantAlleleFraction, Is.EqualTo(0.007547169811320755).Within(1e-12),
                "Numeric IMAF value of 4/530.");
            Assert.That(result.DetectedVariantCount, Is.EqualTo(2),
                "Loci 1 and 2 have alt reads => detected = 2 => positive.");
        });
    }

    // M7 — panel-level Poisson p = 1 - e^(-n*f*m). n=1000, f=0.001 (IMAF), m=16 => lambda=16.
    [Test]
    public void DetectMRD_PanelDetectionProbability_PoissonM16()
    {
        // 16 markers each with VAF exactly 0.001 (1 alt / 1000 total) => IMAF = 0.001, m = 16.
        var panel = new List<OncologyAnalyzer.TumorMarker>();
        for (int i = 0; i < 16; i++)
        {
            panel.Add(Marker(1, 1000));
        }

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel, genomeEquivalents: 1000);

        Assert.Multiple(() =>
        {
            Assert.That(result.IntegratedMutantAlleleFraction, Is.EqualTo(0.001).Within(1e-12),
                "IMAF = 16/16000 = 0.001.");
            // 1 - e^(-16) = 0.9999998874648253 (Signatera white paper Fig 2: p = 1 - e^(-nfm)).
            Assert.That(result.DetectionProbability, Is.EqualTo(1.0 - Math.Exp(-16.0)).Within(1e-12),
                "Panel Poisson p = 1 - e^(-n*f*m) with n=1000, f=0.001, m=16 => lambda=16.");
            Assert.That(result.DetectionProbability, Is.EqualTo(0.9999998874648253).Within(1e-12),
                "Numeric value of 1 - e^(-16).");
        });
    }

    // S1 — parameterized threshold: with positivityThreshold=1, a single detected variant is positive.
    [Test]
    public void DetectMRD_ThresholdOne_SingleDetectedIsPositive()
    {
        var panel = new[] { Marker(4, 200), Marker(0, 200), Marker(0, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel, positivityThreshold: 1);

        Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
            "1 detected >= threshold 1 => positive (parameterized threshold).");
    }

    // S2 — stricter threshold=3: 2 detected is below it => negative.
    [Test]
    public void DetectMRD_ThresholdThree_TwoDetectedIsNegative()
    {
        var panel = new[] { Marker(5, 200), Marker(3, 200), Marker(0, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel, positivityThreshold: 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(2), "Two markers detected.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "2 detected < threshold 3 => negative.");
        });
    }

    // S3 — minSupportingReads=3: a locus with only 2 alt reads is NOT counted as detected.
    [Test]
    public void DetectMRD_MinSupportingReadsThree_LowSupportNotDetected()
    {
        var panel = new[] { Marker(5, 200), Marker(2, 200), Marker(3, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel, minSupportingReads: 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.DetectedVariantCount, Is.EqualTo(2),
                "Only the loci with >= 3 alt reads (5 and 3) count; the 2-read locus is excluded.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "2 detected >= 2 => positive.");
        });
    }

    // C1 — null panel => ArgumentNullException.
    [Test]
    public void DetectMRD_NullPanel_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectMRD(null!),
            "A null marker panel is invalid.");
    }

    // C2 — empty panel => ArgumentException (nothing to interrogate).
    [Test]
    public void DetectMRD_EmptyPanel_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectMRD(Array.Empty<OncologyAnalyzer.TumorMarker>()),
            "An empty panel cannot be assessed for MRD.");
    }

    // C3 — positivityThreshold < 1 => ArgumentOutOfRangeException.
    [Test]
    public void DetectMRD_InvalidThreshold_Throws()
    {
        var panel = new[] { Marker(5, 200), Marker(3, 200) };

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectMRD(panel, positivityThreshold: 0),
            "Positivity threshold must be at least 1.");
    }

    // C5 — minSupportingReads < 1 => ArgumentOutOfRangeException (contract: r_min >= 1).
    [Test]
    public void DetectMRD_InvalidMinSupportingReads_Throws()
    {
        var panel = new[] { Marker(5, 200), Marker(3, 200) };

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectMRD(panel, minSupportingReads: 0),
            "Minimum supporting reads must be at least 1 (MRD_Detection contract §3.1).");
    }

    // C6 — genomeEquivalents < 0 => ArgumentOutOfRangeException (contract: n >= 0).
    [Test]
    public void DetectMRD_NegativeGenomeEquivalents_Throws()
    {
        var panel = new[] { Marker(5, 200), Marker(3, 200) };

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectMRD(panel, genomeEquivalents: -1),
            "Genome equivalents n cannot be negative (MRD_Detection contract §3.1).");
    }

    // M7b — documented edge case: with the default genomeEquivalents = 0, lambda = n*f*m = 0 so the
    // panel Poisson detection probability p = 1 - e^(-0) = 0 (MRD_Detection.md §6.1, §3.1 default n=0).
    [Test]
    public void DetectMRD_DefaultGenomeEquivalents_DetectionProbabilityZero()
    {
        var panel = new[] { Marker(3, 200), Marker(2, 200) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel);

        Assert.That(result.DetectionProbability, Is.EqualTo(0.0).Within(1e-12),
            "n = 0 (default) => lambda = 0 => p = 1 - e^0 = 0 (no informative sampling).");
    }

    // M3b — documented edge case: all total reads = 0 => IMAF = 0 and (with n given) p = 0
    // (MRD_Detection.md §6.1: "All total reads = 0 => IMAF = 0, p = 0").
    [Test]
    public void DetectMRD_AllTotalReadsZero_ImafAndProbabilityZero()
    {
        var panel = new[] { Marker(0, 0), Marker(0, 0), Marker(0, 0) };

        OncologyAnalyzer.MrdResult result = OncologyAnalyzer.DetectMRD(panel, genomeEquivalents: 1000);

        Assert.Multiple(() =>
        {
            Assert.That(result.IntegratedMutantAlleleFraction, Is.EqualTo(0.0).Within(1e-12),
                "Sigma(total) = 0 => IMAF defined as 0 (no informative reads).");
            Assert.That(result.DetectionProbability, Is.EqualTo(0.0).Within(1e-12),
                "f = IMAF = 0 => lambda = 0 => p = 0.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "0 detected < 2 => MRD-negative.");
        });
    }

    #endregion

    #region IsVariantDetected

    // Direct coverage of the per-locus presence call (MRD_Detection.md §2.2: detected iff a_i >= r_min).
    // Boundary at the default cutoff r_min = 1: exactly 1 alt read counts as detected, 0 does not.
    [Test]
    public void IsVariantDetected_DefaultCutoff_OneAltReadIsDetected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsVariantDetected(Marker(1, 200)), Is.True,
                "1 alt read >= default r_min = 1 => detected.");
            Assert.That(OncologyAnalyzer.IsVariantDetected(Marker(0, 200)), Is.False,
                "0 alt reads < r_min = 1 => not detected.");
        });
    }

    // Boundary at a custom cutoff r_min = 3: exactly 3 detected, 2 not detected (Wan 2020 background rule).
    [Test]
    public void IsVariantDetected_CustomCutoff_BoundaryExactlyAtThreshold()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsVariantDetected(Marker(3, 200), minSupportingReads: 3), Is.True,
                "alt = 3 >= r_min = 3 => detected (boundary inclusive).");
            Assert.That(OncologyAnalyzer.IsVariantDetected(Marker(2, 200), minSupportingReads: 3), Is.False,
                "alt = 2 < r_min = 3 => not detected.");
        });
    }

    // minSupportingReads < 1 is invalid for the per-locus call as well.
    [Test]
    public void IsVariantDetected_InvalidMinSupportingReads_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.IsVariantDetected(Marker(1, 200), minSupportingReads: 0),
            "Minimum supporting reads must be at least 1.");
    }

    #endregion

    #region TrackVariantsOverTime

    // M8 — longitudinal: detected counts [0,1,2,3] => status [neg,neg,pos,pos]; first positive at index 2.
    [Test]
    public void TrackVariantsOverTime_RisingSignal_FirstPositiveAtIndexTwo()
    {
        var t0 = new[] { Marker(0, 200), Marker(0, 200), Marker(0, 200) }; // 0 detected
        var t1 = new[] { Marker(4, 200), Marker(0, 200), Marker(0, 200) }; // 1 detected
        var t2 = new[] { Marker(4, 200), Marker(3, 200), Marker(0, 200) }; // 2 detected
        var t3 = new[] { Marker(4, 200), Marker(3, 200), Marker(2, 200) }; // 3 detected
        var timeline = new[] { t0, t1, t2, t3 };

        OncologyAnalyzer.MrdLongitudinalResult result = OncologyAnalyzer.TrackVariantsOverTime(timeline);

        Assert.Multiple(() =>
        {
            Assert.That(result.Timepoints.Count, Is.EqualTo(4), "Four timepoints in order.");
            Assert.That(result.Timepoints[0].Result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "t0: 0 detected => negative.");
            Assert.That(result.Timepoints[1].Result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Negative),
                "t1: 1 detected < 2 => negative.");
            Assert.That(result.Timepoints[2].Result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "t2: 2 detected => positive.");
            Assert.That(result.Timepoints[3].Result.Status, Is.EqualTo(OncologyAnalyzer.MrdStatus.Positive),
                "t3: 3 detected => positive.");
            Assert.That(result.FirstPositiveIndex, Is.EqualTo(2),
                "Earliest MRD-positive timepoint is index 2.");
        });
    }

    // S4 — all-negative timeline => FirstPositiveIndex = -1.
    [Test]
    public void TrackVariantsOverTime_AllNegative_FirstPositiveMinusOne()
    {
        var t0 = new[] { Marker(0, 200), Marker(0, 200) };
        var t1 = new[] { Marker(3, 200), Marker(0, 200) }; // only 1 detected => negative
        var timeline = new[] { t0, t1 };

        OncologyAnalyzer.MrdLongitudinalResult result = OncologyAnalyzer.TrackVariantsOverTime(timeline);

        Assert.Multiple(() =>
        {
            Assert.That(result.FirstPositiveIndex, Is.EqualTo(-1),
                "No timepoint reaches >= 2 detected => no positive => -1.");
            Assert.That(result.Timepoints[1].TimepointIndex, Is.EqualTo(1),
                "Timepoint indices preserve input order.");
        });
    }

    // C4 — null timepoints => ArgumentNullException.
    [Test]
    public void TrackVariantsOverTime_NullTimepoints_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.TrackVariantsOverTime(null!),
            "A null timepoint series is invalid.");
    }

    #endregion
}
