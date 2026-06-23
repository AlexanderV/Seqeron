// ONCO-MRD-001 — Minimal/Molecular Residual Disease Detection
// Evidence: docs/Evidence/ONCO-MRD-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MRD-001.md
// Source: Reinert et al. (2019). JAMA Oncol 5(8):1124-1131 (tumour-informed MRD; ≥2 of tracked SNVs => positive).
//         Natera Signatera white paper (2020): 16 tracked SNVs; "at least two" => ctDNA-positive; p = 1 - e^(-nfm).
//         Wan et al. (2020). Sci Transl Med 12(548):eaaz8084 (INVAR integrated mutant allele fraction / IMAF).
//         Quoted positivity rule: PMC9265001 Table 1.
//         INVAR2 reference impl (nrlab-CRUK/INVAR2): R/shared/detectionFunctions.R (calc_log_likelihood,
//         estimate_p_EM, calc_likelihood_ratio) and R/4_detection/generalisedLikelihoodRatioTest.R
//         (calculateIMAFv2) — exact GLRT / background-subtraction / AF-weighting formulas for EstimateInvarSignal.

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

    #region IntegratedMutantAlleleFractionV2 (INVAR background-subtracted, depth-weighted)

    // Helper: an INVAR locus (alt, total, tumourAF, background).
    private static OncologyAnalyzer.InvarLocus IL(int alt, int total, double af, double bg) =>
        new(alt, total, af, bg);

    // M9 — INVAR2 calculateIMAFv2: per-locus bs = max(0, VAF - background), then depth-weighted mean.
    // Two loci, both above background; hand-derived:
    //  locus1 VAF = 50/1000 = 0.05, bg = 0.01 => bs = 0.04; depth 1000
    //  locus2 VAF = 20/1000 = 0.02, bg = 0.01 => bs = 0.01; depth 1000
    //  IMAFv2 = (0.04*1000 + 0.01*1000) / 2000 = 50/2000 = 0.025
    [Test]
    public void IntegratedMutantAlleleFractionV2_TwoLociAboveBackground_DepthWeightedSubtractedMean()
    {
        var loci = new[] { IL(50, 1000, 0.5, 0.01), IL(20, 1000, 0.5, 0.01) };

        double imafV2 = OncologyAnalyzer.IntegratedMutantAlleleFractionV2(loci);

        Assert.That(imafV2, Is.EqualTo(0.025).Within(1e-12),
            "IMAFv2 = weighted.mean(max(0, VAF - bg), depth) = (0.04*1000 + 0.01*1000)/2000 = 0.025 (INVAR2).");
    }

    // M10 — a locus whose VAF is at/below background contributes 0 (pmax(0, .) clamps the subtraction).
    //  locus1 VAF = 0.05, bg = 0.01 => bs = 0.04 (depth 1000)
    //  locus2 VAF = 0.005, bg = 0.01 => bs = max(0, -0.005) = 0 (depth 1000)
    //  IMAFv2 = (0.04*1000 + 0*1000)/2000 = 40/2000 = 0.02
    [Test]
    public void IntegratedMutantAlleleFractionV2_LocusBelowBackground_ContributesZero()
    {
        var loci = new[] { IL(50, 1000, 0.5, 0.01), IL(5, 1000, 0.5, 0.01) };

        double imafV2 = OncologyAnalyzer.IntegratedMutantAlleleFractionV2(loci);

        Assert.That(imafV2, Is.EqualTo(0.02).Within(1e-12),
            "Below-background locus subtracts to 0 (pmax) => IMAFv2 = 0.04*1000/2000 = 0.02 (INVAR2).");
    }

    // M11 — pure background everywhere => every bs = 0 => IMAFv2 = 0 (background subtraction removes noise).
    [Test]
    public void IntegratedMutantAlleleFractionV2_PureBackground_Zero()
    {
        // Each locus has VAF exactly equal to background (1/1000 = 0.001 == bg).
        var loci = new[] { IL(1, 1000, 0.5, 0.001), IL(1, 1000, 0.5, 0.001), IL(1, 1000, 0.5, 0.001) };

        double imafV2 = OncologyAnalyzer.IntegratedMutantAlleleFractionV2(loci);

        Assert.That(imafV2, Is.EqualTo(0.0).Within(1e-12),
            "VAF == background at every locus => max(0, VAF-bg) = 0 => IMAFv2 = 0.");
    }

    // C7 — zero-coverage loci contribute no weight; all-zero-coverage => IMAFv2 = 0.
    [Test]
    public void IntegratedMutantAlleleFractionV2_NoCoverage_Zero()
    {
        var loci = new[] { IL(0, 0, 0.5, 0.01) };

        double imafV2 = OncologyAnalyzer.IntegratedMutantAlleleFractionV2(loci);

        Assert.That(imafV2, Is.EqualTo(0.0).Within(1e-12),
            "Total depth 0 => no weight => IMAFv2 = 0.");
    }

    // C8 — null loci => ArgumentNullException.
    [Test]
    public void IntegratedMutantAlleleFractionV2_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.IntegratedMutantAlleleFractionV2(null!),
            "Null loci collection is invalid.");
    }

    #endregion

    #region EstimateInvarSignal (INVAR GLRT: background subtraction + AF weighting)

    // Builds n identical loci with the given per-locus mutant-read count, depth, tumour AF and background.
    private static OncologyAnalyzer.InvarLocus[] UniformLoci(int n, int alt, int total, double af, double bg)
    {
        var loci = new OncologyAnalyzer.InvarLocus[n];
        for (int i = 0; i < n; i++)
        {
            loci[i] = IL(alt, total, af, bg);
        }

        return loci;
    }

    // M12 — Background subtraction removes pure noise: mutant reads only at the background rate
    // (M = 1 per 1000 reads, e = 0.001) => EM p-hat ~ 0 and LR ~ 0 => NOT detected.
    // Reference (INVAR2 formulas, Evidence dataset inj=0): p-hat ~ 3.3e-5, LR ~ -0.0001.
    [Test]
    public void EstimateInvarSignal_PureBackground_NotDetectedAndZeroFraction()
    {
        var loci = UniformLoci(n: 50, alt: 1, total: 1000, af: 0.4, bg: 0.001);

        OncologyAnalyzer.InvarSignalResult result = OncologyAnalyzer.EstimateInvarSignal(loci);

        Assert.Multiple(() =>
        {
            Assert.That(result.EstimatedTumorFraction, Is.EqualTo(0.0).Within(1e-3),
                "Pure-background reads => EM ctDNA fraction p-hat ~ 0 (INVAR2 estimate_p_EM).");
            Assert.That(result.LikelihoodRatio, Is.EqualTo(0.0).Within(1e-2),
                "Pure background => GLRT LR ~ 0 (no ctDNA evidence above background).");
            Assert.That(result.Detected, Is.False,
                "LR ~ 0 (and barely any signal) => not detected with default threshold 0 only via near-zero LR.");
        });
    }

    // M13 — Injected signal is recovered: tumour AF 0.4, background 0.001, depth 1000, injected p = 0.01
    // => expected mutant rate q = 0.01*g + 0.001*0.99 with g = 0.4*0.999+0.6*0.001 => M = round(q*1000) = 5.
    // Reference (Evidence dataset inj=0.01): p-hat ~ 0.01002, LR ~ 4.06.
    [Test]
    public void EstimateInvarSignal_InjectedOnePercent_RecoversFractionAndDetects()
    {
        var loci = UniformLoci(n: 50, alt: 5, total: 1000, af: 0.4, bg: 0.001);

        OncologyAnalyzer.InvarSignalResult result = OncologyAnalyzer.EstimateInvarSignal(loci);

        Assert.Multiple(() =>
        {
            Assert.That(result.EstimatedTumorFraction, Is.EqualTo(0.01002).Within(5e-4),
                "EM recovers the injected ctDNA fraction ~ 0.010 (INVAR2 estimate_p_EM reference).");
            Assert.That(result.LikelihoodRatio, Is.EqualTo(4.06).Within(0.05),
                "GLRT statistic ~ 4.06 for injected p = 0.01 (INVAR2 calc_likelihood_ratio reference).");
            Assert.That(result.Detected, Is.True,
                "LR > 0 with mutant reads present => detected at default threshold.");
            Assert.That(result.LocusCount, Is.EqualTo(50), "All 50 loci are informative (AF > 0).");
        });
    }

    // M14 — Recovery at higher injection: injected p = 0.05 (M = 21) => p-hat ~ 0.0501, LR ~ 44.14.
    [Test]
    public void EstimateInvarSignal_InjectedFivePercent_RecoversFraction()
    {
        var loci = UniformLoci(n: 50, alt: 21, total: 1000, af: 0.4, bg: 0.001);

        OncologyAnalyzer.InvarSignalResult result = OncologyAnalyzer.EstimateInvarSignal(loci);

        Assert.Multiple(() =>
        {
            Assert.That(result.EstimatedTumorFraction, Is.EqualTo(0.0501).Within(1e-3),
                "EM recovers injected ctDNA fraction ~ 0.050 (INVAR2 reference).");
            Assert.That(result.LikelihoodRatio, Is.EqualTo(44.14).Within(0.3),
                "GLRT statistic ~ 44.1 for injected p = 0.05 (INVAR2 reference).");
            Assert.That(result.Detected, Is.True, "Strong signal => detected.");
        });
    }

    // M15 — Monotonicity: more injected signal => strictly larger GLRT statistic (INVAR2 reference table).
    // M/locus for inj {0, 0.005, 0.01, 0.02, 0.05} at AF 0.4, e 0.001, depth 1000 => {1, 5, 5, 9, 21}.
    [Test]
    public void EstimateInvarSignal_RisingSignal_LikelihoodRatioMonotoneIncreasing()
    {
        int[] mutantReads = { 1, 5, 5, 9, 21 };
        var lrs = new double[mutantReads.Length];
        for (int i = 0; i < mutantReads.Length; i++)
        {
            var loci = UniformLoci(n: 50, alt: mutantReads[i], total: 1000, af: 0.4, bg: 0.001);
            lrs[i] = OncologyAnalyzer.EstimateInvarSignal(loci).LikelihoodRatio;
        }

        Assert.Multiple(() =>
        {
            for (int i = 1; i < lrs.Length; i++)
            {
                Assert.That(lrs[i], Is.GreaterThanOrEqualTo(lrs[i - 1] - 1e-9),
                    $"LR must be non-decreasing as injected signal rises (step {i}): {lrs[i - 1]} -> {lrs[i]}.");
            }

            // Pinned reference endpoints: pure background ~ 0, strongest ~ 44 (INVAR2 dataset).
            Assert.That(lrs[0], Is.EqualTo(0.0).Within(1e-2), "Pure background LR ~ 0.");
            Assert.That(lrs[^1], Is.EqualTo(44.14).Within(0.3), "Strongest signal LR ~ 44.1.");
        });
    }

    // M16 — AF weighting boosts sensitivity vs flat pooling on a low-signal mixture.
    // 20 high-AF loci (0.5) + 20 low-AF loci (0.05), depth 2000, e 0.002, injected p = 0.008.
    // Mutant reads per locus from its TRUE AF: high-AF M = round(q_high*2000), low-AF M = round(q_low*2000).
    // Weighted model uses true per-locus AF; "unweighted" replaces every AF by the panel mean AF (0.275).
    // Reference (Evidence AF-weighting dataset): weighted LR ~ 2.66 > unweighted LR ~ 1.91.
    [Test]
    public void EstimateInvarSignal_AfWeighting_HigherLikelihoodRatioThanFlatPooling()
    {
        const int depth = 2000;
        const double e = 0.002;
        const double injected = 0.008;
        const double highAf = 0.5;
        const double lowAf = 0.05;

        double G(double af) => (af * (1 - e)) + ((1 - af) * e);
        int MutantReads(double af)
        {
            double q = (injected * G(af)) + (e * (1 - injected));
            return (int)Math.Round(q * depth);
        }

        int mHigh = MutantReads(highAf);
        int mLow = MutantReads(lowAf);
        double meanAf = ((20 * highAf) + (20 * lowAf)) / 40.0; // 0.275

        var weighted = new List<OncologyAnalyzer.InvarLocus>();
        var flat = new List<OncologyAnalyzer.InvarLocus>();
        for (int i = 0; i < 20; i++)
        {
            weighted.Add(IL(mHigh, depth, highAf, e));
            flat.Add(IL(mHigh, depth, meanAf, e));
        }

        for (int i = 0; i < 20; i++)
        {
            weighted.Add(IL(mLow, depth, lowAf, e));
            flat.Add(IL(mLow, depth, meanAf, e));
        }

        double weightedLr = OncologyAnalyzer.EstimateInvarSignal(weighted).LikelihoodRatio;
        double flatLr = OncologyAnalyzer.EstimateInvarSignal(flat).LikelihoodRatio;

        Assert.Multiple(() =>
        {
            Assert.That(weightedLr, Is.GreaterThan(flatLr),
                "Per-locus AF weighting concentrates signal at high-SNR loci => larger GLRT than flat pooling.");
            Assert.That(weightedLr, Is.EqualTo(2.66).Within(0.1),
                "Weighted GLRT ~ 2.66 (INVAR2 AF-weighting reference).");
            Assert.That(flatLr, Is.EqualTo(1.91).Within(0.1),
                "Flat-pooled GLRT ~ 1.91 (INVAR2 AF-weighting reference).");
        });
    }

    // S5 — detectionThreshold gates the call: a weak signal whose LR is below the threshold is NOT detected.
    [Test]
    public void EstimateInvarSignal_HighDetectionThreshold_WeakSignalNotDetected()
    {
        // Injected p = 0.005 => LR ~ 1.30 (Evidence dataset). A threshold of 5 is not reached.
        var loci = UniformLoci(n: 50, alt: 5, total: 1000, af: 0.4, bg: 0.001);

        OncologyAnalyzer.InvarSignalResult low = OncologyAnalyzer.EstimateInvarSignal(loci, detectionThreshold: 0.0);
        OncologyAnalyzer.InvarSignalResult high = OncologyAnalyzer.EstimateInvarSignal(loci, detectionThreshold: 5.0);

        Assert.Multiple(() =>
        {
            Assert.That(low.Detected, Is.True, "At threshold 0 the positive LR is detected.");
            Assert.That(high.Detected, Is.False,
                "LR ~ 4.06 < threshold 5 => not detected (specificity knob).");
        });
    }

    // S6 — zero background is floored to 1/depth so the likelihood stays finite (INVAR2 doMain guard);
    // a clear signal with e = 0 is still recovered and detected.
    [Test]
    public void EstimateInvarSignal_ZeroBackground_FiniteAndDetects()
    {
        var loci = UniformLoci(n: 50, alt: 10, total: 1000, af: 0.4, bg: 0.0);

        OncologyAnalyzer.InvarSignalResult result = OncologyAnalyzer.EstimateInvarSignal(loci);

        Assert.Multiple(() =>
        {
            Assert.That(double.IsFinite(result.LikelihoodRatio), Is.True,
                "Zero background is floored to 1/depth => log-likelihood finite.");
            Assert.That(result.EstimatedTumorFraction, Is.GreaterThan(0.0),
                "Clear mutant signal => positive ctDNA fraction estimate.");
            Assert.That(result.Detected, Is.True, "Clear signal => detected.");
        });
    }

    // C9 — null loci => ArgumentNullException.
    [Test]
    public void EstimateInvarSignal_NullLoci_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.EstimateInvarSignal(null!),
            "Null loci collection is invalid.");
    }

    // C10 — no informative locus (every tumour AF = 0) => ArgumentException.
    [Test]
    public void EstimateInvarSignal_NoInformativeLocus_Throws()
    {
        var loci = new[] { IL(5, 1000, 0.0, 0.001), IL(3, 1000, 0.0, 0.001) };

        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.EstimateInvarSignal(loci),
            "All tumour AF = 0 => no informative locus to estimate signal.");
    }

    // C11 — negative detection threshold => ArgumentOutOfRangeException.
    [Test]
    public void EstimateInvarSignal_NegativeThreshold_Throws()
    {
        var loci = UniformLoci(n: 5, alt: 5, total: 1000, af: 0.4, bg: 0.001);

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.EstimateInvarSignal(loci, detectionThreshold: -1.0),
            "Detection threshold cannot be negative.");
    }

    // C12 — out-of-range tumour AF (> 1) => ArgumentOutOfRangeException.
    [Test]
    public void EstimateInvarSignal_TumourAfAboveOne_Throws()
    {
        var loci = new[] { IL(5, 1000, 1.5, 0.001) };

        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateInvarSignal(loci),
            "Tumour allele fraction must be in [0, 1].");
    }

    // C13 — out-of-range background rate (>= 1) => ArgumentOutOfRangeException.
    [Test]
    public void EstimateInvarSignal_BackgroundRateAtOne_Throws()
    {
        var loci = new[] { IL(5, 1000, 0.4, 1.0) };

        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateInvarSignal(loci),
            "Background error rate must be in [0, 1).");
    }

    #endregion
}
