// ONCO-CNA-002 — Focal Amplification Detection
// Evidence: docs/Evidence/ONCO-CNA-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CNA-002.md
// Source: Mermel CH et al. (2011). GISTIC2.0. Genome Biology 12:R41.
//         https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
//         GISTIC2 docs broad_len_cutoff=0.98, t_amp=0.1; NCBI Gene oncogene arms.

using System.Globalization;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberArmSegment;
using Thresholds = Seqeron.Genomics.Oncology.OncologyAnalyzer.FocalAmplificationThresholds;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_DetectFocalAmplifications_Tests
{
    // Arm of fixed length 1,000,000 bp so that segment-length fraction is read directly.
    private const long Arm = 1_000_000;

    private static Segment Seg(string arm, long start, long end, double log2) =>
        new(arm, start, end, Arm, log2);

    #region DetectFocalAmplifications

    // M1 — 17q, length 0.50 of arm (< 0.98) and log2 1.0 (> t_amp 0.1) ⇒ focal amplification.
    // Source: Mermel 2011 (<98% ⇒ focal); GISTIC2 t_amp=0.1.
    [Test]
    public void DetectFocalAmplifications_FocalHighAmp_Reported()
    {
        var segments = new[] { Seg("17q", 100_000, 600_000, 1.0) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Has.Count.EqualTo(1),
            "A 0.50-of-arm, log2=1.0 segment is amplified (>0.1) and focal (<0.98), so it must be reported.");
        Assert.That(result[0].Arm, Is.EqualTo("17q"),
            "The reported segment must be the input 17q segment.");
    }

    // M2 — 8q occupying 0.99 of arm (> 0.98) is arm-level even at log2 1.5, so NOT focal.
    // Source: Mermel 2011 — events occupying >98% of an arm are arm-level.
    [Test]
    public void DetectFocalAmplifications_WholeArm_NotReported()
    {
        var segments = new[] { Seg("8q", 0, 990_000, 1.5) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Is.Empty,
            "A 0.99-of-arm segment is arm-level (>0.98), so it must NOT be reported as focal even when highly amplified.");
    }

    // M3 — exactly 0.98 of arm is the cutoff; focal test is strict < 0.98 ⇒ NOT focal.
    // Source: GISTIC2 broad_len_cutoff=0.98; paper "more than 98% ⇒ arm-level".
    [Test]
    public void DetectFocalAmplifications_ExactlyCutoff_NotReported()
    {
        var segments = new[] { Seg("11q", 0, 980_000, 1.0) }; // 980,000/1,000,000 = 0.98 exactly

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Is.Empty,
            "A segment whose length equals exactly 0.98 of the arm is arm-level (focal test is strictly < 0.98).");
    }

    // M4 — log2 0.05 does not exceed t_amp 0.1 ⇒ not amplified, excluded though focal in length.
    // Source: GISTIC2 t_amp=0.1.
    [Test]
    public void DetectFocalAmplifications_LowAmplitude_NotReported()
    {
        var segments = new[] { Seg("7p", 0, 300_000, 0.05) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Is.Empty,
            "log2=0.05 is below t_amp=0.1, so the segment is not amplified and must be excluded despite being focal.");
    }

    // M5 — single-copy gain log2(3/2)=0.585 (> 0.1) at 0.10-of-arm ⇒ focal amplification.
    // Source: CNVkit log2(3/2)=0.585; GISTIC2 t_amp=0.1; Mermel 2011 focal.
    [Test]
    public void DetectFocalAmplifications_JustAboveAmpCutoff_Reported()
    {
        const double SingleCopyGainLog2 = 0.585; // log2(3/2), CNVkit
        var segments = new[] { Seg("12q", 0, 100_000, SingleCopyGainLog2) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Has.Count.EqualTo(1),
            "A single-copy-gain log2 (0.585) exceeds t_amp=0.1 and 0.10-of-arm is focal, so it must be reported.");
    }

    // M4b — log2 exactly at t_amp (0.1) is NOT amplified: GISTIC2 t_amp is "above this positive
    // value", so the amplitude test is strictly greater-than (boundary excluded).
    // Source: GISTIC2 docs t_amp — "Regions with a copy number gain ABOVE this positive value are
    // considered amplified." (https://broadinstitute.github.io/gistic2/)
    [Test]
    public void DetectFocalAmplifications_Log2ExactlyAtTamp_NotReported()
    {
        // 0.10-of-arm (focal) but log2 == t_amp exactly: amplitude test is strict > 0.1 ⇒ not amplified.
        var segments = new[] { Seg("17q", 0, 100_000, OncologyAnalyzer.DefaultAmplificationLog2Threshold) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.That(result, Is.Empty,
            "log2 exactly equal to t_amp (0.1) is not 'above' the threshold, so the segment is not amplified.");
    }

    // M11 — mixed list: only the two focal amplifications survive, in input order.
    // Source: INV-03 (subset, order-preserving).
    [Test]
    public void DetectFocalAmplifications_MixedList_PreservesFocalSubsetInOrder()
    {
        var segments = new[]
        {
            Seg("17q", 0, 500_000, 1.0),   // focal amp  -> keep (index 0)
            Seg("8q", 0, 990_000, 1.5),    // arm-level  -> drop
            Seg("7p", 0, 300_000, 0.05),   // low amp    -> drop
            Seg("11q", 0, 200_000, 0.7),   // focal amp  -> keep (index 3)
        };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(2),
                "Only the two focal amplifications must remain (arm-level and low-amp are filtered out).");
            Assert.That(result[0].Arm, Is.EqualTo("17q"),
                "First kept segment must be the 17q focal amplification (input order preserved).");
            Assert.That(result[1].Arm, Is.EqualTo("11q"),
                "Second kept segment must be the 11q focal amplification (input order preserved).");
        });
    }

    // S1 — custom t_amp=0.3: a log2=0.2 segment is now below threshold ⇒ not reported.
    // Source: GISTIC2 t_amp is a parameter; override changes the amplitude gate.
    [Test]
    public void DetectFocalAmplifications_CustomTamp_BelowThreshold_NotReported()
    {
        var thresholds = new Thresholds(0.3, OncologyAnalyzer.DefaultBroadLengthCutoff);
        var segments = new[] { Seg("17q", 0, 200_000, 0.2) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds);

        Assert.That(result, Is.Empty,
            "With t_amp raised to 0.3, a log2=0.2 segment is below threshold and must not be reported.");
    }

    // S1b — custom broad_len_cutoff: raising it to 0.999 makes a 0.99-of-arm segment focal again.
    // Source: GISTIC2 broad_len_cutoff is a parameter (fraction of arm); the focal test is L/A < cutoff.
    [Test]
    public void DetectFocalAmplifications_CustomBroadLengthCutoff_AdmitsLongerSegment()
    {
        var thresholds = new Thresholds(OncologyAnalyzer.DefaultAmplificationLog2Threshold, 0.999);
        var segments = new[] { Seg("8q", 0, 990_000, 1.5) }; // 0.99 of arm; arm-level under default 0.98

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds);

        Assert.That(result, Has.Count.EqualTo(1),
            "With broad_len_cutoff raised to 0.999, a 0.99-of-arm amplified segment is focal (0.99 < 0.999).");
    }

    // IsFocalAmplification (public predicate) — direct coverage of the conjunction it computes.
    // Source: Mermel 2011 focal (L/A < 0.98) AND GISTIC2 t_amp (log2 > 0.1).
    [Test]
    public void IsFocalAmplification_Predicate_FocalAndAmplified_True()
    {
        var seg = Seg("17q", 100_000, 600_000, 1.0); // 0.50 of arm, log2 1.0

        Assert.That(OncologyAnalyzer.IsFocalAmplification(seg, Thresholds.Default), Is.True,
            "A 0.50-of-arm (focal), log2=1.0 (amplified) segment satisfies both predicates.");
    }

    [Test]
    public void IsFocalAmplification_Predicate_ArmLevel_False()
    {
        var seg = Seg("8q", 0, 990_000, 1.5); // 0.99 of arm ⇒ not focal

        Assert.That(OncologyAnalyzer.IsFocalAmplification(seg, Thresholds.Default), Is.False,
            "A 0.99-of-arm segment is arm-level (not focal), so the predicate is false even when amplified.");
    }

    [Test]
    public void IsFocalAmplification_Predicate_NotAmplified_False()
    {
        var seg = Seg("7p", 0, 300_000, 0.05); // focal length but log2 0.05 < t_amp

        Assert.That(OncologyAnalyzer.IsFocalAmplification(seg, Thresholds.Default), Is.False,
            "A focal-length but low-amplitude (log2=0.05) segment is not amplified, so the predicate is false.");
    }

    [Test]
    public void IsFocalAmplification_Predicate_InvalidSegment_Throws()
    {
        var seg = new Segment("17q", 0, 100, 0, 1.0); // non-positive arm length

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.IsFocalAmplification(seg, Thresholds.Default),
            "A segment with non-positive arm length must throw ArgumentException from the predicate.");
    }

    // C1 — null segments ⇒ ArgumentNullException.
    [Test]
    public void DetectFocalAmplifications_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectFocalAmplifications(null!),
            "Null segments must throw ArgumentNullException.");
    }

    // C2 — empty input ⇒ empty output.
    [Test]
    public void DetectFocalAmplifications_Empty_ReturnsEmpty()
    {
        var result = OncologyAnalyzer.DetectFocalAmplifications(Array.Empty<Segment>());

        Assert.That(result, Is.Empty, "Empty input must yield an empty result.");
    }

    // C4 — non-positive arm length ⇒ ArgumentException.
    [Test]
    public void DetectFocalAmplifications_NonPositiveArmLength_Throws()
    {
        var segments = new[] { new Segment("17q", 0, 100, 0, 1.0) };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectFocalAmplifications(segments),
            "A segment with non-positive arm length must throw ArgumentException.");
    }

    // C4b — End <= Start ⇒ ArgumentException.
    [Test]
    public void DetectFocalAmplifications_EndNotAfterStart_Throws()
    {
        var segments = new[] { new Segment("17q", 500, 500, Arm, 1.0) };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectFocalAmplifications(segments),
            "A segment with End <= Start must throw ArgumentException.");
    }

    // C5 — thresholds outside the GISTIC2 reference-implementation ranges are rejected.
    // Source: GISTIC2 source/gp_gistic2_from_seg.m — numeric_arg(a,'ta',0.1,[0,Inf]),
    // numeric_arg(a,'brlen',0.98,[0 2]); non-numeric (NaN) values throw 'snp:badarg:nonnumeric'.
    [TestCase(-0.1, 0.98, TestName = "DetectFocalAmplifications_NegativeTamp_Throws")]
    [TestCase(double.NaN, 0.98, TestName = "DetectFocalAmplifications_NaNTamp_Throws")]
    [TestCase(0.1, double.NaN, TestName = "DetectFocalAmplifications_NaNBroadLenCutoff_Throws")]
    [TestCase(0.1, -0.01, TestName = "DetectFocalAmplifications_NegativeBroadLenCutoff_Throws")]
    [TestCase(0.1, 2.01, TestName = "DetectFocalAmplifications_BroadLenCutoffAbove2_Throws")]
    public void DetectFocalAmplifications_ThresholdsOutsideGistic2Range_Throw(double tAmp, double cutoff)
    {
        var thresholds = new Thresholds(tAmp, cutoff);
        var segments = new[] { Seg("17q", 0, 100_000, 1.0) };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds),
                "Thresholds outside GISTIC2 ranges (t_amp in [0,Inf], broad_len_cutoff in [0,2]) must be rejected.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.IsFocalAmplification(segments[0], thresholds),
                "The single-segment predicate must reject the same out-of-range thresholds.");
            // Validation happens even for empty input (no silent acceptance of a bad parameter).
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DetectFocalAmplifications(Array.Empty<Segment>(), thresholds),
                "Out-of-range thresholds must be rejected even when there are no segments.");
        });
    }

    // C5b — range endpoints are allowed: t_amp = 0 (any positive gain amplified), broad_len_cutoff = 2
    // (upper bound of GISTIC2 -brlen). A 0.99-of-arm log2 0.05 segment is then a focal amplification.
    // Source: GISTIC2 gp_gistic2_from_seg.m ranges [0,Inf] and [0 2] are inclusive (val < lo || val > hi rejects).
    [Test]
    public void DetectFocalAmplifications_Gistic2RangeEndpoints_Accepted()
    {
        var thresholds = new Thresholds(0.0, 2.0);
        var segments = new[] { Seg("8q", 0, 990_000, 0.05) };

        var result = OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds);

        Assert.That(result, Has.Count.EqualTo(1),
            "With t_amp=0 and broad_len_cutoff=2, log2 0.05 > 0 and 0.99 < 2, so the segment is reported.");
    }

    // C6 — a NaN log2 ratio is a no-call: NaN is not "above" t_amp, so the segment is not amplified.
    // Source: GISTIC2 amplitude test (value compared against t_amp; NaN comparisons are false).
    [Test]
    public void DetectFocalAmplifications_NaNLog2_NotReported()
    {
        var segments = new[] { Seg("17q", 0, 100_000, double.NaN) };

        Assert.That(OncologyAnalyzer.DetectFocalAmplifications(segments), Is.Empty,
            "A NaN (no-call) log2 ratio is not above t_amp, so the segment must not be reported.");
    }

    #endregion

    #region IdentifyAmplifiedOncogenes

    // M6 — 17q focal amp ⇒ ERBB2. Source: NCBI Gene 2064 (17q12).
    [Test]
    public void IdentifyAmplifiedOncogenes_Arm17q_ReturnsErbb2()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("17q", 0, 100_000, 1.0) });

        Assert.That(genes, Is.EquivalentTo(new[] { "ERBB2" }),
            "A focal amplification on 17q maps to ERBB2 (NCBI Gene 17q12).");
    }

    // M7 — 8q focal amp ⇒ MYC. Source: NCBI Gene 4609 (8q24.21).
    [Test]
    public void IdentifyAmplifiedOncogenes_Arm8q_ReturnsMyc()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("8q", 0, 100_000, 1.0) });

        Assert.That(genes, Is.EquivalentTo(new[] { "MYC" }),
            "A focal amplification on 8q maps to MYC (NCBI Gene 8q24.21).");
    }

    // M8 — 7p focal amp ⇒ EGFR. Source: NCBI Gene 1956 (7p11.2).
    [Test]
    public void IdentifyAmplifiedOncogenes_Arm7p_ReturnsEgfr()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("7p", 0, 100_000, 1.0) });

        Assert.That(genes, Is.EquivalentTo(new[] { "EGFR" }),
            "A focal amplification on 7p maps to EGFR (NCBI Gene 7p11.2).");
    }

    // M9 — 11q focal amp ⇒ CCND1. Source: NCBI Gene 595 (11q13.3).
    [Test]
    public void IdentifyAmplifiedOncogenes_Arm11q_ReturnsCcnd1()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("11q", 0, 100_000, 1.0) });

        Assert.That(genes, Is.EquivalentTo(new[] { "CCND1" }),
            "A focal amplification on 11q maps to CCND1 (NCBI Gene 11q13.3).");
    }

    // M10 — 12q focal amp ⇒ both MDM2 and CDK4. Source: NCBI Gene 4193 (12q15), 1019 (12q14.1).
    [Test]
    public void IdentifyAmplifiedOncogenes_Arm12q_ReturnsMdm2AndCdk4()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("12q", 0, 100_000, 1.0) });

        Assert.That(genes, Is.EquivalentTo(new[] { "MDM2", "CDK4" }),
            "A focal amplification on 12q maps to both MDM2 (12q15) and CDK4 (12q14.1) per NCBI Gene.");
    }

    // M12 — only focal amplifications feed the mapper: a low-amp segment (filtered) yields no genes.
    // Source: INV-04.
    [Test]
    public void IdentifyAmplifiedOncogenes_NonAmplifiedArm_NotMapped()
    {
        // DetectFocalAmplifications drops the low-amp 17q segment, so its arm is never mapped.
        var focal = OncologyAnalyzer.DetectFocalAmplifications(new[] { Seg("17q", 0, 100_000, 0.05) });
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(focal);

        Assert.That(genes, Is.Empty,
            "A non-amplified 17q segment is filtered out, so ERBB2 must NOT be reported.");
    }

    // S2 — arm with no panel oncogene (5q) ⇒ empty.
    [Test]
    public void IdentifyAmplifiedOncogenes_ArmWithoutPanelGene_Empty()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { Seg("5q", 0, 100_000, 1.0) });

        Assert.That(genes, Is.Empty, "An amplification on 5q maps to no panel oncogene.");
    }

    // INV-4 — several amplifications on the same/mixed-case arms report each gene once, in panel order.
    // Source: NCBI Gene loci (ERBB2 17q12, CCND1 11q13.3); arm labels matched case-insensitively.
    [Test]
    public void IdentifyAmplifiedOncogenes_DuplicateAndMixedCaseArms_DistinctInPanelOrder()
    {
        var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[]
        {
            Seg("11q", 0, 100_000, 1.0),
            Seg("17Q", 0, 100_000, 1.0),
            Seg("17q", 200_000, 300_000, 2.0),
        });

        Assert.That(genes, Is.EqualTo(new[] { "ERBB2", "CCND1" }),
            "ERBB2 (17q) and CCND1 (11q) must each be reported once, in panel order (ERBB2 before CCND1).");
    }

    // C3 — null amplifications ⇒ ArgumentNullException.
    [Test]
    public void IdentifyAmplifiedOncogenes_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.IdentifyAmplifiedOncogenes(null!),
            "Null amplifications must throw ArgumentNullException.");
    }

    #endregion

    #region Marker-unit arm fraction (F48, GISTIC2 normalize_by_arm_length norm_type = 1)

    // Values locked from GNU Octave running the original GISTIC2 normalize_by_arm_length.m (broadinstitute/gistic2
    // master 26c590bd) on a toy marker map (chr1 p = 10 markers, q = 40 markers; norm_type = 1, ref_length = 2):
    // 3/10 → 0.29999999999999999, 39/40 → 0.97499999999999998, 40/40 → 1, 10/10 → 1.
    private static Segment MSeg(string arm, long start, long end, double log2, int markers, int armMarkers) =>
        new(arm, start, end, Arm, log2) { MarkerCount = markers, ArmMarkerCount = armMarkers };

    [TestCase(3, 10, 0.29999999999999999)]
    [TestCase(39, 40, 0.97499999999999998)]
    [TestCase(40, 40, 1.0)]
    [TestCase(10, 10, 1.0)]
    public void ArmFraction_MarkerCounts_MatchesGistic2Octave(int markers, int armMarkers, double expected)
    {
        var seg = MSeg("1q", 0, 500_000, 1.0, markers, armMarkers);
        Assert.That(seg.ArmFraction, Is.EqualTo(expected),
            "With marker counts the arm fraction is markers ÷ arm markers (GISTIC2 default norm_type = 1).");
    }

    // 0.99 of the arm in bp (broad) but 39/40 = 0.975 markers (< 0.98 ⇒ focal under GISTIC2's default units).
    [Test]
    public void DetectFocalAmplifications_MarkerUnits_FocalWhereBpIsBroad()
    {
        var bp = Seg("17q", 0, 990_000, 1.0);
        var markers = MSeg("17q", 0, 990_000, 1.0, 39, 40);

        Assert.That(OncologyAnalyzer.DetectFocalAmplifications(new[] { bp }), Is.Empty,
            "Without marker counts the bp fraction 0.99 ≥ 0.98 is broad.");
        Assert.That(OncologyAnalyzer.DetectFocalAmplifications(new[] { markers }), Has.Count.EqualTo(1),
            "With marker counts 39/40 = 0.975 < 0.98 is focal (GISTIC2 marker units).");
    }

    // 0.5 of the arm in bp (focal) but all 40/40 arm markers (fraction 1 ≥ 0.98 ⇒ broad).
    [Test]
    public void DetectFocalAmplifications_MarkerUnits_BroadWhereBpIsFocal()
    {
        var markers = MSeg("8q", 0, 500_000, 1.0, 40, 40);
        Assert.That(OncologyAnalyzer.IsFocalAmplification(markers, Thresholds.Default), Is.False,
            "40/40 markers = whole arm in GISTIC2 units ⇒ arm-level, not focal.");
    }

    [Test]
    public void ArmFraction_NoMarkerCounts_StaysBp()
    {
        Assert.That(Seg("8q", 0, 250_000, 1.0).ArmFraction, Is.EqualTo(0.25),
            "Absent marker counts the bp fraction is unchanged.");
    }

    [TestCase(5, null)]
    [TestCase(null, 5)]
    [TestCase(0, 5)]
    [TestCase(6, 5)]
    public void IsFocalAmplification_InvalidMarkerCounts_Throws(int? markers, int? armMarkers)
    {
        var seg = new Segment("8q", 0, 500_000, Arm, 1.0) { MarkerCount = markers, ArmMarkerCount = armMarkers };
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IsFocalAmplification(seg, Thresholds.Default),
            "Marker counts must be both-or-neither with 1 ≤ MarkerCount ≤ ArmMarkerCount.");
    }

    #endregion

    #region IdentifyAmplifiedOncogenes — locus overlap (F49, GISTIC2 genes_at partial_hits = 1)

    // Expected gene lists locked from GNU Octave running the original GISTIC2 genes_at.m (partial_hits default)
    // on the 12-gene GRCh38 panel (GENCODE v22 coordinates from GISTIC2 refgenes/Gencode.v22.170324).
    [TestCase("17", 39_730_426, 39_800_000, new[] { "ERBB2" })]          // touches ERBB2 end
    [TestCase("17", 39_730_427, 39_800_000, new string[0])]             // one base past ERBB2 end
    [TestCase("chr17", 39_600_000, 39_687_914, new[] { "ERBB2" })]      // touches ERBB2 start
    [TestCase("17", 39_600_000, 39_687_913, new string[0])]
    [TestCase("12", 57_756_013, 68_808_172, new[] { "MDM2", "CDK4" })]  // panel order
    [TestCase("12", 57_756_014, 68_808_171, new string[0])]
    [TestCase("8", 127_736_000, 127_740_000, new[] { "MYC" })]          // inside the gene
    public void IdentifyAmplifiedOncogenes_LocusOverlap_MatchesGistic2GenesAt(
        string chromosome, long start, long end, string[] expected)
    {
        var result = OncologyAnalyzer.IdentifyAmplifiedOncogenes(
            new[] { new OncologyAnalyzer.CopyNumberRegion(chromosome, start, end) });
        Assert.That(result, Is.EqualTo(expected),
            "Genes are reported iff gene.start ≤ region.end ∧ gene.end ≥ region.start (closed intervals).");
    }

    // Arm-level would report ERBB2 for any 17q amplification; locus overlap does not.
    [Test]
    public void IdentifyAmplifiedOncogenes_LocusOverlap_FocalElsewhereOnArm_NotErbb2()
    {
        var region = new OncologyAnalyzer.CopyNumberRegion("17", 60_000_000, 61_000_000);
        Assert.That(OncologyAnalyzer.IdentifyAmplifiedOncogenes(new[] { region }), Is.Empty,
            "A 17q amplification away from the ERBB2 locus must not report ERBB2.");
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_LocusOverlap_CustomPanelAndDistinct()
    {
        var panel = new[]
        {
            new OncologyAnalyzer.GeneLocus("G2", "X", 200, 300),
            new OncologyAnalyzer.GeneLocus("G1", "X", 100, 150),
        };
        var regions = new[]
        {
            new OncologyAnalyzer.CopyNumberRegion("chrX", 150, 150),
            new OncologyAnalyzer.CopyNumberRegion("x", 120, 250),
        };
        Assert.That(OncologyAnalyzer.IdentifyAmplifiedOncogenes(regions, panel), Is.EqualTo(new[] { "G2", "G1" }),
            "Each gene once, in panel order; 'chrX'/'x' match 'X'.");
    }

    [Test]
    public void IdentifyAmplifiedOncogenes_LocusOverlap_InvalidInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            OncologyAnalyzer.IdentifyAmplifiedOncogenes((IEnumerable<OncologyAnalyzer.CopyNumberRegion>)null!));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IdentifyAmplifiedOncogenes(
            new[] { new OncologyAnalyzer.CopyNumberRegion("17", 200, 199) }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IdentifyAmplifiedOncogenes(
            new[] { new OncologyAnalyzer.CopyNumberRegion("", 1, 2) }));
    }

    [Test]
    public void DefaultOncogeneLoci_Gistic2Gencode22Coordinates()
    {
        Assert.That(OncologyAnalyzer.DefaultOncogeneLoci, Is.EqualTo(new[]
        {
            new OncologyAnalyzer.GeneLocus("ERBB2", "17", 39_687_914, 39_730_426),
            new OncologyAnalyzer.GeneLocus("MYC", "8", 127_735_434, 127_741_434),
            new OncologyAnalyzer.GeneLocus("EGFR", "7", 55_019_021, 55_256_620),
            new OncologyAnalyzer.GeneLocus("CCND1", "11", 69_641_087, 69_654_474),
            new OncologyAnalyzer.GeneLocus("MDM2", "12", 68_808_172, 68_850_686),
            new OncologyAnalyzer.GeneLocus("CDK4", "12", 57_747_727, 57_756_013),
        }));
    }

    #endregion

    #region F55 — GetChromosomeArmLengths (UCSC cytoBand, GISTIC2 arm split)

    // FIN-B24 WP26 / F55. Values computed in Python from gistic2 26c590bd refgenes/hg38.UCSC.add_mir.160920/cytoBand.txt
    // (hg38) and the cyto struct of support/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat (hg19, Octave load);
    // all 96 arm lengths also equal GISTIC2 normalize_by_arm_length.m chrarms{1}.length run in Octave on the same cyto.
    [TestCase(OncologyAnalyzer.ReferenceGenome.GRCh38, 1_030_800_000L, 2_057_469_832L, 123_400_000L, 125_556_422L, 121_700_000L, 125_100_000L)]
    [TestCase(OncologyAnalyzer.ReferenceGenome.GRCh37, 1_040_600_000L, 2_055_077_412L, 125_000_000L, 124_250_621L, 121_500_000L, 128_900_000L)]
    public void GetChromosomeArmLengths_SumsAndChr1_MatchUcscCytoBand(
        OncologyAnalyzer.ReferenceGenome genome, long sumP, long sumQ, long chr1P, long chr1Q, long acenStart, long acenEnd)
    {
        var arms = OncologyAnalyzer.GetChromosomeArmLengths(genome);
        Assert.That(arms.Select(a => a.Chromosome), Is.EqualTo(Enumerable.Range(1, 22).Select(i => i.ToString(CultureInfo.InvariantCulture)).Append("X").Append("Y")));
        Assert.That(arms.Sum(a => a.PArmLength), Is.EqualTo(sumP), "Σ p-arm (Python, cytoBand)");
        Assert.That(arms.Sum(a => a.QArmLength), Is.EqualTo(sumQ), "Σ q-arm (Python, cytoBand)");
        Assert.That((arms[0].PArmLength, arms[0].QArmLength, arms[0].CentromereStart, arms[0].CentromereEnd),
            Is.EqualTo((chr1P, chr1Q, acenStart, acenEnd)), "chr1 acen rows");
        // Chromosome ends = UCSC chrom.sizes (existing GetAutosomeLengths table).
        Assert.That(arms.Take(22).Select(a => a.Length), Is.EqualTo(OncologyAnalyzer.GetAutosomeLengths(genome)));
        Assert.That(arms.All(a => a.CentromereStart < a.PArmEnd && a.PArmEnd < a.CentromereEnd && a.CentromereEnd < a.Length), Is.True);
    }

    // Acrocentric / sex chromosomes (GISTIC2 split = p/q acen boundary): spot values from the cytoBand acen rows.
    [Test]
    public void GetChromosomeArmLengths_SpotValues_MatchCytoBandAcenRows()
    {
        var g38 = OncologyAnalyzer.GetChromosomeArmLengths(OncologyAnalyzer.ReferenceGenome.GRCh38);
        var g37 = OncologyAnalyzer.GetChromosomeArmLengths(OncologyAnalyzer.ReferenceGenome.GRCh37);
        Assert.Multiple(() =>
        {
            Assert.That(g38[12], Is.EqualTo(new OncologyAnalyzer.ChromosomeArms("13", 16_500_000L, 17_700_000L, 18_900_000L, 114_364_328L)));
            Assert.That(g38[21], Is.EqualTo(new OncologyAnalyzer.ChromosomeArms("22", 13_700_000L, 15_000_000L, 17_400_000L, 50_818_468L)));
            Assert.That(g38[22], Is.EqualTo(new OncologyAnalyzer.ChromosomeArms("X", 58_100_000L, 61_000_000L, 63_800_000L, 156_040_895L)));
            Assert.That(g38[23], Is.EqualTo(new OncologyAnalyzer.ChromosomeArms("Y", 10_300_000L, 10_400_000L, 10_600_000L, 57_227_415L)));
            Assert.That(g37[8], Is.EqualTo(new OncologyAnalyzer.ChromosomeArms("9", 47_300_000L, 49_000_000L, 50_700_000L, 141_213_431L)));
            Assert.That(g37[23].QArmLength, Is.EqualTo(46_873_566L));
            Assert.That(() => OncologyAnalyzer.GetChromosomeArmLengths((OncologyAnalyzer.ReferenceGenome)99),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        });
    }

    #endregion
}
