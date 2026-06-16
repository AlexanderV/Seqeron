// ONCO-CNA-002 — Focal Amplification Detection
// Evidence: docs/Evidence/ONCO-CNA-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CNA-002.md
// Source: Mermel CH et al. (2011). GISTIC2.0. Genome Biology 12:R41.
//         https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
//         GISTIC2 docs broad_len_cutoff=0.98, t_amp=0.1; NCBI Gene oncogene arms.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberArmSegment;
using Thresholds = Seqeron.Genomics.Oncology.OncologyAnalyzer.FocalAmplificationThresholds;

namespace Seqeron.Genomics.Tests;

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

    // C3 — null amplifications ⇒ ArgumentNullException.
    [Test]
    public void IdentifyAmplifiedOncogenes_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.IdentifyAmplifiedOncogenes(null!),
            "Null amplifications must throw ArgumentNullException.");
    }

    #endregion
}
