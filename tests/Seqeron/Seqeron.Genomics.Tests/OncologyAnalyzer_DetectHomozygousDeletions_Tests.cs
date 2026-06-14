// ONCO-CNA-003 — Homozygous (Deep) Deletion Detection
// Evidence: docs/Evidence/ONCO-CNA-003-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CNA-003.md
// Source: Cheng J et al. (2017). Pan-cancer homozygous deletions. Nat Commun 8:1221.
//         https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/  (homozygous = zero copies of both alleles, total CN 0)
//         cBioPortal discrete CNA: -2 Deep Deletion = homozygous; -1 shallow = heterozygous.
//         https://docs.cbioportal.org/file-formats/  https://docs.cbioportal.org/user-guide/faq/
//         CNVkit absolute_threshold integer-CN (defaults -1.1,-0.25,0.2,0.7); NCBI Gene arms.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberArmSegment;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_DetectHomozygousDeletions_Tests
{
    // Fixed arm length; arm fraction is irrelevant to deletion calling (only log2/CN matter).
    private const long Arm = 1_000_000;

    private static Segment Seg(string arm, double log2) => new(arm, 0, 1_000, Arm, log2);

    #region DetectHomozygousDeletions

    // M1 — log2 = -2.0 ⇒ CN 0 (DeepDeletion) ⇒ homozygous deletion.
    // Source: Cheng 2017 (total CN 0); cBioPortal -2; CNVkit (<= -1.1 ⇒ CN 0).
    [Test]
    public void DetectHomozygousDeletions_DeepDeletionSegment_IsReported()
    {
        var segments = new[] { Seg("9p", -2.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "log2 = -2.0 classifies to integer CN 0 (total CN 0 = both alleles lost), so it is a homozygous deletion.");
            Assert.That(result[0].Arm, Is.EqualTo("9p"),
                "The reported segment must be the input CN-0 segment.");
        });
    }

    // M2 — log2 = -0.5 ⇒ CN 1 (single-copy / heterozygous loss) ⇒ NOT homozygous.
    // Source: cBioPortal -1 = shallow/heterozygous (not -2); Cheng 2017 (one allele remains).
    [Test]
    public void DetectHomozygousDeletions_SingleCopyLoss_NotReported()
    {
        var segments = new[] { Seg("3p", -0.5) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "log2 = -0.5 is integer CN 1 (heterozygous single-copy loss), which is NOT a homozygous deletion.");
    }

    // M3 — log2 = 0.0 ⇒ CN 2 (diploid) ⇒ NOT reported.
    // Source: cBioPortal 0 = diploid.
    [Test]
    public void DetectHomozygousDeletions_NeutralSegment_NotReported()
    {
        var segments = new[] { Seg("10q", 0.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "A copy-number-neutral segment (CN 2) must not be reported as a homozygous deletion.");
    }

    // M4 — gain (log2 0.5 ⇒ CN 3) and amplification (log2 1.0 ⇒ CN >=4) ⇒ NOT reported.
    // Source: cBioPortal 1 = gain, 2 = amplification (neither is a deletion).
    [Test]
    public void DetectHomozygousDeletions_GainAndAmplification_NotReported()
    {
        var segments = new[] { Seg("8q", 0.5), Seg("17q", 1.0) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "Gain (CN 3) and amplification (CN >=4) are increases, not homozygous deletions.");
    }

    // M5 — mixed [CN0, CN1, CN0, CN2] ⇒ only the two CN-0 segments, in input order.
    // Source: INV-1 (CN 0 only); INV-3 (order-preserving filter).
    [Test]
    public void DetectHomozygousDeletions_MixedSet_ReturnsCn0InOrder()
    {
        var segments = new[]
        {
            Seg("9p", -2.0),  // CN 0
            Seg("3p", -0.5),  // CN 1
            Seg("10q", -1.5), // CN 0
            Seg("1q", 0.0),   // CN 2
        };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result.Select(s => s.Arm), Is.EqualTo(new[] { "9p", "10q" }),
                "Only the two CN-0 segments are reported, preserving input order.");
            Assert.That(result, Has.Count.EqualTo(2),
                "Exactly the two homozygous-deletion segments are reported.");
        });
    }

    // M6 — log2 exactly -1.1 (the deletion cutoff) ⇒ CN 0 (<= cutoff) ⇒ reported.
    // Source: CNVkit assigns CN by "less than or equal to each threshold in sequence".
    [Test]
    public void DetectHomozygousDeletions_Log2AtDeletionCutoff_IsReported()
    {
        var segments = new[] { Seg("9p", -1.1) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Has.Count.EqualTo(1),
            "log2 exactly at the deletion cutoff -1.1 is <= the cutoff, so CN 0 (homozygous).");
    }

    // M7 — log2 just above the cutoff (-1.0999) ⇒ CN 1 ⇒ NOT reported.
    // Source: CNVkit threshold boundary (> -1.1 ⇒ CN 1).
    [Test]
    public void DetectHomozygousDeletions_Log2JustAboveCutoff_NotReported()
    {
        var segments = new[] { Seg("9p", -1.0999) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "log2 = -1.0999 is above the -1.1 cutoff, so CN 1 (heterozygous), not homozygous.");
    }

    // M14 — custom thresholds: raising the deletion cutoff to -0.4 makes log2 = -0.5 ⇒ CN 0.
    // Source: CNVkit thresholds are parameters of the calling.
    [Test]
    public void DetectHomozygousDeletions_CustomThresholds_ShiftsCn0Boundary()
    {
        var segments = new[] { Seg("10q", -0.5) };
        var custom = new[] { -0.4, -0.25, 0.2, 0.7 }; // deletion cutoff now -0.4

        var withDefault = OncologyAnalyzer.DetectHomozygousDeletions(segments);
        var withCustom = OncologyAnalyzer.DetectHomozygousDeletions(segments, custom);

        Assert.Multiple(() =>
        {
            Assert.That(withDefault, Is.Empty,
                "Under default thresholds log2 = -0.5 is CN 1 (not homozygous).");
            Assert.That(withCustom, Has.Count.EqualTo(1),
                "With deletion cutoff -0.4, log2 = -0.5 (<= -0.4) is CN 0 (homozygous).");
        });
    }

    // S3 — ploidy parameter feeds CNVkit n = ploidy*2^log2; a triploid no-call/neutral is not homozygous.
    // Source: CNVkit n = ploidy*2^log2; NaN no-call returns neutral reference CN (rounded ploidy).
    [Test]
    public void DetectHomozygousDeletions_TriploidPloidy_RespectsBoundary()
    {
        // log2 = -2.0 is <= -1.1 (deletion cutoff is on log2, not ploidy-scaled) ⇒ CN 0 regardless of ploidy.
        var deep = new[] { Seg("9p", -2.0) };
        // NaN under triploid is a neutral no-call (CN = round(3) = 3), not homozygous.
        var noCall = new[] { Seg("9p", double.NaN) };

        var deepResult = OncologyAnalyzer.DetectHomozygousDeletions(deep, thresholds: null, ploidy: 3.0);
        var noCallResult = OncologyAnalyzer.DetectHomozygousDeletions(noCall, thresholds: null, ploidy: 3.0);

        Assert.Multiple(() =>
        {
            Assert.That(deepResult, Has.Count.EqualTo(1),
                "log2 <= -1.1 is CN 0 by the log2 threshold, independent of ploidy.");
            Assert.That(noCallResult, Is.Empty,
                "A NaN log2 is a neutral no-call (CN = rounded ploidy), not a homozygous deletion.");
        });
    }

    // C1 — empty input ⇒ empty result.
    [Test]
    public void DetectHomozygousDeletions_EmptyInput_ReturnsEmpty()
    {
        var result = OncologyAnalyzer.DetectHomozygousDeletions(Array.Empty<Segment>());

        Assert.That(result, Is.Empty, "An empty segment set yields no homozygous deletions.");
    }

    // C2 — null input ⇒ ArgumentNullException.
    [Test]
    public void DetectHomozygousDeletions_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(null!),
            "Null segments must throw ArgumentNullException.");
    }

    // C3 — invalid segment (End <= Start) ⇒ ArgumentException.
    [Test]
    public void DetectHomozygousDeletions_InvalidSegment_Throws()
    {
        var bad = new[] { new Segment("9p", 1_000, 1_000, Arm, -2.0) }; // End == Start

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectHomozygousDeletions(bad),
            "A segment with End <= Start must throw ArgumentException.");
    }

    // C5 — NaN log2 ⇒ neutral no-call ⇒ NOT reported.
    // Source: CNVkit NaN log2 is a no-call returning the neutral reference CN.
    [Test]
    public void DetectHomozygousDeletions_NaNLog2_NotReported()
    {
        var segments = new[] { Seg("9p", double.NaN) };

        var result = OncologyAnalyzer.DetectHomozygousDeletions(segments);

        Assert.That(result, Is.Empty,
            "A NaN log2 is a neutral no-call (CN = ploidy), not a homozygous deletion.");
    }

    #endregion

    #region IsHomozygousDeletion

    // M8 — predicate: CN-0 segment true, CN-1 segment false.
    // Source: INV-1 (CN 0 = homozygous); INV-2 (CN 1 heterozygous, not).
    [Test]
    public void IsHomozygousDeletion_Cn0AndCn1_TrueThenFalse()
    {
        var cn0 = Seg("9p", -2.0);
        var cn1 = Seg("3p", -0.5);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(cn0), Is.True,
                "log2 = -2.0 (CN 0) is a homozygous deletion.");
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(cn1), Is.False,
                "log2 = -0.5 (CN 1) is a heterozygous loss, not homozygous.");
        });
    }

    #endregion

    #region IdentifyDeletedTumorSuppressors

    // M9 — 17p ⇒ TP53. Source: NCBI Gene TP53 17p13.1.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm17p_MapsTp53()
    {
        var dels = new[] { Seg("17p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "TP53" }),
            "A homozygous deletion on 17p maps to TP53 (NCBI Gene 17p13.1).");
    }

    // M10 — 13q ⇒ RB1 and BRCA2 (both on 13q), in panel order RB1 before BRCA2.
    // Source: NCBI Gene RB1 13q14.2, BRCA2 13q13.1.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm13q_MapsRb1AndBrca2()
    {
        var dels = new[] { Seg("13q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "RB1", "BRCA2" }),
            "13q carries both RB1 (13q14.2) and BRCA2 (13q13.1); both reported in panel order.");
    }

    // M11 — 9p ⇒ CDKN2A. Source: NCBI Gene CDKN2A 9p21.3.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm9p_MapsCdkn2a()
    {
        var dels = new[] { Seg("9p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "CDKN2A" }),
            "A homozygous deletion on 9p maps to CDKN2A (NCBI Gene 9p21.3).");
    }

    // M12 — 10q ⇒ PTEN. Source: NCBI Gene PTEN 10q23.31.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm10q_MapsPten()
    {
        var dels = new[] { Seg("10q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "PTEN" }),
            "A homozygous deletion on 10q maps to PTEN (NCBI Gene 10q23.31).");
    }

    // M13 — 17q ⇒ BRCA1. Source: NCBI Gene BRCA1 17q21.31.
    [Test]
    public void IdentifyDeletedTumorSuppressors_Arm17q_MapsBrca1()
    {
        var dels = new[] { Seg("17q", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "BRCA1" }),
            "A homozygous deletion on 17q maps to BRCA1 (NCBI Gene 17q21.31).");
    }

    // S1 — non-panel arm (1p) ⇒ no gene.
    [Test]
    public void IdentifyDeletedTumorSuppressors_NonPanelArm_ReturnsEmpty()
    {
        var dels = new[] { Seg("1p", -2.0) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.Empty,
            "An arm with no panel tumour suppressor (1p) maps to no gene.");
    }

    // S2 — two deletions on 13q ⇒ RB1 and BRCA2 each once.
    [Test]
    public void IdentifyDeletedTumorSuppressors_DuplicateArm_GenesOnce()
    {
        var dels = new[] { Seg("13q", -2.0), Seg("13q", -1.5) };

        var genes = OncologyAnalyzer.IdentifyDeletedTumorSuppressors(dels);

        Assert.That(genes, Is.EqualTo(new[] { "RB1", "BRCA2" }),
            "Repeated deletions on the same arm report each gene once.");
    }

    // C4 — null ⇒ ArgumentNullException.
    [Test]
    public void IdentifyDeletedTumorSuppressors_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.IdentifyDeletedTumorSuppressors(null!),
            "Null deletions must throw ArgumentNullException.");
    }

    #endregion
}
