// ONCO-HRD-001 — Homologous Recombination Deficiency (HRD) composite score
// Evidence: docs/Evidence/ONCO-HRD-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-HRD-001.md
// Source: Telli ML et al. (2016). Clin Cancer Res 22(15):3764–3773.
//         HRD score = unweighted sum of LOH + TAI + LST; HR deficiency defined as HRD score >= 42.

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.AlleleSpecificSegment;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_CalculateHRDScore_Tests
{
    #region CalculateHRDScore

    // M1 — Telli 2016: HRD = unweighted sum LOH+TAI+LST. 20+15+12 = 47.
    [Test]
    public void CalculateHRDScore_StandardComponents_ReturnsUnweightedSum()
    {
        int score = OncologyAnalyzer.CalculateHRDScore(20, 15, 12);

        Assert.That(score, Is.EqualTo(47),
            "Telli 2016: HRD score is the unweighted sum LOH+TAI+LST = 20+15+12 = 47 (a weighted/averaged form would not give 47).");
    }

    // M2 — Telli 2016 sum on a low-signal triple: 5+4+3 = 12.
    [Test]
    public void CalculateHRDScore_SmallComponents_ReturnsUnweightedSum()
    {
        int score = OncologyAnalyzer.CalculateHRDScore(5, 4, 3);

        Assert.That(score, Is.EqualTo(12),
            "HRD = 5+4+3 = 12; the unweighted sum of the three genomic-scar counts (Telli 2016).");
    }

    // C1 — INV-2: the unweighted sum is order-independent (commutative).
    [Test]
    public void CalculateHRDScore_ComponentOrderPermuted_ReturnsSameSum()
    {
        int a = OncologyAnalyzer.CalculateHRDScore(20, 15, 12);
        int b = OncologyAnalyzer.CalculateHRDScore(12, 20, 15);
        int c = OncologyAnalyzer.CalculateHRDScore(15, 12, 20);

        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(47), "Permutation (20,15,12) sums to 47.");
            Assert.That(b, Is.EqualTo(47), "Permutation (12,20,15) sums to 47 — unweighted sum is commutative (Telli 2016).");
            Assert.That(c, Is.EqualTo(47), "Permutation (15,12,20) sums to 47 — unweighted sum is commutative (Telli 2016).");
        });
    }

    // S2 — component counts are non-negative event counts; a negative count is invalid.
    [Test]
    public void CalculateHRDScore_NegativeComponent_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateHRDScore(-1, 0, 0),
                "Negative LOH count is invalid: components are non-negative event counts (Abkevich 2012).");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateHRDScore(0, -1, 0),
                "Negative TAI count is invalid (Birkbak 2012).");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateHRDScore(0, 0, -1),
                "Negative LST count is invalid (Popova 2012).");
        });
    }

    #endregion

    #region ClassifyHRDStatus

    // M3 — Telli 2016: HR deficiency defined as HRD score >= 42; cutoff is inclusive at 42.
    [Test]
    public void ClassifyHRDStatus_ExactCutoff_ReturnsHrdHigh()
    {
        OncologyAnalyzer.HrdStatus status = OncologyAnalyzer.ClassifyHRDStatus(42);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh),
            "Telli 2016 cutoff is inclusive: score of exactly 42 is HRD-high (a strict >42 rule would wrongly call this negative).");
    }

    // M4 — one below the cutoff is HRD-negative.
    [Test]
    public void ClassifyHRDStatus_OneBelowCutoff_ReturnsHrdNegative()
    {
        OncologyAnalyzer.HrdStatus status = OncologyAnalyzer.ClassifyHRDStatus(41);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
            "Score 41 < 42 is HRD-negative (Telli 2016); confirms the boundary is at 42, not 41.");
    }

    // M5 — well above the cutoff is HRD-high.
    [Test]
    public void ClassifyHRDStatus_WellAboveCutoff_ReturnsHrdHigh()
    {
        OncologyAnalyzer.HrdStatus status = OncologyAnalyzer.ClassifyHRDStatus(100);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh),
            "Score 100 >= 42 is HRD-high (Telli 2016).");
    }

    // S4 — zero score is HRD-negative (well-defined low signal).
    [Test]
    public void ClassifyHRDStatus_ZeroScore_ReturnsHrdNegative()
    {
        OncologyAnalyzer.HrdStatus status = OncologyAnalyzer.ClassifyHRDStatus(0);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
            "Score 0 < 42 is HRD-negative (Telli 2016).");
    }

    // S3 — a negative score is invalid (score is a sum of non-negative counts).
    [Test]
    public void ClassifyHRDStatus_NegativeScore_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyHRDStatus(-1),
            "A negative HRD score is impossible (sum of non-negative counts) and must be rejected.");
    }

    [Test]
    public void HrdHighScoreThreshold_MatchesTelli2016()
    {
        Assert.That(OncologyAnalyzer.HrdHighScoreThreshold, Is.EqualTo(42),
            "The HRD-high cutoff constant must be 42 per Telli 2016 / myChoice CDx.");
    }

    #endregion

    #region DetectHRD

    // M6 — components summing to exactly 42 → HRD-high (boundary via end-to-end path).
    [Test]
    public void DetectHRD_ComponentsSummingTo42_IsHrdHigh()
    {
        OncologyAnalyzer.HrdResult result =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(14, 14, 14));

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(42), "14+14+14 = 42 (unweighted sum, Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh),
                "Score 42 is at the inclusive cutoff → HRD-high (Telli 2016).");
        });
    }

    // M7 — components summing to 41 → HRD-negative (just below boundary).
    [Test]
    public void DetectHRD_ComponentsSummingTo41_IsHrdNegative()
    {
        OncologyAnalyzer.HrdResult result =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(14, 13, 14));

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(41), "14+13+14 = 41 (Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 41 < 42 → HRD-negative (Telli 2016).");
        });
    }

    // M8 — end-to-end high case preserves components, score, and status.
    [Test]
    public void DetectHRD_HighComponents_ReturnsScoreAndHrdHigh()
    {
        var components = new OncologyAnalyzer.HrdComponents(20, 15, 12);

        OncologyAnalyzer.HrdResult result = OncologyAnalyzer.DetectHRD(components);

        Assert.Multiple(() =>
        {
            Assert.That(result.Components, Is.EqualTo(components), "Input components are preserved in the result.");
            Assert.That(result.Score, Is.EqualTo(47), "20+15+12 = 47 (Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh),
                "Score 47 >= 42 → HRD-high (Telli 2016).");
        });
    }

    // M9 — end-to-end negative case.
    [Test]
    public void DetectHRD_LowComponents_ReturnsScoreAndHrdNegative()
    {
        OncologyAnalyzer.HrdResult result =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(5, 4, 3));

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(12), "5+4+3 = 12 (Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 12 < 42 → HRD-negative (Telli 2016).");
        });
    }

    // S1 — near-diploid / low-signal tumour: all components zero.
    [Test]
    public void DetectHRD_NearDiploidZeroComponents_IsHrdNegative()
    {
        OncologyAnalyzer.HrdResult result =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(0, 0, 0));

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(0), "Near-diploid tumour with no scars sums to 0.");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 0 < 42 → HRD-negative (low-signal edge case, Telli 2016).");
        });
    }

    [Test]
    public void DetectHRD_NegativeComponent_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(-1, 0, 0)),
            "DetectHRD must reject negative component counts (delegates validation to CalculateHRDScore).");
    }

    #endregion

    #region DetectHRD from allele-specific segments (LOH derived end-to-end)

    // The ONCO-LOH-001 Evidence dataset: scarHRD calc.hrd yields HRD-LOH = 1 (only the chr1 20 Mb LOH
    // region qualifies; the other segments fail size / homdel / het / whole-chromosome rules). This is the
    // cross-checked LOH count that the segment-driven DetectHRD overload must derive.
    private static List<Segment> LohEvidenceDataset() => new()
    {
        new Segment("1", 0, 20_000_000, 1, 0),          // 20 Mb LOH, chr1 also has a het → counted
        new Segment("1", 20_000_000, 60_000_000, 1, 1), // het retained → not LOH
        new Segment("2", 0, 10_000_000, 2, 0),          // LOH but length ≤ 15 Mb → not counted
        new Segment("3", 0, 16_000_000, 1, 0),          // whole chr3 (all minor=0) → excluded
        new Segment("4", 0, 30_000_000, 0, 0),          // homozygous deletion (major=0) → not LOH
        new Segment("5", 0, 15_000_000, 1, 0),          // length exactly 15 Mb → not > 15 Mb
        new Segment("5", 15_000_000, 50_000_000, 1, 1), // het retained → not LOH
    };

    // M10 — the segment-driven overload derives LOH=1 from the scarHRD-verified dataset and sums it with
    // caller-supplied TAI/LST. 1 (derived LOH) + 25 + 16 = 42 → HRD-high (Telli 2016 inclusive cutoff).
    [Test]
    public void DetectHRD_FromSegments_DerivesLohAndSumsToScore()
    {
        OncologyAnalyzer.HrdResult result = OncologyAnalyzer.DetectHRD(LohEvidenceDataset(), tai: 25, lst: 16);

        Assert.Multiple(() =>
        {
            Assert.That(result.Components.Loh, Is.EqualTo(1),
                "LOH must be DERIVED from the segments (scarHRD calc.hrd on the ONCO-LOH-001 dataset = 1), not supplied; a wrong derivation would not give 1.");
            Assert.That(result.Components.Tai, Is.EqualTo(25), "TAI is the caller-supplied value (25), unchanged.");
            Assert.That(result.Components.Lst, Is.EqualTo(16), "LST is the caller-supplied value (16), unchanged.");
            Assert.That(result.Score, Is.EqualTo(42), "Combined HRD = derived 1 + 25 + 16 = 42 (unweighted sum, Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdHigh),
                "Score 42 is at the inclusive cutoff → HRD-high (Telli 2016).");
        });
    }

    // M11 — the derived-LOH path must agree with the standalone DetectLOH derivation for the same segments,
    // and the overload must equal the components-overload fed that derived LOH (consistency of the two paths).
    [Test]
    public void DetectHRD_FromSegments_MatchesDetectLohPlusComponentsPath()
    {
        List<Segment> segments = LohEvidenceDataset();
        int derivedLoh = OncologyAnalyzer.DetectLOH(segments).Score;

        OncologyAnalyzer.HrdResult viaSegments = OncologyAnalyzer.DetectHRD(segments, tai: 4, lst: 3);
        OncologyAnalyzer.HrdResult viaComponents =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(derivedLoh, 4, 3));

        Assert.Multiple(() =>
        {
            Assert.That(derivedLoh, Is.EqualTo(1), "Standalone DetectLOH derives 1 on the scarHRD dataset.");
            Assert.That(viaSegments, Is.EqualTo(viaComponents),
                "The segment overload must equal the components overload fed the same derived LOH (1+4+3=8 → HRD-negative).");
            Assert.That(viaSegments.Score, Is.EqualTo(8), "Derived 1 + 4 + 3 = 8 (Telli 2016).");
            Assert.That(viaSegments.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 8 < 42 → HRD-negative (Telli 2016).");
        });
    }

    // M12 — no LOH segments → derived LOH = 0; the score is then exactly TAI + LST.
    [Test]
    public void DetectHRD_FromSegmentsWithNoLoh_DerivesZeroLoh()
    {
        var segments = new[]
        {
            new Segment("7", 0, 40_000_000, 1, 1),  // balanced het, not LOH
            new Segment("8", 0, 30_000_000, 2, 2),  // balanced, not LOH
        };

        OncologyAnalyzer.HrdResult result = OncologyAnalyzer.DetectHRD(segments, tai: 10, lst: 5);

        Assert.Multiple(() =>
        {
            Assert.That(result.Components.Loh, Is.EqualTo(0), "No segment is LOH → derived HRD-LOH = 0.");
            Assert.That(result.Score, Is.EqualTo(15), "0 (derived) + 10 + 5 = 15 (Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 15 < 42 → HRD-negative (Telli 2016).");
        });
    }

    // S5 — null segments rejected by the segment-driven overload.
    [Test]
    public void DetectHRD_FromSegments_NullSegments_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectHRD(null!, tai: 0, lst: 0),
            "The segment-driven overload must reject null segments before deriving LOH.");
    }

    // S6 — negative caller-supplied TAI/LST rejected (delegates to CalculateHRDScore validation).
    [Test]
    public void DetectHRD_FromSegments_NegativeTaiOrLst_Throws()
    {
        var segments = new[] { new Segment("1", 0, 20_000_000, 1, 0), new Segment("1", 20_000_000, 60_000_000, 1, 1) };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DetectHRD(segments, tai: -1, lst: 0),
                "Negative caller-supplied TAI must be rejected (Birkbak 2012: counts are non-negative).");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DetectHRD(segments, tai: 0, lst: -1),
                "Negative caller-supplied LST must be rejected (Popova 2012: counts are non-negative).");
        });
    }

    #endregion
}
