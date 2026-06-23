// ONCO-HRD-001 — Homologous Recombination Deficiency (HRD) composite score
// Evidence: docs/Evidence/ONCO-HRD-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-HRD-001.md
// Source: Telli ML et al. (2016). Clin Cancer Res 22(15):3764–3773.
//         HRD score = unweighted sum of LOH + TAI + LST; HR deficiency defined as HRD score >= 42.
//         HRD-TAI: Birkbak NJ et al. (2012). Cancer Discov 2(4):366. HRD-LST: Popova T et al. (2012).
//         Cancer Res 72(21):5454. TAI/LST derivation + centromere coords: scarHRD calc.ai_new / calc.lst
//         (https://github.com/sztup/scarHRD) and UCSC cytoBand acen (hg38/hg19), retrieved 2026-06-23.

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

    #region CalculateHrdTaiScore (TAI derived from segments — Birkbak 2012 / scarHRD calc.ai_new)

    // GRCh38 chr1 centromere (UCSC cytoBand acen) = [121_700_000, 125_100_000].
    private const long Chr1CentromereStart = 121_700_000L;
    private const long Chr1CentromereEnd = 125_100_000L;

    // M13 — TAI counts a p-telomeric imbalance (first segment, end before centromere start) AND a
    // q-telomeric imbalance (last segment, start after centromere end) on the same chromosome → 2.
    // scarHRD calc.ai_new: AI==2 at a chromosome end on one arm relative to the centromere → AI<-1.
    [Test]
    public void CalculateHrdTaiScore_BothTelomericArmsImbalanced_CountsTwo()
    {
        var segments = new[]
        {
            new Segment("1", 0, 50_000_000, 2, 1),                       // first, imbalanced, end < cen start → p-telomeric
            new Segment("1", 50_000_000, Chr1CentromereStart, 1, 1),     // balanced interior
            new Segment("1", 130_000_000, 248_956_422, 2, 1),            // last, imbalanced, start > cen end → q-telomeric
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(2),
            "Both the p-arm-terminal and q-arm-terminal imbalanced segments touch a sub-telomere without crossing the centromere → TAI = 2 (Birkbak 2012; scarHRD calc.ai_new).");
    }

    // M14 — an interstitial imbalance (neither first nor last) is NOT telomeric → not counted.
    [Test]
    public void CalculateHrdTaiScore_InterstitialImbalance_NotCounted()
    {
        var segments = new[]
        {
            new Segment("1", 0, 50_000_000, 1, 1),                  // first, balanced
            new Segment("1", 50_000_000, 110_000_000, 3, 1),        // interior, imbalanced (interstitial AI=2)
            new Segment("1", 130_000_000, 248_956_422, 1, 1),       // last, balanced
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(0),
            "An imbalanced interstitial segment (neither chromosome end) is interstitial AI, not telomeric → TAI = 0 (scarHRD calc.ai_new).");
    }

    // M15 — a first imbalanced segment whose END crosses the centromere start is NOT telomeric.
    // scarHRD condition is strict: sample.chrom.seg[1,4] < chrominfo[i,2].
    [Test]
    public void CalculateHrdTaiScore_FirstSegmentCrossingCentromere_NotCounted()
    {
        var segments = new[]
        {
            new Segment("1", 0, 130_000_000, 2, 1),                 // first, imbalanced, END 130M > cen start 121.7M → crosses
            new Segment("1", 130_000_000, 248_956_422, 1, 1),       // last, balanced
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(0),
            "A terminal imbalanced segment that crosses the centromere (end ≥ centromere start) is not telomeric → TAI = 0 (Birkbak 2012: must not cross the centromere).");
    }

    // M16 — a single imbalanced segment spanning the chromosome is whole-chromosome AI (scarHRD AI=3), not telomeric.
    [Test]
    public void CalculateHrdTaiScore_SingleImbalancedSegment_WholeChromosomeNotCounted()
    {
        var segments = new[] { new Segment("1", 0, 248_956_422, 2, 1) }; // one segment, imbalanced

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(0),
            "A single imbalanced segment is whole-chromosome AI (AI=3 in scarHRD), never telomeric → TAI = 0.");
    }

    // M17 — sub-1 Mb segments are dropped before TAI assignment (scarHRD min.size = 1e6); a < 1 Mb terminal
    // imbalanced fragment does not produce a telomeric event.
    [Test]
    public void CalculateHrdTaiScore_SubMegabaseTerminalSegment_Dropped()
    {
        var segments = new[]
        {
            new Segment("1", 0, 500_000, 2, 1),                     // < 1 Mb imbalanced terminal fragment → dropped
            new Segment("1", 500_000, 248_956_422, 1, 1),           // balanced, spans the rest
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(0),
            "The < 1 Mb terminal fragment is removed before AI assignment (scarHRD min.size=1e6), leaving one balanced segment → TAI = 0.");
    }

    // S7 — only the q-telomeric side imbalanced → exactly 1.
    [Test]
    public void CalculateHrdTaiScore_OnlyQArmTelomericImbalanced_CountsOne()
    {
        var segments = new[]
        {
            new Segment("1", 0, 50_000_000, 1, 1),                  // first, balanced
            new Segment("1", 130_000_000, 248_956_422, 2, 1),       // last, imbalanced, start > cen end → q-telomeric
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(1),
            "Only the q-arm-terminal segment is imbalanced and clears the centromere end → TAI = 1.");
    }

    // S8 — sex chromosomes are excluded from TAI (centromere table is autosome-only).
    [Test]
    public void CalculateHrdTaiScore_SexChromosome_Excluded()
    {
        var segments = new[]
        {
            new Segment("X", 0, 50_000_000, 2, 1),
            new Segment("X", 60_000_000, 156_040_895, 2, 1),
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);

        Assert.That(tai, Is.EqualTo(0),
            "Segments on chrX are excluded from TAI — the centromere table is autosome-only (scarHRD restricts to autosomes).");
    }

    // S9 — GRCh37 uses the hg19 centromere table; chr1 q-arm telomeric event clears the hg19 centromere end.
    [Test]
    public void CalculateHrdTaiScore_GRCh37CentromereTable_Used()
    {
        // hg19 chr1 centromere end = 128_900_000; a last segment starting at 130M is q-telomeric under hg19.
        var segments = new[]
        {
            new Segment("1", 0, 50_000_000, 1, 1),
            new Segment("1", 130_000_000, 249_250_621, 2, 1),
        };

        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments, OncologyAnalyzer.ReferenceGenome.GRCh37);

        Assert.That(tai, Is.EqualTo(1),
            "Under GRCh37 (hg19 chr1 centromere end 128.9 Mb) the last imbalanced segment starting at 130 Mb is q-telomeric → TAI = 1.");
    }

    [Test]
    public void CalculateHrdTaiScore_NullSegments_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateHrdTaiScore(null!),
            "TAI derivation must reject null segments.");
    }

    [Test]
    public void CalculateHrdTaiScore_EmptySegments_ReturnsZero()
    {
        int tai = OncologyAnalyzer.CalculateHrdTaiScore(Array.Empty<Segment>());

        Assert.That(tai, Is.EqualTo(0), "No segments → no telomeric AI → TAI = 0.");
    }

    #endregion

    #region CalculateHrdLstScore (LST derived from segments — Popova 2012 / scarHRD calc.lst)

    // M18 — two adjacent ≥10 Mb segments on the p-arm separated by a < 3 Mb gap, with different allele state
    // (so they are not merged) → 1 LST. scarHRD calc.lst: both flagged large, gap < 3 Mb → one break.
    [Test]
    public void CalculateHrdLstScore_TwoAdjacentLargeArmSegments_CountsOne()
    {
        var segments = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),                  // p-arm, 40 Mb, large
            new Segment("1", 40_000_000, 80_000_000, 1, 1),         // p-arm, 40 Mb, large, different state → not merged
        };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(1),
            "Two adjacent ≥10 Mb p-arm segments of different allele state with a < 3 Mb gap form one large-scale state transition → LST = 1 (Popova 2012; scarHRD calc.lst).");
    }

    // M19 — when one of the two neighbours is < 10 Mb the break is NOT counted (only one side is large).
    [Test]
    public void CalculateHrdLstScore_OneNeighbourBelow10Mb_NotCounted()
    {
        var segments = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),                  // 40 Mb, large
            new Segment("1", 40_000_000, 45_000_000, 1, 1),         // 5 Mb (≥3 Mb so not smoothed) but < 10 Mb → not large
            new Segment("1", 45_000_000, 110_000_000, 3, 1),        // 65 Mb large, but its left neighbour is not large
        };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(0),
            "Neither adjacent pair has BOTH sides ≥10 Mb (the 5 Mb middle segment breaks both pairs) → LST = 0 (Popova 2012: both regions must be ≥10 Mb).");
    }

    // M20 — 3 Mb smoothing removes a short interstitial segment so the two flanking ≥10 Mb segments become
    // adjacent and count as one LST. scarHRD iterative while(<3e6) removal + re-merge.
    [Test]
    public void CalculateHrdLstScore_ShortSegmentSmoothed_ExposesTransition()
    {
        var segments = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),                  // 40 Mb large
            new Segment("1", 40_000_000, 42_000_000, 4, 0),         // 2 Mb (< 3 Mb) → smoothed away
            new Segment("1", 42_000_000, 90_000_000, 1, 1),         // 48 Mb large, different state from the first
        };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(1),
            "The 2 Mb middle segment is smoothed out (< 3 Mb), leaving two adjacent ≥10 Mb p-arm segments → LST = 1 (scarHRD 3 Mb smoothing then break count).");
    }

    // M21 — a chromosome with fewer than 2 segments yields no transition (scarHRD: nrow < 2 → next).
    [Test]
    public void CalculateHrdLstScore_SingleSegment_NotCounted()
    {
        var segments = new[] { new Segment("1", 0, 80_000_000, 2, 1) };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(0), "A single segment cannot form a transition → LST = 0 (scarHRD skips chromosomes with < 2 segments).");
    }

    // M22 — a q-arm transition is counted: two adjacent ≥10 Mb segments past the centromere end.
    [Test]
    public void CalculateHrdLstScore_QArmTransition_CountsOne()
    {
        var segments = new[]
        {
            new Segment("1", Chr1CentromereEnd, 180_000_000, 2, 1),       // q-arm, ~55 Mb large
            new Segment("1", 180_000_000, 248_956_422, 1, 1),            // q-arm, ~69 Mb large, different state
        };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(1),
            "Two adjacent ≥10 Mb q-arm segments of different state with a < 3 Mb gap → one q-arm LST = 1 (scarHRD calc.lst q.arm block).");
    }

    // S10 — sex chromosomes are excluded from LST (scarHRD removes chr23/24/X/Y).
    [Test]
    public void CalculateHrdLstScore_SexChromosome_Excluded()
    {
        var segments = new[]
        {
            new Segment("X", 0, 40_000_000, 2, 1),
            new Segment("X", 40_000_000, 80_000_000, 1, 1),
        };

        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        Assert.That(lst, Is.EqualTo(0),
            "chrX segments are excluded from LST (scarHRD calc.lst removes chr23/24/X/Y) → LST = 0.");
    }

    [Test]
    public void CalculateHrdLstScore_NullSegments_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateHrdLstScore(null!),
            "LST derivation must reject null segments.");
    }

    [Test]
    public void CalculateHrdLstScore_EmptySegments_ReturnsZero()
    {
        int lst = OncologyAnalyzer.CalculateHrdLstScore(Array.Empty<Segment>());

        Assert.That(lst, Is.EqualTo(0), "No segments → no transitions → LST = 0.");
    }

    // P1 — property/invariant (INV-6/INV-7): TAI and LST are order-independent — both sort segments
    // per chromosome before scoring, so permuting the input must not change either count.
    [Test]
    public void CalculateHrdTaiAndLstScore_InputOrderPermuted_ReturnsSameCounts()
    {
        var inOrder = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),
            new Segment("1", 40_000_000, Chr1CentromereStart, 1, 1),
            new Segment("1", 130_000_000, 248_956_422, 2, 1),
        };
        var permuted = new[] { inOrder[2], inOrder[0], inOrder[1] };

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateHrdTaiScore(permuted),
                Is.EqualTo(OncologyAnalyzer.CalculateHrdTaiScore(inOrder)),
                "TAI must be order-independent (segments are sorted per chromosome before scoring).");
            Assert.That(OncologyAnalyzer.CalculateHrdLstScore(permuted),
                Is.EqualTo(OncologyAnalyzer.CalculateHrdLstScore(inOrder)),
                "LST must be order-independent (segments are sorted per chromosome before scoring).");
        });
    }

    #endregion

    #region DetectHRD from segments — all three components derived (scarHRD sum_HRD0)

    // M23 — the all-derived overload computes LOH, TAI and LST from one segment set and sums them.
    // chr1 here yields: TAI = 2 (both terminal imbalances clear the centromere) and LST = 1 (two adjacent
    // ≥10 Mb p-arm segments of different state). LOH = 0 (no minor==0 region > 15 Mb that isn't whole-chr).
    [Test]
    public void DetectHRD_AllDerivedFromSegments_SumsThreeComponents()
    {
        var segments = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),                  // p-arm large, imbalanced, terminal (end < cen start) → p-telomeric TAI; p-arm LST side
            new Segment("1", 40_000_000, Chr1CentromereStart, 1, 1),// p-arm large, balanced, adjacent → forms the LST pair with the first
            new Segment("1", 130_000_000, 248_956_422, 2, 1),       // q-arm, imbalanced, terminal (start > cen end) → q-telomeric TAI
        };

        OncologyAnalyzer.HrdResult result = OncologyAnalyzer.DetectHRD(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result.Components.Loh, Is.EqualTo(0), "No qualifying LOH region (no minor==0 region > 15 Mb) → derived LOH = 0.");
            Assert.That(result.Components.Tai, Is.EqualTo(2),
                "Derived TAI = 2: first segment p-telomeric (end < centromere start) and last segment q-telomeric (start > centromere end).");
            Assert.That(result.Components.Lst, Is.EqualTo(1),
                "Derived LST = 1: the two adjacent ≥10 Mb p-arm segments of different allele state.");
            Assert.That(result.Score, Is.EqualTo(3), "Unweighted sum 0 + 2 + 1 = 3 (scarHRD sum_HRD0; Telli 2016).");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.HrdStatus.HrdNegative),
                "Score 3 < 42 → HRD-negative (Telli 2016).");
        });
    }

    // M24 — the all-derived overload must equal the components overload fed the three standalone derivations
    // (consistency: DetectHRD(segments) == DetectHRD(HrdComponents(DetectLOH, CalculateHrdTaiScore, CalculateHrdLstScore))).
    [Test]
    public void DetectHRD_AllDerived_MatchesStandaloneComponentDerivations()
    {
        var segments = new[]
        {
            new Segment("1", 0, 40_000_000, 2, 1),
            new Segment("1", 40_000_000, Chr1CentromereStart, 1, 1),
            new Segment("1", 130_000_000, 248_956_422, 2, 1),
        };

        int loh = OncologyAnalyzer.DetectLOH(segments).Score;
        int tai = OncologyAnalyzer.CalculateHrdTaiScore(segments);
        int lst = OncologyAnalyzer.CalculateHrdLstScore(segments);

        OncologyAnalyzer.HrdResult viaSegments = OncologyAnalyzer.DetectHRD(segments);
        OncologyAnalyzer.HrdResult viaComponents =
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(loh, tai, lst));

        Assert.That(viaSegments, Is.EqualTo(viaComponents),
            "The all-derived overload must equal the components overload fed the three standalone derivations (LOH/TAI/LST).");
    }

    [Test]
    public void DetectHRD_AllDerived_NullSegments_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectHRD((IEnumerable<Segment>)null!),
            "The all-derived overload must reject null segments.");
    }

    #endregion
}
