// ONCO-LOH-001 — Loss of Heterozygosity (HRD-LOH) detection
// Evidence: docs/Evidence/ONCO-LOH-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-LOH-001.md
// Source: Abkevich V et al. (2012). Br J Cancer 107(10):1776–1782 (PMID 23047548).
//         scarHRD calc.hrd: LOH = minor==0 & major!=0; length > 15e6; whole-chromosome LOH excluded.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.AlleleSpecificSegment;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_DetectLOH_Tests
{
    // The Evidence synthetic dataset (Evidence §Test Datasets); expected HRD-LOH score = 1.
    private static List<Segment> EvidenceDataset() => new()
    {
        new Segment("1", 0, 20_000_000, 1, 0),          // 20 Mb LOH, chr1 also has a het → counted
        new Segment("1", 20_000_000, 60_000_000, 1, 1), // het retained → not LOH
        new Segment("2", 0, 10_000_000, 2, 0),          // LOH but length ≤ 15 Mb → not counted
        new Segment("3", 0, 16_000_000, 1, 0),          // whole chr3 (all minor=0) → excluded
        new Segment("4", 0, 30_000_000, 0, 0),          // homozygous deletion (major=0) → not LOH
        new Segment("5", 0, 15_000_000, 1, 0),          // length exactly 15 Mb → not > 15 Mb
        new Segment("5", 15_000_000, 50_000_000, 1, 1), // het retained → not LOH
    };

    #region DetectLOH / CalculateHrdLohScore

    // M1 — scarHRD calc.hrd on the full Evidence dataset: only chr1's 20 Mb LOH qualifies.
    [Test]
    public void DetectLOH_EvidenceDataset_ScoreIsOne()
    {
        OncologyAnalyzer.LohResult result = OncologyAnalyzer.DetectLOH(EvidenceDataset());

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(1),
                "Evidence dataset: only the chr1 20 Mb LOH region (minor=0, major=1, >15 Mb, not whole-chr) is counted; the other 6 segments each fail one rule (size / homdel / het / whole-chr).");
            Assert.That(result.Regions, Has.Count.EqualTo(1),
                "Regions list must match the score (one qualifying region).");
            Assert.That(result.Regions[0].Chromosome, Is.EqualTo("1"),
                "The single counted region is on chr1.");
            Assert.That(result.Regions[0].Length, Is.EqualTo(20_000_000L),
                "The counted region spans 0..20,000,000 = 20 Mb (length = end - start).");
        });
    }

    // M2 — INV-03/INV-04: a >15 Mb LOH segment on a chromosome that also has a non-LOH segment is counted.
    [Test]
    public void DetectLOH_LohSegmentOver15Mb_IsCounted()
    {
        var segments = new[]
        {
            new Segment("7", 0, 18_000_000, 2, 0),          // 18 Mb LOH
            new Segment("7", 18_000_000, 40_000_000, 1, 1), // het (keeps chr7 from being whole-chr LOH)
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(1),
            "An 18 Mb segment with minor=0, major=2 on a chromosome that is not entirely LOH is one HRD-LOH region.");
    }

    // M3 — INV-04: strict boundary. A segment of length exactly 15,000,000 bp is NOT counted (> not >=).
    [Test]
    public void DetectLOH_SegmentExactly15Mb_IsNotCounted()
    {
        var segments = new[]
        {
            new Segment("8", 0, 15_000_000, 1, 0),          // length exactly 15 Mb
            new Segment("8", 15_000_000, 30_000_000, 1, 1), // het, prevents whole-chr LOH
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(0),
            "scarHRD filter is strict (end-start > 15e6); a region of exactly 15,000,000 bp is excluded. A '>=' implementation would wrongly return 1.");
    }

    // M4 — INV-03: homozygous deletion (minor=0 AND major=0) is NOT loss of heterozygosity.
    [Test]
    public void DetectLOH_HomozygousDeletion_IsNotLoh()
    {
        var segments = new[]
        {
            new Segment("9", 0, 30_000_000, 0, 0),          // both alleles lost → homdel, not LOH
            new Segment("9", 30_000_000, 50_000_000, 1, 1), // het, prevents whole-chr LOH
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(0),
            "scarHRD requires major != 0 (segSamp[,nA] != 0); a 30 Mb homozygous deletion is not counted even though >15 Mb.");
    }

    // M5 — INV-03: a heterozygous-retained segment (minor != 0) is not LOH regardless of size.
    [Test]
    public void DetectLOH_HeterozygousRetained_IsNotLoh()
    {
        var segments = new[]
        {
            new Segment("10", 0, 40_000_000, 1, 1),         // 40 Mb but het retained
            new Segment("10", 40_000_000, 60_000_000, 2, 0),// LOH but kept small to avoid >15 Mb confound
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(1),
            "Only the minor=0 segment (40M..60M = 20 Mb) is LOH; the 40 Mb het segment (minor=1) is never LOH.");
    }

    // M6 — INV-05: whole-chromosome LOH (every segment minor=0) is excluded (Abkevich '< whole chromosome').
    [Test]
    public void DetectLOH_WholeChromosomeLoh_IsExcluded()
    {
        var segments = new[]
        {
            new Segment("11", 0, 16_000_000, 1, 0),          // all of chr11 is LOH...
            new Segment("11", 16_000_000, 60_000_000, 2, 0), // ...so chr11 is whole-chromosome LOH
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(0),
            "scarHRD chrDel: a chromosome where every segment has minor=0 is whole-chromosome LOH and excluded, even though each segment is >15 Mb LOH.");
    }

    // M7 — INV-05: the same LOH segment IS counted once a non-LOH segment makes the chromosome partial.
    [Test]
    public void DetectLOH_LohWithNonLohSegmentOnSameChromosome_IsCounted()
    {
        var segments = new[]
        {
            new Segment("11", 0, 16_000_000, 1, 0),          // 16 Mb LOH
            new Segment("11", 16_000_000, 60_000_000, 1, 1), // het → chr11 no longer whole-chr LOH
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(1),
            "Adding a het segment removes chr11 from chrDel, so the 16 Mb LOH region is now counted (contrast with M6).");
    }

    // M2b — copy-neutral LOH: minor=0 with total CN = 2 (major=2) is still LOH (scarHRD nB==0 & nA!=0,
    // independent of total copy number). A wrong impl keying on total<2 ("copy loss only") would miss this.
    [Test]
    public void DetectLOH_CopyNeutralLoh_IsCounted()
    {
        var segments = new[]
        {
            new Segment("6", 0, 20_000_000, 2, 0),          // 20 Mb copy-NEUTRAL LOH (total CN = 2, minor = 0)
            new Segment("6", 20_000_000, 40_000_000, 1, 1), // het → chr6 not whole-chromosome LOH
        };

        int score = OncologyAnalyzer.CalculateHrdLohScore(segments);

        Assert.That(score, Is.EqualTo(1),
            "scarHRD LOH = (minor==0 & major!=0), independent of total copy number; a 20 Mb copy-neutral LOH (major=2, minor=0) is one HRD-LOH region.");
    }

    #endregion

    #region CalculateLOHFraction

    // M8 — LOH fraction = LOH length / total covered length: 20M / 60M = 0.3333...
    [Test]
    public void CalculateLOHFraction_PartialChromosome_ReturnsLengthWeightedFraction()
    {
        var segments = new[]
        {
            new Segment("1", 0, 20_000_000, 1, 0),          // 20 Mb LOH
            new Segment("1", 20_000_000, 60_000_000, 1, 1), // 40 Mb het
        };

        double fraction = OncologyAnalyzer.CalculateLOHFraction(segments, "1");

        Assert.That(fraction, Is.EqualTo(20_000_000.0 / 60_000_000.0).Within(1e-10),
            "LOH fraction = LOH length (20 Mb) / total covered length (60 Mb) = 1/3.");
    }

    // M9 — chromosome with no LOH → fraction 0.
    [Test]
    public void CalculateLOHFraction_NoLoh_ReturnsZero()
    {
        var segments = new[] { new Segment("2", 0, 50_000_000, 1, 1) };

        double fraction = OncologyAnalyzer.CalculateLOHFraction(segments, "2");

        Assert.That(fraction, Is.EqualTo(0.0).Within(1e-10),
            "A single het segment contributes no LOH length → fraction 0.");
    }

    // M10 — chromosome fully under LOH → fraction 1.
    [Test]
    public void CalculateLOHFraction_FullLoh_ReturnsOne()
    {
        var segments = new[] { new Segment("3", 0, 40_000_000, 1, 0) };

        double fraction = OncologyAnalyzer.CalculateLOHFraction(segments, "3");

        Assert.That(fraction, Is.EqualTo(1.0).Within(1e-10),
            "A single 40 Mb LOH segment covers the whole reported length → fraction 1 (INV-02 upper bound).");
    }

    // INV-02 — fraction is bounded within [0,1] across a mixed chromosome.
    [Test]
    public void CalculateLOHFraction_MixedChromosome_IsWithinUnitInterval()
    {
        var segments = EvidenceDataset();

        double f1 = OncologyAnalyzer.CalculateLOHFraction(segments, "1");
        double f5 = OncologyAnalyzer.CalculateLOHFraction(segments, "5");

        Assert.Multiple(() =>
        {
            Assert.That(f1, Is.InRange(0.0, 1.0), "INV-02: chr1 LOH fraction must lie in [0,1].");
            Assert.That(f5, Is.InRange(0.0, 1.0), "INV-02: chr5 LOH fraction must lie in [0,1].");
            // chr1: 20M LOH / (20M+40M) = 1/3; chr5: 15M LOH / (15M+35M) = 0.3.
            Assert.That(f1, Is.EqualTo(1.0 / 3.0).Within(1e-10), "chr1 LOH fraction = 20M/60M.");
            Assert.That(f5, Is.EqualTo(15_000_000.0 / 50_000_000.0).Within(1e-10), "chr5 LOH fraction = 15M/50M = 0.3.");
        });
    }

    #endregion

    #region Edge cases and invariants

    // S1 — oncoscanR merge: two adjacent 8 Mb LOH segments merge to 16 Mb (>15 Mb) → counted once.
    [Test]
    public void DetectLOH_AdjacentLohSegmentsMergeAcross15Mb_CountedAsOne()
    {
        var segments = new[]
        {
            new Segment("12", 0, 8_000_000, 1, 0),          // 8 Mb LOH
            new Segment("12", 8_000_000, 16_000_000, 1, 0), // adjacent 8 Mb LOH → merge to 16 Mb
            new Segment("12", 16_000_000, 40_000_000, 1, 1),// het → chr12 not whole-chr LOH
        };

        OncologyAnalyzer.LohResult result = OncologyAnalyzer.DetectLOH(segments);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(1),
                "Two adjacent 8 Mb LOH pieces merge into one 16 Mb LOH region (>15 Mb) per oncoscanR; without merging each 8 Mb piece would be <15 Mb and score 0.");
            Assert.That(result.Regions[0].Length, Is.EqualTo(16_000_000L),
                "The merged region spans 0..16,000,000.");
        });
    }

    // S2 — empty input → score 0, fraction 0.
    [Test]
    public void DetectLOH_EmptyInput_ScoreZero()
    {
        OncologyAnalyzer.LohResult result = OncologyAnalyzer.DetectLOH(Array.Empty<Segment>());

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(0), "No segments → no LOH regions.");
            Assert.That(result.Regions, Is.Empty, "No regions returned for empty input.");
            Assert.That(OncologyAnalyzer.CalculateLOHFraction(Array.Empty<Segment>(), "1"),
                Is.EqualTo(0.0).Within(1e-10), "No covered length → fraction 0.");
        });
    }

    // S3 — null input → ArgumentNullException.
    [Test]
    public void DetectLOH_NullInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.DetectLOH(null!), NUnit.Framework.Throws.ArgumentNullException,
                "DetectLOH rejects null segments.");
            Assert.That(() => OncologyAnalyzer.CalculateLOHFraction(null!, "1"), NUnit.Framework.Throws.ArgumentNullException,
                "CalculateLOHFraction rejects null segments.");
            Assert.That(() => OncologyAnalyzer.CalculateLOHFraction(Array.Empty<Segment>(), null!),
                NUnit.Framework.Throws.ArgumentNullException, "CalculateLOHFraction rejects a null chromosome.");
        });
    }

    // S4 — invalid segment (End <= Start) → ArgumentException.
    [Test]
    public void DetectLOH_NonPositiveLength_Throws()
    {
        var bad = new[] { new Segment("1", 100, 100, 1, 0) };

        Assert.That(() => OncologyAnalyzer.DetectLOH(bad), NUnit.Framework.Throws.ArgumentException,
            "A segment with End <= Start has non-positive length and is invalid.");
    }

    // S4b — negative copy number → ArgumentException.
    [Test]
    public void DetectLOH_NegativeCopyNumber_Throws()
    {
        var bad = new[] { new Segment("1", 0, 20_000_000, -1, 0) };

        Assert.That(() => OncologyAnalyzer.DetectLOH(bad), NUnit.Framework.Throws.ArgumentException,
            "Copy numbers must be non-negative.");
    }

    // S5 — fraction for a chromosome absent from the input → 0.
    [Test]
    public void CalculateLOHFraction_UnknownChromosome_ReturnsZero()
    {
        var segments = new[] { new Segment("1", 0, 20_000_000, 1, 0) };

        double fraction = OncologyAnalyzer.CalculateLOHFraction(segments, "22");

        Assert.That(fraction, Is.EqualTo(0.0).Within(1e-10),
            "Chromosome 22 has no covered length in the input → fraction 0.");
    }

    // C1 — INV-06: the HRD-LOH score is independent of input segment order.
    [Test]
    public void DetectLOH_ShuffledInput_ScoreUnchanged()
    {
        var ordered = EvidenceDataset();
        var shuffled = ordered.AsEnumerable().Reverse().ToList();

        int orderedScore = OncologyAnalyzer.CalculateHrdLohScore(ordered);
        int shuffledScore = OncologyAnalyzer.CalculateHrdLohScore(shuffled);

        Assert.Multiple(() =>
        {
            Assert.That(orderedScore, Is.EqualTo(1), "Baseline score is 1 (Evidence dataset).");
            Assert.That(shuffledScore, Is.EqualTo(orderedScore),
                "INV-06: per-chromosome aggregation makes the count order-independent.");
        });
    }

    #endregion
}
