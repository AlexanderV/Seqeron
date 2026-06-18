// ONCO-PLOIDY-001 — Tumor Ploidy Estimation + Whole-Genome-Doubling detection
// Evidence: docs/Evidence/ONCO-PLOIDY-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-PLOIDY-001.md
// Source: Patchwork (Genome Biology, PMC4053982) — ploidy = length-weighted mean total CN;
//         Van Loo P et al. (2010) PNAS 107(39):16910–16915 (ASCAT, n-scale ploidy);
//         Bielski CM et al. (2018) Nat Genet 50:1189–1195 / facets-suite is_genome_doubled (PMID 30013179).

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.AlleleSpecificSegment;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_EstimatePloidy_Tests
{
    private const double Tolerance = 1e-10;

    #region EstimatePloidy

    // M1 — Patchwork weighted mean: CN 2/4/3 over lengths 100/100/50 Mb → 750M/250M = 3.0.
    [Test]
    public void EstimatePloidy_WorkedExample_ReturnsThree()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 100_000_000, 1, 1), // total CN 2
            new("2", 0, 100_000_000, 2, 2), // total CN 4
            new("3", 0,  50_000_000, 2, 1), // total CN 3
        };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        // ψ = Σ(CN·L)/Σ(L) = (2·100 + 4·100 + 3·50)·1e6 / 250e6 = 750/250 = 3.0 (Patchwork).
        Assert.That(ploidy, Is.EqualTo(3.0).Within(Tolerance),
            "Length-weighted mean of total CN must be 750M/250M = 3.0 (Patchwork PloidyTum definition).");
    }

    // M2 — pure diploid genome (all 1:1) → ψ = 2.0 exactly (n-scale 2n baseline).
    [Test]
    public void EstimatePloidy_PureDiploid_ReturnsTwo()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 50_000_000, 1, 1),
            new("2", 0, 90_000_000, 1, 1),
            new("3", 0, 12_345_678, 1, 1),
        };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        Assert.That(ploidy, Is.EqualTo(2.0).Within(Tolerance),
            "A genome composed entirely of 1:1 segments is diploid: ψ = 2.0 on the n-scale (ASCAT/Patchwork).");
    }

    // M3 — length weighting: long 1:1 (300 Mb) + short 2:2 (10 Mb) must weight toward 2, not 3.
    [Test]
    public void EstimatePloidy_LongDiploidShortAmplified_IsLengthWeighted()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 300_000_000, 1, 1), // total CN 2, long
            new("2", 0,  10_000_000, 2, 2), // total CN 4, short
        };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        // ψ = (2·300 + 4·10)·1e6 / 310e6 = 640/310 ≈ 2.0645; a plain mean would give 3.0.
        Assert.That(ploidy, Is.EqualTo(640.0 / 310.0).Within(Tolerance),
            "Ploidy must be weighted by segment length (640/310 ≈ 2.0645), not an unweighted mean of 3.0.");
    }

    // M4 — single segment (2:1, total 3) → ψ = 3.0.
    [Test]
    public void EstimatePloidy_SingleSegment_ReturnsItsTotalCopyNumber()
    {
        var segments = new List<Segment> { new("1", 0, 75_000_000, 2, 1) };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        Assert.That(ploidy, Is.EqualTo(3.0).Within(Tolerance),
            "The weighted mean over one segment equals that segment's total CN (2+1 = 3).");
    }

    // M5 — empty segment set → undefined (Σ L = 0) → ArgumentException.
    [Test]
    public void EstimatePloidy_EmptySegments_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.EstimatePloidy(new List<Segment>()),
            "Ploidy is undefined for an empty genome (Σ length = 0); the method must reject it.");
    }

    // M6 — segment with End ≤ Start (Length ≤ 0) → ArgumentException.
    [Test]
    public void EstimatePloidy_NonPositiveLength_Throws()
    {
        var segments = new List<Segment> { new("1", 100, 100, 1, 1) }; // End == Start

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.EstimatePloidy(segments),
            "A segment with End ≤ Start has non-positive length and is invalid input.");
    }

    // M7 — negative copy number → ArgumentException.
    [Test]
    public void EstimatePloidy_NegativeCopyNumber_Throws()
    {
        var segments = new List<Segment> { new("1", 0, 1_000_000, -1, 1) };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.EstimatePloidy(segments),
            "Copy numbers must be non-negative; a negative value is invalid input.");
    }

    // M5/M13 — null input → ArgumentNullException (guard contract).
    [Test]
    public void EstimatePloidy_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.EstimatePloidy(null!),
            "Null segments must raise ArgumentNullException.");
    }

    // S2 — a homozygous-deletion (0:0) segment contributes zeros to the weighted mean.
    [Test]
    public void EstimatePloidy_WithHomozygousDeletionSegment_IncludesZeros()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 40_000_000, 0, 0), // total CN 0 (homozygous deletion)
            new("2", 0, 40_000_000, 2, 2), // total CN 4
        };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        // ψ = (0·40 + 4·40)/80 = 160/80 = 2.0; the CN-0 region is counted with its length.
        Assert.That(ploidy, Is.EqualTo(2.0).Within(Tolerance),
            "A CN-0 segment is included in the weighted mean: (0·40+4·40)/80 = 2.0.");
    }

    // C1 — near-triploid genome: ψ exceeds the >2.7n aneuploidy direction (Van Loo et al.).
    [Test]
    public void EstimatePloidy_NearTriploidGenome_ExceedsAneuploidyDirection()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 90_000_000, 2, 1), // total CN 3
            new("2", 0, 90_000_000, 2, 1), // total CN 3
            new("3", 0, 10_000_000, 2, 2), // total CN 4
        };

        double ploidy = OncologyAnalyzer.EstimatePloidy(segments);

        // ψ = (3·90 + 3·90 + 4·10)·1e6 / 190e6 = 580/190 ≈ 3.0526 (> 2.7n aneuploid).
        Assert.That(ploidy, Is.EqualTo(580.0 / 190.0).Within(Tolerance),
            "Near-triploid genome ploidy = 580/190 ≈ 3.053, above the >2.7n aneuploidy direction (Van Loo et al.).");
    }

    #endregion

    #region DetectWholeGenomeDoubling

    // M8 — 60% of length at major CN ≥ 2 → 0.60 > 0.5 → true (facets-suite).
    [Test]
    public void DetectWholeGenomeDoubling_SixtyPercentElevated_ReturnsTrue()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 60_000_000, 2, 0), // major 2 → elevated, 60 Mb
            new("2", 0, 40_000_000, 1, 1), // major 1 → not elevated, 40 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "0.60 of the genome at major CN ≥ 2 exceeds 0.5 → whole-genome doubled (facets-suite).");
    }

    // M9 — exactly 50% elevated → NOT > 0.5 → false (strict threshold).
    [Test]
    public void DetectWholeGenomeDoubling_ExactlyHalfElevated_ReturnsFalse()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 50_000_000, 2, 1), // major 2 → elevated, 50 Mb
            new("2", 0, 50_000_000, 1, 1), // major 1 → not elevated, 50 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "Exactly 0.50 is not strictly greater than 0.5; the genome is NOT doubled (facets-suite strict >).");
    }

    // M10 — 40% elevated → 0.40 ≤ 0.5 → false.
    [Test]
    public void DetectWholeGenomeDoubling_FortyPercentElevated_ReturnsFalse()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 40_000_000, 2, 2), // major 2 → elevated, 40 Mb
            new("2", 0, 60_000_000, 1, 1), // major 1 → not elevated, 60 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "0.40 of the genome at major CN ≥ 2 is below 0.5 → not doubled.");
    }

    // M11 — all 1:1 (total CN 2, major CN 1) → not doubled: WGD uses MAJOR, not total CN.
    [Test]
    public void DetectWholeGenomeDoubling_AllBalancedDiploid_ReturnsFalse()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 100_000_000, 1, 1),
            new("2", 0, 100_000_000, 1, 1),
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "A balanced diploid genome has major CN 1 (total 2); WGD uses major CN ≥ 2, so it is NOT doubled.");
    }

    // M12 — every segment has major CN ≥ 2 → fraction 1.0 > 0.5 → true.
    [Test]
    public void DetectWholeGenomeDoubling_AllMajorElevated_ReturnsTrue()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 30_000_000, 2, 0), // major 2 (LOH)
            new("2", 0, 70_000_000, 2, 2), // major 2
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "Fraction at major CN ≥ 2 is 1.0 > 0.5 → whole-genome doubled.");
    }

    // S1 — 2:0 LOH segments (major 2, minor 0) over >50% count as elevated.
    [Test]
    public void DetectWholeGenomeDoubling_LohSegments_CountAsElevated()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 70_000_000, 2, 0), // major 2 LOH → elevated, 70 Mb
            new("2", 0, 30_000_000, 1, 1), // major 1 → not elevated, 30 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "LOH segments with major CN 2 are elevated; 0.70 > 0.5 → doubled (WGD keys on major CN, not heterozygosity).");
    }

    // M13 — empty set → ArgumentException (fraction undefined).
    [Test]
    public void DetectWholeGenomeDoubling_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoubling(new List<Segment>()),
            "The genome fraction is undefined for an empty segment set; the method must reject it.");
    }

    // M13 — invalid segment length → ArgumentException (shared validation).
    [Test]
    public void DetectWholeGenomeDoubling_NonPositiveLength_Throws()
    {
        var segments = new List<Segment> { new("1", 200, 100, 2, 2) }; // End < Start

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoubling(segments),
            "A segment with End ≤ Start is invalid input.");
    }

    // M13 — negative copy number → ArgumentException.
    [Test]
    public void DetectWholeGenomeDoubling_NegativeCopyNumber_Throws()
    {
        var segments = new List<Segment> { new("1", 0, 1_000_000, 2, -1) };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoubling(segments),
            "Negative copy numbers are invalid input.");
    }

    // M13 — null input → ArgumentNullException.
    [Test]
    public void DetectWholeGenomeDoubling_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoubling(null!),
            "Null segments must raise ArgumentNullException.");
    }

    #endregion
}
