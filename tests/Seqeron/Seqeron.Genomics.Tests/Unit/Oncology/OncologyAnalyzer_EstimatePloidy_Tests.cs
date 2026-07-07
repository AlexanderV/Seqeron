// ONCO-PLOIDY-001 — Tumor Ploidy Estimation + Whole-Genome-Doubling detection
// Evidence: docs/Evidence/ONCO-PLOIDY-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-PLOIDY-001.md
// Source: Patchwork (Genome Biology, PMC4053982) — ploidy = length-weighted mean total CN;
//         Van Loo P et al. (2010) PNAS 107(39):16910–16915 (ASCAT, n-scale ploidy);
//         Bielski CM et al. (2018) Nat Genet 50:1189–1195 / facets-suite is_genome_doubled (PMID 30013179);
//         UCSC hg38.chrom.sizes / hg19.chrom.sizes (reference autosomal genome lengths, retrieved 2026-06-22).

using Segment = Seqeron.Genomics.Oncology.OncologyAnalyzer.AlleleSpecificSegment;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

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

    #region DetectWholeGenomeDoubling (reference chromosome-size table)

    // Authoritative autosomal genome lengths (Σ chr1–22) computed from UCSC *.chrom.sizes.
    // GRCh38 (hg38.chrom.sizes): 2,875,001,522 bp. GRCh37 (hg19.chrom.sizes): 2,881,033,286 bp.
    private const long GRCh38AutosomalGenome = 2_875_001_522L;
    private const long GRCh37AutosomalGenome = 2_881_033_286L;

    // M14 — embedded GRCh38 autosome sizes match UCSC hg38.chrom.sizes exactly (entry 0 = chr1 … 21 = chr22).
    [Test]
    public void GetAutosomeLengths_GRCh38_MatchesUcscChromSizes()
    {
        // Source: UCSC hg38.chrom.sizes (retrieved 2026-06-22); cross-verified Ensembl GRCh38.p14.
        long[] expected =
        {
            248_956_422L, 242_193_529L, 198_295_559L, 190_214_555L, 181_538_259L,
            170_805_979L, 159_345_973L, 145_138_636L, 138_394_717L, 133_797_422L,
            135_086_622L, 133_275_309L, 114_364_328L, 107_043_718L, 101_991_189L,
            90_338_345L,  83_257_441L,  80_373_285L,  58_617_616L,  64_444_167L,
            46_709_983L,  50_818_468L,
        };

        IReadOnlyList<long> actual = OncologyAnalyzer.GetAutosomeLengths(OncologyAnalyzer.ReferenceGenome.GRCh38);

        Assert.That(actual, Is.EqualTo(expected),
            "Embedded GRCh38 autosome lengths (chr1–22) must equal the authoritative UCSC hg38.chrom.sizes values exactly.");
    }

    // M15 — embedded GRCh37 autosome sizes match UCSC hg19.chrom.sizes exactly.
    [Test]
    public void GetAutosomeLengths_GRCh37_MatchesUcscChromSizes()
    {
        // Source: UCSC hg19.chrom.sizes (retrieved 2026-06-22).
        long[] expected =
        {
            249_250_621L, 243_199_373L, 198_022_430L, 191_154_276L, 180_915_260L,
            171_115_067L, 159_138_663L, 146_364_022L, 141_213_431L, 135_534_747L,
            135_006_516L, 133_851_895L, 115_169_878L, 107_349_540L, 102_531_392L,
            90_354_753L,  81_195_210L,  78_077_248L,  59_128_983L,  63_025_520L,
            48_129_895L,  51_304_566L,
        };

        IReadOnlyList<long> actual = OncologyAnalyzer.GetAutosomeLengths(OncologyAnalyzer.ReferenceGenome.GRCh37);

        Assert.That(actual, Is.EqualTo(expected),
            "Embedded GRCh37 autosome lengths (chr1–22) must equal the authoritative UCSC hg19.chrom.sizes values exactly.");
    }

    // M16 — autosomal genome length = Σ chr1–22 from the embedded table (the WGD denominator).
    [Test]
    public void GetAutosomalGenomeLength_BothBuilds_MatchSummedChromSizes()
    {
        long hg38 = OncologyAnalyzer.GetAutosomalGenomeLength(OncologyAnalyzer.ReferenceGenome.GRCh38);
        long hg19 = OncologyAnalyzer.GetAutosomalGenomeLength(OncologyAnalyzer.ReferenceGenome.GRCh37);

        Assert.Multiple(() =>
        {
            // facets-suite: autosomal_genome = sum(chrom_info$size[chr %in% 1:22]).
            Assert.That(hg38, Is.EqualTo(GRCh38AutosomalGenome),
                "GRCh38 autosomal genome must equal Σ(chr1–22) = 2,875,001,522 bp (UCSC hg38.chrom.sizes).");
            Assert.That(hg19, Is.EqualTo(GRCh37AutosomalGenome),
                "GRCh37 autosomal genome must equal Σ(chr1–22) = 2,881,033,286 bp (UCSC hg19.chrom.sizes).");
        });
    }

    // M8 — just over half of the GRCh38 autosomal genome at major CN ≥ 2 → > 0.5 → true.
    [Test]
    public void DetectWholeGenomeDoubling_JustOverHalfOfReferenceGenome_ReturnsTrue()
    {
        // half = 1,437,500,761 bp; +1 bp tips the fraction strictly above 0.5.
        long elevated = (GRCh38AutosomalGenome / 2) + 1;
        var segments = new List<Segment>
        {
            new("1", 0, elevated, 2, 0), // major 2 → elevated, autosome
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "elevated/Σ(chr1–22) > 0.5 (just over half the GRCh38 autosomal genome) → whole-genome doubled.");
    }

    // M9 — exactly half (the larger floor side) of the reference genome → NOT > 0.5 → false (strict).
    [Test]
    public void DetectWholeGenomeDoubling_ExactlyHalfOfReferenceGenome_ReturnsFalse()
    {
        // GRCh38 autosomal genome is even (…522): half is exact, and 0.5 is NOT > 0.5.
        long half = GRCh38AutosomalGenome / 2; // 1,437,500,761
        var segments = new List<Segment>
        {
            new("1", 0, half, 2, 1), // major 2 → elevated, exactly half the reference genome
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "elevated = exactly half the GRCh38 autosomal genome → fraction 0.5, not strictly > 0.5 → not doubled.");
    }

    // M10 — just under half of the reference genome → < 0.5 → false.
    [Test]
    public void DetectWholeGenomeDoubling_JustUnderHalfOfReferenceGenome_ReturnsFalse()
    {
        long elevated = (GRCh38AutosomalGenome / 2) - 1;
        var segments = new List<Segment>
        {
            new("1", 0, elevated, 2, 2), // major 2 → elevated, just under half
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "elevated just under half the GRCh38 autosomal genome → fraction < 0.5 → not doubled.");
    }

    // M11 — all 1:1 (major CN 1) → numerator 0 → not doubled: WGD uses MAJOR, not total CN.
    [Test]
    public void DetectWholeGenomeDoubling_AllBalancedDiploid_ReturnsFalse()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 200_000_000, 1, 1),
            new("2", 0, 200_000_000, 1, 1),
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "Balanced diploid segments have major CN 1; numerator is 0 → not doubled (WGD keys on major CN ≥ 2).");
    }

    // M12 — supplied-segment bias is gone: a tiny fully-amplified region is NOT WGD against the true genome.
    [Test]
    public void DetectWholeGenomeDoubling_SmallFullyAmplifiedRegion_IsNotDoubledAgainstReferenceGenome()
    {
        // Old supplied-length behaviour: 100 Mb all major≥2 → fraction 1.0 → WGD. Reference table fixes this:
        // 100 Mb / 2.875 Gb ≈ 0.035 ≪ 0.5 → NOT WGD (segments do not tile the genome).
        var segments = new List<Segment>
        {
            new("1", 0, 100_000_000, 2, 2), // major 2, only 100 Mb interrogated
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "100 Mb at major CN ≥ 2 is ~3.5% of the GRCh38 autosomal genome → not doubled; the reference denominator removes the supplied-segment bias.");
    }

    // S1 — 2:0 LOH segments (major 2, minor 0) count as elevated against the reference genome.
    [Test]
    public void DetectWholeGenomeDoubling_LohSegmentsOverHalfReference_CountAsElevated()
    {
        long elevated = (GRCh38AutosomalGenome / 2) + 1;
        var segments = new List<Segment>
        {
            new("1", 0, elevated, 2, 0), // major 2 LOH → elevated
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "LOH segments with major CN 2 are elevated; just over half the reference genome → doubled.");
    }

    // S3 — sex-chromosome / non-autosomal segments are excluded from the numerator (facets-suite chrom %in% 1:22).
    [Test]
    public void DetectWholeGenomeDoubling_SexChromosomeSegments_ExcludedFromNumerator()
    {
        long overHalf = (GRCh38AutosomalGenome / 2) + 1;
        var segments = new List<Segment>
        {
            new("chrX", 0, overHalf, 2, 0), // amplified but NOT an autosome → ignored
            new("chrY", 0, overHalf, 2, 0), // amplified but NOT an autosome → ignored
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.False,
            "chrX/chrY amplifications are excluded from the autosomal WGD numerator (facets-suite chrom %in% 1:22) → not doubled.");
    }

    // S4 — "chr"-prefixed autosomes are recognised identically to bare numbers.
    [Test]
    public void DetectWholeGenomeDoubling_ChrPrefixedAutosomes_AreRecognised()
    {
        long overHalf = (GRCh38AutosomalGenome / 2) + 1;
        var segments = new List<Segment>
        {
            new("chr7", 0, overHalf, 2, 0), // "chr"-prefixed autosome, major 2 → elevated
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(segments);

        Assert.That(wgd, Is.True,
            "A 'chr7' segment is an autosome; just over half the reference genome at major CN ≥ 2 → doubled.");
    }

    // M17 — GRCh37 selector uses the hg19 denominator (different boundary from GRCh38).
    [Test]
    public void DetectWholeGenomeDoubling_GRCh37Selector_UsesHg19Denominator()
    {
        // This length is just over half of GRCh37 but NOT over half of GRCh38, so the two builds disagree.
        long overHalfHg19 = (GRCh37AutosomalGenome / 2) + 1; // ≈ 1,440,516,644 bp
        var segments = new List<Segment> { new("1", 0, overHalfHg19, 2, 0) };

        bool wgdHg19 = OncologyAnalyzer.DetectWholeGenomeDoubling(segments, OncologyAnalyzer.ReferenceGenome.GRCh37);
        bool wgdHg38 = OncologyAnalyzer.DetectWholeGenomeDoubling(segments, OncologyAnalyzer.ReferenceGenome.GRCh38);

        Assert.Multiple(() =>
        {
            Assert.That(wgdHg19, Is.True,
                "Just over half the GRCh37 autosomal genome at major CN ≥ 2 → doubled under the GRCh37 denominator.");
            Assert.That(wgdHg38, Is.True,
                "The same length also exceeds half of GRCh38 (smaller autosomal genome) → doubled under GRCh38 too.");
        });
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

    // M18 — empty set against the reference table → fraction 0/genome = 0 → not doubled (no exception).
    [Test]
    public void DetectWholeGenomeDoubling_EmptySegments_ReturnsFalse()
    {
        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoubling(new List<Segment>());

        Assert.That(wgd, Is.False,
            "With a fixed reference denominator, an empty segment set gives numerator 0 → fraction 0 → not doubled.");
    }

    #endregion

    #region DetectWholeGenomeDoublingFromSuppliedLength (legacy supplied-segment denominator)

    // L1 — 60% of supplied length at major CN ≥ 2 → 0.60 > 0.5 → true (legacy behaviour preserved).
    [Test]
    public void DetectWholeGenomeDoublingFromSuppliedLength_SixtyPercentElevated_ReturnsTrue()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 60_000_000, 2, 0), // major 2 → elevated, 60 Mb
            new("2", 0, 40_000_000, 1, 1), // major 1 → not elevated, 40 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(segments);

        Assert.That(wgd, Is.True,
            "0.60 of the supplied length at major CN ≥ 2 exceeds 0.5 → doubled (legacy supplied-length denominator).");
    }

    // L2 — exactly 50% of supplied length → NOT > 0.5 → false (strict threshold).
    [Test]
    public void DetectWholeGenomeDoublingFromSuppliedLength_ExactlyHalf_ReturnsFalse()
    {
        var segments = new List<Segment>
        {
            new("1", 0, 50_000_000, 2, 1), // major 2 → elevated, 50 Mb
            new("2", 0, 50_000_000, 1, 1), // major 1 → not elevated, 50 Mb
        };

        bool wgd = OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(segments);

        Assert.That(wgd, Is.False,
            "Exactly 0.50 of supplied length is not strictly greater than 0.5 → not doubled.");
    }

    // L3 — empty set → ArgumentException (supplied-length denominator undefined).
    [Test]
    public void DetectWholeGenomeDoublingFromSuppliedLength_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(new List<Segment>()),
            "The supplied-length fraction is undefined for an empty segment set; the legacy overload must reject it.");
    }

    // L4 — null input → ArgumentNullException.
    [Test]
    public void DetectWholeGenomeDoublingFromSuppliedLength_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(null!),
            "Null segments must raise ArgumentNullException.");
    }

    #endregion
}
