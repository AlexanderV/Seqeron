using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for ChromosomeAnalyzer.cs (checklist 04 rows 48-52:
/// CHROM-TELO/CENT/KARYO/ANEU/SYNT-001).
///
/// The canonical suite under-covered the deterministic scoring/classification helpers
/// (arm-ratio classes, G-band stains, heterochromatin detection, telomere measurement,
/// synteny/rearrangement/aneuploidy callers). These pin the published rules with exact
/// values so relational/arithmetic/logical mutants diverge (baseline 61.4%).
/// </summary>
[TestFixture]
public class ChromosomeAnalyzerMutationTests
{
    // ── Arm ratio: p / q, with edge guards ───────────────────────────────────────────

    [Test]
    public void CalculateArmRatio_IsPArmOverQArm()
    {
        // p = 30, q = 100 − 30 = 70 → 30/70.
        ChromosomeAnalyzer.CalculateArmRatio(30, 100).Should().BeApproximately(30.0 / 70.0, 1e-9);
    }

    [Test]
    [TestCase(0, 100)]     // non-positive centromere
    [TestCase(30, 0)]      // non-positive length
    [TestCase(100, 100)]   // q arm zero
    public void CalculateArmRatio_DegenerateInputs_ReturnZero(int centromere, int length)
    {
        ChromosomeAnalyzer.CalculateArmRatio(centromere, length).Should().Be(0);
    }

    // ── End-to-end: CalculateArmRatio → ClassifyChromosomeByArmRatio per Levan (1964) ──
    // Hand-computed: p = centromerePosition, q = length − p, r = long/short = max(p,q)/min(p,q).
    //   centromere 50/100 → p=q=50 → r=1.00 → metacentric
    //   centromere 40/100 → p=40,q=60 → r=1.50 → metacentric (1.5 < 1.7)
    //   centromere 33/100 → p=33,q=67 → r=2.03 → submetacentric (1.7 < 2.03 ≤ 3.0)
    //   centromere 20/100 → p=20,q=80 → r=4.00 → subtelocentric (3.0 < 4.0 < 7.0)
    //   centromere 10/100 → p=10,q=90 → r=9.00 → acrocentric (≥ 7.0)
    [Test]
    [TestCase(50, 100, "Metacentric")]
    [TestCase(40, 100, "Metacentric")]
    [TestCase(33, 100, "Submetacentric")]
    [TestCase(20, 100, "Subtelocentric")]
    [TestCase(10, 100, "Acrocentric")]
    public void ArmRatioPipeline_ClassifiesPerLevan(int centromere, int length, string expected)
    {
        double ratio = ChromosomeAnalyzer.CalculateArmRatio(centromere, length);
        ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(ratio).Should().Be(expected);
    }

    // ── Levan arm-ratio classification: boundary table ────────────────────────────────
    // Levan, Fredga & Sandberg (1964) Hereditas 52(2):201-220. r = long/short (≥1):
    //   metacentric 1.0..1.7, submetacentric 1.7..3.0, subtelocentric 3.0..7.0,
    //   acrocentric ≥7.0; telocentric = single-arm (one arm absent).
    // ClassifyChromosomeByArmRatio normalises any p/q or q/p input to r ≥ 1, so each row
    // below is given as both forms where the p/q form is the reciprocal. Expected values
    // are hand-derived from the thresholds above (NOT read back from the code).

    [Test]
    // r = q/p form (≥ 1)
    [TestCase(1.0, "Metacentric")]      // perfectly symmetric
    [TestCase(1.43, "Metacentric")]     // r<1.7
    [TestCase(1.7, "Metacentric")]      // boundary 1.7 → metacentric (lower category)
    [TestCase(1.71, "Submetacentric")]  // just above 1.7
    [TestCase(2.0, "Submetacentric")]
    [TestCase(3.0, "Submetacentric")]   // boundary 3.0 → submetacentric
    [TestCase(3.01, "Subtelocentric")]  // just above 3.0
    [TestCase(5.0, "Subtelocentric")]
    [TestCase(6.99, "Subtelocentric")]
    [TestCase(7.0, "Acrocentric")]      // boundary 7.0 → acrocentric
    [TestCase(7.01, "Acrocentric")]
    [TestCase(21.0, "Acrocentric")]
    // p/q form (≤ 1) — same chromosomes, reciprocal ratio; must classify identically
    [TestCase(0.7, "Metacentric")]      // q/p = 1.43
    [TestCase(0.5, "Submetacentric")]   // q/p = 2.0
    [TestCase(0.3, "Subtelocentric")]   // q/p = 3.33  (was mis-called Acrocentric before fix)
    [TestCase(0.2, "Subtelocentric")]   // q/p = 5.0   (was mis-called Acrocentric before fix)
    [TestCase(0.1, "Acrocentric")]      // q/p = 10.0  (was mis-called Telocentric before fix)
    public void ClassifyChromosomeByArmRatio_BoundaryTable(double armRatio, string expected)
    {
        ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(armRatio).Should().Be(expected);
    }

    [Test]
    [TestCase(0.0)]    // degenerate: no short arm measured
    [TestCase(-1.0)]   // invalid
    public void ClassifyChromosomeByArmRatio_DegenerateSingleArm_IsTelocentric(double armRatio)
    {
        // Telocentric is the single-arm case (centromere at one end), not a finite ratio.
        ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(armRatio).Should().Be("Telocentric");
    }

    // ── G-band stain thresholds, gene-density, coordinates, naming, arm switch ────────

    [Test]
    public void PredictGBands_StainGeneDensityCoordinatesAndArmSwitch_AreExact()
    {
        // 3 bands of 10 bp: GC 0% (dark), 40% (medium), 100% (light).
        string seq = "AAAAAAAAAA" + "GCGCAAAAAA" + "GCGCGCGCGC"; // 30 bp
        var bands = ChromosomeAnalyzer.PredictGBands("chr1", seq, bandSize: 10).ToList();

        bands.Should().HaveCount(3);

        bands[0].Stain.Should().Be("gpos100");       // GC 0 < 0.37 → dark
        bands[0].GcContent.Should().BeApproximately(0.0, 1e-9);
        bands[0].GeneDensity.Should().BeApproximately(0.0, 1e-9); // gc × 2
        bands[0].Start.Should().Be(0);
        bands[0].End.Should().Be(9);
        bands[0].Name.Should().Be("chr1p1");

        bands[1].Stain.Should().Be("gpos50");        // 0.37 ≤ 0.40 < 0.45 → medium
        bands[1].GcContent.Should().BeApproximately(0.4, 1e-9);
        bands[1].GeneDensity.Should().BeApproximately(0.8, 1e-9);
        bands[1].Name.Should().Be("chr1p2");

        bands[2].Stain.Should().Be("gneg");          // 1.0 ≥ 0.45 → light
        bands[2].GeneDensity.Should().BeApproximately(2.0, 1e-9);
        bands[2].Name.Should().Be("chr1q1");         // arm switched to q at midpoint
    }

    // ── Telomere measurement (TTAGGG / CCCTAA repeats) ────────────────────────────────

    [Test]
    public void AnalyzeTelomeres_PerfectRepeats_MeasuredExactly()
    {
        string fivePrime = string.Concat(Enumerable.Repeat("CCCTAA", 100));  // revcomp(TTAGGG) × 100
        string filler = new string('A', 60);
        string threePrime = string.Concat(Enumerable.Repeat("TTAGGG", 100)); // × 100
        string seq = fivePrime + filler + threePrime; // 1260 bp

        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chrX", seq);

        r.Has5PrimeTelomere.Should().BeTrue();
        r.Has3PrimeTelomere.Should().BeTrue();
        r.TelomereLength5Prime.Should().Be(600);   // 100 × 6
        r.TelomereLength3Prime.Should().Be(600);
        r.RepeatPurity5Prime.Should().BeApproximately(1.0, 1e-9);
        r.RepeatPurity3Prime.Should().BeApproximately(1.0, 1e-9);
        r.IsCriticallyShort.Should().BeTrue("600 < criticalLength 3000");
    }

    [Test]
    public void AnalyzeTelomeres_NoRepeats_NoTelomere()
    {
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chrX", new string('A', 2000));
        r.Has5PrimeTelomere.Should().BeFalse();
        r.Has3PrimeTelomere.Should().BeFalse();
        r.TelomereLength5Prime.Should().Be(0);
        r.TelomereLength3Prime.Should().Be(0);
    }

    [Test]
    [TestCase(2.0, 1.0, 7000.0, 14000.0)]
    [TestCase(0.5, 1.0, 7000.0, 3500.0)]
    public void EstimateTelomereLengthFromTSRatio_IsProportional(double ts, double refRatio, double refLen, double expected)
    {
        ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(ts, refRatio, refLen)
            .Should().BeApproximately(expected, 1e-9);
    }

    [Test]
    public void EstimateCellDivisions_IsLostOverLossPerDivision()
    {
        // (15000 − 10000) / 50 = 100.
        ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength(10000).Should().BeApproximately(100.0, 1e-9);
        // Non-positive loss → 0.
        ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength(10000, lossPerDivision: 0).Should().Be(0);
    }

    // ── Heterochromatin: repeat-content threshold + region accumulation ───────────────

    [Test]
    public void FindHeterochromatinRegions_FullyRepetitive_YieldsOneSpanningRegion()
    {
        string seq = string.Concat(Enumerable.Repeat("ACGT", 45)); // 180 bp tandem
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions(seq, windowSize: 60).ToList();

        regions.Should().ContainSingle();
        regions[0].Start.Should().Be(0);
        regions[0].End.Should().Be(seq.Length - 1);
        regions[0].Type.Should().Be("Constitutive");
    }

    [Test]
    public void FindHeterochromatinRegions_NoRepeats_YieldsNothing()
    {
        // 'N' k-mers are excluded → repeat content 0 → below threshold everywhere.
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions(new string('N', 180), windowSize: 60).ToList();
        regions.Should().BeEmpty();
    }

    // ── Synteny blocks: collinear runs, strand, gene count, coordinates ───────────────

    private static (string, int, int, string, string, int, int, string) Pair(
        int s1, int e1, int s2, int e2) => ("chrA", s1, e1, $"g{s1}", "chrB", s2, e2, $"h{s2}");

    [Test]
    public void FindSyntenyBlocks_ForwardCollinearRun_OneBlockPlusStrand()
    {
        var pairs = new[] { Pair(0, 10, 0, 10), Pair(20, 30, 20, 30), Pair(40, 50, 40, 50), Pair(60, 70, 60, 70) };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs).ToList();

        blocks.Should().ContainSingle();
        blocks[0].GeneCount.Should().Be(4);
        blocks[0].Species1Start.Should().Be(0);
        blocks[0].Species1End.Should().Be(70);
        blocks[0].Strand.Should().Be('+');
    }

    [Test]
    public void FindSyntenyBlocks_ReverseCollinearRun_MinusStrand()
    {
        var pairs = new[] { Pair(0, 10, 60, 70), Pair(20, 30, 40, 50), Pair(40, 50, 20, 30), Pair(60, 70, 0, 10) };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs).ToList();

        blocks.Should().ContainSingle();
        blocks[0].Strand.Should().Be('-');
        blocks[0].GeneCount.Should().Be(4);
    }

    [Test]
    public void FindSyntenyBlocks_FewerThanMinGenes_Empty()
    {
        var pairs = new[] { Pair(0, 10, 0, 10), Pair(20, 30, 20, 30) }; // 2 < default minGenes 3
        ChromosomeAnalyzer.FindSyntenyBlocks(pairs).Should().BeEmpty();
    }

    // ── Rearrangement detection: inversion vs translocation ───────────────────────────

    private static ChromosomeAnalyzer.SyntenyBlock Block(string c1, int s1, int e1, string c2, char strand) =>
        new(c1, s1, e1, c2, s1, e1, strand, 5, 1.0);

    [Test]
    public void DetectRearrangements_OppositeStrandsSameTarget_IsInversion()
    {
        var blocks = new[] { Block("chrA", 0, 100, "chrB", '+'), Block("chrA", 200, 300, "chrB", '-') };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        events.Should().Contain(e => e.Type == "Inversion");
        var inv = events.First(e => e.Type == "Inversion");
        inv.Position1.Should().Be(100);          // current.Species1End
        inv.Position2.Should().Be(200);          // next.Species1Start
        inv.Size.Should().Be(100);               // 200 − 100
    }

    [Test]
    public void DetectRearrangements_DifferentTargetChromosome_IsTranslocation()
    {
        var blocks = new[] { Block("chrA", 0, 100, "chrB", '+'), Block("chrA", 200, 300, "chrC", '+') };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        events.Should().Contain(e => e.Type == "Translocation");
        events.First(e => e.Type == "Translocation").Chromosome2.Should().Be("chrC");
    }

    // ── Aneuploidy: copy number = round(2^logRatio × 2), logRatio = log2(mean/median) ──

    [Test]
    public void DetectAneuploidy_DoubleDepthBin_IsCopyNumberFour()
    {
        var depth = new[] { ("chr1", 0, 200.0), ("chr1", 1000, 200.0) };
        var states = ChromosomeAnalyzer.DetectAneuploidy(depth, medianDepth: 100).ToList();

        states.Should().ContainSingle();
        states[0].CopyNumber.Should().Be(4);                 // round(2^1 × 2)
        states[0].LogRatio.Should().BeApproximately(1.0, 1e-9); // log2(200/100)
        states[0].Confidence.Should().BeApproximately(1.0, 1e-9);
    }

    [Test]
    public void DetectAneuploidy_NonPositiveMedian_Empty()
    {
        var depth = new[] { ("chr1", 0, 200.0) };
        ChromosomeAnalyzer.DetectAneuploidy(depth, medianDepth: 0).Should().BeEmpty();
    }

    // ── Karyotype: "_N" copy-suffix grouping + aneuploidy terms ───────────────────────

    [Test]
    public void AnalyzeKaryotype_GroupsCopiesByBaseName_AndNamesAneuploidy()
    {
        // Three copies of chr1 (named chr1_1/_2/_3) must group to one base "chr1" → Trisomy.
        // If the "_N" suffix were not stripped, each would be its own monosomic group.
        var chroms = new[]
        {
            ("chr1_1", 1000L, false), ("chr1_2", 1000L, false), ("chr1_3", 1000L, false),
        };
        var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, expectedPloidyLevel: 2);

        k.HasAneuploidy.Should().BeTrue();
        k.Abnormalities.Should().ContainSingle().Which.Should().Be("Trisomy chr1");
    }

    [Test]
    public void AnalyzeKaryotype_CorrectPloidy_NoAneuploidy()
    {
        var chroms = new[] { ("chr1_1", 1000L, false), ("chr1_2", 1000L, false) };
        var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, expectedPloidyLevel: 2);

        k.HasAneuploidy.Should().BeFalse();
        k.Abnormalities.Should().BeEmpty();
        k.AutosomeCount.Should().Be(2);
    }

    // ── Telomere isCritical boundary (strict '<' against criticalLength) ──────────────

    [Test]
    public void AnalyzeTelomeres_FivePrimeOnly_CriticalUsesStrictLessThan()
    {
        // 5' telomere only (3' end is plain filler). length5 = 600.
        string seq = string.Concat(Enumerable.Repeat("CCCTAA", 100)) + new string('A', 2000);

        // criticalLength == length5 → NOT critical (strict '<').
        ChromosomeAnalyzer.AnalyzeTelomeres("chrX", seq, criticalLength: 600)
            .IsCriticallyShort.Should().BeFalse();
        // criticalLength one above → critical.
        ChromosomeAnalyzer.AnalyzeTelomeres("chrX", seq, criticalLength: 601)
            .IsCriticallyShort.Should().BeTrue();
    }

    // ── Rearrangements: deletion and duplication branches ─────────────────────────────

    private static ChromosomeAnalyzer.SyntenyBlock B2(
        string c1, int s1, int e1, string c2, int s2, int e2, char strand) =>
        new(c1, s1, e1, c2, s2, e2, strand, 5, 1.0);

    [Test]
    public void DetectRearrangements_AsymmetricGapSameStrand_IsDeletion()
    {
        // gap1 = 300−100 = 200; gap2 = 110−100 = 10; gap1 > 2·gap2 → deletion in species 2.
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 300, 400, "chrB", 110, 210, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        events.Should().Contain(e => e.Type == "Deletion");
        var del = events.First(e => e.Type == "Deletion");
        del.Size.Should().Be(190);        // gap1 − max(0, gap2) = 200 − 10
        del.Position1.Should().Be(100);
    }

    [Test]
    public void DetectRearrangements_OverlappingSpecies1DifferentTarget_IsDuplication()
    {
        var blocks = new[]
        {
            B2("chrA", 0, 200, "chrB", 0, 200, '+'),
            B2("chrA", 100, 300, "chrC", 500, 700, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        events.Should().Contain(e => e.Type == "Duplication");
        var dup = events.First(e => e.Type == "Duplication");
        dup.Position1.Should().Be(100);   // overlapStart = max(0, 100)
        dup.Size.Should().Be(100);        // overlapEnd(200) − overlapStart(100)
        dup.Chromosome2.Should().Be("chrC");
    }

    // ── Synteny: a gap exceeding maxGap splits a collinear run into two blocks ─────────

    [Test]
    public void FindSyntenyBlocks_GapExceedingMaxGap_SplitsBlock()
    {
        // maxGap default 10 → threshold 10·1e6 = 1e7. Two collinear triples separated by ~5e7.
        var pairs = new[]
        {
            Pair(0, 1_000_000, 0, 1_000_000),
            Pair(2_000_000, 3_000_000, 2_000_000, 3_000_000),
            Pair(4_000_000, 5_000_000, 4_000_000, 5_000_000),
            Pair(60_000_000, 61_000_000, 60_000_000, 61_000_000),
            Pair(62_000_000, 63_000_000, 62_000_000, 63_000_000),
            Pair(64_000_000, 65_000_000, 64_000_000, 65_000_000),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs).ToList();

        blocks.Should().HaveCount(2, "the >1e7 gap breaks collinearity into two blocks");
        blocks.Should().OnlyContain(b => b.GeneCount == 3);
    }

    // ── Heterochromatin DetermineHeterochromatinType: position-based classification ───

    [Test]
    public void FindHeterochromatinRegions_RegionNearStart_IsTelomeric()
    {
        // Repetitive block at the very start, then a long N tail → region midpoint position < 0.05.
        string seq = string.Concat(Enumerable.Repeat("ACGT", 15)) + new string('N', 3000); // 3060 bp
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions(seq, windowSize: 60).ToList();

        regions.Should().NotBeEmpty();
        regions[0].Start.Should().Be(0);
        regions[0].End.Should().Be(29, "region closes at the first low-repeat window i=30 → end i−1");
        regions[0].Type.Should().Be("Telomeric");
    }

    [Test]
    public void FindHeterochromatinRegions_RegionNearMiddle_IsCentromeric()
    {
        // Repetitive block centered in the sequence → region midpoint position in (0.45, 0.55).
        string seq = new string('N', 1500) + string.Concat(Enumerable.Repeat("ACGT", 15)) + new string('N', 1500);
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions(seq, windowSize: 60).ToList();

        regions.Should().ContainSingle();
        regions[0].Type.Should().Be("Centromeric");
    }

    [Test]
    public void FindHeterochromatinRegions_RegionOffCentre_IsConstitutive()
    {
        // Repetitive block at ~20% position → neither telomeric (<0.05/>0.95) nor centromeric (0.45-0.55).
        string seq = new string('N', 600) + string.Concat(Enumerable.Repeat("ACGT", 15)) + new string('N', 2400);
        var regions = ChromosomeAnalyzer.FindHeterochromatinRegions(seq, windowSize: 60).ToList();

        regions.Should().NotBeEmpty();
        regions[0].Type.Should().Be("Constitutive");
    }

    // ── AnalyzeCentromere: deterministic detection (windowSize < 1000 ⇒ GcVariability ≡ 0,
    //    so score == repeat content and the first maximal-repeat window is the centromere) ──

    private const string Diverse40 = "TTGCAACTGGAATCCGTACATGGCAATTCGGATCAACGTT"; // non-repetitive 40-mer

    [Test]
    public void AnalyzeCentromere_RepetitiveBlockBetweenDiverseFlanks_DetectedExactly()
    {
        string seq = Diverse40 + string.Concat(Enumerable.Repeat("ACGT", 30)) + Diverse40; // 200 bp

        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr1", seq, windowSize: 40, minAlphaSatelliteContent: 0.3);

        // First fully-repetitive window begins at index 40 (end 40+40=80).
        r.Start.Should().Be(40);
        r.End.Should().Be(80);
        r.Length.Should().Be(40);
        r.AlphaSatelliteContent.Should().BeApproximately(1.0, 1e-9, "repeat content of a perfect tandem window is 1.0");
        // arm ratio q/p = 140/60 ≈ 2.33 → Submetacentric (Levan), not acrocentric.
        r.CentromereType.Should().Be("Submetacentric");
        r.IsAcrocentric.Should().BeFalse();
    }

    [Test]
    public void AnalyzeCentromere_LargeWindow_ExtendsRightAcrossRepeatBlock()
    {
        // windowSize 120 ⇒ half-window 60 ≥ kmer*2, so the boundary-extension loop runs.
        // Repeat block [120,480) → detected window starts at 120, then extends right to 480.
        // Flank120 is non-repetitive (no repeated 15-mer ⇒ repeat content 0), so it neither
        // seeds nor extends the centromere.
        const string Flank120 =
            "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCG";
        string seq = Flank120                                                // 120 bp diverse prefix
                   + string.Concat(Enumerable.Repeat("ACGT", 90))           // 360 bp tandem
                   + Flank120;                                               // 120 bp diverse suffix
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr1", seq, windowSize: 120, minAlphaSatelliteContent: 0.3);

        r.Start.Should().Be(120, "first maximal-repeat window begins where the tandem starts");
        r.End.Should().Be(480, "right-extension walks 60-bp half-windows to the end of the repeat block");
        r.Length.Should().Be(360);
        r.AlphaSatelliteContent.Should().BeApproximately(1.0, 1e-9);
        r.CentromereType.Should().Be("Metacentric"); // arm ratio 1.0
    }

    [Test]
    public void AnalyzeCentromere_NoRepetitiveContent_ReturnsNoCentromere()
    {
        // All-diverse sequence: repeat content never exceeds minAlphaSatelliteContent → no centromere.
        string seq = string.Concat(Enumerable.Repeat(Diverse40, 5)); // 200 bp, non-repetitive
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr1", seq, windowSize: 40, minAlphaSatelliteContent: 0.95);

        r.Start.Should().BeNull();
        r.End.Should().BeNull();
        r.Length.Should().Be(0);
        r.CentromereType.Should().Be("Unknown");
    }

    // ── G-band stain thresholds at exact GC boundaries (37% dark/medium, 45% medium/light) ──

    [Test]
    public void PredictGBands_ExactThresholdGc_UsesStrictLessThan()
    {
        // 100-bp bands with EXACTLY 37% and 45% GC (achievable only at 100-bp granularity).
        string gc37 = new string('G', 37) + new string('A', 63); // GC = 0.37
        string gc45 = new string('G', 45) + new string('A', 55); // GC = 0.45
        string seq = gc37 + gc45;

        var bands = ChromosomeAnalyzer.PredictGBands("chr1", seq, bandSize: 100).ToList();

        // gcContent 0.37 is NOT < 0.37 → medium (gpos50), not dark.
        bands[0].Stain.Should().Be("gpos50");
        // gcContent 0.45 is NOT < 0.45 → light (gneg), not medium.
        bands[1].Stain.Should().Be("gneg");
    }

    [Test]
    public void PredictGBands_AllN_UsesDefaultHalfGc()
    {
        // No valid bases → total == 0 → default GC 0.5 → light band (gneg), not GC 0 (dark).
        var bands = ChromosomeAnalyzer.PredictGBands("chr1", new string('N', 50), bandSize: 100).ToList();
        bands.Should().ContainSingle();
        bands[0].Stain.Should().Be("gneg");
        bands[0].GcContent.Should().BeApproximately(0.5, 1e-9);
    }

    // ── Rearrangement negative cases: pin the deletion/duplication GUARDS ─────────────

    [Test]
    public void DetectRearrangements_SymmetricGap_IsNotDeletion()
    {
        // gap1 (200) is NOT greater than 2·gap2 (2·150=300) → no deletion (pins gap1 > 2·gap2).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 300, 400, "chrB", 250, 350, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().NotContain(e => e.Type == "Deletion");
    }

    [Test]
    public void DetectRearrangements_OppositeStrands_ProduceNoDeletion()
    {
        // Same target chromosome but opposite strands → inversion, never a deletion
        // (pins the `sameChr && sameStrand` guard against the `||` mutant).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 300, 400, "chrB", 500, 600, '-'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().Contain(e => e.Type == "Inversion");
        events.Should().NotContain(e => e.Type == "Deletion");
    }

    [Test]
    public void DetectRearrangements_OverlappingSameTarget_IsNotDuplication()
    {
        // Overlapping species-1 regions mapping to the SAME species-2 location → not a duplication
        // (pins the sameTarget guard: identical Chr/Start/End).
        var blocks = new[]
        {
            B2("chrA", 0, 200, "chrB", 0, 200, '+'),
            B2("chrA", 100, 300, "chrB", 0, 200, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().NotContain(e => e.Type == "Duplication");
    }

    [Test]
    public void DetectRearrangements_OppositeStrandWithDeletionGaps_StillNoDeletion()
    {
        // Opposite strands but gaps that WOULD satisfy deletion (gap1=200 > 2·gap2=20).
        // The `sameChr && sameStrand` guard must block deletion (pins the `&&` → `||` mutant).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 300, 400, "chrB", 110, 210, '-'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().Contain(e => e.Type == "Inversion");
        events.Should().NotContain(e => e.Type == "Deletion");
    }

    [Test]
    public void DetectRearrangements_GapExactlyTwiceSpecies2Gap_IsNotDeletion()
    {
        // gap1 (300) == 2·gap2 (2·150) → NOT strictly greater → no deletion (pins `>` vs `>=`).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 400, 500, "chrB", 250, 350, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().NotContain(e => e.Type == "Deletion");
    }

    [Test]
    public void DetectRearrangements_ZeroSpecies2Gap_IsDeletion()
    {
        // gap2 == 0 (adjacent in species 2) with gap1 > 0 → deletion (pins `gap2 >= 0` vs `> 0`).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 0, 100, '+'),
            B2("chrA", 300, 400, "chrB", 100, 200, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().Contain(e => e.Type == "Deletion");
    }

    // ── Whole-chromosome aneuploidy from copy-number states ───────────────────────────

    private static ChromosomeAnalyzer.CopyNumberState CN(string chr, int cn) =>
        new(chr, 0, 1, cn, 0.0, 1.0);

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_DominantCopyNumberThree_IsTrisomy()
    {
        var states = new[] { CN("chr21", 3), CN("chr21", 3), CN("chr21", 3), CN("chr21", 3), CN("chr21", 2) };
        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();

        calls.Should().ContainSingle();
        calls[0].Chromosome.Should().Be("chr21");
        calls[0].CopyNumber.Should().Be(3);
        calls[0].Type.Should().Be("Trisomy");
    }

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_DiploidOrSubThreshold_NoCall()
    {
        // All copy-number 2 → not aneuploid.
        ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(new[] { CN("chr1", 2), CN("chr1", 2) })
            .Should().BeEmpty();
        // Dominant fraction below minFraction (0.8): 3 of 5 are CN3 (0.6) → no call.
        var mixed = new[] { CN("chr1", 3), CN("chr1", 3), CN("chr1", 3), CN("chr1", 2), CN("chr1", 2) };
        ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(mixed).Should().BeEmpty();
    }

    // ── Synteny: a species-2 gap exceeding maxGap (with small species-1 gap) splits the run ──

    [Test]
    public void DetectRearrangements_MinusStrandDeletion_UsesMinusStrandGapFormula()
    {
        // Both blocks on '-' strand, same target. Deletion gap2 = current.Start2 − next.End2 = 300−295 = 5;
        // gap1 = 200 > 2·5 → deletion. Pins the '-'-branch subtraction (vs the '+' formula / addition mutant).
        var blocks = new[]
        {
            B2("chrA", 0, 100, "chrB", 300, 400, '-'),
            B2("chrA", 300, 400, "chrB", 285, 295, '-'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        events.Should().Contain(e => e.Type == "Deletion");
        events.First(e => e.Type == "Deletion").Size.Should().Be(195); // gap1 − max(0, gap2) = 200 − 5
    }

    [Test]
    public void DetectRearrangements_PartialTargetMatchOverlap_IsDuplication()
    {
        // Overlapping species-1 regions; same target chr+start but DIFFERENT end → not "sameTarget"
        // → duplication. Pins the 3-way `&&` sameTarget guard against the `||` mutant.
        var blocks = new[]
        {
            B2("chrA", 0, 200, "chrB", 0, 200, '+'),
            B2("chrA", 100, 300, "chrB", 0, 250, '+'),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        events.Should().Contain(e => e.Type == "Duplication");
    }

    [Test]
    public void FindSyntenyBlocks_Species2StartEqualsPrevEnd_IsReverseStrand()
    {
        // curr.Start2 == prev.End2 at every step → strictly-greater forward test is FALSE → '-' strand.
        // Pins `curr.Start2 > prev.End2` against the `>=` mutant (which would flip the strand to '+').
        var pairs = new[] { Pair(0, 10, 0, 10), Pair(20, 30, 10, 20), Pair(40, 50, 20, 30) };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs).ToList();

        blocks.Should().ContainSingle();
        blocks[0].Strand.Should().Be('-');
    }

    [Test]
    public void FindSyntenyBlocks_Species2GapExceedingMaxGap_BreaksCollinearity()
    {
        // First three are collinear; the fourth has a small Species1 gap but a huge Species2 gap,
        // so collinearity (gap1 ≤ T AND gap2 ≤ T) must fail → only the leading triple is a block.
        var pairs = new[]
        {
            Pair(0, 10, 0, 10),
            Pair(20, 30, 20, 30),
            Pair(40, 50, 40, 50),
            Pair(60, 70, 50_000_000, 50_000_010),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs).ToList();

        blocks.Should().ContainSingle();
        blocks[0].GeneCount.Should().Be(3, "the large Species2 gap excludes the 4th gene");
    }

    // ── Alpha-satellite / HOR detection (rows 257 CHROM-ALPHASAT-001, 258 CHROM-HOR-001) ──
    //
    // McNulty & Sullivan 2018 (PMC6121732): a HOR is a block of N distinct ~171 bp monomers tiled with high
    // inter-copy identity; alpha satellite is AT-rich (>50%) and tandemly periodic; the CENP-B box is the
    // 17-bp IUPAC consensus YTTCGTTGGAARCGGGA (Masumoto et al. 1989).

    // A valid realisation of the CENP-B consensus (Y->C, R->A), 17 bp.
    private const string CenpBBox = "CTTCGTTGGAAACGGGA";

    [Test]
    public void DetectHigherOrderRepeat_MonomerLengthBoundary_OneIsValidZeroThrows()
    {
        Action zero = () => ChromosomeAnalyzer.DetectHigherOrderRepeat("AAAA", monomerLength: 0);
        zero.Should().Throw<ArgumentOutOfRangeException>("monomerLength < 1 is rejected");

        ChromosomeAnalyzer.HorResult r = default;
        Action one = () => r = ChromosomeAnalyzer.DetectHigherOrderRepeat("AAAA", monomerLength: 1);
        one.Should().NotThrow("monomerLength = 1 is the smallest legal monomer length");
        r.MonomersPerUnit.Should().Be(1, "four identical 1-bp monomers form a period-1 (non-HOR) array");
        r.HasHigherOrderStructure.Should().BeFalse();
    }

    [Test]
    public void DetectHigherOrderRepeat_PerfectPeriodTwoArray_ReportsPeriodTwoStructure()
    {
        // Two distinct monomers (AAAA, GGGG) tiled three times: adjacent monomers are 0% identical (k=1
        // rejected) but monomers two apart are identical (k=2 accepted) -> period 2.
        string array = string.Concat(Enumerable.Repeat("AAAA" + "GGGG", 3)); // 6 monomers of length 4

        var r = ChromosomeAnalyzer.DetectHigherOrderRepeat(array, monomerLength: 4);

        r.HasHigherOrderStructure.Should().BeTrue("a period-2 HOR is a multi-monomer structure");
        r.MonomersPerUnit.Should().Be(2, "smallest inter-HOR-identical period = 2 monomers");
        r.MonomerCount.Should().Be(6, "24 bp / 4 bp = 6 monomers");
        r.HorCopyNumber.Should().Be(3, "6 monomers / period 2 = 3 HOR copies");
        r.HorUnitLengthBp.Should().Be(8, "unit length = period(2) * monomerLength(4)");
        r.MeanInterHorIdentity.Should().BeApproximately(100.0, 1e-9,
            "monomers one unit apart (AAAA/AAAA, GGGG/GGGG) are 100% identical");
        r.MeanIntraHorIdentity.Should().BeApproximately(0.0, 1e-9,
            "the two monomers WITHIN a unit (AAAA vs GGGG) share no aligned identity");
    }

    [Test]
    public void DetectAlphaSatellite_ExactlyMinimumLength_IsAnalysedNotSkipped()
    {
        // Minimum analysable length = monomer(171) + tolerance(5) + 1 = 177 bp. Exactly 177 must be analysed.
        string atRich = string.Concat(Enumerable.Repeat("AT", 89))[..177]; // 177 bp, all A/T

        var r = ChromosomeAnalyzer.DetectAlphaSatellite(atRich);

        atRich.Length.Should().Be(177, "exactly the minimum analysable length");
        r.AtContent.Should().BeApproximately(1.0, 1e-9,
            "an all-A/T array has AT content 1.0 — proving it was analysed, not skipped as too short");
    }

    [Test]
    public void DetectAlphaSatellite_GcRichPeriodicArray_IsNotAlphaSatellite()
    {
        string gcRich = string.Concat(Enumerable.Repeat("GC", 100)); // 200 bp, periodic but AT-poor

        var r = ChromosomeAnalyzer.DetectAlphaSatellite(gcRich);

        r.AtContent.Should().BeApproximately(0.0, 1e-9, "no A/T bases");
        r.IsAlphaSatellite.Should().BeFalse(
            "periodicity alone is insufficient — alpha satellite must also be AT-rich (>50%)");
    }

    [Test]
    public void DetectAlphaSatellite_CenpBBoxAtLastPosition_IsCounted()
    {
        // 160 bp background + 17 bp box = 177 bp; box starts at offset 160 = 177 - 17 (the LAST scan offset).
        string seq = new string('A', 160) + CenpBBox;

        seq.Length.Should().Be(177);
        ChromosomeAnalyzer.DetectAlphaSatellite(seq).CenpBBoxCount.Should().Be(1,
            "the box at the final scan offset must still be found (inclusive scan bound)");
        ChromosomeAnalyzer.DetectAlphaSatellite(new string('A', 177)).CenpBBoxCount.Should().Be(0,
            "no CENP-B box in a poly-A array");
    }

    [Test]
    public void FindCenpBBoxes_MatchesConsensusExactlyAndRejectsFixedPositionMismatch()
    {
        string withBox = "AAAAA" + CenpBBox + "AAAAA"; // box starts at index 5
        ChromosomeAnalyzer.FindCenpBBoxes(withBox).Should().Equal(new[] { 5 },
            "the box is found at its true 0-based start");

        // Corrupt a fixed (non-degenerate) consensus position: consensus[4] = 'G' must match exactly.
        string broken = "AAAAA" + "CTTCATTGGAAACGGGA" + "AAAAA"; // index-4 G -> A
        ChromosomeAnalyzer.FindCenpBBoxes(broken).Should().BeEmpty(
            "a mismatch at a fixed consensus position (G at index 4) breaks the box");
    }
}
