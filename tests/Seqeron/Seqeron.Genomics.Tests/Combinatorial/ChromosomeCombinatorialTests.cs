namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Chromosome area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Chromosome")]
public class ChromosomeCombinatorialTests
{
    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);
    private static string Repeat(string unit, int times) => string.Concat(Enumerable.Repeat(unit, times));

    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-TELO-001 — Telomere detection (Chromosome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 48.
    // Spec: tests/TestSpecs/CHROM-TELO-001.md (canonical AnalyzeTelomeres).
    // Dimensions: repeatMotif(2) × minRepeats(3) × seqLen(3). Grid 2×3×3 = 18.
    //
    // Model (Blackburn 1991; telomere repeats): a telomere is a tandem array of a short motif
    // at a chromosome end — the TTAGGG repeat on the 3′ strand and its reverse complement
    // (CCCTAA) on the 5′ strand. AnalyzeTelomeres measures the contiguous repeat run at each
    // end and flags a telomere when its length reaches minTelomereLength.
    //
    // The combinatorial point: motif, the minimum-length threshold and the chromosome length
    // interact — the planted run length is measured exactly per motif and end, and the
    // presence flag flips precisely when the run reaches the threshold (minRepeats × motifLen).
    // ═══════════════════════════════════════════════════════════════════════

    private const int PlantedRepeats = 100;

    [Test, Combinatorial]
    public void ChromTelo_DetectsRunAtThreshold_PerMotifAndEnd(
        [Values("TTAGGG", "TTTAGGG")] string motif,
        [Values(50, 100, 150)] int minRepeats,
        [Values(3000, 6000, 10000)] int seqLen)
    {
        int teloLen = PlantedRepeats * motif.Length;
        string telo5 = Repeat(RevComp(motif), PlantedRepeats); // 5′ strand carries the reverse complement
        string telo3 = Repeat(motif, PlantedRepeats);
        string middle = new string('A', seqLen - 2 * teloLen);
        string seq = telo5 + middle + telo3;

        var r = AnalyzeTelo(seq, motif, minRepeats);

        r.TelomereLength3Prime.Should().Be(teloLen, "the contiguous 3′ run is measured exactly");
        r.TelomereLength5Prime.Should().Be(teloLen, "the contiguous 5′ run is measured exactly");
        r.RepeatPurity3Prime.Should().BeApproximately(1.0, 1e-9, "perfect repeats are 100% pure");

        bool expectDetected = PlantedRepeats >= minRepeats;
        r.Has3PrimeTelomere.Should().Be(expectDetected, "presence flips at the minTelomereLength threshold");
        r.Has5PrimeTelomere.Should().Be(expectDetected);
    }

    private static ChromosomeAnalyzer.TelomereResult AnalyzeTelo(string seq, string motif, int minRepeats) =>
        ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, motif,
            searchLength: Math.Max(10000, seq.Length),
            minTelomereLength: minRepeats * motif.Length);

    /// <summary>
    /// Interaction witness: a chromosome with no terminal repeat array has no telomere at
    /// either end, whatever the motif.
    /// </summary>
    [Test]
    public void ChromTelo_NonTelomericEnds_NoTelomere()
    {
        string seq = new string('A', 4000) + new string('C', 4000);
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, "TTAGGG", minTelomereLength: 300);

        r.Has5PrimeTelomere.Should().BeFalse();
        r.Has3PrimeTelomere.Should().BeFalse();
        r.TelomereLength3Prime.Should().Be(0);
    }

    /// <summary>
    /// Interaction witness: the motif matters — a TTAGGG array is not recognised when the
    /// caller searches for the 7-mer TTTAGGG motif (frame-shifted, sub-threshold similarity).
    /// </summary>
    [Test]
    public void ChromTelo_MotifMustMatch()
    {
        string seq = new string('A', 2000) + Repeat("TTAGGG", 100); // human telomere at the 3′ end
        ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, "TTAGGG", minTelomereLength: 300)
            .Has3PrimeTelomere.Should().BeTrue("the matching motif detects it");
        ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, "TTTAGGG", minTelomereLength: 300)
            .TelomereLength3Prime.Should().BeLessThan(300, "a different motif does not register the array");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-CENT-001 — Centromere detection (Chromosome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 49.
    // Spec: tests/TestSpecs/CHROM-CENT-001.md (canonical AnalyzeCentromere).
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model (alpha-satellite centromeres): a centromere is a long tandem-satellite array of
    // high repeat content. AnalyzeCentromere slides a window, scoring repeatContent·(1−GC
    // variability), and calls a centromere where repeat content exceeds minAlphaSatelliteContent.
    //
    // The combinatorial point: window size, the satellite-content threshold and chromosome
    // length interact during the scan, yet a strong planted satellite array is detected and
    // localised to that array across all window sizes and lengths, for every threshold it clears.
    // ═══════════════════════════════════════════════════════════════════════

    private const int CentFlank = 600;

    [Test, Combinatorial]
    public void ChromCent_DetectsAndLocalisesSatelliteArray(
        [Values(200, 400, 800)] int windowSize,
        [Values(0.3, 0.5, 0.7)] double threshold,
        [Values(3000, 6000, 12000)] int seqLen)
    {
        int blockLen = seqLen - 2 * CentFlank;
        string unit = DiverseDna(21, 0xA1u);
        string block = Repeat(unit, blockLen / unit.Length + 1)[..blockLen];
        string seq = DiverseDna(CentFlank, 0x11u) + block + DiverseDna(CentFlank, 0x22u);

        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize, threshold);

        r.Start.Should().NotBeNull("a strong satellite array is a centromere");
        r.End.Should().NotBeNull();
        r.Length.Should().BeGreaterThan(0);
        r.AlphaSatelliteContent.Should().BeGreaterThan(0);

        r.Start!.Value.Should().BeGreaterThanOrEqualTo(0);
        r.End!.Value.Should().BeLessThanOrEqualTo(seqLen);
        (r.Start!.Value < seqLen - CentFlank && r.End!.Value > CentFlank)
            .Should().BeTrue("the detected region overlaps the planted satellite block");
    }

    /// <summary>
    /// Interaction witness: a non-repetitive chromosome has no centromere (repeat content is
    /// below any reasonable satellite threshold).
    /// </summary>
    [Test]
    public void ChromCent_UniqueSequence_NoCentromere()
    {
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", DiverseDna(4000, 0x33u), windowSize: 400,
            minAlphaSatelliteContent: 0.3);
        r.Start.Should().BeNull("unique sequence is not a centromere");
    }

    /// <summary>
    /// Interaction witness: an unreachable threshold suppresses detection even for a strong
    /// satellite array — the content gate is genuinely applied.
    /// </summary>
    [Test]
    public void ChromCent_ThresholdAboveContent_NoCentromere()
    {
        string unit = DiverseDna(21, 0xA1u);
        string seq = DiverseDna(CentFlank, 0x11u) + Repeat(unit, 100) + DiverseDna(CentFlank, 0x22u);
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize: 400, minAlphaSatelliteContent: 1.5);
        r.Start.Should().BeNull("no region exceeds an impossible threshold");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-ANEU-001 — Aneuploidy detection (Chromosome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 51.
    // Spec: tests/TestSpecs/CHROM-ANEU-001.md (canonical DetectAneuploidy +
    //       IdentifyWholeChromosomeAneuploidy).
    // Dimensions: nChrom(3) × depth(3) × threshold(3). Grid 3×3×3 = 27.
    //
    // Model (read-depth copy-number): bin copy number is round(2·depth/medianDepth) clamped
    // to [0,10] (Wikipedia: Aneuploidy). A whole chromosome is called aneuploid when one
    // non-disomic copy number dominates at least minFraction of its bins.
    //
    // The combinatorial point: chromosome count, the per-chromosome depth (copy level) and the
    // calling fraction interact — per-bin copy number follows the depth ratio exactly, and the
    // whole-chromosome call fires only when the dominant non-disomic level clears the threshold.
    // ═══════════════════════════════════════════════════════════════════════

    private const double MedianDepth = 30.0;
    private const int AneuBinSize = 1000;
    private const double AneuFraction = 0.7; // 7 of 10 bins carry the aneuploid level

    [Test, Combinatorial]
    public void ChromAneu_BinCopyNumberAndWholeChromCall(
        [Values(1, 2, 3)] int nChrom,
        [Values(1, 2, 3)] int copyLevel,   // 1=monosomy, 2=disomy, 3=trisomy
        [Values(0.5, 0.8, 1.0)] double threshold)
    {
        var depthData = new List<(string, int, double)>();
        for (int c = 0; c < nChrom; c++)
            for (int b = 0; b < 10; b++)
            {
                double depth = copyLevel != 2 && b < 7 ? copyLevel / 2.0 * MedianDepth : MedianDepth;
                depthData.Add(($"chr{c}", b * AneuBinSize, depth));
            }

        var states = ChromosomeAnalyzer.DetectAneuploidy(depthData, MedianDepth, AneuBinSize).ToList();

        states.Should().HaveCount(nChrom * 10);
        foreach (var s in states)
        {
            int bin = s.Start / AneuBinSize;
            int expectedCn = copyLevel != 2 && bin < 7 ? copyLevel : 2;
            s.CopyNumber.Should().Be(expectedCn, "bin CN follows round(2·depth/median)");
            s.Confidence.Should().BeInRange(0.0, 1.0);
        }

        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, threshold).ToList();
        bool shouldCall = copyLevel != 2 && AneuFraction >= threshold;
        calls.Should().HaveCount(shouldCall ? nChrom : 0,
            "a chromosome is called iff its dominant non-disomic level clears minFraction");
        if (shouldCall)
            calls.Should().OnlyContain(call => call.CopyNumber == copyLevel
                && call.Type == (copyLevel == 1 ? "Monosomy" : "Trisomy"));
    }

    /// <summary>
    /// Interaction witness: a uniform 1.5× chromosome is trisomy (CN 3 in every bin), called at
    /// any fraction; a uniform 1× (diploid-depth) chromosome is disomy and never called.
    /// </summary>
    [Test]
    public void ChromAneu_UniformTrisomy_Called_DisomyNot()
    {
        var trisomy = Enumerable.Range(0, 10).Select(b => ("chrT", b * AneuBinSize, 1.5 * MedianDepth)).Cast<(string, int, double)>();
        var trisomyStates = ChromosomeAnalyzer.DetectAneuploidy(trisomy, MedianDepth, AneuBinSize).ToList();
        trisomyStates.Should().OnlyContain(s => s.CopyNumber == 3);
        ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(trisomyStates, 0.8)
            .Should().ContainSingle(c => c.Type == "Trisomy");

        var disomy = Enumerable.Range(0, 10).Select(b => ("chrD", b * AneuBinSize, MedianDepth)).Cast<(string, int, double)>();
        var disomyStates = ChromosomeAnalyzer.DetectAneuploidy(disomy, MedianDepth, AneuBinSize).ToList();
        disomyStates.Should().OnlyContain(s => s.CopyNumber == 2);
        ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(disomyStates, 0.8).Should().BeEmpty();
    }

    /// <summary>
    /// Interaction witness: copy number is clamped to 10 for extreme depth, and empty / zero-median
    /// inputs yield no states (no division by zero).
    /// </summary>
    [Test]
    public void ChromAneu_Clamp_AndDegenerateInputs()
    {
        var deep = new[] { ("chr", 0, 1000.0) };
        ChromosomeAnalyzer.DetectAneuploidy(deep, MedianDepth, AneuBinSize).First().CopyNumber
            .Should().Be(10, "copy number is clamped to a maximum of 10");

        ChromosomeAnalyzer.DetectAneuploidy(Array.Empty<(string, int, double)>(), MedianDepth, AneuBinSize).Should().BeEmpty();
        ChromosomeAnalyzer.DetectAneuploidy(deep, 0.0, AneuBinSize).Should().BeEmpty("zero median is rejected");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-SYNT-001 — Synteny block detection (Chromosome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 52.
    // Spec: tests/TestSpecs/CHROM-SYNT-001.md (canonical FindSyntenyBlocks).
    // Dimensions: nGenes(3) × minBlockSize(3) × nChroms(3). Grid 3×3×3 = 27.
    //
    // Model (comparative genomics): a synteny block is a maximal collinear run of orthologous
    // genes between two genomes (same strand orientation, gaps within maxGap). FindSyntenyBlocks
    // groups ortholog pairs by chromosome pair and emits a block per collinear run of ≥ minGenes.
    //
    // The combinatorial point: gene count, the minimum-block size and the number of chromosome
    // pairs interact — a perfectly collinear run of G genes on each of nChroms chromosome pairs
    // yields nChroms blocks of size G exactly when G ≥ minGenes, and none otherwise.
    // ═══════════════════════════════════════════════════════════════════════

    private static List<(string, int, int, string, string, int, int, string)> Orthologs(int nChroms, int nGenes, bool forward = true)
    {
        var list = new List<(string, int, int, string, string, int, int, string)>();
        for (int c = 0; c < nChroms; c++)
            for (int g = 0; g < nGenes; g++)
            {
                int s1 = g * 1000;
                int s2 = (forward ? g : nGenes - 1 - g) * 1000;
                list.Add(($"A{c}", s1, s1 + 500, $"g1_{c}_{g}", $"B{c}", s2, s2 + 500, $"g2_{c}_{g}"));
            }
        return list;
    }

    [Test, Combinatorial]
    public void ChromSynt_CollinearRuns_PerChromosomePair(
        [Values(2, 3, 5)] int nGenes,
        [Values(2, 3, 4)] int minBlockSize,
        [Values(1, 2, 3)] int nChroms)
    {
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(Orthologs(nChroms, nGenes), minBlockSize, maxGap: 10).ToList();

        int expected = nGenes >= minBlockSize ? nChroms : 0;
        blocks.Should().HaveCount(expected, "one block per chromosome pair iff the run meets minGenes");

        if (expected > 0)
        {
            blocks.Should().OnlyContain(b => b.GeneCount == nGenes, "the whole collinear run is one block");
            blocks.Should().OnlyContain(b => b.Strand == '+', "increasing positions are forward synteny");
            blocks.Should().OnlyContain(b => b.Species1Start <= b.Species1End && b.Species2Start <= b.Species2End);
            blocks.Select(b => b.Species1Chromosome)
                .Should().BeEquivalentTo(Enumerable.Range(0, nChroms).Select(c => $"A{c}"));
        }
    }

    /// <summary>
    /// Interaction witness: a run whose orthologs run in the opposite order on the second genome
    /// is an inverted (reverse-strand) synteny block.
    /// </summary>
    [Test]
    public void ChromSynt_InvertedRun_IsReverseStrand()
    {
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(Orthologs(1, 5, forward: false), minGenes: 3, maxGap: 10).ToList();
        blocks.Should().ContainSingle();
        blocks[0].Strand.Should().Be('-', "decreasing second-genome order is an inversion");
        blocks[0].GeneCount.Should().Be(5);
    }

    /// <summary>
    /// Worked example: four collinear genes on one chromosome pair form a single forward block
    /// spanning the first gene's start to the last gene's end.
    /// </summary>
    [Test]
    public void ChromSynt_WorkedExample_SingleForwardBlock()
    {
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(Orthologs(1, 4), minGenes: 3, maxGap: 10).ToList();
        blocks.Should().ContainSingle();
        var b = blocks[0];
        b.GeneCount.Should().Be(4);
        b.Strand.Should().Be('+');
        b.Species1Chromosome.Should().Be("A0");
        b.Species1Start.Should().Be(0);
        b.Species1End.Should().Be(3 * 1000 + 500);
    }
}
