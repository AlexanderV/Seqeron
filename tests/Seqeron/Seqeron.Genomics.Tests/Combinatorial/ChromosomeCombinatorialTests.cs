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
}
