namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Annotation area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Annotation")]
public class AnnotationCombinatorialTests
{
    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-ORF-001 — Open-reading-frame detection (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 28.
    // Spec: tests/TestSpecs/ANNOT-ORF-001.md (canonical GenomeAnnotator.FindOrfs).
    // Dimensions: minLen(3) × frame(3: fwd/rev/both) × startCodon(2: ATG/alt). Grid 3×3×2 = 18.
    //
    // Model (Wikipedia ORF; Rosalind): an ORF is a maximal run start→stop in one of the
    // six frames; FindOrfs reports its amino-acid length as the codons from the start up
    // to (excluding) the stop, filters by minLength, recognises the canonical ATG plus the
    // alternative bacterial starts GTG/TTG, and — with requireStartCodon — emits only
    // runs that begin at a start codon AND terminate at a stop. searchBothStrands governs
    // whether the reverse complement is scanned; reverse hits carry IsReverseComplement.
    //
    // The combinatorial point: minLength, strand and start codon interact. A length-9 ORF
    // is detected for minLength ≤ 9 and filtered above it; it is found on whichever strand
    // it sits (with the correct IsReverseComplement flag); and both ATG- and GTG-started
    // ORFs are recognised. The grid embeds one clean ORF (no internal start/stop) in
    // start/stop-free poly-C filler so detection is exact and unambiguous.
    // ═══════════════════════════════════════════════════════════════════════

    public enum Orientation { Forward, Reverse, Both }

    private const int OrfAaLength = 9;                       // ATG/GTG + 8×AAA, stop excluded
    private static string MakeOrf(string start) => start + string.Concat(Enumerable.Repeat("AAA", 8)) + "TAA";
    private static readonly string Filler = new('C', 12);   // no start/stop codon on either strand

    [Test, Combinatorial]
    public void AnnotOrf_Detection_AcrossMinLenStrandAndStartCodon(
        [Values(5, 9, 15)] int minLength,
        [Values(Orientation.Forward, Orientation.Reverse, Orientation.Both)] Orientation orientation,
        [Values("ATG", "GTG")] string startCodon)
    {
        string orf = MakeOrf(startCodon);
        string template = orientation switch
        {
            Orientation.Forward => Filler + orf + Filler,
            Orientation.Reverse => Filler + RevComp(orf) + Filler,
            _ => Filler + orf + Filler + RevComp(orf) + Filler,
        };

        var orfs = GenomeAnnotator
            .FindOrfs(template, minLength, searchBothStrands: true, requireStartCodon: true)
            .ToList();

        bool detected = minLength <= OrfAaLength;
        var forward = orfs.Where(o => o.Sequence == orf && !o.IsReverseComplement).ToList();
        var reverse = orfs.Where(o => o.Sequence == orf && o.IsReverseComplement).ToList();

        int expFwd = orientation != Orientation.Reverse && detected ? 1 : 0;
        int expRev = orientation != Orientation.Forward && detected ? 1 : 0;
        forward.Should().HaveCount(expFwd, $"forward ORF presence under {orientation}, minLength {minLength}");
        reverse.Should().HaveCount(expRev, $"reverse ORF presence under {orientation}, minLength {minLength}");

        // Every reported match obeys the ORF invariants (start/stop/in-frame).
        foreach (var o in forward.Concat(reverse))
        {
            o.Sequence[..3].Should().BeOneOf("ATG", "GTG", "TTG");
            o.Sequence[^3..].Should().BeOneOf("TAA", "TAG", "TGA");
            (o.Sequence.Length % 3).Should().Be(0);
        }
    }

    /// <summary>
    /// Interaction witness: searchBothStrands is the strand gate. A reverse-strand ORF is
    /// invisible to a forward-only search but found once both strands are scanned.
    /// </summary>
    [Test]
    public void AnnotOrf_ForwardOnlySearch_IgnoresReverseStrandOrf()
    {
        string orf = MakeOrf("ATG");
        string template = Filler + RevComp(orf) + Filler; // ORF only on the reverse strand

        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeFalse("forward-only search cannot see a reverse ORF");
        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: true, requireStartCodon: true)
            .Any(o => o.Sequence == orf && o.IsReverseComplement).Should().BeTrue("both-strand search finds it");
    }

    /// <summary>
    /// Interaction witness: minLength is an inclusive amino-acid threshold — an ORF exactly
    /// at the threshold is kept, one amino acid longer than the threshold-minus-one is the
    /// boundary. — spec M06/M06b.
    /// </summary>
    [Test]
    public void AnnotOrf_MinLength_IsInclusiveThreshold()
    {
        string orf = MakeOrf("ATG"); // 9 aa
        string template = Filler + orf + Filler;

        GenomeAnnotator.FindOrfs(template, minLength: OrfAaLength, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeTrue("length == minLength is included");
        GenomeAnnotator.FindOrfs(template, minLength: OrfAaLength + 1, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeFalse("length < minLength is excluded");
    }

    /// <summary>
    /// Interaction witness: the alternative start TTG is also recognised, and requiring a
    /// start codon suppresses a start-less sense run that ends in a stop.
    /// </summary>
    [Test]
    public void AnnotOrf_StartCodonPolicy_AltStartRecognised_RequireStartFilters()
    {
        string ttgOrf = MakeOrf("TTG");
        string template = Filler + ttgOrf + Filler;
        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == ttgOrf).Should().BeTrue("TTG is an alternative start codon");

        // A frame with no start codon: AAA…AAA TAA. requireStartCodon must yield nothing.
        string startless = string.Concat(Enumerable.Repeat("AAA", 9)) + "TAA";
        GenomeAnnotator.FindOrfs(startless, minLength: 1, searchBothStrands: false, requireStartCodon: true)
            .Should().BeEmpty("with requireStartCodon a run lacking a start codon is not an ORF");
    }

    /// <summary>
    /// Worked example: the textbook minimal ORF ATG-AAA-TAA is one frame-1 forward ORF
    /// spanning [0,9) translating to M-K (stop trimmed from the protein).
    /// </summary>
    [Test]
    public void AnnotOrf_MinimalOrf_WorkedExample()
    {
        var orfs = GenomeAnnotator.FindOrfs("ATGAAATAA", minLength: 1, searchBothStrands: false, requireStartCodon: true).ToList();

        orfs.Should().ContainSingle();
        var o = orfs[0];
        o.Start.Should().Be(0);
        o.End.Should().Be(9);
        o.Frame.Should().Be(1);
        o.IsReverseComplement.Should().BeFalse();
        o.Sequence.Should().Be("ATGAAATAA");
        o.ProteinSequence.Should().StartWith("MK");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-GENE-001 — Gene prediction with ribosome-binding sites (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 29.
    // Spec: tests/TestSpecs/ANNOT-GENE-001.md (canonical PredictGenes +
    //       FindRibosomeBindingSites).
    // Dimensions: minOrfLen(3) × rbsWindow(3) × scoring(2). Grid 3×3×2 = 18.
    //
    // Model (prokaryotic ORF-based gene finding; Shine & Dalgarno 1975): PredictGenes
    // turns every qualifying ORF (≥ minOrfLength aa, both strands) into a CDS gene with a
    // sequential id and a protein_length attribute = (End−Start)/3 − 1 (stop excluded).
    // FindRibosomeBindingSites locates the Shine-Dalgarno motif (AGGAGG…) upstream of a
    // forward ORF start at an aligned spacing within [minDistance, maxDistance] (Chen 1994:
    // functional 4–15 nt); the both-strand variant additionally reports reverse-strand SDs.
    //
    // The combinatorial point: minOrfLen gates gene calling, the SD spacing (rbsWindow)
    // gates RBS detection, and the scan mode (forward-only vs both-strand "scoring") agrees
    // on forward hits — three orthogonal knobs each asserted per cell. The construct plants
    // one 40-aa forward gene preceded by AGGAGG at a chosen spacing in start/stop-free filler.
    // ═══════════════════════════════════════════════════════════════════════

    public enum ScanMode { ForwardOnly, BothStrands }

    private const int GeneAaLength = 40;                       // ATG + 39×AAA, stop excluded
    private const int GeneNtSpan = (GeneAaLength + 1) * 3;     // ORF DNA incl. stop = 123
    private const int SdPosition = 12;                         // AGGAGG starts after 12 nt of filler
    private static string MakeGeneOrf() => "ATG" + string.Concat(Enumerable.Repeat("AAA", 39)) + "TAA";

    private static string BuildGeneTemplate(int sdSpacing) =>
        new string('C', SdPosition) + "AGGAGG" + new string('C', sdSpacing) + MakeGeneOrf() + new string('C', 12);

    [Test, Combinatorial]
    public void AnnotGene_PredictionAndRbs_AcrossMinOrfLenSpacingAndScanMode(
        [Values(20, 40, 50)] int minOrfLen,
        [Values(6, 15, 30)] int sdSpacing,
        [Values(ScanMode.ForwardOnly, ScanMode.BothStrands)] ScanMode scanMode)
    {
        string template = BuildGeneTemplate(sdSpacing);

        // ── Gene prediction depends only on minOrfLen vs the 40-aa gene. ──
        var genes = GenomeAnnotator.PredictGenes(template, minOrfLen, "gene").ToList();
        if (minOrfLen <= GeneAaLength)
        {
            genes.Should().ContainSingle("the 40-aa ORF is called iff minOrfLen ≤ 40");
            var g = genes[0];
            g.Type.Should().Be("CDS");
            g.Strand.Should().Be('+');
            g.Start.Should().BeLessThan(g.End);
            (g.End - g.Start).Should().Be(GeneNtSpan);
            g.GeneId.Should().MatchRegex(@"^gene_\d{4}$");
            int proteinLength = int.Parse(g.Attributes["protein_length"]);
            proteinLength.Should().Be((g.End - g.Start) / 3 - 1, "protein_length excludes the stop codon");
            proteinLength.Should().Be(GeneAaLength);
        }
        else
        {
            genes.Should().BeEmpty("a 40-aa ORF is filtered when minOrfLen = 50");
        }

        // ── RBS detection depends only on the SD spacing vs [4,15]; both scan modes agree on '+'. ──
        IEnumerable<int> rbsPositions = scanMode == ScanMode.ForwardOnly
            ? GenomeAnnotator.FindRibosomeBindingSites(template, upstreamWindow: 40, minDistance: 4, maxDistance: 15)
                .Select(h => h.position)
            : GenomeAnnotator.FindRibosomeBindingSitesBothStrands(template, upstreamWindow: 40, minDistance: 4, maxDistance: 15)
                .Where(h => h.strand == '+').Select(h => h.position);

        bool sdDetected = sdSpacing is >= 4 and <= 15;
        rbsPositions.Contains(SdPosition).Should().Be(sdDetected,
            $"AGGAGG at spacing {sdSpacing} is detected iff 4 ≤ spacing ≤ 15");
    }

    /// <summary>
    /// Interaction witness: the scan-mode "scoring" axis is real — a reverse-strand
    /// Shine-Dalgarno site is reported by the both-strand scan (strand '-') but is invisible
    /// to the forward-only scan.
    /// </summary>
    [Test]
    public void AnnotGene_BothStrandScan_FindsReverseStrandRbs_ForwardOnlyMissesIt()
    {
        string forwardConstruct = "AGGAGG" + new string('C', 6) + MakeGeneOrf();
        string template = new string('C', 12) + RevComp(forwardConstruct) + new string('C', 12);

        var both = GenomeAnnotator
            .FindRibosomeBindingSitesBothStrands(template, upstreamWindow: 40, minDistance: 4, maxDistance: 15)
            .ToList();
        both.Should().Contain(h => h.strand == '-', "the reverse-strand SD is found by the both-strand scan");

        int reverseSd = both.First(h => h.strand == '-').position;
        var forward = GenomeAnnotator
            .FindRibosomeBindingSites(template, upstreamWindow: 40, minDistance: 4, maxDistance: 15)
            .Select(h => h.position).ToHashSet();
        forward.Should().NotContain(reverseSd, "the forward-only scan cannot see a reverse-strand SD");
    }

    /// <summary>
    /// Worked example: two stacked forward ORFs are called as gene_0001 / gene_0002 with
    /// ascending ids and CDS type. — spec gene-id pattern invariant.
    /// </summary>
    [Test]
    public void AnnotGene_MultipleGenes_NumberedSequentially()
    {
        string two = new string('C', 9) + MakeGeneOrf() + new string('C', 9) + MakeGeneOrf() + new string('C', 9);
        var genes = GenomeAnnotator.PredictGenes(two, minOrfLength: 20, prefix: "gene")
            .Where(g => g.Strand == '+').OrderBy(g => g.Start).ToList();

        genes.Should().HaveCountGreaterThanOrEqualTo(2);
        genes[0].GeneId.Should().Be("gene_0001");
        genes[1].GeneId.Should().Be("gene_0002");
        genes.Should().OnlyContain(g => g.Type == "CDS");
    }
}
