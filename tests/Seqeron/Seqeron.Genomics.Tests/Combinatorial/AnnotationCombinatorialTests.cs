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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-PROM-001 — Bacterial promoter-motif detection (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 30.
    // Spec: tests/TestSpecs/ANNOT-PROM-001.md (canonical FindPromoterMotifs).
    // Dimensions: threshold(3) × windowSize(3) × motifSet(2). Grid 3×3×2 = 18.
    //
    // Model (Pribnow 1975; Harley & Reynolds 1987): FindPromoterMotifs scans for the −35
    // box (TTGACA) and −10/Pribnow box (TATAAT) and their prefix/suffix variants, scoring
    // each by summed E. coli position-occurrence probabilities normalised to the consensus
    // (e.g. −35 TTGAC=0.855, TGACA=0.815, TTGA=0.710; −10 ATAAT=0.813, TATAA=0.801,
    // TATA=0.665). The method takes only a sequence, so the checklist axes map to caller
    // operations: motifSet = box type filter, threshold = a score cutoff over the variant
    // table, windowSize = the length of the searched substring.
    //
    // The combinatorial point: threshold and motifSet INTERACT because the two boxes have
    // different variant-score spectra — at cutoff 0.84 the −35 box admits two variants
    // (1.000, 0.855) but the −10 box admits only the consensus (its next variant is 0.813).
    // windowSize gates whether the planted box lies inside the searched region at all.
    // ═══════════════════════════════════════════════════════════════════════

    public enum BoxType { Minus35, Minus10 }

    private const int PromBoxStart = 10; // box occupies [10,16) in 'A' filler (no spurious motifs)

    /// <summary>Distinct variant count of a full-consensus box passing the score cutoff (from the documented table).</summary>
    private static int ExpectedPassingVariants(BoxType box, double threshold) => box == BoxType.Minus35
        ? (threshold > 0.99 ? 1 : threshold > 0.855 ? 1 : threshold > 0.815 ? 2 : threshold > 0.71 ? 3 : 4)
        : (threshold > 0.99 ? 1 : threshold > 0.813 ? 1 : threshold > 0.801 ? 2 : threshold > 0.665 ? 3 : 4);

    [Test, Combinatorial]
    public void AnnotProm_ScoreThresholdAndWindow_GateBoxVariants(
        [Values(0.99, 0.84, 0.70)] double threshold,
        [Values(8, 16, 30)] int windowSize,
        [Values(BoxType.Minus35, BoxType.Minus10)] BoxType box)
    {
        string consensus = box == BoxType.Minus35 ? "TTGACA" : "TATAAT";
        string typeName = box == BoxType.Minus35 ? "-35 box" : "-10 box";
        string template = new string('A', PromBoxStart) + consensus + new string('A', 20);
        string searched = template[..windowSize];

        var hits = GenomeAnnotator.FindPromoterMotifs(searched).Where(h => h.type == typeName).ToList();
        var passing = hits.Where(h => h.score >= threshold).Select(h => h.sequence).Distinct().ToList();

        bool boxInWindow = windowSize >= PromBoxStart + consensus.Length;
        if (boxInWindow)
        {
            hits.Should().Contain(h => h.sequence == consensus && h.score == 1.0, "the full consensus scores 1.0");
            passing.Should().HaveCount(ExpectedPassingVariants(box, threshold),
                $"{typeName} variants passing cutoff {threshold}");
        }
        else
        {
            hits.Should().BeEmpty("the box lies outside the searched window");
        }
    }

    /// <summary>
    /// Interaction witness: at the same score cutoff (0.84) the two boxes admit different
    /// numbers of variants — the −35 spectrum has a 0.855 variant above the cutoff while the
    /// −10 spectrum's next variant (0.813) falls below it.
    /// </summary>
    [Test]
    public void AnnotProm_ThresholdAndBoxType_Interact()
    {
        string m35 = new string('A', 10) + "TTGACA" + new string('A', 10);
        string m10 = new string('A', 10) + "TATAAT" + new string('A', 10);

        int n35 = GenomeAnnotator.FindPromoterMotifs(m35)
            .Where(h => h.type == "-35 box" && h.score >= 0.84).Select(h => h.sequence).Distinct().Count();
        int n10 = GenomeAnnotator.FindPromoterMotifs(m10)
            .Where(h => h.type == "-10 box" && h.score >= 0.84).Select(h => h.sequence).Distinct().Count();

        n35.Should().Be(2, "TTGACA(1.000) and TTGAC(0.855) clear 0.84");
        n10.Should().Be(1, "only TATAAT(1.000) clears 0.84; ATAAT(0.813) does not");
    }

    /// <summary>
    /// Worked example: a canonical E. coli promoter (−35 box, 17-bp spacer, −10 box) yields
    /// both consensus boxes with the right type, position and score 1.0, and the documented
    /// variant scores match the source table.
    /// </summary>
    [Test]
    public void AnnotProm_CanonicalPromoter_WorkedExample()
    {
        string promoter = new string('A', 5) + "TTGACA" + new string('G', 17) + "TATAAT" + new string('A', 5);
        var hits = GenomeAnnotator.FindPromoterMotifs(promoter).ToList();

        hits.Should().Contain(h => h.type == "-35 box" && h.sequence == "TTGACA" && h.position == 5 && h.score == 1.0);
        hits.Should().Contain(h => h.type == "-10 box" && h.sequence == "TATAAT" && h.position == 28 && h.score == 1.0);

        hits.First(h => h.sequence == "TTGAC").score.Should().BeApproximately(0.855, 1e-9);
        hits.First(h => h.sequence == "ATAAT").score.Should().BeApproximately(0.813, 1e-9);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-GFF-001 — GFF3 feature I/O (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 31.
    // Spec: tests/TestSpecs/ANNOT-GFF-001.md (canonical ParseGff3 / ToGff3).
    // Dimensions: featureType(3) × strand(3) × phase(3). Grid 3×3×3 = 27.
    //
    // Model (Sequence Ontology GFF3 v1.26): a feature line is nine tab-separated columns;
    // column 3 is the SO type, column 7 the strand (+/−/.), column 8 the phase (0/1/2 or .).
    // ParseGff3 maps these to GenomicFeature.Type/Strand/Phase, defaulting phase "." to null.
    //
    // The combinatorial point: the three columns are parsed independently and must not
    // bleed into one another — every (type, strand, phase) triple must round-trip into the
    // corresponding GenomicFeature fields with the other columns (coordinates, id) intact.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void AnnotGff_ParsesTypeStrandPhase_Independently(
        [Values("gene", "mRNA", "CDS")] string featureType,
        [Values("+", "-", ".")] string strand,
        [Values("0", "1", "2")] string phase)
    {
        string line = $"chr1\t.\t{featureType}\t100\t200\t.\t{strand}\t{phase}\tID=f1;Name=demo";

        var features = GenomeAnnotator.ParseGff3(new[] { line }).ToList();

        features.Should().ContainSingle();
        var f = features[0];
        f.Type.Should().Be(featureType, "column 3 is the feature type");
        f.Strand.Should().Be(strand[0], "column 7 is the strand");
        f.Phase.Should().Be(int.Parse(phase), "column 8 is the phase");
        f.Start.Should().Be(100);
        f.End.Should().Be(200);
        f.Score.Should().BeNull("a '.' score column parses to null");
        f.FeatureId.Should().Be("f1", "the ID attribute drives the feature id");
        f.Attributes["Name"].Should().Be("demo");
    }

    /// <summary>
    /// Interaction witness: phase is REQUIRED for CDS — ToGff3 emits phase "0" for a CDS but
    /// "." for any other type, so a round-trip yields Phase 0 only for the CDS feature.
    /// — GFF3 v1.26 NOTE 4.
    /// </summary>
    [Test]
    public void AnnotGff_PhaseColumn_DependsOnFeatureType()
    {
        var cds = new GenomeAnnotator.GeneAnnotation("g1", 99, 200, '+', "CDS", "prot",
            new Dictionary<string, string>());

        string cdsLine = GenomeAnnotator.ToGff3(new[] { cds }, "chr1").Skip(1).First();
        var parsedCds = GenomeAnnotator.ParseGff3(new[] { cdsLine }).Single();
        parsedCds.Phase.Should().Be(0, "CDS carries phase 0");
        parsedCds.Start.Should().Be(100, "ToGff3 writes 1-based Start (ann.Start + 1)");
        parsedCds.Strand.Should().Be('+');

        var gene = cds with { Type = "gene" };
        string geneLine = GenomeAnnotator.ToGff3(new[] { gene }, "chr1").Skip(1).First();
        GenomeAnnotator.ParseGff3(new[] { geneLine }).Single().Phase.Should().BeNull("non-CDS phase is '.'");
    }

    /// <summary>
    /// Worked example: reserved attribute characters (';', '=') survive a ToGff3 → ParseGff3
    /// round-trip via GFF3 percent-encoding/decoding. — GFF3 v1.26 column-9 encoding.
    /// </summary>
    [Test]
    public void AnnotGff_AttributeEncoding_RoundTrips()
    {
        var ann = new GenomeAnnotator.GeneAnnotation("id;x=y", 0, 9, '+', "gene", "a=b;c",
            new Dictionary<string, string>());

        string line = GenomeAnnotator.ToGff3(new[] { ann }, "chr1").Skip(1).First();
        line.Should().NotContain("id;x=y", "reserved characters are percent-encoded on output");

        var f = GenomeAnnotator.ParseGff3(new[] { line }).Single();
        f.FeatureId.Should().Be("id;x=y", "ID is decoded on input");
        f.Attributes["product"].Should().Be("a=b;c");
    }
}
