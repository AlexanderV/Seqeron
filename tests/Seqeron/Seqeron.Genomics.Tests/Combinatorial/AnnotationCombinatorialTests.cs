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
    // Unit: ANNOT-PROM-001 — Bacterial σ70 promoter boxes (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 30.
    // Spec: tests/TestSpecs/ANNOT-PROM-001.md (canonical FindPromoterMotifs → MotifFinder.FindSigma70Promoters).
    // Dimensions: plantedMismatches(3) × spacer(4) × mutatedBox(2). Grid 3×4×2 = 24.
    //
    // Model (Harley & Reynolds 1987): a σ70 promoter is a −35 hexamer ≈ TTGACA and a −10 hexamer ≈ TATAAT separated
    // by a 15–21 bp spacer; FindPromoterMotifs reports both boxes of each pair within the mismatch limits, score
    // 1 − d/6. The planted pair has one box mutated by plantedMismatches; the limit is 1 mismatch per box.
    //
    // The combinatorial point: the mismatch limit and the spacer window INTERACT — a pair is reported only when
    // BOTH hold (d ≤ 1 AND 15 ≤ spacer ≤ 21); failing either removes BOTH boxes (no unpaired box is reported).
    // ═══════════════════════════════════════════════════════════════════════

    public enum BoxType { Minus35, Minus10 }

    private static readonly string[] Mutated35 = { "TTGACA", "TTGACG", "TTGAGG" }; // d = 0, 1, 2
    private static readonly string[] Mutated10 = { "TATAAT", "TATAAG", "TATAGG" }; // d = 0, 1, 2

    [Test, Combinatorial]
    public void AnnotProm_MismatchLimitAndSpacer_GatePair(
        [Values(0, 1, 2)] int plantedMismatches,
        [Values(14, 15, 21, 22)] int spacer,
        [Values(BoxType.Minus35, BoxType.Minus10)] BoxType box)
    {
        string m35 = box == BoxType.Minus35 ? Mutated35[plantedMismatches] : "TTGACA";
        string m10 = box == BoxType.Minus10 ? Mutated10[plantedMismatches] : "TATAAT";
        string seq = "GGGG" + m35 + new string('G', spacer) + m10 + "GGGG";

        var hits = GenomeAnnotator.FindPromoterMotifs(seq, 1, 1, 15, 21).ToList();

        bool expected = plantedMismatches <= 1 && spacer is >= 15 and <= 21;
        if (expected)
        {
            double s35 = box == BoxType.Minus35 ? (6 - plantedMismatches) / 6.0 : 1.0;
            double s10 = box == BoxType.Minus10 ? (6 - plantedMismatches) / 6.0 : 1.0;
            hits.Should().Equal((4, "-35 box", m35, s35), (10 + spacer, "-10 box", m10, s10));
        }
        else
        {
            hits.Should().BeEmpty("the pair fails the mismatch limit or the 15–21 bp spacer window");
        }
    }

    /// <summary>
    /// Interaction witness: the same 2-mismatch −35 box is reported under the default limit (2) and dropped under
    /// limit 1, while a spacer of 22 drops the consensus pair under either limit.
    /// </summary>
    [Test]
    public void AnnotProm_LimitAndSpacer_Interact()
    {
        string twoMm = "GGGG" + "TTGAGG" + new string('G', 17) + "TATAAT" + "GGGG";
        GenomeAnnotator.FindPromoterMotifs(twoMm).Should().HaveCount(2);
        GenomeAnnotator.FindPromoterMotifs(twoMm, 1, 1, 15, 21).Should().BeEmpty();

        string far = "GGGG" + "TTGACA" + new string('G', 22) + "TATAAT" + "GGGG";
        GenomeAnnotator.FindPromoterMotifs(far).Should().BeEmpty();
        GenomeAnnotator.FindPromoterMotifs(far, 1, 1, 15, 21).Should().BeEmpty();
    }

    /// <summary>
    /// Worked example: a canonical E. coli promoter (−35 box, 17-bp spacer, −10 box) yields exactly both consensus
    /// boxes at score 1.0.
    /// </summary>
    [Test]
    public void AnnotProm_CanonicalPromoter_WorkedExample()
    {
        string promoter = new string('A', 5) + "TTGACA" + new string('G', 17) + "TATAAT" + new string('A', 5);
        var hits = GenomeAnnotator.FindPromoterMotifs(promoter).ToList();

        hits.Should().Contain((5, "-35 box", "TTGACA", 1.0));
        hits.Should().Contain((28, "-10 box", "TATAAT", 1.0));
        hits.Should().OnlyContain(h => h.score >= 4.0 / 6);
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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-CODING-001 — Coding-potential score (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 216.
    // Spec: tests/TestSpecs/ANNOT-CODING-001.md (canonical GenomeAnnotator.CalculateCodingPotential).
    // ADVANCED §10.
    // Dimensions: method(2) × seqLen(3). Grid 2×3 = 6 (full, exhaustive ⊇ pairwise).
    //
    // Model (CPAT, Wang 2013): the coding potential is the mean log-ratio of coding to non-coding
    // hexamer frequencies over the frame-0 hexamers — positive for coding-like sequence, negative for
    // non-coding-like sequence.
    //
    // Axis mapping (documented): method → the sequence regime {CodingLike, NoncodingLike}; seqLen →
    // length. Engineered hexamer tables make AAAAAA coding-favoured and TTTTTT non-coding-favoured.
    // The combinatorial point: the score is the log-ratio mean (here ±ln 9), with the correct sign for
    // each regime at every length.
    // ═══════════════════════════════════════════════════════════════════════

    public enum CodingRegime { CodingLike, NoncodingLike }

    [Test, Combinatorial]
    public void AnnotCoding_LogRatioSign_AcrossRegimeAndLength(
        [Values(CodingRegime.CodingLike, CodingRegime.NoncodingLike)] CodingRegime method,
        [Values(12, 30, 60)] int seqLen)
    {
        var coding = new Dictionary<string, double> { ["AAAAAA"] = 0.9, ["TTTTTT"] = 0.1 };
        var noncoding = new Dictionary<string, double> { ["AAAAAA"] = 0.1, ["TTTTTT"] = 0.9 };
        string seq = new string(method == CodingRegime.CodingLike ? 'A' : 'T', seqLen);

        double score = GenomeAnnotator.CalculateCodingPotential(seq, coding, noncoding);
        double expected = method == CodingRegime.CodingLike ? Math.Log(0.9 / 0.1) : Math.Log(0.1 / 0.9);

        score.Should().BeApproximately(expected, 1e-9, "coding potential = mean log(coding/noncoding) hexamer ratio");
        if (method == CodingRegime.CodingLike) score.Should().BeGreaterThan(0);
        else score.Should().BeLessThan(0);
    }

    /// <summary>
    /// Interaction witness — a coding-favoured sequence scores above a non-coding-favoured one, and a
    /// sequence shorter than one hexamer scores 0.
    /// </summary>
    [Test]
    public void AnnotCoding_CodingScoresAboveNoncoding()
    {
        var coding = new Dictionary<string, double> { ["AAAAAA"] = 0.9, ["TTTTTT"] = 0.1 };
        var noncoding = new Dictionary<string, double> { ["AAAAAA"] = 0.1, ["TTTTTT"] = 0.9 };

        GenomeAnnotator.CalculateCodingPotential(new string('A', 30), coding, noncoding)
            .Should().BeGreaterThan(GenomeAnnotator.CalculateCodingPotential(new string('T', 30), coding, noncoding));
        GenomeAnnotator.CalculateCodingPotential("AAA", coding, noncoding).Should().Be(0, "shorter than a hexamer ⇒ 0");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-REPEAT-001 — Repetitive-element detection (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 218.
    // Spec: tests/TestSpecs/ANNOT-REPEAT-001.md (canonical FindRepetitiveElements / ClassifyRepeat).
    // ADVANCED §10.
    // Dimensions: minLen(3) × repeatType(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (RepeatMasker; Wikipedia tandem/inverted repeats): FindRepetitiveElements reports tandem
    // repeats (head-to-tail copies of a motif) and inverted repeats (W·G·W̄ᴿ) whose total span ≥
    // minRepeatLength; ClassifyRepeat assigns a RepeatMasker-style class.
    //
    // Axis mapping (documented): minLen → minRepeatLength; repeatType → the planted element kind
    // {Tandem, Inverted, Both}. Engineered constructs span 24 bp. The combinatorial point: the planted
    // element is reported with its correct type exactly when its span (24) ≥ minLen, and every reported
    // element meets the length floor.
    // ═══════════════════════════════════════════════════════════════════════

    public enum RepeatKind { Tandem, Inverted, Both }

    [Test, Combinatorial]
    public void AnnotRepeat_DetectsPlantedElement_AcrossMinLenAndType(
        [Values(10, 20, 40)] int minLen,
        [Values(RepeatKind.Tandem, RepeatKind.Inverted, RepeatKind.Both)] RepeatKind repeatType)
    {
        const string tandem = "ACGTAC" + "ACGTAC" + "ACGTAC" + "ACGTAC"; // unit ACGTAC ×4 = 24 bp
        string arm = "GATTACAGGC"; // 10 bp
        string inverted = arm + "TT" + RevComp(arm); // W·gap·W̄ᴿ, span 22

        string seq = repeatType switch
        {
            RepeatKind.Tandem => "CC" + tandem + "CC",
            RepeatKind.Inverted => "CC" + inverted + "CC",
            _ => "CC" + tandem + "CCGGCC" + inverted + "CC",
        };

        var elements = GenomeAnnotator.FindRepetitiveElements(seq, minLen, minCopies: 2).ToList();
        elements.Should().OnlyContain(e => e.end - e.start + 1 >= minLen, "every element meets the length floor");
        elements.Should().OnlyContain(e => e.type == "tandem_repeat" || e.type == "inverted_repeat");

        // Tandem gate is on total span (24 bp); inverted gate is on the arm length (10 bp).
        bool tandemFits = 24 >= minLen, invertedFits = arm.Length >= minLen;
        if (repeatType is RepeatKind.Tandem or RepeatKind.Both)
            elements.Any(e => e.type == "tandem_repeat").Should().Be(tandemFits, "tandem found iff its 24-bp span ≥ minLen");
        if (repeatType is RepeatKind.Inverted or RepeatKind.Both)
            elements.Any(e => e.type == "inverted_repeat").Should().Be(invertedFits, "inverted found iff its arm length ≥ minLen");
    }

    /// <summary>
    /// Interaction witness — ClassifyRepeat: a library element wins by name, a short microsatellite is
    /// a Simple_repeat, and an unmatched non-repeat is Unknown.
    /// </summary>
    [Test]
    public void AnnotRepeat_ClassifyByLibraryThenSimpleThenUnknown()
    {
        var db = new Dictionary<string, string> { ["GGGCCCAAATTT"] = "LINE/L1" };

        GenomeAnnotator.ClassifyRepeat("AAGGGCCCAAATTTAA", db).Should().Be("LINE/L1", "a library element match wins");
        GenomeAnnotator.ClassifyRepeat("ATATATATAT", db).Should().Be("Simple_repeat", "a 2-bp microsatellite is a simple repeat");
        GenomeAnnotator.ClassifyRepeat("GATTACAGCTAG", db).Should().Be("Unknown", "no library match and not a simple repeat");
    }
}
