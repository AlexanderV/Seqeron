namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Splicing area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Splicing")]
public class SplicingCombinatorialTests
{
    private const int SiteOffset = 30;

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-DONOR-001 — 5′ splice-site (donor) prediction (Splicing)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 77.
    // Spec: tests/TestSpecs/SPLICE-DONOR-001.md (canonical FindDonorSites).
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model (Shapiro & Senapathy 1987): the U2 5′ splice site is the consensus MAG|GURAGU with an
    // invariant GU at the intron start; FindDonorSites scans for GU and scores the consensus window.
    // The scoring window is fixed internally, so windowSize maps to the searched-substring length;
    // threshold is the minimum score, seqLen the total length.
    //
    // The combinatorial point: search window, score threshold and length interact — a planted
    // consensus donor is reported exactly when it lies inside the searched window AND its score
    // clears the threshold.
    // ═══════════════════════════════════════════════════════════════════════

    private static string DonorSequence(int seqLen) =>
        new string('A', SiteOffset) + "GTAAGT" + new string('A', seqLen - SiteOffset - 6); // GU consensus at offset 30

    [Test, Combinatorial]
    public void SpliceDonor_DetectionGatedByWindowAndThreshold(
        [Values(20, 40, 80)] int windowSize,
        [Values(0.3, 0.6, 0.9)] double threshold,
        [Values(90, 130, 200)] int seqLen)
    {
        string seq = DonorSequence(seqLen);
        double plantedScore = SpliceSitePredictor.FindDonorSites(seq, 0.0)
            .First(s => s.Position == SiteOffset).Score;

        var hits = SpliceSitePredictor.FindDonorSites(seq[..windowSize], threshold).ToList();

        bool expected = windowSize >= SiteOffset + 6 && plantedScore >= threshold;
        hits.Any(s => s.Position == SiteOffset && s.Type == SpliceSitePredictor.SpliceSiteType.Donor)
            .Should().Be(expected, "the donor is reported iff within the window and above the threshold");

        hits.Should().OnlyContain(s => s.Score >= threshold, "every reported site clears the threshold");
    }

    /// <summary>
    /// Interaction witness: a strong consensus donor (GTAAGT) scores higher than a degenerate one
    /// (GTxxxx), so raising the threshold drops the weak site first.
    /// </summary>
    [Test]
    public void SpliceDonor_ConsensusScoresHigherThanDegenerate()
    {
        double consensus = SpliceSitePredictor.FindDonorSites(new string('A', 30) + "GTAAGT" + new string('A', 30), 0.0)
            .First(s => s.Position == 30).Score;
        double degenerate = SpliceSitePredictor.FindDonorSites(new string('A', 30) + "GTCCCC" + new string('A', 30), 0.0)
            .First(s => s.Position == 30).Score;
        consensus.Should().BeGreaterThan(degenerate, "the GURAGU consensus scores above a non-consensus GU");
    }

    /// <summary>
    /// Interaction witness: every reported donor begins with the invariant GU dinucleotide, and a
    /// GU-free sequence yields none.
    /// </summary>
    [Test]
    public void SpliceDonor_RequiresGuDinucleotide()
    {
        var hits = SpliceSitePredictor.FindDonorSites(DonorSequence(90), 0.0).ToList();
        hits.Should().NotBeEmpty();
        hits.Should().OnlyContain(s => s.Motif.Contains("GU"), "donor motifs contain the GU intron start");

        SpliceSitePredictor.FindDonorSites(new string('A', 60), 0.0).Should().BeEmpty("no GU ⇒ no donor");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-ACCEPTOR-001 — 3′ splice-site (acceptor) prediction (Splicing)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 78.
    // Spec: tests/TestSpecs/SPLICE-ACCEPTOR-001.md (canonical FindAcceptorSites).
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model (Shapiro & Senapathy 1987): the 3′ splice site is a polypyrimidine tract followed by
    // the invariant AG at the intron end; FindAcceptorSites scans AG and scores the upstream
    // pyrimidine context. As with the donor, windowSize maps to the searched-substring length,
    // threshold is the minimum score, seqLen the total length.
    //
    // The combinatorial point: search window, threshold and length interact — a planted acceptor
    // (poly-C tract + AG) is reported exactly when its AG is scannable within the window and its
    // score clears the threshold.
    // ═══════════════════════════════════════════════════════════════════════

    private static string AcceptorSequence(int seqLen) =>
        new string('C', 30) + "AG" + new string('G', seqLen - 32); // polypyrimidine tract + AG, A at offset 30, G at 31

    [Test, Combinatorial]
    public void SpliceAcceptor_DetectionGatedByWindowAndThreshold(
        [Values(25, 40, 80)] int windowSize,
        [Values(0.3, 0.6, 0.9)] double threshold,
        [Values(90, 130, 200)] int seqLen)
    {
        string seq = AcceptorSequence(seqLen);
        double plantedScore = SpliceSitePredictor.FindAcceptorSites(seq, 0.0)
            .First(s => s.Position == 31).Score;

        var hits = SpliceSitePredictor.FindAcceptorSites(seq[..windowSize], threshold).ToList();

        bool expected = windowSize >= 32 && plantedScore >= threshold;
        hits.Any(s => s.Position == 31 && s.Type == SpliceSitePredictor.SpliceSiteType.Acceptor)
            .Should().Be(expected, "the acceptor is reported iff scannable in the window and above threshold");

        hits.Should().OnlyContain(s => s.Score >= threshold);
    }

    /// <summary>
    /// Interaction witness: a strong polypyrimidine tract upstream of the AG scores higher than a
    /// purine-rich (poly-A) upstream — the pyrimidine context drives the score.
    /// </summary>
    [Test]
    public void SpliceAcceptor_PyrimidineTract_RaisesScore()
    {
        double strong = SpliceSitePredictor.FindAcceptorSites(new string('C', 30) + "AG" + new string('G', 30), 0.0)
            .First(s => s.Position == 31).Score;
        double weak = SpliceSitePredictor.FindAcceptorSites(new string('A', 30) + "AG" + new string('G', 30), 0.0)
            .First(s => s.Position == 31).Score;
        strong.Should().BeGreaterThan(weak, "a polypyrimidine tract is favourable for the acceptor");
    }

    /// <summary>
    /// Interaction witness: every acceptor ends in the invariant AG and an AG-free sequence yields none.
    /// </summary>
    [Test]
    public void SpliceAcceptor_RequiresAg()
    {
        SpliceSitePredictor.FindAcceptorSites(AcceptorSequence(90), 0.0)
            .Should().OnlyContain(s => s.Motif.Contains("AG"), "acceptor motifs contain the AG intron end");
        SpliceSitePredictor.FindAcceptorSites(new string('C', 60), 0.0).Should().BeEmpty("no AG ⇒ no acceptor");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-PREDICT-001 — Intron / gene-structure prediction (Splicing)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 79.
    // Spec: tests/TestSpecs/SPLICE-PREDICT-001.md (canonical PredictIntrons /
    //       PredictGeneStructure).
    // Dimensions: minIntron(3) × maxIntron(3) × minExon(3) × scoring(2). Grid 3×3×3×2 = 54.
    //
    // Model: an intron is a donor→acceptor span whose length lies in [minIntron, maxIntron] and
    // whose combined splice-site score clears the threshold; the gene structure derives exons
    // (≥ minExon) between non-overlapping introns, with the spliced product = concatenated exons.
    //
    // The combinatorial point: the two length bounds, the exon-length floor and the score
    // threshold jointly constrain the output — the grid asserts every predicted intron/exon obeys
    // all four bounds for all 54 parameter combinations, and the spliced sequence is consistent.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string SpliceGene =
        new string('A', 40) + ("GTAAGT" + new string('A', 40) + new string('C', 20) + "AG") +
        new string('A', 40) + ("GTAAGT" + new string('A', 40) + new string('C', 20) + "AG") +
        new string('A', 40);

    [Test, Combinatorial]
    public void SplicePredict_OutputObeysAllBounds(
        [Values(40, 60, 70)] int minIntron,
        [Values(80, 120, 200)] int maxIntron,
        [Values(20, 30, 50)] int minExon,
        [Values(0.3, 0.6)] double scoring)
    {
        var introns = SpliceSitePredictor.PredictIntrons(SpliceGene, minIntron, maxIntron, scoring).ToList();
        introns.Should().OnlyContain(i => i.Length >= minIntron && i.Length <= maxIntron, "intron length within bounds");
        introns.Should().OnlyContain(i => i.Score >= scoring, "intron score clears the threshold");
        introns.Should().OnlyContain(i => i.Length == i.End - i.Start + 1, "length matches the span");

        var gs = SpliceSitePredictor.PredictGeneStructure(SpliceGene, minExon, minIntron, scoring);
        gs.Exons.Should().OnlyContain(e => e.Length >= minExon, "exons meet the minimum length");
        gs.Introns.Should().OnlyContain(i => i.Length >= minIntron && i.Score >= scoring);
        gs.SplicedSequence.Should().Be(string.Concat(gs.Exons.Select(e => e.Sequence)),
            "the spliced product is the concatenation of the reported exons");
    }

    /// <summary>
    /// Interaction witness: each length / score bound is monotone — relaxing maxIntron, lowering
    /// minIntron, or lowering the score threshold can only admit more introns (superset).
    /// </summary>
    [Test]
    public void SplicePredict_BoundsAreMonotone()
    {
        HashSet<(int, int)> Spans(int minI, int maxI, double sc) =>
            SpliceSitePredictor.PredictIntrons(SpliceGene, minI, maxI, sc).Select(i => (i.Start, i.End)).ToHashSet();

        Spans(40, 80, 0.3).Should().BeSubsetOf(Spans(40, 200, 0.3), "larger maxIntron admits more");
        Spans(70, 200, 0.3).Should().BeSubsetOf(Spans(40, 200, 0.3), "smaller minIntron admits more");
        Spans(40, 200, 0.6).Should().BeSubsetOf(Spans(40, 200, 0.3), "lower score threshold admits more");
    }

    /// <summary>
    /// Interaction witness: the engineered gene yields at least one intron under permissive
    /// settings (the contract grid is not vacuous), and predicted introns start at a donor GU.
    /// </summary>
    [Test]
    public void SplicePredict_NonVacuous_AndIntronsStartAtDonor()
    {
        var introns = SpliceSitePredictor.PredictIntrons(SpliceGene, 40, 200, 0.3).ToList();
        introns.Should().NotBeEmpty("the engineered donor/acceptor pairs form introns");
        introns.Should().OnlyContain(i => i.DonorSite.Type == SpliceSitePredictor.SpliceSiteType.Donor
                                          || i.DonorSite.Type == SpliceSitePredictor.SpliceSiteType.U12Donor);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-MAXENT3-001 — MaxEntScan 3' acceptor score (Splicing)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 250.
    // Spec: tests/TestSpecs/SPLICE-MAXENT3-001.md (SpliceSitePredictor.ScoreAcceptorMaxEnt). ADVANCED §10.
    //
    // Sources: Yeo & Burge (2004) J Comput Biol 11:377 (maximum-entropy splice-site model, score3).
    //
    // Model: score3 scores a 23-nt acceptor window (20 intron + 3 exon) whose intron/exon boundary is the
    //   invariant AG dinucleotide; the published reference window scores 2.89 bits.
    //
    // Dimensions: windowOffset(3) × dinucleotide(2). Grid 3×2 = 6 (exhaustive).
    //
    // The combinatorial point: the upstream polypyrimidine context (windowOffset) and the boundary
    // dinucleotide interact. In every context the canonical AG acceptor scores strictly higher than the
    // same window with a disrupted (non-AG) boundary — the model's defining splice signal.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] AcceptorUpstreams =
    {
        "ttccaaacgaacttttgt", // reference (score3 = 2.89 at canonical AG)
        "tctctctctctctctcct", // strong polypyrimidine tract
        "ttttccttttccttttgt", // alternative pyrimidine context
    };

    private static string AcceptorWindow(int offset, string dinucleotide) =>
        AcceptorUpstreams[offset] + dinucleotide + "gga"; // 18 + 2 + 3 = 23 nt

    [Test, Combinatorial]
    public void AcceptorMaxEnt_CanonicalAgBeatsDisrupted_AcrossOffsetAndDinucleotide(
        [Values(0, 1, 2)] int windowOffset,
        [Values("AG", "cc")] string dinucleotide)
    {
        string window = AcceptorWindow(windowOffset, dinucleotide);
        double score = SpliceSitePredictor.ScoreAcceptorMaxEnt(window);

        double.IsNaN(score).Should().BeFalse("score3 is finite for a 23-nt window");
        SpliceSitePredictor.ScoreAcceptorMaxEnt(window).Should().Be(score, "score3 is deterministic");

        double canonical = SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorWindow(windowOffset, "AG"));
        double disrupted = SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorWindow(windowOffset, "cc"));
        canonical.Should().BeGreaterThan(disrupted, "the invariant AG dinucleotide is required for a strong 3' acceptor");
    }

    /// <summary>Worked-example anchor: the published reference window scores 2.89 bits (Yeo &amp; Burge 2004).</summary>
    [Test]
    public void AcceptorMaxEnt_ReferenceWindow_Reproduces289()
    {
        Math.Round(SpliceSitePredictor.ScoreAcceptorMaxEnt("ttccaaacgaacttttgtAGgga"), 2)
            .Should().Be(2.89, "MaxEntScan score3 of the canonical documented window is 2.89 bits");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-MAXENT5-001 — MaxEntScan 5' donor score (Splicing)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 251.
    // Spec: tests/TestSpecs/SPLICE-MAXENT5-001.md (SpliceSitePredictor.ScoreDonorMaxEnt). ADVANCED §10.
    //
    // Sources: Yeo & Burge (2004) J Comput Biol 11:377 (maximum-entropy splice-site model, score5).
    //
    // Model: score5 scores a 9-nt donor window (3 exon + 6 intron) whose exon/intron boundary is the
    //   invariant GT dinucleotide; the published reference window cagGTAAGT scores 10.86 bits.
    //
    // Dimensions: windowOffset(3) × dinucleotide(2). Grid 3×2 = 6 (exhaustive).
    //
    // The combinatorial point: the exonic context (windowOffset) and the boundary dinucleotide interact.
    // In every context the canonical GT donor scores strictly higher than the same window with a
    // disrupted (non-GT) boundary.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] DonorExons = { "cag", "caa", "gag" };

    private static string DonorWindow(int offset, string dinucleotide) =>
        DonorExons[offset] + dinucleotide + "AAGT"; // 3 + 2 + 4 = 9 nt

    [Test, Combinatorial]
    public void DonorMaxEnt_CanonicalGtBeatsDisrupted_AcrossOffsetAndDinucleotide(
        [Values(0, 1, 2)] int windowOffset,
        [Values("GT", "cc")] string dinucleotide)
    {
        string window = DonorWindow(windowOffset, dinucleotide);
        double score = SpliceSitePredictor.ScoreDonorMaxEnt(window);

        double.IsNaN(score).Should().BeFalse("score5 is finite for a 9-nt window");
        SpliceSitePredictor.ScoreDonorMaxEnt(window).Should().Be(score, "score5 is deterministic");

        double canonical = SpliceSitePredictor.ScoreDonorMaxEnt(DonorWindow(windowOffset, "GT"));
        double disrupted = SpliceSitePredictor.ScoreDonorMaxEnt(DonorWindow(windowOffset, "cc"));
        canonical.Should().BeGreaterThan(disrupted, "the invariant GT dinucleotide is required for a strong 5' donor");
    }

    /// <summary>Worked-example anchor: the published reference window scores 10.86 bits (Yeo &amp; Burge 2004).</summary>
    [Test]
    public void DonorMaxEnt_ReferenceWindow_Reproduces1086()
    {
        Math.Round(SpliceSitePredictor.ScoreDonorMaxEnt("cagGTAAGT"), 2)
            .Should().Be(10.86, "MaxEntScan score5 of the canonical documented window is 10.86 bits");
    }
}
