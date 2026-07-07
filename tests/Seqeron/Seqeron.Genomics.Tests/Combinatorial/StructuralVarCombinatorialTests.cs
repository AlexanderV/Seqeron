namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the StructuralVar area (StructuralVariantAnalyzer,
/// Seqeron.Genomics.Annotation).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("StructuralVar")]
public class StructuralVarCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SV-BREAKPOINT-001 — Split-read breakpoint clustering (StructuralVar)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 201.
    // Spec: tests/TestSpecs/SV-BREAKPOINT-001.md (canonical FindBreakpoints / RefineBreakpoint). ADVANCED §10.
    // Dimensions: splitReads(3) × minMapQ(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (ClipCrop, Suzuki 2011; SoftSearch, Hart 2013): split-read junctions are clustered within a
    // tolerance window and a breakpoint is reported only when its cluster has ≥ minSupport reads;
    // RefineBreakpoint returns the consensus (modal) junction coordinate.
    //
    // Axis mapping (documented — SplitRead carries no mapping quality): splitReads → the number of
    // reads supporting the planted junction; minMapQ → the minSupport threshold. The combinatorial
    // point: the breakpoint is reported exactly when support ≥ minSupport, with SupportingReads equal
    // to the cluster size and the refined coordinate equal to the planted junction.
    // ═══════════════════════════════════════════════════════════════════════

    private static StructuralVariantAnalyzer.SplitRead Split(string id, int junction) =>
        new(id, "chr1", junction - 30, junction, 30, "ACGT");

    [Test, Combinatorial]
    public void SvBreakpoint_ReportedWhenSupportMet_AcrossReadsAndMinSupport(
        [Values(2, 3, 4)] int splitReads,
        [Values(2, 3, 4)] int minSupport)
    {
        const int junction = 10_000;
        var reads = Enumerable.Range(0, splitReads).Select(i => Split($"r{i}", junction)).ToList();

        var breakpoints = StructuralVariantAnalyzer.FindBreakpoints(reads, clusterTolerance: 5, minSupport: minSupport).ToList();

        bool expectReported = splitReads >= minSupport;
        (breakpoints.Count != 0).Should().Be(expectReported, "a breakpoint is reported iff its support clears minSupport");

        if (expectReported)
        {
            var bp = breakpoints.Should().ContainSingle().Subject;
            bp.Position1.Should().Be(junction, "all reads share the junction coordinate");
            bp.SupportingReads.Should().Be(splitReads, "support equals the cluster size");
        }

        // RefineBreakpoint returns the consensus (modal) junction regardless of the support threshold.
        StructuralVariantAnalyzer.RefineBreakpoint("chr1", junction - 5, junction + 5, reads)
            .Should().Be(junction, "the modal junction is the consensus coordinate");
    }

    /// <summary>
    /// Interaction witness — junctions farther apart than the tolerance form separate breakpoints, and
    /// a region with no supporting read refines to null.
    /// </summary>
    [Test]
    public void SvBreakpoint_ToleranceSeparatesClusters_AndEmptyRegionIsNull()
    {
        var reads = new[]
        {
            Split("a", 1000), Split("b", 1002),     // within tolerance 5 ⇒ one cluster
            Split("c", 2000), Split("d", 2001),     // a second cluster far away
        };

        StructuralVariantAnalyzer.FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2)
            .Should().HaveCount(2, "two tolerance-separated junction clusters ⇒ two breakpoints");

        StructuralVariantAnalyzer.RefineBreakpoint("chr1", 5000, 6000, reads)
            .Should().BeNull("no split read's junction lies in the region");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SV-CNV-001 — Read-depth copy-number calling (StructuralVar)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 202.
    // Spec: tests/TestSpecs/SV-CNV-001.md (canonical DetectCNV / SegmentCopyNumber). ADVANCED §10.
    // Dimensions: binSize(3) × coverageRatio(3) × ploidy(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Yoon 2009 read-depth; CNVkit): a window's log2 ratio = log2(depth/reference) and its
    // integer copy number is CN = round(2 · 2^log2) = round(2 · ratio) for a diploid baseline.
    //
    // Axis mapping (documented): binSize → the depth window size; coverageRatio → depth/reference
    // (0.5/1.0/2.0 ⇒ CN 1/2/4); ploidy → the two CNV entry points (DetectCNV from depth vs
    // SegmentCopyNumber from log2 ratios), which must agree on the copy number. The combinatorial point:
    // both callers convert the coverage ratio to the same integer copy number at every bin size.
    // ═══════════════════════════════════════════════════════════════════════

    public enum CnvEntry { FromDepth, FromLogRatio }

    [Test, Combinatorial]
    public void SvCnv_CopyNumberFromCoverageRatio_AcrossBinSizeRatioEntry(
        [Values(50, 100, 200)] int binSize,
        [Values(0.5, 1.0, 2.0)] double coverageRatio,
        [Values(CnvEntry.FromDepth, CnvEntry.FromLogRatio)] CnvEntry ploidy)
    {
        int expectedCn = (int)Math.Round(2 * coverageRatio); // 0.5→1, 1.0→2, 2.0→4
        double log2 = Math.Log2(coverageRatio);

        if (ploidy == CnvEntry.FromDepth)
        {
            const int baseDepth = 100;
            int windowDepth = (int)Math.Round(baseDepth * coverageRatio);
            var depth = Enumerable.Repeat(windowDepth, binSize * 4).ToList();

            var segments = StructuralVariantAnalyzer.DetectCNV(depth, windowSize: binSize, referenceDepth: baseDepth).ToList();
            segments.Should().HaveCount(4, "four full windows");
            segments.Should().OnlyContain(s => s.CopyNumber == expectedCn, "CN = round(2·ratio)");
            segments.Should().OnlyContain(s => Math.Abs(s.LogRatio - log2) < 1e-9, "log2 ratio = log2(depth/reference)");
        }
        else
        {
            var seg = StructuralVariantAnalyzer.SegmentCopyNumber(Enumerable.Repeat(log2, 4)).Single();
            seg.CopyNumber.Should().Be(expectedCn, "the log-ratio caller agrees on CN = round(2·2^log2)");
            seg.ProbeCount.Should().Be(4, "four equal-CN windows merge into one segment");
        }
    }

    /// <summary>
    /// Interaction witness — copy number is monotone in coverage: a deletion (ratio 0.5) calls fewer
    /// copies than a normal region (1.0), which calls fewer than a duplication (2.0).
    /// </summary>
    [Test]
    public void SvCnv_CopyNumberMonotoneInCoverage()
    {
        int Cn(double ratio) => StructuralVariantAnalyzer.SegmentCopyNumber(new[] { Math.Log2(ratio) }).Single().CopyNumber;

        Cn(0.5).Should().BeLessThan(Cn(1.0));
        Cn(1.0).Should().BeLessThan(Cn(2.0));
        Cn(1.0).Should().Be(2, "a copy-neutral region is diploid");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SV-DETECT-001 — Paired-end SV detection & classification (StructuralVar)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 203.
    // Spec: tests/TestSpecs/SV-DETECT-001.md (canonical ClassifySV / DetectSVs). ADVANCED §10.
    // Dimensions: svType(3) × minSize(3) × support(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Medvedev 2009 PEM signatures; BreakDancer): a discordant read-pair signature classifies
    // by orientation/span — FR over-long ⇒ deletion, RF ⇒ duplication, same-strand ⇒ inversion,
    // interchromosomal ⇒ translocation; DetectSVs clusters discordant pairs and reports a cluster only
    // when its support ≥ minSupport.
    //
    // Axis mapping (documented): svType → the planted PEM signature; minSize → the number of supporting
    // read pairs planted; support → minSupport. The combinatorial point: ClassifySV returns the planted
    // type, and DetectSVs reports the SV (with that type and support count) exactly when the planted
    // support clears minSupport.
    // ═══════════════════════════════════════════════════════════════════════

    public enum PlantedSv { Deletion, Duplication, Inversion }

    private static (string, string, int, char, string, int, char, int) Pair(PlantedSv sv, int i)
    {
        int p1 = 1000 + i * 10;
        return sv switch
        {
            // FR, span ≫ insert size ⇒ deletion.
            PlantedSv.Deletion => ($"r{i}", "chr1", p1, '+', "chr1", p1 + 2000, '-', 2000),
            // RF (outward) ⇒ duplication.
            PlantedSv.Duplication => ($"r{i}", "chr1", p1, '-', "chr1", p1 + 400, '+', 400),
            // Same strand ⇒ inversion.
            _ => ($"r{i}", "chr1", p1, '+', "chr1", p1 + 400, '+', 400),
        };
    }

    [Test, Combinatorial]
    public void SvDetect_ClassifyAndReport_AcrossTypeSupportPlanted(
        [Values(PlantedSv.Deletion, PlantedSv.Duplication, PlantedSv.Inversion)] PlantedSv svType,
        [Values(2, 3, 4)] int nPairs,
        [Values(2, 3, 4)] int minSupport)
    {
        var pairs = Enumerable.Range(0, nPairs).Select(i => Pair(svType, i)).ToList();
        var expectedType = svType switch
        {
            PlantedSv.Deletion => StructuralVariantAnalyzer.SVType.Deletion,
            PlantedSv.Duplication => StructuralVariantAnalyzer.SVType.Duplication,
            _ => StructuralVariantAnalyzer.SVType.Inversion,
        };

        // ClassifySV maps the signature to the expected SV type.
        var sig = new StructuralVariantAnalyzer.ReadPairSignature("r", "chr1",
            pairs[0].Item3, pairs[0].Item4, "chr1", pairs[0].Item6, pairs[0].Item7, pairs[0].Item8, true);
        StructuralVariantAnalyzer.ClassifySV(sig).Should().Be(expectedType, "the PEM signature classifies the SV type");

        var svs = StructuralVariantAnalyzer.DetectSVs(pairs, minSupport: minSupport).ToList();

        bool expectReported = nPairs >= minSupport;
        (svs.Count != 0).Should().Be(expectReported, "a cluster is reported iff its support clears minSupport");
        if (expectReported)
        {
            var sv = svs.Should().ContainSingle().Subject;
            sv.Type.Should().Be(expectedType);
            sv.SupportingReads.Should().Be(nPairs, "support equals the cluster size");
        }
    }

    /// <summary>
    /// Interaction witness — the PEM signature rules: interchromosomal ⇒ translocation, a too-small
    /// FR span ⇒ insertion, and a normal FR pair is concordant (not discordant).
    /// </summary>
    [Test]
    public void SvDetect_SignatureRules()
    {
        StructuralVariantAnalyzer.SVType Classify(string c1, int p1, char s1, string c2, int p2, char s2, int ins) =>
            StructuralVariantAnalyzer.ClassifySV(new StructuralVariantAnalyzer.ReadPairSignature("r", c1, p1, s1, c2, p2, s2, ins, true));

        Classify("chr1", 100, '+', "chr2", 100, '-', 400).Should().Be(StructuralVariantAnalyzer.SVType.Translocation);
        Classify("chr1", 100, '+', "chr1", 150, '-', 50).Should().Be(StructuralVariantAnalyzer.SVType.Insertion, "span ≪ insert size");
        Classify("chr1", 100, '+', "chr1", 500, '-', 2000).Should().Be(StructuralVariantAnalyzer.SVType.Deletion, "span ≫ insert size");

        // A normal FR pair within the insert-size window is concordant ⇒ not flagged discordant.
        StructuralVariantAnalyzer.FindDiscordantPairs(new[] { ("r", "chr1", 100, '+', "chr1", 500, '-', 400) })
            .Should().BeEmpty("a normal FR pair is concordant");
    }
}
