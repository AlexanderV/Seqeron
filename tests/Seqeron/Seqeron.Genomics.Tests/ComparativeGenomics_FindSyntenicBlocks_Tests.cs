// COMPGEN-SYNTENY-001 — Synteny / Collinearity Block Detection
// Evidence: docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-SYNTENY-001.md
// Source: Wang Y et al. (2012). MCScanX. Nucleic Acids Research 40(7):e49.
//         https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336
//         MatchScore=50, GapPenalty=-1, NumberofGaps<25, min 5 collinear pairs / score 250.

using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ComparativeGenomics_FindSyntenicBlocks_Tests
{
    #region Helpers

    private static ComparativeGenomics.Gene G1(int idx, char strand = '+')
        => new($"g{idx}", "genome1", idx * 100, idx * 100 + 50, strand);

    private static ComparativeGenomics.Gene G2(int idx, char strand = '+')
        => new($"h{idx}", "genome2", idx * 100, idx * 100 + 50, strand);

    // Build n genome-1 genes g0..g(n-1) in order.
    private static List<ComparativeGenomics.Gene> Genome1(int n)
        => Enumerable.Range(0, n).Select(i => G1(i)).ToList();

    // Build genome-2 genes from a list of target indices (defines genome-2 order/positions).
    private static List<ComparativeGenomics.Gene> Genome2(IEnumerable<int> targetOrder)
        => targetOrder.Select((_, pos) => new ComparativeGenomics.Gene(
               $"h{pos}", "genome2", pos * 100, pos * 100 + 50, '+')).ToList();

    #endregion

    #region FindSyntenicBlocks — MUST Tests

    // M1 — Five adjacent forward anchors -> one forward block of 5 (score 5*50=250).
    // Source: MCScanX MatchScore=50, >=5 collinear pairs / score 250 (PMC3326336).
    [Test]
    public void FindSyntenicBlocks_FiveAdjacentForwardAnchors_ReturnsSingleForwardBlockOfFive()
    {
        // Arrange: g0..g4 in genome1; h0..h4 in genome2; ortholog gi -> hi (identical order).
        var genome1 = Genome1(5);
        var genome2 = Enumerable.Range(0, 5)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // Assert: exactly one forward block of 5 anchors (score 250 == threshold).
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "5 adjacent collinear anchors form exactly one block");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(5), "all 5 anchored pairs belong to the block");
            Assert.That(blocks[0].IsInverted, Is.False, "same genome-2 order => forward block");
        });
    }

    // M2 — Five anchors with reversed genome-2 order -> one inverted block.
    // Source: MCScanX sorts in both transcriptional directions (PMC3326336).
    [Test]
    public void FindSyntenicBlocks_FiveReversedAnchors_ReturnsInvertedBlock()
    {
        // Arrange: g0..g4; genome2 has h0..h4 but ortholog gi -> h(4-i) (reverse order).
        var genome1 = Genome1(5);
        var genome2 = Enumerable.Range(0, 5)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '-'))
            .ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{4 - i}");

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "reversed collinear anchors form one block");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(5), "all 5 reversed pairs belong to the block");
            Assert.That(blocks[0].IsInverted, Is.True, "decreasing genome-2 order => inverted block");
        });
    }

    // M3 — Four adjacent anchors (score 200) are below the 250 / 5-pair threshold -> no block.
    // Source: MCScanX requires >=5 collinear pairs (PMC3326336).
    [Test]
    public void FindSyntenicBlocks_FourAdjacentAnchors_ReturnsNoBlock()
    {
        // Arrange
        var genome1 = Genome1(4);
        var genome2 = Enumerable.Range(0, 4)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = Enumerable.Range(0, 4).ToDictionary(i => $"g{i}", i => $"h{i}");

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // Assert: 4 * 50 = 200 < 250 => not reported.
        Assert.That(blocks, Is.Empty, "4 collinear pairs score 200 < 250 and are not a reported block");
    }

    // M4 — A gap of >= maxGap intervening genes breaks the chain into sub-runs below threshold.
    // Source: MCScanX NumberofGaps < 25 (PMC3326336). Here maxGap=5 for a small input.
    [Test]
    public void FindSyntenicBlocks_GapExceedsMaxGap_BreaksChainAndReportsNothing()
    {
        // Arrange: genome1 g0..g5 (6 anchors); genome2 positions: 0,1,2 then jump to 50,51,52.
        // Gap between the two runs is |50-2|-1 = 47 intervening genes >> maxGap.
        var genome1 = Genome1(6);
        var positions = new[] { 0, 1, 2, 50, 51, 52 };
        var genome2 = Enumerable.Range(0, 60)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = new Dictionary<string, string>();
        for (int i = 0; i < positions.Length; i++)
            map[$"g{i}"] = $"h{positions[i]}";

        // Act: maxGap small so each run (3 anchors) stays separate and below the 5-pair minimum.
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map, maxGap: 5).ToList();

        // Assert
        Assert.That(blocks, Is.Empty, "gap of 47 genes splits anchors into two 3-pair runs, both < 5 pairs");
    }

    // M5 — Empty genome -> empty result, no exception.
    [Test]
    public void FindSyntenicBlocks_EmptyGenome1_ReturnsEmpty()
    {
        var genome1 = new List<ComparativeGenomics.Gene>();
        var genome2 = Genome2(new[] { 0, 1, 2, 3, 4 });
        var map = new Dictionary<string, string>();

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        Assert.That(blocks, Is.Empty, "no genes in genome1 => no anchors => no blocks");
    }

    // M6 — No orthologs -> empty result.
    [Test]
    public void FindSyntenicBlocks_NoOrthologs_ReturnsEmpty()
    {
        var genome1 = Genome1(5);
        var genome2 = Genome2(new[] { 0, 1, 2, 3, 4 });
        var map = new Dictionary<string, string>();

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        Assert.That(blocks, Is.Empty, "empty ortholog map => no anchors => no blocks");
    }

    // M7 — Two separated 5-anchor forward runs -> two non-overlapping blocks.
    // Source: MCScanX reports non-overlapping chains (PMC3326336).
    [Test]
    public void FindSyntenicBlocks_TwoSeparatedRuns_ReturnsTwoBlocks()
    {
        // Arrange: g0..g4 map to h0..h4 (run A); g5..g9 map to h20..h24 (run B).
        // The gap between runs (|20-4|-1 = 15) exceeds maxGap=5, so they are distinct blocks.
        var genome1 = Genome1(10);
        var genome2 = Enumerable.Range(0, 30)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = new Dictionary<string, string>();
        for (int i = 0; i < 5; i++) map[$"g{i}"] = $"h{i}";
        for (int i = 5; i < 10; i++) map[$"g{i}"] = $"h{15 + i}"; // h20..h24

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map, maxGap: 5).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(2), "two separated 5-anchor runs => two blocks");
            Assert.That(blocks.Sum(b => b.GeneCount), Is.EqualTo(10), "each anchor belongs to exactly one block");
            Assert.That(blocks.All(b => !b.IsInverted), Is.True, "both runs are forward");
        });
    }

    // M8 — Block coordinates span the parent genes' min/max and satisfy Start<=End.
    // Source: coordinate contract (INV-05).
    [Test]
    public void FindSyntenicBlocks_ForwardBlock_HasCoordinatesWithinParentGeneSpan()
    {
        // Arrange: g0..g4 have Start = i*100, End = i*100+50.
        var genome1 = Genome1(5);
        var genome2 = Enumerable.Range(0, 5)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        // Act
        var block = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).Single();

        // Assert: Start1 = g0.Start = 0; End1 = g4.End = 450; Start2 = 0; End2 = 450.
        Assert.Multiple(() =>
        {
            Assert.That(block.Start1, Is.EqualTo(0), "block starts at first gene Start");
            Assert.That(block.End1, Is.EqualTo(450), "block ends at last gene End (4*100+50)");
            Assert.That(block.Start2, Is.EqualTo(0), "genome-2 span starts at first target gene");
            Assert.That(block.End2, Is.EqualTo(450), "genome-2 span ends at last target gene");
            Assert.That(block.Start1, Is.LessThanOrEqualTo(block.End1), "Start1 <= End1");
            Assert.That(block.Start2, Is.LessThanOrEqualTo(block.End2), "Start2 <= End2");
        });
    }

    #endregion

    #region FindSyntenicBlocks — SHOULD Tests

    // S1 — Six anchors with one 1-gene gap still score >= 250 and form one block of 6.
    // Source: GapPenalty=-1; score = 6*50 - 1*1 = 299 >= 250 (PMC3326336).
    [Test]
    public void FindSyntenicBlocks_SixAnchorsWithOneGap_ReturnsBlockOfSix()
    {
        // Arrange: genome-2 positions 0,1,2,4,5,6 (one intervening gene between index 2 and 4).
        var genome1 = Genome1(6);
        var positions = new[] { 0, 1, 2, 4, 5, 6 };
        var genome2 = Enumerable.Range(0, 10)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = new Dictionary<string, string>();
        for (int i = 0; i < positions.Length; i++) map[$"g{i}"] = $"h{positions[i]}";

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "gapped but consistent chain forms one block");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(6), "all 6 anchors chain (score 6*50-1=299>=250)");
            Assert.That(blocks[0].IsInverted, Is.False, "increasing order => forward");
        });
    }

    // S2 — Null genome1 -> ArgumentNullException.
    [Test]
    public void FindSyntenicBlocks_NullGenome1_Throws()
    {
        var genome2 = Genome2(new[] { 0, 1, 2 });
        var map = new Dictionary<string, string>();

        Assert.Throws<ArgumentNullException>(
            () => ComparativeGenomics.FindSyntenicBlocks(null!, genome2, map).ToList(),
            "null genome1 must throw ArgumentNullException");
    }

    // S3 — An ortholog whose target gene is absent in genome2 is skipped, not a crash.
    [Test]
    public void FindSyntenicBlocks_OrthologTargetAbsent_SkipsAnchor()
    {
        // Arrange: 5 valid anchors plus one dangling ortholog (g5 -> hX not present).
        var genome1 = Genome1(6);
        var genome2 = Enumerable.Range(0, 5)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");
        map["g5"] = "hDoesNotExist";

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // Assert: dangling anchor ignored; the 5 valid anchors still form one block.
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "dangling ortholog is skipped, valid anchors still chain");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(5), "only the 5 resolvable anchors are counted");
        });
    }

    #endregion

    #region FindSyntenicBlocks — COULD Tests

    // C1 — VisualizeSynteny renders one descriptive line per block (delegate smoke test).
    [Test]
    public void VisualizeSynteny_OneForwardBlock_RendersDescriptiveLine()
    {
        var genome1 = Genome1(5);
        var genome2 = Enumerable.Range(0, 5)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        var text = ComparativeGenomics.VisualizeSynteny(blocks);

        Assert.Multiple(() =>
        {
            Assert.That(text, Does.Contain("genome1"), "rendering names genome1");
            Assert.That(text, Does.Contain("genome2"), "rendering names genome2");
            Assert.That(text, Does.Contain("forward"), "forward block labelled forward");
            Assert.That(text, Does.Contain("5 genes"), "rendering states the gene count");
        });
    }

    // C2 — Property: every reported block has GeneCount >= 5 and coordinates within parent bounds.
    // O(n^2)-class invariant property test (INV-01, INV-05).
    [Test]
    public void FindSyntenicBlocks_AllReportedBlocks_SatisfyMinSizeAndCoordinateInvariants()
    {
        // Arrange: a deterministic mixed input — one forward 5-run and one short 3-run.
        var genome1 = Genome1(8);
        var genome2 = Enumerable.Range(0, 30)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();
        var map = new Dictionary<string, string>();
        for (int i = 0; i < 5; i++) map[$"g{i}"] = $"h{i}";          // forward 5-run
        for (int i = 5; i < 8; i++) map[$"g{i}"] = $"h{20 + i}";     // 3-run, below threshold

        // Act
        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map, maxGap: 5).ToList();

        // Assert: invariants hold for all reported blocks.
        Assert.That(blocks, Is.Not.Empty, "at least the 5-run block is reported");
        Assert.Multiple(() =>
        {
            foreach (var b in blocks)
            {
                Assert.That(b.GeneCount, Is.GreaterThanOrEqualTo(5), "INV-01: reported block has >= 5 anchors");
                Assert.That(b.Start1, Is.LessThanOrEqualTo(b.End1), "INV-05: Start1 <= End1");
                Assert.That(b.Start2, Is.LessThanOrEqualTo(b.End2), "INV-05: Start2 <= End2");
            }
        });
    }

    #endregion
}
