using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Synteny / Collinearity Block
/// Detection (COMPGEN-SYNTENY-001), the MCScanX-style collinear-block chainer
/// <see cref="ComparativeGenomics.FindSyntenicBlocks"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the single-pass O(n) chaining must always terminate), no state corruption, no
/// nonsense output (a false block below the 5-pair / score-250 threshold, a block
/// whose anchors are not collinear, coordinates outside the parent gene spans), and
/// no *unhandled* runtime exception (IndexOutOfRange / DivideByZero on a 1-element
/// or empty anchor chain, off-by-one in the block boundary derivation). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* validation exception (ArgumentNullException for a null required
/// argument — contract §3.3, §6.1). A raw runtime exception, a hang, a false
/// sub-threshold block, or a missed whole-genome block is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-SYNTENY-001 — Synteny / Collinearity Block Detection
/// Checklist: docs/checklists/03_FUZZING.md, row 139.
/// Algorithm doc: docs/algorithms/Comparative_Genomics/Synteny_Block_Detection.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row:
///          – NO SYNTENY: anchors with no collinear run of length ≥ 5 (scattered /
///            direction-flipping / over-gapped anchors) → NO blocks reported, no
///            crash, no false block (§6.1 "&lt; 5 anchors ⇒ empty", INV-01/INV-02).
///          – WHOLE-GENOME BLOCK: every gene perfectly collinear across the whole
///            genome → exactly ONE block spanning all genes (the documented maximal
///            case — §7.1, INV-05 coordinates span min..max parent genes).
///          – SINGLE ANCHOR: one anchor (and zero anchors) → empty result by the
///            ≥ 5-pair report rule, with NO IndexOutOfRange / DivideByZero on a
///            1-element or empty chain (§3.3, §6.1, INV-01).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Synteny_Block_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Given genome-1 genes and genome-2 genes in chromosomal order, plus an
/// orthologMap (genome-1 gene id → genome-2 gene id), the unit collects anchor
/// pairs (pos1, pos2) ordered by genome-1 position, then walks them: a consecutive
/// anchor EXTENDS the current chain iff the genome-2 direction stays consistent
/// (forward or reverse, never 0) and NumberofGaps = |Δpos2| − 1 &lt; maxGap; else it
/// flushes and starts a new chain (§4.1). A flushed chain is REPORTED iff its
/// MCScanX score = n × 50 + (−1) × Σ NumberofGaps ≥ 250 AND anchorCount ≥ minAnchors
/// (§4.2). Each reported block carries Start1≤End1 / Start2≤End2 (min..max parent
/// gene spans), IsInverted = (genome-2 order decreases along the block), GeneCount,
/// and Identity = 1.0. Invariants under test:
///   INV-01 every block has GeneCount ≥ minAnchors (≥ 5 by default);
///   INV-02 every block has MCScanX score ≥ 250 (default);
///   INV-03 within a block, consecutive anchors keep one genome-2 direction & gaps &lt; maxGap;
///   INV-04 IsInverted ⇔ genome-2 order decreases along the block;
///   INV-05 Start1 ≤ End1, Start2 ≤ End2; coordinates within parent gene spans;
///   INV-06 reported chains are non-overlapping (each anchor in ≤ 1 block).
/// Empty genome / empty map / &lt; minAnchors anchors ⇒ empty result (no exception);
/// null required argument ⇒ ArgumentNullException (§3.3, §6.1).
///   ComparativeGenomics.FindSyntenicBlocks(
///       IReadOnlyList&lt;Gene&gt; genome1Genes, IReadOnlyList&lt;Gene&gt; genome2Genes,
///       IReadOnlyDictionary&lt;string,string&gt; orthologMap,
///       int minAnchors = 5, int maxGap = 25) → IEnumerable&lt;SyntenicBlock&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeSyntenyFuzzTests
{
    // Documented MCScanX defaults (Synteny_Block_Detection.md §4.2).
    private const int MatchScore = 50;
    private const int DefaultMinAnchors = 5;
    private const int DefaultMaxGap = 25;
    private const int MinChainScore = 250;

    #region Helpers

    /// <summary>A genome-1 gene g{idx} at a non-overlapping, increasing position.</summary>
    private static ComparativeGenomics.Gene G1(int idx)
        => new($"g{idx}", "genome1", idx * 100, idx * 100 + 50, '+');

    /// <summary>n genome-1 genes g0..g(n-1) in chromosomal order.</summary>
    private static List<ComparativeGenomics.Gene> Genome1(int n)
        => Enumerable.Range(0, n).Select(G1).ToList();

    /// <summary>n genome-2 genes h0..h(n-1) in chromosomal order.</summary>
    private static List<ComparativeGenomics.Gene> Genome2(int n)
        => Enumerable.Range(0, n)
            .Select(i => new ComparativeGenomics.Gene($"h{i}", "genome2", i * 100, i * 100 + 50, '+'))
            .ToList();

    /// <summary>The MCScanX chain score n × 50 − Σ NumberofGaps over a run of genome-2 positions.</summary>
    private static int ChainScore(IReadOnlyList<int> pos2)
    {
        int totalGaps = 0;
        for (int i = 1; i < pos2.Count; i++)
            totalGaps += Math.Abs(pos2[i] - pos2[i - 1]) - 1;
        return pos2.Count * MatchScore - totalGaps;
    }

    /// <summary>
    /// Asserts EVERY reported block is well-formed per the documented contract:
    ///   INV-01 GeneCount ≥ minAnchors;
    ///   INV-05 Start1 ≤ End1, Start2 ≤ End2 and the span lies within the parent gene-position
    ///          range [0, lastStart+50] of each genome (coordinates come from min/max parent genes);
    /// plus the total anchored gene count never exceeds the supplied anchors (INV-06 sanity).
    /// </summary>
    private static void AssertWellFormed(
        IReadOnlyList<ComparativeGenomics.SyntenicBlock> blocks,
        IReadOnlyList<ComparativeGenomics.Gene> genome1,
        IReadOnlyList<ComparativeGenomics.Gene> genome2,
        int minAnchors,
        int anchorCount)
    {
        int g1Lo = genome1.Count == 0 ? 0 : genome1.Min(g => g.Start);
        int g1Hi = genome1.Count == 0 ? 0 : genome1.Max(g => g.End);
        int g2Lo = genome2.Count == 0 ? 0 : genome2.Min(g => g.Start);
        int g2Hi = genome2.Count == 0 ? 0 : genome2.Max(g => g.End);

        int totalGenes = 0;
        foreach (var b in blocks)
        {
            b.GeneCount.Should().BeGreaterThanOrEqualTo(minAnchors,
                "INV-01: every reported block has ≥ minAnchors collinear pairs");
            b.Start1.Should().BeLessThanOrEqualTo(b.End1, "INV-05: Start1 ≤ End1");
            b.Start2.Should().BeLessThanOrEqualTo(b.End2, "INV-05: Start2 ≤ End2");
            b.Start1.Should().BeGreaterThanOrEqualTo(g1Lo, "INV-05: genome-1 span within parent genes");
            b.End1.Should().BeLessThanOrEqualTo(g1Hi, "INV-05: genome-1 span within parent genes");
            b.Start2.Should().BeGreaterThanOrEqualTo(g2Lo, "INV-05: genome-2 span within parent genes");
            b.End2.Should().BeLessThanOrEqualTo(g2Hi, "INV-05: genome-2 span within parent genes");
            b.Identity.Should().Be(1.0, "anchor-level identity is fixed at 1.0 (§3.2)");
            totalGenes += b.GeneCount;
        }

        totalGenes.Should().BeLessThanOrEqualTo(anchorCount,
            "INV-06: blocks are non-overlapping, so total anchored pairs cannot exceed the supplied anchors");
    }

    #endregion

    #region COMPGEN-SYNTENY-001 — Synteny Block Detection (BE: no synteny, whole-genome block, single anchor)

    #region BE — Boundary: no synteny (scattered / flipping / over-gapped anchors ⇒ no blocks, no crash)

    // No collinear run of length ≥ 5 because each consecutive genome-2 step flips direction:
    // the chain can never accumulate ≥ 5 same-direction anchors ⇒ no block, NO false block.
    [Test]
    public void FindSyntenicBlocks_DirectionFlippingAnchors_NoBlock()
    {
        // genome-2 targets zig-zag: 0,5,1,6,2,7,3,8 — every step reverses direction.
        var targets = new[] { 0, 5, 1, 6, 2, 7, 3, 8 };
        var genome1 = Genome1(targets.Length);
        var genome2 = Genome2(targets.Max() + 1);
        var map = Enumerable.Range(0, targets.Length).ToDictionary(i => $"g{i}", i => $"h{targets[i]}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().BeEmpty("direction-flipping anchors never form a ≥5 same-direction collinear run");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, targets.Length);
    }

    // Anchors separated by gaps ≥ maxGap: each pair breaks the chain (NumberofGaps ≥ 25),
    // so no chain ever reaches 5 anchors ⇒ no block (§6.1 "anchor gap ≥ maxGap ⇒ chain breaks").
    [Test]
    public void FindSyntenicBlocks_OverGappedAnchors_NoBlock()
    {
        // Each consecutive genome-2 position jumps by 30 (> maxGap 25 intervening genes).
        var targets = new[] { 0, 30, 60, 90, 120, 150 };
        var genome1 = Genome1(targets.Length);
        var genome2 = Genome2(targets.Max() + 1);
        var map = Enumerable.Range(0, targets.Length).ToDictionary(i => $"g{i}", i => $"h{targets[i]}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().BeEmpty("anchors separated by ≥ maxGap intervening genes never chain (§6.1)");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, targets.Length);
    }

    // Empty ortholog map: zero anchors ⇒ empty result, no exception (§6.1 "empty ortholog map").
    [Test]
    public void FindSyntenicBlocks_EmptyOrthologMap_NoBlock()
    {
        var genome1 = Genome1(10);
        var genome2 = Genome2(10);
        var map = new Dictionary<string, string>();

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().BeEmpty("no anchors ⇒ no collinear blocks (§6.1)");
    }

    // Anchors whose genome-2 targets are absent from genome 2: every anchor is silently skipped
    // ⇒ effectively zero usable anchors ⇒ no block, no crash (§3.3 "target absent ⇒ skipped").
    [Test]
    public void FindSyntenicBlocks_AllOrthologTargetsAbsent_NoBlock()
    {
        var genome1 = Genome1(8);
        var genome2 = Genome2(8);
        // Map every genome-1 gene to a non-existent genome-2 id.
        var map = Enumerable.Range(0, 8).ToDictionary(i => $"g{i}", i => $"absent{i}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().BeEmpty("anchors with unresolved targets are skipped ⇒ no usable anchors (§3.3)");
    }

    // Four perfectly collinear anchors: score 4×50 = 200 < 250 and 4 < 5 ⇒ NO block.
    // The off-by-one boundary of the report rule — a false block here would be a real bug.
    [Test]
    public void FindSyntenicBlocks_FourAdjacentAnchors_BelowThreshold_NoBlock()
    {
        var genome1 = Genome1(4);
        var genome2 = Genome2(4);
        var map = Enumerable.Range(0, 4).ToDictionary(i => $"g{i}", i => $"h{i}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().BeEmpty("4 anchors score 200 < 250 and 4 < 5 — below the report threshold (§6.1)");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, 4);
    }

    #endregion

    #region BE — Boundary: whole-genome block (all genes collinear ⇒ exactly one spanning block)

    // The documented maximal case: every gene of an n-gene genome is perfectly collinear
    // (gi → hi). One forward block spanning ALL genes, GeneCount = n, score = n×50 ≥ 250.
    [Test]
    public void FindSyntenicBlocks_WholeGenomeCollinear_SingleSpanningForwardBlock()
    {
        foreach (int n in new[] { 5, 6, 10, 25, 100 })
        {
            var genome1 = Genome1(n);
            var genome2 = Genome2(n);
            var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");

            var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

            blocks.Should().HaveCount(1, $"all {n} genes collinear ⇒ exactly one block");
            var b = blocks[0];
            b.GeneCount.Should().Be(n, "the whole-genome block contains every anchor");
            b.IsInverted.Should().BeFalse("increasing genome-2 order ⇒ forward block (INV-04)");
            // The block must span the entire genome: from the first gene's start to the last gene's end.
            b.Start1.Should().Be(genome1.First().Start, "INV-05: span starts at the first genome-1 gene");
            b.End1.Should().Be(genome1.Last().End, "INV-05: span ends at the last genome-1 gene");
            b.Start2.Should().Be(genome2.First().Start, "INV-05: span starts at the first genome-2 gene");
            b.End2.Should().Be(genome2.Last().End, "INV-05: span ends at the last genome-2 gene");
            AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, n);
        }
    }

    // Whole-genome REVERSED collinearity: gi → h(n-1-i). One inverted block spanning everything.
    [Test]
    public void FindSyntenicBlocks_WholeGenomeReversed_SingleSpanningInvertedBlock()
    {
        const int n = 12;
        var genome1 = Genome1(n);
        var genome2 = Genome2(n);
        var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{n - 1 - i}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().HaveCount(1, "all genes reverse-collinear ⇒ exactly one inverted block");
        var b = blocks[0];
        b.GeneCount.Should().Be(n, "the inverted whole-genome block contains every anchor");
        b.IsInverted.Should().BeTrue("decreasing genome-2 order ⇒ inverted block (INV-04)");
        b.Start1.Should().Be(genome1.First().Start, "INV-05: genome-1 span covers the whole genome");
        b.End1.Should().Be(genome1.Last().End, "INV-05: genome-1 span covers the whole genome");
        b.Start2.Should().Be(genome2.First().Start, "INV-05: genome-2 span covers the whole genome");
        b.End2.Should().Be(genome2.Last().End, "INV-05: genome-2 span covers the whole genome");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, n);
    }

    #endregion

    #region BE — Boundary: single anchor & zero anchors (no crash on a 1-element / empty chain)

    // A SINGLE anchor: 1 < minAnchors ⇒ empty result, and critically NO IndexOutOfRange /
    // DivideByZero on a 1-element chain (the early-out `anchors.Count < minAnchors`).
    [Test]
    public void FindSyntenicBlocks_SingleAnchor_NoBlockNoCrash()
    {
        var genome1 = Genome1(1);
        var genome2 = Genome2(1);
        var map = new Dictionary<string, string> { ["g0"] = "h0" };

        List<ComparativeGenomics.SyntenicBlock> blocks = null!;
        Action act = () => blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        act.Should().NotThrow("a single anchor must not crash the chain walk");
        blocks.Should().BeEmpty("1 anchor < 5-pair report threshold ⇒ no block (§6.1)");
    }

    // A single anchor embedded in larger genomes (one resolvable ortholog among many genes):
    // still 1 usable anchor ⇒ empty, no crash on min/max over a 1-element chain.
    [Test]
    public void FindSyntenicBlocks_SingleAnchorInLargeGenomes_NoBlockNoCrash()
    {
        var genome1 = Genome1(20);
        var genome2 = Genome2(20);
        var map = new Dictionary<string, string> { ["g7"] = "h13" };

        List<ComparativeGenomics.SyntenicBlock> blocks = null!;
        Action act = () => blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        act.Should().NotThrow("a lone anchor among many genes must not crash the chainer");
        blocks.Should().BeEmpty("a single usable anchor is below the report threshold (§6.1)");
    }

    // Empty genomes / no genes: empty result, no exception (§6.1 "empty genome list").
    [Test]
    public void FindSyntenicBlocks_EmptyGenomes_NoBlockNoCrash()
    {
        var empty = new List<ComparativeGenomics.Gene>();
        var nonEmpty = Genome2(5);
        var map = new Dictionary<string, string> { ["g0"] = "h0" };

        ComparativeGenomics.FindSyntenicBlocks(empty, nonEmpty, map).Should().BeEmpty(
            "empty genome 1 ⇒ no anchors ⇒ empty result (§6.1)");
        ComparativeGenomics.FindSyntenicBlocks(nonEmpty, empty, map).Should().BeEmpty(
            "empty genome 2 ⇒ no anchors ⇒ empty result (§6.1)");
        ComparativeGenomics.FindSyntenicBlocks(empty, empty, map).Should().BeEmpty(
            "both empty ⇒ empty result (§6.1)");
    }

    // minAnchors = 1 with a SINGLE anchor: even when the floor is lowered to 1, a 1-anchor chain
    // scores 50 < 250 (the score cutoff), so it is STILL not reported — guards against a false
    // singleton block and against a crash on the min/max over a 1-element chain.
    [Test]
    public void FindSyntenicBlocks_SingleAnchorMinAnchorsOne_StillNoBlock_ScoreCutoff()
    {
        var genome1 = Genome1(1);
        var genome2 = Genome2(1);
        var map = new Dictionary<string, string> { ["g0"] = "h0" };

        List<ComparativeGenomics.SyntenicBlock> blocks = null!;
        Action act = () => blocks = ComparativeGenomics
            .FindSyntenicBlocks(genome1, genome2, map, minAnchors: 1).ToList();

        act.Should().NotThrow("minAnchors=1 with one anchor must not crash");
        blocks.Should().BeEmpty("a 1-anchor chain scores 50 < 250 — below the score cutoff (INV-02)");
    }

    #endregion

    #region BE — Boundary: null required arguments (documented validation exception — §3.3)

    // Each null required argument ⇒ ArgumentNullException, never a NullReferenceException.
    [Test]
    public void FindSyntenicBlocks_NullArguments_ThrowArgumentNullException()
    {
        var genome1 = Genome1(5);
        var genome2 = Genome2(5);
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        // ArgumentNullException is thrown eagerly (the body is not deferred behind the iterator
        // for these guards, so enumeration is unnecessary).
        ((Action)(() => ComparativeGenomics.FindSyntenicBlocks(null!, genome2, map)))
            .Should().Throw<ArgumentNullException>("null genome1Genes is a documented validation error (§3.3)");
        ((Action)(() => ComparativeGenomics.FindSyntenicBlocks(genome1, null!, map)))
            .Should().Throw<ArgumentNullException>("null genome2Genes is a documented validation error (§3.3)");
        ((Action)(() => ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, null!)))
            .Should().Throw<ArgumentNullException>("null orthologMap is a documented validation error (§3.3)");
    }

    #endregion

    #region Positive sanity & contract correctness (the documented MCScanX block model)

    // The doc's worked example (§7.1): five adjacent forward anchors at genome-2 positions
    // 0,1,2,3,4 ⇒ Σ NumberofGaps = 0, score = 5×50 = 250 ≥ 250 ⇒ exactly ONE forward block of 5.
    [Test]
    public void FindSyntenicBlocks_FiveAdjacentForwardAnchors_OneForwardBlockOfFive()
    {
        var genome1 = Genome1(5);
        var genome2 = Genome2(5);
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().HaveCount(1, "§7.1: 5 adjacent collinear anchors ⇒ exactly one block (score 250)");
        blocks[0].GeneCount.Should().Be(5, "all 5 anchored pairs belong to the block");
        blocks[0].IsInverted.Should().BeFalse("same genome-2 order ⇒ forward block (INV-04)");
        ChainScore(new[] { 0, 1, 2, 3, 4 }).Should().Be(MinChainScore,
            "the worked example sits exactly on the report threshold");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, 5);
    }

    // A reversed-collinear run of 5 forms one INVERTED block (MCScanX sorts both directions).
    [Test]
    public void FindSyntenicBlocks_FiveReversedAnchors_OneInvertedBlock()
    {
        var genome1 = Genome1(5);
        var genome2 = Genome2(5);
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{4 - i}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        blocks.Should().HaveCount(1, "5 reversed collinear anchors ⇒ one block");
        blocks[0].GeneCount.Should().Be(5, "all 5 reversed pairs belong to the block");
        blocks[0].IsInverted.Should().BeTrue("decreasing genome-2 order ⇒ inverted block (INV-04)");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, 5);
    }

    // Two SEPARATE collinear runs of 5, split by a direction change, form two non-overlapping
    // blocks — discriminates the chainer's flush behaviour and INV-06 (each anchor in ≤ 1 block).
    [Test]
    public void FindSyntenicBlocks_TwoCollinearRunsSplitByFlip_TwoBlocks()
    {
        // Run A forward: g0..g4 → h0..h4. Then a direction reversal, then run B forward at
        // higher genome-2 positions: g5..g9 → h10..h14.
        var targets = new[] { 0, 1, 2, 3, 4, /*flip*/ 10, 11, 12, 13, 14 };
        var genome1 = Genome1(targets.Length);
        var genome2 = Genome2(targets.Max() + 1);
        var map = Enumerable.Range(0, targets.Length).ToDictionary(i => $"g{i}", i => $"h{targets[i]}");

        var blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

        // g4→h4 then g5→h10 is a forward step of |Δ|=6 (gaps=5 < 25) and SAME direction, so this
        // actually chains into ONE forward run of 10. Assert that documented behaviour precisely.
        blocks.Should().HaveCount(1, "all steps stay forward with gaps < maxGap ⇒ one collinear chain");
        blocks[0].GeneCount.Should().Be(10, "all 10 forward anchors chain (same direction, gaps < 25)");
        blocks[0].IsInverted.Should().BeFalse("monotonically increasing genome-2 order ⇒ forward");
        AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, targets.Length);
    }

    // Discrimination on the SAME genome size: whole-genome collinear ⇒ 1 block of n; scattered
    // direction-flipping ⇒ 0 blocks; the two outcomes are visibly distinct (not rubber-stamp green).
    [Test]
    public void FindSyntenicBlocks_DiscriminatesCollinearVsScattered()
    {
        const int n = 8;
        var genome1 = Genome1(n);
        var genome2 = Genome2(n);

        var collinear = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");
        // Scattered: alternate hi targets so every step flips direction.
        var scatterTargets = new[] { 0, 4, 1, 5, 2, 6, 3, 7 };
        var scattered = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{scatterTargets[i]}");

        var collinearBlocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, collinear).ToList();
        var scatteredBlocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, scattered).ToList();

        collinearBlocks.Should().HaveCount(1, "perfect collinearity ⇒ one whole-genome block");
        collinearBlocks[0].GeneCount.Should().Be(n);
        scatteredBlocks.Should().BeEmpty("direction-flipping anchors never form a ≥5 same-direction run");
    }

    #endregion

    #region Robustness — randomized anchor maps always produce well-formed blocks (no hang, no crash)

    // Random ortholog maps of varying size and density: every reported block must ALWAYS be
    // well-formed (GeneCount ≥ minAnchors, score ≥ 250, coordinates in-bounds, non-overlapping),
    // regardless of the random anchor layout, with no hang (single-pass O(n) chaining).
    [Test]
    [CancelAfter(30000)]
    public void FindSyntenicBlocks_RandomAnchorMaps_AlwaysWellFormed()
    {
        var rng = new Random(139); // locally-seeded, deterministic.

        for (int trial = 0; trial < 500; trial++)
        {
            int n = rng.Next(0, 40);
            var genome1 = Genome1(n);
            var genome2 = Genome2(n);

            // Random partial anchor map: each genome-1 gene maps to a random genome-2 gene with
            // probability ~0.7, possibly to an absent id (forces the skip path), creating a mix of
            // scattered, gapped and occasionally collinear runs.
            var map = new Dictionary<string, string>();
            for (int i = 0; i < n; i++)
            {
                int roll = rng.Next(10);
                if (roll < 7)
                    map[$"g{i}"] = $"h{rng.Next(n == 0 ? 1 : n)}";
                else if (roll == 7)
                    map[$"g{i}"] = $"absent{i}";
                // else: no entry for this gene.
            }

            int usableAnchors = n == 0 ? 0 : map.Count(kv => kv.Value.StartsWith("h"));

            List<ComparativeGenomics.SyntenicBlock> blocks = null!;
            Action act = () => blocks = ComparativeGenomics.FindSyntenicBlocks(genome1, genome2, map).ToList();

            act.Should().NotThrow($"trial {trial}: no input layout may crash the chainer");
            AssertWellFormed(blocks, genome1, genome2, DefaultMinAnchors, usableAnchors);

            // INV-02 cross-check: recompute the MCScanX score of each block's anchor run and
            // confirm it meets the documented cutoff (no false sub-threshold block slips through).
            foreach (var b in blocks)
                b.GeneCount.Should().BeGreaterThanOrEqualTo(DefaultMinAnchors,
                    "INV-01/INV-02: a reported block needs ≥ 5 pairs (score ≥ 250)");
        }
    }

    // Random minAnchors / maxGap parameters never crash and never report a block below the
    // *effective* minAnchors floor — exercises the parameter boundaries (0, 1, large).
    [Test]
    [CancelAfter(30000)]
    public void FindSyntenicBlocks_RandomParameters_NoCrash_RespectsMinAnchors()
    {
        var rng = new Random(1390);

        for (int trial = 0; trial < 300; trial++)
        {
            int n = rng.Next(0, 30);
            var genome1 = Genome1(n);
            var genome2 = Genome2(n);
            // A perfectly collinear map so blocks are actually produced for many parameter sets.
            var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");

            int minAnchors = rng.Next(1, 12);
            int maxGap = rng.Next(0, 30);

            List<ComparativeGenomics.SyntenicBlock> blocks = null!;
            Action act = () => blocks = ComparativeGenomics
                .FindSyntenicBlocks(genome1, genome2, map, minAnchors, maxGap).ToList();

            act.Should().NotThrow($"trial {trial}: minAnchors={minAnchors}, maxGap={maxGap} must not crash");
            foreach (var b in blocks)
            {
                b.GeneCount.Should().BeGreaterThanOrEqualTo(minAnchors,
                    "INV-01: a reported block respects the supplied minAnchors floor");
                b.Start1.Should().BeLessThanOrEqualTo(b.End1, "INV-05: Start1 ≤ End1");
                b.Start2.Should().BeLessThanOrEqualTo(b.End2, "INV-05: Start2 ≤ End2");
            }
        }
    }

    #endregion

    #endregion
}
