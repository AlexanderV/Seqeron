using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for paired-end-mapping structural-variant detection — SV-DETECT-001.
/// The unit under test is the find-discordant → cluster → classify pipeline in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs:
///   • <see cref="StructuralVariantAnalyzer.DetectSVs"/> — canonical entry point that
///     flags discordant read pairs, clusters nearby supporters and emits one
///     <see cref="StructuralVariantAnalyzer.StructuralVariant"/> per cluster meeting
///     the minimum read-pair support;
///   • <see cref="StructuralVariantAnalyzer.ClassifySV"/> — maps one read-pair PEM
///     signature to an <see cref="StructuralVariantAnalyzer.SVType"/>;
///   • <see cref="StructuralVariantAnalyzer.FindDiscordantPairs"/> — the μ ± c·σ span
///     cutoff and FR-orientation anomaly test (exercised via DetectSVs and directly).
///
/// SCOPE. This file is scoped strictly to SV-DETECT-001 — SV typing and calling from
/// paired-end *mapping signatures* (DEL/DUP/INV/INS/TRA, size, position). Split-read
/// breakpoint localization (SV-BREAKPOINT-001) and read-depth CNV calling
/// (SV-CNV-001) are SEPARATE units handled by their own fuzz fixtures and are not
/// exercised here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts that the code
/// NEVER fails in an undisciplined way: no hang, no nonsense output, no *unhandled*
/// runtime exception (e.g. an overflow on a MaxInt position/insert size, a Min/Max
/// over an empty cluster, a NaN/Infinity leaking into Quality/Length), and no silent
/// corruption of the support count or SV type. Every input must resolve to EITHER a
/// well-defined, theory-correct value OR a documented, intentional outcome
/// (ArgumentNullException for null readPairs). The headline boundaries for THIS unit
/// (checklist row 203: "identical genomes, overlapping SVs, empty") are:
///   • identical genomes → a donor with NO rearrangement maps as all-concordant FR
///     pairs inside μ ± c·σ ⇒ FindDiscordantPairs flags NONE ⇒ DetectSVs emits NO SV
///     (§2.2, §3.3, §6.1, INV-05);
///   • overlapping SVs → many discordant supporters within clusterDistance collapse to
///     ONE SV; supporters > clusterDistance apart split into separate SVs; the SV type
///     is the cluster representative's signature, support = cluster size (§4.1, INV-06);
///   • empty input → an empty, non-null result; no Min/Max-of-empty crash (§6.1);
///   • 0 / -1 / MaxInt positions & insert sizes → no overflow into a negative Length or
///     a mis-signed span comparison (Boundary Exploitation MaxInt; §2.2).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SV-DETECT-001 — SV detection from paired-end signatures (StructuralVar)
/// Checklist: docs/checklists/03_FUZZING.md, row 203.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 203): "identical genomes, overlapping SVs, empty".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (SV_Detection.md)
/// (docs/algorithms/StructuralVar/SV_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • discordant-by-span iff insert size strictly outside [μ−c·σ, μ+c·σ]; bounds are
///     inclusive ⇒ a span exactly on a bound is CONCORDANT      (§2.2, §6.1, INV-05)
///   • only FR (strand1 '+', strand2 '-') is concordant orientation; FF/RR and RF are
///     discordant                                               (§2.2 ¶ after table, §4.1[1])
///   • classification order (first match wins): chr1≠chr2 → Translocation; same-strand
///     → Inversion; RF ('-','+') → Duplication; FR & span>μ+c·σ → Deletion; FR &
///     span<μ−c·σ → Insertion; else → ComplexRearrangement      (§2.2, §4.2)
///   • inter-chromosomal precedence over orientation (ASM-02)   (§2.3, INV-01)
///   • DetectSVs emits an SV per cluster iff supporters ≥ minSupport (BreakDancer -r,
///     default 2)                                               (§3.1, §4.1[3], INV-06)
///   • a cluster's Start = min Position1, End = max Position2, Length = |End−Start|,
///     SupportingReads = cluster size                           (§3.2, §4)
///   • defaults: μ 400, σ 50, c 3, clusterDistance 500, minSupport 2  (§3.1)
///   • empty readPairs ⇒ empty result; null readPairs ⇒ ArgumentNullException (§3.3, §6.1)
///   • worked example: three FR pairs span 5000 (>550) same chr ⇒ ONE Deletion,
///     SupportingReads 3                                        (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class StructuralVariantDetectionFuzzTests
{
    // Documented library/gate defaults (§3.1), mirrored locally so the test owns the
    // span bounds and the support gate independently of the source constants. With
    // μ=400, σ=50, c=3 the concordant span window is the inclusive [250, 550].
    private const int Mean = 400;            // expected insert size μ
    private const int Sd = 50;               // insert-size σ
    private const double Cutoff = 3.0;       // anomaly cutoff c (BreakDancer -c)
    private const int LowerBound = 250;      // μ − c·σ = 400 − 150
    private const int UpperBound = 550;      // μ + c·σ = 400 + 150
    private const int DefaultMinSupport = 2; // BreakDancer -r
    private const int DefaultClusterDistance = 500;

    // Builds a read-pair tuple in the documented field order
    // (ReadId, Chr1, Pos1, Strand1, Chr2, Pos2, Strand2, InsertSize) — §3.1.
    private static (string, string, int, char, string, int, char, int) Pair(
        string readId, string chr1, int pos1, char strand1,
        string chr2, int pos2, char strand2, int insertSize)
        => (readId, chr1, pos1, strand1, chr2, pos2, strand2, insertSize);

    // A concordant forward-reverse pair: same chromosome, '+'/'-', span IN [250,550].
    // This is what an *identical* (un-rearranged) donor genome maps as (§2.1, §2.2).
    private static (string, string, int, char, string, int, char, int) Concordant(
        string readId, string chrom, int pos, int span)
        => Pair(readId, chrom, pos, '+', chrom, pos + span, '-', span);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural invariants on EVERY accepted SV set, regardless
    // of input. This is what stops a fuzz test from rubber-stamping a green call:
    // support must equal a real cluster size ≥ the gate, Length/Quality must be finite
    // and non-negative, the type must be a defined enum member, and the SV count can
    // never exceed the input pair count (each pair joins at most one reported cluster).
    private static void AssertWellFormed(
        IReadOnlyList<StructuralVariant> svs, int inputPairCount, int minSupport)
    {
        svs.Count.Should().BeLessThanOrEqualTo(
            inputPairCount, "each pair joins at most one reported cluster (§4)");

        foreach (var sv in svs)
        {
            sv.SupportingReads.Should().BeGreaterThanOrEqualTo(
                minSupport, "every emitted SV passed the minimum-support gate (INV-06)");
            sv.SupportingReads.Should().BeGreaterThan(0, "a cluster has at least one member");

            sv.Length.Should().BeGreaterThanOrEqualTo(0, "Length = |End − Start| is non-negative (§3.2)");
            Enum.IsDefined(typeof(SVType), sv.Type).Should().BeTrue("type is a defined SVType (§3.2)");

            double.IsNaN(sv.Quality).Should().BeFalse("quality is a finite score");
            double.IsInfinity(sv.Quality).Should().BeFalse("quality is a finite score");
        }
    }

    #region SV-DETECT-001 — Positive sanity (documented worked example)

    [Test]
    public void DetectSVs_DocumentedWorkedExample_OneDeletionSupportThree()
    {
        // Docs §7.1: three FR pairs on chr1, span 5000 > μ+c·σ = 550 ⇒ all three are
        // discordant deletion-signature pairs; their Position1 {1000,1010,1020} lie
        // within clusterDistance 500 ⇒ ONE cluster, classified Deletion (FR, span
        // larger). Start = min Pos1 = 1000, End = max Pos2 = 6020, Length = 5020,
        // SupportingReads = 3. Hand-checked from the doc, not the code.
        var pairs = new[]
        {
            Pair("r1", "chr1", 1000, '+', "chr1", 6000, '-', 5000),
            Pair("r2", "chr1", 1010, '+', "chr1", 6010, '-', 5000),
            Pair("r3", "chr1", 1020, '+', "chr1", 6020, '-', 5000),
        };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff,
            clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport).ToList();

        svs.Should().HaveCount(1, "three deletion-signature pairs within clusterDistance ⇒ one SV (§7.1)");
        svs[0].Type.Should().Be(SVType.Deletion, "same chr, FR, span 5000 > 550 ⇒ Deletion (§2.2)");
        svs[0].SupportingReads.Should().Be(3, "cluster size = 3 (§7.1, INV-06)");
        svs[0].Start.Should().Be(1000, "min Position1 over the cluster (§3.2)");
        svs[0].End.Should().Be(6020, "max Position2 over the cluster (§3.2)");
        svs[0].Length.Should().Be(5020, "|End − Start| = 6020 − 1000 (§3.2)");
        svs[0].Chromosome.Should().Be("chr1");
        AssertWellFormed(svs, inputPairCount: 3, minSupport: DefaultMinSupport);
    }

    [Test]
    public void ClassifySV_EachSignature_MapsToItsDocumentedType()
    {
        // Each PEM signature, hand-derived from §2.2 / §4.2, in classification order.
        // chr1 ≠ chr2 ⇒ Translocation (precedence over orientation, INV-01).
        ClassifySV(new ReadPairSignature("t", "chr1", 0, '+', "chr2", 0, '-', Mean, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Translocation, "inter-chromosomal (§2.2, INV-01)");

        // Inter-chr AND same strand still ⇒ Translocation (ASM-02 precedence, M5).
        ClassifySV(new ReadPairSignature("t2", "chr1", 0, '+', "chr2", 0, '+', Mean, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Translocation, "chr difference precedes orientation (ASM-02)");

        // Same chr, same strand (FF) ⇒ Inversion.
        ClassifySV(new ReadPairSignature("i", "chr1", 0, '+', "chr1", 0, '+', Mean, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Inversion, "same-strand intra-chr ⇒ Inversion (§2.2)");

        // Same chr, RF (everted: '-','+') ⇒ Duplication.
        ClassifySV(new ReadPairSignature("d", "chr1", 0, '-', "chr1", 0, '+', Mean, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Duplication, "RF everted ⇒ tandem Duplication (§2.2)");

        // Same chr, FR, span > 550 ⇒ Deletion.
        ClassifySV(new ReadPairSignature("del", "chr1", 0, '+', "chr1", 0, '-', 5000, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Deletion, "FR span larger than insert ⇒ Deletion (§2.2)");

        // Same chr, FR, span < 250 ⇒ Insertion.
        ClassifySV(new ReadPairSignature("ins", "chr1", 0, '+', "chr1", 0, '-', 100, true),
            Mean, Sd, Cutoff).Should().Be(SVType.Insertion, "FR span smaller than insert ⇒ Insertion (§2.2)");
    }

    #endregion

    #region SV-DETECT-001 — BE boundary: identical genomes (no SVs called)

    [Test]
    public void FindDiscordantPairs_IdenticalGenome_AllConcordant_FlagsNone()
    {
        // "identical genomes" boundary: an un-rearranged donor maps as concordant FR
        // pairs whose span sits inside [250,550]. None is interchromosomal, none is
        // out-of-span, none is non-FR ⇒ FindDiscordantPairs flags ZERO (§2.2, INV-05).
        var rng = new Random(203_010);
        var pairs = Enumerable.Range(0, 50)
            .Select(i => Concordant($"c{i}", "chr1", pos: rng.Next(0, 100_000), span: rng.Next(LowerBound, UpperBound + 1)))
            .ToList();

        FindDiscordantPairs(pairs, Mean, Sd, Cutoff)
            .Should().BeEmpty("every pair is concordant FR within μ ± c·σ ⇒ no anomaly (§2.2, INV-05)");
    }

    [Test]
    public void DetectSVs_IdenticalGenome_CallsNoStructuralVariants()
    {
        // The same un-rearranged donor through the full pipeline ⇒ no discordant pairs
        // to cluster ⇒ NO SV emitted. This is the headline "identical genomes" case:
        // a variant caller must report nothing when the donor equals the reference
        // (§3.3 empty result of empty discordant set; §6.1).
        var pairs = Enumerable.Range(0, 30)
            .Select(i => Concordant($"c{i}", "chr1", pos: 1000 + i * 10, span: Mean))
            .ToList();

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff,
            clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport).ToList();

        svs.Should().BeEmpty("an identical (un-rearranged) genome yields no discordant signatures ⇒ no SV");
        AssertWellFormed(svs, inputPairCount: pairs.Count, minSupport: DefaultMinSupport);
    }

    [Test]
    public void FindDiscordantPairs_SpanExactlyOnBounds_IsConcordant_OneBeyond_IsDiscordant()
    {
        // §6.1 / INV-05: the bound is inclusive — span exactly μ ± c·σ (250 or 550) is
        // concordant; one unit beyond (249 / 551) is discordant. This pins that an
        // identical-genome span sitting ON the boundary is NOT mis-flagged.
        FindDiscordantPairs(new[] { Concordant("lo", "chr1", 1000, LowerBound) }, Mean, Sd, Cutoff)
            .Should().BeEmpty("span = μ − c·σ = 250 is on the inclusive bound ⇒ concordant (INV-05)");
        FindDiscordantPairs(new[] { Concordant("hi", "chr1", 1000, UpperBound) }, Mean, Sd, Cutoff)
            .Should().BeEmpty("span = μ + c·σ = 550 is on the inclusive bound ⇒ concordant (INV-05)");

        FindDiscordantPairs(new[] { Concordant("below", "chr1", 1000, LowerBound - 1) }, Mean, Sd, Cutoff)
            .Should().ContainSingle("span 249 is strictly below the lower bound ⇒ discordant (INV-05)");
        FindDiscordantPairs(new[] { Concordant("above", "chr1", 1000, UpperBound + 1) }, Mean, Sd, Cutoff)
            .Should().ContainSingle("span 551 is strictly above the upper bound ⇒ discordant (INV-05)");
    }

    #endregion

    #region SV-DETECT-001 — BE boundary: overlapping SVs (clustering)

    [Test]
    public void DetectSVs_OverlappingDeletionSupporters_CollapseToOneSv()
    {
        // "overlapping SVs" boundary: several deletion-signature pairs whose breakpoints
        // overlap (all within clusterDistance 500) are evidence for ONE event, so they
        // collapse into a single Deletion with support = the count (§4.1 clustering,
        // INV-06). Independently computed: Pos1 {1000,1100,1200,1300} all ≤ 500 apart
        // pairwise-adjacent ⇒ one cluster of 4.
        var pairs = new[]
        {
            Pair("a", "chr1", 1000, '+', "chr1", 6000, '-', 5000),
            Pair("b", "chr1", 1100, '+', "chr1", 6100, '-', 5000),
            Pair("c", "chr1", 1200, '+', "chr1", 6200, '-', 5000),
            Pair("d", "chr1", 1300, '+', "chr1", 6300, '-', 5000),
        };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff,
            clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport).ToList();

        svs.Should().HaveCount(1, "overlapping supporters within clusterDistance ⇒ one SV (§4.1)");
        svs[0].Type.Should().Be(SVType.Deletion);
        svs[0].SupportingReads.Should().Be(4, "all four overlapping pairs support one event (INV-06)");
        svs[0].Start.Should().Be(1000, "min Position1 (§3.2)");
        svs[0].End.Should().Be(6300, "max Position2 (§3.2)");
        AssertWellFormed(svs, inputPairCount: 4, minSupport: DefaultMinSupport);
    }

    [Test]
    public void DetectSVs_TwoDistantSvGroups_SplitIntoSeparateSvs()
    {
        // Two deletion-signature groups whose positions are FAR apart (gap ≫
        // clusterDistance) must NOT merge — they are two distinct events. Group A at
        // ~1000, group B at ~50000; the 49000 gap > 500 ⇒ two clusters, each support 2
        // (§4.1; the clustering split rule mirrors the breakpoint sibling).
        var pairs = new[]
        {
            Pair("a1", "chr1", 1000, '+', "chr1", 6000, '-', 5000),
            Pair("a2", "chr1", 1010, '+', "chr1", 6010, '-', 5000),
            Pair("b1", "chr1", 50000, '+', "chr1", 56000, '-', 6000),
            Pair("b2", "chr1", 50010, '+', "chr1", 56010, '-', 6000),
        };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff,
            clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport).ToList();

        svs.Should().HaveCount(2, "two groups beyond clusterDistance ⇒ two separate SVs (§4.1)");
        svs.Should().OnlyContain(sv => sv.Type == SVType.Deletion && sv.SupportingReads == 2);
        AssertWellFormed(svs, inputPairCount: 4, minSupport: DefaultMinSupport);
    }

    [Test]
    public void DetectSVs_OverlappingCluster_BelowMinSupport_EmitsNothing()
    {
        // A single overlapping pair (cluster size 1) is below the default minSupport 2
        // ⇒ NO SV, even though the pair is a valid deletion signature (INV-06 gate). The
        // support gate, not the signature, suppresses it.
        var pairs = new[] { Pair("solo", "chr1", 1000, '+', "chr1", 6000, '-', 5000) };

        DetectSVs(pairs, Mean, Sd, Cutoff, clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport)
            .Should().BeEmpty("one deletion-signature pair is below minSupport 2 ⇒ no SV (INV-06)");

        // The same pair with minSupport relaxed to 1 IS reported — confirms the gate,
        // not the classification, was the cause.
        DetectSVs(pairs, Mean, Sd, Cutoff, clusterDistance: DefaultClusterDistance, minSupport: 1)
            .Should().ContainSingle().Which.Type.Should().Be(SVType.Deletion);
    }

    #endregion

    #region SV-DETECT-001 — BE boundary: empty / null input

    [Test]
    public void DetectSVs_EmptyInput_ReturnsEmptyNonNull()
    {
        // §6.1: no pairs to classify ⇒ empty result. The Min/Max-over-cluster hazard
        // must not surface — the empty-discordant short-circuit guards it.
        DetectSVs(Array.Empty<(string, string, int, char, string, int, char, int)>(),
                Mean, Sd, Cutoff, clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport)
            .Should().NotBeNull().And.BeEmpty("empty readPairs ⇒ empty result (§6.1)");

        FindDiscordantPairs(Array.Empty<(string, string, int, char, string, int, char, int)>(), Mean, Sd, Cutoff)
            .Should().NotBeNull().And.BeEmpty("empty readPairs ⇒ no discordant pairs (§6.1)");
    }

    [Test]
    public void DetectSVs_NullInput_ThrowsArgumentNullException()
    {
        // §3.3: null readPairs ⇒ ArgumentNullException, thrown eagerly (the public
        // method validates before deferring to the iterator).
        Action detect = () => DetectSVs(null!, Mean, Sd, Cutoff);
        detect.Should().Throw<ArgumentNullException>();

        Action find = () => FindDiscordantPairs(null!, Mean, Sd, Cutoff).ToList();
        find.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region SV-DETECT-001 — BE boundary: extreme coordinates / insert sizes (0 / -1 / MaxInt)

    [Test]
    public void DetectSVs_ZeroAndNegativePositions_NoCrash_WellFormed()
    {
        // 0 / -1 positions (BE): coordinates at and below the origin must not break the
        // Min/Max/Length arithmetic. Two FR deletion-signature pairs at pos 0 and -1
        // (within clusterDistance) cluster into one Deletion; Length = |maxPos2 −
        // minPos1| stays non-negative.
        var pairs = new[]
        {
            Pair("z", "chr1", 0, '+', "chr1", 5000, '-', 5000),
            Pair("n", "chr1", -1, '+', "chr1", 4999, '-', 5000),
        };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff,
            clusterDistance: DefaultClusterDistance, minSupport: DefaultMinSupport).ToList();

        svs.Should().HaveCount(1, "two overlapping deletion pairs near the origin ⇒ one SV");
        svs[0].Type.Should().Be(SVType.Deletion);
        svs[0].Start.Should().Be(-1, "min Position1 (§3.2)");
        AssertWellFormed(svs, inputPairCount: 2, minSupport: DefaultMinSupport);
    }

    [Test]
    public void ClassifySV_MaxIntInsertSize_ClassifiesAsDeletion_NoOverflow()
    {
        // MaxInt insert size (BE): an FR same-chr pair with span int.MaxValue is far
        // above the upper bound ⇒ Deletion. The span comparison is against a double
        // bound (no int overflow); assert the documented signature, not a throw (§2.2).
        var pair = new ReadPairSignature("big", "chr1", 0, '+', "chr1", 0, '-', int.MaxValue, true);
        ClassifySV(pair, Mean, Sd, Cutoff)
            .Should().Be(SVType.Deletion, "span int.MaxValue ≫ μ + c·σ ⇒ Deletion (§2.2)");
    }

    [Test]
    [CancelAfter(30_000)]
    public void DetectSVs_RandomizedBoundarySweep_NeverThrows_WellFormed()
    {
        // Randomized BE sweep over the degenerate shapes: empty / singleton pair sets,
        // identical-genome (all-concordant) batches, overlapping vs distant clusters,
        // all four discordant signatures, and extreme/negative/MaxInt coordinates and
        // insert sizes. Every accepted result must satisfy the documented invariants and
        // the call must never throw or hang (CancelAfter 30s). The expected SV count and
        // per-cluster support are recomputed INDEPENDENTLY from §2.2/§4.1 (chr-grouped
        // discordant filter → linear cluster sweep → support gate), so a wrong
        // classification, cutoff or cluster rule would be caught.
        var rng = new Random(203_777);
        for (int trial = 0; trial < 3000; trial++)
        {
            int n = rng.Next(0, 14);                  // includes empty / singleton
            int minSupport = rng.Next(1, 4);
            int clusterDistance = rng.Next(0, 800);   // includes the 0 boundary
            int contigs = rng.Next(1, 3);

            var pairs = new List<(string, string, int, char, string, int, char, int)>(n);
            for (int i = 0; i < n; i++)
            {
                string chr1 = $"chr{rng.Next(0, contigs)}";
                // Mostly intra-chr (so DEL/DUP/INV/INS dominate), occasionally inter-chr.
                bool interChr = rng.Next(0, 6) == 0;
                string chr2 = interChr ? $"chr{rng.Next(0, contigs) + 100}" : chr1;

                char s1 = rng.Next(0, 2) == 0 ? '+' : '-';
                char s2 = rng.Next(0, 2) == 0 ? '+' : '-';

                int pos1 = rng.Next(0, 5) switch
                {
                    0 => 0,
                    1 => -rng.Next(0, 50),
                    2 => int.MaxValue - rng.Next(0, 6),
                    3 => 1000,                          // identical-collision bait
                    _ => rng.Next(0, 100_000),
                };

                int insertSize = rng.Next(0, 6) switch
                {
                    0 => Mean,                          // concordant span (identical genome)
                    1 => rng.Next(LowerBound, UpperBound + 1), // concordant span
                    2 => int.MaxValue,                  // deletion span (MaxInt)
                    3 => 0,                             // insertion span (BE 0)
                    4 => rng.Next(0, LowerBound),       // insertion span
                    _ => rng.Next(UpperBound + 1, 20_000), // deletion span
                };

                long p2 = (long)pos1 + insertSize;
                int pos2 = p2 > int.MaxValue ? int.MaxValue : (int)p2;
                pairs.Add(Pair($"r{i}", chr1, pos1, s1, chr2, pos2, s2, insertSize));
            }

            List<StructuralVariant> svs = null!;
            Action act = () => svs = DetectSVs(pairs, Mean, Sd, Cutoff, clusterDistance, minSupport).ToList();
            act.Should().NotThrow($"trial {trial}: degenerate read-pair set must not throw");

            AssertWellFormed(svs, pairs.Count, minSupport);

            // Independent spec recomputation (§2.2 discordant filter; §4.1 cluster sweep).
            var discordant = pairs
                .Where(p =>
                {
                    var (id, c1, x1, st1, c2, x2, st2, ins) = p;
                    return c1 != c2 || ins < LowerBound || ins > UpperBound || !(st1 == '+' && st2 == '-');
                })
                .Select(p => (p.Item2, p.Item3, p.Item5, p.Item6)) // (Chr1, Pos1, Chr2, Pos2)
                .OrderBy(t => t.Item1, StringComparer.CurrentCulture).ThenBy(t => t.Item2)
                .ToList();

            int expectedCount = 0;
            var expectedSupports = new List<int>();
            if (discordant.Count > 0)
            {
                int size = 1;
                for (int k = 1; k <= discordant.Count; k++)
                {
                    bool same = k < discordant.Count &&
                        discordant[k].Item1 == discordant[k - 1].Item1 &&     // Chr1
                        discordant[k].Item3 == discordant[k - 1].Item3 &&     // Chr2
                        Math.Abs((long)discordant[k].Item2 - discordant[k - 1].Item2) <= clusterDistance &&
                        Math.Abs((long)discordant[k].Item4 - discordant[k - 1].Item4) <= clusterDistance;
                    if (same)
                    {
                        size++;
                    }
                    else
                    {
                        if (size >= minSupport) { expectedCount++; expectedSupports.Add(size); }
                        size = 1;
                    }
                }
            }

            svs.Count.Should().Be(expectedCount,
                $"trial {trial}: emitted SV count must match the spec clustering (§2.2/§4.1)");
            svs.Select(s => s.SupportingReads).OrderBy(s => s)
                .Should().BeEquivalentTo(expectedSupports.OrderBy(s => s),
                    $"trial {trial}: per-cluster support must match cluster sizes (INV-06)");
        }
    }

    #endregion
}
