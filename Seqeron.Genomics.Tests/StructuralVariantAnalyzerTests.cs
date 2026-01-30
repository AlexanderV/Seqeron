using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class StructuralVariantAnalyzerTests
{
    #region Discordant Read Pair Tests

    [Test]
    public void FindDiscordantPairs_LargeInsertSize_MarksAsDiscordant()
    {
        var readPairs = new List<(string, string, int, char, string, int, char, int)>
        {
            ("read1", "chr1", 1000, '+', "chr1", 1400, '-', 400), // Normal - proper FR orientation
            ("read2", "chr1", 1000, '+', "chr1", 10000, '-', 9000) // Large insert
        };

        var discordant = StructuralVariantAnalyzer.FindDiscordantPairs(
            readPairs, expectedInsertSize: 400, insertSizeStdDev: 50).ToList();

        Assert.That(discordant, Has.Count.EqualTo(1));
        Assert.That(discordant[0].ReadId, Is.EqualTo("read2"));
        Assert.That(discordant[0].IsDiscordant, Is.True);
    }

    [Test]
    public void FindDiscordantPairs_Interchromosomal_MarksAsDiscordant()
    {
        var readPairs = new List<(string, string, int, char, string, int, char, int)>
        {
            ("read1", "chr1", 1000, '+', "chr2", 5000, '-', 0)
        };

        var discordant = StructuralVariantAnalyzer.FindDiscordantPairs(readPairs).ToList();

        Assert.That(discordant, Has.Count.EqualTo(1));
        Assert.That(discordant[0].Chromosome1, Is.EqualTo("chr1"));
        Assert.That(discordant[0].Chromosome2, Is.EqualTo("chr2"));
    }

    [Test]
    public void FindDiscordantPairs_AbnormalOrientation_MarksAsDiscordant()
    {
        var readPairs = new List<(string, string, int, char, string, int, char, int)>
        {
            ("read1", "chr1", 1000, '+', "chr1", 1500, '+', 500) // Same strand
        };

        var discordant = StructuralVariantAnalyzer.FindDiscordantPairs(readPairs).ToList();

        Assert.That(discordant, Has.Count.EqualTo(1));
    }

    [Test]
    public void FindDiscordantPairs_NormalPairs_NotReturned()
    {
        var readPairs = new List<(string, string, int, char, string, int, char, int)>
        {
            ("read1", "chr1", 1000, '+', "chr1", 1400, '-', 400)
        };

        var discordant = StructuralVariantAnalyzer.FindDiscordantPairs(
            readPairs, expectedInsertSize: 400, insertSizeStdDev: 50).ToList();

        Assert.That(discordant, Is.Empty);
    }

    [Test]
    public void ClusterDiscordantPairs_MultiplePairs_CreatesSV()
    {
        var pairs = new List<StructuralVariantAnalyzer.ReadPairSignature>
        {
            new("r1", "chr1", 1000, '+', "chr1", 10000, '-', 9000, true),
            new("r2", "chr1", 1050, '+', "chr1", 10050, '-', 9000, true),
            new("r3", "chr1", 1100, '+', "chr1", 10100, '-', 9000, true)
        };

        var svs = StructuralVariantAnalyzer.ClusterDiscordantPairs(pairs, minSupport: 3).ToList();

        Assert.That(svs, Has.Count.EqualTo(1));
        Assert.That(svs[0].SupportingReads, Is.EqualTo(3));
    }

    [Test]
    public void ClusterDiscordantPairs_InsufficientSupport_NoSV()
    {
        var pairs = new List<StructuralVariantAnalyzer.ReadPairSignature>
        {
            new("r1", "chr1", 1000, '+', "chr1", 10000, '-', 9000, true),
            new("r2", "chr1", 1050, '+', "chr1", 10050, '-', 9000, true)
        };

        var svs = StructuralVariantAnalyzer.ClusterDiscordantPairs(pairs, minSupport: 5).ToList();

        Assert.That(svs, Is.Empty);
    }

    [Test]
    public void ClusterDiscordantPairs_DifferentChromosomes_SeparateClusters()
    {
        var pairs = new List<StructuralVariantAnalyzer.ReadPairSignature>
        {
            new("r1", "chr1", 1000, '+', "chr1", 10000, '-', 9000, true),
            new("r2", "chr1", 1050, '+', "chr1", 10050, '-', 9000, true),
            new("r3", "chr1", 1100, '+', "chr1", 10100, '-', 9000, true),
            new("r4", "chr2", 1000, '+', "chr2", 10000, '-', 9000, true),
            new("r5", "chr2", 1050, '+', "chr2", 10050, '-', 9000, true),
            new("r6", "chr2", 1100, '+', "chr2", 10100, '-', 9000, true)
        };

        var svs = StructuralVariantAnalyzer.ClusterDiscordantPairs(pairs, minSupport: 3).ToList();

        Assert.That(svs, Has.Count.EqualTo(2));
    }

    #endregion

    #region Split Read Tests

    [Test]
    public void FindSplitReads_SoftClippedRead_Detected()
    {
        var alignments = new List<(string, string, int, string, string)>
        {
            ("read1", "chr1", 1000, "50S100M", "ACGT" + new string('N', 146))
        };

        var splits = StructuralVariantAnalyzer.FindSplitReads(alignments, minClipLength: 20).ToList();

        Assert.That(splits, Has.Count.EqualTo(1));
        Assert.That(splits[0].ClipLength, Is.EqualTo(50));
    }

    [Test]
    public void FindSplitReads_ShortClip_NotDetected()
    {
        var alignments = new List<(string, string, int, string, string)>
        {
            ("read1", "chr1", 1000, "10S100M", "ACGTACGTAC" + new string('N', 100))
        };

        var splits = StructuralVariantAnalyzer.FindSplitReads(alignments, minClipLength: 20).ToList();

        Assert.That(splits, Is.Empty);
    }

    [Test]
    public void FindSplitReads_RightClip_Detected()
    {
        var alignments = new List<(string, string, int, string, string)>
        {
            ("read1", "chr1", 1000, "100M50S", new string('N', 100) + "ACGT" + new string('N', 46))
        };

        var splits = StructuralVariantAnalyzer.FindSplitReads(alignments, minClipLength: 20).ToList();

        Assert.That(splits, Has.Count.EqualTo(1));
    }

    [Test]
    public void ClusterSplitReads_MultipleSplits_CreatesBreakpoint()
    {
        var splits = new List<StructuralVariantAnalyzer.SplitRead>
        {
            new("r1", "chr1", 1000, 5000, 50, "ACGT"),
            new("r2", "chr1", 1005, 5005, 50, "ACGT"),
            new("r3", "chr1", 1010, 5010, 50, "ACGT")
        };

        var breakpoints = StructuralVariantAnalyzer.ClusterSplitReads(splits, minSupport: 2).ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(1));
        Assert.That(breakpoints[0].SupportingReads, Is.EqualTo(3));
    }

    [Test]
    public void ClusterSplitReads_DistantReads_SeparateBreakpoints()
    {
        var splits = new List<StructuralVariantAnalyzer.SplitRead>
        {
            new("r1", "chr1", 1000, 5000, 50, "ACGT"),
            new("r2", "chr1", 1005, 5005, 50, "ACGT"),
            new("r3", "chr1", 2000, 6000, 50, "ACGT"),
            new("r4", "chr1", 2005, 6005, 50, "ACGT")
        };

        var breakpoints = StructuralVariantAnalyzer.ClusterSplitReads(
            splits, clusterDistance: 10, minSupport: 2).ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(2));
    }

    #endregion

    #region Copy Number Analysis Tests

    [Test]
    public void SegmentCopyNumber_ConstantValues_SingleSegment()
    {
        var probes = Enumerable.Range(0, 20)
            .Select(i => ("chr1", i * 1000, 0.0, 0.5))
            .ToList();

        var segments = StructuralVariantAnalyzer.SegmentCopyNumber(probes, minProbes: 5).ToList();

        Assert.That(segments, Has.Count.EqualTo(1));
        Assert.That(segments[0].CopyNumber, Is.EqualTo(2));
    }

    [Test]
    public void SegmentCopyNumber_DeletionRegion_DetectsCNVs()
    {
        var probes = new List<(string, int, double, double)>();

        // Normal region
        for (int i = 0; i < 10; i++)
            probes.Add(("chr1", i * 1000, 0.0, 0.5));

        // Deletion region
        for (int i = 10; i < 20; i++)
            probes.Add(("chr1", i * 1000, -1.0, 0.0));

        // Normal region
        for (int i = 20; i < 30; i++)
            probes.Add(("chr1", i * 1000, 0.0, 0.5));

        var segments = StructuralVariantAnalyzer.SegmentCopyNumber(
            probes, changeThreshold: 0.5, minProbes: 5).ToList();

        Assert.That(segments.Count, Is.GreaterThanOrEqualTo(2));
        Assert.That(segments.Any(s => s.CopyNumber < 2), Is.True);
    }

    [Test]
    public void IdentifyCNVs_DeletionSegment_ReturnsCNV()
    {
        var segments = new List<StructuralVariantAnalyzer.CopyNumberSegment>
        {
            new("chr1", 0, 100000, -1.0, 1, 0.0, 50)
        };

        var cnvs = StructuralVariantAnalyzer.IdentifyCNVs(segments, minLength: 10000).ToList();

        Assert.That(cnvs, Has.Count.EqualTo(1));
        Assert.That(cnvs[0].Type, Is.EqualTo(StructuralVariantAnalyzer.SVType.Deletion));
    }

    [Test]
    public void IdentifyCNVs_DuplicationSegment_ReturnsCNV()
    {
        var segments = new List<StructuralVariantAnalyzer.CopyNumberSegment>
        {
            new("chr1", 0, 100000, 0.58, 3, 0.33, 50) // Log2(3/2) ≈ 0.58
        };

        var cnvs = StructuralVariantAnalyzer.IdentifyCNVs(segments, minLength: 10000).ToList();

        Assert.That(cnvs, Has.Count.EqualTo(1));
        Assert.That(cnvs[0].Type, Is.EqualTo(StructuralVariantAnalyzer.SVType.Duplication));
    }

    [Test]
    public void IdentifyCNVs_NormalCopyNumber_NoCNV()
    {
        var segments = new List<StructuralVariantAnalyzer.CopyNumberSegment>
        {
            new("chr1", 0, 100000, 0.0, 2, 0.5, 50)
        };

        var cnvs = StructuralVariantAnalyzer.IdentifyCNVs(segments).ToList();

        Assert.That(cnvs, Is.Empty);
    }

    [Test]
    public void IdentifyCNVs_ShortSegment_Filtered()
    {
        var segments = new List<StructuralVariantAnalyzer.CopyNumberSegment>
        {
            new("chr1", 0, 1000, -1.0, 1, 0.0, 5)
        };

        var cnvs = StructuralVariantAnalyzer.IdentifyCNVs(segments, minLength: 10000).ToList();

        Assert.That(cnvs, Is.Empty);
    }

    #endregion

    #region SV Merging and Filtering Tests

    [Test]
    public void MergeOverlappingSVs_OverlappingDeletions_Merged()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 5000, StructuralVariantAnalyzer.SVType.Deletion, 4000, 50, 5, null),
            new("SV2", "chr1", 2000, 6000, StructuralVariantAnalyzer.SVType.Deletion, 4000, 60, 6, null)
        };

        var merged = StructuralVariantAnalyzer.MergeOverlappingSVs(svs, overlapFraction: 0.3).ToList();

        Assert.That(merged, Has.Count.EqualTo(1));
        Assert.That(merged[0].Start, Is.EqualTo(1000));
        Assert.That(merged[0].End, Is.EqualTo(6000));
    }

    [Test]
    public void MergeOverlappingSVs_DifferentTypes_NotMerged()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 5000, StructuralVariantAnalyzer.SVType.Deletion, 4000, 50, 5, null),
            new("SV2", "chr1", 2000, 6000, StructuralVariantAnalyzer.SVType.Duplication, 4000, 60, 6, null)
        };

        var merged = StructuralVariantAnalyzer.MergeOverlappingSVs(svs).ToList();

        Assert.That(merged, Has.Count.EqualTo(2));
    }

    [Test]
    public void MergeOverlappingSVs_NonOverlapping_NotMerged()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null),
            new("SV2", "chr1", 10000, 11000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 60, 6, null)
        };

        var merged = StructuralVariantAnalyzer.MergeOverlappingSVs(svs).ToList();

        Assert.That(merged, Has.Count.EqualTo(2));
    }

    [Test]
    public void FilterSVs_ByQuality_FiltersLowQuality()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 10, 5, null),
            new("SV2", "chr1", 5000, 6000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null)
        };

        var filtered = StructuralVariantAnalyzer.FilterSVs(svs, minQuality: 20).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("SV2"));
    }

    [Test]
    public void FilterSVs_BySupport_FiltersLowSupport()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 1, null),
            new("SV2", "chr1", 5000, 6000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 10, null)
        };

        var filtered = StructuralVariantAnalyzer.FilterSVs(svs, minSupport: 5).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
    }

    [Test]
    public void FilterSVs_ByLength_FiltersShort()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 1010, StructuralVariantAnalyzer.SVType.Deletion, 10, 50, 5, null),
            new("SV2", "chr1", 5000, 6000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null)
        };

        var filtered = StructuralVariantAnalyzer.FilterSVs(svs, minLength: 100).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
    }

    #endregion

    #region SV Annotation Tests

    [Test]
    public void AnnotateSVs_OverlapsGene_ReturnsAnnotation()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 5000, StructuralVariantAnalyzer.SVType.Deletion, 4000, 50, 5, null)
        };

        var genes = new List<(string, string, int, int, IReadOnlyList<(int, int)>)>
        {
            ("GENE1", "chr1", 2000, 4000, new List<(int, int)> { (2000, 2500), (3000, 3500) })
        };

        var annotations = StructuralVariantAnalyzer.AnnotateSVs(svs, genes).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].AffectedGenes, Contains.Item("GENE1"));
        Assert.That(annotations[0].AffectedExons.Count, Is.GreaterThan(0));
        Assert.That(annotations[0].FunctionalImpact, Is.EqualTo("HIGH"));
    }

    [Test]
    public void AnnotateSVs_NoOverlap_NoAffectedGenes()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null)
        };

        var genes = new List<(string, string, int, int, IReadOnlyList<(int, int)>)>
        {
            ("GENE1", "chr1", 5000, 6000, new List<(int, int)> { (5000, 5500) })
        };

        var annotations = StructuralVariantAnalyzer.AnnotateSVs(svs, genes).ToList();

        Assert.That(annotations[0].AffectedGenes, Is.Empty);
        Assert.That(annotations[0].FunctionalImpact, Is.EqualTo("LOW"));
    }

    [Test]
    public void AnnotateSVs_IntronicDeletion_ModifierImpact()
    {
        var svs = new List<StructuralVariantAnalyzer.StructuralVariant>
        {
            new("SV1", "chr1", 2600, 2900, StructuralVariantAnalyzer.SVType.Deletion, 300, 50, 5, null)
        };

        var genes = new List<(string, string, int, int, IReadOnlyList<(int, int)>)>
        {
            ("GENE1", "chr1", 2000, 4000, new List<(int, int)> { (2000, 2500), (3000, 3500) })
        };

        var annotations = StructuralVariantAnalyzer.AnnotateSVs(svs, genes).ToList();

        Assert.That(annotations[0].AffectedGenes, Contains.Item("GENE1"));
        Assert.That(annotations[0].AffectedExons, Is.Empty);
        Assert.That(annotations[0].FunctionalImpact, Is.EqualTo("MODIFIER"));
    }

    #endregion

    #region SV Genotyping Tests

    [Test]
    public void GenotypeSV_HighAltFraction_HomozygousAlt()
    {
        var sv = new StructuralVariantAnalyzer.StructuralVariant(
            "SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null);

        var (genotype, quality) = StructuralVariantAnalyzer.GenotypeSV(sv, refReads: 2, altReads: 50, totalReads: 52);

        Assert.That(genotype, Is.EqualTo("1/1"));
        Assert.That(quality, Is.GreaterThan(0));
    }

    [Test]
    public void GenotypeSV_LowAltFraction_HomozygousRef()
    {
        var sv = new StructuralVariantAnalyzer.StructuralVariant(
            "SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null);

        var (genotype, quality) = StructuralVariantAnalyzer.GenotypeSV(sv, refReads: 50, altReads: 2, totalReads: 52);

        Assert.That(genotype, Is.EqualTo("0/0"));
    }

    [Test]
    public void GenotypeSV_MidAltFraction_Heterozygous()
    {
        var sv = new StructuralVariantAnalyzer.StructuralVariant(
            "SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null);

        var (genotype, quality) = StructuralVariantAnalyzer.GenotypeSV(sv, refReads: 25, altReads: 25, totalReads: 50);

        Assert.That(genotype, Is.EqualTo("0/1"));
    }

    [Test]
    public void GenotypeSV_NoReads_MissingGenotype()
    {
        var sv = new StructuralVariantAnalyzer.StructuralVariant(
            "SV1", "chr1", 1000, 2000, StructuralVariantAnalyzer.SVType.Deletion, 1000, 50, 5, null);

        var (genotype, quality) = StructuralVariantAnalyzer.GenotypeSV(sv, refReads: 0, altReads: 0, totalReads: 0);

        Assert.That(genotype, Is.EqualTo("./."));
        Assert.That(quality, Is.EqualTo(0));
    }

    #endregion

    #region Breakpoint Assembly Tests

    [Test]
    public void AssembleBreakpointSequence_MultipleSplits_ReturnsLongest()
    {
        var splits = new List<StructuralVariantAnalyzer.SplitRead>
        {
            new("r1", "chr1", 1000, 5000, 20, "ACGTACGTACGTACGTACGT"),
            new("r2", "chr1", 1005, 5005, 50, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
            new("r3", "chr1", 1010, 5010, 30, "ACGTACGTACGTACGTACGTACGTACGTAC")
        };

        var seq = StructuralVariantAnalyzer.AssembleBreakpointSequence(splits);

        Assert.That(seq, Is.Not.Null);
        Assert.That(seq!.Length, Is.EqualTo(50));
    }

    [Test]
    public void AssembleBreakpointSequence_EmptyInput_ReturnsNull()
    {
        var seq = StructuralVariantAnalyzer.AssembleBreakpointSequence(
            new List<StructuralVariantAnalyzer.SplitRead>());

        Assert.That(seq, Is.Null);
    }

    [Test]
    public void FindMicrohomology_ExactMatch_ReturnsLength()
    {
        string leftFlank = "ACGTACGTACGT";
        string rightFlank = "ACGTNNNNNNNN";

        var (length, seq) = StructuralVariantAnalyzer.FindMicrohomology(leftFlank, rightFlank);

        Assert.That(length, Is.EqualTo(4)); // ACGT
        Assert.That(seq, Is.EqualTo("ACGT"));
    }

    [Test]
    public void FindMicrohomology_NoMatch_ReturnsZero()
    {
        string leftFlank = "AAAAAAAAAA";
        string rightFlank = "TTTTTTTTTT";

        var (length, seq) = StructuralVariantAnalyzer.FindMicrohomology(leftFlank, rightFlank);

        Assert.That(length, Is.EqualTo(0));
        Assert.That(seq, Is.Empty);
    }

    [Test]
    public void FindMicrohomology_EmptyInput_ReturnsZero()
    {
        var (length, seq) = StructuralVariantAnalyzer.FindMicrohomology("", "ACGT");

        Assert.That(length, Is.EqualTo(0));
    }

    #endregion

    #region Edge Cases Tests

    [Test]
    public void FindDiscordantPairs_EmptyInput_ReturnsEmpty()
    {
        var discordant = StructuralVariantAnalyzer.FindDiscordantPairs(
            new List<(string, string, int, char, string, int, char, int)>()).ToList();

        Assert.That(discordant, Is.Empty);
    }

    [Test]
    public void ClusterDiscordantPairs_EmptyInput_ReturnsEmpty()
    {
        var svs = StructuralVariantAnalyzer.ClusterDiscordantPairs(
            new List<StructuralVariantAnalyzer.ReadPairSignature>()).ToList();

        Assert.That(svs, Is.Empty);
    }

    [Test]
    public void SegmentCopyNumber_EmptyInput_ReturnsEmpty()
    {
        var segments = StructuralVariantAnalyzer.SegmentCopyNumber(
            new List<(string, int, double, double)>()).ToList();

        Assert.That(segments, Is.Empty);
    }

    [Test]
    public void MergeOverlappingSVs_EmptyInput_ReturnsEmpty()
    {
        var merged = StructuralVariantAnalyzer.MergeOverlappingSVs(
            new List<StructuralVariantAnalyzer.StructuralVariant>()).ToList();

        Assert.That(merged, Is.Empty);
    }

    [Test]
    public void ClusterSplitReads_EmptyInput_ReturnsEmpty()
    {
        var breakpoints = StructuralVariantAnalyzer.ClusterSplitReads(
            new List<StructuralVariantAnalyzer.SplitRead>()).ToList();

        Assert.That(breakpoints, Is.Empty);
    }

    #endregion
}
