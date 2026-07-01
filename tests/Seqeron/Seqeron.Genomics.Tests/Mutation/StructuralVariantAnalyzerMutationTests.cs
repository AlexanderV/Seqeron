using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;
using SVType = Seqeron.Genomics.Annotation.StructuralVariantAnalyzer.SVType;
using StructuralVariant = Seqeron.Genomics.Annotation.StructuralVariantAnalyzer.StructuralVariant;
using ReadPairSignature = Seqeron.Genomics.Annotation.StructuralVariantAnalyzer.ReadPairSignature;
using SplitRead = Seqeron.Genomics.Annotation.StructuralVariantAnalyzer.SplitRead;
using CopyNumberSegment = Seqeron.Genomics.Annotation.StructuralVariantAnalyzer.CopyNumberSegment;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// SV-DETECT/BREAKPOINT/CNV-001 mutation killers: exact-value / boundary tests for the deterministic
/// SV helpers the canonical fixtures only smoke-tested — PEM signature classification (Medvedev 2009),
/// SAM CIGAR soft-clip parsing + aligned-length, reciprocal-overlap merging, quality/support filtering,
/// read-depth CNV calling from segments, gene/exon annotation + impact, VAF genotyping, and breakpoint
/// microhomology.
/// </summary>
[TestFixture]
public class StructuralVariantAnalyzerMutationTests
{
    private const double Tol = 1e-9;

    #region ClassifySV — PEM signatures (Medvedev et al. 2009)

    private static ReadPairSignature Pair(string chr1, string chr2, char s1, char s2, int insertSize)
        => new("r", chr1, 1000, s1, chr2, 2000, s2, insertSize, true);

    [Test]
    public void ClassifySV_InterchromosomalIsTranslocation()
        => Assert.That(ClassifySV(Pair("chr1", "chr2", '+', '-', 400)), Is.EqualTo(SVType.Translocation));

    [Test]
    public void ClassifySV_SameStrandIsInversion()
        => Assert.That(ClassifySV(Pair("chr1", "chr1", '+', '+', 400)), Is.EqualTo(SVType.Inversion));

    [Test]
    public void ClassifySV_ReverseForwardIsDuplication()
        => Assert.That(ClassifySV(Pair("chr1", "chr1", '-', '+', 400)), Is.EqualTo(SVType.Duplication));

    // FR mates (concordant orientation), default bounds = 400 ± 3·50 = [250, 550].
    [TestCase(600, SVType.Deletion)]              // span > upper ⇒ deletion
    [TestCase(200, SVType.Insertion)]             // span < lower ⇒ insertion
    [TestCase(551, SVType.Deletion)]              // just above upper
    [TestCase(550, SVType.ComplexRearrangement)]  // == upper ⇒ NOT deletion (kills > → >=)
    [TestCase(250, SVType.ComplexRearrangement)]  // == lower ⇒ NOT insertion (kills < → <=)
    [TestCase(300, SVType.ComplexRearrangement)]  // inside ⇒ complex (kills upper→lower / cutoff/sd mutants)
    [TestCase(400, SVType.ComplexRearrangement)]  // inside (kills lower→expected+cutoff·sd)
    [TestCase(500, SVType.ComplexRearrangement)]  // inside (kills upper→cutoff/sd)
    public void ClassifySV_ForwardReverseSpanThresholds(int insertSize, SVType expected)
        => Assert.That(ClassifySV(Pair("chr1", "chr1", '+', '-', insertSize)), Is.EqualTo(expected));

    #endregion

    #region GenotypeSV — VAF-based genotype + quality

    [TestCase(20, 1, 100, "0/0", 60.0)]   // AF 0.01 < 0.1 ⇒ hom-ref, q = refReads·3
    [TestCase(1, 30, 31, "1/1", 90.0)]    // AF 0.97 > 0.9 ⇒ hom-alt, q = altReads·3
    [TestCase(5, 5, 10, "0/1", 20.0)]     // AF 0.5 ∈ [0.3,0.7] ⇒ het, q = (ref+alt)·2
    [TestCase(7, 3, 10, "0/1", 20.0)]     // AF 0.3 boundary ⇒ het (kills >= → >)
    [TestCase(3, 7, 10, "0/1", 20.0)]     // AF 0.7 boundary ⇒ het (kills <= → <)
    [TestCase(8, 2, 10, "0/1", 15.0)]     // AF 0.2 ⇒ likely-het else branch, q = (ref+alt)·1.5
    [TestCase(90, 10, 100, "0/1", 99.0)]  // AF 0.1 boundary ⇒ NOT hom-ref (kills < → <=)
    [TestCase(10, 90, 100, "0/1", 99.0)]  // AF 0.9 boundary ⇒ NOT hom-alt (kills > → >=)
    public void GenotypeSV_VafThresholds(int refReads, int altReads, int total, string gt, double q)
    {
        var sv = MakeSV("SV1", "chr1", 100, 200, SVType.Deletion, 100, 60, 5);
        var (genotype, quality) = GenotypeSV(sv, refReads, altReads, total);
        Assert.That(genotype, Is.EqualTo(gt));
        Assert.That(quality, Is.EqualTo(q).Within(Tol));
    }

    [Test]
    public void GenotypeSV_NoReadsIsNoCall()
    {
        var (genotype, quality) = GenotypeSV(MakeSV("SV1", "chr1", 100, 200, SVType.Deletion, 100, 60, 5), 5, 5, 0);
        Assert.That(genotype, Is.EqualTo("./."));
        Assert.That(quality, Is.EqualTo(0).Within(Tol));
    }

    #endregion

    #region FilterSVs — inclusive quality/support/length boundaries

    [Test]
    public void FilterSVs_KeepsVariantsExactlyAtEachBoundary()
    {
        var svs = new[]
        {
            MakeSV("Q",   "chr1", 0, 100, SVType.Deletion, 100, 20, 5),           // quality == minQuality (20)
            MakeSV("S",   "chr1", 0, 100, SVType.Deletion, 100, 50, 2),           // support == minSupport (2)
            MakeSV("L",   "chr1", 0, 100, SVType.Deletion, 50, 50, 5),            // length == minLength (50)
            MakeSV("MAX", "chr1", 0, 100, SVType.Deletion, 100_000_000, 50, 5),   // length == maxLength
            MakeSV("LowQ","chr1", 0, 100, SVType.Deletion, 100, 19, 5),           // quality below min ⇒ dropped
        };
        var kept = FilterSVs(svs).Select(v => v.Id).ToList();
        Assert.That(kept, Is.EquivalentTo(new[] { "Q", "S", "L", "MAX" }));
    }

    #endregion

    #region FindSplitReads / CIGAR parsing

    [Test]
    public void FindSplitReads_LeftSoftClipUsesPrimaryPositionAndPrefixSequence()
    {
        string seq = new string('C', 30) + new string('A', 70);
        var sr = FindSplitReads(new[] { ("r1", "chr1", 100, "30S70M", seq) }).Single();
        Assert.That(sr.ClipLength, Is.EqualTo(30));
        Assert.That(sr.SupplementaryPosition, Is.EqualTo(100));          // left clip ⇒ supp = primary
        Assert.That(sr.ClippedSequence, Is.EqualTo(new string('C', 30))); // prefix (kills isLeft pos==0 mutant)
    }

    [Test]
    public void FindSplitReads_RightSoftClipUsesAlignedEndAndSuffixSequence()
    {
        string seq = new string('A', 70) + new string('G', 30);
        var sr = FindSplitReads(new[] { ("r2", "chr1", 100, "70M30S", seq) }).Single();
        Assert.That(sr.ClipLength, Is.EqualTo(30));
        Assert.That(sr.SupplementaryPosition, Is.EqualTo(170));          // 100 + aligned(70)
        Assert.That(sr.ClippedSequence, Is.EqualTo(new string('G', 30))); // suffix
    }

    [TestCase("50M10D30S", 160)] // aligned = 50 + 10 (D consumes reference)
    [TestCase("50M10N30S", 160)] // aligned = 50 + 10 (N consumes reference)
    public void FindSplitReads_DeletionAndSkipConsumeReferenceInAlignedLength(string cigar, int suppPos)
    {
        string seq = new string('A', 50) + new string('G', 30);
        var sr = FindSplitReads(new[] { ("r", "chr1", 100, cigar, seq) }).Single();
        Assert.That(sr.SupplementaryPosition, Is.EqualTo(suppPos));
    }

    [Test]
    public void FindSplitReads_ClipExactlyAtMinLengthIsReported()
    {
        string seq = new string('C', 20) + new string('A', 80);
        var reads = FindSplitReads(new[] { ("r", "chr1", 100, "20S80M", seq) }, minClipLength: 20).ToList();
        Assert.That(reads, Has.Count.EqualTo(1)); // 20 >= 20 (kills >= → >)
    }

    [Test]
    public void FindSplitReads_ClipBelowMinLengthIsDropped()
    {
        string seq = new string('C', 19) + new string('A', 81);
        Assert.That(FindSplitReads(new[] { ("r", "chr1", 100, "19S81M", seq) }, minClipLength: 20), Is.Empty);
    }

    #endregion

    #region MergeOverlappingSVs / reciprocal overlap

    [Test]
    public void MergeOverlappingSVs_MergesAtExactlyHalfReciprocalOverlap()
    {
        var a = MakeSV("A", "chr1", 100, 200, SVType.Deletion, 100, 60, 3, inserted: "AAA");
        var b = MakeSV("B", "chr1", 150, 250, SVType.Deletion, 100, 80, 4, inserted: "GGG");
        // overlap = 50 / min(100,100) = 0.5 == fraction ⇒ merge (kills >= → >).
        var merged = MergeOverlappingSVs(new[] { a, b }, overlapFraction: 0.5).ToList();

        Assert.That(merged, Has.Count.EqualTo(1));
        Assert.That(merged[0].Start, Is.EqualTo(100));
        Assert.That(merged[0].End, Is.EqualTo(250));
        Assert.That(merged[0].Length, Is.EqualTo(150));            // max(End) − min(Start)
        Assert.That(merged[0].SupportingReads, Is.EqualTo(7));     // 3 + 4
        Assert.That(merged[0].Quality, Is.EqualTo(80).Within(Tol)); // max(60,80)
        Assert.That(merged[0].InsertedSequence, Is.EqualTo("AAA")); // current ?? next (kills coalescing mutants)
    }

    [Test]
    public void MergeOverlappingSVs_DoesNotMergeBelowOverlapFraction()
    {
        var a = MakeSV("A", "chr1", 100, 200, SVType.Deletion, 100, 60, 3);
        var b = MakeSV("B", "chr1", 190, 290, SVType.Deletion, 100, 80, 4);
        // overlap = 10 / 100 = 0.1 < 0.5 ⇒ kept separate (kills overlapLen arithmetic mutants).
        Assert.That(MergeOverlappingSVs(new[] { a, b }, overlapFraction: 0.5).Count(), Is.EqualTo(2));
    }

    #endregion

    #region IdentifyCNVs

    [Test]
    public void IdentifyCNVs_CallsLongNonDiploidSegmentsWithTypeLengthQuality()
    {
        var segments = new[]
        {
            Seg("chr1", 0, 20000, -1.0, 1),      // deletion, length 20000, q = |−1|·50 = 50
            Seg("chr1", 0, 20000, 0.0, 2),       // diploid ⇒ skipped
            Seg("chr1", 5000, 25000, 0.585, 3),  // duplication, length 20000
            Seg("chr1", 5000, 12000, -1.0, 1),   // length 7000 < 10000 ⇒ skipped
        };
        var cnvs = IdentifyCNVs(segments).ToList();

        Assert.That(cnvs, Has.Count.EqualTo(2));
        var del = cnvs.Single(c => c.Type == SVType.Deletion);
        Assert.That(del.Length, Is.EqualTo(20000));        // End − Start (kills End+Start)
        Assert.That(del.Quality, Is.EqualTo(50).Within(Tol)); // |LogRatio|·50 (kills /50)
        Assert.That(cnvs.Any(c => c.Type == SVType.Duplication), Is.True);
    }

    #endregion

    #region AnnotateSVs — gene/exon overlap + impact

    [Test]
    public void AnnotateSVs_FlagsAffectedExonsAndHighImpact()
    {
        var sv = MakeSV("SV1", "chr1", 150, 160, SVType.Deletion, 10, 60, 5);
        var genes = new (string, string, int, int, IReadOnlyList<(int, int)>)[]
        {
            ("g1", "chr1", 100, 200, new[] { (100, 120), (155, 165), (180, 200) }), // exon2 hit
            ("g2", "chr1", 300, 400, new[] { (300, 320) }),                          // no overlap
        };
        var ann = AnnotateSVs(new[] { sv }, genes).Single();

        Assert.That(ann.AffectedGenes, Is.EquivalentTo(new[] { "g1" }));      // g2 excluded (kills overlap boundary)
        Assert.That(ann.AffectedExons, Is.EquivalentTo(new[] { "g1:exon2" })); // 1-based index (kills i+1 → i-1)
        Assert.That(ann.FunctionalImpact, Is.EqualTo("HIGH"));                 // deletion hitting exon
        Assert.That(ann.IsPathogenic, Is.True);                               // HIGH ⇒ pathogenic (kills || → &&)
    }

    [Test]
    public void AnnotateSVs_GeneOverlapWithoutExonIsModifierAndNotPathogenic()
    {
        var sv = MakeSV("SV2", "chr1", 150, 160, SVType.Deletion, 10, 60, 5);
        var genes = new (string, string, int, int, IReadOnlyList<(int, int)>)[]
        {
            ("g3", "chr1", 140, 400, new[] { (300, 320) }), // overlaps gene body, misses every exon
        };
        var ann = AnnotateSVs(new[] { sv }, genes).Single();
        Assert.That(ann.AffectedGenes, Is.EquivalentTo(new[] { "g3" }));
        Assert.That(ann.AffectedExons, Is.Empty);
        Assert.That(ann.FunctionalImpact, Is.EqualTo("MODIFIER"));
        Assert.That(ann.IsPathogenic, Is.False);
    }

    #endregion

    #region FindMicrohomology + CreateBreakpoint quality

    [Test]
    public void FindMicrohomology_FindsSingleBaseHomology()
    {
        // Longest suffix of left == prefix of right is the single base "T" (kills len >= 1 → len > 1).
        var (length, seq) = FindMicrohomology("ACGT", "TGCA");
        Assert.That(length, Is.EqualTo(1));
        Assert.That(seq, Is.EqualTo("T"));
    }

    [Test]
    public void FindMicrohomology_EmptyFlankYieldsNone()
        => Assert.That(FindMicrohomology("", "ACGT"), Is.EqualTo((0, "")));

    [Test]
    public void FindBreakpoints_ClustersJunctionsWithMeanPositionAndQuality()
    {
        var reads = new[]
        {
            new SplitRead("r1", "chr1", 90, 1000, 25, "AAAA"),
            new SplitRead("r2", "chr1", 92, 1002, 25, "CCCC"),
        };
        var bp = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2).Single();
        Assert.That(bp.SupportingReads, Is.EqualTo(2));
        Assert.That(bp.Position1, Is.EqualTo(1001));                 // round(mean(1000,1002))
        Assert.That(bp.Quality, Is.EqualTo(30).Within(Tol));        // min(2·15, 100)
    }

    #endregion

    private static StructuralVariant MakeSV(
        string id, string chr, int start, int end, SVType type, int length, double quality, int support,
        string? inserted = null)
        => new(id, chr, start, end, type, length, quality, support, inserted);

    private static CopyNumberSegment Seg(string chr, int start, int end, double logR, int cn)
        => new(chr, start, end, logR, cn, 0.5, 10);
}
