using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Annotation.VariantAnnotator;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;
using Transcript = Seqeron.Genomics.Annotation.VariantAnnotator.Transcript;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// VARIANT-ANNOT-001 mutation killers (batch 5): the reverse-strand ('-') consequence dispatch
/// (upstream/downstream/splice-donor/splice-acceptor/splice-region/UTR/start-lost), the exact
/// upstream(5 kb)/downstream(500 bp) flank distances of OverlapsOrNear / DetermineConsequence, and
/// the CDS-position accumulation (forward + reverse) that drives the protein position and codon
/// number. Conforms to Ensembl VEP <c>VariationEffect.pm</c> consequence predicates: on the minus
/// strand the splice donor sits at the exon 5' (genomic-low) edge and the acceptor at the 3'
/// (genomic-high) edge, and the start codon coincides with the genomic-high CDS bound (CdsEnd).
/// </summary>
[TestFixture]
public class VariantAnnotator_MutationKillers5_Tests
{
    private const double Tol = 1e-9;

    // exon1 100..150, intron 151..199, exon2 200..250 on the MINUS strand.
    // CDS 120..240 ⇒ coding exons (120,150) and (200,240); start codon is at CdsEnd (240) on '-'.
    // ⇒ 3'UTR 100..119 (genomic-low end), 5'UTR 241..250 (genomic-high end).
    private static Transcript MinusTwoExon() => new(
        "ENST_M", "ENSG_M", "GENE_M", "chr1",
        100, 250, '-',
        new List<(int, int)> { (100, 150), (200, 250) },
        new List<(int, int)> { (120, 150), (200, 240) },
        120, 240);

    private static readonly string RefSeq = new string('A', 300);

    private static ConsequenceType MinusConsequence(int pos, string r, string a, VariantType t) =>
        AnnotateVariant(new Variant("chr1", pos, r, a, t), new[] { MinusTwoExon() }, RefSeq).Single().Consequence;

    #region Minus-strand consequence dispatch

    [Test]
    public void Minus_UpstreamGeneVariant_PastGenomicHighEnd()
        // '-' strand: 5' is the genomic-high side. A variant just past End (250) within 5 kb is upstream.
        => Assert.That(MinusConsequence(300, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.UpstreamGeneVariant));

    [Test]
    public void Minus_DownstreamGeneVariant_BeforeGenomicLowEnd()
        // '-' strand: 3' is the genomic-low side. A variant just before Start (100) within 500 bp is downstream.
        => Assert.That(MinusConsequence(60, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.DownstreamGeneVariant));

    [Test]
    public void Minus_SpliceAcceptor_AtGenomicHighExonEdge()
        // '-' acceptor = 2 bp into the intron at the exon 3'/genomic-high edge (exon1 end 150 → 151..152).
        => Assert.That(MinusConsequence(151, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceAcceptorVariant));

    [Test]
    public void Minus_SpliceDonor_AtGenomicLowExonEdge()
        // '-' donor = 2 bp into the intron at the exon 5'/genomic-low edge (exon2 start 200 → 198..199).
        => Assert.That(MinusConsequence(199, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceDonorVariant));

    [Test]
    public void Minus_SpliceRegion_3To8BpIntoIntron()
        // 3..8 bp into the intron from exon1's genomic-high edge (153..158).
        => Assert.That(MinusConsequence(155, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceRegionVariant));

    [Test]
    public void Minus_ThreePrimeUtr_BelowCdsAtGenomicLowEnd()
        // In exon1 but below CdsStart (120); on '-' the genomic-low UTR is the 3' UTR.
        => Assert.That(MinusConsequence(110, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.ThreePrimeUtrVariant));

    [Test]
    public void Minus_FivePrimeUtr_AboveCdsAtGenomicHighEnd()
        // In exon2 but above CdsEnd (240); on '-' the genomic-high UTR is the 5' UTR.
        => Assert.That(MinusConsequence(245, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.FivePrimeUtrVariant));

    [Test]
    public void Minus_StartLost_AtCdsEndStartCodon()
        // The ATG start codon lies at the genomic-high CDS bound (CdsEnd 240) on '-'; position 240 is within [238,240].
        => Assert.That(MinusConsequence(240, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.StartLost));

    [Test]
    public void Minus_CodingMissense_DefaultWhenNotStartCodon()
        => Assert.That(MinusConsequence(210, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.MissenseVariant));

    #endregion

    #region Exact flank distances (upstream 5 kb / downstream 500 bp)

    // Single-exon '+' transcript far from the origin so flank boundaries use positive coordinates.
    private static Transcript PlusFar() => new(
        "ENST_F", "ENSG_F", "GENE_F", "chr1",
        10000, 10100, '+',
        new List<(int, int)> { (10000, 10100) },
        new List<(int, int)>(), null, null);

    private static Transcript MinusFar() => new(
        "ENST_MF", "ENSG_MF", "GENE_MF", "chr1",
        10000, 10100, '-',
        new List<(int, int)> { (10000, 10100) },
        new List<(int, int)>(), null, null);

    [Test]
    public void Plus_UpstreamBoundary_Exactly5000BpIsUpstreamNotIntergenic()
    {
        // Start - varEnd == 5000 ⇒ still upstream_gene_variant (inclusive ≤ 5000), and the variant
        // must survive OverlapsOrNear's identical 5 kb upstream flank.
        var c = AnnotateVariant(new Variant("chr1", 5000, "A", "G", VariantType.SNV), new[] { PlusFar() }).Single().Consequence;
        Assert.That(c, Is.EqualTo(ConsequenceType.UpstreamGeneVariant));
    }

    [Test]
    public void Plus_DownstreamBoundary_Exactly500BpIsDownstreamNotIntergenic()
    {
        // varStart - End == 500 ⇒ still downstream_gene_variant (inclusive ≤ 500).
        var c = AnnotateVariant(new Variant("chr1", 10600, "A", "G", VariantType.SNV), new[] { PlusFar() }).Single().Consequence;
        Assert.That(c, Is.EqualTo(ConsequenceType.DownstreamGeneVariant));
    }

    [Test]
    public void Minus_DownstreamBoundary_Exactly500BpIsDownstreamNotIntergenic()
    {
        // '-' downstream: Start - varEnd == 500 ⇒ still downstream_gene_variant.
        var c = AnnotateVariant(new Variant("chr1", 9500, "A", "G", VariantType.SNV), new[] { MinusFar() }).Single().Consequence;
        Assert.That(c, Is.EqualTo(ConsequenceType.DownstreamGeneVariant));
    }

    #endregion

    #region CDS-position / protein-position accumulation

    // exon1 100..150, exon2 200..250; CDS 120..240 ('+'). Coding exons (120,150),(200,240).
    private static Transcript PlusTwoExon() => new(
        "ENST_P2", "ENSG_P2", "GENE_P2", "chr1",
        100, 250, '+',
        new List<(int, int)> { (100, 150), (200, 250) },
        new List<(int, int)> { (120, 150), (200, 240) },
        120, 240);

    [Test]
    public void Plus_CdsPosition_FirstCodingExon()
    {
        // Forward CDS: pos 135 in coding exon1 (120..150) ⇒ cds = (135-120)+1 = 16, protein = (16-1)/3+1 = 6.
        var a = AnnotateVariant(new Variant("chr1", 135, "A", "G", VariantType.SNV), new[] { PlusTwoExon() }, RefSeq).Single();
        Assert.That(a.CdsPosition, Is.EqualTo(16));
        Assert.That(a.ProteinPosition, Is.EqualTo(6));
        Assert.That(a.AminoAcidChange, Is.EqualTo("p.X6Y"));
    }

    [Test]
    public void Plus_CdsPosition_SecondCodingExonAccumulatesFirstExonLength()
    {
        // exon1 contributes (150-120+1)=31 bases; pos 210 in exon2 ⇒ cds = 31+(210-200)+1 = 42, protein 14.
        var a = AnnotateVariant(new Variant("chr1", 210, "A", "G", VariantType.SNV), new[] { PlusTwoExon() }, RefSeq).Single();
        Assert.That(a.CdsPosition, Is.EqualTo(42));
        Assert.That(a.ProteinPosition, Is.EqualTo(14));
        Assert.That(a.AminoAcidChange, Is.EqualTo("p.X14Y"));
        // Missense ⇒ SIFT/PolyPhen from the simplified biochemical model on the (X→Y) placeholder pair:
        // 'X' is in no biochemical group and 'Y' is aromatic ⇒ dissimilar ⇒ similarity 0.3.
        Assert.That(a.SiftScore!.Value, Is.EqualTo(0.3).Within(Tol));
        Assert.That(a.PolyphenScore!.Value, Is.EqualTo(0.7).Within(Tol));
    }

    [Test]
    public void Minus_CdsPosition_FirstCodingExonFromGenomicHighEnd()
    {
        // Reverse CDS walks from the genomic-high coding exon: pos 210 in (200,240) ⇒ cds = (240-210)+1 = 31, protein 11.
        var a = AnnotateVariant(new Variant("chr1", 210, "A", "G", VariantType.SNV), new[] { MinusTwoExon() }, RefSeq).Single();
        Assert.That(a.CdsPosition, Is.EqualTo(31));
        Assert.That(a.ProteinPosition, Is.EqualTo(11));
        Assert.That(a.AminoAcidChange, Is.EqualTo("p.X11Y"));
    }

    [Test]
    public void Minus_CdsPosition_SecondCodingExonAccumulatesGenomicHighExonLength()
    {
        // Genomic-high coding exon (200,240) contributes (240-200+1)=41; pos 145 in (120,150)
        // ⇒ cds = 41+(150-145)+1 = 47, protein = (47-1)/3+1 = 16.
        var a = AnnotateVariant(new Variant("chr1", 145, "A", "G", VariantType.SNV), new[] { MinusTwoExon() }, RefSeq).Single();
        Assert.That(a.CdsPosition, Is.EqualTo(47));
        Assert.That(a.ProteinPosition, Is.EqualTo(16));
    }

    #endregion
}
