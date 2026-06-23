using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Annotation.VariantAnnotator;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;
using Transcript = Seqeron.Genomics.Annotation.VariantAnnotator.Transcript;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// VARIANT-ANNOT-001 mutation killers (batch 2): drives the AnnotateVariant consequence-dispatch
/// pipeline (intergenic / upstream / downstream / coding) and the ACMG-style PredictPathogenicity
/// point accumulation and classification.
/// </summary>
[TestFixture]
public class VariantAnnotator_MutationKillers2_Tests
{
    private const double Tol = 1e-9;

    // Single-exon '+' transcript: exon 90..130, CDS 100..120; referenceSequence[0] == genomic 100.
    private const string RefWindow = "ATGCAAGAATTATAANAAGGG";
    private const int RefStart = 100;

    private static Transcript CodingTranscript() => new(
        "ENST_TEST", "ENSG_TEST", "GENE_TEST", "chr1",
        90, 130, '+',
        new List<(int, int)> { (90, 130) },
        new List<(int, int)> { (100, 120) },
        100, 120);

    #region AnnotateVariant dispatch

    [Test]
    public void AnnotateVariant_MissenseInCds_HasModerateImpactAndScores()
    {
        // GAA(Glu, codon3 @106-108) base2 A>T → GTA(Val): missense.
        var v = new Variant("chr1", 107, "A", "T", VariantType.SNV);
        var a = AnnotateVariant(v, new[] { CodingTranscript() }, RefWindow).Single();

        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.MissenseVariant));
        Assert.That(a.Impact, Is.EqualTo(ImpactLevel.Moderate));
        Assert.That(a.AminoAcidChange, Is.Not.Null); // AnnotateVariant aligns refSeq from pos 1
        Assert.That(a.TranscriptId, Is.EqualTo("ENST_TEST"));
        Assert.That(a.SiftScore, Is.Not.Null);     // missense ⇒ SIFT/PolyPhen computed
        Assert.That(a.PolyphenScore, Is.Not.Null);
    }

    [Test]
    public void AnnotateVariant_NoOverlappingTranscript_IsIntergenic()
    {
        var v = new Variant("chr1", 100000, "A", "G", VariantType.SNV);
        var a = AnnotateVariant(v, new[] { CodingTranscript() }, RefWindow).Single();
        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.IntergenicVariant));
        Assert.That(a.Impact, Is.EqualTo(ImpactLevel.Modifier));
        Assert.That(a.TranscriptId, Is.EqualTo(""));
    }

    [Test]
    public void AnnotateVariant_JustUpstream_IsUpstreamGeneVariant()
    {
        // '+' strand, variant before transcript start 90 but within 5 kb ⇒ upstream_gene_variant.
        var v = new Variant("chr1", 88, "A", "G", VariantType.SNV);
        var a = AnnotateVariant(v, new[] { CodingTranscript() }, RefWindow).Single();
        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.UpstreamGeneVariant));
    }

    [Test]
    public void AnnotateVariant_JustDownstream_IsDownstreamGeneVariant()
    {
        // '+' strand, variant past transcript end 130 but within 500 bp ⇒ downstream_gene_variant.
        var v = new Variant("chr1", 135, "A", "G", VariantType.SNV);
        var a = AnnotateVariant(v, new[] { CodingTranscript() }, RefWindow).Single();
        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.DownstreamGeneVariant));
    }

    #endregion

    #region PredictPathogenicity (ACMG point accumulation)

    private static VariantAnnotation Annotation(
        ConsequenceType consequence, ImpactLevel impact,
        double? sift = null, double? polyphen = null) =>
        new(new Variant("chr1", 100, "A", "G", VariantType.SNV),
            "ENST_TEST", "ENSG_TEST", "GENE_TEST",
            consequence, impact, null, null, null, null,
            sift, polyphen, null, null, null);

    [Test]
    public void PredictPathogenicity_HighImpactClinvarPathogenicLof_IsPathogenicActionable()
    {
        // PVS1(8) + PP5 ClinVar pathogenic(4) + PS3 LOF(4) = 16 ⇒ Pathogenic.
        var ann = Annotation(ConsequenceType.StopGained, ImpactLevel.High);
        var p = PredictPathogenicity(ann, inClinvar: true, clinvarSignificance: "Pathogenic",
            functionalEvidence: new[] { "LOF functional assay" });

        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.Pathogenic));
        Assert.That(p.IsActionable, Is.True);
    }

    [Test]
    public void PredictPathogenicity_HighImpactOnly_IsLikelyPathogenic()
    {
        // PVS1(8) alone ⇒ net 8 ⇒ Likely Pathogenic (still actionable).
        var p = PredictPathogenicity(Annotation(ConsequenceType.FrameshiftVariant, ImpactLevel.High));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.LikelyPathogenic));
        Assert.That(p.IsActionable, Is.True);
    }

    [Test]
    public void PredictPathogenicity_CommonSynonymousClinvarBenign_IsBenign()
    {
        // BA1 MAF>5%(4) + BP7 synonymous(1) + BP6 ClinVar benign(4) = benign 9 ⇒ Benign.
        var ann = Annotation(ConsequenceType.SynonymousVariant, ImpactLevel.Low);
        var p = PredictPathogenicity(ann, populationFrequency: 0.1, inClinvar: true,
            clinvarSignificance: "Benign");
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.Benign));
        Assert.That(p.IsActionable, Is.False);
    }

    [Test]
    public void PredictPathogenicity_ModeratelyCommon_IsLikelyBenign()
    {
        // BS1 MAF>1%(2) only ⇒ net -2 ⇒ Likely Benign.
        var p = PredictPathogenicity(Annotation(ConsequenceType.MissenseVariant, ImpactLevel.Moderate),
            populationFrequency: 0.02);
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.LikelyBenign));
    }

    [Test]
    public void PredictPathogenicity_NoEvidence_IsUncertainWithBaselineConfidence()
    {
        var p = PredictPathogenicity(Annotation(ConsequenceType.MissenseVariant, ImpactLevel.Moderate));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.UncertainSignificance));
        Assert.That(p.ConfidenceScore, Is.EqualTo(0.5).Within(Tol)); // total points 0 ⇒ 0.5
    }

    #endregion
}
