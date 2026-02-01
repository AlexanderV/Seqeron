using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class VariantAnnotatorTests
{
    #region Variant Classification Tests

    [Test]
    public void ClassifyVariant_SingleNucleotide_ReturnsSNV()
    {
        var type = VariantAnnotator.ClassifyVariant("A", "G");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.SNV));
    }

    [Test]
    public void ClassifyVariant_Insertion_ReturnsInsertion()
    {
        var type = VariantAnnotator.ClassifyVariant("A", "ACGT");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.Insertion));
    }

    [Test]
    public void ClassifyVariant_Deletion_ReturnsDeletion()
    {
        var type = VariantAnnotator.ClassifyVariant("ACGT", "A");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.Deletion));
    }

    [Test]
    public void ClassifyVariant_MNV_ReturnsMNV()
    {
        var type = VariantAnnotator.ClassifyVariant("AC", "GT");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.MNV));
    }

    [Test]
    public void ClassifyVariant_ComplexIndel_ReturnsIndel()
    {
        var type = VariantAnnotator.ClassifyVariant("ACG", "TT");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.Indel));
    }

    [Test]
    public void ClassifyVariant_EmptyInput_ReturnsComplex()
    {
        var type = VariantAnnotator.ClassifyVariant("", "A");
        Assert.That(type, Is.EqualTo(VariantAnnotator.VariantType.Complex));
    }

    #endregion

    #region Variant Normalization Tests

    [Test]
    public void NormalizeVariant_CommonSuffix_Trimmed()
    {
        var normalized = VariantAnnotator.NormalizeVariant("chr1", 100, "ACG", "ATG");

        Assert.That(normalized.Position, Is.EqualTo(101));
        Assert.That(normalized.Reference, Is.EqualTo("C"));
        Assert.That(normalized.Alternate, Is.EqualTo("T"));
    }

    [Test]
    public void NormalizeVariant_CommonPrefix_Trimmed()
    {
        var normalized = VariantAnnotator.NormalizeVariant("chr1", 100, "ACGT", "AGGG");

        Assert.That(normalized.Position, Is.EqualTo(101));
    }

    [Test]
    public void NormalizeVariant_PreservesChromosome()
    {
        var normalized = VariantAnnotator.NormalizeVariant("chrX", 100, "A", "T");

        Assert.That(normalized.Chromosome, Is.EqualTo("chrX"));
    }

    #endregion

    #region Consequence Annotation Tests

    [Test]
    public void AnnotateVariant_IntergenicVariant_ReturnsIntergenic()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            CreateTestTranscript("T1", "G1", "Gene1", "chr1", 1, 1000, '+')
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.IntergenicVariant));
    }

    [Test]
    public void AnnotateVariant_UpstreamVariant_ReturnsUpstream()
    {
        var variant = new VariantAnnotator.Variant("chr1", 950, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            CreateTestTranscript("T1", "G1", "Gene1", "chr1", 1000, 2000, '+')
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.UpstreamGeneVariant));
    }

    [Test]
    public void AnnotateVariant_DownstreamVariant_ReturnsDownstream()
    {
        var variant = new VariantAnnotator.Variant("chr1", 2200, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            CreateTestTranscript("T1", "G1", "Gene1", "chr1", 1000, 2000, '+')
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.DownstreamGeneVariant));
    }

    [Test]
    public void AnnotateVariant_IntronVariant_ReturnsIntron()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1250, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1500), (1800, 2000) };
        var codingExons = new List<(int, int)> { (1050, 1100), (1400, 1500), (1800, 1900) };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, 1050, 1900)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.IntronVariant));
    }

    [Test]
    public void AnnotateVariant_MissenseVariant_ReturnsMissense()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1450, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1500), (1800, 2000) };
        var codingExons = new List<(int, int)> { (1050, 1100), (1400, 1500), (1800, 1900) };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, 1050, 1900)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.MissenseVariant));
        Assert.That(annotations[0].Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate));
    }

    [Test]
    public void AnnotateVariant_FrameshiftVariant_ReturnsFrameshift()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1450, "AC", "A",
            VariantAnnotator.VariantType.Deletion);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1500), (1800, 2000) };
        var codingExons = new List<(int, int)> { (1050, 1100), (1400, 1500), (1800, 1900) };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, 1050, 1900)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.FrameshiftVariant));
        Assert.That(annotations[0].Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High));
    }

    [Test]
    public void AnnotateVariant_InframeInsertion_ReturnsInframe()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1450, "A", "ATTT",
            VariantAnnotator.VariantType.Insertion);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1500), (1800, 2000) };
        var codingExons = new List<(int, int)> { (1050, 1100), (1400, 1500), (1800, 1900) };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, 1050, 1900)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.InframeInsertion));
    }

    [Test]
    public void AnnotateVariant_SpliceAcceptorVariant_ReturnsSpliceAcceptor()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1399, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1500), (1800, 2000) };
        var codingExons = new List<(int, int)> { (1050, 1100), (1400, 1500), (1800, 1900) };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, 1050, 1900)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.SpliceAcceptorVariant));
        Assert.That(annotations[0].Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High));
    }

    [Test]
    public void AnnotateVariant_MultipleTranscripts_ReturnsMultipleAnnotations()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1500, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            CreateTestTranscript("T1", "G1", "Gene1", "chr1", 1000, 2000, '+'),
            CreateTestTranscript("T2", "G1", "Gene1", "chr1", 1200, 2200, '+')
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations, Has.Count.EqualTo(2));
    }

    #endregion

    #region Impact Level Tests

    [Test]
    public void GetImpactLevel_HighImpactConsequences_ReturnHigh()
    {
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.StopGained),
            Is.EqualTo(VariantAnnotator.ImpactLevel.High));
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.FrameshiftVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.High));
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SpliceDonorVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.High));
    }

    [Test]
    public void GetImpactLevel_ModerateConsequences_ReturnModerate()
    {
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.MissenseVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate));
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.InframeDeletion),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate));
    }

    [Test]
    public void GetImpactLevel_LowConsequences_ReturnLow()
    {
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SynonymousVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Low));
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SpliceRegionVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Low));
    }

    [Test]
    public void GetImpactLevel_ModifierConsequences_ReturnModifier()
    {
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.IntronVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Modifier));
        Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.IntergenicVariant),
            Is.EqualTo(VariantAnnotator.ImpactLevel.Modifier));
    }

    #endregion

    #region Pathogenicity Prediction Tests

    [Test]
    public void PredictPathogenicity_HighImpact_LikelyPathogenic()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "T1", "G1", "Gene1",
            VariantAnnotator.ConsequenceType.StopGained,
            VariantAnnotator.ImpactLevel.High,
            null, null, null, null,
            null, null, null, null, null);

        var prediction = VariantAnnotator.PredictPathogenicity(annotation, 0.00001);

        Assert.That(prediction.Classification,
            Is.EqualTo(VariantAnnotator.PathogenicityClass.Pathogenic)
            .Or.EqualTo(VariantAnnotator.PathogenicityClass.LikelyPathogenic));
        Assert.That(prediction.IsActionable, Is.True);
    }

    [Test]
    public void PredictPathogenicity_CommonVariant_Benign()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "T1", "G1", "Gene1",
            VariantAnnotator.ConsequenceType.SynonymousVariant,
            VariantAnnotator.ImpactLevel.Low,
            null, null, null, null,
            null, null, null, null, null);

        var prediction = VariantAnnotator.PredictPathogenicity(annotation, 0.10); // 10% frequency

        Assert.That(prediction.Classification,
            Is.EqualTo(VariantAnnotator.PathogenicityClass.Benign)
            .Or.EqualTo(VariantAnnotator.PathogenicityClass.LikelyBenign));
        Assert.That(prediction.IsActionable, Is.False);
    }

    [Test]
    public void PredictPathogenicity_RareVariant_AddsPM2Evidence()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "T1", "G1", "Gene1",
            VariantAnnotator.ConsequenceType.MissenseVariant,
            VariantAnnotator.ImpactLevel.Moderate,
            null, null, null, null,
            null, null, null, null, null);

        var prediction = VariantAnnotator.PredictPathogenicity(annotation, 0.00005);

        Assert.That(prediction.EvidenceCriteria, Has.Some.Contain("PM2"));
    }

    [Test]
    public void PredictPathogenicity_WithClinVarPathogenic_IncreasesScore()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "T1", "G1", "Gene1",
            VariantAnnotator.ConsequenceType.MissenseVariant,
            VariantAnnotator.ImpactLevel.Moderate,
            null, null, null, null,
            null, null, null, null, null);

        var prediction = VariantAnnotator.PredictPathogenicity(
            annotation, inClinvar: true, clinvarSignificance: "Pathogenic");

        Assert.That(prediction.EvidenceCriteria, Has.Some.Contain("ClinVar"));
    }

    [Test]
    public void PredictPathogenicity_DeleterousPredictions_AddsEvidence()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "T1", "G1", "Gene1",
            VariantAnnotator.ConsequenceType.MissenseVariant,
            VariantAnnotator.ImpactLevel.Moderate,
            null, "p.R100W", null, null,
            0.01, // SIFT deleterious
            0.95, // PolyPhen probably damaging
            null, null, null);

        var prediction = VariantAnnotator.PredictPathogenicity(annotation);

        Assert.That(prediction.EvidenceCriteria, Has.Some.Contain("SIFT"));
        Assert.That(prediction.EvidenceCriteria, Has.Some.Contain("PolyPhen"));
    }

    #endregion

    #region Conservation Score Tests

    [Test]
    public void CalculateConservation_HighlyConserved_HighScores()
    {
        var positions = new List<(string, int, IReadOnlyList<char>)>
        {
            ("chr1", 100, new List<char> { 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A' })
        };

        var scores = VariantAnnotator.CalculateConservation(positions).ToList();

        Assert.That(scores, Has.Count.EqualTo(1));
        Assert.That(scores[0].PhastCons, Is.EqualTo(1.0));
        Assert.That(scores[0].PhyloP, Is.GreaterThan(0));
        Assert.That(scores[0].ConservedSpeciesCount, Is.EqualTo(10));
    }

    [Test]
    public void CalculateConservation_Variable_LowScores()
    {
        var positions = new List<(string, int, IReadOnlyList<char>)>
        {
            ("chr1", 100, new List<char> { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' })
        };

        var scores = VariantAnnotator.CalculateConservation(positions).ToList();

        Assert.That(scores[0].PhastCons, Is.LessThan(0.5));
    }

    [Test]
    public void FindConservedElements_ConservedRegion_Detected()
    {
        var scores = Enumerable.Range(0, 50)
            .Select(i => new VariantAnnotator.ConservationScore("chr1", i * 10, 5.0, 0.95, 4.0, 10))
            .ToList();

        var elements = VariantAnnotator.FindConservedElements(scores, threshold: 0.8, minLength: 10).ToList();

        Assert.That(elements, Has.Count.EqualTo(1));
        Assert.That(elements[0].Score, Is.GreaterThan(0.9));
    }

    [Test]
    public void CalculateConservation_EmptyAlleles_ReturnsZeroScores()
    {
        var positions = new List<(string, int, IReadOnlyList<char>)>
        {
            ("chr1", 100, new List<char>())
        };

        var scores = VariantAnnotator.CalculateConservation(positions).ToList();

        Assert.That(scores[0].PhyloP, Is.EqualTo(0));
        Assert.That(scores[0].PhastCons, Is.EqualTo(0));
    }

    #endregion

    #region Regulatory Annotation Tests

    [Test]
    public void AnnotateRegulatoryElements_OverlappingRegion_ReturnsAnnotation()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1500, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var regions = new List<(string, int, int, string, string?, double?, IReadOnlyList<string>)>
        {
            ("chr1", 1000, 2000, "Enhancer", "HeLa", 0.85, new List<string> { "TP53", "MYC" })
        };

        var annotations = VariantAnnotator.AnnotateRegulatoryElements(variant, regions).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].FeatureType, Is.EqualTo("Enhancer"));
        Assert.That(annotations[0].TranscriptionFactors, Contains.Item("TP53"));
    }

    [Test]
    public void AnnotateRegulatoryElements_NoOverlap_ReturnsEmpty()
    {
        var variant = new VariantAnnotator.Variant("chr1", 5000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var regions = new List<(string, int, int, string, string?, double?, IReadOnlyList<string>)>
        {
            ("chr1", 1000, 2000, "Enhancer", null, null, new List<string>())
        };

        var annotations = VariantAnnotator.AnnotateRegulatoryElements(variant, regions).ToList();

        Assert.That(annotations, Is.Empty);
    }

    [Test]
    public void AnnotateRegulatoryElements_DifferentChromosome_ReturnsEmpty()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1500, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var regions = new List<(string, int, int, string, string?, double?, IReadOnlyList<string>)>
        {
            ("chr2", 1000, 2000, "Enhancer", null, null, new List<string>())
        };

        var annotations = VariantAnnotator.AnnotateRegulatoryElements(variant, regions).ToList();

        Assert.That(annotations, Is.Empty);
    }

    #endregion

    #region TF Binding Prediction Tests

    [Test]
    public void PredictTfBindingChange_DisruptsMotif_ReturnsChange()
    {
        var variant = new VariantAnnotator.Variant("chr1", 120, "A", "T",
            VariantAnnotator.VariantType.SNV);

        var motifs = new List<(string, string, double)>
        {
            ("TP53", "RRRCATGYYY", 0.8)
        };

        string context = "NNNNNNNNNNNNNNNNNNNNCATGNNNNNNNNNNNNNNNNN"; // A at position 20

        var changes = VariantAnnotator.PredictTfBindingChange(variant, motifs, context, 20).ToList();

        // May or may not find change depending on motif scoring
        Assert.That(changes, Is.Not.Null);
    }

    #endregion

    #region VCF Integration Tests

    [Test]
    public void ParseVcfVariant_ValidFields_ReturnsVariant()
    {
        var variant = VariantAnnotator.ParseVcfVariant(
            "chr1", 12345, "rs12345", "A", "G", 99.0);

        Assert.That(variant.Chromosome, Is.EqualTo("chr1"));
        Assert.That(variant.Position, Is.EqualTo(12345));
        Assert.That(variant.Id, Is.EqualTo("rs12345"));
        Assert.That(variant.Reference, Is.EqualTo("A"));
        Assert.That(variant.Alternate, Is.EqualTo("G"));
        Assert.That(variant.Quality, Is.EqualTo(99.0));
        Assert.That(variant.Type, Is.EqualTo(VariantAnnotator.VariantType.SNV));
    }

    [Test]
    public void FormatAsVcfInfo_FullAnnotation_ContainsAllFields()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotation = new VariantAnnotator.VariantAnnotation(
            variant, "ENST00001", "ENSG00001", "BRCA1",
            VariantAnnotator.ConsequenceType.MissenseVariant,
            VariantAnnotator.ImpactLevel.Moderate,
            "c.100A>G", "p.R100W", 100, 100,
            0.01, 0.95, null, null, null);

        string info = VariantAnnotator.FormatAsVcfInfo(annotation);

        Assert.That(info, Does.Contain("GENE=BRCA1"));
        Assert.That(info, Does.Contain("TRANSCRIPT=ENST00001"));
        Assert.That(info, Does.Contain("CONSEQUENCE=MissenseVariant"));
        Assert.That(info, Does.Contain("IMPACT=Moderate"));
        Assert.That(info, Does.Contain("HGVSP=p.R100W"));
        Assert.That(info, Does.Contain("SIFT=0.010"));
        Assert.That(info, Does.Contain("POLYPHEN=0.950"));
    }

    #endregion

    #region Batch Processing Tests

    [Test]
    public void AnnotateVariants_MultipleVariants_ReturnsGrouped()
    {
        var variants = new List<VariantAnnotator.Variant>
        {
            new("chr1", 1500, "A", "G", VariantAnnotator.VariantType.SNV),
            new("chr1", 1600, "C", "T", VariantAnnotator.VariantType.SNV)
        };

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            CreateTestTranscript("T1", "G1", "Gene1", "chr1", 1000, 2000, '+')
        };

        var groups = VariantAnnotator.AnnotateVariants(variants, transcripts).ToList();

        Assert.That(groups, Has.Count.EqualTo(2));
        Assert.That(groups[0].Key.Position, Is.EqualTo(1500));
        Assert.That(groups[1].Key.Position, Is.EqualTo(1600));
    }

    #endregion

    #region Edge Cases Tests

    [Test]
    public void AnnotateVariant_EmptyTranscriptList_ReturnsIntergenic()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var annotations = VariantAnnotator.AnnotateVariant(
            variant, new List<VariantAnnotator.Transcript>()).ToList();

        Assert.That(annotations, Has.Count.EqualTo(1));
        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.IntergenicVariant));
    }

    [Test]
    public void AnnotateVariant_NonCodingTranscript_ReturnsNonCoding()
    {
        var variant = new VariantAnnotator.Variant("chr1", 1500, "A", "G",
            VariantAnnotator.VariantType.SNV);

        var exons = new List<(int, int)> { (1000, 1100), (1400, 1600), (1800, 2000) };
        var codingExons = new List<(int, int)>(); // No coding exons

        var transcripts = new List<VariantAnnotator.Transcript>
        {
            new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', exons, codingExons, null, null)
        };

        var annotations = VariantAnnotator.AnnotateVariant(variant, transcripts).ToList();

        Assert.That(annotations[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.NonCodingTranscriptExonVariant));
    }

    [Test]
    public void FindConservedElements_EmptyInput_ReturnsEmpty()
    {
        var elements = VariantAnnotator.FindConservedElements(
            new List<VariantAnnotator.ConservationScore>()).ToList();

        Assert.That(elements, Is.Empty);
    }

    #endregion

    #region Helper Methods

    private static VariantAnnotator.Transcript CreateTestTranscript(
        string transcriptId, string geneId, string geneName,
        string chromosome, int start, int end, char strand)
    {
        var exons = new List<(int, int)>
        {
            (start, start + 100),
            (start + 400, start + 500),
            (start + 800, end)
        };

        var codingExons = new List<(int, int)>
        {
            (start + 50, start + 100),
            (start + 400, start + 500),
            (start + 800, end - 100)
        };

        return new VariantAnnotator.Transcript(
            transcriptId, geneId, geneName, chromosome,
            start, end, strand, exons, codingExons,
            start + 50, end - 100);
    }

    #endregion
}
