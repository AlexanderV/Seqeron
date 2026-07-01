using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FormatVcfInfoTests
{
    private static AnnotatorVariantDto Snv() =>
        new("chr17", 43_045_712, "G", "A", "SNV", null, null);

    private static VariantAnnotationDto FullAnnotation() =>
        new(
            Variant: Snv(),
            TranscriptId: "ENST1",
            GeneId: "ENSG1",
            GeneName: "BRCA1",
            Consequence: "MissenseVariant",
            Impact: "Moderate",
            CodonChange: "c.509G>A",
            AminoAcidChange: "p.Arg170His",
            ProteinPosition: 170,
            CdsPosition: 509,
            SiftScore: 0.02,
            PolyphenScore: 0.95,
            CaddScore: null,
            ExistingVariation: null,
            PopulationFrequencies: null);

    [Test]
    public void FormatVcfInfo_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FormatVcfInfo(FullAnnotation()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FormatVcfInfo(null!));
        // Unknown consequence/impact enum strings are rejected by the DTO -> domain mapping.
        var badConsequence = FullAnnotation() with { Consequence = "NotAConsequence" };
        Assert.Throws<ArgumentException>(() => AnnotationTools.FormatVcfInfo(badConsequence));
    }

    [Test]
    public void FormatVcfInfo_Binding_InvokesSuccessfully()
    {
        // VariantAnnotator.FormatAsVcfInfo joins GENE/TRANSCRIPT/CONSEQUENCE/IMPACT and, when present,
        // HGVSP (AminoAcidChange), HGVSC (CodonChange), SIFT (F3), POLYPHEN (F3) with ';'.
        var result = AnnotationTools.FormatVcfInfo(FullAnnotation());

        Assert.That(result.Info, Is.EqualTo(
            "GENE=BRCA1;TRANSCRIPT=ENST1;CONSEQUENCE=MissenseVariant;IMPACT=Moderate;" +
            "HGVSP=p.Arg170His;HGVSC=c.509G>A;SIFT=0.020;POLYPHEN=0.950"));
    }

    [Test]
    public void FormatVcfInfo_OptionalFieldsOmitted_WhenNull()
    {
        // With no HGVS/SIFT/POLYPHEN, only the four mandatory fields are emitted.
        var minimal = FullAnnotation() with
        {
            CodonChange = null,
            AminoAcidChange = null,
            SiftScore = null,
            PolyphenScore = null
        };
        var result = AnnotationTools.FormatVcfInfo(minimal);

        Assert.That(result.Info, Is.EqualTo(
            "GENE=BRCA1;TRANSCRIPT=ENST1;CONSEQUENCE=MissenseVariant;IMPACT=Moderate"));
    }
}
