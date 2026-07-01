using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictPathogenicityTests
{
    private static AnnotatorVariantDto Snv() => new("chr1", 1000, "A", "G", "SNV", null, null);

    private static VariantAnnotationDto Annotation(string consequence, string impact) => new(
        Snv(), "T1", "G1", "Gene1", consequence, impact,
        null, null, null, null, null, null, null, null, null);

    [Test]
    public void PredictPathogenicity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictPathogenicity(Annotation("StopGained", "High"), 0.00001));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PredictPathogenicity(null!));
    }

    [Test]
    public void PredictPathogenicity_HighImpactRare_IsPathogenicAndActionable()
    {
        // Mirrors VariantAnnotatorTests.PredictPathogenicity_HighImpact_LikelyPathogenic.
        var result = AnnotationTools.PredictPathogenicity(Annotation("StopGained", "High"), populationFrequency: 0.00001);

        Assert.Multiple(() =>
        {
            Assert.That(result.Prediction.Classification, Is.EqualTo("Pathogenic").Or.EqualTo("LikelyPathogenic"));
            Assert.That(result.Prediction.IsActionable, Is.True);
            Assert.That(result.Prediction.EvidenceCriteria, Has.Some.Contains("PVS1"));
        });
    }

    [Test]
    public void PredictPathogenicity_CommonSynonymous_IsBenignNotActionable()
    {
        // Mirrors PredictPathogenicity_CommonVariant_Benign: synonymous + 10% AF.
        var result = AnnotationTools.PredictPathogenicity(Annotation("SynonymousVariant", "Low"), populationFrequency: 0.10);

        Assert.Multiple(() =>
        {
            Assert.That(result.Prediction.Classification, Is.EqualTo("Benign").Or.EqualTo("LikelyBenign"));
            Assert.That(result.Prediction.IsActionable, Is.False);
        });
    }
}
