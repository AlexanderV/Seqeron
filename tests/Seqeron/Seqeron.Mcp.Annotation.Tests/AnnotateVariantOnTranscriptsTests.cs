using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnnotateVariantOnTranscriptsTests
{
    // Exons/coding-exons mirror the authoritative VariantAnnotatorTests fixtures
    // (Seqeron.Genomics.Tests/Unit/Annotation/VariantAnnotatorTests.cs).
    private static readonly IReadOnlyList<GenomicIntervalDto> Exons = new List<GenomicIntervalDto>
    {
        new(1000, 1100), new(1400, 1500), new(1800, 2000)
    };

    private static readonly IReadOnlyList<GenomicIntervalDto> CodingExons = new List<GenomicIntervalDto>
    {
        new(1050, 1100), new(1400, 1500), new(1800, 1900)
    };

    private static TranscriptDto Transcript() =>
        new("T1", "G1", "Gene1", "chr1", 1000, 2000, '+', Exons, CodingExons, 1050, 1900);

    private static AnnotatorVariantDto Snv(int position) =>
        new("chr1", position, "A", "G", "SNV", null, null);

    [Test]
    public void AnnotateVariantOnTranscripts_Schema_ValidatesCorrectly()
    {
        var ts = new List<TranscriptDto> { Transcript() };
        Assert.DoesNotThrow(() => AnnotationTools.AnnotateVariantOnTranscripts(Snv(1450), ts));

        Assert.Throws<ArgumentNullException>(
            () => AnnotationTools.AnnotateVariantOnTranscripts(null!, ts));
        // empty reference
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.AnnotateVariantOnTranscripts(new("chr1", 1450, "", "G", "SNV", null, null), ts));
        // null / empty transcripts
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.AnnotateVariantOnTranscripts(Snv(1450), null!));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.AnnotateVariantOnTranscripts(Snv(1450), new List<TranscriptDto>()));
        // unknown variant type surfaces as ArgumentException from FromDto
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.AnnotateVariantOnTranscripts(new("chr1", 1450, "A", "G", "Bogus", null, null), ts));
    }

    [Test]
    public void AnnotateVariantOnTranscripts_MissenseVariant_ReturnsMissenseModerate()
    {
        // chr1:1450 A>G inside coding exon [1400,1500] -> MissenseVariant / Moderate.
        var ts = new List<TranscriptDto> { Transcript() };

        var result = AnnotationTools.AnnotateVariantOnTranscripts(Snv(1450), ts);

        Assert.That(result.Annotations, Has.Count.EqualTo(1));
        var a = result.Annotations[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.TranscriptId, Is.EqualTo("T1"));
            Assert.That(a.GeneId, Is.EqualTo("G1"));
            Assert.That(a.GeneName, Is.EqualTo("Gene1"));
            Assert.That(a.Consequence, Is.EqualTo("MissenseVariant"));
            Assert.That(a.Impact, Is.EqualTo("Moderate"));
        });
    }

    [Test]
    public void AnnotateVariantOnTranscripts_IntronVariant_ReturnsIntronModifier()
    {
        // chr1:1250 A>G between exons [1000,1100] and [1400,1500] -> IntronVariant / Modifier.
        var ts = new List<TranscriptDto> { Transcript() };

        var a = AnnotationTools.AnnotateVariantOnTranscripts(Snv(1250), ts).Annotations[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.Consequence, Is.EqualTo("IntronVariant"));
            Assert.That(a.Impact, Is.EqualTo("Modifier"));
        });
    }

    [Test]
    public void AnnotateVariantOnTranscripts_NoRelevantTranscript_ReturnsIntergenic()
    {
        // chr1:1000000 far from transcript [1000,2000] -> single IntergenicVariant / Modifier.
        var ts = new List<TranscriptDto> { Transcript() };

        var result = AnnotationTools.AnnotateVariantOnTranscripts(Snv(1000000), ts);

        Assert.That(result.Annotations, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Annotations[0].Consequence, Is.EqualTo("IntergenicVariant"));
            Assert.That(result.Annotations[0].Impact, Is.EqualTo("Modifier"));
        });
    }
}
