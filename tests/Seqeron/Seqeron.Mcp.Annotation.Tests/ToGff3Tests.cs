using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ToGff3Tests
{
    private static GeneAnnotationDto Gene() => new(
        GeneId: "gene1",
        Start: 99,   // 0-based internal
        End: 500,
        Strand: '+',
        Type: "CDS",
        Product: "hypothetical protein",
        Attributes: new Dictionary<string, string> { ["frame"] = "1" });

    [Test]
    public void ToGff3_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ToGff3(new[] { Gene() }, "chr1"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ToGff3(Array.Empty<GeneAnnotationDto>(), "chr1"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ToGff3(null!, "chr1"));
    }

    [Test]
    public void ToGff3_EmitsVersionHeaderThenNineColumns()
    {
        // Mirrors GenomeAnnotator_GFF3_Tests.ToGff3_GeneratesValidTabDelimitedOutput.
        var result = AnnotationTools.ToGff3(new[] { Gene() }, "chr1");

        Assert.That(result.Lines, Has.Count.EqualTo(2));
        Assert.That(result.Lines[0], Is.EqualTo("##gff-version 3"));

        var fields = result.Lines[1].Split('\t');
        Assert.That(fields, Has.Length.EqualTo(9));
        Assert.Multiple(() =>
        {
            Assert.That(fields[0], Is.EqualTo("chr1"));      // seqid
            Assert.That(fields[1], Is.EqualTo("."));         // source
            Assert.That(fields[2], Is.EqualTo("CDS"));       // type
            Assert.That(fields[3], Is.EqualTo("100"));       // start = 99 + 1 (1-based)
            Assert.That(fields[4], Is.EqualTo("500"));       // end
            Assert.That(fields[5], Is.EqualTo("."));         // score
            Assert.That(fields[6], Is.EqualTo("+"));         // strand
            Assert.That(fields[7], Is.EqualTo("0"));         // phase (CDS)
            Assert.That(fields[8], Does.Contain("ID=gene1"));
            Assert.That(fields[8], Does.Contain("product=hypothetical protein"));
        });
    }

    [Test]
    public void ToGff3_RoundTripsThroughParseGff3()
    {
        var lines = AnnotationTools.ToGff3(new[] { Gene() }, "chr1").Lines;
        var parsed = AnnotationTools.ParseGff3(string.Join("\n", lines));

        Assert.That(parsed.Features, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(parsed.Features[0].FeatureId, Is.EqualTo("gene1"));
            Assert.That(parsed.Features[0].Type, Is.EqualTo("CDS"));
            Assert.That(parsed.Features[0].Start, Is.EqualTo(100)); // 1-based as serialized
            Assert.That(parsed.Features[0].End, Is.EqualTo(500));
            Assert.That(parsed.Features[0].Strand, Is.EqualTo('+'));
        });
    }

    [Test]
    public void ToGff3_SourceScoreAttributeValuesAndSeqid_Delegated()
    {
        // Mirrors GenomeAnnotator_GFF3_Tests.ToGff3_AttributeValues_CommaSeparatedAndRoundTrip /
        // ToGff3_SeqidOutsideAllowedSet_PercentEncoded (B11 F10).
        var exon = new GeneAnnotationDto("e1", 1049, 1500, '+', "exon", "", new Dictionary<string, string>(),
            Source: "src", Score: 0.5,
            AttributeValues: new Dictionary<string, List<string>> { ["Parent"] = new() { "mRNA00001", "mRNA00002" } });

        var line = AnnotationTools.ToGff3(new[] { exon }, "chr 1").Lines[1];

        Assert.That(line, Is.EqualTo("chr%201\tsrc\texon\t1050\t1500\t0.5\t+\t.\tID=e1;Parent=mRNA00001,mRNA00002"));
    }

    [Test]
    public void ToGff3_InvalidStrand_Throws()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ToGff3(new[] { Gene() with { Strand = 'x' } }, "chr1"));
    }
}
