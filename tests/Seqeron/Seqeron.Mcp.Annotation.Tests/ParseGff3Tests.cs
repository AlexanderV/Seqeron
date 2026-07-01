using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ParseGff3Tests
{
    [Test]
    public void ParseGff3_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ParseGff3("chr1\tENSEMBL\tCDS\t1000\t2000\t95.5\t-\t2\tID=cds1"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ParseGff3(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ParseGff3(null!));
    }

    [Test]
    public void ParseGff3_AllColumns_ExtractsCorrectly()
    {
        // Mirrors GenomeAnnotator_GFF3_Tests.ParseGff3_AllColumns_ExtractsCorrectly.
        var result = AnnotationTools.ParseGff3("chr1\tENSEMBL\tCDS\t1000\t2000\t95.5\t-\t2\tID=cds1;Name=TestCDS");

        Assert.That(result.Features, Has.Count.EqualTo(1));
        var f = result.Features[0];
        Assert.Multiple(() =>
        {
            Assert.That(f.FeatureId, Is.EqualTo("cds1"));
            Assert.That(f.Type, Is.EqualTo("CDS"));
            Assert.That(f.Start, Is.EqualTo(1000));
            Assert.That(f.End, Is.EqualTo(2000));
            Assert.That(f.Score, Is.EqualTo(95.5).Within(0.01));
            Assert.That(f.Strand, Is.EqualTo('-'));
            Assert.That(f.Phase, Is.EqualTo(2));
            Assert.That(f.Attributes["Name"], Is.EqualTo("TestCDS"));
        });
    }

    [Test]
    public void ParseGff3_UndefinedScoreAndPhase_ReturnNull()
    {
        // Mirrors ParseGff3_UndefinedScore_ReturnsNull / ParseGff3_UndefinedPhase_ReturnsNull.
        var result = AnnotationTools.ParseGff3("chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1");

        Assert.That(result.Features, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Features[0].Score, Is.Null);
            Assert.That(result.Features[0].Phase, Is.Null);
        });
    }

    [Test]
    public void ParseGff3_SkipsCommentsAndDirectives()
    {
        // Mirrors ParseGff3_SkipsComments / ParseGff3_SkipsDirectives.
        string gff = "##gff-version 3\n# a comment\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1";
        var result = AnnotationTools.ParseGff3(gff);

        Assert.That(result.Features, Has.Count.EqualTo(1));
        Assert.That(result.Features[0].FeatureId, Is.EqualTo("gene1"));
    }
}
