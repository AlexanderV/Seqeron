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

    [Test]
    public void ParseGff3_MultiValuedAttributeSeqidSource_DelegatedFields()
    {
        // Mirrors GenomeAnnotator_GFF3_Tests.ParseGff3_SpecCanonicalGene_MultipleParentageIsValueList (B11 F8).
        string gff = "##gff-version 3\r\nctg123\tsrc%09x\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002;Note=a%2Cb\r\n##FASTA\r\n>ctg123\r\nACGT";
        var f = AnnotationTools.ParseGff3(gff).Features.Single();

        Assert.Multiple(() =>
        {
            Assert.That(f.Seqid, Is.EqualTo("ctg123"));
            Assert.That(f.Source, Is.EqualTo("src\tx"));
            Assert.That(f.AttributeValues!["Parent"], Is.EqualTo(new[] { "mRNA00001", "mRNA00002" }));
            Assert.That(f.AttributeValues!["Note"], Is.EqualTo(new[] { "a,b" }));
            Assert.That(f.Attributes["Note"], Is.EqualTo("a,b"));
        });
    }

    [Test]
    public void ParseGff3_SpecInvalidLine_Skipped()
    {
        // Mirrors GenomeAnnotator_GFF3_Tests.ParseGff3_SpecInvalidLine_Skipped (B11 F9): strand 'Z'.
        var result = AnnotationTools.ParseGff3("chr1\tsrc\tgene\t1\t100\t.\tZ\t.\tID=bad\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=ok");
        Assert.That(result.Features.Select(x => x.FeatureId), Is.EqualTo(new[] { "ok" }));
    }
}
