using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class BedParseTests
{
    [Test]
    public void BedParse_Schema_ValidatesCorrectly()
    {
        var validBed = "chr1\t100\t200";
        Assert.DoesNotThrow(() => ParsersTools.BedParse(validBed));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedParse(null!));
    }

    [Test]
    public void BedParse_Binding_ParsesBed3()
    {
        var bed = "chr1\t100\t200\nchr2\t300\t500";
        var result = ParsersTools.BedParse(bed);

        Assert.That(result.Count, Is.EqualTo(2));
        Assert.That(result.Records[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(result.Records[0].ChromStart, Is.EqualTo(100));
        Assert.That(result.Records[0].ChromEnd, Is.EqualTo(200));
        Assert.That(result.Records[0].Length, Is.EqualTo(100));
        Assert.That(result.Records[1].Chrom, Is.EqualTo("chr2"));
        Assert.That(result.Records[1].Length, Is.EqualTo(200));
    }

    [Test]
    public void BedParse_Binding_ParsesBed6()
    {
        var bed = "chr1\t100\t200\tgene1\t500\t+";
        var result = ParsersTools.BedParse(bed, "bed6");

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Records[0].Name, Is.EqualTo("gene1"));
        Assert.That(result.Records[0].Score, Is.EqualTo(500));
        Assert.That(result.Records[0].Strand, Is.EqualTo("+"));
    }

    [Test]
    public void BedParse_Format_InvalidFormatThrows()
    {
        var bed = "chr1\t100\t200";
        Assert.Throws<ArgumentException>(() => ParsersTools.BedParse(bed, "invalid"));
    }

    [Test]
    public void BedParse_Binding_SkipsCommentLines()
    {
        var bed = "#comment\nchr1\t100\t200\ntrack name=test\nchr2\t300\t400";
        var result = ParsersTools.BedParse(bed);

        Assert.That(result.Count, Is.EqualTo(2));
        Assert.That(result.Records[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(result.Records[1].Chrom, Is.EqualTo("chr2"));
    }
}
