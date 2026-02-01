using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class BedFilterTests
{
    private const string TestBed = "chr1\t100\t200\tgene1\t500\t+\nchr1\t300\t400\tgene2\t300\t-\nchr2\t100\t500\tgene3\t800\t+";

    [Test]
    public void BedFilter_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.BedFilter(TestBed));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedFilter(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedFilter(null!));
    }

    [Test]
    public void BedFilter_Binding_FiltersByChrom()
    {
        var result = ParsersTools.BedFilter(TestBed, chrom: "chr1");

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.TotalCount, Is.EqualTo(3));
        Assert.That(result.Records.All(r => r.Chrom == "chr1"), Is.True);
    }

    [Test]
    public void BedFilter_Binding_FiltersByRegion()
    {
        var result = ParsersTools.BedFilter(TestBed, chrom: "chr1", regionStart: 150, regionEnd: 350);

        Assert.That(result.PassedCount, Is.EqualTo(2)); // gene1 (100-200) and gene2 (300-400) overlap with 150-350
        Assert.That(result.Records[0].Name, Is.EqualTo("gene1"));
        Assert.That(result.Records[1].Name, Is.EqualTo("gene2"));
    }

    [Test]
    public void BedFilter_Binding_FiltersByStrand()
    {
        var result = ParsersTools.BedFilter(TestBed, strand: "+");

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.Records.All(r => r.Strand == "+"), Is.True);
    }

    [Test]
    public void BedFilter_Binding_FiltersByLength()
    {
        var result = ParsersTools.BedFilter(TestBed, minLength: 200);

        Assert.That(result.PassedCount, Is.EqualTo(1)); // Only gene3 (400bp)
        Assert.That(result.Records[0].Name, Is.EqualTo("gene3"));
    }

    [Test]
    public void BedFilter_Binding_FiltersByScore()
    {
        var result = ParsersTools.BedFilter(TestBed, minScore: 400);

        Assert.That(result.PassedCount, Is.EqualTo(2)); // gene1 (500) and gene3 (800)
    }

    [Test]
    public void BedFilter_Binding_CombinesFilters()
    {
        var result = ParsersTools.BedFilter(TestBed, chrom: "chr1", strand: "+");

        Assert.That(result.PassedCount, Is.EqualTo(1)); // Only gene1
        Assert.That(result.Records[0].Name, Is.EqualTo("gene1"));
    }

    [Test]
    public void BedFilter_Strand_InvalidStrandThrows()
    {
        Assert.Throws<ArgumentException>(() => ParsersTools.BedFilter(TestBed, strand: "x"));
    }
}
