using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class BedMergeTests
{
    [Test]
    public void BedMerge_Schema_ValidatesCorrectly()
    {
        var validBed = "chr1\t100\t200";
        Assert.DoesNotThrow(() => ParsersTools.BedMerge(validBed));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedMerge(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedMerge(null!));
    }

    [Test]
    public void BedMerge_Binding_MergesOverlapping()
    {
        var bed = "chr1\t100\t200\nchr1\t150\t250\nchr1\t400\t500";
        var result = ParsersTools.BedMerge(bed);

        Assert.That(result.OriginalCount, Is.EqualTo(3));
        Assert.That(result.MergedCount, Is.EqualTo(2)); // First two merge, third separate
        Assert.That(result.Records[0].ChromStart, Is.EqualTo(100));
        Assert.That(result.Records[0].ChromEnd, Is.EqualTo(250)); // Merged
        Assert.That(result.Records[1].ChromStart, Is.EqualTo(400));
    }

    [Test]
    public void BedMerge_Binding_KeepsSeparateChromosomes()
    {
        var bed = "chr1\t100\t200\nchr2\t100\t200";
        var result = ParsersTools.BedMerge(bed);

        Assert.That(result.MergedCount, Is.EqualTo(2)); // Different chroms, no merge
    }

    [Test]
    public void BedMerge_Binding_MergesAdjacent()
    {
        var bed = "chr1\t100\t200\nchr1\t200\t300"; // Adjacent (end == start)
        var result = ParsersTools.BedMerge(bed);

        Assert.That(result.MergedCount, Is.EqualTo(1));
        Assert.That(result.Records[0].ChromStart, Is.EqualTo(100));
        Assert.That(result.Records[0].ChromEnd, Is.EqualTo(300));
    }
}
