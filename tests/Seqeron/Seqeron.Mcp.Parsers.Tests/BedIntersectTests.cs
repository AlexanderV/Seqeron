using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class BedIntersectTests
{
    [Test]
    public void BedIntersect_Schema_ValidatesCorrectly()
    {
        var bedA = "chr1\t100\t200";
        var bedB = "chr1\t150\t250";
        Assert.DoesNotThrow(() => ParsersTools.BedIntersect(bedA, bedB));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedIntersect("", bedB));
        Assert.Throws<ArgumentException>(() => ParsersTools.BedIntersect(bedA, ""));
    }

    [Test]
    public void BedIntersect_Binding_FindsOverlap()
    {
        var bedA = "chr1\t100\t200";
        var bedB = "chr1\t150\t250";
        var result = ParsersTools.BedIntersect(bedA, bedB);

        Assert.That(result.IntersectionCount, Is.EqualTo(1));
        Assert.That(result.Records[0].ChromStart, Is.EqualTo(150)); // Intersection start
        Assert.That(result.Records[0].ChromEnd, Is.EqualTo(200));   // Intersection end
    }

    [Test]
    public void BedIntersect_Binding_NoOverlap()
    {
        var bedA = "chr1\t100\t200";
        var bedB = "chr1\t300\t400";
        var result = ParsersTools.BedIntersect(bedA, bedB);

        Assert.That(result.IntersectionCount, Is.EqualTo(0));
    }

    [Test]
    public void BedIntersect_Binding_DifferentChromosomes()
    {
        var bedA = "chr1\t100\t200";
        var bedB = "chr2\t100\t200";
        var result = ParsersTools.BedIntersect(bedA, bedB);

        Assert.That(result.IntersectionCount, Is.EqualTo(0));
    }

    [Test]
    public void BedIntersect_Binding_MultipleOverlaps()
    {
        var bedA = "chr1\t100\t300";
        var bedB = "chr1\t150\t200\nchr1\t250\t350";
        var result = ParsersTools.BedIntersect(bedA, bedB);

        Assert.That(result.IntersectionCount, Is.EqualTo(2));
    }
}
