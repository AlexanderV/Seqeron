using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class IsWithinCodingRegionTests
{
    // AUG start at index 0; downstream positions tested for reading frame.
    private const string Seq = "AUGAAAAAAAAA";

    [Test]
    public void IsWithinCodingRegion_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.IsWithinCodingRegion(Seq, 6));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsWithinCodingRegion("", 0));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsWithinCodingRegion(null!, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.IsWithinCodingRegion(Seq, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.IsWithinCodingRegion(Seq, 100));
    }

    [Test]
    public void IsWithinCodingRegion_Binding_InvokesSuccessfully()
    {
        // SpliceSitePredictor.IsWithinCodingRegion: finds an upstream AUG, returns (position - i) % 3 == frame.
        // AUG at i=0; position 6 in frame 0: (6-0)%3 == 0 -> true.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.IsWithinCodingRegion(Seq, 6, frame: 0).IsCoding, Is.True);
            // position 5 in frame 0: (5-0)%3 == 2 != 0 -> false.
            Assert.That(AnnotationTools.IsWithinCodingRegion(Seq, 5, frame: 0).IsCoding, Is.False);
            // position 5 in frame 2: (5-0)%3 == 2 == 2 -> true.
            Assert.That(AnnotationTools.IsWithinCodingRegion(Seq, 5, frame: 2).IsCoding, Is.True);
        });
    }

    [Test]
    public void IsWithinCodingRegion_NoUpstreamStart_ReturnsFalse()
    {
        // No AUG upstream of the position -> not coding.
        Assert.That(AnnotationTools.IsWithinCodingRegion("CCCCCCCCC", 6, frame: 0).IsCoding, Is.False);
    }
}
