using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CanPairTests
{
    [Test]
    public void CanPair_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CanPair("A", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CanPair("", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CanPair("AU", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CanPair("A", "UU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CanPair("A", null!));
    }

    [Test]
    public void CanPair_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.CanPair: Watson-Crick A-U / G-C plus G-U wobble; DNA T treated as U.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.CanPair("A", "U").CanPair, Is.True);   // Watson-Crick
            Assert.That(AnnotationTools.CanPair("G", "C").CanPair, Is.True);   // Watson-Crick
            Assert.That(AnnotationTools.CanPair("G", "U").CanPair, Is.True);   // wobble
            Assert.That(AnnotationTools.CanPair("A", "T").CanPair, Is.True);   // T normalised to U
            Assert.That(AnnotationTools.CanPair("A", "G").CanPair, Is.False);  // non-pairing
            Assert.That(AnnotationTools.CanPair("A", "C").CanPair, Is.False);  // non-pairing
        });
    }
}
