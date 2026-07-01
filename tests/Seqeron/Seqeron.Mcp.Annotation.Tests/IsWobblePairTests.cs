using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class IsWobblePairTests
{
    [Test]
    public void IsWobblePair_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.IsWobblePair("G", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsWobblePair("", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsWobblePair("GG", "U"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsWobblePair("G", null!));
    }

    [Test]
    public void IsWobblePair_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.IsWobblePair: only G-U / U-G (DNA T treated as U).
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.IsWobblePair("G", "U").IsWobble, Is.True);
            Assert.That(AnnotationTools.IsWobblePair("U", "G").IsWobble, Is.True);
            Assert.That(AnnotationTools.IsWobblePair("G", "T").IsWobble, Is.True);   // T -> U
            Assert.That(AnnotationTools.IsWobblePair("A", "U").IsWobble, Is.False);  // Watson-Crick, not wobble
            Assert.That(AnnotationTools.IsWobblePair("G", "C").IsWobble, Is.False);  // Watson-Crick, not wobble
        });
    }
}
