using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class Log2TransformTests
{
    [Test]
    public void Log2Transform_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.Log2Transform(new List<double> { 1.0 }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.Log2Transform(null!));
    }

    [Test]
    public void Log2Transform_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.Log2Transform: log2(x + pseudocount). Default pseudocount 1.
        // [0,1,3,7] -> log2([1,2,4,8]) = [0,1,2,3].
        var result = AnnotationTools.Log2Transform(new List<double> { 0.0, 1.0, 3.0, 7.0 });
        Assert.That(result.Transformed, Is.EqualTo(new[] { 0.0, 1.0, 2.0, 3.0 }).Within(1e-9));
    }

    [Test]
    public void Log2Transform_CustomPseudocount()
    {
        // pseudocount 3: [1,5] -> log2([4,8]) = [2,3].
        var result = AnnotationTools.Log2Transform(new List<double> { 1.0, 5.0 }, pseudocount: 3.0);
        Assert.That(result.Transformed, Is.EqualTo(new[] { 2.0, 3.0 }).Within(1e-9));
    }
}
