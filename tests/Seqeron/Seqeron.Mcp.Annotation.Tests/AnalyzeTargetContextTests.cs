using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnalyzeTargetContextTests
{
    [Test]
    public void AnalyzeTargetContext_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.AnalyzeTargetContext("AAAAAAAAAAAAAAAAAAAA", 8, 11));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AnalyzeTargetContext("", 0, 1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AnalyzeTargetContext(null!, 0, 1));
        // start < 0, end past sequence, start > end all out of range
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AnalyzeTargetContext("AAAA", -1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AnalyzeTargetContext("AAAA", 0, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AnalyzeTargetContext("AAAA", 3, 1));
    }

    [Test]
    public void AnalyzeTargetContext_Binding_InvokesSuccessfully()
    {
        // 20-nt all-A window. auContent = 20/20 = 1.0.
        // nearStart: 8 < 20*0.15 (=3.0) -> false; nearEnd: 11 > 20*0.85 (=17.0) -> false.
        // contextScore = 1.0*0.5 + 0.3 (mid-transcript bonus) = 0.8.
        var au = AnnotationTools.AnalyzeTargetContext("AAAAAAAAAAAAAAAAAAAA", 8, 11, 30);
        Assert.Multiple(() =>
        {
            Assert.That(au.AuContent, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(au.NearStart, Is.False);
            Assert.That(au.NearEnd, Is.False);
            Assert.That(au.ContextScore, Is.EqualTo(0.8).Within(1e-9));
        });

        // 20-nt all-GC window: no A/U -> auContent = 0; mid-transcript bonus only -> 0.3.
        var gc = AnnotationTools.AnalyzeTargetContext("GCGCGCGCGCGCGCGCGCGC", 8, 11, 30);
        Assert.Multiple(() =>
        {
            Assert.That(gc.AuContent, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(gc.ContextScore, Is.EqualTo(0.3).Within(1e-9));
        });

        // Site at the very start: nearStart true (0 < 3.0).
        var start = AnnotationTools.AnalyzeTargetContext("AAAAAAAAAAAAAAAAAAAA", 0, 1, 30);
        Assert.That(start.NearStart, Is.True);
    }
}
