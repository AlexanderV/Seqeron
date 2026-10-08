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
        // Grimson (2007) / TargetScan context (MiRnaAnalyzer.AnalyzeTargetContext):
        // 40-nt poly(A), site 20..23: weighted local AU = 1.0; NearStart = (21 < 15) false;
        // NearEnd = 23 > 34 false; d5 = 20, d3 = 16 ⇒ EndProximity = 1 − 16/18;
        // ContextScore = 0.5·1.0 + 0.5·(1 − 16/18).
        var au = AnnotationTools.AnalyzeTargetContext(new string('A', 40), 20, 23, 30);
        Assert.Multiple(() =>
        {
            Assert.That(au.AuContent, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(au.NearStart, Is.False);
            Assert.That(au.NearEnd, Is.False);
            Assert.That(au.ContextScore, Is.EqualTo(0.5 + 0.5 * (1 - 16.0 / 18.0)).Within(1e-9));
        });

        // No A/U ⇒ AuContent = 0; only the positional term remains.
        var gc = AnnotationTools.AnalyzeTargetContext(string.Concat(Enumerable.Repeat("GC", 20)), 20, 23, 30);
        Assert.Multiple(() =>
        {
            Assert.That(gc.AuContent, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(gc.ContextScore, Is.EqualTo(0.5 * (1 - 16.0 / 18.0)).Within(1e-9));
        });

        // Site within the first 15 nt of the 3'UTR (TargetScan MIN_DIST_TO_CDS): nearStart, score 0.
        var start = AnnotationTools.AnalyzeTargetContext("AAAAAAAAAAAAAAAAAAAA", 0, 1, 30);
        Assert.That(start.NearStart, Is.True);
        Assert.That(start.ContextScore, Is.EqualTo(0.0));
    }
}
