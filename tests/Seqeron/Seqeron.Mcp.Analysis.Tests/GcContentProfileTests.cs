using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>gc_content_profile</c> MCP tool.
/// Expected values computed from GC = (G+C)/(A+T+U+G+C)*100 per window, NOT the
/// wrapper output.
/// </summary>
[TestFixture]
public class GcContentProfileTests
{
    [Test]
    public void GcContentProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GcContentProfile("GGCC", 2, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GcContentProfile("", 2, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GcContentProfile(null!, 2, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GcContentProfile("GGCC", 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GcContentProfile("GGCC", 2, 0));
    }

    [Test]
    public void GcContentProfile_Binding_InvokesSuccessfully()
    {
        // "GGCC" window 2 step 1 -> GG,GC,CC all 100%.
        var rich = AnalysisTools.GcContentProfile("GGCC", 2, 1).Values;
        Assert.Multiple(() =>
        {
            Assert.That(rich, Has.Length.EqualTo(3));
            Assert.That(rich, Is.All.EqualTo(100.0).Within(1e-9));
        });

        // "ATGC" window 4 -> GC=2/4 = 50%.
        var half = AnalysisTools.GcContentProfile("ATGC", 4, 1).Values;
        Assert.Multiple(() =>
        {
            Assert.That(half, Has.Length.EqualTo(1));
            Assert.That(half[0], Is.EqualTo(50.0).Within(1e-9));
        });

        // "ATAT" window 2 step 2 -> two AT windows, 0%.
        var none = AnalysisTools.GcContentProfile("ATAT", 2, 2).Values;
        Assert.Multiple(() =>
        {
            Assert.That(none, Has.Length.EqualTo(2));
            Assert.That(none, Is.All.EqualTo(0.0).Within(1e-9));
        });
    }
}
