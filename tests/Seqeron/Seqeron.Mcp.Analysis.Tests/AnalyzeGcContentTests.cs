using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>analyze_gc_content</c> MCP tool.
/// Expected values are taken from the algorithm's own unit tests
/// (Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_AnalyzeGcContent_Tests.cs,
/// spec SEQ-GC-ANALYSIS-001), NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class AnalyzeGcContentTests
{
    [Test]
    public void AnalyzeGcContent_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.AnalyzeGcContent("GGGCCAT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeGcContent(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeGcContent(null!));
        // Non-DNA input must be rejected by the DNA guard.
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeGcContent("ACGTX"));
    }

    [Test]
    public void AnalyzeGcContent_Binding_InvokesSuccessfully()
    {
        // "GGGCCAT": G=3,C=2,A=1,T=1,n=7.
        // GC% = (3+2)/7*100 = 71.42857142857143; GC skew = (3-2)/5 = 0.2; AT skew = (1-1)/2 = 0.0.
        var r = AnalysisTools.AnalyzeGcContent("GGGCCAT");
        Assert.Multiple(() =>
        {
            Assert.That(r.OverallGcContent, Is.EqualTo(5.0 / 7.0 * 100.0).Within(1e-10));
            Assert.That(r.OverallGcSkew, Is.EqualTo(0.2).Within(1e-10));
            Assert.That(r.OverallAtSkew, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(r.SequenceLength, Is.EqualTo(7));
            // len < default window (1000): profiles empty, variances 0.
            Assert.That(r.WindowedGcSkew, Is.Empty);
            Assert.That(r.WindowedGcContent, Is.Empty);
            Assert.That(r.GcSkewVariance, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(r.GcContentVariance, Is.EqualTo(0.0).Within(1e-10));
        });
    }

    [Test]
    public void AnalyzeGcContent_WindowedPopulationVariance()
    {
        // "GGCC" window 2 step 2 -> windows GG(+1,100%), CC(-1,100%).
        // GcSkewVariance = ((1-0)^2+(-1-0)^2)/2 = 1.0 (population, /N not /N-1).
        // GcContentVariance = variance of {100,100} = 0.0.
        var r = AnalysisTools.AnalyzeGcContent("GGCC", windowSize: 2, stepSize: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.WindowedGcSkew, Has.Length.EqualTo(2));
            Assert.That(r.GcSkewVariance, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(r.GcContentVariance, Is.EqualTo(0.0).Within(1e-10));
            Assert.That(r.WindowedGcContent[0].WindowStart, Is.EqualTo(0));
            Assert.That(r.WindowedGcContent[0].WindowEnd, Is.EqualTo(1));
            Assert.That(r.WindowedGcContent[0].Position, Is.EqualTo(1));
        });
    }

    // fraction=true forwards to Core AnalyzeGcContent(fraction): GC content in [0,1] (Biopython 1.88
    // gc_fraction("GGGCCAT") = 0.7142857142857143); skew values are unaffected.
    [Test]
    public void AnalyzeGcContent_Fraction_ReportsFractions()
    {
        var r = AnalysisTools.AnalyzeGcContent("GGGCCAT", windowSize: 4, stepSize: 3, fraction: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.OverallGcContent, Is.EqualTo(5.0 / 7.0).Within(1e-12));
            Assert.That(r.OverallGcSkew, Is.EqualTo(0.2).Within(1e-12));
            // windows GGGC (gc_fraction 1.0), CCAT (0.5)
            Assert.That(r.WindowedGcContent.Select(p => p.GcContent), Is.EqualTo(new[] { 1.0, 0.5 }).Within(1e-12));
            Assert.That(r.GcContentVariance, Is.EqualTo(0.0625).Within(1e-12));
            // Biopython GC_skew("GGGC")=0.5, GC_skew("CCAT")=-1.0 → numpy.var = 0.5625 (unaffected by fraction).
            Assert.That(r.GcSkewVariance, Is.EqualTo(0.5625).Within(1e-12));
        });
    }
}
