using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>detect_aneuploidy</c>. GenomeAssemblyAnalyzer... ChromosomeAnalyzer.DetectAneuploidy
/// bins depth by position; per bin: logRatio = log2(mean/median), copyNumber = round(2^logRatio * 2).
/// A bin at exactly the median depth gives logRatio 0, copyNumber 2, confidence 1.
/// </summary>
[TestFixture]
public class DetectAneuploidyTests
{
    [Test]
    public void DetectAneuploidy_Schema_ValidatesCorrectly()
    {
        var depth = new List<DepthSample> { new("chr1", 0, 30.0) };
        Assert.DoesNotThrow(() => ChromosomeTools.DetectAneuploidy(depth, 30.0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectAneuploidy(null!, 30.0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectAneuploidy(depth, 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectAneuploidy(depth, 30.0, binSize: 0));
    }

    [Test]
    public void DetectAneuploidy_Binding_InvokesSuccessfully()
    {
        // All positions in bin 0, depth == median 30 -> diploid.
        var depth = new List<DepthSample>();
        for (int p = 0; p < 5; p++) depth.Add(new DepthSample("chr1", p, 30.0));

        var result = ChromosomeTools.DetectAneuploidy(depth, 30.0, 1000000);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var s = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(s.Chromosome, Is.EqualTo("chr1"));
            Assert.That(s.Start, Is.EqualTo(0));
            Assert.That(s.End, Is.EqualTo(999999));
            Assert.That(s.CopyNumber, Is.EqualTo(2));
            Assert.That(s.LogRatio, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(s.Confidence, Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
