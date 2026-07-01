using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>predict_g_bands</c>. ChromosomeAnalyzer.PredictGBands scores each band by GC:
/// GC &lt; darkThreshold -> gpos100, &lt; lightThreshold -> gpos50, else gneg. Bands are named
/// "{chr}{arm}{n}", switching from p to q at the sequence midpoint.
/// </summary>
[TestFixture]
public class PredictGBandsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void PredictGBands_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.PredictGBands("chr1", "ACGTACGT", 4));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.PredictGBands("chr1", ""));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.PredictGBands("", "ACGT"));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.PredictGBands("chr1", "ACGT", bandSize: 0));
    }

    [Test]
    public void PredictGBands_Binding_InvokesSuccessfully()
    {
        // 10 bp AT (GC 0) + 10 bp GC (GC 1), bandSize 10 -> two bands.
        var result = ChromosomeTools.PredictGBands("chr1", Rep("AT", 5) + Rep("GC", 5), 10);

        Assert.That(result.Items, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            // AT-rich band -> dark (gpos100), p arm.
            Assert.That(result.Items[0].Name, Is.EqualTo("chr1p1"));
            Assert.That(result.Items[0].Stain, Is.EqualTo("gpos100"));
            Assert.That(result.Items[0].GcContent, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(result.Items[0].Start, Is.EqualTo(0));
            Assert.That(result.Items[0].End, Is.EqualTo(9));
            // GC-rich band -> light (gneg), q arm after midpoint switch.
            Assert.That(result.Items[1].Name, Is.EqualTo("chr1q1"));
            Assert.That(result.Items[1].Stain, Is.EqualTo("gneg"));
            Assert.That(result.Items[1].GcContent, Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
