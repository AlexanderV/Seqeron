using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;
using SourceGa = Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>gap_distribution</c>. GenomeAssemblyAnalyzer.AnalyzeGapDistribution reports count,
/// mean/median/max gap length and per-type counts (median = lengths[count/2] after ascending sort).
/// </summary>
[TestFixture]
public class GapDistributionTests
{
    [Test]
    public void GapDistribution_Schema_ValidatesCorrectly()
    {
        var gaps = new List<SourceGa.GapInfo> { new("s", 0, 4, 5, "Short") };
        Assert.DoesNotThrow(() => ChromosomeTools.GapDistribution(gaps));
        Assert.DoesNotThrow(() => ChromosomeTools.GapDistribution(new List<SourceGa.GapInfo>()));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.GapDistribution(null!));
    }

    [Test]
    public void GapDistribution_Binding_InvokesSuccessfully()
    {
        // Lengths 5 (Short) and 100 (Long): mean=52.5, median=lengths[1]=100, max=100.
        var gaps = new List<SourceGa.GapInfo> { new("s", 0, 4, 5, "Short"), new("s", 10, 109, 100, "Long") };
        var gd = ChromosomeTools.GapDistribution(gaps);

        Assert.Multiple(() =>
        {
            Assert.That(gd.Count, Is.EqualTo(2));
            Assert.That(gd.MeanLength, Is.EqualTo(52.5).Within(1e-9));
            Assert.That(gd.MedianLength, Is.EqualTo(100).Within(1e-9));
            Assert.That(gd.MaxLength, Is.EqualTo(100));
            Assert.That(gd.TypeCounts["Short"], Is.EqualTo(1));
            Assert.That(gd.TypeCounts["Long"], Is.EqualTo(1));
        });
    }
}
