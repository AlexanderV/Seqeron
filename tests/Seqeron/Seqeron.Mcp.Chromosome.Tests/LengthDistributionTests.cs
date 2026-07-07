using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>length_distribution</c>. GenomeAssemblyAnalyzer.CalculateLengthDistribution buckets
/// each length into the first "&lt;bin" it falls under (default bins 100,500,1000,...,1000000),
/// or the "&gt;=maxBin" overflow bucket.
/// </summary>
[TestFixture]
public class LengthDistributionTests
{
    [Test]
    public void LengthDistribution_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.LengthDistribution(new List<int> { 100, 500 }));
        Assert.DoesNotThrow(() => ChromosomeTools.LengthDistribution(new List<int>()));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.LengthDistribution(null!));
    }

    [Test]
    public void LengthDistribution_Binding_InvokesSuccessfully()
    {
        // 50 -> "<100", 150 -> "<500", 600 -> "<1000".
        var dist = ChromosomeTools.LengthDistribution(new List<int> { 50, 150, 600 }).Distribution;

        Assert.Multiple(() =>
        {
            Assert.That(dist["<100"], Is.EqualTo(1));
            Assert.That(dist["<500"], Is.EqualTo(1));
            Assert.That(dist["<1000"], Is.EqualTo(1));
            Assert.That(dist["<5000"], Is.EqualTo(0));
            Assert.That(dist[">=1000000"], Is.EqualTo(0));
        });
    }

    [Test]
    public void LengthDistribution_OverflowBucket()
    {
        // 2,000,000 exceeds the largest default bin (1,000,000) -> overflow bucket.
        var dist = ChromosomeTools.LengthDistribution(new List<int> { 2_000_000 }).Distribution;
        Assert.That(dist[">=1000000"], Is.EqualTo(1));
    }
}
