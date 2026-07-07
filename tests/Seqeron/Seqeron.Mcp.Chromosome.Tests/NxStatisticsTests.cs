using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>nx_statistics</c>. GenomeAssemblyAnalyzer.CalculateNx: Nx is the length of the
/// shortest sequence whose cumulative (descending) coverage first reaches the threshold%; Lx is
/// how many sequences that took (Miller, Koren &amp; Sutton 2010).
/// </summary>
[TestFixture]
public class NxStatisticsTests
{
    [Test]
    public void NxStatistics_Schema_ValidatesCorrectly()
    {
        var lengths = new List<int> { 100, 50 };
        Assert.DoesNotThrow(() => ChromosomeTools.NxStatistics(lengths, 150, 50));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.NxStatistics(null!, 150, 50));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.NxStatistics(lengths, -1, 50));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.NxStatistics(lengths, 150, 101));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.NxStatistics(lengths, 150, -1));
    }

    [Test]
    public void NxStatistics_Binding_InvokesSuccessfully()
    {
        // Sorted desc [100,90,80,70,60,50,40,30,20,10], total = 550.
        // 50% = 275: 100+90+80=270 (<275), +70=340 (>=275) -> N50=70, L50=4, cumulative=340.
        var lengths = new List<int> { 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 };
        var nx = ChromosomeTools.NxStatistics(lengths, 550, 50);

        Assert.Multiple(() =>
        {
            Assert.That(nx.Threshold, Is.EqualTo(50));
            Assert.That(nx.Nx, Is.EqualTo(70));
            Assert.That(nx.Lx, Is.EqualTo(4));
            Assert.That(nx.CumulativeLength, Is.EqualTo(340));
        });
    }
}
