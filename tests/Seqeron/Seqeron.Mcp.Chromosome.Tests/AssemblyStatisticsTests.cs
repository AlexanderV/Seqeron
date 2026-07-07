using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>assembly_statistics</c>. Expected values from GenomeAssemblyAnalyzer.CalculateStatistics
/// (N50/L50 per Miller, Koren &amp; Sutton 2010).
/// </summary>
[TestFixture]
public class AssemblyStatisticsTests
{
    [Test]
    public void AssemblyStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.AssemblyStatistics(new List<NamedSequence> { new("s", "ACGT") }));
        Assert.DoesNotThrow(() => ChromosomeTools.AssemblyStatistics(new List<NamedSequence>()));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssemblyStatistics(null!));
    }

    [Test]
    public void AssemblyStatistics_Binding_InvokesSuccessfully()
    {
        // s1 = "ACGTACGTGG" (10 bp, GC=6), s2 = "ACGTAC" (6 bp, GC=3) -> 16 bp, 9 GC / 16 = 0.5625.
        var seqs = new List<NamedSequence> { new("s1", "ACGTACGTGG"), new("s2", "ACGTAC") };
        var stats = ChromosomeTools.AssemblyStatistics(seqs);

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalSequences, Is.EqualTo(2));
            Assert.That(stats.TotalLength, Is.EqualTo(16));
            Assert.That(stats.LargestContig, Is.EqualTo(10));
            Assert.That(stats.SmallestContig, Is.EqualTo(6));
            Assert.That(stats.MeanLength, Is.EqualTo(8.0).Within(1e-9));
            Assert.That(stats.GcContent, Is.EqualTo(0.5625).Within(1e-9));
            // Sorted desc [10, 6]: cumulative >= 50% of 16 (=8) at first -> N50=10, L50=1.
            Assert.That(stats.N50, Is.EqualTo(10));
            Assert.That(stats.L50, Is.EqualTo(1));
        });
    }
}
