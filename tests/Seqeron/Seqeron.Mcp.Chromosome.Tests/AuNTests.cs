using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>au_n</c>. GenomeAssemblyAnalyzer.CalculateAuN = sum(l^2) / sum(l).
/// </summary>
[TestFixture]
public class AuNTests
{
    [Test]
    public void AuN_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.AuN(new List<int> { 100, 50 }));
        Assert.DoesNotThrow(() => ChromosomeTools.AuN(new List<int>()));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AuN(null!));
    }

    [Test]
    public void AuN_Binding_InvokesSuccessfully()
    {
        // sum(l^2) for 100..10 = 38500; sum(l)=550; auN = 38500/550 = 70.
        var lengths = new List<int> { 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 };
        Assert.That(ChromosomeTools.AuN(lengths).AuN, Is.EqualTo(70.0).Within(1e-9));

        // Single contig: auN equals its length.
        Assert.That(ChromosomeTools.AuN(new List<int> { 123 }).AuN, Is.EqualTo(123.0).Within(1e-9));
        // Empty -> 0.
        Assert.That(ChromosomeTools.AuN(new List<int>()).AuN, Is.EqualTo(0.0));
    }

    [Test]
    public void AuN_GenomeSize_ComputesAuNG()
    {
        // A3-1: auNG = sum(l^2) / G (QUAST au_metric with reference_length): 22000 / 1000 = 22.
        var lengths = new List<int> { 100, 80, 60, 40, 20 };
        Assert.Multiple(() =>
        {
            Assert.That(ChromosomeTools.AuN(lengths, 1000).AuN, Is.EqualTo(22.0).Within(1e-12));
            Assert.That(ChromosomeTools.AuN(lengths, null).AuN, Is.EqualTo(22000.0 / 300.0).Within(1e-12), "omitted genomeSize = auN");
            Assert.Throws<ArgumentException>(() => ChromosomeTools.AuN(lengths, 0));
            Assert.Throws<ArgumentException>(() => ChromosomeTools.AuN(lengths, -5));
        });
    }
}
