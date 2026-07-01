using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CodonUsageStatisticsTests
{
    [Test]
    public void CodonUsageStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.codon_usage_statistics("ATGTTTAAA"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_usage_statistics(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_usage_statistics(null!));
    }

    [Test]
    public void CodonUsageStatistics_Binding_InvokesSuccessfully()
    {
        // Codons ATG(M), TTT(F), AAA(K):
        //   TotalCodons = 3.
        //   GC1: A,T,A -> 0/3 = 0%.
        //   GC2: T,T,A -> 0/3 = 0%.
        //   GC3: G,T,A -> 1/3 = 33.333%.
        //   OverallGc = (0+0+33.333)/3 = 11.111%.
        //   GC3s: synonymous-at-third codons are TTT(Phe) and AAA(Lys); pos3 = T,A -> 0/2 = 0%.
        //   RSCU[TTT] = observed 1 / expected (1/2) = 2.0.
        var stats = MolToolsTools.codon_usage_statistics("ATGTTTAAA");

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalCodons, Is.EqualTo(3));
            Assert.That(stats.CodonCounts["ATG"], Is.EqualTo(1));
            Assert.That(stats.CodonCounts["TTT"], Is.EqualTo(1));
            Assert.That(stats.CodonCounts["AAA"], Is.EqualTo(1));
            Assert.That(stats.Gc1, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(stats.Gc2, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(stats.Gc3, Is.EqualTo(100.0 / 3.0).Within(1e-9));
            Assert.That(stats.Gc3s, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(stats.OverallGc, Is.EqualTo((100.0 / 3.0) / 3.0).Within(1e-9));
            Assert.That(stats.Rscu["TTT"], Is.EqualTo(2.0).Within(1e-9));
        });
    }
}
