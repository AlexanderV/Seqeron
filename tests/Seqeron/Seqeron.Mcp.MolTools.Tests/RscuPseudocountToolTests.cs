using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

// CODON-RSCU-001 (audit A2-9): optional `pseudocount` on the MCP `rscu` tool
// (codonbias 0.5.0 RelativeSynonymousCodonUsage(pseudocount) / cubar est_rscu(pseudo_cnt) formula).
[TestFixture]
public class RscuPseudocountToolTests
{
    [Test]
    public void Rscu_Pseudocount_HandDerived()
    {
        // Phe TTT×2 TTC×1, c = 1: TTT = 2·3/(3+2) = 1.2, TTC = 2·2/5 = 0.8; absent Met → 1.
        var rscu = MolToolsTools.rscu("TTTTTTTTC", pseudocount: 1).Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(rscu["TTT"], Is.EqualTo(1.2).Within(1e-12));
            Assert.That(rscu["TTC"], Is.EqualTo(0.8).Within(1e-12));
            Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(1e-12));
        });
    }

    [Test]
    public void Rscu_Pseudocount_DelegatesToCore_WithGeneticCode()
    {
        const string seq = "ATAATGATGTGATGGTGGAGATAA";
        var expected = CodonUsageAnalyzer.CalculateRscu(seq, GeneticCode.GetByTableNumber(2), 0.5);
        var actual = MolToolsTools.rscu(seq, genetic_code: 2, pseudocount: 0.5).Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(actual, Is.EqualTo(expected));
            Assert.That(actual["ATA"], Is.EqualTo(0.75).Within(1e-12), "Met(ATA,ATG) table 2: 2·1.5/4");
        });
    }

    [Test]
    public void Rscu_DefaultPseudocount_IsZero()
    {
        var rscu = MolToolsTools.rscu("TTTTTTTTC").Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.rscu("TTTTTTTTC", pseudocount: 0).Rscu, Is.EqualTo(rscu));
            Assert.That(rscu["ATG"], Is.EqualTo(0.0), "absent family 0 (CodonW) without a pseudocount");
        });
    }

    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void Rscu_InvalidPseudocount_Throws(double pc) =>
        Assert.Throws<ArgumentException>(() => MolToolsTools.rscu("TTT", pseudocount: pc));
}
