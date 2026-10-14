using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

// CODON-RSCU-001 (audit A3-4): optional `family_level` on the MCP `rscu` tool
// (cubar est_rscu(level = "subfam"): subfam = amino acid + first two codon bases).
[TestFixture]
public class RscuFamilyLevelToolTests
{
    [Test]
    public void Rscu_Subfamily_HandDerived()
    {
        // CTG×3 CTA×1 TTA×1: Leu_CT = 4·3/4 = 3, 4·1/4 = 1; Leu_TT = 2·1/1 = 2. Amino-acid level: CTG 6·3/5 = 3.6.
        var sub = MolToolsTools.rscu("CTGCTGCTGCTATTA", family_level: RscuFamilyLevel.Subfamily).Rscu;
        var aa = MolToolsTools.rscu("CTGCTGCTGCTATTA").Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(sub["CTG"], Is.EqualTo(3.0).Within(1e-12));
            Assert.That(sub["CTA"], Is.EqualTo(1.0).Within(1e-12));
            Assert.That(sub["TTA"], Is.EqualTo(2.0).Within(1e-12));
            Assert.That(aa["CTG"], Is.EqualTo(3.6).Within(1e-12), "default stays the amino-acid family");
        });
    }

    [Test]
    public void Rscu_Subfamily_DelegatesToCore_WithGeneticCodeAndPseudocount()
    {
        const string seq = "ATAATGATGTGATGGTGGAGATAATCAAGC";
        var expected = CodonUsageAnalyzer.CalculateRscu(seq, GeneticCode.GetByTableNumber(2), 1.0,
            RscuFamilyLevel.Subfamily);
        var actual = MolToolsTools.rscu(seq, genetic_code: 2, pseudocount: 1, family_level: RscuFamilyLevel.Subfamily).Rscu;
        Assert.That(actual, Is.EqualTo(expected));
    }

    [Test]
    public void Rscu_UnknownFamilyLevel_Throws() =>
        Assert.Throws<ArgumentException>(() => MolToolsTools.rscu("TTT", family_level: (RscuFamilyLevel)9));
}
