using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Pins the two published on-target efficacy scorers exposed by the MolTools server. The expected
/// values are the reference ones locked by Seqeron.Genomics.Tests (Doench 2014: the two worked
/// examples distributed with CRISPOR's doenchScore.py; Rule Set 2: a row of the verified Azimuth
/// oracle on which the reference and Microsoft's own fixture agree).
/// </summary>
[TestFixture]
public class CalculateOnTargetScoreTests
{
    private const string Doench2014Example1 = "TATAGCTGCGATCTGAGGTAGGGAGGGACC"; // reference 0.713089368437
    private const string Doench2014Example2 = "TCCGCACCTGTCACGGTCGGGGCTTGGCGC"; // reference 0.0189838463593
    private const string Oracle30Mer = "AAAAAAAAAAAAAAAGCTAACAGCAGGAGT"; // nopos_oracle ref_score 0.529037

    [Test]
    public void CalculateOnTargetDoench2014_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_doench2014(Doench2014Example1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014("ACGTACGT"));
        // 30-mer without an NGG PAM at offsets 25-26.
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_doench2014("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    }

    [Test]
    public void CalculateOnTargetDoench2014_Binding_ReproducesTheReferenceExamples()
    {
        Assert.That(MolToolsTools.calculate_on_target_doench2014(Doench2014Example1).Score,
            Is.EqualTo(71.3089368437).Within(1e-4));
        Assert.That(MolToolsTools.calculate_on_target_doench2014(Doench2014Example2).Score,
            Is.EqualTo(1.89838463593).Within(1e-4));
    }

    [Test]
    public void CalculateOnTargetRuleSet2_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer));
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, 25.0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2("ACGTACGT"));
        // The gene-context arguments must be supplied together.
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, null));
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, null, 25.0));
    }

    [Test]
    public void CalculateOnTargetRuleSet2_Binding_ReproducesTheVerifiedOracleRow()
    {
        Assert.That(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer).Score,
            Is.EqualTo(0.529037).Within(1e-5));

        // The gene-context (full) model is a different model and must give a different value.
        Assert.That(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, 25.0).Score,
            Is.Not.EqualTo(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer).Score));
    }
}
