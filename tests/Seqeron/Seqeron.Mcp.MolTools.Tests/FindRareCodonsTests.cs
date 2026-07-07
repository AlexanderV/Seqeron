using NUnit.Framework;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindRareCodonsTests
{
    private static CodonUsageTableInput PheTable() =>
        new(CodonFrequencies: new Dictionary<string, double> { ["UUU"] = 0.8, ["UUC"] = 0.1 });

    [Test]
    public void FindRareCodons_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_rare_codons("TTTTTC", PheTable()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_rare_codons("", PheTable()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_rare_codons(null!, PheTable()));
        Assert.Throws<ArgumentNullException>(() => MolToolsTools.find_rare_codons("TTTTTC", null!));
    }

    [Test]
    public void FindRareCodons_Binding_InvokesSuccessfully()
    {
        // Table: UUU=0.8 (common), UUC=0.1 (< 0.15 threshold, rare).
        // "TTTTTC" -> UUU (pos 0, kept), UUC (pos 3, rare).
        var rare = MolToolsTools.find_rare_codons("TTTTTC", PheTable()).RareCodons;

        Assert.Multiple(() =>
        {
            Assert.That(rare, Has.Count.EqualTo(1));
            Assert.That(rare[0].Position, Is.EqualTo(3));
            Assert.That(rare[0].Codon, Is.EqualTo("UUC"));
            Assert.That(rare[0].AminoAcid, Is.EqualTo("F"));
            Assert.That(rare[0].Frequency, Is.EqualTo(0.1).Within(1e-9));
        });

        // Raising the threshold above 0.8 makes even UUU rare (both codons reported).
        var both = MolToolsTools.find_rare_codons("TTTTTC", PheTable(), threshold: 0.9).RareCodons;
        Assert.That(both, Has.Count.EqualTo(2));
    }
}
