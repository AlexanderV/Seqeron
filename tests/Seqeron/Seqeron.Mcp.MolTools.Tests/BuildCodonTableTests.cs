using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class BuildCodonTableTests
{
    [Test]
    public void BuildCodonTable_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.build_codon_table("TTTTTTTTC", "MyOrg"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.build_codon_table("", "MyOrg"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.build_codon_table(null!, "MyOrg"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.build_codon_table("TTT", ""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.build_codon_table("TTT", null!));
    }

    [Test]
    public void BuildCodonTable_Binding_InvokesSuccessfully()
    {
        // Reference "TTTTTTTTC" -> codons TTT,TTT,TTC -> RNA UUU,UUU,UUC.
        // UUU and UUC both code Phe(F). Phe total = 3.
        //   freq[UUU] = 2/3, freq[UUC] = 1/3.
        var result = MolToolsTools.build_codon_table("TTTTTTTTC", "MyOrg");

        Assert.Multiple(() =>
        {
            Assert.That(result.OrganismName, Is.EqualTo("MyOrg"));
            Assert.That(result.CodonFrequencies["UUU"], Is.EqualTo(2.0 / 3.0).Within(1e-9));
            Assert.That(result.CodonFrequencies["UUC"], Is.EqualTo(1.0 / 3.0).Within(1e-9));
            // Frequencies are computed on the RNA alphabet (T replaced by U).
            Assert.That(result.CodonFrequencies.ContainsKey("TTT"), Is.False);
            // Genetic-code map assigns Phe to both codons.
            Assert.That(result.CodonToAminoAcid["UUU"], Is.EqualTo("F"));
            Assert.That(result.CodonToAminoAcid["UUC"], Is.EqualTo("F"));
        });
    }
}
