using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CodonAdaptationIndexTests
{
    private static Dictionary<string, double> PheRscu() =>
        new() { ["TTT"] = 1.0, ["TTC"] = 0.5 };

    [Test]
    public void CodonAdaptationIndex_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.codon_adaptation_index("TTT", PheRscu()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_adaptation_index("", PheRscu()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_adaptation_index(null!, PheRscu()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.codon_adaptation_index("TTT", null!));
    }

    [Test]
    public void CodonAdaptationIndex_Binding_InvokesSuccessfully()
    {
        var rscu = PheRscu(); // Phe: max RSCU = 1.0 => w[TTT]=1.0, w[TTC]=0.5.

        // "TTC" -> w = 0.5 -> CAI = 0.5.
        var single = MolToolsTools.codon_adaptation_index("TTC", rscu);
        // "TTTTTC" -> {1.0, 0.5} -> CAI = sqrt(0.5) = 0.70710678...
        var pair = MolToolsTools.codon_adaptation_index("TTTTTC", rscu);

        Assert.Multiple(() =>
        {
            Assert.That(single.Cai, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(pair.Cai, Is.EqualTo(Math.Sqrt(0.5)).Within(1e-9));
        });
    }
}
