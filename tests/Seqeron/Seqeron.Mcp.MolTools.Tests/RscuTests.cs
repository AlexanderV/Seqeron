using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RscuTests
{
    [Test]
    public void Rscu_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.rscu("TTTTTTTTC"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.rscu(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.rscu(null!));
    }

    [Test]
    public void Rscu_Binding_InvokesSuccessfully()
    {
        // "TTTTTTTTC" -> TTT,TTT,TTC. Phe family {TTT,TTC}, total 3, expected = 3/2 = 1.5.
        //   RSCU[TTT] = 2 / 1.5 = 1.333...
        //   RSCU[TTC] = 1 / 1.5 = 0.666...
        var rscu = MolToolsTools.rscu("TTTTTTTTC").Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(rscu["TTT"], Is.EqualTo(2.0 / 1.5).Within(1e-9));
            Assert.That(rscu["TTC"], Is.EqualTo(1.0 / 1.5).Within(1e-9));
        });

        // Equal usage -> RSCU 1.0 each. "TTTTTC" -> TTT,TTC, total 2, expected 1.0.
        var equal = MolToolsTools.rscu("TTTTTC").Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(equal["TTT"], Is.EqualTo(1.0).Within(1e-9));
            Assert.That(equal["TTC"], Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
