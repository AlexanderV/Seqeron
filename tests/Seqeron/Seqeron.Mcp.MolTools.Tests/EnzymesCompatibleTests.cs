using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EnzymesCompatibleTests
{
    [Test]
    public void EnzymesCompatible_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.enzymes_compatible("EcoRV", "SmaI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_compatible("", "SmaI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_compatible("EcoRV", null!));
    }

    [Test]
    public void EnzymesCompatible_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Two blunt cutters -> compatible.
            Assert.That(MolToolsTools.enzymes_compatible("EcoRV", "SmaI").Compatible, Is.True);

            // BamHI and BglII share a 5' GATC overhang -> compatible.
            Assert.That(MolToolsTools.enzymes_compatible("BamHI", "BglII").Compatible, Is.True);

            // EcoRI (AATT overhang) vs BamHI (GATC overhang) -> not compatible.
            Assert.That(MolToolsTools.enzymes_compatible("EcoRI", "BamHI").Compatible, Is.False);

            // Unknown enzyme name -> false (no exception).
            Assert.That(MolToolsTools.enzymes_compatible("NotAnEnzyme", "EcoRI").Compatible, Is.False);
        });
    }
}
