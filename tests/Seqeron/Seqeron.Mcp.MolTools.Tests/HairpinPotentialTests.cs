using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class HairpinPotentialTests
{
    [Test]
    public void HairpinPotential_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.hairpin_potential("GGGGAAACCCC"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.hairpin_potential(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.hairpin_potential(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.hairpin_potential("GGGGAAACCCC", 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.hairpin_potential("GGGGAAACCCC", 4, -1));
    }

    [Test]
    public void HairpinPotential_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // GGGG (stem) + AAA (loop) + CCCC (complementary stem) -> hairpin.
            Assert.That(MolToolsTools.hairpin_potential("GGGGAAACCCC").HasHairpin, Is.True);

            // No self-complementary stem -> no hairpin.
            Assert.That(MolToolsTools.hairpin_potential("AAAAAAAAAAA").HasHairpin, Is.False);

            // Shorter than 2*stem + loop (= 11 for defaults) -> false.
            Assert.That(MolToolsTools.hairpin_potential("GGGGCCCC").HasHairpin, Is.False);
        });
    }
}
