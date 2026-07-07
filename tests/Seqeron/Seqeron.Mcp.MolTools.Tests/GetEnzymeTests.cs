using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class GetEnzymeTests
{
    [Test]
    public void GetEnzyme_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.get_enzyme("EcoRI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme("   "));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme(null!));
    }

    [Test]
    public void GetEnzyme_Binding_InvokesSuccessfully()
    {
        // EcoRI: GAATTC, cut forward 1 / reverse 5, Escherichia coli.
        var ecoRI = MolToolsTools.get_enzyme("EcoRI").Enzyme;
        Assert.Multiple(() =>
        {
            Assert.That(ecoRI, Is.Not.Null);
            Assert.That(ecoRI!.Name, Is.EqualTo("EcoRI"));
            Assert.That(ecoRI.RecognitionSequence, Is.EqualTo("GAATTC"));
            Assert.That(ecoRI.CutPositionForward, Is.EqualTo(1));
            Assert.That(ecoRI.CutPositionReverse, Is.EqualTo(5));
            Assert.That(ecoRI.Organism, Is.EqualTo("Escherichia coli"));
        });

        // Lookup is case-insensitive.
        Assert.That(MolToolsTools.get_enzyme("ecori").Enzyme, Is.Not.Null);

        // Unknown enzyme -> null (documented, not an error).
        Assert.That(MolToolsTools.get_enzyme("NotAnEnzyme").Enzyme, Is.Null);
    }
}
