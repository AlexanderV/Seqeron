using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class MolecularWeightProteinTests
{
    [Test]
    public void MolecularWeightProtein_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.MolecularWeightProtein("MAEGEITTFT"));
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightProtein(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightProtein(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightProtein("MAEGJ")); // J is invalid
    }

    [Test]
    public void MolecularWeightProtein_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.MolecularWeightProtein("MAEGEITTFT");
        Assert.That(result.MolecularWeight, Is.GreaterThan(0));
        Assert.That(result.Unit, Is.EqualTo("Da"));

        // Single amino acid - just check it returns a positive value
        var single = SequenceTools.MolecularWeightProtein("M");
        Assert.That(single.MolecularWeight, Is.GreaterThan(0));
    }
}
