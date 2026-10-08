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
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightProtein("MAEG1")); // '1' is not an amino-acid code (J = Xle is a valid IUPAC code)
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
