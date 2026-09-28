using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class HydrophobicityTests
{
    [Test]
    public void Hydrophobicity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.Hydrophobicity("MAEGEITTFT"));
        Assert.Throws<ArgumentException>(() => SequenceTools.Hydrophobicity(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.Hydrophobicity(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.Hydrophobicity("MAEG1")); // '1' is not an amino-acid code (J = Xle is a valid IUPAC code)
    }

    [Test]
    public void Hydrophobicity_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.Hydrophobicity("MAEGEITTFT");
        Assert.That(result.Gravy, Is.TypeOf<double>());

        // Test hydrophobic protein (ILV-rich)
        var hydrophobicResult = SequenceTools.Hydrophobicity("ILVILVILVV");
        Assert.That(hydrophobicResult.Gravy, Is.GreaterThan(0));

        // Test hydrophilic protein (DERK-rich)
        var hydrophilicResult = SequenceTools.Hydrophobicity("DERKDERKDE");
        Assert.That(hydrophilicResult.Gravy, Is.LessThan(0));
    }
}
