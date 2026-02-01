using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class AminoAcidCompositionTests
{
    [Test]
    public void AminoAcidComposition_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.AminoAcidComposition("MAEGEITTFT"));
        Assert.Throws<ArgumentException>(() => SequenceTools.AminoAcidComposition(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.AminoAcidComposition(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.AminoAcidComposition("MAEGJ")); // J is invalid
    }

    [Test]
    public void AminoAcidComposition_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.AminoAcidComposition("MAEGEITTFT");
        Assert.That(result.Length, Is.EqualTo(10));
        Assert.That(result.Counts, Contains.Key("M"));
        Assert.That(result.Counts["M"], Is.EqualTo(1));
        Assert.That(result.Counts["T"], Is.EqualTo(3));
        Assert.That(result.MolecularWeight, Is.GreaterThan(0));
        Assert.That(result.IsoelectricPoint, Is.InRange(0, 14));
    }
}
