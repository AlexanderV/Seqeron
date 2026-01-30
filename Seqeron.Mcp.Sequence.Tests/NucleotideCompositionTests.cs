using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class NucleotideCompositionTests
{
    [Test]
    public void NucleotideComposition_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.NucleotideComposition("ATGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.NucleotideComposition(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.NucleotideComposition(null!));
    }

    [Test]
    public void NucleotideComposition_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.NucleotideComposition("ATGCATGC");
        Assert.That(result.Length, Is.EqualTo(8));
        Assert.That(result.A, Is.EqualTo(2));
        Assert.That(result.T, Is.EqualTo(2));
        Assert.That(result.G, Is.EqualTo(2));
        Assert.That(result.C, Is.EqualTo(2));
        Assert.That(result.GcContent, Is.EqualTo(0.5));

        var rna = SequenceTools.NucleotideComposition("AUGCAUGC");
        Assert.That(rna.U, Is.EqualTo(2));
        Assert.That(rna.T, Is.EqualTo(0));
    }
}
