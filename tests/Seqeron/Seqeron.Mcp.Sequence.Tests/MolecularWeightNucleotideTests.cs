using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class MolecularWeightNucleotideTests
{
    [Test]
    public void MolecularWeightNucleotide_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.MolecularWeightNucleotide("ATGC", true));
        Assert.DoesNotThrow(() => SequenceTools.MolecularWeightNucleotide("AUGC", false));
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightNucleotide("", true));
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightNucleotide(null!, true));
    }

    [Test]
    public void MolecularWeightNucleotide_Binding_InvokesSuccessfully()
    {
        var dna = SequenceTools.MolecularWeightNucleotide("ATGC", true);
        Assert.That(dna.MolecularWeight, Is.GreaterThan(0));
        Assert.That(dna.Unit, Is.EqualTo("Da"));
        Assert.That(dna.SequenceType, Is.EqualTo("DNA"));

        var rna = SequenceTools.MolecularWeightNucleotide("AUGC", false);
        Assert.That(rna.MolecularWeight, Is.GreaterThan(0));
        Assert.That(rna.SequenceType, Is.EqualTo("RNA"));
    }
}
