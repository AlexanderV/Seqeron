using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class DnaValidateTests
{
    [Test]
    public void DnaValidate_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.DnaValidate("ATGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.DnaValidate(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.DnaValidate(null!));
    }

    [Test]
    public void DnaValidate_Binding_InvokesSuccessfully()
    {
        var valid = SequenceTools.DnaValidate("ATGCATGC");
        Assert.That(valid.Valid, Is.True);
        Assert.That(valid.Length, Is.EqualTo(8));
        Assert.That(valid.Error, Is.Null);

        var invalid = SequenceTools.DnaValidate("ATGXATGC");
        Assert.That(invalid.Valid, Is.False);
        Assert.That(invalid.Length, Is.EqualTo(8));
        Assert.That(invalid.Error, Does.Contain("X"));
    }

    // iupac=true delegates to Core SequenceExtensions.IndexOfInvalidIupacDna: the alphabet is Biopython 1.88
    // IUPACData.ambiguous_dna_letters "GATCRYWSMKHBVDN" (scikit-bio 0.7.4 DNA("ACGTN") also valid, DNA("ACGTU") ValueError).
    [Test]
    public void DnaValidate_Iupac_AcceptsAmbiguityCodes()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.DnaValidate("ACGTN").Valid, Is.False, "strict default unchanged");
            Assert.That(SequenceTools.DnaValidate("ACGTN").Error, Is.EqualTo("Invalid nucleotide 'N' at position 4"));
            Assert.That(SequenceTools.DnaValidate("ACGTRYSWKMBDHVN", iupac: true).Valid, Is.True);
            Assert.That(SequenceTools.DnaValidate("acgtn", iupac: true).Valid, Is.True);
            var u = SequenceTools.DnaValidate("ACGNU", iupac: true);
            Assert.That(u.Valid, Is.False);
            Assert.That(u.Error, Is.EqualTo("Invalid nucleotide 'U' at position 4"));
            Assert.That(SequenceTools.DnaValidate("ACX", iupac: true).Error, Is.EqualTo("Invalid nucleotide 'X' at position 2"));
            // Empty input stays rejected (MCP null/empty convention) in both modes.
            Assert.Throws<ArgumentException>(() => SequenceTools.DnaValidate("", iupac: true));
        });
    }
}
