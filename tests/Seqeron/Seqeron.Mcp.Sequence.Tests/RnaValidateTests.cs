using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class RnaValidateTests
{
    [Test]
    public void RnaValidate_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.RnaValidate("AUGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.RnaValidate(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.RnaValidate(null!));
    }

    [Test]
    public void RnaValidate_Binding_InvokesSuccessfully()
    {
        var valid = SequenceTools.RnaValidate("AUGCAUGC");
        Assert.That(valid.Valid, Is.True);
        Assert.That(valid.Length, Is.EqualTo(8));
        Assert.That(valid.Error, Is.Null);

        var invalid = SequenceTools.RnaValidate("AUGTATGC"); // T is invalid in RNA
        Assert.That(invalid.Valid, Is.False);
        Assert.That(invalid.Error, Does.Contain("T"));
    }

    // iupac=true delegates to Core SequenceExtensions.IndexOfInvalidIupacRna: Biopython 1.88
    // IUPACData.ambiguous_rna_letters "GAUCRYWSMKHBVDN" (scikit-bio 0.7.4 RNA("ACGUN") valid, RNA("ACGTN") ValueError).
    [Test]
    public void RnaValidate_Iupac_AcceptsAmbiguityCodes()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.RnaValidate("ACGUN").Valid, Is.False, "strict default unchanged");
            Assert.That(SequenceTools.RnaValidate("ACGURYSWKMBDHVN", iupac: true).Valid, Is.True);
            var t = SequenceTools.RnaValidate("ACGTN", iupac: true);
            Assert.That(t.Valid, Is.False);
            Assert.That(t.Error, Is.EqualTo("Invalid nucleotide 'T' at position 3"));
            Assert.Throws<ArgumentException>(() => SequenceTools.RnaValidate("", iupac: true));
        });
    }
}
