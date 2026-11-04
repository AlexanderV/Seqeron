using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ProteinValidateTests
{
    [Test]
    public void ProteinValidate_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ProteinValidate("MAEGEITTFT"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ProteinValidate(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ProteinValidate(null!));
    }

    [Test]
    public void ProteinValidate_Binding_InvokesSuccessfully()
    {
        var valid = SequenceTools.ProteinValidate("MAEGEITTFT");
        Assert.That(valid.Valid, Is.True);
        Assert.That(valid.Length, Is.EqualTo(10));
        Assert.That(valid.Error, Is.Null);

        var invalid = SequenceTools.ProteinValidate("MAEGEITTF1"); // '1' is not an amino-acid code
        Assert.That(invalid.Valid, Is.False);
        Assert.That(invalid.Error, Does.Contain("1"));

        // IUPAC ambiguity codes B (Asx), Z (Glx), J (Xle) are valid (Biopython extended_protein_letters).
        Assert.That(SequenceTools.ProteinValidate("MAEGBZJ").Valid, Is.True);

        var withStop = SequenceTools.ProteinValidate("MAEG*");
        Assert.That(withStop.Valid, Is.True); // * is valid stop codon
    }

    // FIN-B02 WP5: the offending character is located with ASCII-only folding (as ProteinSequence,
    // F31), so a non-ASCII letter that ToUpperInvariant would fold onto a valid residue
    // (U+017F 'ſ' -> 'S', U+0131 'ı' -> 'I') is reported, not a generic "Invalid sequence".
    [TestCase("MAſEG", "\u017f", 2)]
    [TestCase("ıMAEG", "\u0131", 0)]
    public void ProteinValidate_NonAsciiLetter_IsReportedVerbatim(string sequence, string ch, int position)
    {
        var result = SequenceTools.ProteinValidate(sequence);
        Assert.Multiple(() =>
        {
            Assert.That(result.Valid, Is.False);
            Assert.That(result.Error, Is.EqualTo($"Invalid amino acid '{ch}' at position {position}"));
        });
    }
}
