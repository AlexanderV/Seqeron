using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class IsValidRnaTests
{
    [Test]
    public void IsValidRna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IsValidRna("AUGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsValidRna(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsValidRna(null!));
    }

    [Test]
    public void IsValidRna_Binding_InvokesSuccessfully()
    {
        // Valid RNA
        var valid = SequenceTools.IsValidRna("AUGCGAUCGAUCG");
        Assert.That(valid.IsValid, Is.True);
        Assert.That(valid.Length, Is.EqualTo(13));

        // Invalid RNA (contains T - DNA nucleotide)
        var invalid = SequenceTools.IsValidRna("ATGC");
        Assert.That(invalid.IsValid, Is.False);

        // Invalid RNA (contains invalid character)
        var invalidChar = SequenceTools.IsValidRna("AUGX");
        Assert.That(invalidChar.IsValid, Is.False);

        // Case insensitive
        var lowercase = SequenceTools.IsValidRna("augc");
        Assert.That(lowercase.IsValid, Is.True);
    }
}
