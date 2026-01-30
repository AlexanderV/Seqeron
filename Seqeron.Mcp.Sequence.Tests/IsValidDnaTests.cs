using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class IsValidDnaTests
{
    [Test]
    public void IsValidDna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IsValidDna("ATGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsValidDna(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsValidDna(null!));
    }

    [Test]
    public void IsValidDna_Binding_InvokesSuccessfully()
    {
        // Valid DNA
        var valid = SequenceTools.IsValidDna("ATGCGATCGATCG");
        Assert.That(valid.IsValid, Is.True);
        Assert.That(valid.Length, Is.EqualTo(13));

        // Invalid DNA (contains U - RNA nucleotide)
        var invalid = SequenceTools.IsValidDna("AUGC");
        Assert.That(invalid.IsValid, Is.False);

        // Invalid DNA (contains invalid character)
        var invalidChar = SequenceTools.IsValidDna("ATGX");
        Assert.That(invalidChar.IsValid, Is.False);

        // Case insensitive
        var lowercase = SequenceTools.IsValidDna("atgc");
        Assert.That(lowercase.IsValid, Is.True);
    }
}
