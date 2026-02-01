using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class IupacMatchesTests
{
    [Test]
    public void IupacMatches_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IupacMatches("A", "R"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatches("", "R"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatches("A", ""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatches(null!, "R"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatches("A", null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatches("AG", "R")); // must be single char
    }

    [Test]
    public void IupacMatches_Binding_InvokesSuccessfully()
    {
        // A matches R (purine)
        var result = SequenceTools.IupacMatches("A", "R");
        Assert.That(result.Matches, Is.True);
        Assert.That(result.Nucleotide, Is.EqualTo("A"));
        Assert.That(result.IupacCode, Is.EqualTo("R"));

        // G matches R (purine)
        var result2 = SequenceTools.IupacMatches("G", "R");
        Assert.That(result2.Matches, Is.True);

        // C does not match R (purine)
        var result3 = SequenceTools.IupacMatches("C", "R");
        Assert.That(result3.Matches, Is.False);

        // Any nucleotide matches N
        var result4 = SequenceTools.IupacMatches("A", "N");
        Assert.That(result4.Matches, Is.True);
    }
}
