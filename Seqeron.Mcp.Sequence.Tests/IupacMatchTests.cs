using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class IupacMatchTests
{
    [Test]
    public void IupacMatch_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IupacMatch("R", "A"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatch("", "A"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatch("A", ""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatch(null!, "A"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatch("A", null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacMatch("AG", "A")); // must be single char
    }

    [Test]
    public void IupacMatch_Binding_InvokesSuccessfully()
    {
        // R (purine = A or G) matches A
        var result = SequenceTools.IupacMatch("R", "A");
        Assert.That(result.Matches, Is.True);
        Assert.That(result.Code1, Is.EqualTo("R"));
        Assert.That(result.Code2, Is.EqualTo("A"));

        // R matches G
        var result2 = SequenceTools.IupacMatch("R", "G");
        Assert.That(result2.Matches, Is.True);

        // R does not match Y (pyrimidine = C or T)
        var result3 = SequenceTools.IupacMatch("R", "Y");
        Assert.That(result3.Matches, Is.False);
    }
}
