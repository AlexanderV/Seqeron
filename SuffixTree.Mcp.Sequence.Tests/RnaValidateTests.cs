using NUnit.Framework;
using SuffixTree.Mcp.Sequence.Tools;

namespace SuffixTree.Mcp.Sequence.Tests;

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
}
