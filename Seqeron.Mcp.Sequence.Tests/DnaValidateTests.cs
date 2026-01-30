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
}
