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

        var invalid = SequenceTools.ProteinValidate("MAEGEITTFJ"); // J is invalid
        Assert.That(invalid.Valid, Is.False);
        Assert.That(invalid.Error, Does.Contain("J"));

        var withStop = SequenceTools.ProteinValidate("MAEG*");
        Assert.That(withStop.Valid, Is.True); // * is valid stop codon
    }
}
