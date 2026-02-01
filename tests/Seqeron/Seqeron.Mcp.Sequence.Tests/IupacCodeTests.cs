using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class IupacCodeTests
{
    [Test]
    public void IupacCode_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IupacCode("AG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacCode(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IupacCode(null!));
    }

    [Test]
    public void IupacCode_Binding_InvokesSuccessfully()
    {
        // AG = Purine = R
        var result = SequenceTools.IupacCode("AG");
        Assert.That(result.Code, Is.EqualTo("R"));
        Assert.That(result.InputBases, Is.EqualTo("AG"));

        // CT = Pyrimidine = Y
        var result2 = SequenceTools.IupacCode("CT");
        Assert.That(result2.Code, Is.EqualTo("Y"));

        // Single base
        var result3 = SequenceTools.IupacCode("A");
        Assert.That(result3.Code, Is.EqualTo("A"));
    }
}
