using NUnit.Framework;
using SuffixTree.Mcp.Sequence.Tools;

namespace SuffixTree.Mcp.Sequence.Tests;

[TestFixture]
public class IsoelectricPointTests
{
    [Test]
    public void IsoelectricPoint_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.IsoelectricPoint("MAEGEITTFT"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("MAEGJ")); // J is invalid
    }

    [Test]
    public void IsoelectricPoint_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.IsoelectricPoint("MAEGEITTFT");
        Assert.That(result.PI, Is.InRange(0, 14));

        // Test acidic protein (more D, E)
        var acidicResult = SequenceTools.IsoelectricPoint("DDDEEE");
        Assert.That(acidicResult.PI, Is.LessThan(7.0));

        // Test basic protein (more K, R)
        var basicResult = SequenceTools.IsoelectricPoint("KKKRRR");
        Assert.That(basicResult.PI, Is.GreaterThan(7.0));
    }
}
