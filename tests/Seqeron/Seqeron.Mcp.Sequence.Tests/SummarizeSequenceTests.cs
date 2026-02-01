using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class SummarizeSequenceTests
{
    [Test]
    public void SummarizeSequence_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.SummarizeSequence("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.SummarizeSequence(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.SummarizeSequence(null!));
    }

    [Test]
    public void SummarizeSequence_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.SummarizeSequence("ATGCGATCGATCG");
        Assert.That(result.Length, Is.EqualTo(13));
        Assert.That(result.GcContent, Is.InRange(0, 100));
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.Complexity, Is.InRange(0, 1));
        Assert.That(result.MeltingTemperature, Is.TypeOf<double>());
        Assert.That(result.Composition, Contains.Key("A"));
        Assert.That(result.Composition, Contains.Key("T"));
        Assert.That(result.Composition, Contains.Key("G"));
        Assert.That(result.Composition, Contains.Key("C"));
    }
}
