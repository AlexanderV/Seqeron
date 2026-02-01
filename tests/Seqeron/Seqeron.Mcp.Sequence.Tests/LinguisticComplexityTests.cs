using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class LinguisticComplexityTests
{
    [Test]
    public void LinguisticComplexity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.LinguisticComplexity("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.LinguisticComplexity(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.LinguisticComplexity(null!));
    }

    [Test]
    public void LinguisticComplexity_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.LinguisticComplexity("ATGCGATCGATCG");
        Assert.That(result.Complexity, Is.TypeOf<double>());
        Assert.That(result.Complexity, Is.InRange(0, 1));

        // Low complexity sequence (homopolymer) has lower complexity
        var lowComplexity = SequenceTools.LinguisticComplexity("AAAAAAAAAA");
        var highComplexity = SequenceTools.LinguisticComplexity("ATGCATGCATGC");
        Assert.That(highComplexity.Complexity, Is.GreaterThan(lowComplexity.Complexity));
    }
}
