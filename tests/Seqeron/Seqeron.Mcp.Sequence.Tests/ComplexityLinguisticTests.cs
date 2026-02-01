using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityLinguisticTests
{
    [Test]
    public void ComplexityLinguistic_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityLinguistic("ATGCGATCGATCG", 10));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityLinguistic("", 10));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityLinguistic(null!, 10));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityLinguistic("ATGC", 0)); // maxWordLength < 1
    }

    [Test]
    public void ComplexityLinguistic_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityLinguistic("ATGCGATCGATCG", 10);
        Assert.That(result.Complexity, Is.GreaterThan(0));
        Assert.That(result.Complexity, Is.LessThanOrEqualTo(1));
        Assert.That(result.MaxWordLength, Is.EqualTo(10));

        // Low complexity (repetitive) sequence has lower linguistic complexity
        var lowComplexity = SequenceTools.ComplexityLinguistic("AAAAAAAAAAAAAAAA", 10);
        var highComplexity = SequenceTools.ComplexityLinguistic("ATGCATGCATGCATGC", 10);
        Assert.That(highComplexity.Complexity, Is.GreaterThan(lowComplexity.Complexity));
    }
}
