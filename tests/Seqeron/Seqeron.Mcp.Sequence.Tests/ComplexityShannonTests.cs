using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityShannonTests
{
    [Test]
    public void ComplexityShannon_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityShannon("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityShannon(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityShannon(null!));
    }

    [Test]
    public void ComplexityShannon_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityShannon("ATGCGATCGATCG");
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.Entropy, Is.LessThanOrEqualTo(2)); // Max for DNA is 2 bits

        // Low complexity (repetitive) sequence has lower entropy
        var lowComplexity = SequenceTools.ComplexityShannon("AAAAAAAAAAAAAAAA");
        var highComplexity = SequenceTools.ComplexityShannon("ATGCATGCATGCATGC");
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }

    /// <summary>
    /// A1-12 / F25: documented alphabet A/C/G/T/U with U counted as T, other symbols ignored.
    /// Reference: scipy 1.x <c>scipy.stats.entropy(counts, base=2)</c> over {A,C,G,T(+U)}:
    /// [1,1,1,1] → 2.0, [3,3] → 1.0, [2,1,1] → 1.5.
    /// </summary>
    [TestCase("ACGU", 2.0)]
    [TestCase("acgu", 2.0)]
    [TestCase("AUAUAU", 1.0)]
    [TestCase("ACGTN", 2.0)]     // N ignored
    [TestCase("AACG-", 1.5)]     // gap ignored
    public void ComplexityShannon_AcgtuAlphabet_UAsT_MatchesScipy(string seq, double expected)
    {
        Assert.That(SequenceTools.ComplexityShannon(seq).Entropy, Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void ComplexityShannon_RnaEqualsDnaTranscript()
    {
        Assert.That(SequenceTools.ComplexityShannon("AUGGCUUAAC").Entropy,
            Is.EqualTo(SequenceTools.ComplexityShannon("ATGGCTTAAC").Entropy));
    }
}
