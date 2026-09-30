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

    // 2026-09 B03 F21: linguistic_complexity now delegates to the canonical SequenceComplexity sum form,
    // so it agrees with complexity_linguistic for the same word length (only the defaults differ: 6 vs 10).
    // Reference: Rosalind LING sample ATTTGGATT → 0.875 (m = N = 9); m = 6 → 29/34 (Python reference).
    [TestCase("ATTTGGATT", 9, 0.875)]
    [TestCase("ATTTGGATT", 6, 29.0 / 34.0)]
    [TestCase("AAAAAAAAAA", 6, 2.0 / 13.0)]
    [TestCase("ATGCATGCATGC", 6, 24.0 / 49.0)]
    public void LinguisticComplexity_AgreesWithComplexityLinguistic_AndSumFormReference(string seq, int m, double expected)
    {
        var a = SequenceTools.LinguisticComplexity(seq, m).Complexity;
        var b = SequenceTools.ComplexityLinguistic(seq, m).Complexity;
        Assert.That(a, Is.EqualTo(b));
        Assert.That(a, Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void LinguisticComplexity_AlphabetSize_FixedAlphabet()
    {
        // Protein, a = 20: MKVLAAGIVGLLLAA, m = 15 -> 9/10 (exact brute force; B04 F37).
        Assert.That(SequenceTools.LinguisticComplexity("MKVLAAGIVGLLLAA", 15, alphabetSize: 20).Complexity, Is.EqualTo(0.9).Within(1e-15));
        Assert.Throws<ArgumentException>(() => SequenceTools.LinguisticComplexity("ACGTN", 6, alphabetSize: 4));
    }
}
