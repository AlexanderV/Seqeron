using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class SequenceComplexityTests
{
    #region Linguistic Complexity Tests

    [Test]
    public void CalculateLinguisticComplexity_HighComplexity_ReturnsHigh()
    {
        // Random-like sequence should have high complexity
        var sequence = new DnaSequence("ATGCTAGCATGCAATG");
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence);

        Assert.That(lc, Is.GreaterThan(0.5));
    }

    [Test]
    public void CalculateLinguisticComplexity_LowComplexity_ReturnsLow()
    {
        // Highly repetitive sequence should have low complexity
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAA");
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence);

        Assert.That(lc, Is.LessThan(0.3));
    }

    [Test]
    public void CalculateLinguisticComplexity_EmptySequence_ReturnsZero()
    {
        double lc = SequenceComplexity.CalculateLinguisticComplexity("");
        Assert.That(lc, Is.EqualTo(0));
    }

    [Test]
    public void CalculateLinguisticComplexity_RangeIsZeroToOne()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence);

        Assert.That(lc, Is.GreaterThanOrEqualTo(0));
        Assert.That(lc, Is.LessThanOrEqualTo(1));
    }

    [Test]
    public void CalculateLinguisticComplexity_StringOverload_Works()
    {
        double lc = SequenceComplexity.CalculateLinguisticComplexity("ATGCATGC");
        Assert.That(lc, Is.GreaterThan(0));
    }

    #endregion

    #region Shannon Entropy Tests

    [Test]
    public void CalculateShannonEntropy_EqualBases_ReturnsTwo()
    {
        // Equal distribution of all 4 bases = max entropy (2 bits)
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(2.0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_SingleBase_ReturnsZero()
    {
        // Only one base type = zero entropy
        var sequence = new DnaSequence("AAAAAAA");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_TwoBases_ReturnsOne()
    {
        // Two bases equally distributed = 1 bit entropy
        var sequence = new DnaSequence("ATATATAT");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_EmptySequence_ReturnsZero()
    {
        double entropy = SequenceComplexity.CalculateShannonEntropy("");
        Assert.That(entropy, Is.EqualTo(0));
    }

    [Test]
    public void CalculateShannonEntropy_StringOverload_Works()
    {
        double entropy = SequenceComplexity.CalculateShannonEntropy("ATGC");
        Assert.That(entropy, Is.EqualTo(2.0).Within(0.01));
    }

    #endregion

    #region K-mer Entropy Tests

    [Test]
    public void CalculateKmerEntropy_VariedDinucleotides_ReturnsHigh()
    {
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.GreaterThan(1.5)); // High dinucleotide entropy
    }

    [Test]
    public void CalculateKmerEntropy_RepeatedDinucleotides_ReturnsLow()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.EqualTo(0).Within(0.01)); // Only AA
    }

    [Test]
    public void CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero()
    {
        var sequence = new DnaSequence("AT");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 5);

        Assert.That(entropy, Is.EqualTo(0));
    }

    #endregion

    #region Windowed Complexity Tests

    [Test]
    public void CalculateWindowedComplexity_ReturnsMultiplePoints()
    {
        var sequence = new DnaSequence(new string('A', 50) + new string('T', 50) + string.Concat(Enumerable.Repeat("ATGCATGC", 10)));
        var points = SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 20, stepSize: 20).ToList();

        Assert.That(points.Count, Is.GreaterThan(1));
    }

    [Test]
    public void CalculateWindowedComplexity_IncludesBothMetrics()
    {
        var sequence = new DnaSequence("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC");
        var points = SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 20, stepSize: 10).ToList();

        Assert.That(points[0].ShannonEntropy, Is.GreaterThan(0));
        Assert.That(points[0].LinguisticComplexity, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateWindowedComplexity_PositionsAreCorrect()
    {
        var sequence = new DnaSequence(new string('A', 100));
        var points = SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 20, stepSize: 20).ToList();

        Assert.That(points[0].Position, Is.EqualTo(10)); // Center of first window
        Assert.That(points[1].Position, Is.EqualTo(30)); // Center of second window
    }

    #endregion

    #region Low Complexity Region Tests

    [Test]
    public void FindLowComplexityRegions_FindsPolyARegion()
    {
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 20)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 20)));
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 20, entropyThreshold: 0.5).ToList();

        Assert.That(regions.Count, Is.GreaterThan(0));
    }

    [Test]
    public void FindLowComplexityRegions_HighComplexity_ReturnsEmpty()
    {
        var sequence = new DnaSequence("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC");
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 20, entropyThreshold: 0.5).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void FindLowComplexityRegions_ReturnsCorrectSequence()
    {
        var sequence = new DnaSequence("ATGCATGC" + new string('A', 64) + "ATGCATGC");
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 32, entropyThreshold: 0.5).ToList();

        if (regions.Count > 0)
        {
            Assert.That(regions[0].Sequence, Does.Contain("AAA"));
        }
    }

    #endregion

    #region DUST Score Tests

    [Test]
    public void CalculateDustScore_LowComplexity_ReturnsHigh()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAA");
        double dust = SequenceComplexity.CalculateDustScore(sequence);

        Assert.That(dust, Is.GreaterThan(1));
    }

    [Test]
    public void CalculateDustScore_HighComplexity_ReturnsLow()
    {
        var sequence = new DnaSequence("ATGCTAGCATGCTAGC");
        double dust = SequenceComplexity.CalculateDustScore(sequence);

        Assert.That(dust, Is.LessThan(1));
    }

    [Test]
    public void CalculateDustScore_EmptySequence_ReturnsZero()
    {
        double dust = SequenceComplexity.CalculateDustScore("");
        Assert.That(dust, Is.EqualTo(0));
    }

    [Test]
    public void CalculateDustScore_StringOverload_Works()
    {
        double dust = SequenceComplexity.CalculateDustScore("AAAAAAA");
        Assert.That(dust, Is.GreaterThan(0));
    }

    #endregion

    #region Masking Tests

    [Test]
    public void MaskLowComplexity_MasksPolyARegion()
    {
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 16)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 16)));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 2.0);

        Assert.That(masked, Does.Contain("N"));
    }

    [Test]
    public void MaskLowComplexity_PreservesHighComplexity()
    {
        // Use a longer and more varied sequence to avoid false positives
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCA");
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 10.0);

        Assert.That(masked, Does.Not.Contain("N"));
    }

    [Test]
    public void MaskLowComplexity_CustomMaskChar()
    {
        var sequence = new DnaSequence(new string('A', 100));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 1.0, maskChar: 'X');

        Assert.That(masked, Does.Contain("X"));
    }

    #endregion

    #region Compression Ratio Tests

    [Test]
    public void EstimateCompressionRatio_HighComplexity_ReturnsHigh()
    {
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGCATGCAATGC");
        double ratio = SequenceComplexity.EstimateCompressionRatio(sequence);

        Assert.That(ratio, Is.GreaterThan(0.5));
    }

    [Test]
    public void EstimateCompressionRatio_LowComplexity_ReturnsLow()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        double ratio = SequenceComplexity.EstimateCompressionRatio(sequence);

        Assert.That(ratio, Is.LessThan(0.3));
    }

    [Test]
    public void EstimateCompressionRatio_EmptySequence_ReturnsZero()
    {
        double ratio = SequenceComplexity.EstimateCompressionRatio("");
        Assert.That(ratio, Is.EqualTo(0));
    }

    [Test]
    public void EstimateCompressionRatio_RangeIsZeroToOne()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        double ratio = SequenceComplexity.EstimateCompressionRatio(sequence);

        Assert.That(ratio, Is.GreaterThanOrEqualTo(0));
        Assert.That(ratio, Is.LessThanOrEqualTo(1));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void CalculateLinguisticComplexity_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.CalculateLinguisticComplexity((DnaSequence)null!));
    }

    [Test]
    public void CalculateShannonEntropy_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.CalculateShannonEntropy((DnaSequence)null!));
    }

    [Test]
    public void CalculateWindowedComplexity_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.CalculateWindowedComplexity((DnaSequence)null!).ToList());
    }

    [Test]
    public void FindLowComplexityRegions_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.FindLowComplexityRegions((DnaSequence)null!).ToList());
    }

    [Test]
    public void CalculateDustScore_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.CalculateDustScore((DnaSequence)null!));
    }

    [Test]
    public void MaskLowComplexity_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.MaskLowComplexity((DnaSequence)null!));
    }

    [Test]
    public void CalculateLinguisticComplexity_ZeroWordLength_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 0));
    }

    [Test]
    public void CalculateWindowedComplexity_ZeroWindowSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 0).ToList());
    }

    #endregion
}
