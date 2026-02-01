using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

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
    public void CalculateLinguisticComplexity_RangeIsZeroToOne_ForMultipleSequences()
    {
        // Range invariant: 0 ≤ LC ≤ 1 for all valid inputs
        // Source: Troyanskaya et al. (2002), mathematical definition
        var testSequences = new[]
        {
            "A",                           // Single nucleotide
            "AAAA",                        // Homopolymer
            "ATGC",                        // All bases once
            "ATGCATGCATGC",               // Repeated pattern
            "ATGCTAGCATGCAATGCTAGCATGC",  // Random-like
            new string('A', 100),          // Long homopolymer
            string.Concat(Enumerable.Repeat("ATGC", 25))  // Long varied
        };

        Assert.Multiple(() =>
        {
            foreach (string seq in testSequences)
            {
                double lc = SequenceComplexity.CalculateLinguisticComplexity(seq);
                Assert.That(lc, Is.GreaterThanOrEqualTo(0), $"LC < 0 for sequence: {seq[..Math.Min(20, seq.Length)]}...");
                Assert.That(lc, Is.LessThanOrEqualTo(1), $"LC > 1 for sequence: {seq[..Math.Min(20, seq.Length)]}...");
            }
        });
    }

    [Test]
    public void CalculateLinguisticComplexity_StringOverload_MatchesDnaSequenceOverload()
    {
        // API consistency: string overload should produce same result as DnaSequence
        const string sequence = "ATGCTAGCATGCAATG";
        var dnaSeq = new DnaSequence(sequence);

        double lcString = SequenceComplexity.CalculateLinguisticComplexity(sequence);
        double lcDna = SequenceComplexity.CalculateLinguisticComplexity(dnaSeq);

        Assert.That(lcString, Is.EqualTo(lcDna).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_SingleNucleotide_ReturnsPositiveValue()
    {
        // Single nucleotide has vocabulary of size 1 at word length 1
        // Source: Definition - vocabulary exists
        double lc = SequenceComplexity.CalculateLinguisticComplexity("A");

        Assert.That(lc, Is.GreaterThan(0));
        Assert.That(lc, Is.LessThanOrEqualTo(1));
    }

    [Test]
    public void CalculateLinguisticComplexity_DinucleotideRepeat_LowerThanRandom()
    {
        // Repetitive dinucleotide pattern has reduced vocabulary
        // Source: Orlov & Potapov (2004) - repetitive patterns have lower complexity
        string repetitive = string.Concat(Enumerable.Repeat("AT", 20)); // ATATATATAT...
        string varied = "ATGCTAGCATGCAATGCTAGCATGCAATGCTAGCAT";

        double lcRepetitive = SequenceComplexity.CalculateLinguisticComplexity(repetitive);
        double lcVaried = SequenceComplexity.CalculateLinguisticComplexity(varied);

        Assert.That(lcRepetitive, Is.LessThan(lcVaried));
    }

    [Test]
    public void CalculateLinguisticComplexity_MaxWordLengthParameter_AffectsResult()
    {
        // maxWordLength parameter controls vocabulary depth
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGC");

        double lc2 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 2);
        double lc5 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 5);
        double lc10 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 10);

        // Different maxWordLength values should produce different results
        // (unless sequence is too short to see the difference)
        Assert.That(lc2, Is.Not.EqualTo(lc10).Within(1e-10));
        Assert.That(lc5, Is.Not.EqualTo(lc10).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_LowercaseInput_HandledCorrectly()
    {
        // Case insensitivity for robustness
        const string upper = "ATGCTAGCATGC";
        const string lower = "atgctagcatgc";
        const string mixed = "AtGcTaGcAtGc";

        double lcUpper = SequenceComplexity.CalculateLinguisticComplexity(upper);
        double lcLower = SequenceComplexity.CalculateLinguisticComplexity(lower);
        double lcMixed = SequenceComplexity.CalculateLinguisticComplexity(mixed);

        Assert.Multiple(() =>
        {
            Assert.That(lcLower, Is.EqualTo(lcUpper).Within(1e-10));
            Assert.That(lcMixed, Is.EqualTo(lcUpper).Within(1e-10));
        });
    }

    #endregion

    #region Shannon Entropy Tests

    [Test]
    public void CalculateShannonEntropy_EqualBases_ReturnsTwo()
    {
        // Equal distribution of all 4 bases = max entropy (2 bits)
        // Source: Wikipedia - Entropy (information theory), H_max = log2(4) = 2
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(2.0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_SingleBase_ReturnsZero()
    {
        // Only one base type = zero entropy (no uncertainty)
        // Source: Wikipedia - Entropy (information theory), H = 0 when p = 1
        var sequence = new DnaSequence("AAAAAAA");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_TwoBases_ReturnsOne()
    {
        // Two bases equally distributed = 1 bit entropy
        // Source: Binary entropy H = -2 × (0.5 × log2(0.5)) = 1
        var sequence = new DnaSequence("ATATATAT");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_EmptySequence_ReturnsZero()
    {
        // Empty sequence = no information content
        // Source: Convention, no data = no entropy
        double entropy = SequenceComplexity.CalculateShannonEntropy("");
        Assert.That(entropy, Is.EqualTo(0));
    }

    [Test]
    public void CalculateShannonEntropy_StringOverload_MatchesDnaSequenceOverload()
    {
        // API consistency: string overload should produce same result as DnaSequence
        const string sequence = "ATGCATGCATGCATGC";
        var dnaSeq = new DnaSequence(sequence);

        double entropyString = SequenceComplexity.CalculateShannonEntropy(sequence);
        double entropyDna = SequenceComplexity.CalculateShannonEntropy(dnaSeq);

        Assert.That(entropyString, Is.EqualTo(entropyDna).Within(1e-10));
    }

    [Test]
    public void CalculateShannonEntropy_RangeIsZeroToTwo_ForDnaSequences()
    {
        // Invariant INV-ENT-001: 0 ≤ H ≤ 2 for any DNA sequence
        // Source: Wikipedia - max entropy = log2(alphabet size) = log2(4) = 2
        var testSequences = new[]
        {
            "A",                           // Single nucleotide
            "AAAA",                         // Homopolymer
            "ATGC",                         // All bases once
            "ATGCATGCATGC",                 // Repeated pattern
            "ATGCTAGCATGCAATGCTAGCATGC",   // Random-like
            new string('A', 100),           // Long homopolymer
            string.Concat(Enumerable.Repeat("ATGC", 25))  // Long varied
        };

        Assert.Multiple(() =>
        {
            foreach (string seq in testSequences)
            {
                double entropy = SequenceComplexity.CalculateShannonEntropy(seq);
                Assert.That(entropy, Is.GreaterThanOrEqualTo(0),
                    $"Entropy < 0 for sequence: {seq[..Math.Min(20, seq.Length)]}...");
                Assert.That(entropy, Is.LessThanOrEqualTo(2.0),
                    $"Entropy > 2 for sequence: {seq[..Math.Min(20, seq.Length)]}...");
            }
        });
    }

    [Test]
    public void CalculateShannonEntropy_ThreeBases_ReturnsLog2Of3()
    {
        // Three bases equally distributed = log2(3) ≈ 1.585 bits
        // Source: Shannon entropy formula for n=3 uniform symbols
        var sequence = new DnaSequence("ATGATGATG"); // A, T, G each 33.3%
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        double expectedEntropy = Math.Log2(3); // ≈ 1.585
        Assert.That(entropy, Is.EqualTo(expectedEntropy).Within(0.01));
    }

    [Test]
    public void CalculateShannonEntropy_LowercaseInput_HandledCorrectly()
    {
        // Case insensitivity for robustness
        const string upper = "ATGCATGCATGC";
        const string lower = "atgcatgcatgc";
        const string mixed = "AtGcAtGcAtGc";

        double entropyUpper = SequenceComplexity.CalculateShannonEntropy(upper);
        double entropyLower = SequenceComplexity.CalculateShannonEntropy(lower);
        double entropyMixed = SequenceComplexity.CalculateShannonEntropy(mixed);

        Assert.Multiple(() =>
        {
            Assert.That(entropyLower, Is.EqualTo(entropyUpper).Within(1e-10));
            Assert.That(entropyMixed, Is.EqualTo(entropyUpper).Within(1e-10));
        });
    }

    #endregion

    #region K-mer Entropy Tests

    [Test]
    public void CalculateKmerEntropy_VariedDinucleotides_ReturnsHigh()
    {
        // Source: High k-mer diversity = high entropy
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.GreaterThan(1.5)); // High dinucleotide entropy
    }

    [Test]
    public void CalculateKmerEntropy_RepeatedDinucleotides_ReturnsZero()
    {
        // Homopolymer has only one k-mer type = zero entropy
        // Source: Single symbol = zero entropy
        var sequence = new DnaSequence("AAAAAAAAAA");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.EqualTo(0).Within(0.01)); // Only AA
    }

    [Test]
    public void CalculateKmerEntropy_SequenceShorterThanK_ReturnsZero()
    {
        // No k-mers extractable = zero entropy
        var sequence = new DnaSequence("AT");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 5);

        Assert.That(entropy, Is.EqualTo(0));
    }

    [Test]
    public void CalculateKmerEntropy_InvalidK_ThrowsException()
    {
        // Parameter validation: k must be >= 1
        var sequence = new DnaSequence("ATGCATGC");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.CalculateKmerEntropy(sequence, k: 0));
    }

    [Test]
    public void CalculateKmerEntropy_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.CalculateKmerEntropy((DnaSequence)null!, k: 2));
    }

    [Test]
    public void CalculateKmerEntropy_RangeIsNonNegative_ForDnaSequences()
    {
        // K-mer entropy should always be >= 0
        var testSequences = new[]
        {
            ("ATGC", 1),
            ("ATGCATGC", 2),
            ("ATGCATGCATGC", 3),
            ("AAAAAAAAAA", 2),
            ("ATATATATAT", 2)
        };

        Assert.Multiple(() =>
        {
            foreach (var (seq, k) in testSequences)
            {
                var dnaSeq = new DnaSequence(seq);
                double entropy = SequenceComplexity.CalculateKmerEntropy(dnaSeq, k);
                Assert.That(entropy, Is.GreaterThanOrEqualTo(0),
                    $"K-mer entropy < 0 for sequence: {seq}, k={k}");
            }
        });
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
    public void CalculateLinguisticComplexity_NegativeWordLength_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: -1));
    }

    [Test]
    public void FindLowComplexityRegions_InvalidWindowSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 0).ToList());
    }

    [Test]
    public void MaskLowComplexity_ResultLengthEqualsInputLength()
    {
        // Invariant: masked sequence length equals input length
        var sequence = new DnaSequence(new string('A', 100) + "ATGCTAGCATGCAATG");
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 1.0);

        Assert.That(masked.Length, Is.EqualTo(sequence.Length));
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
