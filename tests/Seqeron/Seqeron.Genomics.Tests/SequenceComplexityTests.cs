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
        // "ATGCTAGCATGCAATG" (N=16, maxWord=10): rich vocabulary
        // Hand-calculated: obs=91, max=103 => LC = 91/103
        // Source: Troyanskaya et al. (2002) summation formula
        var sequence = new DnaSequence("ATGCTAGCATGCAATG");
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence);

        Assert.That(lc, Is.EqualTo(91.0 / 103.0).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_LowComplexity_ReturnsLow()
    {
        // Homopolymer "AAAAAAAAAAAAAAAA" (N=16, maxWord=10):
        // Each word length i has exactly 1 observed word; V_max = min(4^i, N-i+1)
        // Hand-calculated: obs=10, max=103 => LC = 10/103
        // Source: Orlov & Potapov (2004)
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAA");
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence);

        Assert.That(lc, Is.EqualTo(10.0 / 103.0).Within(1e-10));
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
    public void CalculateLinguisticComplexity_SingleNucleotide_ReturnsOne()
    {
        // "A" (N=1): i=1: obs=1, V_max=min(4,1)=1 => LC = 1/1 = 1.0
        // Vocabulary is saturated: the only possible word is observed
        // Source: Troyanskaya et al. (2002) formula definition
        double lc = SequenceComplexity.CalculateLinguisticComplexity("A");

        Assert.That(lc, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateLinguisticComplexity_DinucleotideRepeat_LowerThanRandom()
    {
        // Repetitive dinucleotide pattern has reduced vocabulary
        // Source: Orlov & Potapov (2004) - repetitive patterns have lower complexity
        string repetitive = string.Concat(Enumerable.Repeat("AT", 20)); // ATATATATAT... (40bp)
        string varied = "ATGCTAGCATGCAATGCTAGCATGCAATGCTAGCAT"; // 36bp

        double lcRepetitive = SequenceComplexity.CalculateLinguisticComplexity(repetitive);
        double lcVaried = SequenceComplexity.CalculateLinguisticComplexity(varied);

        // Dinucleotide repeat has very limited vocabulary at all word lengths
        Assert.That(lcRepetitive, Is.LessThan(0.1), "Dinucleotide repeat should have very low complexity");
        Assert.That(lcVaried, Is.GreaterThan(0.4), "Varied sequence should have moderate-high complexity");
        Assert.That(lcRepetitive, Is.LessThan(lcVaried));
    }

    [Test]
    public void CalculateLinguisticComplexity_MaxWordLengthParameter_AffectsResult()
    {
        // maxWordLength parameter controls vocabulary depth
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGC");

        double lc1 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 1);
        double lc2 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 2);
        double lc5 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 5);
        double lc10 = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 10);

        // maxWordLength=1: only unigrams, obs=4/max=4 => LC=1.0
        Assert.That(lc1, Is.EqualTo(1.0));
        // Different maxWordLength values produce different results
        Assert.That(lc2, Is.Not.EqualTo(lc10).Within(1e-10));
        Assert.That(lc5, Is.Not.EqualTo(lc10).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_WikipediaExample_MatchesHandCalculation()
    {
        // Wikipedia "Linguistic sequence complexity" — example: ACGGGAAGCTGATTCCA (N=17)
        // Troyanskaya summation formula: LC = Σ observed / Σ possible
        // i=1: obs=4, max=min(4,17)=4;  i=2: obs=14, max=min(16,16)=16
        // i=3: obs=15, max=min(64,15)=15; i=4: obs=14, max=min(256,14)=14
        // Total: obs=47, max=49 => LC = 47/49
        // Wikipedia U-values: U1=4/4, U2=14/16, U3=15/15, U4=14/14
        const string wikiSequence = "ACGGGAAGCTGATTCCA";

        double lc = SequenceComplexity.CalculateLinguisticComplexity(wikiSequence, maxWordLength: 4);

        Assert.That(lc, Is.EqualTo(47.0 / 49.0).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_WikipediaDinucleotideRepeat_MatchesHandCalculation()
    {
        // Wikipedia: ACACACACACACACACA (N=17) — dinucleotide repeat
        // Wikipedia states: U1=2/4 (only A,C); U2=2/16 (only AC,CA)
        // Troyanskaya summation: obs=2 at every word length → obs=20, max=112
        // LC = 20/112 = 5/28
        const string dinucRepeat = "ACACACACACACACACA";

        double lc = SequenceComplexity.CalculateLinguisticComplexity(dinucRepeat, maxWordLength: 10);

        Assert.That(lc, Is.EqualTo(5.0 / 28.0).Within(1e-10));
    }

    [Test]
    public void CalculateLinguisticComplexity_MaximalComplexity_ReturnsOne()
    {
        // "ATGC" (N=4, maxWord=10): all positions have unique words at every length
        // i=1: obs=4/4; i=2: obs=3/3; i=3: obs=2/2; i=4: obs=1/1
        // Total: obs=10, max=10 => LC = 1.0 (maximum complexity)
        // Source: Troyanskaya et al. (2002) — saturated vocabulary
        double lc = SequenceComplexity.CalculateLinguisticComplexity("ATGC");

        Assert.That(lc, Is.EqualTo(1.0));
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
    public void CalculateKmerEntropy_VariedDinucleotides_ReturnsExact()
    {
        // "ATGCATGCATGCATGC" k=2: 15 dinucleotides, counts: AT=4, TG=4, GC=4, CA=3
        // H = -(3×(4/15)×log₂(4/15) + (3/15)×log₂(3/15))
        // Source: Shannon entropy formula applied to k-mer frequency distribution
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        double expected = -(3.0 * (4.0 / 15) * Math.Log2(4.0 / 15) + (3.0 / 15) * Math.Log2(3.0 / 15));
        Assert.That(entropy, Is.EqualTo(expected).Within(1e-10));
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
    public void CalculateKmerEntropy_RangeIsNonNegativeAndBounded_ForDnaSequences()
    {
        // K-mer entropy should be 0 ≤ H ≤ log2(4^k) for DNA sequences
        // Source: Shannon entropy maximum = log2(number of possible symbols)
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
                double maxEntropy = Math.Log2(Math.Pow(4, k));
                Assert.That(entropy, Is.GreaterThanOrEqualTo(0),
                    $"K-mer entropy < 0 for sequence: {seq}, k={k}");
                Assert.That(entropy, Is.LessThanOrEqualTo(maxEntropy),
                    $"K-mer entropy > log2(4^{k})={maxEntropy} for sequence: {seq}, k={k}");
            }
        });
    }

    [Test]
    public void CalculateKmerEntropy_UniformDinucleotides_ReturnsLog2Of3()
    {
        // "ATCG" k=2: 3 unique dinucleotides (AT, TC, CG), each appearing once
        // H = log2(3) ≈ 1.585 (maximum entropy for 3 symbols)
        // Source: Shannon entropy for uniform distribution
        var sequence = new DnaSequence("ATCG");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.EqualTo(Math.Log2(3)).Within(1e-10));
    }

    #endregion

    #region Windowed Complexity Tests

    [Test]
    public void CalculateWindowedComplexity_ReturnsCorrectPointCount()
    {
        // 50A + 50T + 80bp = 180 total; floor((180-20)/20)+1 = 9 windows
        var sequence = new DnaSequence(new string('A', 50) + new string('T', 50) + string.Concat(Enumerable.Repeat("ATGCATGC", 10)));
        var points = SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 20, stepSize: 20).ToList();

        Assert.That(points.Count, Is.EqualTo(9));
    }

    [Test]
    public void CalculateWindowedComplexity_IncludesBothMetrics_ExactValues()
    {
        // 68bp of repeated ATGC: each 20bp window has equal base distribution
        // Shannon entropy of 20bp "ATGCATGCATGCATGCATGC" = 2.0 (4 bases equally distributed)
        // LC computed with maxWordLength=min(6,20)=6 for each window
        var sequence = new DnaSequence("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC");
        var points = SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 20, stepSize: 10).ToList();

        Assert.That(points[0].ShannonEntropy, Is.EqualTo(2.0).Within(1e-10));
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
        // 80bp (ATGC×20) + 64A + 80bp (ATGC×20) = 224bp total
        // The poly-A stretch has entropy=0, well below threshold=0.5
        // Exactly 1 low-complexity region should be detected, starting near the poly-A
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 20)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 20)));
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 20, entropyThreshold: 0.5).ToList();

        Assert.That(regions.Count, Is.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(79));
        Assert.That(regions[0].End, Is.EqualTo(146));
        Assert.That(regions[0].MinEntropy, Is.EqualTo(0));
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
        // "ATGCATGC" (8bp) + 64A + "ATGCATGC" (8bp) = 80bp total
        // With window=32, threshold=0.5: region starts at pos 6, ends at 75, length=70
        // MinEntropy=0 (pure homopolymer windows)
        var sequence = new DnaSequence("ATGCATGC" + new string('A', 64) + "ATGCATGC");
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 32, entropyThreshold: 0.5).ToList();

        Assert.That(regions.Count, Is.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(6));
        Assert.That(regions[0].End, Is.EqualTo(75));
        Assert.That(regions[0].Length, Is.EqualTo(70));
        Assert.That(regions[0].MinEntropy, Is.EqualTo(0));
    }

    #endregion

    #region DUST Score Tests

    [Test]
    public void CalculateDustScore_LowComplexity_ReturnsHigh()
    {
        // "AAAAAAAAAAAAAAAAAA" (N=18): 16 AAA triplets
        // score = 16×15/2 = 120, DUST = 120/(16-1) = 8.0
        // Source: Morgulis et al. (2006) symmetric DUST formula
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAA");
        double dust = SequenceComplexity.CalculateDustScore(sequence);

        Assert.That(dust, Is.EqualTo(8.0).Within(1e-10));
    }

    [Test]
    public void CalculateDustScore_HighComplexity_ReturnsLow()
    {
        // "ATGCTAGCATGCTAGC" (N=16): diverse triplets, low pair counts
        // Hand-calculated: score=6/(14-1) = 6/13
        // Source: Morgulis et al. (2006)
        var sequence = new DnaSequence("ATGCTAGCATGCTAGC");
        double dust = SequenceComplexity.CalculateDustScore(sequence);

        Assert.That(dust, Is.EqualTo(6.0 / 13.0).Within(1e-10));
    }

    [Test]
    public void CalculateDustScore_EmptySequence_ReturnsZero()
    {
        double dust = SequenceComplexity.CalculateDustScore("");
        Assert.That(dust, Is.EqualTo(0));
    }

    [Test]
    public void CalculateDustScore_SequenceShorterThanWordSize_ReturnsZero()
    {
        // "AT" (N=2) < wordSize=3: no triplets extractable → 0
        double dust = SequenceComplexity.CalculateDustScore("AT");
        Assert.That(dust, Is.EqualTo(0));
    }

    [Test]
    public void CalculateDustScore_StringOverload_ReturnsExact()
    {
        // "AAAAAAA" (N=7): 5 triplets, all AAA
        // score = 5×4/2 = 10, DUST = 10/(5-1) = 2.5
        // Source: Morgulis et al. (2006) symmetric formula
        double dust = SequenceComplexity.CalculateDustScore("AAAAAAA");
        Assert.That(dust, Is.EqualTo(2.5).Within(1e-10));
    }

    #endregion

    #region Masking Tests

    [Test]
    public void MaskLowComplexity_MasksLowComplexityWindows()
    {
        // ATGC×16 (64bp) + A×64 + ATGC×16 (64bp) = 192bp total, window=64, threshold=2.0
        // ATGC×16 window DUST ≈ 7.4 (4 recurring triplets), A×64 DUST = 31.0
        // All windows exceed threshold=2.0, so entire sequence is masked
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 16)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 16)));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 2.0);

        Assert.That(masked.Length, Is.EqualTo(192));
        Assert.That(masked.Count(c => c == 'N'), Is.EqualTo(192));
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
        // 100A, window=64, threshold=1.0: DUST(A×64) = 31.0 >> 1.0
        // All positions covered by at least one window are masked with 'X'
        var sequence = new DnaSequence(new string('A', 100));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 1.0, maskChar: 'X');

        Assert.That(masked.Length, Is.EqualTo(100));
        Assert.That(masked.Count(c => c == 'X'), Is.EqualTo(100));
    }

    #endregion

    #region Compression Ratio Tests

    [Test]
    public void EstimateCompressionRatio_HighComplexity_ReturnsExact()
    {
        // "ATGCTAGCATGCAATGCTAGCATGCAATGC" (N=30): diverse vocabulary
        // unique substrings = 112, expected = 216 → ratio = 112/216 = 14/27
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGCATGCAATGC");
        double ratio = SequenceComplexity.EstimateCompressionRatio(sequence);

        Assert.That(ratio, Is.EqualTo(14.0 / 27.0).Within(1e-10));
    }

    [Test]
    public void EstimateCompressionRatio_LowComplexity_ReturnsExact()
    {
        // 31×A: only 1 unique substring per length → unique=10, expected=224
        // ratio = 10/224 = 5/112
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        double ratio = SequenceComplexity.EstimateCompressionRatio(sequence);

        Assert.That(ratio, Is.EqualTo(5.0 / 112.0).Within(1e-10));
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
    public void EstimateCompressionRatio_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceComplexity.EstimateCompressionRatio((DnaSequence)null!));
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

    [Test]
    public void CalculateWindowedComplexity_ZeroStepSize_ThrowsException()
    {
        var sequence = new DnaSequence("ATGCATGCATGCATGCATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceComplexity.CalculateWindowedComplexity(sequence, windowSize: 10, stepSize: 0).ToList());
    }

    [Test]
    public void MaskLowComplexity_ShortSequence_PreservesOriginal()
    {
        // When sequence length < windowSize, no windows are processed → original returned
        var sequence = new DnaSequence("ATGC");
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 0.0);

        Assert.That(masked, Is.EqualTo("ATGC"));
    }

    #endregion
}
