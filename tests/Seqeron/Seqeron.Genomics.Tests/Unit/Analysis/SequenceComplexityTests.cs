namespace Seqeron.Genomics.Tests.Unit.Analysis;

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

    [TestCase(9)]
    [TestCase(int.MaxValue)]
    public void CalculateLinguisticComplexity_RosalindLingSample_AllWordLengths_Returns0875(int maxWordLength)
    {
        // Rosalind LING sample: lc("ATTTGGATT") = 0.875 = sub(s)/m(4,9) = 35/40,
        // i.e. Troyanskaya et al. (2002) LC = A(s)/M(s) summed over ALL lengths 1..N.
        // With maxWordLength ≥ N the Orlov & Potapov (2004) truncated sum equals it.
        // Reference: 2026-09 review Python recomputation (fractions) = 7/8.
        double lc = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence("ATTTGGATT"), maxWordLength);

        Assert.That(lc, Is.EqualTo(0.875).Within(1e-12));
    }

    [TestCase("ATGCTAGCATGCAATG", 28.0 / 31.0)]
    [TestCase("AAAAAAAAAAAAAAAA", 4.0 / 31.0)]
    [TestCase("ACACACACACACACACA", 33.0 / 140.0)]
    [TestCase("ACGGGAAGCTGATTCCA", 69.0 / 70.0)]
    public void CalculateLinguisticComplexity_TroyanskayaFullLength_MatchesReference(string sequence, double expected)
    {
        // Troyanskaya et al. (2002): LC = Σ_{l=1..N} A_l / Σ_{l=1..N} min(4^l, N−l+1).
        // Expected values from the 2026-09 review Python reference (exact fractions).
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: sequence.Length);

        Assert.That(lc, Is.EqualTo(expected).Within(1e-12));
    }

    [TestCase("AAACCCGGGTTT", 51.0 / 55.0)]
    [TestCase("AACCGGTTACGT", 52.0 / 55.0)]
    [TestCase("ACGTACGTACGT", 28.0 / 55.0)]
    [TestCase("AAAACCCCGGGG", 9.0 / 11.0)]
    [TestCase("AAAAAACCCCCC", 3.0 / 5.0)]
    [TestCase("AAAAAAAAAACC", 4.0 / 11.0)]
    public void CalculateLinguisticComplexity_OrlovSumForm_UniversalmotifSequences_MatchesReference(string sequence, double expected)
    {
        // Orlov & Potapov (2004) CL = Σ_{i=1..m} V_i / Σ_{i=1..m} min(4^i, N−i+1), m = 7.
        // The per-length (V_i, V_max,i) used by the Python reference reproduce the product form
        // of universalmotif::calc_complexity(method = "Trifonov", max word size 7) to 4 dp
        // (0.6364, 0.7273, 0.01231, 0.2386, 0.0227, 0.0011), independently confirming V_max,i.
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength: 7);

        Assert.That(lc, Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void CalculateLinguisticComplexity_WordLengthsBeyond4Pow31_NoOverflow()
    {
        // m ≥ 32 makes 4^i exceed long.MaxValue; V_max,i must still be N−i+1.
        // Python reference (exact): m = N = 40 → 749/761.
        const string sequence = "ACGTTGCAAGGCTTACCGATGCATCGGATCCTAGGCTAAC";

        double lc = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(sequence), maxWordLength: 40);

        Assert.That(lc, Is.EqualTo(749.0 / 761.0).Within(1e-12));
    }

    [TestCase(13, 782.0 / 1209.0)]
    [TestCase(20, 1405.0 / 1937.0)]
    [TestCase(50, 1971.0 / 2251.0)]
    [TestCase(120, 6427.0 / 6987.0)]
    public void CalculateLinguisticComplexity_SuffixTreePath_RepeatRichSequence_MatchesReference(int maxWordLength, double expected)
    {
        // m > 12 counts V_i from the suffix tree (Troyanskaya et al. 2002). Sequence contains an
        // exact 20-nt repeat, a (CAG)10 microsatellite and an A/C-only tail so internal nodes,
        // multi-length edges and leaf edges are all exercised. Expected: Python reference (exact).
        const string sequence =
            "GCTAAAGACAATTACATAACATACACGTCACAGCAGCAGCAGCAGCAGCAGCAGCAGCAGGCTAAAGACAATTACATAACC" +
            "AAACAAAACCCCCCCAAAACCCCCAACACACCAACCCCC";

        double lcDna = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(sequence), maxWordLength);
        double lcString = SequenceComplexity.CalculateLinguisticComplexity(sequence.ToLowerInvariant(), maxWordLength);

        Assert.Multiple(() =>
        {
            Assert.That(lcDna, Is.EqualTo(expected).Within(1e-12));
            Assert.That(lcString, Is.EqualTo(expected).Within(1e-12));
        });
    }

    [Test]
    public void CalculateLinguisticComplexity_FullLengthLongHomopolymer_ExactAndLinearTime()
    {
        // Troyanskaya all-length LC of A^N: V_i = 1 for every i, so LC = N / Σ_i min(4^i, N−i+1).
        // N = 200,000 would need ~2·10^10 substring characters by direct enumeration; the suffix-tree
        // path is linear. Expected denominator computed in closed form below.
        const int n = 200_000;
        long possible = 0;
        for (int i = 1; i <= n; i++)
            possible += i < 16 ? Math.Min(1L << (2 * i), n - i + 1) : n - i + 1;

        double lc = SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(new string('A', n)), int.MaxValue);

        Assert.That(lc, Is.EqualTo((double)n / possible).Within(1e-15));
    }

    // Alphabet size a of M = Σ min(a^i, N − i + 1) (Troyanskaya et al. 2002; Rosalind LING "alphabet of size a"):
    // {A,C,G,T/U} extended by any other symbol present. Expected values: Python brute force (exact Fractions,
    // scratch lc_ref.py) — 3000 random cases over ACGT/ACGTN/ACGU/ACGTU/IUPAC alphabets, m 1..70, 0 mismatches.
    // Before the fix ACGTN gave 15/14 > 1 (V_max assumed 4^i).
    [TestCase("ACGTN", 5, 1.0)]                                      // 15/15
    [TestCase("acgtn", 10, 1.0)]                                     // upper-cased, m clamped to N
    [TestCase("ATGCATGCNN", 10, 22.0 / 25.0)]                        // 44/50, a = 5
    [TestCase("NNNNNNNN", 8, 8.0 / 33.0)]                            // a = 5 ({N} ∪ ACGT)
    [TestCase("ACGTNNNNACGTNNNNACGTRYACGTNNNN", 6, 69.0 / 142.0)]    // hash path, a = 7
    [TestCase("ACGTNNNNACGTNNNNACGTRYACGTNNNN", 30, 345.0 / 442.0)]  // suffix-tree path (m > 12), a = 7
    [TestCase("ACGUACGUAAUU", 12, 6.0 / 7.0)]                        // RNA: a = 4 (U replaces T)
    [TestCase("ACGTUACGTU", 10, 0.8)]                                // T and U both present: a = 5
    public void CalculateLinguisticComplexity_NonAcgtSymbols_AlphabetExtended_MatchesBruteForce(
        string sequence, int maxWordLength, double expected)
    {
        double lc = SequenceComplexity.CalculateLinguisticComplexity(sequence, maxWordLength);

        Assert.That(lc, Is.EqualTo(expected).Within(1e-15));
    }

    [Test]
    public void CalculateLinguisticComplexity_RnaEqualsDnaCounterpart()
    {
        // U plays the role of T (a = 4), so the RNA and DNA spellings have identical LC.
        Assert.That(SequenceComplexity.CalculateLinguisticComplexity("ACGUACGUAAUU", 12),
            Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity("ACGTACGTAATT", 12)));
    }

    [Test]
    public void CalculateLinguisticComplexity_LargeAlphabet_NeverExceedsOne_NoOverflow()
    {
        // 300 distinct caseless CJK symbols: every substring is distinct, so LC = 1 exactly for the hash and
        // suffix-tree paths; a^i with a = 304 would overflow long by i = 8 without the saturation guard.
        string s = new(Enumerable.Range(0x4E00, 300).Select(c => (char)c).ToArray());

        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(s, 10), Is.EqualTo(1.0));
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(s, 300), Is.EqualTo(1.0));
        });
    }

    #endregion

    #region Shannon Entropy Tests

    [Test]
    public void CalculateShannonEntropy_EqualBases_ReturnsTwo()
    {
        // Equal distribution of all 4 bases = max entropy (2 bits)
        // Source: Wikipedia - Entropy (information theory), H_max = log₂(4) = 2
        // Hand calc: p=0.25 exact in float, log₂(0.25)=−2 exact → H = 2.0 exact
        var sequence = new DnaSequence("ATGCATGCATGCATGC");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(2.0));
    }

    [Test]
    public void CalculateShannonEntropy_SingleBase_ReturnsZero()
    {
        // Only one base type = zero entropy (no uncertainty)
        // Source: Wikipedia - Entropy (information theory), H = 0 when p = 1
        // Hand calc: p=1.0, log₂(1.0)=0 → H = 0.0 exact
        var sequence = new DnaSequence("AAAAAAA");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateShannonEntropy_TwoBases_ReturnsOne()
    {
        // Two bases equally distributed = 1 bit entropy
        // Source: Binary entropy H = −2 × (0.5 × log₂(0.5)) = 1
        // Hand calc: p=0.5 exact in float, log₂(0.5)=−1 exact → H = 1.0 exact
        var sequence = new DnaSequence("ATATATAT");
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        Assert.That(entropy, Is.EqualTo(1.0));
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
        // Both paths feed same uppercase string to same CalculateShannonEntropyCore → bitwise identical
        const string sequence = "ATGCATGCATGCATGC";
        var dnaSeq = new DnaSequence(sequence);

        double entropyString = SequenceComplexity.CalculateShannonEntropy(sequence);
        double entropyDna = SequenceComplexity.CalculateShannonEntropy(dnaSeq);

        Assert.That(entropyString, Is.EqualTo(entropyDna));
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
        // Three bases equally distributed = log₂(3) ≈ 1.58496 bits
        // Source: Shannon entropy formula for n=3 uniform symbols
        // Hand calc: A=T=G=3/9 → H = −3×(⅓×log₂⅓) = log₂(3)
        var sequence = new DnaSequence("ATGATGATG"); // A, T, G each 33.3%
        double entropy = SequenceComplexity.CalculateShannonEntropy(sequence);

        double expectedEntropy = Math.Log2(3); // ≈ 1.58496
        Assert.That(entropy, Is.EqualTo(expectedEntropy).Within(1e-10));
    }

    [Test]
    public void CalculateShannonEntropy_LowercaseInput_HandledCorrectly()
    {
        // Case insensitivity: ToUpperInvariant() produces same string → bitwise identical result
        const string upper = "ATGCATGCATGC";
        const string lower = "atgcatgcatgc";
        const string mixed = "AtGcAtGcAtGc";

        double entropyUpper = SequenceComplexity.CalculateShannonEntropy(upper);
        double entropyLower = SequenceComplexity.CalculateShannonEntropy(lower);
        double entropyMixed = SequenceComplexity.CalculateShannonEntropy(mixed);

        Assert.Multiple(() =>
        {
            Assert.That(entropyLower, Is.EqualTo(entropyUpper));
            Assert.That(entropyMixed, Is.EqualTo(entropyUpper));
        });
    }

    [Test]
    public void CalculateShannonEntropy_NonDnaCharacters_Ignored()
    {
        // Alphabet is {A,T,G,C}: non-ATGC chars excluded from both numerator and denominator
        // Source: Shannon entropy computed over fixed DNA alphabet
        Assert.Multiple(() =>
        {
            // "ATGCNN" → only ATGC counted (4 bases, p=0.25 each) → H = 2.0 exact
            double entropy = SequenceComplexity.CalculateShannonEntropy("ATGCNN");
            Assert.That(entropy, Is.EqualTo(2.0));

            // All non-ATGC → total=0 → H = 0.0
            double entropyAllN = SequenceComplexity.CalculateShannonEntropy("NNNNNN");
            Assert.That(entropyAllN, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void CalculateShannonEntropy_RnaUracil_CountedAsFourthNucleotide()
    {
        // RNA U is the same nucleotide class as DNA T (IUPAC-IUB 1970), alphabet {A,C,G,T/U}.
        // Reference: scipy.stats.entropy(counts, base=2) with U→T:
        //   "ACGU" → [1,1,1,1] → 2.0; "AAUU" → [2,0,0,2] → 1.0; "acgu" → 2.0;
        //   "ACGUN" → 2.0 (N excluded); "GGGGCCCAU" → [1,3,4,1] → 1.7527152789797047.
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateShannonEntropy("ACGU"), Is.EqualTo(2.0));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("AAUU"), Is.EqualTo(1.0));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("acgu"), Is.EqualTo(2.0));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("ACGUN"), Is.EqualTo(2.0));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("GGGGCCCAU"),
                Is.EqualTo(1.7527152789797047).Within(1e-12));
            // RNA and its DNA equivalent give the same entropy.
            Assert.That(SequenceComplexity.CalculateShannonEntropy("GGGGCCCAU"),
                Is.EqualTo(SequenceComplexity.CalculateShannonEntropy("GGGGCCCAT")));
        });
    }

    [Test]
    public void CalculateShannonEntropy_MatchesScipyReference()
    {
        // Reference: scipy.stats.entropy([nA,nC,nG,nT], base=2), non-ACGT/U excluded.
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateShannonEntropy("AAAAAAAAAAAAACGT"),
                Is.EqualTo(0.9933927290103627).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("AACGTTTGCA"),
                Is.EqualTo(1.970950594454669).Within(1e-12));
            Assert.That(SequenceComplexity.CalculateShannonEntropy("ACGTNNRYacgt"), Is.EqualTo(2.0));
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
        // Homopolymer has only one k-mer type "AA" = zero entropy
        // Source: Shannon (1948) — single symbol: p=1, H = −1×log₂(1) = 0 exact
        var sequence = new DnaSequence("AAAAAAAAAA");
        double entropy = SequenceComplexity.CalculateKmerEntropy(sequence, k: 2);

        Assert.That(entropy, Is.EqualTo(0.0)); // Only AA
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

    // Windowed complexity tests moved to canonical file
    // SequenceComplexity_CalculateWindowedComplexity_Tests.cs (SEQ-COMPLEX-WINDOW-001).

    #region Low Complexity Region Tests

    [Test]
    public void FindLowComplexityRegions_FindsPolyARegion()
    {
        // 80bp (ATGC×20) + 64A + 80bp (ATGC×20) = 224bp total, w=20, threshold 0.5.
        // Flagged windows (H < 0.5): starts 79 ("C"+19A, H=0.286) .. 126 (19A+"T", H=0.286);
        // start 127 (18A+"TG", H=0.569) is not flagged. Region = union of flagged windows
        // = [79, 126+19] = 79..145 (BBDuk maskLowEntropy window-union rule; Python reference
        // scipy.stats.entropy(base=2) per window + bit-mask union gives (79, 145, 67, 0.0)).
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 20)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 20)));
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 20, entropyThreshold: 0.5).ToList();

        Assert.That(regions.Count, Is.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(79));
        Assert.That(regions[0].End, Is.EqualTo(145));
        Assert.That(regions[0].Length, Is.EqualTo(67));
        Assert.That(regions[0].MinEntropy, Is.EqualTo(0));
        Assert.That(regions[0].Sequence, Is.EqualTo("C" + new string('A', 65) + "T"));
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
        // "ATGCATGC" (8bp) + 64A + "ATGCATGC" (8bp) = 80bp total, w=32, threshold 0.5.
        // Flagged windows: starts 6..43 (last flagged window 43..74 = 29A+"ATG"); region =
        // union = 6..74, length 69 (Python reference: (6, 74, 69, 0.0)). MinEntropy=0 (poly-A windows).
        var sequence = new DnaSequence("ATGCATGC" + new string('A', 64) + "ATGCATGC");
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 32, entropyThreshold: 0.5).ToList();

        Assert.That(regions.Count, Is.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(6));
        Assert.That(regions[0].End, Is.EqualTo(74));
        Assert.That(regions[0].Length, Is.EqualTo(69));
        Assert.That(regions[0].MinEntropy, Is.EqualTo(0));
        Assert.That(regions[0].Sequence, Is.EqualTo("GC" + new string('A', 65) + "TG"));
    }

    [Test]
    public void FindLowComplexityRegions_InvalidArguments_ThrowEagerlyBeforeEnumeration()
    {
        // Validation must happen at call time, not on first MoveNext of the lazy iterator.
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentNullException>(() => SequenceComplexity.FindLowComplexityRegions((DnaSequence)null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 0));
    }

    [Test]
    public void FindLowComplexityRegions_OverlappingFlaggedWindows_MergeIntoOneRegion()
    {
        // w=8, threshold 0.6. Flagged windows (H < 0.6): starts 1,2,3,4 (end 11) and 7,8 (ends 14,15);
        // windows 5,6 have H = H(2/8,6/8) = 0.811 and are not flagged, but window 7 overlaps the
        // union 1..11, so the masked positions form ONE run 1..15 (BBDuk bit-mask union).
        // MinEntropy = H(1/8,7/8) = 0.5435644431995964. Python reference: (1, 15, 15, 0.5435644431995964).
        // (Pre-fix code emitted two overlapping regions 1..12 and 7..16.)
        var sequence = new DnaSequence("CAAAAACAAAAACAAACAAA");
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 8, entropyThreshold: 0.6).ToList();

        Assert.That(regions.Count, Is.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(1));
        Assert.That(regions[0].End, Is.EqualTo(15));
        Assert.That(regions[0].Length, Is.EqualTo(15));
        Assert.That(regions[0].MinEntropy, Is.EqualTo(0.5435644431995964).Within(1e-12));
        Assert.That(regions[0].Sequence, Is.EqualTo("AAAAACAAAAACAAA"));
    }

    [Test]
    public void FindLowComplexityRegions_TwoSeparatedTracts_ReturnsDisjointRegionsIncludingTrailing()
    {
        // A10 + (ACGT)×3 + C10 (32 bp), w=8, threshold 1.0: flagged windows 0..4 (union 0..11)
        // and 21..24 (union 21..31, trailing to sequence end). Python reference:
        // [(0, 11, 12, 0.0), (21, 31, 11, 0.0)].
        var sequence = new DnaSequence(new string('A', 10) + "ACGTACGTACGT" + new string('C', 10));
        var regions = SequenceComplexity.FindLowComplexityRegions(sequence, windowSize: 8, entropyThreshold: 1.0).ToList();

        Assert.That(regions.Select(r => (r.Start, r.End, r.Length)),
            Is.EqualTo(new[] { (0, 11, 12), (21, 31, 11) }));
        Assert.That(regions.Select(r => r.MinEntropy), Is.EqualTo(new[] { 0.0, 0.0 }));
        Assert.That(regions.Select(r => r.Sequence), Is.EqualTo(new[] { "AAAAAAAAAAAC", "TCCCCCCCCCC" }));
    }

    #endregion

    #region DUST Score Tests

    [Test]
    public void CalculateDustScore_LowComplexity_ReturnsHigh()
    {
        // "AAAAAAAAAAAAAAAAAA" (L=18): ℓ = 16 AAA triplets
        // Σ = 16×15/2 = 120, DUST = 120/(ℓ-1) = 120/15 = 8.0
        // Source: Morgulis et al. (2006); NCBI symdust / lh3/sdust normalise by ℓ-1
        var sequence = new DnaSequence("AAAAAAAAAAAAAAAAAA");
        double dust = SequenceComplexity.CalculateDustScore(sequence);

        Assert.That(dust, Is.EqualTo(8.0).Within(1e-10));
    }

    [Test]
    public void CalculateDustScore_HighComplexity_ReturnsLow()
    {
        // "ATGCTAGCATGCTAGC" (L=16): 14 triplets, ATG,TGC,GCT,CTA,TAG,AGC each ×2 (GCA,CAT ×1)
        // Σ = 6×(2·1/2) = 6, DUST = 6/(ℓ-1) = 6/13
        // Source: Morgulis et al. (2006); lh3/sdust (ℓ-1 normaliser)
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
        // "AAAAAAA" (L=7): 5 triplets, all AAA
        // Σ = 5×4/2 = 10, DUST = 10/(ℓ-1) = 10/4 = 2.5 (sdust -t 20 masks a 7-A run)
        // Source: Morgulis et al. (2006); lh3/sdust (ℓ-1 normaliser)
        double dust = SequenceComplexity.CalculateDustScore("AAAAAAA");
        Assert.That(dust, Is.EqualTo(2.5).Within(1e-10));
    }

    #endregion

    #region Masking Tests

    [Test]
    public void MaskLowComplexity_MasksLowComplexityWindows()
    {
        // ATGC×16 (64bp) + A×64 + ATGC×16 (64bp) = 192bp total, window=64, threshold=2.0
        // Reference: lh3/sdust -w 64 -t 20 reports the single interval [0,192) ⇒ all masked
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat("ATGC", 16)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 16)));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 2.0);

        Assert.That(masked.Length, Is.EqualTo(192));
        Assert.That(masked.Count(c => c == 'N'), Is.EqualTo(192));
    }

    [Test]
    public void MaskLowComplexity_PreservesHighComplexity()
    {
        // Reference: lh3/sdust -w 64 -t 100 reports no interval for this 78-bp sequence
        var sequence = new DnaSequence("ATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCAATGCTAGCATGCA");
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 10.0);

        Assert.That(masked, Does.Not.Contain("N"));
    }

    [Test]
    public void MaskLowComplexity_CustomMaskChar()
    {
        // 100A, window=64, threshold=1.0: reference lh3/sdust -w 64 -t 10 ⇒ [0,100) masked with 'X'
        var sequence = new DnaSequence(new string('A', 100));
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 1.0, maskChar: 'X');

        Assert.That(masked.Length, Is.EqualTo(100));
        Assert.That(masked.Count(c => c == 'X'), Is.EqualTo(100));
    }

    #endregion

    #region Compression Ratio Tests

    // Canonical Lempel–Ziv / compression-ratio tests live in
    // SequenceComplexity_EstimateCompressionRatio_Tests.cs (SEQ-COMPLEX-COMPRESS-001).
    // The earlier heuristic-output assertions were removed when the implementation
    // was corrected to the source-backed Lempel–Ziv (1976) complexity.

    [Test]
    public void EstimateCompressionRatio_EmptySequence_ReturnsZero()
    {
        double ratio = SequenceComplexity.EstimateCompressionRatio("");
        Assert.That(ratio, Is.EqualTo(0));
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
    public void MaskLowComplexity_ShortSequence_PreservesOriginal()
    {
        // ATGC: two distinct triplets, raw score 0 is not > 0.0 ⇒ nothing masked (sdust -t 0: no output)
        var sequence = new DnaSequence("ATGC");
        string masked = SequenceComplexity.MaskLowComplexity(sequence, windowSize: 64, threshold: 0.0);

        Assert.That(masked, Is.EqualTo("ATGC"));
    }

    #endregion
}
