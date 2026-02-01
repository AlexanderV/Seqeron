using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for KMER-FREQ-001: K-mer Frequency Analysis.
    /// 
    /// Canonical methods: GetKmerSpectrum, GetKmerFrequencies, CalculateKmerEntropy
    /// 
    /// Evidence:
    /// - Wikipedia: K-mer spectrum, frequency distribution as genomic signature
    /// - Wikipedia: Shannon entropy formula H = -Σ p log₂(p)
    /// - Shannon (1948): Maximum entropy when equiprobable, zero entropy for certainty
    /// - Rosalind KMER: K-mer composition problem
    /// </summary>
    [TestFixture]
    public class KmerAnalyzer_Frequency_Tests
    {
        private const double Tolerance = 0.0001;

        #region GetKmerFrequencies - Frequency Sum Invariant (M1)

        /// <summary>
        /// M1: Sum of all k-mer frequencies equals 1.0.
        /// Evidence: Mathematical definition of probability distribution.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_StandardSequence_SumsToOne()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("ACGTACGT", 2);
            double total = frequencies.Values.Sum();

            Assert.That(total, Is.EqualTo(1.0).Within(Tolerance),
                "Sum of all k-mer frequencies must equal 1.0");
        }

        /// <summary>
        /// M1: Frequency sum invariant holds for homopolymer.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_Homopolymer_SumsToOne()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("AAAAAAA", 3);
            double total = frequencies.Values.Sum();

            Assert.That(total, Is.EqualTo(1.0).Within(Tolerance));
        }

        /// <summary>
        /// M1: Frequency sum invariant holds for various k values.
        /// </summary>
        [Test]
        [TestCase("ACGTACGTACGT", 1)]
        [TestCase("ACGTACGTACGT", 2)]
        [TestCase("ACGTACGTACGT", 3)]
        [TestCase("ACGTACGTACGT", 4)]
        public void GetKmerFrequencies_VariousK_SumsToOne(string sequence, int k)
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies(sequence, k);

            if (frequencies.Count > 0)
            {
                double total = frequencies.Values.Sum();
                Assert.That(total, Is.EqualTo(1.0).Within(Tolerance),
                    $"Frequency sum must equal 1.0 for k={k}");
            }
        }

        #endregion

        #region GetKmerFrequencies - Edge Cases (M2)

        /// <summary>
        /// M2: Empty sequence returns empty dictionary.
        /// Evidence: Wikipedia pseudocode - no k-mers when L=0.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_EmptySequence_ReturnsEmptyDictionary()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("", 2);
            Assert.That(frequencies, Is.Empty);
        }

        /// <summary>
        /// M2: k > sequence length returns empty dictionary.
        /// Evidence: L - k + 1 ≤ 0 when k > L.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_KGreaterThanLength_ReturnsEmptyDictionary()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("ACG", 4);
            Assert.That(frequencies, Is.Empty);
        }

        #endregion

        #region GetKmerFrequencies - Calculation Correctness (M3)

        /// <summary>
        /// M3: Single k-mer type has frequency 1.0.
        /// Evidence: Mathematical definition - only k-mer gets 100% probability.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_SingleKmerType_HasFrequencyOne()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("AAA", 2);

            Assert.Multiple(() =>
            {
                Assert.That(frequencies, Has.Count.EqualTo(1));
                Assert.That(frequencies["AA"], Is.EqualTo(1.0).Within(Tolerance));
            });
        }

        /// <summary>
        /// S1: Individual frequency correctness for mixed sequence.
        /// </summary>
        [Test]
        public void GetKmerFrequencies_MixedSequence_CorrectFrequencies()
        {
            // AAACGT: 2-mers are AA, AA, AC, CG, GT (5 total)
            // AA appears 2 times = 2/5 = 0.4
            // AC, CG, GT each appear 1 time = 1/5 = 0.2 each
            var frequencies = KmerAnalyzer.GetKmerFrequencies("AAACGT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(frequencies["AA"], Is.EqualTo(0.4).Within(Tolerance));
                Assert.That(frequencies["AC"], Is.EqualTo(0.2).Within(Tolerance));
                Assert.That(frequencies["CG"], Is.EqualTo(0.2).Within(Tolerance));
                Assert.That(frequencies["GT"], Is.EqualTo(0.2).Within(Tolerance));
            });
        }

        #endregion

        #region GetKmerSpectrum - Spectrum Correctness (M4)

        /// <summary>
        /// M4: Spectrum correctly maps multiplicity to count.
        /// Evidence: Wikipedia K-mer - "k-mer spectrum shows multiplicity vs count"
        /// 
        /// For "ACGTACGT" with k=4:
        /// - 4-mers: ACGT(pos 0), CGTA(pos 1), GTAC(pos 2), TACG(pos 3), ACGT(pos 4)
        /// - ACGT appears 2 times, CGTA/GTAC/TACG each appear 1 time
        /// - Spectrum: {1: 3, 2: 1}
        /// </summary>
        [Test]
        public void GetKmerSpectrum_StandardSequence_ReturnsCorrectSpectrum()
        {
            var spectrum = KmerAnalyzer.GetKmerSpectrum("ACGTACGT", 4);

            Assert.Multiple(() =>
            {
                Assert.That(spectrum[1], Is.EqualTo(3),
                    "3 k-mers should appear exactly once");
                Assert.That(spectrum[2], Is.EqualTo(1),
                    "1 k-mer should appear exactly twice");
            });
        }

        /// <summary>
        /// M4: Homopolymer has single k-mer appearing multiple times.
        /// </summary>
        [Test]
        public void GetKmerSpectrum_Homopolymer_SingleEntryWithHighCount()
        {
            // "AAAAA" with k=2: AA appears 4 times (positions 0,1,2,3)
            var spectrum = KmerAnalyzer.GetKmerSpectrum("AAAAA", 2);

            Assert.Multiple(() =>
            {
                Assert.That(spectrum, Has.Count.EqualTo(1));
                Assert.That(spectrum[4], Is.EqualTo(1),
                    "One k-mer (AA) should appear 4 times");
            });
        }

        /// <summary>
        /// S2: Spectrum with multiple multiplicities.
        /// </summary>
        [Test]
        public void GetKmerSpectrum_MultipleMultiplicities_AllCaptured()
        {
            // Create sequence with known k-mer frequencies
            // "AABCABC" doesn't work for DNA, use "AACGAACG" with k=2
            // AA(2), AC(2), CG(2), GA(1), AA(counted above)
            // Actually: AA, AC, CG, GA, AA, AC, CG = AA:2, AC:2, CG:2, GA:1
            var spectrum = KmerAnalyzer.GetKmerSpectrum("AACGAACG", 2);

            Assert.Multiple(() =>
            {
                Assert.That(spectrum.ContainsKey(1), Is.True, "Should have k-mers appearing once");
                Assert.That(spectrum.ContainsKey(2), Is.True, "Should have k-mers appearing twice");
            });
        }

        #endregion

        #region GetKmerSpectrum - Spectrum Total Invariant (M5)

        /// <summary>
        /// M5: Sum of (multiplicity × count) equals L - k + 1.
        /// Evidence: Mathematical property - total k-mers must be conserved.
        /// </summary>
        [Test]
        [TestCase("ACGTACGT", 2)]
        [TestCase("ACGTACGT", 3)]
        [TestCase("ACGTACGT", 4)]
        [TestCase("AAAAAAAAAA", 3)]
        public void GetKmerSpectrum_TotalInvariant_SumEqualsLMinusKPlusOne(string sequence, int k)
        {
            var spectrum = KmerAnalyzer.GetKmerSpectrum(sequence, k);
            int expectedTotal = sequence.Length - k + 1;

            // Sum of (multiplicity × number_of_kmers_with_that_multiplicity)
            int actualTotal = spectrum.Sum(kvp => kvp.Key * kvp.Value);

            Assert.That(actualTotal, Is.EqualTo(expectedTotal),
                $"Spectrum total must equal L - k + 1 = {expectedTotal}");
        }

        #endregion

        #region GetKmerSpectrum - Edge Cases

        /// <summary>
        /// Edge case: Empty sequence returns empty spectrum.
        /// </summary>
        [Test]
        public void GetKmerSpectrum_EmptySequence_ReturnsEmptyDictionary()
        {
            var spectrum = KmerAnalyzer.GetKmerSpectrum("", 4);
            Assert.That(spectrum, Is.Empty);
        }

        /// <summary>
        /// Edge case: k > sequence length returns empty spectrum.
        /// </summary>
        [Test]
        public void GetKmerSpectrum_KGreaterThanLength_ReturnsEmptyDictionary()
        {
            var spectrum = KmerAnalyzer.GetKmerSpectrum("ACG", 5);
            Assert.That(spectrum, Is.Empty);
        }

        #endregion

        #region CalculateKmerEntropy - Zero Entropy (M6)

        /// <summary>
        /// M6: Homopolymer has zero entropy (only one k-mer type).
        /// Evidence: Shannon (1948) - "When entropy is zero, there is no uncertainty"
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_Homopolymer_ReturnsZero()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("AAAA", 2);
            Assert.That(entropy, Is.EqualTo(0.0).Within(Tolerance));
        }

        /// <summary>
        /// M6: Single possible k-mer (k = sequence length) has zero entropy.
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_SingleKmer_ReturnsZero()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACGT", 4);
            Assert.That(entropy, Is.EqualTo(0.0).Within(Tolerance),
                "Single k-mer has complete certainty, thus zero entropy");
        }

        #endregion

        #region CalculateKmerEntropy - Maximum Entropy (M7)

        /// <summary>
        /// M7: Uniform distribution yields maximum entropy = log₂(n).
        /// Evidence: Shannon (1948) - "Maximum uncertainty when all outcomes are equally likely"
        /// 
        /// For "ACGT" with k=1: 4 unique bases, each appears once.
        /// H = -4 × (0.25 × log₂(0.25)) = -4 × (0.25 × -2) = 2.0 bits
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_UniformDistribution_ReturnsMaxEntropy()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACGT", 1);
            double expectedMaxEntropy = Math.Log2(4); // 2.0 bits

            Assert.That(entropy, Is.EqualTo(expectedMaxEntropy).Within(Tolerance),
                "Uniform distribution of 4 symbols should have entropy = log₂(4) = 2.0 bits");
        }

        /// <summary>
        /// M7 variant: All distinct k-mers have maximum entropy.
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_AllDistinctKmers_ReturnsLogOfCount()
        {
            // "ACGT" with k=2: AC, CG, GT (3 distinct 2-mers, each appears once)
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACGT", 2);
            int uniqueKmers = 3;
            double expectedMaxEntropy = Math.Log2(uniqueKmers);

            Assert.That(entropy, Is.EqualTo(expectedMaxEntropy).Within(Tolerance));
        }

        #endregion

        #region CalculateKmerEntropy - Edge Cases (M8)

        /// <summary>
        /// M8: Empty sequence returns 0.0 entropy.
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_EmptySequence_ReturnsZero()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("", 2);
            Assert.That(entropy, Is.EqualTo(0.0).Within(Tolerance));
        }

        /// <summary>
        /// M8: k > sequence length returns 0.0 entropy.
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_KGreaterThanLength_ReturnsZero()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACG", 5);
            Assert.That(entropy, Is.EqualTo(0.0).Within(Tolerance));
        }

        #endregion

        #region CalculateKmerEntropy - Bounds Invariant (M9)

        /// <summary>
        /// M9: Entropy is always between 0 and log₂(unique_kmer_count).
        /// Evidence: Shannon (1948) - entropy bounds property.
        /// </summary>
        [Test]
        [TestCase("ACGTACGT", 2)]
        [TestCase("ACGTACGTACGT", 3)]
        [TestCase("AAACCCGGGTTT", 2)]
        [TestCase("ATATATAT", 2)]
        public void CalculateKmerEntropy_BoundsInvariant_EntropWithinValidRange(string sequence, int k)
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies(sequence, k);

            if (frequencies.Count == 0)
            {
                return; // Skip if no k-mers
            }

            double entropy = KmerAnalyzer.CalculateKmerEntropy(sequence, k);
            int uniqueKmers = frequencies.Count;
            double maxEntropy = Math.Log2(uniqueKmers);

            Assert.Multiple(() =>
            {
                Assert.That(entropy, Is.GreaterThanOrEqualTo(0.0),
                    "Entropy must be non-negative");
                Assert.That(entropy, Is.LessThanOrEqualTo(maxEntropy + Tolerance),
                    $"Entropy must not exceed log₂({uniqueKmers}) = {maxEntropy}");
            });
        }

        #endregion

        #region Should Tests - Additional Coverage

        /// <summary>
        /// S3: Intermediate entropy for non-uniform distribution.
        /// </summary>
        [Test]
        public void CalculateKmerEntropy_NonUniformDistribution_IntermediateEntropy()
        {
            // Sequence with some repeated k-mers should have 0 < entropy < max
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACGTACGT", 2);
            var frequencies = KmerAnalyzer.GetKmerFrequencies("ACGTACGT", 2);
            double maxEntropy = Math.Log2(frequencies.Count);

            Assert.Multiple(() =>
            {
                Assert.That(entropy, Is.GreaterThan(0.0));
                Assert.That(entropy, Is.LessThan(maxEntropy));
            });
        }

        /// <summary>
        /// S4: Case insensitivity - mixed case produces same result as uppercase.
        /// </summary>
        [Test]
        public void KmerFrequencyMethods_MixedCase_SameAsUppercase()
        {
            string lowercase = "acgtacgt";
            string uppercase = "ACGTACGT";
            int k = 2;

            var freqLower = KmerAnalyzer.GetKmerFrequencies(lowercase, k);
            var freqUpper = KmerAnalyzer.GetKmerFrequencies(uppercase, k);

            var specLower = KmerAnalyzer.GetKmerSpectrum(lowercase, k);
            var specUpper = KmerAnalyzer.GetKmerSpectrum(uppercase, k);

            double entropyLower = KmerAnalyzer.CalculateKmerEntropy(lowercase, k);
            double entropyUpper = KmerAnalyzer.CalculateKmerEntropy(uppercase, k);

            Assert.Multiple(() =>
            {
                Assert.That(freqLower.Count, Is.EqualTo(freqUpper.Count));
                Assert.That(specLower.Count, Is.EqualTo(specUpper.Count));
                Assert.That(entropyLower, Is.EqualTo(entropyUpper).Within(Tolerance));
            });
        }

        #endregion

        #region Could Tests - Extended Coverage

        /// <summary>
        /// C2: Methods work correctly across various k values.
        /// </summary>
        [Test]
        public void KmerFrequencyMethods_VariousKValues_AllMethodsWork()
        {
            const string sequence = "ACGTACGTACGTACGT";

            Assert.Multiple(() =>
            {
                for (int k = 1; k <= 5; k++)
                {
                    var frequencies = KmerAnalyzer.GetKmerFrequencies(sequence, k);
                    var spectrum = KmerAnalyzer.GetKmerSpectrum(sequence, k);
                    double entropy = KmerAnalyzer.CalculateKmerEntropy(sequence, k);

                    Assert.That(frequencies.Count, Is.GreaterThan(0), $"k={k}: should have frequencies");
                    Assert.That(spectrum.Count, Is.GreaterThan(0), $"k={k}: should have spectrum");
                    Assert.That(entropy, Is.GreaterThanOrEqualTo(0), $"k={k}: entropy should be non-negative");
                }
            });
        }

        #endregion
    }
}
