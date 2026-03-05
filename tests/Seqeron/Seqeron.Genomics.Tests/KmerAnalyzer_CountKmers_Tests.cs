using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for KMER-COUNT-001: K-mer Counting.
    /// 
    /// Canonical methods: KmerAnalyzer.CountKmers, CountKmersSpan, CountKmersBothStrands
    /// 
    /// Evidence:
    /// - Wikipedia: K-mer definition, L − k + 1 formula, algorithm pseudocode
    /// - Rosalind KMER: K-mer composition problem with sample dataset
    /// - Rosalind BA1E: Clump finding for k-mer counting validation
    /// </summary>
    [TestFixture]
    public class KmerAnalyzer_CountKmers_Tests
    {
        #region Edge Cases - Empty and Boundary Conditions

        /// <summary>
        /// M1: Empty sequence returns empty dictionary.
        /// Evidence: Wikipedia pseudocode - loop from 0 to L-k+1 yields nothing when L=0.
        /// </summary>
        [Test]
        public void CountKmers_EmptySequence_ReturnsEmptyDictionary()
        {
            var counts = KmerAnalyzer.CountKmers("", 4);
            Assert.That(counts, Is.Empty);
        }

        /// <summary>
        /// M1 variant: Null sequence returns empty dictionary.
        /// Evidence: Defensive programming for null inputs.
        /// </summary>
        [Test]
        public void CountKmers_NullSequence_ReturnsEmptyDictionary()
        {
            string? nullSequence = null;
            var counts = KmerAnalyzer.CountKmers(nullSequence!, 4);
            Assert.That(counts, Is.Empty);
        }

        /// <summary>
        /// M2: k > sequence length returns empty dictionary.
        /// Evidence: Wikipedia formula L − k + 1 becomes negative/zero.
        /// </summary>
        [Test]
        public void CountKmers_KLargerThanSequence_ReturnsEmptyDictionary()
        {
            var counts = KmerAnalyzer.CountKmers("ACG", 4);
            Assert.That(counts, Is.Empty);
        }

        /// <summary>
        /// M3: k ≤ 0 throws ArgumentOutOfRangeException.
        /// Evidence: k must be positive for valid k-mer definition.
        /// </summary>
        [Test]
        [TestCase(0)]
        [TestCase(-1)]
        [TestCase(-10)]
        public void CountKmers_InvalidK_ThrowsArgumentOutOfRangeException(int k)
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                KmerAnalyzer.CountKmers("ACGT", k));
        }

        /// <summary>
        /// S4: k = sequence length yields single k-mer with count 1.
        /// Evidence: L − k + 1 = 1 when k = L.
        /// </summary>
        [Test]
        public void CountKmers_KEqualSequenceLength_ReturnsSingleKmer()
        {
            var counts = KmerAnalyzer.CountKmers("ACGT", 4);

            Assert.That(counts, Has.Count.EqualTo(1));
            Assert.That(counts["ACGT"], Is.EqualTo(1));
        }

        #endregion

        #region Invariant: Total Count = L - k + 1

        /// <summary>
        /// M5: Total count invariant - sum of all k-mer counts equals L − k + 1.
        /// Evidence: Wikipedia - "a sequence of length L will have L − k + 1 k-mers".
        /// </summary>
        [Test]
        [TestCase("ACGT", 2)]            // 4 - 2 + 1 = 3
        [TestCase("ACGTACGT", 4)]        // 8 - 4 + 1 = 5
        [TestCase("AAAAAAAAAA", 3)]      // 10 - 3 + 1 = 8
        [TestCase("ATCGATCGATCG", 1)]    // 12 - 1 + 1 = 12
        public void CountKmers_TotalCountInvariant_SumEqualsLMinusKPlusOne(string sequence, int k)
        {
            var counts = KmerAnalyzer.CountKmers(sequence, k);
            int expectedTotal = sequence.Length - k + 1;
            int actualTotal = counts.Values.Sum();

            Assert.That(actualTotal, Is.EqualTo(expectedTotal),
                $"Sum of k-mer counts should equal L - k + 1 = {expectedTotal}");
        }

        /// <summary>
        /// M5 property: Invariant holds for various k values on same sequence.
        /// Evidence: Wikipedia formula applies universally.
        /// </summary>
        [Test]
        public void CountKmers_TotalCountInvariant_HoldsForAllValidK()
        {
            const string sequence = "ACGTACGTACGT";
            int L = sequence.Length;

            Assert.Multiple(() =>
            {
                for (int k = 1; k <= L; k++)
                {
                    var counts = KmerAnalyzer.CountKmers(sequence, k);
                    int expected = L - k + 1;
                    int actual = counts.Values.Sum();

                    Assert.That(actual, Is.EqualTo(expected),
                        $"Failed for k={k}: expected {expected}, got {actual}");
                }
            });
        }

        #endregion

        #region Counting Correctness

        /// <summary>
        /// M8: Distinct k-mers counted correctly.
        /// Evidence: Rosalind KMER problem - distinct k-mers identified.
        /// </summary>
        [Test]
        public void CountKmers_SimpleSequence_CountsDistinctKmersCorrectly()
        {
            var counts = KmerAnalyzer.CountKmers("ACGTACGT", 4);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(4)); // ACGT, CGTA, GTAC, TACG
                Assert.That(counts["ACGT"], Is.EqualTo(2));
                Assert.That(counts["CGTA"], Is.EqualTo(1));
                Assert.That(counts["GTAC"], Is.EqualTo(1));
                Assert.That(counts["TACG"], Is.EqualTo(1));
            });
        }

        /// <summary>
        /// M6: Homopolymer sequence - single k-mer with count L − k + 1.
        /// Evidence: Wikipedia - all same bases produce repeated identical k-mers.
        /// </summary>
        [Test]
        [TestCase("AAAA", 2, "AA", 3)]
        [TestCase("TTTTTT", 3, "TTT", 4)]
        [TestCase("GGGGGGGG", 4, "GGGG", 5)]
        [TestCase("CCCC", 4, "CCCC", 1)]
        public void CountKmers_Homopolymer_SingleKmerWithCorrectCount(
            string sequence, int k, string expectedKmer, int expectedCount)
        {
            var counts = KmerAnalyzer.CountKmers(sequence, k);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(1));
                Assert.That(counts.ContainsKey(expectedKmer), Is.True);
                Assert.That(counts[expectedKmer], Is.EqualTo(expectedCount));
            });
        }

        /// <summary>
        /// M9: Overlapping k-mers counted correctly.
        /// Evidence: Wikipedia sliding window - overlapping substrings extracted.
        /// </summary>
        [Test]
        public void CountKmers_OverlappingKmers_AllCounted()
        {
            // "AAACGT" has 2-mers: AA, AA, AC, CG, GT
            var counts = KmerAnalyzer.CountKmers("AAACGT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts["AA"], Is.EqualTo(2)); // Overlapping AAs
                Assert.That(counts["AC"], Is.EqualTo(1));
                Assert.That(counts["CG"], Is.EqualTo(1));
                Assert.That(counts["GT"], Is.EqualTo(1));
            });
        }

        /// <summary>
        /// S3: k = 1 counts individual nucleotides.
        /// Evidence: Wikipedia - k=1 gives "monomers" (individual bases).
        /// </summary>
        [Test]
        public void CountKmers_KEqualsOne_CountsNucleotides()
        {
            var counts = KmerAnalyzer.CountKmers("AACGT", 1);

            Assert.Multiple(() =>
            {
                Assert.That(counts["A"], Is.EqualTo(2));
                Assert.That(counts["C"], Is.EqualTo(1));
                Assert.That(counts["G"], Is.EqualTo(1));
                Assert.That(counts["T"], Is.EqualTo(1));
            });
        }

        #endregion

        #region Case Sensitivity

        /// <summary>
        /// M7: Case-insensitive counting - lowercase normalized to uppercase.
        /// Evidence: Wikipedia/algorithm norm — k-mers are case-insensitive.
        /// Input: "acgt", k=2 → 3 two-mers: AC, CG, GT (all uppercase, each count=1).
        /// </summary>
        [Test]
        public void CountKmers_LowercaseSequence_NormalizedToUppercase()
        {
            var counts = KmerAnalyzer.CountKmers("acgt", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(3));
                Assert.That(counts.Values.Sum(), Is.EqualTo(3));
                Assert.That(counts["AC"], Is.EqualTo(1));
                Assert.That(counts["CG"], Is.EqualTo(1));
                Assert.That(counts["GT"], Is.EqualTo(1));
                Assert.That(counts.ContainsKey("ac"), Is.False);
            });
        }

        /// <summary>
        /// S1: Mixed case input counted as same k-mer.
        /// Evidence: DNA sequences should be case-insensitive.
        /// Input: "AcGtACgt", k=4 → 5 four-mers: ACGT(×2), CGTA, GTAC, TACG.
        /// </summary>
        [Test]
        public void CountKmers_MixedCase_TreatedAsSameKmer()
        {
            var counts = KmerAnalyzer.CountKmers("AcGtACgt", 4);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(4));
                Assert.That(counts.Values.Sum(), Is.EqualTo(5));
                Assert.That(counts["ACGT"], Is.EqualTo(2));
                Assert.That(counts["CGTA"], Is.EqualTo(1));
                Assert.That(counts["GTAC"], Is.EqualTo(1));
                Assert.That(counts["TACG"], Is.EqualTo(1));
            });
        }

        #endregion

        #region Non-Standard Characters

        /// <summary>
        /// S2: Non-DNA characters (like N) are counted as-is.
        /// Evidence: Genomic data often contains N for unknown bases.
        /// Input: "ACNGT", k=2 → 4 two-mers: AC, CN, NG, GT (each count=1).
        /// </summary>
        [Test]
        public void CountKmers_WithAmbiguousBase_CountedAsIs()
        {
            var counts = KmerAnalyzer.CountKmers("ACNGT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(4));
                Assert.That(counts.Values.Sum(), Is.EqualTo(4));
                Assert.That(counts["AC"], Is.EqualTo(1));
                Assert.That(counts["CN"], Is.EqualTo(1));
                Assert.That(counts["NG"], Is.EqualTo(1));
                Assert.That(counts["GT"], Is.EqualTo(1));
            });
        }

        #endregion

        #region Span-Based API (CountKmersSpan)

        /// <summary>
        /// M10: CountKmersSpan produces same results as CountKmers.
        /// Evidence: API consistency - both should implement same algorithm.
        /// </summary>
        [Test]
        public void CountKmersSpan_ProducesSameResultAsCountKmers()
        {
            const string sequence = "ACGTACGTACGT";
            const int k = 4;

            var stringCounts = KmerAnalyzer.CountKmers(sequence, k);
            var spanCounts = sequence.AsSpan().CountKmersSpan(k);

            Assert.Multiple(() =>
            {
                Assert.That(spanCounts.Count, Is.EqualTo(stringCounts.Count));
                foreach (var kvp in stringCounts)
                {
                    Assert.That(spanCounts.ContainsKey(kvp.Key), Is.True,
                        $"Span result missing key: {kvp.Key}");
                    Assert.That(spanCounts[kvp.Key], Is.EqualTo(kvp.Value),
                        $"Count mismatch for {kvp.Key}");
                }
            });
        }

        /// <summary>
        /// M10 variant: Span handles edge cases same as string version.
        /// Evidence: API consistency.
        /// </summary>
        [Test]
        public void CountKmersSpan_EmptySpan_ReturnsEmptyDictionary()
        {
            var counts = ReadOnlySpan<char>.Empty.CountKmersSpan(4);
            Assert.That(counts, Is.Empty);
        }

        /// <summary>
        /// Span-based counting with k > length returns empty.
        /// </summary>
        [Test]
        public void CountKmersSpan_KLargerThanSpan_ReturnsEmptyDictionary()
        {
            ReadOnlySpan<char> span = "ACG".AsSpan();
            var counts = span.CountKmersSpan(4);
            Assert.That(counts, Is.Empty);
        }

        /// <summary>
        /// CountKmersSpan must throw for k ≤ 0, matching CountKmers behavior.
        /// Evidence: k must be positive for valid k-mer definition (Wikipedia).
        /// </summary>
        [Test]
        [TestCase(0)]
        [TestCase(-1)]
        public void CountKmersSpan_InvalidK_ThrowsArgumentOutOfRangeException(int k)
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                ReadOnlySpan<char> span = "ACGT".AsSpan();
                span.CountKmersSpan(k);
            });
        }

        /// <summary>
        /// Mutation killer: k == sequence.Length should return exactly one k-mer 
        /// (boundary for Length &lt; k guard — must NOT be empty).
        /// </summary>
        [Test]
        public void CountKmersSpan_KEqualToLength_ReturnsSingleKmer()
        {
            ReadOnlySpan<char> span = "ACGT".AsSpan();
            var counts = span.CountKmersSpan(4);
            Assert.That(counts, Has.Count.EqualTo(1));
            Assert.That(counts["ACGT"], Is.EqualTo(1));
        }

        #endregion

        #region Both Strands (CountKmersBothStrands)

        /// <summary>
        /// M11: CountKmersBothStrands combines forward and reverse complement counts.
        /// Evidence: DNA double-helix - both strands are biologically relevant.
        /// </summary>
        [Test]
        public void CountKmersBothStrands_CombinesForwardAndReverseComplement()
        {
            // ACGT reverse complement is ACGT (palindromic)
            // Forward 2-mers: AC, CG, GT
            // RevComp 2-mers: AC, CG, GT (ACGT → ACGT)
            var dna = new DnaSequence("ACGT");
            var counts = KmerAnalyzer.CountKmersBothStrands(dna, 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts["AC"], Is.EqualTo(2)); // 1 forward + 1 revcomp
                Assert.That(counts["CG"], Is.EqualTo(2));
                Assert.That(counts["GT"], Is.EqualTo(2));
            });
        }

        /// <summary>
        /// M11 variant: Non-palindromic sequence adds different k-mers from reverse complement.
        /// </summary>
        [Test]
        public void CountKmersBothStrands_NonPalindromicSequence_AddsNewKmers()
        {
            // "AAA" forward: AA (count 2)
            // "AAA" reverse complement: TTT → TT (count 2)
            var dna = new DnaSequence("AAA");
            var counts = KmerAnalyzer.CountKmersBothStrands(dna, 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts["AA"], Is.EqualTo(2));
                Assert.That(counts["TT"], Is.EqualTo(2));
            });
        }

        /// <summary>
        /// M11: Total invariant - both strands contribute L − k + 1 each.
        /// </summary>
        [Test]
        public void CountKmersBothStrands_TotalCountInvariant()
        {
            var dna = new DnaSequence("ACGTACGT");
            int k = 3;
            var counts = KmerAnalyzer.CountKmersBothStrands(dna, k);

            // Each strand contributes L - k + 1 = 8 - 3 + 1 = 6
            int expectedTotal = 2 * (dna.Sequence.Length - k + 1);
            int actualTotal = counts.Values.Sum();

            Assert.That(actualTotal, Is.EqualTo(expectedTotal));
        }

        #endregion

        #region DnaSequence Wrapper

        /// <summary>
        /// Smoke test: DnaSequence overload delegates correctly.
        /// Evidence: Wrapper should produce same result as string version.
        /// </summary>
        [Test]
        public void CountKmers_DnaSequence_DelegatesToStringVersion()
        {
            var dna = new DnaSequence("ACGTACGT");
            var dnaResult = KmerAnalyzer.CountKmers(dna, 4);
            var stringResult = KmerAnalyzer.CountKmers(dna.Sequence, 4);

            Assert.That(dnaResult, Is.EqualTo(stringResult));
        }

        #endregion

        #region Wikipedia Evidence Tests

        /// <summary>
        /// Wikipedia exact example: "The sequence ATGG has two 3-mers: ATG and TGG."
        /// Source: https://en.wikipedia.org/wiki/K-mer (lead diagram caption)
        /// </summary>
        [Test]
        public void CountKmers_WikipediaExample_ATGG_TwoThreeMers()
        {
            var counts = KmerAnalyzer.CountKmers("ATGG", 3);

            Assert.Multiple(() =>
            {
                Assert.That(counts, Has.Count.EqualTo(2));
                Assert.That(counts["ATG"], Is.EqualTo(1));
                Assert.That(counts["TGG"], Is.EqualTo(1));
            });
        }

        /// <summary>
        /// Wikipedia k-mer table: sequence "GTAGAGCTGT" (length 10).
        /// k=2 yields 9 total k-mers (L−k+1 = 10−2+1 = 9), 7 unique.
        /// Listed in table: GT, TA, AG, GA, AG, GC, CT, TG, GT.
        /// Source: https://en.wikipedia.org/wiki/K-mer#Introduction
        /// </summary>
        [Test]
        public void CountKmers_WikipediaTable_GTAGAGCTGT_TwoMers()
        {
            var counts = KmerAnalyzer.CountKmers("GTAGAGCTGT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts.Values.Sum(), Is.EqualTo(9), "Total = L − k + 1 = 9");
                Assert.That(counts, Has.Count.EqualTo(7), "7 unique 2-mers");
                Assert.That(counts["GT"], Is.EqualTo(2));
                Assert.That(counts["TA"], Is.EqualTo(1));
                Assert.That(counts["AG"], Is.EqualTo(2));
                Assert.That(counts["GA"], Is.EqualTo(1));
                Assert.That(counts["GC"], Is.EqualTo(1));
                Assert.That(counts["CT"], Is.EqualTo(1));
                Assert.That(counts["TG"], Is.EqualTo(1));
            });
        }

        /// <summary>
        /// Wikipedia k-mer table: GTAGAGCTGT with various k values.
        /// Validates L − k + 1 k-mers for each k from 1 to 10.
        /// Source: https://en.wikipedia.org/wiki/K-mer#Introduction
        /// </summary>
        [Test]
        [TestCase(1, 10)]  // 10 - 1 + 1 = 10
        [TestCase(2, 9)]   // GT,TA,AG,GA,AG,GC,CT,TG,GT
        [TestCase(3, 8)]   // GTA,TAG,AGA,GAG,AGC,GCT,CTG,TGT
        [TestCase(4, 7)]   // GTAG,TAGA,AGAG,GAGC,AGCT,GCTG,CTGT
        [TestCase(5, 6)]
        [TestCase(6, 5)]
        [TestCase(7, 4)]
        [TestCase(8, 3)]
        [TestCase(9, 2)]
        [TestCase(10, 1)]  // GTAGAGCTGT
        public void CountKmers_WikipediaTable_GTAGAGCTGT_TotalKmersPerK(int k, int expectedTotal)
        {
            var counts = KmerAnalyzer.CountKmers("GTAGAGCTGT", k);
            Assert.That(counts.Values.Sum(), Is.EqualTo(expectedTotal));
        }

        /// <summary>
        /// Wikipedia k-mer table: GTAGAGCTGT k=4 yields 7 unique four-mers, each with count 1.
        /// Table: GTAG, TAGA, AGAG, GAGC, AGCT, GCTG, CTGT.
        /// Source: https://en.wikipedia.org/wiki/K-mer#Introduction
        /// </summary>
        [Test]
        public void CountKmers_WikipediaTable_GTAGAGCTGT_FourMers()
        {
            var counts = KmerAnalyzer.CountKmers("GTAGAGCTGT", 4);

            Assert.Multiple(() =>
            {
                Assert.That(counts.Values.Sum(), Is.EqualTo(7), "Total = L − k + 1 = 7");
                Assert.That(counts, Has.Count.EqualTo(7), "All 7 four-mers are unique");
                Assert.That(counts["GTAG"], Is.EqualTo(1));
                Assert.That(counts["TAGA"], Is.EqualTo(1));
                Assert.That(counts["AGAG"], Is.EqualTo(1));
                Assert.That(counts["GAGC"], Is.EqualTo(1));
                Assert.That(counts["AGCT"], Is.EqualTo(1));
                Assert.That(counts["GCTG"], Is.EqualTo(1));
                Assert.That(counts["CTGT"], Is.EqualTo(1));
            });
        }

        #endregion

        #region Rosalind KMER Evidence Tests

        /// <summary>
        /// Rosalind KMER problem: 4-mer composition of sample dataset.
        /// Validates specific k-mer counts from the known sample output (256 values).
        /// Source: https://rosalind.info/problems/kmer/ (Sample Dataset / Sample Output)
        /// </summary>
        [Test]
        public void CountKmers_RosalindKmerSample_SpecificCounts()
        {
            const string rosalindSequence =
                "CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG" +
                "CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT" +
                "TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA" +
                "AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG" +
                "GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA" +
                "CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA" +
                "CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG";

            var counts = KmerAnalyzer.CountKmers(rosalindSequence, 4);

            Assert.Multiple(() =>
            {
                // Total invariant: L − k + 1 = 415 − 4 + 1 = 412
                Assert.That(counts.Values.Sum(), Is.EqualTo(412),
                    "Rosalind: total 4-mers = 412");

                // Specific counts from Rosalind sample output (lexicographic positions verified)
                Assert.That(counts["AAAA"], Is.EqualTo(4), "Rosalind position 0");
                Assert.That(counts["AACT"], Is.EqualTo(5), "Rosalind position 7");
                Assert.That(counts["AGCT"], Is.EqualTo(3), "Rosalind position 39");
                Assert.That(counts["GCTG"], Is.EqualTo(5), "Rosalind position 158");
                Assert.That(counts["TATA"], Is.EqualTo(2), "Rosalind position 204");
                Assert.That(counts["TTTT"], Is.EqualTo(1), "Rosalind position 255");

                // Zero-count k-mer should not appear in dictionary
                Assert.That(counts.ContainsKey("CCCC"), Is.False,
                    "Rosalind position 85 = 0, should not be in dictionary");
            });
        }

        /// <summary>
        /// Rosalind KMER: exact unique 4-mer count = 209 (from 256 possible, 47 absent).
        /// Source: https://rosalind.info/problems/kmer/ (Sample Output: 47 zeros in 256 values)
        /// </summary>
        [Test]
        public void CountKmers_RosalindKmerSample_ExactUniqueCount()
        {
            const string rosalindSequence =
                "CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG" +
                "CCCTGACTACCGAACCAGTTGTGAGTACTCAACTGGGTGAGAGTGCAGTCCCTATTGAGT" +
                "TTCCGAGACTCACCGGGATTTTCGATCCAGCCTCAGTCCAGTCTTGTGGCCAACTCACCA" +
                "AATGACGTTGGAATATCCCTGTCTAGCTCACGCAGTACTTAGTAAGAGGTCGCTGCAGCG" +
                "GGGCAAGGAGATCGGAAAATGTGCTCTATATGCGACTAAAGCTCCTAACTTACACGTAGA" +
                "CTTGCCCGTGTTAAAAACTCGGCTCACATGCTGTCTGCGGCTGGCTGTATACAGTATCTA" +
                "CCTAATACCCTTCAGTTCGCCGCACAAAAGCTGGGAGTTACCGCGGAAATCACAG";

            var counts = KmerAnalyzer.CountKmers(rosalindSequence, 4);

            // Rosalind sample output has 256 values, 47 are zero → 209 unique k-mers
            Assert.That(counts.Count, Is.EqualTo(209),
                "Rosalind sample: 209 unique 4-mers (256 possible minus 47 absent)");
        }

        #endregion

        #region Case Normalization Regression Tests

        /// <summary>
        /// Regression: CountKmersSpan must normalize to uppercase (same as CountKmers).
        /// Deviation found: CountKmersSpan previously did NOT call ToUpperInvariant().
        /// Evidence: Wikipedia/algorithm doc — k-mers are case-insensitive.
        /// </summary>
        [Test]
        public void CountKmersSpan_LowercaseInput_NormalizesToUppercase()
        {
            ReadOnlySpan<char> span = "acgtacgt".AsSpan();
            var counts = span.CountKmersSpan(4);

            Assert.Multiple(() =>
            {
                Assert.That(counts.ContainsKey("ACGT"), Is.True, "Lowercase should normalize to ACGT");
                Assert.That(counts["ACGT"], Is.EqualTo(2));
                Assert.That(counts.ContainsKey("acgt"), Is.False, "Lowercase key must not exist");
            });
        }

        /// <summary>
        /// Regression: CountKmersSpan mixed case must match CountKmers result.
        /// This test would have caught the original deviation.
        /// </summary>
        [Test]
        public void CountKmersSpan_MixedCase_MatchesCountKmers()
        {
            const string sequence = "AcGtAcGt";
            var stringCounts = KmerAnalyzer.CountKmers(sequence, 3);
            var spanCounts = sequence.AsSpan().CountKmersSpan(3);

            Assert.Multiple(() =>
            {
                Assert.That(spanCounts.Count, Is.EqualTo(stringCounts.Count),
                    "Span and string must produce same number of unique k-mers");
                foreach (var kvp in stringCounts)
                {
                    Assert.That(spanCounts.ContainsKey(kvp.Key), Is.True,
                        $"Span result missing key: {kvp.Key}");
                    Assert.That(spanCounts[kvp.Key], Is.EqualTo(kvp.Value),
                        $"Count mismatch for {kvp.Key}");
                }
            });
        }

        /// <summary>
        /// Regression: CancellationToken overload must normalize to uppercase.
        /// Deviation found: previously used sequence.AsSpan() without ToUpperInvariant().
        /// </summary>
        [Test]
        public void CountKmers_CancellationOverload_NormalizesCase()
        {
            var cts = new CancellationTokenSource();
            var counts = KmerAnalyzer.CountKmers("acgtacgt", 4, cts.Token);

            Assert.Multiple(() =>
            {
                Assert.That(counts.ContainsKey("ACGT"), Is.True, "Lowercase should normalize to ACGT");
                Assert.That(counts["ACGT"], Is.EqualTo(2));
                Assert.That(counts.ContainsKey("acgt"), Is.False, "Lowercase key must not exist");
            });
        }

        #endregion
    }
}
