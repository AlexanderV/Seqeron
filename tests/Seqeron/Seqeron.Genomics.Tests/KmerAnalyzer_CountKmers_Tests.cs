using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

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
        /// Evidence: Implementation-specific - sequences may have mixed case.
        /// </summary>
        [Test]
        public void CountKmers_LowercaseSequence_NormalizedToUppercase()
        {
            var counts = KmerAnalyzer.CountKmers("acgt", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts.ContainsKey("AC"), Is.True);
                Assert.That(counts.ContainsKey("CG"), Is.True);
                Assert.That(counts.ContainsKey("GT"), Is.True);
                // Lowercase keys should not exist
                Assert.That(counts.ContainsKey("ac"), Is.False);
            });
        }

        /// <summary>
        /// S1: Mixed case input counted as same k-mer.
        /// Evidence: DNA sequences should be case-insensitive.
        /// </summary>
        [Test]
        public void CountKmers_MixedCase_TreatedAsSameKmer()
        {
            var counts = KmerAnalyzer.CountKmers("AcGtACgt", 4);

            // Both "AcGt" and "ACgt" should normalize to "ACGT"
            Assert.That(counts["ACGT"], Is.EqualTo(2));
        }

        #endregion

        #region Non-Standard Characters

        /// <summary>
        /// S2: Non-DNA characters (like N) are counted as-is.
        /// Evidence: Genomic data often contains N for unknown bases.
        /// </summary>
        [Test]
        public void CountKmers_WithAmbiguousBase_CountedAsIs()
        {
            var counts = KmerAnalyzer.CountKmers("ACNGT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts.ContainsKey("AC"), Is.True);
                Assert.That(counts.ContainsKey("CN"), Is.True);
                Assert.That(counts.ContainsKey("NG"), Is.True);
                Assert.That(counts.ContainsKey("GT"), Is.True);
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

        #region Real-World Scenarios

        /// <summary>
        /// Biological example: Finding TATA box occurrences.
        /// Evidence: Wikipedia - TATA box is a common promoter motif.
        /// </summary>
        [Test]
        public void CountKmers_PromoterAnalysis_FindsTataBox()
        {
            const string promoter = "GCGCGCTATAAAAGGGGCTATAAAAATTT";
            var counts = KmerAnalyzer.CountKmers(promoter, 4);

            Assert.That(counts["TATA"], Is.EqualTo(2));
        }

        /// <summary>
        /// Rosalind-style: Verify k-mer composition for known sequence.
        /// Evidence: Rosalind KMER problem - lexicographic k-mer counting.
        /// </summary>
        [Test]
        public void CountKmers_RosalindStyle_CorrectComposition()
        {
            // Simple verification using known 2-mer counts
            var counts = KmerAnalyzer.CountKmers("AACGTT", 2);

            Assert.Multiple(() =>
            {
                Assert.That(counts["AA"], Is.EqualTo(1));
                Assert.That(counts["AC"], Is.EqualTo(1));
                Assert.That(counts["CG"], Is.EqualTo(1));
                Assert.That(counts["GT"], Is.EqualTo(1));
                Assert.That(counts["TT"], Is.EqualTo(1));
            });
        }

        #endregion
    }
}
