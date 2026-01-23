using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    /// <summary>
    /// Tests for K-mer search operations: FindMostFrequentKmers, FindUniqueKmers, FindClumps.
    /// 
    /// Test Unit: KMER-FIND-001
    /// Evidence: Wikipedia (K-mer), Rosalind BA1B (frequent words), Rosalind BA1E (clump finding)
    /// </summary>
    [TestFixture]
    public class KmerAnalyzer_Find_Tests
    {
        #region FindMostFrequentKmers - Must Tests

        /// <summary>
        /// M1: Rosalind BA1B sample dataset.
        /// Source: https://rosalind.info/problems/ba1b/
        /// Sequence "ACGTTGCATGTCGCATGATGCATGAGAGCT" with k=4 should return CATG and GCAT.
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_RosalindBA1B_Sample_ReturnsCatgAndGcat()
        {
            // Arrange
            const string sequence = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
            const int k = 4;

            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers(sequence, k).ToList();

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result, Does.Contain("CATG"), "CATG should be a most frequent 4-mer");
                Assert.That(result, Does.Contain("GCAT"), "GCAT should be a most frequent 4-mer");
            });
        }

        /// <summary>
        /// M2: When multiple k-mers share maximum count, all should be returned.
        /// Source: Rosalind BA1B definition: "All most frequent k-mers"
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_TiedMaxCount_ReturnsAllTied()
        {
            // Arrange - ACGT appears twice in ACGTACGT, all others appear once
            const string sequence = "ACGTACGT";
            const int k = 4;

            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers(sequence, k).ToList();

            // Assert - only ACGT appears twice
            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1), "Only one 4-mer appears twice");
                Assert.That(result, Does.Contain("ACGT"), "ACGT is the most frequent");
            });
        }

        /// <summary>
        /// M3: Single most frequent k-mer.
        /// Source: Rosalind BA1B definition
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_SingleMostFrequent_ReturnsIt()
        {
            // Arrange - AA appears twice in AAACGT, others appear once
            const string sequence = "AAACGT";
            const int k = 2;

            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers(sequence, k).ToList();

            // Assert
            Assert.That(result, Does.Contain("AA"), "AA appears twice, should be most frequent");
        }

        /// <summary>
        /// M4: Empty sequence returns empty collection.
        /// Source: Implementation contract (edge case)
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_EmptySequence_ReturnsEmpty()
        {
            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers("", 4).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// M5: When k > sequence length, no k-mers exist.
        /// Source: Wikipedia pseudocode: L - k + 1 k-mers; if L &lt; k, no k-mers
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_KGreaterThanSequenceLength_ReturnsEmpty()
        {
            // Arrange
            const string sequence = "ACG";
            const int k = 5;

            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers(sequence, k).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        #endregion

        #region FindMostFrequentKmers - Should Tests

        /// <summary>
        /// S1: Case insensitivity.
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_MixedCase_TreatedAsUppercase()
        {
            // Arrange
            const string sequenceLower = "aaacgt";
            const string sequenceUpper = "AAACGT";
            const int k = 2;

            // Act
            var resultLower = KmerAnalyzer.FindMostFrequentKmers(sequenceLower, k).ToList();
            var resultUpper = KmerAnalyzer.FindMostFrequentKmers(sequenceUpper, k).ToList();

            // Assert
            Assert.That(resultLower, Is.EquivalentTo(resultUpper));
        }

        /// <summary>
        /// S2: Single character sequence with k=1.
        /// </summary>
        [Test]
        public void FindMostFrequentKmers_SingleCharacter_ReturnsThatCharacter()
        {
            // Arrange
            const string sequence = "A";
            const int k = 1;

            // Act
            var result = KmerAnalyzer.FindMostFrequentKmers(sequence, k).ToList();

            // Assert
            Assert.That(result, Is.EquivalentTo(new[] { "A" }));
        }

        #endregion

        #region FindUniqueKmers - Must Tests

        /// <summary>
        /// M6: Basic uniqueness detection.
        /// Source: Wikipedia (K-mer): k-mers appearing exactly once
        /// In ACGTACGT with k=4: ACGT appears 2×, CGTA/GTAC/TACG appear 1×.
        /// </summary>
        [Test]
        public void FindUniqueKmers_ReturnsKmersAppearingOnce()
        {
            // Arrange
            const string sequence = "ACGTACGT";
            const int k = 4;

            // Act
            var result = KmerAnalyzer.FindUniqueKmers(sequence, k).ToList();

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result, Does.Contain("CGTA"), "CGTA appears once");
                Assert.That(result, Does.Contain("GTAC"), "GTAC appears once");
                Assert.That(result, Does.Contain("TACG"), "TACG appears once");
                Assert.That(result, Does.Not.Contain("ACGT"), "ACGT appears twice, not unique");
            });
        }

        /// <summary>
        /// M7: When no k-mer repeats, all are unique.
        /// Source: Mathematical invariant
        /// </summary>
        [Test]
        public void FindUniqueKmers_AllUnique_ReturnsAll()
        {
            // Arrange - "ACGT" has 3 distinct 2-mers, each appearing once
            const string sequence = "ACGT";
            const int k = 2;

            // Act
            var result = KmerAnalyzer.FindUniqueKmers(sequence, k).ToList();

            // Assert
            Assert.That(result, Has.Count.EqualTo(3), "All 3 k-mers are unique");
        }

        /// <summary>
        /// M8: Homopolymer has no unique k-mers.
        /// Source: Edge case
        /// </summary>
        [Test]
        public void FindUniqueKmers_Homopolymer_ReturnsEmpty()
        {
            // Arrange - "AAAA" only has "AA" which appears 3 times
            const string sequence = "AAAA";
            const int k = 2;

            // Act
            var result = KmerAnalyzer.FindUniqueKmers(sequence, k).ToList();

            // Assert
            Assert.That(result, Is.Empty, "No k-mer appears exactly once");
        }

        #endregion

        #region FindUniqueKmers - Should Tests

        /// <summary>
        /// S3: Empty sequence returns empty.
        /// </summary>
        [Test]
        public void FindUniqueKmers_EmptySequence_ReturnsEmpty()
        {
            // Act
            var result = KmerAnalyzer.FindUniqueKmers("", 3).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        #endregion

        #region FindClumps - Must Tests

        /// <summary>
        /// M9: Rosalind BA1E sample dataset.
        /// Source: https://rosalind.info/problems/ba1e/
        /// Sequence with k=5, L=75, t=4 should find CGACA, GAAGA, AATGT.
        /// </summary>
        [Test]
        public void FindClumps_RosalindBA1E_Sample_ReturnsExpectedClumps()
        {
            // Arrange - Rosalind BA1E sample
            const string sequence = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC";
            const int k = 5;
            const int windowSize = 75;
            const int minOccurrences = 4;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result, Does.Contain("CGACA"), "CGACA should form a (75,4)-clump");
                Assert.That(result, Does.Contain("GAAGA"), "GAAGA should form a (75,4)-clump");
                Assert.That(result, Does.Contain("AATGT"), "AATGT should form a (75,4)-clump");
            });
        }

        /// <summary>
        /// M10: Simple clump detection.
        /// Source: Rosalind BA1E definition
        /// "AAAAA" has "AAA" appearing 3 times in a 5-bp window.
        /// </summary>
        [Test]
        public void FindClumps_SimpleClump_Found()
        {
            // Arrange
            const string sequence = "AAAAA";
            const int k = 3;
            const int windowSize = 5;
            const int minOccurrences = 3;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.That(result, Does.Contain("AAA"), "AAA forms a (5,3)-clump");
        }

        /// <summary>
        /// M11: No clump found when threshold not met.
        /// Source: Rosalind BA1E definition
        /// </summary>
        [Test]
        public void FindClumps_NoClump_ReturnsEmpty()
        {
            // Arrange - looking for 3 occurrences but none exist
            const string sequence = "ACGT";
            const int k = 2;
            const int windowSize = 4;
            const int minOccurrences = 3;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// M12: Invalid parameters return empty.
        /// Source: Implementation contract
        /// </summary>
        [Test]
        public void FindClumps_KGreaterThanWindow_ReturnsEmpty()
        {
            // Arrange
            const string sequence = "ACGTACGT";
            const int k = 5;
            const int windowSize = 4;
            const int minOccurrences = 2;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// M12b: Empty sequence returns empty.
        /// Source: Implementation contract
        /// </summary>
        [Test]
        public void FindClumps_EmptySequence_ReturnsEmpty()
        {
            // Act
            var result = KmerAnalyzer.FindClumps("", 3, 5, 2).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        /// <summary>
        /// M12c: Window larger than sequence returns empty.
        /// Source: Implementation contract
        /// </summary>
        [Test]
        public void FindClumps_WindowLargerThanSequence_ReturnsEmpty()
        {
            // Arrange
            const string sequence = "ACGT";
            const int k = 2;
            const int windowSize = 10;
            const int minOccurrences = 2;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.That(result, Is.Empty);
        }

        #endregion

        #region FindClumps - Should Tests

        /// <summary>
        /// S4: Clump at boundary.
        /// </summary>
        [Test]
        public void FindClumps_ClumpAtStart_Found()
        {
            // Arrange - AAA appears 3 times at start within window of 5
            const string sequence = "AAAAACGT";
            const int k = 2;
            const int windowSize = 4;
            const int minOccurrences = 3;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.That(result, Does.Contain("AA"), "AA appears 3 times in window at start");
        }

        /// <summary>
        /// S5: Multiple distinct clumps.
        /// </summary>
        [Test]
        public void FindClumps_MultipleDistinctClumps_ReturnsAll()
        {
            // Arrange - Both AA and CC form clumps
            const string sequence = "AAAAACCCCC";
            const int k = 2;
            const int windowSize = 5;
            const int minOccurrences = 3;

            // Act
            var result = KmerAnalyzer.FindClumps(sequence, k, windowSize, minOccurrences).ToList();

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result, Does.Contain("AA"), "AA forms a clump");
                Assert.That(result, Does.Contain("CC"), "CC forms a clump");
            });
        }

        #endregion
    }
}
