using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for auxiliary K-mer analysis methods not covered by dedicated test units.
    /// 
    /// Covered by dedicated test files:
    /// - KMER-COUNT-001: KmerAnalyzer_CountKmers_Tests.cs
    /// - KMER-FREQ-001: KmerAnalyzer_Frequency_Tests.cs
    /// - KMER-FIND-001: KmerAnalyzer_Find_Tests.cs
    /// 
    /// This file contains tests for auxiliary methods:
    /// - KmerDistance
    /// - FindKmersWithMinCount
    /// - GenerateAllKmers
    /// - FindKmerPositions
    /// - AnalyzeKmers
    /// </summary>
    [TestFixture]
    public class KmerAnalyzerTests
    {
        #region K-mer Distance

        [Test]
        public void KmerDistance_IdenticalSequences_ReturnsZero()
        {
            double distance = KmerAnalyzer.KmerDistance("ACGT", "ACGT", 2);
            Assert.That(distance, Is.EqualTo(0.0).Within(0.0001));
        }

        [Test]
        public void KmerDistance_DifferentSequences_ReturnsPositive()
        {
            double distance = KmerAnalyzer.KmerDistance("AAAA", "TTTT", 2);
            Assert.That(distance, Is.GreaterThan(0));
        }

        [Test]
        public void KmerDistance_SimilarSequences_SmallDistance()
        {
            double distance1 = KmerAnalyzer.KmerDistance("ACGT", "ACGA", 2);
            double distance2 = KmerAnalyzer.KmerDistance("ACGT", "TTTT", 2);

            Assert.That(distance1, Is.LessThan(distance2));
        }

        #endregion

        #region K-mers With Min Count

        [Test]
        public void FindKmersWithMinCount_FiltersCorrectly()
        {
            var kmers = KmerAnalyzer.FindKmersWithMinCount("ACGTACGT", 4, 2).ToList();

            Assert.That(kmers, Has.Count.EqualTo(1));
            Assert.That(kmers[0].Kmer, Is.EqualTo("ACGT"));
            Assert.That(kmers[0].Count, Is.EqualTo(2));
        }

        [Test]
        public void FindKmersWithMinCount_OrdersByCountDescending()
        {
            var kmers = KmerAnalyzer.FindKmersWithMinCount("AAAACGT", 2, 1).ToList();

            Assert.That(kmers[0].Count, Is.GreaterThanOrEqualTo(kmers[^1].Count));
        }

        #endregion

        #region Generate All K-mers

        [Test]
        public void GenerateAllKmers_Dna_GeneratesCorrectCount()
        {
            var kmers = KmerAnalyzer.GenerateAllKmers(2, "ACGT").ToList();

            // 4^2 = 16 possible 2-mers
            Assert.That(kmers, Has.Count.EqualTo(16));
        }

        [Test]
        public void GenerateAllKmers_K3_Generates64()
        {
            var kmers = KmerAnalyzer.GenerateAllKmers(3, "ACGT").ToList();

            // 4^3 = 64 possible 3-mers
            Assert.That(kmers, Has.Count.EqualTo(64));
        }

        [Test]
        public void GenerateAllKmers_InvalidK_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                KmerAnalyzer.GenerateAllKmers(0).ToList());
        }

        [Test]
        public void GenerateAllKmers_EmptyAlphabet_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() =>
                KmerAnalyzer.GenerateAllKmers(2, "").ToList());
        }

        #endregion

        #region Find K-mer Positions

        [Test]
        public void FindKmerPositions_MultipleOccurrences_ReturnsAll()
        {
            var positions = KmerAnalyzer.FindKmerPositions("ACGTACGT", "ACGT").ToList();

            Assert.That(positions, Has.Count.EqualTo(2));
            Assert.That(positions[0], Is.EqualTo(0));
            Assert.That(positions[1], Is.EqualTo(4));
        }

        [Test]
        public void FindKmerPositions_NotFound_ReturnsEmpty()
        {
            var positions = KmerAnalyzer.FindKmerPositions("ACGT", "TTTT").ToList();
            Assert.That(positions, Is.Empty);
        }

        #endregion

        // Note: CountKmersBothStrands tests moved to KmerAnalyzer_CountKmers_Tests.cs (KMER-COUNT-001)

        #region Analyze K-mers

        [Test]
        public void AnalyzeKmers_ReturnsCorrectStatistics()
        {
            var stats = KmerAnalyzer.AnalyzeKmers("ACGTACGT", 4);

            Assert.That(stats.TotalKmers, Is.EqualTo(5));  // 5 positions for 4-mers in 8-char sequence
            Assert.That(stats.UniqueKmers, Is.EqualTo(4)); // ACGT, CGTA, GTAC, TACG
            Assert.That(stats.MaxCount, Is.EqualTo(2));    // ACGT appears twice
            Assert.That(stats.MinCount, Is.EqualTo(1));    // Others appear once
        }

        [Test]
        public void AnalyzeKmers_EmptySequence_ReturnsZeros()
        {
            var stats = KmerAnalyzer.AnalyzeKmers("", 4);

            Assert.That(stats.TotalKmers, Is.EqualTo(0));
            Assert.That(stats.UniqueKmers, Is.EqualTo(0));
            Assert.That(stats.Entropy, Is.EqualTo(0));
        }

        #endregion

        #region Real-World Cases

        [Test]
        public void CountKmers_PromotorAnalysis()
        {
            // Looking for TATA box (TATAAAA) - simplified as TATA
            string promotor = "GCGCGCTATAAAAGGGGCTATAAAAATTT";
            var counts = KmerAnalyzer.CountKmers(promotor, 4);

            Assert.That(counts["TATA"], Is.EqualTo(2));
        }

        #endregion
    }
}
