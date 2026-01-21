using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    [TestFixture]
    public class KmerAnalyzerTests
    {
        #region Count K-mers

        [Test]
        public void CountKmers_SimpleSequence_CountsCorrectly()
        {
            var counts = KmerAnalyzer.CountKmers("ACGTACGT", 4);

            Assert.That(counts["ACGT"], Is.EqualTo(2));
            Assert.That(counts["CGTA"], Is.EqualTo(1));
            Assert.That(counts["GTAC"], Is.EqualTo(1));
            Assert.That(counts["TACG"], Is.EqualTo(1));
        }

        [Test]
        public void CountKmers_AllSame_SingleCount()
        {
            var counts = KmerAnalyzer.CountKmers("AAAA", 2);

            Assert.That(counts.Count, Is.EqualTo(1));
            Assert.That(counts["AA"], Is.EqualTo(3));
        }

        [Test]
        public void CountKmers_EmptySequence_ReturnsEmptyDictionary()
        {
            var counts = KmerAnalyzer.CountKmers("", 4);
            Assert.That(counts, Is.Empty);
        }

        [Test]
        public void CountKmers_KLargerThanSequence_ReturnsEmptyDictionary()
        {
            var counts = KmerAnalyzer.CountKmers("ACG", 4);
            Assert.That(counts, Is.Empty);
        }

        [Test]
        public void CountKmers_InvalidK_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmers("ACGT", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmers("ACGT", -1));
        }

        [Test]
        public void CountKmers_DnaSequence_Works()
        {
            var dna = new DnaSequence("ACGTACGT");
            var counts = KmerAnalyzer.CountKmers(dna, 4);

            Assert.That(counts["ACGT"], Is.EqualTo(2));
        }

        #endregion

        #region K-mer Spectrum

        [Test]
        public void GetKmerSpectrum_ReturnsFrequencyDistribution()
        {
            // ACGTACGT: ACGT appears 2 times, others appear 1 time
            var spectrum = KmerAnalyzer.GetKmerSpectrum("ACGTACGT", 4);

            Assert.That(spectrum[1], Is.EqualTo(3)); // 3 k-mers appear once
            Assert.That(spectrum[2], Is.EqualTo(1)); // 1 k-mer appears twice
        }

        #endregion

        #region Most Frequent K-mers

        [Test]
        public void FindMostFrequentKmers_SingleMostFrequent_ReturnsIt()
        {
            var mostFrequent = KmerAnalyzer.FindMostFrequentKmers("AAACGT", 2).ToList();

            Assert.That(mostFrequent, Does.Contain("AA"));
        }

        [Test]
        public void FindMostFrequentKmers_MultipleMostFrequent_ReturnsAll()
        {
            // In ACGTACGT, ACGT appears twice, but so do multiple overlapping sequences
            var mostFrequent = KmerAnalyzer.FindMostFrequentKmers("ACGTACGT", 4).ToList();

            Assert.That(mostFrequent, Does.Contain("ACGT"));
        }

        [Test]
        public void FindMostFrequentKmers_EmptySequence_ReturnsEmpty()
        {
            var mostFrequent = KmerAnalyzer.FindMostFrequentKmers("", 4).ToList();
            Assert.That(mostFrequent, Is.Empty);
        }

        #endregion

        #region K-mer Frequencies

        [Test]
        public void GetKmerFrequencies_SumsToOne()
        {
            var frequencies = KmerAnalyzer.GetKmerFrequencies("ACGTACGT", 2);
            double total = frequencies.Values.Sum();

            Assert.That(total, Is.EqualTo(1.0).Within(0.0001));
        }

        [Test]
        public void GetKmerFrequencies_CalculatesCorrectly()
        {
            // AAA: AA appears twice out of 2 possible 2-mers
            var frequencies = KmerAnalyzer.GetKmerFrequencies("AAA", 2);

            Assert.That(frequencies["AA"], Is.EqualTo(1.0));
        }

        #endregion

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

        #region Unique K-mers

        [Test]
        public void FindUniqueKmers_ReturnsKmersAppearingOnce()
        {
            // In ACGTACGT: ACGT appears twice, others appear once
            var unique = KmerAnalyzer.FindUniqueKmers("ACGTACGT", 4).ToList();

            Assert.That(unique, Does.Contain("CGTA"));
            Assert.That(unique, Does.Contain("GTAC"));
            Assert.That(unique, Does.Contain("TACG"));
            Assert.That(unique, Does.Not.Contain("ACGT"));
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

        #region K-mer Entropy

        [Test]
        public void CalculateKmerEntropy_UniformDistribution_HighEntropy()
        {
            // All unique k-mers = high entropy
            double entropy = KmerAnalyzer.CalculateKmerEntropy("ACGT", 1);

            // Maximum entropy for 4 symbols = log2(4) = 2
            Assert.That(entropy, Is.EqualTo(2.0).Within(0.0001));
        }

        [Test]
        public void CalculateKmerEntropy_SingleRepeated_ZeroEntropy()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("AAAA", 2);

            // Only one k-mer, so entropy is 0
            Assert.That(entropy, Is.EqualTo(0.0));
        }

        [Test]
        public void CalculateKmerEntropy_EmptySequence_ReturnsZero()
        {
            double entropy = KmerAnalyzer.CalculateKmerEntropy("", 2);
            Assert.That(entropy, Is.EqualTo(0.0));
        }

        #endregion

        #region Find Clumps

        [Test]
        public void FindClumps_SimpleClump_Found()
        {
            // AAAAA has AAA appearing 3 times in a window of 5
            var clumps = KmerAnalyzer.FindClumps("AAAAA", 3, 5, 3).ToList();

            Assert.That(clumps, Does.Contain("AAA"));
        }

        [Test]
        public void FindClumps_NoClump_ReturnsEmpty()
        {
            // Looking for 3 occurrences in window of 4, but sequence too short
            var clumps = KmerAnalyzer.FindClumps("ACGT", 2, 4, 3).ToList();

            Assert.That(clumps, Is.Empty);
        }

        [Test]
        public void FindClumps_InvalidParameters_ReturnsEmpty()
        {
            var clumps = KmerAnalyzer.FindClumps("ACGT", 5, 4, 2).ToList(); // k > window
            Assert.That(clumps, Is.Empty);
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

        #region Count K-mers Both Strands

        [Test]
        public void CountKmersBothStrands_CountsBothStrands()
        {
            var dna = new DnaSequence("ACGT");
            var counts = KmerAnalyzer.CountKmersBothStrands(dna, 2);

            // Forward: AC, CG, GT
            // Reverse (ACGT is palindrome): AC, CG, GT
            // Each appears twice total
            Assert.That(counts["AC"], Is.EqualTo(2));
            Assert.That(counts["CG"], Is.EqualTo(2));
            Assert.That(counts["GT"], Is.EqualTo(2));
        }

        #endregion

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

        [Test]
        public void FindClumps_RepetitiveRegion()
        {
            // Simulating a region with repetitive elements
            string sequence = "AAAAACGTAAAAACGT";
            var clumps = KmerAnalyzer.FindClumps(sequence, 4, 8, 2).ToList();

            // AAAA should form a clump
            Assert.That(clumps, Does.Contain("AAAA"));
        }

        #endregion
    }
}
