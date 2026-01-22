using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    /// <summary>
    /// Test suite for ApproximateMatcher utility methods.
    /// 
    /// PAT-APPROX-001 (Hamming Distance) tests: ApproximateMatcher_HammingDistance_Tests.cs
    /// PAT-APPROX-002 (Edit Distance) tests: ApproximateMatcher_EditDistance_Tests.cs
    /// 
    /// This file contains tests for utility methods:
    /// - FindBestMatch
    /// - CountApproximateOccurrences
    /// - FindFrequentKmersWithMismatches
    /// </summary>
    [TestFixture]
    public class ApproximateMatcherTests
    {
        #region Find Best Match (uses Hamming internally - utility method)

        [Test]
        public void FindBestMatch_ExactMatch_ReturnsZeroDistance()
        {
            var result = ApproximateMatcher.FindBestMatch("ACGTACGT", "ACGT");

            Assert.That(result, Is.Not.Null);
            Assert.That(result!.Value.Distance, Is.EqualTo(0));
            Assert.That(result!.Value.IsExact, Is.True);
        }

        [Test]
        public void FindBestMatch_NoExactMatch_ReturnsBest()
        {
            // TTTTTTTT vs ACGT: best match is TTTT with 3 mismatches (A→T, C→T, G→T)
            var result = ApproximateMatcher.FindBestMatch("TTTTTTTT", "ACGT");

            Assert.That(result, Is.Not.Null);
            Assert.That(result!.Value.Distance, Is.EqualTo(3)); // T matches T at position 3
        }

        [Test]
        public void FindBestMatch_EmptySequence_ReturnsNull()
        {
            var result = ApproximateMatcher.FindBestMatch("", "ACGT");
            Assert.That(result, Is.Null);
        }

        [Test]
        public void FindBestMatch_PatternTooLong_ReturnsNull()
        {
            var result = ApproximateMatcher.FindBestMatch("AC", "ACGT");
            Assert.That(result, Is.Null);
        }

        #endregion

        #region Count Approximate Occurrences (delegates to FindWithMismatches - smoke tests)

        [Test]
        [Description("Smoke test: CountApproximateOccurrences delegates to FindWithMismatches")]
        public void CountApproximateOccurrences_ExactMatches_CountsCorrectly()
        {
            int count = ApproximateMatcher.CountApproximateOccurrences("ACGTACGT", "ACGT", 0);
            Assert.That(count, Is.EqualTo(2));
        }

        [Test]
        public void CountApproximateOccurrences_WithMismatches_CountsAll()
        {
            // ACGT appears exactly, ACGG would match with 1 mismatch
            int count = ApproximateMatcher.CountApproximateOccurrences("ACGTACGT", "ACGG", 1);
            Assert.That(count, Is.EqualTo(2));
        }

        #endregion

        #region Frequent K-mers with Mismatches

        [Test]
        public void FindFrequentKmersWithMismatches_SimpleCase_FindsMostFrequent()
        {
            // In ACGT, 2-mers are AC, CG, GT
            // With 1 mismatch, AC matches CC, AA, AG, TC, etc.
            var result = ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 2, 0).ToList();

            Assert.That(result, Has.Count.EqualTo(3)); // AC, CG, GT each appear once
        }

        [Test]
        public void FindFrequentKmersWithMismatches_RepeatSequence_FindsRepeated()
        {
            // AAAAAA has 3 occurrences of AAAA
            var result = ApproximateMatcher.FindFrequentKmersWithMismatches("AAAAAA", 4, 0).ToList();

            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Kmer, Is.EqualTo("AAAA"));
            Assert.That(result[0].Count, Is.EqualTo(3));
        }

        [Test]
        public void FindFrequentKmersWithMismatches_InvalidK_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 0, 1).ToList());
        }

        [Test]
        public void FindFrequentKmersWithMismatches_NegativeD_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 2, -1).ToList());
        }

        #endregion

        // Note: Real-world use case tests for FindWithMismatches are in
        // ApproximateMatcher_HammingDistance_Tests.cs (PAT-APPROX-001)
    }
}
