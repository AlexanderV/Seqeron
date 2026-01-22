using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    /// <summary>
    /// Test suite for ApproximateMatcher methods NOT covered by PAT-APPROX-001.
    /// 
    /// PAT-APPROX-001 (Hamming Distance) tests are in: ApproximateMatcher_HammingDistance_Tests.cs
    /// PAT-APPROX-002 (Edit Distance) tests are in this file (Edit Distance region)
    /// 
    /// This file also contains tests for utility methods:
    /// - FindBestMatch
    /// - CountApproximateOccurrences
    /// - FindFrequentKmersWithMismatches
    /// </summary>
    [TestFixture]
    public class ApproximateMatcherTests
    {
        #region Edit Distance (PAT-APPROX-002 scope)

        [Test]
        public void EditDistance_IdenticalStrings_ReturnsZero()
        {
            int distance = ApproximateMatcher.EditDistance("ACGT", "ACGT");
            Assert.That(distance, Is.EqualTo(0));
        }

        [Test]
        public void EditDistance_OneSubstitution_ReturnsOne()
        {
            int distance = ApproximateMatcher.EditDistance("ACGT", "ACGG");
            Assert.That(distance, Is.EqualTo(1));
        }

        [Test]
        public void EditDistance_OneInsertion_ReturnsOne()
        {
            int distance = ApproximateMatcher.EditDistance("ACGT", "ACGGT");
            Assert.That(distance, Is.EqualTo(1));
        }

        [Test]
        public void EditDistance_OneDeletion_ReturnsOne()
        {
            int distance = ApproximateMatcher.EditDistance("ACGT", "ACT");
            Assert.That(distance, Is.EqualTo(1));
        }

        [Test]
        public void EditDistance_EmptyAndNonEmpty_ReturnsLength()
        {
            Assert.That(ApproximateMatcher.EditDistance("", "ACGT"), Is.EqualTo(4));
            Assert.That(ApproximateMatcher.EditDistance("ACGT", ""), Is.EqualTo(4));
        }

        [Test]
        public void EditDistance_ComplexCase_CalculatesCorrectly()
        {
            // GATTACA → GCATGCU
            // Requires: G→G, A→C, T→A, T→T, A→G, C→C, insert U
            int distance = ApproximateMatcher.EditDistance("GATTACA", "GCATGCU");
            Assert.That(distance, Is.EqualTo(4));
        }

        [Test]
        public void EditDistance_CaseInsensitive()
        {
            int distance = ApproximateMatcher.EditDistance("acgt", "ACGT");
            Assert.That(distance, Is.EqualTo(0));
        }

        #endregion

        #region Find With Edits (PAT-APPROX-002 scope)

        [Test]
        public void FindWithEdits_ExactMatch_Found()
        {
            var matches = ApproximateMatcher.FindWithEdits("ACGTACGT", "ACGT", 0).ToList();

            Assert.That(matches, Has.Count.EqualTo(2));
        }

        [Test]
        public void FindWithEdits_WithInsertion_Found()
        {
            // Looking for ACT in ACGT (G is inserted)
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACT", 1).ToList();

            Assert.That(matches.Any(m => m.Distance == 1), Is.True);
        }

        [Test]
        public void FindWithEdits_WithDeletion_Found()
        {
            // Looking for ACGGT in ACGT (G is deleted in sequence)
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACG", 1).ToList();

            Assert.That(matches.Any(m => m.Distance == 0), Is.True);
        }

        #endregion

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
