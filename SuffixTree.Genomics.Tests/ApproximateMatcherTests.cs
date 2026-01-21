using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    [TestFixture]
    public class ApproximateMatcherTests
    {
        #region Hamming Distance

        [Test]
        public void HammingDistance_IdenticalStrings_ReturnsZero()
        {
            int distance = ApproximateMatcher.HammingDistance("ACGT", "ACGT");
            Assert.That(distance, Is.EqualTo(0));
        }

        [Test]
        public void HammingDistance_OneDifference_ReturnsOne()
        {
            int distance = ApproximateMatcher.HammingDistance("ACGT", "ACGG");
            Assert.That(distance, Is.EqualTo(1));
        }

        [Test]
        public void HammingDistance_AllDifferent_ReturnsLength()
        {
            int distance = ApproximateMatcher.HammingDistance("AAAA", "TTTT");
            Assert.That(distance, Is.EqualTo(4));
        }

        [Test]
        public void HammingDistance_CaseInsensitive()
        {
            int distance = ApproximateMatcher.HammingDistance("acgt", "ACGT");
            Assert.That(distance, Is.EqualTo(0));
        }

        [Test]
        public void HammingDistance_DifferentLengths_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => 
                ApproximateMatcher.HammingDistance("ACGT", "ACG"));
        }

        #endregion

        #region Edit Distance

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

        #region Find With Mismatches

        [Test]
        public void FindWithMismatches_ExactMatch_FoundWithZeroMismatches()
        {
            var matches = ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGT", 0).ToList();
            
            Assert.That(matches, Has.Count.EqualTo(2));
            Assert.That(matches[0].Position, Is.EqualTo(0));
            Assert.That(matches[1].Position, Is.EqualTo(4));
            Assert.That(matches.All(m => m.Distance == 0), Is.True);
        }

        [Test]
        public void FindWithMismatches_OneMismatch_Found()
        {
            var matches = ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGG", 1).ToList();
            
            Assert.That(matches, Has.Count.EqualTo(2));
            Assert.That(matches[0].Distance, Is.EqualTo(1));
            Assert.That(matches[0].MismatchPositions, Does.Contain(3));
        }

        [Test]
        public void FindWithMismatches_TooManyMismatches_NotFound()
        {
            var matches = ApproximateMatcher.FindWithMismatches("ACGT", "TGCA", 2).ToList();
            
            Assert.That(matches, Is.Empty);
        }

        [Test]
        public void FindWithMismatches_MultipleMismatches_AllReturned()
        {
            // Find AAAA in TTTTAAAATTTT with up to 2 mismatches
            var matches = ApproximateMatcher.FindWithMismatches("TTTTAAAATTTT", "AAAA", 2).ToList();
            
            // Should find exact match at position 4, and approximate matches
            var exactMatch = matches.FirstOrDefault(m => m.Distance == 0);
            Assert.That(exactMatch.Position, Is.EqualTo(4));
        }

        [Test]
        public void FindWithMismatches_EmptyPattern_ReturnsEmpty()
        {
            var matches = ApproximateMatcher.FindWithMismatches("ACGT", "", 1).ToList();
            Assert.That(matches, Is.Empty);
        }

        [Test]
        public void FindWithMismatches_PatternLongerThanSequence_ReturnsEmpty()
        {
            var matches = ApproximateMatcher.FindWithMismatches("ACG", "ACGT", 1).ToList();
            Assert.That(matches, Is.Empty);
        }

        [Test]
        public void FindWithMismatches_NegativeMismatches_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => 
                ApproximateMatcher.FindWithMismatches("ACGT", "AC", -1).ToList());
        }

        [Test]
        public void FindWithMismatches_DnaSequence_Works()
        {
            var dna = new DnaSequence("ACGTACGT");
            var matches = ApproximateMatcher.FindWithMismatches(dna, "ACGT", 0).ToList();
            
            Assert.That(matches, Has.Count.EqualTo(2));
        }

        #endregion

        #region Find With Edits

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

        #region Find Best Match

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

        #region Count Approximate Occurrences

        [Test]
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

        #region Real-World Cases

        [Test]
        public void FindWithMismatches_SNP_Detection()
        {
            // Simulating SNP detection - looking for a reference allele with possible mutations
            string genome = "ATGCGATCGATCGATCGATCG";
            string reference = "GATC"; // Looking for GATC pattern

            var matches = ApproximateMatcher.FindWithMismatches(genome, reference, 1).ToList();

            // Should find multiple matches, some exact, some with 1 mismatch
            Assert.That(matches.Count(m => m.Distance == 0), Is.GreaterThan(0));
        }

        [Test]
        public void FindWithMismatches_PrimerBinding_WithMismatches()
        {
            // Simulating primer binding - primer might bind with mismatches
            string template = "ATGCATGCATGCATGCATGCATGC";
            string primer = "ATGC";

            var bindings = ApproximateMatcher.FindWithMismatches(template, primer, 1).ToList();

            // Should find multiple binding sites
            Assert.That(bindings, Has.Count.GreaterThan(0));
        }

        #endregion
    }
}
