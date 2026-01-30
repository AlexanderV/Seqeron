using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    [TestFixture]
    public class GenomicAnalyzerTests
    {
        #region Repeat Finding

        [Test]
        public void FindLongestRepeat_SimpleRepeat_FindsIt()
        {
            var dna = new DnaSequence("ACGTACGT");
            var repeat = GenomicAnalyzer.FindLongestRepeat(dna);

            Assert.That(repeat.Sequence, Is.EqualTo("ACGT"));
            Assert.That(repeat.Count, Is.EqualTo(2));
        }

        [Test]
        public void FindLongestRepeat_NoRepeat_ReturnsEmpty()
        {
            var dna = new DnaSequence("ACGT");
            var repeat = GenomicAnalyzer.FindLongestRepeat(dna);

            Assert.That(repeat.IsEmpty, Is.True);
        }

        // NOTE: FindTandemRepeats tests moved to GenomicAnalyzer_TandemRepeat_Tests.cs
        // as part of REP-TANDEM-001 consolidation with comprehensive coverage.

        [Test]
        public void FindRepeats_MultipleRepeats_FindsAll()
        {
            var dna = new DnaSequence("ACGTACGTTTTTACGT");
            var repeats = GenomicAnalyzer.FindRepeats(dna, minLength: 3).ToList();

            Assert.That(repeats.Any(r => r.Sequence == "ACGT"), Is.True);
        }

        #endregion

        #region Motif Finding

        [Test]
        public void FindMotif_ExistingMotif_FindsAllOccurrences()
        {
            var dna = new DnaSequence("ACGTACGTACGT");
            var positions = GenomicAnalyzer.FindMotif(dna, "ACGT");

            Assert.That(positions, Has.Count.EqualTo(3));
            Assert.That(positions, Does.Contain(0));
            Assert.That(positions, Does.Contain(4));
            Assert.That(positions, Does.Contain(8));
        }

        [Test]
        public void FindMotif_NonExistingMotif_ReturnsEmpty()
        {
            var dna = new DnaSequence("ACGTACGT");
            var positions = GenomicAnalyzer.FindMotif(dna, "TTTT");

            Assert.That(positions, Is.Empty);
        }

        [Test]
        public void FindMotif_CaseInsensitive_Works()
        {
            var dna = new DnaSequence("ACGTACGT");
            var positions = GenomicAnalyzer.FindMotif(dna, "acgt");

            Assert.That(positions, Has.Count.EqualTo(2));
        }

        // Note: Comprehensive Palindrome Detection tests are in RepeatFinder_Palindrome_Tests.cs
        // See Test Unit REP-PALIN-001
        // These are minimal smoke tests for GenomicAnalyzer.FindPalindromes API

        [Test]
        [Description("Smoke test: GenomicAnalyzer.FindPalindromes finds EcoRI site")]
        public void FindPalindromes_EcoRI_SmokeTest()
        {
            // EcoRI site: GAATTC (is a palindrome in DNA sense)
            var dna = new DnaSequence("AAAGAATTCAAA");
            var palindromes = GenomicAnalyzer.FindPalindromes(dna, minLength: 6, maxLength: 6).ToList();

            Assert.That(palindromes, Has.Count.EqualTo(1));
            Assert.That(palindromes[0].Sequence, Is.EqualTo("GAATTC"));
        }

        [Test]
        [Description("Smoke test: GenomicAnalyzer.FindPalindromes finds multiple sites")]
        public void FindPalindromes_MultipleSites_SmokeTest()
        {
            // Two restriction sites
            var dna = new DnaSequence("GAATTCAAAAGAATTC");
            var palindromes = GenomicAnalyzer.FindPalindromes(dna, minLength: 6, maxLength: 6).ToList();

            Assert.That(palindromes, Has.Count.EqualTo(2));
        }

        [Test]
        public void FindKnownMotifs_MultipleMotifs_FindsAll()
        {
            var dna = new DnaSequence("ATGAAATTTGGGCCC");
            var motifs = new[] { "ATG", "AAA", "GGG", "XXX" };

            var found = GenomicAnalyzer.FindKnownMotifs(dna, motifs);

            Assert.That(found.ContainsKey("ATG"), Is.True);
            Assert.That(found.ContainsKey("AAA"), Is.True);
            Assert.That(found.ContainsKey("GGG"), Is.True);
            Assert.That(found.ContainsKey("XXX"), Is.False);
        }

        #endregion

        #region Sequence Comparison

        [Test]
        public void FindLongestCommonRegion_CommonSubsequence_FindsIt()
        {
            var dna1 = new DnaSequence("AAACGTCGTAAA");
            var dna2 = new DnaSequence("TTTCGTCGTTTT");

            var common = GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);

            Assert.That(common.Sequence, Is.EqualTo("CGTCGT"));
        }

        [Test]
        public void FindLongestCommonRegion_NoCommon_ReturnsEmpty()
        {
            var dna1 = new DnaSequence("AAAA");
            var dna2 = new DnaSequence("TTTT");

            var common = GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);

            Assert.That(common.IsEmpty, Is.True);
        }

        [Test]
        public void CalculateSimilarity_IdenticalSequences_Returns100()
        {
            var dna1 = new DnaSequence("ACGTACGT");
            var dna2 = new DnaSequence("ACGTACGT");

            double similarity = GenomicAnalyzer.CalculateSimilarity(dna1, dna2);

            Assert.That(similarity, Is.EqualTo(100.0));
        }

        [Test]
        public void CalculateSimilarity_CompletelyDifferent_ReturnsLow()
        {
            var dna1 = new DnaSequence("AAAAAAAA");
            var dna2 = new DnaSequence("TTTTTTTT");

            double similarity = GenomicAnalyzer.CalculateSimilarity(dna1, dna2);

            Assert.That(similarity, Is.LessThan(50.0));
        }

        #endregion

        #region Open Reading Frames

        [Test]
        public void FindOpenReadingFrames_SimpleOrf_FindsIt()
        {
            // ATG (start) + 30 codons + TAA (stop) = 99bp ORF
            string seq = "ATG" + new string('A', 93) + "TAA";
            var dna = new DnaSequence(seq);

            var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 90).ToList();

            Assert.That(orfs, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(orfs.Any(o => o.Sequence.StartsWith("ATG")), Is.True);
        }

        [Test]
        public void FindOpenReadingFrames_MultipleFrames_FindsAll()
        {
            // ORF in frame 1
            string seq = "ATGAAAAAATAA" + "ATGCCCCCCTAGATGGGGGGGTGA";
            var dna = new DnaSequence(seq);

            var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 9).ToList();

            // Should find ORFs starting with ATG
            Assert.That(orfs.All(o => o.Sequence.StartsWith("ATG")), Is.True);
        }

        [Test]
        public void FindOpenReadingFrames_NoOrf_ReturnsEmpty()
        {
            var dna = new DnaSequence("TTTTTTTTT");

            var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 30).ToList();

            Assert.That(orfs, Is.Empty);
        }

        #endregion
    }
}
