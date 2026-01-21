using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    [TestFixture]
    public class DnaSequenceTests
    {
        #region Construction

        [Test]
        public void Constructor_ValidSequence_CreatesSequence()
        {
            var dna = new DnaSequence("ACGT");
            Assert.That(dna.Sequence, Is.EqualTo("ACGT"));
        }

        [Test]
        public void Constructor_LowercaseSequence_NormalizesToUppercase()
        {
            var dna = new DnaSequence("acgt");
            Assert.That(dna.Sequence, Is.EqualTo("ACGT"));
        }

        [Test]
        public void Constructor_MixedCase_NormalizesToUppercase()
        {
            var dna = new DnaSequence("AcGt");
            Assert.That(dna.Sequence, Is.EqualTo("ACGT"));
        }

        [Test]
        public void Constructor_InvalidNucleotide_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() => new DnaSequence("ACGX"));
        }

        [Test]
        public void Constructor_EmptySequence_CreatesEmpty()
        {
            var dna = new DnaSequence("");
            Assert.That(dna.Sequence, Is.Empty);
            Assert.That(dna.Length, Is.EqualTo(0));
        }

        [Test]
        public void TryCreate_ValidSequence_ReturnsTrue()
        {
            bool result = DnaSequence.TryCreate("ACGT", out var dna);
            Assert.That(result, Is.True);
            Assert.That(dna!.Sequence, Is.EqualTo("ACGT"));
        }

        [Test]
        public void TryCreate_InvalidSequence_ReturnsFalse()
        {
            bool result = DnaSequence.TryCreate("ACGX", out var dna);
            Assert.That(result, Is.False);
            Assert.That(dna, Is.Null);
        }

        #endregion

        #region Complement

        [Test]
        public void Complement_ReturnsCorrectComplement()
        {
            var dna = new DnaSequence("ACGT");
            var complement = dna.Complement();
            Assert.That(complement.Sequence, Is.EqualTo("TGCA"));
        }

        [Test]
        public void Complement_LongerSequence_ReturnsCorrectComplement()
        {
            var dna = new DnaSequence("AATTCCGG");
            var complement = dna.Complement();
            Assert.That(complement.Sequence, Is.EqualTo("TTAAGGCC"));
        }

        #endregion

        #region Reverse Complement

        [Test]
        public void ReverseComplement_ReturnsCorrectReverseComplement()
        {
            var dna = new DnaSequence("ACGT");
            var revComp = dna.ReverseComplement();
            Assert.That(revComp.Sequence, Is.EqualTo("ACGT")); // ACGT is its own reverse complement!
        }

        [Test]
        public void ReverseComplement_AsymmetricSequence_Works()
        {
            var dna = new DnaSequence("AACG");
            var revComp = dna.ReverseComplement();
            Assert.That(revComp.Sequence, Is.EqualTo("CGTT"));
        }

        [Test]
        public void ReverseComplement_EcoRI_Site()
        {
            // EcoRI recognition site: GAATTC
            var dna = new DnaSequence("GAATTC");
            var revComp = dna.ReverseComplement();
            Assert.That(revComp.Sequence, Is.EqualTo("GAATTC")); // Palindrome!
        }

        #endregion

        #region GC Content

        [Test]
        public void GcContent_AllGC_Returns100()
        {
            var dna = new DnaSequence("GCGCGC");
            Assert.That(dna.GcContent(), Is.EqualTo(100.0));
        }

        [Test]
        public void GcContent_NoGC_Returns0()
        {
            var dna = new DnaSequence("ATATAT");
            Assert.That(dna.GcContent(), Is.EqualTo(0.0));
        }

        [Test]
        public void GcContent_HalfGC_Returns50()
        {
            var dna = new DnaSequence("ACGT");
            Assert.That(dna.GcContent(), Is.EqualTo(50.0));
        }

        [Test]
        public void GcContent_Empty_Returns0()
        {
            var dna = new DnaSequence("");
            Assert.That(dna.GcContent(), Is.EqualTo(0.0));
        }

        #endregion

        #region Transcription

        [Test]
        public void Transcribe_ReplacesThymineWithUracil()
        {
            var dna = new DnaSequence("ATGC");
            string rna = dna.Transcribe();
            Assert.That(rna, Is.EqualTo("AUGC"));
        }

        [Test]
        public void Transcribe_NoThymine_UnchangedExceptT()
        {
            var dna = new DnaSequence("ACGACG");
            string rna = dna.Transcribe();
            Assert.That(rna, Is.EqualTo("ACGACG"));
        }

        #endregion

        #region Subsequence

        [Test]
        public void Subsequence_ReturnsCorrectSubsequence()
        {
            var dna = new DnaSequence("ACGTACGT");
            var sub = dna.Subsequence(2, 4);
            Assert.That(sub.Sequence, Is.EqualTo("GTAC"));
        }

        [Test]
        public void Indexer_ReturnsCorrectNucleotide()
        {
            var dna = new DnaSequence("ACGT");
            Assert.That(dna[0], Is.EqualTo('A'));
            Assert.That(dna[1], Is.EqualTo('C'));
            Assert.That(dna[2], Is.EqualTo('G'));
            Assert.That(dna[3], Is.EqualTo('T'));
        }

        #endregion

        #region Equality

        [Test]
        public void Equals_SameSequence_ReturnsTrue()
        {
            var dna1 = new DnaSequence("ACGT");
            var dna2 = new DnaSequence("ACGT");
            Assert.That(dna1.Equals(dna2), Is.True);
        }

        [Test]
        public void Equals_DifferentSequence_ReturnsFalse()
        {
            var dna1 = new DnaSequence("ACGT");
            var dna2 = new DnaSequence("TGCA");
            Assert.That(dna1.Equals(dna2), Is.False);
        }

        #endregion
    }
}
