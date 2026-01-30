using NUnit.Framework;

namespace Seqeron.Genomics.Tests
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

        // Note: Comprehensive validation tests are in SequenceExtensions_SequenceValidation_Tests.cs (SEQ-VALID-001)
        // Factory/constructor tests here verify delegation to validation logic

        [Test]
        [Description("Smoke test: Constructor throws on invalid nucleotide")]
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
        [Description("Smoke test: TryCreate returns true for valid sequence")]
        public void TryCreate_ValidSequence_ReturnsTrue()
        {
            bool result = DnaSequence.TryCreate("ACGT", out var dna);
            Assert.That(result, Is.True);
            Assert.That(dna!.Sequence, Is.EqualTo("ACGT"));
        }

        [Test]
        [Description("Smoke test: TryCreate returns false for invalid sequence")]
        public void TryCreate_InvalidSequence_ReturnsFalse()
        {
            bool result = DnaSequence.TryCreate("ACGX", out var dna);
            Assert.That(result, Is.False);
            Assert.That(dna, Is.Null);
        }

        #endregion

        #region Complement

        [Test]
        [Description("Smoke test: DnaSequence.Complement() delegates to SequenceExtensions.GetComplementBase")]
        public void Complement_DelegatesToGetComplementBase_SmokeTest()
        {
            var dna = new DnaSequence("ACGT");
            var complement = dna.Complement();
            Assert.That(complement.Sequence, Is.EqualTo("TGCA"));
        }

        // Comprehensive complement tests are in SequenceExtensions_Complement_Tests.cs

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

        // Note: GC Content detailed tests are in SequenceExtensions_CalculateGcContent_Tests.cs
        // DnaSequence.GcContent() delegates to SequenceExtensions.CalculateGcContentFast()

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
