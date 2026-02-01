using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    [TestFixture]
    public class RnaSequenceTests
    {
        #region Construction

        [Test]
        public void Constructor_ValidSequence_CreatesSequence()
        {
            var rna = new RnaSequence("ACGU");
            Assert.That(rna.Sequence, Is.EqualTo("ACGU"));
        }

        [Test]
        public void Constructor_LowercaseSequence_NormalizesToUppercase()
        {
            var rna = new RnaSequence("acgu");
            Assert.That(rna.Sequence, Is.EqualTo("ACGU"));
        }

        [Test]
        public void Constructor_MixedCase_NormalizesToUppercase()
        {
            var rna = new RnaSequence("AcGu");
            Assert.That(rna.Sequence, Is.EqualTo("ACGU"));
        }

        [Test]
        public void Constructor_InvalidNucleotide_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() => new RnaSequence("ACGX"));
        }

        [Test]
        public void Constructor_ThymineInvalid_ThrowsArgumentException()
        {
            // T is DNA, not RNA
            Assert.Throws<ArgumentException>(() => new RnaSequence("ACGT"));
        }

        [Test]
        public void Constructor_EmptySequence_CreatesEmpty()
        {
            var rna = new RnaSequence("");
            Assert.That(rna.Sequence, Is.Empty);
            Assert.That(rna.Length, Is.EqualTo(0));
        }

        [Test]
        public void TryCreate_ValidSequence_ReturnsTrue()
        {
            bool result = RnaSequence.TryCreate("ACGU", out var rna);
            Assert.That(result, Is.True);
            Assert.That(rna!.Sequence, Is.EqualTo("ACGU"));
        }

        [Test]
        public void TryCreate_InvalidSequence_ReturnsFalse()
        {
            bool result = RnaSequence.TryCreate("ACGX", out var rna);
            Assert.That(result, Is.False);
            Assert.That(rna, Is.Null);
        }

        #endregion

        #region Complement

        [Test]
        [Description("Smoke test: RnaSequence.Complement() delegates to SequenceExtensions.GetRnaComplementBase")]
        public void Complement_DelegatesToGetRnaComplementBase_SmokeTest()
        {
            var rna = new RnaSequence("ACGU");
            var complement = rna.Complement();
            Assert.That(complement.Sequence, Is.EqualTo("UGCA"));
        }

        // Comprehensive RNA complement tests are in SequenceExtensions_Complement_Tests.cs

        #endregion

        #region Reverse Complement

        [Test]
        public void ReverseComplement_ReturnsCorrectReverseComplement()
        {
            var rna = new RnaSequence("ACGU");
            var revComp = rna.ReverseComplement();
            Assert.That(revComp.Sequence, Is.EqualTo("ACGU")); // ACGU is its own reverse complement!
        }

        [Test]
        public void ReverseComplement_AsymmetricSequence_Works()
        {
            var rna = new RnaSequence("AACG");
            var revComp = rna.ReverseComplement();
            Assert.That(revComp.Sequence, Is.EqualTo("CGUU"));
        }

        #endregion

        // Note: GC Content detailed tests are in SequenceExtensions_CalculateGcContent_Tests.cs
        // RnaSequence.GcContent() delegates to SequenceExtensions.CalculateGcContentFast()

        #region AU Content

        [Test]
        public void AuContent_AllAU_Returns100()
        {
            var rna = new RnaSequence("AUAUAU");
            Assert.That(rna.AuContent(), Is.EqualTo(100.0));
        }

        [Test]
        public void AuContent_NoAU_Returns0()
        {
            var rna = new RnaSequence("GCGCGC");
            Assert.That(rna.AuContent(), Is.EqualTo(0.0));
        }

        [Test]
        public void AuContent_HalfAU_Returns50()
        {
            var rna = new RnaSequence("ACGU");
            Assert.That(rna.AuContent(), Is.EqualTo(50.0));
        }

        #endregion

        #region Reverse Transcription

        [Test]
        public void ReverseTranscribe_ReplacesUracilWithThymine()
        {
            var rna = new RnaSequence("AUGC");
            DnaSequence dna = rna.ReverseTranscribe();
            Assert.That(dna.Sequence, Is.EqualTo("ATGC"));
        }

        [Test]
        public void ReverseTranscribe_AllNucleotides_Works()
        {
            var rna = new RnaSequence("ACGU");
            DnaSequence dna = rna.ReverseTranscribe();
            Assert.That(dna.Sequence, Is.EqualTo("ACGT"));
        }

        #endregion

        #region From DNA

        [Test]
        public void FromDna_TranscribesDnaToRna()
        {
            var dna = new DnaSequence("ATGC");
            var rna = RnaSequence.FromDna(dna);
            Assert.That(rna.Sequence, Is.EqualTo("AUGC"));
        }

        #endregion

        #region Codons

        [Test]
        public void GetCodons_Frame0_ReturnsCodons()
        {
            var rna = new RnaSequence("AUGGCUUAA");
            var codons = rna.GetCodons(0).ToList();
            Assert.That(codons, Has.Count.EqualTo(3));
            Assert.That(codons[0], Is.EqualTo("AUG"));
            Assert.That(codons[1], Is.EqualTo("GCU"));
            Assert.That(codons[2], Is.EqualTo("UAA"));
        }

        [Test]
        public void GetCodons_Frame1_ReturnsCodons()
        {
            var rna = new RnaSequence("GAUGGCUUAA");
            var codons = rna.GetCodons(1).ToList();
            Assert.That(codons, Has.Count.EqualTo(3));
            Assert.That(codons[0], Is.EqualTo("AUG"));
            Assert.That(codons[1], Is.EqualTo("GCU"));
            Assert.That(codons[2], Is.EqualTo("UAA"));
        }

        [Test]
        public void GetCodons_InvalidFrame_ThrowsException()
        {
            var rna = new RnaSequence("AUGGCUUAA");
            Assert.Throws<ArgumentOutOfRangeException>(() => rna.GetCodons(3).ToList());
        }

        #endregion

        #region Subsequence

        [Test]
        public void Subsequence_ReturnsCorrectSubsequence()
        {
            var rna = new RnaSequence("ACGUACGU");
            var sub = rna.Subsequence(2, 4);
            Assert.That(sub.Sequence, Is.EqualTo("GUAC"));
        }

        #endregion

        #region Indexer

        [Test]
        public void Indexer_ReturnsCorrectNucleotide()
        {
            var rna = new RnaSequence("ACGU");
            Assert.That(rna[0], Is.EqualTo('A'));
            Assert.That(rna[1], Is.EqualTo('C'));
            Assert.That(rna[2], Is.EqualTo('G'));
            Assert.That(rna[3], Is.EqualTo('U'));
        }

        #endregion

        #region Equality

        [Test]
        public void Equals_SameSequence_ReturnsTrue()
        {
            var rna1 = new RnaSequence("ACGU");
            var rna2 = new RnaSequence("ACGU");
            Assert.That(rna1.Equals(rna2), Is.True);
        }

        [Test]
        public void Equals_DifferentSequence_ReturnsFalse()
        {
            var rna1 = new RnaSequence("ACGU");
            var rna2 = new RnaSequence("UGCA");
            Assert.That(rna1.Equals(rna2), Is.False);
        }

        [Test]
        public void GetHashCode_SameSequence_SameHash()
        {
            var rna1 = new RnaSequence("ACGU");
            var rna2 = new RnaSequence("ACGU");
            Assert.That(rna1.GetHashCode(), Is.EqualTo(rna2.GetHashCode()));
        }

        #endregion

        #region ToString

        [Test]
        public void ToString_ReturnsSequence()
        {
            var rna = new RnaSequence("ACGU");
            Assert.That(rna.ToString(), Is.EqualTo("ACGU"));
        }

        #endregion

        #region SuffixTree Integration

        [Test]
        public void SuffixTree_IsCreatedLazily()
        {
            var rna = new RnaSequence("ACGUACGU");
            var tree = rna.SuffixTree;
            Assert.That(tree, Is.Not.Null);
            Assert.That(tree.Contains("ACGU"), Is.True);
            Assert.That(tree.Contains("GUAC"), Is.True);
        }

        #endregion
    }
}
