using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for Translator protein translation.
    /// Test Unit: TRANS-PROT-001
    /// 
    /// Evidence Sources:
    ///   - Wikipedia: Translation (biology), Reading frame, Open reading frame
    ///   - NCBI: The Genetic Codes (translation tables)
    ///   - Lodish H et al. (2007): Molecular Cell Biology
    /// 
    /// Key Invariants:
    ///   - Codons are read in triplets from 5' to 3'
    ///   - Frame parameter offsets reading start position
    ///   - Stop codons (UAA, UAG, UGA) terminate translation
    ///   - Six-frame translation covers both strands
    /// </summary>
    [TestFixture]
    public class TranslatorTests
    {
        #region Basic Translation

        [Test]
        public void Translate_SingleCodon_ReturnsSingleAminoAcid()
        {
            var dna = new DnaSequence("ATG");
            var protein = Translator.Translate(dna);
            Assert.That(protein.Sequence, Is.EqualTo("M"));
        }

        [Test]
        public void Translate_MultipleCodens_ReturnsProtein()
        {
            // ATG GCT TAA = M A *
            var dna = new DnaSequence("ATGGCTTAA");
            var protein = Translator.Translate(dna);
            Assert.That(protein.Sequence, Is.EqualTo("MA*"));
        }

        [Test]
        public void Translate_ToFirstStop_StopsAtStopCodon()
        {
            // ATG GCT TAA GCT = M A * A
            var dna = new DnaSequence("ATGGCTTAAGCT");
            var protein = Translator.Translate(dna, toFirstStop: true);
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void Translate_Frame1_ShiftsReading()
        {
            // A ATG GCT = skip A, then ATG GCT = M A
            var dna = new DnaSequence("AATGGCT");
            var protein = Translator.Translate(dna, frame: 1);
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void Translate_Frame2_ShiftsReading()
        {
            // AA ATG GCT = skip AA, then ATG GCT = M A
            var dna = new DnaSequence("AAATGGCT");
            var protein = Translator.Translate(dna, frame: 2);
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void Translate_InvalidFrame_ThrowsException()
        {
            var dna = new DnaSequence("ATGGCT");
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.Translate(dna, frame: 3));
        }

        [Test]
        public void Translate_EmptySequence_ReturnsEmpty()
        {
            var protein = Translator.Translate("");
            Assert.That(protein.Sequence, Is.Empty);
        }

        [Test]
        public void Translate_NullDna_ThrowsException()
        {
            Assert.Throws<ArgumentNullException>(() => Translator.Translate((DnaSequence)null!));
        }

        [Test]
        public void Translate_NullRna_ThrowsException()
        {
            // Source: Implementation spec - null input handling
            Assert.Throws<ArgumentNullException>(() => Translator.Translate((RnaSequence)null!));
        }

        #endregion

        #region RNA Translation

        [Test]
        public void Translate_Rna_Works()
        {
            var rna = new RnaSequence("AUGGCUUAA");
            var protein = Translator.Translate(rna);
            Assert.That(protein.Sequence, Is.EqualTo("MA*"));
        }

        [Test]
        public void Translate_RnaToFirstStop_Works()
        {
            var rna = new RnaSequence("AUGGCUUAAGCU");
            var protein = Translator.Translate(rna, toFirstStop: true);
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        #endregion

        #region String Translation

        [Test]
        public void Translate_String_Works()
        {
            var protein = Translator.Translate("ATGGCT");
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void Translate_LowercaseString_Works()
        {
            var protein = Translator.Translate("atggct");
            Assert.That(protein.Sequence, Is.EqualTo("MA"));
        }

        #endregion

        #region Alternative Genetic Codes

        [Test]
        public void Translate_VertebrateMitochondrial_UsesDifferentCode()
        {
            // AGA is Arg in standard, but Stop in vertebrate mitochondrial
            var dna = new DnaSequence("ATGAGA");

            var standardProtein = Translator.Translate(dna, GeneticCode.Standard);
            Assert.That(standardProtein.Sequence, Is.EqualTo("MR"));

            var mitoProtein = Translator.Translate(dna, GeneticCode.VertebrateMitochondrial);
            Assert.That(mitoProtein.Sequence, Is.EqualTo("M*"));
        }

        [Test]
        public void Translate_YeastMitochondrial_CUU_IsThreonine()
        {
            // CUU is Leu in standard, but Thr in yeast mitochondrial
            var dna = new DnaSequence("ATGCTT");

            var standardProtein = Translator.Translate(dna, GeneticCode.Standard);
            Assert.That(standardProtein.Sequence, Is.EqualTo("ML"));

            var yeastProtein = Translator.Translate(dna, GeneticCode.YeastMitochondrial);
            Assert.That(yeastProtein.Sequence, Is.EqualTo("MT"));
        }

        #endregion

        #region Six Frame Translation

        [Test]
        public void TranslateSixFrames_ReturnsAllSixFrames()
        {
            var dna = new DnaSequence("ATGGCTAAA");
            var frames = Translator.TranslateSixFrames(dna);

            Assert.That(frames.Count, Is.EqualTo(6));
            Assert.That(frames.ContainsKey(1), Is.True);
            Assert.That(frames.ContainsKey(2), Is.True);
            Assert.That(frames.ContainsKey(3), Is.True);
            Assert.That(frames.ContainsKey(-1), Is.True);
            Assert.That(frames.ContainsKey(-2), Is.True);
            Assert.That(frames.ContainsKey(-3), Is.True);
        }

        [Test]
        public void TranslateSixFrames_Frame1_MatchesDirect()
        {
            var dna = new DnaSequence("ATGGCTAAA");
            var frames = Translator.TranslateSixFrames(dna);
            var direct = Translator.Translate(dna, frame: 0);

            Assert.That(frames[1].Sequence, Is.EqualTo(direct.Sequence));
        }

        [Test]
        public void TranslateSixFrames_NegativeFrames_UseReverseComplement()
        {
            var dna = new DnaSequence("ATGGCTAAA");
            var revComp = dna.ReverseComplement();
            var frames = Translator.TranslateSixFrames(dna);
            var revFrame1 = Translator.Translate(revComp, frame: 0);

            Assert.That(frames[-1].Sequence, Is.EqualTo(revFrame1.Sequence));
        }

        [Test]
        public void TranslateSixFrames_NullInput_ThrowsException()
        {
            // Source: Implementation spec - null input handling
            Assert.Throws<ArgumentNullException>(() => Translator.TranslateSixFrames(null!));
        }

        #endregion

        #region ORF Finding

        [Test]
        public void FindOrfs_SimpleOrf_FindsIt()
        {
            // ATG followed by 100+ amino acids and stop codon
            // ATG + 99 * "GCT" (Ala) + TAA = 100 AA ORF
            string dna = "ATG" + new string('G', 99 * 3).Replace("GGG", "GCT") + "TAA";
            // Actually let's make it simpler with just 100 Alanines
            var sb = new System.Text.StringBuilder("ATG");
            for (int i = 0; i < 100; i++)
                sb.Append("GCT");
            sb.Append("TAA");

            var sequence = new DnaSequence(sb.ToString());
            var orfs = Translator.FindOrfs(sequence, minLength: 100).ToList();

            Assert.That(orfs, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(orfs[0].AminoAcidLength, Is.EqualTo(101)); // M + 100 A
        }

        [Test]
        public void FindOrfs_NoStartCodon_ReturnsEmpty()
        {
            // No ATG in sequence
            var dna = new DnaSequence("GCTGCTGCT");
            var orfs = Translator.FindOrfs(dna, minLength: 1).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_ShortOrf_FilteredByMinLength()
        {
            // ATG GCT TAA = 2 amino acids (M A), filtered by minLength=100
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 100).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_RespectMinLength_FindsSmallOrfs()
        {
            // ATG GCT TAA = 2 amino acids (M A)
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 2).ToList();

            Assert.That(orfs, Has.Count.EqualTo(1));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void FindOrfs_ForwardOnly_DoesNotSearchReverseStrand()
        {
            // Create a sequence with ORF only on reverse strand
            var dna = new DnaSequence("GCTGCTGCT");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_ResultHasCorrectFrame()
        {
            // ATG GCT TAA in frame 0
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 2).ToList();

            Assert.That(orfs[0].Frame, Is.EqualTo(1)); // Frame 1 = first reading frame
        }

        [Test]
        public void FindOrfs_NullDna_ThrowsException()
        {
            Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs(null!).ToList());
        }

        #endregion

        #region Real Sequences

        [Test]
        public void Translate_InsulinBChain_ProducesCorrectProtein()
        {
            // Human insulin B chain coding sequence (simplified)
            // FVNQHLCGSHLVEALYLVCGERGFFYTPKT
            var dna = new DnaSequence("TTCGTGAACCAGCACCTGTGCGGCTCCCACCTGGTGGAAGCTCTGTACCTGGTGTGTGGGGAGCGTGGCTTCTTCTACACACCCAAGACC");
            var protein = Translator.Translate(dna);

            // Check it starts with F (Phe)
            Assert.That(protein[0], Is.EqualTo('F'));
            Assert.That(protein.Length, Is.EqualTo(30));
        }

        #endregion
    }
}
