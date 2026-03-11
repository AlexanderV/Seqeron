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

        [Test]
        public void Translate_SequenceShorterThan3_ReturnsEmpty()
        {
            // Less than one complete codon — no amino acid produced
            // Source: Wikipedia - codons are read in triplets
            var protein = Translator.Translate("AT");
            Assert.That(protein.Sequence, Is.Empty);
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
        public void Translate_DnaString_ConvertsTToU()
        {
            // DNA string with T must produce identical result to RNA with U
            // Source: Wikipedia - DNA T corresponds to RNA U in the genetic code
            var dnaResult = Translator.Translate("ATGGCTTAA");
            var rnaResult = Translator.Translate(new RnaSequence("AUGGCUUAA"));
            Assert.That(dnaResult.Sequence, Is.EqualTo("MA*"));
            Assert.That(dnaResult.Sequence, Is.EqualTo(rnaResult.Sequence));
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
            // All 3 negative frames must match direct translation of reverse complement
            // Source: Wikipedia Reading frame - 6 frames from double-stranded DNA
            var dna = new DnaSequence("ATGGCTAAA");
            var revComp = dna.ReverseComplement();
            var frames = Translator.TranslateSixFrames(dna);

            for (int f = 0; f < 3; f++)
            {
                var expected = Translator.Translate(revComp, frame: f);
                Assert.That(frames[-(f + 1)].Sequence, Is.EqualTo(expected.Sequence),
                    $"Frame -{f + 1} should match reverse complement frame {f}");
            }
        }

        [Test]
        public void TranslateSixFrames_NullInput_ThrowsException()
        {
            // Source: Implementation spec - null input handling
            Assert.Throws<ArgumentNullException>(() => Translator.TranslateSixFrames(null!));
        }

        [Test]
        public void TranslateSixFrames_EmptySequence_ReturnsEmptyFrames()
        {
            var dna = new DnaSequence("");
            var frames = Translator.TranslateSixFrames(dna);

            Assert.That(frames.Count, Is.EqualTo(6));
            Assert.That(frames.Values.All(p => p.Sequence == ""), Is.True);
        }

        #endregion

        #region ORF Finding (Smoke Tests — canonical tests in GenomeAnnotator_ORF_Tests.cs)

        [Test]
        public void FindOrfs_SimpleOrf_FindsIt()
        {
            // ATG GCT TTC TAA = M A F * → protein "MAF" (3 amino acids)
            // Source: Wikipedia ORF - start codon to stop codon
            var dna = new DnaSequence("ATGGCTTTCTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.That(orfs, Has.Count.EqualTo(1));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MAF"));
        }

        [Test]
        public void FindOrfs_NoStartCodon_ReturnsEmpty()
        {
            // No ATG/TTG/CTG (start codons) in any reading frame → no ORF
            // Source: Wikipedia ORF - requires start codon
            var dna = new DnaSequence("GCCGCCGCCTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_RespectMinLength_FindsSmallOrfs()
        {
            // ATG GCT TAA = 2 amino acids (M A), meets minLength=2
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 2, searchBothStrands: false).ToList();

            Assert.That(orfs, Has.Count.EqualTo(1));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MA"));
        }

        [Test]
        public void FindOrfs_ShortOrf_FilteredByMinLength()
        {
            // ATG GCT TAA = protein "MA" (2 aa) < minLength 5 → filtered out
            // Source: Wikipedia ORF - short ORFs excluded by length threshold
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 5, searchBothStrands: false).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_ForwardOnly_DoesNotSearchReverseStrand()
        {
            // TTAGCCGCCCAT has no start codons on forward strand (any frame)
            // Reverse complement = ATGGGCGGCTAA which has ATG...TAA ORF
            // With searchBothStrands=false, reverse-strand ORF must not appear
            var dna = new DnaSequence("TTAGCCGCCCAT");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.That(orfs, Is.Empty);
        }

        [Test]
        public void FindOrfs_BothStrands_SearchesReverseComplement()
        {
            // Forward: TTA GCC GCC CAT → no ATG → no ORFs
            // Reverse complement: ATG GGC GGC TAA → M G G * → ORF found
            // Source: Wikipedia ORF - six-frame search covers both strands
            var dna = new DnaSequence("TTAGCCGCCCAT");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: true).ToList();

            Assert.That(orfs, Has.Count.EqualTo(1));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MGG"));
            Assert.That(orfs[0].Frame, Is.Negative, "ORF should be on reverse strand");
        }

        [Test]
        public void FindOrfs_OrfResult_HasCorrectPositions()
        {
            // ATG GCT TAA at positions 0-8 in frame 0 (reported as frame 1)
            // Source: Wikipedia ORF - ORF spans from start to stop codon
            var dna = new DnaSequence("ATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(orfs, Has.Count.EqualTo(1));
                Assert.That(orfs[0].StartPosition, Is.EqualTo(0));
                Assert.That(orfs[0].EndPosition, Is.EqualTo(8));
                Assert.That(orfs[0].Frame, Is.EqualTo(1));
                Assert.That(orfs[0].NucleotideLength, Is.EqualTo(9));
                Assert.That(orfs[0].AminoAcidLength, Is.EqualTo(2));
            });
        }

        [Test]
        public void FindOrfs_NullDna_ThrowsException()
        {
            Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs(null!).ToList());
        }

        [Test]
        public void FindOrfs_MultipleOrfs_FindsAll()
        {
            // Two ORFs in same frame: ATG GCT TAA | ATG GCT TAA
            // Source: Wikipedia ORF - multiple ORFs can exist in same sequence
            var dna = new DnaSequence("ATGGCTTAAATGGCTTAA");
            var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

            Assert.That(orfs, Has.Count.EqualTo(2));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MA"));
            Assert.That(orfs[1].Protein.Sequence, Is.EqualTo("MA"));
        }

        #endregion

        #region Real Sequences

        [Test]
        public void Translate_InsulinBChain_ProducesCorrectProtein()
        {
            // Human insulin B chain coding sequence
            // Source: UniProt P01308, positions 25-54 of preproinsulin
            // DNA from NCBI RefSeq NM_000207.3
            var dna = new DnaSequence("TTCGTGAACCAGCACCTGTGCGGCTCCCACCTGGTGGAAGCTCTGTACCTGGTGTGTGGGGAGCGTGGCTTCTTCTACACACCCAAGACC");
            var protein = Translator.Translate(dna);

            Assert.That(protein.Sequence, Is.EqualTo("FVNQHLCGSHLVEALYLVCGERGFFYTPKT"));
        }

        #endregion
    }
}
