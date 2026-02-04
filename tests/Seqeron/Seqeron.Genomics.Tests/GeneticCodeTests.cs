using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for GeneticCode codon translation.
    /// Test Unit: TRANS-CODON-001
    /// 
    /// Evidence Sources:
    ///   - Wikipedia: Genetic code, Start codon, Stop codon
    ///   - NCBI: The Genetic Codes (translation tables 1-33)
    ///   - Crick FH (1968): The origin of the genetic code
    /// </summary>
    [TestFixture]
    public class GeneticCodeTests
    {
        #region Standard Genetic Code

        /// <summary>
        /// Verifies GeneticCode metadata.
        /// Source: Implementation specification
        /// </summary>
        [Test]
        public void Standard_HasCorrectName()
        {
            Assert.That(GeneticCode.Standard.Name, Is.EqualTo("Standard"));
            Assert.That(GeneticCode.Standard.TableNumber, Is.EqualTo(1));
        }

        /// <summary>
        /// Verifies that standard genetic code has all 64 codons (4³).
        /// Source: NCBI Translation Table 1
        /// </summary>
        [Test]
        public void Standard_Has64Codons()
        {
            Assert.That(GeneticCode.Standard.CodonTable.Count, Is.EqualTo(64));
        }

        /// <summary>
        /// Verifies the three standard stop codons: UAA (ochre), UAG (amber), UGA (opal).
        /// Source: Wikipedia (Stop codon), NCBI Table 1
        /// </summary>
        [Test]
        public void Standard_HasThreeStopCodons()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.StopCodons.Count, Is.EqualTo(3));
                Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UAA"), "Ochre stop codon");
                Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UAG"), "Amber stop codon");
                Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UGA"), "Opal stop codon");
            });
        }

        /// <summary>
        /// Standard code has AUG as the only designated start codon.
        /// Source: NCBI Table 1, Wikipedia (Start codon)
        /// </summary>
        [Test]
        public void Standard_HasOneStartCodon()
        {
            Assert.That(GeneticCode.Standard.StartCodons.Count, Is.EqualTo(1));
            Assert.That(GeneticCode.Standard.StartCodons, Does.Contain("AUG"));
        }

        #endregion

        #region Translate

        /// <summary>
        /// AUG (start codon) translates to Methionine in all contexts.
        /// Source: NCBI Table 1, Wikipedia (Start codon)
        /// </summary>
        [Test]
        public void Translate_AUG_ReturnsMethionine()
        {
            char aa = GeneticCode.Standard.Translate("AUG");
            Assert.That(aa, Is.EqualTo('M'));
        }

        /// <summary>
        /// DNA codons (T instead of U) are automatically normalized.
        /// Source: Implementation specification (T→U conversion)
        /// </summary>
        [Test]
        public void Translate_DnaCodon_Works()
        {
            // ATG should be converted to AUG internally
            char aa = GeneticCode.Standard.Translate("ATG");
            Assert.That(aa, Is.EqualTo('M'));
        }

        /// <summary>
        /// Codons are case-insensitive.
        /// Source: Implementation specification
        /// </summary>
        [Test]
        public void Translate_LowercaseCodon_Works()
        {
            char aa = GeneticCode.Standard.Translate("aug");
            Assert.That(aa, Is.EqualTo('M'));
        }

        /// <summary>
        /// All three stop codons translate to '*'.
        /// Source: NCBI format convention, Wikipedia (Stop codon)
        /// </summary>
        [Test]
        public void Translate_AllStopCodons_ReturnsAsterisk()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.Translate("UAA"), Is.EqualTo('*'), "Ochre");
                Assert.That(GeneticCode.Standard.Translate("UAG"), Is.EqualTo('*'), "Amber");
                Assert.That(GeneticCode.Standard.Translate("UGA"), Is.EqualTo('*'), "Opal");
            });
        }

        /// <summary>
        /// Codons must be exactly 3 characters.
        /// Source: Definition of codon (triplet code)
        /// </summary>
        [Test]
        public void Translate_InvalidCodonLength_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate("AU"));
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate("AUGC"));
        }

        /// <summary>
        /// Empty codon is invalid.
        /// Source: Definition of codon
        /// </summary>
        [Test]
        public void Translate_EmptyCodon_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate(""));
        }

        /// <summary>
        /// Null codon is invalid.
        /// Source: Implementation specification
        /// </summary>
        [Test]
        public void Translate_NullCodon_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate(null!));
        }

        /// <summary>
        /// All 64 codons produce valid single-letter amino acid codes.
        /// Source: NCBI Table 1 (20 amino acids + stop)
        /// </summary>
        [Test]
        public void Translate_AllCodons_ProduceValidAminoAcids()
        {
            var validAa = "ACDEFGHIKLMNPQRSTVWY*";
            foreach (var codon in GeneticCode.Standard.CodonTable.Keys)
            {
                char aa = GeneticCode.Standard.Translate(codon);
                Assert.That(validAa, Does.Contain(aa.ToString()), $"Codon {codon} produced invalid AA '{aa}'");
            }
        }

        /// <summary>
        /// Complete verification of all 64 standard codons against NCBI Table 1.
        /// Source: NCBI Translation Table 1 (Standard)
        /// </summary>
        [Test]
        public void Translate_CompleteStandardCodonTable_MatchesNcbi()
        {
            // Expected amino acids for all 64 codons (NCBI Table 1)
            // Order: UUU, UUC, UUA, UUG, CUU, CUC, CUA, CUG, AUU, AUC, AUA, AUG, GUU, GUC, GUA, GUG
            //        UCU, UCC, UCA, UCG, CCU, CCC, CCA, CCG, ACU, ACC, ACA, ACG, GCU, GCC, GCA, GCG
            //        UAU, UAC, UAA, UAG, CAU, CAC, CAA, CAG, AAU, AAC, AAA, AAG, GAU, GAC, GAA, GAG
            //        UGU, UGC, UGA, UGG, CGU, CGC, CGA, CGG, AGU, AGC, AGA, AGG, GGU, GGC, GGA, GGG
            var expectedTable = new Dictionary<string, char>
            {
                // First codon position U
                ["UUU"] = 'F', ["UUC"] = 'F', ["UUA"] = 'L', ["UUG"] = 'L',
                ["UCU"] = 'S', ["UCC"] = 'S', ["UCA"] = 'S', ["UCG"] = 'S',
                ["UAU"] = 'Y', ["UAC"] = 'Y', ["UAA"] = '*', ["UAG"] = '*',
                ["UGU"] = 'C', ["UGC"] = 'C', ["UGA"] = '*', ["UGG"] = 'W',
                
                // First codon position C
                ["CUU"] = 'L', ["CUC"] = 'L', ["CUA"] = 'L', ["CUG"] = 'L',
                ["CCU"] = 'P', ["CCC"] = 'P', ["CCA"] = 'P', ["CCG"] = 'P',
                ["CAU"] = 'H', ["CAC"] = 'H', ["CAA"] = 'Q', ["CAG"] = 'Q',
                ["CGU"] = 'R', ["CGC"] = 'R', ["CGA"] = 'R', ["CGG"] = 'R',
                
                // First codon position A
                ["AUU"] = 'I', ["AUC"] = 'I', ["AUA"] = 'I', ["AUG"] = 'M',
                ["ACU"] = 'T', ["ACC"] = 'T', ["ACA"] = 'T', ["ACG"] = 'T',
                ["AAU"] = 'N', ["AAC"] = 'N', ["AAA"] = 'K', ["AAG"] = 'K',
                ["AGU"] = 'S', ["AGC"] = 'S', ["AGA"] = 'R', ["AGG"] = 'R',
                
                // First codon position G
                ["GUU"] = 'V', ["GUC"] = 'V', ["GUA"] = 'V', ["GUG"] = 'V',
                ["GCU"] = 'A', ["GCC"] = 'A', ["GCA"] = 'A', ["GCG"] = 'A',
                ["GAU"] = 'D', ["GAC"] = 'D', ["GAA"] = 'E', ["GAG"] = 'E',
                ["GGU"] = 'G', ["GGC"] = 'G', ["GGA"] = 'G', ["GGG"] = 'G'
            };

            Assert.That(expectedTable.Count, Is.EqualTo(64), "Test data completeness");

            foreach (var (codon, expectedAa) in expectedTable)
            {
                var actualAa = GeneticCode.Standard.Translate(codon);
                Assert.That(actualAa, Is.EqualTo(expectedAa), 
                    $"Codon {codon}: expected {expectedAa}, got {actualAa}");
            }
        }

        /// <summary>
        /// Mixed case codons work correctly.
        /// Source: Implementation specification (case-insensitive)
        /// </summary>
        [Test]
        public void Translate_MixedCaseCodon_Works()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.Translate("AuG"), Is.EqualTo('M'));
                Assert.That(GeneticCode.Standard.Translate("uAa"), Is.EqualTo('*'));
                Assert.That(GeneticCode.Standard.Translate("Uuu"), Is.EqualTo('F'));
            });
        }

        #endregion

        #region IsStartCodon / IsStopCodon

        /// <summary>
        /// AUG is the universal start codon (RNA format).
        /// Source: NCBI Table 1, Wikipedia (Start codon)
        /// </summary>
        [Test]
        public void IsStartCodon_AUG_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStartCodon("AUG"), Is.True);
        }

        /// <summary>
        /// ATG (DNA format) is recognized as start codon.
        /// Source: T→U normalization
        /// </summary>
        [Test]
        public void IsStartCodon_ATG_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStartCodon("ATG"), Is.True);
        }

        /// <summary>
        /// Non-start codons return false.
        /// Source: NCBI Table 1 (only AUG is start in standard code)
        /// </summary>
        [Test]
        public void IsStartCodon_NonStartCodon_ReturnsFalse()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.IsStartCodon("UUU"), Is.False);
                Assert.That(GeneticCode.Standard.IsStartCodon("GUG"), Is.False, "GUG is not start in standard code");
                Assert.That(GeneticCode.Standard.IsStartCodon("UUG"), Is.False, "UUG is not start in standard code");
            });
        }

        /// <summary>
        /// Invalid input returns false (not exception).
        /// Source: Implementation specification (defensive behavior)
        /// </summary>
        [Test]
        public void IsStartCodon_InvalidInput_ReturnsFalse()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.IsStartCodon("AU"), Is.False);
                Assert.That(GeneticCode.Standard.IsStartCodon(""), Is.False);
                Assert.That(GeneticCode.Standard.IsStartCodon(null!), Is.False);
            });
        }

        /// <summary>
        /// UAA (ochre) is a stop codon.
        /// Source: Wikipedia (Stop codon), NCBI Table 1
        /// </summary>
        [Test]
        public void IsStopCodon_UAA_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStopCodon("UAA"), Is.True);
        }

        /// <summary>
        /// All three standard stop codons are recognized.
        /// Source: NCBI Table 1, Wikipedia (Stop codon)
        /// </summary>
        [Test]
        public void IsStopCodon_AllThreeStandard_ReturnsTrue()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.IsStopCodon("UAA"), Is.True, "Ochre");
                Assert.That(GeneticCode.Standard.IsStopCodon("UAG"), Is.True, "Amber");
                Assert.That(GeneticCode.Standard.IsStopCodon("UGA"), Is.True, "Opal");
            });
        }

        /// <summary>
        /// DNA format stop codons are recognized.
        /// Source: T→U normalization
        /// </summary>
        [Test]
        public void IsStopCodon_DnaFormat_ReturnsTrue()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.IsStopCodon("TAA"), Is.True);
                Assert.That(GeneticCode.Standard.IsStopCodon("TAG"), Is.True);
                Assert.That(GeneticCode.Standard.IsStopCodon("TGA"), Is.True);
            });
        }

        /// <summary>
        /// Non-stop codons return false.
        /// Source: NCBI Table 1
        /// </summary>
        [Test]
        public void IsStopCodon_NonStopCodon_ReturnsFalse()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.Standard.IsStopCodon("AUG"), Is.False);
                Assert.That(GeneticCode.Standard.IsStopCodon("UUU"), Is.False);
            });
        }

        #endregion

        #region GetCodonsForAminoAcid (Degeneracy)

        /// <summary>
        /// Methionine has only one codon (AUG) - no degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Methionine_ReturnsOneCodon()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('M').ToList();
            Assert.That(codons, Has.Count.EqualTo(1));
            Assert.That(codons[0], Is.EqualTo("AUG"));
        }

        /// <summary>
        /// Leucine has six codons (UUA, UUG, CUU, CUC, CUA, CUG) - 6-fold degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Leucine_ReturnsSixCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('L').ToList();
            Assert.That(codons, Has.Count.EqualTo(6));
            Assert.That(codons, Does.Contain("UUA"));
            Assert.That(codons, Does.Contain("UUG"));
            Assert.That(codons, Does.Contain("CUU"));
            Assert.That(codons, Does.Contain("CUC"));
            Assert.That(codons, Does.Contain("CUA"));
            Assert.That(codons, Does.Contain("CUG"));
        }

        /// <summary>
        /// Serine has six codons (UCU, UCC, UCA, UCG, AGU, AGC) - 6-fold degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Serine_ReturnsSixCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('S').ToList();
            Assert.That(codons, Has.Count.EqualTo(6));
        }

        /// <summary>
        /// Arginine has six codons - 6-fold degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Arginine_ReturnsSixCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('R').ToList();
            Assert.That(codons, Has.Count.EqualTo(6));
            Assert.That(codons, Does.Contain("AGA"));
            Assert.That(codons, Does.Contain("AGG"));
        }

        /// <summary>
        /// Tryptophan has only one codon (UGG) - no degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Tryptophan_ReturnsOneCodon()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('W').ToList();
            Assert.That(codons, Has.Count.EqualTo(1));
            Assert.That(codons[0], Is.EqualTo("UGG"));
        }

        /// <summary>
        /// Isoleucine has three codons (AUU, AUC, AUA) - 3-fold degeneracy.
        /// Source: NCBI Table 1 degeneracy
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Isoleucine_ReturnsThreeCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('I').ToList();
            Assert.That(codons, Has.Count.EqualTo(3));
        }

        /// <summary>
        /// Stop codons (*) return the three standard stops.
        /// Source: NCBI Table 1
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_Stop_ReturnsThreeCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('*').ToList();
            Assert.That(codons, Has.Count.EqualTo(3));
            Assert.That(codons, Does.Contain("UAA"));
            Assert.That(codons, Does.Contain("UAG"));
            Assert.That(codons, Does.Contain("UGA"));
        }

        /// <summary>
        /// Lowercase amino acid input works (case-insensitive).
        /// Source: Implementation specification
        /// </summary>
        [Test]
        public void GetCodonsForAminoAcid_LowercaseInput_Works()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('m').ToList();
            Assert.That(codons, Has.Count.EqualTo(1));
            Assert.That(codons[0], Is.EqualTo("AUG"));
        }

        #endregion

        #region Alternative Genetic Codes - Vertebrate Mitochondrial (Table 2)

        /// <summary>
        /// In vertebrate mitochondrial code, UGA encodes Tryptophan (not Stop).
        /// Source: NCBI Table 2, Wikipedia (Genetic code - variations)
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_UGA_IsTryptophan()
        {
            char aa = GeneticCode.VertebrateMitochondrial.Translate("UGA");
            Assert.That(aa, Is.EqualTo('W'));
        }

        /// <summary>
        /// In vertebrate mitochondrial code, AGA is a stop codon (not Arginine).
        /// Source: NCBI Table 2
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_AGA_IsStopCodon()
        {
            Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("AGA"), Is.True);
            Assert.That(GeneticCode.VertebrateMitochondrial.Translate("AGA"), Is.EqualTo('*'));
        }

        /// <summary>
        /// In vertebrate mitochondrial code, AGG is also a stop codon (not Arginine).
        /// Source: NCBI Table 2
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_AGG_IsStopCodon()
        {
            Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("AGG"), Is.True);
            Assert.That(GeneticCode.VertebrateMitochondrial.Translate("AGG"), Is.EqualTo('*'));
        }

        /// <summary>
        /// In vertebrate mitochondrial code, AUA encodes Methionine (not Isoleucine).
        /// Source: NCBI Table 2
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_AUA_IsMethionine()
        {
            char aa = GeneticCode.VertebrateMitochondrial.Translate("AUA");
            Assert.That(aa, Is.EqualTo('M'));
        }

        /// <summary>
        /// Vertebrate mitochondrial code has four start codons: AUG, AUA, AUU, AUC.
        /// Source: NCBI Table 2
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_HasFourStartCodons()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.VertebrateMitochondrial.StartCodons.Count, Is.EqualTo(4));
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("AUG"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("AUA"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("AUU"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStartCodon("AUC"), Is.True);
            });
        }

        /// <summary>
        /// Vertebrate mitochondrial code has four stop codons: UAA, UAG, AGA, AGG.
        /// Source: NCBI Table 2
        /// </summary>
        [Test]
        public void VertebrateMitochondrial_HasFourStopCodons()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.VertebrateMitochondrial.StopCodons.Count, Is.EqualTo(4));
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("UAA"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("UAG"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("AGA"), Is.True);
                Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("AGG"), Is.True);
            });
        }

        #endregion

        #region Alternative Genetic Codes - Yeast Mitochondrial (Table 3)

        /// <summary>
        /// In yeast mitochondrial code, CUU encodes Threonine (not Leucine).
        /// Source: NCBI Table 3
        /// </summary>
        [Test]
        public void YeastMitochondrial_CUU_IsThreonine()
        {
            char aa = GeneticCode.YeastMitochondrial.Translate("CUU");
            Assert.That(aa, Is.EqualTo('T'));
        }

        /// <summary>
        /// All four CUx codons encode Threonine in yeast mitochondrial code.
        /// Source: NCBI Table 3
        /// </summary>
        [Test]
        public void YeastMitochondrial_AllCUxCodons_AreThreonine()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.YeastMitochondrial.Translate("CUU"), Is.EqualTo('T'));
                Assert.That(GeneticCode.YeastMitochondrial.Translate("CUC"), Is.EqualTo('T'));
                Assert.That(GeneticCode.YeastMitochondrial.Translate("CUA"), Is.EqualTo('T'));
                Assert.That(GeneticCode.YeastMitochondrial.Translate("CUG"), Is.EqualTo('T'));
            });
        }

        /// <summary>
        /// In yeast mitochondrial code, UGA encodes Tryptophan (not Stop).
        /// Source: NCBI Table 3
        /// </summary>
        [Test]
        public void YeastMitochondrial_UGA_IsTryptophan()
        {
            char aa = GeneticCode.YeastMitochondrial.Translate("UGA");
            Assert.That(aa, Is.EqualTo('W'));
            Assert.That(GeneticCode.YeastMitochondrial.IsStopCodon("UGA"), Is.False);
        }

        /// <summary>
        /// Yeast mitochondrial code has only two stop codons: UAA, UAG.
        /// Source: NCBI Table 3
        /// </summary>
        [Test]
        public void YeastMitochondrial_HasTwoStopCodons()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.YeastMitochondrial.StopCodons.Count, Is.EqualTo(2));
                Assert.That(GeneticCode.YeastMitochondrial.IsStopCodon("UAA"), Is.True);
                Assert.That(GeneticCode.YeastMitochondrial.IsStopCodon("UAG"), Is.True);
            });
        }

        #endregion

        #region Alternative Genetic Codes - Bacterial/Plastid (Table 11)

        /// <summary>
        /// Bacterial/Plastid code uses alternative start codons: AUG, GUG, UUG.
        /// Source: NCBI Table 11
        /// </summary>
        [Test]
        public void BacterialPlastid_HasAlternativeStartCodons()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("AUG"), Is.True);
                Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("GUG"), Is.True);
                Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("UUG"), Is.True);
            });
        }

        /// <summary>
        /// Bacterial/Plastid code has same codon table as standard, just different starts.
        /// Source: NCBI Table 11
        /// </summary>
        [Test]
        public void BacterialPlastid_CodonTable_SameAsStandard()
        {
            // Verify a few key translations match standard
            Assert.Multiple(() =>
            {
                Assert.That(GeneticCode.BacterialPlastid.Translate("AUG"), Is.EqualTo('M'));
                Assert.That(GeneticCode.BacterialPlastid.Translate("UGA"), Is.EqualTo('*'));
                Assert.That(GeneticCode.BacterialPlastid.Translate("AGA"), Is.EqualTo('R'));
            });
        }

        #endregion

        #region GetByTableNumber

        /// <summary>
        /// NCBI Table 1 (Standard) is retrieved by table number.
        /// Source: NCBI Genetic Codes tables numbering
        /// </summary>
        [Test]
        public void GetByTableNumber_1_ReturnsStandard()
        {
            var code = GeneticCode.GetByTableNumber(1);
            Assert.Multiple(() =>
            {
                Assert.That(code, Is.EqualTo(GeneticCode.Standard));
                Assert.That(code.TableNumber, Is.EqualTo(1));
                Assert.That(code.Name, Is.EqualTo("Standard"));
            });
        }

        /// <summary>
        /// NCBI Table 2 (Vertebrate Mitochondrial) is retrieved by table number.
        /// Source: NCBI Genetic Codes tables numbering
        /// </summary>
        [Test]
        public void GetByTableNumber_2_ReturnsVertebrateMitochondrial()
        {
            var code = GeneticCode.GetByTableNumber(2);
            Assert.Multiple(() =>
            {
                Assert.That(code, Is.EqualTo(GeneticCode.VertebrateMitochondrial));
                Assert.That(code.TableNumber, Is.EqualTo(2));
                Assert.That(code.Name, Is.EqualTo("Vertebrate Mitochondrial"));
            });
        }

        /// <summary>
        /// NCBI Table 3 (Yeast Mitochondrial) is retrieved by table number.
        /// Source: NCBI Genetic Codes tables numbering
        /// </summary>
        [Test]
        public void GetByTableNumber_3_ReturnsYeastMitochondrial()
        {
            var code = GeneticCode.GetByTableNumber(3);
            Assert.Multiple(() =>
            {
                Assert.That(code, Is.EqualTo(GeneticCode.YeastMitochondrial));
                Assert.That(code.TableNumber, Is.EqualTo(3));
                Assert.That(code.Name, Is.EqualTo("Yeast Mitochondrial"));
            });
        }

        /// <summary>
        /// NCBI Table 11 (Bacterial, Archaeal and Plant Plastid) is retrieved by table number.
        /// Source: NCBI Genetic Codes tables numbering
        /// </summary>
        [Test]
        public void GetByTableNumber_11_ReturnsBacterialPlastid()
        {
            var code = GeneticCode.GetByTableNumber(11);
            Assert.Multiple(() =>
            {
                Assert.That(code, Is.EqualTo(GeneticCode.BacterialPlastid));
                Assert.That(code.TableNumber, Is.EqualTo(11));
                Assert.That(code.Name, Is.EqualTo("Bacterial, Archaeal and Plant Plastid"));
            });
        }

        /// <summary>
        /// Invalid table number throws ArgumentException.
        /// Source: Implementation specification (explicit error handling)
        /// </summary>
        [Test]
        public void GetByTableNumber_Invalid_ThrowsException()
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentException>(() => GeneticCode.GetByTableNumber(0));
                Assert.Throws<ArgumentException>(() => GeneticCode.GetByTableNumber(99));
                Assert.Throws<ArgumentException>(() => GeneticCode.GetByTableNumber(-1));
            });
        }

        /// <summary>
        /// All supported table numbers return valid genetic codes.
        /// Source: NCBI Genetic Codes
        /// </summary>
        [Test]
        public void GetByTableNumber_AllSupported_ReturnsValidCodes()
        {
            int[] supportedTables = { 1, 2, 3, 11 };
            
            foreach (int tableNum in supportedTables)
            {
                var code = GeneticCode.GetByTableNumber(tableNum);
                Assert.That(code, Is.Not.Null, $"Table {tableNum} should return a valid code");
                Assert.That(code.TableNumber, Is.EqualTo(tableNum), $"Table number should match");
            }
        }

        #endregion
    }
}
