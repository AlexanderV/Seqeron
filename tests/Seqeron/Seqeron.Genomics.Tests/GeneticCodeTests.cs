using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    [TestFixture]
    public class GeneticCodeTests
    {
        #region Standard Genetic Code

        [Test]
        public void Standard_HasCorrectName()
        {
            Assert.That(GeneticCode.Standard.Name, Is.EqualTo("Standard"));
            Assert.That(GeneticCode.Standard.TableNumber, Is.EqualTo(1));
        }

        [Test]
        public void Standard_Has64Codons()
        {
            Assert.That(GeneticCode.Standard.CodonTable.Count, Is.EqualTo(64));
        }

        [Test]
        public void Standard_HasThreeStopCodons()
        {
            Assert.That(GeneticCode.Standard.StopCodons.Count, Is.EqualTo(3));
            Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UAA"));
            Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UAG"));
            Assert.That(GeneticCode.Standard.StopCodons, Does.Contain("UGA"));
        }

        [Test]
        public void Standard_HasOneStartCodon()
        {
            Assert.That(GeneticCode.Standard.StartCodons.Count, Is.EqualTo(1));
            Assert.That(GeneticCode.Standard.StartCodons, Does.Contain("AUG"));
        }

        #endregion

        #region Translate

        [Test]
        public void Translate_AUG_ReturnsMethionine()
        {
            char aa = GeneticCode.Standard.Translate("AUG");
            Assert.That(aa, Is.EqualTo('M'));
        }

        [Test]
        public void Translate_DnaCodon_Works()
        {
            // ATG should be converted to AUG internally
            char aa = GeneticCode.Standard.Translate("ATG");
            Assert.That(aa, Is.EqualTo('M'));
        }

        [Test]
        public void Translate_LowercaseCodon_Works()
        {
            char aa = GeneticCode.Standard.Translate("aug");
            Assert.That(aa, Is.EqualTo('M'));
        }

        [Test]
        public void Translate_StopCodon_ReturnsAsterisk()
        {
            Assert.That(GeneticCode.Standard.Translate("UAA"), Is.EqualTo('*'));
            Assert.That(GeneticCode.Standard.Translate("UAG"), Is.EqualTo('*'));
            Assert.That(GeneticCode.Standard.Translate("UGA"), Is.EqualTo('*'));
        }

        [Test]
        public void Translate_InvalidCodonLength_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate("AU"));
            Assert.Throws<ArgumentException>(() => GeneticCode.Standard.Translate("AUGC"));
        }

        [Test]
        public void Translate_AllCodons_ProduceValidAminoAcids()
        {
            var validAa = "ACDEFGHIKLMNPQRSTVWY*";
            foreach (var codon in GeneticCode.Standard.CodonTable.Keys)
            {
                char aa = GeneticCode.Standard.Translate(codon);
                Assert.That(validAa, Does.Contain(aa.ToString()));
            }
        }

        #endregion

        #region IsStartCodon / IsStopCodon

        [Test]
        public void IsStartCodon_AUG_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStartCodon("AUG"), Is.True);
        }

        [Test]
        public void IsStartCodon_ATG_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStartCodon("ATG"), Is.True);
        }

        [Test]
        public void IsStartCodon_Other_ReturnsFalse()
        {
            Assert.That(GeneticCode.Standard.IsStartCodon("UUU"), Is.False);
        }

        [Test]
        public void IsStopCodon_UAA_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStopCodon("UAA"), Is.True);
        }

        [Test]
        public void IsStopCodon_TAA_ReturnsTrue()
        {
            Assert.That(GeneticCode.Standard.IsStopCodon("TAA"), Is.True);
        }

        [Test]
        public void IsStopCodon_Other_ReturnsFalse()
        {
            Assert.That(GeneticCode.Standard.IsStopCodon("AUG"), Is.False);
        }

        #endregion

        #region GetCodonsForAminoAcid

        [Test]
        public void GetCodonsForAminoAcid_Methionine_ReturnsOneCodon()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('M').ToList();
            Assert.That(codons, Has.Count.EqualTo(1));
            Assert.That(codons[0], Is.EqualTo("AUG"));
        }

        [Test]
        public void GetCodonsForAminoAcid_Leucine_ReturnsSixCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('L').ToList();
            Assert.That(codons, Has.Count.EqualTo(6));
        }

        [Test]
        public void GetCodonsForAminoAcid_Serine_ReturnsSixCodons()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('S').ToList();
            Assert.That(codons, Has.Count.EqualTo(6));
        }

        [Test]
        public void GetCodonsForAminoAcid_Tryptophan_ReturnsOneCodon()
        {
            var codons = GeneticCode.Standard.GetCodonsForAminoAcid('W').ToList();
            Assert.That(codons, Has.Count.EqualTo(1));
            Assert.That(codons[0], Is.EqualTo("UGG"));
        }

        #endregion

        #region Alternative Genetic Codes

        [Test]
        public void VertebrateMitochondrial_UGA_IsTryptophan()
        {
            char aa = GeneticCode.VertebrateMitochondrial.Translate("UGA");
            Assert.That(aa, Is.EqualTo('W'));
        }

        [Test]
        public void VertebrateMitochondrial_AGA_IsStopCodon()
        {
            Assert.That(GeneticCode.VertebrateMitochondrial.IsStopCodon("AGA"), Is.True);
            Assert.That(GeneticCode.VertebrateMitochondrial.Translate("AGA"), Is.EqualTo('*'));
        }

        [Test]
        public void VertebrateMitochondrial_AUA_IsMethionine()
        {
            char aa = GeneticCode.VertebrateMitochondrial.Translate("AUA");
            Assert.That(aa, Is.EqualTo('M'));
        }

        [Test]
        public void YeastMitochondrial_CUU_IsThreonine()
        {
            char aa = GeneticCode.YeastMitochondrial.Translate("CUU");
            Assert.That(aa, Is.EqualTo('T'));
        }

        [Test]
        public void BacterialPlastid_HasAlternativeStartCodons()
        {
            Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("AUG"), Is.True);
            Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("GUG"), Is.True);
            Assert.That(GeneticCode.BacterialPlastid.IsStartCodon("UUG"), Is.True);
        }

        #endregion

        #region GetByTableNumber

        [Test]
        public void GetByTableNumber_1_ReturnsStandard()
        {
            var code = GeneticCode.GetByTableNumber(1);
            Assert.That(code, Is.EqualTo(GeneticCode.Standard));
        }

        [Test]
        public void GetByTableNumber_2_ReturnsVertebrateMitochondrial()
        {
            var code = GeneticCode.GetByTableNumber(2);
            Assert.That(code, Is.EqualTo(GeneticCode.VertebrateMitochondrial));
        }

        [Test]
        public void GetByTableNumber_Invalid_ThrowsException()
        {
            Assert.Throws<ArgumentException>(() => GeneticCode.GetByTableNumber(99));
        }

        #endregion
    }
}
