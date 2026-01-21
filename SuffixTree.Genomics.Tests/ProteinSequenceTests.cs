using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    [TestFixture]
    public class ProteinSequenceTests
    {
        #region Construction

        [Test]
        public void Constructor_ValidSequence_CreatesSequence()
        {
            var protein = new ProteinSequence("ACDEFGHIKLMNPQRSTVWY");
            Assert.That(protein.Sequence, Is.EqualTo("ACDEFGHIKLMNPQRSTVWY"));
        }

        [Test]
        public void Constructor_LowercaseSequence_NormalizesToUppercase()
        {
            var protein = new ProteinSequence("acdefg");
            Assert.That(protein.Sequence, Is.EqualTo("ACDEFG"));
        }

        [Test]
        public void Constructor_WithStopCodon_Succeeds()
        {
            var protein = new ProteinSequence("MKVL*");
            Assert.That(protein.Sequence, Is.EqualTo("MKVL*"));
        }

        [Test]
        public void Constructor_WithUnknownAminoAcid_Succeeds()
        {
            var protein = new ProteinSequence("MXVLX");
            Assert.That(protein.Sequence, Is.EqualTo("MXVLX"));
        }

        [Test]
        public void Constructor_InvalidAminoAcid_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() => new ProteinSequence("MKVLZ"));
        }

        [Test]
        public void Constructor_EmptySequence_CreatesEmpty()
        {
            var protein = new ProteinSequence("");
            Assert.That(protein.Sequence, Is.Empty);
            Assert.That(protein.Length, Is.EqualTo(0));
        }

        [Test]
        public void TryCreate_ValidSequence_ReturnsTrue()
        {
            bool result = ProteinSequence.TryCreate("MKVL", out var protein);
            Assert.That(result, Is.True);
            Assert.That(protein!.Sequence, Is.EqualTo("MKVL"));
        }

        [Test]
        public void TryCreate_InvalidSequence_ReturnsFalse()
        {
            bool result = ProteinSequence.TryCreate("MKVLZ", out var protein);
            Assert.That(result, Is.False);
            Assert.That(protein, Is.Null);
        }

        #endregion

        #region Molecular Weight

        [Test]
        public void MolecularWeight_SingleAminoAcid_ReturnsWeight()
        {
            var protein = new ProteinSequence("M");
            // Methionine weight: 149.21
            Assert.That(protein.MolecularWeight(), Is.EqualTo(149.21));
        }

        [Test]
        public void MolecularWeight_Dipeptide_CalculatesCorrectly()
        {
            var protein = new ProteinSequence("MK");
            // M (149.21) + K (146.19) - H2O (18.015) = 277.385 ≈ 277.38
            Assert.That(protein.MolecularWeight(), Is.EqualTo(277.39).Within(0.1));
        }

        [Test]
        public void MolecularWeight_Empty_ReturnsZero()
        {
            var protein = new ProteinSequence("");
            Assert.That(protein.MolecularWeight(), Is.EqualTo(0));
        }

        [Test]
        public void MolecularWeight_Insulin_ApproximatelyCorrect()
        {
            // Human insulin A chain: GIVEQCCTSICSLYQLENYCN (21 residues)
            var insulinA = new ProteinSequence("GIVEQCCTSICSLYQLENYCN");
            // Theoretical MW is approximately 2383 Da
            double mw = insulinA.MolecularWeight();
            Assert.That(mw, Is.GreaterThan(2000).And.LessThan(2600));
        }

        #endregion

        #region Isoelectric Point

        [Test]
        public void IsoelectricPoint_BasicProtein_HighpI()
        {
            // Protein with many basic residues (K, R)
            var protein = new ProteinSequence("KKKKRRRR");
            double pI = protein.IsoelectricPoint();
            Assert.That(pI, Is.GreaterThan(10));
        }

        [Test]
        public void IsoelectricPoint_AcidicProtein_LowpI()
        {
            // Protein with many acidic residues (D, E)
            var protein = new ProteinSequence("DDDDEEEE");
            double pI = protein.IsoelectricPoint();
            Assert.That(pI, Is.LessThan(4));
        }

        [Test]
        public void IsoelectricPoint_NeutralProtein_NeutralRange()
        {
            // Mostly neutral amino acids
            var protein = new ProteinSequence("GGGAAAVVV");
            double pI = protein.IsoelectricPoint();
            Assert.That(pI, Is.InRange(5, 7));
        }

        [Test]
        public void IsoelectricPoint_Empty_ReturnsZero()
        {
            var protein = new ProteinSequence("");
            Assert.That(protein.IsoelectricPoint(), Is.EqualTo(0));
        }

        #endregion

        #region GRAVY (Hydropathicity)

        [Test]
        public void Gravy_HydrophobicProtein_Positive()
        {
            // Hydrophobic amino acids (I, L, V, F)
            var protein = new ProteinSequence("ILVFILVF");
            Assert.That(protein.Gravy(), Is.GreaterThan(0));
        }

        [Test]
        public void Gravy_HydrophilicProtein_Negative()
        {
            // Hydrophilic amino acids (K, R, D, E)
            var protein = new ProteinSequence("KRDEKRDE");
            Assert.That(protein.Gravy(), Is.LessThan(0));
        }

        [Test]
        public void Gravy_Empty_ReturnsZero()
        {
            var protein = new ProteinSequence("");
            Assert.That(protein.Gravy(), Is.EqualTo(0));
        }

        #endregion

        #region Amino Acid Composition

        [Test]
        public void AminoAcidComposition_CountsCorrectly()
        {
            var protein = new ProteinSequence("MMMKKLLAA");
            var composition = protein.AminoAcidComposition();

            Assert.That(composition['M'], Is.EqualTo(3));
            Assert.That(composition['K'], Is.EqualTo(2));
            Assert.That(composition['L'], Is.EqualTo(2));
            Assert.That(composition['A'], Is.EqualTo(2));
        }

        [Test]
        public void AminoAcidComposition_AllDifferent_AllOnes()
        {
            var protein = new ProteinSequence("MKVL");
            var composition = protein.AminoAcidComposition();

            Assert.That(composition.Values, Is.All.EqualTo(1));
            Assert.That(composition.Count, Is.EqualTo(4));
        }

        #endregion

        #region Type Percentage

        [Test]
        public void TypePercentage_AllNonpolar_Returns100()
        {
            var protein = new ProteinSequence("AAGG");  // A and G are nonpolar
            Assert.That(protein.TypePercentage(AminoAcidType.Nonpolar), Is.EqualTo(100));
        }

        [Test]
        public void TypePercentage_HalfBasic_Returns50()
        {
            var protein = new ProteinSequence("KKAA");  // K is basic, A is nonpolar
            Assert.That(protein.TypePercentage(AminoAcidType.Basic), Is.EqualTo(50));
        }

        [Test]
        public void TypePercentage_NoAcidic_Returns0()
        {
            var protein = new ProteinSequence("MKVL");
            Assert.That(protein.TypePercentage(AminoAcidType.Acidic), Is.EqualTo(0));
        }

        #endregion

        #region Three Letter Code

        [Test]
        public void ToThreeLetterCode_SingleAminoAcid_ReturnsCode()
        {
            var protein = new ProteinSequence("M");
            Assert.That(protein.ToThreeLetterCode(), Is.EqualTo("Met"));
        }

        [Test]
        public void ToThreeLetterCode_MultipleAminoAcids_ReturnsDashSeparated()
        {
            var protein = new ProteinSequence("MKV");
            Assert.That(protein.ToThreeLetterCode(), Is.EqualTo("Met-Lys-Val"));
        }

        [Test]
        public void ToThreeLetterCode_WithStopCodon_IncludesTer()
        {
            var protein = new ProteinSequence("MK*");
            Assert.That(protein.ToThreeLetterCode(), Is.EqualTo("Met-Lys-Ter"));
        }

        [Test]
        public void ToThreeLetterCode_WithUnknown_IncludesXaa()
        {
            var protein = new ProteinSequence("MXK");
            Assert.That(protein.ToThreeLetterCode(), Is.EqualTo("Met-Xaa-Lys"));
        }

        #endregion

        #region Find Motif

        [Test]
        public void FindMotif_FoundOnce_ReturnsPosition()
        {
            var protein = new ProteinSequence("MKVLLCDE");
            var positions = protein.FindMotif("LLC").ToList();
            Assert.That(positions, Has.Count.EqualTo(1));
            Assert.That(positions[0], Is.EqualTo(3));
        }

        [Test]
        public void FindMotif_FoundMultiple_ReturnsAllPositions()
        {
            var protein = new ProteinSequence("MKVMKVMKV");
            var positions = protein.FindMotif("MKV").ToList();
            Assert.That(positions, Has.Count.EqualTo(3));
            Assert.That(positions, Is.EqualTo(new[] { 0, 3, 6 }));
        }

        [Test]
        public void FindMotif_NotFound_ReturnsEmpty()
        {
            var protein = new ProteinSequence("MKVLLCDE");
            var positions = protein.FindMotif("ZZZ").ToList();
            Assert.That(positions, Is.Empty);
        }

        [Test]
        public void FindMotif_EmptyPattern_ReturnsEmpty()
        {
            var protein = new ProteinSequence("MKVL");
            var positions = protein.FindMotif("").ToList();
            Assert.That(positions, Is.Empty);
        }

        #endregion

        #region Subsequence

        [Test]
        public void Subsequence_ReturnsCorrectSubsequence()
        {
            var protein = new ProteinSequence("MKVLACDE");
            var sub = protein.Subsequence(2, 4);
            Assert.That(sub.Sequence, Is.EqualTo("VLAC"));
        }

        #endregion

        #region Indexer

        [Test]
        public void Indexer_ReturnsCorrectAminoAcid()
        {
            var protein = new ProteinSequence("MKVL");
            Assert.That(protein[0], Is.EqualTo('M'));
            Assert.That(protein[1], Is.EqualTo('K'));
            Assert.That(protein[2], Is.EqualTo('V'));
            Assert.That(protein[3], Is.EqualTo('L'));
        }

        #endregion

        #region Equality

        [Test]
        public void Equals_SameSequence_ReturnsTrue()
        {
            var p1 = new ProteinSequence("MKVL");
            var p2 = new ProteinSequence("MKVL");
            Assert.That(p1.Equals(p2), Is.True);
        }

        [Test]
        public void Equals_DifferentSequence_ReturnsFalse()
        {
            var p1 = new ProteinSequence("MKVL");
            var p2 = new ProteinSequence("MKVK");
            Assert.That(p1.Equals(p2), Is.False);
        }

        [Test]
        public void GetHashCode_SameSequence_SameHash()
        {
            var p1 = new ProteinSequence("MKVL");
            var p2 = new ProteinSequence("MKVL");
            Assert.That(p1.GetHashCode(), Is.EqualTo(p2.GetHashCode()));
        }

        #endregion

        #region Static Properties

        [Test]
        public void StandardAminoAcids_Contains20AminoAcids()
        {
            Assert.That(ProteinSequence.StandardAminoAcids.Count, Is.EqualTo(20));
        }

        [Test]
        public void ValidCharacters_Contains22Characters()
        {
            // 20 standard + * + X
            Assert.That(ProteinSequence.ValidCharacters.Count, Is.EqualTo(22));
        }

        [Test]
        public void Properties_Contains20Entries()
        {
            Assert.That(ProteinSequence.Properties.Count, Is.EqualTo(20));
        }

        [Test]
        public void Properties_AllHaveValidData()
        {
            foreach (var kvp in ProteinSequence.Properties)
            {
                Assert.That(kvp.Value.Name, Is.Not.Empty);
                Assert.That(kvp.Value.ThreeLetterCode, Has.Length.EqualTo(3));
                Assert.That(kvp.Value.MolecularWeight, Is.GreaterThan(50));
            }
        }

        #endregion

        #region SuffixTree Integration

        [Test]
        public void SuffixTree_IsCreatedLazily()
        {
            var protein = new ProteinSequence("MKVLMKVL");
            var tree = protein.SuffixTree;
            Assert.That(tree, Is.Not.Null);
            Assert.That(tree.Contains("MKVL"), Is.True);
            Assert.That(tree.Contains("VLMK"), Is.True);
        }

        #endregion
    }
}
