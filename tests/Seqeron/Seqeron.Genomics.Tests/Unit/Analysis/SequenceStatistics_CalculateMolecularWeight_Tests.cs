// SEQ-MW-001 — Molecular Weight Calculation
// Evidence: docs/Evidence/SEQ-MW-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-MW-001.md
// Source: Expasy Compute pI/Mw (Gasteiger et al. 2005); Biopython Bio.SeqUtils.molecular_weight
//         and Bio.Data.IUPACData (average mass tables, water = 18.0153 Da).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateMolecularWeight_Tests
{
    // Expected values are exact average-mass derivations (Biopython IUPACData; water = 18.0153 Da).
    private const double Tolerance = 1e-4;

    #region CalculateMolecularWeight (protein)

    // M1 — protein "AGC": 89.0932 + 75.0666 + 121.1582 - 2*18.0153 = 249.2874 Da
    // Evidence: Biopython molecular_weight("AGC","protein") -> 249.29; Expasy Mw formula.
    [Test]
    public void CalculateMolecularWeight_TripeptideAGC_ReturnsExactAverageMass()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("AGC");

        Assert.That(mw, Is.EqualTo(249.2874).Within(Tolerance),
            "Protein Mw = sum of free amino-acid average masses minus (n-1) waters; Biopython docstring 249.29");
    }

    // M4 — single Gly: zero peptide bonds, so result is the free amino-acid mass 75.0666 Da
    // Evidence: Expasy formula (n-1=0); Biopython protein_weights["G"] = 75.0666.
    [Test]
    public void CalculateMolecularWeight_SingleGlycine_ReturnsFreeAminoAcidMass()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("G");

        Assert.That(mw, Is.EqualTo(75.0666).Within(Tolerance),
            "One residue forms zero peptide bonds, so Mw equals the free amino-acid mass of glycine");
    }

    // M7 — empty protein -> 0
    [Test]
    public void CalculateMolecularWeight_EmptyString_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("");

        Assert.That(mw, Is.EqualTo(0), "No monomers means molecular weight 0 (INV-01)");
    }

    // M8 — null protein -> 0
    [Test]
    public void CalculateMolecularWeight_Null_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight(null!);

        Assert.That(mw, Is.EqualTo(0), "Null input returns 0 without throwing (INV-01)");
    }

    // S1 — case-insensitive
    [Test]
    public void CalculateMolecularWeight_LowercaseInput_EqualsUppercase()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("agc");

        Assert.That(mw, Is.EqualTo(249.2874).Within(Tolerance),
            "Input is upper-cased before lookup, so lowercase yields the same Mw as 'AGC' (INV-04)");
    }

    // S3 — unknown symbol skipped (no mass, no bond): "AG*C" == "AGC"
    // Evidence: ASSUMPTION-02 (deviation from Biopython reject-on-unknown).
    [Test]
    public void CalculateMolecularWeight_UnknownSymbol_IsSkipped()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("AG*C");

        Assert.That(mw, Is.EqualTo(249.2874).Within(Tolerance),
            "Unrecognized '*' contributes no mass and no peptide bond, so Mw equals that of 'AGC'");
    }

    // C1 — bond-count invariant (protein): "AG" = mA + mG - one water
    // Evidence: INV-03; 89.0932 + 75.0666 - 18.0153 = 146.1445 Da.
    [Test]
    public void CalculateMolecularWeight_Dipeptide_RemovesExactlyOneWater()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("AG");

        Assert.That(mw, Is.EqualTo(146.1445).Within(Tolerance),
            "Two residues form exactly one peptide bond, removing one water (INV-03)");
    }

    // Branch: no recognized monomers -> 0 (distinct from null/empty short-circuit).
    // Evidence: spec §3.3 "A sequence containing no recognized monomers returns 0".
    [Test]
    public void CalculateMolecularWeight_AllUnknownSymbols_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("***");

        Assert.That(mw, Is.EqualTo(0),
            "A sequence with no recognized amino acids contributes no mass and returns 0 (spec §3.3)");
    }

    #endregion

    #region CalculateNucleotideMolecularWeight (DNA/RNA)

    // M2 — DNA "AGC": 331.2218 + 347.2212 + 307.1971 - 2*18.0153 = 949.6095 Da
    // Evidence: Biopython molecular_weight("AGC","DNA") -> 949.61.
    [Test]
    public void CalculateNucleotideMolecularWeight_DnaAGC_ReturnsExactAverageMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: true);

        Assert.That(mw, Is.EqualTo(949.6095).Within(Tolerance),
            "DNA Mw = sum of monophosphate masses minus (n-1) waters; Biopython docstring 949.61");
    }

    // M3 — RNA "AGC": 347.2212 + 363.2206 + 323.1965 - 2*18.0153 = 997.6077 Da
    // Evidence: Biopython molecular_weight("AGC","RNA") -> 997.61.
    [Test]
    public void CalculateNucleotideMolecularWeight_RnaAGC_ReturnsExactAverageMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("AGC", isDna: false);

        Assert.That(mw, Is.EqualTo(997.6077).Within(Tolerance),
            "RNA Mw = sum of monophosphate masses minus (n-1) waters; Biopython docstring 997.61");
    }

    // M5 — single DNA A: zero phosphodiester bonds -> monophosphate mass 331.2218 Da
    // Evidence: Biopython unambiguous_dna_weights["A"] = 331.2218.
    [Test]
    public void CalculateNucleotideMolecularWeight_SingleDnaA_ReturnsMonophosphateMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("A", isDna: true);

        Assert.That(mw, Is.EqualTo(331.2218).Within(Tolerance),
            "One nucleotide forms zero bonds, so Mw equals the dAMP monophosphate mass");
    }

    // M6 — single RNA A: monophosphate mass 347.2212 Da
    // Evidence: Biopython unambiguous_rna_weights["A"] = 347.2212.
    [Test]
    public void CalculateNucleotideMolecularWeight_SingleRnaA_ReturnsMonophosphateMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("A", isDna: false);

        Assert.That(mw, Is.EqualTo(347.2212).Within(Tolerance),
            "One ribonucleotide forms zero bonds, so Mw equals the AMP monophosphate mass");
    }

    // M9 — empty nucleotide -> 0
    [Test]
    public void CalculateNucleotideMolecularWeight_EmptyString_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("", isDna: true);

        Assert.That(mw, Is.EqualTo(0), "No monomers means molecular weight 0 (INV-01)");
    }

    // M10 — null nucleotide -> 0
    [Test]
    public void CalculateNucleotideMolecularWeight_Null_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight(null!, isDna: true);

        Assert.That(mw, Is.EqualTo(0), "Null input returns 0 without throwing (INV-01)");
    }

    // S2 — case-insensitive (DNA)
    [Test]
    public void CalculateNucleotideMolecularWeight_LowercaseInput_EqualsUppercase()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("agc", isDna: true);

        Assert.That(mw, Is.EqualTo(949.6095).Within(Tolerance),
            "Input is upper-cased before lookup, so lowercase yields the same Mw as 'AGC' (INV-04)");
    }

    // S4 — unknown nucleotide symbol skipped: "AG*C" == "AGC" (DNA)
    // Evidence: ASSUMPTION-02 (deviation from Biopython reject-on-unknown).
    [Test]
    public void CalculateNucleotideMolecularWeight_UnknownSymbol_IsSkipped()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("AG*C", isDna: true);

        Assert.That(mw, Is.EqualTo(949.6095).Within(Tolerance),
            "Unrecognized '*' contributes no mass and no bond, so Mw equals that of DNA 'AGC'");
    }

    // C2 — bond-count invariant (DNA): "AG" = mA + mG - one water
    // Evidence: INV-04; 331.2218 + 347.2212 - 18.0153 = 660.4277 Da.
    [Test]
    public void CalculateNucleotideMolecularWeight_Dinucleotide_RemovesExactlyOneWater()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("AG", isDna: true);

        Assert.That(mw, Is.EqualTo(660.4277).Within(Tolerance),
            "Two nucleotides form exactly one phosphodiester bond, removing one water (INV-04)");
    }

    // Full DNA alphabet, exercises the T table entry (untested by AGC).
    // Evidence: Biopython unambiguous_dna_weights {A 331.2218, C 307.1971, G 347.2212, T 322.2085};
    //           331.2218 + 307.1971 + 347.2212 + 322.2085 - 3*18.0153 = 1253.8027 Da.
    [Test]
    public void CalculateNucleotideMolecularWeight_DnaACGT_ReturnsExactAverageMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("ACGT", isDna: true);

        Assert.That(mw, Is.EqualTo(1253.8027).Within(Tolerance),
            "Sum of all four DNA monophosphate masses minus three waters (Biopython tables)");
    }

    // Full RNA alphabet, exercises the U table entry (untested by AGC).
    // Evidence: Biopython unambiguous_rna_weights {A 347.2212, C 323.1965, G 363.2206, U 324.1813};
    //           347.2212 + 323.1965 + 363.2206 + 324.1813 - 3*18.0153 = 1303.7737 Da.
    [Test]
    public void CalculateNucleotideMolecularWeight_RnaACGU_ReturnsExactAverageMass()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("ACGU", isDna: false);

        Assert.That(mw, Is.EqualTo(1303.7737).Within(Tolerance),
            "Sum of all four RNA monophosphate masses minus three waters (Biopython tables)");
    }

    // Branch: no recognized monomers -> 0 (distinct from null/empty short-circuit).
    // Evidence: spec §3.3 "A sequence containing no recognized monomers returns 0".
    [Test]
    public void CalculateNucleotideMolecularWeight_AllUnknownSymbols_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("***", isDna: true);

        Assert.That(mw, Is.EqualTo(0),
            "A sequence with no recognized nucleotides contributes no mass and returns 0 (spec §3.3)");
    }

    #endregion
}
