using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// Review 2026-09, B02 F32 (TRANS-PROT-001; B03 request R4): the protein tools validate through
/// <c>ProteinSequence.TryCreate</c>, which now accepts selenocysteine (U) and pyrrolysine (O)
/// (Biopython 1.88 <c>IUPACData.extended_protein_letters</c> = "ACDEFGHIKLMNPQRSTVWYBXZJUO").
/// Locked values: Biopython <c>molecular_weight(s, "protein")</c> MKUAO 735.7975, MSKU 514.4988;
/// EMBOSS 6.6 iep pI MSKU 9.1500 (U neutral, = MSK); Biopython <c>gravy("MKA")</c> −0.06666666666666665
/// (U/O skipped; Biopython raises KeyError on them).
/// </summary>
[TestFixture]
public class ProteinSecPylTests
{
    [TestCase("MKUAO")]
    [TestCase("mkuao")]
    [TestCase("U")]
    [TestCase("O")]
    public void ProteinValidate_SecPyl_Valid(string sequence)
    {
        var result = SequenceTools.ProteinValidate(sequence);
        Assert.That(result.Valid, Is.True);
        Assert.That(result.Error, Is.Null);
        Assert.That(result.Length, Is.EqualTo(sequence.Length));
    }

    [TestCase("MKUAO", 735.7975)]
    [TestCase("MSKU", 514.4988)]
    [TestCase("MSKO", 601.759)]
    public void MolecularWeightProtein_SecPyl_MatchesBiopython(string sequence, double expected)
    {
        Assert.That(SequenceTools.MolecularWeightProtein(sequence).MolecularWeight, Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void IsoelectricPoint_SecPyl_MatchesEmbossIep()
    {
        Assert.That(SequenceTools.IsoelectricPoint("MSKU").PI, Is.EqualTo(9.15));
        Assert.That(SequenceTools.IsoelectricPoint("ACDUEFGOHK").PI, Is.EqualTo(5.43));
    }

    [Test]
    public void Hydrophobicity_SecPyl_Skipped()
    {
        Assert.That(SequenceTools.Hydrophobicity("MKUAO").Gravy, Is.EqualTo(-0.06666666666666665).Within(1e-12));
    }

    [Test]
    public void AminoAcidComposition_SecPyl_CountedAndWeighed()
    {
        var result = SequenceTools.AminoAcidComposition("MKUAO");
        Assert.Multiple(() =>
        {
            Assert.That(result.Length, Is.EqualTo(5));
            Assert.That(result.Counts["U"], Is.EqualTo(1));
            Assert.That(result.Counts["O"], Is.EqualTo(1));
            Assert.That(result.MolecularWeight, Is.EqualTo(735.7975).Within(1e-9));
            Assert.That(result.IsoelectricPoint, Is.EqualTo(9.15));
        });
    }
}
