using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Review 2026-09, B02 F32 (TRANS-PROT-001; B03 request R4): selenocysteine (U, Sec) and pyrrolysine
/// (O, Pyl) are members of Biopython 1.88 <c>IUPACData.extended_protein_letters</c>
/// ("ACDEFGHIKLMNPQRSTVWYBXZJUO") with three-letter codes Sec/Pyl
/// (<c>protein_letters_1to3_extended</c>; <c>seq3("MKUAO*X")</c> = "MetLysSecAlaPylTerXaa"), but
/// <see cref="ProteinSequence"/> rejected them.
/// Reference values (all reproduced, see B02.md F32):
/// <list type="bullet">
/// <item>MW: Biopython <c>molecular_weight(s, "protein")</c> = <c>ProteinAnalysis(s).molecular_weight()</c>:
///   MKUAO 735.7975, MSKU 514.4988, MSKO 601.759, ACDUEFGOHK 1293.3102, U 168.0532, O 255.3134
///   (EMBOSS 6.6 pepstats: 735.80 / 514.50 / 601.76 / 1293.32).</item>
/// <item>pI (EMBOSS Epk default): EMBOSS 6.6 iep treats U/O as neutral — MSKU = MSKO = MSK = 9.1500,
///   MKUAO = MKA = 9.1500, ACDUEFGOHK = ACDEFGHK = 5.4305.</item>
/// <item>pI / charge (Bjellqvist): Biopython <c>IsoelectricPoint</c> treats U/O as neutral —
///   pi(MKUAO) = pi(MKA) = 8.500108…, charge_at_pH(7.0) = 0.49935568654234397 for both.</item>
/// <item>GRAVY: Biopython <c>ProteinAnalysis("MKUAO").gravy()</c> raises <c>KeyError: 'U'</c> (no
///   Kyte–Doolittle value); Seqeron skips residues absent from the scale exactly as for X/B/Z/J, so
///   GRAVY(MKUAO) = Biopython gravy("MKA") = −0.06666666666666665.</item>
/// </list>
/// </summary>
[TestFixture]
public class ProteinSequence_SecPyl_Tests
{
    [TestCase("U")]
    [TestCase("O")]
    [TestCase("MKUAO")]
    [TestCase("mkuao")]
    [TestCase("ACDUEFGOHK")]
    [Description("F32: U (Sec) and O (Pyl) are accepted (Biopython extended_protein_letters)")]
    public void Constructor_SecPyl_Accepted(string input)
    {
        Assert.That(ProteinSequence.TryCreate(input, out var protein), Is.True);
        Assert.That(protein!.Sequence, Is.EqualTo(input.ToUpperInvariant()));
    }

    [Test]
    [Description("F32: three-letter codes as Biopython seq3 / protein_letters_1to3_extended")]
    public void ToThreeLetterCode_SecPyl_AsBiopython()
    {
        // Biopython seq3("MKUAO*X") = "MetLysSecAlaPylTerXaa".
        Assert.That(new ProteinSequence("MKUAO*X").ToThreeLetterCode(),
            Is.EqualTo("Met-Lys-Sec-Ala-Pyl-Ter-Xaa"));
        Assert.That(ProteinSequence.NonStandardAminoAcids['U'], Is.EqualTo("Sec"));
        Assert.That(ProteinSequence.NonStandardAminoAcids['O'], Is.EqualTo("Pyl"));
    }

    [TestCase("U", 168.0532)]
    [TestCase("O", 255.3134)]
    [TestCase("MKUAO", 735.7975)]
    [TestCase("MSKU", 514.4988)]
    [TestCase("MSKO", 601.759)]
    [TestCase("ACDUEFGOHK", 1293.3102)]
    [Description("F32: MW = Biopython molecular_weight(seq, 'protein') incl. Sec/Pyl masses")]
    public void MolecularWeight_SecPyl_MatchesBiopython(string input, double expected)
    {
        Assert.That(new ProteinSequence(input).MolecularWeight(), Is.EqualTo(expected).Within(1e-9));
    }

    [TestCase("MSKU", 9.15)]
    [TestCase("MSKO", 9.15)]
    [TestCase("MKUAO", 9.15)]
    [TestCase("ACDUEFGOHK", 5.43)]
    [Description("F32: EMBOSS iep pI — U/O neutral (iep: 9.1500 / 9.1500 / 9.1500 / 5.4305)")]
    public void IsoelectricPoint_SecPyl_MatchesEmbossIep(string input, double expected)
    {
        Assert.That(new ProteinSequence(input).IsoelectricPoint(), Is.EqualTo(expected));
    }

    [Test]
    [Description("F32: U/O neutral on both pK scales — identical to the sequence without them")]
    public void NetChargeAndPi_SecPyl_NeutralLikeReferences()
    {
        Assert.Multiple(() =>
        {
            // EMBOSS scale: MKUAO behaves as MKA.
            Assert.That(ProteinPhysicochemistry.NetCharge("MKUAO", 7.0),
                Is.EqualTo(ProteinPhysicochemistry.NetCharge("MKA", 7.0)).Within(1e-12));
            // Bjellqvist scale: Biopython IsoelectricPoint("MKUAO").charge_at_pH(7.0) / pi().
            Assert.That(ProteinPhysicochemistry.NetCharge("MKUAO", 7.0, ProteinPkaScale.Bjellqvist),
                Is.EqualTo(0.49935568654234397).Within(1e-9));
            Assert.That(ProteinPhysicochemistry.IsoelectricPoint("MKUAO", ProteinPkaScale.Bjellqvist),
                Is.EqualTo(8.50));
        });
    }

    [Test]
    [Description("F32: GRAVY skips U/O (no Kyte-Doolittle value) like X/B/Z/J")]
    public void Gravy_SecPyl_SkippedLikeAmbiguityCodes()
    {
        // Biopython ProteinAnalysis("MKUAO").gravy() raises KeyError 'U'; gravy("MKA") = -0.06666666666666665.
        Assert.That(new ProteinSequence("MKUAO").Gravy(), Is.EqualTo(-0.06666666666666665).Within(1e-12));
    }

    [Test]
    [Description("F32: composition counts U/O; F47: TypePercentage per EMBOSS pepstats, where O (Pyl) is Basic (Eamino.dat)")]
    public void CompositionAndTypePercentage_SecPyl()
    {
        var protein = new ProteinSequence("MKUAOU");
        var composition = protein.AminoAcidComposition();
        Assert.Multiple(() =>
        {
            Assert.That(composition['U'], Is.EqualTo(2));
            Assert.That(composition['O'], Is.EqualTo(1));
            Assert.That(composition.Values.Sum(), Is.EqualTo(6));
            // EMBOSS 6.6.0 pepstats "MKUAOU": Basic 2 / 33.333 (K + O); Tiny 3 / 50.000 (A + U×2).
            Assert.That(protein.TypePercentage(AminoAcidType.Basic), Is.EqualTo(33.333).Within(5e-4));
            Assert.That(protein.TypePercentage(AminoAcidType.Tiny), Is.EqualTo(50.000).Within(5e-4));
        });
    }

    [Test]
    [Description("F32: standard-table translation never emits U/O (UGA/UAG → '*')")]
    public void Translator_StopCodonsNotRecoded()
    {
        Assert.That(Translator.Translate(new DnaSequence("ATGTGATAG")).Sequence, Is.EqualTo("M**"));
    }
}
