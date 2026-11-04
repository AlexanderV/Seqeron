// SEQ-MW-001 / SEQ-PI-001 / SEQ-HYDRO-001 / SEQ-AA-COMP-001 / SEQ-CODON-FREQ-001 — routing to the
// canonical implementations (review 2026-09, FIN-B03 WP1: audit rows A1, A2; B02 request / R16).
// Evidence: docs/Validation/review-2026-09/B03.md, F22.
// Reference runs (2026-10-14): EMBOSS 6.6.0 pepstats / iep and Biopython 1.88.
//   pepstats "MSKſ" -> "from 1 to 3", Residues = 3, MW 364.46, pI 9.1500 (identical to "MSK");
//   iep "MSKſ" = iep "MSK" = 9.1500; pepstats "AıI" -> "from 1 to 2" (ı dropped);
//   Biopython ProteinAnalysis("MSK").molecular_weight() = 364.4609, ProteinAnalysis("AI").gravy() = 3.15.
//   (Biopython upper-cases with Python's Unicode str.upper(), so it reads ſ as S and ı as I —
//   MW("MSKſ") = 451.5382, gravy("AıI") = 3.6; the library follows EMBOSS / B02 TRANS-PROT-001:
//   ASCII-only case folding, non-ASCII letters are non-residues.)

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_ProteinPhysicochemistryRouting_Tests
{
    private const char LongS = 'ſ';     // ſ LATIN SMALL LETTER LONG S (ToUpperInvariant -> 'S')
    private const char DotlessI = 'ı';  // ı LATIN SMALL LETTER DOTLESS I (ToUpperInvariant -> 'I')

    [Test]
    public void MolecularWeight_NonAsciiLongS_IsSkippedLikeEmboss()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("MSK" + LongS);

        Assert.Multiple(() =>
        {
            Assert.That(mw, Is.EqualTo(SequenceStatistics.CalculateMolecularWeight("MSK")),
                "ſ is not a residue (EMBOSS pepstats: Residues = 3)");
            Assert.That(mw, Is.EqualTo(364.4609).Within(1e-4), "Biopython MW(\"MSK\") = 364.4609; pepstats 364.46");
        });
    }

    [Test]
    public void IsoelectricPoint_NonAsciiLongS_EqualsEmbossIep()
    {
        Assert.That(SequenceStatistics.CalculateIsoelectricPoint("MSK" + LongS), Is.EqualTo(9.15),
            "EMBOSS iep \"MSKſ\" = iep \"MSK\" = 9.1500");
    }

    [Test]
    public void Gravy_NonAsciiDotlessI_IsNotIsoleucine()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateHydrophobicity("A" + DotlessI + "I"), Is.EqualTo(3.15).Within(1e-12),
                "ı dropped (pepstats 'from 1 to 2'): GRAVY = Biopython gravy(\"AI\") = 3.15");
            Assert.That(SequenceStatistics.CalculateHydrophobicity(DotlessI.ToString()), Is.EqualTo(0.0),
                "no recognised residue -> 0");
        });
    }

    [Test]
    public void HydrophobicityProfile_NonAsciiLetters_HaveNoScaleValue_LowercaseAsciiFolds()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile(DotlessI + "a" + LongS + "i", windowSize: 1).ToArray();

        Assert.That(profile, Is.EqualTo(new[] { 0.0, 1.8, 0.0, 4.5 }),
            "ı/ſ contribute 0 (unknown residue); a -> A (1.8), i -> I (4.5)");
    }

    [Test]
    public void AminoAcidComposition_NonAsciiLetters_AreNotResidues()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("MSK" + LongS + DotlessI);

        Assert.Multiple(() =>
        {
            Assert.That(c.Length, Is.EqualTo(3), "pepstats \"MSKſ\": Residues = 3");
            Assert.That(c.Counts, Is.EquivalentTo(new Dictionary<char, int> { ['M'] = 1, ['S'] = 1, ['K'] = 1 }));
            Assert.That(c.MolecularWeight, Is.EqualTo(364.4609).Within(1e-4));
        });
    }

    /// <summary>
    /// Random differential (2000 cases, fixed seed): the public SequenceStatistics protein API is
    /// bit-identical to the canonical ProteinPhysicochemistry on inputs mixing upper/lower case,
    /// ambiguity codes (B/Z/X/U/O), junk and the non-ASCII letters ſ/ı.
    /// </summary>
    [Test]
    public void ProteinApi_RandomInputs_BitIdenticalToProteinPhysicochemistry()
    {
        const string alphabet = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwyBZXUOJbzxuo*-1 ſı";
        var rng = new Random(20261014);
        for (int iter = 0; iter < 2000; iter++)
        {
            int len = rng.Next(1, 60);
            var chars = new char[len];
            for (int i = 0; i < len; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
            string s = new(chars);
            double pH = rng.NextDouble() * 14.0;

            Assert.That(SequenceStatistics.CalculateMolecularWeight(s), Is.EqualTo(ProteinPhysicochemistry.MolecularWeight(s)), s);
            Assert.That(SequenceStatistics.CalculateHydrophobicity(s), Is.EqualTo(ProteinPhysicochemistry.Gravy(s)), s);
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint(s),
                Is.EqualTo(ProteinPhysicochemistry.IsoelectricPoint(s, ProteinPkaScale.Emboss)), s);
            Assert.That(SequenceStatistics.CalculateIsoelectricPoint(s, SequenceStatistics.PkaScale.Bjellqvist),
                Is.EqualTo(ProteinPhysicochemistry.IsoelectricPoint(s, ProteinPkaScale.Bjellqvist)), s);
            Assert.That(SequenceStatistics.CalculateNetCharge(s, pH, SequenceStatistics.PkaScale.Emboss),
                Is.EqualTo(ProteinPhysicochemistry.NetCharge(s, pH, ProteinPkaScale.Emboss)), s);
            Assert.That(SequenceStatistics.CalculateNetCharge(s, pH, SequenceStatistics.PkaScale.Bjellqvist),
                Is.EqualTo(ProteinPhysicochemistry.NetCharge(s, pH, ProteinPkaScale.Bjellqvist)), s);
        }
    }

    [Test]
    public void PkaScale_UndefinedValue_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceStatistics.CalculateIsoelectricPoint("ACD", (SequenceStatistics.PkaScale)99));
    }

    /// <summary>
    /// A2 (R16): CalculateCodonFrequencies is Translator.CountCodons normalised by the total
    /// (U→T, ASCII fold, non-ACGT triplets skipped without frame shift, trailing partial dropped).
    /// </summary>
    [Test]
    public void CodonFrequencies_RandomInputs_EqualTranslatorCountCodonsOverTotal()
    {
        const string alphabet = "ACGTUacgtuNRY-ſ";
        var rng = new Random(16);
        for (int iter = 0; iter < 500; iter++)
        {
            int len = rng.Next(0, 40);
            var chars = new char[len];
            for (int i = 0; i < len; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
            string s = new(chars);
            int frame = rng.Next(0, 4);

            var counts = Seqeron.Genomics.Core.Translator.CountCodons(s, frame);
            int total = counts.Values.Sum();
            var expected = counts.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);

            Assert.That(SequenceStatistics.CalculateCodonFrequencies(s, frame), Is.EquivalentTo(expected), $"{s} / frame {frame}");
        }
    }
}
