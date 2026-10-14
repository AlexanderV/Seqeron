namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Review 2026-09, B02 F47 (A1-15): <see cref="ProteinSequence.TypePercentage"/> /
/// <see cref="ProteinSequence.TypeCount"/> follow the EMBOSS <c>pepstats</c>
/// "Property / Residues / Number / Mole%" table. Class membership from EMBOSS <c>Eamino.dat</c>
/// (U as C, O as K, J as I/L; X and '*' in no class); B/Z split to D/N and E/Q by
/// <c>embIepCompC</c>; denominator = pepstats "Residues" (one terminal '*' dropped by the reader).
/// Every expected value is the output of EMBOSS 6.6.0 <c>pepstats</c> (printed to 3 decimals).
/// </summary>
[TestFixture]
public class ProteinSequence_PepstatsClasses_Tests
{
    // Swiss-Prot LACI_ECOLI (P03023), 360 aa.
    private const string LacIEcoli =
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVATSSLALHAPSQIVAAIKSRADQLG" +
        "ASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIA" +
        "LLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDD" +
        "TEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLAPNTQTASPRALADSLMQLARQVSRLESGQ";

    private const string InsulinA = "GIVEQCCTSICSLYQLENYCN";

    // B×3 → D2 N1, Z×4 → E2 Q2; X, J, U, O; internal and terminal '*'. pepstats: Residues = 15.
    private const string Mixed = "MBBBZZZZXJUO*KH*";

    private static readonly AminoAcidType[] Order =
    {
        AminoAcidType.Tiny, AminoAcidType.Small, AminoAcidType.Aliphatic, AminoAcidType.Aromatic,
        AminoAcidType.Nonpolar, AminoAcidType.Polar, AminoAcidType.Charged, AminoAcidType.Basic,
        AminoAcidType.Acidic,
    };

    // Columns: Tiny, Small, Aliphatic, Aromatic, Non-polar, Polar, Charged, Basic, Acidic.
    private static IEnumerable<TestCaseData> PepstatsTables()
    {
        yield return new TestCaseData(InsulinA,
            new[] { 8, 11, 5, 2, 12, 9, 2, 0, 2 },
            new[] { 38.095, 52.381, 23.810, 9.524, 57.143, 42.857, 9.524, 0.000, 9.524 })
            .SetName("Pepstats_InsulinA");
        yield return new TestCaseData(LacIEcoli,
            new[] { 120, 197, 137, 21, 200, 160, 69, 37, 32 },
            new[] { 33.333, 54.722, 38.056, 5.833, 55.556, 44.444, 19.167, 10.278, 8.889 })
            .SetName("Pepstats_LacI");
        yield return new TestCaseData(Mixed,
            new[] { 1, 4, 1, 1, 3, 10, 7, 3, 4 },
            new[] { 6.667, 26.667, 6.667, 6.667, 20.000, 66.667, 46.667, 20.000, 26.667 })
            .SetName("Pepstats_AmbiguousSecPylStops");
        yield return new TestCaseData("MKUAOU",
            new[] { 3, 3, 1, 0, 4, 2, 2, 2, 0 },
            new[] { 50.000, 50.000, 16.667, 0.000, 66.667, 33.333, 33.333, 33.333, 0.000 })
            .SetName("Pepstats_SecPyl");
        yield return new TestCaseData("ACDEFBZXJUOKRH*",
            new[] { 3, 5, 2, 2, 5, 8, 8, 4, 4 },
            new[] { 21.429, 35.714, 14.286, 14.286, 35.714, 57.143, 57.143, 28.571, 28.571 })
            .SetName("Pepstats_SingleBZ_TerminalStop");
        yield return new TestCaseData("MKVL",
            new[] { 0, 1, 2, 0, 3, 1, 1, 1, 0 },
            new[] { 0.000, 25.000, 50.000, 0.000, 75.000, 25.000, 25.000, 25.000, 0.000 })
            .SetName("Pepstats_MKVL");
    }

    [TestCaseSource(nameof(PepstatsTables))]
    public void ResidueClasses_MatchPepstats(string sequence, int[] numbers, double[] molePercents)
    {
        var protein = new ProteinSequence(sequence);
        Assert.Multiple(() =>
        {
            for (int i = 0; i < Order.Length; i++)
            {
                Assert.That(protein.TypeCount(Order[i]), Is.EqualTo(numbers[i]), $"{Order[i]} Number");
                Assert.That(protein.TypePercentage(Order[i]), Is.EqualTo(molePercents[i]).Within(5e-4), $"{Order[i]} Mole%");
            }
        });
    }

    [TestCase("MKRHDE*", AminoAcidType.Basic, 50.000)]   // pepstats Residues = 6 (terminal '*' dropped)
    [TestCase("AK**", AminoAcidType.Tiny, 33.333)]       // Residues = 3 (only one terminal '*' dropped)
    [TestCase("*AK", AminoAcidType.Polar, 33.333)]       // Residues = 3 (leading '*' counted)
    [TestCase("AK*HD*GG", AminoAcidType.Tiny, 37.500)]   // Residues = 8 (internal '*' counted)
    [TestCase("AAXXHK", AminoAcidType.Aromatic, 16.667)] // Residues = 6 (X counted, in no class)
    public void Denominator_MatchesPepstatsResidues(string sequence, AminoAcidType type, double expected)
    {
        Assert.That(new ProteinSequence(sequence).TypePercentage(type), Is.EqualTo(expected).Within(5e-4));
    }

    [TestCase("BBZZAK", AminoAcidType.Polar, 5, 83.333)]    // B2 → D1 N1, Z2 → E1 Q1
    [TestCase("BBZZAK", AminoAcidType.Charged, 3, 50.000)]
    [TestCase("BBBZZZZ", AminoAcidType.Small, 3, 42.857)]   // B3 → D2 N1
    [TestCase("BBBZZZZ", AminoAcidType.Acidic, 4, 57.143)]  // D2 + E2
    public void AmbiguityCodes_SplitAsEmbossEmbIepCompC(string sequence, AminoAcidType type, int number, double expected)
    {
        var protein = new ProteinSequence(sequence);
        Assert.Multiple(() =>
        {
            Assert.That(protein.TypeCount(type), Is.EqualTo(number));
            Assert.That(protein.TypePercentage(type), Is.EqualTo(expected).Within(5e-4));
        });
    }

    [Test]
    public void TypePercentage_IsUnrounded()
    {
        // pepstats prints 38.095; the exact value 800/21 is returned (the old method rounded to 2 dp).
        Assert.That(new ProteinSequence(InsulinA).TypePercentage(AminoAcidType.Tiny),
            Is.EqualTo(800.0 / 21).Within(1e-12));
    }

    [TestCase("")]
    [TestCase("*")]
    public void TypePercentage_NoResidues_ReturnsZero(string sequence)
    {
        Assert.That(new ProteinSequence(sequence).TypePercentage(AminoAcidType.Polar), Is.EqualTo(0));
    }
}
