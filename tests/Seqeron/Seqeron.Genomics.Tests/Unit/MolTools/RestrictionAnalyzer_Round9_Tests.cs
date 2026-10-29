namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 9 (the confirmation audit of round 8) on the restriction slice: B08 A9-1, A9-2, A9-3 - test strength of the
/// congruence rule (<c>ReverseReadingCutsAlike</c> compares modulo the length on a circle only) and of the pattern bound. No defect was
/// found in the code of round 8; each row kills a surviving mutant the auditor built. Expected values are exact arithmetic.
/// TestSpec: TestSpecs/RESTR-FIND-001.md ("Round 9").
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Round9_Tests
{
    private const string EcoRiLocus = "AAAGAATTCAAA";

    // ------------------------------------------------------------------------------------------------
    // A9-1: "linear stays absolute". GAATTC at 0 of GAATTCAAAAAA (n = 12), offsets (-3, -3): the forward reading cuts at -3, outside the molecule, so
    // its record does not exist; the reverse reading cuts at 0 + 6 + 3 = 9, a real bond. -3 and 9 are congruent modulo 12, but on a linear molecule
    // they are different places, so the reverse record is a cleavage of its own: digest 0+9 / 9+3, TotalSites 1. (A modulo rule applied to a linear
    // molecule took the single real cleavage for a duplicate of the missing forward one and lost it.)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void LinearMoleculeShorterThanTheOffsets_TheReverseReadingIsItsOwnCleavage()
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", -3, -3, Organism: "");
        var dna = new DnaSequence("GAATTCAAAAAA");

        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, new[] { enzyme })
            .Select(f => (f.StartPosition, f.Length, f.Sequence)).ToList();
        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, new[] { enzyme });

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Is.EqualTo(new[] { (0, 9, "GAATTCAAA"), (9, 3, "AAA") }));
            Assert.That(map.TotalSites, Is.EqualTo(1));
            Assert.That(map.UniqueCutters, Is.EqualTo(new[] { "X" }));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A9-2: the two-pair rule with offsets that are congruent but not equal. GAATTC at 3 of a 12-mer; the reverse reading of pair (f, r) cuts the top
    // strand at 3 + 6 - r and the bottom at 3 + 6 - f. A two-sided enzyme is one cleavage per window when the two readings cut the same SET of bonds.
    // ------------------------------------------------------------------------------------------------

    public static IEnumerable<TestCaseData> TwoPairRows()
    {
        // (1,5)+(7,11): forward top cuts 4 and 10; reverse top cuts 9-5 = 4 and 9-11 = -2 = 10 (mod 12): the same set, pair by pair
        yield return new TestCaseData(1, 5, 7, 11, new[] { (4, 6, "AATTCA"), (10, 6, "AAAAAG") }).SetName("1_5__7_11");
        // (13,17)+(7,11): 13 = 1 and 17 = 5 (mod 12): the same enzyme
        yield return new TestCaseData(13, 17, 7, 11, new[] { (4, 6, "AATTCA"), (10, 6, "AAAAAG") }).SetName("13_17__7_11");
        // (1,3)+(15,17): forward top cuts 4 and 3+15 = 18 = 6; reverse top cuts 9-3 = 6 and 9-17 = -8 = 4: the same set, crosswise (the mirror)
        yield return new TestCaseData(1, 3, 15, 17, new[] { (4, 2, "AA"), (6, 10, "TTCAAAAAAG") }).SetName("1_3__15_17");
    }

    [TestCaseSource(nameof(TwoPairRows))]
    public void CircularTwoSidedEnzyme_ReadingsThatCutTheSameBondsModuloTheLength_AreOneCleavage(
        int f1, int r1, int f2, int r2, (int Start, int Length, string Sequence)[] expected)
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", f1, r1, Organism: "") { SecondCut = (f2, r2) };
        var dna = new DnaSequence(EcoRiLocus);

        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme });
        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme })
            .OrderBy(f => f.StartPosition).Select(f => (f.StartPosition, f.Length, f.Sequence)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(1), "one cleavage per window: the reverse reading cuts the same two bonds");
            Assert.That(map.UniqueCutters, Is.EqualTo(new[] { "X" }));
            Assert.That(fragments, Is.EqualTo(expected.OrderBy(e => e.Start).ToList()));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A9-3: the length bound comes before the IUPAC scan (no work, and no echo of the whole pattern, for an oversized one); AreCompatible checks the
    // size limits before the pattern, as FindSites does
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void FindPatternOccurrences_OversizedPattern_IsRefusedBeforeTheAlphabetIsScanned()
    {
        // 64 valid codes then a bad character: the length is the reported reason, and the message does not echo the 2 000-character pattern
        string pattern = new string('N', 2000) + "X";

        var e = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", pattern));

        Assert.Multiple(() =>
        {
            Assert.That(e!.Message, Does.Contain("2001").And.Contain("longest supported"));
            Assert.That(e.Message, Does.Not.Contain("IUPAC"));
            Assert.That(e.Message.Length, Is.LessThan(300), "the pattern is not echoed");
        });
    }

    [Test]
    public void AreCompatible_NameLimitIsReported_BeforeABadPattern()
    {
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        var longNameBadPattern = new RestrictionEnzyme(new string('Y', 65), "GAAXXC", 1, 5, Organism: "");

        var e1 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(longNameBadPattern, ecoRi));
        var e2 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(ecoRi, longNameBadPattern));

        Assert.Multiple(() =>
        {
            Assert.That(e1!.Message, Does.Contain("longest supported name"));
            Assert.That(e2!.Message, Does.Contain("longest supported name"));
        });
    }
}
