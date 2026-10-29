namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 8 (the confirmation audit of round 7) on the restriction slice: B08 A8-1, A8-2, A8-3, A8-7, A8-8.
/// Each test failed on the code the auditor read (9b454d19) or, for the mutation rows, on the mutant the auditor built;
/// expected values are exact arithmetic or the documented rules, not the output of the member under test.
/// TestSpec: TestSpecs/RESTR-FIND-001.md, RESTR-DIGEST-001.md ("Round 8").
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Round8_Tests
{
    private const string EcoRiLocus = "AAAGAATTCAAA";

    // ------------------------------------------------------------------------------------------------
    // A8-1: the second cut pair of ReverseReadingCutsAlike (mutant "len as int" survived: the second-pair terms wrap again).
    // GAATTC at 3 of a 12-mer, cuts (5, 1) with the second pair (Min+5, Min+1) (stagger 4 for both), 2^31 = 8 (mod 12):
    //   forward reading top cuts   3 + 5 = 8   and 3 + (Min+5) = 3 + 4 + 5 = 12 = 0
    //   reverse reading top cuts   9 - 1 = 8    and 9 - (Min+1) = 8 + 8 = 16 = 4
    // The first pair is alike (5 = 6 - 1), the second is not, so the union of the boundaries is {0, 4, 8}.
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void CircularDigestAndMap_TwoSidedCustomEnzymeAtTheIntExtremes_CutAtTheUnionOfBothReadings()
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", 5, 1, Organism: "") { SecondCut = (int.MinValue + 5, int.MinValue + 1) };
        var dna = new DnaSequence(EcoRiLocus);

        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme })
            .OrderBy(f => f.StartPosition).Select(f => (f.StartPosition, f.Length, f.Sequence)).ToList();
        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme });

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Is.EqualTo(new[] { (0, 4, "AAAG"), (4, 4, "AATT"), (8, 4, "CAAA") }));
            Assert.That(map.TotalSites, Is.EqualTo(2));
            Assert.That(map.UniqueCutters, Is.Empty);
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A8-2: AreCompatible validates both definitions even when the ends decide the answer first (blunt x blunt, opposite polarity)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void AreCompatible_RefusesABadDefinition_EvenWhenTheEndsDecideTheAnswerFirst()
    {
        var eco = RestrictionAnalyzer.GetEnzyme("EcoRV")!;     // blunt
        var pst = RestrictionAnalyzer.GetEnzyme("PstI")!;      // 3' overhang
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;   // 5' overhang
        var badBlunt = new RestrictionEnzyme("BadBlunt", "GAAXXC", 3, 3, Organism: "");
        var badFivePrime = new RestrictionEnzyme("BadFive", "GAAXXC", 1, 5, Organism: "");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(badBlunt, eco), "blunt x blunt: True without the check");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(eco, badBlunt), "blunt x blunt, swapped");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(badFivePrime, pst), "5' x 3': False without the check");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(pst, badFivePrime), "3' x 5', swapped");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(badBlunt, ecoRi), "blunt x 5'");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A8-3: which check comes first, and which parameter is blamed
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void EnzymeLimits_AreCheckedBeforeThePattern_AndBlameTheCallersParameter()
    {
        // 65 bases that are also not IUPAC: the size limit is the reported reason (no work is spent on the characters of an oversized pattern)
        var oversized = new RestrictionEnzyme("Big", new string('X', 65), 1, 5, Organism: "");
        var longName = new RestrictionEnzyme(new string('Y', 65), "GAAXXC", 1, 5, Organism: "");
        var dna = new DnaSequence(EcoRiLocus);
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.Multiple(() =>
        {
            var e1 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, oversized, MoleculeTopology.Linear));
            Assert.That(e1!.ParamName, Is.EqualTo("enzyme"));
            Assert.That(e1.Message, Does.Contain("longest supported"));

            var e2 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, new[] { oversized }));
            Assert.That(e2!.ParamName, Is.EqualTo("enzymes"));
            Assert.That(e2.Message, Does.Contain("longest supported"));

            var e3 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetOverhangSequence(oversized));
            Assert.That(e3!.ParamName, Is.EqualTo("enzyme"));

            var e4 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(oversized, ecoRi));
            Assert.That(e4!.ParamName, Is.EqualTo("enzyme1"));
            var e5 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(ecoRi, oversized));
            Assert.That(e5!.ParamName, Is.EqualTo("enzyme2"));

            var e6 = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, longName, MoleculeTopology.Linear));
            Assert.That(e6!.Message, Does.Contain("longest supported name"), "the name limit is reported before the bad pattern");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A8-7: the pattern of the raw scan is bounded like an enzyme's recognition sequence (CPU is n * m)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void FindPatternOccurrences_PatternBeyondTheLimit_IsRefusedAtTheCall()
    {
        var tooLong = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", new string('N', 65)));
        var huge = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences(new string('A', 100), new string('N', 200_000)));

        Assert.Multiple(() =>
        {
            Assert.That(tooLong!.ParamName, Is.EqualTo("recognitionPattern"));
            Assert.That(tooLong.Message, Does.Contain("65").And.Contain("64"));
            Assert.That(huge!.ParamName, Is.EqualTo("recognitionPattern"));
        });
    }

    [Test]
    public void FindPatternOccurrences_PatternAtTheLimit_Works()
    {
        string pattern = new('N', RestrictionAnalyzer.MaxRecognitionLength);

        Assert.That(RestrictionAnalyzer.FindPatternOccurrences(new string('A', 100), pattern).Count(), Is.EqualTo(100 - 64 + 1));
    }

    // ------------------------------------------------------------------------------------------------
    // A8-8: on a circle two cuts are the same bond when congruent modulo the length, so the reverse reading of a site that cuts the
    // same bonds is one cleavage whichever representative of the offsets is given. GAATTC at 3 of a 12-mer; (1, 5) cuts top 4 / bottom 8
    // on both readings (9 - 5 = 4, 9 - 1 = 8); (13, 17), (-11, -7), (25, 29) are congruent to it modulo 12.
    // ------------------------------------------------------------------------------------------------

    [TestCase(1, 5)]
    [TestCase(13, 17)]
    [TestCase(-11, -7)]
    [TestCase(25, 29)]
    [TestCase(1 + 12 * 100_000, 5 + 12 * 100_000)]
    public void CircularMap_OffsetsCongruentModuloTheLength_CountOneCleavage_AsTheDigestDoes(int forward, int reverse)
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", forward, reverse, Organism: "");
        var dna = new DnaSequence(EcoRiLocus);

        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme });
        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1), "one cut");
            Assert.That(fragments[0].StartPosition, Is.EqualTo(4));
            Assert.That(map.TotalSites, Is.EqualTo(1), "one cleavage: the reverse reading cuts the same two bonds");
            Assert.That(map.UniqueCutters, Is.EqualTo(new[] { "X" }));
        });
    }

    [TestCase(1, 3)]
    [TestCase(13, 15)]
    public void CircularMap_AnEnzymeWhoseReverseReadingCutsOtherBonds_StillCountsTwoCleavages(int forward, int reverse)
    {
        // (1, 3): forward top 4 / bottom 6; reverse top 9 - 3 = 6 / bottom 9 - 1 = 8 - different bonds. (13, 15) is congruent to it.
        var enzyme = new RestrictionEnzyme("X", "GAATTC", forward, reverse, Organism: "");

        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(EcoRiLocus), MoleculeTopology.Circular, new[] { enzyme });

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2));
            Assert.That(map.UniqueCutters, Is.Empty);
        });
    }

    [Test]
    public void LinearMap_OffsetsAreAbsolute_ThereIsNoModulo()
    {
        // On a linear molecule only the plain relation counts: (1, 5) is one cleavage; (13, 17) puts every cut outside the 12-mer (no site at all).
        var alike = new RestrictionEnzyme("A", "GAATTC", 1, 5, Organism: "");
        var outside = new RestrictionEnzyme("B", "GAATTC", 13, 17, Organism: "");
        var dna = new DnaSequence(EcoRiLocus);

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, new[] { alike }).TotalSites, Is.EqualTo(1));
            Assert.That(RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Linear, new[] { outside }).TotalSites, Is.EqualTo(0));
        });
    }
}
