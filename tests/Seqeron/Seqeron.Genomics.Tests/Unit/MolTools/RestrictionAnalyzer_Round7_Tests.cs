namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 7 (the confirmation audit of round 6) on the restriction slice: B08 A7-2, A7-3, A7-4 and O7-1.
/// Each test failed on the code the auditor read (cbcecfd3) or, for the mutation rows, on the mutant the auditor built;
/// expected values are exact arithmetic or the documented rules, not the output of the member under test.
/// TestSpec: TestSpecs/RESTR-FIND-001.md, RESTR-DIGEST-001.md ("Round 7").
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Round7_Tests
{
    private const string EcoRiLocus = "AAAGAATTCAAA";

    // ------------------------------------------------------------------------------------------------
    // A7-2: ReverseReadingCutsAlike compared the cut offsets in int. For an offset within `len` of int.MinValue, `len - rev` wraps and
    // can equal `fwd`, so the reverse reading of the site was taken for "the same cleavage" and dropped from the digest and the map.
    // Expected cuts are exact arithmetic on GAATTC at 3 of a 12-mer (site length 6), 2^31 = 8 (mod 12), checked outside the library:
    //   forward reading: top cut = 3 + fwd;           reverse reading: top cut = 3 + 6 - rev        (all modulo 12)
    // ------------------------------------------------------------------------------------------------

    public static IEnumerable<TestCaseData> WrapRows()
    {
        // fwd = Min+4, rev = Min+2: top cuts (3 + 4 + 4) = 11 and (9 - (Min+2)) = 7 + 2^31 = 7 + 8 = 3 (mod 12)
        yield return new TestCaseData(int.MinValue + 4, int.MinValue + 2, new[] { (3, 8, "GAATTCAA"), (11, 4, "AAAA") }).SetName("Min+4_Min+2");
        // blunt pair (Min+3, Min+3): top cuts (3 + 4 + 3) = 10 and (6 - Min) = 6 + 8 = 14 = 2 (mod 12)
        yield return new TestCaseData(int.MinValue + 3, int.MinValue + 3, new[] { (2, 8, "AGAATTCA"), (10, 4, "AAAA") }).SetName("Min+3_Min+3");
        // the control on the other side of the int range, where len - rev does not wrap: Max = 7 (mod 12); top cuts 3 + (Max-3) = 7 and 9 - (Max-1) = 3
        yield return new TestCaseData(int.MaxValue - 3, int.MaxValue - 1, new[] { (3, 4, "GAAT"), (7, 8, "TCAAAAAA") }).SetName("Max-3_Max-1");
    }

    [TestCaseSource(nameof(WrapRows))]
    public void CircularDigest_CustomEnzymeAtTheIntExtremes_CutsAtBothReadings(int forward, int reverse, (int Start, int Length, string Sequence)[] expected)
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", forward, reverse, Organism: "");
        var dna = new DnaSequence(EcoRiLocus);

        var fragments = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme })
            .OrderBy(f => f.StartPosition).Select(f => (f.StartPosition, f.Length, f.Sequence)).ToList();
        var rawFragments = RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Circular, new[] { enzyme }, null)
            .OrderBy(f => f.StartPosition).Select(f => (f.StartPosition, f.Length, f.Sequence)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Is.EqualTo(expected.OrderBy(e => e.Start).ToList()), "DnaSequence overload");
            Assert.That(rawFragments, Is.EqualTo(fragments), "raw-string overload");
        });
    }

    [TestCaseSource(nameof(WrapRows))]
    public void CircularMap_AndPartialDigest_CustomEnzymeAtTheIntExtremes_CountTwoCleavages(int forward, int reverse, (int Start, int Length, string Sequence)[] expected)
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", forward, reverse, Organism: "");
        var dna = new DnaSequence(EcoRiLocus);

        var map = RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme });
        var partial = RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Circular, (int?)null, new[] { enzyme }).ToList();
        var summary = RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Circular, new[] { enzyme });

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2), "two cleavages: the reverse reading cuts other bonds");
            Assert.That(map.UniqueCutters, Is.Empty, "an enzyme that cuts twice is not a unique cutter");
            Assert.That(partial, Has.Count.EqualTo(5), "k = 2 cuts on a circle: k^2 + 1 partial-digest fragments");
            Assert.That(summary.TotalFragments, Is.EqualTo(expected.Length));
            Assert.That(summary.FragmentSizes.OrderBy(x => x), Is.EqualTo(expected.Select(e => e.Length).OrderBy(x => x)));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A7-3: the other two fields of a caller-supplied enzyme that multiply the output - the recognition sequence and the name - are bounded
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void EnzymeLimits_ArePinnedLiterally()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.MaxCutStagger, Is.EqualTo(64));
            Assert.That(RestrictionAnalyzer.MaxRecognitionLength, Is.EqualTo(64));
            Assert.That(RestrictionAnalyzer.MaxEnzymeNameLength, Is.EqualTo(64));
            Assert.That(RestrictionAnalyzer.RebaseEnzymes.Values.Max(e => e.RecognitionSequence.Length), Is.EqualTo(15), "REBASE's longest site");
            Assert.That(RestrictionAnalyzer.RebaseEnzymes.Values.Max(e => e.Name.Length), Is.EqualTo(9), "REBASE's longest name");
        });
    }

    private static void AssertRefusedByEveryMemberThatTakesAnEnzyme(RestrictionEnzyme enzyme, string label)
    {
        var dna = new DnaSequence("AAAGAATTCAAAGAATTCAAAGGGGGGGGGGGG");
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Circular), $"{label}: FindSites circular");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Linear), $"{label}: FindSites linear");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme }), $"{label}: Digest");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna.Sequence, MoleculeTopology.Linear, new[] { enzyme }, null), $"{label}: Digest raw");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Circular, new[] { enzyme }), $"{label}: GetDigestSummary");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme }), $"{label}: CreateMap");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Circular, (int?)null, new[] { enzyme }), $"{label}: PartialDigest");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetOverhangSequence(enzyme), $"{label}: GetOverhangSequence");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(enzyme, ecoRi), $"{label}: AreCompatible first");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(ecoRi, enzyme), $"{label}: AreCompatible second");
        });
    }

    [Test]
    public void CustomEnzyme_RecognitionSequenceBeyondTheLimit_IsRefusedAtTheCall_ByEveryMember()
    {
        AssertRefusedByEveryMemberThatTakesAnEnzyme(
            new RestrictionEnzyme("Long", new string('N', RestrictionAnalyzer.MaxRecognitionLength + 1), 1, 5, Organism: ""), "65 N");
        AssertRefusedByEveryMemberThatTakesAnEnzyme(
            new RestrictionEnzyme("Huge", new string('N', 200_000), 1, 5, Organism: ""), "200 000 N (the request that exhausted memory)");
    }

    [Test]
    public void CustomEnzyme_NameBeyondTheLimit_IsRefusedAtTheCall_ByEveryMember()
    {
        AssertRefusedByEveryMemberThatTakesAnEnzyme(
            new RestrictionEnzyme(new string('X', RestrictionAnalyzer.MaxEnzymeNameLength + 1), "GAATTC", 1, 5, Organism: ""), "65-character name");
    }

    [Test]
    public void CustomEnzyme_RecognitionSequenceAndNameAtTheLimit_AreAccepted()
    {
        var enzyme = new RestrictionEnzyme(
            new string('X', RestrictionAnalyzer.MaxEnzymeNameLength),
            new string('N', RestrictionAnalyzer.MaxRecognitionLength), 32, 32, Organism: "");
        var dna = new DnaSequence(new string('A', 100));

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Linear).Count(), Is.EqualTo(2 * (100 - 64 + 1)), "an N pattern matches every window, on both strands");
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(enzyme), Is.Empty, "a blunt enzyme");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // O7-1: GetOverhangSequence(RestrictionEnzyme) and AreCompatible(RestrictionEnzyme, RestrictionEnzyme) now check the definition as FindSites does
    // ------------------------------------------------------------------------------------------------

    [TestCase("")]
    [TestCase("GAAXXC")]
    [TestCase("GAAUUC")]
    [TestCase("GAA TTC")]
    public void GetOverhangSequence_AndAreCompatible_RefuseAnEmptyOrNonIupacRecognitionSequence(string site)
    {
        var enzyme = new RestrictionEnzyme("Bad", site, 1, 5, Organism: "");
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetOverhangSequence(enzyme));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(enzyme, ecoRi));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(ecoRi, enzyme));
        });
    }

    [Test]
    public void GetOverhangSequence_AndAreCompatible_ReadALowerCaseRecognitionSequenceAsUpperCase()
    {
        var lower = new RestrictionEnzyme("EcoRIlower", "gaattc", 1, 5, Organism: "");
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(lower), Is.EqualTo("AATT"), "was \"aatt\"");
            Assert.That(RestrictionAnalyzer.AreCompatible(lower, ecoRi), Is.True, "was an ArgumentOutOfRangeException from the IUPAC table");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A7-4: the 3' branch of the overhang read-out with a base that has no complement (the surviving mutant M10)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void Digest_ThreePrimeDistalOverhangWindowHoldingAUracil_KeepsTheUracil()
    {
        // GAATTC at 2 of AAGAATTCAUCGAAAAAAA (n = 19). Cuts (10, 6): top strand at 2 + 10 = 12, bottom at 2 + 6 = 8, so the 3' overhang
        // window is bases 8..11 = A U C G. The left fragment ends with it as it reads (top strand); the right fragment carries it on the bottom
        // strand, i.e. reverse-complemented: reverse GCUA, complement C G U(no complement, kept) T = CGUT.
        var enzyme = new RestrictionEnzyme("ThreePrimeDistal", "GAATTC", 10, 6, Organism: "");

        var digest = RestrictionAnalyzer.Digest("AAGAATTCAUCGAAAAAAA", MoleculeTopology.Linear, new[] { enzyme }, null).ToList();

        Assert.That(digest, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(digest[0].RightOverhang, Is.EqualTo("AUCG"));
            Assert.That(digest[1].LeftOverhang, Is.EqualTo("CGUT"));
            Assert.That(digest[0].RightOverhangType, Is.EqualTo(OverhangType.ThreePrime));
        });
    }
}
