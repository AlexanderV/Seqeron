namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinderTests
{
    #region Exact Motif Finding Tests

    [Test]
    public void FindExactMotif_SingleOccurrence_FindsIt()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        var positions = MotifFinder.FindExactMotif(sequence, "TGCA").ToList();

        Assert.That(positions.Count, Is.EqualTo(2));
        Assert.That(positions, Does.Contain(1));
        Assert.That(positions, Does.Contain(5));
    }

    [Test]
    public void FindExactMotif_NoOccurrence_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAAAAA");
        var positions = MotifFinder.FindExactMotif(sequence, "TGCA").ToList();

        Assert.That(positions, Is.Empty);
    }

    [Test]
    public void FindExactMotif_OverlappingMatches_FindsAll()
    {
        var sequence = new DnaSequence("AAAA");
        var positions = MotifFinder.FindExactMotif(sequence, "AA").ToList();

        Assert.That(positions.Count, Is.EqualTo(3));
    }

    [Test]
    public void FindExactMotif_EmptyMotif_ReturnsEmpty()
    {
        var sequence = new DnaSequence("ATGC");
        var positions = MotifFinder.FindExactMotif(sequence, "").ToList();

        Assert.That(positions, Is.Empty);
    }

    [Test]
    public void FindExactMotif_CaseInsensitive()
    {
        var sequence = new DnaSequence("ATGCATGC");
        var positions = MotifFinder.FindExactMotif(sequence, "atgc").ToList();

        Assert.That(positions.Count, Is.EqualTo(2));
    }

    #endregion

    #region Degenerate Motif Finding Tests (Smoke - comprehensive tests in IupacMotifMatchingTests)

    // NOTE: Comprehensive IUPAC degenerate motif tests are in IupacMotifMatchingTests.cs (PAT-IUPAC-001)

    [Test]
    [Description("Smoke: Y (pyrimidine) pattern works via MotifFinder API")]
    public void FindDegenerateMotif_PyrimidineY_MatchesCT()
    {
        var sequence = new DnaSequence("CATTAT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "YAT").ToList();

        Assert.That(matches.Count, Is.EqualTo(2)); // CAT and TAT
    }

    [Test]
    [Description("Smoke: N (any) pattern works via MotifFinder API")]
    public void FindDegenerateMotif_AnyN_MatchesAll()
    {
        var sequence = new DnaSequence("ATGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "NNG").ToList();

        Assert.That(matches.Count, Is.EqualTo(1));
        Assert.That(matches[0].MatchedSequence, Is.EqualTo("ATG"));
    }

    #endregion

    #region PWM Tests (Smoke - comprehensive tests in MotifFinder_PWM_Tests)

    // NOTE: Comprehensive PWM tests are in MotifFinder_PWM_Tests.cs (PAT-PWM-001)
    // These tests are retained as smoke tests for MotifFinder API verification.

    [Test]
    [Description("Smoke: CreatePwm returns valid matrix")]
    public void CreatePwm_Smoke_ReturnsValidMatrix()
    {
        var sequences = new[] { "ATGC" };
        var pwm = MotifFinder.CreatePwm(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(pwm.Length, Is.EqualTo(4));
            Assert.That(pwm.Consensus, Is.EqualTo("ATGC"));
        });
    }

    [Test]
    [Description("Smoke: ScanWithPwm finds trained sequence")]
    public void ScanWithPwm_Smoke_FindsMatch()
    {
        var sequences = new[] { "ATGC", "ATGC", "ATGC" };
        var pwm = MotifFinder.CreatePwm(sequences);

        var sequence = new DnaSequence("AAAATGCAAA");
        var matches = MotifFinder.ScanWithPwm(sequence, pwm, threshold: 0).ToList();

        Assert.That(matches.Any(m => m.MatchedSequence == "ATGC"));
    }

    #endregion

    // NOTE: GenerateConsensus tests are in the canonical file
    // MotifFinder_GenerateConsensus_Tests.cs (MOTIF-GENERATE-001).

    // NOTE: DiscoverMotifs tests are in the canonical file
    // MotifFinder_DiscoverMotifs_Tests.cs (MOTIF-DISCOVER-001).

    #region Shared Motif Tests

    [Test]
    public void FindSharedMotifs_FindsCommonKmer()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("TGCATGCA"),
            new DnaSequence("GCATGCAT")
        };

        var shared = MotifFinder.FindSharedMotifs(sequences, k: 4, minSequences: 3).ToList();

        Assert.That(shared.Count, Is.GreaterThan(0));
    }

    [Test]
    public void FindSharedMotifs_ReturnsPrevalence()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var shared = MotifFinder.FindSharedMotifs(sequences, k: 4, minSequences: 2).ToList();

        Assert.That(shared.First().Prevalence, Is.EqualTo(1.0));
    }

    [Test]
    public void FindSharedMotifs_FiltersNotShared()
    {
        var sequences = new[]
        {
            new DnaSequence("AAAA"),
            new DnaSequence("TTTT")
        };

        var shared = MotifFinder.FindSharedMotifs(sequences, k: 4, minSequences: 2).ToList();

        Assert.That(shared, Is.Empty);
    }

    #endregion

    // NOTE: FindRegulatoryElements / KnownMotifs tests moved to the canonical
    // MotifFinder_FindRegulatoryElements_Tests.cs (MOTIF-REGULATORY-001).

    #region Edge Cases

    [Test]
    public void FindExactMotif_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindExactMotif(null!, "ATG").ToList());
    }

    // NOTE: Canonical null test for FindDegenerateMotif is in IupacMotifMatchingTests.cs (PAT-IUPAC-001)

    // NOTE: PWM null tests moved to MotifFinder_PWM_Tests.cs (PAT-PWM-001)

    // NOTE: DiscoverMotifs null/k-range tests moved to the canonical
    // MotifFinder_DiscoverMotifs_Tests.cs (MOTIF-DISCOVER-001).

    [Test]
    public void FindSharedMotifs_NullSequences_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindSharedMotifs(null!).ToList());
    }

    #endregion
}
