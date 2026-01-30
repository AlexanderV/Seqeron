using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

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
    // These tests are retained as smoke tests for MotifFinder API verification.

    [Test]
    [Description("Smoke: R (purine) pattern works via MotifFinder API")]
    public void FindDegenerateMotif_PurineR_MatchesAG()
    {
        var sequence = new DnaSequence("ATGCATGC");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "RTG").ToList();

        Assert.That(matches.Count, Is.EqualTo(2)); // ATG matches at 0 and 4
    }

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

    [Test]
    [Description("Smoke: W (weak) pattern works via MotifFinder API")]
    public void FindDegenerateMotif_WeakW_MatchesAT()
    {
        var sequence = new DnaSequence("ATATAT");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "WWW").ToList();

        Assert.That(matches.Count, Is.EqualTo(4)); // All AT combinations
    }

    [Test]
    [Description("Smoke: MatchedSequence property works correctly")]
    public void FindDegenerateMotif_ReturnsMatchedSequence()
    {
        var sequence = new DnaSequence("CAGCTG");
        var matches = MotifFinder.FindDegenerateMotif(sequence, "CANNTG").ToList();

        Assert.That(matches.Count, Is.EqualTo(1));
        Assert.That(matches[0].MatchedSequence, Is.EqualTo("CAGCTG"));
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

    #region Consensus Sequence Tests

    [Test]
    public void GenerateConsensus_IdenticalSequences_ReturnsSame()
    {
        var sequences = new[] { "ATGC", "ATGC", "ATGC" };
        string consensus = MotifFinder.GenerateConsensus(sequences);

        Assert.That(consensus, Is.EqualTo("ATGC"));
    }

    [Test]
    public void GenerateConsensus_MixedBases_ReturnsIupac()
    {
        var sequences = new[] { "ATGC", "GTGC" }; // First position: A and G = R
        string consensus = MotifFinder.GenerateConsensus(sequences);

        Assert.That(consensus[0], Is.EqualTo('R').Or.EqualTo('A').Or.EqualTo('G'));
    }

    [Test]
    public void GenerateConsensus_Empty_ReturnsEmpty()
    {
        var sequences = Array.Empty<string>();
        string consensus = MotifFinder.GenerateConsensus(sequences);

        Assert.That(consensus, Is.Empty);
    }

    [Test]
    public void GenerateConsensus_AllDifferent_ReturnsMostCommon()
    {
        // When bases are equally distributed, the algorithm returns the most common one
        // (alphabetically first if tied: A > C > G > T)
        var sequences = new[] { "AAAA", "TTTT", "GGGG", "CCCC" };
        string consensus = MotifFinder.GenerateConsensus(sequences);

        // With 25% threshold, all bases are exactly at threshold, so it falls back to most common
        Assert.That(consensus.Length, Is.EqualTo(4));
        Assert.That(consensus, Does.Match("^[ACGTN]+$")); // Valid DNA chars
    }

    #endregion

    #region Motif Discovery Tests

    [Test]
    public void DiscoverMotifs_FindsRepeatedKmer()
    {
        var sequence = new DnaSequence("ATGCATGCATGC");
        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        Assert.That(motifs.Any(m => m.Sequence == "ATGC"));
    }

    [Test]
    public void DiscoverMotifs_ReturnsPositions()
    {
        var sequence = new DnaSequence("ATGCATGC");
        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        var atgcMotif = motifs.First(m => m.Sequence == "ATGC");
        Assert.That(atgcMotif.Positions, Does.Contain(0));
        Assert.That(atgcMotif.Positions, Does.Contain(4));
    }

    [Test]
    public void DiscoverMotifs_CalculatesEnrichment()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 3, minCount: 1).ToList();

        var aaaMotif = motifs.First(m => m.Sequence == "AAA");
        Assert.That(aaaMotif.Enrichment, Is.GreaterThan(1));
    }

    [Test]
    public void DiscoverMotifs_FiltersByMinCount()
    {
        var sequence = new DnaSequence("ATGCAAAA");
        var motifs = MotifFinder.DiscoverMotifs(sequence, k: 4, minCount: 2).ToList();

        Assert.That(motifs.All(m => m.Count >= 2));
    }

    #endregion

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

    #region Regulatory Element Tests

    [Test]
    public void FindRegulatoryElements_FindsTataBox()
    {
        var sequence = new DnaSequence("ATATATAAAGCATGC");
        var elements = MotifFinder.FindRegulatoryElements(sequence).ToList();

        Assert.That(elements.Any(e => e.Name == "TATA Box"));
    }

    [Test]
    public void FindRegulatoryElements_FindsPolyASignal()
    {
        var sequence = new DnaSequence("ATGCAATAAAGCATGC");
        var elements = MotifFinder.FindRegulatoryElements(sequence).ToList();

        Assert.That(elements.Any(e => e.Name == "Poly(A) Signal"));
    }

    [Test]
    public void FindRegulatoryElements_FindsEBox()
    {
        var sequence = new DnaSequence("CAGCTG");
        var elements = MotifFinder.FindRegulatoryElements(sequence).ToList();

        Assert.That(elements.Any(e => e.Name == "E-box"));
    }

    [Test]
    public void FindRegulatoryElements_ReturnsDescription()
    {
        var sequence = new DnaSequence("TATAAA");
        var elements = MotifFinder.FindRegulatoryElements(sequence).ToList();

        Assert.That(elements.First().Description, Does.Contain("promoter"));
    }

    [Test]
    public void KnownMotifs_ContainsExpectedPatterns()
    {
        Assert.That(MotifFinder.KnownMotifs.TataBox, Is.EqualTo("TATAAA"));
        Assert.That(MotifFinder.KnownMotifs.ShineDalgarno, Is.EqualTo("AGGAGG"));
        Assert.That(MotifFinder.KnownMotifs.PolyASignal, Is.EqualTo("AATAAA"));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void FindExactMotif_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindExactMotif(null!, "ATG").ToList());
    }

    // NOTE: Canonical null test for FindDegenerateMotif is in IupacMotifMatchingTests.cs (PAT-IUPAC-001)
    [Test]
    [Description("Smoke: Null sequence throws ArgumentNullException")]
    public void FindDegenerateMotif_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindDegenerateMotif(null!, "ATG").ToList());
    }

    // NOTE: PWM null tests moved to MotifFinder_PWM_Tests.cs (PAT-PWM-001)

    [Test]
    public void DiscoverMotifs_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.DiscoverMotifs(null!).ToList());
    }

    [Test]
    public void DiscoverMotifs_ZeroK_ThrowsException()
    {
        var sequence = new DnaSequence("ATGC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            MotifFinder.DiscoverMotifs(sequence, k: 0).ToList());
    }

    [Test]
    public void FindSharedMotifs_NullSequences_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindSharedMotifs(null!).ToList());
    }

    [Test]
    public void FindRegulatoryElements_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.FindRegulatoryElements(null!).ToList());
    }

    [Test]
    public void GenerateConsensus_NullSequences_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.GenerateConsensus(null!));
    }

    #endregion
}
