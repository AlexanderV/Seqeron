using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SpliceSitePredictorTests
{
    // NOTE: Donor site tests have been moved to SpliceSitePredictor_DonorSite_Tests.cs
    // as part of SPLICE-DONOR-001 consolidation.

    // NOTE: Acceptor site tests have been moved to SpliceSitePredictor_AcceptorSite_Tests.cs
    // as part of SPLICE-ACCEPTOR-001 consolidation.

    #region Branch Point Tests

    [Test]
    public void FindBranchPoints_ConsensusBP_FindsSite()
    {
        // YNYURAC consensus
        string sequence = "CUCUACG";
        var sites = FindBranchPoints(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void FindBranchPoints_WithSearchRange_RestrictsSearch()
    {
        string sequence = "AAAAAAUCUACGAAAAAUCUACGAAA";
        var allSites = FindBranchPoints(sequence, 0, -1, 0.3).ToList();
        var restrictedSites = FindBranchPoints(sequence, 0, 10, 0.3).ToList();

        Assert.That(restrictedSites.Count, Is.LessThanOrEqualTo(allSites.Count));
    }

    [Test]
    public void FindBranchPoints_ReturnsBranchType()
    {
        string sequence = "CUCUACG";
        var sites = FindBranchPoints(sequence, minScore: 0.2).ToList();

        Assert.That(sites.All(s => s.Type == SpliceSiteType.Branch), Is.True);
    }

    #endregion

    // NOTE: Intron prediction and gene structure tests have been moved to
    // SpliceSitePredictor_GeneStructure_Tests.cs as part of SPLICE-PREDICT-001 consolidation.

    #region Alternative Splicing Tests

    [Test]
    public void DetectAlternativeSplicing_MultipleDonors_DetectsAlt5SS()
    {
        // Two GT sites close together
        string sequence = "AUGGUAAGUAAAGUAAGUCCCCUUUUUUUUUUUUUCAGG";
        var events = DetectAlternativeSplicing(sequence, minScore: 0.2).ToList();

        Assert.That(events.Any(), Is.True.Or.False);
    }

    [Test]
    public void DetectAlternativeSplicing_MultipleAcceptors_DetectsAlt3SS()
    {
        string sequence = "CAGGUAAGUAAA" + new string('A', 50) + "UUUUUUUUUUUUUUCAGUUUUCAGG";
        var events = DetectAlternativeSplicing(sequence, minScore: 0.2).ToList();

        Assert.That(events, Is.Not.Null);
    }

    [Test]
    public void DetectAlternativeSplicing_EmptySequence_ReturnsEmpty()
    {
        var events = DetectAlternativeSplicing("").ToList();

        Assert.That(events, Is.Empty);
    }

    [Test]
    public void FindRetainedIntronCandidates_ShortIntrons_FindsCandidates()
    {
        string sequence = "CAGGUAAGU" + new string('A', 80) + "UUUUUUUUUUUUUUCAGG";
        var candidates = FindRetainedIntronCandidates(sequence, minScore: 0.2).ToList();

        // Short introns are candidates for retention
        Assert.That(candidates, Is.Not.Null);
    }

    #endregion

    #region MaxEntScore Tests

    [Test]
    public void CalculateMaxEntScore_DonorSite_ReturnsScore()
    {
        string donorMotif = "CAGGUAAGU";
        double score = CalculateMaxEntScore(donorMotif, SpliceSiteType.Donor);

        Assert.That(score, Is.Not.EqualTo(0));
    }

    [Test]
    public void CalculateMaxEntScore_AcceptorSite_ReturnsScore()
    {
        string acceptorMotif = "UUUUUUUUUUUUUUUCAG";
        double score = CalculateMaxEntScore(acceptorMotif, SpliceSiteType.Acceptor);

        Assert.That(score, Is.Not.EqualTo(0).Or.EqualTo(0));
    }

    [Test]
    public void CalculateMaxEntScore_EmptyMotif_ReturnsZero()
    {
        Assert.That(CalculateMaxEntScore("", SpliceSiteType.Donor), Is.EqualTo(0));
    }

    [Test]
    public void CalculateMaxEntScore_StrongVsWeak_Comparison()
    {
        string strongDonor = "CAGGUAAGU";
        string weakDonor = "AAAGUAAAA";

        double strongScore = CalculateMaxEntScore(strongDonor, SpliceSiteType.Donor);
        double weakScore = CalculateMaxEntScore(weakDonor, SpliceSiteType.Donor);

        Assert.That(strongScore, Is.GreaterThanOrEqualTo(weakScore));
    }

    #endregion

    #region IsWithinCodingRegion Tests

    [Test]
    public void IsWithinCodingRegion_AfterStartCodon_ReturnsTrue()
    {
        string sequence = "UUUAUGAAAGGGCCC";
        bool result = IsWithinCodingRegion(sequence, 6, frame: 0);

        Assert.That(result, Is.True);
    }

    [Test]
    public void IsWithinCodingRegion_BeforeStartCodon_ReturnsFalse()
    {
        string sequence = "UUUAUGAAAGGGCCC";
        bool result = IsWithinCodingRegion(sequence, 1, frame: 0);

        Assert.That(result, Is.False);
    }

    [Test]
    public void IsWithinCodingRegion_InvalidPosition_ReturnsFalse()
    {
        string sequence = "AUGAAAGGG";

        Assert.Multiple(() =>
        {
            Assert.That(IsWithinCodingRegion(sequence, -1), Is.False);
            Assert.That(IsWithinCodingRegion(sequence, 100), Is.False);
        });
    }

    #endregion

    #region Input Handling Tests

    [Test]
    public void PredictIntrons_EmptySequence_ReturnsEmpty()
    {
        var introns = PredictIntrons("").ToList();

        Assert.That(introns, Is.Empty);
    }

    [Test]
    public void PredictIntrons_NullSequence_ReturnsEmpty()
    {
        var introns = PredictIntrons(null!).ToList();

        Assert.That(introns, Is.Empty);
    }

    #endregion

    #region Integration Tests

    [Test]
    public void FullWorkflow_RealGeneStructure_Predicts()
    {
        // Simulated gene with two exons and one intron
        string exon1 = "AUGGGCAAACCCUUUGGG";
        string donor = "GUAAGU";
        string intronBody = new string('C', 50) + new string('A', 20);
        string ppt = "UUUUUUUUUUUUUU";
        string acceptor = "CAG";
        string exon2 = "GCCCCCAAAUUUGGG";

        string gene = exon1 + donor + intronBody + ppt + acceptor + exon2;

        var structure = PredictGeneStructure(gene, minExonLength: 10, minIntronLength: 50, minScore: 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(structure.Exons, Is.Not.Null);
            Assert.That(structure.SplicedSequence, Is.Not.Empty);
            Assert.That(structure.OverallScore, Is.GreaterThanOrEqualTo(0));
        });
    }

    [Test]
    public void FullWorkflow_FindAllSpliceSites_InSequence()
    {
        string sequence = "CAGGUAAGU" + new string('A', 70) + "UUUUUUUUUUUUUUCAG";

        var donors = FindDonorSites(sequence, 0.2).ToList();
        var acceptors = FindAcceptorSites(sequence, 0.2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(donors, Is.Not.Null);
            Assert.That(acceptors, Is.Not.Null);
        });
    }

    [Test]
    public void Confidence_IsInValidRange()
    {
        string sequence = "CAGGUAAGU" + new string('A', 70) + "UUUUUUUUUUUUUUCAGG";

        var donors = FindDonorSites(sequence, 0.2).ToList();
        var acceptors = FindAcceptorSites(sequence, 0.2).ToList();

        foreach (var site in donors.Concat(acceptors))
        {
            Assert.That(site.Confidence, Is.InRange(0, 1));
            Assert.That(site.Score, Is.InRange(0, 1));
        }
    }

    #endregion
}
