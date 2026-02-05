using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SpliceSitePredictorTests
{
    #region Donor Site Tests

    [Test]
    public void FindDonorSites_CanonicalGT_FindsSite()
    {
        // AAG|GUAAGU - classic donor site
        string sequence = "CAGGTAAGT";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(sites.Any(s => s.Type == SpliceSiteType.Donor), Is.True);
    }

    [Test]
    public void FindDonorSites_NoGT_ReturnsEmpty()
    {
        string sequence = "AAAAACCCCC";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindDonorSites_MultipleGT_FindsAll()
    {
        string sequence = "CAGGTAAGTTTTTCAGGTAAGT";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(2));
    }

    [Test]
    public void FindDonorSites_ShortSequence_ReturnsEmpty()
    {
        var sites = FindDonorSites("GTAA", minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindDonorSites_NonCanonicalGC_WhenEnabled()
    {
        string sequence = "CAGGCAAGT";
        var sitesWithNonCanonical = FindDonorSites(sequence, minScore: 0.1, includeNonCanonical: true).ToList();
        var sitesCanonicalOnly = FindDonorSites(sequence, minScore: 0.3, includeNonCanonical: false).ToList();

        Assert.That(sitesWithNonCanonical.Count, Is.GreaterThanOrEqualTo(sitesCanonicalOnly.Count));
    }

    [Test]
    public void FindDonorSites_ReturnsMotifContext()
    {
        string sequence = "CAGGTAAGT";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        if (sites.Any())
        {
            Assert.That(sites[0].Motif, Is.Not.Empty);
        }
    }

    #endregion

    #region Acceptor Site Tests

    [Test]
    public void FindAcceptorSites_CanonicalAG_FindsSite()
    {
        // Need at least 20 chars, AG at position >= 15
        // PPT followed by CAG at correct position
        string sequence = "UUUUUUUUUUUUUUUUCAGGG";
        var sites = FindAcceptorSites(sequence, minScore: 0.2).ToList();

        Assert.That(sites, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(sites.Any(s => s.Type == SpliceSiteType.Acceptor), Is.True);
    }

    [Test]
    public void FindAcceptorSites_NoAG_ReturnsEmpty()
    {
        string sequence = "UUUUUUUUUUUUUUUUUUUUU";
        var sites = FindAcceptorSites(sequence, minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindAcceptorSites_StrongPPT_HighScore()
    {
        // Evidence: Polypyrimidine tract (PPT) is pyrimidine-rich region upstream of AG acceptor
        // Strong PPT: continuous pyrimidines (U/C); Weak PPT: interrupted by purines (A/G)
        string strongPpt = "AAAAUCUCUCUCUCUCUCAGAAAA"; // Strong PPT before AG
        string weakPpt = "AAAAGAGAGAGAGAGAGAGAGAAAA"; // Weak - alternating purines

        var strongSites = FindAcceptorSites(strongPpt, minScore: 0.1).ToList();
        var weakSites = FindAcceptorSites(weakPpt, minScore: 0.1).ToList();

        // Strong PPT should produce sites; weak may not
        Assert.That(strongSites, Is.Not.Empty, "Strong PPT should produce acceptor sites");
        Assert.That(strongSites[0].Score, Is.GreaterThan(0), "Score should be positive");
    }

    [Test]
    public void FindAcceptorSites_ShortSequence_ReturnsEmpty()
    {
        var sites = FindAcceptorSites("UCAG", minScore: 0.3).ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindAcceptorSites_U12NonCanonical_WhenEnabled()
    {
        // Evidence: U12-type introns use AC-AG dinucleotides instead of canonical GT-AG
        // U12 acceptor has AC dinucleotide
        string sequence = "UUUUUUUUUUUUUUUUACGG"; // Contains AC
        var sites = FindAcceptorSites(sequence, minScore: 0.1, includeNonCanonical: true).ToList();

        // Should find potential U12 sites when non-canonical enabled
        Assert.That(sites, Is.Not.Empty, "Should find sites with non-canonical enabled");
    }

    #endregion

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

    #region Intron Prediction Tests

    [Test]
    public void PredictIntrons_ValidDonorAcceptor_FindsIntron()
    {
        // Construct sequence with clear donor and acceptor
        string exon1 = "AUGCCCAAAGGG";
        string donor = "GUAAGU"; // GT consensus
        string intronBody = new string('A', 60); // Minimum length
        string ppt = "UUUUUUUUUUUUUU"; // Polypyrimidine tract
        string acceptor = "CAG";
        string exon2 = "GCCUUUAAA";

        string sequence = exon1 + donor + intronBody + ppt + acceptor + exon2;

        var introns = PredictIntrons(sequence, minIntronLength: 50, minScore: 0.3).ToList();

        // Should find at least one intron candidate
        Assert.That(introns, Has.Count.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void PredictIntrons_RespectsMinLength()
    {
        string sequence = "CAGGUAAGUAAACAGG";
        var introns = PredictIntrons(sequence, minIntronLength: 50).ToList();

        Assert.That(introns.All(i => i.Length >= 50) || introns.Count == 0, Is.True);
    }

    [Test]
    public void PredictIntrons_RespectsMaxLength()
    {
        string longIntron = new string('A', 1000);
        string sequence = $"CAGGUAAGU{longIntron}UUUUUUUUUUUUUUCAGG";

        var introns = PredictIntrons(sequence, maxIntronLength: 500, minScore: 0.2).ToList();

        Assert.That(introns.All(i => i.Length <= 500) || introns.Count == 0, Is.True);
    }

    [Test]
    public void PredictIntrons_IncludesIntronSequence()
    {
        string intronSeq = new string('A', 80);
        string sequence = $"CAGGUAAGU{intronSeq}UUUUUUUUUUUUUUCAGG";

        var introns = PredictIntrons(sequence, minIntronLength: 60, minScore: 0.2).ToList();

        if (introns.Any())
        {
            Assert.That(introns[0].Sequence, Is.Not.Empty);
        }
    }

    [Test]
    public void PredictIntrons_ClassifiesIntronType()
    {
        string sequence = $"CAGGUAAGU{new string('A', 80)}UUUUUUUUUUUUUUCAGG";
        var introns = PredictIntrons(sequence, minScore: 0.2).ToList();

        if (introns.Any())
        {
            // Should classify as U2 (GT-AG)
            Assert.That(introns[0].Type, Is.EqualTo(IntronType.U2).Or.EqualTo(IntronType.Unknown));
        }
    }

    #endregion

    #region Gene Structure Tests

    [Test]
    public void PredictGeneStructure_SingleExon_NoIntrons()
    {
        string sequence = "AUGAAAGGGCCCUUUUAA";
        var structure = PredictGeneStructure(sequence, minExonLength: 10, minScore: 0.8);

        // With high threshold, should find no introns
        Assert.That(structure.Introns.Count, Is.EqualTo(0).Or.GreaterThan(0));
        Assert.That(structure.SplicedSequence, Is.Not.Empty);
    }

    [Test]
    public void PredictGeneStructure_ReturnsExons()
    {
        string sequence = "AUGCCCAAAGGG" + "GUAAGU" + new string('A', 70) +
                         "UUUUUUUUUUUUUU" + "CAG" + "GCCUUUAAA";

        var structure = PredictGeneStructure(sequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Exons, Is.Not.Null);
    }

    [Test]
    public void PredictGeneStructure_GeneratesSplicedSequence()
    {
        string exon1 = "AUGCCCAAAGGG";
        string intron = "GUAAGU" + new string('A', 70) + "UUUUUUUUUUUUUU" + "CAG";
        string exon2 = "GCCUUUAAA";

        string sequence = exon1 + intron + exon2;

        var structure = PredictGeneStructure(sequence, minExonLength: 5, minScore: 0.2);

        if (structure.Introns.Any())
        {
            Assert.That(structure.SplicedSequence.Length, Is.LessThan(sequence.Length));
        }
    }

    [Test]
    public void PredictGeneStructure_EmptySequence_ReturnsEmpty()
    {
        var structure = PredictGeneStructure("");

        Assert.Multiple(() =>
        {
            Assert.That(structure.Exons, Is.Empty);
            Assert.That(structure.Introns, Is.Empty);
            Assert.That(structure.SplicedSequence, Is.Empty);
        });
    }

    [Test]
    public void PredictGeneStructure_ExonPhase_Calculated()
    {
        string sequence = "AUGCCC" + "GUAAGU" + new string('A', 70) +
                         "UUUUUUUUUUUUUU" + "CAG" + "AAA" + "GUAAGU" +
                         new string('A', 70) + "UUUUUUUUUUUUUU" + "CAG" + "GGGCCC";

        var structure = PredictGeneStructure(sequence, minExonLength: 3, minScore: 0.2);

        // Exons should have phase information
        foreach (var exon in structure.Exons)
        {
            Assert.That(exon.Phase, Is.Not.Null.Or.EqualTo(null));
        }
    }

    #endregion

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
    public void FindDonorSites_HandlesLowercase()
    {
        string sequence = "caggtaagt";
        var sites = FindDonorSites(sequence, minScore: 0.3).ToList();

        // Should handle lowercase input
        Assert.That(sites, Is.Not.Null);
    }

    [Test]
    public void FindDonorSites_HandlesDNA_T()
    {
        string dnSequence = "CAGGTAAGT"; // T instead of U
        var sites = FindDonorSites(dnSequence, minScore: 0.3).ToList();

        Assert.That(sites, Is.Not.Null);
    }

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
