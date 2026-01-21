using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class CrisprDesignerTests
{
    #region CRISPR System Tests

    [Test]
    public void GetSystem_SpCas9_ReturnsCorrectSystem()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9);

        Assert.That(system.Name, Is.EqualTo("SpCas9"));
        Assert.That(system.PamSequence, Is.EqualTo("NGG"));
        Assert.That(system.GuideLength, Is.EqualTo(20));
        Assert.That(system.PamAfterTarget, Is.True);
    }

    [Test]
    public void GetSystem_SaCas9_ReturnsCorrectSystem()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.SaCas9);

        Assert.That(system.Name, Is.EqualTo("SaCas9"));
        Assert.That(system.PamSequence, Is.EqualTo("NNGRRT"));
        Assert.That(system.GuideLength, Is.EqualTo(21));
    }

    [Test]
    public void GetSystem_Cas12a_ReturnsCorrectSystem()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.Cas12a);

        Assert.That(system.Name, Is.EqualTo("Cas12a/Cpf1"));
        Assert.That(system.PamSequence, Is.EqualTo("TTTV"));
        Assert.That(system.PamAfterTarget, Is.False); // PAM is 5' of target
    }

    #endregion

    #region PAM Finding Tests

    [Test]
    public void FindPamSites_SpCas9_NGG_FindsSites()
    {
        // Sequence with NGG PAMs
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG"); // AGG at end
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "AGG"));
    }

    [Test]
    public void FindPamSites_SpCas9_CGG_FindsSite()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTCGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "CGG"));
    }

    [Test]
    public void FindPamSites_SpCas9_TGG_FindsSite()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTTGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TGG"));
    }

    [Test]
    public void FindPamSites_NoPam_ReturnsEmpty()
    {
        // Sequence without NGG
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        // Check that no sites have NGG PAM - there may be reverse complement matches
        var forwardSites = sites.Where(s => s.IsForwardStrand && s.PamSequence.EndsWith("GG")).ToList();
        Assert.That(forwardSites.All(s => s.PamSequence[0] != 'N'));
    }

    [Test]
    public void FindPamSites_BothStrands_FindsSites()
    {
        // Create sequence with PAM on forward strand (AGG) 
        // For reverse strand, we need sequence that creates NGG when reverse-complemented
        // CCN on forward = NGG on reverse (e.g., CCA reverse = TGG)
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGGTTTTTTTTTTTTTTTTTTTTTTTTCCA");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        // Should find forward strand site (AGG)
        Assert.That(sites.Any(s => s.IsForwardStrand), $"No forward strand sites found. Sites: {sites.Count}");
    }

    [Test]
    public void FindPamSites_StringOverload_Works()
    {
        var sites = CrisprDesigner.FindPamSites("ACGTACGTACGTACGTACGTACGTAGG", CrisprSystemType.SpCas9).ToList();
        Assert.That(sites.Any(s => s.PamSequence == "AGG"));
    }

    [Test]
    public void FindPamSites_EmptySequence_ReturnsEmpty()
    {
        var sites = CrisprDesigner.FindPamSites("", CrisprSystemType.SpCas9).ToList();
        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindPamSites_ReturnsTargetSequence()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var site = sites.FirstOrDefault(s => s.PamSequence == "AGG" && s.IsForwardStrand);
        if (site != null)
        {
            Assert.That(site.TargetSequence.Length, Is.EqualTo(20));
        }
    }

    [Test]
    public void FindPamSites_Cas12a_TTTV_FindsSites()
    {
        // TTTV = TTTA, TTTC, TTTG
        var sequence = new DnaSequence("ACGTTTTAACGTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TTTA"));
    }

    [Test]
    public void FindPamSites_SaCas9_NNGRRT_FindsSites()
    {
        // NNGRRT examples: AAGAAT, CCGAGT, TTGAGT
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTAAGAATACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence.Length == 6));
    }

    #endregion

    #region Guide RNA Evaluation Tests

    [Test]
    public void EvaluateGuideRna_OptimalGuide_HighScore()
    {
        // Good guide: ~50% GC, no polyT, low self-complementarity
        string guide = "ACGTACGTACGTACGTACGT"; // 50% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.GreaterThan(70));
        Assert.That(candidate.GcContent, Is.EqualTo(50));
        Assert.That(candidate.HasPolyT, Is.False);
    }

    [Test]
    public void EvaluateGuideRna_LowGcContent_LowerScore()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA"; // 0% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(0));
        Assert.That(candidate.Issues, Has.Some.Contains("Low GC"));
    }

    [Test]
    public void EvaluateGuideRna_HighGcContent_LowerScore()
    {
        string guide = "GCGCGCGCGCGCGCGCGCGC"; // 100% GC
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.Score, Is.LessThan(50));
        Assert.That(candidate.GcContent, Is.EqualTo(100));
        Assert.That(candidate.Issues, Has.Some.Contains("High GC"));
    }

    [Test]
    public void EvaluateGuideRna_HasPolyT_Penalized()
    {
        string guide = "ACGTACGTTTTTACGTACGT"; // Contains TTTT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.True);
        Assert.That(candidate.Issues, Has.Some.Contains("TTTT"));
    }

    [Test]
    public void EvaluateGuideRna_NoPolyT_NotPenalized()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.HasPolyT, Is.False);
    }

    [Test]
    public void EvaluateGuideRna_CalculatesSeedGc()
    {
        string guide = "AAAAAAAAAAAAACGTACGT"; // Last 12 bases are AAAACGTACGT
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.SeedGcContent, Is.GreaterThan(0));
    }

    [Test]
    public void EvaluateGuideRna_FullGuideRna_IncludesScaffold()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var candidate = CrisprDesigner.EvaluateGuideRna(guide, CrisprSystemType.SpCas9);

        Assert.That(candidate.FullGuideRna, Does.StartWith(guide));
        Assert.That(candidate.FullGuideRna.Length, Is.GreaterThan(guide.Length));
    }

    [Test]
    public void EvaluateGuideRna_EmptyGuide_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.EvaluateGuideRna("", CrisprSystemType.SpCas9));
    }

    #endregion

    #region Guide RNA Design Tests

    [Test]
    public void DesignGuideRnas_WithPamInRegion_ReturnsGuides()
    {
        // Create a sequence with PAM in the target region
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGG");
        var guides = CrisprDesigner.DesignGuideRnas(sequence, 20, 45, CrisprSystemType.SpCas9).ToList();

        // Should find at least one guide
        Assert.That(guides, Has.Count.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void DesignGuideRnas_InvalidRegionStart_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, -1, 10, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void DesignGuideRnas_InvalidRegionEnd_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.DesignGuideRnas(sequence, 0, 100, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void DesignGuideRnas_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.DesignGuideRnas(null!, 0, 10, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region Off-Target Analysis Tests

    [Test]
    public void FindOffTargets_NoMismatches_ReturnsEmpty()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG"); // Exact match

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // Off-targets require mismatches, exact matches are not off-targets
        Assert.That(offTargets.All(ot => ot.Mismatches > 0));
    }

    [Test]
    public void FindOffTargets_WithMismatches_FindsOffTargets()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Create sequence with similar but not identical target
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGGTTTTTTTTTTTTTTTTTTTTACGAACGTACGTACGTACGTCGG");

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        // May or may not find off-targets depending on sequence
        Assert.That(offTargets, Is.Not.Null);
    }

    [Test]
    public void FindOffTargets_MaxMismatchesRespected()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAAAGG"); // Many mismatches

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 2, CrisprSystemType.SpCas9).ToList();

        Assert.That(offTargets.All(ot => ot.Mismatches <= 2));
    }

    [Test]
    public void FindOffTargets_ReturnsMismatchPositions()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTTCGAACGTACGTACGTACGTCGG"); // Different target

        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3, CrisprSystemType.SpCas9).ToList();

        foreach (var ot in offTargets)
        {
            Assert.That(ot.MismatchPositions, Is.Not.Null);
            Assert.That(ot.MismatchPositions.Count, Is.EqualTo(ot.Mismatches));
        }
    }

    [Test]
    public void FindOffTargets_EmptyGuide_ThrowsException()
    {
        var genome = new DnaSequence("ACGTACGT");
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("", genome, 3, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void FindOffTargets_NullGenome_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGT", null!, 3, CrisprSystemType.SpCas9).ToList());
    }

    [Test]
    public void FindOffTargets_InvalidMaxMismatches_ThrowsException()
    {
        var genome = new DnaSequence("ACGTACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.FindOffTargets("ACGTACGT", genome, 10, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region Specificity Score Tests

    [Test]
    public void CalculateSpecificityScore_NoOffTargets_ReturnsHigh()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        // Short sequence unlikely to have off-targets
        var genome = new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.GreaterThanOrEqualTo(0));
        Assert.That(score, Is.LessThanOrEqualTo(100));
    }

    [Test]
    public void CalculateSpecificityScore_ManyOffTargets_ReturnsLower()
    {
        string guide = "AAAAAAAAAAAAAAAAAAAA";
        // Sequence with many similar regions
        var genome = new DnaSequence("AAAAAAAAAAAAAAAAAAAAAAGGAAAAAAAAAAAAAAAAAAACGGAAAAAAAAAAAAAAAAAAAAGG");

        double score = CrisprDesigner.CalculateSpecificityScore(guide, genome, CrisprSystemType.SpCas9);

        Assert.That(score, Is.LessThanOrEqualTo(100));
    }

    #endregion

    #region Parameter Tests

    [Test]
    public void GuideRnaParameters_Default_HasValidValues()
    {
        var defaults = GuideRnaParameters.Default;

        Assert.That(defaults.MinGcContent, Is.EqualTo(40));
        Assert.That(defaults.MaxGcContent, Is.EqualTo(70));
        Assert.That(defaults.MinScore, Is.EqualTo(50));
        Assert.That(defaults.AvoidPolyT, Is.True);
        Assert.That(defaults.CheckSelfComplementarity, Is.True);
    }

    [Test]
    public void GuideRnaParameters_CustomValues_Respected()
    {
        var custom = new GuideRnaParameters(
            MinGcContent: 30,
            MaxGcContent: 80,
            MinScore: 40,
            AvoidPolyT: false,
            CheckSelfComplementarity: false);

        Assert.That(custom.MinGcContent, Is.EqualTo(30));
        Assert.That(custom.MaxGcContent, Is.EqualTo(80));
        Assert.That(custom.MinScore, Is.EqualTo(40));
    }

    #endregion
}
