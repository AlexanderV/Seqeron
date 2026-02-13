namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for CRISPR guide design: PAM finding, guide RNA design, off-target scoring.
///
/// Test Units: CRISPR-PAM-001, CRISPR-GUIDE-001, CRISPR-OFF-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class CrisprProperties
{
    // A synthetic sequence with known NGG PAMs
    private const string TestSequence = "ACGTACGTACGTACGTACGTGGG" +
                                        "TTTTTTTTTTTTTTTTTTTTTGG" +
                                        "ACGTACGTACGTACGTACGTACGT";

    // -- CRISPR-PAM-001 --

    /// <summary>
    /// All found PAM sites contain the expected PAM motif.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PamSites_ContainExpectedMotif()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9).ToList();

        foreach (var site in sites)
            Assert.That(site.PamSequence, Does.Match("[ACGT]GG"),
                $"PAM '{site.PamSequence}' at {site.Position} must match NGG");
    }

    /// <summary>
    /// PAM positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PamSites_Positions_WithinBounds()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = CrisprDesigner.FindPamSites(dna).ToList();

        foreach (var site in sites)
        {
            Assert.That(site.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(site.Position, Is.LessThan(TestSequence.Length));
        }
    }

    /// <summary>
    /// Target sequence length matches the system's guide length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PamSites_TargetLength_MatchesSystem()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9);
        var dna = new DnaSequence(TestSequence);
        var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9).ToList();

        foreach (var site in sites.Where(s => s.TargetSequence.Length > 0))
            Assert.That(site.TargetSequence.Length, Is.EqualTo(system.GuideLength),
                $"Target '{site.TargetSequence}' length should be {system.GuideLength}");
    }

    // -- CRISPR-GUIDE-001 --

    /// <summary>
    /// Guide RNA GC content is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_GcContent_InRange()
    {
        var longSeq = new DnaSequence(string.Concat(Enumerable.Repeat("ACGTACGTACGTACGTACGTACGT", 5)));
        var guides = CrisprDesigner.DesignGuideRnas(longSeq, 0, longSeq.Length - 1).ToList();

        foreach (var g in guides)
            Assert.That(g.GcContent, Is.InRange(0.0, 1.0),
                $"GC content {g.GcContent} out of range for guide at {g.Position}");
    }

    /// <summary>
    /// Guide RNA score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_Score_InRange()
    {
        var longSeq = new DnaSequence(string.Concat(Enumerable.Repeat("ACGTACGTACGTACGTACGTACGT", 5)));
        var guides = CrisprDesigner.DesignGuideRnas(longSeq, 0, longSeq.Length - 1).ToList();

        foreach (var g in guides)
            Assert.That(g.Score, Is.InRange(0.0, 1.0),
                $"Score {g.Score} out of range for guide at {g.Position}");
    }

    // -- CRISPR-OFF-001 --

    /// <summary>
    /// Off-target mismatches is ≤ maxMismatches.
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_Mismatches_WithinLimit()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(TestSequence + guide + "GGGTTTTTTTTTT");
        int maxMismatches = 2;
        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, maxMismatches).ToList();

        foreach (var ot in offTargets)
            Assert.That(ot.Mismatches, Is.LessThanOrEqualTo(maxMismatches),
                $"Off-target at {ot.Position} has {ot.Mismatches} mismatches");
    }

    /// <summary>
    /// Off-target score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_Score_InRange()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(TestSequence + guide + "GGGTTTTTTTTTT");
        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3).ToList();

        foreach (var ot in offTargets)
            Assert.That(ot.OffTargetScore, Is.InRange(0.0, 1.0),
                $"Off-target score {ot.OffTargetScore} out of range");
    }
}
