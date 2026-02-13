namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for splice site prediction: donor/acceptor sites, gene structure.
///
/// Test Units: SPLICE-DONOR-001, SPLICE-ACCEPTOR-001, SPLICE-PREDICT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class SplicingProperties
{
    // Sequence with canonical GT donor and AG acceptor sites
    // Exon1-GT...intron...AG-Exon2
    private const string TestSequence =
        "ATGATGAAAGCCGCCATGGCG" +       // exon1 (21bp)
        "GTAAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCAG" + // intron (~85bp with GT...AG)
        "GCGATGAAAGCCGCCATGGCG";        // exon2 (21bp)

    // -- SPLICE-DONOR-001 --

    /// <summary>
    /// Donor sites have GT/GC motif.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DonorSites_HaveCanonicalMotif()
    {
        var donors = SpliceSitePredictor.FindDonorSites(TestSequence, minScore: 0.1).ToList();

        foreach (var d in donors)
            Assert.That(d.Motif.ToUpperInvariant(), Does.Contain("GT").Or.Contain("GU").Or.Contain("GC"),
                $"Donor at {d.Position} has motif '{d.Motif}'");
    }

    /// <summary>
    /// Donor site scores are in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void DonorSites_Score_InRange()
    {
        var donors = SpliceSitePredictor.FindDonorSites(TestSequence, minScore: 0.0).ToList();

        foreach (var d in donors)
            Assert.That(d.Score, Is.InRange(0.0, 1.0),
                $"Donor score {d.Score} at {d.Position} out of range");
    }

    /// <summary>
    /// Donor site positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DonorSites_Positions_WithinBounds()
    {
        var donors = SpliceSitePredictor.FindDonorSites(TestSequence, minScore: 0.0).ToList();

        foreach (var d in donors)
        {
            Assert.That(d.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(d.Position, Is.LessThan(TestSequence.Length));
        }
    }

    // -- SPLICE-ACCEPTOR-001 --

    /// <summary>
    /// Acceptor sites have AG motif.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AcceptorSites_HaveCanonicalMotif()
    {
        var acceptors = SpliceSitePredictor.FindAcceptorSites(TestSequence, minScore: 0.1).ToList();

        foreach (var a in acceptors)
            Assert.That(a.Motif, Does.EndWith("G").Or.Contain("AG"),
                $"Acceptor at {a.Position} has motif '{a.Motif}'");
    }

    /// <summary>
    /// Acceptor site scores are in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void AcceptorSites_Score_InRange()
    {
        var acceptors = SpliceSitePredictor.FindAcceptorSites(TestSequence, minScore: 0.0).ToList();

        foreach (var a in acceptors)
            Assert.That(a.Score, Is.InRange(0.0, 1.0),
                $"Acceptor score {a.Score} at {a.Position} out of range");
    }

    /// <summary>
    /// Acceptor site positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AcceptorSites_Positions_WithinBounds()
    {
        var acceptors = SpliceSitePredictor.FindAcceptorSites(TestSequence, minScore: 0.0).ToList();

        foreach (var a in acceptors)
        {
            Assert.That(a.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(a.Position, Is.LessThan(TestSequence.Length));
        }
    }

    // -- SPLICE-PREDICT-001 --

    /// <summary>
    /// PredictGeneStructure exon count ≥ 1 (at minimum a single-exon gene).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_HasAtLeastOneExon()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(TestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.Exons.Count, Is.GreaterThanOrEqualTo(1),
            "Gene structure should have at least one exon");
    }

    /// <summary>
    /// Exon and intron counts are consistent: introns = exons - 1 (if any introns).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_IntronCount_IsExonsMinusOne()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(TestSequence, minExonLength: 10, minIntronLength: 30);

        if (structure.Introns.Count > 0)
            Assert.That(structure.Introns.Count, Is.EqualTo(structure.Exons.Count - 1),
                $"Expected {structure.Exons.Count - 1} introns for {structure.Exons.Count} exons");
    }

    /// <summary>
    /// SplicedSequence is shorter than or equal to original (introns removed).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_SplicedSequence_ShorterOrEqual()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(TestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.SplicedSequence.Length, Is.LessThanOrEqualTo(TestSequence.Length),
            "Spliced sequence should be ≤ original length");
    }

    /// <summary>
    /// Gene structure overall score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_OverallScore_InRange()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(TestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.OverallScore, Is.InRange(0.0, 1.0),
            $"Overall score {structure.OverallScore} out of range");
    }
}
