namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for restriction enzyme analysis: site finding and digestion.
///
/// Test Units: RESTR-FIND-001, RESTR-DIGEST-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class RestrictionProperties
{
    // Sequence containing known EcoRI (GAATTC) and BamHI (GGATCC) sites
    private const string TestSequence =
        "AAAGAATTCAAAGGATCCAAAGAATTCAAAGGATCCAAA";

    // -- RESTR-FIND-001 --

    /// <summary>
    /// Found restriction sites contain the recognized sequence matching the enzyme.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_RecognizedSequence_MatchesEnzyme()
    {
        var dna = new DnaSequence(TestSequence);
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI");
        Assert.That(enzyme, Is.Not.Null, "EcoRI not found in enzyme database");

        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();

        foreach (var site in sites)
            Assert.That(site.RecognizedSequence, Does.Contain(enzyme!.RecognitionSequence).IgnoreCase,
                $"Site at {site.Position} recognized '{site.RecognizedSequence}' doesn't match '{enzyme.RecognitionSequence}'");
    }

    /// <summary>
    /// All site positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_Positions_WithinBounds()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();

        foreach (var site in sites)
        {
            Assert.That(site.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(site.Position, Is.LessThan(TestSequence.Length));
        }
    }

    /// <summary>
    /// FindSites for non-existing recognition site returns empty.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSites_NoSite_ReturnsEmpty()
    {
        // Sequence with no HindIII (AAGCTT) site
        var dna = new DnaSequence("ACGTACGTACGTACGT");
        var sites = RestrictionAnalyzer.FindSites(dna, "HindIII").ToList();
        Assert.That(sites, Is.Empty);
    }

    /// <summary>
    /// Enzymes dictionary is not empty and all entries have recognition sequences.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Enzymes_AllHaveRecognitionSequence()
    {
        var enzymes = RestrictionAnalyzer.Enzymes;
        Assert.That(enzymes.Count, Is.GreaterThan(0), "Enzyme database should not be empty");

        foreach (var kvp in enzymes)
        {
            Assert.That(kvp.Value.RecognitionSequence, Is.Not.Null.And.Not.Empty,
                $"Enzyme {kvp.Key} has no recognition sequence");
            Assert.That(kvp.Value.RecognitionLength, Is.GreaterThan(0),
                $"Enzyme {kvp.Key} has zero recognition length");
        }
    }

    // -- RESTR-DIGEST-001 --

    /// <summary>
    /// Digest fragment lengths sum to original sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_FragmentLengths_SumToSequenceLength()
    {
        var dna = new DnaSequence(TestSequence);
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        int totalLength = fragments.Sum(f => f.Length);
        Assert.That(totalLength, Is.EqualTo(TestSequence.Length),
            $"Fragment lengths sum {totalLength} ≠ sequence length {TestSequence.Length}");
    }

    /// <summary>
    /// Number of digest fragments is ≥ 2 when sites exist.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_FragmentCount_AtLeastTwo()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI").ToList();
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        Assert.That(sites.Count, Is.GreaterThan(0), "Should find at least one EcoRI site");
        Assert.That(fragments.Count, Is.GreaterThanOrEqualTo(2),
            $"Expected >= 2 fragments for {sites.Count} sites, got {fragments.Count}");
    }

    /// <summary>
    /// DigestSummary total matches fragment count.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DigestSummary_TotalFragments_MatchesDigest()
    {
        var dna = new DnaSequence(TestSequence);
        var summary = RestrictionAnalyzer.GetDigestSummary(dna, "EcoRI");
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();

        Assert.That(summary.TotalFragments, Is.EqualTo(fragments.Count));
    }

    /// <summary>
    /// Multi-enzyme digest produces more fragments than single enzyme.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_MultiEnzyme_MoreFragments()
    {
        var dna = new DnaSequence(TestSequence);
        var singleFragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
        var multiFragments = RestrictionAnalyzer.Digest(dna, "EcoRI", "BamHI").ToList();

        Assert.That(multiFragments.Count, Is.GreaterThanOrEqualTo(singleFragments.Count),
            "Multi-enzyme digest should produce >= fragments than single enzyme");
    }

    /// <summary>
    /// All fragment lengths are positive.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Digest_AllFragments_HavePositiveLength()
    {
        var dna = new DnaSequence(TestSequence);
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI", "BamHI").ToList();

        foreach (var f in fragments)
            Assert.That(f.Length, Is.GreaterThan(0), $"Fragment {f.FragmentNumber} has non-positive length {f.Length}");
    }
}
