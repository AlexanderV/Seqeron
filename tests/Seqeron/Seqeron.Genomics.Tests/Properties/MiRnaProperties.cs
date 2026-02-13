namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for miRNA analysis: seed extraction, target prediction, pre-miRNA hairpins.
///
/// Test Units: MIRNA-SEED-001, MIRNA-TARGET-001, MIRNA-PRECURSOR-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class MiRnaProperties
{
    // -- MIRNA-SEED-001 --

    /// <summary>
    /// Seed sequence length is 7 (positions 2-8 of mature miRNA).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GetSeedSequence_HasLength7()
    {
        string miRna = "UAGCUUAUCAGACUGAUGUUGA";
        string seed = MiRnaAnalyzer.GetSeedSequence(miRna);

        Assert.That(seed.Length, Is.EqualTo(7),
            $"Seed '{seed}' should have length 7");
    }

    /// <summary>
    /// Seed is a substring of the original miRNA.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GetSeedSequence_IsSubstringOfMiRna()
    {
        string miRna = "UAGCUUAUCAGACUGAUGUUGA";
        string seed = MiRnaAnalyzer.GetSeedSequence(miRna);

        Assert.That(miRna, Does.Contain(seed),
            $"Seed '{seed}' should be substring of '{miRna}'");
    }

    /// <summary>
    /// CreateMiRna preserves name and sequence.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CreateMiRna_PreservesNameAndSequence()
    {
        string name = "hsa-miR-21-5p";
        string seq = "UAGCUUAUCAGACUGAUGUUGA";
        var mirna = MiRnaAnalyzer.CreateMiRna(name, seq);

        Assert.That(mirna.Name, Is.EqualTo(name));
        Assert.That(mirna.Sequence, Is.EqualTo(seq));
        Assert.That(mirna.SeedSequence, Is.Not.Null.And.Not.Empty);
    }

    /// <summary>
    /// Same family miRNAs share seed region.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CompareSeedRegions_IdenticalSeeds_AllMatches()
    {
        var mirna1 = MiRnaAnalyzer.CreateMiRna("miR-A", "UAGCUUAUCAGACUGAUGUUGA");
        var mirna2 = MiRnaAnalyzer.CreateMiRna("miR-B", "UAGCUUAUCAGACUGAUGUUGA");
        var comparison = MiRnaAnalyzer.CompareSeedRegions(mirna1, mirna2);

        Assert.That(comparison.Mismatches, Is.EqualTo(0));
        Assert.That(comparison.IsSameFamily, Is.True);
    }

    // -- MIRNA-TARGET-001 --

    /// <summary>
    /// Target site positions are within mRNA bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindTargetSites_Positions_WithinBounds()
    {
        string mRna = "AUGCCAUUUUAGCUUAUCAGACAACUAUGAAUCCAAUUAGCUUAUCAGACAACUAUUU";
        var mirna = MiRnaAnalyzer.CreateMiRna("miR-test", "UAGCUUAUCAGACUGAUGUUGA");
        var sites = MiRnaAnalyzer.FindTargetSites(mRna, mirna).ToList();

        foreach (var site in sites)
        {
            Assert.That(site.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(site.End, Is.LessThanOrEqualTo(mRna.Length));
        }
    }

    /// <summary>
    /// Target site score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindTargetSites_Score_InRange()
    {
        string mRna = "AUGCCAUUUUAGCUUAUCAGACAACUAUGAAUCCAAUUAGCUUAUCAGACAACUAUUU";
        var mirna = MiRnaAnalyzer.CreateMiRna("miR-test", "UAGCUUAUCAGACUGAUGUUGA");
        var sites = MiRnaAnalyzer.FindTargetSites(mRna, mirna).ToList();

        foreach (var site in sites)
            Assert.That(site.Score, Is.InRange(0.0, 1.0),
                $"Score {site.Score} out of range at {site.Start}");
    }

    // -- MIRNA-PRECURSOR-001 --

    /// <summary>
    /// Pre-miRNA hairpin has structure with balanced parentheses.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_Structure_BalancedParentheses()
    {
        // A sequence long enough to contain potential hairpins
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
        {
            int opens = hp.Structure.Count(c => c == '(');
            int closes = hp.Structure.Count(c => c == ')');
            Assert.That(opens, Is.EqualTo(closes),
                $"Unbalanced structure: {opens} opens vs {closes} closes");
        }
    }

    /// <summary>
    /// Pre-miRNA start/end are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_Bounds_WithinSequence()
    {
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
        {
            Assert.That(hp.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(hp.End, Is.LessThanOrEqualTo(sequence.Length));
            Assert.That(hp.End, Is.GreaterThan(hp.Start));
        }
    }

    /// <summary>
    /// Pre-miRNA structure length equals sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPreMiRnaHairpins_StructureLength_EqualsSequenceLength()
    {
        string sequence = string.Concat(Enumerable.Repeat("GCGCUUUUGCGC", 20));
        var hairpins = MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength: 30).ToList();

        foreach (var hp in hairpins)
            Assert.That(hp.Structure.Length, Is.EqualTo(hp.Sequence.Length),
                $"Structure length {hp.Structure.Length} ≠ sequence length {hp.Sequence.Length}");
    }
}
