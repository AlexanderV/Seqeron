namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for miRNA analysis: target prediction, pre-miRNA hairpins.
///
/// Test Units: MIRNA-TARGET-001, MIRNA-PRECURSOR-001
/// MIRNA-SEED-001 property tests removed — consolidated into canonical MiRnaAnalyzer_SeedAnalysis_Tests.cs
/// (3 duplicates of M-003/M-007/M-009; 1 weak: IsSubstringOfMiRna can't distinguish extraction position)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class MiRnaProperties
{
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
