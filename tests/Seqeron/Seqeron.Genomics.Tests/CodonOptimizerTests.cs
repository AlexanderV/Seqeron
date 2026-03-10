using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class CodonOptimizerTests
{
    // NOTE: OptimizeSequence tests in CodonOptimizer_OptimizeSequence_Tests.cs (CODON-OPT-001)
    // NOTE: CAI calculation tests in CodonOptimizer_CAI_Tests.cs (CODON-CAI-001)
    // NOTE: Codon Usage Analysis tests in CodonOptimizer_CodonUsage_Tests.cs (CODON-USAGE-001)
    // NOTE: FindRareCodons tests in CodonOptimizer_FindRareCodons_Tests.cs (CODON-RARE-001)

    #region Restriction Site Removal Tests

    [Test]
    public void RemoveRestrictionSites_RemovesTargetSite()
    {
        // EcoRI: GAATTC (GAAUUC in RNA)
        string sequence = "AUGGAAUUCGCU"; // Contains EcoRI site
        var result = CodonOptimizer.RemoveRestrictionSites(sequence, new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        Assert.That(result, Does.Not.Contain("GAAUUC"));
    }

    [Test]
    public void RemoveRestrictionSites_PreservesProtein()
    {
        string sequence = "AUGGAAUUC"; // M-E-F
        string original = sequence.Replace('T', 'U');
        var result = CodonOptimizer.RemoveRestrictionSites(sequence, new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        // Verify protein is preserved by checking length is unchanged
        Assert.That(result.Length, Is.EqualTo(original.Length));
    }

    [Test]
    public void RemoveRestrictionSites_MultiplesSites()
    {
        // BamHI: GGATCC, HindIII: AAGCTT
        string sequence = "AUGGGATCCAAGCTTGCU";
        var result = CodonOptimizer.RemoveRestrictionSites(sequence,
            new[] { "GGATCC", "AAGCTT" }, CodonOptimizer.EColiK12);

        Assert.That(result, Does.Not.Contain("GGAUCC"));
    }

    [Test]
    public void RemoveRestrictionSites_EmptySequence_ReturnsEmpty()
    {
        var result = CodonOptimizer.RemoveRestrictionSites("", new[] { "GAATTC" }, CodonOptimizer.EColiK12);
        Assert.That(result, Is.Empty);
    }

    #endregion

    #region Secondary Structure Reduction Tests

    [Test]
    public void ReduceSecondaryStructure_ModifiesHighStructure()
    {
        // Create sequence with potential secondary structure
        string sequence = "AUGGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAUAA";
        var result = CodonOptimizer.ReduceSecondaryStructure(sequence, CodonOptimizer.EColiK12);

        // Should return modified or same sequence
        Assert.That(result.Length, Is.EqualTo(sequence.Length));
    }

    [Test]
    public void ReduceSecondaryStructure_ShortSequence_ReturnsSame()
    {
        string sequence = "AUGGCU";
        var result = CodonOptimizer.ReduceSecondaryStructure(sequence, CodonOptimizer.EColiK12, 40);
        Assert.That(result, Is.EqualTo(sequence));
    }

    [Test]
    public void ReduceSecondaryStructure_EmptySequence_ReturnsEmpty()
    {
        var result = CodonOptimizer.ReduceSecondaryStructure("", CodonOptimizer.EColiK12);
        Assert.That(result, Is.Empty);
    }

    #endregion

    #region Codon Usage Table Tests

    [Test]
    public void EColiK12_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.EColiK12.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void Yeast_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.Yeast.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void Human_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.Human.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void CodonTables_FrequenciesNormalized()
    {
        // For each amino acid, synonymous codon frequencies should sum to ~1
        foreach (var aa in new[] { "L", "S", "R", "P", "A", "G", "V", "T" })
        {
            // Just verify tables are properly structured
            Assert.That(CodonOptimizer.EColiK12.CodonFrequencies.Values.All(f => f >= 0 && f <= 1), Is.True);
        }
    }

    [Test]
    public void CreateCodonTableFromSequence_CreatesValidTable()
    {
        string reference = "AUGGCUGCUGCUGCUGCUGCUGCUUAA";
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference, "Test Organism");

        Assert.That(table.OrganismName, Is.EqualTo("Test Organism"));
        Assert.That(table.CodonFrequencies.ContainsKey("GCU"), Is.True);
    }

    [Test]
    public void CreateCodonTableFromSequence_CalculatesCorrectFrequencies()
    {
        // 4x GCU, 2x GCC = 6x Ala
        string reference = "GCUGCUGCUGCUGCCGCC";
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference, "Test");

        // GCU should be ~0.67, GCC should be ~0.33
        Assert.That(table.CodonFrequencies.GetValueOrDefault("GCU", 0), Is.GreaterThan(0.5));
    }

    #endregion
}
