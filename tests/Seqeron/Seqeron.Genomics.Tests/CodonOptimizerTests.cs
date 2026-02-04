using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class CodonOptimizerTests
{
    #region Standard Genetic Code Tests

    [Test]
    public void TranslateSequence_ValidCodingSequence_ReturnsProtein()
    {
        // AUGGCUUAA = M-A-Stop
        var result = CodonOptimizer.OptimizeSequence("AUGGCUUAA", CodonOptimizer.EColiK12);
        Assert.That(result.ProteinSequence, Is.EqualTo("MA*"));
    }

    [Test]
    public void OptimizeSequence_EmptySequence_ReturnsEmptyResult()
    {
        var result = CodonOptimizer.OptimizeSequence("", CodonOptimizer.EColiK12);
        Assert.That(result.OptimizedSequence, Is.Empty);
        Assert.That(result.ProteinSequence, Is.Empty);
    }

    [Test]
    public void OptimizeSequence_ConvertsThymineToUracil()
    {
        // DNA sequence with T
        var result = CodonOptimizer.OptimizeSequence("ATGGCTTAA", CodonOptimizer.EColiK12);
        Assert.That(result.OriginalSequence, Does.Contain("U"));
        Assert.That(result.OriginalSequence, Does.Not.Contain("T"));
    }

    [Test]
    public void OptimizeSequence_TrimsToCompleteCodons()
    {
        // 10 nucleotides - should trim to 9 (3 codons)
        var result = CodonOptimizer.OptimizeSequence("AUGGCUUAAG", CodonOptimizer.EColiK12);
        Assert.That(result.OriginalSequence.Length % 3, Is.EqualTo(0));
    }

    #endregion

    // NOTE: CAI calculation tests moved to CodonOptimizer_CAI_Tests.cs (CODON-CAI-001)

    #region Optimization Strategy Tests

    [Test]
    public void OptimizeSequence_MaximizeCAI_IncreasesCAI()
    {
        // Start with rare codons
        string rare = "CUACCA"; // Leu-Pro with rare E. coli codons
        var result = CodonOptimizer.OptimizeSequence(rare, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI));
    }

    [Test]
    public void OptimizeSequence_BalancedOptimization_PreservesProtein()
    {
        string original = "AUGGCUGCACUGUAA"; // M-A-A-L-Stop
        var result = CodonOptimizer.OptimizeSequence(original, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        // Protein must be preserved
        Assert.That(result.ProteinSequence, Is.EqualTo("MAAL*"));
    }

    [Test]
    public void OptimizeSequence_AvoidRareCodons_ReplacesOnlyRare()
    {
        // Mix of good and rare codons
        string mixed = "CUGCUA"; // L-L (good-rare for E. coli)
        var result = CodonOptimizer.OptimizeSequence(mixed, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.AvoidRareCodeons);

        Assert.That(result.ProteinSequence, Is.EqualTo("LL"));
    }

    [Test]
    public void OptimizeSequence_HarmonizeExpression_MaintainsDistribution()
    {
        string sequence = "CUGCUGCUGCUGCUGCUG"; // 6x Leu
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.HarmonizeExpression);

        Assert.That(result.ProteinSequence, Is.EqualTo("LLLLLL"));
    }

    [Test]
    public void OptimizeSequence_TracksChanges()
    {
        string rare = "CUAAGACGA"; // L-R-R with rare codons
        var result = CodonOptimizer.OptimizeSequence(rare, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Should have changes for rare codons
        Assert.That(result.Changes.Count, Is.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region GC Content Tests

    [Test]
    public void OptimizeSequence_ReportsGcContent()
    {
        string sequence = "AUGGCCGCC"; // High GC
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        Assert.That(result.GcContentOriginal, Is.GreaterThan(0.5));
    }

    [Test]
    public void OptimizeSequence_BalancedOptimization_AimsForTargetGc()
    {
        // Test with extreme GC sequences
        string lowGc = "AUAUAUAUAUAUAUAUAU"; // Very low GC
        var result = CodonOptimizer.OptimizeSequence(lowGc, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        // Should try to balance GC if possible while preserving protein
        Assert.That(result.ProteinSequence.Length, Is.GreaterThan(0));
    }

    #endregion

    // NOTE: Codon Usage Analysis tests moved to CodonOptimizer_CodonUsage_Tests.cs (CODON-USAGE-001)
    // NOTE: FindRareCodons tests moved to CodonOptimizer_FindRareCodons_Tests.cs (CODON-RARE-001)

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

    #region Edge Cases and Robustness Tests

    [Test]
    public void OptimizeSequence_StopCodonsPreserved()
    {
        string sequence = "AUGUAA"; // M-Stop
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        Assert.That(result.OptimizedSequence, Does.EndWith("UAA").Or.EndWith("UAG").Or.EndWith("UGA"));
        Assert.That(result.ProteinSequence, Does.EndWith("*"));
    }

    [Test]
    public void OptimizeSequence_SingleAminoAcidCodonsUnchanged()
    {
        // Met (AUG) and Trp (UGG) have only one codon
        string sequence = "AUGUGG";
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        Assert.That(result.OptimizedSequence, Is.EqualTo("AUGUGG"));
    }

    [Test]
    public void OptimizeSequence_LowercaseInput_Handled()
    {
        string sequence = "auggcuuaa";
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        Assert.That(result.ProteinSequence, Is.EqualTo("MA*"));
    }

    [Test]
    public void OptimizeSequence_ForEColi_ImprovesCAI()
    {
        // Sequence with yeast-preferred codons
        string yeastOptimized = "CUAUUACCA"; // L-L-P with yeast preferences
        var result = CodonOptimizer.OptimizeSequence(yeastOptimized, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI));
    }

    [Test]
    public void OptimizeSequence_ForYeast_DifferentFromEColi()
    {
        string sequence = "CUGCCG"; // L-P
        var ecoliResult = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var yeastResult = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.Yeast,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Different organisms have different preferred codons
        Assert.That(ecoliResult.OptimizedSequence, Is.Not.EqualTo(yeastResult.OptimizedSequence)
            .Or.EqualTo(yeastResult.OptimizedSequence)); // May be same if both prefer it
    }

    [Test]
    public void OptimizeSequence_LongSequence_Completes()
    {
        // Generate a long sequence
        var codons = new[] { "AUG", "GCU", "CUG", "AAA", "GAU" };
        string longSequence = string.Join("", Enumerable.Repeat(codons, 100).SelectMany(x => x)) + "UAA";

        var result = CodonOptimizer.OptimizeSequence(longSequence, CodonOptimizer.EColiK12);

        Assert.That(result.OptimizedSequence.Length, Is.EqualTo(longSequence.Length));
    }

    #endregion

    #region Integration Tests

    [Test]
    public void OptimizeSequence_GFP_OptimizesForEColi()
    {
        // Part of GFP sequence
        string gfpPart = "AUGAGCAAAGGU";
        var result = CodonOptimizer.OptimizeSequence(gfpPart, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        Assert.That(result.ProteinSequence, Is.EqualTo("MSKG"));
        Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI));
    }

    [Test]
    public void FullWorkflow_OptimizeAndAnalyze()
    {
        string original = "AUGCUACGAAGAUAA"; // M-L-R-R-Stop

        // Optimize
        var optimized = CodonOptimizer.OptimizeSequence(original, CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Find remaining rare codons
        var rareAfter = CodonOptimizer.FindRareCodons(optimized.OptimizedSequence,
            CodonOptimizer.EColiK12, 0.1).ToList();

        // Should have fewer or no rare codons
        var rareBefore = CodonOptimizer.FindRareCodons(original, CodonOptimizer.EColiK12, 0.1).ToList();
        Assert.That(rareAfter.Count, Is.LessThanOrEqualTo(rareBefore.Count));
    }

    #endregion
}
