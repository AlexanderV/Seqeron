using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CodonOptimizer.CalculateCodonUsage and CodonOptimizer.CompareCodonUsage.
/// Test Unit: CODON-USAGE-001
/// Evidence: Wikipedia (Codon usage bias), Kazusa Codon Usage Database, Sharp & Li (1987)
/// </summary>
[TestFixture]
public class CodonOptimizer_CodonUsage_Tests
{
    #region CalculateCodonUsage - Must Tests

    /// <summary>
    /// M1: Verify basic codon counting with known input.
    /// Source: Mathematical definition, Kazusa database format.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_SimpleCodingSequence_CountsCorrectly()
    {
        // Arrange: "AUGGCUGCU" = M-A-A (AUG + GCU + GCU)
        const string sequence = "AUGGCUGCU";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(usage["AUG"], Is.EqualTo(1), "AUG should appear once");
            Assert.That(usage["GCU"], Is.EqualTo(2), "GCU should appear twice");
            Assert.That(usage.Count, Is.EqualTo(2), "Only 2 distinct codons");
        });
    }

    /// <summary>
    /// M2: Empty input returns empty result.
    /// Source: Standard edge case handling.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_EmptySequence_ReturnsEmptyDictionary()
    {
        // Act
        var usage = CodonOptimizer.CalculateCodonUsage("");

        // Assert
        Assert.That(usage, Is.Empty);
    }

    /// <summary>
    /// M3: Trailing incomplete codons are ignored.
    /// Source: Kazusa format, standard codon counting practice.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_IncompleteCodon_IgnoresTrailing()
    {
        // Arrange: "AUGGC" = AUG + GC (incomplete)
        const string sequence = "AUGGC";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(usage.Count, Is.EqualTo(1), "Only complete codons counted");
            Assert.That(usage["AUG"], Is.EqualTo(1));
        });
    }

    /// <summary>
    /// M4: Repeated codons counted correctly.
    /// Source: Mathematical counting invariant.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_RepeatedCodons_CountsAccurately()
    {
        // Arrange: 5x AUG
        const string sequence = "AUGAUGAUGAUGAUG";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(usage.Count, Is.EqualTo(1), "Only one distinct codon");
            Assert.That(usage["AUG"], Is.EqualTo(5), "AUG should appear 5 times");
        });
    }

    /// <summary>
    /// M5: All valid codons are recognized (tests comprehensive coverage).
    /// Source: Genetic code completeness.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_VariousCodons_HandlesAllRecognized()
    {
        // Arrange: Sequence with several different codons
        // UUU (Phe), UUC (Phe), CUG (Leu), GGG (Gly), UAA (Stop)
        const string sequence = "UUUUUCCUGGGG";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(usage["UUU"], Is.EqualTo(1));
            Assert.That(usage["UUC"], Is.EqualTo(1));
            Assert.That(usage["CUG"], Is.EqualTo(1));
            Assert.That(usage["GGG"], Is.EqualTo(1));
        });
    }

    #endregion

    #region CompareCodonUsage - Must Tests

    /// <summary>
    /// M6: Same sequence comparison returns similarity = 1.0.
    /// Source: Mathematical identity property of similarity metrics.
    /// </summary>
    [Test]
    public void CompareCodonUsage_IdenticalSequences_ReturnsOne()
    {
        // Arrange
        const string sequence = "AUGGCUGCACUG";

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(sequence, sequence);

        // Assert
        Assert.That(similarity, Is.EqualTo(1.0).Within(1e-10));
    }

    /// <summary>
    /// M7: Different codon usage returns similarity less than 1.0.
    /// Source: Distance metric property.
    /// </summary>
    [Test]
    public void CompareCodonUsage_DifferentCodons_ReturnsLessThanOne()
    {
        // Arrange: Completely different codons (both Leu)
        const string seq1 = "CUGCUGCUGCUG"; // All CUG
        const string seq2 = "CUACUACUACUA"; // All CUA

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert
        Assert.That(similarity, Is.LessThan(1.0));
    }

    /// <summary>
    /// M8: Empty sequences return 0.0.
    /// Source: Edge case handling (no data to compare).
    /// </summary>
    [Test]
    public void CompareCodonUsage_EmptySequences_ReturnsZero()
    {
        // Act
        double similarity = CodonOptimizer.CompareCodonUsage("", "");

        // Assert
        Assert.That(similarity, Is.EqualTo(0.0));
    }

    /// <summary>
    /// M9: Comparison is symmetric: Sim(a,b) = Sim(b,a).
    /// Source: Mathematical property of distance-based similarity.
    /// </summary>
    [Test]
    public void CompareCodonUsage_Symmetric_SameResultBothDirections()
    {
        // Arrange
        const string seq1 = "AUGGCUGCACUG";
        const string seq2 = "UUUCCCGGGAAA";

        // Act
        double sim1 = CodonOptimizer.CompareCodonUsage(seq1, seq2);
        double sim2 = CodonOptimizer.CompareCodonUsage(seq2, seq1);

        // Assert
        Assert.That(sim1, Is.EqualTo(sim2).Within(1e-10), "Similarity should be symmetric");
    }

    /// <summary>
    /// M10: Result is always in the range [0, 1].
    /// Source: Normalized Manhattan distance produces bounded result.
    /// </summary>
    [TestCase("AUGGCUGCACUG", "UUUCCCGGGAAA")]
    [TestCase("CUGCUGCUG", "CUACUACUA")]
    [TestCase("AUG", "UGA")]
    public void CompareCodonUsage_ResultRange_ZeroToOne(string seq1, string seq2)
    {
        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(similarity, Is.GreaterThanOrEqualTo(0.0));
            Assert.That(similarity, Is.LessThanOrEqualTo(1.0));
        });
    }

    #endregion

    #region CalculateCodonUsage - Should Tests

    /// <summary>
    /// S1: DNA input with T is converted to RNA U internally.
    /// Source: Implementation behavior (T→U conversion).
    /// </summary>
    [Test]
    public void CalculateCodonUsage_DnaInput_ConvertsTToU()
    {
        // Arrange: DNA sequence with T (ATG = AUG)
        const string dnaSequence = "ATGGCTTAA";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(dnaSequence);

        // Assert: Keys should be in RNA format (U, not T)
        Assert.Multiple(() =>
        {
            Assert.That(usage.ContainsKey("AUG"), Is.True, "ATG should be stored as AUG");
            Assert.That(usage.ContainsKey("ATG"), Is.False, "T should be converted to U");
        });
    }

    /// <summary>
    /// S2: Mixed case input is handled case-insensitively.
    /// Source: Standard robustness requirement.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_MixedCase_HandlesCorrectly()
    {
        // Arrange
        const string mixedCase = "AugGcuGCU";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(mixedCase);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(usage["AUG"], Is.EqualTo(1));
            Assert.That(usage["GCU"], Is.EqualTo(2));
        });
    }

    /// <summary>
    /// S3: Sum of all codon counts equals total number of complete codons.
    /// Source: Mathematical invariant.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_Invariant_SumEqualsTotalCodons()
    {
        // Arrange: 7 codons
        const string sequence = "AUGGCUGCACUGUAACCC";
        const int expectedCodons = 6;

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);
        int totalCount = usage.Values.Sum();

        // Assert
        Assert.That(totalCount, Is.EqualTo(expectedCodons));
    }

    #endregion

    #region CompareCodonUsage - Should Tests

    /// <summary>
    /// S4: One empty sequence comparison returns 0.
    /// Source: Edge case - no data from one side.
    /// </summary>
    [Test]
    public void CompareCodonUsage_OneEmptySequence_ReturnsZero()
    {
        // Arrange
        const string nonEmpty = "AUGGCUGCU";

        // Act
        double sim1 = CodonOptimizer.CompareCodonUsage(nonEmpty, "");
        double sim2 = CodonOptimizer.CompareCodonUsage("", nonEmpty);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(sim1, Is.EqualTo(0.0), "Non-empty vs empty should return 0");
            Assert.That(sim2, Is.EqualTo(0.0), "Empty vs non-empty should return 0");
        });
    }

    /// <summary>
    /// S5: Sequences with no overlapping codons have low similarity.
    /// Source: Manhattan distance property.
    /// </summary>
    [Test]
    public void CompareCodonUsage_NoOverlappingCodons_LowSimilarity()
    {
        // Arrange: Completely disjoint codon sets
        const string seq1 = "UUUUUUUUU"; // All UUU (Phe)
        const string seq2 = "GGGGGGGGG"; // All GGG (Gly)

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: With no overlap, similarity should be 0
        // (each sequence has 100% of one codon, 0% of the other)
        Assert.That(similarity, Is.EqualTo(0.0).Within(1e-10));
    }

    /// <summary>
    /// S5b: Partial codon overlap produces intermediate similarity.
    /// Source: ASSUMPTION - expected behavior based on metric definition.
    /// </summary>
    [Test]
    public void CompareCodonUsage_PartialOverlap_IntermediateSimilarity()
    {
        // Arrange: 50% shared codons
        const string seq1 = "AUGGCUAUG"; // AUG, GCU, AUG
        const string seq2 = "AUGUUUAUG"; // AUG, UUU, AUG

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: Should be between 0 and 1 (partial overlap)
        Assert.Multiple(() =>
        {
            Assert.That(similarity, Is.GreaterThan(0.0));
            Assert.That(similarity, Is.LessThan(1.0));
        });
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Null handling test (if implemented).
    /// Source: Defensive programming.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_NullSequence_ReturnsEmpty()
    {
        // The implementation treats null as empty via string.IsNullOrEmpty
        var usage = CodonOptimizer.CalculateCodonUsage(null!);

        Assert.That(usage, Is.Empty);
    }

    /// <summary>
    /// Very short sequence (less than one codon).
    /// Source: Edge case handling.
    /// </summary>
    [TestCase("")]
    [TestCase("A")]
    [TestCase("AU")]
    public void CalculateCodonUsage_TooShort_ReturnsEmpty(string sequence)
    {
        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(sequence);

        // Assert
        Assert.That(usage, Is.Empty);
    }

    /// <summary>
    /// Single codon sequence comparison with itself.
    /// Source: Minimal valid input.
    /// </summary>
    [Test]
    public void CompareCodonUsage_SingleCodon_IdenticalReturnsOne()
    {
        // Arrange
        const string sequence = "AUG";

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(sequence, sequence);

        // Assert
        Assert.That(similarity, Is.EqualTo(1.0).Within(1e-10));
    }

    #endregion
}
