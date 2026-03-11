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
    /// M5: All 64 standard RNA codons are recognized and counted.
    /// Source: Genetic code completeness — 64 codons (61 coding + 3 stop).
    /// </summary>
    [Test]
    public void CalculateCodonUsage_AllSixtyFourCodons_HandlesAll()
    {
        // Arrange: all 64 standard codons in codon-table order, each once
        const string allCodons =
            "UUUUUCUUAUUG" + // UUU UUC UUA UUG
            "CUUCUCCUACUG" + // CUU CUC CUA CUG
            "AUUAUCAUAAUG" + // AUU AUC AUA AUG
            "GUUGUCGUAGUG" + // GUU GUC GUA GUG
            "UCUUCCUCAUCG" + // UCU UCC UCA UCG
            "CCUCCCCCACCG" + // CCU CCC CCA CCG
            "ACUACCACAACG" + // ACU ACC ACA ACG
            "GCUGCCGCAGCG" + // GCU GCC GCA GCG
            "UAUUACUAAUAG" + // UAU UAC UAA UAG
            "CAUCACCAACAG" + // CAU CAC CAA CAG
            "AAUAACAAAAAG" + // AAU AAC AAA AAG
            "GAUGACGAAGAG" + // GAU GAC GAA GAG
            "UGUUGCUGAUGG" + // UGU UGC UGA UGG
            "CGUCGCCGACGG" + // CGU CGC CGA CGG
            "AGUAGCAGAAGG" + // AGU AGC AGA AGG
            "GGUGGCGGAGGG";  // GGU GGC GGA GGG

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(allCodons);

        // Assert: all 64 codons present, each exactly once
        Assert.Multiple(() =>
        {
            Assert.That(usage.Count, Is.EqualTo(64), "All 64 codons should be distinct keys");
            Assert.That(usage.Values.Sum(), Is.EqualTo(64), "Total count = 64");
            Assert.That(usage.Values, Has.All.EqualTo(1), "Each codon appears exactly once");
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
    /// Source: TVD formula. Exact calculation:
    ///   seq1: f(CUG)=3/4, f(CUA)=1/4
    ///   seq2: f(CUA)=3/4, f(CUG)=1/4
    ///   Σ|f₁-f₂| = |3/4-1/4| + |1/4-3/4| = 1/2 + 1/2 = 1
    ///   Similarity = 1 - 1/2 = 0.5
    /// </summary>
    [Test]
    public void CompareCodonUsage_DifferentCodons_ReturnsLessThanOne()
    {
        // Arrange: Partially overlapping Leu codons
        const string seq1 = "CUGCUGCUGCUA"; // 3×CUG + 1×CUA
        const string seq2 = "CUACUACUACUG"; // 3×CUA + 1×CUG

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: exact value from TVD formula
        Assert.That(similarity, Is.EqualTo(0.5).Within(1e-10));
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
    /// Source: TVD symmetry (|x-y| = |y-x|). Exact derivation:
    ///   seq1: f(AUG)=1/2, f(CCC)=1/4, f(UUU)=1/4
    ///   seq2: f(AUG)=1/4, f(UUU)=1/2, f(CCC)=1/4
    ///   Σ|f₁-f₂| = 1/4 + 0 + 1/4 = 1/2 → Similarity = 3/4
    /// </summary>
    [Test]
    public void CompareCodonUsage_Symmetric_SameResultBothDirections()
    {
        // Arrange: partially overlapping 3-codon distributions
        const string seq1 = "AUGAUGCCCUUU"; // AUG×2, CCC×1, UUU×1
        const string seq2 = "AUGUUUUUUCCC"; // AUG×1, UUU×2, CCC×1

        // Act
        double sim1 = CodonOptimizer.CompareCodonUsage(seq1, seq2);
        double sim2 = CodonOptimizer.CompareCodonUsage(seq2, seq1);

        // Assert: symmetric AND exact value
        Assert.Multiple(() =>
        {
            Assert.That(sim1, Is.EqualTo(sim2).Within(1e-10), "Similarity must be symmetric");
            Assert.That(sim1, Is.EqualTo(0.75).Within(1e-10), "TVD similarity = 3/4");
        });
    }

    /// <summary>
    /// M10: Results span [0, 1] with exact TVD values.
    /// Source: TVD formula derivations:
    ///   Case 1: disjoint (AUG vs CCC) → Σ=2 → sim=0.0
    ///   Case 2: f₁(AUG)=1, f₂(AUG)=1/4,f₂(CCC)=3/4 → Σ=3/2 → sim=1/4
    ///   Case 3: f₁(AUG)=3/4,f₁(CCC)=1/4, f₂(AUG)=1/2,f₂(CCC)=1/2 → Σ=1/2 → sim=3/4
    /// </summary>
    [TestCase("AUG", "CCC", 0.0)]
    [TestCase("AUGAUGAUGAUG", "AUGCCCCCCCCC", 0.25)]
    [TestCase("AUGAUGAUGCCC", "AUGAUGCCCCCC", 0.75)]
    public void CompareCodonUsage_ResultRange_ZeroToOne(string seq1, string seq2, double expected)
    {
        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: exact TVD value
        Assert.That(similarity, Is.EqualTo(expected).Within(1e-10));
    }

    #endregion

    #region CalculateCodonUsage - Should Tests

    /// <summary>
    /// S1: DNA input with T is converted to RNA U internally.
    /// Source: Biological equivalence of DNA T and RNA U.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_DnaInput_ConvertsTToU()
    {
        // Arrange: DNA "ATGGCTTAA" → RNA "AUGGCUUAA" → codons AUG, GCU, UAA
        const string dnaSequence = "ATGGCTTAA";

        // Act
        var usage = CodonOptimizer.CalculateCodonUsage(dnaSequence);

        // Assert: all codons in RNA (U) form, complete dictionary
        Assert.Multiple(() =>
        {
            Assert.That(usage.Count, Is.EqualTo(3), "3 distinct codons");
            Assert.That(usage["AUG"], Is.EqualTo(1), "ATG → AUG");
            Assert.That(usage["GCU"], Is.EqualTo(1), "GCT → GCU");
            Assert.That(usage["UAA"], Is.EqualTo(1), "TAA → UAA");
            Assert.That(usage.Keys.Any(k => k.Contains('T')), Is.False, "No T in any key");
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
        // Arrange: 18 chars = 6 complete codons (AUG, GCU, GCA, CUG, UAA, CCC)
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
    /// S5: Sequences with no overlapping codons yield zero similarity.
    /// Source: TVD formula. For disjoint distributions P and Q:
    ///   Σ|P(c)-Q(c)| = Σ P(c) + Σ Q(c) = 1 + 1 = 2
    ///   Similarity = 1 - 2/2 = 0
    /// </summary>
    [Test]
    public void CompareCodonUsage_NoOverlappingCodons_ZeroSimilarity()
    {
        // Arrange: Completely disjoint codon sets
        const string seq1 = "UUUUUUUUU"; // All UUU (Phe)
        const string seq2 = "GGGGGGGGG"; // All GGG (Gly)

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: disjoint distributions → similarity = 0 by TVD formula
        Assert.That(similarity, Is.EqualTo(0.0).Within(1e-10));
    }

    /// <summary>
    /// S6: Partial codon overlap produces intermediate similarity.
    /// Source: TVD formula. Exact calculation:
    ///   seq1: f(AUG)=2/3, f(GCU)=1/3
    ///   seq2: f(AUG)=2/3, f(UUU)=1/3
    ///   Σ|f₁-f₂| = |2/3-2/3| + |1/3-0| + |0-1/3| = 0 + 1/3 + 1/3 = 2/3
    ///   Similarity = 1 - (2/3)/2 = 1 - 1/3 = 2/3
    /// </summary>
    [Test]
    public void CompareCodonUsage_PartialOverlap_IntermediateSimilarity()
    {
        // Arrange: 2/3 codons shared (AUG), 1/3 differ (GCU vs UUU)
        const string seq1 = "AUGGCUAUG"; // AUG×2, GCU×1
        const string seq2 = "AUGUUUAUG"; // AUG×2, UUU×1

        // Act
        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);

        // Assert: exact value from TVD formula
        Assert.That(similarity, Is.EqualTo(2.0 / 3.0).Within(1e-10));
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
