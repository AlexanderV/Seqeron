using NUnit.Framework;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CodonOptimizer.CalculateCAI method.
/// Test Unit: CODON-CAI-001
/// Evidence: Sharp & Li (1987), Wikipedia CAI
/// </summary>
[TestFixture]
public class CodonOptimizer_CAI_Tests
{
    #region Edge Case Tests (Must)

    [Test]
    public void CalculateCAI_EmptySequence_ReturnsZero()
    {
        // Arrange & Act
        double cai = CodonOptimizer.CalculateCAI("", CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(0));
    }

    [Test]
    public void CalculateCAI_NullSequence_ReturnsZero()
    {
        // Arrange & Act
        double cai = CodonOptimizer.CalculateCAI(null!, CodonOptimizer.EColiK12);

        // Assert - Implementation returns 0 for null (safe handling)
        Assert.That(cai, Is.EqualTo(0));
    }

    #endregion

    #region Single-Codon Amino Acid Tests (Must)

    [Test]
    public void CalculateCAI_SingleMetCodon_ReturnsOne()
    {
        // Arrange - AUG is the only codon for Methionine
        // w = 1.0/1.0 = 1.0, CAI = 1.0^(1/1) = 1.0
        const string sequence = "AUG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateCAI_SingleTrpCodon_ReturnsOne()
    {
        // Arrange - UGG is the only codon for Tryptophan
        // w = 1.0/1.0 = 1.0, CAI = 1.0^(1/1) = 1.0
        const string sequence = "UGG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateCAI_MetAndTrp_ReturnsOne()
    {
        // Arrange - Both have unique codons with w=1.0
        // CAI = (1.0 × 1.0)^(1/2) = 1.0
        const string sequence = "AUGUGG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    #endregion

    #region Optimal Codon Tests (Must)

    [Test]
    public void CalculateCAI_AllOptimalCodonsEColi_ReturnsOne()
    {
        // Arrange - All optimal codons for E. coli K12
        // CUG (Leu, 0.47/0.47=1.0), CCG (Pro, 0.49/0.49=1.0), ACC (Thr, 0.40/0.40=1.0)
        // CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
        const string sequence = "CUGCCGACC";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void CalculateCAI_OptimalCodonsWithMet_ReturnsOne()
    {
        // Arrange - AUG (Met) + optimal codons
        // AUG(1.0), CUG(1.0), ACC(1.0) → CAI = 1.0
        const string sequence = "AUGCUGACC";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.01));
    }

    #endregion

    #region Rare Codon Tests (Must)

    [Test]
    public void CalculateCAI_RareCodonsEColi_ReturnsLow()
    {
        // Arrange - CUA (Leu rare, 0.04/0.47≈0.085)
        const string sequence = "CUACUACUA";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Geometric mean of ~0.085 values
        Assert.That(cai, Is.LessThan(0.15));
    }

    [Test]
    public void CalculateCAI_RareArginineCodonsEColi_ReturnsLow()
    {
        // Arrange - AGA, AGG are rare in E. coli (0.07, 0.04)
        // w_AGA = 0.07/0.36 ≈ 0.19, w_AGG = 0.04/0.36 ≈ 0.11
        const string sequence = "AGAAGG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.LessThan(0.25));
    }

    [Test]
    public void CalculateCAI_MixedOptimalAndRare_IntermediateValue()
    {
        // Arrange - Mix of optimal (CUG) and rare (CUA)
        const string sequence = "CUGCUA";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Between rare-only and optimal-only
        Assert.That(cai, Is.GreaterThan(0.1).And.LessThan(0.8));
    }

    #endregion

    #region Range Invariant Tests (Must)

    [Test]
    public void CalculateCAI_AnyValidSequence_RangeIsZeroToOne()
    {
        // Arrange
        string[] sequences =
        {
            "AUGGCUUAA",      // Short with stop
            "CUGCCGACC",       // Optimal
            "CUACUACUA",       // Rare
            "AUGUGGCUGACC",    // Mixed
            "AUGAAAGGGCCC"     // Random
        };

        // Act & Assert
        foreach (string sequence in sequences)
        {
            double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
            Assert.That(cai, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
                $"Sequence '{sequence}' produced CAI={cai} outside [0,1]");
        }
    }

    [Test]
    [TestCase("AUGGCU")]
    [TestCase("CUGCUGCUGCUGCUG")]
    [TestCase("AGAAGAAGAAGAAGA")]
    public void CalculateCAI_Parameterized_AlwaysInRange(string sequence)
    {
        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(cai, Is.GreaterThanOrEqualTo(0));
            Assert.That(cai, Is.LessThanOrEqualTo(1));
        });
    }

    #endregion

    #region Organism Specificity Tests (Must)

    [Test]
    public void CalculateCAI_SameSequence_DifferentOrganisms_DifferentResults()
    {
        // Arrange - CUG is preferred in E. coli but less so in Yeast
        const string sequence = "CUGCCGACC";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);

        // Assert - Different organisms have different preferences
        Assert.That(ecoliCai, Is.Not.EqualTo(yeastCai).Within(0.05));
    }

    [Test]
    public void CalculateCAI_YeastPreferredCodons_HigherInYeast()
    {
        // Arrange - UUA is preferred in yeast (0.28) but rare in E. coli (0.14)
        const string sequence = "UUAUUAUUA";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);

        // Assert - Should be higher in yeast
        Assert.That(yeastCai, Is.GreaterThan(ecoliCai));
    }

    [Test]
    public void CalculateCAI_AllThreeOrganismTables_ProduceValidResults()
    {
        // Arrange
        const string sequence = "AUGCUGCCGACC";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);
        double humanCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Human);

        // Assert - All should be valid
        Assert.Multiple(() =>
        {
            Assert.That(ecoliCai, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
            Assert.That(yeastCai, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
            Assert.That(humanCai, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
        });
    }

    #endregion

    #region Input Format Tests (Must)

    [Test]
    public void CalculateCAI_DnaInputWithThymine_ConvertsToUracil()
    {
        // Arrange - Same sequence in DNA (T) and RNA (U) format
        const string dnaSequence = "ATGCTG";
        const string rnaSequence = "AUGCUG";

        // Act
        double dnaCai = CodonOptimizer.CalculateCAI(dnaSequence, CodonOptimizer.EColiK12);
        double rnaCai = CodonOptimizer.CalculateCAI(rnaSequence, CodonOptimizer.EColiK12);

        // Assert - Should produce identical results
        Assert.That(dnaCai, Is.EqualTo(rnaCai).Within(0.001));
    }

    [Test]
    public void CalculateCAI_LowercaseInput_HandledCorrectly()
    {
        // Arrange
        const string lowercase = "augcug";
        const string uppercase = "AUGCUG";

        // Act
        double lowerCai = CodonOptimizer.CalculateCAI(lowercase, CodonOptimizer.EColiK12);
        double upperCai = CodonOptimizer.CalculateCAI(uppercase, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(lowerCai, Is.EqualTo(upperCai).Within(0.001));
    }

    [Test]
    public void CalculateCAI_MixedCaseInput_HandledCorrectly()
    {
        // Arrange
        const string mixed = "AuGcUg";
        const string standard = "AUGCUG";

        // Act
        double mixedCai = CodonOptimizer.CalculateCAI(mixed, CodonOptimizer.EColiK12);
        double standardCai = CodonOptimizer.CalculateCAI(standard, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(mixedCai, Is.EqualTo(standardCai).Within(0.001));
    }

    #endregion

    #region Stop Codon Tests (Must)

    [Test]
    public void CalculateCAI_SequenceWithStopCodon_ExcludesStopFromCalculation()
    {
        // Arrange - AUG + stop. Only AUG should be counted.
        const string withStop = "AUGUAA";
        const string withoutStop = "AUG";

        // Act
        double withStopCai = CodonOptimizer.CalculateCAI(withStop, CodonOptimizer.EColiK12);
        double withoutStopCai = CodonOptimizer.CalculateCAI(withoutStop, CodonOptimizer.EColiK12);

        // Assert - Should be same since stop is excluded
        Assert.That(withStopCai, Is.EqualTo(withoutStopCai).Within(0.001));
    }

    [Test]
    public void CalculateCAI_OnlyStopCodons_ReturnsZero()
    {
        // Arrange - All stop codons
        const string sequence = "UAAUAGUGA";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - No non-stop codons to evaluate
        Assert.That(cai, Is.EqualTo(0));
    }

    [Test]
    public void CalculateCAI_StopCodonInMiddle_ExcludedFromCalculation()
    {
        // Arrange - CUG + UAA + CCG (stop in middle)
        const string sequence = "CUGUAACCG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Stop excluded, CUG and CCG both optimal → high CAI
        Assert.That(cai, Is.GreaterThan(0.9));
    }

    #endregion

    #region Geometric Mean Property Tests (Must)

    [Test]
    public void CalculateCAI_SingleRareCodon_SignificantlyLowersCAI()
    {
        // Arrange - Geometric mean is sensitive to low values
        const string allOptimal = "CUGCUGCUGCUGCUG"; // 5 optimal Leu
        const string oneRare = "CUGCUGCUACUGCUG";   // 4 optimal + 1 rare

        // Act
        double optimalCai = CodonOptimizer.CalculateCAI(allOptimal, CodonOptimizer.EColiK12);
        double oneRareCai = CodonOptimizer.CalculateCAI(oneRare, CodonOptimizer.EColiK12);

        // Assert - Single rare codon drops CAI significantly
        Assert.That(oneRareCai, Is.LessThan(optimalCai * 0.75),
            "Single rare codon should significantly lower CAI due to geometric mean");
    }

    [Test]
    public void CalculateCAI_MoreRareCodons_LowerCAI()
    {
        // Arrange
        const string oneRare = "CUGCUGCUACUGCUG";   // 1 rare
        const string twoRare = "CUGCUACUACUGCUG";   // 2 rare
        const string threeRare = "CUGCUACUACUACUG"; // 3 rare

        // Act
        double oneRareCai = CodonOptimizer.CalculateCAI(oneRare, CodonOptimizer.EColiK12);
        double twoRareCai = CodonOptimizer.CalculateCAI(twoRare, CodonOptimizer.EColiK12);
        double threeRareCai = CodonOptimizer.CalculateCAI(threeRare, CodonOptimizer.EColiK12);

        // Assert - More rare codons → lower CAI
        Assert.That(oneRareCai, Is.GreaterThan(twoRareCai));
        Assert.That(twoRareCai, Is.GreaterThan(threeRareCai));
    }

    #endregion

    #region Hand-Calculated Verification Tests (Must)

    [Test]
    public void CalculateCAI_HandCalculatedRareCodons_MatchesExpected()
    {
        // Arrange - CUA (Leu rare) + ACU (Thr suboptimal)
        // CUA: w = 0.04/0.47 ≈ 0.085
        // ACU: w = 0.19/0.40 ≈ 0.475
        // CAI = (0.085 × 0.475)^(1/2) ≈ 0.20
        const string sequence = "CUAACU";
        const double expectedCai = 0.20;

        // Act
        double actualCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Within reasonable tolerance for frequency variations
        Assert.That(actualCai, Is.EqualTo(expectedCai).Within(0.10));
    }

    #endregion

    #region Performance Tests (Should)

    [Test]
    public void CalculateCAI_LongSequence_CompletesInReasonableTime()
    {
        // Arrange - 1000 codons (3000 nucleotides)
        string longSequence = string.Concat(Enumerable.Repeat("AUGCUG", 500));

        // Act
        var sw = System.Diagnostics.Stopwatch.StartNew();
        double cai = CodonOptimizer.CalculateCAI(longSequence, CodonOptimizer.EColiK12);
        sw.Stop();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(cai, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
            Assert.That(sw.ElapsedMilliseconds, Is.LessThan(1000), "Should complete in under 1 second");
        });
    }

    #endregion

    #region Incomplete Codon Tests (Should)

    [Test]
    public void CalculateCAI_IncompleteFinalCodon_IgnoredCorrectly()
    {
        // Arrange - AUG + incomplete (only complete codons counted)
        const string sequence = "AUGC"; // AUG + C (incomplete)

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Should be 1.0 (only AUG counted)
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateCAI_TwoIncompleteBases_IgnoredCorrectly()
    {
        // Arrange - CUG + incomplete
        const string sequence = "CUGCC"; // CUG + CC (incomplete)

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Should be 1.0 (only CUG counted, which is optimal)
        Assert.That(cai, Is.EqualTo(1.0).Within(0.01));
    }

    #endregion
}
