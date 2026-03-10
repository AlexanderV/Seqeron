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
        // Arrange - Both single-codon amino acids: AUG (Met, w=1.0) + UGG (Trp, w=1.0)
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
        // Arrange - All optimal codons for E. coli K12 (Kazusa MG1655, species=316407)
        // CUG (Leu, 0.50/0.50=1.0), CCG (Pro, 0.53/0.53=1.0), ACC (Thr, 0.44/0.44=1.0)
        // CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
        const string sequence = "CUGCCGACC";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateCAI_OptimalCodonsWithMet_ReturnsOne()
    {
        // Arrange - AUG (Met, w=1.0) + CUG (Leu optimal, w=1.0) + ACC (Thr optimal, w=1.0)
        // CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
        const string sequence = "AUGCUGACC";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    #endregion

    #region Rare Codon Tests (Must)

    [Test]
    public void CalculateCAI_RareCodonsEColi_ReturnsLow()
    {
        // Arrange - CUA (Leu rare): w = 0.04/0.50 = 0.08
        // All 3 codons identical → CAI = 0.08
        // Source: Kazusa MG1655 (species=316407)
        const string sequence = "CUACUACUA";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Exact: geometric mean of (0.08, 0.08, 0.08) = 0.08
        Assert.That(cai, Is.EqualTo(0.08).Within(0.005));
    }

    [Test]
    public void CalculateCAI_RareArginineCodonsEColi_MatchesHandCalculated()
    {
        // Arrange - AGA, AGG are rare in E. coli K12
        // w_AGA = 0.04/0.40 = 0.10, w_AGG = 0.02/0.40 = 0.05
        // CAI = (0.10 × 0.05)^(1/2) = 0.07071
        // Source: Kazusa MG1655 (species=316407)
        const string sequence = "AGAAGG";
        const double expectedCai = 0.07071;

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(expectedCai).Within(0.005));
    }

    [Test]
    public void CalculateCAI_MixedOptimalAndRare_MatchesHandCalculated()
    {
        // Arrange - CUG (optimal Leu, w=1.0) + CUA (rare Leu, w=0.04/0.50=0.08)
        // CAI = (1.0 × 0.08)^(1/2) = 0.28284
        // Source: Kazusa MG1655 (species=316407)
        const string sequence = "CUGCUA";
        const double expectedCai = 0.28284;

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(cai, Is.EqualTo(expectedCai).Within(0.005));
    }

    #endregion

    #region Range Invariant Tests (Must)

    [Test]
    public void CalculateCAI_AnyValidSequence_RangeIsZeroToOne()
    {
        // Arrange - diverse sequences: optimal, rare, mixed, with stop
        string[] sequences =
        {
            "AUGGCUUAA",      // Short with stop
            "CUGCCGACC",       // Optimal
            "CUACUACUA",       // Rare
            "AUGUGGCUGACC",    // Mixed
            "AUGAAAGGGCCC",    // Random
            "AGAAGAAGAAGAAGA", // Rare arginine
            "AUGGCU"            // DNA-like
        };

        // Act & Assert
        foreach (string sequence in sequences)
        {
            double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
            Assert.That(cai, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
                $"Sequence '{sequence}' produced CAI={cai} outside [0,1]");
        }
    }

    #endregion

    #region Organism Specificity Tests (Must)

    [Test]
    public void CalculateCAI_SameSequence_DifferentOrganisms_DifferentResults()
    {
        // Arrange - CUG-CCG-ACC: optimal for E. coli, suboptimal for yeast
        // Yeast: CUG(0.11/0.29=0.3793), CCG(0.12/0.42=0.2857), ACC(0.22/0.35=0.6286)
        //   CAI = exp((ln(0.3793)+ln(0.2857)+ln(0.6286))/3) = 0.4085
        // Source: Kazusa species=316407 and species=4932
        const string sequence = "CUGCCGACC";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);

        // Assert - Exact hand-calculated values
        Assert.That(ecoliCai, Is.EqualTo(1.0).Within(0.001));
        Assert.That(yeastCai, Is.EqualTo(0.4085).Within(0.005));
    }

    [Test]
    public void CalculateCAI_YeastPreferredCodons_HigherInYeast()
    {
        // Arrange - UUA: E. coli w=0.13/0.50=0.26; Yeast w=0.28/0.29=0.9655
        // E. coli CAI = 0.26; Yeast CAI = 0.9655
        // Source: Kazusa species=316407 and species=4932
        const string sequence = "UUAUUAUUA";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);

        // Assert - Exact values from Kazusa-verified tables
        Assert.That(ecoliCai, Is.EqualTo(0.26).Within(0.005));
        Assert.That(yeastCai, Is.EqualTo(0.9655).Within(0.005));
    }

    [Test]
    public void CalculateCAI_AllThreeOrganismTables_MatchHandCalculated()
    {
        // Arrange - AUG-CUG-CCG-ACC
        // E. coli: all w=1.0 → CAI=1.0
        // Yeast: AUG(1.0), CUG(0.11/0.29=0.3793), CCG(0.12/0.42=0.2857), ACC(0.22/0.35=0.6286)
        //   CAI = exp((0 + ln(0.3793) + ln(0.2857) + ln(0.6286))/4) = 0.5109
        // Human: AUG(1.0), CUG(0.40/0.40=1.0), CCG(0.11/0.32=0.34375), ACC(0.36/0.36=1.0)
        //   CAI = exp((0 + 0 + ln(0.34375) + 0)/4) = 0.7656
        // Source: Kazusa species=316407, 4932, 9606
        const string sequence = "AUGCUGCCGACC";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double yeastCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Yeast);
        double humanCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Human);

        // Assert - Exact values from Kazusa-verified tables
        Assert.Multiple(() =>
        {
            Assert.That(ecoliCai, Is.EqualTo(1.0).Within(0.001));
            Assert.That(yeastCai, Is.EqualTo(0.5109).Within(0.005));
            Assert.That(humanCai, Is.EqualTo(0.7656).Within(0.005));
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

        // Assert - Both should produce CAI=1.0 (AUG w=1.0, CUG w=1.0)
        Assert.That(lowerCai, Is.EqualTo(1.0).Within(0.001));
        Assert.That(upperCai, Is.EqualTo(1.0).Within(0.001));
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
        // Arrange - CUG + UAA + CCG (stop in middle excluded)
        // Only CUG (w=1.0) and CCG (w=1.0) counted → CAI = 1.0
        const string sequence = "CUGUAACCG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Exact: both non-stop codons are optimal
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    #endregion

    #region Geometric Mean Property Tests (Must)

    [Test]
    public void CalculateCAI_SingleRareCodon_SignificantlyLowersCAI()
    {
        // Arrange - Geometric mean is highly sensitive to low values (Sharp & Li 1987)
        // 5×CUG: all w=1.0, CAI=1.0
        // 4×CUG+1×CUA: w=[1,1,1,1,0.08]
        //   CAI = exp((4×ln(1)+ln(0.08))/5) = exp(-2.52573/5) = 0.60360
        const string allOptimal = "CUGCUGCUGCUGCUG";
        const string oneRare = "CUGCUGCUACUGCUG";

        // Act
        double optimalCai = CodonOptimizer.CalculateCAI(allOptimal, CodonOptimizer.EColiK12);
        double oneRareCai = CodonOptimizer.CalculateCAI(oneRare, CodonOptimizer.EColiK12);

        // Assert - Exact hand-calculated values
        Assert.That(optimalCai, Is.EqualTo(1.0).Within(0.001));
        Assert.That(oneRareCai, Is.EqualTo(0.6036).Within(0.005));
    }

    [Test]
    public void CalculateCAI_MoreRareCodons_LowerCAI()
    {
        // Arrange - Monotonicity: more rare codons → lower CAI
        // CUA: w=0.04/0.50=0.08; CUG: w=1.0
        // 1 rare: exp((4×0+ln(0.08))/5) = 0.60360
        // 2 rare: exp((3×0+2×ln(0.08))/5) = 0.36434
        // 3 rare: exp((2×0+3×ln(0.08))/5) = 0.21952
        const string oneRare = "CUGCUGCUACUGCUG";
        const string twoRare = "CUGCUACUACUGCUG";
        const string threeRare = "CUGCUACUACUACUG";

        // Act
        double oneRareCai = CodonOptimizer.CalculateCAI(oneRare, CodonOptimizer.EColiK12);
        double twoRareCai = CodonOptimizer.CalculateCAI(twoRare, CodonOptimizer.EColiK12);
        double threeRareCai = CodonOptimizer.CalculateCAI(threeRare, CodonOptimizer.EColiK12);

        // Assert - Exact hand-calculated values + monotonicity
        Assert.Multiple(() =>
        {
            Assert.That(oneRareCai, Is.EqualTo(0.6036).Within(0.005));
            Assert.That(twoRareCai, Is.EqualTo(0.3643).Within(0.005));
            Assert.That(threeRareCai, Is.EqualTo(0.2195).Within(0.005));
            Assert.That(oneRareCai, Is.GreaterThan(twoRareCai));
            Assert.That(twoRareCai, Is.GreaterThan(threeRareCai));
        });
    }

    #endregion

    #region Hand-Calculated Verification Tests (Must)

    [Test]
    public void CalculateCAI_HandCalculatedRareCodons_MatchesExpected()
    {
        // Arrange - CUA (Leu rare) + ACU (Thr suboptimal)
        // CUA: w = 0.04/0.50 = 0.08
        // ACU: w = 0.16/0.44 = 0.36364
        // CAI = (0.08 × 0.36364)^(1/2) = 0.17056
        // Source: Kazusa MG1655 (species=316407), verified against https://www.kazusa.or.jp/codon/
        const string sequence = "CUAACU";
        const double expectedCai = 0.17056;

        // Act
        double actualCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);

        // Assert - Tight tolerance: values derived directly from Kazusa-verified table
        Assert.That(actualCai, Is.EqualTo(expectedCai).Within(0.005));
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
        Assert.That(cai, Is.EqualTo(1.0).Within(0.001));
    }

    #endregion
}
