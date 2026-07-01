using NUnit.Framework;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

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

    #region Single-Codon Amino Acid Exclusion (Sharp & Li 1987 / Jansen 2003)

    // Source: Jansen, Bauer & Stadler (2003), Nucleic Acids Research — "An Improved
    // Implementation of the Codon Adaptation Index" (PMC2684136), quoting Sharp & Li (1987):
    // "The original paper proposing CAI (Sharp and Li, 1987) specifically stated that codon
    // families containing a single codon (e.g. AUG and UGG in the standard genetic code)
    // should be excluded in computing CAI" because "their corresponding w value will always
    // be 1 regardless of codon usage bias of the gene."
    // The exclusion is opt-in (excludeSingleCodonAminoAcids: true); default is unchanged.

    [Test]
    public void CalculateCAI_DefaultMode_IncludesSingleCodonAminoAcids()
    {
        // Arrange - default behaviour (excludeSingleCodonAminoAcids omitted) must be UNCHANGED:
        // AUG (Met, w=1.0) + UGG (Trp, w=1.0) → CAI = (1×1)^(1/2) = 1.0
        const string sequence = "AUGUGG";

        // Act
        double defaultCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double explicitInclude = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: false);

        // Assert - default == explicit-false == 1.0 (historical inclusive behaviour preserved)
        Assert.Multiple(() =>
        {
            Assert.That(defaultCai, Is.EqualTo(1.0).Within(1e-10),
                "Default mode must still include Met/Trp with w=1.0 (CAI=1.0)");
            Assert.That(explicitInclude, Is.EqualTo(defaultCai).Within(1e-10),
                "excludeSingleCodonAminoAcids:false must equal the default");
        });
    }

    [Test]
    public void CalculateCAI_ExcludeMode_AllSingleCodonAA_ReturnsZero()
    {
        // Arrange - AUG (Met) + UGG (Trp): both are single-codon AAs. When excluded, NO codons
        // remain in the geometric mean (L=0) → returns 0 (no codons to evaluate convention).
        // Per Sharp & Li 1987: a gene of only Met/Trp must not yield an inflated CAI of 1.0.
        const string sequence = "AUGUGG";

        // Act
        double cai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: true);

        // Assert - all codons excluded → 0 (contrast with default 1.0 proving the exclusion fired)
        Assert.That(cai, Is.EqualTo(0).Within(1e-10),
            "With Met/Trp excluded, an all-Met/Trp sequence has no scored codons → CAI=0");
    }

    [Test]
    public void CalculateCAI_ExcludeMode_DropsMetFromGeometricMean()
    {
        // Arrange - AUG (Met, single-codon) + CUA + CUA (Leu rare, w=0.04/0.50=0.08 each).
        // Inclusive (default): CAI = exp((ln1 + ln0.08 + ln0.08)/3) = 0.18566355334451112
        // Exclusive: Met dropped → geometric mean over the two CUA only:
        //   CAI = exp((ln0.08 + ln0.08)/2) = 0.08 exactly
        // Source: Kazusa E. coli K12 (species=316407); Sharp & Li (1987) exclusion rule.
        const string sequence = "AUGCUACUA";

        // Act
        double inclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double exclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: true);

        // Assert - exact independently-derived values
        Assert.Multiple(() =>
        {
            Assert.That(inclusive, Is.EqualTo(0.18566355334451112).Within(1e-10),
                "Default mode keeps Met (w=1.0) in the L=3 geometric mean");
            Assert.That(exclusive, Is.EqualTo(0.08).Within(1e-10),
                "Excluding Met leaves only the two CUA codons → geometric mean = 0.08");
        });
    }

    [Test]
    public void CalculateCAI_ExcludeMode_DropsBothMetAndTrp_ScoresOnlyRemainder()
    {
        // Arrange - AUG (Met) + UGG (Trp) + CUA (Leu rare, w=0.08). Both single-codon AAs excluded,
        // leaving exactly one scored codon CUA → CAI = exp(ln0.08 / 1) = 0.08 exactly.
        // Inclusive (default): exp((ln1 + ln1 + ln0.08)/3) = 0.43088693800637673.
        const string sequence = "AUGUGGCUA";

        // Act
        double inclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double exclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: true);

        // Assert - exact independently-derived values
        Assert.Multiple(() =>
        {
            Assert.That(inclusive, Is.EqualTo(0.43088693800637673).Within(1e-10),
                "Default mode keeps Met and Trp (w=1.0 each) in the L=3 geometric mean");
            Assert.That(exclusive, Is.EqualTo(0.08).Within(1e-10),
                "Excluding Met and Trp leaves only CUA → CAI=0.08");
        });
    }

    [Test]
    public void CalculateCAI_ExcludeMode_NoSingleCodonAA_UnchangedFromDefault()
    {
        // Arrange - CUGCUA contains no Met/Trp, so exclusion must not alter the result.
        // CUG (w=1.0) + CUA (w=0.08) → CAI = (1×0.08)^(1/2) = 0.28284271247461906
        const string sequence = "CUGCUA";

        // Act
        double inclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: false);
        double exclusive = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12, excludeSingleCodonAminoAcids: true);

        // Assert - identical when the sequence has no single-codon amino acids
        Assert.Multiple(() =>
        {
            Assert.That(inclusive, Is.EqualTo(0.28284271247461906).Within(1e-10));
            Assert.That(exclusive, Is.EqualTo(inclusive).Within(1e-10),
                "With no Met/Trp present, the exclusion flag must not change CAI");
        });
    }

    #endregion

    #region Zero-Frequency / No-Data Reference Codon Tests (Must)

    // These cover the implementation's handling of partial reference tables (gaps), which
    // Sharp & Li (1987) did not encounter (they used complete reference sets). They lock the
    // two documented fallbacks so they cannot regress:
    //   (a) codon present-in-sequence but ABSENT from the table while another synonymous codon
    //       IS present (maxFreq > 0, f = 0): w is clamped to 1e-6 (avoids ln(0) = -inf).
    //   (b) amino acid with NO frequency data at all (maxFreq <= 0): w = NaN -> codon skipped.

    [Test]
    public void CalculateCAI_AbsentCodonWithPresentSynonym_ClampsWeightToEpsilon()
    {
        // Arrange - custom table built from a reference of only "CUG" (Leu): CUG freq = 1.0,
        // all other Leu codons absent. Score "CUACUG": CUA is absent (f=0) but Leu's maxFreq
        // is 1.0 (CUG present) -> w_CUA = max(0/1.0, 1e-6) = 1e-6; w_CUG = 1.0.
        // CAI = exp((ln(1e-6) + ln(1.0)) / 2) = 0.001 exactly.
        var partial = CodonOptimizer.CreateCodonTableFromSequence("CUG", "partial-Leu");

        // Act
        double cai = CodonOptimizer.CalculateCAI("CUACUG", partial);

        // Assert - clamp keeps the value finite and bounded, not 0 and not NaN.
        Assert.That(cai, Is.EqualTo(0.001).Within(1e-9),
            "Absent codon with a present synonym must clamp w to 1e-6, not collapse to 0 or NaN");
    }

    [Test]
    public void CalculateCAI_AllCodonsAbsentFromFamily_ClampsToEpsilon()
    {
        // Arrange - only CUA in the sequence, table has CUG=1.0 (CUA absent).
        // Single codon, w clamped to 1e-6 -> CAI = exp(ln(1e-6)/1) = 1e-6.
        var partial = CodonOptimizer.CreateCodonTableFromSequence("CUG", "partial-Leu");

        // Act
        double cai = CodonOptimizer.CalculateCAI("CUA", partial);

        // Assert
        Assert.That(cai, Is.EqualTo(1e-6).Within(1e-12),
            "A lone absent codon (synonym present) must yield the 1e-6 clamp, never ln(0)");
    }

    [Test]
    public void CalculateCAI_AminoAcidWithNoFrequencyData_IsSkipped()
    {
        // Arrange - table built from only "CUG" (Leu) has NO Phe data. Score "UUUCUG":
        // UUU (Phe) has maxFreq = 0 in this table -> w = NaN -> skipped (not counted in L).
        // Only CUG (w = 1.0) remains -> CAI = 1.0.
        var partial = CodonOptimizer.CreateCodonTableFromSequence("CUG", "partial-Leu");

        // Act
        double cai = CodonOptimizer.CalculateCAI("UUUCUG", partial);

        // Assert - the no-data amino acid is dropped, not treated as w=0.
        Assert.That(cai, Is.EqualTo(1.0).Within(1e-12),
            "An amino acid with no table frequency data must be skipped, leaving only CUG (w=1.0)");
    }

    [Test]
    public void CalculateCAI_AllCodonsHaveNoFrequencyData_ReturnsZero()
    {
        // Arrange - table built from only "CUG" (Leu); sequence is all Phe (UUU), which has
        // no data in this table -> every codon NaN-skipped -> count=0 -> returns 0.
        var partial = CodonOptimizer.CreateCodonTableFromSequence("CUG", "partial-Leu");

        // Act
        double cai = CodonOptimizer.CalculateCAI("UUUUUU", partial);

        // Assert
        Assert.That(cai, Is.EqualTo(0),
            "When no codon has table data, count=0 and CAI is 0 by the no-codons convention");
    }

    #endregion
}
