namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for CODON-OPT-001: Sequence Optimization
/// Canonical method: CodonOptimizer.OptimizeSequence(...)
/// 
/// Evidence: Wikipedia (Codon usage bias, Codon Adaptation Index),
/// Sharp &amp; Li (1987), Plotkin &amp; Kudla (2011)
/// </summary>
[TestFixture]
[Category("Codon")]
[Category("CODON-OPT-001")]
public class CodonOptimizer_OptimizeSequence_Tests
{
    #region MUST Tests - Protein Preservation

    [Test]
    [Description("M1: Protein must be preserved across all optimization strategies")]
    public void OptimizeSequence_PreservesProtein_AllStrategies()
    {
        // Arrange - sequence with multiple amino acids: M-A-L-R-Stop
        const string sequence = "AUGGCUCUAAGAUAA";
        var strategies = new[]
        {
            CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            CodonOptimizer.OptimizationStrategy.HarmonizeExpression,
            CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
            CodonOptimizer.OptimizationStrategy.AvoidRareCodeons
        };

        foreach (var strategy in strategies)
        {
            // Act
            var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12, strategy);

            // Assert
            Assert.Multiple(() =>
            {
                Assert.That(result.ProteinSequence, Is.EqualTo("MALR*"),
                    $"Protein must be preserved for strategy {strategy}");
                Assert.That(result.OptimizedSequence.Length, Is.EqualTo(result.OriginalSequence.Length),
                    $"Sequence length must be preserved for strategy {strategy}");
            });
        }
    }

    [Test]
    [Description("M1: Protein preservation for longer sequences")]
    public void OptimizeSequence_LongerSequence_PreservesProtein()
    {
        // Arrange - GFP-like sequence start: M-S-K-G-E-E-L-F
        const string sequence = "AUGAGCAAAGGTGAAGAACUGUUC";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - verify protein is unchanged
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("MSKGEELF"));
            Assert.That(result.OptimizedSequence.Length % 3, Is.EqualTo(0));
        });
    }

    #endregion

    #region MUST Tests - Empty and Edge Cases

    [Test]
    [Description("M2: Empty sequence returns empty result with CAI=0")]
    public void OptimizeSequence_EmptySequence_ReturnsEmptyResult()
    {
        // Act
        var result = CodonOptimizer.OptimizeSequence("", CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.Empty);
            Assert.That(result.ProteinSequence, Is.Empty);
            Assert.That(result.OriginalCAI, Is.EqualTo(0));
            Assert.That(result.OptimizedCAI, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("M3: DNA thymine (T) is converted to RNA uracil (U)")]
    public void OptimizeSequence_DnaInput_ConvertsThymineToUracil()
    {
        // Arrange - DNA sequence with T
        const string dnaSequence = "ATGGCTTAA";

        // Act
        var result = CodonOptimizer.OptimizeSequence(dnaSequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalSequence, Does.Contain("U"));
            Assert.That(result.OriginalSequence, Does.Not.Contain("T"));
            Assert.That(result.OptimizedSequence, Does.Not.Contain("T"));
            Assert.That(result.ProteinSequence, Is.EqualTo("MA*"));
        });
    }

    [Test]
    [Description("M4: Incomplete codons are trimmed to complete triplets")]
    public void OptimizeSequence_IncompleteCodons_TrimsToComplete()
    {
        // Arrange - 10 nucleotides (not divisible by 3)
        const string sequence = "AUGGCUUAAG"; // 10 nt -> should trim to 9

        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalSequence.Length % 3, Is.EqualTo(0));
            Assert.That(result.OriginalSequence.Length, Is.EqualTo(9));
            Assert.That(result.OptimizedSequence.Length % 3, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("M4: Single extra nucleotide is trimmed")]
    public void OptimizeSequence_SingleExtraNucleotide_Trimmed()
    {
        // Arrange - 7 nucleotides
        const string sequence = "AUGGCUA"; // 7 nt -> should trim to 6

        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(result.OriginalSequence.Length, Is.EqualTo(6));
    }

    #endregion

    #region MUST Tests - CAI Behavior

    [Test]
    [Description("M5: MaximizeCAI strategy increases or maintains CAI")]
    public void OptimizeSequence_MaximizeCAI_IncreasesOrMaintainsCAI()
    {
        // Arrange - sequence with rare E. coli codons (CUA=Leu rare, AGA=Arg rare)
        const string rareCodeSequence = "CUAAGACGA"; // L-R-R with rare codons

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            rareCodeSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - optimized CAI should be >= original
        Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI),
            "MaximizeCAI should improve or maintain CAI");
    }

    [Test]
    [Description("M5: MaximizeCAI on already optimal sequence maintains CAI=1.0")]
    public void OptimizeSequence_AlreadyOptimal_MaintainsHighCAI()
    {
        // Arrange - E. coli optimal codons (Kazusa species=316407):
        // CUG = best Leu (0.50/0.50=1.0), CCG = best Pro (0.53/0.53=1.0), ACC = best Thr (0.44/0.44=1.0)
        const string optimalSequence = "CUGCCGACC";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            optimalSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - all codons already optimal → CAI = exp((ln(1)+ln(1)+ln(1))/3) = 1.0
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalCAI, Is.EqualTo(1.0).Within(0.001),
                "All-optimal E. coli codons should have CAI = 1.0");
            Assert.That(result.OptimizedCAI, Is.EqualTo(1.0).Within(0.001));
            Assert.That(result.ChangedCodons, Is.EqualTo(0),
                "No codons should be changed when all are already optimal");
        });
    }

    #endregion

    #region MUST Tests - Special Codons

    [Test]
    [Description("M6: Methionine (AUG) is the only codon - cannot be changed")]
    public void OptimizeSequence_MethionineCodon_Unchanged()
    {
        // Arrange - only methionine
        const string sequence = "AUG";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUG"));
            Assert.That(result.ProteinSequence, Is.EqualTo("M"));
            Assert.That(result.OptimizedCAI, Is.EqualTo(0),
                "AUG (Met) is not scored in CAI (Sharp & Li 1987 / CodonW): no scored codon → 0");
        });
    }

    [Test]
    [Description("M6: Tryptophan (UGG) is the only codon - cannot be changed")]
    public void OptimizeSequence_TryptophanCodon_Unchanged()
    {
        // Arrange - Met + Trp
        const string sequence = "AUGUGG";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGUGG"));
            Assert.That(result.ProteinSequence, Is.EqualTo("MW"));
        });
    }

    [Test]
    [Description("M7: Stop codons are preserved (not optimized away)")]
    public void OptimizeSequence_StopCodons_Preserved()
    {
        // Arrange - sequence ending with stop codon
        const string sequence = "AUGUAA"; // M-Stop (UAA)

        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Does.EndWith("*"));
            Assert.That(result.OptimizedSequence,
                Does.EndWith("UAA").Or.EndWith("UAG").Or.EndWith("UGA"),
                "Stop codon must remain a stop codon");
        });
    }

    [Test]
    [Description("M7: All three stop codons are preserved correctly")]
    [TestCase("AUGUAA", "M*")] // UAA stop
    [TestCase("AUGUAG", "M*")] // UAG stop
    [TestCase("AUGUGA", "M*")] // UGA stop
    public void OptimizeSequence_AllStopCodons_CorrectProtein(string sequence, string expectedProtein)
    {
        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.That(result.ProteinSequence, Is.EqualTo(expectedProtein));
    }

    #endregion

    #region MUST Tests - Organism Specificity

    [Test]
    [Description("M8: Different organisms produce different optimizations")]
    public void OptimizeSequence_DifferentOrganisms_DifferentResults()
    {
        // Arrange - CUG(Leu): E.coli best (0.50), Yeast low (0.11, best=UUG 0.29)
        //           AGA(Arg): Yeast best (0.48), E.coli rare (0.04, best=CGC 0.40)
        const string sequence = "CUGAGA"; // L-R

        // Act
        var ecoliResult = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var yeastResult = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.Yeast, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - same protein, different codon choices per Kazusa data
        Assert.Multiple(() =>
        {
            Assert.That(ecoliResult.ProteinSequence, Is.EqualTo("LR"));
            Assert.That(yeastResult.ProteinSequence, Is.EqualTo("LR"));
            // E. coli: CUG stays (best Leu, 0.50), AGA→CGC (best Arg, 0.40)
            Assert.That(ecoliResult.OptimizedSequence, Is.EqualTo("CUGCGC"),
                "E. coli prefers CUG(Leu) and CGC(Arg) per Kazusa 316407");
            // Yeast: CUG→UUG (best Leu, 0.29), AGA stays (best Arg, 0.48)
            Assert.That(yeastResult.OptimizedSequence, Is.EqualTo("UUGAGA"),
                "Yeast prefers UUG(Leu) and AGA(Arg) per Kazusa 4932");
        });
    }

    [Test]
    [Description("M8: Human optimization differs from E. coli")]
    public void OptimizeSequence_HumanVsEcoli_DifferentCAI()
    {
        // Arrange - CUGCCGACC = L-P-T
        // E. coli (Kazusa 316407): CUG=0.50/0.50=1.0, CCG=0.53/0.53=1.0, ACC=0.44/0.44=1.0 → CAI=1.0
        // Human (Kazusa 9606):     CUG=0.40/0.40=1.0, CCG=0.11/0.32=0.34, ACC=0.36/0.36=1.0 → CAI≈0.700
        const string sequence = "CUGCCGACC";

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double humanCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Human);

        // Assert - verified against hand-computed values from Kazusa data
        Assert.Multiple(() =>
        {
            Assert.That(ecoliCai, Is.EqualTo(1.0).Within(0.001),
                "All codons are optimal for E. coli");
            Assert.That(humanCai, Is.EqualTo(0.700).Within(0.01),
                "exp((ln(1)+ln(0.11/0.32)+ln(1))/3) ≈ 0.700");
        });
    }

    #endregion

    #region MUST Tests - Input Handling

    [Test]
    [Description("M9: Lowercase input is handled correctly")]
    public void OptimizeSequence_LowercaseInput_HandledCorrectly()
    {
        // Arrange
        const string lowercaseSequence = "auggcuuaa";

        // Act
        var result = CodonOptimizer.OptimizeSequence(lowercaseSequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("MA*"));
            Assert.That(result.OriginalSequence, Is.EqualTo("AUGGCUUAA"));
        });
    }

    [Test]
    [Description("M10: OptimizationResult has all fields populated correctly")]
    public void OptimizeSequence_ReturnsValidOptimizationResult()
    {
        // Arrange - AUGGCUUAA = M-A-Stop, default BalancedOptimization
        // GCU(A, 0.16) → GCG(A, 0.36 = best above 0.15 threshold)
        const string sequence = "AUGGCUUAA";

        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert - all values hand-computed from Kazusa E. coli K12 (species=316407)
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalSequence, Is.EqualTo("AUGGCUUAA"));
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGGCGUAA"),
                "GCU→GCG (best Ala codon above 0.15 threshold)");
            Assert.That(result.ProteinSequence, Is.EqualTo("MA*"));
            Assert.That(result.OriginalCAI, Is.EqualTo(0.16 / 0.36).Within(1e-12),
                "AUG not scored (Sharp & Li / CodonW); wi(GCU)=0.16/0.36 → CAI=0.444 (corrected 2026-09 from 0.667)");
            Assert.That(result.OptimizedCAI, Is.EqualTo(1.0).Within(0.001),
                "wi(GCG)=0.36/0.36=1.0 → CAI=1.0");
            Assert.That(result.GcContentOriginal, Is.EqualTo(3.0 / 9.0).Within(0.001));
            Assert.That(result.GcContentOptimized, Is.EqualTo(4.0 / 9.0).Within(0.001));
            Assert.That(result.ChangedCodons, Is.EqualTo(1));
            Assert.That(result.Changes, Has.Count.EqualTo(1));
            Assert.That(result.Changes[0], Is.EqualTo((3, "GCU", "GCG")));
        });
    }

    #endregion

    #region SHOULD Tests - Strategy Specifics

    [Test]
    [Description("S1: AvoidRareCodons only replaces codons below threshold")]
    public void OptimizeSequence_AvoidRareCodons_OnlyReplacesRare()
    {
        // Arrange - CUG is common (0.50), CUA is rare (0.04) in E. coli (Kazusa species=316407)
        const string mixedSequence = "CUGCUA"; // L-L with good-rare for E. coli

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            mixedSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.AvoidRareCodeons);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("LL"));
            // CUG (0.50) is above threshold → must stay unchanged
            Assert.That(result.OptimizedSequence.Substring(0, 3), Is.EqualTo("CUG"),
                "Common codon CUG (freq=0.50) must remain unchanged");
            // CUA (0.04) is below threshold → must be replaced with a non-rare Leu codon
            Assert.That(result.OptimizedSequence.Substring(3, 3), Is.Not.EqualTo("CUA"),
                "Rare codon CUA (freq=0.04) must be replaced");
            Assert.That(result.ChangedCodons, Is.EqualTo(1),
                "Exactly one codon (CUA) should be changed");
        });
    }

    [Test]
    [Description("S2: BalancedOptimization targets GC content within 40-60% range")]
    public void OptimizeSequence_BalancedOptimization_ConsidersGcContent()
    {
        // Arrange - AUGGCCGCC = M-A-A → GC = 7/9 = 0.778 (above 60% target)
        // Phase 1: GCC→GCG (best Ala), GCC→GCG → AUGGCGGCG (GC still 7/9)
        // Phase 2: BalanceGcContent reduces GC by substituting lower-GC synonymous codons
        const string sequence = "AUGGCCGCC";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.GcContentOriginal, Is.EqualTo(7.0 / 9.0).Within(0.001),
                "Original GC = 7/9 ≈ 0.778");
            Assert.That(result.GcContentOptimized, Is.InRange(0.40, 0.60),
                "GC content must be in target range after balancing");
            Assert.That(result.ProteinSequence, Is.EqualTo("MAA"),
                "Protein must be preserved");
        });
    }

    [Test]
    [Description("S3: Changes list accurately tracks modifications")]
    public void OptimizeSequence_TracksChanges_Correctly()
    {
        // Arrange - rare codons that should be changed
        const string rareSequence = "CUAAGACGA"; // L-R-R with rare codons

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            rareSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Changes, Is.Not.Null);
            Assert.That(result.ChangedCodons, Is.EqualTo(result.Changes.Count));
            // Each change should have valid position
            foreach (var change in result.Changes)
            {
                Assert.That(change.Position % 3, Is.EqualTo(0),
                    "Change position should be at codon boundary");
                Assert.That(change.Original.Length, Is.EqualTo(3));
                Assert.That(change.Optimized.Length, Is.EqualTo(3));
            }
        });
    }

    [Test]
    [Description("S4: Long sequences complete in reasonable time")]
    public void OptimizeSequence_LongSequence_Completes()
    {
        // Arrange - generate long sequence (500 codons)
        var codons = new[] { "AUG", "GCU", "CUG", "AAA", "GAU" };
        string longSequence = string.Join("", Enumerable.Repeat(codons, 100).SelectMany(x => x)) + "UAA";

        // Act
        var result = CodonOptimizer.OptimizeSequence(longSequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence.Length, Is.EqualTo(longSequence.Length));
            Assert.That(result.ProteinSequence.Length, Is.EqualTo(501)); // 500 + stop
        });
    }

    [Test]
    [Description("S5: GFP-like sequence optimizes correctly")]
    public void OptimizeSequence_GfpSequence_OptimizesCorrectly()
    {
        // Arrange - part of GFP: MSKG
        const string gfpPart = "AUGAGCAAAGGU";

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            gfpPart,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("MSKG"));
            Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI));
        });
    }

    #endregion

    #region COULD Tests - Advanced Features

    [Test]
    [Description("C1: HarmonizeExpression maintains codon distribution")]
    public void OptimizeSequence_HarmonizeExpression_MaintainsDistribution()
    {
        // Arrange - repetitive sequence
        const string sequence = "CUGCUGCUGCUGCUGCUG"; // 6x Leu

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.HarmonizeExpression);

        // Assert - protein preserved, CAI valid (stochastic strategy)
        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("LLLLLL"));
            Assert.That(result.OptimizedCAI, Is.GreaterThan(0).And.LessThanOrEqualTo(1),
                "Weighted random selection must produce valid codons");
        });
    }

    #endregion

    #region Invariant Tests

    [Test]
    [Description("Invariant: Optimized sequence encodes same protein for all inputs")]
    [TestCase("AUGGCUUAA")]           // M-A-Stop
    [TestCase("CUGCUACUUCUCCUA")]     // All Leu synonyms
    [TestCase("CGUCGCCGACGGAGAAGG")]  // All Arg synonyms
    [TestCase("AUGUAA")]              // Start-Stop
    [TestCase("AUGUGGUAA")]           // M-W-Stop (unique codons)
    public void OptimizeSequence_Invariant_ProteinPreserved(string sequence)
    {
        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Translate original and optimized independently
        var originalProtein = TranslateSequence(result.OriginalSequence);
        var optimizedProtein = TranslateSequence(result.OptimizedSequence);

        // Assert
        Assert.That(optimizedProtein, Is.EqualTo(originalProtein),
            $"Optimized sequence must encode same protein. Original: {result.OriginalSequence}, Optimized: {result.OptimizedSequence}");
    }

    [Test]
    [Description("Invariant: CAI is always in valid range (0,1]")]
    [TestCase("AUGGCUUAA")]
    [TestCase("CUGCUGCUGCUGCUG")]
    [TestCase("AGAAGAAGAAGA")]  // Rare codons
    public void OptimizeSequence_Invariant_CAIInValidRange(string sequence)
    {
        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalCAI, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
            Assert.That(result.OptimizedCAI, Is.GreaterThan(0).And.LessThanOrEqualTo(1));
        });
    }

    [Test]
    [Description("Invariant: CAI matches hand-computed Sharp & Li (1987) formula using Kazusa data")]
    public void OptimizeSequence_CAI_MatchesSharpLiFormula()
    {
        // Hand-computed CAI for CUAAGACGA in E. coli (Kazusa species=316407):
        //   CUA(L): freq=0.04, max(L)=CUG=0.50, wi=0.04/0.50=0.08
        //   AGA(R): freq=0.04, max(R)=CGC=0.40, wi=0.04/0.40=0.10
        //   CGA(R): freq=0.06, max(R)=CGC=0.40, wi=0.06/0.40=0.15
        //   CAI = exp((ln(0.08)+ln(0.10)+ln(0.15))/3)
        //       = exp((-2.5257+-2.3026+-1.8971)/3) = exp(-2.2418) ≈ 0.1063
        const string rareSequence = "CUAAGACGA";
        double expectedCAI = System.Math.Exp(
            (System.Math.Log(0.04 / 0.50) + System.Math.Log(0.04 / 0.40) + System.Math.Log(0.06 / 0.40)) / 3.0);

        double actualCAI = CodonOptimizer.CalculateCAI(rareSequence, CodonOptimizer.EColiK12);

        Assert.That(actualCAI, Is.EqualTo(expectedCAI).Within(0.001),
            "CAI must match Sharp & Li (1987) geometric mean formula");
    }

    [Test]
    [Description("Invariant: MaximizeCAI turns rare codons into optimal → CAI=1.0")]
    public void OptimizeSequence_MaximizeCAI_ProducesOptimalCAI()
    {
        // Arrange - all rare E. coli codons (Kazusa species=316407)
        const string rareSequence = "CUAAGACGA"; // L-R-R, original CAI ≈ 0.106

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            rareSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - MaximizeCAI selects most frequent codon for each AA → wi=1.0 for all → CAI=1.0
        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedCAI, Is.EqualTo(1.0).Within(0.001),
                "MaximizeCAI should produce CAI=1.0 (all codons optimal)");
            Assert.That(result.OriginalCAI, Is.LessThan(0.2),
                "Original rare codons should have very low CAI");
            // Verify specific codon replacements per Kazusa data:
            // CUA(L) → CUG(L, freq=0.50 = max)
            // AGA(R) → CGC(R, freq=0.40 = max)
            // CGA(R) → CGC(R, freq=0.40 = max)
            Assert.That(result.OptimizedSequence.Substring(0, 3), Is.EqualTo("CUG"),
                "CUA should become CUG (best E. coli Leu codon, freq=0.50)");
            Assert.That(result.OptimizedSequence.Substring(3, 3), Is.EqualTo("CGC"),
                "AGA should become CGC (best E. coli Arg codon, freq=0.40)");
            Assert.That(result.OptimizedSequence.Substring(6, 3), Is.EqualTo("CGC"),
                "CGA should become CGC (best E. coli Arg codon, freq=0.40)");
        });
    }

    #endregion

    #region Reference cross-checks (review 2026-09)

    [Test]
    [Description("MaximizeCAI reproduces DNA Chisel CodonOptimize(method='use_best_codon') codon for codon")]
    public void OptimizeSequence_MaximizeCAI_MatchesDnaChiselUseBestCodon()
    {
        // DNA Chisel 3.2.16, DnaOptimizationProblem(constraints=[EnforceTranslation()],
        // objectives=[CodonOptimize(method="use_best_codon",
        //   codon_usage_table=python_codon_tables.get_codons_table('e_coli_316407'))]).optimize()
        // → ATGAGCAAAGGCGAAGAACTGTTTACCGGCGTGGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAAC
        const string sequence = "ATGAGCAAAGGTGAAGAACTGTTCACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTTAAC";
        const string dnaChisel = "ATGAGCAAAGGCGAAGAACTGTTTACCGGCGTGGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAAC";

        var result = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo(dnaChisel.Replace('T', 'U')));
            Assert.That(result.OptimizedCAI, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(TranslateSequence(result.OptimizedSequence),
                Is.EqualTo(TranslateSequence(result.OriginalSequence)));
        });
    }

    [Test]
    [Description("HarmonizeExpression is deterministic and allocates codons by largest remainder (DNA Chisel match_codon_usage optimum)")]
    public void OptimizeSequence_HarmonizeExpression_MatchesTargetUsageDeterministically()
    {
        // 10 Leu codons, E. coli Leu frequencies UUA .13 UUG .13 CUU .10 CUC .10 CUA .04 CUG .50.
        // Largest-remainder rounding of 10×f: floors 1,1,1,1,0,5 (sum 9); the single remaining
        // seat goes to the largest fraction (CUA, 0.4) → UUA1 UUG1 CUU1 CUC1 CUA1 CUG5.
        // The five positions that already hold CUG keep it, the rest take the other codons in
        // NCBI order. DNA Chisel's own match_codon_usage score for this allocation on the
        // 69-nt test gene equals the score its randomised optimizer reaches (−11.77).
        const string sequence = "CUGCUGCUGCUGCUGCUGCUGCUGCUGCUG";

        var first = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.HarmonizeExpression);
        var second = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.HarmonizeExpression);

        Assert.Multiple(() =>
        {
            Assert.That(first.OptimizedSequence, Is.EqualTo("CUGCUGCUGCUGCUGUUAUUGCUUCUCCUA"));
            Assert.That(second.OptimizedSequence, Is.EqualTo(first.OptimizedSequence),
                "the harmonization must be deterministic (it used weighted-random sampling before 2026-09)");
            Assert.That(TranslateSequence(first.OptimizedSequence), Is.EqualTo("LLLLLLLLLL"));
        });
    }

    [Test]
    [Description("BalancedOptimization only makes swaps that move GC toward the target window")]
    public void OptimizeSequence_BalancedOptimization_MakesNoNeutralSwaps()
    {
        // AUGGCCGCC → phase 1 AUGGCGGCG (GC 7/9 = 0.778 > 0.60). Ala codons below GCG: GCU and
        // GCA both have 2 GC, GCC has 3. One swap gives 6/9 (still high), the second 5/9 = 0.556
        // (inside [0.40, 0.60]) and the pass stops — GCA before GCU because it is the more
        // frequent E. coli codon (0.21 vs 0.16).
        var result = CodonOptimizer.OptimizeSequence(
            "AUGGCCGCC", CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGGCAGCA"));
            Assert.That(result.GcContentOptimized, Is.EqualTo(5.0 / 9.0).Within(1e-12));
            Assert.That(result.ChangedCodons, Is.EqualTo(2));
        });
    }

    [Test]
    [Description("AvoidRareCodeons falls back to the best synonymous codon when no synonym reaches the threshold")]
    public void OptimizeSequence_AvoidRareCodons_NoSynonymAboveThreshold_UsesBestCodon()
    {
        // Threshold 0.9: every Leu codon is "rare" (max is CUG at 0.50). The rare codon must
        // still be replaced by the most frequent synonym (DNA Chisel use_best_codon), not left
        // in place as before 2026-09.
        var result = CodonOptimizer.OptimizeSequence(
            "CUA", CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons,
            rareCodonThreshold: 0.9);

        Assert.That(result.OptimizedSequence, Is.EqualTo("CUG"));
    }

    [Test]
    [Description("Ambiguous and invalid triplets: IUPAC codons resolve through GeneticCode, unusable triplets are preserved")]
    [TestCase("AUGGCNUAA", "AUGGCGUAA", "MA*")]   // GCN is unambiguously Ala → best Ala codon
    [TestCase("AUGNNNGCU", "AUGNNNGCG", "MXA")]   // NNN has no unique amino acid → untouched
    [TestCase("AUG#!?GCU", "AUG#!?GCG", "MXA")]   // not a codon at all → untouched, frame kept
    public void OptimizeSequence_AmbiguousCodons_HandledPerGeneticCode(
        string sequence, string expected, string expectedProtein)
    {
        var result = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo(expected));
            Assert.That(result.ProteinSequence, Is.EqualTo(expectedProtein));
        });
    }

    [Test]
    [Description("Frequency ties are broken deterministically by NCBI codon order")]
    public void OptimizeSequence_TiedBestCodons_PicksFirstInNcbiOrder()
    {
        // Human Arg (Kazusa 9606) is the only tie in the three presets: AGA and AGG both 0.21.
        // NCBI table-1 codon order for Arg is CGU CGC CGA CGG AGA AGG, so AGA wins — the same
        // situation where Biopython CodonAdaptationIndex.optimize only warns and picks one.
        var result = CodonOptimizer.OptimizeSequence(
            "CGA", CodonOptimizer.Human, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        Assert.That(result.OptimizedSequence, Is.EqualTo("AGA"));
    }

    #endregion

    #region Helper Methods

    // Canonical translation (no private genetic-code copy in the test either).
    private static string TranslateSequence(string rnaSequence) =>
        string.Concat(Enumerable.Range(0, rnaSequence.Length / 3)
            .Select(i => GeneticCode.Standard.Translate(rnaSequence.Substring(i * 3, 3))));

    #endregion
}
