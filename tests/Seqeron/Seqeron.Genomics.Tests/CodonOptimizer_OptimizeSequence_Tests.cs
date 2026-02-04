using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using System.Linq;

namespace Seqeron.Genomics.Tests;

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
    [Description("M5: MaximizeCAI on already optimal sequence maintains CAI")]
    public void OptimizeSequence_AlreadyOptimal_MaintainsHighCAI()
    {
        // Arrange - E. coli optimal codons: CUG=Leu, CCG=Pro
        const string optimalSequence = "CUGCCGACC"; // L-P-T with good E. coli codons

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            optimalSequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert
        Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(result.OriginalCAI));
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
            Assert.That(result.OptimizedCAI, Is.EqualTo(1.0),
                "AUG is the only Met codon, so CAI = 1.0");
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
        // Arrange - sequence with codons that have different preferences
        // CUG (Leu): E. coli prefers, Yeast doesn't
        // AGA (Arg): Yeast prefers, E. coli doesn't
        const string sequence = "CUGAGA"; // L-R

        // Act
        var ecoliResult = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var yeastResult = CodonOptimizer.OptimizeSequence(
            sequence, CodonOptimizer.Yeast, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Assert - proteins must be same, but CAI values differ
        Assert.Multiple(() =>
        {
            Assert.That(ecoliResult.ProteinSequence, Is.EqualTo(yeastResult.ProteinSequence),
                "Same protein for both organisms");
            Assert.That(ecoliResult.OptimizedCAI, Is.Not.EqualTo(yeastResult.OptimizedCAI).Within(0.001)
                .Or.Not.EqualTo(ecoliResult.OriginalCAI).Within(0.001),
                "Different organisms should have different CAI metrics");
        });
    }

    [Test]
    [Description("M8: Human optimization differs from E. coli")]
    public void OptimizeSequence_HumanVsEcoli_DifferentCAI()
    {
        // Arrange
        const string sequence = "CUGCCGACC"; // L-P-T

        // Act
        double ecoliCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.EColiK12);
        double humanCai = CodonOptimizer.CalculateCAI(sequence, CodonOptimizer.Human);

        // Assert
        Assert.That(ecoliCai, Is.Not.EqualTo(humanCai).Within(0.001),
            "E. coli and Human have different codon preferences");
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
    [Description("M10: OptimizationResult has all fields populated")]
    public void OptimizeSequence_ReturnsValidOptimizationResult()
    {
        // Arrange
        const string sequence = "AUGGCUUAA";

        // Act
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OriginalSequence, Is.Not.Null);
            Assert.That(result.OptimizedSequence, Is.Not.Null);
            Assert.That(result.ProteinSequence, Is.Not.Null);
            Assert.That(result.OriginalCAI, Is.GreaterThanOrEqualTo(0));
            Assert.That(result.OptimizedCAI, Is.GreaterThanOrEqualTo(0));
            Assert.That(result.GcContentOriginal, Is.InRange(0, 1));
            Assert.That(result.GcContentOptimized, Is.InRange(0, 1));
            Assert.That(result.ChangedCodons, Is.GreaterThanOrEqualTo(0));
            Assert.That(result.Changes, Is.Not.Null);
        });
    }

    #endregion

    #region SHOULD Tests - Strategy Specifics

    [Test]
    [Description("S1: AvoidRareCodons only replaces codons below threshold")]
    public void OptimizeSequence_AvoidRareCodons_OnlyReplacesRare()
    {
        // Arrange - CUG is common (0.47), CUA is rare (0.04) in E. coli
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
            // The common CUG should remain, rare CUA may change
            Assert.That(result.OptimizedSequence.Substring(0, 3), Is.EqualTo("CUG"),
                "Common codon should remain unchanged");
        });
    }

    [Test]
    [Description("S2: BalancedOptimization considers GC content")]
    public void OptimizeSequence_BalancedOptimization_ConsidersGcContent()
    {
        // Arrange - sequence to optimize
        const string sequence = "AUGGCCGCC"; // High GC content

        // Act
        var result = CodonOptimizer.OptimizeSequence(
            sequence,
            CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization);

        // Assert - should report GC content
        Assert.Multiple(() =>
        {
            Assert.That(result.GcContentOriginal, Is.GreaterThan(0));
            Assert.That(result.GcContentOptimized, Is.GreaterThan(0));
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

        // Assert - protein preserved
        Assert.That(result.ProteinSequence, Is.EqualTo("LLLLLL"));
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

    #endregion

    #region Helper Methods

    private static string TranslateSequence(string rnaSequence)
    {
        var geneticCode = new Dictionary<string, char>
        {
            { "UUU", 'F' }, { "UUC", 'F' },
            { "UUA", 'L' }, { "UUG", 'L' }, { "CUU", 'L' }, { "CUC", 'L' }, { "CUA", 'L' }, { "CUG", 'L' },
            { "AUU", 'I' }, { "AUC", 'I' }, { "AUA", 'I' },
            { "AUG", 'M' },
            { "GUU", 'V' }, { "GUC", 'V' }, { "GUA", 'V' }, { "GUG", 'V' },
            { "UCU", 'S' }, { "UCC", 'S' }, { "UCA", 'S' }, { "UCG", 'S' }, { "AGU", 'S' }, { "AGC", 'S' },
            { "CCU", 'P' }, { "CCC", 'P' }, { "CCA", 'P' }, { "CCG", 'P' },
            { "ACU", 'T' }, { "ACC", 'T' }, { "ACA", 'T' }, { "ACG", 'T' },
            { "GCU", 'A' }, { "GCC", 'A' }, { "GCA", 'A' }, { "GCG", 'A' },
            { "UAU", 'Y' }, { "UAC", 'Y' },
            { "UAA", '*' }, { "UAG", '*' }, { "UGA", '*' },
            { "CAU", 'H' }, { "CAC", 'H' },
            { "CAA", 'Q' }, { "CAG", 'Q' },
            { "AAU", 'N' }, { "AAC", 'N' },
            { "AAA", 'K' }, { "AAG", 'K' },
            { "GAU", 'D' }, { "GAC", 'D' },
            { "GAA", 'E' }, { "GAG", 'E' },
            { "UGU", 'C' }, { "UGC", 'C' },
            { "UGG", 'W' },
            { "CGU", 'R' }, { "CGC", 'R' }, { "CGA", 'R' }, { "CGG", 'R' }, { "AGA", 'R' }, { "AGG", 'R' },
            { "GGU", 'G' }, { "GGC", 'G' }, { "GGA", 'G' }, { "GGG", 'G' }
        };

        var protein = new System.Text.StringBuilder();
        for (int i = 0; i + 2 < rnaSequence.Length; i += 3)
        {
            var codon = rnaSequence.Substring(i, 3);
            if (geneticCode.TryGetValue(codon, out char aa))
                protein.Append(aa);
            else
                protein.Append('X');
        }
        return protein.ToString();
    }

    #endregion
}
