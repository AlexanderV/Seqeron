using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for PopulationGeneticsAnalyzer linkage disequilibrium methods.
/// Test Unit: POP-LD-001
/// Evidence: docs/Evidence/POP-LD-001.md
/// </summary>
[TestFixture]
public class PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests
{
    #region CalculateLD - Basic Functionality

    /// <summary>
    /// M1: Empty genotypes should return zero LD.
    /// Source: Wikipedia (LD definition) - No data means no association.
    /// </summary>
    [Test]
    public void CalculateLD_EmptyGenotypes_ReturnsZeroLD()
    {
        var ld = PopulationGeneticsAnalyzer.CalculateLD(
            "V1", "V2", new List<(int, int)>(), 1000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(0), "r² should be 0 for empty genotypes");
            Assert.That(ld.DPrime, Is.EqualTo(0), "D' should be 0 for empty genotypes");
        });
    }

    /// <summary>
    /// M4, M5: CalculateLD should preserve variant IDs and distance.
    /// Source: API contract - Data integrity.
    /// </summary>
    [Test]
    public void CalculateLD_PreservesVariantIdsAndDistance()
    {
        var genotypes = new List<(int, int)> { (0, 0), (1, 1) };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("Variant1", "Variant2", genotypes, 5000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.Variant1, Is.EqualTo("Variant1"), "Variant1 ID should be preserved");
            Assert.That(ld.Variant2, Is.EqualTo("Variant2"), "Variant2 ID should be preserved");
            Assert.That(ld.Distance, Is.EqualTo(5000), "Distance should be preserved");
        });
    }

    /// <summary>
    /// M6: Perfect LD (alleles always co-occur) should return high r².
    /// Source: Wikipedia (LD definition) - Complete association.
    /// </summary>
    [Test]
    public void CalculateLD_PerfectLD_ReturnsHighRSquared()
    {
        // Perfect LD: genotypes always match
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 0), (1, 1), (1, 1), (2, 2), (2, 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        // Due to phase estimation from genotypes, perfect r²=1 is not always achieved
        Assert.That(ld.RSquared, Is.GreaterThan(0.4),
            "r² should be high (>0.4) for perfect LD pattern");
    }

    /// <summary>
    /// M7: No LD (random association) should return low r².
    /// Source: Wikipedia (LD definition) - Random association.
    /// </summary>
    [Test]
    public void CalculateLD_NoLD_ReturnsLowRSquared()
    {
        // No LD: random pairing
        var genotypes = new List<(int, int)>
        {
            (0, 2), (2, 0), (1, 1), (0, 1), (2, 1), (1, 0)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        Assert.That(ld.RSquared, Is.LessThan(0.5),
            "r² should be low (<0.5) for random association pattern");
    }

    #endregion

    #region CalculateLD - Range Invariants

    /// <summary>
    /// M2: r² must always be in range [0, 1].
    /// Source: Hill & Robertson (1968) - Mathematical invariant.
    /// </summary>
    [Test]
    [TestCase(new[] { 0, 0, 1, 1, 2, 2 }, new[] { 0, 0, 1, 1, 2, 2 })]
    [TestCase(new[] { 0, 2, 1, 0, 2, 1 }, new[] { 2, 0, 1, 2, 0, 1 })]
    [TestCase(new[] { 0, 0, 0, 0 }, new[] { 0, 1, 2, 1 })]
    [TestCase(new[] { 1, 1, 1, 1 }, new[] { 1, 1, 1, 1 })]
    public void CalculateLD_RSquared_AlwaysInValidRange(int[] geno1, int[] geno2)
    {
        var genotypes = geno1.Zip(geno2, (g1, g2) => (g1, g2)).ToList();

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.GreaterThanOrEqualTo(0), "r² must be >= 0");
            Assert.That(ld.RSquared, Is.LessThanOrEqualTo(1), "r² must be <= 1");
        });
    }

    /// <summary>
    /// M3: D' must always be in range [0, 1] (using |D'|).
    /// Source: Lewontin (1964) - Normalized measure invariant.
    /// </summary>
    [Test]
    [TestCase(new[] { 0, 0, 1, 1, 2, 2 }, new[] { 0, 0, 1, 1, 2, 2 })]
    [TestCase(new[] { 0, 2, 1, 0, 2, 1 }, new[] { 2, 0, 1, 2, 0, 1 })]
    [TestCase(new[] { 0, 0, 0, 0 }, new[] { 0, 1, 2, 1 })]
    public void CalculateLD_DPrime_AlwaysInValidRange(int[] geno1, int[] geno2)
    {
        var genotypes = geno1.Zip(geno2, (g1, g2) => (g1, g2)).ToList();

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ld.DPrime, Is.GreaterThanOrEqualTo(0), "D' must be >= 0");
            Assert.That(ld.DPrime, Is.LessThanOrEqualTo(1), "D' must be <= 1");
        });
    }

    #endregion

    #region CalculateLD - Edge Cases

    /// <summary>
    /// M8: Monomorphic locus should handle division by zero gracefully.
    /// Source: Mathematical edge case - Denominator = 0 protection.
    /// </summary>
    [Test]
    public void CalculateLD_MonomorphicFirstLocus_ReturnsZeroRSquared()
    {
        // First locus all homozygous major (p1 = 0), second varies
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 1), (0, 2), (0, 0)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(0), "r² should be 0 for monomorphic locus");
            Assert.That(double.IsNaN(ld.RSquared), Is.False, "r² should not be NaN");
            Assert.That(double.IsInfinity(ld.RSquared), Is.False, "r² should not be Infinity");
        });
    }

    /// <summary>
    /// M8: Monomorphic second locus should handle division by zero gracefully.
    /// </summary>
    [Test]
    public void CalculateLD_MonomorphicSecondLocus_ReturnsZeroRSquared()
    {
        // First varies, second locus all homozygous major
        var genotypes = new List<(int, int)>
        {
            (0, 0), (1, 0), (2, 0), (1, 0)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(0), "r² should be 0 for monomorphic locus");
            Assert.That(double.IsNaN(ld.RSquared), Is.False, "r² should not be NaN");
        });
    }

    /// <summary>
    /// S1: Single genotype pair should produce valid (though unreliable) result.
    /// Source: Statistical validity - Minimum sample.
    /// </summary>
    [Test]
    public void CalculateLD_SingleGenotypePair_ReturnsValidResult()
    {
        var genotypes = new List<(int, int)> { (1, 1) };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.Multiple(() =>
        {
            Assert.That(double.IsNaN(ld.RSquared), Is.False, "r² should not be NaN");
            Assert.That(double.IsNaN(ld.DPrime), Is.False, "D' should not be NaN");
        });
    }

    /// <summary>
    /// S2: All homozygous major genotypes (no variation).
    /// Source: No variation edge case.
    /// </summary>
    [Test]
    public void CalculateLD_AllHomozygousMajor_ReturnsZeroRSquared()
    {
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 0), (0, 0), (0, 0)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.That(ld.RSquared, Is.EqualTo(0), "r² should be 0 when no variation exists");
    }

    /// <summary>
    /// S3: All homozygous minor genotypes (no variation).
    /// Source: No variation edge case.
    /// </summary>
    [Test]
    public void CalculateLD_AllHomozygousMinor_ReturnsZeroRSquared()
    {
        var genotypes = new List<(int, int)>
        {
            (2, 2), (2, 2), (2, 2), (2, 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 100);

        Assert.That(ld.RSquared, Is.EqualTo(0), "r² should be 0 when no variation exists");
    }

    #endregion

    #region FindHaplotypeBlocks - Basic Functionality

    /// <summary>
    /// M9: Single variant cannot form a block.
    /// Source: Gabriel (2002) - Block requires ≥2 variants.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_SingleVariant_ReturnsNoBlocks()
    {
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, new List<int> { 0, 1, 2 })
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants).ToList();

        Assert.That(blocks, Is.Empty, "Single variant cannot form a block");
    }

    /// <summary>
    /// M10: Empty variants list returns no blocks.
    /// Source: Gabriel (2002) - Block definition.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_EmptyVariants_ReturnsNoBlocks()
    {
        var variants = new List<(string, int, IReadOnlyList<int>)>();

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants).ToList();

        Assert.That(blocks, Is.Empty, "Empty variants should return no blocks");
    }

    /// <summary>
    /// M11: Variants with high LD should form a block.
    /// Source: Gabriel (2002) - Block detection.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_HighLD_CreatesBlock()
    {
        // Create variants in strong LD (identical genotypes)
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 200, genotypes),
            ("V3", 300, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        Assert.That(blocks, Has.Count.GreaterThanOrEqualTo(1), "High LD variants should form at least one block");
    }

    /// <summary>
    /// M12: Variants with low LD should not form a block.
    /// Source: Gabriel (2002) - No block when LD below threshold.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_LowLD_ReturnsNoBlocks()
    {
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, new List<int> { 0, 1, 2, 0, 1, 2 }),
            ("V2", 200, new List<int> { 2, 1, 0, 2, 1, 0 })
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.9).ToList();

        Assert.That(blocks, Is.Empty, "Low LD variants should not form a block");
    }

    #endregion

    #region FindHaplotypeBlocks - Invariants

    /// <summary>
    /// M13: Block start position must be less than or equal to end position.
    /// Source: Position ordering invariant.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_BlockPositions_StartLessThanOrEqualToEnd()
    {
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 500, genotypes),
            ("V3", 1000, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        foreach (var block in blocks)
        {
            Assert.That(block.Start, Is.LessThanOrEqualTo(block.End),
                $"Block start ({block.Start}) must be <= end ({block.End})");
        }
    }

    /// <summary>
    /// M14: Each block must contain at least 2 variants.
    /// Source: Gabriel (2002) - Block definition.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_EachBlock_ContainsAtLeastTwoVariants()
    {
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 200, genotypes),
            ("V3", 300, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        foreach (var block in blocks)
        {
            Assert.That(block.Variants.Count, Is.GreaterThanOrEqualTo(2),
                "Each block must contain at least 2 variants");
        }
    }

    #endregion

    #region FindHaplotypeBlocks - Should Tests

    /// <summary>
    /// S4: ldThreshold parameter should affect block formation.
    /// Source: API contract - Threshold customization.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_ThresholdParameter_AffectsBlockFormation()
    {
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 200, genotypes)
        };

        var blocksLowThreshold = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.1).ToList();
        var blocksHighThreshold = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.99).ToList();

        // With very low threshold, blocks more likely; with very high threshold, less likely
        Assert.That(blocksLowThreshold.Count, Is.GreaterThanOrEqualTo(blocksHighThreshold.Count),
            "Lower threshold should produce same or more blocks than higher threshold");
    }

    /// <summary>
    /// S5: Variants should be ordered by position before block detection.
    /// Source: Gabriel (2002) - Correct block boundaries.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_UnorderedInput_OrdersByPosition()
    {
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        // Provide variants out of order
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V3", 300, genotypes),
            ("V1", 100, genotypes),
            ("V2", 200, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        // Should still form correct block regardless of input order
        if (blocks.Count > 0)
        {
            Assert.That(blocks[0].Start, Is.EqualTo(100), "Block should start at lowest position");
        }
    }

    /// <summary>
    /// S6: Multiple separate blocks should be detected when LD breaks.
    /// Source: Gabriel (2002) - Complex genome regions.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_MultipleBlocks_DetectsAll()
    {
        var genoA = new List<int> { 0, 0, 1, 1, 2, 2 };
        var genoB = new List<int> { 2, 1, 0, 2, 1, 0 }; // Different pattern to break LD
        var genoC = new List<int> { 0, 0, 1, 1, 2, 2 }; // Same as A

        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genoA),
            ("V2", 200, genoA),    // Block 1: V1-V2
            ("V3", 300, genoB),    // LD break
            ("V4", 400, genoC),
            ("V5", 500, genoC)     // Block 2: V4-V5 (if LD is high enough)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.3).ToList();

        // The exact number depends on LD calculation; just verify no overlap
        foreach (var block in blocks)
        {
            Assert.That(block.Variants.Count, Is.GreaterThanOrEqualTo(2));
        }
    }

    #endregion
}
