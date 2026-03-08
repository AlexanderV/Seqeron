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
    /// M6: Perfect LD (identical genotype vectors) should return r² = 1.0 and D' = 1.0.
    /// Source: Wikipedia (LD as a correlation) — squared Pearson correlation of identical
    /// vectors is exactly 1.0. D' is clamped to 1.0 per Lewontin (1964).
    /// Math: X₁=X₂=[0,0,1,1,2,2], mean=1, Cov=Var=2/3 → r²=1.0.
    ///       D=Cov/2=1/3, p₁=p₂=0.5, D_max=0.25, |D|/D_max=4/3 → clamped to 1.0.
    /// </summary>
    [Test]
    public void CalculateLD_PerfectLD_ReturnsExactValues()
    {
        // Perfect LD: genotypes always match → Pearson correlation = 1.0 → r² = 1.0
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 0), (1, 1), (1, 1), (2, 2), (2, 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(1.0).Within(1e-10),
                "r² must be 1.0 for identical genotype vectors (perfect LD)");
            Assert.That(ld.DPrime, Is.EqualTo(1.0).Within(1e-10),
                "D' must be 1.0 for perfect LD (Lewontin 1964, clamped)");
        });
    }

    /// <summary>
    /// M7: No LD (independent alleles) should return r² = 0 and D' = 0.
    /// Source: Wikipedia (LD definition) — when alleles are independent, D = 0 and r² = 0.
    /// Data: 3×3 balanced design where every X₂ value appears equally within each X₁ level,
    /// giving Cov(X₁,X₂) = 0 exactly → D = 0 → D' = 0.
    /// </summary>
    [Test]
    public void CalculateLD_NoLD_ReturnsZeroValues()
    {
        // Balanced 3×3 design: each (X₁, X₂) pair appears once → Cov = 0, r² = 0
        var genotypes = new List<(int, int)>
        {
            (0, 0), (0, 1), (0, 2),
            (1, 0), (1, 1), (1, 2),
            (2, 0), (2, 1), (2, 2)
        };

        var ld = PopulationGeneticsAnalyzer.CalculateLD("V1", "V2", genotypes, 1000);

        Assert.Multiple(() =>
        {
            Assert.That(ld.RSquared, Is.EqualTo(0.0).Within(1e-10),
                "r² must be 0 for perfectly balanced (independent) genotypes");
            Assert.That(ld.DPrime, Is.EqualTo(0.0).Within(1e-10),
                "D' must be 0 when Cov = 0 (no linkage disequilibrium)");
        });
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
    /// Source: Gabriel (2002) — consecutive variants with r² ≥ threshold form a block.
    /// Identical genotype vectors have r² = 1.0 ≥ 0.7.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_HighLD_CreatesBlock()
    {
        // Identical genotypes → r² = 1.0 between all adjacent pairs
        var genotypes = new List<int> { 0, 0, 1, 1, 2, 2 };
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genotypes),
            ("V2", 200, genotypes),
            ("V3", 300, genotypes)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

        Assert.That(blocks, Has.Count.EqualTo(1), "Identical genotypes (r²=1.0) must form exactly one block");
    }

    /// <summary>
    /// M12: Variants with low LD should not form a block.
    /// Source: Gabriel (2002) — no block when r² < threshold.
    /// Data: balanced design gives r² = 0 exactly, which is < 0.7.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_LowLD_ReturnsNoBlocks()
    {
        // Balanced design: every genotype value of V2 appears equally within each V1 level
        // → Cov = 0, r² = 0 < 0.7 → no block
        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, new List<int> { 0, 0, 0, 1, 1, 1, 2, 2, 2 }),
            ("V2", 200, new List<int> { 0, 1, 2, 0, 1, 2, 0, 1, 2 })
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

        Assert.That(blocks, Is.Empty, "Uncorrelated variants (r²=0) must not form a block");
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

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

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

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

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

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

        // Identical genotypes → r² = 1.0 between all adjacent pairs → one block
        Assert.That(blocks, Has.Count.EqualTo(1),
            "Identical genotypes should produce r²=1.0 and form one block");
        Assert.That(blocks[0].Start, Is.EqualTo(100), "Block should start at lowest position");
        Assert.That(blocks[0].End, Is.EqualTo(300), "Block should end at highest position");
        Assert.That(blocks[0].Variants, Has.Count.EqualTo(3), "Block should contain all 3 variants");
    }

    /// <summary>
    /// S6: Multiple separate blocks should be detected when LD breaks.
    /// Source: Gabriel (2002) — consecutive high-LD regions separated by low-LD pair.
    /// V1-V2 have r²=1.0; V2→V3 breaks LD; V4-V5 have r²=1.0.
    /// </summary>
    [Test]
    public void FindHaplotypeBlocks_MultipleBlocks_DetectsAll()
    {
        var genoA = new List<int> { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        // Balanced design against genoA → r² = 0 between V2 and V3
        var genoB = new List<int> { 0, 1, 2, 0, 1, 2, 0, 1, 2 };

        var variants = new List<(string, int, IReadOnlyList<int>)>
        {
            ("V1", 100, genoA),
            ("V2", 200, genoA),    // Block 1: V1-V2 (r²=1.0)
            ("V3", 300, genoB),    // LD break (r²=0 with V2)
            ("V4", 400, genoA),
            ("V5", 500, genoA)     // Block 2: V4-V5 (r²=1.0)
        };

        var blocks = PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold: 0.7).ToList();

        Assert.That(blocks, Has.Count.EqualTo(2), "Should detect two separate blocks");
        Assert.That(blocks[0].Variants, Has.Count.EqualTo(2), "First block: V1-V2");
        Assert.That(blocks[1].Variants, Has.Count.EqualTo(2), "Second block: V4-V5");
    }

    #endregion


}
