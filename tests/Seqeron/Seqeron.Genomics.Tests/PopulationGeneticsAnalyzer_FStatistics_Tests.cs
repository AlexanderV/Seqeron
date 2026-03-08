using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for F-Statistics (Fst, Fis, Fit) in PopulationGeneticsAnalyzer.
/// Test Unit: POP-FST-001
/// 
/// Implementation uses Wright's variance-based Fst: Fst = σ²_S / p̄(1-p̄)
/// Source: Wright (1965) Evolution 19:395-420
/// 
/// Evidence Sources:
/// - Wikipedia: Fixation index (https://en.wikipedia.org/wiki/Fixation_index)
/// - Wikipedia: F-statistics (https://en.wikipedia.org/wiki/F-statistics)
/// - Wright (1965) Evolution 19:395-420
/// - Cavalli-Sforza et al. (1994) History and Geography of Human Genes
/// </summary>
[TestFixture]
public class PopulationGeneticsAnalyzer_FStatistics_Tests
{
    #region CalculateFst Tests

    /// <summary>
    /// Identical populations should have Fst = 0 (complete panmixia).
    /// Source: Wikipedia - "A zero value implies complete panmixia"
    /// </summary>
    [Test]
    public void CalculateFst_IdenticalPopulations_ReturnsZero()
    {
        var pop1 = new List<(double, int)> { (0.5, 100), (0.3, 100) };
        var pop2 = new List<(double, int)> { (0.5, 100), (0.3, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(0).Within(0.001));
    }

    /// <summary>
    /// Fixed differences (p1=1.0, p2=0.0) must yield Fst = 1.0 (complete differentiation).
    /// Source: Wikipedia Fixation_index - "A value of one implies... complete differentiation"
    /// 
    /// Mathematical proof for Wright's variance Fst with equal sample sizes:
    ///   pBar = 0.5, variance = 0.25, het = 0.25, Fst = 0.25/0.25 = 1.0
    /// </summary>
    [Test]
    public void CalculateFst_FixedDifferences_ReturnsOne()
    {
        var pop1 = new List<(double, int)> { (1.0, 100) };
        var pop2 = new List<(double, int)> { (0.0, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(1.0).Within(1e-10));
    }

    /// <summary>
    /// Fst must always be in range [0, 1].
    /// Source: Wikipedia - "values range from 0 to 1"
    /// </summary>
    [Test]
    public void CalculateFst_ValueRange_BetweenZeroAndOne()
    {
        var testCases = new[]
        {
            (new List<(double, int)> { (0.5, 100) }, new List<(double, int)> { (0.5, 100) }),
            (new List<(double, int)> { (0.9, 100) }, new List<(double, int)> { (0.1, 100) }),
            (new List<(double, int)> { (1.0, 100) }, new List<(double, int)> { (0.0, 100) }),
            (new List<(double, int)> { (0.7, 50), (0.3, 50) }, new List<(double, int)> { (0.4, 80), (0.6, 80) }),
        };

        foreach (var (pop1, pop2) in testCases)
        {
            double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

            Assert.Multiple(() =>
            {
                Assert.That(fst, Is.GreaterThanOrEqualTo(0), "Fst must be >= 0");
                Assert.That(fst, Is.LessThanOrEqualTo(1), "Fst must be <= 1");
            });
        }
    }

    /// <summary>
    /// Empty populations should return Fst = 0 (graceful handling).
    /// Source: Implementation contract - undefined case returns 0.
    /// </summary>
    [Test]
    public void CalculateFst_EmptyPopulations_ReturnsZero()
    {
        double fst = PopulationGeneticsAnalyzer.CalculateFst(
            new List<(double, int)>(),
            new List<(double, int)>());

        Assert.That(fst, Is.EqualTo(0));
    }

    /// <summary>
    /// Single locus with hand-calculated exact Fst.
    /// pop1=(0.8, 100), pop2=(0.2, 100):
    ///   pBar = 0.5, var = (100*0.09 + 100*0.09)/200 = 0.09, het = 0.25
    ///   Fst = 0.09/0.25 = 0.36
    /// </summary>
    [Test]
    public void CalculateFst_SingleLocus_ExactValue()
    {
        var pop1 = new List<(double, int)> { (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.2, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(0.36).Within(1e-10));
    }

    /// <summary>
    /// Multi-locus hand-calculated Fst.
    /// pop1=[(0.9,100),(0.8,100)], pop2=[(0.1,100),(0.2,100)]:
    ///   Locus 1: pBar=0.5, var=0.16, het=0.25
    ///   Locus 2: pBar=0.5, var=0.09, het=0.25
    ///   Fst = (0.16+0.09)/(0.25+0.25) = 0.50
    /// </summary>
    [Test]
    public void CalculateFst_MultiLocus_ExactValue()
    {
        var pop1 = new List<(double, int)> { (0.9, 100), (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.1, 100), (0.2, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(0.50).Within(1e-10));
    }

    /// <summary>
    /// Unequal sample sizes: larger populations contribute more to the weighted mean.
    /// Source: Wright (1965) - variance weighted by subpopulation sizes.
    ///
    /// pop1=(0.5, 1000), pop2=(0.9, 10):
    ///   pBar = 509/1010, var ≈ 0.001568, het ≈ 0.249984
    ///   Fst ≈ 0.006274
    ///
    /// With equal sizes (n1=n2=100, same freqs): pBar=0.7, Fst = 4/21 ≈ 0.1905
    /// The unequal-size Fst is much smaller because the large pop (p=0.5) dominates,
    /// pulling pBar close to 0.5 and reducing weighted variance.
    /// </summary>
    [Test]
    public void CalculateFst_UnequalSampleSizes_WeightedCalculation()
    {
        var pop1 = new List<(double, int)> { (0.5, 1000) };
        var pop2 = new List<(double, int)> { (0.9, 10) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        // Hand-calculated: pBar=509/1010, Fst ≈ 0.006274
        Assert.That(fst, Is.EqualTo(0.006274288358450095).Within(1e-10));

        // Compare with equal-size case to demonstrate weighting effect
        var pop1Equal = new List<(double, int)> { (0.5, 100) };
        var pop2Equal = new List<(double, int)> { (0.9, 100) };
        double fstEqual = PopulationGeneticsAnalyzer.CalculateFst(pop1Equal, pop2Equal);

        // Equal sizes: pBar=0.7, Fst = 4/21 ≈ 0.1905 — much larger
        Assert.That(fstEqual, Is.EqualTo(4.0 / 21.0).Within(1e-10));
        Assert.That(fst, Is.LessThan(fstEqual),
            "Larger population weight should pull Fst toward 0 when large pop has intermediate freq");
    }

    #endregion

    #region CalculatePairwiseFst Tests

    /// <summary>
    /// Pairwise Fst matrix should have correct dimensions.
    /// </summary>
    [Test]
    public void CalculatePairwiseFst_ThreePopulations_CorrectDimensions()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100) }),
            ("Pop2", new List<(double, int)> { (0.6, 100) }),
            ("Pop3", new List<(double, int)> { (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        Assert.Multiple(() =>
        {
            Assert.That(matrix.GetLength(0), Is.EqualTo(3));
            Assert.That(matrix.GetLength(1), Is.EqualTo(3));
        });
    }

    /// <summary>
    /// Diagonal of pairwise Fst matrix must be zero (self-comparison).
    /// Source: Mathematical property - Fst(pop, pop) = 0.
    /// </summary>
    [Test]
    public void CalculatePairwiseFst_DiagonalIsZero()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100) }),
            ("Pop2", new List<(double, int)> { (0.6, 100) }),
            ("Pop3", new List<(double, int)> { (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        Assert.Multiple(() =>
        {
            Assert.That(matrix[0, 0], Is.EqualTo(0), "Diagonal[0,0] must be 0");
            Assert.That(matrix[1, 1], Is.EqualTo(0), "Diagonal[1,1] must be 0");
            Assert.That(matrix[2, 2], Is.EqualTo(0), "Diagonal[2,2] must be 0");
        });
    }

    /// <summary>
    /// Pairwise Fst matrix must be symmetric: Fst(i,j) = Fst(j,i).
    /// Source: Mathematical property.
    /// </summary>
    [Test]
    public void CalculatePairwiseFst_SymmetricMatrix()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100) }),
            ("Pop2", new List<(double, int)> { (0.6, 100) }),
            ("Pop3", new List<(double, int)> { (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        Assert.Multiple(() =>
        {
            Assert.That(matrix[0, 1], Is.EqualTo(matrix[1, 0]).Within(1e-10), "Fst[0,1] == Fst[1,0]");
            Assert.That(matrix[0, 2], Is.EqualTo(matrix[2, 0]).Within(1e-10), "Fst[0,2] == Fst[2,0]");
            Assert.That(matrix[1, 2], Is.EqualTo(matrix[2, 1]).Within(1e-10), "Fst[1,2] == Fst[2,1]");
        });
    }

    /// <summary>
    /// All pairwise Fst values must be in valid range [0, 1].
    /// </summary>
    [Test]
    public void CalculatePairwiseFst_AllValuesInRange()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100), (0.3, 100) }),
            ("Pop2", new List<(double, int)> { (0.7, 100), (0.5, 100) }),
            ("Pop3", new List<(double, int)> { (0.1, 100), (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                Assert.That(matrix[i, j], Is.InRange(0.0, 1.0),
                    $"Fst[{i},{j}] must be in [0,1]");
            }
        }
    }

    #endregion

    #region CalculateFStatistics Tests

    /// <summary>
    /// CalculateFStatistics should return all three components with correct values.
    ///
    /// Data: [(20,50,25,50,0.4,0.5), (30,50,15,50,0.5,0.3)]
    ///
    /// Locus 1: pBar=0.45, hetObs=45, hetExp=49, hetTotal=49.5
    /// Locus 2: pBar=0.4,  hetObs=45, hetExp=46, hetTotal=48
    /// hi=90/200=0.45,  hs=95/200=0.475,  ht=97.5/200=0.4875
    /// Fis = 1 - hi/hs = 1/19,  Fit = 1 - hi/ht = 1/13,  Fst = 1 - hs/ht = 1/39
    /// </summary>
    [Test]
    public void CalculateFStatistics_ReturnsAllComponents()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 50, 25, 50, 0.4, 0.5),
            (30, 50, 15, 50, 0.5, 0.3)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Population1, Is.EqualTo("Pop1"));
            Assert.That(stats.Population2, Is.EqualTo("Pop2"));
            Assert.That(stats.Fis, Is.EqualTo(1.0 / 19.0).Within(1e-10), "Fis = 1/19");
            Assert.That(stats.Fit, Is.EqualTo(1.0 / 13.0).Within(1e-10), "Fit = 1/13");
            Assert.That(stats.Fst, Is.EqualTo(1.0 / 39.0).Within(1e-10), "Fst = 1/39");
        });
    }

    /// <summary>
    /// Population names should be preserved in the result.
    /// </summary>
    [Test]
    public void CalculateFStatistics_PopulationNamesPreserved()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (10, 50, 10, 50, 0.3, 0.3)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("European", "Asian", data);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Population1, Is.EqualTo("European"));
            Assert.That(stats.Population2, Is.EqualTo("Asian"));
        });
    }

    /// <summary>
    /// Empty data should return zero F-statistics values.
    /// </summary>
    [Test]
    public void CalculateFStatistics_EmptyData_ReturnsZeroValues()
    {
        var data = new List<(int, int, int, int, double, double)>();

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Fst, Is.EqualTo(0));
            Assert.That(stats.Fis, Is.EqualTo(0));
            Assert.That(stats.Fit, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// Fst should be in range [0, 1].
    /// Fis and Fit can be negative (excess heterozygosity) but typically in [-1, 1].
    /// Source: Wikipedia - Fst range [0,1], Fis/Fit can be negative.
    /// </summary>
    [Test]
    public void CalculateFStatistics_ComponentsInValidRange()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 50, 25, 50, 0.4, 0.5),
            (30, 50, 15, 50, 0.5, 0.3),
            (25, 50, 20, 50, 0.6, 0.4)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Fst, Is.InRange(0.0, 1.0), "Fst must be in [0,1]");
            Assert.That(stats.Fis, Is.InRange(-1.0, 1.0), "Fis must be in [-1,1]");
            Assert.That(stats.Fit, Is.InRange(-1.0, 1.0), "Fit must be in [-1,1]");
        });
    }

    /// <summary>
    /// Verify Wright's partition identity: (1-Fit) = (1-Fis)(1-Fst).
    /// Source: Wikipedia F-statistics §Partition - "(1-F_IT) = (1-F_IS)(1-F_ST)"
    /// 
    /// This is an algebraic identity (not an approximation) for the heterozygosity-based
    /// computation used here: (H_I/H_S)(H_S/H_T) = H_I/H_T.
    /// </summary>
    [Test]
    public void CalculateFStatistics_PartitionRelationship_ExactIdentity()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 100, 30, 100, 0.4, 0.6),
            (35, 100, 25, 100, 0.5, 0.4),
            (15, 100, 40, 100, 0.3, 0.7)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        // (1 - Fit) = (1 - Fis)(1 - Fst)  — exact algebraic identity
        double lhs = 1 - stats.Fit;
        double rhs = (1 - stats.Fis) * (1 - stats.Fst);

        Assert.That(lhs, Is.EqualTo(rhs).Within(1e-10),
            "F-statistics partition formula must hold exactly for heterozygosity-based computation");
    }

    /// <summary>
    /// Hand-calculated F-statistics with exact expected values.
    /// 
    /// Data: [(20, 100, 30, 100, 0.4, 0.6), (35, 100, 25, 100, 0.5, 0.4), (15, 100, 40, 100, 0.3, 0.7)]
    /// 
    /// Locus 1: pBar=0.5, hetObs=50, hetExp=96, hetTotal=100
    /// Locus 2: pBar=0.45, hetObs=60, hetExp=98, hetTotal=99
    /// Locus 3: pBar=0.5, hetObs=55, hetExp=84, hetTotal=100
    /// 
    /// hi = 165/600 = 0.275,  hs = 278/600 ≈ 0.46333,  ht = 299/600 ≈ 0.49833
    /// Fis = 1 - hi/hs ≈ 0.40647
    /// Fit = 1 - hi/ht ≈ 0.44832
    /// Fst = 1 - hs/ht ≈ 0.07020
    /// 
    /// Source: Heterozygosity-based definitions from Wikipedia F-statistics §Definitions.
    /// </summary>
    [Test]
    public void CalculateFStatistics_HandCalculated_ExactValues()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 100, 30, 100, 0.4, 0.6),
            (35, 100, 25, 100, 0.5, 0.4),
            (15, 100, 40, 100, 0.3, 0.7)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        // Hand-calculated: hi=165/600, hs=278/600, ht=299/600
        double expectedFis = 1.0 - (165.0 / 278.0); // 1 - hi/hs
        double expectedFit = 1.0 - (165.0 / 299.0); // 1 - hi/ht
        double expectedFst = 1.0 - (278.0 / 299.0); // 1 - hs/ht

        Assert.Multiple(() =>
        {
            Assert.That(stats.Fis, Is.EqualTo(expectedFis).Within(1e-10), "Fis hand-calculated");
            Assert.That(stats.Fit, Is.EqualTo(expectedFit).Within(1e-10), "Fit hand-calculated");
            Assert.That(stats.Fst, Is.EqualTo(expectedFst).Within(1e-10), "Fst hand-calculated");
        });
    }

    #endregion

    #region Reference Data Validation Tests - Published Fst Values

    /// <summary>
    /// Multi-locus Fst with moderate differentiation.
    /// Uses exact binary fractions to avoid floating-point representation errors.
    ///
    /// pop1 = [(0.25,200), (0.5,200), (0.125,200)]
    /// pop2 = [(0.5,200),  (0.625,200), (0.375,200)]
    ///
    /// Locus 1: pBar=0.375, var=0.015625, het=0.234375
    /// Locus 2: pBar=0.5625, var=0.00390625, het=0.24609375
    /// Locus 3: pBar=0.25, var=0.015625, het=0.1875
    /// sum(var)=9/256, sum(het)=171/256, Fst = 9/171 = 1/19 ≈ 0.0526
    /// Wright scale: moderate differentiation (0.05 – 0.15).
    /// </summary>
    [Test]
    public void CalculateFst_MultiLocusModerate_ExactValue()
    {
        var pop1 = new List<(double, int)> { (0.25, 200), (0.5, 200), (0.125, 200) };
        var pop2 = new List<(double, int)> { (0.5, 200), (0.625, 200), (0.375, 200) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(1.0 / 19.0).Within(1e-10), "Fst = 1/19");
    }

    /// <summary>
    /// Validates Fst with exact values at opposite ends of Wright's interpretation scale.
    /// Source: Hartl &amp; Clark (2007); Wright (1978) scale.
    ///
    /// Little: pop1=[(0.50,100),(0.48,100)], pop2=[(0.52,100),(0.50,100)]
    ///   Locus 1: pBar=0.51, var=0.0001, het=0.2499
    ///   Locus 2: pBar=0.49, var=0.0001, het=0.2499
    ///   Fst = 0.0002/0.4998 = 1/2499 ≈ 0.000400 (&lt; 0.05 = little)
    ///
    /// Very great: pop1=[(0.80,100),(0.85,100)], pop2=[(0.30,100),(0.25,100)]
    ///   Locus 1: pBar=0.55, var=0.0625, het=0.2475
    ///   Locus 2: pBar=0.55, var=0.09, het=0.2475
    ///   Fst = 0.1525/0.495 = 61/198 ≈ 0.3081 (&gt; 0.25 = very great)
    /// </summary>
    [Test]
    public void CalculateFst_WrightInterpretationScale_ExactValues()
    {
        var pop1Similar = new List<(double, int)> { (0.50, 100), (0.48, 100) };
        var pop2Similar = new List<(double, int)> { (0.52, 100), (0.50, 100) };
        double fstLittle = PopulationGeneticsAnalyzer.CalculateFst(pop1Similar, pop2Similar);

        var pop1Different = new List<(double, int)> { (0.80, 100), (0.85, 100) };
        var pop2Different = new List<(double, int)> { (0.30, 100), (0.25, 100) };
        double fstGreat = PopulationGeneticsAnalyzer.CalculateFst(pop1Different, pop2Different);

        Assert.Multiple(() =>
        {
            Assert.That(fstLittle, Is.EqualTo(1.0 / 2499.0).Within(1e-10),
                "Fst = 1/2499 (Wright: little differentiation)");
            Assert.That(fstGreat, Is.EqualTo(61.0 / 198.0).Within(1e-10),
                "Fst = 61/198 (Wright: very great differentiation)");
        });
    }

    /// <summary>
    /// Fixed alleles at multiple loci must yield Fst = 1.0 exactly.
    /// Source: Wikipedia Fixation_index - "A value of one implies complete differentiation"
    /// 
    /// For Wright's variance Fst, each locus with p1=1,p2=0 contributes:
    ///   pBar=0.5, variance=0.25, het=0.25 → ratio=1.0 per locus and in aggregate.
    /// </summary>
    [Test]
    public void CalculateFst_FixedAlleles_MultiLoci_ReturnsOne()
    {
        var pop1Fixed = new List<(double, int)>
        {
            (1.0, 100), (1.0, 100), (1.0, 100), (1.0, 100), (1.0, 100)
        };
        var pop2Fixed = new List<(double, int)>
        {
            (0.0, 100), (0.0, 100), (0.0, 100), (0.0, 100), (0.0, 100)
        };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1Fixed, pop2Fixed);

        Assert.That(fst, Is.EqualTo(1.0).Within(1e-10),
            "Completely fixed alleles must yield Fst = 1.0 (Wright's variance Fst)");
    }

    /// <summary>
    /// Validates Fst with exact values at three differentiation levels + monotonicity.
    /// Source: Wright (1951) island model; Fst increases with allele frequency divergence.
    ///
    /// Similar (p=0.48): pBar=0.49, Fst = 1/2499
    /// Moderate (p=0.35): pBar=0.425, Fst = 9/391
    /// Different (p=0.15): pBar=0.325, Fst = 49/351
    /// </summary>
    [Test]
    public void CalculateFst_IslandModelConsistency_ExactValuesAndMonotonic()
    {
        var popA = new List<(double, int)> { (0.5, 100) };

        var popB_similar = new List<(double, int)> { (0.48, 100) };
        var popB_moderate = new List<(double, int)> { (0.35, 100) };
        var popB_different = new List<(double, int)> { (0.15, 100) };

        double fst_similar = PopulationGeneticsAnalyzer.CalculateFst(popA, popB_similar);
        double fst_moderate = PopulationGeneticsAnalyzer.CalculateFst(popA, popB_moderate);
        double fst_different = PopulationGeneticsAnalyzer.CalculateFst(popA, popB_different);

        Assert.Multiple(() =>
        {
            Assert.That(fst_similar, Is.EqualTo(1.0 / 2499.0).Within(1e-10), "Fst similar = 1/2499");
            Assert.That(fst_moderate, Is.EqualTo(9.0 / 391.0).Within(1e-10), "Fst moderate = 9/391");
            Assert.That(fst_different, Is.EqualTo(49.0 / 351.0).Within(1e-10), "Fst different = 49/351");

            Assert.That(fst_similar, Is.LessThan(fst_moderate),
                "Less differentiation should yield lower Fst");
            Assert.That(fst_moderate, Is.LessThan(fst_different),
                "More differentiation should yield higher Fst");
        });
    }

    #endregion

    #region Missing Coverage Tests

    /// <summary>
    /// All monomorphic sites (no polymorphism) should give Fst = 0.
    /// Source: Evidence §5.1 — "All monomorphic sites → Fst = 0 (no variation)".
    /// Both populations share the same allele frequency at every locus,
    /// so σ²_S = 0 and Fst = 0/het = 0.
    /// </summary>
    [Test]
    public void CalculateFst_MonomorphicSites_ReturnsZero()
    {
        var pop1 = new List<(double, int)> { (0.3, 100), (0.6, 100), (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.3, 100), (0.6, 100), (0.8, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.EqualTo(0.0).Within(1e-10));
    }

    /// <summary>
    /// Both populations fixed for the same allele: pBar = 0 or 1, denominator = 0.
    /// Source: Evidence §5.2 — "Division by zero (pBar = 0 or 1) → Return 0".
    /// </summary>
    [Test]
    public void CalculateFst_BothFixedSameAllele_ReturnsZero()
    {
        // pBar = 0
        var pop1Zero = new List<(double, int)> { (0.0, 100) };
        var pop2Zero = new List<(double, int)> { (0.0, 100) };
        double fstZero = PopulationGeneticsAnalyzer.CalculateFst(pop1Zero, pop2Zero);

        // pBar = 1
        var pop1One = new List<(double, int)> { (1.0, 100) };
        var pop2One = new List<(double, int)> { (1.0, 100) };
        double fstOne = PopulationGeneticsAnalyzer.CalculateFst(pop1One, pop2One);

        Assert.Multiple(() =>
        {
            Assert.That(fstZero, Is.EqualTo(0), "Both fixed at 0: pBar=0 → denominator=0 → Fst=0");
            Assert.That(fstOne, Is.EqualTo(0), "Both fixed at 1: pBar=1 → denominator=0 → Fst=0");
        });
    }

    /// <summary>
    /// Pairwise Fst matrix with exact hand-calculated cell values.
    /// Pop1(0.5), Pop2(0.6), Pop3(0.9) — single locus, equal sizes.
    ///
    /// Fst[0,1]: pBar=0.55, var=0.0025, het=0.2475, Fst = 1/99
    /// Fst[0,2]: pBar=0.70, var=0.04, het=0.21, Fst = 4/21
    /// Fst[1,2]: pBar=0.75, var=0.0225, het=0.1875, Fst = 3/25
    /// </summary>
    [Test]
    public void CalculatePairwiseFst_ExactCellValues()
    {
        var populations = new List<(string, IReadOnlyList<(double, int)>)>
        {
            ("Pop1", new List<(double, int)> { (0.5, 100) }),
            ("Pop2", new List<(double, int)> { (0.6, 100) }),
            ("Pop3", new List<(double, int)> { (0.9, 100) })
        };

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(populations);

        Assert.Multiple(() =>
        {
            Assert.That(matrix[0, 1], Is.EqualTo(1.0 / 99.0).Within(1e-10), "Fst[Pop1,Pop2] = 1/99");
            Assert.That(matrix[0, 2], Is.EqualTo(4.0 / 21.0).Within(1e-10), "Fst[Pop1,Pop3] = 4/21");
            Assert.That(matrix[1, 2], Is.EqualTo(3.0 / 25.0).Within(1e-10), "Fst[Pop2,Pop3] = 3/25");
        });
    }

    /// <summary>
    /// F-statistics with excess heterozygosity: Fis must be negative.
    /// Source: Wikipedia — "Fis can be negative with excess heterozygosity".
    ///
    /// Pop1: p=0.3, HetObs=60/100 (expected=42 from HWE)
    /// Pop2: p=0.7, HetObs=80/100 (expected=42 from HWE)
    ///
    /// hi = 140/200 = 7/10
    /// hs = 84/200 = 21/50
    /// ht = 100/200 = 1/2
    /// Fis = 1 - (7/10)/(21/50) = -2/3
    /// Fit = 1 - (7/10)/(1/2) = -2/5
    /// Fst = 1 - (21/50)/(1/2) = 4/25
    /// Partition: (1-Fit)(=7/5) = (1-Fis)(5/3) × (1-Fst)(21/25) = 7/5 ✓
    /// </summary>
    [Test]
    public void CalculateFStatistics_ExcessHeterozygosity_NegativeFis()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (60, 100, 80, 100, 0.3, 0.7)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Fis, Is.EqualTo(-2.0 / 3.0).Within(1e-10), "Fis = -2/3 (excess het)");
            Assert.That(stats.Fit, Is.EqualTo(-2.0 / 5.0).Within(1e-10), "Fit = -2/5");
            Assert.That(stats.Fst, Is.EqualTo(4.0 / 25.0).Within(1e-10), "Fst = 4/25");

            // Verify partition identity still holds
            double lhs = 1 - stats.Fit;
            double rhs = (1 - stats.Fis) * (1 - stats.Fst);
            Assert.That(lhs, Is.EqualTo(rhs).Within(1e-10), "Partition identity (1-Fit)=(1-Fis)(1-Fst)");
        });
    }

    #endregion
}
