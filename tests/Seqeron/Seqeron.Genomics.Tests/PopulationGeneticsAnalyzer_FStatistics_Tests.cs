using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for F-Statistics (Fst, Fis, Fit) in PopulationGeneticsAnalyzer.
/// Test Unit: POP-FST-001
/// 
/// Evidence Sources:
/// - Wikipedia: Fixation index (https://en.wikipedia.org/wiki/Fixation_index)
/// - Wikipedia: F-statistics (https://en.wikipedia.org/wiki/F-statistics)
/// - Weir & Cockerham (1984) Evolution 38:1358-1370
/// - Wright (1965) Evolution 19:395-420
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
    /// Different populations should yield positive Fst.
    /// Source: Mathematical definition - variance > 0 for different frequencies.
    /// </summary>
    [Test]
    public void CalculateFst_DifferentPopulations_ReturnsPositive()
    {
        var pop1 = new List<(double, int)> { (0.9, 100), (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.1, 100), (0.2, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.GreaterThan(0));
    }

    /// <summary>
    /// Fixed differences (p1=1.0, p2=0.0) should yield high Fst.
    /// Source: Wikipedia - "A value of one implies... complete differentiation"
    /// Note: Actual value depends on estimator; we expect > 0.5.
    /// </summary>
    [Test]
    public void CalculateFst_FixedDifferences_ReturnsHighFst()
    {
        var pop1 = new List<(double, int)> { (1.0, 100) };
        var pop2 = new List<(double, int)> { (0.0, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.That(fst, Is.GreaterThan(0.5));
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
    /// Single locus should produce valid Fst calculation.
    /// Source: Mathematical - algorithm works for n >= 1 loci.
    /// </summary>
    [Test]
    public void CalculateFst_SingleLocus_ValidResult()
    {
        var pop1 = new List<(double, int)> { (0.8, 100) };
        var pop2 = new List<(double, int)> { (0.2, 100) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        Assert.Multiple(() =>
        {
            Assert.That(fst, Is.GreaterThan(0), "Different frequencies should yield positive Fst");
            Assert.That(fst, Is.LessThanOrEqualTo(1), "Fst must be <= 1");
        });
    }

    /// <summary>
    /// Unequal sample sizes should be weighted appropriately.
    /// Source: Weir & Cockerham (1984) - sample size weighting.
    /// </summary>
    [Test]
    public void CalculateFst_UnequalSampleSizes_WeightedCalculation()
    {
        // Large sample with freq 0.5, small sample with freq 0.9
        var pop1 = new List<(double, int)> { (0.5, 1000) };
        var pop2 = new List<(double, int)> { (0.9, 10) };

        double fst = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);

        // Should be valid and positive
        Assert.Multiple(() =>
        {
            Assert.That(fst, Is.GreaterThan(0));
            Assert.That(fst, Is.LessThanOrEqualTo(1));
        });
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
    /// CalculateFStatistics should return all three components (Fis, Fit, Fst).
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
            // Fst should be in valid range; Fis/Fit can be negative
            Assert.That(stats.Fst, Is.GreaterThanOrEqualTo(0));
            Assert.That(double.IsFinite(stats.Fis), Is.True, "Fis should be finite");
            Assert.That(double.IsFinite(stats.Fit), Is.True, "Fit should be finite");
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
    /// Verify the F-statistics partition relationship: (1-Fit) ≈ (1-Fis)(1-Fst).
    /// Source: Wright (1965) - partition formula.
    /// Note: Relationship may not be exact due to estimation method.
    /// </summary>
    [Test]
    public void CalculateFStatistics_PartitionRelationship()
    {
        var data = new List<(int, int, int, int, double, double)>
        {
            (20, 100, 30, 100, 0.4, 0.6),
            (35, 100, 25, 100, 0.5, 0.4),
            (15, 100, 40, 100, 0.3, 0.7)
        };

        var stats = PopulationGeneticsAnalyzer.CalculateFStatistics("Pop1", "Pop2", data);

        // (1 - Fit) = (1 - Fis)(1 - Fst)
        double lhs = 1 - stats.Fit;
        double rhs = (1 - stats.Fis) * (1 - stats.Fst);

        // Allow some tolerance due to estimation variance
        Assert.That(lhs, Is.EqualTo(rhs).Within(0.1),
            "F-statistics should approximately satisfy partition formula");
    }

    #endregion
}
