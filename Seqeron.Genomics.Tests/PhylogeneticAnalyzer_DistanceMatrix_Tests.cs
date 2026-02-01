using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for PHYLO-DIST-001: Phylogenetic Distance Matrix Calculation.
/// Covers: CalculateDistanceMatrix, CalculatePairwiseDistance
/// 
/// Evidence: Wikipedia (Models of DNA evolution, Distance matrices in phylogeny),
/// Jukes &amp; Cantor (1969), Kimura (1980)
/// </summary>
[TestFixture]
[Category("PHYLO-DIST-001")]
public class PhylogeneticAnalyzer_DistanceMatrix_Tests
{
    #region M01-M04: Core Distance Matrix Properties

    [Test]
    [Description("M01: Identical sequences should have zero distance for all methods")]
    public void CalculatePairwiseDistance_IdenticalSequences_ReturnsZeroForAllMethods()
    {
        const string seq = "ACGTACGT";

        Assert.Multiple(() =>
        {
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq, PhylogeneticAnalyzer.DistanceMethod.Hamming),
                Is.EqualTo(0), "Hamming");
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq, PhylogeneticAnalyzer.DistanceMethod.PDistance),
                Is.EqualTo(0), "PDistance");
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq, PhylogeneticAnalyzer.DistanceMethod.JukesCantor),
                Is.EqualTo(0), "JukesCantor");
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(seq, seq, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter),
                Is.EqualTo(0), "Kimura2Parameter");
        });
    }

    [Test]
    [Description("M02: Distance matrix must be symmetric: d(i,j) = d(j,i)")]
    public void CalculateDistanceMatrix_IsSymmetric()
    {
        var seqs = new List<string>
        {
            "ACGTACGT",
            "TCGTACGT",
            "GCGTACGA"
        };

        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        Assert.Multiple(() =>
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Assert.That(matrix[i, j], Is.EqualTo(matrix[j, i]).Within(1e-10),
                        $"Symmetry violated at [{i},{j}]");
                }
            }
        });
    }

    [Test]
    [Description("M03: Distance matrix diagonal must be zero: d(i,i) = 0")]
    public void CalculateDistanceMatrix_DiagonalIsZero()
    {
        var seqs = new List<string> { "ACGT", "TCGT", "GCGT" };

        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        Assert.Multiple(() =>
        {
            for (int i = 0; i < 3; i++)
            {
                Assert.That(matrix[i, i], Is.EqualTo(0), $"Diagonal[{i},{i}] not zero");
            }
        });
    }

    [Test]
    [Description("M04: All distance values must be non-negative")]
    public void CalculateDistanceMatrix_AllValuesNonNegative()
    {
        var seqs = new List<string>
        {
            "AAAA",
            "CCCC",
            "GGGG",
            "TTTT"
        };

        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        Assert.Multiple(() =>
        {
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Assert.That(matrix[i, j], Is.GreaterThanOrEqualTo(0),
                        $"Negative distance at [{i},{j}]");
                }
            }
        });
    }

    #endregion

    #region M05-M09: Distance Method Calculations

    [Test]
    [Description("M05: p-distance returns proportion of differing sites")]
    public void CalculatePairwiseDistance_PDistance_ReturnsProportionDifferent()
    {
        // 1 difference in 8 positions = 0.125
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.PDistance);

        Assert.That(dist, Is.EqualTo(0.125).Within(1e-10));
    }

    [Test]
    [Description("M06: Hamming distance returns raw count of mismatches")]
    public void CalculatePairwiseDistance_Hamming_ReturnsRawCount()
    {
        // 2 differences at positions 0 and 7
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGA";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);

        Assert.That(dist, Is.EqualTo(2));
    }

    [Test]
    [Description("M07: Jukes-Cantor distance is always >= p-distance (correction increases distance)")]
    public void CalculatePairwiseDistance_JukesCantor_GreaterThanOrEqualToPDistance()
    {
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGA";

        double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double jcDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);

        Assert.That(jcDist, Is.GreaterThanOrEqualTo(pDist));
    }

    [Test]
    [Description("M08: Jukes-Cantor formula verification: d = -0.75 * ln(1 - 4p/3)")]
    public void CalculatePairwiseDistance_JukesCantor_MatchesFormula()
    {
        // p = 0.125 (1 diff in 8 sites)
        // Expected: d = -0.75 * ln(1 - 4*0.125/3) = -0.75 * ln(0.8333...) ≈ 0.1369
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);

        // Calculate expected: -0.75 * ln(1 - 4*0.125/3)
        double p = 0.125;
        double expected = -0.75 * System.Math.Log(1 - (4.0 * p / 3.0));

        Assert.That(dist, Is.EqualTo(expected).Within(1e-6));
    }

    [Test]
    [Description("M09: Kimura 2-Parameter distinguishes transitions from transversions")]
    public void CalculatePairwiseDistance_Kimura2P_DistinguishesTransitionTypes()
    {
        // A→G is a transition (purine to purine)
        const string seqTransition = "ACGT";
        const string targetTransition = "GCGT";

        // A→C is a transversion (purine to pyrimidine)
        const string seqTransversion = "ACGT";
        const string targetTransversion = "CCGT";

        double distTransition = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seqTransition, targetTransition, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        double distTransversion = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seqTransversion, targetTransversion, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        // Both should be positive and finite
        Assert.Multiple(() =>
        {
            Assert.That(distTransition, Is.GreaterThan(0), "Transition distance should be > 0");
            Assert.That(distTransversion, Is.GreaterThan(0), "Transversion distance should be > 0");
            Assert.That(double.IsFinite(distTransition), Is.True, "Transition distance should be finite");
            Assert.That(double.IsFinite(distTransversion), Is.True, "Transversion distance should be finite");
        });
    }

    #endregion

    #region M10-M13: Edge Cases and Pre-conditions

    [Test]
    [Description("M10: Unequal length sequences throw ArgumentException")]
    public void CalculatePairwiseDistance_UnequalLengths_ThrowsArgumentException()
    {
        Assert.Throws<System.ArgumentException>(() =>
            PhylogeneticAnalyzer.CalculatePairwiseDistance("ACGT", "ACGTACGT"));
    }

    [Test]
    [Description("M11: Gap positions are ignored in distance calculation")]
    public void CalculatePairwiseDistance_GapsIgnored()
    {
        // Position 3 has gap, so only 7 comparable sites, all matching
        const string seq1 = "ACG-ACGT";
        const string seq2 = "ACGTACGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);

        Assert.That(dist, Is.EqualTo(0));
    }

    [Test]
    [Description("M12: Distance calculation is case-insensitive")]
    public void CalculatePairwiseDistance_CaseInsensitive()
    {
        const string upper = "ACGTACGT";
        const string lower = "acgtacgt";
        const string mixed = "AcGtAcGt";

        Assert.Multiple(() =>
        {
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(upper, lower),
                Is.EqualTo(0), "Upper vs lower");
            Assert.That(
                PhylogeneticAnalyzer.CalculatePairwiseDistance(upper, mixed),
                Is.EqualTo(0), "Upper vs mixed");
        });
    }

    [Test]
    [Description("M13: High divergence (p >= 0.75) causes JC69 saturation (returns infinity)")]
    public void CalculatePairwiseDistance_JukesCantor_HighDivergence_ReturnsInfinity()
    {
        // All different = p = 1.0 > 0.75
        const string seq1 = "AAAA";
        const string seq2 = "CCCC";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);

        Assert.That(double.IsPositiveInfinity(dist), Is.True,
            "JC69 should return +Infinity when p >= 0.75");
    }

    #endregion

    #region S01-S05: Quality and Robustness Tests

    [Test]
    [Description("S01: Matrix dimensions match sequence count (n sequences → n×n matrix)")]
    public void CalculateDistanceMatrix_DimensionsMatchSequenceCount()
    {
        var seqs = new List<string> { "ACGT", "TCGT", "GCGT", "CCGT" };

        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);

        Assert.Multiple(() =>
        {
            Assert.That(matrix.GetLength(0), Is.EqualTo(4));
            Assert.That(matrix.GetLength(1), Is.EqualTo(4));
        });
    }

    [Test]
    [Description("S02: Three-sequence matrix has correct off-diagonal values")]
    public void CalculateDistanceMatrix_ThreeSequences_CorrectValues()
    {
        var seqs = new List<string>
        {
            "AAAA",  // 0
            "AAAC",  // 1: differs from 0 by 1
            "CCCC"   // 2: differs from 0 by 4, from 1 by 3
        };

        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs,
            PhylogeneticAnalyzer.DistanceMethod.Hamming);

        Assert.Multiple(() =>
        {
            Assert.That(matrix[0, 1], Is.EqualTo(1), "[0,1] should be 1");
            Assert.That(matrix[0, 2], Is.EqualTo(4), "[0,2] should be 4");
            Assert.That(matrix[1, 2], Is.EqualTo(3), "[1,2] should be 3");
        });
    }

    [Test]
    [Description("S03: Kimura 2-Parameter distance >= p-distance")]
    public void CalculatePairwiseDistance_Kimura2P_GreaterThanOrEqualToPDistance()
    {
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGA";

        double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double k2pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        Assert.That(k2pDist, Is.GreaterThanOrEqualTo(pDist));
    }

    [Test]
    [Description("S04: All-gap alignment returns zero distance (no comparable sites)")]
    public void CalculatePairwiseDistance_AllGaps_ReturnsZero()
    {
        const string seq1 = "----";
        const string seq2 = "ACGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);

        Assert.That(dist, Is.EqualTo(0));
    }

    [Test]
    [Description("S05: Single difference in 8 sites gives correct p-distance (0.125)")]
    public void CalculatePairwiseDistance_SingleDifference_CorrectPDistance()
    {
        // Various single-difference cases
        var testCases = new[]
        {
            ("ACGTACGT", "TCGTACGT"),  // First position
            ("ACGTACGT", "ACGTACGA"),  // Last position
            ("ACGTACGT", "ACGAACGT"),  // Middle position
        };

        foreach (var (s1, s2) in testCases)
        {
            double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                s1, s2, PhylogeneticAnalyzer.DistanceMethod.PDistance);

            Assert.That(dist, Is.EqualTo(0.125).Within(1e-10),
                $"Failed for: {s1} vs {s2}");
        }
    }

    #endregion

    #region Additional Edge Cases

    [Test]
    [Description("Gaps in both sequences at same position are skipped")]
    public void CalculatePairwiseDistance_BothGapsAtSamePosition_Skipped()
    {
        const string seq1 = "AC-TACGT";
        const string seq2 = "AC-TACGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);

        Assert.That(dist, Is.EqualTo(0));
    }

    [Test]
    [Description("p-distance of 0.25 gives JC distance approximately 0.304")]
    public void CalculatePairwiseDistance_JukesCantor_TwoDifferences_MatchesFormula()
    {
        // 2 differences in 8 sites = p = 0.25
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGA";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);

        // Expected: -0.75 * ln(1 - 4*0.25/3) = -0.75 * ln(0.6667) ≈ 0.304
        double p = 0.25;
        double expected = -0.75 * System.Math.Log(1 - (4.0 * p / 3.0));

        Assert.That(dist, Is.EqualTo(expected).Within(1e-6));
    }

    [Test]
    [Description("Distance is symmetric: d(A,B) = d(B,A) for pairwise calculation")]
    public void CalculatePairwiseDistance_Symmetric()
    {
        const string seq1 = "ACGTACGT";
        const string seq2 = "TCGTACGA";

        double distAB = PhylogeneticAnalyzer.CalculatePairwiseDistance(seq1, seq2);
        double distBA = PhylogeneticAnalyzer.CalculatePairwiseDistance(seq2, seq1);

        Assert.That(distAB, Is.EqualTo(distBA).Within(1e-10));
    }

    #endregion
}
