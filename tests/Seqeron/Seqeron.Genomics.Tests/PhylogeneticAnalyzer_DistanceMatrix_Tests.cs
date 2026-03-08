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
        // Sequences with small divergence → finite JC distances
        var seqs = new List<string>
        {
            "ACGTACGT",
            "TCGTACGT",
            "ACGTACGA",
            "ACGTCCGT"
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
                    Assert.That(double.IsNaN(matrix[i, j]), Is.False,
                        $"NaN distance at [{i},{j}]");
                }
            }
            // Verify off-diagonal values are finite and positive
            Assert.That(matrix[0, 1], Is.GreaterThan(0).And.LessThan(double.PositiveInfinity),
                "Off-diagonal must be finite positive for non-saturated sequences");
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

        // Cross-check with independently calculated reference value (Evidence §2.3)
        // p = 1/8 → d = -3/4 × ln(5/6) ≈ 0.13674
        Assert.That(dist, Is.EqualTo(0.13674).Within(1e-4),
            "JC69 reference: p=1/8 → d ≈ 0.13674");
    }

    [Test]
    [Description("M09: Kimura 2-Parameter distinguishes transitions from transversions (K80, Kimura 1980)")]
    public void CalculatePairwiseDistance_Kimura2P_DistinguishesTransitionTypes()
    {
        // A→G is a transition (purine to purine)
        // Wikipedia: transitions are A↔G (purine↔purine) or C↔T (pyrimidine↔pyrimidine)
        const string seqTransition = "ACGT";
        const string targetTransition = "GCGT"; // 1 transition at pos 0: S=1/4, V=0

        // A→C is a transversion (purine to pyrimidine)
        const string seqTransversion = "ACGT";
        const string targetTransversion = "CCGT"; // 1 transversion at pos 0: S=0, V=1/4

        double distTransition = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seqTransition, targetTransition, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        double distTransversion = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seqTransversion, targetTransversion, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        // K2P formula (Wikipedia): K = -1/2 * ln((1-2p-q) * sqrt(1-2q))
        // Pure transition (S=0.25, V=0): d = -0.5 * ln((1-0.5-0) * sqrt(1-0)) = -0.5 * ln(0.5)
        double expectedTransition = -0.5 * System.Math.Log(0.5);
        // Pure transversion (S=0, V=0.25): d = -0.5 * ln((1-0-0.25) * sqrt(1-0.5))
        double expectedTransversion = -0.5 * System.Math.Log(0.75 * System.Math.Sqrt(0.5));

        Assert.Multiple(() =>
        {
            Assert.That(distTransition, Is.EqualTo(expectedTransition).Within(1e-10),
                "K2P transition distance must match formula: -0.5 * ln((1-2S-V) * sqrt(1-2V))");
            Assert.That(distTransversion, Is.EqualTo(expectedTransversion).Within(1e-10),
                "K2P transversion distance must match formula");
            Assert.That(distTransition, Is.Not.EqualTo(distTransversion).Within(1e-10),
                "K2P must produce different distances for transitions vs transversions");

            // Cross-check with independently calculated reference values (Evidence §5.3)
            // Pure transition (S=1/4, V=0): -1/2 × ln(1/2) ≈ 0.34657
            Assert.That(distTransition, Is.EqualTo(0.34657).Within(1e-4),
                "K2P reference: pure transition (S=0.25, V=0) → d ≈ 0.34657");
            // Pure transversion (S=0, V=1/4): -1/2 × ln(3/4 × √(1/2)) ≈ 0.31713
            Assert.That(distTransversion, Is.EqualTo(0.31713).Within(1e-4),
                "K2P reference: pure transversion (S=0, V=0.25) → d ≈ 0.31713");
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

    [Test]
    [Description("M14: Ambiguous bases (N, R, Y, etc.) are skipped like gaps")]
    public void CalculatePairwiseDistance_AmbiguousBases_SkippedLikeGaps()
    {
        // seq1: ACGTACGT — 8 standard bases
        // seq2: NCGTACGT — position 0 has 'N', should be skipped
        // Comparable sites: positions 1-7 (7 sites), all matching => 0 differences
        const string seq1 = "ACGTACGT";
        const string seq2 = "NCGTACGT";

        Assert.Multiple(() =>
        {
            double hamming = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);
            Assert.That(hamming, Is.EqualTo(0),
                "N position should be skipped, not counted as mismatch");

            double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.PDistance);
            Assert.That(pDist, Is.EqualTo(0),
                "p-distance should be 0 when only N differs");
        });

        // Also verify with IUPAC codes: R (purine), Y (pyrimidine)
        const string seq3 = "RCGTACGT";
        const string seq4 = "ACGTACYT";

        Assert.Multiple(() =>
        {
            double d1 = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                seq1, seq3, PhylogeneticAnalyzer.DistanceMethod.Hamming);
            Assert.That(d1, Is.EqualTo(0), "R position should be skipped");

            double d2 = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                seq1, seq4, PhylogeneticAnalyzer.DistanceMethod.Hamming);
            Assert.That(d2, Is.EqualTo(0), "Y position should be skipped");
        });
    }

    [Test]
    [Description("M15: Null sequences throw ArgumentNullException")]
    public void CalculatePairwiseDistance_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentNullException>(() =>
                PhylogeneticAnalyzer.CalculatePairwiseDistance(null!, "ACGT"));
            Assert.Throws<System.ArgumentNullException>(() =>
                PhylogeneticAnalyzer.CalculatePairwiseDistance("ACGT", null!));
        });
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

    #region S06-S07 and Additional Verification

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
    [Description("S06: K2P formula with mixed transitions and transversions (Kimura 1980)")]
    public void CalculatePairwiseDistance_Kimura2P_MixedChanges_MatchesFormula()
    {
        // ACGTACGT vs GCGTTCGT:
        // pos 0: A→G (transition, purine→purine), pos 4: A→T (transversion, purine→pyrimidine)
        // 8 comparable sites, S = 1/8, V = 1/8
        const string seq1 = "ACGTACGT";
        const string seq2 = "GCGTTCGT";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        // K2P formula (Wikipedia, Kimura 1980): K = -1/2 * ln((1-2S-V) * sqrt(1-2V))
        double s = 1.0 / 8;
        double v = 1.0 / 8;
        double expected = -0.5 * System.Math.Log((1 - 2 * s - v) * System.Math.Sqrt(1 - 2 * v));

        Assert.That(dist, Is.EqualTo(expected).Within(1e-10));
    }

    [Test]
    [Description("S07: K2P saturation: high transversion proportion (V≥0.5) returns infinity")]
    public void CalculatePairwiseDistance_Kimura2P_HighDivergence_ReturnsInfinity()
    {
        // AAAA vs CCCC: A→C is transversion (purine→pyrimidine), V=1.0
        // K2P: 1-2V = 1-2 = -1 ≤ 0 → formula undefined → +∞
        const string seq1 = "AAAA";
        const string seq2 = "CCCC";

        double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        Assert.That(double.IsPositiveInfinity(dist), Is.True,
            "K2P should return +Infinity when formula arguments ≤ 0 (saturation)");
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

    #region C01-C03: Extended Scenarios

    [Test]
    [Description("C01: Distance matrix computation scales for 100 sequences")]
    public void CalculateDistanceMatrix_100Sequences_CompletesSuccessfully()
    {
        // Generate 100 random 200bp DNA sequences with ~10% divergence
        var random = new System.Random(42);
        const string bases = "ACGT";
        var reference = new string(Enumerable.Range(0, 200)
            .Select(_ => bases[random.Next(4)]).ToArray());

        var seqs = new List<string> { reference };
        for (int i = 1; i < 100; i++)
        {
            var chars = reference.ToCharArray();
            for (int j = 0; j < chars.Length; j++)
            {
                if (random.NextDouble() < 0.10)
                    chars[j] = bases[random.Next(4)];
            }
            seqs.Add(new string(chars));
        }

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(seqs);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(matrix.GetLength(0), Is.EqualTo(100));
            Assert.That(matrix.GetLength(1), Is.EqualTo(100));
            // Should complete in reasonable time (< 5 seconds even on slow hardware)
            Assert.That(sw.ElapsedMilliseconds, Is.LessThan(5000),
                "100-sequence distance matrix should compute in < 5s");
        });
    }

    [Test]
    [Description("C02: JC >= p-distance and K2P >= p-distance for all test pairs")]
    public void CalculatePairwiseDistance_AllMethods_ConsistentOrdering()
    {
        // Multiple pairs with varying divergence, all below JC saturation
        var pairs = new[]
        {
            ("ACGTACGT", "TCGTACGT"),   // 1 diff, transversion (p=1/8)
            ("ACGTACGT", "GCGTACGT"),   // 1 diff, transition (p=1/8)
            ("ACGTACGT", "TCGTACGA"),   // 2 diffs (p=2/8)
            ("ACGTACGT", "TCGAACGA"),   // 3 diffs (p=3/8)
        };

        foreach (var (s1, s2) in pairs)
        {
            double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                s1, s2, PhylogeneticAnalyzer.DistanceMethod.PDistance);
            double jcDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                s1, s2, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
            double k2pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
                s1, s2, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

            Assert.Multiple(() =>
            {
                Assert.That(jcDist, Is.GreaterThanOrEqualTo(pDist),
                    $"JC >= p-distance failed for: {s1} vs {s2}");
                Assert.That(k2pDist, Is.GreaterThanOrEqualTo(pDist),
                    $"K2P >= p-distance failed for: {s1} vs {s2}");
            });
        }
    }

    [Test]
    [Description("C03: Complex gap patterns are handled correctly")]
    public void CalculatePairwiseDistance_MixedGaps_CorrectHandling()
    {
        // Multiple gaps in various positions
        // seq1: A - G T - C G T  (gaps at pos 1, 4)
        // seq2: A C G - A C - T  (gaps at pos 3, 6)
        // Comparable positions: 0(A=A), 2(G=G), 5(C=C), 7(T=T)
        // Positions 1,3,4,6 have at least one gap → skipped
        // 4 comparable sites, 0 differences
        const string seq1 = "A-GT-CGT";
        const string seq2 = "ACG-AC-T";

        double hamming = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.Hamming);
        double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq1, seq2, PhylogeneticAnalyzer.DistanceMethod.PDistance);

        Assert.Multiple(() =>
        {
            Assert.That(hamming, Is.EqualTo(0), "No differences at comparable sites");
            Assert.That(pDist, Is.EqualTo(0), "p-distance should be 0");
        });

        // With differences at comparable sites
        // seq3: A - G T - C G T  (gaps at pos 1, 4)
        // seq4: T C G - A C - T  (gaps at pos 3, 6; differs at pos 0: A→T)
        // Comparable positions: 0(A≠T), 2(G=G), 5(C=C), 7(T=T)
        // 4 comparable sites, 1 difference → p = 0.25
        const string seq3 = "A-GT-CGT";
        const string seq4 = "TCG-AC-T";

        double hamming2 = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq3, seq4, PhylogeneticAnalyzer.DistanceMethod.Hamming);
        double pDist2 = PhylogeneticAnalyzer.CalculatePairwiseDistance(
            seq3, seq4, PhylogeneticAnalyzer.DistanceMethod.PDistance);

        Assert.Multiple(() =>
        {
            Assert.That(hamming2, Is.EqualTo(1), "1 difference at comparable sites");
            Assert.That(pDist2, Is.EqualTo(0.25).Within(1e-10), "p = 1/4 = 0.25");
        });
    }

    #endregion
}
