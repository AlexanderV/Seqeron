using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for Position Weight Matrix (PWM) functionality.
/// Test Unit: PAT-PWM-001
/// 
/// Evidence sources:
/// - Wikipedia: Position weight matrix (https://en.wikipedia.org/wiki/Position_weight_matrix)
/// - Kel et al. (2003): MATCH algorithm (https://pmc.ncbi.nlm.nih.gov/articles/PMC169193/)
/// - Rosalind: Consensus and Profile (https://rosalind.info/problems/cons/)
/// - Nishida et al. (2008): Pseudocounts for TF binding sites
/// </summary>
[TestFixture]
[Category("PAT-PWM-001")]
[Description("Position Weight Matrix construction and scanning tests")]
public class MotifFinder_PWM_Tests
{
    #region CreatePwm Construction Tests

    [Test]
    [Description("Single sequence creates valid PWM with length equal to sequence length")]
    public void CreatePwm_SingleSequence_CreatesValidMatrix()
    {
        // Arrange
        var sequences = new[] { "ATGC" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Length, Is.EqualTo(4), "PWM length should match sequence length");
            Assert.That(pwm.Consensus, Is.EqualTo("ATGC"), "Single sequence should be its own consensus");
            Assert.That(pwm.Matrix, Is.Not.Null, "Matrix should be initialized");
            Assert.That(pwm.Matrix.GetLength(0), Is.EqualTo(4), "Matrix should have 4 rows (A,C,G,T)");
            Assert.That(pwm.Matrix.GetLength(1), Is.EqualTo(4), "Matrix should have 4 columns");
        });
    }

    [Test]
    [Description("Multiple identical sequences produce same consensus as single sequence")]
    public void CreatePwm_MultipleIdenticalSequences_CreatesSameConsensus()
    {
        // Arrange
        var sequences = new[] { "ATGC", "ATGC", "ATGC" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert
        Assert.That(pwm.Consensus, Is.EqualTo("ATGC"));
    }

    [Test]
    [Description("Mixed sequences derive consensus from most common base at each position (Wikipedia)")]
    public void CreatePwm_MixedSequences_ConsensusFollowsMaxRule()
    {
        // Arrange - positions designed so consensus is clear
        // Position 0: A=3, T=1 → A
        // Position 1: T=4 → T
        // Position 2: G=3, C=1 → G
        // Position 3: C=4 → C
        var sequences = new[]
        {
            "ATGC",
            "ATGC",
            "ATGC",
            "TTCC"  // Different at positions 0 and 2
        };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert
        Assert.That(pwm.Consensus, Is.EqualTo("ATGC"),
            "Consensus should reflect most common base at each position");
    }

    [Test]
    [Description("Rosalind CONS problem test case - canonical bioinformatics dataset")]
    public void CreatePwm_RosalindCONS_TestCase()
    {
        // Arrange - Rosalind CONS problem sample dataset
        // Source: https://rosalind.info/problems/cons/
        var sequences = new[]
        {
            "ATCCAGCT",
            "GGGCAACT",
            "ATGGATCT",
            "AAGCAACC",
            "TTGGAACT",
            "ATGCCATT",
            "ATGGCACT"
        };
        // Expected consensus: ATGCAACT
        // Profile: A: 5 1 0 0 5 5 0 0
        //          C: 0 0 1 4 2 0 6 1
        //          G: 1 1 6 3 0 1 0 0
        //          T: 1 5 0 0 0 1 1 6

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Length, Is.EqualTo(8), "PWM length should be 8");
            Assert.That(pwm.Consensus, Is.EqualTo("ATGCAACT"),
                "Consensus should match Rosalind expected output");
        });
    }

    [Test]
    [Description("Empty sequence collection throws ArgumentException (Wikipedia: requires input)")]
    public void CreatePwm_EmptySequences_ThrowsArgumentException()
    {
        // Act & Assert
        var ex = Assert.Throws<ArgumentException>(() =>
            MotifFinder.CreatePwm(Array.Empty<string>()));

        Assert.That(ex!.Message, Does.Contain("sequence").IgnoreCase);
    }

    [Test]
    [Description("Unequal length sequences throw ArgumentException (Wikipedia: same length required)")]
    public void CreatePwm_UnequalLengths_ThrowsArgumentException()
    {
        // Arrange
        var sequences = new[] { "ATGC", "ATG" };

        // Act & Assert
        var ex = Assert.Throws<ArgumentException>(() =>
            MotifFinder.CreatePwm(sequences));

        Assert.That(ex!.Message, Does.Contain("length").IgnoreCase);
    }

    [Test]
    [Description("Null sequences collection throws ArgumentNullException")]
    public void CreatePwm_NullSequences_ThrowsArgumentNullException()
    {
        // Act & Assert
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.CreatePwm(null!));
    }

    [Test]
    [Description("Input is normalized to uppercase (CreatePwm calls ToUpperInvariant)")]
    public void CreatePwm_LowercaseInput_NormalizesToUppercase()
    {
        // Arrange
        var sequences = new[] { "atgc", "ATGC" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert
        Assert.That(pwm.Consensus, Is.EqualTo("ATGC"),
            "Lowercase input should be normalized to uppercase");
    }

    #endregion

    #region PWM Properties Tests

    [Test]
    [Description("PWM Length property equals input sequence length")]
    public void Pwm_Length_MatchesInputSequenceLength()
    {
        // Arrange & Act
        var pwm6 = MotifFinder.CreatePwm(new[] { "ATGCAT" });
        var pwm10 = MotifFinder.CreatePwm(new[] { "ATGCATGCAT" });

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pwm6.Length, Is.EqualTo(6));
            Assert.That(pwm10.Length, Is.EqualTo(10));
        });
    }

    [Test]
    [Description("MaxScore > MinScore for non-uniform matrix (M6); MaxScore >= MinScore always")]
    public void Pwm_MaxScore_GreaterThanMinScore_ForNonUniform()
    {
        // Arrange - non-uniform inputs where MaxScore must be strictly greater
        var nonUniformCases = new[]
        {
            new[] { "ATGC" },                           // Single — each position has one dominant base
            new[] { "AAAA" },                           // All-A — A dominates, others are pseudocount-only
            new[] { "ATGC", "GCTA" },                   // Mixed
            new[] { "ATAT", "TATA", "ATAT", "TATA" }    // Alternating — 2 bases dominate each position
        };

        // Act & Assert — strict inequality for non-uniform distribution
        foreach (var sequences in nonUniformCases)
        {
            var pwm = MotifFinder.CreatePwm(sequences);
            Assert.That(pwm.MaxScore, Is.GreaterThan(pwm.MinScore),
                $"MaxScore should be > MinScore for non-uniform: {string.Join(",", sequences)}");
        }
    }

    [Test]
    [Description("Log-odds scoring: perfect match position gets highest score (Wikipedia formula)")]
    public void Pwm_LogOdds_PerfectMatchGetsMaxPositionalScore()
    {
        // Arrange - single sequence means position 0 has only 'A'
        var sequences = new[] { "AAAA" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert - A score at position 0 should be max for that position
        double aScore = pwm.Matrix[0, 0];  // Row 0 = A
        double cScore = pwm.Matrix[1, 0];  // Row 1 = C
        double gScore = pwm.Matrix[2, 0];  // Row 2 = G
        double tScore = pwm.Matrix[3, 0];  // Row 3 = T

        Assert.That(aScore, Is.GreaterThan(cScore), "A should score higher than C at position 0");
        Assert.That(aScore, Is.GreaterThan(gScore), "A should score higher than G at position 0");
        Assert.That(aScore, Is.GreaterThan(tScore), "A should score higher than T at position 0");
    }

    [Test]
    [Description("Pseudocount prevents infinite negative scores for unseen bases (Nishida 2008)")]
    public void Pwm_Pseudocount_PreventsInfiniteScores()
    {
        // Arrange - with single sequence, C,G,T are "unseen" at position 0
        var sequences = new[] { "AAAA" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert - all scores should be finite (no -∞)
        for (int pos = 0; pos < pwm.Length; pos++)
        {
            for (int baseIdx = 0; baseIdx < 4; baseIdx++)
            {
                Assert.That(double.IsFinite(pwm.Matrix[baseIdx, pos]),
                    $"Score at base {baseIdx}, position {pos} should be finite");
            }
        }
    }

    [Test]
    [Description("PWM log-odds formula: log2((count + pseudocount) / (N + 4*pseudocount) / background) — Wikipedia")]
    public void Pwm_LogOddsFormula_NumericalVerification()
    {
        // Arrange — single sequence "AG", pseudocount=0.25, background=0.25
        // Hand-calculated using Wikipedia formula:
        //   N=1, denom = 1 + 4*0.25 = 2.0
        //   freq(A,0) = (1+0.25)/2 = 0.625  → log2(0.625/0.25) = log2(2.5) ≈ 1.32193
        //   freq(C,0) = (0+0.25)/2 = 0.125  → log2(0.125/0.25) = log2(0.5)  = -1.0
        //   freq(G,0) = 0.125                → -1.0
        //   freq(T,0) = 0.125                → -1.0
        //   freq(A,1) = 0.125                → -1.0
        //   freq(G,1) = 0.625                → log2(2.5) ≈ 1.32193
        var pwm = MotifFinder.CreatePwm(new[] { "AG" }, pseudocount: 0.25);

        double expectedHigh = Math.Log2(2.5);   // ≈ 1.32193
        double expectedLow = Math.Log2(0.5);    // = -1.0

        Assert.Multiple(() =>
        {
            // Position 0: A observed
            Assert.That(pwm.Matrix[0, 0], Is.EqualTo(expectedHigh).Within(1e-10), "A at pos 0");
            Assert.That(pwm.Matrix[1, 0], Is.EqualTo(expectedLow).Within(1e-10), "C at pos 0");
            Assert.That(pwm.Matrix[2, 0], Is.EqualTo(expectedLow).Within(1e-10), "G at pos 0");
            Assert.That(pwm.Matrix[3, 0], Is.EqualTo(expectedLow).Within(1e-10), "T at pos 0");

            // Position 1: G observed
            Assert.That(pwm.Matrix[0, 1], Is.EqualTo(expectedLow).Within(1e-10), "A at pos 1");
            Assert.That(pwm.Matrix[1, 1], Is.EqualTo(expectedLow).Within(1e-10), "C at pos 1");
            Assert.That(pwm.Matrix[2, 1], Is.EqualTo(expectedHigh).Within(1e-10), "G at pos 1");
            Assert.That(pwm.Matrix[3, 1], Is.EqualTo(expectedLow).Within(1e-10), "T at pos 1");

            // MaxScore = 2 * log2(2.5), MinScore = 2 * log2(0.5)
            Assert.That(pwm.MaxScore, Is.EqualTo(2 * expectedHigh).Within(1e-10), "MaxScore");
            Assert.That(pwm.MinScore, Is.EqualTo(2 * expectedLow).Within(1e-10), "MinScore");

            Assert.That(pwm.Consensus, Is.EqualTo("AG"), "Consensus");
        });
    }

    [Test]
    [Description("Wikipedia PPM example: 10 sequences produce correct consensus (Wikipedia)")]
    public void CreatePwm_WikipediaExample_ConsensusAndScoresMatchSource()
    {
        // Wikipedia Position Weight Matrix article — 10 sequences of length 9
        // Source: https://en.wikipedia.org/wiki/Position_weight_matrix
        var sequences = new[]
        {
            "GAGGTAAAC",
            "TCCGTAAGT",
            "CAGGTTGGA",
            "ACAGTCAGT",
            "TAGGTCATT",
            "TAGGTACTG",
            "ATGGTAACT",
            "CAGGTATAC",
            "TGTGTGAGT",
            "AAGGTAAGT"
        };

        // Wikipedia PPM (without pseudocounts):
        //   A: 0.3 0.6 0.1 0.0 0.0 0.6 0.7 0.2 0.1
        //   C: 0.2 0.2 0.1 0.0 0.0 0.2 0.1 0.1 0.2
        //   G: 0.1 0.1 0.7 1.0 0.0 0.1 0.1 0.5 0.1
        //   T: 0.4 0.1 0.1 0.0 1.0 0.1 0.1 0.2 0.6
        //
        // Position 4 is all G, position 5 is all T → highest scores there

        var pwm = MotifFinder.CreatePwm(sequences, pseudocount: 0.25);

        // Verify consensus: most frequent base at each position
        // Pos 1:T(4)>A(3)>C(2)>G(1) → T; Pos 2:A(6); Pos 3:G(7); Pos 4:G(10);
        // Pos 5:T(10); Pos 6:A(6); Pos 7:A(7); Pos 8:G(5); Pos 9:T(6)
        Assert.That(pwm.Consensus, Is.EqualTo("TAGGTAAGT"),
            "Consensus should match most frequent base at each position");

        // Verify specific log-odds values for position 4 (all G, index 3)
        // N=10, pseudocount=0.25, denom = 10 + 4*0.25 = 11
        // G: (10+0.25)/11 = 10.25/11 → log2(10.25/11 / 0.25) = log2(10.25/2.75) ≈ 1.89849
        // A,C,T: (0+0.25)/11 = 0.25/11 → log2(0.25/11 / 0.25) = log2(1/11) ≈ -3.45943
        double expectedG4 = Math.Log2(10.25 / 11.0 / 0.25);
        double expectedOther4 = Math.Log2(0.25 / 11.0 / 0.25);

        Assert.Multiple(() =>
        {
            Assert.That(pwm.Matrix[2, 3], Is.EqualTo(expectedG4).Within(1e-10),
                "G at position 4 (all-G column)");
            Assert.That(pwm.Matrix[0, 3], Is.EqualTo(expectedOther4).Within(1e-10),
                "A at position 4 (unseen)");
            Assert.That(pwm.Matrix[1, 3], Is.EqualTo(expectedOther4).Within(1e-10),
                "C at position 4 (unseen)");
            Assert.That(pwm.Matrix[3, 3], Is.EqualTo(expectedOther4).Within(1e-10),
                "T at position 4 (unseen)");

            Assert.That(pwm.Length, Is.EqualTo(9));
        });
    }

    [Test]
    [Description("Scanning score equals sum of positional log-odds (Wikipedia scoring rule)")]
    public void ScanWithPwm_ScoreEqualsSumOfLogOdds()
    {
        // Arrange — build PWM from single sequence, scan with that sequence
        // Score should equal MaxScore (sum of best log-odds at each position)
        var pwm = MotifFinder.CreatePwm(new[] { "AG" }, pseudocount: 0.25);
        var sequence = new DnaSequence("AG");

        // Act
        var matches = MotifFinder.ScanWithPwm(sequence, pwm, threshold: double.MinValue).ToList();

        // Assert
        Assert.That(matches.Count, Is.EqualTo(1));
        Assert.That(matches[0].Score, Is.EqualTo(pwm.MaxScore).Within(1e-10),
            "Perfect match score should equal MaxScore (Wikipedia: sum of log-odds)");
    }

    #endregion

    #region ScanWithPwm Tests

    [Test]
    [Description("ScanWithPwm finds sequence used to train PWM")]
    public void ScanWithPwm_FindsTrainedSequence()
    {
        // Arrange
        var trainingSeq = "ATGC";
        var pwm = MotifFinder.CreatePwm(new[] { trainingSeq, trainingSeq, trainingSeq });
        var targetSequence = new DnaSequence("AAAATGCAAA");

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: 0).ToList();

        // Assert
        Assert.That(matches.Any(m => m.MatchedSequence == "ATGC"),
            "Should find the trained sequence in target");
    }

    [Test]
    [Description("ScanWithPwm returns correct positions for matches")]
    public void ScanWithPwm_ReturnsCorrectPositions()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "ATGC" });
        var targetSequence = new DnaSequence("ATGCATGC");  // ATGC at positions 0 and 4

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: 0).ToList();

        // Assert — ATGC appears at positions 0 and 4, returned in position order (S3)
        var atgcMatches = matches.Where(m => m.MatchedSequence == "ATGC").ToList();
        Assert.Multiple(() =>
        {
            Assert.That(atgcMatches.Count, Is.EqualTo(2), "Should find ATGC at exactly positions 0 and 4");
            Assert.That(atgcMatches.Select(m => m.Position), Is.EqualTo(new[] { 0, 4 }),
                "Matches must be returned in position order");
        });
    }

    [Test]
    [Description("MatchedSequence property contains correct substring from target")]
    public void ScanWithPwm_ReturnsMatchedSequence()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "AAAA" });
        var targetSequence = new DnaSequence("TTTTAAAATTTT");

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm).ToList();

        // Assert - find the AAAA match
        var aaMatches = matches.Where(m => m.MatchedSequence == "AAAA").ToList();
        Assert.That(aaMatches, Is.Not.Empty, "Should find AAAA match");
        Assert.That(aaMatches.First().Position, Is.EqualTo(4), "AAAA is at position 4");
    }

    [Test]
    [Description("Match scores are within valid range [MinScore, MaxScore]")]
    public void ScanWithPwm_ScoreWithinValidRange()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "ATGC" });
        var targetSequence = new DnaSequence("ATGCTTTTGCTA");

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: double.MinValue).ToList();

        // Assert
        foreach (var match in matches)
        {
            Assert.That(match.Score, Is.LessThanOrEqualTo(pwm.MaxScore),
                $"Score {match.Score} exceeds MaxScore {pwm.MaxScore}");
            Assert.That(match.Score, Is.GreaterThanOrEqualTo(pwm.MinScore),
                $"Score {match.Score} below MinScore {pwm.MinScore}");
        }
    }

    [Test]
    [Description("Threshold filters results correctly - only scores >= threshold returned (M9)")]
    public void ScanWithPwm_ThresholdFiltersResults()
    {
        // Arrange — PWM trained on "AAAA", scanned against mixed sequence
        var pwm = MotifFinder.CreatePwm(new[] { "AAAA" });
        var targetSequence = new DnaSequence("AAAATTTTCCCCGGGG");

        // Get all matches to establish baseline
        var allMatches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: double.MinValue).ToList();
        Assert.That(allMatches.Count, Is.GreaterThan(0), "Precondition: must have matches at MinValue threshold");

        // Use a threshold that definitely excludes some but not all: the max score minus epsilon
        // Position 0 (AAAA) is perfect match with MaxScore; others score lower
        double threshold = pwm.MaxScore - 0.001;

        // Act
        var filteredMatches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold).ToList();

        // Assert — verify actual filtering occurred AND all returned scores are >= threshold
        Assert.Multiple(() =>
        {
            Assert.That(filteredMatches.Count, Is.LessThan(allMatches.Count),
                "Filtering must actually remove some matches");
            Assert.That(filteredMatches.Count, Is.GreaterThan(0),
                "At least the perfect match should survive");
            Assert.That(filteredMatches.All(m => m.Score >= threshold), Is.True,
                "All returned matches must have score >= threshold");
        });
    }

    [Test]
    [Description("High threshold (MaxScore) returns only perfect matches (S6)")]
    public void ScanWithPwm_HighThreshold_ReturnsFewerMatches()
    {
        // Arrange — ATGC appears at positions 0 and 12; non-ATGC windows in between
        var pwm = MotifFinder.CreatePwm(new[] { "ATGC" });
        var targetSequence = new DnaSequence("ATGCTTTTGCTAATGC");

        // Act
        var lowThresholdMatches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: double.MinValue).ToList();
        var highThresholdMatches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: pwm.MaxScore).ToList();

        // Assert — strict fewer: high threshold filters out non-perfect matches
        Assert.Multiple(() =>
        {
            Assert.That(lowThresholdMatches.Count, Is.GreaterThan(highThresholdMatches.Count),
                "Low threshold must return more matches than MaxScore threshold");
            Assert.That(highThresholdMatches.Count, Is.EqualTo(2),
                "Only positions 0 and 12 are exact ATGC matches");
            Assert.That(highThresholdMatches.Select(m => m.Position), Is.EqualTo(new[] { 0, 12 }));
        });
    }

    [Test]
    [Description("Sequence shorter than PWM returns empty results")]
    public void ScanWithPwm_SequenceShorterThanPwm_ReturnsEmpty()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "ATGCATGCAT" });  // Length 10
        var targetSequence = new DnaSequence("ATGC");  // Length 4

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm).ToList();

        // Assert
        Assert.That(matches, Is.Empty, "Should return no matches when sequence is shorter than PWM");
    }

    [Test]
    [Description("Non-ACGT characters in training sequences throw ArgumentException (IUPAC-IUB: only A,C,G,T)")]
    public void CreatePwm_NonAcgtCharacters_ThrowsArgumentException()
    {
        // Arrange — N is an IUPAC ambiguity code, not a valid base for PWM training
        var sequences = new[] { "ATNG" };

        // Act & Assert
        var ex = Assert.Throws<ArgumentException>(() => MotifFinder.CreatePwm(sequences));
        Assert.That(ex!.Message, Does.Contain("Invalid character"));
    }

    [Test]
    [Description("Null sequence throws ArgumentNullException")]
    public void ScanWithPwm_NullSequence_ThrowsArgumentNullException()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "ATGC" });

        // Act & Assert
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.ScanWithPwm(null!, pwm).ToList());
    }

    [Test]
    [Description("Null PWM throws ArgumentNullException")]
    public void ScanWithPwm_NullPwm_ThrowsArgumentNullException()
    {
        // Arrange
        var sequence = new DnaSequence("ATGC");

        // Act & Assert
        Assert.Throws<ArgumentNullException>(() =>
            MotifFinder.ScanWithPwm(sequence, null!).ToList());
    }

    #endregion

    #region Invariant Tests

    [Test]
    [Description("All PWM invariants hold for valid input")]
    public void Pwm_AllInvariants_HoldForValidInput()
    {
        // Arrange
        var sequences = new[]
        {
            "ATGCATGC",
            "ATGCATGC",
            "GCTAGCTA",
            "TATATATA"
        };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert - all invariants
        Assert.Multiple(() =>
        {
            // Invariant 1: PWM Length equals input sequence length
            Assert.That(pwm.Length, Is.EqualTo(8), "PWM.Length = sequence length");

            // Invariant 2: Consensus length equals PWM length
            Assert.That(pwm.Consensus.Length, Is.EqualTo(pwm.Length), "Consensus.Length = PWM.Length");

            // Invariant 3: MaxScore >= MinScore
            Assert.That(pwm.MaxScore, Is.GreaterThanOrEqualTo(pwm.MinScore), "MaxScore >= MinScore");

            // Invariant 4: Matrix dimensions are 4 x Length
            Assert.That(pwm.Matrix.GetLength(0), Is.EqualTo(4), "Matrix rows = 4");
            Assert.That(pwm.Matrix.GetLength(1), Is.EqualTo(pwm.Length), "Matrix cols = Length");

            // Invariant 5: All scores are finite
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < pwm.Length; j++)
                {
                    Assert.That(double.IsFinite(pwm.Matrix[i, j]),
                        $"Matrix[{i},{j}] should be finite");
                }
            }

            // Invariant 6: Consensus only contains valid bases
            Assert.That(pwm.Consensus, Does.Match("^[ACGT]+$"), "Consensus contains only A,C,G,T");
        });
    }

    #endregion

    #region Additional Edge Cases

    [Test]
    [Description("PWM from uniform sequences: exact log-odds values (N=3, pseudocount=0.25)")]
    public void CreatePwm_UniformSequence_ExactLogOddsValues()
    {
        // Arrange — 3 identical "AAAA" sequences
        // N=3, p=0.25, denom=3+4*0.25=4.0
        // freq(A) = (3+0.25)/4 = 0.8125 → log2(0.8125/0.25) = log2(3.25) ≈ 1.70044
        // freq(C,G,T) = (0+0.25)/4 = 0.0625 → log2(0.0625/0.25) = log2(0.25) = -2.0
        var sequences = new[] { "AAAA", "AAAA", "AAAA" };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        double expectedHigh = Math.Log2(3.25);  // ≈ 1.70044
        double expectedLow = -2.0;               // log2(0.25)

        // Assert — exact values at every position
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Consensus, Is.EqualTo("AAAA"));
            for (int pos = 0; pos < 4; pos++)
            {
                Assert.That(pwm.Matrix[0, pos], Is.EqualTo(expectedHigh).Within(1e-10),
                    $"A at pos {pos}");
                Assert.That(pwm.Matrix[1, pos], Is.EqualTo(expectedLow).Within(1e-10),
                    $"C at pos {pos}");
                Assert.That(pwm.Matrix[2, pos], Is.EqualTo(expectedLow).Within(1e-10),
                    $"G at pos {pos}");
                Assert.That(pwm.Matrix[3, pos], Is.EqualTo(expectedLow).Within(1e-10),
                    $"T at pos {pos}");
            }
            Assert.That(pwm.MaxScore, Is.EqualTo(4 * expectedHigh).Within(1e-10), "MaxScore");
            Assert.That(pwm.MinScore, Is.EqualTo(4 * expectedLow).Within(1e-10), "MinScore");
        });
    }

    [Test]
    [Description("Maximum entropy: equal frequency at all positions → all log-odds = 0")]
    public void CreatePwm_MaximumEntropy_AllLogOddsZero()
    {
        // Arrange — each column has exactly one A, C, G, T (maximum diversity)
        // N=4, p=0.25, denom = 4 + 4*0.25 = 5.0
        // Each base: freq = (1+0.25)/5 = 0.25 → log2(0.25/0.25) = 0.0
        var sequences = new[]
        {
            "ACGT",
            "CGTA",
            "GTAC",
            "TACG"
        };

        // Act
        var pwm = MotifFinder.CreatePwm(sequences);

        // Assert — all log-odds must be exactly 0, MaxScore = MinScore = 0
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Length, Is.EqualTo(4));
            for (int pos = 0; pos < 4; pos++)
            {
                for (int baseIdx = 0; baseIdx < 4; baseIdx++)
                {
                    Assert.That(pwm.Matrix[baseIdx, pos], Is.EqualTo(0.0).Within(1e-10),
                        $"Base {baseIdx} at pos {pos} should be 0 (maximum entropy)");
                }
            }
            Assert.That(pwm.MaxScore, Is.EqualTo(0.0).Within(1e-10), "MaxScore = 0 at max entropy");
            Assert.That(pwm.MinScore, Is.EqualTo(0.0).Within(1e-10), "MinScore = 0 at max entropy");
        });
    }

    [Test]
    [Description("Multiple matches at same score are all returned")]
    public void ScanWithPwm_MultipleMatchesSameScore_AllReturned()
    {
        // Arrange
        var pwm = MotifFinder.CreatePwm(new[] { "AA" });
        var targetSequence = new DnaSequence("AAAAAA");  // AA at 0,1,2,3,4

        // Act
        var matches = MotifFinder.ScanWithPwm(targetSequence, pwm, threshold: 0).ToList();

        // Assert - should find all overlapping matches
        Assert.That(matches.Count, Is.EqualTo(5), "Should find AA at all 5 possible positions");
    }

    [Test]
    [Description("Pseudocount = 0 produces -∞ for unseen bases (Wikipedia: risk of zero probabilities)")]
    public void CreatePwm_ZeroPseudocount_ProducesNegativeInfinity()
    {
        // Arrange — single sequence "A", pseudocount=0
        // freq(A,0) = (1+0)/(1+0) = 1.0 → log2(1.0/0.25) = log2(4) = 2.0
        // freq(C,0) = 0/(1) = 0.0 → log2(0/0.25) = log2(0) = -∞
        var pwm = MotifFinder.CreatePwm(new[] { "A" }, pseudocount: 0);

        // Assert — verify exact behavior documented in Wikipedia edge case
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Matrix[0, 0], Is.EqualTo(2.0).Within(1e-10),
                "Observed base: log2(1/0.25) = 2.0");
            Assert.That(double.IsNegativeInfinity(pwm.Matrix[1, 0]),
                "C unseen with pseudocount=0 → -∞");
            Assert.That(double.IsNegativeInfinity(pwm.Matrix[2, 0]),
                "G unseen with pseudocount=0 → -∞");
            Assert.That(double.IsNegativeInfinity(pwm.Matrix[3, 0]),
                "T unseen with pseudocount=0 → -∞");
        });
    }

    [Test]
    [Description("Threshold boundary: score exactly equal to threshold is included (>= semantics)")]
    public void ScanWithPwm_ThresholdBoundary_ExactScoreIncluded()
    {
        // Arrange — build PWM from "AG", scan "AG" (only one window, score = MaxScore)
        var pwm = MotifFinder.CreatePwm(new[] { "AG" }, pseudocount: 0.25);
        var sequence = new DnaSequence("AG");

        // Use MaxScore as threshold — the perfect match score equals MaxScore exactly
        double threshold = pwm.MaxScore;

        // Act
        var matches = MotifFinder.ScanWithPwm(sequence, pwm, threshold).ToList();

        // Assert — score == threshold must be included (>= not >)
        Assert.Multiple(() =>
        {
            Assert.That(matches.Count, Is.EqualTo(1),
                "Exact score == threshold must be included (>= semantics)");
            Assert.That(matches[0].Score, Is.EqualTo(threshold).Within(1e-10),
                "Match score should equal threshold exactly");
            Assert.That(matches[0].Position, Is.EqualTo(0));
        });
    }

    #endregion
}
