using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Test suite for PAT-APPROX-001: Approximate Matching (Hamming Distance).
/// 
/// Covers:
/// - ApproximateMatcher.HammingDistance (canonical)
/// - ApproximateMatcher.FindWithMismatches (canonical)
/// - SequenceExtensions.HammingDistance span API (smoke only - tested here for completeness)
/// 
/// Evidence sources:
/// - Wikipedia: Hamming distance (https://en.wikipedia.org/wiki/Hamming_distance)
/// - Rosalind: Counting Point Mutations (https://rosalind.info/problems/hamm/)
/// - Gusfield (1997): Algorithms on Strings, Trees and Sequences
/// </summary>
[TestFixture]
[Category("PAT-APPROX-001")]
[Category("Pattern Matching")]
[Category("Hamming Distance")]
public class ApproximateMatcher_HammingDistance_Tests
{
    #region HammingDistance - Basic Functionality

    [Test]
    [Description("Identical strings have Hamming distance 0 (identity property)")]
    public void HammingDistance_IdenticalStrings_ReturnsZero()
    {
        // Evidence: Wikipedia - "d(s,t) = 0 if and only if s = t"
        Assert.That(ApproximateMatcher.HammingDistance("ACGT", "ACGT"), Is.EqualTo(0));
    }

    [Test]
    [Description("Strings differing in one position have distance 1")]
    public void HammingDistance_OneDifference_ReturnsOne()
    {
        // Last character differs: T vs G
        Assert.That(ApproximateMatcher.HammingDistance("ACGT", "ACGG"), Is.EqualTo(1));
    }

    [Test]
    [Description("Completely different strings have distance equal to length")]
    public void HammingDistance_AllDifferent_ReturnsLength()
    {
        Assert.That(ApproximateMatcher.HammingDistance("AAAA", "TTTT"), Is.EqualTo(4));
    }

    [Test]
    [Description("Hamming distance is case-insensitive for DNA sequences")]
    public void HammingDistance_CaseInsensitive_ReturnsZero()
    {
        Assert.That(ApproximateMatcher.HammingDistance("acgt", "ACGT"), Is.EqualTo(0));
    }

    [Test]
    [Description("Empty strings have Hamming distance 0")]
    public void HammingDistance_EmptyStrings_ReturnsZero()
    {
        Assert.That(ApproximateMatcher.HammingDistance("", ""), Is.EqualTo(0));
    }

    #endregion

    #region HammingDistance - Error Handling

    [Test]
    [Description("Different length strings throw ArgumentException (Hamming constraint)")]
    public void HammingDistance_DifferentLengths_ThrowsArgumentException()
    {
        // Evidence: Wikipedia - "Hamming distance between two equal-length strings"
        var ex = Assert.Throws<ArgumentException>(() =>
            ApproximateMatcher.HammingDistance("ACGT", "ACG"));

        Assert.That(ex!.Message, Does.Contain("equal length"));
    }

    [Test]
    [Description("Null first string throws ArgumentNullException")]
    public void HammingDistance_NullFirstString_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            ApproximateMatcher.HammingDistance(null!, "ACGT"));
    }

    [Test]
    [Description("Null second string throws ArgumentNullException")]
    public void HammingDistance_NullSecondString_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            ApproximateMatcher.HammingDistance("ACGT", null!));
    }

    #endregion

    #region HammingDistance - Mathematical Properties (Wikipedia)

    [Test]
    [Description("Hamming distance is symmetric: d(s,t) = d(t,s)")]
    public void HammingDistance_Symmetry_BothDirectionsEqual()
    {
        // Evidence: Wikipedia - "Symmetry"
        string s1 = "ABCDEF";
        string s2 = "AXYZEF";

        int d1 = ApproximateMatcher.HammingDistance(s1, s2);
        int d2 = ApproximateMatcher.HammingDistance(s2, s1);

        Assert.That(d1, Is.EqualTo(d2));
        Assert.That(d1, Is.EqualTo(3)); // B≠X, C≠Y, D≠Z
    }

    [Test]
    [Description("Triangle inequality: d(a,c) ≤ d(a,b) + d(b,c)")]
    public void HammingDistance_TriangleInequality_Holds()
    {
        // Evidence: Wikipedia - "satisfies the triangle inequality"
        string a = "AAAA";
        string b = "AABB";
        string c = "BBBB";

        int dAB = ApproximateMatcher.HammingDistance(a, b);
        int dBC = ApproximateMatcher.HammingDistance(b, c);
        int dAC = ApproximateMatcher.HammingDistance(a, c);

        Assert.Multiple(() =>
        {
            Assert.That(dAB, Is.EqualTo(2), "d(a,b)");
            Assert.That(dBC, Is.EqualTo(2), "d(b,c)");
            Assert.That(dAC, Is.EqualTo(4), "d(a,c)");
            Assert.That(dAC, Is.LessThanOrEqualTo(dAB + dBC), "Triangle inequality");
        });
    }

    #endregion

    #region HammingDistance - Evidence-Based Test Cases

    [Test]
    [Description("Rosalind HAMM problem canonical test case")]
    public void HammingDistance_RosalindHamm_ReturnsExpectedDistance()
    {
        // Evidence: https://rosalind.info/problems/hamm/
        // Sample Dataset:
        // GAGCCTACTAACGGGAT
        // CATCGTAATGACGGCCT
        // Sample Output: 7

        string s = "GAGCCTACTAACGGGAT";
        string t = "CATCGTAATGACGGCCT";

        int distance = ApproximateMatcher.HammingDistance(s, t);

        Assert.That(distance, Is.EqualTo(7));
    }

    [Test]
    [Description("Wikipedia example: karolin vs kathrin")]
    public void HammingDistance_WikipediaExample1_ReturnsThree()
    {
        // Evidence: Wikipedia - "karolin" and "kathrin" is 3
        Assert.That(ApproximateMatcher.HammingDistance("karolin", "kathrin"), Is.EqualTo(3));
    }

    [Test]
    [Description("Wikipedia example: karolin vs kerstin")]
    public void HammingDistance_WikipediaExample2_ReturnsThree()
    {
        // Evidence: Wikipedia - "karolin" and "kerstin" is 3
        Assert.That(ApproximateMatcher.HammingDistance("karolin", "kerstin"), Is.EqualTo(3));
    }

    [Test]
    [Description("Wikipedia example: kathrin vs kerstin")]
    public void HammingDistance_WikipediaExample3_ReturnsFour()
    {
        // Evidence: Wikipedia - "kathrin" and "kerstin" is 4
        Assert.That(ApproximateMatcher.HammingDistance("kathrin", "kerstin"), Is.EqualTo(4));
    }

    #endregion

    #region HammingDistance - Span API

    [Test]
    [Description("Span-based Hamming distance matches string-based implementation")]
    public void HammingDistance_SpanApi_MatchesStringApi()
    {
        ReadOnlySpan<char> s1 = "ACGTACGT".AsSpan();
        ReadOnlySpan<char> s2 = "ACTTACGT".AsSpan();

        int spanResult = s1.HammingDistance(s2);
        int stringResult = ApproximateMatcher.HammingDistance("ACGTACGT", "ACTTACGT");

        Assert.That(spanResult, Is.EqualTo(stringResult));
        Assert.That(spanResult, Is.EqualTo(1)); // G≠T at position 2
    }

    [Test]
    [Description("Span API throws for unequal lengths")]
    public void HammingDistance_SpanApi_UnequalLengths_ThrowsArgumentException()
    {
        // Note: Cannot use ReadOnlySpan in lambda, so we use TestDelegate pattern
        Assert.Throws<ArgumentException>(HammingDistanceWithUnequalSpans);

        static void HammingDistanceWithUnequalSpans()
        {
            ReadOnlySpan<char> s1 = "ACGT".AsSpan();
            ReadOnlySpan<char> s2 = "ACG".AsSpan();
            _ = s1.HammingDistance(s2);
        }
    }

    #endregion

    #region FindWithMismatches - Basic Functionality

    [Test]
    [Description("Zero mismatches returns only exact matches")]
    public void FindWithMismatches_ZeroMismatches_ReturnsExactMatchesOnly()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGT", 0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(2));
            Assert.That(matches[0].Position, Is.EqualTo(0));
            Assert.That(matches[1].Position, Is.EqualTo(4));
            Assert.That(matches.All(m => m.Distance == 0), Is.True, "All matches should be exact");
            Assert.That(matches.All(m => m.IsExact), Is.True, "IsExact should be true");
        });
    }

    [Test]
    [Description("One mismatch allowed finds approximate matches")]
    public void FindWithMismatches_OneMismatch_FindsApproximateMatches()
    {
        // Search for "ACGG" in "ACGTACGT" with 1 mismatch
        // Position 0: ACGT vs ACGG → distance 1 (T≠G)
        // Position 4: ACGT vs ACGG → distance 1 (T≠G)
        var matches = ApproximateMatcher.FindWithMismatches("ACGTACGT", "ACGG", 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(2));
            Assert.That(matches[0].Distance, Is.EqualTo(1));
            Assert.That(matches[0].MismatchPositions, Does.Contain(3), "Mismatch at position 3 of pattern");
        });
    }

    [Test]
    [Description("Pattern with too many mismatches not found")]
    public void FindWithMismatches_TooManyMismatches_ReturnsEmpty()
    {
        // "ACGT" vs "TGCA" has distance 4 (reverse but not complement)
        var matches = ApproximateMatcher.FindWithMismatches("ACGT", "TGCA", 2).ToList();

        Assert.That(matches, Is.Empty, "Distance 4 > maxMismatches 2");
    }

    [Test]
    [Description("maxMismatches equal to pattern length matches all windows")]
    public void FindWithMismatches_MaxMismatchesEqualsPatternLength_MatchesAllWindows()
    {
        // Pattern "AB" (length 2), maxMismatches = 2 → all positions match
        var matches = ApproximateMatcher.FindWithMismatches("XXXX", "AB", 2).ToList();

        Assert.That(matches, Has.Count.EqualTo(3), "Positions 0, 1, 2 all match");
    }

    [Test]
    [Description("Pattern equals entire sequence with exact match")]
    public void FindWithMismatches_PatternEqualsSequence_ReturnsPositionZero()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACGT", "ACGT", 0).ToList();

        Assert.That(matches, Has.Count.EqualTo(1));
        Assert.That(matches[0].Position, Is.EqualTo(0));
        Assert.That(matches[0].IsExact, Is.True);
    }

    [Test]
    [Description("Pattern equals sequence with mismatches allowed")]
    public void FindWithMismatches_PatternEqualsSequence_WithMismatches_ReturnsMatch()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACGT", "AXXX", 3).ToList();

        Assert.That(matches, Has.Count.EqualTo(1));
        Assert.That(matches[0].Position, Is.EqualTo(0));
        Assert.That(matches[0].Distance, Is.EqualTo(3));
    }

    #endregion

    #region FindWithMismatches - Edge Cases

    [Test]
    [Description("Empty pattern returns empty collection")]
    public void FindWithMismatches_EmptyPattern_ReturnsEmpty()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACGT", "", 1).ToList();
        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("Empty sequence returns empty collection")]
    public void FindWithMismatches_EmptySequence_ReturnsEmpty()
    {
        var matches = ApproximateMatcher.FindWithMismatches("", "AC", 1).ToList();
        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("Pattern longer than sequence returns empty collection")]
    public void FindWithMismatches_PatternLongerThanSequence_ReturnsEmpty()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACG", "ACGT", 1).ToList();
        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("Null sequence returns empty collection")]
    public void FindWithMismatches_NullSequence_ReturnsEmpty()
    {
        var matches = ApproximateMatcher.FindWithMismatches((string)null!, "AC", 1).ToList();
        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("Null pattern returns empty collection")]
    public void FindWithMismatches_NullPattern_ReturnsEmpty()
    {
        var matches = ApproximateMatcher.FindWithMismatches("ACGT", null!, 1).ToList();
        Assert.That(matches, Is.Empty);
    }

    [Test]
    [Description("Negative maxMismatches throws ArgumentOutOfRangeException")]
    public void FindWithMismatches_NegativeMismatches_ThrowsArgumentOutOfRangeException()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            ApproximateMatcher.FindWithMismatches("ACGT", "AC", -1).ToList());
    }

    #endregion

    #region FindWithMismatches - Invariants

    [Test]
    [Description("INV-7: All result distances are ≤ maxMismatches")]
    public void FindWithMismatches_Invariant_AllDistancesWithinLimit()
    {
        const int maxMismatches = 2;
        var matches = ApproximateMatcher.FindWithMismatches(
            "ACGTACGTACGTACGT", "ACGT", maxMismatches).ToList();

        Assert.That(matches.All(m => m.Distance <= maxMismatches), Is.True,
            "All match distances must be ≤ maxMismatches");
    }

    [Test]
    [Description("INV-8: Distance equals actual Hamming distance of matched sequence")]
    public void FindWithMismatches_Invariant_DistanceMatchesActualHammingDistance()
    {
        string pattern = "ACGT";
        var matches = ApproximateMatcher.FindWithMismatches("AXGTAYGT", pattern, 2).ToList();

        foreach (var match in matches)
        {
            int actualDistance = ApproximateMatcher.HammingDistance(match.MatchedSequence, pattern);
            Assert.That(match.Distance, Is.EqualTo(actualDistance),
                $"Distance mismatch at position {match.Position}");
        }
    }

    [Test]
    [Description("INV-9: All result positions are valid")]
    public void FindWithMismatches_Invariant_AllPositionsValid()
    {
        string sequence = "ACGTACGTACGT";
        string pattern = "ACGT";
        var matches = ApproximateMatcher.FindWithMismatches(sequence, pattern, 1).ToList();

        int maxValidPosition = sequence.Length - pattern.Length;

        Assert.That(matches.All(m => m.Position >= 0 && m.Position <= maxValidPosition), Is.True,
            "All positions must be in valid range");
    }

    [Test]
    [Description("MismatchPositions correctly identify differing positions")]
    public void FindWithMismatches_MismatchPositions_AreCorrect()
    {
        // Pattern "ACGT" vs window "AXGX" → mismatches at positions 1 and 3
        var matches = ApproximateMatcher.FindWithMismatches("AXGX", "ACGT", 2).ToList();

        Assert.That(matches, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(matches[0].MismatchPositions, Does.Contain(1), "Position 1: C≠X");
            Assert.That(matches[0].MismatchPositions, Does.Contain(3), "Position 3: T≠X");
            Assert.That(matches[0].MismatchPositions.Count, Is.EqualTo(2));
        });
    }

    #endregion

    #region FindWithMismatches - DnaSequence Overload

    [Test]
    [Description("DnaSequence overload works correctly")]
    public void FindWithMismatches_DnaSequence_Works()
    {
        var dna = new DnaSequence("ACGTACGT");
        var matches = ApproximateMatcher.FindWithMismatches(dna, "ACGT", 0).ToList();

        Assert.That(matches, Has.Count.EqualTo(2));
    }

    [Test]
    [Description("DnaSequence overload produces same results as string overload")]
    public void FindWithMismatches_DnaSequence_MatchesStringOverload()
    {
        // Use valid DNA sequence only (ACGT)
        string seq = "ACGTTTACGT";
        var dna = new DnaSequence(seq);

        var stringMatches = ApproximateMatcher.FindWithMismatches(seq, "ACGT", 1).ToList();
        var dnaMatches = ApproximateMatcher.FindWithMismatches(dna, "ACGT", 1).ToList();

        Assert.That(dnaMatches.Count, Is.EqualTo(stringMatches.Count));
        for (int i = 0; i < stringMatches.Count; i++)
        {
            Assert.That(dnaMatches[i].Position, Is.EqualTo(stringMatches[i].Position));
            Assert.That(dnaMatches[i].Distance, Is.EqualTo(stringMatches[i].Distance));
        }
    }

    #endregion

    #region FindWithMismatches - Cancellation Support

    [Test]
    [Description("Cancellation token is respected")]
    public void FindWithMismatches_WithCancellation_CompletesNormally()
    {
        using var cts = new CancellationTokenSource();
        string sequence = new('A', 1000);

        var matches = ApproximateMatcher.FindWithMismatches(sequence, "AAAA", 0, cts.Token).ToList();

        Assert.That(matches.Count, Is.GreaterThan(0));
    }

    [Test]
    [Description("Cancelled token throws OperationCanceledException")]
    public void FindWithMismatches_WithCancelledToken_Throws()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        string sequence = new('A', 1000);

        Assert.Throws<OperationCanceledException>(() =>
            ApproximateMatcher.FindWithMismatches(sequence, "AAAA", 0, cts.Token).ToList());
    }

    #endregion

    #region Real-World Bioinformatics Use Cases

    [Test]
    [Description("SNP detection: finding reference allele with single mutation")]
    public void FindWithMismatches_SnpDetection_FindsMutatedSites()
    {
        // Simulating SNP detection: looking for GATC with one possible mutation
        string genome = "ATGCGATCGATCGATCGATCG";
        string reference = "GATC";

        var matches = ApproximateMatcher.FindWithMismatches(genome, reference, 1).ToList();

        // Verify we find exact matches
        var exactMatches = matches.Where(m => m.Distance == 0).ToList();
        Assert.That(exactMatches.Count, Is.GreaterThan(0), "Should find exact GATC matches");

        // All matches should have distance 0 or 1
        Assert.That(matches.All(m => m.Distance <= 1), Is.True);
    }

    [Test]
    [Description("Primer binding: finding binding sites with allowed mismatches")]
    public void FindWithMismatches_PrimerBinding_FindsBindingSites()
    {
        string template = "ATGCATGCATGCATGCATGCATGC";
        string primer = "ATGC";

        var bindings = ApproximateMatcher.FindWithMismatches(template, primer, 1).ToList();

        Assert.That(bindings.Count, Is.GreaterThan(0), "Should find binding sites");

        // Verify exact bindings exist
        var exactBindings = bindings.Where(b => b.IsExact).ToList();
        Assert.That(exactBindings.Count, Is.GreaterThan(0), "Should have exact binding sites");
    }

    #endregion
}
