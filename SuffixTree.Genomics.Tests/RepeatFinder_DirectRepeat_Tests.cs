using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for RepeatFinder.FindDirectRepeats (Test Unit REP-DIRECT-001).
/// 
/// Direct repeats are identical sequences appearing multiple times in the same orientation.
/// Example: 5' TTACG------TTACG 3' where ------ is the spacing region.
/// 
/// Sources:
/// - Wikipedia: Direct repeat, Repeated sequence (DNA)
/// - Ussery et al. (2009): Computing for Comparative Microbial Genomics
/// - Richard (2021): PMC8145212 - Trinucleotide repeat expansions
/// </summary>
[TestFixture]
public class RepeatFinder_DirectRepeat_Tests
{
    #region MUST Tests - Core Algorithm

    /// <summary>
    /// M1: Core algorithm detects identical sequences at two positions.
    /// Evidence: Wikipedia - Direct repeat definition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SimpleRepeat_FindsMatchingPair()
    {
        // Arrange: "ACGTA" appears at position 0 and position 9 with spacing of 4
        var sequence = new DnaSequence("ACGTATTTTACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(results[0].RepeatSequence, Is.EqualTo("ACGTA"));
            Assert.That(results[0].FirstPosition, Is.EqualTo(0));
            Assert.That(results[0].SecondPosition, Is.EqualTo(9));
            Assert.That(results[0].Length, Is.EqualTo(5));
            Assert.That(results[0].Spacing, Is.EqualTo(4));
        });
    }

    /// <summary>
    /// M2: Adjacent repeats (minSpacing=0) should be detected.
    /// Evidence: Wikipedia - tandem direct repeats.
    /// </summary>
    [Test]
    public void FindDirectRepeats_AdjacentRepeats_WithZeroSpacing_Found()
    {
        // Arrange: "ACGTA" immediately followed by "ACGTA"
        var sequence = new DnaSequence("ACGTAACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        // Assert
        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
        var repeat = results.First(r => r.RepeatSequence == "ACGTA");
        Assert.Multiple(() =>
        {
            Assert.That(repeat.FirstPosition, Is.EqualTo(0));
            Assert.That(repeat.SecondPosition, Is.EqualTo(5));
            Assert.That(repeat.Spacing, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M3: Sequence without repeated patterns returns empty.
    /// </summary>
    [Test]
    public void FindDirectRepeats_NoRepeats_ReturnsEmpty()
    {
        // Arrange: No 5+ bp pattern repeats
        var sequence = new DnaSequence("ACGTACGT");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M4: Empty input returns empty enumerable.
    /// </summary>
    [Test]
    public void FindDirectRepeats_EmptySequence_ReturnsEmpty()
    {
        // Act
        var results = RepeatFinder.FindDirectRepeats("", 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    #endregion

    #region MUST Tests - Parameter Validation

    /// <summary>
    /// M5: Null DnaSequence throws ArgumentNullException.
    /// </summary>
    [Test]
    public void FindDirectRepeats_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindDirectRepeats((DnaSequence)null!, 5, 10, 1).ToList());
    }

    /// <summary>
    /// M6: minLength less than 2 throws exception.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinLengthTooSmall_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindDirectRepeats(sequence, 1, 10, 1).ToList());
    }

    /// <summary>
    /// M7: maxLength less than minLength throws exception.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MaxLengthLessThanMinLength_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindDirectRepeats(sequence, 10, 5, 1).ToList());
    }

    #endregion

    #region MUST Tests - Invariants

    /// <summary>
    /// M8: Spacing = SecondPosition - FirstPosition - Length.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SpacingCalculation_MatchesFormula()
    {
        // Arrange: Various repeats with known spacing
        var sequence = new DnaSequence("ACGTATTTTACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert: Verify invariant for all results
        foreach (var result in results)
        {
            int expectedSpacing = result.SecondPosition - result.FirstPosition - result.Length;
            Assert.That(result.Spacing, Is.EqualTo(expectedSpacing),
                $"Invariant violated for repeat at {result.FirstPosition}->{result.SecondPosition}");
        }
    }

    /// <summary>
    /// M9: FirstPosition is always less than SecondPosition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_FirstPosition_AlwaysLessThanSecondPosition()
    {
        var sequence = new DnaSequence("ACGTAACGTATTTTACGTA");

        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        foreach (var result in results)
        {
            Assert.That(result.FirstPosition, Is.LessThan(result.SecondPosition),
                $"FirstPosition ({result.FirstPosition}) should be < SecondPosition ({result.SecondPosition})");
        }
    }

    /// <summary>
    /// M10: RepeatSequence equals actual substring at FirstPosition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_RepeatSequence_MatchesSubstringAtFirstPosition()
    {
        var sequenceStr = "ACGTATTTTACGTA";
        var sequence = new DnaSequence(sequenceStr);

        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        foreach (var result in results)
        {
            string actual = sequenceStr.Substring(result.FirstPosition, result.Length);
            Assert.That(result.RepeatSequence, Is.EqualTo(actual),
                $"RepeatSequence should match substring at FirstPosition");
        }
    }

    #endregion

    #region MUST Tests - Filter Thresholds

    /// <summary>
    /// M11: Only repeats with Length >= minLength are returned.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinLength_RespectsThreshold()
    {
        // Arrange: "ACGT" (4bp) and "ACGTA" (5bp) both repeat
        var sequence = new DnaSequence("ACGTAACGTAACGT");

        // Act: minLength=5 should exclude 4bp repeats
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        // Assert
        foreach (var result in results)
        {
            Assert.That(result.Length, Is.GreaterThanOrEqualTo(5));
        }
    }

    /// <summary>
    /// M12: Only repeats with Length <= maxLength are returned.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MaxLength_RespectsThreshold()
    {
        // Arrange: Long repeat that exceeds maxLength
        var repeat = "ACGTACGTAC"; // 10bp
        var sequence = new DnaSequence(repeat + "TTTT" + repeat);

        // Act: maxLength=8 should exclude 10bp repeats
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 8, 1).ToList();

        // Assert
        foreach (var result in results)
        {
            Assert.That(result.Length, Is.LessThanOrEqualTo(8));
        }
    }

    /// <summary>
    /// M13: Only repeats with Spacing >= minSpacing are returned.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinSpacing_RespectsThreshold()
    {
        // Arrange: Adjacent repeats with spacing=0
        var sequence = new DnaSequence("ACGTAACGTA");

        // Act: minSpacing=1 should exclude adjacent repeats
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert: Should be empty because only adjacent repeat exists
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M14: Sequence shorter than 2×minLength returns empty.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SequenceTooShort_ReturnsEmpty()
    {
        // Arrange: 8bp sequence, minLength=5 requires at least 10bp (5+5)
        var sequence = new DnaSequence("ACGTACGT");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    #endregion

    #region SHOULD Tests

    /// <summary>
    /// S1: Sequence with multiple repeat occurrences finds pairs.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MultipleOccurrences_FindsAllPairs()
    {
        // Arrange: "ACGTA" appears 3 times at positions 0, 7, 14
        var sequence = new DnaSequence("ACGTATTACGTATTACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 5, 1).ToList();

        // Assert: Should find pairs with the expected repeat
        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(results.Any(r => r.RepeatSequence == "ACGTA"), Is.True,
            "Should find ACGTA as a direct repeat");
    }

    /// <summary>
    /// S2: String overload produces consistent results with DnaSequence overload.
    /// </summary>
    [Test]
    public void FindDirectRepeats_StringOverload_MatchesDnaSequenceOverload()
    {
        const string seq = "ACGTAACGTA";
        var dnaSequence = new DnaSequence(seq);

        var stringResults = RepeatFinder.FindDirectRepeats(seq, 5, 10, 0).ToList();
        var dnaResults = RepeatFinder.FindDirectRepeats(dnaSequence, 5, 10, 0).ToList();

        Assert.That(stringResults, Has.Count.EqualTo(dnaResults.Count));
        for (int i = 0; i < stringResults.Count; i++)
        {
            Assert.That(stringResults[i], Is.EqualTo(dnaResults[i]));
        }
    }

    /// <summary>
    /// S3: Lowercase input is handled correctly (case-insensitive).
    /// </summary>
    [Test]
    public void FindDirectRepeats_LowercaseInput_HandledCorrectly()
    {
        // Arrange: lowercase sequence
        var results = RepeatFinder.FindDirectRepeats("acgtaacgta", 5, 10, 0).ToList();

        // Assert: Should find repeat (normalized to uppercase)
        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(results[0].RepeatSequence, Is.EqualTo("ACGTA"));
    }

    /// <summary>
    /// S4: Repeats with large intervening region are detected.
    /// Evidence: Wikipedia - interspersed repeats.
    /// </summary>
    [Test]
    public void FindDirectRepeats_LongSpacing_Detected()
    {
        // Arrange: 20bp spacing between repeats, using alternating pattern
        var spacing = "ATGATGATGATGATGATGAT"; // 20bp alternating pattern
        var repeat = "CCCGGGCCC"; // 9bp unique repeat (not found in spacing)
        var sequence = new DnaSequence($"{repeat}{spacing}{repeat}");

        // Act: Use exact repeat length to avoid subpattern matches
        var results = RepeatFinder.FindDirectRepeats(sequence, 9, 9, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Spacing, Is.EqualTo(20));
        Assert.That(results[0].RepeatSequence, Is.EqualTo(repeat));
    }

    /// <summary>
    /// S5: CAG trinucleotide repeat (Huntington's disease related).
    /// Evidence: Richard (2021) - trinucleotide repeat disorders.
    /// </summary>
    [Test]
    public void FindDirectRepeats_BiologicalRepeat_TrinucleotideCAG()
    {
        // Arrange: CAG repeat with spacing (not tandem)
        var sequence = new DnaSequence("CAGCAGTTTTTTCAGCAG");

        // Act: Look for 6bp pattern (CAGCAG)
        var results = RepeatFinder.FindDirectRepeats(sequence, 6, 6, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatSequence, Is.EqualTo("CAGCAG"));
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// C1: Overlapping pattern positions handled correctly.
    /// </summary>
    [Test]
    public void FindDirectRepeats_OverlappingPatterns_AllReported()
    {
        // Arrange: "AAAA" can be found at positions 0,1,2,3,4 in "AAAAAAAA"
        var sequence = new DnaSequence("AAAAAATTTTAAAAAA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 4, 6, 1).ToList();

        // Assert: Multiple overlapping matches should be found
        Assert.That(results, Has.Count.GreaterThan(0));
    }

    #endregion
}
