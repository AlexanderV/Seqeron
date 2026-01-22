using NUnit.Framework;
using SuffixTree.Genomics;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for RepeatFinder.FindPalindromes (Test Unit REP-PALIN-001).
/// 
/// A DNA palindrome is a sequence that equals its reverse complement.
/// This is distinct from textual palindromes which read the same forwards and backwards.
/// 
/// Evidence sources:
/// - Wikipedia: Palindromic sequence (https://en.wikipedia.org/wiki/Palindromic_sequence)
/// - Wikipedia: Restriction enzyme (https://en.wikipedia.org/wiki/Restriction_enzyme)
/// - Rosalind REVP problem (https://rosalind.info/problems/revp/)
/// </summary>
[TestFixture]
[Category("RepeatFinder")]
[Category("Palindrome")]
public class RepeatFinder_Palindrome_Tests
{
    #region MUST Tests - Restriction Enzyme Recognition Sites

    /// <summary>
    /// M1: EcoRI recognition site (GAATTC) is a canonical 6bp palindrome.
    /// Source: Wikipedia - Restriction enzyme
    /// </summary>
    [Test]
    [Description("MUST-01: EcoRI recognition site GAATTC is detected as palindrome")]
    public void FindPalindromes_EcoRI_SixBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("AAAGAATTCAAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GAATTC"));
            Assert.That(results[0].Position, Is.EqualTo(3));
            Assert.That(results[0].Length, Is.EqualTo(6));
        });
    }

    /// <summary>
    /// M2: HindIII recognition site (AAGCTT) is a canonical 6bp palindrome.
    /// Source: Wikipedia - Restriction enzyme
    /// </summary>
    [Test]
    [Description("MUST-02: HindIII recognition site AAGCTT is detected as palindrome")]
    public void FindPalindromes_HindIII_SixBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("CCCAAGCTTCCC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("AAGCTT"));
            Assert.That(results[0].Position, Is.EqualTo(3));
            Assert.That(results[0].Length, Is.EqualTo(6));
        });
    }

    /// <summary>
    /// M3: BamHI recognition site (GGATCC) is a canonical 6bp palindrome.
    /// Source: Wikipedia - Restriction enzyme
    /// </summary>
    [Test]
    [Description("MUST-03: BamHI recognition site GGATCC is detected as palindrome")]
    public void FindPalindromes_BamHI_SixBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("TTTGGATCCTTT");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GGATCC"));
            Assert.That(results[0].Position, Is.EqualTo(3));
            Assert.That(results[0].Length, Is.EqualTo(6));
        });
    }

    #endregion

    #region MUST Tests - Four Base Palindromes

    /// <summary>
    /// M4: GCGC is a 4bp palindrome (reverse complement of GCGC is GCGC).
    /// Source: Rosalind REVP
    /// </summary>
    [Test]
    [Description("MUST-04: 4bp palindrome GCGC is detected")]
    public void FindPalindromes_GCGC_FourBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("AAGCGCAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 4).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GCGC"));
            Assert.That(results[0].Position, Is.EqualTo(2));
        });
    }

    /// <summary>
    /// M5: ATAT is a 4bp palindrome (reverse complement of ATAT is ATAT).
    /// Source: Standard test case
    /// </summary>
    [Test]
    [Description("MUST-05: 4bp palindrome ATAT is detected")]
    public void FindPalindromes_ATAT_FourBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("GGATATTCC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 4).ToList();

        // Assert
        Assert.That(results.Any(r => r.Sequence == "ATAT"), Is.True, "ATAT should be detected");
    }

    #endregion

    #region MUST Tests - Invariants

    /// <summary>
    /// M6: Core invariant - every detected palindrome must equal its reverse complement.
    /// Source: Wikipedia - Palindromic sequence definition
    /// </summary>
    [Test]
    [Description("MUST-06: Every palindrome equals its reverse complement")]
    public void FindPalindromes_AllResults_EqualTheirReverseComplement()
    {
        // Arrange - sequence with multiple known palindromes
        var sequence = new DnaSequence("GAATTCAAAGCGCAAAGGATCC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 12).ToList();

        // Assert
        Assert.That(results, Is.Not.Empty, "Should find at least one palindrome");

        foreach (var result in results)
        {
            var reverseComplement = DnaSequence.GetReverseComplementString(result.Sequence);
            Assert.That(result.Sequence, Is.EqualTo(reverseComplement),
                $"Sequence {result.Sequence} at position {result.Position} should equal its reverse complement");
        }
    }

    /// <summary>
    /// M14: Positions are 0-based indices.
    /// Source: Implementation contract
    /// </summary>
    [Test]
    [Description("MUST-14: Position is 0-based index")]
    public void FindPalindromes_Position_IsZeroBasedIndex()
    {
        // Arrange - palindrome at start of sequence
        var sequence = new DnaSequence("GAATTC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.That(results[0].Position, Is.EqualTo(0), "First position should be 0");
    }

    /// <summary>
    /// M15: Length property equals Sequence.Length.
    /// Source: Invariant
    /// </summary>
    [Test]
    [Description("MUST-15: Length property equals Sequence.Length")]
    public void FindPalindromes_Length_MatchesSequenceLength()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTCAAAGCGGCCGCAAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 12).ToList();

        // Assert
        Assert.That(results, Is.Not.Empty);

        foreach (var result in results)
        {
            Assert.That(result.Length, Is.EqualTo(result.Sequence.Length),
                $"Length {result.Length} should equal Sequence.Length {result.Sequence.Length}");
        }
    }

    #endregion

    #region MUST Tests - Edge Cases

    /// <summary>
    /// M7: Sequence without palindromes returns empty.
    /// Source: Standard edge case
    /// </summary>
    [Test]
    [Description("MUST-07: Non-palindromic sequence returns empty")]
    public void FindPalindromes_NoPalindromes_ReturnsEmpty()
    {
        // Arrange - poly-A has no palindromes (revcomp is poly-T)
        var sequence = new DnaSequence("AAAAAAAAAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 12).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M8: Empty sequence returns empty collection.
    /// Source: Standard boundary
    /// </summary>
    [Test]
    [Description("MUST-08: Empty sequence returns empty")]
    public void FindPalindromes_EmptySequence_ReturnsEmpty()
    {
        // Act
        var results = RepeatFinder.FindPalindromes("", minLength: 4, maxLength: 12).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M9: Null sequence throws ArgumentNullException.
    /// Source: Implementation contract
    /// </summary>
    [Test]
    [Description("MUST-09: Null sequence throws ArgumentNullException")]
    public void FindPalindromes_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindPalindromes((DnaSequence)null!, minLength: 4, maxLength: 12).ToList());
    }

    #endregion

    #region MUST Tests - Parameter Validation

    /// <summary>
    /// M10: Odd minLength throws ArgumentOutOfRangeException.
    /// Source: Implementation constraint (palindromes must be even length)
    /// </summary>
    [Test]
    [Description("MUST-10: Odd minLength throws ArgumentOutOfRangeException")]
    public void FindPalindromes_OddMinLength_ThrowsArgumentOutOfRangeException()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTC");

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindPalindromes(sequence, minLength: 5, maxLength: 12).ToList());
    }

    /// <summary>
    /// M11: minLength less than 4 throws ArgumentOutOfRangeException.
    /// Source: Implementation constraint
    /// </summary>
    [Test]
    [Description("MUST-11: minLength < 4 throws ArgumentOutOfRangeException")]
    public void FindPalindromes_MinLengthLessThan4_ThrowsArgumentOutOfRangeException()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTC");

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindPalindromes(sequence, minLength: 2, maxLength: 12).ToList());
    }

    /// <summary>
    /// M12: maxLength less than minLength throws ArgumentOutOfRangeException.
    /// Source: Implementation contract
    /// </summary>
    [Test]
    [Description("MUST-12: maxLength < minLength throws ArgumentOutOfRangeException")]
    public void FindPalindromes_MaxLengthLessThanMinLength_ThrowsArgumentOutOfRangeException()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTC");

        // Act & Assert
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindPalindromes(sequence, minLength: 8, maxLength: 4).ToList());
    }

    #endregion

    #region MUST Tests - Multiple Palindromes

    /// <summary>
    /// M13: Multiple palindromes in sequence are all detected.
    /// Source: Rosalind REVP
    /// </summary>
    [Test]
    [Description("MUST-13: Multiple palindromes are all detected")]
    public void FindPalindromes_MultiplePalindromes_FindsAll()
    {
        // Arrange - two EcoRI sites
        var sequence = new DnaSequence("GAATTCAAAGAATTC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(2));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[1].Position, Is.EqualTo(9));
        });
    }

    #endregion

    #region MUST Tests - Length Thresholds

    /// <summary>
    /// M16: Only palindromes with length >= minLength are returned.
    /// Source: Algorithm specification
    /// </summary>
    [Test]
    [Description("MUST-16: Only palindromes >= minLength are returned")]
    public void FindPalindromes_MinLength_RespectsThreshold()
    {
        // Arrange - sequence with 4bp and 6bp palindromes
        // GCGC is 4bp, GAATTC is 6bp
        var sequence = new DnaSequence("GCGCAAAGAATTC");

        // Act - only look for 6bp+
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 12).ToList();

        // Assert
        Assert.That(results.All(r => r.Length >= 6), Is.True);
        Assert.That(results.Any(r => r.Sequence == "GAATTC"), Is.True);
        Assert.That(results.Any(r => r.Sequence == "GCGC"), Is.False, "4bp palindrome should not be included");
    }

    /// <summary>
    /// M17: Only palindromes with length <= maxLength are returned.
    /// Source: Algorithm specification
    /// </summary>
    [Test]
    [Description("MUST-17: Only palindromes <= maxLength are returned")]
    public void FindPalindromes_MaxLength_RespectsThreshold()
    {
        // Arrange - sequence with 4bp, 6bp, 8bp palindromes possible
        var sequence = new DnaSequence("AAAGCGGCCGCAAA"); // NotI site is 8bp

        // Act - only look for 4-6bp
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 6).ToList();

        // Assert
        Assert.That(results.All(r => r.Length <= 6), Is.True, "All results should be <= 6bp");
        Assert.That(results.Any(r => r.Length == 8), Is.False, "8bp palindrome should not be included");
    }

    #endregion

    #region SHOULD Tests - API Consistency

    /// <summary>
    /// S1: String overload produces same results as DnaSequence overload.
    /// Source: Implementation contract
    /// </summary>
    [Test]
    [Description("SHOULD-01: String overload matches DnaSequence overload")]
    public void FindPalindromes_StringOverload_MatchesDnaSequenceOverload()
    {
        // Arrange
        var sequenceString = "GAATTC";
        var sequenceDna = new DnaSequence(sequenceString);

        // Act
        var stringResults = RepeatFinder.FindPalindromes(sequenceString, 6, 6).ToList();
        var dnaResults = RepeatFinder.FindPalindromes(sequenceDna, 6, 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(stringResults, Has.Count.EqualTo(dnaResults.Count));
            Assert.That(stringResults[0].Sequence, Is.EqualTo(dnaResults[0].Sequence));
            Assert.That(stringResults[0].Position, Is.EqualTo(dnaResults[0].Position));
        });
    }

    /// <summary>
    /// S2: Lowercase input is handled correctly (case insensitive).
    /// Source: Implementation robustness
    /// </summary>
    [Test]
    [Description("SHOULD-02: Lowercase input handled correctly")]
    public void FindPalindromes_LowercaseInput_HandledCorrectly()
    {
        // Act
        var results = RepeatFinder.FindPalindromes("gaattc", minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GAATTC"));
        });
    }

    #endregion

    #region SHOULD Tests - Rosalind Validation

    /// <summary>
    /// S3: Validate against Rosalind REVP sample dataset.
    /// Source: Rosalind REVP (https://rosalind.info/problems/revp/)
    /// 
    /// Input: TCAATGCATGCGGGTCTATATGCAT
    /// Expected palindromes at (1-based position, length):
    /// 4,6 / 5,4 / 6,6 / 7,4 / 17,4 / 18,4 / 20,6 / 21,4
    /// 
    /// Converting to 0-based: 3,6 / 4,4 / 5,6 / 6,4 / 16,4 / 17,4 / 19,6 / 20,4
    /// </summary>
    [Test]
    [Description("SHOULD-03: Rosalind REVP sample dataset validation")]
    public void FindPalindromes_RosalindSample_CorrectOutput()
    {
        // Arrange
        var sequence = new DnaSequence("TCAATGCATGCGGGTCTATATGCAT");

        // Expected palindromes (0-based position, length)
        var expected = new HashSet<(int Position, int Length)>
        {
            (3, 6),   // ATGCAT at position 4 (1-based) = 3 (0-based)
            (4, 4),   // TGCA
            (5, 6),   // GCATGC
            (6, 4),   // CATG
            (16, 4),  // ATAT
            (17, 4),  // TATA
            (19, 6),  // ATGCAT
            (20, 4),  // TGCA
        };

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 12)
            .Select(r => (r.Position, r.Length))
            .ToHashSet();

        // Assert
        foreach (var exp in expected)
        {
            Assert.That(results.Contains(exp), Is.True,
                $"Expected palindrome at position {exp.Position}, length {exp.Length}");
        }
    }

    #endregion

    #region SHOULD Tests - Additional Restriction Sites

    /// <summary>
    /// S4: NotI recognition site (GCGGCCGC) is an 8bp palindrome.
    /// Source: Wikipedia - Restriction enzyme
    /// </summary>
    [Test]
    [Description("SHOULD-04: NotI 8bp palindrome GCGGCCGC is detected")]
    public void FindPalindromes_NotI_EightBasePalindromeDetected()
    {
        // Arrange
        var sequence = new DnaSequence("AAAGCGGCCGCAAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 8, maxLength: 8).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GCGGCCGC"));
        });
    }

    #endregion

    #region SHOULD Tests - Edge Cases

    /// <summary>
    /// S5: Overlapping palindromes at different lengths are both detected.
    /// Source: Edge case
    /// </summary>
    [Test]
    [Description("SHOULD-05: Overlapping palindromes at different lengths detected")]
    public void FindPalindromes_OverlappingAtDifferentLengths_BothDetected()
    {
        // Arrange - ATGCAT contains both 4bp (TGCA) and 6bp (ATGCAT) palindromes
        var sequence = new DnaSequence("ATGCAT");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 4, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results.Any(r => r.Length == 4), Is.True, "Should find 4bp palindrome");
            Assert.That(results.Any(r => r.Length == 6), Is.True, "Should find 6bp palindrome");
        });
    }

    /// <summary>
    /// S6: When entire sequence is a palindrome, it is returned.
    /// Source: Edge case
    /// </summary>
    [Test]
    [Description("SHOULD-06: Entire palindromic sequence is returned")]
    public void FindPalindromes_EntireSequenceIsPalindrome_Returned()
    {
        // Arrange - GAATTC is a 6bp palindrome
        var sequence = new DnaSequence("GAATTC");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[0].Sequence, Is.EqualTo("GAATTC"));
        });
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// C1: 12bp palindrome at maximum default length is detected.
    /// Source: Algorithm range
    /// </summary>
    [Test]
    [Description("COULD-01: 12bp palindrome detected")]
    public void FindPalindromes_TwelveBasePalindrome_Detected()
    {
        // Arrange - construct a 12bp palindrome: GAATTCGAATTC is NOT a palindrome
        // Need: seq == revcomp(seq)
        // ACGTACGTACGT -> revcomp = ACGTACGTACGT? No, revcomp(ACGT) = ACGT
        // Let's use: GAATTCAATTCG - check if palindrome
        // revcomp(GAATTCAATTCG) = CGAATTGAATTC - not equal
        // Try: GAATTAATTAATTC - 14bp, too long
        // CGAATTAATTCG - revcomp = CGAATTAATTCG ✓
        var sequence = new DnaSequence("AAACGAATTAATTCGAAA");

        // Act
        var results = RepeatFinder.FindPalindromes(sequence, minLength: 12, maxLength: 12).ToList();

        // Assert
        Assert.That(results.Any(r => r.Length == 12 && r.Sequence == "CGAATTAATTCG"), Is.True);
    }

    /// <summary>
    /// C2: GenomicAnalyzer.FindPalindromes smoke test.
    /// Source: Implementation comparison
    /// </summary>
    [Test]
    [Description("COULD-02: GenomicAnalyzer.FindPalindromes smoke test")]
    public void GenomicAnalyzer_FindPalindromes_SmokeTest()
    {
        // Arrange
        var sequence = new DnaSequence("AAAGAATTCAAA");

        // Act
        var results = GenomicAnalyzer.FindPalindromes(sequence, minLength: 6, maxLength: 6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("GAATTC"));
            Assert.That(results[0].Position, Is.EqualTo(3));
        });
    }

    #endregion
}
