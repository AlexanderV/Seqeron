using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Comprehensive tests for DNA/RNA reverse complement operations.
/// Test Unit: SEQ-REVCOMP-001
/// Evidence: Wikipedia Complementarity, Nucleic Acid Sequence (IUPAC NC-IUB 1984), Biopython Bio.Seq
/// </summary>
[TestFixture]
public class SequenceExtensions_ReverseComplement_Tests
{
    #region TryGetReverseComplement - Basic Functionality (MUST-01)

    [Test]
    [Description("MUST-01: Basic reverse complement for palindromic sequence ACGT")]
    public void TryGetReverseComplement_PalindromicSequence_ReturnsSameSequence()
    {
        // ACGT is a biological palindrome: complement is TGCA, reversed is ACGT
        char[] destination = new char[4];
        bool success = "ACGT".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("ACGT"));
    }

    [Test]
    [Description("MUST-01: Asymmetric sequence produces correct reverse complement")]
    public void TryGetReverseComplement_AsymmetricSequence_ReturnsCorrectResult()
    {
        // AACG: complement is TTGC, reversed is CGTT
        char[] destination = new char[4];
        bool success = "AACG".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("CGTT"));
    }

    [Test]
    [Description("MUST-01: Wikipedia example — complementary sequence to TTAC is GTAA")]
    public void TryGetReverseComplement_WikipediaExample_ReturnsGTAA()
    {
        // Evidence: Wikipedia Nucleic acid sequence: "complementary sequence to TTAC is GTAA"
        char[] destination = new char[4];
        bool success = "TTAC".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("GTAA"));
    }

    [Test]
    [Description("MUST-01: Longer sequence reverse complement")]
    public void TryGetReverseComplement_LongerSequence_ReturnsCorrectResult()
    {
        // AATTCCGG → complement TTAAGGCC → reversed CCGGAATT
        char[] destination = new char[8];
        bool success = "AATTCCGG".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("CCGGAATT"));
    }

    #endregion

    #region TryGetReverseComplement - Biological Palindromes (MUST-06)

    [TestCase("GAATTC", Description = "MUST-06: EcoRI recognition site")]
    [TestCase("GGATCC", Description = "MUST-06: BamHI recognition site")]
    [TestCase("AAGCTT", Description = "MUST-06: HindIII recognition site")]
    public void TryGetReverseComplement_BiologicalPalindrome_IsOwnReverseComplement(string palindrome)
    {
        // Evidence: Wikipedia — restriction enzyme palindromes
        char[] destination = new char[palindrome.Length];
        bool success = palindrome.AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(palindrome));
    }

    #endregion

    #region TryGetReverseComplement - Involution Property (MUST-02)

    [Test]
    [Description("MUST-02: ReverseComplement(ReverseComplement(x)) = x for multiple sequences")]
    public void TryGetReverseComplement_Involution_MultipleSequences()
    {
        string[] sequences = { "A", "AT", "ATG", "ATGC", "ACGT", "AACGTTAA", "ATGCATGC", "GGGGCCCC" };

        foreach (string seq in sequences)
        {
            char[] first = new char[seq.Length];
            char[] second = new char[seq.Length];

            seq.AsSpan().TryGetReverseComplement(first);
            ((ReadOnlySpan<char>)first).TryGetReverseComplement(second);

            Assert.That(new string(second), Is.EqualTo(seq),
                $"Involution failed for sequence: {seq}");
        }
    }

    #endregion

    #region TryGetReverseComplement - Edge Cases (MUST-03, MUST-04, MUST-05, SHOULD-04)

    [Test]
    [Description("MUST-03: Empty sequence returns true, destination untouched — Biopython: reverse_complement('') → ''")]
    public void TryGetReverseComplement_EmptySequence_ReturnsTrueAndBufferUntouched()
    {
        char[] destination = new char[4];
        Array.Fill(destination, 'X');

        bool success = ReadOnlySpan<char>.Empty.TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(destination[0], Is.EqualTo('X'), "Buffer should be untouched after empty input");
    }

    [Test]
    [Description("MUST-03: Empty source with empty destination succeeds")]
    public void TryGetReverseComplement_EmptySourceAndDestination_ReturnsTrue()
    {
        bool success = ReadOnlySpan<char>.Empty.TryGetReverseComplement(Span<char>.Empty);

        Assert.That(success, Is.True);
    }

    [Test]
    [Description("MUST-04: All single nucleotides return their complements")]
    public void TryGetReverseComplement_AllSingleBases_ReturnComplements()
    {
        // Evidence: Watson-Crick base pairing (Wikipedia, IUPAC)
        var testCases = new (string input, char expected)[]
        {
            ("A", 'T'),
            ("T", 'A'),
            ("G", 'C'),
            ("C", 'G'),
            ("U", 'A')  // RNA: U pairs with A (IUPAC)
        };

        foreach (var (input, expected) in testCases)
        {
            char[] destination = new char[1];
            bool success = input.AsSpan().TryGetReverseComplement(destination);

            Assert.That(success, Is.True, $"Failed for input: {input}");
            Assert.That(destination[0], Is.EqualTo(expected),
                $"Single base {input} should reverse complement to {expected}");
        }
    }

    [Test]
    [Description("MUST-05: Destination too small returns false, no partial writes")]
    public void TryGetReverseComplement_DestinationTooSmall_ReturnsFalseAndBufferUntouched()
    {
        char[] destination = new char[2];
        Array.Fill(destination, 'X');

        bool success = "ACGT".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.False);
        Assert.That(destination[0], Is.EqualTo('X'), "Buffer must not be partially written");
        Assert.That(destination[1], Is.EqualTo('X'), "Buffer must not be partially written");
    }

    [Test]
    [Description("MUST-05: Empty destination with non-empty source returns false")]
    public void TryGetReverseComplement_EmptyDestinationNonEmptySource_ReturnsFalse()
    {
        bool success = "ACGT".AsSpan().TryGetReverseComplement(Span<char>.Empty);

        Assert.That(success, Is.False);
    }

    [Test]
    [Description("SHOULD-04: Destination larger than source only writes source.Length chars")]
    public void TryGetReverseComplement_DestinationLarger_WritesCorrectly()
    {
        char[] buffer = new char[10];
        Array.Fill(buffer, 'X');

        bool success = "AT".AsSpan().TryGetReverseComplement(buffer);

        Assert.That(success, Is.True);
        Assert.That(buffer[0], Is.EqualTo('A'));  // AT reversed: TA, complement: AT
        Assert.That(buffer[1], Is.EqualTo('T'));
        Assert.That(buffer[2], Is.EqualTo('X'), "Should not overwrite beyond source length");
    }

    #endregion

    #region TryGetReverseComplement - Case Insensitivity (MUST-07)

    [TestCase("acgt", "ACGT", Description = "MUST-07: All lowercase palindromic")]
    [TestCase("AcGt", "ACGT", Description = "MUST-07: Mixed case palindromic")]
    [TestCase("aacg", "CGTT", Description = "MUST-07: Lowercase asymmetric")]
    public void TryGetReverseComplement_CaseInsensitive_ReturnsUppercase(string input, string expected)
    {
        // Evidence: IUPAC notation is uppercase; DnaSequence normalizes to uppercase
        char[] destination = new char[input.Length];
        bool success = input.AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(expected));
    }

    #endregion

    #region TryGetReverseComplement - RNA Support (MUST-08)

    [Test]
    [Description("MUST-08: RNA ACGU reverse complement via DNA rules — U→A, G→C, C→G, A→T → ACGT")]
    public void TryGetReverseComplement_RnaSequence_UsesDnaComplementRules()
    {
        // GetComplementBase is DNA-centric: A→T (not A→U)
        char[] destination = new char[4];
        bool success = "ACGU".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("ACGT"));
    }

    [Test]
    [Description("MUST-08: RNA asymmetric sequence — AACU reverse complement → AGTT")]
    public void TryGetReverseComplement_RnaAsymmetric_ReturnsCorrectResult()
    {
        // AACU: U→A, C→G, A→T, A→T → AGTT
        char[] destination = new char[4];
        bool success = "AACU".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("AGTT"));
    }

    #endregion

    #region TryGetReverseComplement - IUPAC Ambiguity Codes (MUST-09)

    [Test]
    [Description("MUST-09: IUPAC ambiguity codes complemented per NC-IUB 1984 table")]
    public void TryGetReverseComplement_IupacAmbiguityCodes_ComplementedCorrectly()
    {
        // Evidence: Wikipedia Nucleic acid notation — IUPAC complement table (NC-IUB 1984)
        var testCases = new (char input, char expected)[]
        {
            ('R', 'Y'), // Purine (A|G) → Pyrimidine (C|T)
            ('Y', 'R'), // Pyrimidine (C|T) → Purine (A|G)
            ('S', 'S'), // Strong (C|G) → Strong (C|G)
            ('W', 'W'), // Weak (A|T) → Weak (A|T)
            ('K', 'M'), // Keto (G|T) → Amino (A|C)
            ('M', 'K'), // Amino (A|C) → Keto (G|T)
            ('B', 'V'), // Not A (C|G|T) → Not T (A|C|G)
            ('D', 'H'), // Not C (A|G|T) → Not G (A|C|T)
            ('H', 'D'), // Not G (A|C|T) → Not C (A|G|T)
            ('V', 'B'), // Not T (A|C|G) → Not A (C|G|T)
            ('N', 'N'), // Any → Any
        };

        foreach (var (input, expected) in testCases)
        {
            char[] destination = new char[1];
            bool success = new ReadOnlySpan<char>(new[] { input }).TryGetReverseComplement(destination);

            Assert.That(success, Is.True, $"Failed for IUPAC code: {input}");
            Assert.That(destination[0], Is.EqualTo(expected),
                $"IUPAC {input} should complement to {expected}");
        }
    }

    [Test]
    [Description("MUST-09: IUPAC involution — complement(complement(x)) = x for all ambiguity codes")]
    public void TryGetReverseComplement_IupacInvolution_HoldsForAllCodes()
    {
        string iupacCodes = "RYSWKMBDHVN";

        foreach (char code in iupacCodes)
        {
            char[] first = new char[1];
            char[] second = new char[1];

            new ReadOnlySpan<char>(new[] { code }).TryGetReverseComplement(first);
            ((ReadOnlySpan<char>)first).TryGetReverseComplement(second);

            Assert.That(second[0], Is.EqualTo(code),
                $"IUPAC involution failed for: {code}");
        }
    }

    [Test]
    [Description("MUST-09: All 11 lowercase IUPAC codes produce uppercase complements")]
    public void TryGetReverseComplement_LowercaseIupac_ReturnsUppercase()
    {
        // All 11 IUPAC ambiguity codes in lowercase
        var testCases = new (char input, char expected)[]
        {
            ('r', 'Y'), ('y', 'R'), ('s', 'S'), ('w', 'W'),
            ('k', 'M'), ('m', 'K'), ('b', 'V'), ('d', 'H'),
            ('h', 'D'), ('v', 'B'), ('n', 'N'),
        };

        foreach (var (input, expected) in testCases)
        {
            char[] destination = new char[1];
            new ReadOnlySpan<char>(new[] { input }).TryGetReverseComplement(destination);

            Assert.That(destination[0], Is.EqualTo(expected),
                $"Lowercase IUPAC '{input}' should complement to uppercase '{expected}'");
        }
    }

    #endregion

    #region TryGetReverseComplement - Gap / Non-IUPAC Pass-Through (SHOULD-05)

    [Test]
    [Description("SHOULD-05: Gap character is preserved and position-reversed — Biopython: gaps pass through")]
    public void TryGetReverseComplement_WithGap_PreservesAndReverses()
    {
        // A-GT: T→A, G→C, -→-, A→T → AC-T
        char[] destination = new char[4];
        bool success = "A-GT".AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("AC-T"));
    }

    #endregion

    #region TryGetReverseComplement - Biopython Cross-Verification

    [Test]
    [Description("Biopython: Seq('CCCCCGATAGNR').reverse_complement() → 'YNCTATCGGGGG'")]
    public void TryGetReverseComplement_BiopythonExample_CCCCCGATAGNR()
    {
        // Evidence: Biopython Bio.Seq docs — reverse_complement() example
        const string input = "CCCCCGATAGNR";
        const string expected = "YNCTATCGGGGG";
        char[] destination = new char[input.Length];

        bool success = input.AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(expected));
    }

    [Test]
    [Description("Biopython: reverse_complement('ACTG-NH') → 'DN-CAGT'")]
    public void TryGetReverseComplement_BiopythonExample_ACTG_NH()
    {
        // Evidence: Biopython Bio.Seq standalone reverse_complement() example
        const string input = "ACTG-NH";
        const string expected = "DN-CAGT";
        char[] destination = new char[input.Length];

        bool success = input.AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(expected));
    }

    [Test]
    [Description("Biopython: complement('ACTG-NH') → 'TGAC-ND' (complement only, no reverse)")]
    public void TryGetComplement_BiopythonExample_ACTG_NH()
    {
        // Evidence: Biopython Bio.Seq complement() example
        const string input = "ACTG-NH";
        const string expected = "TGAC-ND";
        char[] destination = new char[input.Length];

        bool success = input.AsSpan().TryGetComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(expected));
    }

    #endregion

    #region TryGetReverseComplement - Longer Sequences (SHOULD-02)

    [Test]
    [Description("SHOULD-02: Asymmetric 100-base sequence produces correct reverse complement")]
    public void TryGetReverseComplement_LongAsymmetricSequence_ReturnsCorrectResult()
    {
        // A(50)C(50) → complement T(50)G(50) → reversed G(50)T(50)
        string source = new string('A', 50) + new string('C', 50);
        string expected = new string('G', 50) + new string('T', 50);
        char[] destination = new char[100];

        bool success = source.AsSpan().TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo(expected));
    }

    [Test]
    [Description("SHOULD-02: Involution holds for asymmetric 100-base sequence")]
    public void TryGetReverseComplement_LongSequence_InvolutionHolds()
    {
        // A(50)C(25)G(25) — NOT its own reverse complement
        string original = new string('A', 50) + new string('C', 25) + new string('G', 25);
        char[] first = new char[original.Length];
        char[] second = new char[original.Length];

        original.AsSpan().TryGetReverseComplement(first);
        ((ReadOnlySpan<char>)first).TryGetReverseComplement(second);

        Assert.That(new string(second), Is.EqualTo(original));
    }

    #endregion

    #region DnaSequence.GetReverseComplementString - Static Helper

    [Test]
    [Description("DnaSequence.GetReverseComplementString produces correct result")]
    public void GetReverseComplementString_BasicSequence_ReturnsCorrectResult()
    {
        string result = DnaSequence.GetReverseComplementString("AACG");

        Assert.That(result, Is.EqualTo("CGTT"));
    }

    [Test]
    [Description("DnaSequence.GetReverseComplementString handles empty string")]
    public void GetReverseComplementString_EmptyString_ReturnsEmpty()
    {
        string result = DnaSequence.GetReverseComplementString("");

        Assert.That(result, Is.Empty);
    }

    [Test]
    [Description("DnaSequence.GetReverseComplementString handles null")]
    public void GetReverseComplementString_Null_ReturnsNull()
    {
        string? result = DnaSequence.GetReverseComplementString(null!);

        Assert.That(result, Is.Null);
    }

    [Test]
    [Description("DnaSequence.GetReverseComplementString palindrome")]
    public void GetReverseComplementString_Palindrome_ReturnsSame()
    {
        string result = DnaSequence.GetReverseComplementString("GAATTC");

        Assert.That(result, Is.EqualTo("GAATTC"));
    }

    #endregion
}
