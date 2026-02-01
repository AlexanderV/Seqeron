using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for sequence validation methods.
/// Test Unit: SEQ-VALID-001
/// Evidence: IUPAC-IUB 1970 standard, Wikipedia Nucleic acid notation, Bioinformatics.org IUPAC codes
/// </summary>
[TestFixture]
public class SequenceExtensions_SequenceValidation_Tests
{
    #region IsValidDna - MUST Tests

    [Test]
    [Description("M1: Empty sequence has no invalid characters (vacuous truth)")]
    public void IsValidDna_EmptySequence_ReturnsTrue()
    {
        // ASSUMPTION: Empty sequence returns true (vacuously valid)
        ReadOnlySpan<char> sequence = ReadOnlySpan<char>.Empty;
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M3: All standard DNA bases are valid (IUPAC 1970)")]
    public void IsValidDna_AllStandardBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "ACGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M5: Lowercase DNA bases are valid (case-insensitive)")]
    public void IsValidDna_LowercaseBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "acgt".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M7: Mixed case DNA bases are valid")]
    public void IsValidDna_MixedCase_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "AcGt".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M8: U (Uracil) is invalid for DNA (IUPAC: U is RNA only)")]
    public void IsValidDna_ContainsUracil_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGU".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("M10: Invalid character X returns false (IUPAC: X not standard)")]
    public void IsValidDna_InvalidCharacterX_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGX".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("M11: Numeric characters are invalid (IUPAC 1970)")]
    public void IsValidDna_NumericCharacter_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACG1".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("M12: Whitespace is invalid (IUPAC 1970)")]
    public void IsValidDna_Whitespace_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "AC GT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("M13: Ambiguity code N is invalid in strict mode")]
    public void IsValidDna_AmbiguityCodeN_ReturnsFalse()
    {
        // Implementation choice: strict validation rejects IUPAC ambiguity codes
        ReadOnlySpan<char> sequence = "ACGN".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("M14: Single valid base A is valid")]
    public void IsValidDna_SingleValidBase_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "A".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M15: Single invalid base X is invalid")]
    public void IsValidDna_SingleInvalidBase_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "X".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    #endregion

    #region IsValidRna - MUST Tests

    [Test]
    [Description("M2: Empty sequence has no invalid characters (vacuous truth)")]
    public void IsValidRna_EmptySequence_ReturnsTrue()
    {
        // ASSUMPTION: Empty sequence returns true (vacuously valid)
        ReadOnlySpan<char> sequence = ReadOnlySpan<char>.Empty;
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("M4: All standard RNA bases are valid (IUPAC 1970)")]
    public void IsValidRna_AllStandardBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "ACGU".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("M6: Lowercase RNA bases are valid (case-insensitive)")]
    public void IsValidRna_LowercaseBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "acgu".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("M9: T (Thymine) is invalid for RNA (IUPAC: T is DNA only)")]
    public void IsValidRna_ContainsThymine_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGT".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.False);
    }

    [Test]
    [Description("RNA: Ambiguity code N is invalid in strict mode")]
    public void IsValidRna_AmbiguityCodeN_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGUN".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.False);
    }

    #endregion

    #region IsValidDna - SHOULD Tests (Boundary & Edge Cases)

    [Test]
    [Description("S1: Long valid DNA sequence validates correctly")]
    public void IsValidDna_LongValidSequence_ReturnsTrue()
    {
        // 1000+ character sequence
        string longSequence = new string('A', 500) + new string('C', 500) + new string('G', 500) + new string('T', 500);
        ReadOnlySpan<char> sequence = longSequence.AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("S2: Invalid character at start position")]
    public void IsValidDna_InvalidAtStart_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "XACGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("S3: Invalid character at end position")]
    public void IsValidDna_InvalidAtEnd_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGTX".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("S4: Invalid character in middle position")]
    public void IsValidDna_InvalidInMiddle_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACXGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("S5: Sequence with all same valid base")]
    public void IsValidDna_AllSameBase_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "AAAA".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("S6: Special characters are invalid")]
    public void IsValidDna_SpecialCharacters_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "AC@T".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("All four valid bases in each position")]
    [TestCase("A")]
    [TestCase("C")]
    [TestCase("G")]
    [TestCase("T")]
    public void IsValidDna_EachValidBase_ReturnsTrue(string baseChar)
    {
        ReadOnlySpan<char> sequence = baseChar.AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("IUPAC ambiguity codes are invalid in strict mode")]
    [TestCase("R")] // Purine
    [TestCase("Y")] // Pyrimidine
    [TestCase("S")] // Strong
    [TestCase("W")] // Weak
    [TestCase("K")] // Keto
    [TestCase("M")] // Amino
    [TestCase("B")] // Not A
    [TestCase("D")] // Not C
    [TestCase("H")] // Not G
    [TestCase("V")] // Not T
    public void IsValidDna_IupacAmbiguityCodes_ReturnsFalse(string ambiguityCode)
    {
        ReadOnlySpan<char> sequence = ambiguityCode.AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    #endregion

    #region IsValidDna - COULD Tests (Additional Coverage)

    [Test]
    [Description("C2: Tab character is invalid whitespace")]
    public void IsValidDna_TabCharacter_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "AC\tGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("C3: Newline is invalid whitespace")]
    public void IsValidDna_NewlineCharacter_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "AC\nGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("Gap character (-) is invalid in strict mode")]
    public void IsValidDna_GapCharacter_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "AC-GT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    #endregion

    #region IsValidRna - SHOULD Tests

    [Test]
    [Description("All four valid RNA bases individually")]
    [TestCase("A")]
    [TestCase("C")]
    [TestCase("G")]
    [TestCase("U")]
    public void IsValidRna_EachValidBase_ReturnsTrue(string baseChar)
    {
        ReadOnlySpan<char> sequence = baseChar.AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("Long valid RNA sequence validates correctly")]
    public void IsValidRna_LongValidSequence_ReturnsTrue()
    {
        string longSequence = new string('A', 500) + new string('C', 500) + new string('G', 500) + new string('U', 500);
        ReadOnlySpan<char> sequence = longSequence.AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    #endregion

    #region Invariant Tests

    [Test]
    [Description("INV-3: Case invariance - uppercase and lowercase produce same result")]
    public void IsValidDna_CaseInvariance_UpperLowerSameResult()
    {
        string upper = "ACGT";
        string lower = "acgt";
        string mixed = "AcGt";

        Assert.Multiple(() =>
        {
            Assert.That(upper.AsSpan().IsValidDna(), Is.EqualTo(lower.AsSpan().IsValidDna()));
            Assert.That(upper.AsSpan().IsValidDna(), Is.EqualTo(mixed.AsSpan().IsValidDna()));
        });
    }

    [Test]
    [Description("INV-4: RNA case invariance")]
    public void IsValidRna_CaseInvariance_UpperLowerSameResult()
    {
        string upper = "ACGU";
        string lower = "acgu";

        Assert.That(upper.AsSpan().IsValidRna(), Is.EqualTo(lower.AsSpan().IsValidRna()));
    }

    #endregion
}
