using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Comprehensive tests for DNA/RNA complement operations.
/// Test Unit: SEQ-COMP-001
/// Evidence: Wikipedia Complementarity, Nucleic Acid Sequence (IUPAC), Biopython Bio.Seq
/// </summary>
[TestFixture]
public class SequenceExtensions_Complement_Tests
{
    #region GetComplementBase - Standard Watson-Crick Pairing

    [Test]
    [Description("MUST-01: Standard Watson-Crick base pairing - A complements to T")]
    public void GetComplementBase_Adenine_ReturnsThymine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('A'), Is.EqualTo('T'));
    }

    [Test]
    [Description("MUST-01: Standard Watson-Crick base pairing - T complements to A")]
    public void GetComplementBase_Thymine_ReturnsAdenine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('T'), Is.EqualTo('A'));
    }

    [Test]
    [Description("MUST-01: Standard Watson-Crick base pairing - G complements to C")]
    public void GetComplementBase_Guanine_ReturnsCytosine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('G'), Is.EqualTo('C'));
    }

    [Test]
    [Description("MUST-01: Standard Watson-Crick base pairing - C complements to G")]
    public void GetComplementBase_Cytosine_ReturnsGuanine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('C'), Is.EqualTo('G'));
    }

    [Test]
    [Description("MUST-01: All standard DNA bases using Assert.Multiple")]
    public void GetComplementBase_AllStandardBases_CorrectComplements()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetComplementBase('A'), Is.EqualTo('T'), "A → T");
            Assert.That(SequenceExtensions.GetComplementBase('T'), Is.EqualTo('A'), "T → A");
            Assert.That(SequenceExtensions.GetComplementBase('G'), Is.EqualTo('C'), "G → C");
            Assert.That(SequenceExtensions.GetComplementBase('C'), Is.EqualTo('G'), "C → G");
        });
    }

    #endregion

    #region GetComplementBase - Case Insensitivity

    [Test]
    [Description("MUST-02: Lowercase 'a' returns uppercase 'T'")]
    public void GetComplementBase_LowercaseA_ReturnsUppercaseT()
    {
        Assert.That(SequenceExtensions.GetComplementBase('a'), Is.EqualTo('T'));
    }

    [Test]
    [Description("MUST-02: Lowercase 't' returns uppercase 'A'")]
    public void GetComplementBase_LowercaseT_ReturnsUppercaseA()
    {
        Assert.That(SequenceExtensions.GetComplementBase('t'), Is.EqualTo('A'));
    }

    [Test]
    [Description("MUST-02: Lowercase 'g' returns uppercase 'C'")]
    public void GetComplementBase_LowercaseG_ReturnsUppercaseC()
    {
        Assert.That(SequenceExtensions.GetComplementBase('g'), Is.EqualTo('C'));
    }

    [Test]
    [Description("MUST-02: Lowercase 'c' returns uppercase 'G'")]
    public void GetComplementBase_LowercaseC_ReturnsUppercaseG()
    {
        Assert.That(SequenceExtensions.GetComplementBase('c'), Is.EqualTo('G'));
    }

    [Test]
    [Description("MUST-02: All lowercase bases return uppercase complements")]
    public void GetComplementBase_AllLowercaseBases_ReturnsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetComplementBase('a'), Is.EqualTo('T'), "a → T");
            Assert.That(SequenceExtensions.GetComplementBase('t'), Is.EqualTo('A'), "t → A");
            Assert.That(SequenceExtensions.GetComplementBase('g'), Is.EqualTo('C'), "g → C");
            Assert.That(SequenceExtensions.GetComplementBase('c'), Is.EqualTo('G'), "c → G");
        });
    }

    #endregion

    #region GetComplementBase - RNA Support (Uracil)

    [Test]
    [Description("MUST-03: DNA complement supports RNA uracil - U complements to A")]
    public void GetComplementBase_Uracil_ReturnsAdenine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('U'), Is.EqualTo('A'));
    }

    [Test]
    [Description("MUST-03: Lowercase uracil also works")]
    public void GetComplementBase_LowercaseUracil_ReturnsAdenine()
    {
        Assert.That(SequenceExtensions.GetComplementBase('u'), Is.EqualTo('A'));
    }

    #endregion

    #region GetComplementBase - Involution Property

    [Test]
    [Description("MUST-04: Complement of complement equals original for A")]
    public void GetComplementBase_ComplementOfComplement_Adenine_ReturnsOriginal()
    {
        char original = 'A';
        char complement = SequenceExtensions.GetComplementBase(original);
        char doubleComplement = SequenceExtensions.GetComplementBase(complement);
        Assert.That(doubleComplement, Is.EqualTo(original));
    }

    [Test]
    [Description("MUST-04: Involution property holds for all standard bases")]
    public void GetComplementBase_InvolutionProperty_AllBases()
    {
        char[] bases = { 'A', 'T', 'G', 'C' };

        foreach (char b in bases)
        {
            char complement = SequenceExtensions.GetComplementBase(b);
            char doubleComplement = SequenceExtensions.GetComplementBase(complement);
            Assert.That(doubleComplement, Is.EqualTo(b), $"Complement(Complement({b})) should equal {b}");
        }
    }

    #endregion

    #region GetComplementBase - Unknown Base Handling

    [Test]
    [Description("MUST-05: Unknown base 'N' returns unchanged")]
    public void GetComplementBase_UnknownN_ReturnsUnchanged()
    {
        Assert.That(SequenceExtensions.GetComplementBase('N'), Is.EqualTo('N'));
    }

    [Test]
    [Description("MUST-05: Unknown base 'X' returns unchanged")]
    public void GetComplementBase_UnknownX_ReturnsUnchanged()
    {
        Assert.That(SequenceExtensions.GetComplementBase('X'), Is.EqualTo('X'));
    }

    [Test]
    [Description("MUST-05: Gap character '-' returns unchanged")]
    public void GetComplementBase_GapCharacter_ReturnsUnchanged()
    {
        Assert.That(SequenceExtensions.GetComplementBase('-'), Is.EqualTo('-'));
    }

    [Test]
    [Description("MUST-05: Various non-nucleotide characters return unchanged")]
    public void GetComplementBase_NonNucleotideCharacters_ReturnUnchanged()
    {
        char[] unknowns = { 'N', 'X', '-', '.', '?', '*' };

        foreach (char c in unknowns)
        {
            Assert.That(SequenceExtensions.GetComplementBase(c), Is.EqualTo(c),
                $"Unknown character '{c}' should return unchanged");
        }
    }

    #endregion

    #region TryGetComplement - Core Functionality

    [Test]
    [Description("MUST-07: TryGetComplement produces correct complement for standard sequence")]
    public void TryGetComplement_StandardSequence_ReturnsCorrectComplement()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        char[] buffer = new char[4];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("TGCA"));
        });
    }

    [Test]
    [Description("MUST-07: TryGetComplement handles longer sequences correctly")]
    public void TryGetComplement_LongerSequence_ReturnsCorrectComplement()
    {
        ReadOnlySpan<char> source = "AATTCCGG".AsSpan();
        char[] buffer = new char[8];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("TTAAGGCC"));
        });
    }

    [Test]
    [Description("MUST-06: TryGetComplement returns false when destination too small")]
    public void TryGetComplement_DestinationTooSmall_ReturnsFalse()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = stackalloc char[2];

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.False);
    }

    [Test]
    [Description("MUST-06: TryGetComplement returns false when destination is empty")]
    public void TryGetComplement_EmptyDestination_NonEmptySource_ReturnsFalse()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = Span<char>.Empty;

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.False);
    }

    #endregion

    #region TryGetComplement - Empty Sequence

    [Test]
    [Description("MUST-08: Empty source sequence returns true with no output")]
    public void TryGetComplement_EmptySource_ReturnsTrue()
    {
        ReadOnlySpan<char> source = ReadOnlySpan<char>.Empty;
        Span<char> destination = stackalloc char[10];

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.True);
    }

    [Test]
    [Description("MUST-08: Empty source with empty destination returns true")]
    public void TryGetComplement_EmptySourceAndDestination_ReturnsTrue()
    {
        ReadOnlySpan<char> source = ReadOnlySpan<char>.Empty;
        Span<char> destination = Span<char>.Empty;

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.True);
    }

    #endregion

    #region TryGetComplement - Buffer Size Edge Cases

    [Test]
    [Description("SHOULD-03: Destination exactly equal to source length succeeds")]
    public void TryGetComplement_DestinationExactSize_Succeeds()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = stackalloc char[4];

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.True);
    }

    [Test]
    [Description("SHOULD-04: Destination larger than source writes only source.Length characters")]
    public void TryGetComplement_DestinationLarger_WritesOnlySourceLength()
    {
        ReadOnlySpan<char> source = "AC".AsSpan();
        char[] buffer = new char[10];
        Array.Fill(buffer, 'X');
        Span<char> destination = buffer.AsSpan();

        bool success = source.TryGetComplement(destination);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(buffer[0], Is.EqualTo('T'), "First char should be complement of A");
            Assert.That(buffer[1], Is.EqualTo('G'), "Second char should be complement of C");
            Assert.That(buffer[2], Is.EqualTo('X'), "Third char should be unchanged");
        });
    }

    [Test]
    [Description("SHOULD-01: Single character sequence works correctly")]
    public void TryGetComplement_SingleCharacter_WorksCorrectly()
    {
        ReadOnlySpan<char> source = "A".AsSpan();
        char[] buffer = new char[1];

        bool success = source.TryGetComplement(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(buffer[0], Is.EqualTo('T'));
        });
    }

    #endregion

    #region TryGetComplement - Mixed Case and Special Characters

    [Test]
    [Description("SHOULD-02: Mixed case sequence produces correct uppercase complement")]
    public void TryGetComplement_MixedCase_ProducesUppercaseComplement()
    {
        ReadOnlySpan<char> source = "AcGt".AsSpan();
        char[] buffer = new char[4];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("TGCA"));
        });
    }

    [Test]
    [Description("SHOULD-05: All same base sequence complements correctly")]
    public void TryGetComplement_AllAdenine_ReturnsAllThymine()
    {
        ReadOnlySpan<char> source = "AAAA".AsSpan();
        char[] buffer = new char[4];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("TTTT"));
        });
    }

    [Test]
    [Description("SHOULD-05: All G sequence complements to all C")]
    public void TryGetComplement_AllGuanine_ReturnsAllCytosine()
    {
        ReadOnlySpan<char> source = "GGGG".AsSpan();
        char[] buffer = new char[4];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("CCCC"));
        });
    }

    [Test]
    [Description("Sequence with unknown bases preserves unknowns")]
    public void TryGetComplement_SequenceWithUnknowns_PreservesUnknowns()
    {
        ReadOnlySpan<char> source = "ACNGT".AsSpan();
        char[] buffer = new char[5];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo("TGNCA"));
        });
    }

    #endregion

    #region GetRnaComplementBase - Standard Watson-Crick Pairing (RNA)

    [Test]
    [Description("MUST-09: RNA A complements to U")]
    public void GetRnaComplementBase_Adenine_ReturnsUracil()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('A'), Is.EqualTo('U'));
    }

    [Test]
    [Description("MUST-09: RNA U complements to A")]
    public void GetRnaComplementBase_Uracil_ReturnsAdenine()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('U'), Is.EqualTo('A'));
    }

    [Test]
    [Description("MUST-09: RNA G complements to C")]
    public void GetRnaComplementBase_Guanine_ReturnsCytosine()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('G'), Is.EqualTo('C'));
    }

    [Test]
    [Description("MUST-09: RNA C complements to G")]
    public void GetRnaComplementBase_Cytosine_ReturnsGuanine()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('C'), Is.EqualTo('G'));
    }

    [Test]
    [Description("MUST-09: All standard RNA bases")]
    public void GetRnaComplementBase_AllStandardBases_CorrectComplements()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('A'), Is.EqualTo('U'), "A → U");
            Assert.That(SequenceExtensions.GetRnaComplementBase('U'), Is.EqualTo('A'), "U → A");
            Assert.That(SequenceExtensions.GetRnaComplementBase('G'), Is.EqualTo('C'), "G → C");
            Assert.That(SequenceExtensions.GetRnaComplementBase('C'), Is.EqualTo('G'), "C → G");
        });
    }

    #endregion

    #region GetRnaComplementBase - Case Insensitivity

    [Test]
    [Description("MUST-09: Lowercase RNA bases return uppercase complements")]
    public void GetRnaComplementBase_LowercaseBases_ReturnsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase('a'), Is.EqualTo('U'), "a → U");
            Assert.That(SequenceExtensions.GetRnaComplementBase('u'), Is.EqualTo('A'), "u → A");
            Assert.That(SequenceExtensions.GetRnaComplementBase('g'), Is.EqualTo('C'), "g → C");
            Assert.That(SequenceExtensions.GetRnaComplementBase('c'), Is.EqualTo('G'), "c → G");
        });
    }

    #endregion

    #region GetRnaComplementBase - Unknown Base Handling

    [Test]
    [Description("MUST-10: Unknown bases return 'N' in RNA context")]
    public void GetRnaComplementBase_UnknownBase_ReturnsN()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('X'), Is.EqualTo('N'));
    }

    [Test]
    [Description("MUST-10: Gap character returns 'N' in RNA context")]
    public void GetRnaComplementBase_GapCharacter_ReturnsN()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('-'), Is.EqualTo('N'));
    }

    [Test]
    [Description("MUST-10: Various unknown characters return 'N'")]
    public void GetRnaComplementBase_VariousUnknowns_ReturnN()
    {
        char[] unknowns = { 'X', '-', '.', '?', 'T' }; // Note: T is not valid RNA

        foreach (char c in unknowns)
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase(c), Is.EqualTo('N'),
                $"Unknown RNA character '{c}' should return 'N'");
        }
    }

    #endregion

    #region GetRnaComplementBase - Involution Property

    [Test]
    [Description("MUST-04: RNA complement involution property")]
    public void GetRnaComplementBase_InvolutionProperty_AllBases()
    {
        char[] bases = { 'A', 'U', 'G', 'C' };

        foreach (char b in bases)
        {
            char complement = SequenceExtensions.GetRnaComplementBase(b);
            char doubleComplement = SequenceExtensions.GetRnaComplementBase(complement);
            Assert.That(doubleComplement, Is.EqualTo(b), $"RNA Complement(Complement({b})) should equal {b}");
        }
    }

    #endregion
}
