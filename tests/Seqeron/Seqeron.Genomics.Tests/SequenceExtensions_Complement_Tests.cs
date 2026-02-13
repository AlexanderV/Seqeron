using System;
using System.Linq;
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
    [Description("MUST-01: A↔T, G↔C — Wikipedia Complementarity; IUPAC table; Biopython")]
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
    [Description("MUST-02: Uppercase output for lowercase input — DnaSequence/RnaSequence normalize to uppercase")]
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
    [Description("MUST-03: U/u complement to A — IUPAC table: U complement = A; Biopython: U treated as T")]
    public void GetComplementBase_Uracil_ReturnsAdenine()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceExtensions.GetComplementBase('U'), Is.EqualTo('A'), "U → A");
            Assert.That(SequenceExtensions.GetComplementBase('u'), Is.EqualTo('A'), "u → A");
        });
    }

    #endregion

    #region GetComplementBase - Involution Property

    [Test]
    [Description("MUST-04: comp(comp(x)) = x for all standard bases — mathematical bijection property")]
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
    [Description("MUST-05: Non-nucleotide characters pass through — Biopython: complement('XYZ') → unknowns unchanged")]
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
    [Description("MUST-07: TryGetComplement produces correct complement — Biopython: complement('ACGT')→'TGCA'")]
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
    [Description("MUST-06: Returns false when destination.Length < source.Length")]
    public void TryGetComplement_DestinationTooSmall_ReturnsFalse()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = stackalloc char[2];

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.False);
    }

    #endregion

    #region TryGetComplement - Empty Sequence

    [Test]
    [Description("MUST-08: Empty source returns true, destination not modified")]
    public void TryGetComplement_EmptySource_ReturnsTrue()
    {
        ReadOnlySpan<char> source = ReadOnlySpan<char>.Empty;
        char[] buffer = new char[4];
        Array.Fill(buffer, 'X');

        bool success = source.TryGetComplement(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(buffer[0], Is.EqualTo('X'), "destination should not be modified");
        });
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

    [TestCase("AAAA", "TTTT")]
    [TestCase("GGGG", "CCCC")]
    [Description("SHOULD-05: Homogeneous sequence complement")]
    public void TryGetComplement_HomogeneousSequence_ComplementsCorrectly(string input, string expected)
    {
        ReadOnlySpan<char> source = input.AsSpan();
        char[] buffer = new char[input.Length];

        bool success = source.TryGetComplement(buffer);
        string result = new string(buffer);

        Assert.Multiple(() =>
        {
            Assert.That(success, Is.True);
            Assert.That(result, Is.EqualTo(expected));
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
    [Description("MUST-09: RNA A↔U, G↔C — Wikipedia Complementarity; Biopython complement_rna()")]
    public void GetRnaComplementBase_AllStandardBases_CorrectComplements()
    {
        // Biopython: complement_rna("AUGC") → "UACG"
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
    [Description("MUST-10: Non-nucleotide characters pass through — Biopython: complement_rna('XYZ') → unknowns unchanged")]
    public void GetRnaComplementBase_NonNucleotideCharacters_ReturnUnchanged()
    {
        char[] unknowns = { 'X', '-', '.', '?', '*' };

        foreach (char c in unknowns)
        {
            Assert.That(SequenceExtensions.GetRnaComplementBase(c), Is.EqualTo(c),
                $"Unknown RNA character '{c}' should return unchanged");
        }
    }

    [Test]
    [Description("MUST-09: T complements to A in RNA context (Biopython: complement_rna('T') → 'A')")]
    public void GetRnaComplementBase_Thymine_ReturnsAdenine()
    {
        Assert.That(SequenceExtensions.GetRnaComplementBase('T'), Is.EqualTo('A'));
        Assert.That(SequenceExtensions.GetRnaComplementBase('t'), Is.EqualTo('A'));
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

    #region Biopython Cross-Verification

    [Test]
    [Description("Biopython: complement_rna('CGAUT') → 'GCUAA' — T→A in RNA context")]
    public void GetRnaComplementBase_BiopythonCrossVerification_Cgaut()
    {
        // Biopython: complement_rna(Seq("CGAUT")) → Seq('GCUAA')
        var input = "CGAUT";
        var expected = "GCUAA";
        var result = new string(input.Select(SequenceExtensions.GetRnaComplementBase).ToArray());
        Assert.That(result, Is.EqualTo(expected));
    }

    [Test]
    [Description("Cross-verify DNA complement: complement('CGAUT') → 'GCTAA' (Biopython: U treated as T)")]
    public void GetComplementBase_BiopythonCrossVerification_CgautDna()
    {
        // Biopython: complement(Seq("CGAUT")) → Seq('GCTAA')  (U→A, output T for DNA)
        // Our GetComplementBase: C→G, G→C, A→T, U→A, T→A
        // Expected: GCTAA
        var input = "CGAUT";
        var expected = "GCTAA";
        var result = new string(input.Select(SequenceExtensions.GetComplementBase).ToArray());
        Assert.That(result, Is.EqualTo(expected));
    }

    #endregion
}
