namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Canonical tests for sequence validation methods.
/// Test Unit: SEQ-VALID-001
/// Evidence: IUPAC-IUB 1970 standard, NC-IUB 1984, Wikipedia Nucleic acid notation,
/// Biopython Bio.Data.IUPACData (unambiguous_dna_letters, unambiguous_rna_letters)
/// </summary>
[TestFixture]
public class SequenceExtensions_SequenceValidation_Tests
{
    #region IsValidDna - MUST Tests

    [Test]
    [Description("M1: Empty sequence is valid (Biopython: zero-length sequences are always defined; vacuous truth)")]
    public void IsValidDna_EmptySequence_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = ReadOnlySpan<char>.Empty;
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M3: All standard DNA bases are valid (IUPAC 1970; Biopython unambiguous_dna_letters=GATC)")]
    public void IsValidDna_AllStandardBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "ACGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("M5: Lowercase DNA bases are valid (Wikipedia: lowercase used in sequence files for soft masking)")]
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
    [Description("M13: Ambiguity code N is invalid (NC-IUB 1984: N represents positional variant, not a nucleotide)")]
    public void IsValidDna_AmbiguityCodeN_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACGN".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.False);
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
    [Description("M2: Empty sequence is valid (Biopython: zero-length sequences are always defined; vacuous truth)")]
    public void IsValidRna_EmptySequence_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = ReadOnlySpan<char>.Empty;
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("M4: All standard RNA bases are valid (IUPAC 1970; Biopython unambiguous_rna_letters=GAUC)")]
    public void IsValidRna_AllStandardBases_ReturnsTrue()
    {
        ReadOnlySpan<char> sequence = "ACGU".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    [Test]
    [Description("M6: Lowercase RNA bases are valid (Wikipedia: lowercase used in sequence files for soft masking)")]
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
    [Description("RNA: Ambiguity code N is invalid (NC-IUB 1984: N represents positional variant, not a nucleotide)")]
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
    [Description("IUPAC ambiguity codes are invalid (NC-IUB 1984: represent positional variants, not nucleotides)")]
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
    [Description("C1: Unicode characters are invalid (non-ASCII)")]
    public void IsValidDna_UnicodeCharacter_ReturnsFalse()
    {
        ReadOnlySpan<char> sequence = "ACG\u65E5".AsSpan(); // ACG日
        Assert.That(sequence.IsValidDna(), Is.False);
    }

    [Test]
    [Description("Gap character (-) is invalid (not a nucleotide)")]
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
    [Description("INV-3: Case invariance — uppercase, lowercase, mixed all return true")]
    public void IsValidDna_CaseInvariance_AllCasesReturnTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That("ACGT".AsSpan().IsValidDna(), Is.True);
            Assert.That("acgt".AsSpan().IsValidDna(), Is.True);
            Assert.That("AcGt".AsSpan().IsValidDna(), Is.True);
        });
    }

    [Test]
    [Description("INV-4: RNA case invariance — uppercase, lowercase, mixed all return true")]
    public void IsValidRna_CaseInvariance_AllCasesReturnTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That("ACGU".AsSpan().IsValidRna(), Is.True);
            Assert.That("acgu".AsSpan().IsValidRna(), Is.True);
            Assert.That("AcGu".AsSpan().IsValidRna(), Is.True);
        });
    }

    #endregion

    #region IndexOfInvalid / constructor consistency (review 2026-09)

    [TestCase("", -1)]
    [TestCase("ACGTacgt", -1)]
    [TestCase("XACGT", 0)]
    [TestCase("ACGTX", 4)]
    [TestCase("ACGU", 3)]
    [TestCase("ACGN", 3)]
    [TestCase("AC GT", 2)]
    [Description("Canonical first-invalid index for DNA (unambiguous alphabet GATC, Biopython IUPACData)")]
    public void IndexOfInvalidDna_ReturnsFirstInvalidPosition(string sequence, int expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(sequence.AsSpan().IndexOfInvalidDna(), Is.EqualTo(expected));
            Assert.That(sequence.AsSpan().IsValidDna(), Is.EqualTo(expected < 0));
            Assert.That(DnaSequence.TryCreate(sequence, out _), Is.EqualTo(expected < 0), "INV-5: TryCreate ⇔ IsValidDna");
        });
    }

    [TestCase("", -1)]
    [TestCase("ACGUacgu", -1)]
    [TestCase("ACGT", 3)]
    [TestCase("ACGN", 3)]
    [TestCase("-ACG", 0)]
    [Description("Canonical first-invalid index for RNA (unambiguous alphabet GAUC, Biopython IUPACData)")]
    public void IndexOfInvalidRna_ReturnsFirstInvalidPosition(string sequence, int expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(sequence.AsSpan().IndexOfInvalidRna(), Is.EqualTo(expected));
            Assert.That(sequence.AsSpan().IsValidRna(), Is.EqualTo(expected < 0));
            Assert.That(RnaSequence.TryCreate(sequence, out _), Is.EqualTo(expected < 0));
        });
    }

    [Test]
    [Description("U+017F 'ſ' upper-cases (invariant) to 'S'; never a DNA/RNA base, and ctor error names the position")]
    public void IsValidDna_LongSNonAscii_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That("ACGſ".AsSpan().IndexOfInvalidDna(), Is.EqualTo(3));
            Assert.That("ACGſ".AsSpan().IsValidRna(), Is.False);
            var ex = Assert.Throws<ArgumentException>(() => new DnaSequence("ACGX"));
            Assert.That(ex!.Message, Does.Contain("position 3"));
        });
    }

    #endregion

    #region IsValidIupacDna / IsValidIupacRna (NC-IUB 1984; Biopython ambiguous_*_letters; scikit-bio)

    // Expected values computed with Biopython 1.88 IUPACData.ambiguous_dna_letters ("GATCRYWSMKHBVDN") /
    // ambiguous_rna_letters ("GAUCRYWSMKHBVDN") and confirmed with scikit-bio 0.7.4 DNA/RNA(lowercase=True)
    // (definite ∪ degenerate chars; gaps excluded).
    [TestCase("", true, true)]
    [TestCase("ACGTNRYSWKMBDHV", true, false)]
    [TestCase("acgtnryswkmbdhv", true, false)]
    [TestCase("ACGT", true, false)]
    [TestCase("ACGU", false, true)]
    [TestCase("ACGUNRYSWKMBDHV", false, true)]
    [TestCase("AC-GT", false, false)]
    [TestCase("AC.GT", false, false)]
    [TestCase("ACGX", false, false)]
    [TestCase("ACG N", false, false)]
    [TestCase("ACGTE", false, false)]
    [TestCase("ACGTI", false, false)]
    [TestCase("ſ", false, false)]
    [TestCase("ACGTſ", false, false)]
    public void IsValidIupac_MatchesBiopythonAndScikitBio(string sequence, bool expectedDna, bool expectedRna)
    {
        Assert.Multiple(() =>
        {
            Assert.That(sequence.AsSpan().IsValidIupacDna(), Is.EqualTo(expectedDna), "DNA");
            Assert.That(sequence.AsSpan().IsValidIupacRna(), Is.EqualTo(expectedRna), "RNA");
        });
    }

    [Test]
    [Description("IupacHelper.IsNucleotideCode accepts exactly the 15 codes of Biopython ambiguous_dna_letters")]
    public void IupacHelper_IsNucleotideCode_ExactlyFifteenCodes()
    {
        var accepted = Enumerable.Range(0, 128).Select(i => (char)i).Where(IupacHelper.IsNucleotideCode);
        Assert.That(string.Concat(accepted.OrderBy(c => c)), Is.EqualTo("ABCDGHKMNRSTVWY"));
    }

    [Test]
    [Description("IupacDnaSequence: U+017F 'ſ' must not be folded into IUPAC 'S' (scikit-bio DNA('ſ', lowercase=True) raises)")]
    public void IupacDnaSequence_LongSNonAscii_IsInvalid()
    {
        var seq = new IupacDnaSequence("ACſ");
        Assert.Multiple(() =>
        {
            Assert.That(seq.IsValid(), Is.False);
            Assert.That(seq.Sequence, Is.EqualTo("ACſ"));
            Assert.That(new IupacDnaSequence("acgtnry").Sequence, Is.EqualTo("ACGTNRY"));
        });
    }

    #endregion
}
