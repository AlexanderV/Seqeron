// PARSE-FASTA-001 — FASTA Parsing (opt-in multi-alphabet mode)
// Evidence: docs/Evidence/PARSE-FASTA-001-Evidence.md
// TestSpec: tests/TestSpecs/PARSE-FASTA-001.md
// Sources:
//   - NC-IUB (1985). "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences."
//     Nucleic Acids Research 13(9):3021–3030. IUPAC nucleotide codes: A C G T U R Y S W K M B D H V N + '-'.
//   - IUPAC amino-acid one-letter codes (bioinformatics.org IUPAC table): 20 standard + B Z J X U O + '*'.
//   - NCBI BLAST topics (https://blast.ncbi.nlm.nih.gov/doc/blast-topics/): FASTA accepts IUPAC
//     degenerate nucleotide codes and amino-acid codes (U=Sec, X=any, *=stop).
using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests.Unit.IO;

/// <summary>
/// Tests for the opt-in <see cref="SequenceAlphabet"/> FASTA parsing overloads.
/// The default DNA-only path is covered by <see cref="FastaParserTests"/>.
/// </summary>
[TestFixture]
public class FastaParser_Alphabet_Tests
{
    #region Parse(string, SequenceAlphabet) — RNA

    // A1 — RNA FASTA with U is preserved verbatim in RNA mode.
    // Evidence: RNA alphabet A C G U; lowercase mapped to uppercase (NCBI/Wikipedia FASTA).
    [Test]
    public void Parse_RnaSequenceWithUracil_PreservesSequence()
    {
        const string fasta = ">rna1 messenger\nAUGCAUGC\nGGUU";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.Rna).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(records, Has.Count.EqualTo(1), "single RNA record expected");
            Assert.That(records[0].Id, Is.EqualTo("rna1"), "id is first whitespace token");
            Assert.That(records[0].Description, Is.EqualTo("messenger"), "description is remainder of header");
            Assert.That(records[0].Sequence, Is.EqualTo("AUGCAUGCGGUU"),
                "multi-line RNA concatenated; U preserved (not coerced to T)");
            Assert.That(records[0].Alphabet, Is.EqualTo(SequenceAlphabet.Rna), "alphabet recorded");
        });
    }

    // A2 — lowercase RNA mapped to uppercase.
    [Test]
    public void Parse_LowercaseRna_NormalizedToUppercase()
    {
        const string fasta = ">r\naugc";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.Rna).ToList();

        Assert.That(records[0].Sequence, Is.EqualTo("AUGC"),
            "lowercase rna mapped to uppercase per FASTA convention");
    }

    // A3 — DNA letter T is out-of-alphabet for RNA → rejected.
    [Test]
    public void Parse_ThymineInRnaMode_Throws()
    {
        const string fasta = ">r\nAUGT";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta, SequenceAlphabet.Rna).ToList(),
            "T is not in the RNA alphabet (A C G U) and must be rejected");
    }

    #endregion

    #region Parse(string, SequenceAlphabet) — Protein

    // A4 — protein FASTA with W/Y/* and ambiguity X is preserved in protein mode.
    // Evidence: IUPAC amino-acid codes include W (Trp), Y (Tyr), X (Xaa), '*' (stop).
    [Test]
    public void Parse_ProteinSequenceWithStopAndAmbiguity_PreservesSequence()
    {
        const string fasta = ">prot1 sample protein\nMWYXBZJUO*";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.Protein).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(records, Has.Count.EqualTo(1), "single protein record expected");
            Assert.That(records[0].Id, Is.EqualTo("prot1"));
            Assert.That(records[0].Description, Is.EqualTo("sample protein"));
            Assert.That(records[0].Sequence, Is.EqualTo("MWYXBZJUO*"),
                "all IUPAC amino-acid + ambiguity (B Z J X) + rare (U O) + stop (*) preserved");
            Assert.That(records[0].Alphabet, Is.EqualTo(SequenceAlphabet.Protein));
        });
    }

    // A5 — out-of-alphabet character in protein mode is rejected.
    // Evidence: '@' is not an IUPAC amino-acid code.
    [Test]
    public void Parse_NonAminoAcidCharInProteinMode_Throws()
    {
        const string fasta = ">p\nMWY@";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta, SequenceAlphabet.Protein).ToList(),
            "'@' is not an IUPAC amino-acid code and must be rejected");
    }

    #endregion

    #region Parse(string, SequenceAlphabet) — IUPAC nucleotide

    // A6 — IUPAC-ambiguous DNA (N/R/Y and others) preserved in IUPAC mode.
    // Evidence: NC-IUB (1985) ambiguity codes R Y S W K M B D H V N + gap '-'.
    [Test]
    public void Parse_IupacAmbiguousNucleotides_PreservesSequence()
    {
        const string fasta = ">amb degenerate\nACGTNRYSWKMBDHV-U";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.IupacNucleotide).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(records, Has.Count.EqualTo(1));
            Assert.That(records[0].Id, Is.EqualTo("amb"));
            Assert.That(records[0].Sequence, Is.EqualTo("ACGTNRYSWKMBDHV-U"),
                "full IUPAC nucleotide code set incl. ambiguity, gap '-' and U preserved verbatim");
            Assert.That(records[0].Alphabet, Is.EqualTo(SequenceAlphabet.IupacNucleotide));
        });
    }

    // A7 — protein-only residue is out-of-alphabet for IUPAC nucleotide → rejected.
    // Evidence: 'E' (Glu) is not in the IUPAC nucleotide code set.
    [Test]
    public void Parse_ProteinResidueInIupacNucleotideMode_Throws()
    {
        const string fasta = ">n\nACGTE";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta, SequenceAlphabet.IupacNucleotide).ToList(),
            "'E' is not an IUPAC nucleotide code and must be rejected");
    }

    #endregion

    #region Default strict-DNA regression (documented behaviour preserved)

    // A8 — RNA (U) is REJECTED by the default DNA-only Parse (documented limitation kept as default).
    [Test]
    public void Parse_RnaSequence_DefaultDnaPath_Throws()
    {
        const string fasta = ">r\nAUGC";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta).ToList(),
            "default Parse is strict DNA-only: U must still be rejected");
    }

    // A9 — protein residues are REJECTED by the default DNA-only Parse.
    [Test]
    public void Parse_ProteinSequence_DefaultDnaPath_Throws()
    {
        const string fasta = ">p\nMWYK";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta).ToList(),
            "default Parse is strict DNA-only: protein residues must still be rejected");
    }

    // A10 — IUPAC ambiguity codes are REJECTED by the default DNA-only Parse.
    [Test]
    public void Parse_IupacAmbiguity_DefaultDnaPath_Throws()
    {
        const string fasta = ">a\nACGTNRY";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta).ToList(),
            "default Parse is strict DNA-only: N/R/Y must still be rejected");
    }

    // A11 — StrictDna alphabet mode matches default: rejects U.
    [Test]
    public void Parse_StrictDnaAlphabet_RejectsUracil()
    {
        const string fasta = ">r\nAUGC";

        Assert.Throws<ArgumentException>(() => FastaParser.Parse(fasta, SequenceAlphabet.StrictDna).ToList(),
            "StrictDna alphabet is A/C/G/T only; U rejected");
    }

    // A12 — StrictDna alphabet mode accepts plain DNA and preserves it.
    [Test]
    public void Parse_StrictDnaAlphabet_AcceptsDna()
    {
        const string fasta = ">d desc\nacgtACGT";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.StrictDna).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(records, Has.Count.EqualTo(1));
            Assert.That(records[0].Sequence, Is.EqualTo("ACGTACGT"),
                "plain DNA accepted and uppercased in StrictDna mode");
            Assert.That(records[0].Alphabet, Is.EqualTo(SequenceAlphabet.StrictDna));
        });
    }

    #endregion

    #region Edge cases

    // A13 — empty input returns no records in alphabet mode.
    [Test]
    public void Parse_EmptyInput_AlphabetMode_ReturnsEmpty()
    {
        var records = FastaParser.Parse("", SequenceAlphabet.Protein).ToList();

        Assert.That(records, Is.Empty, "empty input yields no records");
    }

    // A14 — multi-record protein FASTA: each record validated and preserved.
    [Test]
    public void Parse_MultiRecordProtein_ReturnsAll()
    {
        const string fasta = ">p1\nMWY\n>p2 second\nGAVL*";

        var records = FastaParser.Parse(fasta, SequenceAlphabet.Protein).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(records, Has.Count.EqualTo(2));
            Assert.That(records[0].Id, Is.EqualTo("p1"));
            Assert.That(records[0].Sequence, Is.EqualTo("MWY"));
            Assert.That(records[1].Id, Is.EqualTo("p2"));
            Assert.That(records[1].Description, Is.EqualTo("second"));
            Assert.That(records[1].Sequence, Is.EqualTo("GAVL*"));
        });
    }

    #endregion
}
