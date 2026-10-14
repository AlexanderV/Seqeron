// TRANS-CODON-001 / TRANS-PROT-001 — MCP translate_dna / translate_rna: optional geneticCode
// (NCBI table) + IUPAC ambiguity codes + opt-in nucleotide X (FIN-B02 WP5, A1-18).
// Reference: Biopython 1.88 Bio.Seq.Seq(s).translate(table=n, to_stop=b); the expected strings
// below are Biopython 1.88 outputs for exactly these inputs.
using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class TranslateGeneticCodeTests
{
    [TestCase("ATGGCNTAR", 1, false, "MA*")]
    [TestCase("ATGGCNTAR", 1, true, "MA")]
    [TestCase("ATGAGATGA", 1, false, "MR*")]
    [TestCase("ATGAGATGA", 2, false, "M*W")]          // table 2: AGA -> stop, TGA -> Trp
    [TestCase("ATAAGATGAGGG", 2, false, "M*WG")]      // table 2: ATA -> Met
    [TestCase("ATAAGATGAGGG", 2, true, "M")]
    [TestCase("atgrayytgsnn", 1, false, "MBLX")]      // RAY -> B (Asx), YTG -> L
    [TestCase("MGNTGACTG", 11, false, "X*L")]
    [TestCase("TTAMGNGAY", 1, false, "LXD")]
    [TestCase("ATGTGA", 4, false, "MW")]
    [TestCase("CTGTAG", 3, false, "T*")]              // table 3: CTG -> Thr
    [TestCase("ATGTAATGA", 6, false, "MQ*")]          // table 6: TAA -> Gln
    [TestCase("ACGNATGAA", 1, false, "TXE")]
    public void TranslateDna_GeneticCodeAndIupac_MatchesBiopython(string seq, int table, bool toStop, string expected)
    {
        Assert.That(SequenceTools.TranslateDna(seq, 0, toStop, table).Protein, Is.EqualTo(expected));
    }

    [TestCase("AUGGCNUAR", 1, false, "MA*")]
    [TestCase("AUGAGAUGA", 2, false, "M*W")]
    [TestCase("AUGAGAUGA", 1, false, "MR*")]
    public void TranslateRna_GeneticCodeAndIupac_MatchesBiopython(string seq, int table, bool toStop, string expected)
    {
        Assert.That(SequenceTools.TranslateRna(seq, 0, toStop, table).Protein, Is.EqualTo(expected));
    }

    [Test]
    public void Translate_DefaultTable_IsStandard()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.TranslateDna("ATGAGATGA").Protein, Is.EqualTo("MR*"));
            Assert.That(SequenceTools.TranslateRna("AUGAGAUGA").Protein, Is.EqualTo("MR*"));
        });
    }

    // Biopython 1.88: Seq("GCXATGGGX").translate() == "AMG"; Seq("ATGTAX").translate() raises
    // TranslationError (X in a possible-stop codon). X is rejected unless allowNucleotideX.
    [Test]
    public void Translate_NucleotideX_OptIn()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.TranslateDna("GCXATGGGX", allowNucleotideX: true).Protein, Is.EqualTo("AMG"));
            Assert.That(SequenceTools.TranslateRna("GCXAUGGGX", allowNucleotideX: true).Protein, Is.EqualTo("AMG"));
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("GCXATGGGX"));
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGTAX", allowNucleotideX: true));
        });
    }

    [Test]
    public void Translate_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            var ex = Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGZZZ"));
            Assert.That(ex!.Message, Does.Contain("'Z' at position 3"));
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("AUGCGA"));     // U is RNA
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("ATGCGA"));     // T is DNA
            // An invalid character outside complete codons is still rejected (trailing partial codon).
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGCGAZ"));
            // Non-ASCII letter: no Unicode folding ('ſ' is not 'S').
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGſCA"));
            // Unknown NCBI table numbers (7, 8, 17-20 unassigned).
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATG", geneticCode: 7));
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("AUG", geneticCode: 0));
            // to_stop with dual-coding stop tables (Biopython ValueError).
            Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGTGA", 0, true, 27));
        });
    }
}
