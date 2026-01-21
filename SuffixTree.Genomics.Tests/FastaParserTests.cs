using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
{
    [TestFixture]
    public class FastaParserTests
    {
        #region Parsing

        [Test]
        public void Parse_SingleSequence_ParsesCorrectly()
        {
            string fasta = @">seq1 Test sequence
ACGTACGT
ACGTACGT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries, Has.Count.EqualTo(1));
            Assert.That(entries[0].Id, Is.EqualTo("seq1"));
            Assert.That(entries[0].Description, Is.EqualTo("Test sequence"));
            Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ACGTACGTACGTACGT"));
        }

        [Test]
        public void Parse_MultipleSequences_ParsesAll()
        {
            string fasta = @">seq1 First
AAAA
>seq2 Second
CCCC
>seq3 Third
GGGG";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries, Has.Count.EqualTo(3));
            Assert.That(entries[0].Id, Is.EqualTo("seq1"));
            Assert.That(entries[1].Id, Is.EqualTo("seq2"));
            Assert.That(entries[2].Id, Is.EqualTo("seq3"));
        }

        [Test]
        public void Parse_NoDescription_ParsesIdOnly()
        {
            string fasta = @">sequence_only
ACGT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Id, Is.EqualTo("sequence_only"));
            Assert.That(entries[0].Description, Is.Null);
        }

        [Test]
        public void Parse_EmptyInput_ReturnsEmpty()
        {
            var entries = FastaParser.Parse("").ToList();
            Assert.That(entries, Is.Empty);
        }

        [Test]
        public void Parse_MultilineSequence_ConcatenatesLines()
        {
            string fasta = @">gene1
AAAA
CCCC
GGGG
TTTT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("AAAACCCCGGGGTTTT"));
        }

        #endregion

        #region Writing

        [Test]
        public void ToFasta_SingleEntry_FormatsCorrectly()
        {
            var entry = new FastaEntry("seq1", "Test", new DnaSequence("ACGTACGT"));
            var entries = new[] { entry };

            string fasta = FastaParser.ToFasta(entries, lineWidth: 80);

            Assert.That(fasta, Does.StartWith(">seq1 Test"));
            Assert.That(fasta, Does.Contain("ACGTACGT"));
        }

        [Test]
        public void ToFasta_LongSequence_WrapsLines()
        {
            var seq = new DnaSequence(new string('A', 100));
            var entry = new FastaEntry("long", null, seq);
            var entries = new[] { entry };

            string fasta = FastaParser.ToFasta(entries, lineWidth: 80);
            var lines = fasta.Split('\n').Where(l => !l.StartsWith(">") && !string.IsNullOrWhiteSpace(l)).ToList();

            Assert.That(lines[0].Trim().Length, Is.EqualTo(80));
            Assert.That(lines[1].Trim().Length, Is.EqualTo(20));
        }

        #endregion

        #region Round Trip

        [Test]
        public void RoundTrip_ParseAndWrite_PreservesData()
        {
            string original = @">seq1 Description here
ACGTACGTACGTACGT
>seq2 Another one
TTTTCCCCGGGGAAAA";

            var entries = FastaParser.Parse(original).ToList();
            string written = FastaParser.ToFasta(entries);
            var reparsed = FastaParser.Parse(written).ToList();

            Assert.That(reparsed, Has.Count.EqualTo(2));
            Assert.That(reparsed[0].Sequence.Sequence, Is.EqualTo(entries[0].Sequence.Sequence));
            Assert.That(reparsed[1].Sequence.Sequence, Is.EqualTo(entries[1].Sequence.Sequence));
        }

        #endregion
    }
}
