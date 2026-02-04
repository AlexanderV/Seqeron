using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Tests for FASTA format parsing and formatting.
    /// Test Unit: PARSE-FASTA-001
    /// Evidence: Wikipedia (FASTA format), NCBI BLAST Help
    /// </summary>
    [TestFixture]
    public class FastaParserTests
    {
        #region Parse - Basic Functionality

        /// <summary>
        /// M1: Single sequence with header and description parses correctly.
        /// Evidence: Wikipedia - "A sequence begins with a greater-than character followed by a description"
        /// </summary>
        [Test]
        public void Parse_SingleSequence_ReturnsCorrectEntry()
        {
            const string fasta = ">seq1 Test sequence\nACGTACGT\nACGTACGT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries, Has.Count.EqualTo(1));
                Assert.That(entries[0].Id, Is.EqualTo("seq1"));
                Assert.That(entries[0].Description, Is.EqualTo("Test sequence"));
                Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ACGTACGTACGTACGT"));
            });
        }

        /// <summary>
        /// M2: Multiple sequences in multi-FASTA format are all returned.
        /// Evidence: Wikipedia - "A multiple-sequence FASTA format would be obtained by concatenating several single-sequence FASTA files"
        /// </summary>
        [Test]
        public void Parse_MultiSequence_ReturnsAllEntries()
        {
            const string fasta = ">seq1 First\nAAAA\n>seq2 Second\nCCCC\n>seq3 Third\nGGGG";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries, Has.Count.EqualTo(3));
                Assert.That(entries[0].Id, Is.EqualTo("seq1"));
                Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("AAAA"));
                Assert.That(entries[1].Id, Is.EqualTo("seq2"));
                Assert.That(entries[1].Sequence.Sequence, Is.EqualTo("CCCC"));
                Assert.That(entries[2].Id, Is.EqualTo("seq3"));
                Assert.That(entries[2].Sequence.Sequence, Is.EqualTo("GGGG"));
            });
        }

        /// <summary>
        /// M3: Multi-line sequence is concatenated into single sequence.
        /// Evidence: Wikipedia - "interleaved, or on multiple lines"
        /// </summary>
        [Test]
        public void Parse_MultilineSequence_ConcatenatesLines()
        {
            const string fasta = ">gene1\nAAAA\nCCCC\nGGGG\nTTTT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("AAAACCCCGGGGTTTT"));
        }

        /// <summary>
        /// M4: Empty input returns empty enumerable.
        /// </summary>
        [Test]
        public void Parse_EmptyInput_ReturnsEmpty()
        {
            var entries = FastaParser.Parse("").ToList();

            Assert.That(entries, Is.Empty);
        }

        /// <summary>
        /// M4b: Whitespace-only input returns empty.
        /// </summary>
        [Test]
        public void Parse_WhitespaceOnlyInput_ReturnsEmpty()
        {
            var entries = FastaParser.Parse("   \n  \n   ").ToList();

            Assert.That(entries, Is.Empty);
        }

        #endregion

        #region Parse - Header Handling

        /// <summary>
        /// M5: Header with description separates ID and description correctly.
        /// Evidence: Wikipedia - "first word of header" is ID, rest is description
        /// </summary>
        [Test]
        public void Parse_HeaderWithDescription_ParsesBoth()
        {
            const string fasta = ">NM_001 Homo sapiens gene X (GeneX), mRNA\nATGC";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries[0].Id, Is.EqualTo("NM_001"));
                Assert.That(entries[0].Description, Is.EqualTo("Homo sapiens gene X (GeneX), mRNA"));
            });
        }

        /// <summary>
        /// M6: Header without description yields null description.
        /// Evidence: Wikipedia - description is optional
        /// </summary>
        [Test]
        public void Parse_HeaderWithoutDescription_ParsesIdOnly()
        {
            const string fasta = ">sequence_only\nACGT";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries[0].Id, Is.EqualTo("sequence_only"));
                Assert.That(entries[0].Description, Is.Null);
            });
        }

        /// <summary>
        /// M12: Special characters in header (pipes, colons) are preserved.
        /// Evidence: NCBI identifiers use pipes and colons (e.g., gi|12345|gb|AAA00000.1|)
        /// </summary>
        [Test]
        public void Parse_SpecialCharsInHeader_PreservesChars()
        {
            const string fasta = ">gi|12345|gb|AAA00000.1| hypothetical protein\nATGC";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries[0].Id, Is.EqualTo("gi|12345|gb|AAA00000.1|"));
                Assert.That(entries[0].Description, Is.EqualTo("hypothetical protein"));
            });
        }

        /// <summary>
        /// S6: Header property combines ID and description correctly.
        /// </summary>
        [Test]
        public void FastaEntry_Header_CombinesIdAndDescription()
        {
            var entryWithDesc = new FastaEntry("seq1", "A description", new DnaSequence("ATGC"));
            var entryNoDesc = new FastaEntry("seq2", null, new DnaSequence("ATGC"));

            Assert.Multiple(() =>
            {
                Assert.That(entryWithDesc.Header, Is.EqualTo("seq1 A description"));
                Assert.That(entryNoDesc.Header, Is.EqualTo("seq2"));
            });
        }

        /// <summary>
        /// S7: ToString returns useful debugging information.
        /// </summary>
        [Test]
        public void FastaEntry_ToString_ReturnsIdAndLength()
        {
            var entry = new FastaEntry("seq1", "Description", new DnaSequence("ATGCATGC"));

            var result = entry.ToString();

            Assert.That(result, Is.EqualTo("seq1 (8 bp)"));
        }

        #endregion

        #region Parse - Whitespace Handling

        /// <summary>
        /// M7: Whitespace within sequence lines is ignored.
        /// Evidence: Wikipedia - "invalid characters would be ignored (including spaces, tabulators)"
        /// </summary>
        [Test]
        public void Parse_WhitespaceInSequence_IgnoresWhitespace()
        {
            const string fasta = ">seq1\n  ATGC  \n  GGCC  ";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ATGCGGCC"));
        }

        /// <summary>
        /// M13: Leading and trailing whitespace in sequence lines is trimmed.
        /// Evidence: NCBI/Wikipedia - whitespace handling
        /// </summary>
        [Test]
        public void Parse_LeadingTrailingWhitespaceInSequence_Trimmed()
        {
            const string fasta = ">seq1\n   AAAA   \n   TTTT   ";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("AAAATTTT"));
        }

        /// <summary>
        /// S5: Blank lines in input are skipped.
        /// Evidence: Common real-world data contains extra blank lines
        /// </summary>
        [Test]
        public void Parse_BlankLinesInInput_SkipsBlankLines()
        {
            const string fasta = ">seq1\nATGC\n\n\n>seq2\nGGCC\n\n";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries, Has.Count.EqualTo(2));
                Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ATGC"));
                Assert.That(entries[1].Sequence.Sequence, Is.EqualTo("GGCC"));
            });
        }

        /// <summary>
        /// Header only without sequence is not yielded.
        /// ASSUMPTION: Implementation-specific behavior - header without sequence is skipped
        /// </summary>
        [Test]
        public void Parse_HeaderWithoutSequence_NotYielded()
        {
            const string fasta = ">empty_seq\n>has_seq\nATGC";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries, Has.Count.EqualTo(1));
                Assert.That(entries[0].Id, Is.EqualTo("has_seq"));
            });
        }

        #endregion

        #region ToFasta - Output Formatting

        /// <summary>
        /// M8: Single entry formats correctly with header and sequence.
        /// </summary>
        [Test]
        public void ToFasta_SingleEntry_FormatsCorrectly()
        {
            var entry = new FastaEntry("seq1", "Test description", new DnaSequence("ACGTACGT"));

            string fasta = FastaParser.ToFasta(new[] { entry }, lineWidth: 80);

            Assert.Multiple(() =>
            {
                Assert.That(fasta, Does.StartWith(">seq1 Test description"));
                Assert.That(fasta, Does.Contain("ACGTACGT"));
            });
        }

        /// <summary>
        /// M9: Long sequence is wrapped at specified line width.
        /// Evidence: Wikipedia - "typically no more than 80 characters in length"
        /// </summary>
        [Test]
        public void ToFasta_LongSequence_WrapsAtLineWidth()
        {
            var seq = new DnaSequence(new string('A', 100));
            var entry = new FastaEntry("long", null, seq);

            string fasta = FastaParser.ToFasta(new[] { entry }, lineWidth: 80);
            var lines = fasta.Split('\n')
                .Where(l => !l.StartsWith(">") && !string.IsNullOrWhiteSpace(l))
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(lines, Has.Count.EqualTo(2));
                Assert.That(lines[0].Trim().Length, Is.EqualTo(80));
                Assert.That(lines[1].Trim().Length, Is.EqualTo(20));
            });
        }

        /// <summary>
        /// M10: Entry without description omits description from header.
        /// </summary>
        [Test]
        public void ToFasta_NoDescription_OmitsDescription()
        {
            var entry = new FastaEntry("seq1", null, new DnaSequence("ATGC"));

            string fasta = FastaParser.ToFasta(new[] { entry });
            var lines = fasta.Split('\n').Select(l => l.TrimEnd('\r')).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(lines[0], Is.EqualTo(">seq1"));
                Assert.That(lines[0], Does.Not.Contain(" "));
            });
        }

        /// <summary>
        /// S1: Custom line width wraps correctly.
        /// </summary>
        [Test]
        public void ToFasta_CustomLineWidth_WrapsCorrectly()
        {
            var seq = new DnaSequence(new string('G', 50));
            var entry = new FastaEntry("test", null, seq);

            string fasta = FastaParser.ToFasta(new[] { entry }, lineWidth: 20);
            var seqLines = fasta.Split('\n')
                .Select(l => l.TrimEnd('\r'))
                .Where(l => !l.StartsWith(">") && !string.IsNullOrWhiteSpace(l))
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(seqLines, Has.Count.EqualTo(3)); // 20 + 20 + 10
                Assert.That(seqLines[0].Length, Is.EqualTo(20));
                Assert.That(seqLines[1].Length, Is.EqualTo(20));
                Assert.That(seqLines[2].Length, Is.EqualTo(10));
            });
        }

        /// <summary>
        /// S2: Multiple entries format as valid multi-FASTA.
        /// </summary>
        [Test]
        public void ToFasta_MultipleEntries_FormatsAll()
        {
            var entries = new[]
            {
                new FastaEntry("seq1", "First", new DnaSequence("AAAA")),
                new FastaEntry("seq2", "Second", new DnaSequence("CCCC"))
            };

            string fasta = FastaParser.ToFasta(entries);

            Assert.Multiple(() =>
            {
                Assert.That(fasta, Does.Contain(">seq1 First"));
                Assert.That(fasta, Does.Contain("AAAA"));
                Assert.That(fasta, Does.Contain(">seq2 Second"));
                Assert.That(fasta, Does.Contain("CCCC"));
            });
        }

        #endregion

        #region Round Trip Integrity

        /// <summary>
        /// M11: Parse → Format → Parse preserves all data.
        /// </summary>
        [Test]
        public void RoundTrip_ParseAndFormat_PreservesData()
        {
            const string original = ">seq1 Description here\nACGTACGTACGTACGT\n>seq2 Another one\nTTTTCCCCGGGGAAAA";

            var entries = FastaParser.Parse(original).ToList();
            string written = FastaParser.ToFasta(entries);
            var reparsed = FastaParser.Parse(written).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(reparsed, Has.Count.EqualTo(entries.Count));
                Assert.That(reparsed[0].Id, Is.EqualTo(entries[0].Id));
                Assert.That(reparsed[0].Description, Is.EqualTo(entries[0].Description));
                Assert.That(reparsed[0].Sequence.Sequence, Is.EqualTo(entries[0].Sequence.Sequence));
                Assert.That(reparsed[1].Id, Is.EqualTo(entries[1].Id));
                Assert.That(reparsed[1].Description, Is.EqualTo(entries[1].Description));
                Assert.That(reparsed[1].Sequence.Sequence, Is.EqualTo(entries[1].Sequence.Sequence));
            });
        }

        /// <summary>
        /// Round trip with special characters in header.
        /// </summary>
        [Test]
        public void RoundTrip_SpecialCharsInHeader_PreservesData()
        {
            const string original = ">gi|12345|ref|NM_001.1| Homo sapiens\nATGCATGC";

            var entries = FastaParser.Parse(original).ToList();
            string written = FastaParser.ToFasta(entries);
            var reparsed = FastaParser.Parse(written).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(reparsed[0].Id, Is.EqualTo(entries[0].Id));
                Assert.That(reparsed[0].Description, Is.EqualTo(entries[0].Description));
            });
        }

        #endregion

        #region File I/O

        /// <summary>
        /// S3: ParseFile reads FASTA file correctly.
        /// </summary>
        [Test]
        public void ParseFile_ValidFile_ReturnsEntries()
        {
            var tempFile = Path.GetTempFileName();
            try
            {
                File.WriteAllText(tempFile, ">seq1 Test\nATGCATGC\n>seq2\nGGGGCCCC");

                var entries = FastaParser.ParseFile(tempFile).ToList();

                Assert.Multiple(() =>
                {
                    Assert.That(entries, Has.Count.EqualTo(2));
                    Assert.That(entries[0].Id, Is.EqualTo("seq1"));
                    Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ATGCATGC"));
                    Assert.That(entries[1].Id, Is.EqualTo("seq2"));
                    Assert.That(entries[1].Sequence.Sequence, Is.EqualTo("GGGGCCCC"));
                });
            }
            finally
            {
                if (File.Exists(tempFile)) File.Delete(tempFile);
            }
        }

        /// <summary>
        /// S4: WriteFile creates valid FASTA file.
        /// </summary>
        [Test]
        public void WriteFile_ValidEntries_CreatesFile()
        {
            var tempFile = Path.GetTempFileName();
            try
            {
                var entries = new[]
                {
                    new FastaEntry("seq1", "Test", new DnaSequence("ATGCATGC"))
                };

                FastaParser.WriteFile(tempFile, entries, lineWidth: 80);

                var content = File.ReadAllText(tempFile);
                Assert.Multiple(() =>
                {
                    Assert.That(content, Does.StartWith(">seq1 Test"));
                    Assert.That(content, Does.Contain("ATGCATGC"));
                });
            }
            finally
            {
                if (File.Exists(tempFile)) File.Delete(tempFile);
            }
        }

        /// <summary>
        /// C1: ParseFileAsync reads file asynchronously.
        /// </summary>
        [Test]
        public async System.Threading.Tasks.Task ParseFileAsync_ValidFile_ReturnsEntries()
        {
            var tempFile = Path.GetTempFileName();
            try
            {
                await File.WriteAllTextAsync(tempFile, ">seq1 Async test\nATGCGGCC");

                var entries = new System.Collections.Generic.List<FastaEntry>();
                await foreach (var entry in FastaParser.ParseFileAsync(tempFile))
                {
                    entries.Add(entry);
                }

                Assert.Multiple(() =>
                {
                    Assert.That(entries, Has.Count.EqualTo(1));
                    Assert.That(entries[0].Id, Is.EqualTo("seq1"));
                    Assert.That(entries[0].Sequence.Sequence, Is.EqualTo("ATGCGGCC"));
                });
            }
            finally
            {
                if (File.Exists(tempFile)) File.Delete(tempFile);
            }
        }

        #endregion

        #region Edge Cases

        /// <summary>
        /// C2: Very long sequence is handled correctly.
        /// </summary>
        [Test]
        public void Parse_VeryLongSequence_HandlesCorrectly()
        {
            var longSeq = new string('A', 10000);
            var fasta = $">long_seq\n{longSeq}";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.That(entries[0].Sequence.Sequence.Length, Is.EqualTo(10000));
        }

        /// <summary>
        /// Tab delimiter in header separates ID and description.
        /// </summary>
        [Test]
        public void Parse_TabInHeader_SeparatesIdAndDescription()
        {
            const string fasta = ">seq1\tDescription with tab\nATGC";

            var entries = FastaParser.Parse(fasta).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(entries[0].Id, Is.EqualTo("seq1"));
                Assert.That(entries[0].Description, Is.EqualTo("Description with tab"));
            });
        }

        /// <summary>
        /// Sequence exactly at line width boundary wraps correctly.
        /// </summary>
        [Test]
        public void ToFasta_SequenceExactlyLineWidth_NoExtraLine()
        {
            var seq = new DnaSequence(new string('T', 80));
            var entry = new FastaEntry("exact", null, seq);

            string fasta = FastaParser.ToFasta(new[] { entry }, lineWidth: 80);
            var seqLines = fasta.Split('\n')
                .Select(l => l.TrimEnd('\r'))
                .Where(l => !l.StartsWith(">") && !string.IsNullOrWhiteSpace(l))
                .ToList();

            Assert.That(seqLines, Has.Count.EqualTo(1));
            Assert.That(seqLines[0].Length, Is.EqualTo(80));
        }

        #endregion
    }
}
