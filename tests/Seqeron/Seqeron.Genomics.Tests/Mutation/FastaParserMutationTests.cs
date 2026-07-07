namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for FastaParser.cs (checklist 04 row 64, PARSE-FASTA-001).
/// The canonical suite left the async file reader and the FASTA line-wrapping writer
/// under-pinned. These exercise ParseFileAsync (record flush guards) and ToFasta with an
/// exact wrapped output.
/// </summary>
[TestFixture]
public class FastaParserMutationTests
{
    // ── ToFasta: line wrapping must not emit a trailing empty line ────────────────────

    [Test]
    public void ToFasta_WrapsAtExactMultiple_NoTrailingBlankLine()
    {
        // 8 bp at width 4 → exactly two full lines; the loop bound (i < length, not <=) must
        // stop before emitting a third empty line.
        var entries = new[] { new FastaEntry("seq1", null, new DnaSequence("ACGTACGT")) };

        string fasta = FastaParser.ToFasta(entries, lineWidth: 4);

        fasta.Should().Be(">seq1\nACGT\nACGT\n");
    }

    [Test]
    public void ToFasta_HeaderWithDescription_RoundTrips()
    {
        var entries = new[] { new FastaEntry("seq1", "desc here", new DnaSequence("ACGTAC")) };
        FastaParser.ToFasta(entries, lineWidth: 80).Should().Be(">seq1 desc here\nACGTAC\n");
    }

    // ── ParseFileAsync: record flush guards (header != null AND sequence non-empty) ───

    [Test]
    public async Task ParseFileAsync_TwoRecords_ParsedExactly()
    {
        string path = Path.GetTempFileName();
        try
        {
            await File.WriteAllTextAsync(path, ">a desc1\nACGT\nACGT\n>b\nGGGG\n");

            var entries = new List<FastaEntry>();
            await foreach (var e in FastaParser.ParseFileAsync(path))
                entries.Add(e);

            entries.Should().HaveCount(2);
            entries[0].Id.Should().Be("a");
            entries[0].Description.Should().Be("desc1");
            entries[0].Sequence.Sequence.Should().Be("ACGTACGT");
            entries[1].Id.Should().Be("b");
            entries[1].Sequence.Sequence.Should().Be("GGGG");
        }
        finally { File.Delete(path); }
    }

    [Test]
    public async Task ParseFileAsync_HeaderWithNoSequence_IsSkipped()
    {
        // ">b" has no sequence before ">c" → the (header != null && Length > 0) guard must NOT
        // flush an empty record (kills the `&&`→`||` and Length>0 boundary mutants).
        string path = Path.GetTempFileName();
        try
        {
            await File.WriteAllTextAsync(path, ">a\nACGT\n>b\n>c\nGGGG\n");

            var entries = new List<FastaEntry>();
            await foreach (var e in FastaParser.ParseFileAsync(path))
                entries.Add(e);

            entries.Select(e => e.Id).Should().Equal(new[] { "a", "c" }, "the empty 'b' record is not emitted");
            entries[0].Sequence.Sequence.Should().Be("ACGT");
            entries[1].Sequence.Sequence.Should().Be("GGGG");
        }
        finally { File.Delete(path); }
    }

    // ── Typed (alphabet-aware) ParseFileAsync(path, alphabet): the newer overload was uncovered ──

    [Test]
    public async Task ParseFileAsyncTyped_TwoRecordsAndEmptyHeaderSkipped()
    {
        string path = Path.GetTempFileName();
        try
        {
            // ">b" has no sequence before ">c": the (header != null && Length > 0) guard must NOT flush it.
            await File.WriteAllTextAsync(path, ">a desc1\nACGT\nACGT\n>b\n>c\nGGGG\n");

            var records = new List<FastaRecord>();
            await foreach (var r in FastaParser.ParseFileAsync(path, SequenceAlphabet.StrictDna))
                records.Add(r);

            records.Select(r => r.Id).Should().Equal(new[] { "a", "c" }, "the empty 'b' record is skipped");
            records[0].Description.Should().Be("desc1");
            records[0].Sequence.Should().Be("ACGTACGT");
            records[0].Alphabet.Should().Be(SequenceAlphabet.StrictDna);
            records[1].Sequence.Should().Be("GGGG");
        }
        finally { File.Delete(path); }
    }

    // ── Sync typed Parse(content, alphabet): empty-header record must be skipped (ParseReaderTyped flush) ──

    [Test]
    public void ParseTyped_EmptyHeaderRecord_IsSkipped()
    {
        var records = FastaParser.Parse(">a\nACGT\n>b\n>c\nGGGG\n", SequenceAlphabet.StrictDna).ToList();

        records.Select(r => r.Id).Should().Equal(new[] { "a", "c" }, "the empty 'b' record is not emitted");
        records[1].Sequence.Should().Be("GGGG");
    }

    // ── Sync untyped Parse(content): empty-header record must be skipped (ParseReader flush guard) ──

    [Test]
    public void ParseUntyped_EmptyHeaderRecord_IsSkipped()
    {
        var entries = FastaParser.Parse(">a\nACGT\n>b\n>c\nGGGG\n").ToList();

        entries.Select(e => e.Id).Should().Equal(new[] { "a", "c" }, "the empty 'b' record is not emitted");
        entries[1].Sequence.Sequence.Should().Be("GGGG");
    }

    // ── FastaRecord.Header: present vs absent description (kills the Description==null branch) ──

    [Test]
    public void FastaRecordHeader_ReflectsDescriptionPresence()
    {
        var withDesc = new FastaRecord("id1", "some description", "ACGT", SequenceAlphabet.StrictDna);
        var noDesc = new FastaRecord("id2", null, "ACGT", SequenceAlphabet.StrictDna);

        withDesc.Header.Should().Be("id1 some description", "header includes the description when present");
        noDesc.Header.Should().Be("id2", "header is just the id when there is no description");
    }
}
