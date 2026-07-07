namespace Seqeron.Genomics.Tests.Algebraic;

using SequenceRecord = SequenceIO.SequenceRecord;
using SequenceFeature = SequenceIO.SequenceFeature;
using Reference = SequenceIO.Reference;

/// <summary>
/// Algebraic-law tests for the FileIO area (GenBank and EMBL flat-file I/O).
///
/// Algebraic testing pins the parse∘serialize round-trip of the GenBank
/// serializer and the validity of a minimal record. GenBank and EMBL both parse
/// into the shared, format-agnostic <see cref="SequenceRecord"/> model, so an
/// EMBL-sourced record can be round-tripped through the GenBank serializer — the
/// only flat-file writer in the library.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 69, 70.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("FileIO")]
public class FlatFileAlgebraicTests
{
    private static SequenceRecord Record(string id, string accession, string description,
        string organism, string sequence) =>
        new(
            Id: id,
            Accession: accession,
            Description: description,
            Sequence: sequence,
            Organism: organism,
            Taxonomy: null,
            Date: null,
            Features: Array.Empty<SequenceFeature>(),
            References: Array.Empty<Reference>(),
            Metadata: new Dictionary<string, string>());

    private static void AssertCoreFieldsEqual(SequenceRecord expected, SequenceRecord actual)
    {
        actual.Id.Should().Be(expected.Id);
        actual.Accession.Should().Be(expected.Accession);
        actual.Description.Should().Be(expected.Description);
        actual.Organism.Should().Be(expected.Organism);
        actual.Sequence.Should().Be(expected.Sequence);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-GENBANK-001 — GenBank flat-file I/O (FileIO)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 69.
    //
    // Model: ToGenBank serializes a SequenceRecord to GenBank flat-file text and
    //        ParseGenBankString parses it back, preserving the record's identity,
    //        accession, definition, organism and sequence.
    //   — docs/algorithms/FileIO; SequenceIO.ToGenBank / ParseGenBankString.
    //
    // Laws (row 69): RT — ParseGenBank(ToGenBank(record)) = record (core fields).
    //                ID — a minimal record (LOCUS + ORIGIN) parses to a valid record.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>RT: a record survives the GenBank serialize→parse cycle.</summary>
    [Test]
    public void GenBank_RoundTrip_ParseOfSerializeIsIdentity()
    {
        var record = Record("TESTSEQ", "ABC12345", "a sample genbank record",
            "Escherichia coli", "ACGTACGTACGTACGTACGT");
        var parsed = SequenceIO.ParseGenBankString(SequenceIO.ToGenBank(record)).Single();
        AssertCoreFieldsEqual(record, parsed);
    }

    /// <summary>ID: a minimal record (LOCUS + ORIGIN only) parses to a valid record.</summary>
    [Test]
    public void GenBank_Identity_MinimalRecordIsValid()
    {
        string minimal =
            "LOCUS       MINID              4 bp    DNA     linear   UNK 01-JAN-2020\n" +
            "ORIGIN\n" +
            "        1 acgt\n" +
            "//\n";
        var parsed = SequenceIO.ParseGenBankString(minimal).Single();
        parsed.Id.Should().Be("MINID");
        parsed.Sequence.Should().Be("ACGT");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-EMBL-001 — EMBL flat-file I/O (FileIO)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 70.
    //
    // Model: EMBL parses into the same SequenceRecord model as GenBank. The library
    //        provides no EMBL serializer, so the round-trip is exercised through the
    //        shared model and the GenBank writer: an EMBL-parsed record, serialized
    //        as GenBank and re-parsed, reproduces its core fields.
    //   — docs/algorithms/FileIO; SequenceIO.ParseEmblString (+ ToGenBank).
    //
    // Laws (row 70): RT — ParseGenBank(ToGenBank(ParseEmbl(text))) preserves the
    //                record (the SequenceRecord is serializable and recoverable).
    //                ID — a minimal EMBL record parses to a valid record.
    // ═══════════════════════════════════════════════════════════════════════

    private const string MinimalEmbl =
        "ID   EMBLSEQ; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\n" +
        "AC   XY98765;\n" +
        "DE   a sample embl record\n" +
        "OS   Escherichia coli\n" +
        "SQ   Sequence 20 BP;\n" +
        "     acgtacgtac gtacgtacgt                                              20\n" +
        "//\n";

    /// <summary>ID: a minimal EMBL record parses to a valid, populated record.</summary>
    [Test]
    public void Embl_Identity_MinimalRecordIsValid()
    {
        var parsed = SequenceIO.ParseEmblString(MinimalEmbl).Single();
        parsed.Id.Should().Be("EMBLSEQ");
        parsed.Accession.Should().Be("XY98765");
        parsed.Description.Should().Be("a sample embl record");
        parsed.Organism.Should().Be("Escherichia coli");
        parsed.Sequence.Should().Be("ACGTACGTACGTACGTACGT");
    }

    /// <summary>
    /// RT: the EMBL-parsed record is serializable and recoverable — round-tripping
    /// it through the GenBank writer/parser preserves its core fields.
    /// </summary>
    [Test]
    public void Embl_RoundTrip_ThroughSharedModelIsIdentity()
    {
        var fromEmbl = SequenceIO.ParseEmblString(MinimalEmbl).Single();
        var roundTripped = SequenceIO.ParseGenBankString(SequenceIO.ToGenBank(fromEmbl)).Single();
        AssertCoreFieldsEqual(fromEmbl, roundTripped);
    }
}
