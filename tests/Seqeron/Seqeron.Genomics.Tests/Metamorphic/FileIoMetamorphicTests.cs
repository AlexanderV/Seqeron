using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the FileIO area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the FORMAT DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PARSE-FASTA-001 — FASTA parsing / serialization (FileIO).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 64.
///
/// API under test (FastaParser.Parse / FastaParser.ToFasta):
///   FASTA is a sequence of records; each record is a '>'-prefixed header line (split into
///   ID = first token and Description = remainder) followed by sequence lines. Parsing
///   concatenates all non-whitespace characters of the sequence lines and ignores blank
///   lines; serialization writes '>'+Header then the sequence hard-wrapped at lineWidth.
///
/// Relations (derived from the format, NOT from output):
///   • COMP (write→parse→write = write): serialization has a canonical normal form, so
///          re-serialising the parse of a serialised record set reproduces the exact text;
///          parsing also round-trips the (ID, Description, Sequence) of every record.
///   • INV  (blank lines ⇒ same records): blank/whitespace-only lines contribute no sequence
///          characters and are not headers, so inserting them anywhere leaves the parsed
///          records unchanged (sequence-line wrapping is likewise irrelevant).
///   • INV  (trailing newline irrelevant): the final record is emitted whether or not the
///          input ends in a newline, so appending or removing trailing newlines is a no-op.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class FileIoMetamorphicTests
{
    #region PARSE-FASTA-001 — Helpers

    /// <summary>Run-order-independent identity of a FASTA record.</summary>
    private static List<(string Id, string? Description, string Sequence)> FastaRecords(string content) =>
        FastaParser.Parse(content)
            .Select(e => (e.Id, e.Description, e.Sequence.Sequence))
            .ToList();

    private static FastaEntry[] SampleEntries() => new[]
    {
        new FastaEntry("seq1", "first sequence", new DnaSequence("ACGTACGTACGT")),
        new FastaEntry("seq2", null, new DnaSequence("TTTTGGGGCCCCAAAA")),
        new FastaEntry("seq3", "chr1 region", new DnaSequence("GGGGCCCCTTTTAAAA")),
    };

    #endregion

    #region PARSE-FASTA-001 COMP — write→parse→write reproduces the text and the records

    [Test]
    [Description("COMP: parsing a serialised record set recovers every record's (ID, Description, Sequence), and re-serialising reproduces the exact text — serialization has a canonical normal form.")]
    public void Fasta_WriteParseWrite_IsIdempotentAndRecordPreserving()
    {
        var entries = SampleEntries();

        string written = FastaParser.ToFasta(entries);
        var parsed = FastaParser.Parse(written).ToList();

        // parse(write(x)) recovers each record exactly.
        parsed.Select(e => (e.Id, e.Description, e.Sequence.Sequence))
            .Should().Equal(entries.Select(e => (e.Id, e.Description, e.Sequence.Sequence)),
                because: "serialising then parsing must round-trip the ID, description and sequence of every record");

        // write(parse(write(x))) == write(x): the serialised form is a fixed point.
        FastaParser.ToFasta(parsed).Should().Be(written,
            because: "ToFasta produces a canonical normal form, so re-serialising the parsed records yields identical text");
    }

    #endregion

    #region PARSE-FASTA-001 INV — blank lines and line-wrapping do not change records

    [Test]
    [Description("INV: blank/whitespace-only lines carry no sequence characters and are not headers, so inserting them (and re-wrapping the sequence) leaves the parsed records unchanged.")]
    public void Fasta_BlankLinesAndRewrapping_PreserveRecords()
    {
        const string compact =
            ">seq1 first sequence\nACGTACGTACGT\n" +
            ">seq2\nTTTTGGGGCCCCAAAA\n";

        // Same records, but with blank lines scattered throughout and the first sequence
        // split across several lines (FASTA permits arbitrary sequence-line wrapping).
        const string withBlanks =
            "\n\n>seq1 first sequence\n\nACGT\nACGT\n\nACGT\n\n" +
            ">seq2\n\n\nTTTTGGGG\nCCCCAAAA\n\n";

        FastaRecords(withBlanks).Should().Equal(FastaRecords(compact),
            because: "blank lines add no sequence characters and sequence-line wrapping is not significant in FASTA");
    }

    #endregion

    #region PARSE-FASTA-001 INV — trailing newlines are irrelevant

    [Test]
    [Description("INV: the final record is emitted regardless of trailing newlines, so adding or removing them does not change the parsed records.")]
    public void Fasta_TrailingNewlines_AreIrrelevant()
    {
        const string noTrailing =
            ">seq1 first sequence\nACGTACGTACGT\n" +
            ">seq2\nTTTTGGGGCCCCAAAA";

        var baseline = FastaRecords(noTrailing);
        baseline.Should().HaveCount(2, because: "the last record must be emitted even without a trailing newline");

        FastaRecords(noTrailing + "\n").Should().Equal(baseline, because: "one trailing newline is a no-op");
        FastaRecords(noTrailing + "\n\n\n").Should().Equal(baseline, because: "multiple trailing newlines are a no-op");
    }

    #endregion
}
