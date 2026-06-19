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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-FASTQ-001 — FASTQ parsing / quality encoding (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 65.
    //
    // API under test (FastqParser.Parse / ToFastqString / Encode/DecodeQualityScores /
    //                 InterleavePairedReads / SplitInterleavedReads):
    //   A FASTQ record is (@header, sequence, +separator, quality). Quality characters encode
    //   Phred scores via a fixed ASCII offset (33 for Sanger/Illumina-1.8+, 64 for legacy
    //   Illumina), so score = ord(char) − offset and char = chr(score + offset).
    //
    // Relations (derived from the format & the offset encoding, NOT from output):
    //   • COMP (round-trip): parsing a serialised record set recovers each record's
    //          ID/description/sequence/quality, and re-serialising reproduces the exact text.
    //   • INV  (offset consistency): decode∘encode is the identity for in-range scores, and
    //          the SAME score encoded under Phred+64 vs Phred+33 differs by exactly 64−33 = 31
    //          in the character code — the offset is the only thing that changes.
    //   • INV  (interleaved order preserved): splitting an interleaved pair stream restores the
    //          two original read lists in their original order (split∘interleave = identity).
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-FASTQ-001 — Helpers

    private static FastqParser.FastqRecord MakeFastqRecord(string id, string description, string sequence, string quality)
    {
        var scores = FastqParser.DecodeQualityScores(quality, FastqParser.QualityEncoding.Phred33);
        return new FastqParser.FastqRecord(id, description, sequence, quality, scores);
    }

    private static string SerializeFastq(IEnumerable<FastqParser.FastqRecord> records) =>
        string.Concat(records.Select(FastqParser.ToFastqString));

    #endregion

    #region PARSE-FASTQ-001 COMP — write→parse→write reproduces the records and text

    [Test]
    [Description("COMP: parsing serialised FASTQ recovers each record's ID/description/sequence/quality, and re-serialising the parsed records reproduces the exact text.")]
    public void Fastq_WriteParseWrite_IsIdempotentAndRecordPreserving()
    {
        // Quality strings contain '#' (ASCII 35 < 64), so Auto-detection resolves to Phred+33.
        var records = new[]
        {
            MakeFastqRecord("read1", "pair 1", "ACGTACGT", "IIIII###"),
            MakeFastqRecord("read2", "",       "TTGGCCAA", "5I5I5I5I"),
        };

        string written = SerializeFastq(records);
        var parsed = FastqParser.Parse(written).ToList();

        parsed.Select(r => (r.Id, r.Description, r.Sequence, r.QualityString))
            .Should().Equal(records.Select(r => (r.Id, r.Description, r.Sequence, r.QualityString)),
                because: "parse must recover the header, sequence and quality of every record");

        SerializeFastq(parsed).Should().Be(written,
            because: "serialisation is a canonical fixed point: write∘parse∘write = write");
    }

    #endregion

    #region PARSE-FASTQ-001 INV — quality offset is consistent

    [Test]
    [Description("INV: decode∘encode is the identity for in-range scores, and encoding the same score under Phred+64 vs Phred+33 shifts the character code by exactly 31 (= 64−33).")]
    public void Fastq_QualityEncoding_OffsetIsConsistent()
    {
        int[] scores = { 0, 2, 20, 30, 40, 60 }; // all within both Phred+33 (≤93) and Phred+64 (≤62) ranges

        foreach (var encoding in new[] { FastqParser.QualityEncoding.Phred33, FastqParser.QualityEncoding.Phred64 })
        {
            string encoded = FastqParser.EncodeQualityScores(scores, encoding);
            FastqParser.DecodeQualityScores(encoded, encoding).Should().Equal(scores,
                because: $"decoding the {encoding} encoding of a score recovers the score exactly");
        }

        string p33 = FastqParser.EncodeQualityScores(scores, FastqParser.QualityEncoding.Phred33);
        string p64 = FastqParser.EncodeQualityScores(scores, FastqParser.QualityEncoding.Phred64);

        for (int i = 0; i < scores.Length; i++)
            (p64[i] - p33[i]).Should().Be(31,
                because: "Phred+64 and Phred+33 differ only by the offset 64−33 = 31 for the same score");
    }

    #endregion

    #region PARSE-FASTQ-001 INV — interleaving order is preserved

    [Test]
    [Description("INV: splitting an interleaved pair stream restores the two original read lists in order — split∘interleave is the identity for equal-length mate lists.")]
    public void Fastq_InterleaveThenSplit_RestoresOriginalOrder()
    {
        var read1 = new[]
        {
            MakeFastqRecord("r1a", "", "AAAA", "IIII"),
            MakeFastqRecord("r1b", "", "CCCC", "5555"),
        };
        var read2 = new[]
        {
            MakeFastqRecord("r2a", "", "GGGG", "####"),
            MakeFastqRecord("r2b", "", "TTTT", "IIII"),
        };

        var interleaved = FastqParser.InterleavePairedReads(read1, read2).ToList();

        interleaved.Select(r => r.Id).Should().Equal(new[] { "r1a", "r2a", "r1b", "r2b" },
            because: "interleaving emits mates in strict alternating order: read1[0], read2[0], read1[1], read2[1]");

        var (split1, split2) = FastqParser.SplitInterleavedReads(interleaved);

        split1.Should().Equal(read1, because: "the even positions reconstruct read1 in original order");
        split2.Should().Equal(read2, because: "the odd positions reconstruct read2 in original order");
    }

    #endregion
}
