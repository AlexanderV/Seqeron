using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-BED-001 — BED parsing / serialization (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 66.
    //
    // API under test (BedParser.Parse / WriteToStream / Sort):
    //   BED describes genomic intervals with 0-based, half-open coordinates [chromStart,
    //   chromEnd). BED6 columns are chrom, start, end, name, score (0–1000), strand. The
    //   feature length is chromEnd − chromStart.
    //
    // Relations (derived from the format & coordinate semantics, NOT from output):
    //   • INV  (sort preserves content): Sort only reorders records, so its output is a
    //          permutation of the input — every record (all fields) is preserved and Sort is
    //          idempotent.
    //   • COMP (round-trip): serialising BED6 records and re-parsing recovers every field, and
    //          re-serialising reproduces the exact text.
    //   • SHIFT (coordinate integrity): translating an interval by a constant δ (start+δ,
    //          end+δ) leaves its length end−start invariant, and the exact coordinates survive
    //          a write→parse round-trip.
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-BED-001 — Helpers

    private static (string Chrom, int Start, int End, string? Name, int? Score, char? Strand) Bed6(BedParser.BedRecord r)
        => (r.Chrom, r.ChromStart, r.ChromEnd, r.Name, r.Score, r.Strand);

    private static string WriteBed(IEnumerable<BedParser.BedRecord> records, BedParser.BedFormat format = BedParser.BedFormat.BED6)
    {
        using var sw = new StringWriter();
        BedParser.WriteToStream(sw, records, format);
        return sw.ToString();
    }

    private static BedParser.BedRecord[] SampleBedRecords() => new[]
    {
        new BedParser.BedRecord("chr1", 100, 200, "geneA", 500, '+'),
        new BedParser.BedRecord("chr2", 50, 75, "geneB", 0, '-'),
        new BedParser.BedRecord("chr1", 300, 350, "geneC", 1000, '.'),
    };

    #endregion

    #region PARSE-BED-001 INV — sorting reorders but never alters record content

    [Test]
    [Description("INV: Sort is a pure reordering — its output is a permutation of the input (every field of every record preserved) and it is idempotent.")]
    public void Bed_Sort_PreservesRecordContentAndIsIdempotent()
    {
        var records = SampleBedRecords();

        var sorted = BedParser.Sort(records).ToList();

        sorted.Select(Bed6).Should().BeEquivalentTo(records.Select(Bed6),
            because: "sorting changes only the order; the multiset of records (all fields) is unchanged");

        BedParser.Sort(sorted).Should().Equal(sorted,
            because: "Sort is idempotent — re-sorting an already-sorted list is a no-op");
    }

    #endregion

    #region PARSE-BED-001 COMP — write→parse→write reproduces the records and text

    [Test]
    [Description("COMP: serialising BED6 records and re-parsing recovers every field, and re-serialising the parsed records reproduces the exact text.")]
    public void Bed_WriteParseWrite_IsIdempotentAndRecordPreserving()
    {
        var records = SampleBedRecords();

        string written = WriteBed(records);
        var parsed = BedParser.Parse(written).ToList();

        parsed.Select(Bed6).Should().Equal(records.Select(Bed6),
            because: "parsing a BED6 serialisation must recover chrom/start/end/name/score/strand of every record");

        WriteBed(parsed).Should().Be(written,
            because: "serialisation is a canonical fixed point: write∘parse∘write = write");
    }

    #endregion

    #region PARSE-BED-001 SHIFT — uniform coordinate shift preserves length and survives round-trip

    [Test]
    [Description("SHIFT: translating an interval by a constant δ preserves its length end−start, and the shifted coordinates survive a write→parse round-trip exactly.")]
    public void Bed_CoordinateShift_PreservesLengthAndRoundTrips()
    {
        var records = SampleBedRecords();
        var originalLengths = records.Select(r => r.Length).ToList();

        foreach (int delta in new[] { 0, 5, 1000 })
        {
            var shifted = records
                .Select(r => r with { ChromStart = r.ChromStart + delta, ChromEnd = r.ChromEnd + delta })
                .ToList();

            shifted.Select(r => r.Length).Should().Equal(originalLengths,
                because: $"adding {delta} to both endpoints translates the interval without resizing it — length end−start is shift-invariant");

            var roundTripped = BedParser.Parse(WriteBed(shifted)).ToList();
            roundTripped.Select(r => (r.ChromStart, r.ChromEnd))
                .Should().Equal(shifted.Select(r => (r.ChromStart, r.ChromEnd)),
                    because: "the integer coordinates must survive serialisation and parsing unchanged");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-VCF-001 — VCF parsing / serialization (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 67.
    //
    // API under test (VcfParser.Parse / WriteToStream):
    //   VCF stores variants as tab-delimited data lines (CHROM POS ID REF ALT QUAL FILTER
    //   INFO …) preceded by '##' meta lines and a '#CHROM' column header. The INFO column is
    //   a ';'-separated set of key=value pairs (or bare flags) — an unordered map.
    //
    // Relations (derived from the format, NOT from output):
    //   • INV  (comment/meta lines): all '#'-prefixed and blank lines are skipped, so inserting
    //          extra '##' meta lines or blank lines anywhere leaves the variant records unchanged.
    //   • COMP (round-trip): serialising variant records and re-parsing recovers every column,
    //          and re-serialising reproduces the exact text.
    //   • INV  (INFO order): INFO is an unordered key→value map, so permuting the INFO subfields
    //          of a data line yields semantically identical records (same INFO map, same columns).
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-VCF-001 — Helpers

    /// <summary>Run-order- and INFO-order-independent identity of a VCF record.</summary>
    private static (string Chrom, int Pos, string Id, string Ref, string Alt, double? Qual, string Filter, string Info)
        VcfKey(VcfParser.VcfRecord r) =>
        (
            r.Chrom, r.Pos, r.Id, r.Ref,
            string.Join(',', r.Alt),
            r.Qual,
            string.Join(';', r.Filter),
            string.Join(';', r.Info.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key}={kv.Value}"))
        );

    private static string WriteVcf(IEnumerable<VcfParser.VcfRecord> records)
    {
        using var sw = new StringWriter();
        VcfParser.WriteToStream(sw, records);
        return sw.ToString();
    }

    private static VcfParser.VcfRecord[] SampleVcfRecords() => new[]
    {
        new VcfParser.VcfRecord("chr1", 100, "rs1", "A", new[] { "G" }, 50.0, new[] { "PASS" },
            new Dictionary<string, string> { ["DP"] = "30", ["AF"] = "0.5" }),
        new VcfParser.VcfRecord("chr2", 200, ".", "C", new[] { "T" }, 99.0, new[] { "PASS" },
            new Dictionary<string, string> { ["DP"] = "12" }),
    };

    #endregion

    #region PARSE-VCF-001 INV — meta/comment and blank lines do not affect variant records

    [Test]
    [Description("INV: '#'-prefixed meta lines and blank lines are skipped, so scattering extra '##' lines and blanks through the file leaves the parsed variant records unchanged.")]
    public void Vcf_CommentAndBlankLines_DoNotAffectRecords()
    {
        const string compact =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30;AF=0.5\n" +
            "chr2\t200\t.\tC\tT\t99\tPASS\tDP=12\n";

        const string withComments =
            "##fileformat=VCFv4.3\n" +
            "##source=unit-test\n" +
            "##contig=<ID=chr1,length=1000>\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "\n" +
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30;AF=0.5\n" +
            "##an=annotation-line-between-records\n" +
            "\n" +
            "chr2\t200\t.\tC\tT\t99\tPASS\tDP=12\n";

        VcfParser.Parse(withComments).Select(VcfKey)
            .Should().Equal(VcfParser.Parse(compact).Select(VcfKey),
                because: "meta and blank lines carry no variant data, so they cannot change the parsed records");
    }

    #endregion

    #region PARSE-VCF-001 COMP — write→parse→write reproduces the records and text

    [Test]
    [Description("COMP: serialising variant records and re-parsing recovers every column, and re-serialising the parsed records reproduces the exact text.")]
    public void Vcf_WriteParseWrite_IsIdempotentAndRecordPreserving()
    {
        var records = SampleVcfRecords();

        string written = WriteVcf(records);
        var parsed = VcfParser.Parse(written).ToList();

        parsed.Select(VcfKey).Should().Equal(records.Select(VcfKey),
            because: "parsing a VCF serialisation must recover chrom/pos/id/ref/alt/qual/filter/info of every variant");

        WriteVcf(parsed).Should().Be(written,
            because: "serialisation is a canonical fixed point: write∘parse∘write = write");
    }

    #endregion

    #region PARSE-VCF-001 INV — INFO subfield order is irrelevant

    [Test]
    [Description("INV: the INFO column is an unordered key→value map, so permuting its subfields yields a record with the same INFO map and identical other columns.")]
    public void Vcf_InfoFieldOrder_IsIrrelevant()
    {
        const string header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        // Same variant, INFO subfields (including a bare SOMATIC flag) in two different orders.
        string orderA = header + "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30;AF=0.5;SOMATIC\n";
        string orderB = header + "chr1\t100\trs1\tA\tG\t50\tPASS\tSOMATIC;AF=0.5;DP=30\n";

        var recordA = VcfParser.Parse(orderA).Single();
        var recordB = VcfParser.Parse(orderB).Single();

        VcfKey(recordB).Should().Be(VcfKey(recordA),
            because: "INFO is a set of key=value pairs, so its serialised order does not change the parsed variant");
        recordB.Info.Should().BeEquivalentTo(recordA.Info,
            because: "the INFO map (DP, AF and the SOMATIC flag) is identical regardless of subfield order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-GFF-001 — GFF3 parsing / serialization (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 68.
    //
    // API under test (GffParser.Parse / WriteToStream):
    //   GFF3 stores features as 9 tab-delimited columns (seqid source type start end score
    //   strand phase attributes) preceded by '##' directives and '#' comments. Column 9 is a
    //   ';'-separated set of key=value attribute pairs — an unordered map.
    //
    // Relations (derived from the format, NOT from output):
    //   • COMP (round-trip): serialising feature records and re-parsing recovers every column,
    //          and re-serialising reproduces the exact text.
    //   • INV  (attribute order): column-9 attributes are an unordered map, so permuting them
    //          yields semantically identical features (same attribute map, same columns).
    //   • INV  (comment lines): '##' directives, '#' comments and blank lines carry no feature
    //          data, so inserting them anywhere leaves the parsed features unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-GFF-001 — Helpers

    /// <summary>Run-order- and attribute-order-independent identity of a GFF feature.</summary>
    private static (string Seqid, string Source, string Type, int Start, int End, double? Score, char Strand, int? Phase, string Attrs)
        GffKey(GffParser.GffRecord r) =>
        (
            r.Seqid, r.Source, r.Type, r.Start, r.End, r.Score, r.Strand, r.Phase,
            string.Join(';', r.Attributes.OrderBy(kv => kv.Key, System.StringComparer.Ordinal).Select(kv => $"{kv.Key}={kv.Value}"))
        );

    private static string WriteGff(IEnumerable<GffParser.GffRecord> records)
    {
        using var sw = new StringWriter();
        GffParser.WriteToStream(sw, records, GffParser.GffFormat.GFF3);
        return sw.ToString();
    }

    private static GffParser.GffRecord[] SampleGffRecords() => new[]
    {
        new GffParser.GffRecord("chr1", "ensembl", "gene", 1000, 5000, null, '+', null,
            new Dictionary<string, string>(System.StringComparer.Ordinal) { ["ID"] = "gene1", ["Name"] = "BRCA1" }),
        new GffParser.GffRecord("chr1", "ensembl", "exon", 1000, 1200, 0.0, '+', 0,
            new Dictionary<string, string>(System.StringComparer.Ordinal) { ["ID"] = "exon1", ["Parent"] = "gene1" }),
    };

    #endregion

    #region PARSE-GFF-001 COMP — write→parse→write reproduces the features and text

    [Test]
    [Description("COMP: serialising GFF3 features and re-parsing recovers every column, and re-serialising the parsed features reproduces the exact text.")]
    public void Gff_WriteParseWrite_IsIdempotentAndRecordPreserving()
    {
        var records = SampleGffRecords();

        string written = WriteGff(records);
        var parsed = GffParser.Parse(written).ToList();

        parsed.Select(GffKey).Should().Equal(records.Select(GffKey),
            because: "parsing a GFF3 serialisation must recover all nine columns of every feature");

        WriteGff(parsed).Should().Be(written,
            because: "serialisation is a canonical fixed point: write∘parse∘write = write");
    }

    #endregion

    #region PARSE-GFF-001 INV — attribute order is irrelevant

    [Test]
    [Description("INV: column-9 attributes are an unordered map, so permuting them yields a feature with the same attribute map and identical other columns.")]
    public void Gff_AttributeOrder_IsIrrelevant()
    {
        const string header = "##gff-version 3\n";

        string orderA = header + "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tID=g1;Name=BRCA;biotype=protein_coding\n";
        string orderB = header + "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tbiotype=protein_coding;Name=BRCA;ID=g1\n";

        var recordA = GffParser.Parse(orderA).Single();
        var recordB = GffParser.Parse(orderB).Single();

        GffKey(recordB).Should().Be(GffKey(recordA),
            because: "the attribute column is a set of key=value pairs, so its order does not change the parsed feature");
        recordB.Attributes.Should().BeEquivalentTo(recordA.Attributes,
            because: "the attribute map (ID, Name, biotype) is identical regardless of order");
    }

    #endregion

    #region PARSE-GFF-001 INV — directives, comments and blank lines do not affect features

    [Test]
    [Description("INV: '##' directives, '#' comments and blank lines are skipped, so scattering them through the file leaves the parsed features unchanged.")]
    public void Gff_CommentAndBlankLines_DoNotAffectFeatures()
    {
        const string compact =
            "##gff-version 3\n" +
            "chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=BRCA1\n" +
            "chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=gene1\n";

        const string withComments =
            "##gff-version 3\n" +
            "##sequence-region chr1 1 100000\n" +
            "# a free-text comment\n" +
            "\n" +
            "chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=BRCA1\n" +
            "# comment between features\n" +
            "\n" +
            "chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=gene1\n";

        GffParser.Parse(withComments).Select(GffKey)
            .Should().Equal(GffParser.Parse(compact).Select(GffKey),
                because: "directives, comments and blank lines carry no feature data, so they cannot change the parsed features");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-GENBANK-001 — GenBank flat-file parsing (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 69.
    //
    // API under test (GenBankParser.Parse):
    //   GenBank is a fixed-column flat file: a LOCUS line, keyword sections (DEFINITION,
    //   ACCESSION, …) and an ORIGIN block that lists the sequence with a left-margin base
    //   position followed by groups of bases. The parser keeps only the letters of the
    //   ORIGIN block (uppercased), discarding positions and layout whitespace.
    //   (There is no GenBank writer, so the canonical encoder below stands in for one.)
    //
    // Relations (derived from the fixed-column format, NOT from output):
    //   • COMP (round-trip identity): encoding known field values into a canonical GenBank
    //          record and parsing recovers those values exactly (LOCUS/DEFINITION/ACCESSION
    //          and the uppercased sequence), and the sequence is a fixed point of
    //          encode-ORIGIN∘parse.
    //   • INV  (sequence whitespace/layout): the parsed sequence is invariant to the ORIGIN
    //          block's base-position numbers, grouping, line width and letter case — only the
    //          letters, uppercased, survive.
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-GENBANK-001 — Helpers

    /// <summary>Builds a canonical GenBank ORIGIN block (left-margin position + space-separated base groups).</summary>
    private static string OriginBlock(string sequence, int perLine = 60, int groupSize = 10, bool upperCase = false)
    {
        string bases = upperCase ? sequence.ToUpperInvariant() : sequence.ToLowerInvariant();
        var sb = new StringBuilder("ORIGIN\n");
        for (int i = 0; i < bases.Length; i += perLine)
        {
            sb.Append((i + 1).ToString().PadLeft(9));
            for (int j = i; j < Math.Min(i + perLine, bases.Length); j += groupSize)
                sb.Append(' ').Append(bases.Substring(j, Math.Min(groupSize, bases.Length - j)));
            sb.Append('\n');
        }
        return sb.ToString();
    }

    /// <summary>Builds a minimal but valid GenBank record around the given ORIGIN block.</summary>
    private static string GenBankRecordText(string locus, string definition, string accession, string originBlock, int length) =>
        $"LOCUS       {locus}     {length} bp    DNA     linear   UNK 01-JAN-2020\n" +
        $"DEFINITION  {definition}\n" +
        $"ACCESSION   {accession}\n" +
        $"VERSION     {accession}.1\n" +
        originBlock +
        "//\n";

    #endregion

    #region PARSE-GENBANK-001 COMP — encode→parse recovers the fields and sequence

    [Test]
    [Description("COMP: encoding known field values into a canonical GenBank record and parsing recovers them exactly, and the sequence is a fixed point of encode-ORIGIN∘parse.")]
    public void GenBank_EncodeParse_RecoversFieldsAndSequence()
    {
        const string locus = "TEST123";
        const string definition = "Synthetic test sequence";
        const string accession = "TEST123";
        const string sequence = "ACGTACGTACGTACGTTTGGCCAA";

        string text = GenBankRecordText(locus, definition, accession, OriginBlock(sequence), sequence.Length);
        var record = GenBankParser.Parse(text).Single();

        record.Locus.Should().Be(locus, because: "the LOCUS name must round-trip");
        record.Definition.Should().Be(definition, because: "the DEFINITION text must round-trip");
        record.Accession.Should().Be(accession, because: "the ACCESSION must round-trip");
        record.Sequence.Should().Be(sequence, because: "the ORIGIN bases (uppercased) must round-trip exactly");

        // Fixed point: re-encoding the parsed sequence and re-parsing yields the same sequence.
        var reparsed = GenBankParser.Parse(
            GenBankRecordText(locus, definition, accession, OriginBlock(record.Sequence), record.Sequence.Length)).Single();
        reparsed.Sequence.Should().Be(record.Sequence, because: "the sequence is a fixed point of encode-ORIGIN∘parse");
    }

    #endregion

    #region PARSE-GENBANK-001 INV — ORIGIN layout/whitespace does not change the sequence

    [Test]
    [Description("INV: the parsed sequence keeps only the ORIGIN letters (uppercased), so it is invariant to base-position numbers, grouping, line width and letter case.")]
    public void GenBank_OriginLayout_DoesNotChangeSequence()
    {
        const string sequence = "ACGTACGTACGTACGTTTGGCCAATTGGCCAA";
        const string locus = "LAY";
        const string accession = "LAY001";
        const string definition = "Layout invariance";

        // The same sequence rendered with different line widths, group sizes and case.
        var layouts = new[]
        {
            OriginBlock(sequence, perLine: 60, groupSize: 10, upperCase: false),
            OriginBlock(sequence, perLine: 30, groupSize: 10, upperCase: false),
            OriginBlock(sequence, perLine: 20, groupSize: 5,  upperCase: false),
            OriginBlock(sequence, perLine: 60, groupSize: 10, upperCase: true),
        };

        foreach (var origin in layouts)
        {
            var record = GenBankParser.Parse(
                GenBankRecordText(locus, definition, accession, origin, sequence.Length)).Single();

            record.Sequence.Should().Be(sequence,
                because: "position numbers, grouping, line width and case are not part of the sequence — only the letters, uppercased, are kept");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PARSE-EMBL-001 — EMBL flat-file parsing (FileIO).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 70.
    //
    // API under test (EmblParser.Parse):
    //   EMBL uses two-letter line prefixes (ID, AC, DE, SQ, …). The SQ block lists the
    //   sequence as space-separated base groups with a right-margin position; the parser
    //   keeps only the letters of the SQ block (uppercased), discarding positions and layout.
    //   (There is no EMBL writer, so the canonical encoder below stands in for one.)
    //
    // Relations (derived from the prefixed-line format, NOT from output):
    //   • COMP (round-trip identity): encoding known field values into a canonical EMBL record
    //          and parsing recovers them exactly (AC/DE and the uppercased sequence), and the
    //          sequence is a fixed point of encode-SQ∘parse.
    //   • INV  (sequence whitespace/layout): the parsed sequence is invariant to the SQ block's
    //          base-position numbers, grouping, line width and letter case.
    // ───────────────────────────────────────────────────────────────────────────

    #region PARSE-EMBL-001 — Helpers

    /// <summary>Builds a canonical EMBL SQ block (5-space margin, space-separated base groups, right-margin position).</summary>
    private static string EmblSqBlock(string sequence, int perLine = 60, int groupSize = 10, bool upperCase = false)
    {
        string bases = upperCase ? sequence.ToUpperInvariant() : sequence.ToLowerInvariant();
        var sb = new StringBuilder($"SQ   Sequence {sequence.Length} BP;\n");
        for (int i = 0; i < bases.Length; i += perLine)
        {
            sb.Append("    ");
            int end = Math.Min(i + perLine, bases.Length);
            for (int j = i; j < end; j += groupSize)
                sb.Append(' ').Append(bases.Substring(j, Math.Min(groupSize, bases.Length - j)));
            sb.Append("      ").Append(end).Append('\n'); // right-margin base position (digits, dropped on parse)
        }
        return sb.ToString();
    }

    /// <summary>Builds a minimal but valid EMBL record around the given SQ block.</summary>
    private static string EmblRecordText(string accession, string description, string sqBlock, int length) =>
        $"ID   {accession}; SV 1; linear; genomic DNA; STD; HUM; {length} BP.\n" +
        "XX\n" +
        $"AC   {accession};\n" +
        "XX\n" +
        $"DE   {description}\n" +
        "XX\n" +
        sqBlock +
        "//\n";

    #endregion

    #region PARSE-EMBL-001 COMP — encode→parse recovers the fields and sequence

    [Test]
    [Description("COMP: encoding known field values into a canonical EMBL record and parsing recovers them exactly, and the sequence is a fixed point of encode-SQ∘parse.")]
    public void Embl_EncodeParse_RecoversFieldsAndSequence()
    {
        const string accession = "TEST123";
        const string description = "Synthetic test sequence";
        const string sequence = "ACGTACGTACGTACGTTTGGCCAA";

        string text = EmblRecordText(accession, description, EmblSqBlock(sequence), sequence.Length);
        var record = EmblParser.Parse(text).Single();

        record.Accession.Should().Be(accession, because: "the AC/ID accession must round-trip");
        record.Description.Should().Be(description, because: "the DE description must round-trip");
        record.Sequence.Should().Be(sequence, because: "the SQ bases (uppercased) must round-trip exactly");

        var reparsed = EmblParser.Parse(
            EmblRecordText(accession, description, EmblSqBlock(record.Sequence), record.Sequence.Length)).Single();
        reparsed.Sequence.Should().Be(record.Sequence, because: "the sequence is a fixed point of encode-SQ∘parse");
    }

    #endregion

    #region PARSE-EMBL-001 INV — SQ layout/whitespace does not change the sequence

    [Test]
    [Description("INV: the parsed sequence keeps only the SQ letters (uppercased), so it is invariant to base-position numbers, grouping, line width and letter case.")]
    public void Embl_SqLayout_DoesNotChangeSequence()
    {
        const string sequence = "ACGTACGTACGTACGTTTGGCCAATTGGCCAA";
        const string accession = "LAY001";
        const string description = "Layout invariance";

        var layouts = new[]
        {
            EmblSqBlock(sequence, perLine: 60, groupSize: 10, upperCase: false),
            EmblSqBlock(sequence, perLine: 30, groupSize: 10, upperCase: false),
            EmblSqBlock(sequence, perLine: 20, groupSize: 5,  upperCase: false),
            EmblSqBlock(sequence, perLine: 60, groupSize: 10, upperCase: true),
        };

        foreach (var sq in layouts)
        {
            var record = EmblParser.Parse(
                EmblRecordText(accession, description, sq, sequence.Length)).Single();

            record.Sequence.Should().Be(sequence,
                because: "position numbers, grouping, line width and case are not part of the sequence — only the letters, uppercased, are kept");
        }
    }

    #endregion
}
