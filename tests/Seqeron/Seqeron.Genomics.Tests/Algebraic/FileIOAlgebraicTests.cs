namespace Seqeron.Genomics.Tests.Algebraic;

using FastqRecord = FastqParser.FastqRecord;
using BedRecord = BedParser.BedRecord;
using VcfRecord = VcfParser.VcfRecord;
using GffRecord = GffParser.GffRecord;

/// <summary>
/// Algebraic-law tests for the FileIO area (FASTA/FASTQ/BED/VCF/GFF round-trips).
///
/// Algebraic testing pins the parse∘serialize round-trip isomorphism (RT) of each
/// file format: writing a set of records and parsing the result reproduces them,
/// and the empty input is the identity (zero records).
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 64–68.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("FileIO")]
public class FileIOAlgebraicTests
{
    private static string Write(System.Action<TextWriter> write)
    {
        var sw = new StringWriter();
        write(sw);
        return sw.ToString();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-FASTA-001 — FASTA I/O (FileIO), checklist row 64.
    // RT: Parse(ToFasta(entries)) = entries.  ID: empty → 0 records.
    //   — FastaParser.ToFasta / Parse.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Fasta_RoundTrip_PreservesEntries()
    {
        var entries = new[]
        {
            new FastaEntry("seq1", "first", new DnaSequence("ACGTACGTACGT")),
            new FastaEntry("seq2", null, new DnaSequence("TTTGGGCCCAAA")),
            new FastaEntry("seq3", "third", new DnaSequence("AAACCCGGGTTT")),
        };
        var parsed = FastaParser.Parse(FastaParser.ToFasta(entries)).ToList();

        parsed.Should().HaveCount(entries.Length);
        for (int i = 0; i < entries.Length; i++)
        {
            parsed[i].Id.Should().Be(entries[i].Id);
            parsed[i].Sequence.Sequence.Should().Be(entries[i].Sequence.Sequence);
        }
    }

    [Test]
    public void Fasta_Identity_EmptyYieldsNoRecords()
    {
        FastaParser.Parse(FastaParser.ToFasta(System.Array.Empty<FastaEntry>())).Should().BeEmpty();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-FASTQ-001 — FASTQ I/O (FileIO), checklist row 65.
    // RT: Parse(write(records)) = records.  ID: empty → 0 records.
    //   — FastqParser.WriteToStream / Parse.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Fastq_RoundTrip_PreservesRecords()
    {
        var records = new[]
        {
            new FastqRecord("r1", "", "ACGTACGT", "IIIIIIII", new[] { 40, 40, 40, 40, 40, 40, 40, 40 }),
            new FastqRecord("r2", "", "GGGGCCCC", "!!!!####", new[] { 0, 0, 0, 0, 2, 2, 2, 2 }),
        };
        string text = Write(w => FastqParser.WriteToStream(w, records));
        var parsed = FastqParser.Parse(text).ToList();

        parsed.Should().HaveCount(records.Length);
        for (int i = 0; i < records.Length; i++)
        {
            parsed[i].Id.Should().Be(records[i].Id);
            parsed[i].Sequence.Should().Be(records[i].Sequence);
            parsed[i].QualityString.Should().Be(records[i].QualityString);
        }
    }

    [Test]
    public void Fastq_Identity_EmptyYieldsNoRecords()
    {
        string text = Write(w => FastqParser.WriteToStream(w, System.Array.Empty<FastqRecord>()));
        FastqParser.Parse(text).Should().BeEmpty();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-BED-001 — BED I/O (FileIO), checklist row 66.
    // RT: Parse(write(regions)) = regions (BED6 fields).  ID: empty → 0 regions.
    //   — BedParser.WriteToStream / Parse.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Bed_RoundTrip_PreservesRecords()
    {
        var records = new[]
        {
            new BedRecord("chr1", 100, 200, "featA", 500, '+'),
            new BedRecord("chr2", 300, 450, "featB", 750, '-'),
        };
        string text = Write(w => BedParser.WriteToStream(w, records, BedParser.BedFormat.BED6));
        var parsed = BedParser.Parse(text, BedParser.BedFormat.BED6).ToList();

        parsed.Should().HaveCount(records.Length);
        for (int i = 0; i < records.Length; i++)
        {
            parsed[i].Chrom.Should().Be(records[i].Chrom);
            parsed[i].ChromStart.Should().Be(records[i].ChromStart);
            parsed[i].ChromEnd.Should().Be(records[i].ChromEnd);
            parsed[i].Name.Should().Be(records[i].Name);
            parsed[i].Score.Should().Be(records[i].Score);
            parsed[i].Strand.Should().Be(records[i].Strand);
        }
    }

    [Test]
    public void Bed_Identity_EmptyYieldsNoRecords()
    {
        string text = Write(w => BedParser.WriteToStream(w, System.Array.Empty<BedRecord>(), BedParser.BedFormat.BED6));
        BedParser.Parse(text, BedParser.BedFormat.BED6).Should().BeEmpty();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-VCF-001 — VCF I/O (FileIO), checklist row 67.
    // RT: Parse(write(variants)) = variants.  ID: header only → 0 variants.
    //   — VcfParser.WriteToStream / Parse.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Vcf_RoundTrip_PreservesRecords()
    {
        var records = new[]
        {
            new VcfRecord("chr1", 1000, "rs1", "A", new[] { "T" }, 60.0, new[] { "PASS" },
                new Dictionary<string, string> { ["DP"] = "30" }),
            new VcfRecord("chr2", 2000, "rs2", "G", new[] { "C" }, 99.0, new[] { "PASS" },
                new Dictionary<string, string> { ["DP"] = "50" }),
        };
        string text = Write(w => VcfParser.WriteToStream(w, records));
        var parsed = VcfParser.Parse(text).ToList();

        parsed.Should().HaveCount(records.Length);
        for (int i = 0; i < records.Length; i++)
        {
            parsed[i].Chrom.Should().Be(records[i].Chrom);
            parsed[i].Pos.Should().Be(records[i].Pos);
            parsed[i].Id.Should().Be(records[i].Id);
            parsed[i].Ref.Should().Be(records[i].Ref);
            parsed[i].Alt.Should().BeEquivalentTo(records[i].Alt);
        }
    }

    [Test]
    public void Vcf_Identity_HeaderOnlyYieldsNoRecords()
    {
        string text = Write(w => VcfParser.WriteToStream(w, System.Array.Empty<VcfRecord>()));
        VcfParser.Parse(text).Should().BeEmpty();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-GFF-001 — GFF I/O (FileIO), checklist row 68.
    // RT: Parse(write(features)) = features.  ID: empty → 0 features.
    //   — GffParser.WriteToStream / Parse.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Gff_RoundTrip_PreservesRecords()
    {
        var records = new[]
        {
            new GffRecord("chr1", "src", "gene", 100, 200, 0.0, '+', 0,
                new Dictionary<string, string> { ["ID"] = "g1" }),
            new GffRecord("chr2", "src", "mRNA", 300, 400, null, '-', null,
                new Dictionary<string, string> { ["ID"] = "m1" }),
        };
        string text = Write(w => GffParser.WriteToStream(w, records, GffParser.GffFormat.GFF3));
        var parsed = GffParser.Parse(text, GffParser.GffFormat.GFF3).ToList();

        parsed.Should().HaveCount(records.Length);
        for (int i = 0; i < records.Length; i++)
        {
            parsed[i].Seqid.Should().Be(records[i].Seqid);
            parsed[i].Source.Should().Be(records[i].Source);
            parsed[i].Type.Should().Be(records[i].Type);
            parsed[i].Start.Should().Be(records[i].Start);
            parsed[i].End.Should().Be(records[i].End);
            parsed[i].Strand.Should().Be(records[i].Strand);
            parsed[i].Attributes.GetValueOrDefault("ID").Should().Be(records[i].Attributes.GetValueOrDefault("ID"));
        }
    }

    [Test]
    public void Gff_Identity_EmptyYieldsNoRecords()
    {
        string text = Write(w => GffParser.WriteToStream(w, System.Array.Empty<GffRecord>(), GffParser.GffFormat.GFF3));
        GffParser.Parse(text, GffParser.GffFormat.GFF3).Should().BeEmpty();
    }
}
