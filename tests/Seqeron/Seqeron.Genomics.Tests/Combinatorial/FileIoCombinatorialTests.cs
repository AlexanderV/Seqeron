namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the FileIO area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("FileIO")]
public class FileIoCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static string Wrap(string seq, int width) =>
        width <= 0
            ? seq
            : string.Join("\n", Enumerable.Range(0, (seq.Length + width - 1) / width)
                .Select(i => seq.Substring(i * width, Math.Min(width, seq.Length - i * width))));

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-FASTA-001 — FASTA parsing (FileIO)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 64.
    // Spec: tests/TestSpecs/PARSE-FASTA-001.md (canonical FastaParser.Parse).
    // Dimensions: multiRecord(2) × lineWrap(3) × headerStyle(3). Grid 2×3×3 = 18.
    //
    // Model (Pearson FASTA): each record is a '>' header line (id = first token, description =
    // the rest) followed by sequence lines; the residues are the concatenation of all sequence
    // lines with whitespace removed — so the parsed sequence is INVARIANT to line wrapping.
    //
    // The combinatorial point: record count, line-wrap width and header style interact, yet the
    // parser recovers every record's id, description and (wrap-independent) sequence in order.
    // ═══════════════════════════════════════════════════════════════════════

    public enum HeaderStyle { IdOnly, IdDescription, PipedAccession }

    private static (string Id, string? Desc, string Header) Header(int r, HeaderStyle style) => style switch
    {
        HeaderStyle.IdOnly => ($"seq{r}", null, $"seq{r}"),
        HeaderStyle.IdDescription => ($"seq{r}", $"description number {r}", $"seq{r} description number {r}"),
        // NCBI piped style has no space ⇒ the whole token is the id (parser splits on whitespace).
        _ => ($"sp|ACC{r}|seq{r}", null, $"sp|ACC{r}|seq{r}"),
    };

    [Test, Combinatorial]
    public void ParseFasta_RecoversRecords_AcrossWrapAndHeaderStyle(
        [Values(1, 3)] int nRecords,
        [Values(10, 25, 0)] int lineWrap, // 0 = single line
        [Values(HeaderStyle.IdOnly, HeaderStyle.IdDescription, HeaderStyle.PipedAccession)] HeaderStyle style)
    {
        var records = Enumerable.Range(0, nRecords)
            .Select(r => (Hdr: Header(r, style), Seq: DiverseDna(40 + r * 7, (uint)(0x100 + r))))
            .ToList();

        string content = string.Concat(records.Select(rec => $">{rec.Hdr.Header}\n{Wrap(rec.Seq, lineWrap)}\n"));

        var parsed = FastaParser.Parse(content).ToList();

        parsed.Should().HaveCount(nRecords, "every record is parsed, in order");
        for (int r = 0; r < nRecords; r++)
        {
            parsed[r].Id.Should().Be(records[r].Hdr.Id, "id is the first header token");
            parsed[r].Description.Should().Be(records[r].Hdr.Desc, "description is the header remainder");
            parsed[r].Sequence.Sequence.Should().Be(records[r].Seq, "line wrapping does not change the residues");
        }
    }

    /// <summary>
    /// Interaction witness: the same record wrapped at different widths (and unwrapped) parses to
    /// the identical sequence — parsing is line-wrap invariant.
    /// </summary>
    [Test]
    public void ParseFasta_WrapInvariance()
    {
        string seq = DiverseDna(120, 7u);
        string Parse(int w) => FastaParser.Parse($">s\n{Wrap(seq, w)}\n").Single().Sequence.Sequence;

        Parse(10).Should().Be(seq);
        Parse(60).Should().Be(seq);
        Parse(0).Should().Be(seq);
    }

    /// <summary>
    /// Interaction witness: a write→parse round-trip preserves ids and sequences for a
    /// multi-record file at any output line width.
    /// </summary>
    [Test]
    public void ParseFasta_RoundTripsThroughToFasta()
    {
        var entries = Enumerable.Range(0, 3)
            .Select(r => new FastaEntry($"id{r}", $"desc {r}", new DnaSequence(DiverseDna(50 + r * 11, (uint)(900 + r))))).ToList();

        string text = FastaParser.ToFasta(entries, lineWidth: 30);
        var parsed = FastaParser.Parse(text).ToList();

        parsed.Select(p => p.Id).Should().Equal(entries.Select(e => e.Id));
        parsed.Select(p => p.Sequence.Sequence).Should().Equal(entries.Select(e => e.Sequence.Sequence));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-FASTQ-001 — FASTQ parsing (FileIO)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 65.
    // Spec: tests/TestSpecs/PARSE-FASTQ-001.md (canonical FastqParser.Parse).
    // Dimensions: qualEncoding(3: Phred33/64/auto) × multiRecord(2) × seqLen(3). Grid 3×2×3 = 18.
    //
    // Model (Cock 2010, FASTQ): a record is four lines (@id, sequence, +, quality); each quality
    // character encodes a Phred score as char−offset (33 for Sanger/Illumina 1.8+, 64 for
    // Illumina 1.3–1.7). The checklist's Solexa encoding is not separately implemented; Auto maps
    // to detection (which resolves Phred33 for low-score data). QualityString length == sequence length.
    //
    // The combinatorial point: encoding offset, record count and read length interact — the
    // parser recovers each record's sequence, raw quality string and decoded scores (= char−offset)
    // for every cell, in order.
    // ═══════════════════════════════════════════════════════════════════════

    public enum FastqEncoding { Phred33, Phred64, Auto }

    [Test, Combinatorial]
    public void ParseFastq_DecodesScores_AcrossEncodingAndLength(
        [Values(FastqEncoding.Phred33, FastqEncoding.Phred64, FastqEncoding.Auto)] FastqEncoding encoding,
        [Values(1, 3)] int nRecords,
        [Values(8, 20, 50)] int seqLen)
    {
        int offset = encoding == FastqEncoding.Phred64 ? 64 : 33; // Auto-detect resolves to 33 for these low scores
        var parseEncoding = encoding switch
        {
            FastqEncoding.Phred33 => FastqParser.QualityEncoding.Phred33,
            FastqEncoding.Phred64 => FastqParser.QualityEncoding.Phred64,
            _ => FastqParser.QualityEncoding.Auto,
        };

        var records = new List<(string Id, string Seq, int[] Scores, string Qual)>();
        for (int r = 0; r < nRecords; r++)
        {
            string seq = DiverseDna(seqLen, (uint)(0x200 + r));
            int[] scores = Enumerable.Range(0, seqLen).Select(i => (i + r) % 40).ToArray(); // 0..39 safe for both offsets
            string qual = new string(scores.Select(s => (char)(s + offset)).ToArray());
            records.Add(($"read{r}", seq, scores, qual));
        }

        string content = string.Concat(records.Select(rec => $"@{rec.Id}\n{rec.Seq}\n+\n{rec.Qual}\n"));
        var parsed = FastqParser.Parse(content, parseEncoding).ToList();

        parsed.Should().HaveCount(nRecords);
        for (int r = 0; r < nRecords; r++)
        {
            parsed[r].Id.Should().Be(records[r].Id);
            parsed[r].Sequence.Should().Be(records[r].Seq);
            parsed[r].QualityString.Should().Be(records[r].Qual);
            parsed[r].QualityScores.Should().Equal(records[r].Scores, "scores decode as char − offset");
            parsed[r].QualityScores.Count.Should().Be(seqLen, "one score per base");
        }
    }

    /// <summary>
    /// Interaction witness: the same quality string decodes to different Phred scores under
    /// Phred33 vs Phred64 (the offsets differ by 31).
    /// </summary>
    [Test]
    public void ParseFastq_EncodingOffsetShiftsScores()
    {
        // 'h' = 104 ⇒ Phred33 score 71, Phred64 score 40.
        const string fastq = "@r\nACGT\n+\nhhhh\n";
        FastqParser.Parse(fastq, FastqParser.QualityEncoding.Phred33).Single().QualityScores
            .Should().AllBeEquivalentTo(104 - 33);
        FastqParser.Parse(fastq, FastqParser.QualityEncoding.Phred64).Single().QualityScores
            .Should().AllBeEquivalentTo(104 - 64);
    }

    /// <summary>
    /// Interaction witness: a quality-score → string → score round-trip is identity for both
    /// offsets.
    /// </summary>
    [Test]
    public void ParseFastq_EncodeDecode_RoundTrips()
    {
        int[] scores = { 0, 10, 20, 30, 40 };
        foreach (var enc in new[] { FastqParser.QualityEncoding.Phred33, FastqParser.QualityEncoding.Phred64 })
        {
            string q = FastqParser.EncodeQualityScores(scores, enc);
            FastqParser.DecodeQualityScores(q, enc).Should().Equal(scores);
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-VCF-001 — VCF parsing (FileIO)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 67.
    // Spec: tests/TestSpecs/PARSE-VCF-001.md (canonical VcfParser.ParseWithHeader).
    // Dimensions: nSamples(3) × nVariants(3) × vcfVersion(2) × genotypeFormat(3). Grid 3×3×2×3 = 54.
    //
    // Model (VCF 4.x spec): a VCF has a ##fileformat line, meta/##FORMAT lines, a #CHROM header
    // naming the samples, then one tab-delimited record per variant. Column 9 (FORMAT) names the
    // per-sample sub-fields (colon-separated), and each sample column supplies their values.
    //
    // The combinatorial point: sample count, variant count, declared version and FORMAT layout
    // interact — the parser recovers the fileformat, sample names, every variant and, for each
    // sample, the FORMAT-keyed values (incl. GT) exactly.
    // ═══════════════════════════════════════════════════════════════════════

    private static string SampleGt(int v, int s) => new[] { "0/0", "0/1", "1/1" }[(v + s) % 3];

    private static string FormatValue(string key, int v, int s) => key switch
    {
        "GT" => SampleGt(v, s),
        "DP" => "30",
        "AD" => "15,15",
        _ => "99", // GQ
    };

    [Test, Combinatorial]
    public void ParseVcf_RecoversSamplesAndGenotypes(
        [Values(1, 2, 3)] int nSamples,
        [Values(1, 2, 3)] int nVariants,
        [Values("VCFv4.1", "VCFv4.3")] string version,
        [Values("GT", "GT:DP", "GT:AD:DP:GQ")] string formatField)
    {
        string[] keys = formatField.Split(':');
        var sb = new System.Text.StringBuilder();
        sb.Append($"##fileformat={version}\n");
        foreach (var k in keys)
            sb.Append($"##FORMAT=<ID={k},Number=1,Type=String,Description=\"{k}\">\n");
        sb.Append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int s = 0; s < nSamples; s++) sb.Append($"\tS{s}");
        sb.Append('\n');

        for (int v = 0; v < nVariants; v++)
        {
            sb.Append($"chr1\t{v + 1}\trs{v}\tA\tG\t50\tPASS\t.\t{formatField}");
            for (int s = 0; s < nSamples; s++)
                sb.Append('\t').Append(string.Join(":", keys.Select(k => FormatValue(k, v, s))));
            sb.Append('\n');
        }

        var (header, recordsEnum) = VcfParser.ParseWithHeader(sb.ToString());
        var records = recordsEnum.ToList();

        header.FileFormat.Should().Be(version);
        header.SampleNames.Should().HaveCount(nSamples);
        records.Should().HaveCount(nVariants);

        for (int v = 0; v < nVariants; v++)
        {
            var rec = records[v];
            rec.Format.Should().Equal(keys, "FORMAT column lists the sub-fields");
            rec.Samples.Should().NotBeNull();
            rec.Samples!.Count.Should().Be(nSamples);
            for (int s = 0; s < nSamples; s++)
            {
                foreach (var k in keys)
                    rec.Samples[s][k].Should().Be(FormatValue(k, v, s), $"sample {s} field {k}");
                VcfParser.GetGenotype(rec, s).Should().Be(SampleGt(v, s));
            }
        }
    }

    /// <summary>
    /// Interaction witness: extending the FORMAT layout adds the new sub-fields to every sample
    /// without disturbing the genotype.
    /// </summary>
    [Test]
    public void ParseVcf_FormatLayout_DrivesSampleFields()
    {
        const string header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n";
        var gtOnly = VcfParser.ParseWithHeader(header + "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n").Records.Single();
        gtOnly.Samples![0].Should().ContainKey("GT").And.NotContainKey("DP");

        var withDp = VcfParser.ParseWithHeader(header + "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP\t1/1:42\n").Records.Single();
        withDp.Samples![0]["GT"].Should().Be("1/1");
        withDp.Samples[0]["DP"].Should().Be("42");
    }

    /// <summary>
    /// Interaction witness: without a #CHROM header (no sample names) the standard fields still
    /// parse but no per-sample data is attached.
    /// </summary>
    [Test]
    public void ParseVcf_NoHeader_ParsesSitesWithoutSamples()
    {
        var rec = VcfParser.Parse("chr1\t100\trs1\tA\tG\t60\tPASS\tDP=10\tGT\t0/1\n").Single();
        rec.Chrom.Should().Be("chr1");
        rec.Pos.Should().Be(100);
        rec.Alt.Should().Equal("G");
        rec.Samples.Should().BeNull("sample columns need a #CHROM header to bind names");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-GENBANK-001 — GenBank flat-file parsing (FileIO)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 69.
    // Spec: tests/TestSpecs/PARSE-GENBANK-001.md (canonical GenBankParser.Parse).
    // Dimensions: nFeatures(3) × hasSequence(2) × division(3). Grid 3×2×3 = 18.
    //
    // Model (GenBank flat file): a record opens with a LOCUS line (name, length, molecule type,
    // topology, division code, date), carries a FEATURES table (each feature a key + location +
    // qualifiers) and an optional ORIGIN sequence block, terminated by "//".
    //
    // The combinatorial point: feature count, sequence presence and division code interact — the
    // parser reads the division from LOCUS, the exact feature count, and attaches the residues
    // only when an ORIGIN block is present.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] GbFeatureKeys = { "gene", "CDS", "misc_feature" };

    private static string BuildGenBank(string name, string division, int nFeatures, bool hasSequence)
    {
        var sb = new System.Text.StringBuilder();
        sb.Append($"LOCUS       {name}      20 bp    DNA     linear   {division} 01-JAN-2024\n");
        sb.Append("DEFINITION  Combinatorial test record.\n");
        sb.Append($"ACCESSION   {name}\n");
        sb.Append($"VERSION     {name}.1\n");
        sb.Append("FEATURES             Location/Qualifiers\n");
        for (int i = 0; i < nFeatures; i++)
        {
            sb.Append("     ").Append(GbFeatureKeys[i % GbFeatureKeys.Length].PadRight(16))
              .Append($"{i * 3 + 1}..{i * 3 + 3}\n");
            sb.Append(new string(' ', 21)).Append($"/note=\"f{i}\"\n");
        }
        if (hasSequence)
        {
            sb.Append("ORIGIN      \n");
            sb.Append("        1 acgtacgtac gtacgtacgt\n");
        }
        sb.Append("//\n");
        return sb.ToString();
    }

    [Test, Combinatorial]
    public void ParseGenBank_ReadsDivisionFeaturesAndSequence(
        [Values(1, 2, 3)] int nFeatures,
        [Values(true, false)] bool hasSequence,
        [Values("BCT", "PRI", "VRL")] string division)
    {
        string text = BuildGenBank("REC001", division, nFeatures, hasSequence);
        var record = GenBankParser.Parse(text).Single();

        record.Locus.Should().Be("REC001");
        record.MoleculeType.Should().Be("DNA");
        record.Topology.Should().Be("linear");
        record.Division.Should().Be(division, "the division code is read from LOCUS");
        record.Features.Should().HaveCount(nFeatures);
        record.Features.Select(f => f.Key)
            .Should().Equal(Enumerable.Range(0, nFeatures).Select(i => GbFeatureKeys[i % GbFeatureKeys.Length]));

        if (hasSequence)
            record.Sequence.Should().Be("ACGTACGTACGTACGTACGT", "ORIGIN residues are concatenated and upper-cased");
        else
            record.Sequence.Should().BeEmpty("no ORIGIN block ⇒ no sequence");
    }

    /// <summary>
    /// Interaction witness: a feature's qualifiers are parsed; the ORIGIN block toggles whether a
    /// sequence is attached without affecting the feature table.
    /// </summary>
    [Test]
    public void ParseGenBank_QualifiersParsed_SequenceOptional()
    {
        var withSeq = GenBankParser.Parse(BuildGenBank("A", "BCT", 2, hasSequence: true)).Single();
        withSeq.Features[0].Qualifiers.Should().ContainKey("note");
        withSeq.Sequence.Should().HaveLength(20);

        var noSeq = GenBankParser.Parse(BuildGenBank("A", "BCT", 2, hasSequence: false)).Single();
        noSeq.Features.Should().HaveCount(2, "feature table is independent of ORIGIN");
        noSeq.Sequence.Should().BeEmpty();
    }
}
