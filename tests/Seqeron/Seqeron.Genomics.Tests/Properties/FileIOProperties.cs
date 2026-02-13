namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for file I/O parsers: FASTQ, BED, VCF, GFF, GenBank, EMBL.
///
/// Test Units: PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001, PARSE-GFF-001,
///             PARSE-GENBANK-001, PARSE-EMBL-001, ANNOT-GFF-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("IO")]
public class FileIOProperties
{
    // -- PARSE-FASTQ-001 --

    /// <summary>
    /// FASTQ round-trip: parse → toFastqString → parse preserves ID and sequence.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_RoundTrip_PreservesIdAndSequence()
    {
        string fastqContent = "@SEQ1 description\nACGTACGT\n+\nIIIIIIII\n@SEQ2\nTTTTAAAA\n+\n!!!!!!!!";
        var records = FastqParser.Parse(fastqContent).ToList();

        foreach (var r in records)
        {
            string serialized = FastqParser.ToFastqString(r);
            var reparsed = FastqParser.Parse(serialized).First();

            Assert.That(reparsed.Id, Is.EqualTo(r.Id), "ID mismatch after round-trip");
            Assert.That(reparsed.Sequence, Is.EqualTo(r.Sequence), "Sequence mismatch after round-trip");
        }
    }

    /// <summary>
    /// Quality string length equals sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_QualityLength_EqualsSequenceLength()
    {
        string fastqContent = "@SEQ1\nACGTACGT\n+\nIIIIIIII\n@SEQ2\nTTTTAAAACC\n+\n!!!!!!!!!I";
        var records = FastqParser.Parse(fastqContent).ToList();

        foreach (var r in records)
            Assert.That(r.QualityString.Length, Is.EqualTo(r.Sequence.Length),
                $"Quality length {r.QualityString.Length} ≠ sequence length {r.Sequence.Length} for {r.Id}");
    }

    /// <summary>
    /// Decoded quality scores are non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_DecodedScores_NonNegative()
    {
        string qualityString = "IIIIIII!@#$%^&";
        var scores = FastqParser.DecodeQualityScores(qualityString);

        foreach (var s in scores)
            Assert.That(s, Is.GreaterThanOrEqualTo(0), $"Quality score {s} is negative");
    }

    /// <summary>
    /// Encode → Decode round-trip preserves scores.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_EncodeDecodeQuality_RoundTrip()
    {
        var originalScores = new[] { 0, 10, 20, 30, 40 };
        string encoded = FastqParser.EncodeQualityScores(originalScores);
        var decoded = FastqParser.DecodeQualityScores(encoded);

        Assert.That(decoded, Is.EqualTo(originalScores));
    }

    /// <summary>
    /// Statistics total reads matches parsed record count.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_Statistics_TotalReadsMatchesParsed()
    {
        string fastqContent = "@SEQ1\nACGT\n+\nIIII\n@SEQ2\nTTTT\n+\n!!!!";
        var records = FastqParser.Parse(fastqContent).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.TotalReads, Is.EqualTo(records.Count));
    }

    // -- PARSE-BED-001 --

    /// <summary>
    /// BED records: ChromEnd > ChromStart (non-empty intervals).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_ChromEnd_GreaterThanStart()
    {
        string bedContent = "chr1\t100\t200\tgene1\nchr1\t300\t500\tgene2\nchr2\t0\t1000\tgene3";
        var records = BedParser.Parse(bedContent).ToList();

        foreach (var r in records)
            Assert.That(r.ChromEnd, Is.GreaterThan(r.ChromStart),
                $"BED record {r.Name}: end {r.ChromEnd} should be > start {r.ChromStart}");
    }

    /// <summary>
    /// BED record length is positive.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_Length_IsPositive()
    {
        string bedContent = "chr1\t100\t200\tgene1\nchr1\t300\t500\tgene2";
        var records = BedParser.Parse(bedContent).ToList();

        foreach (var r in records)
            Assert.That(r.Length, Is.GreaterThan(0));
    }

    /// <summary>
    /// MergeOverlapping produces non-overlapping intervals.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_MergeOverlapping_NoOverlaps()
    {
        string bedContent = "chr1\t100\t300\nchr1\t200\t400\nchr1\t500\t600";
        var records = BedParser.Parse(bedContent).ToList();
        var merged = BedParser.MergeOverlapping(records).ToList();

        for (int i = 1; i < merged.Count; i++)
        {
            if (merged[i].Chrom == merged[i - 1].Chrom)
                Assert.That(merged[i].ChromStart, Is.GreaterThanOrEqualTo(merged[i - 1].ChromEnd),
                    "Merged intervals should not overlap");
        }
    }

    // -- PARSE-VCF-001 --

    /// <summary>
    /// VCF records have non-empty Ref allele.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_RefAllele_NotEmpty()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tAC\tA\t40\tPASS\t.";
        var records = VcfParser.Parse(vcfContent).ToList();

        foreach (var r in records)
            Assert.That(r.Ref, Is.Not.Null.And.Not.Empty,
                $"VCF record at {r.Chrom}:{r.Pos} has empty Ref");
    }

    /// <summary>
    /// SNP has ref and alt length of 1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_IsSNP_SingleBaseRefAndAlt()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tG\t30\tPASS\t.";
        var record = VcfParser.Parse(vcfContent).First();

        Assert.That(VcfParser.IsSNP(record), Is.True);
        Assert.That(record.Ref.Length, Is.EqualTo(1));
        Assert.That(record.Alt[0].Length, Is.EqualTo(1));
    }

    /// <summary>
    /// ClassifyVariant returns a valid VariantType.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_ClassifyVariant_ReturnsValidType()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tAC\tA\t40\tPASS\t.";
        var records = VcfParser.Parse(vcfContent).ToList();

        foreach (var r in records)
        {
            var vtype = VcfParser.ClassifyVariant(r);
            Assert.That(vtype, Is.TypeOf<VcfParser.VariantType>());
        }
    }

    // -- PARSE-GFF-001 + ANNOT-GFF-001 --

    /// <summary>
    /// GFF records: End ≥ Start.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_End_GreaterOrEqualStart()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1\n" +
                            "chr1\ttest\texon\t300\t500\t.\t+\t.\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();

        foreach (var r in records)
            Assert.That(r.End, Is.GreaterThanOrEqualTo(r.Start),
                $"GFF record {r.Type}: end {r.End} < start {r.Start}");
    }

    /// <summary>
    /// FilterByType returns only features of the requested type.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_FilterByType_ReturnsOnlyRequestedType()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1\n" +
                            "chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();
        var exons = GffParser.FilterByType(records, "exon").ToList();

        foreach (var e in exons)
            Assert.That(e.Type, Is.EqualTo("exon"));
    }

    /// <summary>
    /// GFF statistics gene count is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_Statistics_GeneCount_NonNegative()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();
        var stats = GffParser.CalculateStatistics(records);

        Assert.That(stats.GeneCount, Is.GreaterThanOrEqualTo(0));
    }

    // -- PARSE-GENBANK-001 --

    /// <summary>
    /// GenBank parsed record has non-empty sequence matching reported length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenBank_SequenceLength_MatchesHeader()
    {
        string genBankContent =
            "LOCUS       TEST           10 bp    DNA     linear   UNK\n" +
            "DEFINITION  Test sequence.\n" +
            "ACCESSION   TEST001\n" +
            "VERSION     TEST001.1\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene            1..10\n" +
            "                     /gene=\"test\"\n" +
            "ORIGIN\n" +
            "        1 acgtacgtac\n" +
            "//\n";
        var records = GenBankParser.Parse(genBankContent).ToList();

        foreach (var r in records)
        {
            Assert.That(r.Sequence, Is.Not.Null.And.Not.Empty, "GenBank sequence should not be empty");
            Assert.That(r.Sequence.Length, Is.EqualTo(r.SequenceLength),
                $"Sequence length {r.Sequence.Length} ≠ reported {r.SequenceLength}");
        }
    }

    // -- PARSE-EMBL-001 --

    /// <summary>
    /// EMBL parsed record has non-empty sequence.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Embl_ParsedRecord_HasSequence()
    {
        string emblContent =
            "ID   TEST001; SV 1; linear; genomic DNA; STD; UNK; 10 BP.\n" +
            "XX\n" +
            "AC   TEST001;\n" +
            "XX\n" +
            "DE   Test sequence.\n" +
            "XX\n" +
            "SQ   Sequence 10 BP;\n" +
            "     acgtacgtac                                                           10\n" +
            "//\n";
        var records = EmblParser.Parse(emblContent).ToList();

        foreach (var r in records)
            Assert.That(r.Sequence, Is.Not.Null.And.Not.Empty, "EMBL sequence should not be empty");
    }
}
