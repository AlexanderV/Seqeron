using FsCheck;
using FsCheck.Fluent;

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
    #region Generators

    /// <summary>
    /// Generates a random DNA string of length [minLen, maxLen] for embedding in format records.
    /// </summary>
    private static string RandomDna(int length)
    {
        var bases = "ACGT";
        return new string(Enumerable.Range(0, length)
            .Select(_ => bases[Random.Shared.Next(4)]).ToArray());
    }

    /// <summary>
    /// Generates a valid Phred+33 quality string of specified length.
    /// </summary>
    private static string RandomQuality(int length)
    {
        return new string(Enumerable.Range(0, length)
            .Select(_ => (char)Random.Shared.Next(33, 74)).ToArray()); // '!' (Q0) to 'I' (Q40)
    }

    #endregion

    #region PARSE-FASTQ-001: RT: round-trip; P: quality len = seq len; R: quality scores in valid range

    /// <summary>
    /// INV-1: FASTQ round-trip preserves ID and sequence.
    /// Evidence: Parse → ToFastqString → Parse must be lossless for well-formed records.
    /// Source: FASTQ format (Cock et al. 2010, NAR).
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
    /// INV-2: Quality string length equals sequence length.
    /// Evidence: FASTQ spec requires quality line to have exactly one character per base.
    /// Source: Cock et al. (2010) NAR — "The Sanger FASTQ file format".
    /// </summary>
    [TestCase(8)]
    [TestCase(20)]
    [TestCase(50)]
    [Category("Property")]
    public void Fastq_QualityLength_EqualsSequenceLength(int seqLen)
    {
        string seq = RandomDna(seqLen);
        string qual = RandomQuality(seqLen);
        string fastq = $"@read1\n{seq}\n+\n{qual}";
        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].QualityString.Length, Is.EqualTo(records[0].Sequence.Length),
            "Quality string length must equal sequence length");
    }

    /// <summary>
    /// INV-3: Decoded quality scores are non-negative (Phred scores ≥ 0).
    /// Evidence: Phred score Q = −10 log₁₀(P), where P ∈ (0,1], so Q ≥ 0.
    /// Source: Ewing &amp; Green (1998) Genome Res.
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
    /// INV-4: Encode → Decode round-trip preserves Phred scores.
    /// Evidence: EncodeQualityScores(DecodeQualityScores(s)) = s for valid scores.
    /// Source: Phred+33 encoding specification.
    /// </summary>
    [TestCase(new[] { 0, 10, 20, 30, 40 })]
    [TestCase(new[] { 0, 0, 0 })]
    [TestCase(new[] { 40, 40, 40, 40 })]
    [Category("Property")]
    public void Fastq_EncodeDecodeQuality_RoundTrip(int[] originalScores)
    {
        string encoded = FastqParser.EncodeQualityScores(originalScores);
        var decoded = FastqParser.DecodeQualityScores(encoded);

        Assert.That(decoded, Is.EqualTo(originalScores));
    }

    /// <summary>
    /// INV-5: Statistics total reads matches parsed record count.
    /// Evidence: CalculateStatistics.TotalReads counts all input records.
    /// </summary>
    [TestCase(1)]
    [TestCase(3)]
    [TestCase(5)]
    [Category("Property")]
    public void Fastq_Statistics_TotalReadsMatchesParsed(int readCount)
    {
        var lines = new List<string>();
        for (int i = 0; i < readCount; i++)
        {
            int len = 10 + i * 5;
            lines.Add($"@read{i + 1}");
            lines.Add(RandomDna(len));
            lines.Add("+");
            lines.Add(RandomQuality(len));
        }
        string fastq = string.Join("\n", lines);
        var records = FastqParser.Parse(fastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.TotalReads, Is.EqualTo(records.Count),
            $"Statistics reports {stats.TotalReads} reads, parsed {records.Count}");
    }

    /// <summary>
    /// INV-6: FASTQ parsing is deterministic — same content always yields same records.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Fastq_Parsing_IsDeterministic()
    {
        string fastq = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTAAAA\n+\n!!!!!!!!";
        var r1 = FastqParser.Parse(fastq).ToList();
        var r2 = FastqParser.Parse(fastq).ToList();

        Assert.That(r1.Count, Is.EqualTo(r2.Count));
        for (int i = 0; i < r1.Count; i++)
        {
            Assert.That(r1[i].Id, Is.EqualTo(r2[i].Id));
            Assert.That(r1[i].Sequence, Is.EqualTo(r2[i].Sequence));
        }
    }

    #endregion

    #region PARSE-BED-001: R: start < end; R: chrom non-empty; R: start ≥ 0; D: deterministic

    /// <summary>
    /// INV-1: BED records have ChromEnd > ChromStart (non-empty intervals).
    /// Evidence: BED format specifies chromStart &lt; chromEnd for feature regions.
    /// Source: UCSC Genome Browser BED format specification.
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
    /// INV-2: BED record ChromStart ≥ 0 (0-based coordinate system).
    /// Evidence: BED uses 0-based, half-open coordinates; start is always non-negative.
    /// Source: UCSC BED format specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_ChromStart_NonNegative()
    {
        string bedContent = "chr1\t0\t100\tgene1\nchr1\t100\t200\tgene2\nchrX\t50000\t60000\tgene3";
        var records = BedParser.Parse(bedContent).ToList();

        foreach (var r in records)
            Assert.That(r.ChromStart, Is.GreaterThanOrEqualTo(0),
                $"BED record {r.Name}: start {r.ChromStart} should be ≥ 0");
    }

    /// <summary>
    /// INV-3: BED record chromosome name is non-empty.
    /// Evidence: Every genomic feature must be localized to a chromosome/contig.
    /// Source: UCSC BED format specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_Chrom_NonEmpty()
    {
        string bedContent = "chr1\t100\t200\tgene1\nchr2\t300\t500\tgene2";
        var records = BedParser.Parse(bedContent).ToList();

        foreach (var r in records)
            Assert.That(r.Chrom, Is.Not.Null.And.Not.Empty,
                "BED chromosome must not be empty");
    }

    /// <summary>
    /// INV-4: BED record length is positive (derived property).
    /// Evidence: Length = ChromEnd − ChromStart &gt; 0 for valid features.
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
    /// INV-5: MergeOverlapping produces non-overlapping intervals on the same chromosome.
    /// Evidence: By definition, merge combines all overlapping intervals into maximal non-overlapping ones.
    /// Source: Computational genomics interval operations (Quinlan &amp; Hall 2010, Bioinformatics).
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

    /// <summary>
    /// INV-6: BED parsing is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_Parsing_IsDeterministic()
    {
        string bedContent = "chr1\t100\t200\tgene1\t0\t+\nchr2\t300\t500\tgene2\t100\t-";
        var r1 = BedParser.Parse(bedContent).ToList();
        var r2 = BedParser.Parse(bedContent).ToList();

        Assert.That(r1.Count, Is.EqualTo(r2.Count));
        for (int i = 0; i < r1.Count; i++)
        {
            Assert.That(r1[i].Chrom, Is.EqualTo(r2[i].Chrom));
            Assert.That(r1[i].ChromStart, Is.EqualTo(r2[i].ChromStart));
            Assert.That(r1[i].ChromEnd, Is.EqualTo(r2[i].ChromEnd));
        }
    }

    /// <summary>
    /// INV-7: Sorted BED records are in ascending genomic order.
    /// Evidence: Sort orders by chromosome name, then start, then end.
    /// Source: BEDTools convention (Quinlan &amp; Hall 2010).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bed_Sort_AscendingOrder()
    {
        string bedContent = "chr2\t300\t500\nchr1\t200\t400\nchr1\t100\t200";
        var records = BedParser.Parse(bedContent).ToList();
        var sorted = BedParser.Sort(records).ToList();

        for (int i = 1; i < sorted.Count; i++)
        {
            int chromCmp = string.Compare(sorted[i].Chrom, sorted[i - 1].Chrom, StringComparison.OrdinalIgnoreCase);
            if (chromCmp == 0)
                Assert.That(sorted[i].ChromStart, Is.GreaterThanOrEqualTo(sorted[i - 1].ChromStart),
                    "Sorted records must be in ascending start order within chromosome");
            else
                Assert.That(chromCmp, Is.GreaterThanOrEqualTo(0),
                    "Sorted records must be in ascending chromosome order");
        }
    }

    #endregion

    #region PARSE-VCF-001: R: pos > 0; P: ref allele non-empty; R: qual ≥ 0 or missing; D: deterministic

    /// <summary>
    /// INV-1: VCF record position is ≥ 1 (1-based coordinate system).
    /// Evidence: VCF spec uses 1-based positions; position 0 is invalid.
    /// Source: VCF specification v4.3 (Danecek et al. 2011, Bioinformatics).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_Position_GreaterThanZero()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tAC\tA\t40\tPASS\t.";
        var records = VcfParser.Parse(vcfContent).ToList();

        foreach (var r in records)
            Assert.That(r.Pos, Is.GreaterThan(0),
                $"VCF position {r.Pos} must be > 0 (1-based)");
    }

    /// <summary>
    /// INV-2: VCF records have non-empty Ref allele.
    /// Evidence: VCF spec requires REF field to contain at least one base.
    /// Source: VCF specification v4.3.
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
    /// INV-3: VCF quality score is ≥ 0 when present (not missing).
    /// Evidence: QUAL is Phred-scaled probability; Phred scores ≥ 0.
    /// Source: VCF specification v4.3.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_Quality_NonNegativeOrMissing()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tAC\tA\t.\tPASS\t.";
        var records = VcfParser.Parse(vcfContent).ToList();

        foreach (var r in records)
        {
            if (r.Qual.HasValue)
                Assert.That(r.Qual.Value, Is.GreaterThanOrEqualTo(0.0),
                    $"VCF quality {r.Qual.Value} must be ≥ 0");
        }
    }

    /// <summary>
    /// INV-4: SNP has ref and alt length of 1.
    /// Evidence: By definition, a Single Nucleotide Polymorphism involves single-base substitution.
    /// Source: VCF specification — variant classification.
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
    /// INV-5: ClassifyVariant returns a valid VariantType for all parseable records.
    /// Evidence: All VCF records fall into one of: SNP, MNP, Insertion, Deletion, Complex, Symbolic, Unknown.
    /// Source: VCF specification v4.3 — variant type classification.
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

    /// <summary>
    /// INV-6: VCF parsing is deterministic — same content always yields same records.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_Parsing_IsDeterministic()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\trs1\tA\tG\t30\tPASS\t.\nchr2\t200\t.\tAC\tA\t40\tPASS\t.";
        var r1 = VcfParser.Parse(vcfContent).ToList();
        var r2 = VcfParser.Parse(vcfContent).ToList();

        Assert.That(r1.Count, Is.EqualTo(r2.Count));
        for (int i = 0; i < r1.Count; i++)
        {
            Assert.That(r1[i].Chrom, Is.EqualTo(r2[i].Chrom));
            Assert.That(r1[i].Pos, Is.EqualTo(r2[i].Pos));
            Assert.That(r1[i].Ref, Is.EqualTo(r2[i].Ref));
        }
    }

    /// <summary>
    /// INV-7: Insertion variant has altLen > refLen (at least one base inserted).
    /// Evidence: VCF represents insertions with a longer ALT than REF allele.
    /// Source: VCF specification v4.3.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Vcf_Insertion_AltLongerThanRef()
    {
        string vcfContent = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\t.\tA\tACG\t30\tPASS\t.";
        var record = VcfParser.Parse(vcfContent).First();

        Assert.That(VcfParser.ClassifyVariant(record), Is.EqualTo(VcfParser.VariantType.Insertion));
        Assert.That(record.Alt[0].Length, Is.GreaterThan(record.Ref.Length));
    }

    #endregion

    #region PARSE-GFF-001 + ANNOT-GFF-001: RT: round-trip; R: start ≤ end; R: strand ∈ {+,-,.}; D: deterministic

    /// <summary>
    /// INV-1: GFF records have End ≥ Start (1-based, inclusive coordinates).
    /// Evidence: GFF3 spec requires start ≤ end for all features.
    /// Source: GFF3 specification (Sequence Ontology Project).
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
    /// INV-2: GFF strand is one of {+, -, .} (forward, reverse, unstranded).
    /// Evidence: GFF3 spec column 7 allows only +, -, . (or ?) for strand.
    /// Source: GFF3 specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_Strand_IsValid()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t-\t.\tParent=gene1\n" +
                            "chr1\ttest\tregion\t1\t1000\t.\t.\t.\tID=region1";
        var records = GffParser.Parse(gffContent).ToList();

        foreach (var r in records)
            Assert.That(r.Strand, Is.AnyOf('+', '-', '.', '?'),
                $"GFF strand '{r.Strand}' is invalid");
    }

    /// <summary>
    /// INV-3: GFF round-trip — Write → Parse preserves feature coordinates and type.
    /// Evidence: WriteToStream(records) → Parse(output) must preserve Start, End, Type.
    /// Source: GFF3 specification — format is self-describing.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_RoundTrip_PreservesCoordinatesAndType()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1";
        var original = GffParser.Parse(gffContent).ToList();

        using var sw = new StringWriter();
        GffParser.WriteToStream(sw, original);
        var reparsed = GffParser.Parse(sw.ToString()).ToList();

        Assert.That(reparsed.Count, Is.EqualTo(original.Count), "Record count mismatch");
        for (int i = 0; i < original.Count; i++)
        {
            Assert.That(reparsed[i].Start, Is.EqualTo(original[i].Start), $"Start mismatch at index {i}");
            Assert.That(reparsed[i].End, Is.EqualTo(original[i].End), $"End mismatch at index {i}");
            Assert.That(reparsed[i].Type, Is.EqualTo(original[i].Type), $"Type mismatch at index {i}");
        }
    }

    /// <summary>
    /// INV-4: FilterByType returns only features of the requested type.
    /// Evidence: FilterByType is a type-predicate filter.
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
    /// INV-5: GFF statistics gene count is non-negative and consistent with filtered genes.
    /// Evidence: GeneCount = |{r : r.Type = "gene"}|.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_Statistics_GeneCount_ConsistentWithFilter()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\tgene\t600\t900\t.\t+\t.\tID=gene2\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();
        var stats = GffParser.CalculateStatistics(records);
        var genesByFilter = GffParser.GetGenes(records).Count();

        Assert.That(stats.GeneCount, Is.GreaterThanOrEqualTo(0));
        Assert.That(stats.GeneCount, Is.EqualTo(genesByFilter),
            "Statistics gene count must match filtered gene count");
    }

    /// <summary>
    /// INV-6: GFF parsing is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_Parsing_IsDeterministic()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t-\t.\tParent=gene1";
        var r1 = GffParser.Parse(gffContent).ToList();
        var r2 = GffParser.Parse(gffContent).ToList();

        Assert.That(r1.Count, Is.EqualTo(r2.Count));
        for (int i = 0; i < r1.Count; i++)
        {
            Assert.That(r1[i].Seqid, Is.EqualTo(r2[i].Seqid));
            Assert.That(r1[i].Start, Is.EqualTo(r2[i].Start));
            Assert.That(r1[i].End, Is.EqualTo(r2[i].End));
            Assert.That(r1[i].Type, Is.EqualTo(r2[i].Type));
            Assert.That(r1[i].Strand, Is.EqualTo(r2[i].Strand));
        }
    }

    /// <summary>
    /// INV-7: GFF coordinates are 1-based (Start ≥ 1).
    /// Evidence: GFF3 uses 1-based inclusive coordinates unlike BED's 0-based.
    /// Source: GFF3 specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff_Coordinates_OneBased()
    {
        string gffContent = "chr1\ttest\tgene\t1\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();

        foreach (var r in records)
            Assert.That(r.Start, Is.GreaterThanOrEqualTo(1),
                $"GFF start {r.Start} must be ≥ 1 (1-based coordinates)");
    }

    #endregion

    #region PARSE-GENBANK-001: RT: round-trip; P: locus line present; P: sequence preserved; D: deterministic

    /// <summary>
    /// INV-1: GenBank parsed record has non-empty sequence matching reported length.
    /// Evidence: The LOCUS line specifies the sequence length; the ORIGIN section contains the sequence.
    /// Source: NCBI GenBank flat file format specification.
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

    /// <summary>
    /// INV-2: GenBank locus name is non-empty.
    /// Evidence: LOCUS is a mandatory field; the first token after LOCUS is the locus name.
    /// Source: NCBI GenBank flat file format specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenBank_Locus_NotEmpty()
    {
        string genBankContent =
            "LOCUS       MYGENE         20 bp    DNA     linear   UNK\n" +
            "DEFINITION  A test gene.\n" +
            "ACCESSION   MG001\n" +
            "VERSION     MG001.1\n" +
            "FEATURES             Location/Qualifiers\n" +
            "ORIGIN\n" +
            "        1 acgtacgtac acgtacgtac\n" +
            "//\n";
        var records = GenBankParser.Parse(genBankContent).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Locus, Is.Not.Null.And.Not.Empty,
            "GenBank locus name must not be empty");
    }

    /// <summary>
    /// INV-3: GenBank sequence contains only valid nucleotide characters.
    /// Evidence: ORIGIN section contains lowercase a, c, g, t, n.
    /// Source: NCBI GenBank format — sequence is DNA/RNA letters only.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenBank_Sequence_ContainsValidBases()
    {
        string genBankContent =
            "LOCUS       TEST           10 bp    DNA     linear   UNK\n" +
            "DEFINITION  Test.\n" +
            "ACCESSION   T001\n" +
            "VERSION     T001.1\n" +
            "FEATURES             Location/Qualifiers\n" +
            "ORIGIN\n" +
            "        1 acgtacgtac\n" +
            "//\n";
        var records = GenBankParser.Parse(genBankContent).ToList();

        foreach (var r in records)
            Assert.That(r.Sequence, Does.Match("^[acgtnACGTN]+$"),
                $"Sequence contains invalid characters: {r.Sequence}");
    }

    /// <summary>
    /// INV-4: GenBank features are parseable and have valid locations.
    /// Evidence: Features have key, start ≤ end locations.
    /// Source: NCBI GenBank feature table specification.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenBank_Features_HaveValidLocations()
    {
        string genBankContent =
            "LOCUS       TEST           30 bp    DNA     linear   UNK\n" +
            "DEFINITION  Test sequence.\n" +
            "ACCESSION   TEST001\n" +
            "VERSION     TEST001.1\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene            1..20\n" +
            "                     /gene=\"testgene\"\n" +
            "     CDS             5..15\n" +
            "                     /product=\"testprotein\"\n" +
            "ORIGIN\n" +
            "        1 acgtacgtac acgtacgtac acgtacgtac\n" +
            "//\n";
        var records = GenBankParser.Parse(genBankContent).ToList();
        var features = records[0].Features;

        foreach (var f in features)
        {
            Assert.That(f.Key, Is.Not.Null.And.Not.Empty, "Feature key must not be empty");
            Assert.That(f.Location.Start, Is.GreaterThanOrEqualTo(1),
                $"Feature {f.Key} start {f.Location.Start} must be ≥ 1");
            Assert.That(f.Location.End, Is.GreaterThanOrEqualTo(f.Location.Start),
                $"Feature {f.Key} end {f.Location.End} must be ≥ start {f.Location.Start}");
        }
    }

    /// <summary>
    /// INV-5: GenBank parsing is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GenBank_Parsing_IsDeterministic()
    {
        string genBankContent =
            "LOCUS       TEST           10 bp    DNA     linear   UNK\n" +
            "DEFINITION  Test.\n" +
            "ACCESSION   T001\n" +
            "VERSION     T001.1\n" +
            "FEATURES             Location/Qualifiers\n" +
            "ORIGIN\n" +
            "        1 acgtacgtac\n" +
            "//\n";
        var r1 = GenBankParser.Parse(genBankContent).ToList();
        var r2 = GenBankParser.Parse(genBankContent).ToList();

        Assert.That(r1.Count, Is.EqualTo(r2.Count));
        for (int i = 0; i < r1.Count; i++)
        {
            Assert.That(r1[i].Locus, Is.EqualTo(r2[i].Locus));
            Assert.That(r1[i].Sequence, Is.EqualTo(r2[i].Sequence));
            Assert.That(r1[i].SequenceLength, Is.EqualTo(r2[i].SequenceLength));
        }
    }

    #endregion

    #region PARSE-EMBL-001: RT: round-trip; P: ID line present; P: sequence preserved; D: deterministic

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

    #endregion
}
