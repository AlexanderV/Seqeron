namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for file I/O parsers.
///
/// Test Units: PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001, PARSE-GFF-001,
///             PARSE-GENBANK-001, PARSE-EMBL-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("IO")]
public class FileIOSnapshotTests
{
    [Test]
    public Task FastqParse_Snapshot()
    {
        string fastqContent = "@SEQ1 Human sample\nACGTACGTACGT\n+\nIIIIIIIIIIII\n" +
                              "@SEQ2 Mouse sample\nTTTTAAAACCCC\n+\n!!!!!!!!!!!!";
        var records = FastqParser.Parse(fastqContent)
            .Select(r => new { r.Id, r.Description, r.Sequence, r.QualityString })
            .ToList();

        return Verify(new { Records = records });
    }

    [Test]
    public Task FastqStatistics_Snapshot()
    {
        string fastqContent = "@SEQ1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n" +
                              "@SEQ2\nTTTTAAAA\n+\n!!!!!!!!";
        var records = FastqParser.Parse(fastqContent).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        return Verify(new
        {
            stats.TotalReads,
            stats.TotalBases,
            stats.MeanReadLength,
            stats.MeanQuality,
            stats.MinReadLength,
            stats.MaxReadLength
        });
    }

    [Test]
    public Task BedParse_Snapshot()
    {
        string bedContent = "chr1\t100\t200\tgene1\t500\t+\n" +
                            "chr1\t300\t500\tgene2\t300\t-\n" +
                            "chr2\t0\t1000\tgene3\t100\t+";
        var records = BedParser.Parse(bedContent)
            .Select(r => new { r.Chrom, r.ChromStart, r.ChromEnd, r.Name, r.Strand, r.Length })
            .ToList();

        return Verify(new { Records = records });
    }

    [Test]
    public Task VcfParse_Snapshot()
    {
        string vcfContent = "##fileformat=VCFv4.2\n" +
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                            "chr1\t100\trs123\tA\tG\t30\tPASS\tDP=50\n" +
                            "chr1\t200\t.\tAC\tA\t40\tPASS\tDP=30";
        var records = VcfParser.Parse(vcfContent)
            .Select(r => new
            {
                r.Chrom,
                r.Pos,
                r.Id,
                r.Ref,
                Alt = string.Join(",", r.Alt),
                Type = VcfParser.ClassifyVariant(r).ToString()
            })
            .ToList();

        return Verify(new { Records = records });
    }

    [Test]
    public Task GffParse_Snapshot()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=TestGene\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1\n" +
                            "chr1\ttest\texon\t300\t500\t.\t+\t.\tParent=gene1";
        var records = GffParser.Parse(gffContent)
            .Select(r => new { r.Seqid, r.Source, r.Type, r.Start, r.End, r.Strand })
            .ToList();

        return Verify(new { Records = records });
    }

    [Test]
    public Task GffStatistics_Snapshot()
    {
        string gffContent = "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n" +
                            "chr1\ttest\texon\t100\t200\t.\t+\t.\tParent=gene1\n" +
                            "chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=gene1";
        var records = GffParser.Parse(gffContent).ToList();
        var stats = GffParser.CalculateStatistics(records);

        return Verify(new
        {
            stats.TotalFeatures,
            stats.GeneCount,
            stats.ExonCount
        });
    }
}
