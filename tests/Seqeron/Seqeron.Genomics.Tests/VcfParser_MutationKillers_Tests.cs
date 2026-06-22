using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;
using static Seqeron.Genomics.IO.VcfParser;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Targeted mutation-killing tests for VcfParser.cs (checklist 04 row 67, PARSE-VCF-001).
/// The canonical suite left the VCFv4 header metadata parsing, the per-sample FORMAT/genotype
/// decoding, the variant-type classification, the allele-frequency / Ti-Tv statistics and the
/// writer under-pinned. These pin the exact VCF spec behaviour so the boundary/logical/
/// null-coalescing mutants diverge.
/// </summary>
[TestFixture]
public class VcfParser_MutationKillers_Tests
{
    private const string Vcf =
        "##fileformat=VCFv4.2\n" +
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
        "##FILTER=<ID=q10,Description=\"Quality below 10\">\n" +
        "##contig=<ID=chr1>\n" +
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n" +
        "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30\tGT:DP\t0/1:15\t1/1:20\n" +
        "chr1\t200\t.\tAG\tA\t.\tq10\tDP=10\tGT:DP\t0/0:5\t0/1:8\n";

    // ── ParseWithHeader: header metadata + data records parsed exactly ─────────────────

    [Test]
    public void ParseWithHeader_ParsesMetadataAndRecordsExactly()
    {
        var (header, recordsEnum) = ParseWithHeader(Vcf);
        var records = recordsEnum.ToList();

        header.FileFormat.Should().Be("VCFv4.2");
        header.SampleNames.Should().Equal("S1", "S2");

        header.InfoFields.Should().ContainSingle();
        header.InfoFields[0].Id.Should().Be("DP");
        header.InfoFields[0].Number.Should().Be("1");
        header.InfoFields[0].Type.Should().Be("Integer");
        header.InfoFields[0].Description.Should().Be("Total Depth"); // embedded-quote/comma split

        header.FormatFields.Should().ContainSingle();
        header.FormatFields[0].Id.Should().Be("GT");

        header.FilterFields.Should().ContainSingle();
        header.FilterFields[0].Id.Should().Be("q10");
        header.FilterFields[0].Description.Should().Be("Quality below 10");

        header.OtherMetadata.Should().ContainKey("contig");
        header.OtherMetadata["contig"].Should().Be("<ID=chr1>");

        records.Should().HaveCount(2);
        var r1 = records[0];
        r1.Chrom.Should().Be("chr1");
        r1.Pos.Should().Be(100);
        r1.Id.Should().Be("rs1");
        r1.Ref.Should().Be("A");
        r1.Alt.Should().Equal("G");
        r1.Qual.Should().Be(50.0);
        r1.Filter.Should().Equal("PASS");
        r1.Info["DP"].Should().Be("30");
        r1.Format.Should().Equal("GT", "DP");
        r1.Samples.Should().HaveCount(2);
        r1.Samples![0]["GT"].Should().Be("0/1");
        r1.Samples![0]["DP"].Should().Be("15");
        r1.Samples![1]["GT"].Should().Be("1/1");

        records[1].Qual.Should().BeNull();   // "." ⇒ null
        records[1].Filter.Should().Equal("q10");
    }

    // ── ClassifyVariant: SNP / MNP / Insertion / Deletion / Complex / Symbolic ─────────

    [Test]
    public void ClassifyVariant_DistinguishesAllTypes()
    {
        VcfRecord V(string r, string alt) => new("c", 1, ".", r, alt.Split(','), null,
            System.Array.Empty<string>(), new Dictionary<string, string>());

        ClassifyVariant(V("A", "G")).Should().Be(VcfParser.VariantType.SNP);
        ClassifyVariant(V("AC", "GT")).Should().Be(VcfParser.VariantType.MNP);
        ClassifyVariant(V("A", "ACGT")).Should().Be(VcfParser.VariantType.Insertion);
        ClassifyVariant(V("AGT", "A")).Should().Be(VcfParser.VariantType.Deletion);
        ClassifyVariant(V("AG", "TCC")).Should().Be(VcfParser.VariantType.Complex);
        ClassifyVariant(V("A", "<DEL>")).Should().Be(VcfParser.VariantType.Symbolic);
        ClassifyVariant(V("A", "*")).Should().Be(VcfParser.VariantType.Symbolic);
    }

    [Test]
    public void GetVariantLength_IsAbsoluteRefAltLengthDifference()
    {
        var del = new VcfRecord("c", 1, ".", "AGT", new[] { "A" }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());
        GetVariantLength(del).Should().Be(2);
    }

    [Test]
    public void ClassifyVariant_AltIndexOutOfRange_IsUnknown()
    {
        var rec = new VcfRecord("c", 1, ".", "A", new[] { "G" }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());
        // altIndex == Alt.Length must be rejected (kills `altIndex >= Alt.Length` → `>`,
        // which would index past the array).
        ClassifyVariant(rec, 1).Should().Be(VcfParser.VariantType.Unknown);
        GetVariantLength(rec, 1).Should().Be(0);
    }

    [Test]
    public void IsIndel_TrueForInsertionsAndDeletionsOnly()
    {
        var del = new VcfRecord("c", 1, ".", "AGT", new[] { "A" }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());
        var snp = new VcfRecord("c", 1, ".", "A", new[] { "G" }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());
        IsIndel(del).Should().BeTrue();
        IsIndel(snp).Should().BeFalse();
    }

    // ── Genotype zygosity ──────────────────────────────────────────────────────────────

    [Test]
    public void GenotypeZygosity_HetHomAltHomRef()
    {
        var r1 = ParseWithHeader(Vcf).Records.First();
        IsHet(r1, 0).Should().BeTrue();      // S1 = 0/1
        IsHomAlt(r1, 0).Should().BeFalse();
        IsHomAlt(r1, 1).Should().BeTrue();   // S2 = 1/1
        IsHomRef(r1, 1).Should().BeFalse();

        var r2 = ParseWithHeader(Vcf).Records.Last();
        IsHomRef(r2, 0).Should().BeTrue();   // S1 = 0/0
        IsHet(r2, 1).Should().BeTrue();      // S2 = 0/1
        // 0/0 is homozygous REFERENCE, not alt (the allele[0] != "0" clause must hold).
        IsHomAlt(r2, 0).Should().BeFalse();
    }

    [Test]
    public void IsHomAlt_RequiresExactlyDiploid()
    {
        // A triploid 1/1/1 call is not classified hom-alt: the `alleles.Length == 2` clause must
        // hold (kills `Length == 2 && [0]==[1]` → `||`, which would accept any all-equal ploidy).
        var fmt = new[] { "GT" };
        var samples = new IReadOnlyDictionary<string, string>[]
        {
            new Dictionary<string, string> { ["GT"] = "1/1/1" },
        };
        var rec = new VcfRecord("c", 1, ".", "A", new[] { "G" }, null, new[] { "PASS" },
            new Dictionary<string, string>(), fmt, samples);

        IsHomAlt(rec, 0).Should().BeFalse();
    }

    [Test]
    public void GetGenotype_SampleIndexEqualToCount_IsNull()
    {
        var r1 = ParseWithHeader(Vcf).Records.First(); // 2 samples
        // index == Samples.Count must be rejected (kills `sampleIndex >= Count` → `>`, which would
        // index past the sample list).
        GetGenotype(r1, 2).Should().BeNull();
    }

    [Test]
    public void CalculateTiTvRatio_NoTransversions_IsNull()
    {
        VcfRecord Snp(string r, string a) => new("c", 1, ".", r, new[] { a }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());

        // Only transitions (A→G, C→T) ⇒ transversions == 0 ⇒ ratio undefined/null (kills
        // `transversions > 0` → `>= 0`, which would divide by zero).
        CalculateTiTvRatio(new[] { Snp("A", "G"), Snp("C", "T") }).Should().BeNull();
    }

    [Test]
    public void FilterByRegion_ExcludesOtherChromosomesInRange()
    {
        var recs = new[]
        {
            new VcfRecord("chr1", 150, ".", "A", new[] { "G" }, null, System.Array.Empty<string>(),
                new Dictionary<string, string>()),
            new VcfRecord("chr2", 150, ".", "A", new[] { "G" }, null, System.Array.Empty<string>(),
                new Dictionary<string, string>()),
        };
        // chr2:150 is within [100,200] positionally but on the wrong chromosome (kills the
        // `Chrom.Equals && Pos in range` → `||` mutant).
        FilterByRegion(recs, "chr1", 100, 200).Select(r => r.Chrom).Should().Equal("chr1");
    }

    [Test]
    public void GetReadDepth_ReadsDpFormatField()
    {
        var r1 = ParseWithHeader(Vcf).Records.First();
        GetReadDepth(r1, 0).Should().Be(15);
        GetReadDepth(r1, 1).Should().Be(20);
    }

    // ── Allele frequency: counts ALT (allele index altIndex+1) over called alleles ─────

    [Test]
    public void CalculateAlleleFrequency_CountsAltAllelesOverTotal()
    {
        var r1 = ParseWithHeader(Vcf).Records.First();
        // S1 0/1 → one ALT of two; S2 1/1 → two ALT of two ⇒ 3 / 4.
        CalculateAlleleFrequency(new[] { r1 }).Should().BeApproximately(0.75, 1e-9);
    }

    [Test]
    public void CalculateAlleleFrequency_SkipsMissingGenotypes()
    {
        // S1 0/1 contributes; S2 ./. is missing and must be excluded entirely (kills the
        // `gt == null || gt.Contains('.')` skip → `&&`, which would count the missing alleles).
        var fmt = new[] { "GT" };
        var samples = new IReadOnlyDictionary<string, string>[]
        {
            new Dictionary<string, string> { ["GT"] = "0/1" },
            new Dictionary<string, string> { ["GT"] = "./." },
        };
        var rec = new VcfRecord("c", 1, ".", "A", new[] { "G" }, null,
            new[] { "PASS" }, new Dictionary<string, string>(), fmt, samples);

        // Only the 0/1 sample counts: 1 ALT of 2 alleles ⇒ 0.5.
        CalculateAlleleFrequency(new[] { rec }).Should().BeApproximately(0.5, 1e-9);
    }

    [Test]
    public void CalculateStatistics_NoQualityRecords_MeanQualityIsNull()
    {
        // No record carries QUAL ⇒ qualities list empty ⇒ MeanQuality null (kills `Count > 0`
        // → `>= 0`, which would average an empty list).
        var recs = new[]
        {
            new VcfRecord("c", 1, ".", "A", new[] { "G" }, null, new[] { "PASS" },
                new Dictionary<string, string>()),
        };
        CalculateStatistics(recs).MeanQuality.Should().BeNull();
    }

    // ── Statistics & Ti/Tv ──────────────────────────────────────────────────────────────

    [Test]
    public void CalculateStatistics_CountsTypesPassingAndMeanQuality()
    {
        var stats = CalculateStatistics(ParseWithHeader(Vcf).Records);

        stats.TotalVariants.Should().Be(2);
        stats.SnpCount.Should().Be(1);
        stats.IndelCount.Should().Be(1);     // the AG→A deletion
        stats.ComplexCount.Should().Be(0);
        stats.PassingCount.Should().Be(1);   // only the PASS record
        stats.MeanQuality.Should().Be(50.0); // only one record carries QUAL
    }

    [Test]
    public void CalculateTiTvRatio_DividesTransitionsByTransversions()
    {
        VcfRecord Snp(string r, string a) => new("c", 1, ".", r, new[] { a }, null,
            System.Array.Empty<string>(), new Dictionary<string, string>());

        // A→G is a transition; A→C is a transversion ⇒ ratio = 1/1 = 1.0.
        CalculateTiTvRatio(new[] { Snp("A", "G"), Snp("A", "C") }).Should().Be(1.0);
    }

    // ── Filters: inclusive POS region and quality threshold ────────────────────────────

    [Test]
    public void FilterByRegion_PositionBoundsAreInclusive()
    {
        var recs = ParseWithHeader(Vcf).Records.ToList();
        // Query [100,200] must keep BOTH boundary positions (kills `>= start` → `> start`
        // and `<= end` → `< end`).
        FilterByRegion(recs, "chr1", 100, 200).Select(r => r.Pos).Should().Equal(100, 200);
    }

    [Test]
    public void FilterByQuality_ThresholdIsInclusive()
    {
        var recs = ParseWithHeader(Vcf).Records.ToList();
        // minQuality 50: the 50.0 record is kept (>= boundary); the QUAL-less record dropped.
        FilterByQuality(recs, 50).Select(r => r.Pos).Should().Equal(100);
    }

    [Test]
    public void FilterPassing_RequiresExactlyPass()
    {
        var recs = ParseWithHeader(Vcf).Records.ToList();
        FilterPassing(recs).Select(r => r.Pos).Should().Equal(100);
    }

    // ── Round-trip write: header + sample columns reconstruct the records ──────────────

    [Test]
    public void WriteThenParse_WithSamples_RoundTripsRecords()
    {
        var (header, recordsEnum) = ParseWithHeader(Vcf);
        var records = recordsEnum.ToList();

        using var sw = new StringWriter();
        WriteToStream(sw, records, header, header.SampleNames.ToArray());
        var text = sw.ToString();

        text.Should().Contain("##fileformat=VCFv4.2");
        // The metadata-emitting loops must run (kills their block-removal mutants):
        text.Should().Contain("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
        text.Should().Contain("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        text.Should().Contain("##FILTER=<ID=q10,Description=\"Quality below 10\">");
        text.Should().Contain("\tFORMAT\tS1\tS2");
        // Exact data line incl. FORMAT + both sample columns (kills the sample-emitting logic):
        text.Should().Contain("chr1\t100\trs1\tA\tG\t50.00\tPASS\tDP=30\tGT:DP\t0/1:15\t1/1:20");

        var (_, reparsedEnum) = ParseWithHeader(text);
        var reparsed = reparsedEnum.ToList();

        reparsed.Should().HaveCount(2);
        reparsed[0].Chrom.Should().Be("chr1");
        reparsed[0].Pos.Should().Be(100);
        reparsed[0].Ref.Should().Be("A");
        reparsed[0].Alt.Should().Equal("G");
        reparsed[0].Samples![0]["GT"].Should().Be("0/1");
        reparsed[0].Samples![1]["GT"].Should().Be("1/1");
        reparsed[1].Pos.Should().Be(200);
        reparsed[1].Samples![0]["GT"].Should().Be("0/0");
    }
}
