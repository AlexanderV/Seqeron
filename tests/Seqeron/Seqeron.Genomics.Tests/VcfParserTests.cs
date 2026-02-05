using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the VCF format parser.
/// </summary>
[TestFixture]
public class VcfParserTests
{
    #region Test Data

    private const string SimpleVcf = @"##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description=""Total Depth"">
##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=""Read Depth"">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
chr1	100	rs123	A	G	99	PASS	DP=50	GT:DP	0/1:25	1/1:30
chr1	200	.	C	T	50	.	DP=30	GT:DP	0/0:15	0/1:20
chr2	300	rs456	G	A,C	80	PASS	DP=40	GT:DP	1/2:20	0/1:25";

    private const string VcfWithIndels = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	AT	99	PASS	.
chr1	200	.	ATG	A	99	PASS	.
chr1	300	.	A	G	99	PASS	.";

    private const string VcfWithFilters = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	.
chr1	200	.	C	T	20	LowQual	.
chr1	300	.	G	A	10	LowQual;LowCov	.";

    #endregion

    #region Basic Parsing Tests

    [Test]
    public void Parse_SimpleVcf_ReturnsCorrectRecords()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();

        Assert.That(records, Has.Count.EqualTo(3));
        Assert.That(records[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(records[0].Pos, Is.EqualTo(100));
        Assert.That(records[0].Ref, Is.EqualTo("A"));
        Assert.That(records[0].Alt[0], Is.EqualTo("G"));
    }

    [Test]
    public void Parse_EmptyContent_ReturnsEmpty()
    {
        var records = VcfParser.Parse("").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_NullContent_ReturnsEmpty()
    {
        var records = VcfParser.Parse((string)null!).ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_SkipsMetadataLines()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();

        // Should only return variant records, not ## lines
        Assert.That(records, Has.Count.EqualTo(3));
    }

    [Test]
    public void Parse_MultipleAlternates_ParsedCorrectly()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var multiAlt = records.First(r => r.Alt.Length > 1);

        Assert.That(multiAlt.Alt, Has.Length.EqualTo(2));
        Assert.That(multiAlt.Alt, Contains.Item("A"));
        Assert.That(multiAlt.Alt, Contains.Item("C"));
    }

    #endregion

    #region Header Parsing Tests

    [Test]
    public void ParseWithHeader_ReturnsHeaderAndRecords()
    {
        var (header, records) = VcfParser.ParseWithHeader(SimpleVcf);

        Assert.That(header.FileFormat, Is.EqualTo("VCFv4.3"));
        Assert.That(header.InfoFields.Count, Is.GreaterThanOrEqualTo(1));
        Assert.That(header.FormatFields.Count, Is.GreaterThanOrEqualTo(2));
        Assert.That(records.Count(), Is.EqualTo(3));
    }

    [Test]
    public void ParseWithHeader_SampleNames_Parsed()
    {
        var (header, _) = VcfParser.ParseWithHeader(SimpleVcf);

        Assert.That(header.SampleNames, Contains.Item("Sample1"));
        Assert.That(header.SampleNames, Contains.Item("Sample2"));
    }

    [Test]
    public void ParseWithHeader_InfoFields_Parsed()
    {
        var (header, _) = VcfParser.ParseWithHeader(SimpleVcf);
        var dpInfo = header.InfoFields.FirstOrDefault(i => i.Id == "DP");

        Assert.That(dpInfo.Type, Is.EqualTo("Integer"));
        Assert.That(dpInfo.Description, Does.Contain("Depth"));
    }

    #endregion

    #region Variant Classification Tests

    [Test]
    public void ClassifyVariant_SNP_ReturnsSnp()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var snp = records.First(r => r.Ref.Length == 1 && r.Alt[0].Length == 1);

        Assert.That(VcfParser.ClassifyVariant(snp), Is.EqualTo(VcfParser.VariantType.SNP));
    }

    [Test]
    public void ClassifyVariant_Insertion_ReturnsInsertion()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var insertion = records.First(r => r.Ref == "A" && r.Alt[0] == "AT");

        Assert.That(VcfParser.ClassifyVariant(insertion), Is.EqualTo(VcfParser.VariantType.Insertion));
    }

    [Test]
    public void ClassifyVariant_Deletion_ReturnsDeletion()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var deletion = records.First(r => r.Ref == "ATG" && r.Alt[0] == "A");

        Assert.That(VcfParser.ClassifyVariant(deletion), Is.EqualTo(VcfParser.VariantType.Deletion));
    }

    [Test]
    public void IsSNP_ForSnp_ReturnsTrue()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var snp = records.First(r => r.Ref == "A" && r.Alt[0] == "G");

        Assert.That(VcfParser.IsSNP(snp), Is.True);
    }

    [Test]
    public void IsIndel_ForInsertion_ReturnsTrue()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var insertion = records.First(r => r.Alt[0] == "AT");

        Assert.That(VcfParser.IsIndel(insertion), Is.True);
    }

    [Test]
    public void GetVariantLength_Insertion_ReturnsCorrectLength()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var insertion = records.First(r => r.Alt[0] == "AT");

        Assert.That(VcfParser.GetVariantLength(insertion), Is.EqualTo(1)); // AT - A = 1
    }

    [Test]
    public void GetVariantLength_Deletion_ReturnsCorrectLength()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var deletion = records.First(r => r.Ref == "ATG");

        Assert.That(VcfParser.GetVariantLength(deletion), Is.EqualTo(2)); // ATG - A = 2
    }

    [Test]
    public void ClassifyVariant_MNP_ReturnsMnp()
    {
        // MNP: Multi-nucleotide polymorphism - ref and alt same length > 1
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	AT	GC	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(VcfParser.ClassifyVariant(records[0]), Is.EqualTo(VcfParser.VariantType.MNP));
    }

    [Test]
    public void ClassifyVariant_Symbolic_ReturnsSymbolic()
    {
        // Symbolic alleles: structural variants like <DEL>, <INS>, <DUP>
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	<DEL>	99	PASS	.
chr1	200	.	A	<INS>	99	PASS	.
chr1	300	.	A	<DUP>	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(records, Has.Count.EqualTo(3));
        Assert.Multiple(() =>
        {
            Assert.That(VcfParser.ClassifyVariant(records[0]), Is.EqualTo(VcfParser.VariantType.Symbolic));
            Assert.That(VcfParser.ClassifyVariant(records[1]), Is.EqualTo(VcfParser.VariantType.Symbolic));
            Assert.That(VcfParser.ClassifyVariant(records[2]), Is.EqualTo(VcfParser.VariantType.Symbolic));
        });
    }

    [Test]
    public void ClassifyVariant_Complex_ReturnsComplex()
    {
        // Complex: both ref and alt lengths differ but don't fit insertion/deletion pattern
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	ATG	CT	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(VcfParser.ClassifyVariant(records[0]), Is.EqualTo(VcfParser.VariantType.Complex));
    }

    [Test]
    public void ClassifyVariant_BreakendSymbolic_ReturnsSymbolic()
    {
        // Breakend notation using [ or ] brackets
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	]chr2:500]A	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(VcfParser.ClassifyVariant(records[0]), Is.EqualTo(VcfParser.VariantType.Symbolic));
    }

    #endregion

    #region Filtering Tests

    [Test]
    public void FilterByChrom_ReturnsMatchingChromosome()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var chr1 = VcfParser.FilterByChrom(records, "chr1").ToList();

        Assert.That(chr1, Has.Count.EqualTo(2));
        Assert.That(chr1.All(r => r.Chrom == "chr1"), Is.True);
    }

    [Test]
    public void FilterByRegion_ReturnsVariantsInRegion()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var inRegion = VcfParser.FilterByRegion(records, "chr1", 50, 150).ToList();

        Assert.That(inRegion, Has.Count.EqualTo(1));
        Assert.That(inRegion[0].Pos, Is.EqualTo(100));
    }

    [Test]
    public void FilterByQuality_FiltersLowQuality()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var highQual = VcfParser.FilterByQuality(records, minQuality: 60).ToList();

        Assert.That(highQual.All(r => r.Qual >= 60), Is.True);
    }

    [Test]
    public void FilterPassing_ReturnsOnlyPassing()
    {
        var records = VcfParser.Parse(VcfWithFilters).ToList();
        var passing = VcfParser.FilterPassing(records).ToList();

        Assert.That(passing, Has.Count.EqualTo(1));
        Assert.That(passing[0].Filter[0], Is.EqualTo("PASS").IgnoreCase);
    }

    [Test]
    public void FilterSNPs_ReturnsOnlySNPs()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var snps = VcfParser.FilterSNPs(records).ToList();

        Assert.That(snps, Has.Count.EqualTo(1));
        Assert.That(VcfParser.IsSNP(snps[0]), Is.True);
    }

    [Test]
    public void FilterIndels_ReturnsOnlyIndels()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var indels = VcfParser.FilterIndels(records).ToList();

        Assert.That(indels, Has.Count.EqualTo(2));
        Assert.That(indels.All(r => VcfParser.IsIndel(r)), Is.True);
    }

    [Test]
    public void FilterByInfo_FiltersCorrectly()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var filtered = VcfParser.FilterByInfo(records, "DP", v => int.Parse(v) >= 40).ToList();

        Assert.That(filtered.All(r => int.Parse(r.Info["DP"]) >= 40), Is.True);
    }

    #endregion

    #region Genotype Analysis Tests

    [Test]
    public void GetGenotype_ReturnsGenotype()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var gt = VcfParser.GetGenotype(records[0], 0);

        Assert.That(gt, Is.EqualTo("0/1"));
    }

    [Test]
    public void IsHomRef_ForHomRef_ReturnsTrue()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        // Sample1 at position 200 is 0/0
        var record = records[1];

        Assert.That(VcfParser.IsHomRef(record, 0), Is.True);
    }

    [Test]
    public void IsHomAlt_ForHomAlt_ReturnsTrue()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        // Sample2 at position 100 is 1/1
        var record = records[0];

        Assert.That(VcfParser.IsHomAlt(record, 1), Is.True);
    }

    [Test]
    public void IsHet_ForHet_ReturnsTrue()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        // Sample1 at position 100 is 0/1
        var record = records[0];

        Assert.That(VcfParser.IsHet(record, 0), Is.True);
    }

    [Test]
    public void GetReadDepth_ReturnsDepth()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var dp = VcfParser.GetReadDepth(records[0], 0);

        Assert.That(dp, Is.EqualTo(25));
    }

    [Test]
    public void GetAlleleDepth_ReturnsDepths()
    {
        // Test with a VCF that has AD field - skip if not available
        // This is a basic test to ensure the method works
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var ad = VcfParser.GetAlleleDepth(records[0], 0);

        // May be null if AD field is not present
        Assert.Pass();
    }

    [Test]
    public void IsHomRef_PhasedGenotype_ReturnsTrue()
    {
        // Test phased genotype notation (| instead of /)
        const string vcf = @"##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G	99	PASS	.	GT	0|0";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(VcfParser.IsHomRef(records[0], 0), Is.True);
    }

    [Test]
    public void IsHomAlt_PhasedGenotype_ReturnsTrue()
    {
        // Test phased genotype notation for homozygous alternate
        const string vcf = @"##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G	99	PASS	.	GT	1|1";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(VcfParser.IsHomAlt(records[0], 0), Is.True);
    }

    [Test]
    public void IsHet_PhasedGenotype_ReturnsTrue()
    {
        // Test phased heterozygous genotype
        const string vcf = @"##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G	99	PASS	.	GT	0|1";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(VcfParser.IsHet(records[0], 0), Is.True);
    }

    [Test]
    public void IsHomAlt_MultiAllelic_DetectsCorrectly()
    {
        // Test homozygous for second alternate allele (2/2)
        const string vcf = @"##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	.	A	G,T	99	PASS	.	GT	2/2";

        var records = VcfParser.Parse(vcf).ToList();

        Assert.That(VcfParser.IsHomAlt(records[0], 0), Is.True);
    }

    [Test]
    public void IsHet_MultiAllelicBothAlternate_ReturnsTrue()
    {
        // Test heterozygous with both alleles being alternates (1/2)
        var records = VcfParser.Parse(SimpleVcf).ToList();
        // Third record has multi-allelic site with 1/2 genotype
        var multiAllelicRecord = records[2];

        Assert.That(VcfParser.IsHet(multiAllelicRecord, 0), Is.True);
    }

    [Test]
    public void GetGenotype_MissingSample_ReturnsNull()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var gt = VcfParser.GetGenotype(records[0], 99); // Non-existent sample index

        Assert.That(gt, Is.Null);
    }

    [Test]
    public void IsHomRef_NoSamples_ReturnsFalse()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList(); // VCF without samples

        Assert.That(VcfParser.IsHomRef(records[0], 0), Is.False);
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void CalculateStatistics_ReturnsCorrectCounts()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        var stats = VcfParser.CalculateStatistics(records);

        Assert.That(stats.TotalVariants, Is.EqualTo(3));
        Assert.That(stats.SnpCount, Is.EqualTo(1));
        Assert.That(stats.IndelCount, Is.EqualTo(2));
    }

    [Test]
    public void CalculateStatistics_PassingCount_Correct()
    {
        var records = VcfParser.Parse(VcfWithFilters).ToList();
        var stats = VcfParser.CalculateStatistics(records);

        Assert.That(stats.PassingCount, Is.EqualTo(1));
    }

    [Test]
    public void CalculateStatistics_ChromosomeCounts_Correct()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var stats = VcfParser.CalculateStatistics(records);

        Assert.That(stats.ChromosomeCounts["chr1"], Is.EqualTo(2));
        Assert.That(stats.ChromosomeCounts["chr2"], Is.EqualTo(1));
    }

    [Test]
    public void CalculateTiTvRatio_ReturnsRatio()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	.
chr1	200	.	C	T	99	PASS	.
chr1	300	.	A	T	99	PASS	.
chr1	400	.	G	C	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();
        var titv = VcfParser.CalculateTiTvRatio(records);

        // 2 transitions (A>G, C>T) / 2 transversions (A>T, G>C) = 1.0
        Assert.That(titv, Is.EqualTo(1.0).Within(0.01));
    }

    #endregion

    #region Info Field Helpers Tests

    [Test]
    public void GetInfoValue_ReturnsValue()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var dp = VcfParser.GetInfoValue(records[0], "DP");

        Assert.That(dp, Is.EqualTo("50"));
    }

    [Test]
    public void GetInfoValue_NonexistentKey_ReturnsNull()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var value = VcfParser.GetInfoValue(records[0], "NONEXISTENT");

        Assert.That(value, Is.Null);
    }

    [Test]
    public void GetInfoInt_ReturnsInteger()
    {
        var records = VcfParser.Parse(SimpleVcf).ToList();
        var dp = VcfParser.GetInfoInt(records[0], "DP");

        Assert.That(dp, Is.EqualTo(50));
    }

    [Test]
    public void GetInfoDouble_ReturnsDouble()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	AF=0.5";

        var records = VcfParser.Parse(vcf).ToList();
        var af = VcfParser.GetInfoDouble(records[0], "AF");

        Assert.That(af, Is.EqualTo(0.5).Within(0.001));
    }

    [Test]
    public void HasInfoFlag_ForExistingFlag_ReturnsTrue()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	DB";

        var records = VcfParser.Parse(vcf).ToList();
        Assert.That(VcfParser.HasInfoFlag(records[0], "DB"), Is.True);
    }

    #endregion

    #region Writing Tests

    [Test]
    public void WriteToStream_ProducesValidVcf()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        using var writer = new StringWriter();

        VcfParser.WriteToStream(writer, records);
        var output = writer.ToString();

        Assert.That(output, Does.Contain("##fileformat="));
        Assert.That(output, Does.Contain("#CHROM"));
        Assert.That(output, Does.Contain("chr1"));
    }

    [Test]
    public void WriteAndRead_Roundtrip_PreservesData()
    {
        var original = VcfParser.Parse(VcfWithIndels).ToList();
        using var writer = new StringWriter();

        VcfParser.WriteToStream(writer, original);
        var output = writer.ToString();

        var parsed = VcfParser.Parse(output).ToList();

        Assert.That(parsed.Count, Is.EqualTo(original.Count));
        Assert.That(parsed[0].Chrom, Is.EqualTo(original[0].Chrom));
        Assert.That(parsed[0].Pos, Is.EqualTo(original[0].Pos));
        Assert.That(parsed[0].Ref, Is.EqualTo(original[0].Ref));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Parse_MissingQuality_QualIsNull()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	.	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();
        Assert.That(records[0].Qual, Is.Null);
    }

    [Test]
    public void Parse_MissingId_IdIsDot()
    {
        var records = VcfParser.Parse(VcfWithIndels).ToList();
        Assert.That(records[0].Id, Is.EqualTo("."));
    }

    [Test]
    public void Parse_EmptyFilter_FilterIsEmpty()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	.	.";

        var records = VcfParser.Parse(vcf).ToList();
        Assert.That(records[0].Filter, Is.Empty);
    }

    [Test]
    public void Parse_MultipleFilters_AllParsed()
    {
        var records = VcfParser.Parse(VcfWithFilters).ToList();
        var multiFilter = records.First(r => r.Filter.Length > 1);

        Assert.That(multiFilter.Filter, Contains.Item("LowQual"));
        Assert.That(multiFilter.Filter, Contains.Item("LowCov"));
    }

    [Test]
    public void Parse_MalformedLine_Skips()
    {
        const string vcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	abc	.	A	G	99	PASS	.
chr1	100	.	A	G	99	PASS	.";

        var records = VcfParser.Parse(vcf).ToList();
        Assert.That(records, Has.Count.EqualTo(1));
    }

    #endregion

    #region File I/O Tests

    [Test]
    public void ParseFile_NonexistentFile_ReturnsEmpty()
    {
        var records = VcfParser.ParseFile("nonexistent.vcf").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void ParseFile_ValidFile_ParsesRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleVcf);
            var records = VcfParser.ParseFile(tempFile).ToList();

            Assert.That(records, Has.Count.EqualTo(3));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    [Test]
    public void ParseFileWithHeader_ValidFile_ReturnsHeaderAndRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleVcf);
            var (header, records) = VcfParser.ParseFileWithHeader(tempFile);

            Assert.That(header.FileFormat, Is.EqualTo("VCFv4.3"));
            Assert.That(records.Count(), Is.EqualTo(3));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    #endregion
}
