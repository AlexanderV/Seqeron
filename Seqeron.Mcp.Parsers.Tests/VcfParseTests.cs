using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class VcfParseTests
{
    private const string TestVcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs1	A	G	30	PASS	DP=10
chr1	200	rs2	AT	A	25	PASS	DP=15
chr2	300	rs3	C	T	20	.	DP=8";

    [Test]
    public void VcfParse_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfParse(TestVcf));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfParse(null!));
    }

    [Test]
    public void VcfParse_Binding_ParsesRecords()
    {
        var result = ParsersTools.VcfParse(TestVcf);

        Assert.That(result.Count, Is.EqualTo(3));
        Assert.That(result.Records[0].Chrom, Is.EqualTo("chr1"));
        Assert.That(result.Records[0].Pos, Is.EqualTo(100));
        Assert.That(result.Records[0].Ref, Is.EqualTo("A"));
        Assert.That(result.Records[0].Alt, Contains.Item("G"));
        Assert.That(result.Records[0].VariantType, Is.EqualTo("SNP"));
    }

    [Test]
    public void VcfParse_Binding_ClassifiesVariantTypes()
    {
        var result = ParsersTools.VcfParse(TestVcf);

        Assert.That(result.Records[0].VariantType, Is.EqualTo("SNP")); // A>G
        Assert.That(result.Records[1].VariantType, Is.EqualTo("Deletion")); // AT>A
        Assert.That(result.Records[2].VariantType, Is.EqualTo("SNP")); // C>T
    }
}

[TestFixture]
public class VcfStatisticsTests
{
    private const string TestVcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr1	200	.	AT	A	25	PASS	.
chr2	300	.	C	T	20	.	.";

    [Test]
    public void VcfStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfStatistics(TestVcf));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfStatistics(""));
    }

    [Test]
    public void VcfStatistics_Binding_CalculatesStats()
    {
        var result = ParsersTools.VcfStatistics(TestVcf);

        Assert.That(result.TotalVariants, Is.EqualTo(3));
        Assert.That(result.SnpCount, Is.EqualTo(2));
        Assert.That(result.IndelCount, Is.EqualTo(1));
        Assert.That(result.PassingCount, Is.EqualTo(3)); // "." filter is also considered passing
        Assert.That(result.ChromosomeCounts["chr1"], Is.EqualTo(2));
        Assert.That(result.ChromosomeCounts["chr2"], Is.EqualTo(1));
    }
}

[TestFixture]
public class VcfFilterTests
{
    private const string TestVcf = @"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	30	PASS	.
chr1	200	.	AT	A	25	PASS	.
chr2	300	.	C	T	15	LowQual	.";

    [Test]
    public void VcfFilter_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfFilter(TestVcf));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfFilter(""));
    }

    [Test]
    public void VcfFilter_Binding_FiltersByType()
    {
        var result = ParsersTools.VcfFilter(TestVcf, variantType: "snp");

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.Records.All(r => r.VariantType == "SNP"), Is.True);
    }

    [Test]
    public void VcfFilter_Binding_FiltersByChrom()
    {
        var result = ParsersTools.VcfFilter(TestVcf, chrom: "chr1");

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.Records.All(r => r.Chrom == "chr1"), Is.True);
    }

    [Test]
    public void VcfFilter_Binding_FiltersByQuality()
    {
        var result = ParsersTools.VcfFilter(TestVcf, minQuality: 20);

        Assert.That(result.PassedCount, Is.EqualTo(2)); // 30 and 25
    }

    [Test]
    public void VcfFilter_Binding_FiltersPassOnly()
    {
        var result = ParsersTools.VcfFilter(TestVcf, passOnly: true);

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.Records.All(r => r.Filter.Contains("PASS")), Is.True);
    }

    [Test]
    public void VcfFilter_Binding_CombinesFilters()
    {
        var result = ParsersTools.VcfFilter(TestVcf, variantType: "snp", chrom: "chr1");

        Assert.That(result.PassedCount, Is.EqualTo(1)); // Only chr1 SNP
    }
}
