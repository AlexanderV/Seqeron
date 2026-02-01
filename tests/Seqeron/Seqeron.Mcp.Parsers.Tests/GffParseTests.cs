using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class GffParseTests
{
    private const string TestGff = @"##gff-version 3
chr1	HAVANA	gene	100	500	.	+	.	ID=gene1;Name=TestGene
chr1	HAVANA	mRNA	100	500	.	+	.	ID=mrna1;Parent=gene1
chr1	HAVANA	exon	100	200	.	+	.	ID=exon1;Parent=mrna1
chr1	HAVANA	exon	300	500	.	+	.	ID=exon2;Parent=mrna1
chr2	HAVANA	gene	1000	2000	50.5	-	.	ID=gene2;Name=OtherGene";

    [Test]
    public void GffParse_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GffParse(TestGff));
        Assert.Throws<ArgumentException>(() => ParsersTools.GffParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.GffParse(null!));
    }

    [Test]
    public void GffParse_Binding_ParsesRecords()
    {
        var result = ParsersTools.GffParse(TestGff);

        Assert.That(result.Count, Is.EqualTo(5));
        Assert.That(result.Records[0].Seqid, Is.EqualTo("chr1"));
        Assert.That(result.Records[0].Source, Is.EqualTo("HAVANA"));
        Assert.That(result.Records[0].Type, Is.EqualTo("gene"));
        Assert.That(result.Records[0].Start, Is.EqualTo(100));
        Assert.That(result.Records[0].End, Is.EqualTo(500));
        Assert.That(result.Records[0].Strand, Is.EqualTo("+"));
        Assert.That(result.Records[0].GeneName, Is.EqualTo("TestGene"));
    }

    [Test]
    public void GffParse_Binding_ParsesAttributes()
    {
        var result = ParsersTools.GffParse(TestGff);

        Assert.That(result.Records[0].Attributes, Contains.Key("ID"));
        Assert.That(result.Records[0].Attributes["ID"], Is.EqualTo("gene1"));
        Assert.That(result.Records[0].Attributes, Contains.Key("Name"));
        Assert.That(result.Records[0].Attributes["Name"], Is.EqualTo("TestGene"));
    }

    [Test]
    public void GffParse_Binding_ParsesScore()
    {
        var result = ParsersTools.GffParse(TestGff);

        Assert.That(result.Records[0].Score, Is.Null); // "." score
        Assert.That(result.Records[4].Score, Is.EqualTo(50.5)); // numeric score
    }

    [Test]
    public void GffParse_Binding_CalculatesLength()
    {
        var result = ParsersTools.GffParse(TestGff);

        Assert.That(result.Records[0].Length, Is.EqualTo(401)); // 500 - 100 + 1
    }

    [Test]
    public void GffParse_Binding_RespectsFormat()
    {
        // Auto-detect GFF3 format
        var result = ParsersTools.GffParse(TestGff, "auto");
        Assert.That(result.Count, Is.EqualTo(5));

        // Explicit GFF3 format
        result = ParsersTools.GffParse(TestGff, "gff3");
        Assert.That(result.Count, Is.EqualTo(5));
    }
}

[TestFixture]
public class GffStatisticsTests
{
    private const string TestGff = @"##gff-version 3
chr1	HAVANA	gene	100	500	.	+	.	ID=gene1
chr1	HAVANA	exon	100	200	.	+	.	ID=exon1
chr1	HAVANA	exon	300	500	.	+	.	ID=exon2
chr2	HAVANA	gene	1000	2000	.	-	.	ID=gene2
chr2	HAVANA	CDS	1100	1900	.	-	0	ID=cds1";

    [Test]
    public void GffStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GffStatistics(TestGff));
        Assert.Throws<ArgumentException>(() => ParsersTools.GffStatistics(""));
    }

    [Test]
    public void GffStatistics_Binding_CalculatesStats()
    {
        var result = ParsersTools.GffStatistics(TestGff);

        Assert.That(result.TotalFeatures, Is.EqualTo(5));
        Assert.That(result.GeneCount, Is.EqualTo(2));
        Assert.That(result.ExonCount, Is.EqualTo(2));
        Assert.That(result.SequenceIds, Contains.Item("chr1"));
        Assert.That(result.SequenceIds, Contains.Item("chr2"));
        Assert.That(result.Sources, Contains.Item("HAVANA"));
    }

    [Test]
    public void GffStatistics_Binding_CountsFeatureTypes()
    {
        var result = ParsersTools.GffStatistics(TestGff);

        Assert.That(result.FeatureTypeCounts["gene"], Is.EqualTo(2));
        Assert.That(result.FeatureTypeCounts["exon"], Is.EqualTo(2));
        Assert.That(result.FeatureTypeCounts["CDS"], Is.EqualTo(1));
    }
}

[TestFixture]
public class GffFilterTests
{
    private const string TestGff = @"##gff-version 3
chr1	HAVANA	gene	100	500	.	+	.	ID=gene1
chr1	HAVANA	exon	100	200	.	+	.	ID=exon1
chr1	HAVANA	exon	300	500	.	+	.	ID=exon2
chr2	HAVANA	gene	1000	2000	.	-	.	ID=gene2
chr2	HAVANA	CDS	1100	1900	.	-	0	ID=cds1";

    [Test]
    public void GffFilter_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.GffFilter(TestGff));
        Assert.Throws<ArgumentException>(() => ParsersTools.GffFilter(""));
    }

    [Test]
    public void GffFilter_Binding_FiltersByType()
    {
        var result = ParsersTools.GffFilter(TestGff, featureType: "gene");

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.Records.All(r => r.Type == "gene"), Is.True);
    }

    [Test]
    public void GffFilter_Binding_FiltersBySeqid()
    {
        var result = ParsersTools.GffFilter(TestGff, seqid: "chr1");

        Assert.That(result.PassedCount, Is.EqualTo(3));
        Assert.That(result.Records.All(r => r.Seqid == "chr1"), Is.True);
    }

    [Test]
    public void GffFilter_Binding_FiltersByRegion()
    {
        var result = ParsersTools.GffFilter(TestGff, seqid: "chr1", regionStart: 150, regionEnd: 400);

        // Should match exon1 (100-200), exon2 (300-500), and gene1 (100-500) if they overlap
        Assert.That(result.PassedCount, Is.EqualTo(3));
    }

    [Test]
    public void GffFilter_Binding_CombinesFilters()
    {
        var result = ParsersTools.GffFilter(TestGff, featureType: "exon", seqid: "chr1");

        Assert.That(result.PassedCount, Is.EqualTo(2)); // Only chr1 exons
        Assert.That(result.Records.All(r => r.Type == "exon" && r.Seqid == "chr1"), Is.True);
    }

    [Test]
    public void GffFilter_Binding_CalculatesPercentage()
    {
        var result = ParsersTools.GffFilter(TestGff, featureType: "gene");

        Assert.That(result.TotalCount, Is.EqualTo(5));
        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.PassedPercentage, Is.EqualTo(40.0));
    }
}
